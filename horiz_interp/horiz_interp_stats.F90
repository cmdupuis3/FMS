module horizontal_interpolator_stats_mod

    use mpp_mod, only : mpp_send, mpp_recv, mpp_sync_self
    use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_npes
    use mpp_mod, only : COMM_TAG_1, COMM_TAG_2

    implicit none
    private

    public :: stats, stats_conservative

contains

    !> These statistics are for bilinear interpolation and spherical regrid.
    subroutine stats ( dat, low, high, avg, miss, missing_value, mask )
        real,    intent(in)  :: dat(:,:)
        real,    intent(out) :: low, high, avg
        integer, intent(out) :: miss
        real,    intent(in), optional :: missing_value
        real,    intent(in), optional :: mask(:,:)

        real :: dsum, npts, buffer_real(3)
        integer :: pe, root_pe, npes, p, buffer_int(2)

        pe      = mpp_pe()
        root_pe = mpp_root_pe()
        npes    = mpp_npes()

        dsum = 0.0
        miss = 0

        if (present(missing_value)) then
            miss = count(dat(:,:) == missing_value)
            low  = minval(dat(:,:), dat(:,:) /= missing_value)
            high = maxval(dat(:,:), dat(:,:) /= missing_value)
            dsum =    sum(dat(:,:), dat(:,:) /= missing_value)
        else if(present(mask)) then
            miss = count(mask(:,:) <= 0.5)
            low  = minval(dat(:,:), mask=mask(:,:) > 0.5)
            high = maxval(dat(:,:), mask=mask(:,:) > 0.5)
            dsum =    sum(dat(:,:), mask=mask(:,:) > 0.5)
        else
            miss = 0
            low  = minval(dat(:,:))
            high = maxval(dat(:,:))
            dsum =    sum(dat(:,:))
        endif
        avg = 0.0

        npts = size(dat(:,:)) - miss
        if(pe == root_pe) then
            do p = 1, npes - 1  ! root_pe receive data from other pe
                ! Force use of "scalar", integer pointer mpp interface
                call mpp_recv(buffer_real(1),glen=3, from_pe=p+root_pe, tag=COMM_TAG_1)
                dsum = dsum + buffer_real(1)
                low  = min(low,  buffer_real(2))
                high = max(high, buffer_real(3))
                call mpp_recv(buffer_int(1), glen=2, from_pe=p+root_pe, tag=COMM_TAG_2)
                miss = miss + buffer_int(1)
                npts = npts + buffer_int(2)
            enddo
            if(npts == 0.) then
                print*, 'Warning: no points is valid'
            else
                avg = dsum/real(npts)
            endif
        else   ! other pe send data to the root_pe.
            buffer_real(1) = dsum
            buffer_real(2) = low
            buffer_real(3) = high
            ! Force use of "scalar", integer pointer mpp interface
            call mpp_send(buffer_real(1), plen=3, to_pe=root_pe, tag=COMM_TAG_1)
            buffer_int(1) = miss
            buffer_int(2) = npts
            call mpp_send(buffer_int(1),  plen=2, to_pe=root_pe, tag=COMM_TAG_2)
        endif

        call mpp_sync_self()

        return

    end subroutine stats

    !> These statistics are for conservative schemes
    subroutine stats_conservative ( dat, area, asum, dsum, wsum, low, high, miss, mask )
        real,    intent(in)  :: dat(:,:), area(:,:)
        real,    intent(out) :: asum, dsum, wsum, low, high
        integer, intent(out) :: miss
        real,    intent(in), optional :: mask(:,:)

        integer :: pe, root_pe, npes, p, buffer_int(1)
        real    :: buffer_real(5)

        pe = mpp_pe()
        root_pe = mpp_root_pe()
        npes = mpp_npes()

        ! sum data, data*area; and find min,max on each pe.

        if (present(mask)) then
            asum = sum(area(:,:))
            dsum = sum(area(:,:)*dat(:,:)*mask(:,:))
            wsum = sum(area(:,:)*mask(:,:))
            miss = count(mask(:,:) <= 0.5)
            low  = minval(dat(:,:),mask=mask(:,:) > 0.5)
            high = maxval(dat(:,:),mask=mask(:,:) > 0.5)
        else
            asum = sum(area(:,:))
            dsum = sum(area(:,:)*dat(:,:))
            wsum = sum(area(:,:))
            miss = 0
            low  = minval(dat(:,:))
            high = maxval(dat(:,:))
        endif

        ! other pe send local min, max, avg to the root pe and
        ! root pe receive these information

        if(pe == root_pe) then
            do p = 1, npes - 1
                ! Force use of "scalar", integer pointer mpp interface
                call mpp_recv(buffer_real(1),glen=5,from_pe=root_pe+p, tag=COMM_TAG_1)
                asum = asum + buffer_real(1)
                dsum = dsum + buffer_real(2)
                wsum = wsum + buffer_real(3)
                low  = min(low, buffer_real(4))
                high = max(high, buffer_real(5))
                call mpp_recv(buffer_int(1),glen=1,from_pe=root_pe+p, tag=COMM_TAG_2)
                miss = miss + buffer_int(1)
            enddo
        else
            buffer_real(1) = asum
            buffer_real(2) = dsum
            buffer_real(3) = wsum
            buffer_real(4) = low
            buffer_real(5) = high
            ! Force use of "scalar", integer pointer mpp interface
            call mpp_send(buffer_real(1),plen=5,to_pe=root_pe, tag=COMM_TAG_1)
            buffer_int(1) = miss
            call mpp_send(buffer_int(1),plen=1,to_pe=root_pe, tag=COMM_TAG_2)
        endif

        call mpp_sync_self()

    end subroutine stats_conservative


end module horizontal_interpolator_stats_mod
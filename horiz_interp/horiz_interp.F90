!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!> \brief Performs spatial interpolation between grids.
!!
!! This module can interpolate data from any logically rectangular grid
!! to any logically rectangular grid. Four interpolation schems are used here:
!! conservative, bilinear, bicubic and inverse of square distance weighted.
!! The four interpolation schemes are implemented seperately in
!! horiz_interp_conserver_mod, horiz_interp_blinear_mod, horiz_interp_bicubic_mod
!! and horiz_interp_spherical_mod. bicubic interpolation requires the source grid
!! is regular lon/lat grid. User can choose the interpolation method in the
!! public interface horiz_interp_new through optional argument interp_method,
!! with acceptable value "conservative", "bilinear", "bicubic" and "spherical".
!! The default value is "conservative". There is an optional mask field for
!! missing input data. An optional output mask field may be used in conjunction with
!! the input mask to show where output data exists.
module horiz_interp_mod

    use fms_mod,                    only: write_version_number, fms_error_handler
    use fms_mod,                    only: check_nml_error
    use mpp_mod,                    only: mpp_error, FATAL, stdout, stdlog, mpp_min
    use mpp_mod,                    only: input_nml_file, WARNING, mpp_pe, mpp_root_pe
    use constants_mod,              only: pi
    use horizontal_interpolator_types_mod,      only: conservative1HZI_t, conservative2HZI_t, bilinearHZI_t, sphericalHZI_t, bicubicHZI_t
    use horizontal_interpolator_types_mod,      only: conservativeHZI_t, baseHZI_t, assignment(=)
    use horizontal_interpolator_conservative_mod,  only: horiz_interp_conserve_init, horiz_interp_conserve
    use horizontal_interpolator_conservative_mod,  only: horiz_interp_conservative_new_1dx1d, horiz_interp_conservative_new_1dx2d
    use horizontal_interpolator_conservative_mod,  only: horiz_interp_conservative_new_2dx1d, horiz_interp_conservative_new_2dx2d
    use horizontal_interpolator_conservative_mod,  only: hzi_delete_conservative1, hzi_delete_conservative2
    use horizontal_interpolator_bilinear_mod,  only: horiz_interp_bilinear_init, horiz_interp_bilinear
    use horizontal_interpolator_bilinear_mod,  only: horiz_interp_bilinear_new, hzi_delete_bilinear
    use horizontal_interpolator_bicubic_mod,   only: horiz_interp_bicubic_init, horiz_interp_bicubic
    use horizontal_interpolator_bicubic_mod,   only: horiz_interp_bicubic_new, hzi_delete_bicubic
    use horizontal_interpolator_spherical_mod, only: horiz_interp_spherical_init, horiz_interp_spherical
    use horizontal_interpolator_spherical_mod, only: horiz_interp_spherical_new, hzi_delete_spherical

    implicit none
    private

    public :: horiz_interp_type, horiz_interp, horiz_interp_new, horiz_interp_del, horiz_interp_init, horiz_interp_end, assignment(=)

! <INTERFACE NAME="horiz_interp_new">
!   <OVERVIEW>
!      Allocates space and initializes a derived-type variable
!      that contains pre-computed interpolation indices and weights.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Allocates space and initializes a derived-type variable
!      that contains pre-computed interpolation indices and weights
!      for improved performance of multiple interpolations between
!      the same grids. This routine does not need to be called if you
!      are doing a single grid-to-grid interpolation.
!   </DESCRIPTION>
!   <IN NAME="lon_in" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians">
!      Longitude (in radians) for source data grid. You can pass 1-D lon_in to
!      represent the geographical longitude of regular lon/lat grid, or just
!      pass geographical longitude(lon_in is 2-D). The grid location may be
!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!   </IN>
!   <IN NAME="lat_in" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians">
!      Latitude (in radians) for source data grid. You can pass 1-D lat_in to
!      represent the geographical latitude of regular lon/lat grid, or just
!      pass geographical latitude(lat_in is 2-D). The grid location may be
!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!   </IN>
!   <IN NAME="lon_out" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians" >
!      Longitude (in radians) for destination data grid. You can pass 1-D lon_out to
!      represent the geographical longitude of regular lon/lat grid, or just
!      pass geographical longitude(lon_out is 2-D). The grid location may be
!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!   </IN>
!   <IN NAME="lat_out" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians" >
!      Latitude (in radians) for destination data grid. You can pass 1-D lat_out to
!      represent the geographical latitude of regular lon/lat grid, or just
!      pass geographical latitude(lat_out is 2-D). The grid location may be
!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!   </IN>
!   <IN NAME="verbose" TYPE="integer">
!      Integer flag that controls the amount of printed output.
!      verbose = 0, no output; = 1, min,max,means; = 2, still more
!   </IN>
!   <IN NAME="interp_method" TYPE="character(len=*)" >
!      interpolation method, = "conservative", using conservation scheme,
!      = "bilinear", using bilinear interpolation, = "spherical",using spherical regrid.
!      = "bicubic", using bicubic interpolation. The default value is "convervative".
!   </IN>
!   <IN NAME = "src_modulo" >
!      Indicate the source data grid is cyclic or not.
!   </IN>
!   <IN NAME = "grid_at_center" >
!      Indicate the data is on the center of grid box or the edge of grid box.
!      When true, the data is on the center of grid box. default vaule is false.
!      This option is only available when interp_method = "bilinear" or "bicubic".
!   </IN>
!   <OUT NAME="Interp" >
!      A derived-type variable containing indices and weights used for subsequent
!      interpolations. To reinitialize this variable for a different grid-to-grid
!      interpolation you must first use the "horiz_interp_del" interface.
!   </OUT>

    interface hzi_new_conservative
        module procedure hzi_new_conservative_1dx1d ! Source grid is 1d, destination grid is 1d
        module procedure hzi_new_conservative_1dx2d ! Source grid is 1d, destination grid is 2d
        module procedure hzi_new_conservative_2dx2d ! Source grid is 2d, destination grid is 2d
        module procedure hzi_new_conservative_2dx1d ! Source grid is 2d, destination grid is 1d
    end interface

    interface hzi_new_bilinear
        module procedure hzi_new_bilinear_1dx1d
        module procedure hzi_new_bilinear_1dx2d
        module procedure hzi_new_bilinear_2dx2d
        module procedure hzi_new_bilinear_2dx1d
    end interface

    interface hzi_new_bilinear_centered
        module procedure hzi_new_bilinear_1dx1d_centered
        module procedure hzi_new_bilinear_1dx2d_centered
    end interface

    interface hzi_new_bicubic
        module procedure hzi_new_bicubic_1dx1d
        module procedure hzi_new_bicubic_1dx2d
    end interface

    interface hzi_new_bicubic_centered
        module procedure hzi_new_bicubic_1dx1d_centered
        module procedure hzi_new_bicubic_1dx2d_centered
    end interface

    interface hzi_new_spherical
        module procedure hzi_new_spherical_1dx1d
        module procedure hzi_new_spherical_1dx2d
        module procedure hzi_new_spherical_2dx2d
        module procedure hzi_new_spherical_2dx1d
    end interface

! </INTERFACE>

! <INTERFACE NAME="horiz_interp">
!
!   <OVERVIEW>
!     Subroutine for performing the horizontal interpolation between two grids.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Subroutine for performing the horizontal interpolation between
!     two grids. There are two forms of this interface.
!     Form A requires first calling horiz_interp_new, while Form B
!     requires no initialization.
!   </DESCRIPTION>

!   <IN NAME="Interp" >
!     Derived-type variable containing interpolation indices and weights.
!     Returned by a previous call to horiz_interp_new.
!   </IN>
!   <IN NAME="data_in">
!      Input data on source grid.
!   </IN>
!   <IN NAME="verbose">
!      flag for the amount of print output.
!               verbose = 0, no output; = 1, min,max,means; = 2, still more
!   </IN>
!   <IN NAME="mask_in">
!      Input mask, must be the same size as the input data. The real value of
!      mask_in must be in the range (0.,1.). Set mask_in=0.0 for data points
!      that should not be used or have missing data. It is Not needed for
!      spherical regrid.
!   </IN>
!   <IN NAME="missing_value" >
!      Use the missing_value to indicate missing data.
!   </IN>
!   <IN NAME="missing_permit">
!      numbers of points allowed to miss for the bilinear interpolation. The value
!      should be between 0 and 3.
!   </IN>
!   <IN NAME="lon_in, lat_in" >
!      longitude and latitude (in radians) of source grid. More explanation can
!      be found in the documentation of horiz_interp_new.
!   </IN>
!   <IN NAME="lon_out, lat_out" >
!      longitude and latitude (in radians) of destination grid. More explanation can
!      be found in the documentation of horiz_interp_new.
!   </IN>
!   <OUT NAME="data_out">
!      Output data on destination grid.
!   </OUT>
!   <OUT NAME="mask_out">
!      Output mask that specifies whether data was computed.
!   </OUT>

!   <ERROR MSG="size of input array incorrect" STATUS="FATAL">
!      The input data array does not match the size of the input grid edges
!      specified. If you are using the initialization interface make sure you
!      have the correct grid size.
!   </ERROR>
!   <ERROR MSG="size of output array incorrect" STATUS="FATAL">
!      The output data array does not match the size of the input grid
!      edges specified. If you are using the initialization interface make
!      sure you have the correct grid size.
!   </ERROR>

    interface horiz_interp
        module procedure horiz_interp_base_2d
        module procedure horiz_interp_base_3d
    end interface
! </INTERFACE>



    interface hzi_delete
        module procedure hzi_delete_conservative1
        module procedure hzi_delete_conservative2
        module procedure hzi_delete_bilinear
        module procedure hzi_delete_spherical
        module procedure hzi_delete_bicubic
    end interface


    ! Include variable "version" to be written to log file.
    #include<file_version.h>
    logical :: module_is_initialized = .FALSE.

contains

    subroutine horiz_interp_init
        integer :: unit, ierr, io

        if(module_is_initialized) return
        call write_version_number("HORIZ_INTERP_MOD", version)

        read (input_nml_file, horiz_interp_nml, iostat=io)
        ierr = check_nml_error(io,'horiz_interp_nml')
        if (mpp_pe() == mpp_root_pe() ) then
            unit = stdlog()
            write (unit, nml=horiz_interp_nml)
        endif

        call horiz_interp_conserve_init
        call horiz_interp_bilinear_init
        call horiz_interp_bicubic_init
        call horiz_interp_spherical_init

        module_is_initialized = .true.

    end subroutine horiz_interp_init

    function hzi_new_conservative_1dx1d (lon_in, lat_in, lon_out, lat_out, verbose) result(Interp)

        type(conservative1HZI_t) :: Interp
        real, intent(in), dimension(:)          :: lon_in , lat_in
        real, intent(in), dimension(:)          :: lon_out, lat_out
        integer, intent(in), optional           :: verbose

        call horiz_interp_init
        call horiz_interp_conservative_new_1dx1d (Interp, lon_in, lat_in, lon_out, lat_out, verbose)

    end function hzi_new_conservative_1dx1d

    function hzi_new_conservative_1dx2d (lon_in, lat_in, lon_out, lat_out, verbose, mask_in, mask_out) result(Interp)

        type(conservative2HZI_t) :: Interp
        real, intent(in),  dimension(:)               :: lon_in , lat_in
        real, intent(in),  dimension(:,:)             :: lon_out, lat_out
        integer, intent(in),                 optional :: verbose
        real, intent(in), dimension(:,:),    optional :: mask_in
        real, intent(out),dimension(:,:),    optional :: mask_out

        call horiz_interp_init

        !--- check to see if the source grid is regular lat-lon grid or not.
        if(is_lat_lon(lon_out, lat_out)) then
            if(present(mask_in)) then
                if ( ANY(mask_in < -.0001) .or. ANY(mask_in > 1.0001)  ) call mpp_error(FATAL, &
                    'horiz_interp::hzi_new_conservative_1dx2d: input mask not between 0,1')
                allocate(Interp%mask_in(size(mask_in,1), size(mask_in,2)) )
                Interp%mask_in = mask_in
            end if
            call horiz_interp_conservative_new_1dx2d (Interp, lon_in, lat_in, lon_out(:,1), lat_out(1,:), verbose=verbose)
        else
            call horiz_interp_conservative_new_1dx2d (Interp, lon_in, lat_in, lon_out, lat_out, verbose=verbose, mask_in=mask_in, mask_out=mask_out)
        end if

    end function hzi_new_conservative_1dx2d

    function hzi_new_conservative_2dx2d (lon_in, lat_in, lon_out, lat_out, verbose, mask_in, mask_out, is_latlon_in, is_latlon_out) result(Interp)

        type(conservative2HZI_t) :: Interp
        real, intent(in), dimension(:,:) :: lon_in , lat_in
        real, intent(in), dimension(:,:) :: lon_out, lat_out
        integer, intent(in),              optional :: verbose
        real, intent(in), dimension(:,:), optional :: mask_in
        real, intent(out),dimension(:,:), optional :: mask_out
        logical, intent(in),              optional :: is_latlon_in, is_latlon_out

        logical :: src_is_latlon, dst_is_latlon

        call horiz_interp_init

        if(PRESENT(is_latlon_in)) then
            src_is_latlon = is_latlon_in
        else
            src_is_latlon = is_lat_lon(lon_in, lat_in)
        end if
        if(PRESENT(is_latlon_out)) then
            dst_is_latlon = is_latlon_out
        else
            dst_is_latlon = is_lat_lon(lon_out, lat_out)
        end if
        if(src_is_latlon .AND. dst_is_latlon) then
            if(present(mask_in)) then
                if ( ANY(mask_in < -.0001) .or. ANY(mask_in > 1.0001)  ) call mpp_error(FATAL, &
                    'horiz_interp_conserve_new_2d(horiz_interp_conserve_mod): input mask not between 0,1')
                allocate(Interp%mask_in(size(mask_in,1), size(mask_in,2)) )
                Interp%mask_in = mask_in
            end if
            call horiz_interp_conservative_new_2dx2d (Interp, lon_in(:,1), lat_in(1,:), lon_out(:,1), lat_out(1,:), verbose=verbose)
        else if(src_is_latlon .AND. .not. dst_is_latlon) then
            call horiz_interp_conservative_new_2dx2d (Interp, lon_in(:,1), lat_in(1,:), lon_out,      lat_out,      verbose=verbose, mask_in=mask_in, mask_out=mask_out)
        else if(.not. src_is_latlon .AND. dst_is_latlon) then
            call horiz_interp_conservative_new_2dx2d (Interp, lon_in,      lat_in,      lon_out(:,1), lat_out(1,:), verbose=verbose, mask_in=mask_in, mask_out=mask_out)
        else
            call horiz_interp_conservative_new_2dx2d (Interp, lon_in,      lat_in,      lon_out,      lat_out,      verbose=verbose, mask_in=mask_in, mask_out=mask_out)
        end if

    end function hzi_new_conservative_2dx2d

    function hzi_new_conservative_2dx1d (lon_in, lat_in, lon_out, lat_out, verbose, mask_in, mask_out, is_latlon_in) result(Interp)

        type(conservative2HZI_t) :: Interp
        real, intent(in),  dimension(:,:)          :: lon_in , lat_in
        real, intent(in),  dimension(:)            :: lon_out, lat_out
        integer, intent(in),              optional :: verbose
        real, intent(in), dimension(:,:), optional :: mask_in
        real, intent(out),dimension(:,:), optional :: mask_out
        logical, intent(in),              optional :: is_latlon_in

        integer                           :: i, j, nlon_out, nlat_out
        real, dimension(:,:), allocatable :: lon_dst, lat_dst
        logical                           :: src_is_latlon

        call horiz_interp_init

        nlon_out = size(lon_out(:)); nlat_out = size(lat_out(:))
        allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
        do i = 1, nlon_out
            lon_dst(i,:) = lon_out(i)
        enddo
        do j = 1, nlat_out
            lat_dst(:,j) = lat_out(j)
        enddo

        if(PRESENT(is_latlon_in)) then
            src_is_latlon = is_latlon_in
        else
            src_is_latlon = is_lat_lon(lon_in, lat_in)
        end if

        if(src_is_latlon) then
            if(present(mask_in)) then
                if ( ANY(mask_in < -.0001) .or. ANY(mask_in > 1.0001)  ) call mpp_error(FATAL, &
                    'horiz_interp_conserve_new_1d_dst(horiz_interp_conserve_mod): input mask not between 0,1')
                allocate(Interp%mask_in(size(mask_in,1), size(mask_in,2)) )
                Interp%mask_in = mask_in
            end if
            call horiz_interp_conservative_new_2dx1d (Interp, lon_in(:,1), lat_in(1,:), lon_out, lat_out, verbose=verbose)
        else
            call horiz_interp_conservative_new_2dx1d (Interp, lon_in,      lat_in,      lon_out, lat_out, verbose=verbose, mask_in=mask_in, mask_out=mask_out)
        end if

        deallocate(lon_dst,lat_dst)

    end function hzi_new_conservative_2dx1d

    function hzi_new_bilinear_1dx1d (lon_in, lat_in, lon_out, lat_out, verbose, src_modulo) result(Interp)

        type(bilinearHZI_t) :: Interp
        real, intent(in), dimension(:) :: lon_in , lat_in
        real, intent(in), dimension(:) :: lon_out, lat_out
        integer, intent(in), optional  :: verbose
        logical, intent(in), optional  :: src_modulo

        real, dimension(:,:), allocatable :: lon_dst, lat_dst
        real, dimension(:),   allocatable :: lon_src_1d, lat_src_1d
        integer                           :: i, j, nlon_in, nlat_in, nlon_out, nlat_out

        call horiz_interp_init

        nlon_in  = size(lon_in(:))-1
        nlat_in  = size(lat_in(:))-1
        nlon_out = size(lon_out(:))-1
        nlat_out = size(lat_out(:))-1
        allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
        allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
        do i = 1, nlon_in
            lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5
        enddo
        do j = 1, nlat_in
            lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5
        enddo
        do i = 1, nlon_out
            lon_dst(i,:) = (lon_out(i) + lon_out(i+1)) * 0.5
        enddo
        do j = 1, nlat_out
            lat_dst(:,j) = (lat_out(j) + lat_out(j+1)) * 0.5
        enddo
        call horiz_interp_bilinear_new (Interp, lon_src_1d, lat_src_1d, lon_dst, lat_dst, verbose, src_modulo)
        deallocate(lon_src_1d, lat_src_1d, lon_dst, lat_dst)

    end function hzi_new_bilinear_1dx1d

    function hzi_new_bilinear_1dx2d (Interp, lon_in, lat_in, lon_out, lat_out, verbose, src_modulo) result(Interp)

        type(bilinearHZI_t) :: Interp
        real, intent(in),  dimension(:)    :: lon_in , lat_in
        real, intent(in),  dimension(:,:)  :: lon_out, lat_out
        integer, intent(in), optional      :: verbose
        logical, intent(in), optional      :: src_modulo

        real, dimension(:),   allocatable :: lon_src_1d, lat_src_1d
        integer                           :: i, j, nlon_in, nlat_in

        call horiz_interp_init

        nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
        allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
        do i = 1, nlon_in
            lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5
        enddo
        do j = 1, nlat_in
            lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5
        enddo
        call horiz_interp_bilinear_new (Interp, lon_src_1d, lat_src_1d, lon_out, lat_out, verbose, src_modulo)
        deallocate(lon_src_1d,lat_src_1d)

    end function hzi_new_bilinear_1dx2d

    function hzi_new_bilinear_2dx2d (lon_in, lat_in, lon_out, lat_out, verbose, src_modulo) result(Interp)

        type(bilinearHZI_t) :: Interp
        real, intent(in), dimension(:,:) :: lon_in , lat_in
        real, intent(in), dimension(:,:) :: lon_out, lat_out
        integer, intent(in), optional    :: verbose
        logical, intent(in), optional    :: src_modulo

        call horiz_interp_init
        call horiz_interp_bilinear_new (Interp, lon_in, lat_in, lon_out, lat_out, verbose, src_modulo)

    end function hzi_new_bilinear_2dx2d

    function hzi_new_bilinear_2dx1d (lon_in, lat_in, lon_out, lat_out, verbose, src_modulo) result(Interp)

        type(bilinearHZI_t) :: Interp
        real, intent(in), dimension(:,:) :: lon_in , lat_in
        real, intent(in), dimension(:)   :: lon_out, lat_out
        integer, intent(in),    optional :: verbose
        logical, intent(in),    optional :: src_modulo

        integer                           :: i, j, nlon_out, nlat_out
        real, dimension(:,:), allocatable :: lon_dst, lat_dst

        call horiz_interp_init

        nlon_out = size(lon_out(:)); nlat_out = size(lat_out(:))
        allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
        do i = 1, nlon_out
            lon_dst(i,:) = lon_out(i)
        enddo
        do j = 1, nlat_out
            lat_dst(:,j) = lat_out(j)
        enddo

        call horiz_interp_bilinear_new (Interp, lon_in, lat_in, lon_dst, lat_dst, verbose, src_modulo)
        deallocate(lon_dst,lat_dst)

    end function hzi_new_bilinear_2dx1d

    function hzi_new_bilinear_1dx1d_centered (lon_in, lat_in, lon_out, lat_out, verbose, src_modulo) result(Interp)

        type(bilinearHZI_t) :: Interp
        real, intent(in), dimension(:) :: lon_in , lat_in
        real, intent(in), dimension(:) :: lon_out, lat_out
        integer, intent(in), optional  :: verbose
        logical, intent(in), optional  :: src_modulo

        real, dimension(:,:), allocatable :: lon_dst, lat_dst
        integer                           :: i, j, nlon_out, nlat_out

        call horiz_interp_init
        nlon_out = size(lon_out(:))
        nlat_out = size(lat_out(:))
        allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
        do i = 1, nlon_out
            lon_dst(i,:) = lon_out(i)
        enddo
        do j = 1, nlat_out
            lat_dst(:,j) = lat_out(j)
        enddo

        call horiz_interp_bilinear_new (Interp, lon_in, lat_in, lon_dst, lat_dst, verbose, src_modulo)
        deallocate(lon_dst, lat_dst)

    end function hzi_new_bilinear_1dx1d_centered

    function hzi_new_bilinear_1dx2d_centered (Interp, lon_in, lat_in, lon_out, lat_out, verbose, src_modulo) result(Interp)

        type(bilinearHZI_t) :: Interp
        real, intent(in),  dimension(:)    :: lon_in , lat_in
        real, intent(in),  dimension(:,:)  :: lon_out, lat_out
        integer, intent(in), optional      :: verbose
        logical, intent(in), optional      :: src_modulo

        call horiz_interp_init
        call horiz_interp_bilinear_new (Interp, lon_in, lat_in, lon_out, lat_out, verbose, src_modulo)

    end function hzi_new_bilinear_1dx2d_centered

    function hzi_new_bicubic_1dx1d (lon_in, lat_in, lon_out, lat_out, verbose) result(Interp)

        type(bicubicHZI_t) :: Interp
        real, intent(in), dimension(:) :: lon_in , lat_in
        real, intent(in), dimension(:) :: lon_out, lat_out
        integer, intent(in), optional  :: verbose

        real, dimension(:), allocatable :: lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d
        integer                         :: i, j, nlon_in, nlat_in, nlon_out, nlat_out

        call horiz_interp_init

        !No need to expand to 2d, horiz_interp_bicubic_new does 1d-1d
        nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
        nlon_out = size(lon_out(:))-1; nlat_out = size(lat_out(:))-1
        allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
        allocate(lon_dst_1d(nlon_out), lat_dst_1d(nlat_out))
        do i = 1, nlon_in
            lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5
        enddo
        do j = 1, nlat_in
            lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5
        enddo
        do i = 1, nlon_out
            lon_dst_1d(i) = (lon_out(i) + lon_out(i+1)) * 0.5
        enddo
        do j = 1, nlat_out
            lat_dst_1d(j) = (lat_out(j) + lat_out(j+1)) * 0.5
        enddo
        call horiz_interp_bicubic_new (Interp, lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d, verbose)
        deallocate(lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d)

    end function hzi_new_bicubic_1dx1d

    function hzi_new_bicubic_1dx2d (Interp, lon_in, lat_in, lon_out, lat_out, verbose) result(Interp)

        type(bicubicHZI_t) :: Interp
        real, intent(in), dimension(:)   :: lon_in , lat_in
        real, intent(in), dimension(:,:) :: lon_out, lat_out
        integer, intent(in), optional    :: verbose

        real, dimension(:),   allocatable :: lon_src_1d, lat_src_1d
        integer                           :: i, j, nlon_in, nlat_in

        call horiz_interp_init

        nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
        allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
        do i = 1, nlon_in
            lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5
        enddo
        do j = 1, nlat_in
            lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5
        enddo
        call horiz_interp_bicubic_new (Interp, lon_src_1d, lat_src_1d, lon_out, lat_out, verbose)
        deallocate(lon_src_1d,lat_src_1d)

    end function hzi_new_bicubic_1dx2d

    function hzi_new_bicubic_1dx1d_centered (lon_in, lat_in, lon_out, lat_out, verbose) result(Interp)

        type(bicubicHZI_t) :: Interp
        real, intent(in), dimension(:) :: lon_in , lat_in
        real, intent(in), dimension(:) :: lon_out, lat_out
        integer, intent(in), optional  :: verbose

        call horiz_interp_init

        !No need to expand to 2d, horiz_interp_bicubic_new does 1d-1d
        call horiz_interp_bicubic_new (Interp, lon_in, lat_in, lon_out, lat_out, verbose)

    end function hzi_new_bicubic_1dx1d_centered

    function hzi_new_bicubic_1dx2d_centered (Interp, lon_in, lat_in, lon_out, lat_out, verbose) result(Interp)

        type(bicubicHZI_t) :: Interp
        real, intent(in), dimension(:)   :: lon_in , lat_in
        real, intent(in), dimension(:,:) :: lon_out, lat_out
        integer, intent(in), optional    :: verbose

        call horiz_interp_init
        call horiz_interp_bicubic_new (Interp, lon_in, lat_in, lon_out, lat_out, verbose)

    end function hzi_new_bicubic_1dx2d_centered

    function hzi_new_spherical_1dx1d (lon_in, lat_in, lon_out, lat_out, num_nbrs, max_dist, src_modulo) result(Interp)

        type(sphericalHZI_t) :: Interp
        real, intent(in), dimension(:) :: lon_in , lat_in
        real, intent(in), dimension(:) :: lon_out, lat_out
        integer, intent(in),  optional :: num_nbrs
        real,    intent(in),  optional :: max_dist
        logical, intent(in),  optional :: src_modulo

        real, dimension(:,:), allocatable :: lon_src, lat_src, lon_dst, lat_dst
        integer                           :: i, j, nlon_in, nlat_in, nlon_out, nlat_out

        call horiz_interp_init

        nlon_in  = size(lon_in(:));   nlat_in  = size(lat_in(:))
        nlon_out  = size(lon_out(:)); nlat_out = size(lat_out(:))
        allocate(lon_src(nlon_in,nlat_in), lat_src(nlon_in,nlat_in))
        allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
        do i = 1, nlon_in
            lon_src(i,:) = lon_in(i)
        enddo
        do j = 1, nlat_in
            lat_src(:,j) = lat_in(j)
        enddo
        do i = 1, nlon_out
            lon_dst(i,:) = lon_out(i)
        enddo
        do j = 1, nlat_out
            lat_dst(:,j) = lat_out(j)
        enddo
        call horiz_interp_spherical_new (Interp, lon_src, lat_src, lon_dst, lat_dst, num_nbrs, max_dist, src_modulo)
        deallocate(lon_src, lat_src, lon_dst, lat_dst)

    end function hzi_new_spherical_1dx1d

    function hzi_new_spherical_1dx2d (Interp, lon_in, lat_in, lon_out, lat_out, num_nbrs, max_dist, src_modulo) result(Interp)

        type(sphericalHZI_t) :: Interp
        real, intent(in), dimension(:)   :: lon_in , lat_in
        real, intent(in), dimension(:,:) :: lon_out, lat_out
        integer, intent(in), optional    :: num_nbrs  ! minimum number of neighbors
        real,    intent(in), optional    :: max_dist
        logical, intent(in), optional    :: src_modulo

        real, dimension(:,:), allocatable :: lon_src, lat_src
        integer                           :: i, j, nlon_in, nlat_in

        call horiz_interp_init

        nlon_in  = size(lon_in(:));  nlat_in  = size(lat_in(:))
        allocate(lon_src(nlon_in,nlat_in), lat_src(nlon_in,nlat_in))
        do i = 1, nlon_in
            lon_src(i,:) = lon_in(i)
        enddo
        do j = 1, nlat_in
            lat_src(:,j) = lat_in(j)
        enddo
        call horiz_interp_spherical_new (Interp, lon_src, lat_src, lon_out, lat_out, num_nbrs, max_dist, src_modulo)
        deallocate(lon_src, lat_src)

    end function hzi_new_spherical_1dx2d

    function hzi_new_spherical_2dx2d (lon_in, lat_in, lon_out, lat_out, num_nbrs, max_dist, src_modulo) result(Interp)

        type(sphericalHZI_t) :: Interp
        real, intent(in), dimension(:,:) :: lon_in , lat_in
        real, intent(in), dimension(:,:) :: lon_out, lat_out
        integer, intent(in), optional    :: num_nbrs
        real,    intent(in), optional    :: max_dist
        logical, intent(in), optional    :: src_modulo

        call horiz_interp_init
        call horiz_interp_spherical_new (Interp, lon_in, lat_in, lon_out, lat_out, num_nbrs, max_dist, src_modulo)

    end function hzi_new_spherical_2dx2d

    function hzi_new_spherical_2dx1d (Interp, lon_in, lat_in, lon_out, lat_out, num_nbrs, max_dist, src_modulo) result(Interp)

        type(sphericalHZI_t) :: Interp
        real, intent(in), dimension(:,:) :: lon_in , lat_in
        real, intent(in), dimension(:)   :: lon_out, lat_out
        integer, intent(in),    optional :: num_nbrs
        real,    intent(in),    optional :: max_dist
        logical, intent(in),    optional :: src_modulo

        integer                           :: i, j, nlon_out, nlat_out
        real, dimension(:,:), allocatable :: lon_dst, lat_dst

        call horiz_interp_init

        nlon_out = size(lon_out(:)); nlat_out = size(lat_out(:))
        allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
        do i = 1, nlon_out
            lon_dst(i,:) = lon_out(i)
        enddo
        do j = 1, nlat_out
            lat_dst(:,j) = lat_out(j)
        enddo

        call horiz_interp_spherical_new ( Interp, lon_in, lat_in, lon_dst, lat_dst, num_nbrs, max_dist, src_modulo)
        deallocate(lon_dst,lat_dst)

    end function hzi_new_spherical_2dx1d



    function horiz_interp_base_2d (Interp, data_in, verbose, mask_in, mask_out, missing_value, missing_permit, new_missing_handle) result(data_out)

        type (baseHZI_t), intent(in)                 :: Interp
        real, intent(in),   dimension(:,:)           :: data_in
        real,               dimension(:,:)           :: data_out
        integer, intent(in),                optional :: verbose
        real, intent(in),   dimension(:,:), optional :: mask_in
        real, intent(out),  dimension(:,:), optional :: mask_out
        real, intent(in),                   optional :: missing_value
        integer, intent(in),                optional :: missing_permit
        logical, intent(in),                optional :: new_missing_handle

        select type(Interp)
        !class is (baseHZI_t)
        type is (conservativeHZI_t)
            call horiz_interp_conservative (Interp, data_in, data_out, verbose, mask_in, mask_out)
        type is (bilinearHZI_t)
            call horiz_interp_bilinear     (Interp, data_in, data_out, verbose, mask_in, mask_out, missing_value, missing_permit, new_missing_handle )
        type is (bicubicHZI_t)
            call horiz_interp_bicubic      (Interp, data_in, data_out, verbose, mask_in, mask_out, missing_value, missing_permit)
        type is (sphericalHZI_t)
            call horiz_interp_spherical    (Interp, data_in, data_out, verbose, mask_in, mask_out, missing_value )
        end select

    end function horiz_interp_base_2d

    function horiz_interp_base_3d (Interp, data_in, verbose, mask_in, mask_out, missing_value, missing_permit) result(data_out)

        type (baseHZI_t), intent(in)                   :: Interp
        real, intent(in),   dimension(:,:,:)           :: data_in
        real,               dimension(:,:,:)           :: data_out
        integer, intent(in),                  optional :: verbose
        real, intent(in),   dimension(:,:,:), optional :: mask_in
        real, intent(out),  dimension(:,:,:), optional :: mask_out
        real, intent(in),                     optional :: missing_value
        integer, intent(in),                  optional :: missing_permit

        integer :: n

        do n = 1, size(data_in,3)
            if (present(mask_in))
                if(present(mask_out)) then
                    call horiz_interp_base_2d ( Interp, data_in(:,:,n), data_out(:,:,n), verbose,
                        mask_in = mask_in(:,:,n), mask_out = mask_out(:,:,n), missing_value = missing_value, missing_permit = missing_permit )
                else
                    call horiz_interp_base_2d ( Interp, data_in(:,:,n), data_out(:,:,n), verbose,
                        mask_in = mask_in(:,:,n),                             missing_value = missing_value, missing_permit = missing_permit )
                endif
            else
                if(present(mask_out)) then
                    call horiz_interp_base_2d ( Interp, data_in(:,:,n), data_out(:,:,n), verbose,
                                                  mask_out = mask_out(:,:,n), missing_value = missing_value, missing_permit = missing_permit )
                else
                    call horiz_interp_base_2d ( Interp, data_in(:,:,n), data_out(:,:,n), verbose,
                                                                              missing_value = missing_value, missing_permit = missing_permit )
                endif
            endif
        enddo

    end function horiz_interp_base_3d


    subroutine horiz_interp_end
        return
    end subroutine horiz_interp_end

    function is_lat_lon(lon, lat)

        real, dimension(:,:), intent(in) :: lon, lat
        logical                          :: is_lat_lon
        integer                          :: i, j, nlon, nlat, num

        is_lat_lon = .true.
        nlon = size(lon,1)
        nlat = size(lon,2)
        do j = 1, nlat
        do i = 2, nlon
            if(lat(i,j) .NE. lat(1,j)) then
                is_lat_lon = .false.
                exit
            end if
        end do
        end do

        if(is_lat_lon) then
        do i = 1, nlon
            do j = 2, nlat
                if(lon(i,j) .NE. lon(i,1)) then
                    is_lat_lon = .false.
                    exit
                end if
            end do
        end do
        end if

        num = 0
        if(is_lat_lon) num = 1
        call mpp_min(num)
        if(num == 1) then
            is_lat_lon = .true.
        else
            is_lat_lon = .false.
        end if

        return

    end function is_lat_lon

end module horiz_interp_mod

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
module horizontal_interpolator_types_mod

    use mpp_mod, only : mpp_error, FATAL

    implicit none
    private

    ! Oh, how I wish Fortran enums worked like in normal languages...
    enum, bind(C)
        enumerator CONSERVATIVE_1
        enumerator CONSERVATIVE_2
        enumerator BILINEAR
        enumerator SPHERICAL
        enumerator BICUBIC
    end enum

    type, abstract :: baseHZI_t
        real,    dimension(:,:),   pointer :: area_src => NULL()        !< area of the source grid
        real,    dimension(:,:),   pointer :: area_dst => NULL()        !< area of the destination grid
        real,    dimension(:,:,:), pointer :: src_dist => NULL()        !< distance between destination grid and neighbor source grid.
        logical, dimension(:,:),   pointer :: found_neighbors => NULL() !< indicate whether destination grid has some source grid around it.
        real                               :: max_src_dist
        integer, dimension(:,:),   pointer :: num_found => NULL()
        integer                            :: nlon_src                  !< size of source grid
        integer                            :: nlat_src                  !< size of source grid
        integer                            :: nlon_dst                  !< size of destination grid
        integer                            :: nlat_dst                  !< size of destination grid
        real,    dimension(:,:),   pointer :: rat_x => NULL()           !< the ratio of coordinates of the dest grid (x_dest - x_src_r) / (x_src_l - x_src_r)
        real,    dimension(:,:),   pointer :: rat_y => NULL()           !< the ratio of coordinates of the dest grid (y_dest - y_src_r) / (y_src_l - y_src_r)
        real,    dimension(:),     pointer :: lon_in => NULL()          !< the coordinates of the source grid
        real,    dimension(:),     pointer :: lat_in => NULL()          !< the coordinates of the source grid
    end type baseHZI_t

    type, extends(baseHZI_t) :: conservative1HZI_t
        real,    dimension(:,:), pointer :: faci => NULL()   !< weights
        real,    dimension(:,:), pointer :: facj => NULL()   !< weights
        integer, dimension(:,:), pointer :: ilon => NULL()   !< indices
        integer, dimension(:,:), pointer :: jlat => NULL()   !< indices
    end type conservative1HZI_t

    type, extends(conservative1HZI_t) :: conservative2HZI_t
        integer                          :: nxgrid                  !< number of exchange grid between src and dst grid.
        integer, dimension(:),   pointer :: i_src => NULL()         !< indices in source grid.
        integer, dimension(:),   pointer :: j_src => NULL()         !< indices in source grid.
        integer, dimension(:),   pointer :: i_dst => NULL()         !< indices in destination grid.
        integer, dimension(:),   pointer :: j_dst => NULL()         !< indices in destination grid.
        real,    dimension(:),   pointer :: area_frac_dst => NULL() !< area fraction in destination grid.
        real,    dimension(:,:), pointer :: mask_in => NULL()
    end type conservative2HZI_t

    type, extends(baseHZI_t) :: bilinearHZI_t
        real,  dimension(:,:,:), pointer :: wti => NULL() !< weights
        real,  dimension(:,:,:), pointer :: wtj => NULL() !< weights
    end type bilinearHZI_t

    type, extends(baseHZI_t) :: sphericalHZI_t
    end type sphericalHZI_t

    type, extends(baseHZI_t) :: bicubicHZI_t
        real,    dimension(:,:,:), pointer :: wti  => NULL() !< derivative weights
        integer, dimension(:,:,:), pointer :: ilon => NULL() !< indices
        integer, dimension(:,:,:), pointer :: jlat => NULL() !< indices
    end type bicubicHZI_t

    interface assignment(=)
        module procedure hzi_type_eq_conservative1
        module procedure hzi_type_eq_conservative2
        module procedure hzi_type_eq_bilinear
        module procedure hzi_type_eq_spherical
        module procedure hzi_type_eq_bicubic
    end interface

    public :: CONSERVATIVE_1,     CONSERVATIVE_2,     BILINEAR,      SPHERICAL,      BICUBIC
    public :: conservative1HZI_t, conservative2HZI_t, bilinearHZI_t, sphericalHZI_t, bicubicHZI_t
    public :: assignment(=)

contains

    subroutine hzi_type_eq(hzi_in, hzi_out)
        type(baseHZI_t), intent(inout) :: hzi_out
        type(baseHZI_t), intent(in)    :: hzi_in

        hzi_out%area_src        => hzi_in%area_src
        hzi_out%area_dst        => hzi_in%area_dst
        hzi_out%src_dist        => hzi_in%src_dist
        hzi_out%found_neighbors => hzi_in%found_neighbors
        hzi_out%max_src_dist    =  hzi_in%max_src_dist
        hzi_out%num_found       => hzi_in%num_found
        hzi_out%nlon_src        =  hzi_in%nlon_src
        hzi_out%nlat_src        =  hzi_in%nlat_src
        hzi_out%nlon_dst        =  hzi_in%nlon_dst
        hzi_out%nlat_dst        =  hzi_in%nlat_dst
        hzi_out%rat_x           => hzi_in%rat_x
        hzi_out%rat_y           => hzi_in%rat_y
        hzi_out%lon_in          => hzi_in%lon_in
        hzi_out%lat_in          => hzi_in%lat_in
    end subroutine hzi_type_eq

    subroutine hzi_type_eq_conservative1(hzi_in, hzi_out)
        type(conservative1HZI_t), intent(inout) :: hzi_out
        type(conservative1HZI_t), intent(in)    :: hzi_in

        call hzi_type_eq(hzi_in, hzi_out)

        hzi_out%faci            => hzi_in%faci
        hzi_out%facj            => hzi_in%facj
        hzi_out%ilon            => hzi_in%ilon
        hzi_out%jlat            => hzi_in%jlat
    end subroutine hzi_type_eq_conservative1

    subroutine hzi_type_eq_conservative2(hzi_in, hzi_out)
        type(conservative2HZI_t), intent(inout) :: hzi_out
        type(conservative2HZI_t), intent(in)    :: hzi_in

        call hzi_type_eq(hzi_in, hzi_out)

        hzi_out%nxgrid          =  hzi_in%nxgrid
        hzi_out%i_src           => hzi_in%i_src
        hzi_out%j_src           => hzi_in%j_src
        hzi_out%i_dst           => hzi_in%i_dst
        hzi_out%j_dst           => hzi_in%j_dst
        hzi_out%area_frac_dst   => hzi_in%area_frac_dst
        !TODO: what about mask? is this a bug?
    end subroutine hzi_type_eq_conservative2

    subroutine hzi_type_eq_bilinear(hzi_in, hzi_out)
        type(bilinearHZI_t), intent(inout) :: hzi_out
        type(bilinearHZI_t), intent(in)    :: hzi_in

        call hzi_type_eq(hzi_in, hzi_out)

        hzi_out%wti             => hzi_in%wti
        hzi_out%wtj             => hzi_in%wtj
    end subroutine hzi_type_eq_bilinear

    subroutine hzi_type_eq_spherical(hzi_in, hzi_out)
        type(sphericalHZI_t), intent(inout) :: hzi_out
        type(sphericalHZI_t), intent(in)    :: hzi_in

        call hzi_type_eq(hzi_in, hzi_out)

    end subroutine hzi_type_eq_spherical

    subroutine hzi_type_eq_bicubic(hzi_in, hzi_out)
        type(bicubicHZI_t), intent(inout) :: hzi_out
        type(bicubicHZI_t), intent(in)    :: hzi_in

        call hzi_type_eq(hzi_in, hzi_out)

        hzi_out%wti             => hzi_in%wti
        hzi_out%ilon            => hzi_in%ilon
        hzi_out%jlat            => hzi_in%jlat
    end subroutine hzi_type_eq_bicubic


end module horizontal_interpolator_types_mod
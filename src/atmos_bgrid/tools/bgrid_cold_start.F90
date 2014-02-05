
module bgrid_cold_start_mod

use mpp_mod, only: input_nml_file 
use    bgrid_horiz_mod, only: horiz_grid_type, get_horiz_grid_bound, TGRID
use     bgrid_halo_mod, only: update_halo, TEMP
use gaussian_topog_mod, only: gaussian_topog_init
use            fms_mod, only: file_exist, open_namelist_file, open_restart_file, &
                              open_ieee32_file, check_nml_error, close_file,     &
                              mpp_pe, mpp_root_pe, set_domain, read_data,        &
                              stdlog, error_mesg, FATAL, write_version_number,   &
                              field_size, uppercase
use      constants_mod, only: GRAV, RDGAS

implicit none
private

public cold_start_resol, cold_start

!-----------------------------------------------------------------------
!----- namelist ------
!
!  nlon = number of grid points along the longitude axis (1st dimension)
!  nlat = number of grid points along the latitude axis (2nd dimension)
!  nlev = number of vertical levels (equally spaced in sigma)
!  pref = initial surface pressure in pascals
!  tref = initial temperature in deg kelvin
!  equal_vert_spacing = specifies whether equal vertical spacing in sigma
!                       should be used (TRUE) or unequal spacing using 
!                       the Smagorinski formula (FALSE).
!  topog_option = how to compute topography
!                  possible values are: FLAT, FILE, GAUSS
!  topog_file   = name of topography restart file in topog_option=FILE
!
!  NOTE: nlon and nlat are for the global compute grid
!

integer :: nlon = 0
integer :: nlat = 0
integer :: nlev = 0
   real :: pref = 1000.e2
   real :: tref = 255.
logical :: equal_vert_spacing = .true.
character(len=8)   :: topog_option = 'FLAT'
character(len=128) :: topog_file = ' '

namelist /bgrid_cold_start_nml/ nlon, nlat, nlev, pref, tref, &
                                equal_vert_spacing, topog_option, &
                                topog_file

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: bgrid_cold_start.F90,v 19.0 2012/01/06 19:54:36 fms Exp $'
character(len=128) :: tag = '$Name: tikal $'

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine cold_start_resol ( nx, ny, nz )

      integer, intent(out) :: nx, ny, nz

      integer :: unit, ierr, io, logunit

!-------------- read namelist --------------

   if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=bgrid_cold_start_nml, iostat=io)
ierr = check_nml_error(io,'bgrid_cold_start_nml')
#else
      unit = open_namelist_file ( )
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=bgrid_cold_start_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'bgrid_cold_start_nml')
      enddo
 10   call close_file (unit)
#endif
   endif

!-------- write version and namelist to log file --------

   call write_version_number (version,tag)
   logunit = stdlog()
   if (mpp_pe() == mpp_root_pe()) write (logunit, nml=bgrid_cold_start_nml)

!------- must specify a resolution -----

   if (nlon == 0 .or. nlat == 0 .or. nlev == 0)  &
   call error_mesg ('bgrid_cold_start_mod', 'resolution not specified', FATAL)

!------- otherwise, return resolution to calling program ------

   nx  = nlon
   ny  = nlat
   nz  = nlev

end subroutine cold_start_resol

!#######################################################################

subroutine cold_start ( Hgrid, eta, peta, fis, res, ps, pssl, u, v, t )

   type(horiz_grid_type), intent(inout) :: Hgrid
   real, intent(out), dimension(:)      :: eta, peta
   real, intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:)     :: fis, res, &
                                                              ps, pssl
   real, intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:)   :: u, v, t
  !real, intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:,:) :: r
   integer :: m
   real    :: rt

!   very simple initial condition
!   -----------------------------
!   uniform pressure (adjusted for topog)
!   isothermal
!   no wind


!--- no hybrid option ---
   peta = 0.0

   eta = compute_sigma (size(eta(:))-1)

   call compute_topog (Hgrid, fis, res)

   rt = RDGAS * tref

     ps   = pref*exp(-fis/rt)
     pssl = ps                ! no eta option
     u    = 0.0
     v    = 0.0
     t    = tref

end subroutine cold_start

!#######################################################################

 subroutine compute_topog (Hgrid, fis, res)

   type(horiz_grid_type), intent(inout) :: Hgrid
   real, intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:) :: fis, res

   integer :: unit, ires, jres, siz(4), nc
   logical :: eof_found
!-----------------------------------------------------------------------

!---------------------------------------------
! OPTION 1: No mountains, flat topography
!---------------------------------------------
  if ( trim(uppercase(topog_option)) == 'FLAT' ) then
       fis = 0.
       res = 1.

!---------------------------------------------
! OPTION 2: gaussian-shaped mountains
!---------------------------------------------
  else if ( trim(uppercase(topog_option(1:5))) == 'GAUSS' ) then

       call gaussian_topog_init ( Hgrid%Tmp%alm(:,Hgrid%Tmp%js), &
                                  Hgrid%Tmp%aph(Hgrid%Tmp%is,:), &
                                  fis(:,:) )
       fis = fis * GRAV                       ! convert topog from meters to geop ht
       where (fis(:,:) < 0.0) fis(:,:) = 0.0  ! only in case of roundoff error
       res = 1.

!---------------------------------------------
! OPTION 3: read topography from restart file
!---------------------------------------------
  else if ( trim(uppercase(topog_option)) == 'FILE' ) then

   if ( file_exist(trim(topog_file)) ) then
      ! netcdf or native?
        nc = len_trim(topog_file)
      !------------------------------------------
      ! read topography from netcdf restart file

        if (topog_file(nc-2:nc) == '.nc') then

          ! read fis
            call field_size (topog_file(1:nc), 'fis', siz )
            if ( siz(1) /= Hgrid%nlon .or. siz(2) /= Hgrid%nlat )     &
                 call error_mesg ('bgrid_cold_start_mod', 'incorrect resolution &
                              &or no fis field in netcdf topography file', FATAL)
            call read_data (topog_file(1:nc), 'fis', fis, Hgrid%Tmp%Domain)
            ! read res (if it exists)
            call field_size (topog_file(1:nc), 'res', siz )
            if ( siz(1) == Hgrid%nlon .and. siz(2) == Hgrid%nlat ) then
                call read_data (topog_file(1:nc), 'res', res, Hgrid%Tmp%Domain)
            else
                res = 1.
            endif
      !------------------------------------------
      ! read topography from native restart file
      
        else

            unit = open_restart_file (topog_file(1:nc), action='read')
            read (unit) ires,jres
            if ( ires /= Hgrid%nlon .or. jres /= Hgrid%nlat )  &
                            call error_mesg ('bgrid_cold_start_mod',  &
                        'incorrect resolution in native topography file', FATAL)
            call set_domain (Hgrid%Tmp%Domain)
            call read_data  (unit, fis)                ! read topog
            call read_data  (unit, res, end=eof_found) ! read res (if it exists)
            if (eof_found) res = 1.
            call close_file (unit)

        endif

   else
      ! restart file does not exist
        call error_mesg ('bgrid_cold_start',  &
                         'restart file for topography does not exist', FATAL)
   endif

!---------------------------------------------
! OPTION 4: error condition
!---------------------------------------------
  else
        call error_mesg ('bgrid_cold_start',  &
                         'invalid value for topog_option', FATAL)
  endif

! halos
  call update_halo ( Hgrid, TEMP, fis )
  call update_halo ( Hgrid, TEMP, res )

!-----------------------------------------------------------------------

 end subroutine compute_topog

!#######################################################################

 function compute_sigma (nlev) result (eta)

   integer, intent(in) :: nlev
   real, dimension(nlev+1) :: eta
   real :: dz, qk
   integer :: k

     eta(1) = 0.0
     eta(nlev+1) = 1.0

   if ( equal_vert_spacing ) then

     ! levels with equal sigma spacing
       dz = 1./nlev
       do k = 2, nlev
         eta(k) = eta(k-1) + dz
       enddo

   else

     ! space levels according to Smagorinski (1965, MWR, pp.727-768)
     ! sigma spacing smaller at the top and bottom of model
       do k = 1, nlev-1
         qk = real(2*k)/real(2*nlev)
         eta(k+1) = qk*qk*(3.-2.*qk)
       enddo

   endif

 end function compute_sigma

!#######################################################################

end module bgrid_cold_start_mod


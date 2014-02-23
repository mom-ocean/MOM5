
module climap_albedo_mod

!=======================================================================
!
!  Routines for obtaining land surface albedos from the CLIMAP data set
!
!=======================================================================

use  horiz_interp_mod, only:  horiz_interp

use           fms_mod, only:  open_ieee32_file, &
                              error_mesg, FATAL, file_exist,  &
                              check_nml_error, open_namelist_file,     &
                              close_file, mpp_pe, mpp_root_pe, &
                              write_version_number, stdlog, read_data, &
                              NOTE, mpp_error
use     constants_mod, only:  PI

implicit none
private

!=======================================================================

    public  get_climap_albedo,     get_climap_glacier
    public  get_climap_albedo_mcm, get_climap_glacier_mcm

!--------------------- version number ----------------------------------

character(len=128) :: version = '$Id: climap_albedo.F90,v 15.0 2007/08/14 04:00:26 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'
logical :: do_init = .true.

!=======================================================================

CONTAINS

!#######################################################################

subroutine get_climap_albedo (lonb, latb, n, albedo)

real   , intent(in) , dimension(:,:) :: lonb, latb
integer, intent(in)                  :: n
real   , intent(out), dimension(:,:) :: albedo
  
integer, parameter :: nlon = 180, nlat = 90

real :: pi, sb, wb, dx, dy
real, dimension(nlon, nlat) :: alb
integer :: unit, nn
character(len=32) :: lvltag
!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif
   do_init = .false.


pi = 4.0*atan(1.0)
sb = -pi/2.0               ! southern boundary of albedo data grid
wb = 0.0                   ! western  boundary of albedo data grid
dx = 2.0*pi/float(nlon)    ! longitudinal grid size of albedo data grid
dy = pi/float(nlat)        ! latitudinal  grid size of albedo data grid

if(n < 10)  write(lvltag, 1201) n
if(n >= 10) write(lvltag, 1202) n
if(n >= 100) write(lvltag, 1203) n
1201 format(i1)
1202 format(i2)
1203 format(i3)
if (file_exist('INPUT/albedo.data.nc') ) then
   if (mpp_pe() == mpp_root_pe()) call mpp_error ('climap_albedo_mod', &
            'Reading NetCDF formatted input data file: INPUT/albedo.data.nc', NOTE)
    call read_data('INPUT/albedo.data.nc', 'albedo_'//trim(lvltag), alb, no_domain=.true.)
    call horiz_interp_wrapper (alb,wb,sb,dx,dy,lonb,latb,albedo)
else
   if (file_exist('INPUT/albedo.data')) then
!    unit = open_file (file='INPUT/albedo.data',  &
!                      form='ieee32', action='read')
      if (mpp_pe() == mpp_root_pe()) call mpp_error ('climap_albedo_mod', &
            'Reading native formatted input data file: INPUT/albedo.data', NOTE)
      unit = open_ieee32_file (file='INPUT/albedo.data', action='read')
      do nn = 1, n
         read (unit,end=99,err=98)  alb
      end do
      call horiz_interp_wrapper (alb,wb,sb,dx,dy,lonb,latb,albedo)
      call close_file (unit)
   else
      call error_mesg ('get_climap_albedos',  &
           'eof, albedo.data file does not exist?', FATAL)
   endif
endif

return

!-----------------------------------------------------------------------
  98  call error_mesg ('get_climap_albedos',  &
                       'error in reading albedo.data file', FATAL)

  99  call error_mesg ('get_climap_albedos',  &
                'unexpected EOR when reading albedo.data file', FATAL)
!-----------------------------------------------------------------------

end subroutine get_climap_albedo

!#######################################################################

subroutine get_climap_glacier (lonb, latb, crit_albedo, glacier)

real, intent(in), dimension(:,:)     :: lonb, latb
real, intent(in)                     :: crit_albedo
logical, intent(out), dimension(:,:) :: glacier

real, dimension(size(lonb,1)-1,size(latb,2)-1) :: alb_dry

!----- read third record (max snow albedo) ------

call get_climap_albedo(lonb, latb, 3, alb_dry)

glacier = .false.
where(alb_dry >= crit_albedo) glacier = .true.

return

end subroutine get_climap_glacier

!#######################################################################

subroutine get_climap_albedo_mcm (lonb, latb, nlong, nlatg, is, ie, js, je, n, &
                                  albedo)

real   , intent(in) , dimension(:,:) :: lonb, latb
integer, intent(in)                  :: nlong, nlatg, is, ie, js, je, n
real   , intent(out), dimension(:,:) :: albedo
  
real :: pi, sb, wb, dx, dy
real, dimension(nlong, nlatg) :: alb
integer :: unit, nn, i, j

!--- write version id to logfile ---
if (do_init) then
   if (mpp_pe() == mpp_root_pe()) &
   call write_version_number(version, tagname)
   do_init = .false.


endif
if (file_exist('INPUT/ss_albedo.data')) then
!    unit = open_file (file='INPUT/ss_albedo.data',form='ieee32', action='read')
    unit = open_ieee32_file (file='INPUT/ss_albedo.data', action='read')
    do nn = 1, n
      read (unit,end=99,err=98)  alb
    end do
    do j=js,je
     do i=is,ie
      albedo(i-is+1,j-js+1) = alb(i,j)
     enddo
    enddo
    call close_file (unit)

else
    call error_mesg ('get_climap_albedo_mcm',  &
           'eof, ss_albedo.data file does not exist?', FATAL)
endif

return

!-----------------------------------------------------------------------
  98  call error_mesg ('get_climap_albedo_mcm',  &
                       'error in reading ss_albedo.data file', FATAL)

  99  call error_mesg ('get_climap_albedo_mcm',  &
                'unexpected EOR when reading ss_albedo.data file', FATAL)
!-----------------------------------------------------------------------

end subroutine get_climap_albedo_mcm

!#######################################################################

subroutine get_climap_glacier_mcm (lonb, latb, nlong, nlatg, is, ie, js, je, &
                                  crit_albedo, glacier)

real   , intent(in), dimension(:,:)     :: lonb, latb
integer, intent(in)                     :: nlong, nlatg, is, ie, js, je
real   , intent(in)                     :: crit_albedo
logical, intent(out), dimension(:,:)    :: glacier

real, dimension(size(lonb,1)-1,size(latb,2)-1) :: alb_dry

!----- read albedo ------

call get_climap_albedo_mcm(lonb, latb, nlong, nlatg, is, ie, js, je, 1, alb_dry)

glacier = .false.
where(alb_dry >= crit_albedo) glacier = .true.

return

end subroutine get_climap_glacier_mcm

!#######################################################################

subroutine horiz_interp_wrapper (data_in, wb, sb, dx, dy, lonb_out, latb_out, data_out)
real   , intent(in),  dimension(:,:) :: data_in, lonb_out, latb_out
real   , intent(in)                  :: wb, sb, dx, dy
real   , intent(out), dimension(:,:) :: data_out
real    :: lonb_in(size(data_in,1)+1), latb_in(size(data_in,2)+1), tpi
integer :: i, j, nx, ny

     nx = size(data_in,1)
     ny = size(data_in,2)
     tpi = 2.*PI

     ! longitude boundaries
     do i = 1, nx+1
       lonb_in(i) = wb + float(i-1)*dx
     enddo
       if (abs(lonb_in(nx+1)-lonb_in(1)-tpi) < epsilon(lonb_in)) &
               lonb_in(nx+1)=lonb_in(1)+tpi

     ! latitude boundaries
     do j = 2, ny
       latb_in(j) = sb + float(j-1)*dy
     enddo
       latb_in(1)    = -0.5*PI
       latb_in(ny+1) =  0.5*PI

     call horiz_interp (data_in, lonb_in, latb_in, &
                                 lonb_out, latb_out, data_out)

end subroutine horiz_interp_wrapper

!#######################################################################

end module climap_albedo_mod

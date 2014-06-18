!FDOC_TAG_GFDL
module microphys_cloud_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   
! </CONTACT>
! <REVIEWER EMAIL="">
!   
! </REVIEWER>
! <OVERVIEW>
!   
! </OVERVIEW>
! <DESCRIPTION>
!   
! </DESCRIPTION>
!

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod, only: input_nml_file
use fms_mod, only: open_namelist_file, mpp_pe,          &
                   mpp_root_pe, stdlog, fms_init,       &
                   write_version_number, file_exist,    &
                   check_nml_error, close_file

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private

public microphys_cloud, microphys_cloud_init, microphys_cloud_end

integer :: idim, jdim, kdim

!namelist:
!real :: diam_liq = 40.0  ! vtp
real :: diam_liq = 33.2

namelist /microphys_cloud_nml/  diam_liq

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: microphys_cloud.F90,v 19.0 2012/01/06 20:19:43 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'

logical ::   module_is_initialized = .false.
contains 

!#######################################################################

! <SUBROUTINE NAME="microphys_cloud">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call microphys_cloud ( zhalf, zfull, diam_liq_out, diam_ice)
!
!  </TEMPLATE>
!  <IN NAME="zhalf" TYPE="real">
! 
!  </IN>
!  <IN NAME="zfull" TYPE="real">
! 
!  </IN>
!  <OUT NAME="diam_liq_out" TYPE="real">
! 
!  </OUT>
!  <OUT NAME="diam_ice" TYPE="real">
! 
!  </OUT>
! </SUBROUTINE>
!
subroutine microphys_cloud ( zhalf, zfull, diam_liq_out, diam_ice)

real,    dimension(:,:,:),    intent(in)        :: zhalf, zfull
real,    dimension(:,:,:),    intent(out)       :: diam_liq_out,  &
                                                   diam_ice

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      zhalf        height asl at half levels [m]
!      zfull        height asl at full levels [m]
!
!   intent(out):
!
!      diam_liq
!      diam_ice
!
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

  idim = size (zfull,1)
  jdim = size (zfull,2)
  kdim = size (zfull,3)

  call get_diam ( zhalf, zfull, diam_ice )
  diam_liq_out = diam_liq

end subroutine microphys_cloud

!#######################################################################

! <SUBROUTINE NAME="get_diam">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call get_diam ( zhalf, zfull, diam )
!
!  </TEMPLATE>
!  <IN NAME="zhalf" TYPE="real">
! 
!  </IN>
!  <IN NAME="zfull" TYPE="real">
! 
!  </IN>
!  <OUT NAME="diam" TYPE="real">
! 
!  </OUT>
! </SUBROUTINE>
!
subroutine get_diam ( zhalf, zfull, diam )

real, dimension(:,:,:), intent(in)  :: zhalf
real, dimension(:,:,:), intent(in)  :: zfull
real, dimension(:,:,:), intent(out) :: diam

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(size(zfull,3)+1) :: relht

real :: slope, zref
integer :: i, j, k

  relht(kdim+1) = 0.0

  do i=1,idim
     do j=1,jdim

       zref = zhalf(i,j,kdim+1) + 9.9e3
       relht(1:kdim) = max(0.0, min(1.0, (zfull(i,j,:) - zref)/3.3e3 ))

        k = kdim+1
        do while ( relht(k) < 0.30 .and. k > 1 )
           k=k-1
           slope = (30.72 - 38.50)/(0.30 - 0.00)
           diam(i,j,k) = 38.50 + (relht(k) - 0.00) * slope
        enddo

        do while ( relht(k) < 0.45 .and. k > 1 )
           k=k-1
           slope = (28.28 - 30.72)/(0.45 - 0.30)
           diam(i,j,k) = 30.72 + (relht(k) - 0.30) * slope
        enddo

        do while ( relht(k) < 0.64 .and. k > 1 )
           k=k-1
           slope = (25.62 - 28.28)/(0.64 - 0.45)
           diam(i,j,k) = 28.28 + (relht(k) - 0.45) * slope
        enddo

        do while ( relht(k) < 0.76 .and. k > 1 )
           k=k-1
           slope = (24.80 - 25.62)/(0.76 - 0.64)
           diam(i,j,k) = 25.62 + (relht(k) - 0.64) * slope
        enddo

        do while ( k > 1 ) 
           k=k-1
!           slope = (13.30 - 24.80)/(1.00 - 0.76)                              !vtp
           slope = (18.60 - 24.80)/(1.00 - 0.76)
           diam(i,j,k) = 24.80 + (relht(k) - 0.76) * slope
        enddo

     enddo
  enddo

  return
end subroutine get_diam

!#######################################################################

! <SUBROUTINE NAME="microphys_cloud_init">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call microphys_cloud_init
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine microphys_cloud_init

integer :: unit, ierr, io, logunit

   if (module_is_initialized ) return
   call fms_init

!-----------------------------------------------------------------------
! read namelist       
!-----------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=microphys_cloud_nml, iostat=io)
  ierr = check_nml_error(io,"microphys_cloud_nml")
#else
  if (file_exist('input.nml')) then
     unit =  open_namelist_file ( )
     ierr=1; do while (ierr /= 0)
     read (unit, nml=microphys_cloud_nml, iostat=io, end=10) 
     ierr = check_nml_error (io, 'microphys_cloud_nml')
     enddo                       
10   call close_file (unit)      
  endif                         
#endif
                                    
!------------------------------------------------------------------------
!  write version number and namelist to logfile.
!------------------------------------------------------------------------

  call write_version_number (version, tagname)
  logunit = stdlog()
  if (mpp_pe() == mpp_root_pe() )                                       &
                                write (logunit, nml=microphys_cloud_nml)

   module_is_initialized = .true.

  return
end subroutine microphys_cloud_init

!#######################################################################
subroutine microphys_cloud_end

   module_is_initialized = .true.

  return
end subroutine microphys_cloud_end

!#######################################################################

end module microphys_cloud_mod


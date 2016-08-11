module ocean_increment_velocity_mod
!
!<CONTACT EMAIL="russell.fiedler@csiro.au">  Russell Fiedler
!</CONTACT>
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov.au"> Paul Sandery
!</CONTACT>
!
!<OVERVIEW>
! Thickness weighted velocity tendency [meter^2/sec^2] from increments.
!</OVERVIEW>
!
!<DESCRIPTION>
!
! This increment module performs incremental update analysis (IUA), 
! an approach used in data assimilation and forecasting to reduce
! spurious perturbations when correcting the model state. 
! IUA involves restoring to analysis increments i.e. differences between model
! and analysis fields rather than actual fields (See Bloom et al., 1996). 
! The user can define the period that IUA is carried out
! and also the fraction of the increment to be restored over that period.
!
! This module applies a general 3D source to velocities. The sources
! can occur at any location and with any distribution in the domain
! An array of velocity tendencies due to the increments is augmented through a
! call to increment_velocity_source.  The array of velocity tendencies must be
! reset to zero between calls.
!
! The user is responsible for providing (and registering) the data on
! the model grid of values towards which the velocitys are being driven.
!</DESCRIPTION>
!
! <REFERENCE>
! S.C. Bloom, L.L. Takacs, A.M. da Silva, and D. Ledvina
! Data Assimilation Using Incremental Analysis Updates
! Monthly Weather Review  Volume 124, Issue 6 (June 1996)
! pages 1256--1271 
! </REFERENCE>
!
!<NAMELIST NAME="ocean_increment_velocity_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  For using this module.  Default use_this_module=.false.
!  </DATA>
!  <DATA NAME="fraction_increment" TYPE="real">
!  For prescribing the fraction of the increment.
!  applied to the restoring period.  Default fraction_increment=1.0
!  </DATA>
!  <DATA NAME="days_to_increment" TYPE="integer">
!  For specifying the amount of days to restore.
!  Default days_to_increment=1
!  </DATA>
!  <DATA NAME="secs_to_increment" TYPE="integer">
!  For specifying the amount of seconds to restore.
!  Default secs_to_increment=0
!  </DATA>
!</NAMELIST>
!
use diag_manager_mod,         only: register_diag_field
use fms_mod,                  only: write_version_number, open_namelist_file, close_file
use fms_mod,                  only: file_exist
use fms_mod,                  only: open_namelist_file, check_nml_error, close_file
use fms_mod,                  only: read_data, lowercase, FATAL, WARNING, stdout, stdlog
use mpp_mod,                  only: input_nml_file, mpp_sum, mpp_error
use time_interp_external_mod, only: init_external_field, time_interp_external
use time_manager_mod,         only: time_type, set_date,get_time
use time_manager_mod,         only: operator( + ), operator( - ), operator( // )
use time_manager_mod,         only: operator( > ), operator( == ), operator( <= )

use ocean_domains_mod,        only: get_local_indices
use ocean_parameters_mod,     only: missing_value 
use ocean_types_mod,          only: ocean_domain_type, ocean_grid_type, ocean_thickness_type
use ocean_types_mod,          only: ocean_velocity_type, ocean_time_type 
use ocean_workspace_mod,      only: wrk1, wrk2
use ocean_util_mod,           only: diagnose_3d_u

implicit none

private

#include <ocean_memory.h>


type ocean_increment_type
   integer :: id                                             ! time_interp_external index
   character(len=32) :: name                                 ! velocity name corresponding to increment
   real, dimension(:,:,:), pointer :: damp_coeff   => NULL() ! 3d inverse forcing rate (velocity units/ sec)
end type ocean_increment_type


type(ocean_increment_type), allocatable, dimension(:) :: Increment
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

public ocean_increment_velocity_init
public ocean_increment_velocity_source

character(len=126)  :: version = '$Id: ocean_increment_velocity.F90,v 20.0 2013/12/14 00:16:00 fms Exp $'
character (len=128) :: tagname = '$Name: tikal $'

! for diagnostics 
logical :: used
integer, dimension(:), allocatable :: id_increment_tend

logical :: module_is_initialized = .false.
logical :: damp_coeff_3d         = .false. 
logical :: use_this_module       = .false. 
real    :: time_scale
real    :: fraction_increment    = 1.0
integer :: days_to_increment     = 1
integer :: secs_to_increment     = 0


namelist /ocean_increment_velocity_nml/ use_this_module, fraction_increment, &
                                        days_to_increment, secs_to_increment

! Time info

integer :: days, secs
integer :: days_end_increment, secs_end_increment

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_increment_velocity_init">
!
! <DESCRIPTION>
! This subroutine is intended to be used to initialize the increments.
! Everything in this subroutine is a user prototype, and should be replacable.
! </DESCRIPTION>
!
subroutine ocean_increment_velocity_init(Grid, Domain, Time)

  type(ocean_grid_type),        intent(in), target :: Grid
  type(ocean_domain_type),      intent(in), target :: Domain
  type(ocean_time_type),        intent(in)         :: Time

  integer :: i, j, k, n
  integer :: ioun, io_status, ierr

  character(len=128) :: name

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_increment_mod (ocean_increment_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  allocate( Increment(2) )

  call write_version_number(version, tagname)
  
  ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_Increment_velocity_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_increment_velocity_nml')
#else
  ioun =  open_namelist_file()
  read (ioun,ocean_Increment_velocity_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_increment_velocity_nml')
  call close_file (ioun)
#endif
  write (stdlogunit,ocean_increment_velocity_nml)
  write (stdoutunit,'(/)') 
  write (stdoutunit,ocean_increment_velocity_nml)


  Dom => Domain
  Grd => Grid

#ifndef MOM_STATIC_ARRAYS    
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif

  do n=1,2
     Increment(n)%id = -1
  enddo


  if(.not. use_this_module) return 

  if ( days_to_increment  < 0  .or. secs_to_increment < 0 .or. &
       days_to_increment==0 .and. secs_to_increment == 0 ) then
      call mpp_error(FATAL,&
      '==>Error: invalid restoring period, ensure days_to_increment + secs_to_increment is greater than zero')
  endif

  time_scale=fraction_increment/(days_to_increment*86400.0 + secs_to_increment)

  call get_time(Time%model_time,secs,days)
  days_end_increment= days + days_to_increment
  secs_end_increment=secs + secs_to_increment


  do n=1,2

     allocate(Increment(n)%damp_coeff(isd:ied,jsd:jed,nk))
     Increment(n)%damp_coeff(:,:,:) = 0.0  

     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              if(Grd%umask(i,j,k) == 0.0) then
                  Increment(n)%damp_coeff(i,j,k) = 0.0
              else
                  Increment(n)%damp_coeff(i,j,k) = time_scale
              endif
           enddo
        enddo
     enddo

     ! read forcing rates (1/sec) 
     if ( n == 1 ) then
         name = 'INPUT/u_increment.nc'
     else
         name = 'INPUT/v_increment.nc'
     endif
     if (file_exist(trim(name)) ) then
         if ( n == 1 ) then
             Increment(n)%id = init_external_field(name,'u_inc',domain=Domain%domain2d)
         else
             Increment(n)%id = init_external_field(name,'v_inc',domain=Domain%domain2d)
         endif
         if (Increment(n)%id < 1) then 
             call mpp_error(FATAL,&
             '==>Error: in ocean_increment_velocity_mod: forcing rates are specified but increment values are not')
         endif
         write(stdoutunit,*) '==> Using increment data specified from file '//trim(name) 
     else
         write(stdoutunit,*) '==> '//trim(name)//' not found.  Increment not being applied '
     endif

  enddo

  ! register diagnostic output
  allocate (id_increment_tend(2))
  id_increment_tend = -1

  do n=1,2
     if(n==1) then
        id_increment_tend(n) = register_diag_field ('ocean_model', 'u_increment_tend',        &
                   Grd%vel_axes_uv(1:3), Time%model_time, 'rho*dzt*tendency due to increment',&
                   'kg/m s^2', missing_value=missing_value, range=(/-1.e10,1.e10/))
     else
        id_increment_tend(n) = register_diag_field ('ocean_model', 'v_increment_tend',        &
                   Grd%vel_axes_uv(1:3), Time%model_time, 'rho*dzt*tendency due to increment',&
                   'kg/m s^2', missing_value=missing_value, range=(/-1.e10,1.e10/))
     endif
  enddo


end subroutine ocean_increment_velocity_init
! </SUBROUTINE> NAME="ocean_increment_velocity_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_increment_velocity_source">
!
! <DESCRIPTION>
! This subroutine calculates thickness weighted and density weighted
! time tendencies of velocitys due to damping by increments.
! </DESCRIPTION>
!
subroutine ocean_increment_velocity_source(Time, Thickness, Velocity)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_velocity_type), intent(inout)    :: Velocity

  integer :: taum1, tau
  integer :: i, j, k, n

  if(.not. use_this_module) return 

  taum1 = Time%taum1
  tau   = Time%tau
  wrk1  = 0.0 
  wrk2  = 0.0

! Halt forcing

  call get_time(Time%model_time,secs,days)
  if (days + secs/86400.0 .ge. days_end_increment + secs_end_increment/86400.0 ) then
    Increment(:)%id=0
  endif

 
  do n=1,2

     if (Increment(n)%id > 0) then

         ! get increment value for current time
         call time_interp_external(Increment(n)%id, Time%model_time, wrk1) 
 
         do k=1,nk
            do j=jsd,jed
               do i=isd,ied  
                  wrk2(i,j,k) = Thickness%rho_dzu(i,j,k,tau)*Increment(n)%damp_coeff(i,j,k) &
                                *wrk1(i,j,k)
                  Velocity%accel(i,j,k,n) = Velocity%accel(i,j,k,n) + wrk2(i,j,k)
               enddo
            enddo
         enddo

     endif

     call diagnose_3d_u(Time, Grd, id_increment_tend(n), wrk2(:,:,:))
  enddo

  return

end subroutine ocean_increment_velocity_source
! </SUBROUTINE> NAME="ocean_increment_velocity_source"


end module ocean_increment_velocity_mod

module ocean_increment_eta_mod
!
!<CONTACT EMAIL="russell.fiedler@csiro.au">  Russell Fiedler
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov.au"> Paul Sandery 
!</CONTACT>
!
!<OVERVIEW>
! Thickness weighted eta tendency [meter^2/sec^2] from increments.
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
! This module applies a general 2D source to eta. The sources
! can occur at any location and with any distribution in the domain
! An array of eta tendencies due to the increments is augmented through a
! call to increment_eta_source.  The array of eta tendencies must be
! reset to zero between calls.
!
! The user is responsible for providing (and registering) the data on
! the model grid of values towards which the etas are being driven.
!</DESCRIPTION>
!
! <REFERENCE>
! S.C. Bloom, L.L. Takacs, A.M. da Silva, and D. Ledvina
! Data Assimilation Using Incremental Analysis Updates
! Monthly Weather Review  Volume 124, Issue 6 (June 1996)
! pages 1256--1271 
! </REFERENCE>
!
!<NAMELIST NAME="ocean_increment_eta_nml">
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
use time_manager_mod,         only: time_type, set_date, get_time
use time_manager_mod,         only: operator( + ), operator( - ), operator( // )
use time_manager_mod,         only: operator( > ), operator( == ), operator( <= )

use ocean_domains_mod,        only: get_local_indices
use ocean_parameters_mod,     only: missing_value 
use ocean_types_mod,          only: ocean_domain_type, ocean_grid_type, ocean_thickness_type
use ocean_types_mod,          only: ocean_External_mode_type, ocean_time_type 
use ocean_workspace_mod,      only: wrk1_2d, wrk2_2d
use ocean_util_mod,           only: diagnose_2d
use ocean_parameters_mod,     only: rho0
implicit none

private

#include <ocean_memory.h>


type ocean_increment_type
   integer :: id                                             ! time_interp_external index
   character(len=32) :: name                                 ! eta name corresponding to increment
   real, dimension(:,:), pointer :: damp_coeff   => NULL() ! 2D inverse forcing rate (eta units/ sec)
end type ocean_increment_type


type(ocean_increment_type),SAVE :: Increment
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

public ocean_increment_eta_init
public ocean_increment_eta_source

character(len=126)  :: version = '$Id: ocean_increment_eta.F90,v 20.0 2013/12/14 00:15:56 fms Exp $'
character (len=128) :: tagname = '$Name: tikal $'

! for diagnostics 
logical :: used
integer :: id_increment_tend
logical :: module_is_initialized = .false.
logical :: use_this_module       = .false. 
real    :: time_scale
real    :: fraction_increment    = 1.0
integer :: days_to_increment     = 1
integer :: secs_to_increment     = 0

namelist /ocean_increment_eta_nml/ use_this_module, fraction_increment, days_to_increment, secs_to_increment

! Time info
integer :: days_end_increment, secs_end_increment
integer :: days, secs

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_increment_eta_init">
!
! <DESCRIPTION>
! This subroutine is intended to be used to initialize the increments.
! Everything in this subroutine is a user prototype, and should be replacable.
! </DESCRIPTION>
!
subroutine ocean_increment_eta_init(Grid, Domain, Time)

  type(ocean_grid_type),        intent(in), target :: Grid
  type(ocean_domain_type),      intent(in), target :: Domain
  type(ocean_time_type),        intent(in)         :: Time

  integer :: i, j
  integer :: ioun, io_status, ierr

  character(len=128) :: name

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 


  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_increment_mod (ocean_increment_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.


  call write_version_number(version, tagname)

  ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_Increment_eta_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_increment_eta_nml')
#else
  ioun =  open_namelist_file()
  read (ioun,ocean_Increment_eta_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_increment_eta_nml')
  call close_file (ioun)
#endif
  write (stdlogunit,ocean_Increment_eta_nml)
  write (stdoutunit,'(/)') 
  write (stdoutunit,ocean_Increment_eta_nml)

  Dom => Domain
  Grd => Grid

#ifndef MOM_STATIC_ARRAYS    
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif

  allocate(Increment%damp_coeff(isd:ied,jsd:jed))
  Increment%id = -1

  if(.not. use_this_module) return 

  if ( days_to_increment  < 0  .or. secs_to_increment < 0 .or. &
  days_to_increment==0 .and. secs_to_increment == 0 ) then

     call mpp_error(FATAL,&
      '==>Error: invalid restoring period, ensure days_to_increment + secs_to_increment is greater than zero')
  endif

  time_scale=fraction_increment/(days_to_increment*86400.0 + secs_to_increment)

  call get_time(Time%model_time,secs,days)
  days_end_increment=days + days_to_increment
  secs_end_increment=secs + secs_to_increment


  ! read forcing rates (1/sec) 
  Increment%damp_coeff(:,:) = 0.0  

  do j=jsd,jed
     do i=isd,ied
        if(Grd%tmask(i,j,1) == 0.0) then 
            Increment%damp_coeff(i,j) = 0.0
        else
            Increment%damp_coeff(i,j) = time_scale
        endif
     enddo
  enddo

  ! read forcing data
  name = 'INPUT/eta_increment.nc'
  if(file_exist(trim(name))) then
      Increment%id = init_external_field(name,'eta_inc',domain=Domain%domain2d)
      if (Increment%id < 1) then 
          call mpp_error(FATAL,&
               '==>Error: in ocean_increment_eta_mod: forcing rates are specified but increment values are not')
      endif
      write(stdoutunit,*) '==> Using increment data specified from file '//trim(name) 
  else
      write(stdoutunit,*) '==> '//trim(name)//' not found.  Increment not being applied '
  endif

  ! register diagnostic output
  id_increment_tend = -1

  id_increment_tend = register_diag_field ('ocean_model', 'eta_t_increment_tend', &
                Grd%tracer_axes(1:2), Time%model_time, 'eta_t tendency due to increment', &
                'm/s', missing_value=missing_value, range=(/-1.e10,1.e10/))


end subroutine ocean_increment_eta_init
! </SUBROUTINE> NAME="ocean_increment_eta_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_increment_eta_source">
!
! <DESCRIPTION>
! This subroutine calculates thickness weighted and density weighted
! time tendencies of etas due to damping by increments.
! </DESCRIPTION>
!
subroutine ocean_increment_eta_source(Time, Ext_mode)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_external_mode_type), intent(inout) :: Ext_mode

  integer :: taum1, tau
  integer :: i, j
  integer :: days,secs

  if(.not. use_this_module) return 

  taum1 = Time%taum1
  tau   = Time%tau
  wrk1_2d  = 0.0 
  wrk2_2d  = 0.0

  ! halt forcing
  call get_time(Time%model_time,secs,days)
  if(days + secs/86400.0 .ge. days_end_increment + secs_end_increment/86400.0 ) then
    Increment%id=0
  endif

  if (Increment%id > 0) then

     ! get increment value for current time
     call time_interp_external(Increment%id, Time%model_time, wrk1_2d) 

     do j=jsd,jed
         do i=isd,ied  
            wrk2_2d(i,j) = Increment%damp_coeff(i,j)    &
                      *wrk1_2d(i,j)
           Ext_mode%source(i,j) = Ext_mode%source(i,j)+rho0*wrk2_2d(i,j) 
         enddo
     enddo


   endif

   call diagnose_2d(Time, Grd, id_increment_tend, wrk2_2d(:,:))

  return

end subroutine ocean_increment_eta_source
! </SUBROUTINE> NAME="increment_eta_source"


end module ocean_increment_eta_mod

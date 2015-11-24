module ocean_sponges_eta_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov.au"> Paul Sandery
!</CONTACT>
!
!<OVERVIEW>
! Weighted eta tendency [meter*meter/sec] from sponges.
!</OVERVIEW>
!
!<DESCRIPTION>
!
! This module applies sponge to eta. The sponges
! can occur at any location and with any distribution in the domain, and
! with any time step and damping rate.  Sponges occur where positive
! inverse restore times occur in the field passed to sponge_init.  An
! array of eta tendencies due to the sponges is augmented through a
! call to sponge_eta_source.  The array of eta tendencies must be
! reset to zero between calls.
!
! Different damping rates can be specified by making
! calls to register_sponge_rate - no sponges are applied to fields for
! which uniformly zero inverse damping rates are set with a call to
! register_sponge_rate.  The value towards which a field is damped is
! set with calls to register_sponge_field; successive calls are used to
! set up linear interpolation of this restore rate.
!
! Sponge data and damping coefficients are 2 dimensional. 
!
! The user is responsible for providing (and registering) the data on
! the model grid of values towards which the etas are being driven.
!
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_sponges_eta_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  For using this module.  Default use_this_module=.false.
!  </DATA> 
!
!  <DATA NAME="damp_coeff_3d" TYPE="logical">
!  For case when damping coefficients are full 3d field of values.
!  Default damp_coeff_3d=.false., which means damping coeffs are 
!  2d horizontal array.   
!  </DATA> 
!
!</NAMELIST>
!
use diag_manager_mod,         only: register_diag_field
use fms_mod,                  only: write_version_number, open_namelist_file, close_file
use fms_mod,                  only: file_exist
use fms_mod,                  only: open_namelist_file, check_nml_error, close_file
use fms_mod,                  only: read_data, lowercase, FATAL, WARNING, stdout, stdlog
use mpp_mod,                  only: input_nml_file, mpp_sum, mpp_error, mpp_max
use mpp_domains_mod,          only: mpp_global_sum
use time_interp_external_mod, only: init_external_field, time_interp_external
use time_manager_mod,         only: time_type, set_date, get_time
use time_manager_mod,         only: operator( + ), operator( - ), operator( // )
use time_manager_mod,         only: operator( > ), operator( == ), operator( <= )

use ocean_domains_mod,        only: get_local_indices
use ocean_parameters_mod,     only: missing_value, rho0 
use ocean_types_mod,          only: ocean_domain_type, ocean_grid_type, ocean_options_type
use ocean_types_mod,          only: ocean_time_type, ocean_External_mode_type
use ocean_workspace_mod,      only: wrk1_2d, wrk2_2d
use ocean_util_mod,           only: diagnose_2d

implicit none

private

#include <ocean_memory.h>

type ocean_sponge_type
   integer :: id                                             ! time_interp_external index
   character(len=32) :: name                                 ! eta name corresponding to sponge
   real, dimension(:,:)  , pointer :: damp_coeff2d => NULL() ! 2d inverse damping rate (eta units/ sec)
end type ocean_sponge_type

type(ocean_sponge_type), allocatable, dimension(:) :: Sponge
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

public ocean_sponges_eta_init
public sponge_eta_source

character(len=126)  :: version = '$Id: ocean_sponges_eta.F90,v 20.0 2013/12/14 00:16:22 fms Exp $'
character (len=128) :: tagname = '$Name: tikal $'

! for diagnostics 
logical :: used
integer, dimension(:), allocatable :: id_sponge_tend
logical :: module_is_initialized = .false.
logical :: use_this_module       = .false. 

! Adaptive restoring
real    :: athresh               = 0.5
real    :: npower                = 1.0
real    :: sdiffo                = 0.5
real    :: lambda                = 0.0083
real    :: taumin                = 720
logical :: use_adaptive_restore  = .false.
logical :: use_sponge_after_init = .false.
logical :: use_normalising       = .false.
logical :: use_hard_thump        = .false.
integer :: days_to_restore       = 1
integer :: secs_to_restore       = 0

integer :: secs_restore
integer :: initial_day, initial_secs

namelist /ocean_sponges_eta_nml/ use_this_module 
namelist /ocean_sponges_eta_ofam_nml/ &
    use_adaptive_restore, use_sponge_after_init, use_normalising, &
    use_hard_thump, athresh, taumin, lambda, npower, days_to_restore, &
    secs_to_restore

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_sponges_eta_init">
!
! <DESCRIPTION>
! This subroutine is intended to be used to initialize the sponges.
! Everything in this subroutine is a user prototype, and should be replacable.
! </DESCRIPTION>
!
subroutine ocean_sponges_eta_init(Grid, Domain, Time, dtime, Ocean_options)

  type(ocean_grid_type),          intent(in), target :: Grid
  type(ocean_domain_type),        intent(in), target :: Domain
  type(ocean_time_type),          intent(in)         :: Time
  real,                           intent(in)         :: dtime
  type(ocean_options_type),       intent(inout)      :: Ocean_options

  integer :: i, j
  integer :: ioun, io_status, ierr
  integer :: secs, days
  real    :: dtimer

  character(len=128) :: name

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_sponges_eta_mod (ocean_sponges_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  allocate( Sponge(1) )

  call write_version_number(version, tagname)

  ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
  read(input_nml_file, nml=ocean_sponges_eta_nml, iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_sponges_eta_nml')
  read(input_nml_file, nml=ocean_sponges_eta_ofam_nml, iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_sponges_eta_ofam_nml')
#else
  ioun =  open_namelist_file()
  read(ioun, ocean_sponges_eta_nml, IOSTAT=io_status)
  ierr = check_nml_error(io_status, 'ocean_sponges_eta_nml')
  rewind(ioun)
  read(ioun, ocean_sponges_eta_ofam_nml, IOSTAT=io_status)
  ierr = check_nml_error(io_status, 'ocean_sponges_eta_nml')
  call close_file (ioun)
#endif
  write (stdlogunit, ocean_sponges_eta_nml)
  write (stdoutunit, '(/)')
  write (stdoutunit, ocean_sponges_eta_nml)

  write (stdlogunit, ocean_sponges_eta_ofam_nml)
  write (stdoutunit, '(/)')
  write (stdoutunit, ocean_sponges_eta_ofam_nml)

  Dom => Domain
  Grd => Grid

#ifndef MOM_STATIC_ARRAYS    
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif

  Sponge(1)%id = -1
  dtimer = 1.0/dtime

  if(use_this_module) then
      write(stdoutunit,*)'==>Note from ocean_sponges_eta_mod: Using this module.'
      Ocean_options%ocean_sponges_eta= 'Used ocean eta sponges.'
  else
      write(stdoutunit,*)'==>Note from ocean_sponges_eta_mod: NOT using this module.'
      Ocean_options%ocean_sponges_eta= 'Did NOT use ocean eta sponges.'
      return
  endif

  ! Set up some initial times
  call get_time(Time%model_time, secs, days)
  initial_day = days
  initial_secs = secs

  !read damping rates for eta
  if (use_adaptive_restore) then
     allocate(Sponge(1)%damp_coeff2d(isd:ied, jsd:jed))
     Sponge(1)%damp_coeff2d(:,:) = 0.0
  else
     name = 'INPUT/eta_sponge_coeff.nc'
     if (file_exist(name)) then
        write(stdoutunit,*) '==> Using sponge damping times specified from file ',name
        allocate(Sponge(1)%damp_coeff2d(isd:ied,jsd:jed))
        Sponge(1)%damp_coeff2d(:,:) = 0.0
        call read_data(name,'coeff',Sponge(1)%damp_coeff2d,domain=Dom%domain2d,timelevel=1)

           do j=jsc,jec
              do i=isc,iec
                 if(Grd%tmask(i,j,1) == 0.0) then
                     Sponge(1)%damp_coeff2d(i,j) = 0.0
                 endif
              enddo
           enddo


        ! modify damping rates to allow restoring to be solved implicitly
        ! note: test values between zero and 4.0e-8 revert to damping rates defined above
           do j=jsc,jec
              do i=isc,iec
                 if (dtime*Sponge(1)%damp_coeff2d(i,j) > 4.0e-8) then
                     if (dtime*Sponge(1)%damp_coeff2d(i,j) > 37.0) then
                         Sponge(1)%damp_coeff2d(i,j) = dtimer
                     else
                         Sponge(1)%damp_coeff2d(i,j) = (1.0 - exp(-dtime*Sponge(1)%damp_coeff2d(i,j))) * dtimer
                     endif
                 else if (dtime*Sponge(1)%damp_coeff2d(i,j) <= 0.0) then
                     Sponge(1)%damp_coeff2d(i,j) = 0.0
                 endif
              enddo
           enddo

     endif  ! endif for fileexist 'INPUT/eta_sponge_coeff.nc'
  endif


  name = 'INPUT/eta_sponge.nc'
  if(file_exist(trim(name))) then
     Sponge(1)%id = init_external_field(name,'eta_t',domain=Domain%domain2d)
     if (Sponge(1)%id < 1) then
        if ( use_adaptive_restore ) then
          call mpp_error(FATAL,&
             '==>Error: in ocean_sponges_eta_mod: adaptive restoring is specified but sponge values are not')
        else
          call mpp_error(FATAL,&
             '==>Error: in ocean_sponges_eta_mod: damping rates are specified but sponge values are not')
        endif
     endif
     write(stdoutunit,*) '==> Using increment data specified from file '//trim(name)
  else
     write(stdoutunit,*) '==> '//trim(name)//' not found.  Increment not being applied '
  endif

  allocate (id_sponge_tend(1)) 
  id_sponge_tend = -1

  id_sponge_tend = register_diag_field ('ocean_model', 'eta_t_sponge_tend', &
     Grd%tracer_axes(1:2), Time%model_time, 'eta_t tendency due to sponge', &
     'm/s', missing_value=missing_value, range=(/-1.e10,1.e10/))


end subroutine ocean_sponges_eta_init
! </SUBROUTINE> NAME="ocean_sponges_init"


!#######################################################################
! <SUBROUTINE NAME="sponge_eta_source">
!
! <DESCRIPTION>
! This subroutine calculates thickness weighted and density weighted
! time tendencies due to damping by sponges or damping through adaptive
! restoring.
! </DESCRIPTION>
!
subroutine sponge_eta_source(Time, Ext_mode)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_external_mode_type), intent(inout) :: Ext_mode

  integer :: taum1, tau
  integer :: i, j
   
  real    :: sdiff, sum_val, adaptive_coeff
  integer :: secs, days, numsecs

  logical :: do_adaptive_restore = .false.  ! Only restore in specified time period.
  logical :: do_normal_restore = .false.    ! Only restore in specified time period.
  logical, save :: first_pass = .true.
  real, save :: num_surface_wet_cells

  if(.not. use_this_module) return 

  taum1 = Time%taum1
  tau   = Time%tau
  wrk1_2d = 0.0
  wrk2_2d = 0.0

  if (use_adaptive_restore) then
    if (first_pass) then
        num_surface_wet_cells = mpp_global_sum(Dom%domain2d, Grd%tmask(:,:,1))
    end if

    secs_restore = days_to_restore * 86400 + secs_to_restore
    sdiff = 0.0
    call get_time(Time%model_time, secs, days)
    numsecs = (days - initial_day) * 86400 + secs - initial_secs

    do_adaptive_restore = (numsecs < secs_restore)
    do_normal_restore = (.not. do_adaptive_restore) .and. use_sponge_after_init
  else
    do_normal_restore = .true.
  end if

  if (Sponge(1)%id > 0) then

    call time_interp_external(Sponge(1)%id, Time%model_time, wrk1_2d)

    if (do_adaptive_restore) then
        if (first_pass) then
           if (use_hard_thump) then
              do j = jsd, jed
                 do i = isd, ied
                       Ext_mode%eta_t(i, j, taum1) = wrk1_2d(i, j)
                       Ext_mode%eta_t(i, j, tau) = wrk1_2d(i, j)
                 end do
              end do
           end if

           wrk2_2d = abs(Ext_mode%eta_t(:,:,taum1) - wrk1_2d(:,:)) &
                     * Grd%tmask(:,:,1)
           sum_val = sum(wrk2_2d(isc:iec, jsc:jec))

           call mpp_sum(sum_val)
           sdiffo = sum_val / num_surface_wet_cells
           write (stdout(), *) 'mean eta_diff', sdiffo

           if (.not. use_normalising) then
              sdiffo = 1.0
           else
             sdiffo = maxval(wrk2_2d)
             call mpp_max(sdiffo)
             sdiffo = 1.01 * sdiffo
           end if

           first_pass = .false.
        end if

        ! Calculate idealised forcing tendency timescale for each timestep and
        ! array element
        do j = jsd, jed
           do i = isd, ied
              sdiff = abs(Ext_mode%eta_t(i, j, taum1) - wrk1_2d(i, j))
              sdiff = max(sdiff, 1e-9) !minimum tolerance
              adaptive_coeff = 1. / (taumin - (lambda * real(secs_restore)) &
                                     / (((sdiff / sdiffo)**npower) &
                                        * log(1 - athresh)))
              wrk2_2d(i,j) =  adaptive_coeff &
                              * (wrk1_2d(i, j) - Ext_mode%eta_t(i, j, taum1))
           end do
        end do

    else if (do_normal_restore) then
        do j = jsc, jec
           do i = isc, iec
              wrk2_2d(i,j) = Sponge(1)%damp_coeff2d(i, j) &
                             * (wrk1_2d(i, j) - Ext_mode%eta_t(i, j, taum1))
           end do
        end do
    end if

    do j = jsc, jec
        do i = isc, iec
            Ext_mode%source(i,j) = Ext_mode%source(i,j) + rho0 * wrk2_2d(i,j)
        end do
    end do
  end if

  call diagnose_2d(Time, Grd, id_sponge_tend(1), wrk2_2d(:,:))

  return

end subroutine sponge_eta_source
! </SUBROUTINE> NAME="sponge_eta_source"


end module ocean_sponges_eta_mod

module ocean_sponges_tracer_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Bonnie Samuels 
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> R.W. Hallberg 
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> M.J. Harrison 
!</CONTACT>
!
!<CONTACT EMAIL="swathi@cmmacs.ernet.in"> P. S. Swathi  
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Thickness weighted tracer tendency [tracer*meter/sec] from sponges.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module applies sponges to tracers. The sponges
! can occur at any location and with any distribution in the domain, and
! with any time step and damping rate.  Sponges occur where positive
! inverse restore times occur in the field passed to sponge_init.  An
! array of tracer tendencies due to the sponges is augmented through a
! call to sponge_tracer_source.  The array of tracer tendencies must be
! reset to zero between calls.
!
! Different damping rates can be specified for each field by making
! calls to register_sponge_rate - no sponges are applied to fields for
! which uniformly zero inverse damping rates are set with a call to
! register_sponge_rate.  The value towards which a field is damped is
! set with calls to register_sponge_field; successive calls are used to
! set up linear interpolation of this restore rate.
!
! Sponge data and damping coefficients are generally 3 dimensional. 
!
! The user is responsible for providing (and registering) the data on
! the model grid of values towards which the tracers are being driven.
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_sponges_tracer_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  For using this module.  Default use_this_module=.false.
!  </DATA> 
!  <DATA NAME="damp_coeff_3d" TYPE="logical">
!  For case when damping coefficients are full 3d field of values.
!  Default damp_coeff_3d=.false., which means damping coeffs are 
!  2d horizontal array.   
!  </DATA> 
!</NAMELIST>
!
use diag_manager_mod,         only: register_diag_field
use fms_mod,                  only: write_version_number, open_namelist_file, close_file
use fms_mod,                  only: file_exist
use fms_mod,                  only: open_namelist_file, check_nml_error, close_file
use fms_mod,                  only: read_data, lowercase, FATAL, WARNING, stdout, stdlog
use mpp_mod,                  only: input_nml_file, mpp_sum, mpp_error, mpp_max
use time_interp_external_mod, only: init_external_field, time_interp_external
use time_manager_mod,         only: time_type, set_date, get_time
use time_manager_mod,         only: operator( + ), operator( - ), operator( // )
use time_manager_mod,         only: operator( > ), operator( == ), operator( <= )

use ocean_domains_mod,        only: get_local_indices
use ocean_parameters_mod,     only: missing_value 
use ocean_types_mod,          only: ocean_domain_type, ocean_grid_type, ocean_thickness_type
use ocean_types_mod,          only: ocean_prog_tracer_type, ocean_options_type, ocean_time_type 
use ocean_workspace_mod,      only: wrk1, wrk2
use ocean_util_mod,           only: diagnose_3d

implicit none

private

#include <ocean_memory.h>

type ocean_sponge_type
   integer :: id                                             ! time_interp_external index
   character(len=32) :: name                                 ! tracer name corresponding to sponge
   real, dimension(:,:,:), pointer :: damp_coeff   => NULL() ! 3d inverse damping rate (tracer units/ sec)
   real, dimension(:,:)  , pointer :: damp_coeff2d => NULL() ! 2d inverse damping rate (tracer units/ sec)
end type ocean_sponge_type


type(ocean_sponge_type), allocatable, dimension(:) :: Sponge
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

public ocean_sponges_tracer_init
public sponge_tracer_source

character(len=126)  :: version = '$Id: ocean_sponges_tracer.F90,v 20.0 2013/12/14 00:16:24 fms Exp $'
character (len=128) :: tagname = '$Name: tikal $'

! for diagnostics 
logical :: used
integer, dimension(:), allocatable :: id_sponge_tend

integer :: num_prog_tracers      = 0
logical :: module_is_initialized = .FALSE.
logical :: damp_coeff_3d         = .false. 
logical :: use_this_module       = .false. 

! Adaptive restoring
real    :: athresh               = 0.5
real    :: npower                = 1.0
real    :: lambda                = 0.0083
real    :: taumin                = 720
logical :: use_adaptive_restore  = .false.
logical :: use_sponge_after_init = .false.
logical :: use_normalising       = .false.
logical :: use_hard_thump        = .false.
logical :: deflate               = .false.
real    :: deflate_fraction      = 0.6
integer :: secs_to_restore       = 0
integer :: days_to_restore       = 1

logical :: limit_temp            = .false.
real    :: limit_temp_min        = -1.8
real    :: limit_temp_restore    = 10800.

logical :: limit_salt            = .false.
real    :: limit_salt_min        = 0.01
real    :: limit_salt_restore    = 3600.

integer :: secs_restore
integer :: initial_day
integer :: initial_secs
real, allocatable :: sdiffo(:)

namelist /ocean_sponges_tracer_nml/ use_this_module, damp_coeff_3d
namelist /ocean_sponges_tracer_ofam_nml/ &
    use_adaptive_restore, use_sponge_after_init, use_normalising, &
    use_hard_thump, athresh, taumin, lambda, npower, days_to_restore, &
    secs_to_restore, deflate, deflate_fraction, &
    limit_temp, limit_temp_min, limit_temp_restore, &
    limit_salt, limit_salt_min, limit_salt_restore

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_sponges_tracer_init">
!
! <DESCRIPTION>
! This subroutine is intended to be used to initialize the tracer sponges.
! Everything in this subroutine is a user prototype, and should be replacable.
! </DESCRIPTION>
!
subroutine ocean_sponges_tracer_init(Grid, Domain, Time, T_prog, dtime, Ocean_options)

  type(ocean_grid_type),        intent(in), target :: Grid
  type(ocean_domain_type),      intent(in), target :: Domain
  type(ocean_time_type),        intent(in)         :: Time
  type(ocean_prog_tracer_type), intent(in)         :: T_prog(:)
  real,                         intent(in)         :: dtime
  type(ocean_options_type),     intent(inout)      :: Ocean_options

  integer :: i, j, k, n
  integer :: ioun, io_status, ierr
  integer :: index_temp
  integer :: secs, days
  real    :: dtimer

  character(len=128) :: name

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_sponges_tracer_mod (ocean_sponges_tracer_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  num_prog_tracers = size(T_prog(:))
  do n=1,num_prog_tracers
     if (trim(T_prog(n)%name) == 'temp') index_temp = n
  enddo

  allocate( Sponge(num_prog_tracers) )

  call write_version_number(version, tagname)

  ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
  read(input_nml_file, nml=ocean_sponges_tracer_nml, iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_sponges_tracer_nml')
  read(input_nml_file, nml=ocean_sponges_tracer_ofam_nml, iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_sponges_tracer_ofam_nml')
#else
  ioun =  open_namelist_file()
  read(ioun, ocean_sponges_tracer_nml, IOSTAT=io_status)
  ierr = check_nml_error(io_status, 'ocean_sponges_tracer_nml')
  rewind(ioun)
  read(ioun, ocean_sponges_tracer_ofam_nml, IOSTAT=io_status)
  ierr = check_nml_error(io_status, 'ocean_sponges_tracer_ofam_nml')
  call close_file(ioun)
#endif
  write(stdlogunit, ocean_sponges_tracer_nml)
  write(stdoutunit, '(/)')
  write(stdoutunit, ocean_sponges_tracer_nml)

  write(stdlogunit, ocean_sponges_tracer_ofam_nml)
  write(stdoutunit, '(/)')
  write(stdoutunit, ocean_sponges_tracer_ofam_nml)

  Dom => Domain
  Grd => Grid

#ifndef MOM_STATIC_ARRAYS    
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif

  do n=1,num_prog_tracers
     Sponge(n)%id = -1
  enddo

  dtimer = 1.0/dtime

  if(use_this_module) then
      write(stdoutunit,*)'==>Note from ocean_sponges_tracer_mod: Using this module.'
      Ocean_options%ocean_sponges_tracer= 'Used ocean tracer sponges.'
  else
      write(stdoutunit,*)'==>Note from ocean_sponges_tracer_mod: NOT using ocean tracer sponges.'
      Ocean_options%ocean_sponges_tracer= 'Did NOT use ocean tracer sponges.'
      return
  endif

  if (use_adaptive_restore) then
    allocate(sdiffo(num_prog_tracers))

    ! Set up some initial times
    call get_time(Time%model_time, secs, days)
    initial_day = days
    initial_secs = secs
  end if

  do n = 1, num_prog_tracers
    if (use_adaptive_restore) then
        allocate(Sponge(n)%damp_coeff(isd:ied,jsd:jed,nk))
        Sponge(n)%damp_coeff=0.0
    else
     ! read damping rates (1/sec)

     name = 'INPUT/'//trim(T_prog(n)%name)//'_sponge_coeff.nc'
     if (file_exist(trim(name))) then
        write(stdoutunit,*) '==> Using sponge damping times specified from file '//trim(name)
        allocate(Sponge(n)%damp_coeff(isd:ied,jsd:jed,nk))
        allocate(Sponge(n)%damp_coeff2d(isd:ied,jsd:jed))
        Sponge(n)%damp_coeff(:,:,:) = 0.0
        Sponge(n)%damp_coeff2d(:,:) = 0.0

        if(damp_coeff_3d) then
            call read_data(name,'coeff',Sponge(n)%damp_coeff,domain=Domain%domain2d,timelevel=1)
        else
            call read_data(name,'coeff',Sponge(n)%damp_coeff2d,domain=Domain%domain2d,timelevel=1)
            do k=1,nk
               do j=jsc,jec
                  do i=isc,iec
                     Sponge(n)%damp_coeff(i,j,k) = Sponge(n)%damp_coeff2d(i,j)
                  enddo
               enddo
            enddo
        endif

        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 if(Grd%tmask(i,j,k) == 0.0) then
                     Sponge(n)%damp_coeff(i,j,k) = 0.0
                 endif
              enddo
           enddo
        enddo

        ! modify damping rates to allow restoring to be solved implicitly
        ! note: test values between zero and 4.0e-8 revert to damping rates defined above

        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 if (dtime*Sponge(n)%damp_coeff(i,j,k) > 4.0e-8) then
                     if (dtime*Sponge(n)%damp_coeff(i,j,k) > 37.0) then
                         Sponge(n)%damp_coeff(i,j,k) = dtimer
                     else
                         Sponge(n)%damp_coeff(i,j,k) = (1.0 - exp(-dtime*Sponge(n)%damp_coeff(i,j,k))) * dtimer
                     endif
                 else if (dtime*Sponge(n)%damp_coeff(i,j,k) <= 0.0) then
                     Sponge(n)%damp_coeff(i,j,k) = 0.0
                 endif
              enddo
           enddo
        enddo
     else
        write(stdoutunit,*) '==> '//trim(name)//' not found.  Sponge not being applied '
        cycle
     endif
    endif

    ! Read restoring data

    name = 'INPUT/'//trim(T_prog(n)%name)//'_sponge.nc'
    Sponge(n)%id = init_external_field(name,T_prog(n)%name,domain=Domain%domain2d)
    if (Sponge(n)%id < 1) then
      call mpp_error(FATAL,&
      '==>Error: in ocean_sponges_tracer_mod: damping rates are specified but sponge values are not')
    endif
    write(stdoutunit,*) '==> Using sponge data specified from file '//trim(name)
  enddo

  ! register diagnostic output
  allocate (id_sponge_tend(num_prog_tracers))
  id_sponge_tend = -1

  do n=1,num_prog_tracers
     if(n==index_temp) then
        id_sponge_tend(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sponge_tend',&
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*cp*heating due to sponge',        &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e10,1.e10/))
     else
        id_sponge_tend(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sponge_tend',&
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*tendency due to sponge',          &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e10,1.e10/))
     endif
  enddo


end subroutine ocean_sponges_tracer_init
! </SUBROUTINE> NAME="ocean_sponges_tracer_init"


!#######################################################################
! <SUBROUTINE NAME="sponge_tracer_source">
!
! <DESCRIPTION>
! This subroutine calculates thickness weighted and density weighted
! time tendencies of tracers due to damping by sponges.
! </DESCRIPTION>
!
subroutine sponge_tracer_source(Time, Thickness, T_prog)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

  integer :: taum1, tau
  integer :: i, j, k, n

  real    :: sdiff, sum_val, adaptive_coeff
  integer :: secs, days, numsecs
  logical :: do_adaptive_restore = .false.  ! Only restore in specified time period.
  logical :: do_normal_restore = .false.  ! Only restore in specified time period.
  logical, save :: first_pass = .true.

  if(.not. use_this_module) return

  taum1 = Time%taum1
  tau   = Time%tau
  wrk1  = 0.0

  if (use_adaptive_restore) then
    secs_restore = days_to_restore * 86400 + secs_to_restore
    sdiff = 0.0
    call get_time(Time%model_time, secs, days)
    numsecs = (days - initial_day) * 86400 + secs - initial_secs

    do_adaptive_restore = (numsecs < secs_restore)
    do_normal_restore = (.not. do_adaptive_restore) .and. use_sponge_after_init
  else
    do_normal_restore = .true.
  endif

  do n=1,size(T_prog(:))

    ! Need to reinitialise wrk2 here due to limiting of temperature.
    wrk2  = 0.0

    if (Sponge(n)%id > 0) then

        ! get sponge value for current time
        call time_interp_external(Sponge(n)%id, Time%model_time, wrk1)
        if (do_adaptive_restore) then
           if (first_pass) then
              if (use_hard_thump) then
                 do k = 1, nk
                    do j = jsd, jed
                       do i = isd, ied
                          T_prog(n)%field(i,j,k,taum1) = wrk1(i,j,k)
                          T_prog(n)%field(i,j,k,tau) = wrk1(i,j,k)
                       end do
                    end do
                 end do
              end if
              wrk2 = abs(T_prog(n)%field(:,:,:,taum1) - wrk1(:,:,:)) &
                     * Grd%tmask(:,:,:)
              sum_val = sum(wrk2(isc:iec,jsc:jec,:))
              call mpp_sum(sum_val)
              sdiffo(n) = sum_val / Grd%wet_t_points
              if (.not. use_normalising) then
                sdiffo(n) = 1.0
              else
                sdiffo(n) = maxval(wrk2)
                call mpp_max(sdiffo(n))
                sdiffo(n) = 1.01 * sdiffo(n)
              endif
           endif

           do k = 1, nk
             do j = jsd, jed
               do i = isd, ied
                 sdiff = abs(T_prog(n)%field(i,j,k,taum1) - wrk1(i,j,k))
                 sdiff = max(sdiff, 1e-9) ! Minimum tolerance
                 adaptive_coeff = 1.0 / (taumin - &
                                         (lambda * real(secs_restore)) &
                                         / (((sdiff / sdiffo(n))**npower) & 
                                            * log(1 - athresh)))
                 wrk2(i,j,k) = Thickness%rho_dzt(i,j,k,tau) * adaptive_coeff &
                                  * (wrk1(i,j,k) - T_prog(n)%field(i,j,k,taum1))
               enddo
             enddo
           enddo
        else if (do_normal_restore) then
            do k = 1, nk
                do j = jsc, jec
                    do i = isc, iec
                        wrk2(i,j,k) = Thickness%rho_dzt(i,j,k,tau) &
                                     * Sponge(n)%damp_coeff(i,j,k) &
                                     * (wrk1(i,j,k) - T_prog(n)%field(i,j,k,taum1))
                    end do
                end do
             end do
        end if

    end if

    ! These tracer relaxation schemes were in the original OFAM and AusCOM
    ! source code, and may be needed by other OFAM-based experiments.
    !
    ! For now, let's keep them in the same place, but they are now turned off
    ! by default.  They may need to be moved or deleted in the future.

    ! Limit to -1.8 deg with 3hour restoring (got other limits here for initialising. Not sure if needed (namelist use?).
    if (limit_temp .and. trim(T_prog(n)%name) == 'temp') then
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                wrk2(i,j,k) = wrk2(i,j,k) + Thickness%rho_dzt(i,j,k,tau) &
                    * max(limit_temp_min - T_prog(n)%field(i,j,k,taum1), 0.0) &
                    / limit_temp_restore
!                wrk2(i,j,k) = wrk2(i,j,k)+Thickness%rho_dzt(i,j,k,tau)*min(43.0-T_prog(n)%field(i,j,k,taum1),0.0)/3600.0
             enddo
          enddo
       enddo
    end if

    if (limit_salt .and. trim(T_prog(n)%name) == 'salt') then
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                wrk2(i,j,k) = wrk2(i,j,k) + Thickness%rho_dzt(i,j,k,tau) &
                    * max(limit_salt_min - T_prog(n)%field(i,j,k,taum1), 0.0) &
                    / limit_salt_restore
             enddo
          enddo
       enddo
    end if

    ! Update tendency
    do k = 1, nk
        do j = jsc, jec
            do i = isc, iec
                T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) &
                                              + wrk2(i,j,k)
            end do
        end do
    end do

    if (id_sponge_tend(n) > 0) call diagnose_3d(Time, Grd, id_sponge_tend(n), &
         T_prog(n)%conversion*wrk2(:,:,:))

  enddo

  first_pass = .false.

  return

end subroutine sponge_tracer_source
! </SUBROUTINE> NAME="sponge_tracer_source"


end module ocean_sponges_tracer_mod

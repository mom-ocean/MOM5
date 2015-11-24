module ocean_sponges_velocity_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov.au"> Paul Sandery
!</CONTACT>
!
!<OVERVIEW>
! Thickness weighted velocity tendency [meter*meter/sec*sec] from sponges.
!</OVERVIEW>
!
!<DESCRIPTION>
!
! This module applies sponges to currents. The sponges
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
! the model grid of values towards which the currents are being driven.
!
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_sponges_velocity_nml">
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
use time_interp_external_mod, only: init_external_field, time_interp_external
use time_manager_mod,         only: time_type, set_date, get_time
use time_manager_mod,         only: operator( + ), operator( - ), operator( // )
use time_manager_mod,         only: operator( > ), operator( == ), operator( <= )

use ocean_domains_mod,        only: get_local_indices
use ocean_parameters_mod,     only: missing_value, rho0 
use ocean_types_mod,          only: ocean_domain_type, ocean_grid_type, ocean_thickness_type
use ocean_types_mod,          only: ocean_time_type, ocean_velocity_type, ocean_options_type
use ocean_workspace_mod,      only: wrk1, wrk2 
use ocean_util_mod,           only: diagnose_3d_u

implicit none

private

#include <ocean_memory.h>

type ocean_sponge_type
   integer :: id                                                 ! time_interp_external index
   character(len=32) :: name                                     ! tracer name corresponding to sponge
   real, dimension(:,:),   pointer :: damp_coeff_u2d   => NULL() ! 3d inverse damping rate (tracer units/ sec)
   real, dimension(:,:),   pointer :: damp_coeff_v2d   => NULL() ! 2d inverse damping rate (tracer units/ sec)
   real, dimension(:,:,:), pointer :: damp_coeff_u     => NULL() ! 3d inverse damping rate for u current (m/ sec^2)
   real, dimension(:,:,:), pointer :: damp_coeff_v     => NULL() ! 3d inverse damping rate for v current (m/ sec^2)
end type ocean_sponge_type

type(ocean_sponge_type), allocatable, dimension(:) :: Sponge_u
type(ocean_sponge_type), allocatable, dimension(:) :: Sponge_v

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()


public ocean_sponges_velocity_init
public sponge_velocity_source

character(len=126)  :: version = '$Id: ocean_sponges_velocity.F90,v 20.0 2013/12/14 00:16:26 fms Exp $'
character (len=128) :: tagname = '$Name: tikal $'

! for diagnostics 
logical :: used
integer, dimension(:), allocatable :: id_sponge_tend
logical :: module_is_initialized = .false.
logical :: damp_coeff_3d         = .false. 
logical :: use_this_module       = .false. 

! Adaptive restoring
real    :: athresh               = 0.5
real    :: npower                = 1.0
real    :: sdiffo_u              = 0.5
real    :: sdiffo_v              = 0.5
real    :: lambda                = 0.0083
real    :: taumin                = 720
logical :: use_adaptive_restore  = .false.
logical :: use_sponge_after_init = .false.
logical :: use_normalising       = .false.
logical :: use_hard_thump        = .false.
integer :: secs_to_restore       = 0
integer :: days_to_restore       = 1

integer :: secs_restore
integer :: days_end_restore
integer :: secs_end_restore
integer :: initial_day, initial_secs

namelist /ocean_sponges_velocity_nml/ use_this_module, damp_coeff_3d
namelist /ocean_sponges_velocity_ofam_nml/ &
    use_adaptive_restore, use_sponge_after_init, use_normalising, &
    use_hard_thump, athresh, taumin, lambda, npower, days_to_restore, &
    secs_to_restore

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_sponges_velocity_init">
!
! <DESCRIPTION>
! This subroutine is intended to be used to initialize the sponges.
! Everything in this subroutine is a user prototype, and should be replacable.
! </DESCRIPTION>
!
subroutine ocean_sponges_velocity_init(Grid, Domain, Time, dtime, Ocean_options)

  type(ocean_grid_type),          intent(in), target :: Grid
  type(ocean_domain_type),        intent(in), target :: Domain
  type(ocean_time_type),          intent(in)         :: Time
  real,                           intent(in)         :: dtime
  type(ocean_options_type),       intent(inout)      :: Ocean_options

  integer :: i, j, k
  integer :: ioun, io_status, ierr
  integer :: secs, days
  real    :: dtimer

  character(len=128) :: name

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_sponges_velocity_mod (ocean_sponges_velocity_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  allocate( Sponge_u(1) )
  allocate( Sponge_v(1) )

  call write_version_number(version, tagname)

  ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
  read(input_nml_file, nml=ocean_sponges_velocity_nml, iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_sponges_velocity_nml')
  read(input_nml_file, nml=ocean_sponges_velocity_ofam_nml, iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_sponges_velocity_ofam_nml')
#else
  ioun = open_namelist_file()
  read(ioun, ocean_sponges_velocity_nml, IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_sponges_velocity_nml')
  rewind(ioun)
  read(ioun, ocean_sponges_velocity_ofam_nml, IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_sponges_velocity_ofam_nml')
  call close_file(ioun)
#endif
  write(stdlogunit, ocean_sponges_velocity_nml)
  write(stdoutunit, '(/)')
  write(stdoutunit, ocean_sponges_velocity_nml)

  write(stdlogunit, ocean_sponges_velocity_ofam_nml)
  write(stdoutunit, '(/)')
  write(stdoutunit, ocean_sponges_velocity_ofam_nml)

  Dom => Domain
  Grd => Grid

#ifndef MOM_STATIC_ARRAYS    
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif

  Sponge_u(1)%id = -1
  Sponge_v(1)%id = -1
  dtimer = 1.0/dtime

  if(use_this_module) then
      write(stdoutunit,*)'==>Note from ocean_sponges_velocity_mod: Using this module.'
      Ocean_options%ocean_sponges_velocity= 'Used ocean velocity sponges.'
  else
      write(stdoutunit,*)'==>Note from ocean_sponges_velocity_mod: NOT using this module: no velocity sponges.'
      Ocean_options%ocean_sponges_velocity= 'Did NOT use ocean velocity sponges.'
      return
  endif

  if (use_adaptive_restore) then
    ! Set up some initial times
    call get_time(Time%model_time, secs, days)
    days_end_restore = days + days_to_restore
    secs_end_restore = secs + secs_to_restore
    initial_day = days
    initial_secs = secs

    allocate(Sponge_u(1)%damp_coeff_u(isd:ied, jsd:jed, nk))
    Sponge_u(1)%damp_coeff_u(:,:,:) = 0.0

    allocate(Sponge_v(1)%damp_coeff_v(isd:ied,jsd:jed,nk))
    Sponge_v(1)%damp_coeff_v(:,:,:) = 0.0
  else
    ! Read damping rates for u
    name = 'INPUT/u_sponge_coeff.nc'
    if (file_exist(name)) then
        write(stdoutunit,*) '==> Using sponge damping times specified from file ',name
        allocate(Sponge_u(1)%damp_coeff_u(isd:ied,jsd:jed,nk))
        allocate(Sponge_u(1)%damp_coeff_u2d(isd:ied,jsd:jed))
        Sponge_u(1)%damp_coeff_u(:,:,:) = 0.0
        Sponge_u(1)%damp_coeff_u2d(:,:) = 0.0

        if(damp_coeff_3d) then
            call read_data(name,'coeff',Sponge_u(1)%damp_coeff_u,domain=Domain%domain2d,timelevel=1)
        else
            call read_data(name,'coeff',Sponge_u(1)%damp_coeff_u2d,domain=Domain%domain2d,timelevel=1)
            do k=1,nk
               do j=jsc,jec
                  do i=isc,iec
                     Sponge_u(1)%damp_coeff_u(i,j,k) = Sponge_u(1)%damp_coeff_u2d(i,j)
                  enddo
               enddo
            enddo
        endif

        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 if(Grd%umask(i,j,k) == 0.0) then
                     Sponge_u(1)%damp_coeff_u(i,j,k) = 0.0
                 endif
              enddo
           enddo
        enddo


       ! modify damping rates to allow restoring to be solved implicitly
       ! note: test values between zero and 4.0e-8 revert to damping rates defined above
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 if (dtime*Sponge_u(1)%damp_coeff_u(i,j,k) > 4.0e-8) then
                     if (dtime*Sponge_u(1)%damp_coeff_u(i,j,k) > 37.0) then
                         Sponge_u(1)%damp_coeff_u(i,j,k) = dtimer
                     else
                         Sponge_u(1)%damp_coeff_u(i,j,k) = (1.0 - exp(-dtime*Sponge_u(1)%damp_coeff_u(i,j,k))) * dtimer
                     endif
                 else if (dtime*Sponge_u(1)%damp_coeff_u(i,j,k) <= 0.0) then
                     Sponge_u(1)%damp_coeff_u(i,j,k) = 0.0
                 endif
              enddo
           enddo
        enddo
    endif  ! endif for if fileexist INPUT/u_sponge_coeff.nc

    ! Read damping rates for v-sponge
    name = 'INPUT/v_sponge_coeff.nc'
    if (file_exist(name)) then
        write(stdoutunit,*) '==> Using sponge damping times specified from file ',name
        allocate(Sponge_v(1)%damp_coeff_v(isd:ied,jsd:jed,nk))
        allocate(Sponge_v(1)%damp_coeff_v2d(isd:ied,jsd:jed))
        Sponge_v(1)%damp_coeff_v(:,:,:) = 0.0
        Sponge_v(1)%damp_coeff_v2d(:,:) = 0.0

        if(damp_coeff_3d) then
            call read_data(name,'coeff',Sponge_v(1)%damp_coeff_v,domain=Domain%domain2d,timelevel=1)
        else
            call read_data(name,'coeff',Sponge_v(1)%damp_coeff_v2d,domain=Domain%domain2d,timelevel=1)
            do k=1,nk
               do j=jsc,jec
                  do i=isc,iec
                     Sponge_v(1)%damp_coeff_v(i,j,k) = Sponge_v(1)%damp_coeff_v2d(i,j)
                  enddo
               enddo
            enddo
        endif

        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 if(Grd%umask(i,j,k) == 0.0) then
                     Sponge_v(1)%damp_coeff_v(i,j,k) = 0.0
                 endif
              enddo
           enddo
        enddo

        ! modify damping rates to allow restoring to be solved implicitly
        ! note: test values between zero and 4.0e-8 revert to damping rates defined above
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 if (dtime*Sponge_v(1)%damp_coeff_v(i,j,k) > 4.0e-8) then
                     if (dtime*Sponge_v(1)%damp_coeff_v(i,j,k) > 37.0) then
                         Sponge_v(1)%damp_coeff_v(i,j,k) = dtimer
                     else
                         Sponge_v(1)%damp_coeff_v(i,j,k) = (1.0 - exp(-dtime*Sponge_v(1)%damp_coeff_v(i,j,k))) * dtimer
                     endif
                 else if (dtime*Sponge_v(1)%damp_coeff_v(i,j,k) <= 0.0) then
                     Sponge_v(1)%damp_coeff_v(i,j,k) = 0.0
                 endif
              enddo
           enddo
        enddo

    endif   ! endif for if fileexist INPUT/v_sponge_coeff.nc
  end if

  ! read forcing data for u-sponge
  name = 'INPUT/u_sponge.nc'
  if (file_exist(trim(name)) ) then
      Sponge_u(1)%id = init_external_field(name,'u',domain=Domain%domain2d)
      if (Sponge_u(1)%id < 1) then
         if ( use_adaptive_restore ) then
            call mpp_error(FATAL,&
               '==>Error: in ocean_sponges_velocity_mod: adaptive forcing is specified but restoring values are not')
         else
            call mpp_error(FATAL,&
               '==>Error: in ocean_sponges_velocity_mod: forcing rates are specified but restoring values are not')
         endif
      endif
      write(stdoutunit,*) '==> Using restoring data specified from file '//trim(name)
  else
      write(stdoutunit,*) '==> '//trim(name)//' not found.  Increment not being applied '
  endif


  ! read forcing data for v-sponge
  name = 'INPUT/v_sponge.nc'
  if (file_exist(trim(name)) ) then
      Sponge_v(1)%id = init_external_field(name,'v',domain=Domain%domain2d)
      if (Sponge_v(1)%id < 1) then
         if ( use_adaptive_restore ) then
            call mpp_error(FATAL,&
               '==>Error: in ocean_sponges_velocity_mod: adaptive forcing is specified but restoring values are not')
         else
            call mpp_error(FATAL,&
               '==>Error: in ocean_sponges_velocity_mod: forcing rates are specified but restoring values are not')
         endif
      endif
      write(stdoutunit,*) '==> Using restoring data specified from file '//trim(name)
  else
      write(stdoutunit,*) '==> '//trim(name)//' not found.  Increment not being applied '
  endif


  ! register diagnostic output
  allocate (id_sponge_tend(2))
  id_sponge_tend = -1

  id_sponge_tend(1) = register_diag_field ('ocean_model', 'u_sponge_tend',       &
       Grd%vel_axes_uv(1:3), Time%model_time, 'rho*dzt*u_tendency due to sponge',&
       '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1.e10,1.e10/))
  id_sponge_tend(2) = register_diag_field ('ocean_model', 'v_sponge_tend',       &
       Grd%vel_axes_uv(1:3), Time%model_time, 'rho*dzt*v_tendency due to sponge',&
       '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1.e10,1.e10/))


end subroutine ocean_sponges_velocity_init
! </SUBROUTINE> NAME="ocean_sponges_velocity_init"



!#######################################################################
! <SUBROUTINE NAME="sponge_velocity_source">
!
! <DESCRIPTION>
! This subroutine calculates thickness weighted and density weighted
! time tendencies due to damping by sponges.
! </DESCRIPTION>
!
subroutine sponge_velocity_source(Time, Thickness, Velocity)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_velocity_type),    intent(inout) :: Velocity

  integer :: secs, days
  integer :: taum1, tau
  integer :: i, j, k

  real    :: sdiff, sum_val, adaptive_coeff
  integer :: numsecs

  logical :: do_adaptive_restore = .false.  ! Only restore in specified time period.
  logical :: do_normal_restore = .false.    !
  logical, save :: first_pass = .true.

  if (.not. use_this_module) return

  taum1 = Time%taum1
  tau   = Time%tau

  if (use_adaptive_restore) then
    secs_restore = days_to_restore * 86400 + secs_to_restore
    sdiff = 0.0
    call get_time(Time%model_time, secs, days)
    numsecs = (days - initial_day) * 86400 + secs - initial_secs

    do_adaptive_restore = (numsecs < secs_restore)
    do_normal_restore = (.not. do_adaptive_restore) .and. use_sponge_after_init

    write(stdout(),*) 'dar,dnr', do_adaptive_restore, do_normal_restore
  else
     do_normal_restore = .true.
  endif

  wrk1 = 0.0
  wrk2 = 0.0
  if (Sponge_u(1)%id > 0) then

      call time_interp_external(Sponge_u(1)%id, Time%model_time, wrk1)

      if ( do_adaptive_restore ) then
         if ( first_pass ) then
            if ( use_hard_thump ) then
               do k = 1,nk
                  do j=jsd,jed
                     do i=isd,ied
                        Velocity%u(i,j,k,1,taum1)=wrk1(i,j,k)
                        Velocity%u(i,j,k,1,tau)=wrk1(i,j,k)
                     enddo
                  enddo
               enddo
            endif
            wrk2=abs(Velocity%u(:,:,:,1,taum1)-wrk1(:,:,:))*Grd%umask(:,:,:)
            sum_val=sum(wrk2(isc:iec,jsc:jec,:))

            call mpp_sum(sum_val)
            sdiffo_u = sum_val/Grd%wet_u_points
            !if ( mpp_pe() == mpp_root_pe() ) write (stdout(),*)'mean U_diff', sdiffo
            write (stdout(),*)'mean U_diff', sdiffo_u
            if ( .not. use_normalising ) then
               sdiffo_u=1.0
            else
                ! Work out max difference + threshhold
                sdiffo_u  = maxval(wrk2)
                call mpp_max(sdiffo_u)
                sdiffo_u  = 1.01*sdiffo_u
             endif
          endif

         ! Calculate idealised forcing tendency timescale for each timestep and array element
         do k=1,nk
            do j=jsd,jed
               do i=isd,ied
                  sdiff=abs(Velocity%u(i,j,k,1,taum1)-wrk1(i,j,k))
                  sdiff = max(sdiff,1e-9) !minimum tolerance
                  adaptive_coeff=1.0/(taumin - (lambda*real(secs_restore))/(((sdiff/sdiffo_v)**npower)*log(1-athresh)))
                  wrk2(i,j,k) = adaptive_coeff*(wrk1(i,j,k) - Velocity%u(i,j,k,1,taum1))
                  Velocity%accel(i,j,k,1) = Velocity%accel(i,j,k,1) + Thickness%rho_dzu(i,j,k,tau)*wrk2(i,j,k)
               enddo
            enddo
         enddo

      elseif ( do_normal_restore ) then
         do k=1,nk
            do j=jsd,jed
               do i=isd,ied
                  wrk2(i,j,k) = Sponge_u(1)%damp_coeff_u(i,j,k)*(wrk1(i,j,k) - &
                                Velocity%u(i,j,k,1,taum1))
                  Velocity%accel(i,j,k,1) = Velocity%accel(i,j,k,1) + &
                                            Thickness%rho_dzu(i,j,k,tau)*wrk2(i,j,k)
               enddo
            enddo
         enddo
      endif

  endif

  call diagnose_3d_u(Time, Grd, id_sponge_tend(1), wrk2(:,:,:))

  wrk1 = 0.0
  wrk2 = 0.0
  if (Sponge_v(1)%id > 0) then

      call time_interp_external(Sponge_v(1)%id, Time%model_time, wrk1)

      if ( do_adaptive_restore ) then
         if ( first_pass ) then
            if ( use_hard_thump ) then
               do k = 1,nk
                  do j=jsd,jed
                     do i=isd,ied
                        Velocity%u(i,j,k,2,taum1)=wrk1(i,j,k)
                        Velocity%u(i,j,k,2,tau)=wrk1(i,j,k)
                     enddo
                  enddo
               enddo
            endif
            wrk2=abs(Velocity%u(:,:,:,2,taum1)-wrk1(:,:,:))*Grd%umask(:,:,:)
            sum_val=sum(wrk2(isc:iec,jsc:jec,:))

            call mpp_sum(sum_val)
            sdiffo_v = sum_val/Grd%wet_u_points

            if ( .not. use_normalising ) then
               sdiffo_v=1.0
            else
                sdiffo_v  = maxval(wrk2)
                call mpp_max(sdiffo_v)
                sdiffo_v  = 1.01*sdiffo_v
             endif
          endif


         do k=1,nk
            do j=jsd,jed
               do i=isd,ied
                  sdiff=abs(Velocity%u(i,j,k,2,taum1)-wrk1(i,j,k))
                  sdiff = max(sdiff,1e-9) !minimum tolerance
                  adaptive_coeff=1.0/(taumin - (lambda*real(secs_restore))/(((sdiff/sdiffo_v)**npower)*log(1-athresh)))
                  wrk2(i,j,k) = adaptive_coeff*(wrk1(i,j,k) - Velocity%u(i,j,k,2,taum1))
                  Velocity%accel(i,j,k,2) = Velocity%accel(i,j,k,2) + Thickness%rho_dzu(i,j,k,tau)*wrk2(i,j,k)
               enddo
            enddo
         enddo
      endif

      if ( do_normal_restore ) then
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk2(i,j,k) = Sponge_v(1)%damp_coeff_v(i,j,k)*(wrk1(i,j,k) - &
                                Velocity%u(i,j,k,2,taum1))
                  Velocity%accel(i,j,k,2) = Velocity%accel(i,j,k,2) + &
                                            Thickness%rho_dzu(i,j,k,tau)*wrk2(i,j,k)
               enddo
            enddo
         enddo
      endif
  endif

  call diagnose_3d_u(Time, Grd, id_sponge_tend(2), wrk2(:,:,:))

  first_pass = .false.

  return

end subroutine sponge_velocity_source
! </SUBROUTINE> NAME="sponge_velocity_source"

end module ocean_sponges_velocity_mod

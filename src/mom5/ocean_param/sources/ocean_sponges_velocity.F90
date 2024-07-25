!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ocean_sponges_velocity_mod
  !
  !<CONTACT EMAIL="p.sandery@bom.gov.au"> Paul Sandery
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
  ! We now allow for the possibility of restoring adaptively according to 
  ! Sandery et al. (2011). The user  may specify parameters as defined in 
  ! in that paper. 
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
  use diag_manager_mod,         only: register_diag_field, send_data
  use fms_mod,                  only: write_version_number, open_namelist_file, close_file
  use fms_mod,                  only: file_exist
  use fms_mod,                  only: open_namelist_file, check_nml_error, close_file
  use fms_mod,                  only: read_data, lowercase, FATAL, WARNING, stdout, stdlog
  use mpp_mod,                  only: mpp_sum, mpp_error, mpp_max
  use time_interp_external_mod, only: init_external_field, time_interp_external
  use time_manager_mod,         only: time_type, set_date, get_time
  use time_manager_mod,         only: operator( + ), operator( - ), operator( // )
  use time_manager_mod,         only: operator( > ), operator( == ), operator( <= )
  use mpp_domains_mod,          only: mpp_update_domains
  use ocean_domains_mod,        only: get_local_indices
  use ocean_parameters_mod,     only: missing_value, rho0 
  use ocean_types_mod,          only: ocean_domain_type, ocean_grid_type, ocean_thickness_type
  use ocean_types_mod,          only: ocean_time_type, ocean_velocity_type, ocean_options_type
  use ocean_workspace_mod,      only: wrk1, wrk2 

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

  character(len=126)  :: version = '$Id: ocean_sponges_velocity.F90,v 1.1.2.3.42.1 2009/10/10 00:42:54 nnz Exp $'
  character (len=128) :: tagname = '$Name: mom4p1_pubrel_dec2009_nnz $'

  ! for diagnostics 
  logical :: used
  integer, dimension(:), allocatable :: id_sponge_tend
  logical :: module_is_initialized = .false.
  logical :: damp_coeff_3d         = .false. 
  logical :: use_this_module       = .false. 

  ! Adaptive restoring

  real    :: athresh               = 0.5
  real    :: npower                = 1.0
  real    :: sdiffo_u                = 0.5
  real    :: sdiffo_v                = 0.5
  real    :: lambda                = 0.0083
  real    :: taumin                = 3600
  logical :: use_adaptive_restore  = .false.
  logical :: use_sponge_after_init = .false.
  logical :: use_normalising         = .false.
  logical :: use_hard_thump        = .false.
  logical :: use_increment          = .false.
  integer :: secs_to_restore        = 0
  integer :: days_to_restore        = 0
  integer :: secs_restore
  integer :: days_end_restore
  integer :: secs_end_restore
  integer :: initial_day, initial_secs
  logical :: limit_u = .false.
  logical :: limit_v = .false.
  real    :: limit_u_max      = 3.0 
  real    :: limit_u_min      = -3.0
  real    :: limit_v_max      = 3.0
  real    :: limit_v_min      = -3.0

  namelist /ocean_sponges_velocity_nml/ use_this_module, damp_coeff_3d
  namelist /ocean_sponges_velocity_OFAM_nml/ use_adaptive_restore,use_sponge_after_init, use_normalising,use_hard_thump,&
       athresh, taumin,lambda,npower,days_to_restore, secs_to_restore, limit_u, limit_v

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

    integer :: i, j, k, n
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

    call write_version_number()

    ! provide for namelist over-ride of default values
    ioun =  open_namelist_file()
    read (ioun,ocean_sponges_velocity_nml,IOSTAT=io_status)
    write (stdlogunit,ocean_sponges_velocity_nml)
    write (stdoutunit,'(/)') 
    write (stdoutunit,ocean_sponges_velocity_nml)
    ierr = check_nml_error(io_status,'ocean_sponges_velocity_nml')
    call close_file (ioun)

    ! provide for namelist over-ride of default values
    ioun =  open_namelist_file()
    read (ioun,ocean_sponges_velocity_OFAM_nml,IOSTAT=io_status)
    write (stdlogunit,ocean_sponges_velocity_OFAM_nml)
    write (stdoutunit,'(/)') 
    write (stdoutunit,ocean_sponges_velocity_OFAM_nml)
    ierr = check_nml_error(io_status,'ocean_sponges_velocity_OFAM_nml')
    call close_file (ioun)

    Dom => Domain
    Grd => Grid

#ifndef MOM4_STATIC_ARRAYS    
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

    !set up some initial times
    call get_time(Time%model_time,secs,days)
    days_end_restore=days + days_to_restore
    secs_end_restore=secs + secs_to_restore
    initial_day = days
    initial_secs = secs

    !read damping rates for u
    if (.not. use_hard_thump  .and. .not. use_adaptive_restore) then
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
          endif !damp_coeff_3d
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
       endif !if fileexist INPUT/u_sponge_coeff.nc
    endif !.not. use_hard_thump .and. .not. use_adaptive_restore

    !read damping rates for v-sponge 
    if (.not. use_hard_thump  .and. .not. use_adaptive_restore) then
       allocate(Sponge_v(1)%damp_coeff_v(isd:ied,jsd:jed,nk))
       Sponge_v(1)%damp_coeff_v(:,:,:) = 0.0
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
          endif! damp_coeff_3d
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
       endif  !fileexist INPUT/v_sponge_coeff.nc
    endif !.not. use_hard_thump .and. .not. use_adaptive_restore

    ! read forcing data for u-sponge 
    name = 'INPUT/u_sponge.nc'
    if (file_exist(trim(name)) ) then
       Sponge_u(1)%id = init_external_field(name,'u',domain=Domain%domain2d)
       if (Sponge_u(1)%id < 1) then
          call mpp_error(FATAL,&
               '==>Error: in ocean_sponges_velocity_mod: Sponge values for u are required but not specified')
       endif
       write(stdoutunit,*) '==> Using restoring data specified from file '//trim(name)
    else
       write(stdoutunit,*) '==> '//trim(name)//' not found. Sponge not being applied '
    endif

    ! read forcing data for v-sponge 
    name = 'INPUT/v_sponge.nc'
    if (file_exist(trim(name)) ) then
       Sponge_v(1)%id = init_external_field(name,'v',domain=Domain%domain2d)
       if (Sponge_v(1)%id < 1) then
          call mpp_error(FATAL,&
               '==>Error: in ocean_sponges_velocity_mod: Sponge values for v are required but not specified')
       endif
       write(stdoutunit,*) '==> Using restoring data specified from file '//trim(name)
    else
       write(stdoutunit,*) '==> '//trim(name)//' not found. Sponge not being applied '
    endif

    ! register diagnostic output
    allocate (id_sponge_tend(2))
    id_sponge_tend = -1

    id_sponge_tend(1) = register_diag_field ('ocean_model', 'u_sponge_tend',       &
         Grd%vel_axes_uv(1:3), Time%model_time, 'u_tendency due to sponge',&
         '(m/s^2)', missing_value=missing_value, range=(/-1.e-2,1.e-2/))
    id_sponge_tend(2) = register_diag_field ('ocean_model', 'v_sponge_tend',       &
         Grd%vel_axes_uv(1:3), Time%model_time, 'v_tendency due to sponge',&
         '(m/s^2)', missing_value=missing_value, range=(/-1.e-2,1.e-2/))


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
    integer :: i, j, k, n

    real    :: sdiff, sum_val, adaptive_coeff
    integer :: numsecs

    logical :: do_adaptive_restore = .false.  ! Only restore in specified time period.
    logical,save :: first_pass = .true.

    if (.not. use_this_module) return 

    taum1 = Time%taum1
    tau   = Time%tau

    wrk1 = 0.0
    wrk2 = 0.0
    sdiff=0.0

    call get_time(Time%model_time,secs,days)
    !-----------------------------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------------------------
    if (Sponge_u(1)%id > 0) then
       call time_interp_external(Sponge_u(1)%id, Time%model_time, wrk1)
       ! filter
       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                if ( wrk1(i,j,k).ge.limit_u_max ) then
                   wrk1(i,j,k)=limit_u_max
                endif
                if ( wrk1(i,j,k).le.limit_u_min ) then
                   wrk1(i,j,k)=limit_u_min
                endif
             enddo
          enddo
       enddo
       !-----------------------------------------------------------------------------------------------------------------
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
             call mpp_update_domains (Velocity%u(:,:,:,1,taum1) ,Dom%domain2d) 
             call mpp_update_domains (Velocity%u(:,:,:,1,tau) ,Dom%domain2d) 
             wrk2=abs(Velocity%u(:,:,:,1,taum1)-wrk1(:,:,:))*Grd%umask(:,:,:)
          endif !use_hard_thump
          sum_val=sum(wrk2(isc:iec,jsc:jec,:))
          call mpp_sum(sum_val)
          sdiffo_u = sum_val/Grd%wet_u_points
          write (stdout(),*)'mean absolute deviation u = ', sdiffo_u
          if ( .not. use_normalising ) then
             sdiffo_u=1.0
          else
             !work out max difference + threshhold
             sdiffo_u  = maxval(wrk2)
             call mpp_max(sdiffo_u)
             sdiffo_u  = 1.01*sdiffo_u
          endif
       endif !first_pass
       !--------------------------------------------------------------------------------------------------------------
       if (use_hard_thump == .false.) then
          if (use_adaptive_restore) then
             secs_restore = days_to_restore*86400 + secs_to_restore
             numsecs = (days-initial_day)*86400 + secs - initial_secs
             do_adaptive_restore = ( numsecs < secs_restore )
          endif
          if ( do_adaptive_restore ) then
             !calculate idealised forcing tendency timescale for each timestep and array element
             do k=1,nk
                do j=jsd,jed
                   do i=isd,ied
                      sdiff=abs(Velocity%u(i,j,k,1,taum1)-wrk1(i,j,k))
                      sdiff = max(sdiff,1e-9) !minimum tolerance
                      adaptive_coeff=1.0/(taumin - (lambda*real(secs_restore))/(((sdiff/sdiffo_u)**npower)*log(1-athresh)))
                      wrk2(i,j,k) = adaptive_coeff*(wrk1(i,j,k) - Velocity%u(i,j,k,1,taum1))
                      Velocity%accel(i,j,k,1) = Velocity%accel(i,j,k,1) + Thickness%rho_dzu(i,j,k,tau)*wrk2(i,j,k)
                   enddo
                enddo
             enddo
          else !normal restore
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
          endif !do_adaptive_restore
       endif !use_hard_thump
       !-----------------------------------------------------------------------------------------------------------------
    endif !Sponge_u
    !--------------------------------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------------------
    if (id_sponge_tend(1) > 0) then 
       used = send_data(id_sponge_tend(1),                       &
            wrk2(:,:,:), Time%model_time, rmask=Grd%umask(:,:,:),&
            is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif
    !--------------------------------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------------------
    wrk1 = 0.0
    wrk2 = 0.0
    if (Sponge_v(1)%id > 0) then
       call time_interp_external(Sponge_v(1)%id, Time%model_time, wrk1)
       ! filter
       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                if ( wrk1(i,j,k).ge.limit_v_max ) then
                   wrk1(i,j,k)=limit_v_max
                endif
                if ( wrk1(i,j,k).le.limit_v_min ) then
                   wrk1(i,j,k)=limit_v_min
                endif
             enddo
          enddo
       enddo
       !--------------------------------------------------------------------------------------------------------------
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
             call mpp_update_domains (Velocity%u(:,:,:,2,taum1) ,Dom%domain2d) 
             call mpp_update_domains (Velocity%u(:,:,:,2,tau) ,Dom%domain2d) 
             wrk2=abs(Velocity%u(:,:,:,2,taum1)-wrk1(:,:,:))*Grd%umask(:,:,:)
          endif
          sum_val=sum(wrk2(isc:iec,jsc:jec,:))
          call mpp_sum(sum_val)
          sdiffo_v = sum_val/Grd%wet_u_points
          write (stdout(),*)'mean absolute deviation v = ', sdiffo_v
          if ( .not. use_normalising ) then
             sdiffo_v=1.0
          else
             sdiffo_v  = maxval(wrk2)
             call mpp_max(sdiffo_v)
             sdiffo_v  = 1.01*sdiffo_v
          endif
       endif !first_pass
       !--------------------------------------------------------------------------------------------------------------
       if (use_hard_thump == .false.) then
          if (use_adaptive_restore) then
             secs_restore = days_to_restore*86400 + secs_to_restore
             numsecs = (days-initial_day)*86400 + secs - initial_secs
             do_adaptive_restore = ( numsecs < secs_restore )
          endif
          if ( do_adaptive_restore ) then
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
          else !normal restore
             do k=1,nk
                do j=jsd,jed
                   do i=isd,ied
                      wrk2(i,j,k) = Sponge_v(1)%damp_coeff_v(i,j,k)*(wrk1(i,j,k) - &
                           Velocity%u(i,j,k,2,taum1))
                      Velocity%accel(i,j,k,2) = Velocity%accel(i,j,k,2) + &
                           Thickness%rho_dzu(i,j,k,tau)*wrk2(i,j,k)
                   enddo
                enddo
             enddo
          endif !do_adaptive_restore
       endif !use_hard_thump
       !-----------------------------------------------------------------------------------------------------------------
    endif !Sponge_v
    !--------------------------------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------------------
    if (id_sponge_tend(2) > 0) then 
       used = send_data(id_sponge_tend(2),                       &
            wrk2(:,:,:), Time%model_time, rmask=Grd%umask(:,:,:),&
            is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif
    !--------------------------------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------------------
    first_pass = .false.

    return

  end subroutine sponge_velocity_source
  ! </SUBROUTINE> NAME="sponge_velocity_source"
end module ocean_sponges_velocity_mod

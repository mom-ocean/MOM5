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
module ocean_sponges_tracer_mod
  !
  !<CONTACT EMAIL="Bonnie.Samuels@noaa.gov"> Bonnie Samuels 
  !</CONTACT>
  !
  !<CONTACT EMAIL="Robert.Hallberg@noaa.gov"> R.W. Hallberg 
  !</CONTACT>
  !
  !<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> M.J. Harrison 
  !</CONTACT>
  !
  !<CONTACT EMAIL="swathi@cmmacs.ernet.in"> P. S. Swathi  
  !</CONTACT>
  !
  !<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R. C. Pacanowski
  !</CONTACT>
  !
  !<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S. M. Griffies 
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
  !
  ! OFAM MODIFICATIONS:
  !
  ! Adaptive nudiging is permitted.
  ! We can limit the length of time for which nudging is applied.
  ! There is an option to restore to a climatology in the bottom 4 layers
  ! We can limit the temp and salt to fix values.

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
  use diag_manager_mod,         only: register_diag_field, send_data
  use fms_mod,                  only: write_version_number, open_namelist_file, close_file
  use fms_mod,                  only: file_exist
  use fms_mod,                  only: open_namelist_file, check_nml_error, close_file
  use fms_mod,                  only: read_data, lowercase, FATAL, WARNING, stdout, stdlog
  use mpp_mod,                  only: mpp_sum, mpp_error, mpp_max, mpp_min,mpp_pe
  use time_interp_external_mod, only: init_external_field, time_interp_external
  use time_manager_mod,         only: time_type, set_date, get_time
  use time_manager_mod,         only: operator( + ), operator( - ), operator( // )
  use time_manager_mod,         only: operator( > ), operator( == ), operator( <= )
  use mpp_domains_mod,          only: mpp_update_domains 
  use ocean_domains_mod,        only: get_local_indices
  use ocean_parameters_mod,     only: missing_value 
  use ocean_types_mod,          only: ocean_domain_type, ocean_grid_type, ocean_thickness_type
  use ocean_types_mod,          only: ocean_prog_tracer_type, ocean_options_type, ocean_time_type 
  use ocean_workspace_mod,      only: wrk1, wrk2

  implicit none

  private

#include <ocean_memory.h>

  type ocean_sponge_type
     integer :: id                                             ! time_interp_external index
     character(len=32) :: name                                 ! tracer name corresponding to sponge
     real, dimension(:,:,:), pointer :: damp_coeff   => NULL() ! 3d inverse damping rate (tracer units/ sec)
     real, dimension(:,:)  , pointer :: damp_coeff2d => NULL() ! 2d inverse damping rate (tracer units/ sec)
  end type ocean_sponge_type

  type ocean_cars_type
     integer :: id                                             ! time_interp_external index
     character(len=32) :: name                                 ! tracer name corresponding to sponge
     real, dimension(:,:,:)  , pointer :: field => NULL()        ! field to be damped
  end type ocean_cars_type



  type(ocean_sponge_type), allocatable, dimension(:) :: Sponge
  type(ocean_cars_type), allocatable, dimension(:) :: Cars
  type(ocean_domain_type), pointer :: Dom => NULL()
  type(ocean_grid_type),   pointer :: Grd => NULL()

  public ocean_sponges_tracer_init
  public sponge_tracer_source

  character(len=126)  :: version = '$Id: ocean_sponges_tracer.F90,v 1.1.2.4.42.1 2009/10/10 00:42:53 nnz Exp $'
  character (len=128) :: tagname = '$Name: mom4p1_pubrel_dec2009_nnz $'

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
  real    :: taumin                = 3600
  logical :: use_adaptive_restore  = .false.
  logical :: use_sponge_after_init = .false.
  logical :: use_normalising       = .false.
  logical :: use_hard_thump        = .false.
  logical :: use_increment         = .false.
  integer :: secs_to_restore       = 0
  integer :: days_to_restore       = 0
  integer :: secs_restore
  integer :: initial_day
  integer :: initial_secs
  real,allocatable :: sdiffo(:)
  logical :: restore_cars_deep = .false.
  logical :: limit_temp = .false.
  logical :: limit_salt = .false.
  real    :: limit_temp_max      = 40.0
  real    :: limit_temp_min      = -1.8
  real    :: limit_salt_max      = 40.0
  real    :: limit_salt_min      = 0.1

  namelist /ocean_sponges_tracer_nml/ use_this_module, damp_coeff_3d
  namelist /ocean_sponges_tracer_OFAM_nml/ use_adaptive_restore,use_sponge_after_init, use_normalising, use_hard_thump, &
       athresh, taumin,lambda,npower,days_to_restore, secs_to_restore, &
       restore_cars_deep, limit_temp, limit_salt, &
       limit_temp_max, limit_temp_min, limit_salt_max, limit_salt_min

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

    type(ocean_grid_type),            intent(in), target :: Grid
    type(ocean_domain_type),       intent(in), target :: Domain
    type(ocean_time_type),           intent(in)           :: Time
    type(ocean_prog_tracer_type), intent(in)           :: T_prog(:)
    real,                                      intent(in)           :: dtime
    type(ocean_options_type),       intent(inout)       :: Ocean_options

    integer :: i, j, k, n
    integer :: ioun, io_status, ierr
    integer :: index_temp
    integer :: secs, days
    real     :: dtimer

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
    allocate( Cars(num_prog_tracers) )

    call write_version_number()

    ! provide for namelist over-ride of default values
    ioun =  open_namelist_file()
    read (ioun,ocean_sponges_tracer_nml,IOSTAT=io_status)
    write (stdlogunit,ocean_sponges_tracer_nml)
    write (stdoutunit,'(/)') 
    write (stdoutunit,ocean_sponges_tracer_nml)
    ierr = check_nml_error(io_status,'ocean_sponges_tracer_nml')
    call close_file (ioun)

    ! provide for namelist over-ride of default values
    ioun =  open_namelist_file()
    read (ioun,ocean_sponges_tracer_OFAM_nml,IOSTAT=io_status)
    write (stdlogunit,ocean_sponges_tracer_OFAM_nml)
    write (stdoutunit,'(/)') 
    write (stdoutunit,ocean_sponges_tracer_OFAM_nml)
    ierr = check_nml_error(io_status,'ocean_sponges_tracer_OFAM_nml')
    call close_file (ioun)

    Dom => Domain
    Grd => Grid

#ifndef MOM4_STATIC_ARRAYS    
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grid%nk
#endif

    do n=1,num_prog_tracers
       Sponge(n)%id = -1
    enddo
    allocate(sdiffo(num_prog_tracers))

    dtimer = 1.0/dtime

    if(use_this_module) then
       write(stdoutunit,*)'==>Note from ocean_sponges_tracer_mod: Using this module.'
       Ocean_options%ocean_sponges_tracer= 'Used ocean tracer sponges.'
    else
       write(stdoutunit,*)'==>Note from ocean_sponges_tracer_mod: NOT using ocean tracer sponges.'
       Ocean_options%ocean_sponges_tracer= 'Did NOT use ocean tracer sponges.'
       return
    endif

    !set up some initial times
    call get_time(Time%model_time,secs,days)
    initial_day = days
    initial_secs = secs

    if ( limit_temp ) then
       write(stdoutunit,*) '==> Preventing freezing by restoring to -1.8C where applicable'
    endif
    if ( limit_salt ) then
       write(stdoutunit,*) '==> Preventing excessive freshening by restoring salinity to 0.1 where applicable'
       write(stdoutunit,*) '==> Preventing excessive saltiness by restoring salinity to 42.5 where applicable'
    endif

    do n=1,num_prog_tracers
       if ( restore_cars_deep ) then
          name = 'INPUT/'//trim(T_prog(n)%name)//'_cars.nc'
          if (file_exist(trim(name))) then
             write(stdoutunit,*) '==> damping '//trim(T_prog(n)%name)//' at depth from file '//trim(name) 
             allocate(Cars(n)%field(isd:ied,jsd:jed,nk))
             call read_data(name,T_prog(n)%name,Cars(n)%field,domain=Domain%domain2d,timelevel=1)
             Cars(n)%id = 1
          else
             Cars(n)%id = 0
          endif
       endif

       ! read damping rates (1/sec) 
       if (.not. use_hard_thump  .and. .not. use_adaptive_restore) then
          allocate(Sponge(n)%damp_coeff(isd:ied,jsd:jed,nk))
          Sponge(n)%damp_coeff=0.0
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
             write(stdoutunit,*) '==> '//trim(name)//' not found.  Sponge damping coefficients not being applied '
          endif !(file_exist(trim(name)))
       endif !.not. use_hard_thump .and. .not. use_adaptive_restore

       ! read restoring data
       name = 'INPUT/'//trim(T_prog(n)%name)//'_sponge.nc'
       if (file_exist(trim(name)) ) then
          Sponge(n)%id = init_external_field(name,T_prog(n)%name,domain=Domain%domain2d)
          if ( Sponge(n)%id < 1) then 
             call mpp_error(FATAL,&
                  '==>Error: in ocean_sponges_tracer_mod: Sponge values are required but not specified')
          endif
          write(stdoutunit,*) '==> Using restoring data specified from file '//trim(name)
       else 
          write(stdoutunit,*) '==> '//trim(name)//' not found. Sponge not being applied '
       endif
    enddo !n=1,num_prog_tracers

    ! register diagnostic output
    allocate (id_sponge_tend(num_prog_tracers))
    id_sponge_tend = -1

    do n=1,num_prog_tracers
       if(n==index_temp) then
          id_sponge_tend(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sponge_tend',&
               Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*cp*heating due to sponge',        &
               trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e6,1.e6/))
       else
          id_sponge_tend(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sponge_tend',&
               Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*tendency due to sponge',          &
               trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e-1,1.e-1/))
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
    logical,save :: first_pass = .true.
    real, parameter :: one_over_year = 1./(365.0*86400.0)

    if(.not. use_this_module) return 

    taum1 = Time%taum1
    tau   = Time%tau
    wrk1  = 0.0 


    sdiff=0.0
    call get_time(Time%model_time,secs,days)

    do n=1,size(T_prog(:))
       ! Need to reinitialise wrk2 here due to limiting of temperature.
       wrk2  = 0.0
       !-----------------------------------------------------------------------------------------------------------------
       if (Sponge(n)%id > 0) then
          ! get sponge value for current time
          call time_interp_external(Sponge(n)%id, Time%model_time, wrk1)
          !--------------------------------------------------------------------------------------------------------------
          if ( first_pass ) then
             if ( use_hard_thump ) then
                do k = 1,nk
                   do j=jsd,jed
                      do i=isd,ied
                         T_prog(n)%field(i,j,k,taum1)=wrk1(i,j,k)
                         T_prog(n)%field(i,j,k,tau)=wrk1(i,j,k)
                      enddo
                   enddo
                enddo
                call mpp_update_domains ( T_prog(n)%field(:,:,:,tau),Dom%domain2d) 
                call mpp_update_domains ( T_prog(n)%field(:,:,:,taum1),Dom%domain2d) 
                wrk2=abs(T_prog(n)%field(:,:,:,taum1)-wrk1(:,:,:))*Grd%tmask(:,:,:)
             endif !use_hard_thump
             sum_val=sum(wrk2(isc:iec,jsc:jec,:))
             call  mpp_sum(sum_val)
             sdiffo(n) = sum_val/Grd%wet_t_points
             write (stdout(),*)'mean absolute deviation ',trim(T_prog(n)%name),' = ',sdiffo(n)
             if ( .not. use_normalising ) then
                sdiffo(n)=1.0
             else
                sdiffo(n)  = maxval(wrk2)
                call mpp_max(sdiffo(n))
                sdiffo(n)  = 1.01*sdiffo(n)
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
                         sdiff=abs(T_prog(n)%field(i,j,k,taum1)-wrk1(i,j,k))
                         sdiff = max(sdiff,1e-9) !minimum tolerance
                         adaptive_coeff = 1.0/(taumin - (lambda*real(secs_restore))/(((sdiff/sdiffo(n))**npower)*log(1-athresh)))
                         wrk2(i,j,k) = Thickness%rho_dzt(i,j,k,tau)* adaptive_coeff    &
                              *(wrk1(i,j,k)-T_prog(n)%field(i,j,k,taum1))
                      enddo
                   enddo
                enddo
             else !normal restore
                do k=1,nk
                   do j=jsd,jed
                      do i=isd,ied  
                         wrk2(i,j,k) = Thickness%rho_dzt(i,j,k,tau)*Sponge(n)%damp_coeff(i,j,k)    &
                              *(wrk1(i,j,k)-T_prog(n)%field(i,j,k,taum1))
                      enddo
                   enddo
                enddo
             endif
          endif !use_hard_thump
          !--------------------------------------------------------------------------------------------------------------
       endif !Sponge(1)%id
       !-----------------------------------------------------------------------------------------------------------------
       if ( restore_cars_deep .and. Cars(n)%id == 1 ) then
          do k=nk-4,nk
             do j=jsd,jed
                do i=isd,ied  
                   wrk2(i,j,k) = Thickness%rho_dzt(i,j,k,tau)*one_over_year    &
                        *(Cars(n)%field(i,j,k)-T_prog(n)%field(i,j,k,taum1))
                enddo
             enddo
          enddo
       endif
       !--------------------------------------------------------------------------------------------------------------
       !this limiting scheme restores to sponge data rather than hard limits
       if (trim(T_prog(n)%name) == 'temp' .and. limit_temp) then
          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   if ( Grd%tmask(i,j,k)==1.0 ) then
                      if ( T_prog(n)%field(i,j,k,taum1).ge.limit_temp_max ) then
                         sdiff=abs(limit_temp_max-T_prog(n)%field(i,j,k,taum1))
                         sdiff = max(sdiff,1e-9) !minimum tolerance
                         adaptive_coeff = 1.0/(taumin-(lambda*real(86400))/(((sdiff/sdiffo(n))**npower)*log(1-athresh)))
                         wrk2(i,j,k) = Grd%tmask(i,j,k)*Thickness%rho_dzt(i,j,k,tau)*adaptive_coeff    &
                              *(limit_temp_max-T_prog(n)%field(i,j,k,taum1))
                      endif
                      if ( T_prog(n)%field(i,j,k,taum1).le.limit_temp_min ) then
                         sdiff=abs(limit_temp_min-T_prog(n)%field(i,j,k,taum1))
                         sdiff = max(sdiff,1e-9) !minimum tolerance
                         adaptive_coeff = 1.0/(taumin-(lambda*real(86400))/(((sdiff/sdiffo(n))**npower)*log(1-athresh)))
                         wrk2(i,j,k) = Grd%tmask(i,j,k)*Thickness%rho_dzt(i,j,k,tau)*adaptive_coeff    &
                              *(limit_temp_min-T_prog(n)%field(i,j,k,taum1))
                      endif
                   endif
                enddo
             enddo
          enddo
       endif !trim
       !--------------------------------------------------------------------------------------------------------------
       if (trim(T_prog(n)%name) == 'salt' .and. limit_salt) then
          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   if ( Grd%tmask(i,j,k)==real(1) ) then
                      if ( T_prog(n)%field(i,j,k,taum1).ge.limit_salt_max ) then
                         sdiff=abs(limit_salt_max-T_prog(n)%field(i,j,k,taum1))
                         sdiff = max(sdiff,1e-9) !minimum tolerance
                         adaptive_coeff = 1.0/(taumin-(lambda*real(86400))/(((sdiff/sdiffo(n))**npower)*log(1-athresh)))
                         wrk2(i,j,k) = Grd%tmask(i,j,k)*Thickness%rho_dzt(i,j,k,tau)*adaptive_coeff    &
                              *(limit_salt_max-T_prog(n)%field(i,j,k,taum1))
                      endif
                      if ( T_prog(n)%field(i,j,k,taum1).le.limit_salt_min ) then
                         sdiff=abs(limit_salt_min-T_prog(n)%field(i,j,k,taum1))
                         sdiff = max(sdiff,1e-9) !minimum tolerance
                         adaptive_coeff = 1.0/(taumin-(lambda*real(86400))/(((sdiff/sdiffo(n))**npower)*log(1-athresh)))
                         wrk2(i,j,k) = Grd%tmask(i,j,k)*Thickness%rho_dzt(i,j,k,tau)*adaptive_coeff    &
                              *(limit_salt_min-T_prog(n)%field(i,j,k,taum1))
                      endif
                   endif
                enddo
             enddo
          enddo
       endif !trim
       !--------------------------------------------------------------------------------------------------------------
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec  
                T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + wrk2(i,j,k)
             enddo
          enddo
       enddo
       if (id_sponge_tend(n) > 0) used = send_data(id_sponge_tend(n),                 &
            T_prog(n)%conversion*wrk2(:,:,:), Time%model_time, rmask=Grd%tmask(:,:,:), &
            is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
       !--------------------------------------------------------------------------------------------------------------
    enddo !do n=1,size(T_prog(:))

    first_pass = .false.

    return

  end subroutine sponge_tracer_source
  ! </SUBROUTINE> NAME="sponge_tracer_source"
end module ocean_sponges_tracer_mod

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
module ocean_sponges_eta_mod
  !
  !<CONTACT EMAIL="p.sandery@bom.gov.au"> Paul Sandery
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
  ! RASF: Added Paul's adaptive restoring scheme
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
  !<NAMELIST NAME="ocean_sponges_eta_OFAM_nml">
  !
  !  <DATA NAME="use_adaptive_restore" TYPE="logical">
  !  Use adaptive restoring scheme.
  !  </DATA> 
  !  <DATA NAME="use_normalising" TYPE="logical">
  !   Normalise????
  !  </DATA> 
  !  <DATA NAME="use_hard_thump" TYPE="logical">
  !  Force surface height to a fixed value on first time step.
  !  </DATA> 
  !  <DATA NAME="athresh" TYPE="real">
  !  ...
  !  </DATA> 
  !  <DATA NAME="taumin" TYPE="real">
  !  Shortest restoring timescale?
  !  </DATA> 
  !  <DATA NAME="lambda" TYPE="real">
  !  ...
  !  </DATA> 
  !  <DATA NAME="npower" TYPE="real">
  !  Power law for restoring. RMSE is raised to this power.
  !  </DATA> 
  !</NAMELIST>
  use diag_manager_mod,         only: register_diag_field, send_data
  use fms_mod,                  only: write_version_number, open_namelist_file, close_file
  use fms_mod,                  only: file_exist
  use fms_mod,                  only: open_namelist_file, check_nml_error, close_file
  use fms_mod,                  only: read_data, lowercase, FATAL, WARNING, stdout, stdlog
  use mpp_mod,                  only: mpp_sum, mpp_error, mpp_max
  use mpp_domains_mod,          only: mpp_global_sum, mpp_update_domains
  use time_interp_external_mod, only: init_external_field, time_interp_external
  use time_manager_mod,         only: time_type, set_date, get_time
  use time_manager_mod,         only: operator( + ), operator( - ), operator( // )
  use time_manager_mod,         only: operator( > ), operator( == ), operator( <= )

  use ocean_domains_mod,        only: get_local_indices
  use ocean_parameters_mod,     only: missing_value, rho0 
  use ocean_types_mod,          only: ocean_domain_type, ocean_grid_type, ocean_options_type
  use ocean_types_mod,          only: ocean_time_type, ocean_External_mode_type
  use ocean_workspace_mod,      only: wrk1_2d, wrk2_2d

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

  character(len=126)  :: version = '$Id: ocean_sponges_eta.F90,v 1.1.2.4.42.1 2009/10/10 00:42:52 nnz Exp $'
  character (len=128) :: tagname = '$Name: mom4p1_pubrel_dec2009_nnz $'

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
  real    :: taumin                = 3600
  logical :: use_adaptive_restore  = .false.
  logical :: use_sponge_after_init = .false.
  logical :: use_normalising       = .false.
  logical :: use_hard_thump        = .false.
  integer :: secs_to_restore       = 0
  integer :: days_to_restore       = 0
  integer :: secs_restore
  integer :: initial_day, initial_secs

  namelist /ocean_sponges_eta_nml/ use_this_module 
  namelist /ocean_sponges_eta_OFAM_nml/ use_adaptive_restore,use_sponge_after_init, use_normalising, use_hard_thump, &
       athresh, taumin,lambda,npower,days_to_restore, secs_to_restore

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
    type(ocean_domain_type),     intent(in), target :: Domain
    type(ocean_time_type),         intent(in)            :: Time
    real,                                    intent(in)           :: dtime
    type(ocean_options_type),     intent(inout)       :: Ocean_options

    integer :: i, j, k, n
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

    call write_version_number()

    ! provide for namelist over-ride of default values
    ioun =  open_namelist_file()
    read (ioun,ocean_sponges_eta_nml,IOSTAT=io_status)
    write (stdlogunit,ocean_sponges_eta_nml)
    write (stdoutunit,'(/)') 
    write (stdoutunit,ocean_sponges_eta_nml)
    ierr = check_nml_error(io_status,'ocean_sponges_eta_nml')
    call close_file (ioun)

    ! provide for namelist over-ride of default values
    ioun =  open_namelist_file()
    read (ioun,ocean_sponges_eta_OFAM_nml,IOSTAT=io_status)
    write (stdlogunit,ocean_sponges_eta_OFAM_nml)
    write (stdoutunit,'(/)') 
    write (stdoutunit,ocean_sponges_eta_OFAM_nml)
    ierr = check_nml_error(io_status,'ocean_sponges_eta_OFAM_nml')
    call close_file (ioun)

    Dom => Domain
    Grd => Grid

#ifndef MOM4_STATIC_ARRAYS    
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

    !set up some initial times
    call get_time(Time%model_time,secs,days)
    initial_day = days
    initial_secs = secs

    !read damping rates for eta
    if (.not. use_hard_thump  .and. .not. use_adaptive_restore) then
       allocate(Sponge(1)%damp_coeff2d(isd:ied,jsd:jed))
       Sponge(1)%damp_coeff2d(:,:) = 0.0
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
          !modify damping rates to allow restoring to be solved implicitly
          !note: test values between zero and 4.0e-8 revert to damping rates defined above
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
       endif !(file_exist(name))
    endif !.not. use_hard_thump .and. .not. use_adaptive_restore

    name = 'INPUT/eta_sponge.nc'
    if(file_exist(trim(name))) then
       Sponge(1)%id = init_external_field(name,'eta_t',domain=Domain%domain2d)
       if (Sponge(1)%id < 1) then
          call mpp_error(FATAL,&
               '==>Error: in ocean_sponges_eta_mod: sponge values are required but not specified')
       endif
       write(stdoutunit,*) '==> Using sponge data specified from file '//trim(name)
    else
       write(stdoutunit,*) '==> '//trim(name)//' not found.  Sponge not being applied '
    endif

    allocate (id_sponge_tend(1)) 
    id_sponge_tend = -1 ! needed?
    id_sponge_tend = register_diag_field ('ocean_model', 'eta_t_sponge_tend', &
         Grd%tracer_axes(1:2), Time%model_time, 'eta_t tendency due to sponge', &
         'm/s', missing_value=missing_value, range=(/-1.e-2,1.e-2/))

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
    type(ocean_time_type),                intent(in)     :: Time
    type(ocean_external_mode_type), intent(inout) :: Ext_mode
    integer :: taum1, tau
    integer :: i, j
    real      :: sdiff=0.0, sum_val, adaptive_coeff
    integer :: secs, days, numsecs
    logical  :: do_adaptive_restore = .false.  ! Only restore in specified time period.
    logical,save :: first_pass = .true.
    real,   save :: num_surface_wet_cells

    if(.not. use_this_module) return 

    taum1 = Time%taum1
    tau   = Time%tau
    wrk1_2d = 0.0
    wrk2_2d = 0.0
    call get_time(Time%model_time,secs,days)
    !-----------------------------------------------------------------------------------------------------------------
    if (Sponge(1)%id > 0) then
       call time_interp_external(Sponge(1)%id, Time%model_time, wrk1_2d)
       if (first_pass) then
          if (use_hard_thump) then
             do j=jsd,jed
                do i=isd,ied
                   Ext_mode%eta_t(i,j,taum1)=wrk1_2d(i,j) ! 6:
                   Ext_mode%eta_t(i,j,tau)=wrk1_2d(i,j)
                enddo
             enddo
             call mpp_update_domains (Ext_mode%eta_t(:,:,taum1) ,Dom%domain2d)
             call mpp_update_domains (Ext_mode%eta_t(:,:,tau) ,Dom%domain2d)
             wrk2_2d=abs(Ext_mode%eta_t(:,:,taum1)-wrk1_2d(:,:))*Grd%tmask(:,:,1) ! 8: 9/10 only get called if H=0 so 8 must mean H=1
          endif !use_hard_thump
          sum_val=sum(wrk2_2d(isc:iec,jsc:jec))
          call mpp_sum(sum_val)
          num_surface_wet_cells = mpp_global_sum(Dom%domain2d,Grd%tmask(:,:,1)) ! 1:
          sdiffo = sum_val/num_surface_wet_cells
          write (stdout(),*)'mean absolute deviation eta ', sdiffo
          if ( use_normalising == .false.) then
             sdiffo=1.0
          else
             sdiffo  = maxval(wrk2_2d)
             call mpp_max(sdiffo)
             sdiffo  = 1.01*sdiffo
          endif
       endif !first_pass
       !--------------------------------------------------------------------------------------------------------------
       if (use_hard_thump == .false.) then
          if (use_adaptive_restore) then
             secs_restore = days_to_restore*86400 + secs_to_restore
             numsecs = (days-initial_day)*86400 + secs - initial_secs
             do_adaptive_restore = ( numsecs < secs_restore )
          endif
          if (do_adaptive_restore) then
             ! calculate idealised forcing tendency timescale for each timestep and array element
             do j=jsd,jed
                do i=isd,ied
                   sdiff=abs(Ext_mode%eta_t(i,j,taum1)-wrk1_2d(i,j))
                   sdiff = max(sdiff,1e-9) ! minimum tolerance
                   adaptive_coeff =1./(taumin - (lambda*real(secs_restore))/(((sdiff/sdiffo)**npower)*log(1-athresh)))
                   wrk2_2d(i,j) =  adaptive_coeff * (wrk1_2d(i,j) - Ext_mode%eta_t(i,j,taum1))
                enddo
             enddo
          else !normal restore
             do j=jsd,jed
                do i=isd,ied
                   wrk2_2d(i,j) = Sponge(1)%damp_coeff2d(i,j)*(wrk1_2d(i,j) - Ext_mode%eta_t(i,j,taum1))
                enddo
             enddo
          endif
       endif !use_hard_thump
       !--------------------------------------------------------------------------------------------------------------
       do j=jsc,jec
          do i=isc,iec
             Ext_mode%source(i,j) = Ext_mode%source(i,j) + rho0*wrk2_2d(i,j)
          enddo
       enddo
    endif !Sponge(1)%id
    !-----------------------------------------------------------------------------------------------------------------
    if (id_sponge_tend(1) > 0) used = send_data(id_sponge_tend(1), &
         wrk2_2d(:,:), Time%model_time, rmask=Grd%tmask(:,:,1),      &
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    first_pass = .false.

    return

  end subroutine sponge_eta_source
  ! </SUBROUTINE> NAME="sponge_eta_source"
end module ocean_sponges_eta_mod

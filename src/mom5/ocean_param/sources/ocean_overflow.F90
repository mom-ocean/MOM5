module ocean_overflow_mod
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="m.bates@student.unsw.edu.au"> Michael Bates 
!</CONTACT>
!
!<OVERVIEW>
! Tracer source from discharging dense shallow water into the abyss
! at the parcel's depth of neutral buoyancy.  
!</OVERVIEW>
!
!<DESCRIPTION>
! Tracer source from discharging dense shallow water into the abyss
! at the parcel's depth of neutral buoyancy.  Follow the approach
! of Campin and Goosse (1999), as well as modifications.  
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Campin and Goosse (1999): Parameterization of density-driven downsloping flow 
! for a coarse-resolution model in z-coordinate", Tellus 51A, pages 412-430
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R. C. Pacanowski, and A. Rosati
! A Technical Guide to MOM4 (2003)
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, Elements of MOM (2012)
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_overflow_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  For using this module.  Default is false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging 
!  </DATA> 
!
!  <DATA NAME="overflow_mu" TYPE="real" UNITS="inverse seconds">
!  Dissipation rate for the bottom friction.  Campin and Goosse 
!  suggest overflow_mu=10^-4
!  </DATA> 
!  <DATA NAME="overflow_delta" TYPE="real" UNITS="dimensionless">
!  Fraction of a grid cell participating in the overflow process. 
!  Campin and Goosse suggest overflow_delta=1/3. 
!  </DATA> 
!  <DATA NAME="overflow_umax" TYPE="real" UNITS="m/s">
!  Maximum downslope speed. 
!  </DATA> 
!  <DATA NAME="do_bitwise_exact_sum" TYPE="logical">
!  Set true to do bitwise exact global sum. When it is false, the global
!  sum will be non-bitwise_exact, but will significantly increase efficiency. 
!  The default value is false.
!  </DATA> 
!
!  <DATA NAME="no_return_flow" TYPE="logical">
!  Set true to remove the Campin and Goose return flow "piping".
!  Default no_return_flow=.false. to recover the standard approach from 
!  Campin and Goose. 
!  </DATA> 
!
!</NAMELIST>
!

use constants_mod,       only: epsln
use diag_manager_mod,    only: register_diag_field, register_static_field
use fms_mod,             only: write_version_number, error_mesg, FATAL, NOTE
use fms_mod,             only: open_namelist_file, check_nml_error, close_file, stdout, stdlog
use mpp_domains_mod,     only: mpp_update_domains, CGRID_NE
use mpp_domains_mod,     only: mpp_global_sum, BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
use mpp_domains_mod,     only: mpp_define_domains, mpp_update_domains
use mpp_domains_mod,     only: cyclic_global_domain, global_data_domain 
use mpp_mod,             only: input_nml_file, mpp_error

use ocean_density_mod,    only: density
use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: missing_value, rho0r, grav
use ocean_parameters_mod, only: TERRAIN_FOLLOWING
use ocean_types_mod,      only: ocean_domain_type, ocean_grid_type, ocean_time_type, ocean_options_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_density_type, ocean_thickness_type
use ocean_util_mod,       only: write_timestamp, diagnose_2d, diagnose_3d
use ocean_tracer_util_mod,only: diagnose_3d_rho
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk4, wrk1_v, wrk2_v

implicit none

private 

public ocean_overflow_init
public overflow
private watermass_diag_init
private watermass_diag


type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

#include <ocean_memory.h>

! for the algorithm
real    :: overflow_factor ! grav*overflow_delta/(rho0*overflow_mu)  (m^4/kg/sec)
integer :: ip(4),jq(4)     ! to label four surrounding cells 
integer :: ijsc(4),ijec(4) ! for shifting do-loop indices to minimize calls to mpp_update_domains

real, dimension(:,:),   allocatable   :: slope_x          ! topog slope (dimensionless) in the i-direction
real, dimension(:,:),   allocatable   :: slope_y          ! topog slope (dimensionless) in the j-direction 
real, dimension(:,:,:), allocatable   :: topog_slope      ! topog slope (dimensionless) for 4-directions 
real, dimension(:,:,:), allocatable   :: topog_step       ! where downslope flow is topgraphically possible 
real, dimension(:,:,:), allocatable   :: flux_int_z       ! for diagnostics 
real, dimension(:,:,:), allocatable   :: source_overflow  ! thickness weighted tendency from overflow  
real, dimension(:,:,:), allocatable   :: overflow_flux    ! tracer flux associated with overflow 
real, dimension(:,:),   allocatable   :: overflow_xflux   ! mass flux in i-direction
real, dimension(:,:),   allocatable   :: overflow_yflux   ! mass flux in j-direction 

! internally set for computing watermass diagnostics
logical :: compute_watermass_diag = .false. 

! for diagnostic manager 
logical :: used
integer :: id_over_slope_x=-1
integer :: id_over_slope_y=-1
integer :: id_topog_step_1=-1
integer :: id_topog_step_2=-1
integer :: id_topog_step_3=-1
integer :: id_topog_step_4=-1
integer :: id_overflow_xflux=-1 
integer :: id_overflow_yflux=-1
integer :: id_tx_trans_overflow=-1
integer :: id_ty_trans_overflow=-1
integer, dimension(:), allocatable :: id_overflow
integer, dimension(:), allocatable :: id_overflow_xflux_int_z
integer, dimension(:), allocatable :: id_overflow_yflux_int_z

integer :: id_neut_rho_overfl          =-1
integer :: id_wdian_rho_overfl         =-1
integer :: id_neut_rho_overfl_on_nrho  =-1
integer :: id_wdian_rho_overfl_on_nrho =-1
integer :: id_tform_rho_overfl_on_nrho =-1

integer :: id_neut_temp_overfl          =-1
integer :: id_wdian_temp_overfl         =-1
integer :: id_neut_temp_overfl_on_nrho  =-1
integer :: id_wdian_temp_overfl_on_nrho =-1
integer :: id_tform_temp_overfl_on_nrho =-1

integer :: id_neut_salt_overfl          =-1
integer :: id_wdian_salt_overfl         =-1
integer :: id_neut_salt_overfl_on_nrho  =-1
integer :: id_wdian_salt_overfl_on_nrho =-1
integer :: id_tform_salt_overfl_on_nrho =-1

! for specifying transport units
! can either be Sv or mks
character(len=32) :: transport_units='Sv' 
character(len=32) :: transport_dims ='Sv (10^9 kg/s)' 
real              :: transport_convert=1.0e-9 

character(len=128) :: version=&
       '=>Using: ocean_overflow.f90 ($Id: ocean_overflow.F90,v 20.0 2013/12/14 00:16:06 fms Exp $)'
character (len=128) :: tagname=&
     '$Name: tikal $'

! number of prognostic tracers
integer :: num_prog_tracers=0

! processor zero writes to unit 6
integer :: unit=6  

! initialization flag 
logical :: module_is_initialized=.false.

! flag for mpp_global_sum
integer :: global_sum_flag      

! time step 
real :: dtime 
real :: p5dtimer

! set from nml
logical :: use_this_module        = .false. ! must be set .true. in nml to enable this scheme
logical :: debug_this_module      = .false. ! for debugging
logical :: do_bitwise_exact_sum   = .false. ! set true to get slower sum that is same for all PEs. 
logical :: no_return_flow         = .false. ! set true to remove the pre-defined return flow "piping"
real    :: overflow_mu            = 1.e-4   ! frictional dissipation rate (sec^-1) at bottom  
real    :: overflow_delta         = 0.3333  ! fraction of a grid cell participating in overflow 
real    :: overflow_umax          = 0.01    ! maximum downslope speed (m/s)

namelist /ocean_overflow_nml/ use_this_module, debug_this_module, overflow_mu,     &
                              overflow_delta, overflow_umax, do_bitwise_exact_sum, &
                              no_return_flow, transport_units 

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_overflow_init">
!
! <DESCRIPTION>
! Initial set up for mixing of tracers into the abyss next to topography.
! </DESCRIPTION>
!
  subroutine ocean_overflow_init(Grid, Domain, Time, Dens, T_prog, Ocean_options, &
                                 vert_coordinate_type, dtim, debug)

    type(ocean_grid_type),        intent(in), target   :: Grid
    type(ocean_domain_type),      intent(in), target   :: Domain
    type(ocean_time_type),        intent(in), target   :: Time
    type(ocean_density_type),     intent(in)           :: Dens
    type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
    type(ocean_options_type),     intent(inout)        :: Ocean_options
    integer,                      intent(in)           :: vert_coordinate_type 
    real,                         intent(in)           :: dtim
    logical,                      intent(in), optional :: debug

    integer :: io_status, ioun, ierr
    integer :: i,j,m,n,kbot

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

    if ( module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_overflow_mod (ocean_overflow_init): module already initialized')
    endif 

    module_is_initialized = .TRUE.

    call write_version_number( version, tagname )
#ifndef MOM_STATIC_ARRAYS
    call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
    nk = Grid%nk
#endif
    Dom => Domain
    Grd => Grid

    num_prog_tracers = size(T_prog(:))
    dtime = dtim 
    p5dtimer = 0.5/dtime

    ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_overflow_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_overflow_nml')
#else
    ioun =  open_namelist_file()
    read (ioun,ocean_overflow_nml,IOSTAT=io_status)
    ierr = check_nml_error(io_status,'ocean_overflow_nml')
    call close_file (ioun)
#endif
    write (stdlogunit,ocean_overflow_nml)
    write (stdoutunit,'(/)') 
    write (stdoutunit,ocean_overflow_nml)

    if(do_bitwise_exact_sum) then
       global_sum_flag = BITWISE_EXACT_SUM
    else
       global_sum_flag = NON_BITWISE_EXACT_SUM
    endif

    if (PRESENT(debug) .and. .not. debug_this_module) then
       debug_this_module = debug
    endif 
    if(debug_this_module) then 
       write(stdoutunit,'(a)') '==>Note: running ocean_overflow_mod with debug_this_module=.true.'  
    endif 
    
    if(transport_units=='Sv') then
        transport_convert=1.0e-9 
        transport_dims   = 'Sv (10^9 kg/s)'
    else
        transport_convert=1.0
        transport_dims   = 'kg/s'
    endif


    if(.not. use_this_module) then 
      call mpp_error(NOTE,&
      '==>From ocean_overflow_mod: NOT using Campin and Goosse overflow scheme.')
      Ocean_options%overflow = 'Did NOT use Campin and Goosse overflow scheme.'
      return
    else 
      if(vert_coordinate_type == TERRAIN_FOLLOWING) then 
          call mpp_error(FATAL, &
          '==>ocean_overflow_mod: this module is NOT for use with TERRAIN_FOLLOWING vert coodinates.')
      endif
      Ocean_options%overflow = 'Used the Campin and Goosse overflow scheme.'

      if(no_return_flow) then
         call mpp_error(NOTE,&
         '==>From ocean_overflow_mod: USING modified Campin and Goosse overflow scheme with no specified return flow.')
      else
         call mpp_error(NOTE,&
         '==>From ocean_overflow_mod: USING original Campin and Goosse overflow scheme with specified return flow.')
      endif

    endif

    if(debug_this_module) then 
      call mpp_error(NOTE,'==>From ocean_overflow_mod: USING debug_this_module')
    endif 

    ! for holding the tracer tendency from the overflow scheme 
    allocate (source_overflow(isd:ied,jsd:jed,nk))
    source_overflow(:,:,:) = 0.0    

    ! compute a common factor that is time independent 
    overflow_factor = grav*overflow_delta*rho0r/(epsln+overflow_mu)

    ! compute topographic slope arrays for the i-slope and j-slope
    ! slopes are centered on the i-face and j-face of tracer cells 
    allocate (slope_x(isd:ied,jsd:jed))
    allocate (slope_y(isd:ied,jsd:jed))
    slope_x = 0.0 
    slope_y = 0.0
    do j=jsc,jec
      do i=isc,iec
        slope_x(i,j) = (Grd%ht(i+1,j)-Grd%ht(i,j))*Grd%dxter(i,j)     
        slope_y(i,j) = (Grd%ht(i,j+1)-Grd%ht(i,j))*Grd%dytnr(i,j)
        slope_x(i,j) = abs(slope_x(i,j))
        slope_y(i,j) = abs(slope_y(i,j))
      enddo
    enddo 
    call mpp_update_domains(slope_x(:,:),Dom%domain2d)
    call mpp_update_domains(slope_y(:,:),Dom%domain2d)

    ! topographic slope for the four surrounding directions 
    allocate (topog_slope(isd:ied,jsd:jed,4))
    topog_slope(:,:,:) = 0.0
    do j=jsc,jec
       do i=isc,iec
          m=1 ; topog_slope(i,j,m) = slope_x(i,j)  
          m=2 ; topog_slope(i,j,m) = slope_y(i,j)  
          m=3 ; topog_slope(i,j,m) = slope_x(i-1,j)
          m=4 ; topog_slope(i,j,m) = slope_y(i,j-1)
       enddo
    enddo

    ! compute directions from an (i,j) point where topography deepens.
    ! these directions may potentially have downslope flow. 
    ! insist that downslope flow occurs only when there are more kmt cells
    ! in the adjacent column. 
    ! also insist that downslope flow does not involve k=1 cells. 
    allocate (topog_step(isd:ied,jsd:jed,4))
    topog_step(:,:,:) = 0.0
    do j=jsc,jec
      do i=isc,iec
         kbot = Grd%kmt(i,j)
         if(kbot > 1) then 
           if(Grd%kmt(i+1,j) > kbot) topog_step(i,j,1)=1.0
           if(Grd%kmt(i,j+1) > kbot) topog_step(i,j,2)=1.0
           if(Grd%kmt(i-1,j) > kbot) topog_step(i,j,3)=1.0
           if(Grd%kmt(i,j-1) > kbot) topog_step(i,j,4)=1.0
         endif 
      enddo
    enddo   
    ! Block out the bipolar fold in order to ensure tracer conservation.
    ! The reason we do so is related to how the algorithm reaches between  
    ! two adjacent columns of tracer points.  When the column straddles 
    ! the bipolar fold, the code logic is not general and so it actually 
    ! attempts to reach to a non-adjacent column.  Shy of generalizing 
    ! the logic, we simply do not consider overflow for points along fold. 
    if(jec+Dom%joff==Dom%jeg) then 
       topog_step(:,jec,:) = 0.0
    endif 

    ! labels for the four tracer cells surrounding a 
    ! central tracer cell (moving counter-clockwise)
    m=1 ; ip(m)=1  ; jq(m)=0  ; ijsc(m) = 1 ; ijec(m) = 0
    m=2 ; ip(m)=0  ; jq(m)=1  ; ijsc(m) = 1 ; ijec(m) = 0
    m=3 ; ip(m)=-1 ; jq(m)=0  ; ijsc(m) = 0 ; ijec(m) = 1
    m=4 ; ip(m)=0  ; jq(m)=-1 ; ijsc(m) = 0 ; ijec(m) = 1 

    ! for overflow flux 
    allocate (overflow_flux(isd:ied,jsd:jed,4))
    overflow_flux(:,:,:) = 0.0

    ! for diagnosing vertically integrated tracer flux
    allocate (flux_int_z(isd:ied,jsd:jed,4))
    flux_int_z = 0.0

    ! for diagnostics with the traditional overflow scheme 
    allocate (overflow_xflux(isd:ied,jsd:jed))
    allocate (overflow_yflux(isd:ied,jsd:jed))
    overflow_xflux(:,:) = 0.0
    overflow_yflux(:,:) = 0.0


    ! register/send diagnostics 
    id_over_slope_x = register_static_field ('ocean_model', 'over_slope_x', Grd%tracer_axes_flux_x(1:2), &
                 '|d(ht)/dx| on T-cell face', 'm/m', missing_value=missing_value, range=(/-1.e9,1.e9/))
    call diagnose_2d(Time, Grd, id_over_slope_x, slope_x(:,:))

    id_over_slope_y = register_static_field ('ocean_model', 'over_slope_y', Grd%tracer_axes_flux_y(1:2), &
                 '|d(ht)/dy| on T-cell face', 'm/m', missing_value=missing_value, range=(/-1.e9,1.e9/))
    call diagnose_2d(Time, Grd, id_over_slope_y, slope_y(:,:))

    id_topog_step_1 = register_static_field ('ocean_model', 'topog_step_1', Grd%tracer_axes(1:2), &
                 'topog_step_1', 'dimensionless', missing_value=missing_value, range=(/-1.0,1.0/))
    call diagnose_2d(Time, Grd, id_topog_step_1, topog_step(:,:,1))

    id_topog_step_2 = register_static_field ('ocean_model', 'topog_step_2', Grd%tracer_axes(1:2), &
                 'topog_step_2', 'dimensionless', missing_value=missing_value, range=(/-1.0,1.0/))
    call diagnose_2d(Time, Grd, id_topog_step_2, topog_step(:,:,2))

    id_topog_step_3 = register_static_field ('ocean_model', 'topog_step_3', Grd%tracer_axes(1:2), &
                 'topog_step_3', 'dimensionless', missing_value=missing_value, range=(/-1.0,1.0/))
    call diagnose_2d(Time, Grd, id_topog_step_3, topog_step(:,:,3))

    id_topog_step_4 = register_static_field ('ocean_model', 'topog_step_4', Grd%tracer_axes(1:2), &
                 'topog_step_4', 'dimensionless', missing_value=missing_value, range=(/-1.0,1.0/))
    call diagnose_2d(Time, Grd, id_topog_step_4, topog_step(:,:,4))


    id_overflow_xflux = register_diag_field ('ocean_model', 'overflow_xflux', Grd%tracer_axes_flux_x(1:2), &
         Time%model_time, 'off-shelf i-transport from overflow', &
        'kg/s', missing_value=missing_value, range=(/-1.e14,1.e14/))

    id_overflow_yflux = register_diag_field ('ocean_model', 'overflow_yflux', Grd%tracer_axes_flux_y(1:2), &
         Time%model_time, 'off-shelf j-transport from overflow', &
        'kg/s', missing_value=missing_value, range=(/-1.e9,1.e9/))

    id_tx_trans_overflow = register_diag_field ('ocean_model', 'tx_trans_overflow',               &
         Grd%tracer_axes_flux_x(1:3), Time%model_time, 'i-directed mass transport from overflow', &
        trim(transport_dims), missing_value=missing_value, range=(/-1.e20,1.e20/))

    id_ty_trans_overflow = register_diag_field ('ocean_model', 'ty_trans_overflow',               &
         Grd%tracer_axes_flux_x(1:3), Time%model_time, 'j-directed mass transport from overflow', &
        trim(transport_dims), missing_value=missing_value, range=(/-1.e20,1.e20/))

    call watermass_diag_init(Time,Dens)

    allocate (id_overflow(num_prog_tracers))
    allocate (id_overflow_xflux_int_z(num_prog_tracers))
    allocate (id_overflow_yflux_int_z(num_prog_tracers))
    id_overflow = -1
    id_overflow_xflux_int_z=-1
    id_overflow_yflux_int_z=-1

    do n=1,num_prog_tracers
       if(T_prog(n)%name == 'temp') then 
           id_overflow(n) = register_diag_field ('ocean_model', &
                'overflow_'//trim(T_prog(n)%name),              &
                Grd%tracer_axes(1:3), Time%model_time,          &
                'cp*overflow*rho*dzt*temp',                     &
                'Watt/m^2', missing_value=missing_value, range=(/-1.e9,1.e9/))
           id_overflow_xflux_int_z(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_xflux_overflow_int_z',                &
              Grd%tracer_axes_flux_x(1:2), Time%model_time,                 &
              'vert-integ cp*overflow_xflux*rho*dzt*dyt*temp',              &
              'Watt', missing_value=missing_value, range=(/-1.e9,1.e9/))
           id_overflow_yflux_int_z(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_yflux_overflow_int_z',                &
              Grd%tracer_axes_flux_y(1:2), Time%model_time,                 &
              'vert-integ cp*overflow_yflux*rho*dzt*dxt*temp',              &
              'Watt', missing_value=missing_value, range=(/-1.e9,1.e9/))
       else
           id_overflow(n) = register_diag_field ('ocean_model',      &
                'overflow_'//trim(T_prog(n)%name),                   &
                Grd%tracer_axes(1:3), Time%model_time,               &
                'overflow*rho*dzt*tracer for '//trim(T_prog(n)%name),&
                trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e9,1.e9/))
           id_overflow_xflux_int_z(n) = register_diag_field ('ocean_model',             &
              trim(T_prog(n)%name)//'_xflux_overflow_int_z',                            &
              Grd%tracer_axes_flux_x(1:2), Time%model_time,                             &
              'vert-integ overflow_xflux*dyt*rho*dzt*tracer for '//trim(T_prog(n)%name),&
              'kg/sec', missing_value=missing_value, range=(/-1.e9,1.e9/))
           id_overflow_yflux_int_z(n) = register_diag_field ('ocean_model',             &
              trim(T_prog(n)%name)//'_yflux_overflow_int_z',                            &
              Grd%tracer_axes_flux_y(1:2), Time%model_time,                             &
              'vert-integ overflow_yflux*dxt*rho*dzt*tracer for '//trim(T_prog(n)%name),&
              'kg/sec', missing_value=missing_value, range=(/-1.e9,1.e9/)) 
       endif
    enddo

    if (debug_this_module) then
       write(stdoutunit,*) ' '
       write(stdoutunit,*) '==Global sums from ocean_overflow_mod== '
       call write_timestamp(Time%model_time)
       write(stdoutunit,'(a,es24.17)') &
       'slope_x        = ',mpp_global_sum(Dom%domain2d,slope_x(:,:), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'slope_y        = ',mpp_global_sum(Dom%domain2d,slope_y(:,:), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_slope(1) = ',mpp_global_sum(Dom%domain2d,topog_slope(:,:,1), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_slope(2) = ',mpp_global_sum(Dom%domain2d,topog_slope(:,:,2), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_slope(3) = ',mpp_global_sum(Dom%domain2d,topog_slope(:,:,3), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_slope(4) = ',mpp_global_sum(Dom%domain2d,topog_slope(:,:,4), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_step(1)  = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,1), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_step(2)  = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,2), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_step(3)  = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,3), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_step(4)  = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,4), global_sum_flag)
    endif

  end subroutine ocean_overflow_init
! </SUBROUTINE> NAME="ocean_overflow_init"


!#######################################################################
! <SUBROUTINE NAME="overflow">
!
! <DESCRIPTION>
! Compute thickness and density weighted tracer source [tracer*rho*m/s]
! due to upstream tracer advection in regions where density-driven
! overflows are favorable. 
!
! The MOM implementation of the Campin and Goosse (1999)
! algorithm is detailed in MOM Elements.
!
! </DESCRIPTION>
!
subroutine overflow (Time, Thickness, T_prog, Dens, index_temp, index_salt)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(inout) :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_density_type),     intent(in)    :: Dens
  integer,                      intent(in)    :: index_temp
  integer,                      intent(in)    :: index_salt

  integer, dimension(isd:ied,jsd:jed,4) :: kup
  integer, dimension(isd:ied,jsd:jed,4) :: kdw

  real, dimension(0:nk)              :: source_check
  real, dimension(0:nk)              :: flux
  real, dimension(nk)                :: source_column
  real, dimension(isd:ied,jsd:jed)   :: tmp

  integer :: tau, taum1
  integer :: i, j, k, m, n
  real    :: temp_so, salt_so, press, density_check
  real    :: overflow_thickness, overflow_speed
  real    :: max_overflow_flux

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if(.not. use_this_module) return 

  if(.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_overflow_mod (overflow): module must be initialized')
  endif 

  tau     = Time%tau
  taum1   = Time%taum1
  wrk2_v  = 0.0

  ! compute grid details for cells participating in downslope flow 
  ! note:  "so" = "shallow ocean" cell 
  ! note:  "do" = cells within the "deep ocean" column 

  kup(:,:,:) = 0
  kdw(:,:,:) = 0
  overflow_flux(:,:,:) = 0.0
  overflow_xflux(:,:)  = 0.0
  overflow_yflux(:,:)  = 0.0

  do m=1,4
    do j=jsc,jec
       do i=isc,iec

           max_overflow_flux = Grd%dat(i,j)*p5dtimer

           ! check to see if downslope flow is topographically possible
           if(topog_step(i,j,m) == 1.0) then  

               ! check to see if density of so-cell > density of adjacent do-cell  
               k = Grd%kmt(i,j) 
               if(Dens%rho(i,j,k,tau) > Dens%rho(i+ip(m),j+jq(m),k,tau)) then   

                   ! k-level of shallow-ocn (so) bottom cell
                   kup(i,j,m) = Grd%kmt(i,j) 

                   ! speed (m/sec) of downslope flow 
                   k=kup(i,j,m) 
                   temp_so     = T_prog(index_temp)%field(i,j,k,tau)         
                   salt_so     = Dens%rho_salinity(i,j,k,tau) 
                   overflow_speed = overflow_factor*topog_slope(i,j,m) &
                                    *(Dens%rho(i,j,k,tau)-Dens%rho(i+ip(m),j+jq(m),k,tau))
                   overflow_speed = min(overflow_speed,overflow_umax)

                   ! kdw = k-level in the do-column where so-cell is neutrally buoyant or at bottom
                   kdw(i,j,m) = kup(i,j,m)
                   do k=kup(i,j,m)+1,Grd%kmt(i+ip(m),j+jq(m))
                      press         = Dens%pressure_at_depth(i+ip(m),j+jq(m),k)
                      density_check = density(salt_so,temp_so,press)
                      if(density_check > Dens%rho(i+ip(m),j+jq(m),k,tau)) then 
                          kdw(i,j,m)=k
                      else 
                          exit
                      endif
                   enddo

                   ! compute m^2/sec flux leaving vertical face of cell  
                   if(m==1 .or. m==3) then 
                       overflow_flux(i,j,m)=overflow_speed*Grd%dyt(i,j)
                   elseif(m==2 .or. m==4) then 
                       overflow_flux(i,j,m)=overflow_speed*Grd%dxt(i,j)
                   endif

                   ! check for too large flux. based on 1D diffusive stability criteria, 
                   ! no more than one-half the mass of a grid cell can be mixed per time 
                   ! step, which means
                   ! overflow_flux(m^2/sec)*rho_dzt*dtime < 0.5*rho_dzt*dat, 
                   ! or equivalently overflow_flux(m^2/sec) < 0.5*dat/dtime
                   overflow_flux(i,j,m) = min(overflow_flux(i,j,m),max_overflow_flux) 

                   ! convert from m^2/sec to (kg/sec) using min rho_dzt
                   k=kup(i,j,m)
                   overflow_thickness = min(Thickness%rho_dzt(i,j,k,tau), &
                                            Thickness%rho_dzt(i+ip(m),j+jq(m),k,tau)) 
                   do k=kup(i,j,m),kdw(i,j,m)
                     overflow_thickness = min(overflow_thickness,Thickness%rho_dzt(i+ip(m),j+jq(m),k,tau))
                   enddo
                   overflow_flux(i,j,m) = overflow_flux(i,j,m)*overflow_thickness

                   ! for diagnostics: mass transport (kg/sec) off the shelf
                   if(m==1) overflow_xflux(i,j) = overflow_xflux(i,j) + overflow_flux(i,j,m)
                   if(m==2) overflow_yflux(i,j) = overflow_yflux(i,j) + overflow_flux(i,j,m)
                   if(m==3) overflow_xflux(i,j) = overflow_xflux(i,j) - overflow_flux(i,j,m)
                   if(m==4) overflow_yflux(i,j) = overflow_yflux(i,j) - overflow_flux(i,j,m)

               endif    ! if-test for density difference
           endif      ! if-test for topog_step(i,j,m)==1.0
        enddo      ! i-end  
     enddo      ! j-end
  enddo      ! m-end for the four surrounding cells 

  call mpp_update_domains(overflow_flux(:,:,:),Dom%domain2d)
  call mpp_update_domains(kup(:,:,:),Dom%domain2d)
  call mpp_update_domains(kdw(:,:,:),Dom%domain2d)

  if (debug_this_module) then
     write(stdoutunit,*) ' '
     write(stdoutunit,*) '==Global sums from ocean_overflow_mod== '
     call write_timestamp(Time%model_time)
     write(stdoutunit,'(a,es24.17)')  'overflow_flux(1)(kg/sec) = ',&
       mpp_global_sum(Dom%domain2d,overflow_flux(:,:,1), global_sum_flag)
     write(stdoutunit,'(a,es24.17)')  'overflow_flux(2)(kg/sec) = ',&
       mpp_global_sum(Dom%domain2d,overflow_flux(:,:,2), global_sum_flag)
     write(stdoutunit,'(a,es24.17)')  'overflow_flux(3)(kg/sec) = ',&
       mpp_global_sum(Dom%domain2d,overflow_flux(:,:,3), global_sum_flag)
     write(stdoutunit,'(a,es24.17)')  'overflow_flux(4)(kg/sec) = ',&
       mpp_global_sum(Dom%domain2d,overflow_flux(:,:,4), global_sum_flag)
     write(stdoutunit,'(a,i22)')      'kup(1)           = ',&
       mpp_global_sum(Dom%domain2d,kup(:,:,1), global_sum_flag)
     write(stdoutunit,'(a,i22)')      'kup(2)           = ',&
       mpp_global_sum(Dom%domain2d,kup(:,:,2), global_sum_flag)
     write(stdoutunit,'(a,i22)')      'kup(3)           = ',&
       mpp_global_sum(Dom%domain2d,kup(:,:,3), global_sum_flag)
     write(stdoutunit,'(a,i22)')      'kup(4)           = ',&
       mpp_global_sum(Dom%domain2d,kup(:,:,4), global_sum_flag)
     write(stdoutunit,'(a,i22)')      'kdw(1)           = ',&
       mpp_global_sum(Dom%domain2d,kdw(:,:,1), global_sum_flag)
     write(stdoutunit,'(a,i22)')      'kdw(2)           = ',&
       mpp_global_sum(Dom%domain2d,kdw(:,:,2), global_sum_flag)
     write(stdoutunit,'(a,i22)')      'kdw(3)           = ',&
       mpp_global_sum(Dom%domain2d,kdw(:,:,3), global_sum_flag)
     write(stdoutunit,'(a,i22)')      'kdw(4)           = ',&
       mpp_global_sum(Dom%domain2d,kdw(:,:,4), global_sum_flag)
  endif

  ! compute fluxes and the source for those 
  ! cells participating in downslope transport 

  ! do the mass first
  if (no_return_flow) then

     source_overflow(:,:,:) = 0.0
     flux_int_z(:,:,:)      = 0.0
     wrk1(:,:,:)=0.0 ; wrk2(:,:,:)=0.0 ; wrk3(:,:,:)=0.0 ; wrk4(:,:,:)=0.0
     flux(:) = 0.0

     ! extended do-loop limits perform extra computations, but
     ! save time by reducing need to call to mpp_update_domains.
     do m=1,4
        do j=jsc-ijsc(m),jec+ijec(m)
           do i=isc-ijsc(m),iec+ijec(m)
              
              ! see if (i,j) is a so-cell participating in downslope flow in the m-direction
              if(kup(i,j,m) > 0) then

                 ! initialize upstream mass flux leaving tracer cells
                 flux(0) = 0.0

                 ! initialize thickness*density weighted mass tendency in column
                 source_column(:) = 0.0

                 ! mass flux horizontally leaving so-cell at k=kup and entering do-cell at k=kdw
                 k=kup(i,j,m)
                 flux(0) = overflow_flux(i,j,m)

                 ! for diagnosing mass transport in 3-d space  
                 if(m==1) wrk2_v(i,j,k,1) = wrk2_v(i,j,k,1) + overflow_flux(i,j,m)
                 if(m==2) wrk2_v(i,j,k,2) = wrk2_v(i,j,k,2) + overflow_flux(i,j,m)
                 if(m==3) wrk2_v(i,j,k,1) = wrk2_v(i,j,k,1) - overflow_flux(i,j,m)
                 if(m==4) wrk2_v(i,j,k,2) = wrk2_v(i,j,k,2) - overflow_flux(i,j,m)

                 ! mass source for so-cell at k=kup [(kg/m^3)*(m/s)]
                 source_overflow(i,j,k) = source_overflow(i,j,k) - Grd%datr(i,j)*flux(0)

                 ! mass source for do-cell at k=kdw  [(kg/m^3)*(m/s)]
                 k = kdw(i,j,m)
                 source_column(k) = Grd%datr(i+ip(m),j+jq(m)) * flux(0)

                 if(m==1) then 
                    wrk1(i,j,:) = source_column(:)
                    wrk2_v(i+ip(m),j+jq(m),k,1) = wrk2_v(i+ip(m),j+jq(m),k,1) + overflow_flux(i,j,m)
                 endif 
                 if(m==2) then 
                    wrk2(i,j,:) = source_column(:)
                    wrk2_v(i+ip(m),j+jq(m),k,2) = wrk2_v(i+ip(m),j+jq(m),k,2) + overflow_flux(i,j,m) 
                 endif 
                 if(m==3) then 
                    wrk3(i,j,:) = source_column(:)
                    wrk2_v(i+ip(m),j+jq(m),k,1) = wrk2_v(i+ip(m),j+jq(m),k,1) - overflow_flux(i,j,m) 
                 endif 
                 if(m==4) then 
                    wrk4(i,j,:) = source_column(:)
                    wrk2_v(i+ip(m),j+jq(m),k,2) = wrk2_v(i+ip(m),j+jq(m),k,2) - overflow_flux(i,j,m)
                 endif 

              endif ! kup(i,j,m) > 0
           enddo ! i-end
        enddo ! j-end
     enddo ! m-end for the four directions

     ! each cell can have a mass source five different ways:
     ! 1:   when cell is the so-cell
     ! 2-5: when cell is one of the do-cells for adjacent so-cells.
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              source_overflow(i,j,k)  = source_overflow(i,j,k)   &
                                        + wrk1(i-1,j,k) &
                                        + wrk2(i,j-1,k) &
                                        + wrk3(i+1,j,k) &
                                        + wrk4(i,j+1,k)
              if (source_overflow(i,j,k) .ne. 0.0) then
                 Thickness%mass_source(i,j,k) = Thickness%mass_source(i,j,k) + source_overflow(i,j,k)
              endif
           enddo
        enddo
     enddo

     call mpp_update_domains(Thickness%mass_source(:,:,:),Dom%domain2d)
     
  endif ! no_return_flow?

  ! now do the tracer
  ! tracer source = convergence of upstream advective fluxes 
  wrk1_v(:,:,:,:) = 0.0
  do n=1,num_prog_tracers 

     source_overflow(:,:,:) = 0.0
     flux_int_z(:,:,:)      = 0.0 
     wrk1=0.0 ; wrk2=0.0 ; wrk3=0.0 ; wrk4=0.0

     ! extended do-loop limits perform extra computations, but  
     ! save time by reducing need to call to mpp_update_domains.
     do m=1,4
        do j=jsc-ijsc(m),jec+ijec(m)
           do i=isc-ijsc(m),iec+ijec(m)

              ! see if (i,j) is a so-cell participating in downslope flow in the m-direction
              if(kup(i,j,m) > 0) then 

                  ! initialize upstream tracer flux leaving tracer cells   
                  flux(:) = 0.0

                  ! initialize thickness*density weighted tendency in column 
                  source_column(:) = 0.0

                  ! flux horizontally leaving so-cell at k=kup and entering do-cell at k=kdw 
                  k=kup(i,j,m) 
                  flux(0) = overflow_flux(i,j,m)*T_prog(n)%field(i,j,k,taum1)

                  if (.not. no_return_flow) then
                     ! flux leaving do-cells
                     do k=kup(i,j,m),kdw(i,j,m)
                        flux(k) = overflow_flux(i,j,m)*T_prog(n)%field(i+ip(m),j+jq(m),k,taum1)
                     enddo

                     ! source for so-cell at k=kup [tracer*(kg/m^3)*(m/s)] 
                     k = kup(i,j,m)
                     source_overflow(i,j,k) = source_overflow(i,j,k) + Grd%datr(i,j)*flux(k) 

                     ! source for do-cell at k=kdw  [tracer*(kg/m^3)*(m/s)]
                     k = kdw(i,j,m)
                     source_column(k) = - Grd%datr(i+ip(m),j+jq(m))*flux(k)

                     ! source for do-cell with kup <= k <= kdw-1 [tracer*(kg/m^3)*(m/s)]
                     if(kdw(i,j,m) > kup(i,j,m)) then 
                        do k = kup(i,j,m),kdw(i,j,m)-1
                           source_column(k) = Grd%datr(i+ip(m),j+jq(m))*(flux(k+1)-flux(k)) 
                        enddo
                     endif

                  endif !.not. no_return_flow

                  ! for diagnosing vertically integrated horizontal flux
                  flux_int_z(i,j,m) = flux_int_z(i,j,m) + flux(0)-flux(kup(i,j,m))
                  
                  ! source for so-cell at k=kup [tracer*(kg/m^3)*(m/s)]
                  k = kup(i,j,m)
                  source_overflow(i,j,k) = source_overflow(i,j,k) - Grd%datr(i,j)*flux(0)

                  ! source for do-cell at k=kdw  [tracer*(kg/m^3)*(m/s)]
                  k = kdw(i,j,m)
                  source_column(k) = source_column(k) + Grd%datr(i+ip(m),j+jq(m))*flux(0)
                  
                  if(m==1) wrk1(i,j,:) = source_column(:) 
                  if(m==2) wrk2(i,j,:) = source_column(:) 
                  if(m==3) wrk3(i,j,:) = source_column(:) 
                  if(m==4) wrk4(i,j,:) = source_column(:)


              endif  ! kup(i,j,m) > 0
           enddo ! i-end 
        enddo  ! j-end
     enddo  ! m-end for the four directions 

     ! each cell can have a source five different ways:
     ! 1:   when cell is the so-cell
     ! 2-5: when cell is one of the do-cells for adjacent so-cells.
     do k=1,nk
        do j=jsc,jec   
           do i=isc,iec
              source_overflow(i,j,k)  = source_overflow(i,j,k)   &
                                        + wrk1(i-1,j,k) & 
                                        + wrk2(i,j-1,k) &
                                        + wrk3(i+1,j,k) &
                                        + wrk4(i,j+1,k)
              if (source_overflow(i,j,k) .ne. 0.0) then
                 T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + source_overflow(i,j,k)
              endif
           enddo
        enddo
     enddo

     ! for some diagnostics 
     if(n==index_temp) then 
         wrk1_v(:,:,:,1) = 0.0
         do k=1,nk
            do j=jsc,jec   
               do i=isc,iec
                  wrk1_v(i,j,k,1) = source_overflow(i,j,k)
               enddo
            enddo
         enddo
     endif
     if(n==index_salt) then 
         wrk1_v(:,:,:,2) = 0.0
         do k=1,nk
            do j=jsc,jec   
               do i=isc,iec
                  wrk1_v(i,j,k,2) = source_overflow(i,j,k)
               enddo
            enddo
         enddo
     endif

     if (debug_this_module) then
        source_check(:) = 0.0
        do k=1,nk
           do j=jsc,jec  
             do i=isc,iec
               tmp(i,j) =  T_prog(n)%conversion*dtime*Grd%dat(i,j)*source_overflow(i,j,k)
             enddo
           enddo    
           source_check(k) = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)
        enddo
        do k=1,nk
           source_check(0) = source_check(0) + source_check(k)
        enddo
        if(T_prog(n)%name == 'temp' ) then 
           write(stdoutunit,'(a,i3,a,es24.17)') &
           'For tracer(',n,'), (temp   ) total overflow source is (J) ', source_check(0)
        elseif(T_prog(n)%name == 'salt' ) then 
           write(stdoutunit,'(a,i3,a,es24.17)') &
           'For tracer(',n,'), (salt   ) total overflow source is (kg)', source_check(0)
        else 
           write(stdoutunit,'(a,i3,a,es24.17)') &
           'For tracer(',n,'), (passive) total overflow source is (kg)', source_check(0)
        endif 
     endif

     ! send some diagnostics that are functions of tracer
     if(id_overflow(n) > 0) then 
        call diagnose_3d(Time, Grd, id_overflow(n), T_prog(n)%conversion*source_overflow(:,:,:))
     endif
     if(id_overflow_xflux_int_z(n) > 0) then 
        call diagnose_2d(Time, Grd, id_overflow_xflux_int_z(n),             &
                T_prog(n)%conversion*(flux_int_z(:,:,1)-flux_int_z(:,:,3)))
     endif 
     if(id_overflow_yflux_int_z(n) > 0) then 
        call diagnose_2d(Time, Grd, id_overflow_yflux_int_z(n),             &
                T_prog(n)%conversion*(flux_int_z(:,:,2)-flux_int_z(:,:,4)))
     endif

  enddo  ! n-end for num_prog_tracers 


  ! diagnostics for tracer independent mass transports 
  call diagnose_2d(Time, Grd, id_overflow_xflux, overflow_xflux(:,:))
  call diagnose_2d(Time, Grd, id_overflow_yflux, overflow_yflux(:,:))

  if (id_tx_trans_overflow > 0) then 
     call diagnose_3d(Time, Grd, id_tx_trans_overflow, transport_convert*wrk2_v(:,:,:,1))
  endif
  if (id_ty_trans_overflow > 0) then 
     call diagnose_3d(Time, Grd, id_ty_trans_overflow, transport_convert*wrk2_v(:,:,:,2))
  endif

  call watermass_diag(Time, Dens, wrk1_v)


end subroutine overflow
! </SUBROUTINE> NAME="overflow"



!#######################################################################
! <SUBROUTINE NAME="watermass_diag_init">
!
! <DESCRIPTION>
! Initialization of watermass diagnostic output files. 
! </DESCRIPTION>
!
subroutine watermass_diag_init(Time, Dens)

  type(ocean_time_type),    intent(in) :: Time
  type(ocean_density_type), intent(in) :: Dens

  integer :: stdoutunit
  stdoutunit=stdout()
  compute_watermass_diag = .false. 

  id_neut_rho_overfl = register_diag_field ('ocean_model', 'neut_rho_overfl',&
       Grd%tracer_axes(1:3), Time%model_time,                                &
       'update of locally referenced potrho from overflow scheme',           &
       '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_overfl > 0) compute_watermass_diag = .true. 

  id_neut_rho_overfl_on_nrho = register_diag_field ('ocean_model',                  &
       'neut_rho_overfl_on_nrho',Dens%neutralrho_axes(1:3), Time%model_time,        &
       'update of locally ref potrho from overflow as binned to neutral rho layers',&
       '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_overfl_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_overfl = register_diag_field ('ocean_model',      &
       'wdian_rho_overfl', Grd%tracer_axes(1:3), Time%model_time,&
       'dianeutral mass transport due to overflow',              &
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_overfl > 0) compute_watermass_diag = .true. 

  id_wdian_rho_overfl_on_nrho = register_diag_field ('ocean_model',                &
       'wdian_rho_overfl_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,     &
       'dianeutral mass transport due to overflow as binned to neutral rho layers',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_overfl_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_overfl_on_nrho = register_diag_field ('ocean_model',           &
       'tform_rho_overfl_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,&
       'watermass transform due to overflow as binned to neutral rho layers', &
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_overfl_on_nrho > 0) compute_watermass_diag = .true. 



  id_neut_temp_overfl = register_diag_field ('ocean_model', 'neut_temp_overfl',&
       Grd%tracer_axes(1:3), Time%model_time,                                  &
       'temp related update of locally referenced potrho from overflow scheme',&
       '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_overfl > 0) compute_watermass_diag = .true. 

  id_neut_temp_overfl_on_nrho = register_diag_field ('ocean_model',                           &
       'neut_temp_overfl_on_nrho',Dens%neutralrho_axes(1:3), Time%model_time,                 &
       'temp related update of locally ref potrho from overflow binned to neutral rho layers',&
       '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_overfl_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_temp_overfl = register_diag_field ('ocean_model',      &
       'wdian_temp_overfl', Grd%tracer_axes(1:3), Time%model_time,&
       'temp related dianeutral mass transport due to overflow',  &
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_overfl > 0) compute_watermass_diag = .true. 

  id_wdian_temp_overfl_on_nrho = register_diag_field ('ocean_model',                         &
       'wdian_temp_overfl_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,              &
       'temp related dianeutral mass transport due to overflow binned to neutral rho layers',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_overfl_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_temp_overfl_on_nrho = register_diag_field ('ocean_model',                      &
       'tform_temp_overfl_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,           &
       'temp related watermass transform due to overflow as binned to neutral rho layers',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_temp_overfl_on_nrho > 0) compute_watermass_diag = .true. 


  id_neut_salt_overfl = register_diag_field ('ocean_model', 'neut_salt_overfl',&
       Grd%tracer_axes(1:3), Time%model_time,                                  &
       'salt related update of locally referenced potrho from overflow scheme',&
       '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_overfl > 0) compute_watermass_diag = .true. 

  id_neut_salt_overfl_on_nrho = register_diag_field ('ocean_model',                           &
       'neut_salt_overfl_on_nrho',Dens%neutralrho_axes(1:3), Time%model_time,                 &
       'salt related update of locally ref potrho from overflow binned to neutral rho layers',&
       '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_overfl_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_salt_overfl = register_diag_field ('ocean_model',      &
       'wdian_salt_overfl', Grd%tracer_axes(1:3), Time%model_time,&
       'salt related dianeutral mass transport due to overflow',  &
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_overfl > 0) compute_watermass_diag = .true. 

  id_wdian_salt_overfl_on_nrho = register_diag_field ('ocean_model',                         &
       'wdian_salt_overfl_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,              &
       'salt related dianeutral mass transport due to overflow binned to neutral rho layers',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_overfl_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_salt_overfl_on_nrho = register_diag_field ('ocean_model',                      &
       'tform_salt_overfl_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,           &
       'salt related watermass transform due to overflow as binned to neutral rho layers',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_salt_overfl_on_nrho > 0) compute_watermass_diag = .true. 



  if(compute_watermass_diag) then 
      write(stdoutunit,'(/a/)') &
      '==>Note: running ocean_overflow_mod w/ compute_watermass_diag=.true.'  
  endif

end subroutine watermass_diag_init
! </SUBROUTINE> NAME="watermass_diag_init"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag">
!
! <DESCRIPTION>
! Diagnose effects from overflow on the watermass transformation.
! </DESCRIPTION>
!
subroutine watermass_diag(Time, Dens, wrk1_v)

  type(ocean_time_type),           intent(in)    :: Time
  type(ocean_density_type),        intent(in)    :: Dens
  real, dimension(isd:,jsd:,:,:),  intent(inout) :: wrk1_v

  integer :: i,j,k

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_overflow (watermass_diag): module needs initialization ')
  endif 

  if (.not. compute_watermass_diag) return

  ! full effects on rho 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Grd%tmask(i,j,k)*(Dens%drhodT(i,j,k)*wrk1_v(i,j,k,1)+Dens%drhodS(i,j,k)*wrk1_v(i,j,k,2))
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)  
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_overfl,wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_overfl, wrk3(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_overfl_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_overfl_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_overfl_on_nrho, wrk4)

  ! temperature effects on rho 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Grd%tmask(i,j,k)*Dens%drhodT(i,j,k)*wrk1_v(i,j,k,1)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)  
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_temp_overfl,wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_overfl, wrk3(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_overfl_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_overfl_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_overfl_on_nrho, wrk4)

  ! salinity effects on rho 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Grd%tmask(i,j,k)*Dens%drhodS(i,j,k)*wrk1_v(i,j,k,2)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)  
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_salt_overfl,wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_overfl, wrk3(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_overfl_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_overfl_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_overfl_on_nrho, wrk4)


end subroutine watermass_diag
! </SUBROUTINE>  NAME="watermass_diag"


end module ocean_overflow_mod

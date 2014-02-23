module ocean_overexchange_mod
#define COMP isc:iec,jsc:jec
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Exchange of tracer properties as dense shallow parcel discharged
! into deeper water to approach the depth of neutral buoyancy.  
!</OVERVIEW>
!
!<DESCRIPTION>
! Exchange of tracer properties as dense shallow parcel is discharged
! into deeper water to approach the parcel's depth of neutral buoyancy.
! This module can be characterized as a mixture of the approach from
! Campin and Goosse (1999) and and dynamically determined xlandmix. 
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
!<NAMELIST NAME="ocean_overexchange_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  For using this module.  Default use_this_module=.false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging 
!  </DATA> 
!  <DATA NAME="do_bitwise_exact_sum" TYPE="logical">
!  Set true to do bitwise exact global sum. When it is false, the global
!  sum will be non-bitwise_exact, but will significantly increase efficiency. 
!  The default value is false.
!  </DATA> 
!
!  <DATA NAME="overexch_npts" TYPE="integer">
!  Number of points used in search for the exchange method. 
!  Default overexch_npts=1.  Note: it is not possible to have 
!  overexch_npts greater than or equal to the computational domain
!  extents, as this would require updates across multiple processors. 
!  Default overexch_npts=1. 
!  </DATA> 
!  <DATA NAME="overexch_weight_far" TYPE="logical">
!  To place more weight on points further from central point.  This may 
!  be done to enhance properties getting downslope.  Default is 
!  overexch_weight_far=.false.
!  </DATA> 
!  <DATA NAME="overexch_width" TYPE="integer">
!  Width of the re-weighting function used to emphasize the points further 
!  along in the search for exchange points.  Default overexch_width=1.
!  </DATA> 
!  <DATA NAME="overexch_stability" TYPE="real">
!  Stability factor for determining the maximum overexch_flux.
!  Default overexch_stability=0.25
!  </DATA> 
!  <DATA NAME="overexch_min_thickness" TYPE="real" UNITS="metre">
!  Minimum bottom cell thickness allowed for use of this scheme.
!  Have found that with very thin cells, the model can become very 
!  unstable.  Default overexch_min_thickness=4.0
!  </DATA> 
!  <DATA NAME="overexch_check_extrema" TYPE="logical">
!  Check to be sure there are no global tracer extrema formed due
!  to the overexch process. Note that this approach DOES NOT 
!  conserve tracer, so it is not generally recommended. 
!  Default overexch_check_extrema=.false.
!  </DATA> 
!
!  <DATA NAME="overflow_mu" TYPE="real" UNITS="sec^-1">
!  Dissipation rate for the bottom friction.  Campin and Goosse 
!  suggest overflow_mu=10^-4
!  </DATA> 
!  <DATA NAME="overflow_delta" TYPE="real" UNITS="dimensionless">
!  Fraction of a grid cell participating in the overflow process. 
!  Campin and Goosse suggest overflow_delta=1/3. 
!  </DATA> 
!  <DATA NAME="overflow_umax" TYPE="real" UNITS="m/s">
!  Maximum downslope speed used for determining the exchange rate. 
!  Default overflow_umax=1.0.
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
use mpp_mod,             only: input_nml_file, mpp_error

use ocean_density_mod,     only: density
use ocean_domains_mod,     only: get_local_indices, set_ocean_domain
use ocean_parameters_mod,  only: missing_value, rho0r, rho0, grav
use ocean_parameters_mod,  only: TERRAIN_FOLLOWING, PRESSURE_BASED
use ocean_types_mod,       only: ocean_domain_type, ocean_grid_type, ocean_time_type, ocean_options_type
use ocean_types_mod,       only: ocean_prog_tracer_type, ocean_density_type, ocean_thickness_type
use ocean_util_mod,        only: write_timestamp, diagnose_2d, diagnose_3d
use ocean_util_mod,        only: write_chksum_3d, write_chksum_2d, write_chksum_2d_int
use ocean_tracer_util_mod, only: diagnose_3d_rho
use ocean_workspace_mod,   only: wrk1, wrk2, wrk3, wrk4, wrk1_v

implicit none

private 

public ocean_overexchange_init
public overexchange
private watermass_diag_init
private watermass_diag

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

#include <ocean_memory.h>

integer, dimension(:,:), allocatable :: ip    ! shifting in i-direction
integer, dimension(:,:), allocatable :: jq    ! shifting in j-direction

real, dimension(:,:),     allocatable   :: slope_x          ! topog slope (dimensionless) in the i-direction
real, dimension(:,:),     allocatable   :: slope_y          ! topog slope (dimensionless) in the j-direction 
real, dimension(:,:,:),   allocatable   :: topog_step       ! where downslope flow is topgraphically possible 
real, dimension(:,:,:),   allocatable   :: source_overexch  ! (rho*dzt*tracer) tendency
real, dimension(:,:,:,:), allocatable   :: overexch_flux    ! mass flux associated with overexch 

! fields with larger halos 
integer, dimension(:,:),     allocatable :: kmt_ex
integer, dimension(:,:,:),   allocatable :: kup_ex
integer, dimension(:,:,:,:), allocatable :: kdw_ex
real, dimension(:,:,:),   allocatable :: topog_slope_ex    
real, dimension(:,:,:),   allocatable :: topog_step_ex
real, dimension(:,:,:),   allocatable :: tracer_ex
real, dimension(:,:,:),   allocatable :: temp_ex
real, dimension(:,:,:),   allocatable :: salt_ex
real, dimension(:,:,:),   allocatable :: press_ex
real, dimension(:,:,:),   allocatable :: rho_dzt_ex
real, dimension(:,:,:),   allocatable :: rho_ex
real, dimension(:,:,:),   allocatable :: rho_salinity_ex
real, dimension(:,:),     allocatable :: dxt_ex
real, dimension(:,:),     allocatable :: dyt_ex
real, dimension(:,:),     allocatable :: datr_ex

! internally set for computing watermass diagnostics
logical :: compute_watermass_diag = .false. 


! domain for overexch_exchange 
type(ocean_domain_type), save :: Overexchange_domain 


! for diagnostic manager 
logical :: used
integer :: id_slope_x      =-1
integer :: id_slope_y      =-1
integer :: id_topog_step_1 =-1
integer :: id_topog_step_2 =-1
integer :: id_topog_step_3 =-1
integer :: id_topog_step_4 =-1

integer :: id_neut_rho_overex          =-1
integer :: id_wdian_rho_overex         =-1
integer :: id_neut_rho_overex_on_nrho  =-1
integer :: id_wdian_rho_overex_on_nrho =-1
integer :: id_tform_rho_overex_on_nrho =-1

integer :: id_neut_temp_overex          =-1
integer :: id_wdian_temp_overex         =-1
integer :: id_neut_temp_overex_on_nrho  =-1
integer :: id_wdian_temp_overex_on_nrho =-1
integer :: id_tform_temp_overex_on_nrho =-1

integer :: id_neut_salt_overex          =-1
integer :: id_wdian_salt_overex         =-1
integer :: id_neut_salt_overex_on_nrho  =-1
integer :: id_wdian_salt_overex_on_nrho =-1
integer :: id_tform_salt_overex_on_nrho =-1

integer, dimension(:), allocatable :: id_overexch


character(len=128) :: version=&
       '=>Using: ocean_overexchange.f90 ($Id: ocean_overexchange.F90,v 20.0 2013/12/14 00:16:04 fms Exp $)'
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

! halo 
integer :: ijhalo=1 

! grav*overflow_delta/(rho0*overflow_mu)  (m^4/kg/sec)
real :: overexch_factor 

! time step 
real :: dtime 

! for setting max overexch_flux 
real :: overexch_stability_rate

! set from nml
logical :: use_this_module        = .false. ! must be set .true. in nml to enable this scheme
logical :: debug_this_module      = .false. ! for debugging
logical :: do_bitwise_exact_sum   = .false. ! set true to get slower sum that is same for all PEs. 
logical :: overexch_weight_far    = .false. ! to place more weight on points further from central point
logical :: overexch_check_extrema = .false. ! to be sure there are no global extrema formed from overexch
integer :: overexch_npts          = 1       ! number points searched for the exchange 
integer :: overexch_width         = 1       ! width for exponential that damps points near to central point
real    :: overflow_mu            = 1.e-4   ! frictional dissipation rate (sec^-1) at bottom  
real    :: overflow_delta         = 0.3333  ! fraction of a grid cell participating in overflow 
real    :: overflow_umax          = 1.0     ! maximum downslope speed (m/s) for setting exchange rate 
real    :: overexch_stability     = 0.25    ! stability factor for setting the max overexch_flux 
real    :: overexch_min_thickness = 4.0     ! metres dzt needed in order for a cell to participate in exchange

namelist /ocean_overexchange_nml/ use_this_module, debug_this_module, overexch_npts, overexch_width, &
                                  overexch_weight_far, overflow_mu, overflow_delta, overflow_umax,   &
                                  do_bitwise_exact_sum, overexch_stability, overexch_min_thickness,  &
                                  overexch_check_extrema

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_overexchange_init">
!
! <DESCRIPTION>
! Initial set up for mixing of tracers into the abyss next to topography.
! </DESCRIPTION>
!
  subroutine ocean_overexchange_init(Grid, Domain, Time, Dens, T_prog, Ocean_options, &
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
      '==>Error from ocean_overexchange_mod (ocean_overexchange_init): module already initialized')
    endif 

    module_is_initialized = .TRUE.

    call write_version_number( version, tagname )
#ifndef MOM_STATIC_ARRAYS
    call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
    nk = Grid%nk
#endif
    Dom => Domain
    Grd => Grid

    dtime = dtim 
    num_prog_tracers = size(T_prog(:))

    ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_overexchange_nml, iostat=io_status) 
    ierr = check_nml_error(io_status,'ocean_overexchange_nml')
#else
    ioun =  open_namelist_file()
    read (ioun,ocean_overexchange_nml,IOSTAT=io_status)
    ierr = check_nml_error(io_status,'ocean_overexchange_nml')
    call close_file (ioun)
#endif
    write (stdlogunit,ocean_overexchange_nml)
    write (stdoutunit,'(/)') 
    write (stdoutunit,ocean_overexchange_nml)

    if(do_bitwise_exact_sum) then
       global_sum_flag = BITWISE_EXACT_SUM
    else
       global_sum_flag = NON_BITWISE_EXACT_SUM
    endif

    if (PRESENT(debug) .and. .not. debug_this_module) then
       debug_this_module = debug
    endif 
    if(debug_this_module) then 
       write(stdoutunit,'(a)') '==>Note: running ocean_overexchange_mod with debug_this_module=.true.'  
    endif 

    if(.not. use_this_module) then 
      call mpp_error(NOTE,&
      '==>From ocean_overexchange_mod: NOT using overflow exchange scheme.')
      Ocean_options%overexchange = 'Did NOT use overflow exchange scheme.'
      return
    else 
      if(vert_coordinate_type == TERRAIN_FOLLOWING) then 
          call mpp_error(FATAL, &
          '==>ocean_overexchange_mod: this module is NOT for use with TERRAIN_FOLLOWING vert coodinates.')
      endif
      Ocean_options%overexchange = 'Used the overflow exchange scheme. This scheme is experimental. User beware!'
      call mpp_error(NOTE,&
      '==>From ocean_overexchange_mod: USING overflow exchange scheme.')
    endif
    if(debug_this_module) then 
      call mpp_error(NOTE,'==>From ocean_overexchange_mod: USING debug_this_module')
    endif 

    ! for setting the maximum overexch_flux 
    overexch_stability_rate = overexch_stability/dtime

    ! allocate for tracer source arrays 
    allocate(source_overexch(isd:ied,jsd:jed,nk))
    source_overexch(:,:,:) = 0.0

    if(overexch_npts<1) then
      call mpp_error(FATAL, &
      '==>ocean_overexchange_mod: overexch_npts < 1 not allowed. In nml, set overexch_npts>=1.')
    endif  
    write(stdoutunit,'(a,i3)')' In ocean_overexchange_mod: overexch_npts = ',overexch_npts  
    write(stdoutunit,'(a)')   ' Be sure this number is smaller than dimensions of computational domain.'

    ! halo size for extended domain and larger arrays 
    ijhalo=overexch_npts

    ! k-level of central point in search 
    allocate(kup_ex(isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,4))
    kup_ex(:,:,:) = 0

    ! k-level of deeper points in horizontally adjacent columns 
    allocate(kdw_ex(isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,ijhalo,4))
    kdw_ex(:,:,:,:) = 0

    ! mass flux setting strength of tracer exchange 
    allocate(overexch_flux(isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,ijhalo,4))
    overexch_flux(:,:,:,:) = 0.0

    ! define extended domain 
    call set_ocean_domain(Overexchange_domain,Grd,xhalo=ijhalo,yhalo=ijhalo,name='overexch',maskmap=Dom%maskmap)

    ! time independent arrays defined over the Overexch_domain region 
    allocate(kmt_ex (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo))
    allocate(dxt_ex (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo))
    allocate(dyt_ex (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo))
    allocate(datr_ex(isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo))
    kmt_ex(:,:)     = 0
    datr_ex(:,:)    = 0.0
    dxt_ex(:,:)     = 0.0
    dyt_ex(:,:)     = 0.0
    kmt_ex (isc:iec,jsc:jec) = Grd%kmt (isc:iec,jsc:jec)
    datr_ex(isc:iec,jsc:jec) = Grd%datr(isc:iec,jsc:jec)
    dxt_ex (isc:iec,jsc:jec) = Grd%dxt (isc:iec,jsc:jec)
    dyt_ex (isc:iec,jsc:jec) = Grd%dyt (isc:iec,jsc:jec)
    call mpp_update_domains (kmt_ex(:,:), Overexchange_domain%domain2d)
    call mpp_update_domains (datr_ex(:,:),Overexchange_domain%domain2d)
    call mpp_update_domains (dxt_ex(:,:), Overexchange_domain%domain2d)
    call mpp_update_domains (dyt_ex(:,:), Overexchange_domain%domain2d)

    ! time dependent arrays defined over the Overexch_domain region 
    allocate(temp_ex         (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,nk))
    allocate(salt_ex         (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,nk))
    allocate(tracer_ex       (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,nk))
    allocate(press_ex        (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,nk))
    allocate(rho_salinity_ex (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,nk))
    allocate(rho_ex          (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,nk))
    allocate(rho_dzt_ex      (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,nk))
    temp_ex(:,:,:)         = 0.0
    salt_ex(:,:,:)         = 0.0
    tracer_ex(:,:,:)       = 0.0
    press_ex(:,:,:)        = 0.0
    rho_salinity_ex(:,:,:) = 0.0
    rho_ex(:,:,:)          = 0.0
    rho_dzt_ex(:,:,:)      = 0.0

    ! compute a common factor that is time independent 
    overexch_factor = grav*overflow_delta*rho0r/(epsln+overflow_mu)

    ! compute topographic slope arrays for the i-slope and j-slope
    ! slopes are centered on the i-face and j-face of tracer cells 
    allocate (slope_x(isd:ied,jsd:jed))
    allocate (slope_y(isd:ied,jsd:jed))
    slope_x = 0.0 
    slope_y = 0.0
    do j=jsc,jec
      do i=isc,iec
        slope_x(i,j) = (Grd%ht(i+1,j)-Grd%ht(i,j))*Grd%dxter(i,j)*Grd%tmask(i,j,1)*Grd%tmask(i+1,j,1)
        slope_y(i,j) = (Grd%ht(i,j+1)-Grd%ht(i,j))*Grd%dytnr(i,j)*Grd%tmask(i,j,1)*Grd%tmask(i,j+1,1)
        slope_x(i,j) = abs(slope_x(i,j))
        slope_y(i,j) = abs(slope_y(i,j))
      enddo
    enddo 
    call mpp_update_domains(slope_x(:,:),Dom%domain2d)
    call mpp_update_domains(slope_y(:,:),Dom%domain2d)

    ! topographic slope for the four surrounding directions 
    allocate (topog_slope_ex(isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,4))
    topog_slope_ex(:,:,:) = 0.0
    do j=jsc,jec
       do i=isc,iec
          m=1 ; topog_slope_ex(i,j,m) = slope_x(i,j)  
          m=2 ; topog_slope_ex(i,j,m) = slope_y(i,j)  
          m=3 ; topog_slope_ex(i,j,m) = slope_x(i-1,j)
          m=4 ; topog_slope_ex(i,j,m) = slope_y(i,j-1)
       enddo
    enddo
    call mpp_update_domains (topog_slope_ex(:,:,:), Overexchange_domain%domain2d)

    ! compute directions from an (i,j) point where topography deepens.
    ! these directions may potentially have downslope exchange. 
    ! insist that downslope exchange occurs only when there are more 
    ! kmt cells in the adjacent column. also insist that downslope 
    ! exchange does not involve k=1 cells. 
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
    ! adjacent columns of tracer points.  When the column straddles 
    ! the bipolar fold, the code logic is not general and so it actually 
    ! attempts to reach to a non-adjacent column.  Shy of generalizing 
    ! the logic, we simply do not consider overexch for points along fold. 
    if(jec+Dom%joff==Dom%jeg) then 
       topog_step(:,jec,:) = 0.0
    endif 

    ! fill the larger array 
    allocate (topog_step_ex(isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,4))
    topog_step_ex(:,:,:) = 0.0
    topog_step_ex(isc:iec,jsc:jec,:) = topog_step(isc:iec,jsc:jec,:)
    call mpp_update_domains (topog_step_ex(:,:,:), Overexchange_domain%domain2d)

    ! labels for the four quadranats surrounding a 
    ! central tracer cell (moving counter-clockwise)
    allocate(ip(0:ijhalo,4))
    allocate(jq(0:ijhalo,4))

    do n=0,ijhalo
        m=1 
        ip(n,m)   = 1+(n-1)
        jq(n,m)   = 0  

        m=2 
        ip(n,m)   = 0  
        jq(n,m)   = 1+(n-1)  

        m=3 
        ip(n,m)   = -1-(n-1) 
        jq(n,m)   = 0  

        m=4 
        ip(n,m)   =  0  
        jq(n,m)   = -1-(n-1)
    enddo

    ! register/send diagnostics 
    id_slope_x = register_static_field ('ocean_model', 'slope_x', Grd%tracer_axes_flux_x(1:2), &
                 '|d(ht)/dx| on T-cell face', 'm/m', missing_value=missing_value, range=(/-1.e9,1.e9/))
    call diagnose_2d(Time, Grd, id_slope_x, slope_x(:,:))

    id_slope_y = register_static_field ('ocean_model', 'slope_y', Grd%tracer_axes_flux_y(1:2), &
                 '|d(ht)/dy| on T-cell face', 'm/m', missing_value=missing_value, range=(/-1.e9,1.e9/))
    call diagnose_2d(Time, Grd, id_slope_y, slope_y(:,:))

    id_topog_step_1 = register_static_field ('ocean_model', 'topog_step_1', Grd%tracer_axes(1:2), &
                 'topog_step_1', 'dimensionless', missing_value=missing_value, range=(/-1.0,1.0/))
    call diagnose_2d(Time, Grd, id_topog_step_1, topog_step(1,:,:))

    id_topog_step_2 = register_static_field ('ocean_model', 'topog_step_2', Grd%tracer_axes(1:2), &
                 'topog_step_2', 'dimensionless', missing_value=missing_value, range=(/-1.0,1.0/))
    call diagnose_2d(Time, Grd, id_topog_step_2, topog_step(2,:,:))

    id_topog_step_3 = register_static_field ('ocean_model', 'topog_step_3', Grd%tracer_axes(1:2), &
                 'topog_step_3', 'dimensionless', missing_value=missing_value, range=(/-1.0,1.0/))
    call diagnose_2d(Time, Grd, id_topog_step_3, topog_step(3,:,:))

    id_topog_step_4 = register_static_field ('ocean_model', 'topog_step_4', Grd%tracer_axes(1:2), &
                 'topog_step_4', 'dimensionless', missing_value=missing_value, range=(/-1.0,1.0/))
    call diagnose_2d(Time, Grd, id_topog_step_4, topog_step(4,:,:))

    call watermass_diag_init(Time,Dens)
   
    allocate (id_overexch(num_prog_tracers))
    id_overexch = -1

    do n=1,num_prog_tracers
       if(T_prog(n)%name == 'temp') then 
           id_overexch(n) = register_diag_field ('ocean_model', &
                'overexch_'//trim(T_prog(n)%name),              &
                Grd%tracer_axes(1:3), Time%model_time,          &
                'cp*overexch*rho*dzt*temp',                     &
                'Watt/m^2', missing_value=missing_value, range=(/-1.e9,1.e9/))
       else
           id_overexch(n) = register_diag_field ('ocean_model', 'overexch_'//trim(T_prog(n)%name), &
                Grd%tracer_axes(1:3), Time%model_time,               &
                'overexch*rho*dzt*tracer for '//trim(T_prog(n)%name),&
                trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e9,1.e9/))
       endif
    enddo

    if (debug_this_module) then
       write(stdoutunit,*) ' '
       write(stdoutunit,*) '==Global sums from initialization of ocean_overexchange_mod== '
       call write_timestamp(Time%model_time)
       write(stdoutunit,'(a,es24.17)') &
       'slope_x          = ',mpp_global_sum(Dom%domain2d,slope_x(:,:), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'slope_y          = ',mpp_global_sum(Dom%domain2d,slope_y(:,:), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_slope_ex(1)= ',mpp_global_sum(Overexchange_domain%domain2d,topog_slope_ex(:,:,1), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_slope_ex(2)= ',mpp_global_sum(Overexchange_domain%domain2d,topog_slope_ex(:,:,2), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_slope_ex(3)= ',mpp_global_sum(Overexchange_domain%domain2d,topog_slope_ex(:,:,3), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_slope_ex(4)= ',mpp_global_sum(Overexchange_domain%domain2d,topog_slope_ex(:,:,4), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_step(1)    = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,1), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_step(2)    = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,2), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_step(3)    = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,3), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_step(4)    = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,4), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_step_ex(1) = ',&
        mpp_global_sum(Overexchange_domain%domain2d,topog_step_ex(:,:,1), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_step_ex(2) = ',&
        mpp_global_sum(Overexchange_domain%domain2d,topog_step_ex(:,:,2), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_step_ex(3) = ',&
        mpp_global_sum(Overexchange_domain%domain2d,topog_step_ex(:,:,3), global_sum_flag)
       write(stdoutunit,'(a,es24.17)') &
       'topog_step_ex(4) = ',&
        mpp_global_sum(Overexchange_domain%domain2d,topog_step_ex(:,:,4), global_sum_flag)
    endif

  end subroutine ocean_overexchange_init
! </SUBROUTINE> NAME="ocean_overexchange_init"


!#######################################################################
! <SUBROUTINE NAME="overexchange">
!
! <DESCRIPTION>
! Compute thickness and density weighted tracer source [tracer*rho*m/s]
! due to exchange of tracer properties in regions where density-driven
! overflows are favorable. Allow for exchanges to occur over horizontally 
! distant points, so long as the dense shallow parcel finds that it 
! will sit on the bottom of the horizontally adjacent columns.  Doing
! so requires a search algorithm, which requires some if-test logic
! as well as extended halos.  Note that the halos cannot be extended
! to larger than the size of the computational domain on a processor.
! This restriction limits the extent that we can search horizontally.
!
! This scheme can be characterized as a dynamical xlandmix based on 
! the scheme of Campin and Goosse.  The rates for the exchange are 
! functions of the topographic slope and the density differences
! between parcels.  
!
! </DESCRIPTION>
!
subroutine overexchange (Time, Thickness, T_prog, Dens, index_temp, index_salt)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_density_type),     intent(in)    :: Dens
  integer,                      intent(in)    :: index_temp
  integer,                      intent(in)    :: index_salt

  integer :: tau, taum1
  integer :: i, j, k, m, n, nt, npts 
  integer :: kp1, km1
  integer :: kmtij, ku, kd
  integer :: iip0, iip1, jjq0, jjq1
  integer :: iip1r, jjq1r
  integer :: nm1

  real    :: temp_so, salt_so, press, density_check
  real    :: overexch_thickness, overexch_speed
  real    :: weight, arg
  real    :: delta, delta_rho(ijhalo)
  real    :: max_overexch_flux
  real    :: tnew

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return 

  if(.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_overexchange_mod (overexchange): module must be initialized')
  endif 

  tau     = Time%tau
  taum1   = Time%taum1

  delta_rho(:) = 0.0

  ! extend some fields to larger domain 
  rho_salinity_ex = 0.0
  rho_ex          = 0.0
  rho_dzt_ex      = 0.0
  press_ex        = 0.0
  temp_ex         = 0.0
  salt_ex         = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           rho_ex(i,j,k)          = Dens%rho(i,j,k,taum1) 
           rho_salinity_ex(i,j,k) = Dens%rho_salinity(i,j,k,taum1) 
           rho_dzt_ex(i,j,k)      = Thickness%rho_dzt(i,j,k,taum1) 
           press_ex(i,j,k)        = Dens%pressure_at_depth(i,j,k)           
           temp_ex(i,j,k)         = T_prog(index_temp)%field(i,j,k,taum1)           
           salt_ex(i,j,k)         = T_prog(index_salt)%field(i,j,k,taum1)           
        enddo
     enddo
  enddo

  call mpp_update_domains(rho_ex(:,:,:),          Overexchange_domain%domain2d, complete=.false.)
  call mpp_update_domains(rho_salinity_ex(:,:,:), Overexchange_domain%domain2d, complete=.false.)
  call mpp_update_domains(rho_dzt_ex(:,:,:),      Overexchange_domain%domain2d, complete=.false.)
  call mpp_update_domains(press_ex(:,:,:),        Overexchange_domain%domain2d, complete=.false.)
  call mpp_update_domains(temp_ex(:,:,:),         Overexchange_domain%domain2d, complete=.false.)
  call mpp_update_domains(salt_ex(:,:,:),         Overexchange_domain%domain2d, complete=.true.)

  ! compute grid details for cells participating in downslope exch 
  ! note:  "so" = "shallow ocean" cell...the central point  
  ! note:  "do" = cells within "deep ocean" columns 
  ! this part of the computation is independent of the tracer 

  kup_ex(:,:,:)          = 0
  kdw_ex(:,:,:,:)        = 0
  overexch_flux(:,:,:,:) = 0.0

  do j=jsc,jec
     do i=isc,iec

        kmtij = kmt_ex(i,j)
        if(kmtij > 1) then   

            temp_so           = temp_ex(i,j,kmtij)
            salt_so           = rho_salinity_ex(i,j,kmtij)
            max_overexch_flux = Grd%dat(i,j)*overexch_stability_rate

            ! 4-directions surrounding each cell 
            do m=1,4

               ! search for density favorable exchanges
               npts=0 
               nloop1_npts: do n=1,overexch_npts

                  iip0 = i+ip(n-1,m)
                  jjq0 = j+jq(n-1,m)
                  iip1 = i+ip(n,m)
                  jjq1 = j+jq(n,m)

                  ! check if downslope exchange is topographically possible
                  if(topog_step_ex(iip0,jjq0,m)==1.0) then  

                      ! check if density of shallow central cell > density of do-cell  
                      if(rho_ex(i,j,kmtij) > rho_ex(iip1,jjq1,kmtij)) then   

                          ! kdw = k-level in do-column where central cell is neutrally buoyant or at the bottom 
                          kdw_ex(i,j,n,m) = 0
                          kloop1 : do k=kmtij+1,kmt_ex(iip1,jjq1)
                             press         = press_ex(iip1,jjq1,k)
                             density_check = density(salt_so,temp_so,press)
                             if(density_check > rho_ex(iip1,jjq1,k)) then 
                                 kdw_ex(i,j,n,m) = k
                                 delta_rho(n)    = density_check-rho_ex(iip1,jjq1,k) 
                             else 
                                 exit kloop1
                             endif
                          enddo kloop1

                          ! check that n-parcel has a k-level greater than the n-1 parcel.
                          ! if not, then do not perform an exchange with the n-column.
                          nm1 = max(n-1,1)
                          if(n > 1) then 
                              if(kdw_ex(i,j,n,m) <= kdw_ex(i,j,nm1,m)) then  
                                  kdw_ex(i,j,n,m) = 0
                              endif
                          endif

                          ! set strength of downslope exchange  
                          if(kdw_ex(i,j,n,m) > 0) then 

                              ! add this cell to the number of cells for exchange 
                              npts=npts+1

                              if(n==1) then 
                                  kup_ex(i,j,m) = kmtij
                                  delta = rho_ex(i,j,kmtij)-rho_ex(iip1,jjq1,kmtij)
                              else 
                                  delta = delta_rho(n)
                              endif

                              ku = kup_ex(i,j,m) 
                              kd = kdw_ex(i,j,n,m) 

                              ! convert overexch_speed (m/s) to overexch_flux (kg/sec) 
                              ! speed (m/sec) setting strength of downslope exchange  
                              overexch_speed = overexch_factor*topog_slope_ex(iip0,jjq0,m)*delta 
                              overexch_speed = min(overexch_speed,overflow_umax)

                              if(m==1 .or. m==3) then 
                                  overexch_flux(i,j,n,m)=overexch_speed*Grd%dyt(i,j)
                              elseif(m==2 .or. m==4) then 
                                  overexch_flux(i,j,n,m)=overexch_speed*Grd%dxt(i,j)
                              endif

                              ! check for too large flux. based on 1D diffusive stability criteria, 
                              ! no more than one-half the mass of a grid cell can be mixed per time 
                              ! step, which means
                              ! overexch_flux(m^2/sec)*rho_dzt*dtime < 0.5*rho_dzt*dat, 
                              ! or equivalently overexch_flux(m^2/sec) < 0.5*dat/dtime
                              overexch_flux(i,j,n,m) = min(overexch_flux(i,j,n,m),max_overexch_flux) 

                              ! convert from m^2/sec to kg/sec using minimum rho_dzt 
                              overexch_thickness     = min(rho_dzt_ex(i,j,ku),rho_dzt_ex(iip0,jjq0,kd)) 
                              overexch_flux(i,j,n,m) = overexch_flux(i,j,n,m)*overexch_thickness

                              ! mask out regions with too small bottom cells
                              if(rho0r*overexch_thickness < overexch_min_thickness) then 
                                  overexch_flux(i,j,n,m) = 0.0
                              endif


                          endif  ! kdw_ex(i,j,n,m) > 0 if-test 

                      endif      ! rho_ex(i,j,k) > rho_ex(iip1,jjq1,k) if-test 
                  endif          ! topog_step_ex(iip0,jjq0,m)==1.0 if-test


                  ! if kdw is not on the bottom, then exit n-loop since finished with search 
                  if(kdw_ex(i,j,n,m) < kmt_ex(iip1,jjq1)) then 
                      exit nloop1_npts
                  endif

               enddo nloop1_npts ! n-loop for search reaching out from central point

               ! place more weight on the farther points to 
               ! encourage tracer exchange further downslope.  
               if(overexch_weight_far .and. npts > 0) then 
                   weight = overexch_flux(i,j,1,m)
                   do n=1,npts
                      arg = float((n-npts)/overexch_width)
                      overexch_flux(i,j,n,m) = weight*exp(arg)
                   enddo
               endif

               ! if no exchange points 
               if(npts==0) then 
                   overexch_flux(i,j,:,m) = 0.0
               endif

            enddo  ! m-loop

        endif ! kmtij > 1

     enddo     ! i-loop
  enddo        ! j-loop


  ! extend arrays to wide halos 
  call mpp_update_domains(kup_ex(:,:,:),          Overexchange_domain%domain2d)
  call mpp_update_domains(kdw_ex(:,:,:,:),        Overexchange_domain%domain2d)
  call mpp_update_domains(overexch_flux(:,:,:,:), Overexchange_domain%domain2d)


  ! compute tracer tendency for cells participating in downslope exchange 
  do nt=1,num_prog_tracers 

     ! place tracer concentration in the wider array tracer_ex 
     if(nt==index_temp) then 
         do k=1,nk
            do j=jsc-ijhalo,jec+ijhalo
               do i=isc-ijhalo,iec+ijhalo
                  tracer_ex(i,j,k) = temp_ex(i,j,k)
               enddo
            enddo
         enddo
     elseif(nt==index_salt) then 
         do k=1,nk
            do j=jsc-ijhalo,jec+ijhalo
               do i=isc-ijhalo,iec+ijhalo
                  tracer_ex(i,j,k) = salt_ex(i,j,k)
               enddo
            enddo
         enddo
     else
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  tracer_ex(i,j,k) = T_prog(nt)%field(i,j,k,taum1)
               enddo
            enddo
         enddo
         call mpp_update_domains (tracer_ex(:,:,:), Overexchange_domain%domain2d)
     endif

     ! initialize array holding tendency for (rho*dzt*tracer) 
     source_overexch(:,:,:) = 0.0

     ! compute tendency, noting that each (i,j,k) cell can generally be part of 
     ! an exchange as either a central shallow cell, or as one of the deep cells. 
     do j=jsc,jec
        do i=isc,iec
           do m=1,4
              do n=1,overexch_npts

                 ! source for shallow cell 
                 iip1 = i+ip(n,m)
                 jjq1 = j+jq(n,m)
                 ku   = kup_ex(i,j,m) 
                 kd   = kdw_ex(i,j,n,m) 
                 if(ku > 0 .and. kd > ku) then 
                     source_overexch(i,j,ku) = source_overexch(i,j,ku) &
                          + overexch_flux(i,j,n,m)*datr_ex(i,j)*(tracer_ex(iip1,jjq1,kd)-tracer_ex(i,j,ku))
                 endif

                 ! source for deep cell  
                 iip1r = i-ip(n,m)
                 jjq1r = j-jq(n,m)
                 ku    = kup_ex(iip1r,jjq1r,m)  
                 kd    = kdw_ex(iip1r,jjq1r,n,m)  
                 if(ku > 0 .and. kd > ku) then 
                     source_overexch(i,j,kd) = source_overexch(i,j,kd) &
                          + overexch_flux(iip1r,jjq1r,n,m)*datr_ex(i,j)*(tracer_ex(iip1r,jjq1r,ku)-tracer_ex(i,j,kd))
                 endif

              enddo ! n-loop 
           enddo    ! m-loop
        enddo       ! i-loop
     enddo          ! j-loop


     ! ensure that this scheme does not create any global extrema
     if(overexch_check_extrema) then 
         do k=1,nk
            kp1=min(nk,k+1)
            km1=max(1,k-1)
            do j=jsc,jec   
               do i=isc,iec
                  if(Grd%tmask(i,j,k) > 0.0) then  
                      tnew = tracer_ex(i,j,k) + source_overexch(i,j,k)*dtime/rho_dzt_ex(i,j,k)
                      if(tnew > T_prog(nt)%max_tracer_limit .or. tnew < T_prog(nt)%min_tracer_limit) then 
                         source_overexch(i,j,k) = 0.0
                      endif 
                  endif
               enddo
            enddo
         enddo
     endif

     ! fill tracer tendency array 
     do k=1,nk
        do j=jsc,jec   
           do i=isc,iec
              T_prog(nt)%th_tendency(i,j,k) = T_prog(nt)%th_tendency(i,j,k) + source_overexch(i,j,k)
           enddo
        enddo
     enddo


     if(id_overexch(nt) > 0) then 
        call diagnose_3d(Time, Grd, id_overexch(nt), T_prog(nt)%conversion*source_overexch(:,:,:))
     endif

     if(nt==index_temp) then 
         wrk1_v(:,:,:,1) = 0.0
         do k=1,nk
            do j=jsc,jec   
               do i=isc,iec
                  wrk1_v(i,j,k,1) = source_overexch(i,j,k)
               enddo
            enddo
         enddo
     endif
     if(nt==index_salt) then 
         wrk1_v(:,:,:,2) = 0.0
         do k=1,nk
            do j=jsc,jec   
               do i=isc,iec
                  wrk1_v(i,j,k,2) = source_overexch(i,j,k)
               enddo
            enddo
         enddo
     endif

  enddo  ! nt-end for num_prog_tracers 

  call watermass_diag(Time, Dens, wrk1_v)

  if (debug_this_module) then
      write(stdoutunit,*) ' '
      write(stdoutunit,*) '==Global sums from ocean_overexchange_mod== '
      call write_timestamp(Time%model_time)

      write(stdoutunit,'(a,es24.17)') &
           'rho_ex     = ',&
           mpp_global_sum(Overexchange_domain%domain2d,rho_ex(:,:,:), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
           'rho_dzt_ex = ',&
           mpp_global_sum(Overexchange_domain%domain2d,rho_dzt_ex(:,:,:), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
           'press_ex   = ',&
           mpp_global_sum(Overexchange_domain%domain2d,press_ex(:,:,:), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
           'temp_ex     = ',&
           mpp_global_sum(Overexchange_domain%domain2d,temp_ex(:,:,:), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
           'salt_ex     = ',&
           mpp_global_sum(Overexchange_domain%domain2d,salt_ex(:,:,:), global_sum_flag)
      write(stdoutunit,'(a,i24)') &
           'kmt_ex     = ',&
           mpp_global_sum(Overexchange_domain%domain2d,kmt_ex(:,:), global_sum_flag)

      do m=1,4
         write(stdoutunit,'(a,i24)') &
              'kup_ex            = ',  &
              mpp_global_sum(Overexchange_domain%domain2d,kup_ex(:,:,m), global_sum_flag)
         do n=1,overexch_npts
            write(stdoutunit,'(a,es24.17)') &
                 'overexch_flux     = ',  &
                 mpp_global_sum(Overexchange_domain%domain2d,overexch_flux(:,:,n,m), global_sum_flag)
            write(stdoutunit,'(a,i24)') &
                 'kdw_ex            = ',  &
                 mpp_global_sum(Overexchange_domain%domain2d,kdw_ex(:,:,n,m), global_sum_flag)
         enddo
      enddo

      write(stdoutunit,*) ' '
      write(stdoutunit,*) '==Global check sums from ocean_overexchange_mod== '
      call write_timestamp(Time%model_time)
      call write_chksum_3d('rho_ex', rho_ex(COMP,:))
      call write_chksum_3d('rho_dzt_ex', rho_dzt_ex(COMP,:))
      call write_chksum_3d('press_ex', press_ex(COMP,:))
      call write_chksum_3d('temp_ex', temp_ex(COMP,:))
      call write_chksum_3d('salt_ex', salt_ex(COMP,:))
      call write_chksum_2d_int('kmt_ex;', kmt_ex(COMP))
      do m=1,4
         call write_chksum_2d_int('kup_ex', kup_ex(COMP,m))
         do n=1,overexch_npts
            call write_chksum_2d_int('kdw_ex', kdw_ex(COMP,n,m))
            call write_chksum_2d('overexch_flux', overexch_flux(COMP,n,m))
         enddo
      enddo

  endif

end subroutine overexchange
! </SUBROUTINE> NAME="overexchange"


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
  
  id_neut_rho_overex = register_diag_field ('ocean_model', 'neut_rho_overex',&
    Grd%tracer_axes(1:3), Time%model_time,                                   &
    'update of locally referenced potrho from overexchange scheme',          &
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_overex > 0) compute_watermass_diag = .true.

  id_neut_rho_overex_on_nrho = register_diag_field ('ocean_model',                   &
    'neut_rho_overex_on_nrho',Dens%neutralrho_axes(1:3), Time%model_time,            &
    'update of locally ref potrho from overexchange as binned to neutral rho layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_overex_on_nrho > 0) compute_watermass_diag = .true.

  id_wdian_rho_overex = register_diag_field ('ocean_model',   &
    'wdian_rho_overex', Grd%tracer_axes(1:3), Time%model_time,&
    'dianeutral mass transport due to overexchange',          &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_overex > 0) compute_watermass_diag = .true.

  id_wdian_rho_overex_on_nrho = register_diag_field ('ocean_model',                 &
    'wdian_rho_overex_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,         &
    'dianeutral mass transport due to overexchange as binned to neutral rho layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_overex_on_nrho > 0) compute_watermass_diag = .true.

  id_tform_rho_overex_on_nrho = register_diag_field ('ocean_model',           &
    'tform_rho_overex_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,   &
    'watermass transform due to overexchange as binned to neutral rho layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_overex_on_nrho > 0) compute_watermass_diag = .true.


  ! temperature contributions 
  id_neut_temp_overex = register_diag_field ('ocean_model', 'neut_temp_overex', &
    Grd%tracer_axes(1:3), Time%model_time,                                      &
    'temp related update of locally referenced potrho from overexchange scheme',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_overex > 0) compute_watermass_diag = .true.

  id_neut_temp_overex_on_nrho = register_diag_field ('ocean_model',                            &
    'neut_temp_overex_on_nrho',Dens%neutralrho_axes(1:3), Time%model_time,                     &
    'temp related update of locally ref potrho from overexchange binned to neutral rho layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_overex_on_nrho > 0) compute_watermass_diag = .true.

  id_wdian_temp_overex = register_diag_field ('ocean_model',     &
    'wdian_temp_overex', Grd%tracer_axes(1:3), Time%model_time,  &
    'temp related dianeutral mass transport due to overexchange',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_overex > 0) compute_watermass_diag = .true.

  id_wdian_temp_overex_on_nrho = register_diag_field ('ocean_model',                          &
    'wdian_temp_overex_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
    'temp related dianeutral mass transport due to overexchange binned to neutral rho layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_overex_on_nrho > 0) compute_watermass_diag = .true.

  id_tform_temp_overex_on_nrho = register_diag_field ('ocean_model',                       &
    'tform_temp_overex_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,               &
    'temp related watermass transform due to overexchange as binned to neutral rho layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_temp_overex_on_nrho > 0) compute_watermass_diag = .true.


  ! salinity contributions 
  id_neut_salt_overex = register_diag_field ('ocean_model', 'neut_salt_overex', &
    Grd%tracer_axes(1:3), Time%model_time,                                      &
    'salt related update of locally referenced potrho from overexchange scheme',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_overex > 0) compute_watermass_diag = .true.

  id_neut_salt_overex_on_nrho = register_diag_field ('ocean_model',                            &
    'neut_salt_overex_on_nrho',Dens%neutralrho_axes(1:3), Time%model_time,                     &
    'salt related update of locally ref potrho from overexchange binned to neutral rho layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_overex_on_nrho > 0) compute_watermass_diag = .true.

  id_wdian_salt_overex = register_diag_field ('ocean_model',     &
    'wdian_salt_overex', Grd%tracer_axes(1:3), Time%model_time,  &
    'salt related dianeutral mass transport due to overexchange',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_overex > 0) compute_watermass_diag = .true.

  id_wdian_salt_overex_on_nrho = register_diag_field ('ocean_model',                          &
    'wdian_salt_overex_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
    'salt related dianeutral mass transport due to overexchange binned to neutral rho layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_overex_on_nrho > 0) compute_watermass_diag = .true.

  id_tform_salt_overex_on_nrho = register_diag_field ('ocean_model',                    &
    'tform_salt_overex_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,            &
    'salt related watermass transform due to overexchange binned to neutral rho layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_salt_overex_on_nrho > 0) compute_watermass_diag = .true.


  if(compute_watermass_diag) then 
      write(stdoutunit,'(/a/)') &
      '==>Note: running ocean_overexchange_mod w/ compute_watermass_diag=.true.'  
  endif

end subroutine watermass_diag_init
! </SUBROUTINE> NAME="watermass_diag_init"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag">
!
! <DESCRIPTION>
! Diagnose effects from overexchange on the watermass transformation.
! </DESCRIPTION>
!
subroutine watermass_diag(Time, Dens, wrk1_v)

  type(ocean_time_type),           intent(in)    :: Time
  type(ocean_density_type),        intent(in)    :: Dens
  real, dimension(isd:,jsd:,:,:),  intent(inout) :: wrk1_v

  integer :: i,j,k

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_overexchange (watermass_diag): module needs initialization ')
  endif 

  if(.not. compute_watermass_diag) return

  ! full rho contributions 
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

  call diagnose_3d(Time, Grd, id_neut_rho_overex, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_overex, wrk3(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_overex_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_overex_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_overex_on_nrho, wrk4)

  ! temp contributions 
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

  call diagnose_3d(Time, Grd, id_neut_temp_overex, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_overex, wrk3(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_overex_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_overex_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_overex_on_nrho, wrk4)

  ! salinity contributions 
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

  call diagnose_3d(Time, Grd, id_neut_salt_overex, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_overex, wrk3(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_overex_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_overex_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_overex_on_nrho, wrk4)

end subroutine watermass_diag
! </SUBROUTINE>  NAME="watermass_diag"


end module ocean_overexchange_mod

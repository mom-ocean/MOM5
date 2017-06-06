module ocean_vert_kpp_mom4p0_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> A. Rosati
!</CONTACT>
!
!<REVIEWER EMAIL="wily@ucar.edu"> Bill Large 
!</REVIEWER>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen Griffies 
!</REVIEWER>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> M.J. Harrison 
!</REVIEWER>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">  Hyun-Chul Lee 
!</REVIEWER>
!
!<OVERVIEW>
! Vertical viscosity and diffusivity according KPP using 
! code from MOM4p0, which is hard-wired for full cell 
! GEOPOTENTIAL vertical coordinate model. It remains part of 
! MOM for legacy purposes.  
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes vertical viscosity and diffusivity according to 
! the K-profile parameterization scheme of Large, McWilliams, and 
! Doney (1994). It computes both local and non-local mixing.
!
! This module contains code that is hard-wired for GEOPOTENTIAL coordinates, 
! and so is NOT generally recommended.  It remains part of MOM for 
! legacy purposes. 
!
! This version of KPP has been implemented only for the Bgrid.
!
! This module also adds mixing due to barotropic tide drag 
! (coastal_tide_mix) and baroclinic tides (int_tide_mix). 
! The barotropic (coastal_tides) and baroclinic (int_tides) mixing 
! schemes are directly analogous to those available in the 
! module ocean_vert_tidal.F90.  However, some of the averaging and 
! smoothing operations differ, and so detailed comparisons will show
! differences. We retain the code here, in KPP, for legacy purposes. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! W.G. Large and J.C. McWilliams and S.C. Doney
! Oceanic vertical mixing: A review and a model with
! a nonlocal boundary layer parameterization
! Reviews of Geophysics (1994) vol 32 pages 363-403
! </REFERENCE>
!
! <REFERENCE>
! Hyun-Chul Lee, A. Rosati, and M.J. Spelman
! Barotropic tidal mixing impact in a coupled climate model:
! ocean condition and meridional overturning circulation 
! in the northern Atlantic
! Ocean Modelling, vol 11, pages 464--477
! </REFERENCE>
!
! <NOTE>
! Original numerical algorithm by Bill Large at NCAR June 6, 1994
! </NOTE>
!
! <NOTE>
! Equation numbers in the code refer to the Large etal paper. 
! </NOTE>
!
! <NOTE>
! Surface fresh water contributes to surface buoyancy via conversion to 
! a locally implied salt flux. 
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_vert_kpp_mom4p0_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Logical switch to enable kpp diffusion.  Default is false. 
!  </DATA> 
!
!  <DATA NAME="shear_instability" TYPE="logical">
!  logical switch for shear instability mixing.
!  Default shear_instability=.true.
!  </DATA> 
!  <DATA NAME="double_diffusion" TYPE="logical">
!  Logical switch for double-diffusive mixing.
!  Default double_diffusion=.true.  
!  </DATA> 
!  <DATA NAME="diff_cbt_iw" UNITS="m^2/sec" TYPE="real">
!  Background vertical diffusivity.  Note that if using Bryan-Lewis as a 
!  background diffusivity, then should set diff_cbt_iw=0.0. 
!  </DATA> 
!  <DATA NAME="visc_cbu_iw" UNITS="m^2/sec" TYPE="real">
!  Background vertical viscosity
!  </DATA> 
!  <DATA NAME="visc_cbu_limit" UNITS="m^2/sec" TYPE="real">
!  Enhanced vertical viscosity due to shear instability 
!  </DATA> 
!  <DATA NAME="diff_cbt_limit" UNITS="m^2/sec" TYPE="real">
!  Enhanced vertical diffusivity due to shear instability 
!  </DATA> 
!  <DATA NAME="visc_con_limit" UNITS="m^2/sec" TYPE="real">
!  Enhanced vertical viscosity in regions of convection
!  </DATA> 
!  <DATA NAME="diff_con_limit" UNITS="m^2/sec" TYPE="real">
!  Enhanced vertical diffusivity in regions of convection
!  </DATA> 
!  <DATA NAME="concv" UNITS="dimensionless" TYPE="real">
!  constant for pure convection (eqn. 23 of Large etal)
!  </DATA> 
!  <DATA NAME="Ricr" UNITS="dimensionless" TYPE="real">
!  Critical bulk Richardson number.  Default from NCAR is 
!  0.3, though this number has a large uncertainty and some
!  find that 1.0 can be of use. 
!  </DATA> 
!  <DATA NAME="non_local_kpp" TYPE="logical">
!  logical switch for enabling the non-local mixing aspect of kpp. 
!  Default is .true. as this is what the original KPP scheme suggests. 
!  </DATA> 
!  <DATA NAME="smooth_blmc" TYPE="logical">
!  Smooth boundary layer diffusitivies to remove grid scale noise.
!  Such noise is apparent in the diagnosed mixed layer depth as well
!  as the SST, especially when running coupled models where forcing 
!  has high temporal frequency. 
!  </DATA> 
!
!  <DATA NAME="coastal_tidal_mix" TYPE="logical">
!  For adding an extra vertical shear associated with tidal mixing.
!  This method has found to be of use for mixing near shelves.  
!  </DATA> 
!  <DATA NAME="p_tide" TYPE="real">
!  The p constant in the Munk-Anderson scheme
!  Default p_tide=-0.25
! </DATA> 
!  <DATA NAME="sigma_tide" TYPE="real">
!  The sigma constant in the Munk-Anderson scheme
!  Default sigma_tide=3.0
! </DATA> 
!
!  <DATA NAME="int_tidal_mix" TYPE="logical">
!  For adding an internal tidal mixing over rough topography.
!  This method has found to be of use for mixing in the rough topography in open ocean.  
!  Default int_tidal_mix=.false.
!  </DATA> 
!  <DATA NAME="int_tide_zeta1" TYPE="real" UNITS="metre">
!  Shallow depth for computation of internal tide.
!  Default int_tide_zeta1=300.0
! </DATA> 
!  <DATA NAME="int_tide_zeta2" TYPE="real" UNITS="metre">
!  Deeper depth for computation of internal tide.
!  Default int_tide_zeta2=1800.0
! </DATA> 
!  <DATA NAME="int_tide_min_depth" TYPE="real" UNITS="metre">
!  Minimum depth for internal tide mixing to be computed. 
!  Default int_tide_min_depth=100.0
! </DATA> 
!  <DATA NAME="int_tide_q" TYPE="real" UNITS="dimensionless">
!  Fraction of internal tide energy locally dissipated. 
!  Default int_tide_q=.33333
! </DATA> 
!  <DATA NAME="int_tide_gamma" TYPE="real" UNITS="dimensionless">
!  Dimensionless efficiency for converting energy dissipation to diffusivity.
!  Default int_tide_gamma=0.2
! </DATA> 
!
!  <DATA NAME="wsfc_combine_runoff_calve" TYPE="logical">
!  For computing wsfc as in the mom4p0d code, where we combine
!  the runoff+calving into a single field called river.  
!  The alternative keeps the fields separate, as would be appropriate
!  for a land model that separately tracks the tracer content in the 
!  calving and runoff. 
!  Default wsfc_combine_runoff_calve=.true., as this will recover
!  the previous behaviour, to the bit. 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: epsln
use diag_manager_mod, only: register_diag_field, register_static_field
use fms_mod,          only: FATAL, NOTE, stdout, stdlog
use fms_mod,          only: write_version_number, open_namelist_file, check_nml_error, close_file
use fms_mod,          only: read_data
use mpp_domains_mod,  only: mpp_update_domains, NUPDATE, EUPDATE
use mpp_mod,          only: input_nml_file, mpp_error

use ocean_density_mod,     only: density, density_delta_z, density_delta_sfc
use ocean_domains_mod,     only: get_local_indices
use ocean_parameters_mod,  only: GEOPOTENTIAL, cp_ocean, rho0, rho0r, grav, missing_value
use ocean_types_mod,       only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,       only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,       only: ocean_velocity_type, ocean_density_type
use ocean_types_mod,       only: ocean_time_type, ocean_time_steps_type, ocean_thickness_type
use ocean_workspace_mod,   only: wrk1, wrk2, wrk3, wrk4, wrk5
use ocean_util_mod,        only: diagnose_2d, diagnose_3d, diagnose_sum
use ocean_tracer_util_mod, only: diagnose_3d_rho


implicit none

private

public ocean_vert_kpp_mom4p0_init
public vert_mix_kpp_mom4p0

private bldepth
private wscale
private ri_iwmix
private ddmix
private blmix_kpp
private enhance
private ri_for_kpp 
private watermass_diag_init
private watermass_diag 

#include <ocean_memory.h>

#ifdef MOM_STATIC_ARRAYS

real, dimension(isd:ied,jsd:jed)                :: bfsfc      ! surface buoyancy forcing    (m^2/s^3)
real, dimension(isd:ied,jsd:jed)                :: ws         ! scalar velocity scale (m/s)
real, dimension(isd:ied,jsd:jed)                :: wm         ! momentum velocity scale (m/s)
real, dimension(isd:ied,jsd:jed)                :: Ustk2      ! magnitude of surface stokes drift velocity ^2 (m^2 / s^2)
real, dimension(isd:ied,jsd:jed)                :: caseA      ! = 1 in case A; =0 in case B
real, dimension(isd:ied,jsd:jed)                :: stable     ! = 1 in stable forcing; =0 in unstable
real, dimension(isd:ied,jsd:jed,3)              :: dkm1       ! boundary layer diff_cbt at kbl-1 level
real, dimension(isd:ied,jsd:jed,nk,3)           :: blmc       ! boundary layer mixing coefficients
real, dimension(isd:ied,jsd:jed)                :: sigma      ! normalized depth (d / hbl)
real, dimension(isd:ied,jsd:jed)                :: rhosfc     ! potential density of sfc layer(kg/m^3)
real, dimension(isd:ied,jsd:jed,nk)             :: talpha     ! -d(rho)/ d(pot.temperature)  (kg/m^3/C)
real, dimension(isd:ied,jsd:jed,nk)             :: sbeta      ! d(rho)/ d(salinity)       (kg/m^3/PSU)
real, dimension(isd:ied,jsd:jed,nk)             :: alphaDT    ! alpha * DT  across interfaces (kg/m^3)
real, dimension(isd:ied,jsd:jed,nk)             :: betaDS     ! beta  * DS  across interfaces (kg/m^3)
real, dimension(isd:ied,jsd:jed)                :: ustar      ! surface friction velocity       (m/s)
real, dimension(isd:ied,jsd:jed)                :: Bo         ! surface turb buoy. forcing  (m^2/s^3)
real, dimension(isd:ied,jsd:jed)                :: Bosol      ! radiative buoy forcing      (m^2/s^3)
real, dimension(isd:ied,jsd:jed,nk)             :: dbloc      ! local delta buoy at interfaces(m/s^2)
real, dimension(isd:ied,jsd:jed,nk)             :: dbsfc      ! delta buoy w/ respect to sfc  (m/s^2)
real, dimension(isd:ied,jsd:jed,nk)             :: dVsq       ! (velocity shear re sfc)^2   (m/s)^2
real, dimension(isd:ied,jsd:jed,2)              :: Rib        ! Bulk Richardson number
real, dimension(isd:ied,jsd:jed,3)              :: gat1
real, dimension(isd:ied,jsd:jed,3)              :: dat1
real, dimension(isd:ied,jsd:jed)                :: tidal_vel_amp   ! tidal amplitude velocity  (m/s) 
real, dimension(isd:ied,jsd:jed,nk)             :: tidal_fric_turb ! tidal friction            (1/s^2) 

real, dimension(isd:ied,jsd:jed)                :: rough_amp          ! roughness amplitude of topography  (m)
real, dimension(isd:ied,jsd:jed)                :: eflux_int_tide_bf  ! energy flux /buoyancy frequency
real, dimension(isd:ied,jsd:jed,nk)             :: f_int_tide         ! vertical distribution function   

type wsfc_type
  real, dimension(isd:ied,jsd:jed)              :: wsfc        ! rho0r*(stf - pme*(t(i,j,k=1)-tpme) - river*(t(i,j,k=1)-triver))
end type wsfc_type
real, dimension(isd:ied,jsd:jed)                :: sw_frac_hbl ! fractional shortwave penetration at base of mixed layer
integer, dimension(isd:ied,jsd:jed)             :: kbl         ! index of first grid level below hbl

real, private, dimension(isd:ied,jsd:jed,nk)    :: riu         ! Richardson number at base of U-cells
real, private, dimension(isd:ied,jsd:jed,nk)    :: rit         ! Richardson number at base of T-cells
real, private, dimension(isd:ied,jsd:jed,nk)    :: riu_tide    ! Richardson number at base of U-cells due to tidal friction
real, private, dimension(isd:ied,jsd:jed,nk)    :: rit_tide    ! Richardson number at base of T-cells due to tidal friction
real, private, dimension(isd:ied,jsd:jed,nk)    :: ghats       ! nonlocal transport (s/m^2)
real, private, dimension(isd:ied,jsd:jed)       :: hblt        ! boundary layer depth with tmask 
real, private, dimension(isd:ied,jsd:jed)       :: hbl         ! boundary layer depth
real, private, dimension(isd:ied,jsd:jed,nk)    :: bf_int_tide ! Buoyancy frequency in t-point


#else

real, dimension(:,:), allocatable      :: bfsfc    ! surface buoyancy forcing    (m^2/s^3)
real, dimension(:,:), allocatable      :: ws       ! scalar velocity scale (m/s)
real, dimension(:,:), allocatable      :: wm       ! momentum velocity scale (m/s)
real, dimension(:,:), allocatable      :: Ustk2    ! Magnitude of surface stokes drift velocity ^2 (m^2/s^2)
real, dimension(:,:), allocatable      :: caseA    ! = 1 in case A; =0 in case B
real, dimension(:,:), allocatable      :: stable   ! = 1 in stable forcing; =0 in unstable
real, dimension(:,:,:), allocatable    :: dkm1     ! boundary layer diff_cbt at kbl-1 level
real, dimension(:,:,:,:), allocatable  :: blmc     ! boundary layer mixing coefficients
real, dimension(:,:), allocatable      :: sigma    ! normalized depth (d / hbl)
real, dimension(:,:), allocatable      :: rhosfc   ! potential density of sfc layer(kg/m^3)
real, dimension(:,:,:), allocatable    :: talpha   ! -d(rho)/ d(pot.temperature)  (kg/m^3/C)
real, dimension(:,:,:), allocatable    :: sbeta    ! d(rho)/ d(salinity)       (kg/m^3/PSU)
real, dimension(:,:,:), allocatable    :: alphaDT  ! alpha * DT  across interfaces (kg/m^3)
real, dimension(:,:,:), allocatable    :: betaDS   ! beta  * DS  across interfaces (kg/m^3)
real, dimension(:,:), allocatable      :: ustar    ! surface friction velocity       (m/s)
real, dimension(:,:), allocatable      :: Bo       ! surface turb buoy. forcing  (m^2/s^3)
real, dimension(:,:), allocatable      :: Bosol    ! radiative buoy forcing      (m^2/s^3)
real, dimension(:,:,:), allocatable    :: dbloc    ! local delta buoy at interfaces(m/s^2)
real, dimension(:,:,:), allocatable    :: dbsfc    ! delta buoy w/ respect to sfc  (m/s^2)
real, dimension(:,:,:), allocatable    :: dVsq     ! (velocity shear re sfc)^2   (m/s)^2
real, dimension(:,:,:), allocatable    :: Rib      ! Bulk Richardson number
real, dimension(:,:,:), allocatable    :: gat1
real, dimension(:,:,:), allocatable    :: dat1

real, dimension(:,:), allocatable      :: tidal_vel_amp   ! tidal amplitude velocity (m/s)
real, dimension(:,:,:), allocatable    :: tidal_fric_turb ! tidal friction           (1/s^2) 

real, dimension(:,:), allocatable      :: rough_amp          ! roughness amplitude of topography  (m)
real, dimension(:,:), allocatable      :: eflux_int_tide_bf  ! energy flux /buoyancy frequency
real, dimension(:,:,:), allocatable    :: f_int_tide         ! vertical distribution function   

type wsfc_type
  real, dimension(:,:), pointer        :: wsfc => NULL()  ! rho0r*(stf - pme*(t(i,j,k=1)-tpme) - river*(t(i,j,k=1)-triver))
end type wsfc_type
real, dimension(:,:), allocatable      :: sw_frac_hbl     ! fractional shortwave penetration at base of mixed layer
integer, dimension(:,:), allocatable   :: kbl             ! index of first grid level below hbl

real, private, dimension(:,:,:), allocatable :: riu         ! Richardson number at base of U-cells
real, private, dimension(:,:,:), allocatable :: rit         ! Richardson number at base of T-cells
real, private, dimension(:,:,:), allocatable :: rit_tide    ! Richardson number at base of T-cells due to tidal friction
real, private, dimension(:,:,:), allocatable :: riu_tide    ! Richardson number at base of U-cells due to tidal friction
real, private, dimension(:,:,:), allocatable :: ghats       ! nonlocal transport (s/m^2)
real, private, dimension(:,:),   allocatable :: hblt        ! boundary layer depth with tmask 
real, private, dimension(:,:),   allocatable :: hbl         ! boundary layer depth
real, private, dimension(:,:,:), allocatable :: bf_int_tide ! Buoyancy frequency in t-point

#endif

type(wsfc_type), dimension(:), allocatable     :: wsfc
real, parameter :: epsilon = 0.1
real, parameter :: vonk    = 0.4
real, parameter :: conc1   = 5.0
real, parameter :: zmin    = -4.e-7  ! m3/s3 limit for lookup table of wm and ws
real, parameter :: zmax    = 0.0     ! m3/s3 limit for lookup table of wm and ws
real, parameter :: umin    = 0.0     ! m/s limit for lookup table of wm and ws
real, parameter :: umax    = 0.04    ! m/s limit for lookup table of wm and ws

real :: Ricr               = 0.3     ! critical bulk Richardson Number
real :: visc_cbu_limit     = 50.0e-4 ! max visc due to shear instability
real :: diff_cbt_limit     = 50.0e-4 ! max diff due to shear instability
real :: visc_con_limit     = 0.1     ! m^2/s. visc due to convective instability
real :: diff_con_limit     = 0.1     ! m^2/s. diff due to convective instability
real :: visc_cbu_iw        = 1.0e-4  ! m^2/s. visc background due to internal waves
real :: diff_cbt_iw        = 0.1e-4  ! m^2/s. diffusivity background due to internal waves
real :: Vtc                          ! non-dimensional coefficient for velocity 
                                     ! scale of turbulant velocity shear        
                                     ! (=function of concv,concs,epsilon,vonk,Ricr)

real    :: cg            ! non-dimensional coefficient for counter-gradient term
real    :: deltaz        ! delta zehat in table
real    :: deltau        ! delta ustar in table
real    :: concv  = 1.8  ! constant for pure convection (eqn. 23)
real    :: Lgam = 1.04     ! adjustment to non-gradient flux (McWilliam & Sullivan 2000)
real    :: Cw_0 = 0.15     ! eq. (13) in Smyth et al (2002)
real    :: l_smyth = 2.0   ! eq. (13) in Smyth et al (2002)
real    :: LTmax = 5.0     ! maximum Langmuir turbulence enhancement factor (langmuirfactor) allowed
real    :: Wstfac = 0.6    ! stability adjustment coefficient, eq. (13) in Smyth et al (2002)

! for vertical coordinate 
integer :: vert_coordinate

real :: rho_cp
real :: inv_rho_cp

! for global area normalization
real :: cellarea_r

logical :: use_this_module           = .false. ! for turning on the kpp scheme 
logical :: shear_instability         = .true.  ! for shear instability mixing
logical :: double_diffusion          = .true.  ! for double-diffusive mixing
logical :: wsfc_combine_runoff_calve = .true.  ! for combining runoff+calving to compute wfsc
real    :: shear_instability_flag    = 1.0     ! set to 1.0 if shear_instability=.true.

! internally set for computing watermass diagnstics
logical :: compute_watermass_diag = .false. 

! for diagnostics 
integer, dimension(:), allocatable :: id_nonlocal(:)
integer, dimension(:), allocatable :: id_nonlocal_on_nrho(:)
logical :: used
integer :: id_diff_cbt_kpp_t   =-1
integer :: id_diff_cbt_kpp_s   =-1
integer :: id_diff_cbt_int     =-1
integer :: id_visc_cbt_int     =-1
integer :: id_diff_cbt_coast   =-1
integer :: id_visc_cbt_coast   =-1
integer :: id_hblt             =-1
integer :: id_tidal_vel_amp    =-1
integer :: id_tidal_fric_turb  =-1
integer :: id_rough_amp        =-1
integer :: id_eflux_int_tide_bf=-1
integer :: id_f_int_tide       =-1
integer :: id_bf_int_tide      =-1

integer  :: id_neut_rho_kpp_nloc          =-1
integer  :: id_wdian_rho_kpp_nloc         =-1
integer  :: id_tform_rho_kpp_nloc         =-1
integer  :: id_neut_rho_kpp_nloc_on_nrho  =-1
integer  :: id_wdian_rho_kpp_nloc_on_nrho =-1
integer  :: id_tform_rho_kpp_nloc_on_nrho =-1

integer :: id_eta_tend_kpp_nloc     =-1
integer :: id_eta_tend_kpp_nloc_glob=-1

integer  :: id_neut_temp_kpp_nloc          =-1
integer  :: id_wdian_temp_kpp_nloc         =-1
integer  :: id_tform_temp_kpp_nloc         =-1
integer  :: id_neut_temp_kpp_nloc_on_nrho  =-1
integer  :: id_wdian_temp_kpp_nloc_on_nrho =-1
integer  :: id_tform_temp_kpp_nloc_on_nrho =-1

integer  :: id_neut_salt_kpp_nloc          =-1
integer  :: id_wdian_salt_kpp_nloc         =-1
integer  :: id_tform_salt_kpp_nloc         =-1
integer  :: id_neut_salt_kpp_nloc_on_nrho  =-1
integer  :: id_wdian_salt_kpp_nloc_on_nrho =-1
integer  :: id_tform_salt_kpp_nloc_on_nrho =-1


logical :: non_local_kpp      = .true.  ! enable/disable non-local term in KPP
logical :: smooth_blmc        = .false. ! smooth boundary layer diffusitivies to remove grid scale noise
logical :: do_langmuir        = .false. ! whether or not calculate langmuir turbulence enhancement factor

logical :: coastal_tidal_mix  = .false. ! add tidal speed (m/s) to Ri to increase coastal mixing
real    :: sigma_tide         = 3.0     ! the sigma constant in the Munk-Anderson scheme
real    :: p_tide             = -0.25   ! the p constant in the Munk-Anderson scheme
real    :: coastal_tidal_flag = 0.0     ! set to unity if coastal_tidal_mix=.true.

logical :: int_tidal_mix      = .false. ! add internal tidal dissipation over rough topograpgy in the open ocean
real    :: int_tide_zeta1     = 300.0   ! shallow depth (m) for computation of internal tide mixing
real    :: int_tide_zeta2     = 1800.0  ! deeper depth (m) for computation of internal tide mixing
real    :: int_tide_min_depth = 100.0   ! minimum depth (m) a region must posses to employ int tide mixing
real    :: int_tide_q         = 0.33333 ! Fraction of internal tide energy locally dissipated.
real    :: int_tide_gamma     = 0.2     ! Dimensionless efficiency for converting energy dissipation to diffusivity.

integer, parameter :: nni = 890         ! number of values for zehat in the look up table
integer, parameter :: nnj = 480         ! number of values for ustar in the look up table
real, dimension(0:nni+1,0:nnj+1) :: wmt ! lookup table for wm, the turbulent velocity scale for momentum
real, dimension(0:nni+1,0:nnj+1) :: wst ! lookup table for ws, the turbulent velocity scale scalars

real :: tracer_timestep = 0

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type), pointer   :: Grd => NULL()

integer :: num_prog_tracers=0, index_temp, index_salt
integer :: num_diag_tracers=0, index_frazil

character(len=256) :: version=&
     '$Id: ocean_vert_kpp_mom4p0.F90,v 20.0 2013/12/14 00:16:42 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

logical :: module_is_initialized = .FALSE.

namelist /ocean_vert_kpp_mom4p0_nml/ use_this_module, shear_instability, double_diffusion,  &
                              diff_cbt_iw, visc_cbu_iw,                             &
                              visc_cbu_limit, diff_cbt_limit,                       &
                              visc_con_limit, diff_con_limit,                       &
                              concv, Ricr, non_local_kpp, smooth_blmc,              &
                              Lgam, Cw_0,l_smyth, LTmax, Wstfac,                    &
                              coastal_tidal_mix, p_tide, sigma_tide,                &
                              int_tidal_mix, int_tide_zeta1, int_tide_zeta2,        &
                              int_tide_min_depth, int_tide_q, int_tide_gamma,       &
                              wsfc_combine_runoff_calve, do_langmuir

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_vert_kpp_mom4p0_init">
!
! <DESCRIPTION>
! Initialization for the KPP vertical mixing scheme
!
!     input:
!       dzt    = thickness of vertical levels (m)                         <BR/>
!       km     = number of vertical levels                                <BR/> 
!       yt     = latitude of grid points (deg)                            <BR/> 
!       jmt    = number of latitudes                                      <BR/>  
!       dtxcel = time step accelerator as a function of level             <BR/>    
!       dtimet = forward time step for tracer diffusion (sec)             <BR/>   
!       dtimeu = forward time step for velcotiy friction (sec)            <BR/>   
!       error  = logical to signal problems                               <BR/>
!       cifdef = array of character strings for listing enabled "ifdefs"  <BR/>
!       ifdmax = size of "cifdef"                                         <BR/> 
!       nifdef = current number of enabled "ifdefs"                       <BR/>   
!       vmixset= logical to determine if a vertical mixing scheme was     <BR/> 
!                chosen
!
!     output:
!       shear_instability = logical switch for shear instability mixing   <BR/>   
!       double_diffusion = logical switch for double-diffusive mixing     <BR/>
!       visc_cbu_limit = visc max due to shear instability  (m**2/sec)    <BR/>
!       diff_cbt_limit = diffusivity ..                     (m**2/sec)    <BR/>
!       visc_cbu_iw  = visc background due to internal waves(m**2/sec)    <BR/>
!       diff_cbt_iw  = diffusivity ..                       (m**2/sec)    <BR/> 
!       visc_con_limit = visc due to convective instability (m**2/sec)    <BR/>
!       diff_con_limit = diffusivity ..                     (m**2/sec)    <BR/>
!       Vtc = non-dimensional constant used in calc. bulk Ri              <BR/>
!       cg  = constant used in calc.nonlocal transport term               <BR/>  
!       wmt = turbulent velocity scale for momentum                       <BR/>  
!       wst = turbulent velocity scale for scaler                         <BR/> 
!       error  = true if some inconsistancy was found                    
!
! </DESCRIPTION>
!
subroutine ocean_vert_kpp_mom4p0_init (Grid, Domain, Time, Time_steps, Dens, T_prog, T_diag, ver_coordinate)
  
  type(ocean_grid_type),        intent(in), target  :: Grid
  type(ocean_domain_type),      intent(in), target  :: Domain
  type(ocean_time_type),        intent(in)          :: Time
  type(ocean_time_steps_type),  intent(in)          :: Time_steps
  type(ocean_density_type),     intent(in)          :: Dens
  type(ocean_prog_tracer_type), intent(in)          :: T_prog(:)
  type(ocean_diag_tracer_type), intent(in)          :: T_diag(:)
  integer,                      intent(in)          :: ver_coordinate
  
  real, parameter :: cstar  = 10.0  ! proportionality coefficient for nonlocal transport
  real, parameter :: conam  = 1.257
  real, parameter :: concm  = 8.380 
  real, parameter :: conc2  = 16.0
  real, parameter :: zetam  = -0.2
  real, parameter :: conas  = -28.86
  real, parameter :: concs  = 98.96
  real, parameter :: conc3  = 16.0
  real, parameter :: zetas  = -1.0

  real :: zehat ! = zeta * ustar**3
  real :: zeta  ! = stability parameter d/L
  real :: usta

  real :: sqrt_cd, A_tidal
  real :: tidal_vel_amp_t, kapa_int_tide, rho_ref_int_tide
  real :: delz_int_tide, active_cells
  real :: zeta_int_tide1, zeta_int_tide2
  real :: upper_int_tide1, lower_int_tide1
  real :: upper_int_tide2, lower_int_tide2
  integer :: i, j, k, kmax, n, ioun, io_status, ierr

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_vert_kpp_mod (ocean_vert_kpp_mom4p0_init): module is already initialized')
  endif 

  module_is_initialized = .TRUE.
  vert_coordinate = ver_coordinate 
  rho_cp     = rho0*cp_ocean 
  inv_rho_cp = 1.0/rho_cp 

  call write_version_number(version, tagname)

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_vert_kpp_mom4p0_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_vert_kpp_mom4p0_nml')
#else
  ioun =  open_namelist_file ()
  read  (ioun, ocean_vert_kpp_mom4p0_nml,iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_vert_kpp_mom4p0_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_vert_kpp_mom4p0_nml)  
  write (stdlogunit, ocean_vert_kpp_mom4p0_nml)

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk
#endif 

  if(use_this_module) then 
    write(stdoutunit,'(/1x,a)') '==> NOTE: USING KPP_mom4p0 vertical mixing scheme.'
    write(stdoutunit,'(1x,a)')  '          This scheme is hard-wired for GEOPOTENTIAL coordinates.'
    write(stdoutunit,'(1x,a)')  '          It is not generally recommended for use, other than for legacy.'
    write(stdoutunit,'(1x,a)')  '==> NOTE: KPP is typically run with penetrative shortwave heating.'
    write(stdoutunit,'(1x,a/)') '==> NOTE: KPP is typically run with a seasonal and/or diurnal cycle.'
  else 
    write(stdoutunit,'(/1x,a)')'==> NOTE: NOT USING KPP_mom4p0 vertical mixing scheme.'
    return
  endif 

  if(vert_coordinate /= GEOPOTENTIAL) then 
    call mpp_error(NOTE,&
    '==>ocean_vert_kpp_mom4p0_mod: This module is hard-wired for full cell GEOPOTENTIAL coordinate. ')
  endif  

  write(stdoutunit,'(/a,f10.2)')'==>From ocean_vert_kpp_mom4p0_mod: time step for vert-frict of (secs)', &
                               Time_steps%dtime_u 
  write(stdoutunit,'(/a,f10.2)')'==>From ocean_vert_kpp_mom4p0_mod: time step for vert-diff  of (secs)', &
                               Time_steps%dtime_t 

  tracer_timestep = Time_steps%dtts

  index_temp=-1;index_salt=-1
  num_prog_tracers = size(T_prog(:))
  do n= 1, num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo
  allocate( wsfc(num_prog_tracers) )

  index_frazil=-1
  num_diag_tracers = size(T_diag(:))
  do n= 1, num_diag_tracers
     if (T_diag(n)%name == 'frazil') index_frazil = n
  enddo

  if (index_temp < 1 .or. index_salt < 1) then 
    call mpp_error(FATAL,'==>Error: ocean_vert_kpp_mom4p0_mod: temp and/or salt not present in tracer array')
  endif 

  Dom => Domain
  Grd => Grid

  cellarea_r = 1.0/(epsln + Grd%tcellsurf)

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk

  allocate (ghats(isd:ied,jsd:jed,nk))
  allocate (riu(isd:ied,jsd:jed,nk))
  allocate (riu_tide(isd:ied,jsd:jed,nk))
  allocate (rit(isd:ied,jsd:jed,nk))
  allocate (rit_tide(isd:ied,jsd:jed,nk))
  allocate (hblt(isd:ied,jsd:jed))
  allocate (hbl(isd:ied,jsd:jed))
  do n=1,num_prog_tracers 
    allocate( wsfc(n)%wsfc(isd:ied,jsd:jed) )
  enddo


  allocate (bfsfc(isd:ied,jsd:jed))         ! surface buoyancy forcing    (m^2/s^3)
  allocate (caseA(isd:ied,jsd:jed))         ! = 1 in case A; =0 in case B
  allocate (stable(isd:ied,jsd:jed))        ! = 1 in stable forcing; =0 in unstable
  allocate (dkm1(isd:ied,jsd:jed,3))        ! boundary layer diff_cbt at kbl-1 level
  allocate (blmc(isd:ied,jsd:jed,nk,3))     ! boundary layer mixing coefficients
  allocate (sigma(isd:ied,jsd:jed))         ! normalized depth (d / hbl)
  allocate (rhosfc(isd:ied,jsd:jed))        ! potential density of sfc layer(g/m^3)
  allocate (talpha(isd:ied,jsd:jed,nk))     ! -d(rho)/ d(pot.temperature)  (g/m^3/C)
  allocate (sbeta(isd:ied,jsd:jed,nk))      ! d(rho)/ d(salinity)       (g/m^3/PSU)
  allocate (alphaDT(isd:ied,jsd:jed,nk))    ! alpha * DT  across interfaces (g/m^3)
  allocate (betaDS(isd:ied,jsd:jed,nk))     ! beta  * DS  across interfaces (g/m^3)
  allocate (ustar(isd:ied,jsd:jed))         ! surface friction velocity       (m/s)
  allocate (Bo(isd:ied,jsd:jed))            ! surface turb buoy. forcing  (m^2/s^3)
  allocate (Bosol(isd:ied,jsd:jed))         ! radiative buoy forcing      (m^2/s^3)
  allocate (dbloc(isd:ied,jsd:jed,nk))      ! local delta buoy at interfaces(m/s^2)
  allocate (dbsfc(isd:ied,jsd:jed,nk))      ! delta buoy w/ respect to sfc  (m/s^2)
  allocate (dVsq(isd:ied,jsd:jed,nk))       ! (velocity shear re sfc)^2   (m/s)^2
  allocate (kbl(isd:ied,jsd:jed))           ! index of first grid level below hbl
  allocate (Rib(isd:ied,jsd:jed,2))         ! Bulk Richardson number
  allocate (wm(isd:ied,jsd:jed))            ! momentum turbulent velocity scales  (m/s)
  allocate (ws(isd:ied,jsd:jed))            ! scalar turbulent velocity scales  (m/s)
  allocate (Ustk2(isd:ied,jsd:jed))         ! Magnitude of surface stokes drift velocity ^2 (m^2/s^2)
  allocate (gat1(isd:ied,jsd:jed,3))
  allocate (dat1(isd:ied,jsd:jed,3))
  allocate(sw_frac_hbl(isd:ied,jsd:jed))
  allocate (tidal_vel_amp(isd:ied,jsd:jed))      ! tidal speed (m/s)
  allocate (tidal_fric_turb(isd:ied,jsd:jed,nk)) ! tidal friction (1/s^2)
  allocate (rough_amp(isd:ied,jsd:jed))          ! roughness amplitude of topography  (m)
  allocate (eflux_int_tide_bf(isd:ied,jsd:jed))  ! energy flux /buoyancy frequency
  allocate (f_int_tide(isd:ied,jsd:jed,nk))      ! vertical distribution function
  allocate (bf_int_tide(isd:ied,jsd:jed,nk))     ! Buoyancy frequency
#endif

  ghats(:,:,:)     = 0.0
  riu(:,:,:)       = 0.0
  riu_tide(:,:,:)  = 0.0
  rit(:,:,:)       = 0.0
  rit_tide(:,:,:)  = 0.0
  hblt(:,:)        = 0.0
  hbl(:,:)         = 0.0
  sw_frac_hbl(:,:) = 0.0
  Ustk2(:,:)       = 0.0

  do n = 1, num_prog_tracers  
    wsfc(n)%wsfc(:,:) = 0.0
  enddo  
  tidal_vel_amp(:,:)     = 0.0
  tidal_fric_turb(:,:,:) = 0.0
  rough_amp(:,:)         = 0.0
  eflux_int_tide_bf(:,:) = 0.0
  f_int_tide(:,:,:)      = 0.0
  bf_int_tide(:,:,:)     = 0.0


!-----------------------------------------------------------------------
!     initialize some constants for kmix subroutines, and initialize
!     for kmix subroutine "wscale" the 2D-lookup table for wm and ws
!     as functions of ustar and zetahat (=vonk*sigma*hbl*bfsfc).
!-----------------------------------------------------------------------
 

!-----------------------------------------------------------------------
! define some non-dimensional constants  (recall epsilon=0.1)
!-----------------------------------------------------------------------

!     Vtc used in eqn. 23
      Vtc     = concv * sqrt(0.2/concs/epsilon) / vonk**2 / Ricr

!     cg = cs in eqn. 20
      cg      = cstar * vonk * (concs * vonk * epsilon)**(1./3.)

!-----------------------------------------------------------------------
! construct the wm and ws lookup tables (eqn. 13 & B1)
!-----------------------------------------------------------------------

      deltaz = (zmax-zmin)/(nni+1) 
      deltau = (umax-umin)/(nnj+1)
      
      do i=0,nni+1
         zehat = deltaz*(i) + zmin
         do j=0,nnj+1
            usta = deltau*(j) + umin
            zeta = zehat/(usta**3+epsln)

            if(zehat >= 0.) then
               wmt(i,j) = vonk*usta/(1.+conc1*zeta)
               wst(i,j) = wmt(i,j)
            else
               if(zeta > zetam) then
                  wmt(i,j) = vonk* usta * (1.-conc2*zeta)**(1./4.)
               else
                  wmt(i,j) = vonk* (conam*usta**3-concm*zehat)**(1./3.)
               endif
               if(zeta > zetas) then
                  wst(i,j) = vonk* usta * (1.-conc3*zeta)**(1./2.)
               else
                  wst(i,j) = vonk* (conas*usta**3-concs*zehat)**(1./3.)
               endif
            endif   
         enddo   
      enddo

      ! set flag for shear instability mixing 
      if (shear_instability) then
          write (stdoutunit,'(a)') 'Computing vertical mixing from shear instability in KPP module.'
          shear_instability_flag=1.0
      else 
          write (stdoutunit,'(a)') 'NOT computing vertical mixing from shear instability in KPP module.'
          shear_instability_flag=0.0
      endif 

!-----------------------------------------------------------------------
! read tidal speed (m/s) from Global Inverse Solution
! TPX06.0 created by OSU (only M2 tides), interpolated 1/4 x 1/4  
! to ocean model grid by Gaussian method
!
! After read data, then compute static tidal frictional turbulence
!-----------------------------------------------------------------------
      if (coastal_tidal_mix) then
          write (stdoutunit,'(a)') 'Computing vertical mixing from barotropic tides within KPP_mom4p0 module.'
          coastal_tidal_flag=1.0
      endif 
      if (int_tidal_mix) then
          write (stdoutunit,'(a)') 'Computing vertical mixing from internal tides within KPP_mom4p0 module.'
      endif 

      if (coastal_tidal_mix .or. int_tidal_mix) then
          call read_data('INPUT/tideamp.nc','tideamp', tidal_vel_amp, Domain%domain2d, 1)
          write (stdoutunit,*) 'Completed read of tidal velocity amplitude'
          call mpp_update_domains (tidal_vel_amp(:,:), Dom%domain2d)
      endif

      ! calculate the tidal frictional turbulence (frictioanl shear term) based 
      ! on the tidal velocity amplitude as a static variable.
      if (coastal_tidal_mix) then
          sqrt_cd = sqrt(2.4e-3)
          do j=jsd,jed
             do i=isd,ied
                A_tidal = sqrt_cd * tidal_vel_amp(i,j)/vonk
                kmax = Grd%kmt(i,j)
                if (kmax > 0) then
                    do k=1,kmax
                       tidal_fric_turb(i,j,k) = 0.5*(A_tidal/(Grd%zw(kmax)-Grd%zt(k)))**2
                    enddo
                endif
             enddo
          enddo

         ! update halo regions to compute tidal frictional turbulence to place onto U-cell 
         call mpp_update_domains (tidal_fric_turb(:,:,:), Dom%domain2d)

      endif

!-----------------------------------------------------------------------
! Parameterization of internal tidal dissipation over rough topography
! Read input data of roughness amplitude of topography, rough_amp (meter unit)
! logical switch is int_tidal_mix, default is .false.
! eflux_int_tide_bf: internal tidal mixing energy flux devided by buoyancy frequency
! kapa_int_tide : wave number of rough topography = 2*3.14/(10 km), (Jayne and St. Laurent, 2001)
! rho_ref_int_tide=1035 , refernce sea water density, (Gill,1982)
! f_int_tide , vertical distribution function
! zeta_int_tide=500 m, vertical e-folding scale
!-----------------------------------------------------------------------
      if (int_tidal_mix) then

          call read_data('INPUT/rough_amp.nc','roughness_length', rough_amp, Domain%domain2d, 1)
          write (stdoutunit,*) 'Completed read of roughness amplitude of topography'

      kapa_int_tide=2.0*3.1415926/(10.*1.e3)
      rho_ref_int_tide=1035.0                 

          do j=jsc,jec
             do i=isc,iec

        active_cells = Grd%umask(i,j,1) + Grd%umask(i-1,j,1) + Grd%umask(i,j-1,1) + Grd%umask(i-1,j-1,1) + epsln
                tidal_vel_amp_t=(tidal_vel_amp(i,j)+tidal_vel_amp(i,j-1)        &
                                +tidal_vel_amp(i-1,j)+tidal_vel_amp(i-1,j-1))/active_cells
                eflux_int_tide_bf(i,j)=0.5*rho_ref_int_tide*kapa_int_tide*rough_amp(i,j)*rough_amp(i,j) &
                                  *tidal_vel_amp_t*tidal_vel_amp_t
             enddo
          enddo

      ! update halo regions to compute tidal frictional turbulence to place onto U-cell 
      call mpp_update_domains (eflux_int_tide_bf(:,:), Dom%domain2d)


          ! define the verical distribution function for energy flux, f_int_tide
          ! this uses static rigid lid geopotential vertical coordinate depths.  
          ! this is strictly inaccurate.  see the ocean_vert_tidal.F90 for better way. 
          do j=jsd,jed
             do i=isd,ied
                kmax = Grd%kmt(i,j)
                if (kmax > 0) then
                    do k=1,kmax
                        delz_int_tide=Grd%zw(kmax)-Grd%zt(k)
                        upper_int_tide1=exp(-1.0*(delz_int_tide/zeta_int_tide1))
                        lower_int_tide1=zeta_int_tide1 &
                        *(1.0-exp(-1.0*(Grd%zw(kmax)/zeta_int_tide1)))
                        upper_int_tide2=exp(-1.0*(delz_int_tide/zeta_int_tide2))
                        lower_int_tide2=zeta_int_tide2 &
                        *(1.0-exp(-1.0*(Grd%zw(kmax)/zeta_int_tide2)))
                        f_int_tide(i,j,k) = 0.5*(upper_int_tide1/lower_int_tide1+upper_int_tide2/lower_int_tide2)
                    enddo
                endif
             enddo
          enddo
      
      endif

!-----------------------------------------------------------------------
! error checks and diagnostics setup 
!-----------------------------------------------------------------------

  if (Time_steps%aidif /= 1.0) then
    call mpp_error(FATAL,&
    '==>Error in ocean_vert_kpp_mom4p0: kpp must use aidif=1 for implicit vertical mixing')
  endif

  ! register diagnostics 

  allocate(id_nonlocal(num_prog_tracers))
  allocate(id_nonlocal_on_nrho(num_prog_tracers))
  id_nonlocal=-1
  id_nonlocal_on_nrho=-1
  do n = 1, num_prog_tracers
     if(n==index_temp) then
        id_nonlocal(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_nonlocal_KPP', &
                     Grd%tracer_axes(1:3), Time%model_time,                                         &
                     'cp*rho*dzt*nonlocal tendency from KPP', trim(T_prog(n)%flux_units),           &
                     missing_value=missing_value, range=(/-1.e10,1.e10/))
        id_nonlocal_on_nrho(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_nonlocal_KPP_on_nrho', &
                     Dens%neutralrho_axes(1:3), Time%model_time,                                         &
                     'cp*rho*dzt*nonlocal tendency from KPP binned to neutral density', trim(T_prog(n)%flux_units), &
                     missing_value=missing_value, range=(/-1.e20,1.e20/))
     else
        id_nonlocal(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_nonlocal_KPP', &
                     Grd%tracer_axes(1:3), Time%model_time,                                         &
                     'rho*dzt*nonlocal tendency from KPP', trim(T_prog(n)%flux_units),              &
                     missing_value=missing_value, range=(/-1.e10,1.e10/))
        id_nonlocal_on_nrho(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_nonlocal_KPP_on_nrho', &
                     Dens%neutralrho_axes(1:3), Time%model_time,                                         &
                     'rho*dzt*nonlocal tendency from KPP binned to neutral density', trim(T_prog(n)%flux_units),              &
                     missing_value=missing_value, range=(/-1.e20,1.e20/))
     endif
  enddo

  id_diff_cbt_kpp_t = register_diag_field('ocean_model','diff_cbt_kpp_t',Grd%tracer_axes(1:3),&
       Time%model_time, 'vert diffusivity from kpp for temp', 'm^2/sec',                      &
       missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_diff_cbt_kpp_s = register_diag_field('ocean_model','diff_cbt_kpp_s',Grd%tracer_axes(1:3),&
       Time%model_time, 'vert diffusivity from kpp for salt', 'm^2/sec',                      &
       missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_diff_cbt_int = register_diag_field('ocean_model','diff_cbt_int',Grd%tracer_axes(1:3), &
       Time%model_time, 'vert diffusivity from kpp internal tide mixing', 'm^2/sec',       &
       missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_visc_cbt_int = register_diag_field('ocean_model','visc_cbt_int',Grd%tracer_axes(1:3),&
       Time%model_time, 'vert viscosity from kpp internal tide mixing', 'm^2/sec',        &
       missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_diff_cbt_coast = register_diag_field('ocean_model','diff_cbt_coast',Grd%tracer_axes(1:3),&
       Time%model_time, 'vert diffusivity from kpp coastal tide mixing', 'm^2/sec',           &
       missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_visc_cbt_coast = register_diag_field('ocean_model','visc_cbt_coast',Grd%tracer_axes(1:3),&
       Time%model_time, 'vert viscosity from kpp coastal tide mixing', 'm^2/sec',             &
       missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_hblt           = register_diag_field('ocean_model','hblt',Grd%tracer_axes(1:2), &
       Time%model_time, 'T-cell boundary layer depth from KPP', 'm',                 &
       missing_value = missing_value, range=(/-1.e5,1.e6/),                          &
       standard_name='ocean_mixed_layer_thickness_defined_by_mixing_scheme')

  id_bf_int_tide    = register_diag_field('ocean_model','bf_int_tide',Grd%tracer_axes(1:3),&
       Time%model_time, 'buoyancy frequency', '1/sec',                                     &
       missing_value = missing_value, range=(/0.0,1.e-1/))

  ! static fields 
  id_tidal_vel_amp  = register_static_field('ocean_model','tidal_vel_amp',Grd%tracer_axes(1:2), &
       'M2 tidal speed from OSU tide model', 'm/s', missing_value=-10.0, range=(/-10.0,1.e6/))

  id_tidal_fric_turb  = register_static_field('ocean_model','tidal_fric_turb',Grd%tracer_axes(1:3), &
       'tidal friction', '1/s^2', missing_value=-10.0, range=(/-10.0,1.e6/))

   
  id_rough_amp  = register_static_field('ocean_model','rough_amp',Grd%tracer_axes(1:2), &
       ' roughness amplitude of topography', 'm', missing_value=-10.0, range=(/-10.0,1.e6/))

  id_eflux_int_tide_bf  = register_static_field('ocean_model','eflux_int_tide_bf',Grd%tracer_axes(1:2), &
       'energy flux devided by buoyancy frequency', '(W/m**2)*s', missing_value=-10.0, range=(/-10.0,1.e9/))

  id_f_int_tide  = register_static_field('ocean_model','f_int_tide',Grd%tracer_axes(1:3), &
       'verical distriibution function for energy flux', ' ', missing_value=-10.0, range=(/-10.0,1.e3/))

  call diagnose_2d(Time, Grd, id_tidal_vel_amp, tidal_vel_amp(:,:))
  call diagnose_3d(Time, Grd, id_tidal_fric_turb, tidal_fric_turb(:,:,:))
  call diagnose_2d(Time, Grd, id_rough_amp, rough_amp(:,:))
  call diagnose_2d(Time, Grd, id_eflux_int_tide_bf, eflux_int_tide_bf(:,:))
  call diagnose_3d(Time, Grd, id_f_int_tide, f_int_tide(:,:,:))

 call watermass_diag_init(Time,Dens)


end subroutine ocean_vert_kpp_mom4p0_init
! </SUBROUTINE> NAME="ocean_vert_kpp_mom4p0_init">



!#######################################################################
! <SUBROUTINE NAME="vert_mix_kpp_mom4p0">
!
! <DESCRIPTION>
! This subroutine computes the vertical diffusivity and viscosity according
! to the KPP scheme of Large etal.  In brief, the scheme does the 
! following:
!
! --Compute interior mixing everywhere:                               
!   interior mixing gets computed at all cell interfaces due to constant
!   internal wave background activity ("visc_cbu_iw" and "diff_cbt_iw").
!   Mixing is enhanced in places of static instability (local Ri < 0).
!   Additionally, mixing can be enhanced by contribution from shear 
!   instability which is a function of the local Ri.
!
! --Double diffusion:
!   Interior mixing can be enhanced by double diffusion due to salt
!   fingering and diffusive convection ("double_diffusion=.true.").
!
! --Boundary layer:
!
!   (A) Boundary layer depth:
!       at every gridpoint the depth of the oceanic boundary layer 
!       ("hbl") gets computed by evaluating bulk richardson numbers.
!
!   (B) Boundary layer mixing:
!       within the boundary layer, above hbl, vertical mixing is 
!       determined by turbulent surface fluxes, and interior mixing at
!       the lower boundary, i.e. at hbl.
!
! NOTE: Use smf_bgrid since this uses the primary smf array read in from 
! the coupler in ocean_core/ocean_sbc.F90 when using the FMS coupler.
!
! </DESCRIPTION>
!
subroutine vert_mix_kpp_mom4p0 (aidif, Time, Thickness, Velocity, T_prog, T_diag, Dens, &
                         swflx, sw_frac_zt, pme, river, visc_cbu, diff_cbt, hblt_depth, do_wave)

  real,                            intent(in)    :: aidif
  type(ocean_time_type),           intent(in)    :: Time
  type(ocean_thickness_type),      intent(in)    :: Thickness
  type(ocean_velocity_type),       intent(in)    :: Velocity
  type(ocean_prog_tracer_type),    intent(inout) :: T_prog(:)
  type(ocean_diag_tracer_type),    intent(in)    :: T_diag(:)
  type(ocean_density_type),        intent(in)    :: Dens
  real, dimension(isd:,jsd:),      intent(in)    :: swflx
  real, dimension(isd:,jsd:,:),    intent(in)    :: sw_frac_zt
  real, dimension(isd:,jsd:),      intent(in)    :: pme
  real, dimension(isd:,jsd:),      intent(in)    :: river
  real, dimension(isd:,jsd:),      intent(inout) :: hblt_depth
  real, dimension(isd:,jsd:,:),    intent(inout) :: visc_cbu
  real, dimension(isd:,jsd:,:,:),  intent(inout) :: diff_cbt
  logical,                         intent(in)    :: do_wave

  real, dimension(isd:ied,jsd:jed,nk) :: dbloc1, dbsfc1
  real, dimension(isd:ied,jsd:jed)    :: frazil
  real                                :: smftu, smftv, active_cells
  integer                             :: i, j, k, n, ki
  integer                             :: tau, taum1

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, '==>Error from vert_mix_kpp_mom4p0: module must be initialized')
  endif 

  if(.not. use_this_module) return 

  tau   = Time%tau
  taum1 = Time%taum1

  visc_cbu(:,:,:) = 0.0

  if(index_frazil > 0) then 
    do j = jsd, jed
       do i = isd, ied
          frazil(i,j) = T_diag(index_frazil)%field(i,j,1)
       enddo
    enddo
  else 
    frazil(:,:) = 0.0
  endif 

!---------assign Ustk2
    if (do_wave) then
       do j = jsd, jed
          do i = isd, ied
             Ustk2(i,j) = Velocity%ustoke(i,j)**2 + Velocity%vstoke(i,j)**2
          enddo
       enddo
    endif

!-----------------------------------------------------------------------
!     compute gradient Ri 
!-----------------------------------------------------------------------

  call ri_for_kpp(Time, Thickness, aidif, Velocity, T_prog(index_temp)%field(:,:,:,taum1), &
                  Dens%rho_salinity(:,:,:,taum1), Dens%pressure_at_depth,                  &
                  Dens%rho(:,:,:,tau), Dens%rho(:,:,:,taum1))

!-----------------------------------------------------------------------
!     compute vertical difference of velocity squared
!-----------------------------------------------------------------------

    do k=1,nk
      do j=jsd,jed
        do i=isd,ied
            dVsq(i,j,k) = (Velocity%u(i,j,1,1,tau) - Velocity%u(i,j,k,1,tau))**2   &
                        + (Velocity%u(i,j,1,2,tau) - Velocity%u(i,j,k,2,tau))**2
        enddo
      enddo
    enddo

!-----------------------------------------------------------------------
!     density related quantities
!     --------------------------
!     based on mom equation of state, which computes normalized
!     density  drho{k,kr} for layer k with respect to reference layer
!     kr, i.e. drho{T(k),S(k), tr(kr), sr(kr), zt(kr)}, as the
!     difference of drho=rho{k,kr}-rhor(kr):
!
!     density of surface layer                                   (kg/m3)
!           rho   = drho{T(1),S(1),zt(1)} + 1. + rhor(zt(1))
!     local buoyancy gradient at km interfaces:                  (m/s2)
!           dbloc = g/rho{k+1,k+1} * [ drho{k,k+1}-drho{k+1,k+1} ]
!     buoyancy difference with respect to "zref",i.e. the surface(m/s2)
!           dbsfc = g * [ drho{1,k}/rho{1,k} - drho{k,k}/rho{k,k} ]
!     thermal expansion coefficient without 1/rho factor       (kg/m3/C)
!           talpha= -d(rho{k,k})/d(T(k))
!     salt expansion coefficient without 1/rho factor        (kg/m3/PSU)
!           sbeta = d(rho{k,k})/d(S(k))
!-----------------------------------------------------------------------

        ! compute thermal and haline expansion coefficients (with factor of rho included).
        ! result from density_derivs is d(rho)/d(theta), yet need minus this for kpp   
        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 talpha(i,j,k) = -Dens%drhodT(i,j,k)
                 sbeta(i,j,k)  =  Dens%drhodS(i,j,k)
              enddo
           enddo
        enddo

        rhosfc(:,:) = Dens%rho(:,:,1,tau)
        dbsfc(:,:,1) = 0.0

        dbloc1(:,:,:) = density_delta_z(Dens%rho(:,:,:,tau), Dens%rho_salinity(:,:,:,tau),   &
                        T_prog(index_temp)%field(:,:,:,tau), Dens%pressure_at_depth(:,:,:))
        dbsfc1(:,:,:) = density_delta_sfc(Dens%rho(:,:,:,tau), Dens%rho_salinity(:,:,:,tau), &
                        T_prog(index_temp)%field(:,:,:,tau),   Dens%pressure_at_depth(:,:,:))
        do k=2,nk
           do j=jsd,jed
              do i=isd,ied
                 dbloc(i,j,k-1) = -grav * dbloc1(i,j,k-1)/(epsln+Dens%rho(i,j,k,tau))
                 dbsfc(i,j,k)   = -grav * dbsfc1(i,j,k-1)/(epsln+Dens%rho(i,j,k,tau))
              enddo
           enddo
        enddo
        dbloc(:,:,nk) = dbloc(:,:,nk-1)

!-----------------------------------------------------------------------
!     kinematic surface fluxes on the "t-grid" :
!     --------------------------------------------------------
!
!     wusfc = kinematic zonal velocity sfc flux on jp       (m2/s2)
!           = -rho0r*smf_bgrid(i,j,1)
!     wvsfc = kinematic merid velocity                      (m2/s2)
!           = -rho0r*smf_bgrid(i,j,2)
!     wsfc(temp) = kinematic temperature                    (C*m/s)
!                = rho0r*(stf(i,j,temp) - pme *(t(i,j,1,temp)-tpme(temp)) 
!                                       -river*(t(i,j,1,temp)-triver(temp))) 
!     wsfc(salt) = kinematic salinity                       (PSU*m/s)
!                = rho0r*(stf(i,j,salt) - pme  *(t(i,j,1,salt)-tpme(salt)) 
!                                       - river*(t(i,j,1,salt)-triver(salt))) 
!     wtsol = kinematic solar temp                          (C*m/s)
!           = -swflx(i,j)/rho_cp 
!
!     friction velocity, turbulent and radiative sfc buoyancy forcing
!     ---------------------------------------------------------------
!
!     ustar(i) = sqrt(sqrt(wusfc(i)**2 + wvsfc(i)**2))      (m/s)
!     Bo(i)    = grav*(talpha(i,j,1)*wsfc(temp)-sbeta(i,j)*wsfc(salt))/rho(i)
!     Bosol(i) = grav*(talpha(i,j,1)*wtsol(i)                         )/rho(i)
!                Both Bo and Bosol have units               (m2/s3)
!-----------------------------------------------------------------------

      ! stf, pme, and river are mass fluxes and so 
      ! require rho0r to convert to kinematic units.
      ! algebraically we have 
      ! rho0*wsfc = stf - pme*(T(1)-tpme) - runoff*(T(1)-trunoff) - calving*(T(1)-tcalving)
      !           = stf - pme*(T(1)-tpme) - (runoff + calving)*T(1) + runoff*trunoff + calving*tcalving
      !           = stf - pme*(T(1)-tpme) - river*T(1) + runoff_tracer_flux + calving_tracer_flux    
        if(wsfc_combine_runoff_calve) then
            ! older method 
            do n=1,num_prog_tracers
               do j=jsd,jed
                  do i=isd,ied
                     wsfc(n)%wsfc(i,j) = (T_prog(n)%stf(i,j) -                          &
                          pme(i,j)  *(T_prog(n)%field(i,j,1,tau)-T_prog(n)%tpme(i,j)) - &
                          river(i,j)*(T_prog(n)%field(i,j,1,tau)-T_prog(n)%triver(i,j)) &
                          )*rho0r*Grd%tmask(i,j,1)
                  enddo
               enddo
            enddo
        else 
            ! newer method 
            do n=1,num_prog_tracers
               do j=jsd,jed
                  do i=isd,ied
                     wsfc(n)%wsfc(i,j) = (T_prog(n)%stf(i,j)                                    &
                          -pme(i,j)  *(T_prog(n)%field(i,j,1,tau)-T_prog(n)%tpme(i,j))          &
                          -river(i,j)*T_prog(n)%field(i,j,1,tau)                                &
                          +T_prog(n)%runoff_tracer_flux(i,j)+T_prog(n)%calving_tracer_flux(i,j) &
                          )*rho0r*Grd%tmask(i,j,1)
                  enddo
               enddo
            enddo
        endif


      do j=jsc,jec
        do i=isc,iec

          ! ustar is needed on the "T-grid".  It is assumed that masking of 
          ! smf over land was performed inside of the ocean_sbc module. 
          ! smf has units of N/m^2 and so we need rho0r to get ustar in m/s.   
          ! swflx has units W/m^2 so needs 1/rho_cp to get to C*m/s units.
          ! these are proper units for buoyancy fluxes.
          ! use smf_bgrid for either MOM_BGRID or MOM_CGRID, since smf_bgrid
          ! is taken straight from the FMS coupler, so is more basic than smf_cgrid. 
          active_cells = Grd%umask(i,j,1)   + Grd%umask(i-1,j,1)   &
                        +Grd%umask(i,j-1,1) + Grd%umask(i-1,j-1,1) + epsln
          smftu = rho0r*(Velocity%smf_bgrid(i,j,1)   + Velocity%smf_bgrid(i-1,j,1)     &
                        +Velocity%smf_bgrid(i,j-1,1) + Velocity%smf_bgrid(i-1,j-1,1))  &
                  /active_cells
          smftv = rho0r*(Velocity%smf_bgrid(i,j,2)   + Velocity%smf_bgrid(i-1,j,2)    &
                        +Velocity%smf_bgrid(i,j-1,2) + Velocity%smf_bgrid(i-1,j-1,2)) &
                   /active_cells
          ustar(i,j) = sqrt( sqrt(smftu**2 + smftv**2) )
          
          Bo(i,j)    = grav * (talpha(i,j,1) * &
                 (wsfc(index_temp)%wsfc(i,j)+frazil(i,j)/(rho_cp*tracer_timestep))  &
                 -sbeta (i,j,1) * wsfc(index_salt)%wsfc(i,j)) &     
                  / (epsln + rhosfc(i,j))                     
          Bosol(i,j) = grav * talpha(i,j,1) * swflx(i,j)*inv_rho_cp &
                       / (epsln + rhosfc(i,j))

        enddo
      enddo

!-----------------------------------------------------------------------
!     compute interior mixing coefficients everywhere, due to constant 
!     internal wave activity, static instability, and local shear 
!     instability.
!-----------------------------------------------------------------------
      call ri_iwmix(Time, Dens%rho(:,:,:,tau), visc_cbu, diff_cbt)


!-----------------------------------------------------------------------
!     add double diffusion
!-----------------------------------------------------------------------

      if (double_diffusion) then
        call ddmix (Time, T_prog, Dens, diff_cbt)  
      endif

!-----------------------------------------------------------------------
! set seafloor values to zero for blmix 
!-----------------------------------------------------------------------
      do ki = 1,nk-1
        do j = jsc,jec
          do i = isc,iec
            visc_cbu(i,j,ki)    = visc_cbu(i,j,ki)  *Grd%tmask(i,j,min(ki+1,nk))
            diff_cbt(i,j,ki,1)  = diff_cbt(i,j,ki,1)*Grd%tmask(i,j,min(ki+1,nk))
            diff_cbt(i,j,ki,2)  = diff_cbt(i,j,ki,2)*Grd%tmask(i,j,min(ki+1,nk))
          enddo
        enddo
      enddo
!
!-----------------------------------------------------------------------
!     boundary layer mixing coefficients: diagnose new b.l. depth
!-----------------------------------------------------------------------

      call bldepth(sw_frac_zt, do_wave) 
 
!-----------------------------------------------------------------------
!     boundary layer diffusivities
!-----------------------------------------------------------------------

      call blmix_kpp(diff_cbt, visc_cbu, do_wave)

!-----------------------------------------------------------------------
!     enhance diffusivity at interface kbl - 1
!-----------------------------------------------------------------------

      call enhance(diff_cbt, visc_cbu) 

!-----------------------------------------------------------------------
!     combine interior and b.l. coefficients and nonlocal term
!-----------------------------------------------------------------------

      ! smooth diffusivities with a 5 point filter
      ! to remove 2 delta x noise
      if (smooth_blmc) then
          
          call mpp_update_domains(blmc,Dom%domain2d)

          do k=1,nk-1
             do j=jsc,jec
                do i=isc,iec
                   active_cells = 4.0*Grd%tmask(i,j,k) + &
                        Grd%tmask(i-1,j,k)  +&
                        Grd%tmask(i+1,j,k)  +&
                        Grd%tmask(i,j-1,k)  +&
                        Grd%tmask(i,j+1,k)
                   if (active_cells == 8.0) then
                       wrk1(i,j,k) = & 
                            (4.0*blmc(i,j,k,2)   +&
                                 blmc(i-1,j,k,2) +&
                                 blmc(i+1,j,k,2) +&
                                 blmc(i,j-1,k,2) +&
                                 blmc(i,j+1,k,2)) / active_cells
                       wrk2(i,j,k) =  &
                            (4.0*blmc(i,j,k,3) +&
                              blmc(i-1,j,k,3)  +&
                              blmc(i+1,j,k,3)  +&
                              blmc(i,j-1,k,3)  +&
                              blmc(i,j+1,k,3)) / active_cells
                   else
                       wrk1(i,j,k) = blmc(i,j,k,2)
                       wrk2(i,j,k) = blmc(i,j,k,3)
                   endif
                enddo
             enddo
          enddo
   
          do k=1,nk-1
             do j=jsc,jec
                do i=isc,iec
                   blmc(i,j,k,2) = wrk1(i,j,k)
                   blmc(i,j,k,3) = wrk2(i,j,k)
                enddo
             enddo
          enddo

      endif

      
      do ki=1,nk-1
         do j=jsc,jec
            do i=isc,iec               
               if (ki < kbl(i,j)) then
                   visc_cbu(i,j,ki)   = blmc(i,j,ki,1)
                   diff_cbt(i,j,ki,1) = blmc(i,j,ki,2)
                   diff_cbt(i,j,ki,2) = blmc(i,j,ki,3)
               else
                   ghats(i,j,ki)=0.0
               endif

               diff_cbt(i,j,ki,1) = diff_cbt(i,j,ki,1)*Grd%tmask(i,j,min(ki+1,nk))
               diff_cbt(i,j,ki,2) = diff_cbt(i,j,ki,2)*Grd%tmask(i,j,min(ki+1,nk))               
            enddo
         enddo
      enddo

!     update halo regions to compute spatial avg of viscosity to place onto U-cell 
      call mpp_update_domains (visc_cbu(:,:,:), Dom%domain2d, flags=NUPDATE+EUPDATE)

      
! compute viscosity on "U-grid" by averaging visc_cbu from "T-grid". note: zero out land values

      do ki=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               visc_cbu(i,j,ki)   = .25*Grd%umask(i,j,min(ki+1,nk))* &
                    (visc_cbu(i,j,ki) + visc_cbu(i+1,j,ki) +&
                    visc_cbu(i,j+1,ki) + &
                    visc_cbu(i+1,j+1,ki) )
            enddo
         enddo
             
      enddo
      
      visc_cbu(:,:,nk)   = 0.0
      diff_cbt(:,:,nk,1) = 0.0
      diff_cbt(:,:,nk,2) = 0.0

!-----------------------------------------------------------------------
!      boundary layer depth
!-----------------------------------------------------------------------

       do j=jsc,jec
         do i=isc,iec
              hblt(i,j)       = hbl(i,j) * Grd%tmask(i,j,1)
              hblt_depth(i,j) = hblt(i,j)
         enddo
       enddo

!-----------------------------------------------------------------------
! add non-local contribution to th_tendency term
! note th_tendency has units of (kg/m^3)*(m/s)*tracer 
!-----------------------------------------------------------------------

       wrk1(:,:,:)  = 0.0 ! initialize work array for diagnostic storage

       if (non_local_kpp) then

           do n=1,num_prog_tracers

              if (n==index_salt) then

                  do k=1,1
                     do j=jsc,jec
                        do i=isc,iec
                           wrk1(i,j,k) = rho0*wsfc(n)%wsfc(i,j)*(-diff_cbt(i,j,k,2)*ghats(i,j,k))
                           T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + wrk1(i,j,k)
                        enddo
                     enddo
                  enddo
                  do k=2,nk
                     do j=jsc,jec
                        do i=isc,iec
                           wrk1(i,j,k) = rho0*wsfc(n)%wsfc(i,j)*(diff_cbt(i,j,k-1,2)*ghats(i,j,k-1) - &
                                                                 diff_cbt(i,j,k,2)*ghats(i,j,k))
                           T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + wrk1(i,j,k)
                        enddo
                     enddo
                  enddo

              else 

                  do k=1,1
                     do j=jsc,jec
                        do i=isc,iec
                           wrk1(i,j,k) = rho0*wsfc(n)%wsfc(i,j)*(-diff_cbt(i,j,k,1)*ghats(i,j,k))
                           T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + wrk1(i,j,k)
                        enddo
                     enddo
                  enddo
                  do k=2,nk
                     do j=jsc,jec
                        do i=isc,iec
                           wrk1(i,j,k) =  rho0*wsfc(n)%wsfc(i,j)*(diff_cbt(i,j,k-1,1)*ghats(i,j,k-1) - &
                                                                  diff_cbt(i,j,k,1)*ghats(i,j,k))
                           T_prog(n)%th_tendency(i,j,k) =  T_prog(n)%th_tendency(i,j,k) + wrk1(i,j,k)
                        enddo
                     enddo
                  enddo

              endif   ! endif for n==index_salt 


              ! diagnostic for the nonlocal term
              if (id_nonlocal(n) > 0) then 
                 call diagnose_3d(Time, Grd, id_nonlocal(n),T_prog(n)%conversion*wrk1(:,:,:))
              endif 
              if (id_nonlocal_on_nrho(n) > 0) then
                 call diagnose_3d_rho(Time, Dens, id_nonlocal_on_nrho(n),T_prog(n)%conversion*wrk1)
              endif

              ! save nonlocal term for later diagnostics of 
              ! evolution of locally ref potrho and dianeutral velocity component 
              if(n==index_temp) then 
                  do k=1,nk
                     do j=jsc,jec
                        do i=isc,iec
                           T_prog(index_temp)%wrk1(i,j,k) = wrk1(i,j,k)
                        enddo
                     enddo
                  enddo
              endif
              if(n==index_salt) then 
                  do k=1,nk
                     do j=jsc,jec
                        do i=isc,iec
                           T_prog(index_salt)%wrk1(i,j,k) = wrk1(i,j,k)
                        enddo
                     enddo
                  enddo
              endif
                   

           enddo   ! enddo for n-loop 

           call watermass_diag(Time, T_prog, Dens)


       endif  ! endif for non_local_kpp

!-----------------------------------------------------------------------
!     send mixing related fields resulting just from kpp
!-----------------------------------------------------------------------

       call diagnose_3d(Time, Grd, id_diff_cbt_kpp_t, diff_cbt(:,:,:,1))
       call diagnose_3d(Time, Grd, id_diff_cbt_kpp_s, diff_cbt(:,:,:,2))
       call diagnose_2d(Time, Grd, id_hblt, hblt(:,:))
       call diagnose_3d(Time, Grd, id_bf_int_tide, bf_int_tide(:,:,:))

end subroutine vert_mix_kpp_mom4p0
! </SUBROUTINE> NAME="vert_mix_kpp_mom4p0"


!#######################################################################
! <SUBROUTINE NAME="bldepth">
!
! <DESCRIPTION>
! The oceanic planetray boundary layer depth, hbl, is determined as
! the shallowest depth where the bulk richardson number is
! equal to the critical value, Ricr.
!
! Bulk Richardson numbers are evaluated by computing velocity and
! buoyancy differences between values at zt(kl) and surface
! reference values.
!
! In this configuration, the reference values are equal to the
! values in the surface layer.  
! When using a very fine vertical grid, these values should be 
! computed as the vertical average of velocity and buoyancy from 
! the surface down to epsilon*zt(kl).
!
! When the bulk richardson number at k exceeds Ricr, hbl is
! linearly interpolated between grid levels zt(k) and zt(k-1).
!
! The water column and the surface forcing are diagnosed for 
! stable/ustable forcing conditions, and where hbl is relative 
! to grid points (caseA), so that conditional branches can be 
! avoided in later subroutines.
!
!  model  
!      real zt(1:nk)              = vertical grid (m)                           <BR/>
!      real dzt(1:nk)             = layer thicknesses (m)                       <BR/>  
!
!  input
!      real dbloc(ij_bounds,nk)   = local delta buoyancy         (m/s^2)        <BR/> 
!      real dbsfc(ij_bounds,nk)   = delta buoyancy w/ respect to sfc(m/s)^2     <BR/> 
!      real ustar(ij_bounds)      = surface friction velocity     (m/s)         <BR/>   
!      real Bo(ij_bounds)         = surface turbulent buoyancy forcing(m^2/s^3) <BR/> 
!      real Bosol(ij_bounds)      = radiative buoyancy forcing (m^2/s^3)        <BR/>  
!      real f(ij_bounds)          = Coriolis parameter            (1/s)         <BR/>
!      integer jwtype(ij_bounds)  = Jerlov water type           (1 to 5)        <BR/>   
!
!  output
!      real hbl(ij_bounds)        ! boundary layer depth              (m)       <BR/>  
!      real bfsfc(ij_bounds) !Bo+radiation absorbed to d=hbf*hbl(m^2/s^3)       <BR/> 
!      real stable(ij_bounds)     ! =1 in stable forcing; =0 unstable           <BR/> 
!      real caseA(ij_bounds)      ! =1 in case A, =0 in case B                  <BR/> 
!      integer kbl(ij_bounds)     ! index of first grid level below hbl         <BR/>
! </DESCRIPTION>
!
subroutine bldepth(sw_frac_zt, do_wave)

  real, intent(in), dimension(isd:,jsd:,:) :: sw_frac_zt   !3-D array of shortwave fract
  logical, intent(in)                      :: do_wave

  real            :: Ritop         ! numerator of bulk Richardson Number
  real            :: bvsq, Vtsq
  real            :: active_cells, dVsq_avg_tcell
  real            :: hekman, hmonob, hlimit
  
  integer         :: ka, ku, kl, ksave
  integer         :: i, j
  integer         :: iwscale_use_hbl_eq_zt

  real, parameter :: cekman = 0.7  ! constant for Ekman depth
  real, parameter :: cmonob = 1.0  ! constant for Monin-Obukhov depth

! find bulk Richardson number at every grid level until > Ricr
!
! note: the reference depth is -epsilon/2.*zt(k), but the reference
!       u,v,t,s values are simply the surface layer values,
!       and not the averaged values from 0 to 2*ref.depth,
!       which is necessary for very fine grids(top layer < 2m thickness)
! note: max values when Ricr never satisfied are
!       kbl(i)=kmt(i,j) and hbl(i)=zt(i,kmt(i,j))
!
!
!     initialize hbl and kbl to bottomed out values

        do j=jsc,jec
          do i=isc,iec
            Rib(i,j,1) = 0.0
            kbl(i,j)   = MAX(Grd%kmt(i,j),2)
            hbl(i,j)   = Grd%zt(kbl(i,j))
          enddo
        enddo

!       indices for array Rib(i,j,k), the bulk Richardson number.     
       
        ka = 1
        ku = 2

        do kl=2,nk

!       compute bfsfc = sw fraction at hbf * zt
        do j=jsc,jec
         do i=isc,iec
            bfsfc(i,j)  = Bo(i,j)  + Bosol(i,j) * (0.0 - sw_frac_zt(i,j,kl))  ! default assumes sw included in stf. 
                                                                              ! otherwise use (1.0-sw_frac)         
         
            stable(i,j) = 0.5 + SIGN( 0.5, bfsfc(i,j) )
            sigma(i,j)  = stable(i,j) * 1. + (1.-stable(i,j)) * epsilon
          enddo
        enddo

!-----------------------------------------------------------------------
!         compute velocity scales at sigma, for hbl = zt(kl):
!-----------------------------------------------------------------------

          iwscale_use_hbl_eq_zt=1
          call wscale ( iwscale_use_hbl_eq_zt, Grd%zt(kl), do_wave)


        do j=jsc,jec
          do i=isc,iec

!-----------------------------------------------------------------------
!           compute the turbulent shear contribution to Rib
!           eqn. (23)
!-----------------------------------------------------------------------

            bvsq =0.5*                                           &
                 ( dbloc(i,j,kl-1) / (Grd%dzw(kl-1))+             &
                   dbloc(i,j,kl  ) / (Grd%dzw(kl)) )
            Vtsq =  Grd%zt(kl) * ws(i,j) * sqrt(abs(bvsq)) * Vtc
 
!-----------------------------------------------------------------------
!           compute bulk Richardson number at new level
!           note: Ritop needs to be zero on land and ocean bottom
!           points so that the following if statement gets triggered
!           correctly. otherwise, hbl might get set to (big) negative
!           values, that might exceed the limit for the "exp" function
!           in "swfrac"
!
!           eqn. (21)
!
!
!           numerator of bulk richardson number on grid levels
!           --------------------------------------------------
!           note: land and ocean bottom values need to be set to zero
!           so that the subroutine "bldepth" works correctly
!-----------------------------------------------------------------------

            Ritop = (Grd%zt(kl)-Grd%zt(1)) * dbsfc(i,j,kl) * Grd%tmask(i,j,kl)

            active_cells   = Grd%umask(i,j,kl) + Grd%umask(i-1,j,kl) + Grd%umask(i,j-1,kl) + Grd%umask(i-1,j-1,kl) + epsln
            dVsq_avg_tcell = (dVsq(i,j,kl) + dVsq(i-1,j,kl) + dVsq(i,j-1,kl) + dVsq(i-1,j-1,kl))/active_cells
            Rib(i,j,ku) = Ritop / ( dVsq_avg_tcell + Vtsq + epsln )

            if((kbl(i,j) == Grd%kmt(i,j)).and.(Rib(i,j,ku) > Ricr)) then

!              linearly interpolate to find hbl where Rib = Ricr

               hbl(i,j) = Grd%zt(kl-1) + (Grd%dzw(kl-1)) * (Ricr - Rib(i,j,ka))  & 
                                     /(Rib(i,j,ku)-Rib(i,j,ka) + epsln)
               kbl(i,j) = kl
            endif

          enddo
        enddo

         ksave = ka
         ka    = ku
         ku    = ksave

      enddo

!-----------------------------------------------------------------------
!     find stability and buoyancy forcing for boundary layer
!-----------------------------------------------------------------------

      do j=jsc,jec
        do i = isc,iec
          ! Linear interpolation of sw_frac_zt to depth of hbl
          sw_frac_hbl(i,j)=(sw_frac_zt(i,j,kbl(i,j)) - sw_frac_zt(i,j,kbl(i,j)-1)) / Grd%dzt(kbl(i,j)) &
              * (hbl(i,j)-Grd%zt(kbl(i,j)))+sw_frac_zt(i,j,kbl(i,j))
          bfsfc(i,j)  = Bo(i,j) + Bosol(i,j) * (0.0 - sw_frac_hbl(i,j))
          stable(i,j) = 0.5 + SIGN( 0.5, bfsfc(i,j) )
          bfsfc(i,j)  = bfsfc(i,j) + stable(i,j) * epsln  ! ensures bfsfc never=0
        enddo
      enddo


!-----------------------------------------------------------------------
!        check hbl limits for hekman or hmonob
!        eqn. (24)
!-----------------------------------------------------------------------

      do j=jsc,jec
        do i = isc,iec
          if (bfsfc(i,j) > 0.0) then
             hekman = cekman * ustar(i,j) / (abs(Grd%f(i,j))+epsln)
             hmonob = cmonob * ustar(i,j)*ustar(i,j)*ustar(i,j)     &
                     /vonk / (bfsfc(i,j)+epsln) 
             hlimit = stable(i,j)     * AMIN1(hekman,hmonob) +      &
                     (stable(i,j)-1.) * (Grd%zt(nk))
             hbl(i,j) = AMIN1(hbl(i,j),hlimit)
             hbl(i,j) = AMAX1(hbl(i,j),Grd%zt(1))
          endif
          kbl(i,j) =Grd%kmt(i,j)
        enddo
      enddo

!-----------------------------------------------------------------------
!     find new kbl
!-----------------------------------------------------------------------

      do kl=2,nk
        do j=jsc,jec
          do i = isc,iec
            if((kbl(i,j) == Grd%kmt(i,j)).and.(Grd%zt(kl) > hbl(i,j))) then
              kbl(i,j) = kl
            endif
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     find stability and buoyancy forcing for final hbl values
!-----------------------------------------------------------------------

      do j=jsc,jec
        do i = isc,iec
          ! Linear interpolation of sw_frac_zt to depth of hbl
          sw_frac_hbl(i,j)=(sw_frac_zt(i,j,kbl(i,j)) - sw_frac_zt(i,j,kbl(i,j)-1)) / Grd%dzt(kbl(i,j)) &
              * (hbl(i,j)-Grd%zt(kbl(i,j)))+sw_frac_zt(i,j,kbl(i,j))
          bfsfc(i,j)  = Bo(i,j) + Bosol(i,j) * (0.0 - sw_frac_hbl(i,j))
          stable(i,j) = 0.5 + SIGN( 0.5, bfsfc(i,j) )
          bfsfc(i,j)  = bfsfc(i,j) + stable(i,j) * epsln 
        enddo
      enddo

!-----------------------------------------------------------------------
!     determine caseA and caseB
!     (if hbl is below (deeper than) the mid point of level kbl
!     then caseA=0  else caseA=1)
!-----------------------------------------------------------------------

      do j=jsc,jec
        do i = isc,iec
          caseA(i,j)  = 0.5 + SIGN( 0.5,Grd%zt(kbl(i,j)) - 0.5*Grd%dzt(kbl(i,j)) - hbl(i,j) )
        enddo
      enddo

end subroutine bldepth
! </SUBROUTINE> NAME="bldepth"



!#######################################################################
! <SUBROUTINE NAME="wscale">
!
! <DESCRIPTION>
! Compute turbulent velocity scales.
! Use a 2D-lookup table for wm and ws as functions of ustar and
! zetahat (=vonk*sigma*hbl*bfsfc).
!
! Note: the lookup table is only used for unstable conditions
! (zehat <= 0), in the stable domain wm (=ws) gets computed
! directly.
!
!  model
!
!  input                                                               <BR/>
!      real sigma(ij_bounds)  = normalized depth (d/hbl)               <BR/>
!      real hbl(ij_bounds)    = boundary layer depth (m)               <BR/>  
!      real ustar(ij_bounds)  = surface friction velocity    (m/s)     <BR/>  
!      real bfsfc(ij_bounds)  = total surface buoyancy flux (m^2/s^3)  <BR/>

!  output                                                               <BR/>  
!      real wm(ij_bounds),ws(ij_bounds) ! turbulent velocity scales at sigma

! local                                                                 <BR/>  
!      real zehat           ! = zeta *  ustar**3
!
! </DESCRIPTION>
!
subroutine wscale(iwscale_use_hbl_eq_zt, zt_kl, do_wave)

  real,    intent(in) :: zt_kl
  integer, intent(in) :: iwscale_use_hbl_eq_zt
  logical, intent(in) :: do_wave

  real                :: zdiff, udiff, zfrac, ufrac, fzfrac
  real                :: wam, wbm, was, wbs, u3, langmuirfactor, Cw_smyth
  real                :: zehat           ! = zeta *  ustar**3
  integer             :: iz, izp1, ju, jup1
  integer             :: i, j

!-----------------------------------------------------------------------
!     use lookup table for zehat < zmax only; otherwise use
!     stable formulae
!-----------------------------------------------------------------------

      do j=jsc,jec
        do i=isc,iec
          if (iwscale_use_hbl_eq_zt == 1) then
            zehat = vonk * sigma(i,j) * zt_kl * bfsfc(i,j)
          else
            zehat = vonk * sigma(i,j) * hbl(i,j) * bfsfc(i,j)
          end if

         if (zehat <= zmax) then

            zdiff  = zehat-zmin
            iz = int( zdiff/deltaz )
            iz = min( iz , nni )
            iz = max( iz , 0  )
            izp1=iz+1
            
            udiff  = ustar(i,j)-umin

            ju = int( min(udiff/deltau,float(nnj)))
            ju = max( ju , 0  )
            jup1=ju+1

            zfrac = zdiff/deltaz - float(iz)
            ufrac = udiff/deltau - float(ju)

            fzfrac= 1.-zfrac
            wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
            wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
            wm(i,j) = (1.-ufrac)* wbm          + ufrac*wam

            was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
            wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
            ws(i,j) = (1.-ufrac)* wbs          + ufrac*was
          else
            u3    = ustar(i,j)*ustar(i,j)*ustar(i,j)
            wm(i,j) = vonk * ustar(i,j) * u3 / ( u3 + conc1*zehat + epsln )
            ws(i,j) = wm(i,j)
          endif
        enddo
     enddo

!----------- if do_wave, add Langmuir turbulence enhancement factor

      if (do_wave .and. do_langmuir) then
         do j=jsc,jec
            do i=isc,iec
               Cw_smyth=Cw_0*(ustar(i,j)*ustar(i,j)*ustar(i,j)/(ustar(i,j)*ustar(i,j)*ustar(i,j) &
                    + Wstfac* 0.4 *bfsfc(i,j)*hbl(i,j) + epsln))**l_smyth
               langmuirfactor=sqrt(1+Cw_smyth*Ustk2(i,j)/(ustar(i,j)*ustar(i,j) + epsln))
               langmuirfactor = max(1.0, langmuirfactor)
               langmuirfactor = min(LTmax, langmuirfactor)
               ws(i,j)=ws(i,j)*langmuirfactor
               wm(i,j)=wm(i,j)*langmuirfactor
            enddo
         enddo
      endif

end subroutine wscale
! </SUBROUTINE> NAME="wscale"



!#######################################################################
! <SUBROUTINE NAME="ri_iwmix">
!
! <DESCRIPTION>
! Compute interior viscosity and diffusivity due 
! to shear instability (dependent on a local richardson number),
! to background internal wave activity, and 
! to static instability (local richardson number < 0).
!
!     inputs:
!
!      nk             = number of vertical levels                               <BR/>
!      visc_cbu_iw    = background "visc_cbu" (m**2/sec) due to internal waves  <BR/>
!      diff_cbt_iw    = background "diff_cbt" (m**2/sec) due to internal waves  <BR/>
!      visc_cbu_limit = largest "visc_cbu" in regions of gravitational          <BR/>
!                       instability (m**2/sec)                                  <BR/>
!      diff_cbt_limit = largest "diff_cbt" in regions of gravitational          <BR/>
!                       instability (m**2/sec)                       
!
!     outputs:
!
!      visc_cbu = viscosity coefficient at bottom of "u" cells (m**2/s)         <BR/>   
!      diff_cbt = diffusion coefficient at bottom of "t" cells (m**2/s)         <BR/>  
! </DESCRIPTION>
!
subroutine ri_iwmix(Time, rho, visc_cbu, diff_cbt)

  type(ocean_time_type),          intent(in)    :: Time
  real, dimension(isd:,jsd:,:),   intent(in)    :: rho
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbu
  real, dimension(isd:,jsd:,:,:), intent(inout) :: diff_cbt
  
  real, parameter :: Riinfty = 0.8  ! local Richardson Number limit for shear instability
  real            :: Rigg, ratio, frit, fcont, frit_tide, rigg_tide
  real            :: q_int_tide, gamma_int_tide
  real            :: kv_up, kv_dn
  real            :: kv_visc_int_tide, kv_diff_int_tide
  integer         :: i, j, k, kmax

  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0

!-----------------------------------------------------------------------
!     diffusion and viscosity coefficients on bottom of "T" cells
!-----------------------------------------------------------------------

      do k=1,nk-1
        do j=jsc,jec
          do i=isc,iec
  
!-----------------------------------------------------------------------
!           evaluate function of Ri for shear instability eqn. (28b&c)
!-----------------------------------------------------------------------

            Rigg  = AMAX1( rit(i,j,k) , 0.0)
            ratio = AMIN1( Rigg/Riinfty , 1.0 )
            frit  = (1. - ratio*ratio)
            frit  = frit*frit*frit

!-----------------------------------------------------------------------
!           evaluate function of Ri by Munk-Anderson scheme
!-----------------------------------------------------------------------

            rigg_tide  = AMAX1( rit_tide(i,j,k) , 0.0)
            frit_tide = (1.0 + sigma_tide*rigg_tide)**p_tide   

!-----------------------------------------------------------------------
!           evaluate function of Ri for convection eqn. (28a)
!-----------------------------------------------------------------------

            fcont  = 0.5 * ( abs(rit(i,j,k)) - rit(i,j,k) )
            fcont  = fcont / (AMAX1(fcont,epsln))
     
!-----------------------------------------------------------------------
!           mixing due to internal wave activity and static instability 
!           eqn. (29).  Note: temporarily have visc_cbu on T-point.
!-----------------------------------------------------------------------

            visc_cbu(i,j,k)       = visc_cbu_iw + fcont * visc_con_limit   
            diff_cbt(i,j,k,1)     = diff_cbt_iw + fcont * diff_con_limit
            diff_cbt(i,j,k,2)     = diff_cbt_iw + fcont * diff_con_limit

!-----------------------------------------------------------------------
!           add contribution due to shear instability and coastal tide mixing
!-----------------------------------------------------------------------

            wrk1(i,j,k) = coastal_tidal_flag*diff_cbt_limit*frit_tide
            wrk2(i,j,k) = coastal_tidal_flag*visc_cbu_limit*frit_tide

            diff_cbt(i,j,k,1) = diff_cbt(i,j,k,1)              &
                 + shear_instability_flag*diff_cbt_limit*frit  &
                 + wrk1(i,j,k)
            diff_cbt(i,j,k,2) = diff_cbt(i,j,k,2)              &
                 + shear_instability_flag*diff_cbt_limit*frit  &
                 + wrk1(i,j,k)
            visc_cbu(i,j,k)   = visc_cbu(i,j,k)                &
                 + shear_instability_flag*visc_cbu_limit*frit  &
                 + wrk2(i,j,k)

          enddo
        enddo
      enddo

      ! internal tidal mixing where water depth is deeper than 100 m.
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0
      if (shear_instability .and. int_tidal_mix) then
          do j=jsd,jed
             do i=isd,ied
                kmax = Grd%kmt(i,j)
                if (kmax > 2 .and. Grd%zw(kmax)>100.0) then
                    do k=1,kmax

                       q_int_tide=0.3333333 
                       gamma_int_tide=0.2
                       kv_up=q_int_tide*gamma_int_tide*eflux_int_tide_bf(i,j)*bf_int_tide(i,j,kmax-1)*f_int_tide(i,j,k)
                       kv_dn=rho(i,j,k)*bf_int_tide(i,j,k)*bf_int_tide(i,j,k)+epsln
                       kv_visc_int_tide=AMIN1(kv_up/kv_dn,visc_cbu_limit)
                       kv_diff_int_tide=AMIN1(kv_up/kv_dn,diff_cbt_limit)

                       visc_cbu(i,j,k)   = visc_cbu(i,j,k)   +  kv_visc_int_tide
                       diff_cbt(i,j,k,1) = diff_cbt(i,j,k,1) +  kv_diff_int_tide
                       diff_cbt(i,j,k,2) = diff_cbt(i,j,k,2) +  kv_diff_int_tide

                       wrk3(i,j,k)       = kv_diff_int_tide
                       wrk4(i,j,k)       = kv_visc_int_tide

                    enddo
                endif
             enddo
          enddo
      endif


      call diagnose_3d(Time, Grd, id_diff_cbt_coast, wrk1(:,:,:))
      call diagnose_3d(Time, Grd, id_visc_cbt_coast, wrk2(:,:,:))

      call diagnose_3d(Time, Grd, id_diff_cbt_int, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, id_visc_cbt_coast, wrk4(:,:,:))

  
end subroutine ri_iwmix
! </SUBROUTINE> NAME="ri_iwmix"


!#######################################################################
! <SUBROUTINE NAME="ddmix">
!
! <DESCRIPTION>
! Rrho dependent interior flux parameterization.
! Add double-diffusion diffusivities to Ri-mix values at blending
! interface and below. salt fingering code modified july 2003
! by stephen.griffies based on NCAR CCSM2.x
!
!     inputs: 
!
!      nk     = number of vertical levels                                 <BR/>
!      real talpha(imt,km,jmw)   ! d(rho)/ d(pot.temperature) (kg/m^3/C)  <BR/> 
!      real sbeta(imt,km,jmw)    ! d(rho)/ d(salinity)     (kg/m^3/PSU)  
!
!     outputs:
!
!      diff_cbt = diffusion coefficient at bottom of "t" cells (m**2/s)
!
!
! local
!
!      real alphaDT(imt,km,jmw)   ! alpha * DT  across interfaces          <BR/>   
!      real betaDS(imt,km,jmw)    ! beta  * DS  across interfaces          <BR/>  
!
! </DESCRIPTION>
!
subroutine ddmix (Time, T_prog, Dens, diff_cbt)

type(ocean_time_type),          intent(in)    :: Time
type(ocean_prog_tracer_type),   intent(in)    :: T_prog(:)
type(ocean_density_type),       intent(in)    :: Dens
real, dimension(isd:,jsd:,:,:), intent(inout) :: diff_cbt

real diffdd   ! double diffusion diffusivity scale
real prandtl  ! prandtl number
real Rrho     ! double-diffusive density ratio

real, parameter :: Rrho0               = 1.9    ! limit for double diffusive density ratio
real, parameter :: dsfmax              = 1.e-4  ! (m^2/s) max diffusivity in case of salt fingering
real, parameter :: viscosity_molecular = 1.5e-6 ! (m^2/s)

integer :: i, j, k, ki, tau

tau = Time%tau

!     for double diffusion
!     --------------------
!     temperature and salt contributions of buoyancy gradients
!     on interfaces

    do k=1,nk-1
      do j=jsc,jec
        do i=isc,iec
            alphaDT(i,j,k) = 0.5 * (talpha(i,j,k) + talpha(i,j,k+1) ) &
                             * (T_prog(index_temp)%field(i,j,k,tau) - T_prog(index_temp)%field(i,j,k+1,tau))
            betaDS (i,j,k) = 0.5 * (sbeta (i,j,k) + sbeta(i,j,k+1)  ) &
                             *(Dens%rho_salinity(i,j,k,tau) - Dens%rho_salinity(i,j,k+1,tau))
        enddo
      enddo
    enddo
    do j=jsc,jec
      do i=isc,iec
        alphaDT(i,j,nk) = 0.0
        betaDS (i,j,nk) = 0.0
      enddo
    enddo

    do ki=1,nk
       do j=jsc,jec
          do i=isc,iec

!-----------------------------------------------------------------------
!           salt fingering case
!           eqn. (31)
!-----------------------------------------------------------------------

             if ((alphaDT(i,j,ki) > betaDS(i,j,ki)) .and. &
                 (betaDS (i,j,ki) > 0.         )) then

                 Rrho   = MIN(alphaDT(i,j,ki) / betaDS(i,j,ki), Rrho0)

!                diffdd = dsfmax*(1.0-((Rrho-1)/(Rrho0-1))**2)**pexp2     ! (very old code) 
!                diffdd = 1.0-((Rrho-1)/(Rrho0-1))**2                     ! (less old code)
                 diffdd = 1.0-((Rrho-1.0)/(Rrho0-1.0))                    ! (new code)
                 diffdd = dsfmax*diffdd*diffdd*diffdd

                 diff_cbt(i,j,ki,1) = diff_cbt(i,j,ki,1) + 0.7*diffdd
                 diff_cbt(i,j,ki,2) = diff_cbt(i,j,ki,2) + diffdd

!            diffusive convection eqn. (32)
             else if ( (alphaDT(i,j,ki) < 0.0) .and.  &
                       (betaDS(i,j,ki)  < 0.0) .and.  &
                       (alphaDT(i,j,ki) > betaDS(i,j,ki)) ) then

                 Rrho    = alphaDT(i,j,ki) / betaDS(i,j,ki) 
                 diffdd  = viscosity_molecular*0.909*exp(4.6*exp(-0.54*(1/Rrho-1)))  

!                eqn. (34)
                 prandtl = 0.15*Rrho
                 if (Rrho > 0.5) prandtl = (1.85-0.85/Rrho)*Rrho

                 diff_cbt(i,j,ki,1) = diff_cbt(i,j,ki,1) + diffdd
                 diff_cbt(i,j,ki,2) = diff_cbt(i,j,ki,2) + prandtl*diffdd

             endif

          enddo
       enddo
    enddo

end subroutine ddmix
! </SUBROUTINE> NAME="ddmix"


!#######################################################################
! <SUBROUTINE NAME="blmix_kpp">
!
! <DESCRIPTION>
! Mixing coefficients within boundary layer depend on surface
! forcing and the magnitude and gradient of interior mixing below
! the boundary layer ("matching").
!
! CAUTION: if mixing bottoms out at hbl = zt(nk) then
! fictitious layer at nk+1 is needed with small but finite width 
! dzt(nk+1) (eg. epsln = 1.e-20).  
!
!     inputs:
!
!      real ustar(ij_bounds)    ! surface friction velocity         (m/s)  <BR/>
!      real bfsfc(ij_bounds)    ! surface buoyancy forcing     (m^2/s^3)   <BR/>
!      real hbl(ij_bounds)      ! boundary layer depth              (m)    <BR/>
!      real stable(ij_bounds)   ! = 1 in stable forcing                    <BR/>
!      real caseA(ij_bounds)    ! = 1 in case A                            <BR/>
!      integer kbl(ij_bounds)   ! index of first grid level below hbl
!
!     outputs:
!
!      visc_cbu = viscosity coefficient at bottom of "u" cells (m**2/s)    <BR/>     
!      diff_cbt = diffusion coefficient at bottom of "t" cells (m**2/s)    <BR/>
!
!      real dkm1(ij_bounds,3)    = boundary layer diff_cbt at kbl-1 level  <BR/> 
!      real blmc(ij_bounds,nk,3) = boundary layer mixing coeff.(m**2/s)    <BR/> 
!      real ghats(ij_bounds,nk)  = nonlocal scalar transport               <BR/> 
!
!    local:
!
!      real gat1(ij_bounds,3)                                               <BR/>    
!      real dat1(ij_bounds,3)                                               <BR/> 
!      real sigma(ij_bounds)              = normalized depth (d / hbl)      <BR/> 
!      real ws(ij_bounds), wm(ij_bounds)  = turbulent velocity scales (m/s) 
!
! </DESCRIPTION>
!
subroutine blmix_kpp(diff_cbt, visc_cbu, do_wave)

  real, dimension(isd:ied,jsd:jed,nk,2) :: diff_cbt
  real, dimension(isd:ied,jsd:jed,nk)   :: visc_cbu
  logical,intent(in)                    :: do_wave

  real     :: zt_kl_dummy
  real     :: delhat, R, dvdzup, dvdzdn
  real     :: viscp, difsp, diftp, visch, difsh, difth, f1
  real     :: sig, a1, a2, a3, Gm, Gs, Gt

  integer  :: iwscale_use_hbl_eq_zt
  integer  :: i, j, ki
  integer  :: kn, knm1, knp1

!-----------------------------------------------------------------------
!       compute velocity scales at hbl. (recall epsilon=0.1)
!-----------------------------------------------------------------------
        do j = jsc,jec
          do i = isc,iec
            sigma(i,j) = stable(i,j) * 1.0 + (1.-stable(i,j)) * epsilon
          enddo
        enddo

        iwscale_use_hbl_eq_zt = 0
        zt_kl_dummy=0.0

        call wscale (iwscale_use_hbl_eq_zt, zt_kl_dummy, do_wave)

      do j=jsc,jec
        do i = isc,iec 
          kn    = ifix(caseA(i,j)+epsln) *(kbl(i,j) -1) +   &
                  (1-ifix(caseA(i,j)+epsln)) * kbl(i,j)
          knm1 = max(kn-1,1)
          knp1 = min(kn+1,nk)

!-----------------------------------------------------------------------
!         find the interior viscosities and derivatives at hbl(i) 
!         eqn. (18)
!-----------------------------------------------------------------------

          delhat = 0.5*Grd%dzt(kn) + Grd%zt(kn) - hbl(i,j)
          R      = 1.0 - delhat / Grd%dzt(kn)
          dvdzup = (visc_cbu(i,j,knm1) - visc_cbu(i,j,kn))/Grd%dzt(kn)
          dvdzdn = (visc_cbu(i,j,kn) - visc_cbu(i,j,knp1))/Grd%dzt(knp1)
          viscp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+         &
                               R  * (dvdzdn + abs(dvdzdn)) )

          dvdzup = (diff_cbt(i,j,knm1,2) - diff_cbt(i,j,kn,2))/Grd%dzt(kn)
          dvdzdn = (diff_cbt(i,j,kn,2)   - diff_cbt(i,j,knp1,2))/Grd%dzt(knp1)
          difsp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+         &
                               R  * (dvdzdn + abs(dvdzdn)) )

          dvdzup = (diff_cbt(i,j,knm1,1) - diff_cbt(i,j,kn,1))/Grd%dzt(kn)
          dvdzdn = (diff_cbt(i,j,kn,1)   - diff_cbt(i,j,knp1,1))/Grd%dzt(knp1)
          diftp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+         &
                               R  * (dvdzdn + abs(dvdzdn)) )

          visch  = visc_cbu(i,j,kn)       + viscp * delhat
          difsh  = diff_cbt(i,j,kn,2) + difsp * delhat
          difth  = diff_cbt(i,j,kn,1)     + diftp * delhat

          f1 = stable(i,j) * conc1 * bfsfc(i,j) / (ustar(i,j)**4+epsln) 

          gat1(i,j,1) = visch / (hbl(i,j)+epsln) / (wm(i,j)+epsln)
          dat1(i,j,1) = -viscp / (wm(i,j)+epsln) + f1 * visch
          dat1(i,j,1) = min(dat1(i,j,1),0.) 

          gat1(i,j,2) = difsh  / (hbl(i,j)+epsln) / (ws(i,j)+epsln)
          dat1(i,j,2) = -difsp / (ws(i,j)+epsln) + f1 * difsh 
          dat1(i,j,2) = min(dat1(i,j,2),0.) 
  
          gat1(i,j,3) = difth /  (hbl(i,j)+epsln) / (ws(i,j)+epsln)
          dat1(i,j,3) = -diftp / (ws(i,j)+epsln) + f1 * difth 
          dat1(i,j,3) = min(dat1(i,j,3),0.) 
        enddo
      enddo

      blmc(:,:,:,:) = 0.0

      do ki=1,nk

!-----------------------------------------------------------------------
!         compute turbulent velocity scales on the interfaces
!-----------------------------------------------------------------------

        do j  = jsc,jec
          do i  = isc,iec
            sig     = (Grd%zt(ki) + 0.5 * Grd%dzt(ki)) / (hbl(i,j)+epsln)
            sigma(i,j) = stable(i,j)*sig                          &
                     + (1.-stable(i,j))*AMIN1(sig,epsilon)
          enddo
        enddo

          iwscale_use_hbl_eq_zt=0
          zt_kl_dummy=0.0

          call wscale(iwscale_use_hbl_eq_zt, zt_kl_dummy, do_wave)

!-----------------------------------------------------------------------
!         compute the dimensionless shape functions at the interfaces
!         eqn. (11)
!-----------------------------------------------------------------------

          
        do j=jsc,jec
          do i = isc,iec
            if (ki < kbl(i,j)) then
              sig = (Grd%zt(ki) + 0.5 * Grd%dzt(ki)) / (hbl(i,j)+epsln)
              a1 = sig - 2.
              a2 = 3.-2.*sig
              a3 = sig - 1.

              Gm = a1 + a2 * gat1(i,j,1) + a3 * dat1(i,j,1)
                   
              Gs = a1 + a2 * gat1(i,j,2) + a3 * dat1(i,j,2)
                  
              Gt = a1 + a2 * gat1(i,j,3) + a3 * dat1(i,j,3)
                  

!-----------------------------------------------------------------------
!             compute boundary layer diffusivities at the interfaces
!             eqn. (10)
!-----------------------------------------------------------------------

              
              blmc(i,j,ki,1) = hbl(i,j) * wm(i,j) * sig * (1.+sig*Gm) 
              blmc(i,j,ki,2) = hbl(i,j) * ws(i,j) * sig * (1.+sig*Gt) 
              blmc(i,j,ki,3) = hbl(i,j) * ws(i,j) * sig * (1.+sig*Gs) 

!-----------------------------------------------------------------------
!             nonlocal transport term = ghats * <ws>o (eqn. 20)
!             To include Langmuir turbulence effects, multiply ghats
!             by a factor of Lgam (McWilliam & Sullivan 2001)
!-----------------------------------------------------------------------
              if (do_wave .and. do_langmuir) then
                 ghats(i,j,ki) = Lgam * (1.-stable(i,j)) * cg    &
                        / (ws(i,j) * hbl(i,j) + epsln)
              else
                 ghats(i,j,ki) = (1.-stable(i,j)) * cg    &
                        / (ws(i,j) * hbl(i,j) + epsln)
              endif
            endif
          enddo
        enddo
      enddo

     
!-----------------------------------------------------------------------
!     find diffusivities at kbl-1 grid level 
!-----------------------------------------------------------------------

      do j=jsc,jec
        do i= isc,iec
          sig      =  Grd%zt(kbl(i,j)-1)  / (hbl(i,j)+epsln)
          sigma(i,j) =stable(i,j) * sig                  &
                   + (1.-stable(i,j)) * MIN(sig,epsilon)
        enddo

        iwscale_use_hbl_eq_zt=0
        zt_kl_dummy=0.0

        call wscale(iwscale_use_hbl_eq_zt, zt_kl_dummy, do_wave)

        do i = isc,iec
          sig = Grd%zt(kbl(i,j)-1) / (hbl(i,j)+epsln)
          a1= sig - 2.
          a2 = 3.-2.*sig
          a3 = sig - 1.
          Gm = a1 + a2 * gat1(i,j,1) + a3 * dat1(i,j,1)
          Gs = a1 + a2 * gat1(i,j,2) + a3 * dat1(i,j,2)
          Gt = a1 + a2 * gat1(i,j,3) + a3 * dat1(i,j,3)
          dkm1(i,j,1) = hbl(i,j) * wm(i,j) * sig * (1. + sig * Gm)
          dkm1(i,j,2) = hbl(i,j) * ws(i,j) * sig * (1. + sig * Gs)
          dkm1(i,j,3) = hbl(i,j) * ws(i,j) * sig * (1. + sig * Gt)
        enddo
      enddo
      
end subroutine blmix_kpp
! </SUBROUTINE> NAME="blmix_kpp"



!#######################################################################
! <SUBROUTINE NAME="enhance">
!
! <DESCRIPTION>
! Enhance the diffusivity at the kbl-.5 interface
!
! input
!      integer kbl(ij_bounds)   =  grid above hbl                      <BR/>
!      real hbl(ij_bounds)      =  boundary layer depth (m)            <BR/>
!      real dkm1(ij_bounds,3)   =  bl diffusivity at kbl-1 grid level  <BR/>
!      real caseA(ij_bounds)    =  1 in caseA, = 0 in case B
!
! input/output
!      real ghats(ij_bounds,nk)  =  nonlocal transport     (s/m**2)     <BR/>
!      modified ghats at kbl(i)-1 interface        
! output
!      real blmc(ij_bounds,nk,3) = enhanced boundary layer mixing coefficient
!
! local
!  real delta                    =  fraction hbl lies beteen zt neighbors
!
! </DESCRIPTION>
subroutine enhance(diff_cbt, visc_cbu)

real, dimension(isd:ied,jsd:jed,nk,2) :: diff_cbt
real, dimension(isd:ied,jsd:jed,nk)   :: visc_cbu

  real    :: delta, dkmp5, dstar
  integer :: i, j, ki

      do ki=1,nk-1
        do j=jsc,jec
          do i = isc,iec

            if(ki == (kbl(i,j) - 1) ) then

              delta = (hbl(i,j)-Grd%zt(ki)) / (Grd%zt(ki+1)-Grd%zt(ki))

              dkmp5 = caseA(i,j) * visc_cbu(i,j,ki)      &
                    + (1.-caseA(i,j)) * blmc(i,j,ki,1)
              dstar = (1.-delta)**2 * dkm1(i,j,1) + delta**2 * dkmp5
              blmc(i,j,ki,1) = (1.-delta) * visc_cbu(i,j,ki)  &
                             + delta * dstar

              ! temperature:
              dkmp5 = caseA(i,j) * diff_cbt(i,j,ki,1)  &
                    + (1.-caseA(i,j)) * blmc(i,j,ki,2)
              dstar = (1.-delta)**2 * dkm1(i,j,3) + delta**2 * dkmp5    
              blmc(i,j,ki,2) = (1.-delta) * diff_cbt(i,j,ki,1)  &
                           + delta * dstar

              ! salinity:   
              dkmp5 = caseA(i,j) * diff_cbt(i,j,ki,2)  &
                    + (1.-caseA(i,j)) * blmc(i,j,ki,3)
              dstar = (1.-delta)**2 * dkm1(i,j,2) + delta**2 * dkmp5
              blmc(i,j,ki,3) = (1.-delta) * diff_cbt(i,j,ki,2)  &
                           + delta * dstar
            
              ghats(i,j,ki) = (1.-caseA(i,j)) * ghats(i,j,ki)

            endif

          enddo
        enddo
      enddo

end subroutine enhance
! </SUBROUTINE> NAME="enhance"


!#######################################################################
! <SUBROUTINE NAME="ri_for_kpp">
!
! <DESCRIPTION>
! Compute Richardson number on tracer and velocity cell bottoms. 
!  rit = richardson number at bottom of T cells <BR/>
!  riu = richardson number at bottom of U cells
! </DESCRIPTION>
!
subroutine ri_for_kpp (Time, Thickness, aidif, Velocity, theta, salinity, &
                       pressure_at_depth, rho, rho_taum1)

  type(ocean_time_type),        intent(in)  :: Time
  type(ocean_thickness_type),   intent(in)  :: Thickness
  real,                         intent(in)  :: aidif
  type(ocean_velocity_type),    intent(in)  :: Velocity
  real, dimension(isd:,jsd:,:), intent(in)  :: theta
  real, dimension(isd:,jsd:,:), intent(in)  :: salinity 
  real, dimension(isd:,jsd:,:), intent(in)  :: pressure_at_depth
  real, dimension(isd:,jsd:,:), intent(in)  :: rho
  real, dimension(isd:,jsd:,:), intent(in)  :: rho_taum1
  
  integer :: i, j, k, tlev, kmax, mr, tau, taum1
  real    :: active_cells, fx, t1, riu_prev, tmp
  real    :: coef1_int_tide, drho_int_tide

  logical :: smooth_richardson_number = .true.
  integer :: num_smoothings = 1 ! for vertical smoothing of Richardson number


  fx   = -0.25*grav*rho0r

  tau   = Time%tau
  taum1 = Time%taum1

  ! compute density difference across bottom of T cells at tau-1

  if (aidif == 1.0) then
    tlev = tau
    wrk1(:,:,:) = density_delta_z(rho(:,:,:), salinity(:,:,:), theta(:,:,:), pressure_at_depth(:,:,:))
  else
    tlev = taum1
    wrk1(:,:,:) = density_delta_z(rho_taum1(:,:,:), salinity(:,:,:), theta(:,:,:), pressure_at_depth(:,:,:))
  endif

  ! Buoyancy frequency for internal tidal mixing at the bottom of t-cell
  ! bf_int_tide, low limit is zero. no negative frequency here.
  if (int_tidal_mix) then
      do k=1,nk-1
         do j=jsd,jed
            do i=isd,ied
               coef1_int_tide=-1.*grav/(rho(i,j,k)+epsln)
               drho_int_tide=AMIN1( wrk1(i,j,k), 0.0)
               bf_int_tide(i,j,k)=SQRT( coef1_int_tide*drho_int_tide/Thickness%dzwt(i,j,k) )
            enddo
         enddo
      enddo
  endif
     


  ! compute richardson numbers on bottom of U cells due to tide
  if(coastal_tidal_mix) then 
      do k=1,nk-1
         do j=jsd,jed-1
            do i=isd,ied-1
               t1 = fx*Thickness%dzwt(i,j,k)
               riu_tide(i,j,k) = t1*Grd%umask(i,j,k+1)                                    &
                    *(wrk1(i,j+1,k) + wrk1(i+1,j+1,k) + wrk1(i,j,k)   + wrk1(i+1,j,k)) /  &
                    ( tidal_fric_turb(i,j,k) + epsln)
            enddo
         enddo
      enddo
  endif

  ! compute richardson numbers on bottom of U cells

  do k=1,nk-1
     do j=jsd,jed-1
        do i=isd,ied-1
           t1 = fx*Thickness%dzwt(i,j,k)
        riu(i,j,k) = t1*Grd%umask(i,j,k+1)*(wrk1(i,j+1,k) + wrk1(i+1,j+1,k) + wrk1(i,j,k)   + wrk1(i+1,j,k)) /&
                ((Velocity%u(i,j,k,1,tlev) - Velocity%u(i,j,k+1,1,tlev))**2 + &
                (Velocity%u(i,j,k,2,tlev) - Velocity%u(i,j,k+1,2,tlev))**2  &
                + epsln)
        enddo
     enddo
  enddo

  if (smooth_richardson_number) then

  ! smooth Richardson number in the vertical using a 1-2-1 filter:

      do mr = 1,num_smoothings
         do j=jsd,jed-1
            do i=isd,ied-1
               riu_prev    =  0.25 * riu(i,j,1)
               kmax = Grd%kmt(i,j)
               if (kmax > 0) then
                   do k=2,kmax-2
                      tmp        =  riu(i,j,k)
                      riu(i,j,k) = riu_prev + 0.5 * riu(i,j,k) + 0.25 * riu(i,j,k+1)
                      riu_prev   =  0.25 * tmp
                   enddo
               endif
            enddo
         enddo
      enddo
  endif

  ! compute richardson numbers on bottom of T cells as average
  ! of four nearest richardson numbers on bottom of U cells.
  ! (do not consider land cells in the average... only active ones)

  do k=1,nk-1
    do j=jsc,jec
      do i=isc,iec
        active_cells = Grd%umask(i,j,k+1) + Grd%umask(i-1,j,k+1) + Grd%umask(i,j-1,k+1) + Grd%umask(i-1,j-1,k+1) + epsln
        rit(i,j,k)   = (riu(i,j,k) + riu(i-1,j,k) + riu(i,j-1,k) + riu(i-1,j-1,k))/active_cells

        rit_tide(i,j,k)   = (riu_tide(i,j,k) + riu_tide(i-1,j,k) + riu_tide(i,j-1,k) + riu_tide(i-1,j-1,k))/active_cells

        ! make sure no static instability exists (one that is not seen
        ! by the Richardson number).  This may happen due to
        ! horizontal averaging used in calculating the Richardson
        ! number.

        if (rit(i,j,k) > 0.0 .and. wrk1(i,j,k) > 0.0) then
          rit(i,j,k) = -10.0
        endif
        if (rit_tide(i,j,k) > 0.0 .and. wrk1(i,j,k) > 0.0) then
          rit_tide(i,j,k) = -10.0
        endif
      enddo
    enddo
  enddo

end subroutine ri_for_kpp
! </SUBROUTINE> NAME="ri_for_kpp"


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


  id_neut_rho_kpp_nloc = register_diag_field ('ocean_model',      & 
    'neut_rho_kpp_nloc',                                          &
    Grd%tracer_axes(1:3), Time%model_time,                        &
    'locally referenced potrho tendency due to KPP nonlocal term',&
    '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e10,1e10/))
  if(id_neut_rho_kpp_nloc > 0) compute_watermass_diag=.true.

  id_wdian_rho_kpp_nloc = register_diag_field ('ocean_model',   &
    'wdian_rho_kpp_nloc', Grd%tracer_axes(1:3), Time%model_time,&
    'dianeutral mass transport due to KPP nonlocal term',       &
    'kg/sec',missing_value=missing_value, range=(/-1e10,1e10/))
  if(id_wdian_rho_kpp_nloc > 0) compute_watermass_diag=.true.

  id_tform_rho_kpp_nloc = register_diag_field ('ocean_model',               & 
    'tform_rho_kpp_nloc', Grd%tracer_axes(1:3), Time%model_time,            &
    'watermass transform due to KPP nonlocal on levels (pre-layer binning)',&
    'kg/sec',missing_value=missing_value, range=(/-1e10,1e10/))
  if(id_tform_rho_kpp_nloc > 0) compute_watermass_diag=.true.

  id_neut_rho_kpp_nloc_on_nrho = register_diag_field ('ocean_model',                    &
    'neut_rho_kpp_nloc_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,            &
    'update of locally ref potrho from KPP nolocal as binned to neutral density layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_rho_kpp_nloc_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_rho_kpp_nloc_on_nrho = register_diag_field ('ocean_model',                   &
    'wdian_rho_kpp_nloc_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,           &
    'dianeutral mass transport due to KPP nonlocal as binned to neutral density layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_rho_kpp_nloc_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_rho_kpp_nloc_on_nrho = register_diag_field ('ocean_model',             &
    'tform_rho_kpp_nloc_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,     &
    'watermass transform due to KPP nonlocal as binned to neutral density layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_rho_kpp_nloc_on_nrho > 0) compute_watermass_diag=.true.

  id_eta_tend_kpp_nloc= -1          
  id_eta_tend_kpp_nloc= register_diag_field ('ocean_model','eta_tend_kpp_nloc',&
       Grd%tracer_axes(1:2), Time%model_time,                                  &
       'non-Bouss steric sea level tendency from kpp_nloc tendency', 'm/s',    &
       missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_kpp_nloc > 0) compute_watermass_diag=.true.

  id_eta_tend_kpp_nloc_glob= -1          
  id_eta_tend_kpp_nloc_glob= register_diag_field ('ocean_model', 'eta_tend_kpp_nloc_glob',&
       Time%model_time,                                                                   &
       'global mean non-bouss steric sea level tendency from kpp_nloc tendency',          &
       'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_kpp_nloc_glob > 0) compute_watermass_diag=.true.


  id_neut_temp_kpp_nloc = register_diag_field ('ocean_model',                  &
    'neut_temp_kpp_nloc',                                                      &
    Grd%tracer_axes(1:3), Time%model_time,                                     &
    'temp related locally referenced potrho tendency due to KPP nonlocal term',&
    '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e10,1e10/))
  if(id_neut_temp_kpp_nloc > 0) compute_watermass_diag=.true.

  id_wdian_temp_kpp_nloc = register_diag_field ('ocean_model',        &
    'wdian_temp_kpp_nloc', Grd%tracer_axes(1:3), Time%model_time,     &
    'temp related dianeutral mass transport due to KPP nonlocal term',&
    'kg/sec',missing_value=missing_value, range=(/-1e10,1e10/))
  if(id_wdian_temp_kpp_nloc > 0) compute_watermass_diag=.true.

  id_tform_temp_kpp_nloc = register_diag_field ('ocean_model',                           &  
    'tform_temp_kpp_nloc', Grd%tracer_axes(1:3), Time%model_time,                        &
    'temp related watermass transform due to KPP nonlocal on levels (pre-layer binning)',&
    'kg/sec',missing_value=missing_value, range=(/-1e10,1e10/))
  if(id_tform_temp_kpp_nloc > 0) compute_watermass_diag=.true.

  id_neut_temp_kpp_nloc_on_nrho = register_diag_field ('ocean_model',                                &
    'neut_temp_kpp_nloc_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                        &
    'temp related update of locally ref potrho from KPP nolocal as binned to neutral density layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_temp_kpp_nloc_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_temp_kpp_nloc_on_nrho = register_diag_field ('ocean_model',                               &
    'wdian_temp_kpp_nloc_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
    'temp related dianeutral mass transport due to KPP nonlocal as binned to neutral density layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_temp_kpp_nloc_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_temp_kpp_nloc_on_nrho = register_diag_field ('ocean_model',                         &
    'tform_temp_kpp_nloc_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                 &
    'temp related watermass transform due to KPP nonlocal as binned to neutral density layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_temp_kpp_nloc_on_nrho > 0) compute_watermass_diag=.true.



  id_neut_salt_kpp_nloc = register_diag_field ('ocean_model',                  &
    'neut_salt_kpp_nloc',                                                      &
    Grd%tracer_axes(1:3), Time%model_time,                                     &
    'salt related locally referenced potrho tendency due to KPP nonlocal term',&
    '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e10,1e10/))
  if(id_neut_salt_kpp_nloc > 0) compute_watermass_diag=.true.

  id_wdian_salt_kpp_nloc = register_diag_field ('ocean_model',        &
    'wdian_salt_kpp_nloc', Grd%tracer_axes(1:3), Time%model_time,     &
    'salt related dianeutral mass transport due to KPP nonlocal term',&
    'kg/sec',missing_value=missing_value, range=(/-1e10,1e10/))
  if(id_wdian_salt_kpp_nloc > 0) compute_watermass_diag=.true.

  id_tform_salt_kpp_nloc = register_diag_field ('ocean_model',                           &  
    'tform_salt_kpp_nloc', Grd%tracer_axes(1:3), Time%model_time,                        &
    'salt related watermass transform due to KPP nonlocal on levels (pre-layer binning)',&
    'kg/sec',missing_value=missing_value, range=(/-1e10,1e10/))
  if(id_tform_salt_kpp_nloc > 0) compute_watermass_diag=.true.

  id_neut_salt_kpp_nloc_on_nrho = register_diag_field ('ocean_model',                                &
    'neut_salt_kpp_nloc_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                        &
    'salt related update of locally ref potrho from KPP nolocal as binned to neutral density layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_salt_kpp_nloc_on_nrho > 0) compute_watermass_diag=.true.

  id_wdian_salt_kpp_nloc_on_nrho = register_diag_field ('ocean_model',                               &
    'wdian_salt_kpp_nloc_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
    'salt related dianeutral mass transport due to KPP nonlocal as binned to neutral density layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_salt_kpp_nloc_on_nrho > 0) compute_watermass_diag=.true.

  id_tform_salt_kpp_nloc_on_nrho = register_diag_field ('ocean_model',                         &
    'tform_salt_kpp_nloc_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                 &
    'salt related watermass transform due to KPP nonlocal as binned to neutral density layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_salt_kpp_nloc_on_nrho > 0) compute_watermass_diag=.true.


  if(compute_watermass_diag) then 
    write(stdoutunit,'(/a/)') &
    '==>Note: running ocean_vert_kpp_mom4p0_mod w/ compute_watermass_diag=.true.'  
  endif 

end subroutine watermass_diag_init
! </SUBROUTINE> NAME="watermass_diag_init"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag">
!
! <DESCRIPTION>
! Diagnose effects from KPP nonlocal on the watermass transformation.
! </DESCRIPTION>
!
subroutine watermass_diag(Time, T_prog, Dens)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_prog_tracer_type),   intent(in) :: T_prog(:)
  type(ocean_density_type),       intent(in) :: Dens

  integer :: i,j,k,tau
  real, dimension(isd:ied,jsd:jed) :: eta_tend

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_vert_kpp (watermass_diag): module needs initialization ')
  endif 

  if(.not. compute_watermass_diag) return 

  tau = Time%tau


  ! rho diagnostics = sum of temp and salt contributions  
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  wrk5(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Dens%drhodT(i,j,k)*T_prog(index_temp)%wrk1(i,j,k) &
                        +Dens%drhodS(i,j,k)*T_prog(index_salt)%wrk1(i,j,k)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
           wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2)  ! for eta_tend 
        enddo
     enddo
  enddo
  call diagnose_3d(Time, Grd, id_neut_rho_kpp_nloc, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_kpp_nloc, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_kpp_nloc, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_kpp_nloc_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_kpp_nloc_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_kpp_nloc_on_nrho, wrk4)

  if(id_eta_tend_kpp_nloc > 0 .or. id_eta_tend_kpp_nloc_glob > 0) then
      eta_tend(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_kpp_nloc, eta_tend(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_kpp_nloc_glob, eta_tend, cellarea_r)
  endif


  ! temp related contribution  
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Dens%drhodT(i,j,k)*T_prog(index_temp)%wrk1(i,j,k) 
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo
  call diagnose_3d(Time, Grd, id_neut_temp_kpp_nloc, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_kpp_nloc, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_temp_kpp_nloc, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_kpp_nloc_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_kpp_nloc_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_kpp_nloc_on_nrho, wrk4)

  ! salt related contribution  
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Dens%drhodS(i,j,k)*T_prog(index_salt)%wrk1(i,j,k)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo
  call diagnose_3d(Time, Grd, id_neut_salt_kpp_nloc, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_kpp_nloc, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_salt_kpp_nloc, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_kpp_nloc_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_kpp_nloc_on_nrho, wrk3)
  call diagnose_3d_rho(TIme, Dens, id_tform_salt_kpp_nloc_on_nrho, wrk4)

end subroutine watermass_diag
! </SUBROUTINE>  NAME="watermass_diag"


end module ocean_vert_kpp_mom4p0_mod

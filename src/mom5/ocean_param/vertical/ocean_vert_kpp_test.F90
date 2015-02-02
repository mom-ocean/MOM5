module ocean_vert_kpp_test_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> A. Rosati
!</CONTACT>
!
!<CONTACT EMAIL="martin.schmidt@io-warnemuende.de"> Martin Schmidt 
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
!<OVERVIEW>
! Vertical viscosity and diffusivity according KPP.  
!
! This module has extra code options to handle regions of extremely fresh water.  
!
! It also has options for both the Cgrid and Bgrid. 
!
! It is undergoing further development in collaboration with NCAR scientists.  
! So it will undergo significant change during 2012.  
!
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes vertical viscosity and diffusivity according to 
! the K-profile parameterization scheme of Large, McWilliams, and 
! Doney (1994). It computes both local and non-local mixing.
! The code has been updated to MOM4p1, so that vertical grid increments
! are suitable for generalized vertical coordinate models.  When run
! as geopotential model, there will be some differences, since the
! MOM4.0 code (available in ocean_vert_kpp_mom4p0.F90) incorrectly
! ignored the free surface undulations affecting the top model grid
! cell thickness.  
!
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
! <REFERENCE>
! Danabasoglu etal (2006) 
! Diurnal coupling in the tropical oceans of CCSM3
! Journal of Climate (2006) vol 19 pages 2347--2365
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
!<NAMELIST NAME="ocean_vert_kpp_test_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Logical switch to enable kpp diffusion.  Default is false. 
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  Logical switch for debugging. Default debug_this_module=.false. 
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
!
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
!
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
!  Default smooth_blmc=.false.
!
!  Warning: This smoother can cause some problems with ghat in regions
!  of zero surface forcing.  To understand details, one needs 
!  the paper of Large et al. Vertical diffusion has the general form 
!  <wx> = K(x_z - ghats)
!  In the surface layer a vertical scale function ws is estimated.
!  We have K ~ ws and ghats ~1/ws. If wind stress is zero the vertical 
!  scale ws becomes zero too. Hence, ghats is very large 
!  (something finite, since it is divided by ws+epsln). Now it may happen,
!  that the bouyancy flux becomes negative (~ -10-30). This enables
!  the nonlocal scheme. Because the mixing coefficient in the 
!  surface boundary layer scales with ws the corresponding
!  time tendency should be of the order (1/ws * ws = finite). However,
!  if smooth_blmc is enabled, it may happen, that from neighbouring
!  points with different mixing depth a finite value for
!  the mixing coefficient leaks in. In this case
!  the tracer time tendency from the nonlocal scheme becomes huge
!  and the model fails. 
!
!  The smoother destroys the consistency between ghats and diff_cbt.
!  In most cases this should not matter, but the example shows,
!  that sudden model failure is possible under otherwise
!  stable and smooth conditions.
! 
!  </DATA> 
!  <DATA NAME="kl_min" TYPE="integer">
!  Lower loop index for finding new kbl. Needed for use with certain 
!  tests of OBC, where kl_min=1 needed, whereas default in original 
!  implementation has kl_min=2.  Default in MOM is kl_min=2. 
!  </DATA> 
!  <DATA NAME="kbl_standard_method" TYPE="logical">
!  For computing kbl as in the MOM4p0d code, which is taken from 
!  the original NCAR scheme.  If false, then will slightly modify 
!  the logic.  The modified logic has been found necessary when running
!  with as few as two grid cells in the vertical.  
!  Default kbl_standard_method=.true.
!  </DATA> 
!  <DATA NAME="limit_with_hekman" TYPE="logical">
!  Limiting the boundary layer depth with the Ekman depth may result in a 
!  shallow boundary layer. In this case the internal values of the vertical 
!  mixing and viscosity coefficients may be large. This results in 
!  unrealistically large non-local vertical mixing
!  Default limit_with_hekman=.true.
!  </DATA> 
!  <DATA NAME="limit_ghats" TYPE="logical">
!  Limits the non-local vertical tracer flux to the value of the tracer 
!  surface flux.   
!  Default limit_ghats=.false.
!  </DATA> 
!  <DATA NAME="hbl_with_rit" TYPE="logical">
!  The default method for determination of the boundary layer depth may fail
!  if the water column is instable (negative Richardson number) below or above 
!  the layer that contains the diagnosed hbl.  
!  With hbl_with_rit=.true. the search for the boundary layer depth is continued 
!  downward in this case even if the bulk Richardson number exceeds the 
!  critical value. This removes a lot of noise from the boundary layer depth. 
!  Default hbl_with_rit=.false.
!  </DATA> 
!  <DATA NAME="radiation_large" TYPE="logical">
!  Remove the shortwave radiation leaving the boundary layer to the ocean interior 
!  (hence, not absorbed in the boundary layer) from non-local vertical heat flux
!  Default radiation_large=.false.
!  </DATA> 
!  <DATA NAME="radiation_zero" TYPE="logical">
!  Remove the all shortwave radiation from non-local vertical heat flux.
!  Default radiation_zero=.false.
!  </DATA> 
!  <DATA NAME="radiation_iow" TYPE="logical">
!  Keep only the shortwave radiation absorbed between the surface and a certain level
!  in non-local vertical heat flux through this level.
!  Default radiation_iow=.false.
!  </DATA> 
!
!  <DATA NAME="bvf_from_below" TYPE="logical">
!  Use BV-freq. at the cell bottom instead of the cell top
!  as in Danabasoglu et al. (2006).
!  Default bvf_from_below=.false., as this will recover 
!  older behaviour.  
!  </DATA> 
!  <DATA NAME="variable_vtc" TYPE="logical">
!  Make vtc dependent on BV-freq. as in Danabasoglu et al. (2006).
!  Default variable_vtc=.false., as this will recover 
!  older behaviour.  
!  </DATA> 
!  <DATA NAME="use_max_shear" TYPE="logical">
!  Use maximum shear instead of 4-point average 
!  (as in Danabasoglu et al. (2006)).
!  Default use_max_shear=.false., as this will recover 
!  legacy behaviour.  
!  </DATA> 
!
!  <DATA NAME="linear_hbl" TYPE="logical">
!  Use linear interpolation to find the position of hbl.
!  If set to false, then use the quadratic interpolation 
!  as in Danabasoglu et al. (2006). The quadratic approach
!  generally yields a slightly deeper surface boundary layer.
!  Default linear_hbl=.true., as this will recover 
!  older behaviour.  
!  </DATA> 

!  <DATA NAME="wsfc_combine_runoff_calve" TYPE="logical">
!  For computing wsfc as in the MOM4p0d code, where we combine
!  the runoff+calving into a single field called river.  
!  The alternative keeps the fields separate, as would be appropriate
!  for a land model that separately tracks the tracer content in the 
!  calving and runoff. 
!  Default wsfc_combine_runoff_calve=.true., as this will recover
!  the previous behaviour, to the bit. 
!  </DATA> 
!
!  <DATA NAME="smooth_ri_kmax_eq_kmu" TYPE="logical">
!  When smoothing the Richardson number, we do so over a vertical 
!  column with max k-levels set by either kmt or kmu.  The proper 
!  approach is kmu, since we are smoothing riu.  But for backwards
!  compatibility, we default to smooth_ri_kmax_eq_kmu=.false. 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: epsln
use diag_manager_mod, only: register_diag_field, register_static_field
use fms_mod,          only: FATAL, NOTE, WARNING, stdout, stdlog
use fms_mod,          only: write_version_number, open_namelist_file, check_nml_error, close_file
use fms_mod,          only: read_data
use mpp_domains_mod,  only: mpp_update_domains, NUPDATE, EUPDATE
use mpp_mod,          only: input_nml_file, mpp_error

use ocean_density_mod,     only: density, density_delta_z, density_delta_sfc
use ocean_domains_mod,     only: get_local_indices
use ocean_parameters_mod,  only: MOM_BGRID, MOM_CGRID 
use ocean_parameters_mod,  only: grav, cp_ocean, missing_value, von_karman, onefourth
use ocean_parameters_mod,  only: ZSIGMA, PSIGMA, rho0r, rho0, PRESSURE_BASED
use ocean_types_mod,       only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,       only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,       only: ocean_velocity_type, ocean_density_type
use ocean_types_mod,       only: ocean_time_type, ocean_time_steps_type, ocean_thickness_type
use ocean_vert_util_mod,   only: ri_for_bgrid, ri_for_cgrid
use ocean_workspace_mod,   only: wrk1, wrk2, wrk3, wrk4, wrk5
use ocean_workspace_mod,   only: wrk1_v
use ocean_util_mod,        only: diagnose_2d, diagnose_3d, diagnose_3d_u, diagnose_sum
use ocean_tracer_util_mod, only: diagnose_3d_rho


implicit none

private

public ocean_vert_kpp_test_init
public vert_mix_kpp_test

private bldepth
private wscale
private ri_iwmix
private ddmix
private blmix_kpp
private enhance
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
real, dimension(isd:ied,jsd:jed,3)              :: Rib        ! Bulk Richardson number
real, dimension(isd:ied,jsd:jed,3)              :: gat1
real, dimension(isd:ied,jsd:jed,3)              :: dat1

type wsfc_type
  real, dimension(isd:ied,jsd:jed)              :: wsfc        ! rho0r*(stf - pme*(t(i,j,k=1)-tpme) - river*(t(i,j,k=1)-triver))
end type wsfc_type
real, dimension(isd:ied,jsd:jed)                :: sw_frac_hbl ! fractional shortwave penetration at base of mixed layer
integer, dimension(isd:ied,jsd:jed)             :: kbl         ! index of first grid level below hbl

real, private, dimension(isd:ied,jsd:jed,nk)    :: riu         ! Richardson number at base of U-cells
real, private, dimension(isd:ied,jsd:jed,nk)    :: rit         ! Richardson number at base of T-cells
real, private, dimension(isd:ied,jsd:jed,nk)    :: ghats       ! nonlocal transport (s/m^2)
real, private, dimension(isd:ied,jsd:jed)       :: hblt        ! boundary layer depth with tmask 
real, private, dimension(isd:ied,jsd:jed)       :: hbl         ! boundary layer depth

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

type wsfc_type
  real, dimension(:,:), pointer        :: wsfc => NULL()  ! rho0r*(stf - pme*(t(i,j,k=1)-tpme) - river*(t(i,j,k=1)-triver))
end type wsfc_type
real, dimension(:,:), allocatable      :: sw_frac_hbl     ! fractional shortwave penetration at base of mixed layer
integer, dimension(:,:), allocatable   :: kbl             ! index of first grid level below hbl

real, private, dimension(:,:,:), allocatable :: riu         ! Richardson number at base of U-cells
real, private, dimension(:,:,:), allocatable :: rit         ! Richardson number at base of T-cells
real, private, dimension(:,:,:), allocatable :: ghats       ! nonlocal transport (s/m^2)
real, private, dimension(:,:),   allocatable :: hblt        ! boundary layer depth with tmask 
real, private, dimension(:,:),   allocatable :: hbl         ! boundary layer depth

#endif

type(wsfc_type), dimension(:), allocatable :: wsfc
real, parameter :: epsilon = 0.1
real, parameter :: conc1   = 5.0
real, parameter :: zmin    = -4.e-7  ! m3/s3 limit for lookup table of wm and ws
real, parameter :: zmax    = 0.0     ! m3/s3 limit for lookup table of wm and ws
real, parameter :: umin    = 0.0     ! m/s limit for lookup table of wm and ws
real, parameter :: umax    = 0.04    ! m/s limit for lookup table of wm and ws
real, parameter :: Riinfty = 0.8     ! local Richardson number limit for shear instability

real :: Ricr               = 0.3     ! critical bulk Richardson Number
real :: visc_cbu_limit     = 50.0e-4 ! max visc due to shear instability
real :: diff_cbt_limit     = 50.0e-4 ! max diff due to shear instability
real :: visc_con_limit     = 0.1     ! m^2/s. visc due to convective instability
real :: diff_con_limit     = 0.1     ! m^2/s. diff due to convective instability
real :: visc_cbu_iw        = 1.0e-4  ! m^2/s. visc background due to internal waves
real :: diff_cbt_iw        = 0.1e-4  ! m^2/s. diffusivity background due to internal waves
real :: Vtc                          ! non-dimensional coefficient for velocity 
                                     ! scale of turbulant velocity shear        
                                     ! (=function of concv,concs,epsilon,von_karman,Ricr)

real :: cg              ! non-dimensional coefficient for counter-gradient term
real :: deltaz          ! delta zehat in table
real :: deltau          ! delta ustar in table
real :: concv    = 1.8  ! constant for pure convection (eqn. 23)
real :: concv_r         ! inverse concv
real :: vtc_flag = 0.0  ! default to the older approach.   
real :: Lgam = 1.04     ! adjustment to non-gradient flux (McWilliam & Sullivan 2000)
real :: Cw_0 = 0.15     ! eq. (13) in Smyth et al (2002)
real :: l_smyth = 2.0   ! eq. (13) in Smyth et al (2002)
real :: LTmax = 5.0     ! maximum Langmuir turbulence enhancement factor (langmuirfactor) allowed
real :: Wstfac = 0.6    ! stability adjustment coefficient, eq. (13) in Smyth et al (2002)

! for global area normalization
real :: cellarea_r

! for vertical coordinate 
integer :: vert_coordinate
integer :: vert_coordinate_class

real :: rho_cp 
real :: inv_rho_cp

! for Bgrid or Cgrid
integer :: horz_grid
logical :: calc_visc_on_cgrid=.false.  ! internally set: calculate viscosity directly on c-grid

integer :: kbl_max  ! helps to limit vertical loops (defaults to nk)
integer :: kl_min=2 
logical :: kbl_standard_method=.true.

logical :: use_this_module   = .false.   ! logical switch for turning on the kpp scheme 
logical :: shear_instability = .true.    ! logical switch for shear instability mixing
logical :: double_diffusion  = .true.    ! logical switch for double-diffusive mixing
logical :: limit_ghats       = .false.   ! logical switch to limit ghats*diff to 1
logical :: limit_with_hekman = .true.    ! logical switch to limit hbl with hekman
logical :: radiation_large   = .false.   ! logical switch to enable short wave radiation in non-local scheme h_gamma=h
logical :: radiation_zero    = .false.   ! logical switch to enable short wave radiation in non-local scheme h_gamma=0
logical :: radiation_iow     = .false.   ! logical switch to enable short wave radiation in non-local scheme IOW-method
logical :: hbl_with_rit      = .false.   ! logical switch to enable hbl correction for negative rit
logical :: use_sbl_bottom_flux = .false. ! logical switch to add flux through the sbl bottom to the non-local term
logical :: wsfc_combine_runoff_calve = .true.  ! for combining runoff+calving to compute wfsc
logical :: bvf_from_below    = .false.   ! Use BV-freq. at the cell bottom instead of the cell top (Danabasoglu et al.) 
logical :: variable_vtc      = .false.   ! Make vtc dependent on BV-freq.
logical :: use_max_shear     = .false.   ! Use maximum shear instead of 4-point average
logical :: linear_hbl        = .true.    ! To use the linear interpolation as Large etal (1994).
                                         ! Set to .false. for quadratic interpolation as in Danabasoglu et al.
logical :: smooth_ri_kmax_eq_kmu=.false. ! to set details for smoothing the richardson number
real    :: shear_instability_flag=1.0    ! set to 1.0 if shear_instability=.true.

! internally set for computing watermass diagnstics
logical :: compute_watermass_diag = .false. 

! for diagnostics 
integer, dimension(:), allocatable :: id_nonlocal(:)
integer, dimension(:), allocatable :: id_ghats(:)
integer, dimension(:), allocatable :: id_wsfc(:)
integer, dimension(:), allocatable :: id_wbot(:)

logical  :: used
integer  :: id_diff_cbt_kpp_t =-1
integer  :: id_diff_cbt_kpp_s =-1
integer  :: id_visc_cbt_kpp   =-1
integer  :: id_visc_cbu_kpp   =-1
integer  :: id_hblt           =-1
integer  :: id_ws             =-1

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


logical :: non_local_kpp = .true.  ! enable/disable non-local term in KPP
logical :: smooth_blmc   = .false. ! smooth boundary layer diffusitivies to remove grid scale noise
logical :: do_langmuir   = .false. ! whether or not calcualate langmuir turbulence enhancement factor

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
     '$Id: ocean_vert_kpp_test.F90,v 20.0 2013/12/14 00:16:46 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

logical :: module_is_initialized = .FALSE.
logical :: debug_this_module     = .FALSE.

namelist /ocean_vert_kpp_test_nml/ use_this_module, shear_instability, double_diffusion,  &
                                   diff_cbt_iw, visc_cbu_iw,                              &
                                   visc_cbu_limit, diff_cbt_limit,                        &
                                   visc_con_limit, diff_con_limit,                        &
                                   concv, Ricr, non_local_kpp, smooth_blmc,               &
                                   Lgam, Cw_0,l_smyth, LTmax, Wstfac,                     &
                                   kl_min, kbl_standard_method, debug_this_module,        &
                                   limit_with_hekman, limit_ghats, hbl_with_rit,          &
                                   radiation_large, radiation_zero, radiation_iow,        &
                                   use_sbl_bottom_flux, wsfc_combine_runoff_calve,        &
			           bvf_from_below, variable_vtc, use_max_shear,           &  
   			           linear_hbl, smooth_ri_kmax_eq_kmu, do_langmuir
                                 

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_vert_kpp_test_init">
!
! <DESCRIPTION>
! Initialization for the KPP vertical mixing scheme
! </DESCRIPTION>
!
subroutine ocean_vert_kpp_test_init (Grid, Domain, Time, Time_steps, Dens, T_prog, T_diag, &
                                     ver_coordinate, ver_coordinate_class, hor_grid)
  
  type(ocean_grid_type),        intent(in), target  :: Grid
  type(ocean_domain_type),      intent(in), target  :: Domain
  type(ocean_time_type),        intent(in)          :: Time
  type(ocean_time_steps_type),  intent(in)          :: Time_steps
  type(ocean_density_type),     intent(in)          :: Dens
  type(ocean_prog_tracer_type), intent(inout)       :: T_prog(:)
  type(ocean_diag_tracer_type), intent(in)          :: T_diag(:)
  integer,                      intent(in)          :: ver_coordinate
  integer,                      intent(in)          :: ver_coordinate_class
  integer,                      intent(in)          :: hor_grid 
  
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

  integer :: i, j, n, ioun, io_status, ierr
  integer :: rib_dim = 3

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 


  if ( module_is_initialized ) then 
     call mpp_error(FATAL, &
      '==>Error from ocean_vert_kpp_test_mod (ocean_vert_kpp_test_init): module is already initialized')
  endif 

  module_is_initialized = .TRUE.
  vert_coordinate       = ver_coordinate 
  vert_coordinate_class = ver_coordinate_class
  rho_cp     = rho0*cp_ocean 
  inv_rho_cp = 1.0/rho_cp 
  horz_grid  = hor_grid

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_vert_kpp_test_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_vert_kpp_test_nml')
#else
  ioun =  open_namelist_file ()
  read  (ioun, ocean_vert_kpp_test_nml,iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_vert_kpp_test_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_vert_kpp_test_nml)  
  write (stdlogunit, ocean_vert_kpp_test_nml)

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk
#endif 
  kbl_max=nk

  if(use_this_module) then 
    write(stdoutunit,'(/1x,a)') '==> NOTE: USING KPP vertical mixing scheme.'
    write(stdoutunit,'(1x,a)')  '==> NOTE: KPP is typically run with penetrative shortwave heating.'
    write(stdoutunit,'(1x,a/)') '==> NOTE: KPP is typically run with a seasonal and/or diurnal cycle.'
  else 
    write(stdoutunit,'(/1x,a)')'==> NOTE: NOT USING KPP vertical mixing scheme.'
    return
  endif 

  if(debug_this_module) then 
    write(stdoutunit,'(/1x,a)')'==> NOTE: Running KPP with debug_this_module=.true. Set vmixing coeffs to constant.'
  endif 

  if(horz_grid==MOM_CGRID) then 
     calc_visc_on_cgrid=.true.
  else
     calc_visc_on_cgrid=.false.
  endif 

  if(vert_coordinate==ZSIGMA .or. vert_coordinate==PSIGMA) then 
    call mpp_error(WARNING,&
    '==>ocean_vert_kpp_test_mod: ZSIGMA or PSIGMA are not fully working w/ KPP. ')
    call mpp_error(WARNING,&
    'Problems exist /w tracer balances not closing AND instabilities w/ nonlocal transport.')
  endif  

  write(stdoutunit,'(/a,f10.2)') &
  '==>Note from ocean_vert_kpp_test_mod: using forward time step for vert-frict of (secs)', &
                               Time_steps%dtime_u 
  write(stdoutunit,'(/a,f10.2)') &
  '==>Note from ocean_vert_kpp_test_mod: using forward time step for vert-diff  of (secs)', &
                               Time_steps%dtime_t 
  if(kbl_standard_method) then
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: adjust kbl to hbl with the standard method.'
  else
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: adjust kbl to hbl with the non-standard method.'
    if (kl_min > 1) write(stdoutunit,'(1x,a,I4)') '==>WARNING from ocean_vert_kpp_test_mod:'// &
                    ' change kl_min to 1 to avoid negative mixing coefficients with low wind stress.'
  endif
  if(use_sbl_bottom_flux) then
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: '// &
    'Add tracer flux through the SBL-bottom to the nonlocal scheme. This is experimental.'  
  endif                               
  if(bvf_from_below) then
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: '// &
    'Find BVF from a forward derivatives.'  
  else
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: '// &
    'Find BVF from a backward derivatives.'  
  endif                               
  if(variable_vtc) then
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: '// &
    'Diagnose hbl with concv that depends on BVF.'  
  else
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: '// &
    'Diagnose hbl with constant concv.'  
  endif                               
  if(use_max_shear) then          
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: '// &
    'Use maximum shear around a tracer cell for diagnostics of hbl.'  
  else
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: '// &
    'Use average shear around a tracer cell for diagnostics of hbl.'  
  endif                               
  if (linear_hbl) then
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: Use linear interpolation to find hbl'
  else
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: Use quadratic interpolation to find hbl'
  endif
  if(limit_ghats) then
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: Limit ghats*diff_cbt to 1.'
  else
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: Do not limit ghats*diff_cbt to 1.'
  endif
  if(limit_with_hekman) then
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: Limit hbl with hekman for stable case.'
  else
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: Do not limit hbl with hekman for stable case.'
  endif
  if(calc_visc_on_cgrid) then  
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: Calculate viscosity at c-grid.'
  else
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: Interpolate viscosity to c-grid.'
  endif
  if(radiation_large) then
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: sw-radiation in non-local scheme h_gamma=h.'
  elseif (radiation_zero) then
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: sw-radiation in non-local scheme with h_gamma=0.'
  elseif (radiation_iow) then
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: sw-radiation in non-local scheme with IOW-method.'
  else
    write(stdoutunit,'(1x,a)') '==> NOTE from ocean_vert_kpp_test_mod: Leave full sw-radiation in non-local surface flux.'
  endif
  
  if(radiation_large .and. radiation_zero) call mpp_error(FATAL,&
      '==>ocean_vert_kpp_test_mod: Do not enable radiation_large and radiation_zero together. ')   
  if(radiation_large .and. radiation_iow) call mpp_error(FATAL,&
      '==>ocean_vert_kpp_test_mod: Do not enable radiation_large and radiation_iow together. ')   
  if(radiation_zero .and.  radiation_iow) call mpp_error(FATAL,&
      '==>ocean_vert_kpp_test_mod: Do not enable radiation_zero and radiation_iow together. ')   

  Dom => Domain
  Grd => Grid

  tracer_timestep = Time_steps%dtts
  index_temp=-1
  index_salt=-1
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
    call mpp_error(FATAL, &
    '==>Error: ocean_vert_kpp_test_mod: temp and/or salt not present in tracer array')
  endif 


  ! for computing Vtsq 
  if (variable_vtc) then
     vtc_flag=1.0   ! as in Danabasoglu etal (2006) 
  else 
     vtc_flag=0.0   ! as in Large etal (1994)
  endif    

  if (linear_hbl) rib_dim=2 
  cellarea_r = 1.0/(epsln + Grd%tcellsurf)

     
#ifndef MOM_STATIC_ARRAYS

  if(radiation_iow) then
    allocate( T_prog(index_temp)%radiation(isd:ied,jsd:jed,nk)) 
    T_prog(index_temp)%radiation(:,:,:) = 0.0
  endif

  allocate (ghats(isd:ied,jsd:jed,nk))
  allocate (riu(isd:ied,jsd:jed,nk))
  allocate (rit(isd:ied,jsd:jed,nk))
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
  allocate (Rib(isd:ied,jsd:jed,rib_dim))   ! Bulk Richardson number
  allocate (wm(isd:ied,jsd:jed))            ! momentum turbulent velocity scales  (m/s)
  allocate (ws(isd:ied,jsd:jed))            ! scalar turbulent velocity scales  (m/s)
  allocate (Ustk2(isd:ied,jsd:jed))         ! Magnitude of surface stokes drift velocity ^2 (m^2/s^2)
  allocate (gat1(isd:ied,jsd:jed,3))
  allocate (dat1(isd:ied,jsd:jed,3))
  allocate(sw_frac_hbl(isd:ied,jsd:jed))

#endif

  kbl(:,:)         = 0
  ghats(:,:,:)     = 0.0
  riu(:,:,:)       = 0.0
  rit(:,:,:)       = 0.0
  hblt(:,:)        = 0.0
  hbl(:,:)         = 0.0
  sw_frac_hbl(:,:) = 0.0
  Ustk2(:,:)       = 0.0

  do n = 1, num_prog_tracers  
    wsfc(n)%wsfc(:,:) = 0.0
  enddo  

!-----------------------------------------------------------------------
!     initialize some constants for kmix subroutines, and initialize
!     for kmix subroutine "wscale" the 2D-lookup table for wm and ws
!     as functions of ustar and zetahat (=von_karman*sigma*hbl*bfsfc).
!-----------------------------------------------------------------------
 

!-----------------------------------------------------------------------
! define some non-dimensional constants  (recall epsilon=0.1)
!-----------------------------------------------------------------------

      concv_r = 1.0/concv

!     Vtc used in eqn. 23 
      Vtc     = concv * sqrt(0.2/concs/epsilon) / von_karman**2 / Ricr

!     cg = cs in eqn. 20
      cg      = cstar * von_karman * (concs * von_karman * epsilon)**(1./3.)

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
               wmt(i,j) = von_karman*usta/(1.+conc1*zeta)
               wst(i,j) = wmt(i,j)
            else
               if(zeta > zetam) then
                  wmt(i,j) = von_karman* usta * (1.-conc2*zeta)**(1./4.)
               else
                  wmt(i,j) = von_karman* (conam*usta**3-concm*zehat)**(1./3.)
               endif
               if(zeta > zetas) then
                  wst(i,j) = von_karman* usta * (1.-conc3*zeta)**(1./2.)
               else
                  wst(i,j) = von_karman* (conas*usta**3-concs*zehat)**(1./3.)
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
! error checks and diagnostics setup 
!-----------------------------------------------------------------------

  if (Time_steps%aidif /= 1.0) then
    call mpp_error(FATAL,&
    '==>Error in ocean_vert_kpp_test_mod (ocean_vert_kpp_test_init): KPP must use aidif=1 for implicit vertical mixing')
  endif

  ! register diagnostics 

  allocate(id_nonlocal(num_prog_tracers))
  allocate(id_ghats(2))
  allocate(id_wsfc(num_prog_tracers))
  allocate(id_wbot(num_prog_tracers))
  id_nonlocal=-1
  do n = 1, num_prog_tracers
     if(n==index_temp) then
        id_nonlocal(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_nonlocal_KPP', &
                     Grd%tracer_axes(1:3), Time%model_time,                                         &
                     'cp*rho*dzt*nonlocal tendency from KPP', trim(T_prog(n)%flux_units),           &
                     missing_value=missing_value, range=(/-1.e10,1.e10/))
        id_wsfc(n)   = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_wsfc_KPP',       &
                     Grd%tracer_axes(1:2), Time%model_time,                                         &
                     'cp*rho*dzt*surface tendency from KPP', trim(T_prog(n)%flux_units),            &
                     missing_value=missing_value, range=(/-1.e10,1.e10/))
     else
        id_nonlocal(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_nonlocal_KPP', &
                     Grd%tracer_axes(1:3), Time%model_time,                                         &
                     'rho*dzt*nonlocal tendency from KPP', trim(T_prog(n)%flux_units),              &
                     missing_value=missing_value, range=(/-1.e10,1.e10/))
        id_wsfc(n)   =  register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_wsfc_KPP',      &
                     Grd%tracer_axes(1:2), Time%model_time,                                         &
                     'rho*dzt*surface tendency from KPP', trim(T_prog(n)%flux_units),               &
                     missing_value=missing_value, range=(/-1.e10,1.e10/))
     endif
     id_wbot(n)   = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_wbot_KPP', &
          Grd%tracer_axes(1:2), Time%model_time,   &
          'tracer flux through sbl-bottom', trim(T_prog(n)%flux_units),    &
          missing_value=missing_value, range=(/-1.e10,1.e10/))
  enddo
  id_ghats(1) = register_diag_field ('ocean_model', 'temp_ghats_KPP', &
       Grd%tracer_axes(1:3), Time%model_time,      &
       'nonlocal term ghats * diff_cbt from KPP', 'none',     &
       missing_value=missing_value, range=(/-1.e10,1.e10/))
  id_ghats(2) = register_diag_field ('ocean_model', 'salt_ghats_KPP', &
       Grd%tracer_axes(1:3), Time%model_time,              &
       'nonlocal term ghats * diff_cbt from KPP', 'none',     &
       missing_value=missing_value, range=(/-1.e10,1.e10/))

  id_diff_cbt_kpp_t = register_diag_field('ocean_model','diff_cbt_kpp_t',         &
       Grd%tracer_axes(1:3),Time%model_time, 'vert diffusivity from kpp for temp',&
       'm^2/sec', missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_diff_cbt_kpp_s = register_diag_field('ocean_model','diff_cbt_kpp_s',          &
       Grd%tracer_axes(1:3), Time%model_time, 'vert diffusivity from kpp for salt',&
       'm^2/sec', missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_visc_cbt_kpp = register_diag_field('ocean_model','visc_cbt_kpp',            &
       Grd%tracer_axes(1:3),Time%model_time, 'vert viscosity from kpp on T-cell',&
       'm^2/sec', missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_visc_cbu_kpp = register_diag_field('ocean_model','visc_cbu_kpp',            &
       Grd%vel_axes_uv(1:3),Time%model_time, 'vert viscosity from kpp on U-cell',&
       'm^2/sec', missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_hblt = register_diag_field('ocean_model','hblt',Grd%tracer_axes(1:2), &
       Time%model_time, 'T-cell boundary layer depth from KPP', 'm',       &
       missing_value = missing_value, range=(/-1.e5,1.e6/),                &
       standard_name='ocean_mixed_layer_thickness_defined_by_mixing_scheme')

  id_ws = register_diag_field('ocean_model','wscale',Grd%tracer_axes(1:3), &
       Time%model_time, 'wscale from KPP', 'm',                            &
       missing_value = missing_value, range=(/-1.e5,1.e6/))

  call watermass_diag_init(Time,Dens)

end subroutine ocean_vert_kpp_test_init
! </SUBROUTINE> NAME="ocean_vert_kpp_test_init">


!#######################################################################
! <SUBROUTINE NAME="vert_mix_kpp_test">
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
subroutine vert_mix_kpp_test (aidif, Time, Thickness, Velocity, T_prog, T_diag, Dens, &
           swflx, sw_frac_zt, pme, river, visc_cbu, visc_cbt, diff_cbt, hblt_depth, do_wave)

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
  real, dimension(isd:,jsd:,:),    intent(inout) :: visc_cbt
  real, dimension(isd:,jsd:,:,:),  intent(inout) :: diff_cbt
  logical,                         intent(in)    :: do_wave

  real, dimension(isd:ied,jsd:jed,nk) :: dbloc1, dbsfc1
  real, dimension(isd:ied,jsd:jed)    :: frazil
  real                                :: smftu, smftv, active_cells
  real                                :: temporary, temporary1
  integer                             :: i, j, k, n, ki, ni, kp
  integer                             :: tau, taum1

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, '==>Error from ocean_vert_kpp_test: module must be initialized')
  endif 

  if(.not. use_this_module) return 

  tau   = Time%tau
  taum1 = Time%taum1

  visc_cbu(:,:,:) = 0.0
  visc_cbt(:,:,:) = 0.0

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
!     compute gradient Ri and dVsq  
!-----------------------------------------------------------------------

  if(horz_grid == MOM_BGRID) then 

      call ri_for_bgrid(Time, Thickness%dzwt, Dens%drhodT, Dens%drhodS,          &
           T_prog(index_temp)%field(:,:,:,taum1), Dens%rho_salinity(:,:,:,taum1),&
           Velocity%u(:,:,:,1,tau), Velocity%u(:,:,:,2,tau), rit, riu)

      ! compute squared vertical difference of velocity
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               wrk2(i,j,k) = (Velocity%u(i,j,1,1,tau) - Velocity%u(i,j,k,1,tau))**2 &
                           + (Velocity%u(i,j,1,2,tau) - Velocity%u(i,j,k,2,tau))**2
            enddo
         enddo
      enddo
      if (use_max_shear) then  ! as in Danabasoglu et al (2006) 
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   dVsq(i,j,k) = max(wrk2(i,j,k),wrk2(i-1,j,k),wrk2(i,j-1,k),wrk2(i-1,j-1,k))
                enddo
             enddo
          enddo
      else   ! as in legacy GFDL simulations 
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   active_cells = Grd%umask(i,j,k) + Grd%umask(i-1,j,k) + Grd%umask(i,j-1,k) + Grd%umask(i-1,j-1,k) + epsln
                   dVsq(i,j,k)  = (wrk2(i,j,k) + wrk2(i-1,j,k) + wrk2(i,j-1,k) + wrk2(i-1,j-1,k))/active_cells
                enddo
             enddo
          enddo
      endif

  else  ! cgrid  

      ! average horizontal C-grid velocity onto T-cell center 
      wrk1_v(:,:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,1) = Grd%tmask(i,j,k)                                                  &
                    *(Velocity%u(i,j,k,tau,1)*Grd%dtw(i,j)+Velocity%u(i-1,j,k,tau,1)*Grd%dte(i,j)) &
                    /(epsln+Grd%dxt(i,j))
               wrk1_v(i,j,k,2) = Grd%tmask(i,j,k)                                                  &
                    *(Velocity%u(i,j,k,tau,2)*Grd%dts(i,j)+Velocity%u(i,j-1,k,tau,2)*Grd%dtn(i,j)) &
                    /(epsln+Grd%dyt(i,j))
            enddo
         enddo
      enddo
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               dVsq(i,j,k) = (wrk1_v(i,j,1,1) - wrk1_v(i,j,k,1))**2 + (wrk1_v(i,j,1,2) - wrk1_v(i,j,k,2))**2
            enddo
         enddo
      enddo

      call ri_for_cgrid(Time, Thickness%dzwt, Dens%drhodT, Dens%drhodS,          &
           T_prog(index_temp)%field(:,:,:,taum1), Dens%rho_salinity(:,:,:,taum1),&
           wrk1_v(:,:,:,1), wrk1_v(:,:,:,2), rit, riu)


  endif  ! endif for Bgrid/Cgrid 


!-----------------------------------------------------------------------
!     density related quantities
!     --------------------------
!     based on mom equation of state, compute normalized
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
!     saline contraction coefficient without 1/rho factor      (kg/m3/PSU)
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

        rhosfc(:,:)  = Dens%rho(:,:,1,tau)
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

        ! ustar is needed on the "T-grid".  It is assumed that masking of 
        ! smf over land was performed inside of the ocean_sbc module. 
        ! smf has units of N/m^2 and so we need rho0r to get ustar in m/s.   
        ! swflx has units W/m^2 so needs 1/rho_cp to get to C*m/s units.
        ! these are proper units for buoyancy fluxes. 
        !
        ! use smf_bgrid for either MOM_BGRID or MOM_CGRID, since smf_bgrid
        ! is taken straight from the FMS coupler, so is more basic than smf_cgrid. 
        do j=jsc,jec
           do i=isc,iec
               active_cells = Grd%umask(i,j,1)   + Grd%umask(i-1,j,1)   &
                             +Grd%umask(i,j-1,1) + Grd%umask(i-1,j-1,1) + epsln
               smftu = rho0r*(Velocity%smf_bgrid(i,j,1)   + Velocity%smf_bgrid(i-1,j,1)     &
                             +Velocity%smf_bgrid(i,j-1,1) + Velocity%smf_bgrid(i-1,j-1,1))  &
                             /active_cells
               smftv = rho0r*(Velocity%smf_bgrid(i,j,2)   + Velocity%smf_bgrid(i-1,j,2)    &
                             +Velocity%smf_bgrid(i,j-1,2) + Velocity%smf_bgrid(i-1,j-1,2)) &
                             /active_cells
               ustar(i,j) = sqrt( sqrt(smftu**2 + smftv**2) )
             
               Bo(i,j) = grav * (talpha(i,j,1) * &
                      (wsfc(index_temp)%wsfc(i,j)+frazil(i,j)/(rho_cp*tracer_timestep))  &
                      -sbeta (i,j,1) * wsfc(index_salt)%wsfc(i,j)) &     
                       / (epsln + rhosfc(i,j))                     
               Bosol(i,j) = grav * talpha(i,j,1) * swflx(i,j)*inv_rho_cp &
                            / (epsln + rhosfc(i,j))
           enddo
        enddo


      ! compute interior mixing coefficients everywhere, due to constant 
      ! internal wave activity, static instability, and local shear instability.
      call ri_iwmix(visc_cbt, diff_cbt)

      if (double_diffusion) then
        call ddmix (Time, T_prog, Dens, diff_cbt)  
      endif


      ! set seafloor values to zero for blmix 
      do ki = 1,nk-1
         do j = jsc,jec
             do i = isc,iec
                visc_cbt(i,j,ki)   = visc_cbt(i,j,ki)  *Grd%tmask(i,j,min(ki+1,nk))
                diff_cbt(i,j,ki,1) = diff_cbt(i,j,ki,1)*Grd%tmask(i,j,min(ki+1,nk))
                diff_cbt(i,j,ki,2) = diff_cbt(i,j,ki,2)*Grd%tmask(i,j,min(ki+1,nk))
             enddo
         enddo
      enddo


      ! boundary layer mixing coefficients: diagnose new bldepth
      call bldepth(Thickness, sw_frac_zt, do_wave) 
 

      !  boundary layer diffusivities
      call blmix_kpp(Thickness, diff_cbt, visc_cbt, do_wave)
      call diagnose_3d(Time, Grd, id_ws, wrk1(:,:,:))


      ! enhance diffusivity and viscosity at interface kbl-1
      call enhance(Thickness, diff_cbt, visc_cbt) 

!-----------------------------------------------------------------------
!     combine interior and b.l. coefficients and nonlocal term
!-----------------------------------------------------------------------

      ! smooth diffusivities with a 5 point filter
      ! to remove 2 delta x noise
      if (smooth_blmc) then
          
          wrk1(:,:,:) = 0.0
          wrk2(:,:,:) = 0.0
          wrk3(:,:,:) = 0.0

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
                            (4.0*blmc(i,j,k,1)   +&
                                 blmc(i-1,j,k,1) +&
                                 blmc(i+1,j,k,1) +&
                                 blmc(i,j-1,k,1) +&
                                 blmc(i,j+1,k,1)) / active_cells
                       wrk2(i,j,k) = & 
                            (4.0*blmc(i,j,k,2)   +&
                                 blmc(i-1,j,k,2) +&
                                 blmc(i+1,j,k,2) +&
                                 blmc(i,j-1,k,2) +&
                                 blmc(i,j+1,k,2)) / active_cells
                       wrk3(i,j,k) =  &
                            (4.0*blmc(i,j,k,3) +&
                              blmc(i-1,j,k,3)  +&
                              blmc(i+1,j,k,3)  +&
                              blmc(i,j-1,k,3)  +&
                              blmc(i,j+1,k,3)) / active_cells
                   else
                       wrk1(i,j,k) = blmc(i,j,k,1)
                       wrk2(i,j,k) = blmc(i,j,k,2)
                       wrk3(i,j,k) = blmc(i,j,k,3)
                   endif
                enddo
             enddo
          enddo
   
          do k=1,nk-1
             do j=jsc,jec
                do i=isc,iec
                   blmc(i,j,k,1) = wrk1(i,j,k)
                   blmc(i,j,k,2) = wrk2(i,j,k)
                   blmc(i,j,k,3) = wrk3(i,j,k)
                enddo
             enddo
          enddo

      endif
      
      do ki=1,nk-1
         do j=jsc,jec
            do i=isc,iec               
               if (ki < kbl(i,j)) then
                   visc_cbt(i,j,ki)   = blmc(i,j,ki,1)
                   diff_cbt(i,j,ki,1) = blmc(i,j,ki,2)
                   diff_cbt(i,j,ki,2) = blmc(i,j,ki,3)
               else
                   ghats(i,j,ki)=0.0
               endif

               visc_cbt(i,j,ki)   = visc_cbt(i,j,ki)*Grd%tmask(i,j,min(ki+1,nk))               
               diff_cbt(i,j,ki,1) = diff_cbt(i,j,ki,1)*Grd%tmask(i,j,min(ki+1,nk))
               diff_cbt(i,j,ki,2) = diff_cbt(i,j,ki,2)*Grd%tmask(i,j,min(ki+1,nk))               
            enddo
         enddo
      enddo

      ! update halo regions to compute spatial avg of visc_cbt to get visc_cbu
      call mpp_update_domains (visc_cbt(:,:,:), Dom%domain2d, flags=NUPDATE+EUPDATE)
      
      ! compute viscosity on "U-grid" by averaging visc_cbt from "T-grid". 
      ! note: zero out land values
      do ki=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               visc_cbu(i,j,ki) = onefourth*Grd%umask(i,j,min(ki+1,nk))* &
              (visc_cbt(i,j,ki)+visc_cbt(i+1,j,ki)+visc_cbt(i,j+1,ki)+visc_cbt(i+1,j+1,ki))
            enddo
         enddo             
      enddo

      visc_cbt(:,:,nk)   = 0.0
      visc_cbu(:,:,nk)   = 0.0
      diff_cbt(:,:,nk,1) = 0.0
      diff_cbt(:,:,nk,2) = 0.0

      if(debug_this_module) then 
          do k=1,nk-1
             do j=jsd,jed
                do i=isd,ied
                   visc_cbu(i,j,k)   = Grd%umask(i,j,k+1)*visc_cbu_iw
                   visc_cbt(i,j,k)   = Grd%tmask(i,j,k+1)*visc_cbu_iw
                   diff_cbt(i,j,k,1) = Grd%tmask(i,j,k+1)*diff_cbt_iw
                   diff_cbt(i,j,k,2) = Grd%tmask(i,j,k+1)*diff_cbt_iw
                enddo
             enddo
          enddo
      endif

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

       if (non_local_kpp) then

           do n=1,num_prog_tracers

              ! initialize work array for diagnostic storage
              T_prog(n)%wrk1(:,:,:) = 0.0
              wrk1(:,:,:) = 0.0
              wrk2(:,:,:) = 0.0

              ni = 1
              if (n==index_salt) ni = 2

              if(use_sbl_bottom_flux) then
                  do j=jsc,jec
                     do i=isc,iec
                        ki = max(kbl(i,j),1)
                        kp = min(kbl(i,j)+1,nk)
                        temporary = Thickness%dzwt(i,j,ki)
                        temporary1 = min(diff_cbt(i,j,ki,ni), temporary*temporary/tracer_timestep)
                        wrk2(i,j,1) =  Grd%tmask(i,j,kp)*temporary1 &
                             * (T_prog(n)%field(i,j,ki,taum1) - T_prog(n)%field(i,j,kp,taum1)) &
                             /(temporary + epsln)
                        wsfc(n)%wsfc(i,j) = wsfc(n)%wsfc(i,j)  + Grd%tmask(i,j,kp)*temporary1  &
                             * (T_prog(n)%field(i,j,ki,taum1) - T_prog(n)%field(i,j,kp,taum1)) &
                             /(temporary + epsln)
                     enddo
                  enddo
                  call diagnose_2d(TIme, Grd, id_wbot(n), wrk2(:,:,1))
              endif


              ! remove solar radiation absorbed from 0 to level k 
              ! requires absorbed radiation in T_prog(n)%radiation
              if (n==index_temp .and. radiation_iow) then    

                  k = 1
                  do j=jsc,jec
                     do i=isc,iec
                        wrk2(i,j,k) = wsfc(n)%wsfc(i,j) + T_prog(n)%radiation(i,j,k) * inv_rho_cp
                     enddo
                  enddo
                  do k=2,kbl_max
                     do j=jsc,jec
                        do i=isc,iec
                           wrk2(i,j,k) = wrk2(i,j,k-1) + T_prog(n)%radiation(i,j,k) * inv_rho_cp
                        enddo
                     enddo
                  enddo

                  if (limit_ghats) then
                      k=1
                      do j=jsc,jec
                         do i=isc,iec
                            T_prog(n)%wrk1(i,j,k) = rho0 * wrk2(i,j,k)*(-min(diff_cbt(i,j,k,ni)*ghats(i,j,k),1.))
                            T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
                         enddo
                      enddo
                      do k=2,kbl_max
                         do j=jsc,jec
                            do i=isc,iec
                               T_prog(n)%wrk1(i,j,k) =  rho0 * ( wrk2(i,j,k-1)*min(diff_cbt(i,j,k-1,ni)*ghats(i,j,k-1),1.) &
                                                               - wrk2(i,j,k)  *min(diff_cbt(i,j,k,ni)  *ghats(i,j,k)  ,1.))
                               T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
                            enddo
                         enddo
                      enddo
                  else 
                      k=1
                      do j=jsc,jec
                         do i=isc,iec
                            T_prog(n)%wrk1(i,j,k) = rho0 * wrk2(i,j,k)*(-diff_cbt(i,j,k,ni)*ghats(i,j,k))
                            T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
                         enddo
                      enddo
                      do k=2,kbl_max
                         do j=jsc,jec
                            do i=isc,iec
                               T_prog(n)%wrk1(i,j,k) =  rho0 * ( wrk2(i,j,k-1) * diff_cbt(i,j,k-1,ni)*ghats(i,j,k-1) &
                                                               - wrk2(i,j,k)   * diff_cbt(i,j,k,ni)  *ghats(i,j,k) )
                               T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
                            enddo
                         enddo
                      enddo
                  endif  ! endif for limit_ghats 


              else ! alternative treatments of radiation 


                  ! remove solar radiation penetrating the boundary layer
                  if (n==index_temp .and. radiation_large) then   
                      do j=jsc,jec
                         do i=isc,iec
                            wsfc(n)%wsfc(i,j) = wsfc(n)%wsfc(i,j) - swflx(i,j) * inv_rho_cp * sw_frac_hbl(i,j)
                         enddo
                      enddo
                  endif


                  ! remove all solar radiation 
                  if (n==index_temp .and. radiation_zero) then    
                      do j=jsc,jec
                         do i=isc,iec
                            wsfc(n)%wsfc(i,j) = wsfc(n)%wsfc(i,j) - swflx(i,j) * inv_rho_cp 
                         enddo
                      enddo
                  endif

                  if (limit_ghats) then
                      k=1
                      do j=jsc,jec
                         do i=isc,iec
                            T_prog(n)%wrk1(i,j,k) = rho0*wsfc(n)%wsfc(i,j)*(-min(diff_cbt(i,j,k,ni)*ghats(i,j,k),1.))
                            T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
                         enddo
                      enddo
                      do k=2,kbl_max
                         do j=jsc,jec
                            do i=isc,iec
                               T_prog(n)%wrk1(i,j,k) = &
                                    rho0*wsfc(n)%wsfc(i,j)*(min(diff_cbt(i,j,k-1,ni)*ghats(i,j,k-1),1.) - &
                                    min(diff_cbt(i,j,k,ni)  *ghats(i,j,k),1.))
                               T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
                            enddo
                         enddo
                      enddo
                  else
                      k=1
                      do j=jsc,jec
                         do i=isc,iec
                            T_prog(n)%wrk1(i,j,k) = rho0*wsfc(n)%wsfc(i,j)*(-diff_cbt(i,j,k,ni)*ghats(i,j,k))
                            T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
                         enddo
                      enddo
                      do k=2,kbl_max
                         do j=jsc,jec
                            do i=isc,iec
                               T_prog(n)%wrk1(i,j,k) = &
                                    rho0*wsfc(n)%wsfc(i,j)*(diff_cbt(i,j,k-1,ni)*ghats(i,j,k-1) - &
                                    diff_cbt(i,j,k,ni)  *ghats(i,j,k))
                               T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
                            enddo
                         enddo
                      enddo
                  endif ! endif for limit_ghats


              endif   ! endif for treatments of solar radiation 


              ! diagnostics 
              if (n==index_temp .or. n==index_salt) then
                  if(id_ghats(ni) > 0) then 
                      if (limit_ghats) then
                         call diagnose_3d(Time, Grd, id_ghats(ni), min(ghats(:,:,:)*diff_cbt(:,:,:,ni),1.))
                      else
                         call diagnose_3d(Time, Grd, id_ghats(n), ghats(:,:,:)*diff_cbt(:,:,:,ni))
                      endif
                  endif
              endif

              if (id_nonlocal(n) > 0) then 
                 call diagnose_3d(Time, Grd, id_nonlocal(n),T_prog(n)%conversion*T_prog(n)%wrk1(:,:,:))
              endif

           enddo   ! enddo for n-loop 

           call watermass_diag(Time, T_prog, Dens)

       endif   ! endif for non_local_kpp


!-----------------------------------------------------------------------
!     send mixing related fields resulting just from kpp
!-----------------------------------------------------------------------

       call diagnose_3d(Time, Grd, id_diff_cbt_kpp_t, diff_cbt(:,:,:,1))
       call diagnose_3d(Time, Grd, id_diff_cbt_kpp_s, diff_cbt(:,:,:,2))
       call diagnose_3d(Time, Grd, id_visc_cbt_kpp, visc_cbt(:,:,:))
       call diagnose_3d_u(Time, Grd, id_visc_cbu_kpp, visc_cbu(:,:,:))
       call diagnose_2d(Time, Grd, id_hblt, hblt(:,:))


end subroutine vert_mix_kpp_test
! </SUBROUTINE> NAME="vert_mix_kpp_test"


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
!
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
!  input
!      real dbloc     = local delta buoyancy         (m/s^2)     
!      real dbsfc     = delta buoyancy w/ respect to sfc(m/s)^2  
!      real ustar     = surface friction velocity     (m/s)      
!      real Bo        = surface turbulent buoyancy forcing(m^2/s^3)
!      real Bosol     = radiative buoyancy forcing (m^2/s^3)       
!      real f         = Coriolis parameter            (1/s)        
!
!  output
!      real hbl        ! boundary layer depth              (m)      
!      real bfsfc      !Bo+radiation absorbed to d=hbf*hbl(m^2/s^3)      
!      real stable     ! =1 in stable forcing; =0 unstable          
!      real caseA      ! =1 in case A, =0 in case B                 
!      integer kbl     ! index of first grid level below hbl        
!
! </DESCRIPTION>
!
subroutine bldepth(Thickness, sw_frac_zt, do_wave)

  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: sw_frac_zt   !3-D array of shortwave fract
  logical,                      intent(in) :: do_wave

  real     :: Ritop         ! numerator of bulk Richardson Number
  real     :: bvfr, Vtsq
  real     :: hekman, hmonob, hlimit 
  real     :: a_co, b_co, c_co, slope_up, sqrt_arg, z_upper, z_up, zkl
  
  integer  :: ka, ku, kl, ksave, iwet, klp1, klm1, klm2
  integer  :: kupper, kup, kdn
  integer  :: i, j
  integer  :: iwscale_use_hbl_eq_zt

  real, parameter :: c2     = 2.0
  real, parameter :: c4     = 4.0
  real, parameter :: c0     = 0.0
  real, parameter :: cekman = 0.7  ! constant for Ekman depth
  real, parameter :: cmonob = 1.0  ! constant for Monin-Obukhov depth
  real, parameter :: vtc_fac = 200.! factor after Danabasoglu et al. (A.3) 

! find bulk Richardson number at every grid level until > Ricr
!
! note: the reference depth is -epsilon/2.*zt(k), but the reference
!       u,v,t,s values are simply the surface layer values,
!       and not the averaged values from 0 to 2*ref.depth,
!       which is necessary for very fine grids(top layer < 2m thickness)
! note: max values when Ricr never satisfied are
!       kbl(i)=kmt(i,j) and hbl(i)=zt(i,kmt(i,j))

!     Count the number of wet points. 
!     This is the amount of cells hbl is to be found. 
      iwet = 0


!      initialize hbl and kbl to bottomed out values
      do j=jsc,jec
        do i=isc,iec
          Rib(i,j,:) = 0.0
          kbl(i,j)   = MAX(Grd%kmt(i,j),2)
          hbl(i,j)   = Thickness%depth_zt(i,j,kbl(i,j))
          iwet       = iwet + min(Grd%kmt(i,j),1)
        enddo
      enddo


      ! Following Large etal (1994), do linear interpolation to find hbl
      if (linear_hbl) then

      ! indices for array Rib(i,j,k), the bulk Richardson number.     
      ka = 1
      ku = 2

      do kl=2,nk
        klp1 = min(kl+1,nk)
          klm1 = kl-1

        ! compute bfsfc = sw fraction at hbf * zt
        do j=jsc,jec
         do i=isc,iec
            bfsfc(i,j)  = Bo(i,j)  + Bosol(i,j) * (0.0 - sw_frac_zt(i,j,kl))  ! default assumes sw included in stf. 
                                                                              ! otherwise use (1.0-sw_frac)         
            stable(i,j) = 0.5 + SIGN( 0.5, bfsfc(i,j) )
            sigma(i,j)  = stable(i,j) * 1. + (1.-stable(i,j)) * epsilon
          enddo
        enddo

        ! compute velocity scales at sigma, for hbl = zt(kl):
        iwscale_use_hbl_eq_zt=1
        call wscale (iwscale_use_hbl_eq_zt, Thickness%depth_zt(:,:,kl), do_wave)

        do j=jsc,jec
          do i=isc,iec

            if((kbl(i,j) == Grd%kmt(i,j))) then

              ! compute the turbulent shear contribution to Rib
              ! eqn. (23)
              if (bvf_from_below) then
                 bvfr = sqrt(abs(0.5*                            &
                      ( dbloc(i,j,kl  ) / Thickness%dzwt(i,j,kl) +  &
                      dbloc(i,j,klp1) / Thickness%dzwt(i,j,klp1) ) ))
              else
                 bvfr = sqrt(abs(0.5*                              &
                      ( dbloc(i,j,klm1) / Thickness%dzwt(i,j,klm1) +  &
                      dbloc(i,j,kl  ) / Thickness%dzwt(i,j,kl) ) ))
              endif
    
              ! to ensure bitwise compatible with earlier code where default was vtc_flag=0.0 
              Vtsq =   Thickness%depth_zt(i,j,kl) * ws(i,j) * bvfr  &
                       * Vtc *(1.0-vtc_flag)                        &
                     + Thickness%depth_zt(i,j,kl) * ws(i,j) * bvfr  &
                       * max(concv, vtc_flag*(concv +.4 - vtc_fac*bvfr)) * (Vtc*concv_r) * vtc_flag 

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
!           numerator of bulk richardson number on grid levels
!           --------------------------------------------------
!           note: land and ocean bottom values need to be set to zero
!           so that the subroutine "bldepth" works correctly
!-----------------------------------------------------------------------
              Ritop = (Thickness%depth_zt(i,j,kl)-Thickness%depth_zt(i,j,1)) * dbsfc(i,j,kl) * Grd%tmask(i,j,kl)
              Rib(i,j,ku) = Ritop / ( dVsq(i,j,kl) + Vtsq + epsln )

              if(Rib(i,j,ku) > Ricr) then

                  if(((rit(i,j,kl-1).lt.0).or.(rit(i,j,kl).lt.0)).and.hbl_with_rit) then  

                      ! Rib(i,j,ku) is not relevant, because locally unstable
                      Rib(i,j,ku) =  Ricr*0.1

                  else

                      ! linearly interpolate to find hbl where Rib = Ricr
                      hbl(i,j) = Thickness%depth_zt(i,j,kl-1) + Thickness%dzwt(i,j,kl-1) * (Ricr - Rib(i,j,ka))  & 
                           /(Rib(i,j,ku)-Rib(i,j,ka) + epsln)
                      kbl(i,j) = kl
                      iwet     = iwet - 1

                  endif

            endif

            endif  !kbl(i,j) == Grd%kmt(i,j), nothing to do otherwise, hbl is found
          enddo
        enddo

        if (iwet.eq.0) cycle  !all hbl at wet points are defined now -> ready
        
         ksave = ka
         ka    = ku
         ku    = ksave

      enddo


      ! perform quadratic interpolation instead of linear interpolation 
      else


        ! indices for array Rib(i,j,k), the bulk Richardson number.     
        kupper = 1
        kup    = 2
        kdn    = 3

        do kl=2,nk
          klp1 = min(kl+1,nk)
          klm1 = kl-1
          klm2 = max(kl-2,1)

          ! compute bfsfc = sw fraction at hbf * zt
          do j=jsc,jec
            do i=isc,iec
               bfsfc(i,j)  = Bo(i,j)  + Bosol(i,j) * (0.0 - sw_frac_zt(i,j,kl))  ! default assumes sw included in stf. 
                                                                                 ! otherwise use (1.0-sw_frac)         
               stable(i,j) = 0.5 + SIGN( 0.5, bfsfc(i,j) )
               sigma(i,j)  = stable(i,j) * 1. + (1.-stable(i,j)) * epsilon
            enddo
          enddo

          
          !  compute velocity scales at sigma, for hbl = zt(kl):
          iwscale_use_hbl_eq_zt=1
          call wscale (iwscale_use_hbl_eq_zt, Thickness%depth_zt(:,:,kl), do_wave)

          do j=jsc,jec
            do i=isc,iec

              if((kbl(i,j) == Grd%kmt(i,j))) then

                ! compute the turbulent shear contribution to Rib
                ! eqn. (23).
                if (bvf_from_below) then
                   bvfr = sqrt(abs(0.5*                              &
                        ( dbloc(i,j,kl  ) / Thickness%dzwt(i,j,kl) +    &
                        dbloc(i,j,klp1) / Thickness%dzwt(i,j,klp1) ) ))
                else
                   bvfr = sqrt(abs(0.5*                              &
                        ( dbloc(i,j,klm1) / Thickness%dzwt(i,j,klm1) +  &
                        dbloc(i,j,kl  ) / Thickness%dzwt(i,j,kl) ) ))
                endif
    
                ! to ensure bitwise compatible with earlier code where default was vtc_flag=0.0 
                Vtsq =   Thickness%depth_zt(i,j,kl) * ws(i,j) * bvfr  &
                         * Vtc *(1.0-vtc_flag)                        &
                       + Thickness%depth_zt(i,j,kl) * ws(i,j) * bvfr  &
                       * max(concv, vtc_flag*(concv +.4 - vtc_fac*bvfr)) * (Vtc*concv_r) * vtc_flag 

                !-----------------------------------------------------------------------
                ! compute bulk Richardson number at new level
                ! note: Ritop needs to be zero on land and ocean bottom
                ! points so that the following if statement gets triggered
                ! correctly. otherwise, hbl might get set to (big) negative
                ! values, that might exceed the limit for the "exp" function
                ! in "swfrac"
                !
                ! eqn. (21)
                !
                ! numerator of bulk richardson number on grid levels
                ! --------------------------------------------------
                ! note: land and ocean bottom values need to be set to zero
                ! so that the subroutine "bldepth" works correctly
                !-----------------------------------------------------------------------
                Ritop = (Thickness%depth_zt(i,j,kl)-Thickness%depth_zt(i,j,1)) * dbsfc(i,j,kl) * Grd%tmask(i,j,kl)

                Rib(i,j,kdn) = Ritop / ( dVsq(i,j,kl) + Vtsq + epsln )

                if(Rib(i,j,kdn) > Ricr) then

                  if(((rit(i,j,kl-1).lt.0).or.(rit(i,j,kl).lt.0)).and.hbl_with_rit) then  

                      ! Rib(i,j,ku) is not relevant, because locally unstable
                      Rib(i,j,kdn) =  Ricr*0.1			    

                  else

                      ! quadratically interpolate to find hbl where Rib = Ricr, adapted from MODULE: vmix_kpp (NCAR)
                      zkl      =   Thickness%depth_zt(i,j,kl)
                      z_up     = - Thickness%depth_zt(i,j,klm1)
                      z_upper  = - Thickness%depth_zt(i,j,klm2)
                      slope_up =  (Rib(i,j,kupper) - Rib(i,j,kup))/(z_up - z_upper + epsln)
                      a_co = (Rib(i,j,kdn) - Rib(i,j,kup) -         &
                              slope_up*(zkl + z_up) )/(z_up + zkl)**2
                      b_co = slope_up + c2 * a_co * z_up
                      c_co = Rib(i,j,kup) + z_up*(a_co*z_up + slope_up) - Ricr
                      sqrt_arg = b_co**2 - c4*a_co*c_co
                      if ( ( abs(b_co) > epsln .and. abs(a_co)/abs(b_co) <= epsln ) .or. sqrt_arg <= c0 ) then

                         hbl(i,j) = -z_up + (z_up + zkl) *    &
                                   (Ricr - Rib(i,j,kup))/               &
                                   (Rib(i,j,kdn) - Rib(i,j,kup))
                      else
                         hbl(i,j) = (-b_co + sqrt(sqrt_arg)) / (c2*a_co)
                      endif
		      
		      kbl(i,j) = kl
                      iwet     = iwet - 1

                  endif
 
                endif

              endif  !kbl(i,j) == Grd%kmt(i,j), nothing to do otherwise, hbl is found
            enddo
          enddo

          if (iwet.eq.0) cycle  !all hbl at wet points are defined now -> ready
         
          ksave   = kupper
          kupper = kup
          kup    = kdn
          kdn    = ksave

        enddo

      endif  ! endif for linear_hbl

!-----------------------------------------------------------------------
!     find stability and buoyancy forcing for boundary layer
!-----------------------------------------------------------------------

      do j=jsc,jec
        do i = isc,iec

          ! Linear interpolation of sw_frac_zt to depth of hbl
          ! (this is inaccurate, since swflux decays exponentially)
          sw_frac_hbl(i,j)=(sw_frac_zt(i,j,kbl(i,j)) - sw_frac_zt(i,j,kbl(i,j)-1)) / Thickness%dzt(i,j,kbl(i,j)) &
              * (hbl(i,j)-Thickness%depth_zt(i,j,kbl(i,j)))+sw_frac_zt(i,j,kbl(i,j))
          bfsfc(i,j)  = Bo(i,j) + Bosol(i,j) * (0.0 - sw_frac_hbl(i,j))
          stable(i,j) = 0.5 + SIGN( 0.5, bfsfc(i,j) )
          bfsfc(i,j)  = bfsfc(i,j) + stable(i,j) * epsln  ! ensures bfsfc never=0

        enddo
      enddo


!-----------------------------------------------------------------------
!        check hbl limits for hekman or hmonob
!        eqn. (24)
!-----------------------------------------------------------------------
      if(limit_with_hekman) then
      do j=jsc,jec
        do i = isc,iec
          if (bfsfc(i,j) > 0.0) then
             hekman = cekman * ustar(i,j) / (abs(Grd%f(i,j))+epsln)
             hmonob = cmonob * ustar(i,j)*ustar(i,j)*ustar(i,j)     &
                     /von_karman / (bfsfc(i,j)+epsln) 
             hlimit = stable(i,j)     * AMIN1(hekman,hmonob) +      &
                     (stable(i,j)-1.) * (Thickness%depth_zt(i,j,nk))
             hbl(i,j) = AMIN1(hbl(i,j),hlimit)
             hbl(i,j) = AMAX1(hbl(i,j),Thickness%depth_zt(i,j,1))
          endif
          kbl(i,j) =Grd%kmt(i,j)
        enddo
      enddo
      else
        do j=jsc,jec
          do i = isc,iec
            if (bfsfc(i,j) > 0.0) then
              hmonob = cmonob * ustar(i,j)*ustar(i,j)*ustar(i,j)     &
                      /von_karman / (bfsfc(i,j)+epsln) 
              hlimit = stable(i,j)      * hmonob +      &
                       (stable(i,j)-1.) * (Thickness%depth_zt(i,j,nk))
              hbl(i,j) = AMIN1(hbl(i,j),hlimit)
              hbl(i,j) = AMAX1(hbl(i,j),Thickness%depth_zt(i,j,1))
            endif
            kbl(i,j) =Grd%kmt(i,j)
          enddo
        enddo
      endif

!-----------------------------------------------------------------------
!     find new kbl
!-----------------------------------------------------------------------
      if(kbl_standard_method .and. vert_coordinate_class/=PRESSURE_BASED) then 
          do kl=kl_min,nk
             do j=jsc,jec
                do i = isc,iec
                   if((kbl(i,j) == Grd%kmt(i,j)).and.(Thickness%depth_zt(i,j,kl) > hbl(i,j))) then
                       kbl(i,j) = kl
                   endif
                enddo
             enddo
          enddo
      elseif(kbl_standard_method .and. vert_coordinate_class==PRESSURE_BASED) then 
          do kl=kl_min,nk
             do j=jsc,jec
                do i = isc,iec
                   if((kbl(i,j) == Grd%kmt(i,j)).and.(Thickness%depth_zt(i,j,kl) > hbl(i,j))) then
                       kbl(i,j) = min(kl,Grd%kmt(i,j))
                   endif
                enddo
             enddo
          enddo
      else 
          do kl=kl_min,nk
             do j=jsc,jec
                do i = isc,iec
                   if((kbl(i,j) == Grd%kmt(i,j)).and.(Thickness%depth_zt(i,j,kl) >= hbl(i,j))) then
                       kbl(i,j) = min(kl,Grd%kmt(i,j))
                   endif
                enddo
             enddo
          enddo
      endif

      kbl_max = maxval(kbl)

!-----------------------------------------------------------------------
!     find stability and buoyancy forcing for final hbl values
!-----------------------------------------------------------------------

      do j=jsc,jec
        do i = isc,iec
          if (kbl(i,j)==0 .or. kbl(i,j)==1) cycle

          ! Linear interpolation of sw_frac_zt to depth of hbl
          ! (this is inaccurate, since swflux decays exponentially)
          sw_frac_hbl(i,j)=(sw_frac_zt(i,j,kbl(i,j)) - sw_frac_zt(i,j,kbl(i,j)-1)) / Thickness%dzt(i,j,kbl(i,j)) &
              * (hbl(i,j)-Thickness%dzt(i,j,kbl(i,j)))+sw_frac_zt(i,j,kbl(i,j))
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
          if (kbl(i,j)==0) cycle
          caseA(i,j)  = 0.5 + &
          SIGN( 0.5,Thickness%depth_zt(i,j,kbl(i,j)) - 0.5*Thickness%dzt(i,j,kbl(i,j)) - hbl(i,j) )
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
! zetahat (=von_karman*sigma*hbl*bfsfc).
!
! Note: the lookup table is only used for unstable conditions
! (zehat <= 0), in the stable domain wm (=ws) gets computed
! directly.
!
! Note: the loop has been doubled to allow NEC compilers for vectorisation.
! Speed gain was observed at the SX-6.
! Later compiler versions may do better.
!
!
!  input                                                               <BR/>
!      real sigma  = normalized depth (d/hbl)               <BR/>
!      real hbl    = boundary layer depth (m)               <BR/>  
!      real ustar  = surface friction velocity    (m/s)     <BR/>  
!      real bfsfc  = total surface buoyancy flux (m^2/s^3)  <BR/>

!  output                                                               <BR/>  
!      real wm,ws ! turbulent velocity scales at sigma

! local                                                                 <BR/>  
!      real zehat           ! = zeta *  ustar**3
!
! </DESCRIPTION>
!
subroutine wscale(iwscale_use_hbl_eq_zt, zt_kl, do_wave)

  integer,                    intent(in) :: iwscale_use_hbl_eq_zt
  real, dimension(isd:,jsd:), intent(in) :: zt_kl
  logical,                    intent(in) :: do_wave

  real                :: zdiff, udiff, zfrac, ufrac, fzfrac
  real                :: wam, wbm, was, wbs, u3, langmuirfactor, Cw_smyth
  real                :: zehat           ! = zeta *  ustar**3
  integer             :: iz, izp1, ju, jup1
  integer             :: i, j

!-----------------------------------------------------------------------
!     use lookup table for zehat < zmax only; otherwise use
!     stable formulae
!-----------------------------------------------------------------------

      if (iwscale_use_hbl_eq_zt == 1) then

        do j=jsc,jec
          do i=isc,iec
            zehat = von_karman * sigma(i,j) * zt_kl(i,j) * bfsfc(i,j)

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
              wm(i,j) = von_karman * ustar(i,j) * u3 / ( u3 + conc1*zehat + epsln )
              ws(i,j) = wm(i,j)
            endif

          enddo
        enddo

      else

        do j=jsc,jec
          do i=isc,iec
            zehat = von_karman * sigma(i,j) * hbl(i,j)   * bfsfc(i,j)

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
              wm(i,j) = von_karman * ustar(i,j) * u3 / ( u3 + conc1*zehat + epsln )
              ws(i,j) = wm(i,j)
            endif

          enddo
        enddo

      endif

!----------- if do_wave, add Langmuir turbulence enhancement factor

      if (do_wave .and. do_langmuir) then
         do j=jsc,jec
            do i=isc,iec
               Cw_smyth=Cw_0*(ustar(i,j)*ustar(i,j)*ustar(i,j)/(ustar(i,j)*ustar(i,j)*ustar(i,j) &
                    + Wstfac*von_karman*bfsfc(i,j)*hbl(i,j) + epsln))**l_smyth
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
! Diffusion and viscosity coefficients are on bottom
! of T-cells.
!
! </DESCRIPTION>
!
subroutine ri_iwmix(visc_cbt, diff_cbt)

  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbt
  real, dimension(isd:,jsd:,:,:), intent(inout) :: diff_cbt
  
  real    :: Rigg, ratio, frit, fcont
  integer :: i, j, k

  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec
  
!-----------------------------------------------------------------------
!           evaluate function of Ri for shear instability eqn. (28b&c)
!           only do so here using rit.
!-----------------------------------------------------------------------

           Rigg  = AMAX1( rit(i,j,k) , 0.0)
           ratio = AMIN1( Rigg/Riinfty , 1.0 )
           frit  = (1. - ratio*ratio)
           frit  = frit*frit*frit

!-----------------------------------------------------------------------
!          evaluate function of Ri for convection eqn. (28a)
!-----------------------------------------------------------------------

           fcont = 0.5 * ( abs(rit(i,j,k)) - rit(i,j,k) )
           fcont = fcont / (AMAX1(fcont,epsln))
     
!-----------------------------------------------------------------------
!           mixing due to internal wave activity and static instability 
!           eqn. (29).  Note: visc_cbt on T-point.
!-----------------------------------------------------------------------

           visc_cbt(i,j,k)   = visc_cbu_iw + fcont * visc_con_limit   
           diff_cbt(i,j,k,1) = diff_cbt_iw + fcont * diff_con_limit
           diff_cbt(i,j,k,2) = diff_cbt_iw + fcont * diff_con_limit

!-----------------------------------------------------------------------
!          add contribution due to shear instability
!-----------------------------------------------------------------------

           diff_cbt(i,j,k,1) = diff_cbt(i,j,k,1)              &
                + shear_instability_flag*diff_cbt_limit*frit 
           diff_cbt(i,j,k,2) = diff_cbt(i,j,k,2)              &
                + shear_instability_flag*diff_cbt_limit*frit 
           visc_cbt(i,j,k)   = visc_cbt(i,j,k)                &
                + shear_instability_flag*visc_cbu_limit*frit 

        enddo
     enddo
  enddo
 
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
!      real talpha   ! d(rho)/ d(pot.temperature) (kg/m^3/C)  <BR/> 
!      real sbeta    ! d(rho)/ d(salinity)     (kg/m^3/PSU)  
!
!      diff_cbt = diffusion coefficient at bottom of "t" cells (m**2/s)
!
! local
!      real alphaDT  ! alpha * DT  across interfaces   
!      real betaDS   ! beta  * DS  across interfaces   
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

                 diffdd = 1.0-((Rrho-1.0)/(Rrho0-1.0)) 
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
!      real ustar    ! surface friction velocity         (m/s)  
!      real bfsfc    ! surface buoyancy forcing     (m^2/s^3)   
!      real hbl      ! boundary layer depth              (m)    
!      real stable   ! = 1 in stable forcing                    
!      real caseA    ! = 1 in case A                            
!      integer kbl   ! index of first grid level below hbl
!
!     outputs:
!
!      visc_cbt = viscosity coefficient at bottom of "t" cells (m**2/s)   
!      diff_cbt = diffusion coefficient at bottom of "t" cells (m**2/s)   
!
!      real dkm1(,3)    = boundary layer diff_cbt at kbl-1 level  
!      real blmc(,nk,3) = boundary layer mixing coeff.(m**2/s)    
!      real ghats(,nk)  = nonlocal scalar transport               
!
!    local:
!
!      real gat1(,3)                                              
!      real dat1(,3)                                              
!      real sigma()              = normalized depth (d / hbl)     
!      real ws(), wm()  = turbulent velocity scales (m/s) 
!
! </DESCRIPTION>
!
subroutine blmix_kpp(Thickness, diff_cbt, visc_cbt, do_wave)

  type(ocean_thickness_type),     intent(in)    :: Thickness
  real, dimension(isd:,jsd:,:,:), intent(inout) :: diff_cbt
  real, dimension(isd:,jsd:,:) ,  intent(inout) :: visc_cbt
  logical,                        intent(in)    :: do_wave

  real, dimension(isd:ied,jsd:jed) :: zt_kl_dummy

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
        zt_kl_dummy(:,:)      = 0.0

        call wscale (iwscale_use_hbl_eq_zt, zt_kl_dummy, do_wave)

      do j=jsc,jec
        do i = isc,iec 
          if (kbl(i,j) .eq. 0) cycle
          kn    = ifix(caseA(i,j)+epsln) *(kbl(i,j) -1) +   &
                  (1-ifix(caseA(i,j)+epsln)) * kbl(i,j)
          knm1 = max(kn-1,1)
          knp1 = min(kn+1,nk)

!-----------------------------------------------------------------------
!         find the interior viscosities and derivatives at hbl(i) 
!         eqn. (18)
!-----------------------------------------------------------------------

          delhat = 0.5*Thickness%dzt(i,j,kn) + Thickness%depth_zt(i,j,kn) - hbl(i,j)
          R      = 1.0 - delhat / Thickness%dzt(i,j,kn)
          dvdzup = (visc_cbt(i,j,knm1) - visc_cbt(i,j,kn))/Thickness%dzt(i,j,kn)
          dvdzdn = (visc_cbt(i,j,kn) - visc_cbt(i,j,knp1))/Thickness%dzt(i,j,knp1)
          viscp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+         &
                               R  * (dvdzdn + abs(dvdzdn)) )

          dvdzup = (diff_cbt(i,j,knm1,2) - diff_cbt(i,j,kn,2))/Thickness%dzt(i,j,kn)
          dvdzdn = (diff_cbt(i,j,kn,2)   - diff_cbt(i,j,knp1,2))/Thickness%dzt(i,j,knp1)
          difsp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+         &
                               R  * (dvdzdn + abs(dvdzdn)) )

          dvdzup = (diff_cbt(i,j,knm1,1) - diff_cbt(i,j,kn,1))/Thickness%dzt(i,j,kn)
          dvdzdn = (diff_cbt(i,j,kn,1)   - diff_cbt(i,j,knp1,1))/Thickness%dzt(i,j,knp1)
          diftp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+         &
                               R  * (dvdzdn + abs(dvdzdn)) )

          visch  = visc_cbt(i,j,kn)       + viscp * delhat
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

      do ki=1,kbl_max !results are needed only for ki<kbl, hence, limiting this loops saves a lot of wscale calls

!-----------------------------------------------------------------------
!         compute turbulent velocity scales on the interfaces
!-----------------------------------------------------------------------

        do j=jsc,jec
          do i=isc,iec
            sig        = (Thickness%depth_zt(i,j,ki) + 0.5 * Thickness%dzt(i,j,ki)) / (hbl(i,j)+epsln)
            sigma(i,j) = stable(i,j)*sig + (1.-stable(i,j))*AMIN1(sig,epsilon)
          enddo
        enddo

        iwscale_use_hbl_eq_zt = 0
        zt_kl_dummy(:,:)      = 0.0

        call wscale(iwscale_use_hbl_eq_zt, zt_kl_dummy, do_wave)

!-----------------------------------------------------------------------
!         compute the dimensionless shape functions at the interfaces
!         eqn. (11)
!-----------------------------------------------------------------------

        if (id_ws > 0) wrk1(:,:,ki) =  ws(:,:)
          
        do j=jsc,jec
          do i = isc,iec
            if (ki < kbl(i,j)) then
              sig = (Thickness%depth_zt(i,j,ki) + 0.5 * Thickness%dzt(i,j,ki)) / (hbl(i,j)+epsln)
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
          if (kbl(i,j)== 0 .or. kbl(i,j)==1) cycle
          sig        =  Thickness%depth_zt(i,j,kbl(i,j)-1)  / (hbl(i,j)+epsln)
          sigma(i,j) = stable(i,j) * sig + (1.-stable(i,j)) * MIN(sig,epsilon)
        enddo
      enddo

      iwscale_use_hbl_eq_zt = 0
      zt_kl_dummy(:,:)      = 0.0

      call wscale(iwscale_use_hbl_eq_zt, zt_kl_dummy, do_wave)

      do j=jsc,jec
        do i = isc,iec
          if (kbl(i,j)==0 .or. kbl(i,j)==1) cycle
          sig = Thickness%depth_zt(i,j,kbl(i,j)-1) / (hbl(i,j)+epsln)
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
!      integer kbl =  grid above hbl                      
!      real hbl    =  boundary layer depth (m)            
!      real dkm1   =  bl diffusivity at kbl-1 grid level  
!      real caseA  =  1 in caseA, = 0 in case B
!
! input/output
!      real ghats =  nonlocal transport     (s/m**2)    
!      modified ghats at kbl(i)-1 interface        
! output
!      real blmc = enhanced boundary layer mixing coefficient
!
! local
!  real delta =  fraction hbl lies beteen zt neighbors
!
! </DESCRIPTION>
subroutine enhance(Thickness, diff_cbt, visc_cbt)

type(ocean_thickness_type),     intent(in) :: Thickness
real, dimension(isd:,jsd:,:,:), intent(in) :: diff_cbt
real, dimension(isd:,jsd:,:),   intent(in) :: visc_cbt

  real    :: delta, dkmp5, dstar
  integer :: i, j, ki

      do ki=1,nk-1
        do j=jsc,jec
          do i = isc,iec

            if(ki == (kbl(i,j) - 1) ) then

              delta =  (hbl(i,j)-Thickness%depth_zt(i,j,ki)) &
                     / (Thickness%depth_zt(i,j,ki+1)-Thickness%depth_zt(i,j,ki))

              dkmp5 = caseA(i,j) * visc_cbt(i,j,ki)      &
                    + (1.-caseA(i,j)) * blmc(i,j,ki,1)
              dstar = (1.-delta)**2 * dkm1(i,j,1) + delta**2 * dkmp5
              blmc(i,j,ki,1) = (1.-delta) * visc_cbt(i,j,ki)  &
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

  id_wdian_rho_kpp_nloc = register_diag_field ('ocean_model',    &
    'wdian_rho_kpp_nloc', Grd%tracer_axes(1:3), Time%model_time, &
    'dianeutral mass transport due to KPP nonlocal term',        &
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
    '==>Note: running ocean_vert_kpp_mod w/ compute_watermass_diag=.true.'  
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
  call diagnose_3d_rho(Time, Dens, id_tform_salt_kpp_nloc_on_nrho, wrk4)

end subroutine watermass_diag
! </SUBROUTINE>  NAME="watermass_diag"


end module ocean_vert_kpp_test_mod

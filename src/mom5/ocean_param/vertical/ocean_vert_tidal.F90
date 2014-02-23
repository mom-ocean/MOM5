module ocean_vert_tidal_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<REVIEWER EMAIL="hsimmons@iarc.uaf.edu"> Harper Simmons 
!</REVIEWER>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Hyun-Chul Lee 
!</REVIEWER>
!
!<OVERVIEW>
! This module computes a vertical diffusivity and vertical 
! viscosity deduced from barotropic and baroclinic tidal 
! dissipation.  Assume Prandtl number unity. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes a vertical diffusivity and vertical 
! viscosity deduced from barotropic and baroclinic tidal 
! dissipation. For the baroclinic dissipation, we follow
! Simmons etal, and for the barotropic dissipation we follow 
! Lee etal. Assume Prandtl number unity. 
!
! This code is more general than that in the ocean_vert_kpp_mom4p0_mod.
! The KPP_mom4p0 code remains part of MOM for legacy purposes. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Simmons, Jayne, St. Laurent, and Weaver, 2004:
! Tidally driven mixing in a numerical model of the ocean 
! general circulation.  Ocean Modelling, vol. 6,
! pages 245-263.
! </REFERENCE>
!
! <REFERENCE>
! Jayne and St. Laurent, 2001:
! Parameterizing tidal dissipation over rough topography.
! Geophysical Research Letters, vol. 28, pages 811-814.
! </REFERENCE>
!
! <REFERENCE>
! Hyun-Chul Lee, A. Rosati, and M.J. Spelman, 2006: 
! Barotropic tidal mixing effects in a coupled climate model:
! ocean conditions in the northern Atlantic
! Ocean Modelling, vol 11, pages 464--477
! </REFERENCE>
!
! <REFERENCE>
! Osborn, T.R., 1980: Estimates of the local rate of vertical diffusion 
! from dissipation measurements.  JPO, vol. 10, pages 83-89.
! </REFERENCE>
!
! <REFERENCE>
! Munk and Anderson, 1948: Notes on a theory of the thermocline. 
! Journal of Marine Research, vol 3. pages 276-295.
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, Elements of MOM (2012)
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_vert_tidal_nml">
!  <DATA NAME="use_this_module=" TYPE="logical">
!  Must be .true. to use this module. Default is false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging purposes.
!  </DATA> 
!
!  <DATA NAME="use_wave_dissipation=" TYPE="logical">
!  Set to .true. for using the Simmons etal scheme for 
!  obtaining a diffusivity and viscosity based on internal 
!  wave breaking.  This is a general form of the KPP 
!  scheme "int_tidal_mix". 
!  Default use_wave_dissipation=.false.
!  </DATA> 
!  <DATA NAME="use_drag_dissipation=" TYPE="logical">
!  Set to .true. for using the Lee etal scheme for 
!  obtaining a diffusivity and viscosity based on drag
!  of barotropic tides on bottom.  This is a general 
!  form of the KPP scheme "coastal_tidal_mix".
!  Default use_drag_dissipation=.false.
!  </DATA> 
!  <DATA NAME="use_leewave_dissipation=" TYPE="logical">
!  Set to .true. for using a prototype Nikurashin scheme for 
!  obtaining a diffusivity and viscosity based on breaking 
!  leewaves. This scheme is not related to tides, but it 
!  is incorporated to the baroclinic tide parameterization scheme
!  as a prototype.  It will be placed into ts own module when 
!  the parameterization matures. 
!  Default use_leewave_dissipation=.false.
!  </DATA> 
!
!  <DATA NAME="read_leewave_dissipation" TYPE="logical">
!  If .true. then read in leewave dissipation from a file. 
!  Default read_leewave_dissipation=.false.
!  </DATA> 
!
!  <DATA NAME="read_wave_dissipation" TYPE="logical">
!  If .true. then read in wave dissipation computed from 
!  Jayne and St.Laurent (2001) tide model (or another model). 
!  Default read_wave_dissipation=.false.
!  </DATA> 
!  <DATA NAME="fixed_wave_dissipation" TYPE="logical">
!  If .true. then fix the wave dissipation from that 
!  read in by the tide model, such as Jayne and St.Laurent (2001).
!  This power dissipation will be employed
!  for computing wave induced mixing.  
!  Default fixed_wave_dissipation=.false.
!  </DATA> 
!
!  <DATA NAME="read_roughness" TYPE="logical">
!  If .true. then read in bottom roughness amplitude h, 
!  where roughness_length = kappa*h^2, with kappa a 
!  representative roughness wavelength and h a  
!  representative topographic amplitude.  This information is 
!  used for the Simmons etal wave dissipation parameterization.  
!  </DATA> 
!  <DATA NAME="reading_roughness_length" TYPE="logical">
!  If .true., then the field in the roughness file is 
!  roughness_length = kappa*h^2, with kappa a 
!  representative roughness wavelength and h a  
!  representative topographic amplitude.  This information is 
!  used for the Simmons etal wave dissipation parameterization.  
!  Default reading_roughness_length=.false.
!  </DATA> 
!  <DATA NAME="reading_roughness_amp" TYPE="logical">
!  If .true., then the field in the roughness file is 
!  roughness_amp=h, where roughness_length=kappa*h^2.
!  This information is used for the Simmons etal wave
!  dissipation parameterization.  
!  Default reading_roughness_amp=.false.
!  </DATA> 
!  <DATA NAME="default_roughness_length" UNITS="m"  TYPE="real">
!  Default value for kappa*h^2 = roughness length for use 
!  in the absence of a roughness length dataset. MOM default
!  is default_roughness_length=25.0m.
!  </DATA> 
!
!  <DATA NAME="read_tide_speed" TYPE="logical">
!  If .true. then read in tidal speed (m/s) from a tidal model.   
!  This information is used for the computing the energy dissipation 
!  from tides.
!  scheme.  
!  </DATA> 
!  <DATA NAME="tide_speed_data_on_t_grid" TYPE="logical">
!  To set the input tide speed data on T-grid, set to true.
!  Otherwise, set to false.
!  Default tide_speed_data_on_t_grid=.true.
!  </DATA> 
!
!  <DATA NAME="roughness_scale" UNITS="m"  TYPE="real">
!  Scale for the roughness that characterizes the roughness
!  affecting the tidal dissipation process. Used for setting 
!  roughness_length via roughness_length = kappa*h^2, with 
!  kappa = 2pi/(roughness_scale) and h=topography amplitude. 
!  Default roughness_scale=1e4 as in Jayne and St. Laurent (2001)
!  </DATA> 
!  <DATA NAME="default_tide_speed" UNITS="m/s"  TYPE="real">
!  Default value for tidal speed for use in the absence of a 
!  value from a tidal model. 
!  </DATA> 
!  <DATA NAME="speed_min" UNITS="m/s"  TYPE="real">
!  For the drag scheme, we set the diffusivity as well as the 
!  Richardson number to zero if the tide speed is less than 
!  speed_min.  This serves two purposes: 1/ to reduce overflows
!  in some of the diagnostics; 2/ to set the drag induced diffusivity
!  to zero in regions where the tide speed is small. Default
!  speed_min=5e-3m/s. 
!  </DATA> 
!
!  <DATA NAME="shelf_depth_cutoff" UNITS="m" TYPE="real">
!  For use in defining a mask for the Simmons scheme, with depths
!  shallower than shelf_depth_cutoff removed from the scheme. 
!  shelf_depth_cutoff=1000m in Simmons etal.  
!  Default shelf_depth_cutoff=-1000m so there is no cutoff. 
!  </DATA> 
!
!  <DATA NAME="decay_scale" UNITS="m"  TYPE="real">
!  In the Simmons etal vertical profile function, the exponential decay 
!  scale is determined by this parameter.  Default = 500m as in Simmons 
!  etal (2004).  This vertical profile determines how to deposit the 
!  internal wave energy within a vertical column. 
!  </DATA> 
!
!  <DATA NAME="tidal_diss_efficiency" UNITS="dimensionless"  TYPE="real">
!  Fraction of barotropic tidal energy that is dissipated locally, as 
!  opposed to that which propagates away.  Default=1/3 as in 
!  Simmons etal (2004).
!  </DATA> 
!
!  <DATA NAME="mixing_efficiency" UNITS="dimensionless"  TYPE="real">
!  Fraction of energy that is dissipated which is converted into dianeutral 
!  diffusion of tracer.  Default=0.2 based on Osborn (1980).
!  </DATA> 
!  <DATA NAME="mixing_efficiency_n2depend" TYPE="logical">
!  Allow for mixing efficiency to be a function of 
!  N^2/(N^2+Omega^2), which is close to unity except in 
!  regions where N is very small. 
!  Default mixing_efficiency_n2depend=.false.
!  </DATA> 
!
!  <DATA NAME="wave_energy_flux_max" UNITS="W/m2"  TYPE="real">
!  The maximum mechanical energy from internal tides that is 
!  provided for mixing.  Default wave_energy_flux_max=0.1Watt/m^2. 
!  </DATA> 
!
!  <DATA NAME="wave_diffusivity_monotonic" TYPE="logical">
!  Enforce a monotonic decay of the wave dissipation diffusivity,
!  with largest values near bottom and smaller as move to shallower
!  water.  This behaviour is not guaranteed in general, since the 
!  division by the buoyancy frequency can give non-monotone diffusivities.  
!  Default wave_diffusivity_monotonic=.true. 
!  </DATA> 
!
!  <DATA NAME="munk_anderson_p" UNITS="dimensionless" TYPE="real">
!  The p constant in the Munk-Anderson scheme employed by Lee etal. 
!  This parameter is minus the "p_tide" parameter in the KPP schemes. 
!  Default munk_anderson_p=0.25
! </DATA> 
!  <DATA NAME="munk_anderson_sigma" UNITS="dimensionless" TYPE="real">
!  The sigma constant in the Munk-Anderson scheme employed by Lee etal. 
!  This parameter is called "sigma_tide" in the KPP schemes. 
!  Default munk_anderson_sigma=3.0
! </DATA> 
!  <DATA NAME="drag_dissipation_use_cdbot" TYPE="logical">
!  For using the cdbot_array computed from ocean_bbc.F90 module.  
!  Default drag_dissipation_use_cdbot=.false., as this is consistent 
!  with earlier simulations.  
! </DATA> 
!  <DATA NAME="bottom_drag_cd" UNITS="dimensionless" TYPE="real">
!  Bottom drag coefficient from Lee etal. Default bottom_drag_cd=2.4e-3 
! </DATA> 
!
!  <DATA NAME="background_diffusivity" UNITS="m^2/s"  TYPE="real">
!  Background vertical diffusivity not accounted for by the tidal schemes
!  nor any other scheme such as KPP.  Default=0.1e-4. 
!  </DATA> 
!  <DATA NAME="background_viscosity" UNITS="m^2/s"  TYPE="real">
!  Background vertical viscosity not accounted for by the tidal schemes
!  nor any other scheme such as KPP.  Default=0.1e-4. 
!  </DATA> 
!  <DATA NAME="max_wave_diffusivity" UNITS="m^2/s"  TYPE="real">
!  Maximum tracer diffusivity deduced from the wave dissipation
!  scheme from Simmons etal. Default = 5.e-3 m^2/sec.
!  </DATA> 
!
!  <DATA NAME="max_drag_diffusivity" UNITS="m^2/s"  TYPE="real">
!  Maximum tracer diffusivity deduced from the drag dissipation scheme
!  from Lee etal. Default = 5.e-3 m^2/sec.
!  </DATA> 
!  <DATA NAME="drag_dissipation_efold" TYPE="logical">
!  For setting an efolding whereby the drag dissipation diffusivity 
!  exponentially decreases as move upward in the water column.    
!  There are good reasons to set this logical to true, as the scheme 
!  can produce unreasonably large diffusivities far from the bottom, if 
!  there are tides in the deep ocean. 
!  Default drag_dissipation_efold=.true. 
!  </DATA> 
!  <DATA NAME="drag_dissipation_tide_period" UNITS="s"  TYPE="real">
!  Characteristic tide period for use in computing efolding depth for
!  the tide drag scheme.  Default = 12*60*60 = 12hours for semi-diurnal tide.
!  </DATA> 
!  <DATA NAME="drag_mask_deep" TYPE="logical">
!  For masking out the deep ocean regions for the drag dissipation
!  scheme.  This scheme is meant to apply only in shallow shelves,
!  so it is physically relevant to mask it out.  We apply a mask as
!  determined by the ratio of the frictional tide depth scale and the 
!  total ocean depth.    
!  Default drag_mask_deep=.true. 
!  </DATA> 
!  <DATA NAME="drag_mask_deep_ratio" TYPE="real">
!  For determining the drag dissipation mask.
!  The mask = 0 in regions where 
!  tide_depth/total_depth < drag_mask_deep_ratio
!  Default drag_mask_deep_ratio=0.1 
!  </DATA> 
!  <DATA NAME="smooth_ri_drag_cgrid" TYPE="logical">
!  For smoothing the raw C-grid Richardson number computed for
!  the drag scheme on the Cgrid. Default smooth_ri_drag_cgrid=.true.
!  </DATA> 
!
!  <DATA NAME="use_legacy_methods" TYPE="logical">
!  To compute all mixing coefficients using legacy methods. 
!  There are good reasons to prefer the newer approaches, which motivates 
!  setting the default use_legacy_methods=.false.
!  </DATA> 
!  <DATA NAME="drhodz_min" UNITS="kg/m^3"  TYPE="real">
!  Minimum absolute value for the drhodz used to compute N^2 and rhoN2.
!  This value is needed in order to regularize the diffusivity computed
!  from the tide mixing schemes. Default is drhodz_min=1e-10, which 
!  is much smaller than the (N^2)min = 10^-8 sec^-2 used by Simmons 
!  etal. There is some sensitivity to the choice of drhodz_min, with 
!  larger values leading to reduced deep diffusivities, due to the 
!  N^-2 dependence in the diffusivity calculation. 
!  </DATA> 
!  <DATA NAME="smooth_bvfreq_bottom" TYPE="logical">
!  For smoothing the buoyancy frequency at the bottom.
!  Default smooth_bvfreq_bottom=.true.
!  </DATA> 
!  <DATA NAME="vel_micom_smooth" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM Laplacian mixing 
!  coefficient used in the Laplacian smoothing of diffusivities. 
!  Default vel_micom_smooth=0.2.
!  </DATA>
!
!  <DATA NAME="smooth_rho_N2" TYPE="logical">
!  For smoothing the rho_N2 field via a 1-2-1 filter in 
!  vertical.  This is useful to produce smoother diffusivities. 
!  Default is smooth_rho_N2=.true.
!  </DATA> 
!  <DATA NAME="num_121_passes" TYPE="integer">
!  Number of passes of 1-2-1 filter in vertical for 
!  smoothing the rho_N2 field. Default num_121_passes=1.
!  </DATA> 
!
!</NAMELIST>

use constants_mod,     only: pi, epsln
use diag_manager_mod,  only: register_diag_field, register_static_field
use fms_mod,           only: write_version_number, open_namelist_file, close_file, check_nml_error
use fms_mod,           only: stdout, stdlog, read_data, NOTE, FATAL, WARNING
use mpp_domains_mod,   only: mpp_update_domains
use mpp_mod,           only: input_nml_file, mpp_error

use ocean_domains_mod,    only: get_local_indices
use ocean_operators_mod,  only: LAP_T
use ocean_parameters_mod, only: MOM_BGRID, MOM_CGRID 
use ocean_parameters_mod, only: missing_value, onehalf, onefourth
use ocean_parameters_mod, only: von_karman, rho0, rho0r, omega_earth, grav
use ocean_types_mod,      only: ocean_time_type, ocean_domain_type, ocean_grid_type, ocean_options_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_thickness_type, ocean_density_type, ocean_velocity_type 
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk1_2d
use ocean_util_mod,       only: diagnose_2d, diagnose_3d, diagnose_2d_u, diagnose_3d_u

implicit none

private

! for diagnostics 
integer :: id_tide_speed_wave    =-1
integer :: id_tide_speed_drag    =-1
integer :: id_tide_speed_mask    =-1
integer :: id_roughness_length   =-1
integer :: id_roughness_amp      =-1
integer :: id_roughness_klevel   =-1
integer :: id_energy_flux        =-1
integer :: id_power_waves        =-1
integer :: id_leewave_dissipation=-1
integer :: id_power_diss_wave    =-1
integer :: id_power_diss_drag    =-1
integer :: id_power_diss_tides   =-1
integer :: id_power_diss_leewave =-1
integer :: id_rinumber_drag      =-1
integer :: id_drag_dissipation   =-1
integer :: id_bvfreq_bottom      =-1
integer :: id_mix_efficiency     =-1
integer :: id_bvfreq             =-1
integer :: id_diff_cbt_wave      =-1
integer :: id_diff_cbt_drag      =-1
integer :: id_diff_cbt_leewave   =-1
integer :: id_diff_cbt_tides     =-1
integer :: id_visc_cbt_wave      =-1
integer :: id_visc_cbt_drag      =-1
integer :: id_visc_cbt_leewave   =-1
integer :: id_visc_cbt_tides     =-1
integer :: id_visc_cbu_wave      =-1
integer :: id_visc_cbu_leewave   =-1
integer :: id_visc_cbu_drag      =-1
integer :: id_visc_cbu_tides     =-1
integer :: id_drag_diss_efold    =-1
integer :: id_tide_diff_cbt_back =-1
integer :: id_tide_visc_cbu_back =-1
logical :: used

#include <ocean_memory.h>

real, private, dimension(:,:),   allocatable :: smooth_lap          !2D array of micom diffusivities (m^2/sec) for smoothing
real, private, dimension(:,:),   allocatable :: roughness_amp       ! roughness amplitude (m) from topography
real, private, dimension(:,:),   allocatable :: roughness_length    ! roughness length (m) from topography
real, private, dimension(:,:),   allocatable :: wave_dissipation    ! wave dissipation (W/m2) from tide model 
real, private, dimension(:,:),   allocatable :: leewave_dissipation ! leewave dissipation (W/m2) 
real, private, dimension(:,:),   allocatable :: tide_speed_t        ! T-cell speed (m/s) from barotropic tide model 
real, private, dimension(:,:),   allocatable :: tide_speed_u        ! U-cell speed (m/s) from barotropic tide model 
real, private, dimension(:,:),   allocatable :: tide_speed_mask     ! U-cell mask for when tide speed=0
real, private, dimension(:,:),   allocatable :: rescaled_speed_u    ! U-cell speed (m/s) for Lee etal calculation of Ri
real, private, dimension(:,:),   allocatable :: rescaled_speed_t    ! T-cell speed (m/s) for Lee etal calculation of Ri
real, private, dimension(:,:),   allocatable :: efold_depth_r       ! T-cell inverse efold depth (1/m) for Lee etal
real, private, dimension(:,:),   allocatable :: energy_flux         ! energy flux (W/m^2) out of ext-tide to int-tide
real, private, dimension(:,:),   allocatable :: wave_term           ! static term in wave energy flux calculation
real, private, dimension(:,:),   allocatable :: bvfreq_bottom       ! buoyancy frequency (sec^-1) at ocean bottom  

real, private, dimension(:,:,:), allocatable :: mix_efficiency      ! dimensionless mixing efficiency 
real, private, dimension(:,:,:), allocatable :: bvfreq              ! buoyancy frequency (sec^-1) 
real, private, dimension(:,:,:), allocatable :: rho_N2              ! rho*squared buoyancy frequency (kg/m^3)*(sec^-2) 
real, private, dimension(:,:,:), allocatable :: drhodT              ! partial rho / partial temperature 
real, private, dimension(:,:,:), allocatable :: drhodS              ! partial rho / partial salinity 
real, private, dimension(:,:,:), allocatable :: diff_drag           ! diffusivity (m^2/sec) from drag mixing scheme 
real, private, dimension(:,:,:), allocatable :: diff_wave           ! diffusivity (m^2/sec) from wave mixing scheme 
real, private, dimension(:,:,:), allocatable :: diff_leewave        ! diffusivity (m^2/sec) from leewave mixing scheme 
real, private, dimension(:,:),   allocatable :: tmask_deep          ! nonzero for points deeper than shelf_depth_cutoff

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

character(len=128)  :: version='$$'
character (len=128) :: tagname = '$Name: tikal $'

public vert_mix_tidal 
public ocean_vert_tidal_init

private vert_mix_wave
private vert_mix_drag_bgrid 
private vert_mix_drag_cgrid 
private vert_mix_drag_legacy 
private compute_bvfreq

integer :: index_temp=-1
integer :: index_salt=-1 

! for Bgrid or Cgrid
integer :: horz_grid

real :: p5rho0
real :: decay_scale_inv
real :: sqrt_cd
real :: von_karman_inv
real :: dtime 
real :: roughness_kappa
real :: omega_earth2 

logical :: module_is_initialized  = .false.

! nml parameters 

logical :: use_this_module             = .false.
logical :: use_legacy_methods          = .false.  
logical :: debug_this_module           = .false. 
logical :: use_wave_dissipation        = .false.
logical :: use_drag_dissipation        = .false.
logical :: use_leewave_dissipation     = .false.
logical :: read_roughness              = .false.
logical :: read_wave_dissipation       = .false. 
logical :: read_leewave_dissipation    = .false. 
logical :: reading_roughness_amp       = .false.
logical :: reading_roughness_length    = .false.
logical :: read_tide_speed             = .false.
logical :: wave_diffusivity_monotonic  = .true.
logical :: tide_speed_data_on_t_grid   = .true.
logical :: fixed_wave_dissipation      = .false.
logical :: mixing_efficiency_n2depend  = .false. 
logical :: drag_dissipation_efold      = .true.
logical :: smooth_bvfreq_bottom        = .true.  
logical :: drag_mask_deep              = .true.
logical :: smooth_ri_drag_cgrid        = .true. 

logical :: smooth_rho_N2            = .true.   ! for smoothing the rho_N2 field in vertical with 1-2-1
integer :: num_121_passes           = 1        ! number of 1-2-1 passes 

real    :: drag_mask_deep_ratio         = 0.1 
real    :: roughness_scale              = 85e3     ! (metre)
real    :: default_roughness_length     = 25.0     ! (metre)
real    :: default_tide_speed           = .01      ! (m/s)
real    :: shelf_depth_cutoff           = -1000.0  ! (metre)
real    :: decay_scale                  = 500.0    ! (metre)
real    :: tidal_diss_efficiency        = 0.33333  ! (dimensionless) from Simmons etal 
real    :: mixing_efficiency            = 0.2      ! (dimensionless) from Osborne
real    :: munk_anderson_p              = 0.25     ! (dimensionless) from Munk and Anderson 
real    :: munk_anderson_sigma          = 3.0      ! (dimensionless) from Munk and Anderson 
real    :: bottom_drag_cd               = 2.4e-3   ! (dimensionless) bottom drag from Lee etal
real    :: background_diffusivity       = 0.1e-4   ! (m^2/sec)
real    :: background_viscosity         = 0.1e-4   ! (m^2/sec)
real    :: max_wave_diffusivity         = 5.0e-3   ! (m^2/sec)
real    :: max_drag_diffusivity         = 5.0e-3   ! (m^2/sec)
real    :: drhodz_min                   = 1.e-10   ! (kg/m^4) minimum abs(drhodz) used to compute N^2 
real    :: speed_min                    = 5.e-3    ! (m/s) below which set a mask=0 for drag mixing diffusivity 
real    :: wave_energy_flux_max         = 0.1      ! (W/m^2) 
real    :: drag_dissipation_tide_period = 43200.   ! seconds   
real    :: vel_micom_smooth             = 0.2      ! m/sec for smoothing 
logical :: drag_dissipation_use_cdbot   = .false.  ! for using cdbot_array from ocean_bbc 


namelist /ocean_vert_tidal_nml/ use_this_module, use_legacy_methods, debug_this_module,         &
                                use_wave_dissipation, use_drag_dissipation,                     & 
                                read_roughness, read_tide_speed,                                &
                                default_roughness_length, default_tide_speed,                   &
                                shelf_depth_cutoff, decay_scale, roughness_scale,               &
                                tidal_diss_efficiency, mixing_efficiency,                       &
                                mixing_efficiency_n2depend,                                     &
                                munk_anderson_p, munk_anderson_sigma,                           &
                                drag_dissipation_efold, drag_dissipation_tide_period,           &
                                drag_mask_deep, drag_mask_deep_ratio,                           &    
                                bottom_drag_cd, drhodz_min, speed_min,                          &
                                background_diffusivity, background_viscosity,                   &
                                max_wave_diffusivity, max_drag_diffusivity,                     &
                                smooth_bvfreq_bottom, vel_micom_smooth,                         & 
                                smooth_rho_N2, num_121_passes, wave_diffusivity_monotonic,      &
                                tide_speed_data_on_t_grid,                                      &
                                reading_roughness_amp, reading_roughness_length,                &
                                read_wave_dissipation, fixed_wave_dissipation,                  &
                                wave_energy_flux_max,                                           &
                                use_leewave_dissipation, read_leewave_dissipation,              &
                                drag_dissipation_use_cdbot               

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_vert_tidal_init">
!
! <DESCRIPTION>
! Initialization for the ocean_vert_tidal module.
! </DESCRIPTION>
  subroutine ocean_vert_tidal_init(Grid, Domain, Time, T_prog, Velocity, Ocean_options, dtime_t, vert_mix_scheme, hor_grid)
  
    type(ocean_grid_type),        intent(in), target :: Grid
    type(ocean_domain_type),      intent(in), target :: Domain
    type(ocean_time_type),        intent(in)         :: Time
    type(ocean_prog_tracer_type), intent(in)         :: T_prog(:)
    type(ocean_velocity_type),    intent(in)         :: Velocity
    type(ocean_options_type),     intent(inout)      :: Ocean_options 
    real,                         intent(in)         :: dtime_t 
    character(len=10),            intent(in)         :: vert_mix_scheme
    integer,                      intent(in)         :: hor_grid

    real    :: active_cells, temporary 
    integer :: unit, io_status, ierr
    integer :: i,j,n
    integer :: roughness_has_been_read=0

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

    if ( module_is_initialized ) then 
      call mpp_error(FATAL,&
      '==>Error from ocean_vert_tidal_mod (ocean_vert_tidal_init) module already initialized')
    endif 

    module_is_initialized = .TRUE.

    call write_version_number( version, tagname )

#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_vert_tidal_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_vert_tidal_nml')
#else
    unit = open_namelist_file()
    read(unit, ocean_vert_tidal_nml,iostat=io_status)
    ierr = check_nml_error(io_status, 'ocean_vert_tidal_nml')
    call close_file(unit)
#endif
    write (stdoutunit,'(/)')
    write (stdoutunit,ocean_vert_tidal_nml)    
    write (stdlogunit,ocean_vert_tidal_nml)

    Dom => Domain
    Grd => Grid
    dtime = dtime_t 

#ifndef MOM_STATIC_ARRAYS
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grid%nk
#endif

    if(use_this_module) then 
      call mpp_error(NOTE, '==>Note: USING ocean_vert_tidal_mod')
    else 
      call mpp_error(NOTE, '==>Note: NOT using ocean_vert_tidal_mod')
      Ocean_options%tidal_wave_mix = 'Did NOT use tidal wave mixing option for vertical mixing.'
      Ocean_options%tidal_drag_mix = 'Did NOT use tidal drag mixing option for vertical mixing.'
      return 
    endif 

    decay_scale_inv = 1.0/decay_scale
    p5rho0          = 0.5*rho0
    sqrt_cd         = sqrt(bottom_drag_cd)
    von_karman_inv  = 1.0/von_karman
    roughness_kappa = 2.0*pi/roughness_scale   
    omega_earth2    = omega_earth**2
    horz_grid       = hor_grid

    do n=1, size(T_prog(:))
       if (T_prog(n)%name == 'temp') index_temp = n
       if (T_prog(n)%name == 'salt') index_salt = n
    enddo

    ! allocate arrays needed for buoyancy frequency 
    allocate (mix_efficiency(isd:ied,jsd:jed,nk))
    allocate (bvfreq_bottom(isd:ied,jsd:jed))
    allocate (bvfreq(isd:ied,jsd:jed,nk))
    allocate (rho_N2(isd:ied,jsd:jed,nk))
    allocate (drhodT(isd:ied,jsd:jed,nk))
    allocate (drhodS(isd:ied,jsd:jed,nk))
    allocate (diff_drag(isd:ied,jsd:jed,nk))
    allocate (diff_wave(isd:ied,jsd:jed,nk))
    allocate (diff_leewave(isd:ied,jsd:jed,nk))
    allocate (tmask_deep(isd:ied,jsd:jed))
    mix_efficiency  = 0.0
    bvfreq_bottom   = 0.0
    bvfreq          = 0.0
    rho_N2          = 0.0
    drhodT          = 0.0
    drhodS          = 0.0
    diff_drag       = 0.0
    diff_wave       = 0.0
    diff_leewave    = 0.0
    tmask_deep(:,:) = Grd%tmask(:,:,1)

    allocate (smooth_lap(isd:ied,jsd:jed))
    allocate (wave_term(isd:ied,jsd:jed))
    allocate (energy_flux(isd:ied,jsd:jed))
    allocate (wave_dissipation(isd:ied,jsd:jed))
    allocate (leewave_dissipation(isd:ied,jsd:jed))
    allocate (roughness_length(isd:ied,jsd:jed))
    allocate (roughness_amp(isd:ied,jsd:jed))
    allocate (tide_speed_t(isd:ied,jsd:jed))
    allocate (tide_speed_u(isd:ied,jsd:jed))
    allocate (tide_speed_mask(isd:ied,jsd:jed))
    allocate (rescaled_speed_u(isd:ied,jsd:jed))
    allocate (rescaled_speed_t(isd:ied,jsd:jed))
    allocate (efold_depth_r(isd:ied,jsd:jed))
    smooth_lap(:,:)          = Grd%tmask(:,:,1)*vel_micom_smooth*2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(Grd%dxt(:,:)+Grd%dyt(:,:))
    roughness_length(:,:)    = Grd%tmask(:,:,1)*default_roughness_length
    roughness_amp(:,:)       = sqrt(Grd%tmask(:,:,1)/roughness_kappa)
    tide_speed_t(:,:)        = Grd%tmask(:,:,1)*default_tide_speed
    tide_speed_u(:,:)        = Grd%umask(:,:,1)*default_tide_speed
    tide_speed_mask(:,:)     = 0.0
    rescaled_speed_u(:,:)    = 0.0
    rescaled_speed_t(:,:)    = 0.0
    efold_depth_r(:,:)       = 0.0
    wave_term(:,:)           = 0.0
    energy_flux(:,:)         = 0.0
    wave_dissipation(:,:)    = 0.0
    leewave_dissipation(:,:) = 0.0

    if(use_wave_dissipation) then 
        write (stdoutunit,'(a)') &
             'Using Simmons etal scheme to compute dia-surface diffusivity and viscosity based on internal wave breaking.'
        Ocean_options%tidal_wave_mix = 'Used tidal wave mixing option for vertical mixing.'
    else
        write(stdoutunit,'(a)') 'NOT using Simmons etal scheme for dia-surface diffusivity and viscosity.'
        Ocean_options%tidal_wave_mix = 'Did NOT use tidal wave mixing option for vertical mixing.'
    endif

    if(use_drag_dissipation) then 
        write (stdoutunit,'(a)') &
        'Using Lee etal scheme to compute dia-surface diffusivity and viscosity based on barotropic tide drag on bottom.'
        Ocean_options%tidal_drag_mix = 'Used tidal drag mixing option for vertical mixing.'
        if(vert_mix_scheme == 'kpp_mom4p0') then 
            write (stdoutunit,'(a)') &
            '===>WARNING: Using kpp_mom4p0, where Lee etal scheme can be enabled. Be sure not to double count!!!'
        endif
    else
        write(stdoutunit,'(a)') 'NOT using Lee etal scheme from ocean_vert_tidal_mod for dia-surface mixing.'
        Ocean_options%tidal_drag_mix = 'Did NOT use tidal drag mixing option for vertical mixing.'
    endif

    if(use_leewave_dissipation) then 
        write (stdoutunit,'(a)') &
             'Using prototype for Nikurashin scheme to compute dia-surface diff and visc from breaking leewaves.'
        Ocean_options%leewave_mix = 'Using prototype for leewave mixing option for vertical mixing.'
    else
        write(stdoutunit,'(a)') 'NOT using Nikurashin scheme for dia-surface diffusivity and viscosity.'
        Ocean_options%leewave_mix = 'Did NOT use prototype breaking leewave paramaterization of vertical mixing.'
    endif

    if(.not. use_drag_dissipation .and. .not. use_wave_dissipation) then
        call mpp_error(WARNING, &
        '==>ocean_vert_tidal: No dissipation mechanism is set to determine dia-surface mixing.')
    endif

    if(.not. use_wave_dissipation .and. use_leewave_dissipation) then
        call mpp_error(WARNING, &
        '==>ocean_vert_tidal: The prototype leewave mixing scheme must run with use_wave_dissipation=.true.')
    endif

    ! read in topographic amplitude ("h" in "kappa*h^2" from Simmons etal 2004) on T-grid 
    if(read_roughness) then 
        if(reading_roughness_length) then 
           roughness_has_been_read=roughness_has_been_read+1
           call read_data('INPUT/roughness_length.nc','roughness_length', roughness_length, Domain%domain2d)
           write (stdoutunit,*) '==>ocean_vert_tidal_mod: Completed read of topographic roughness length on T-grid.'
           call mpp_update_domains(roughness_length(:,:), Dom%domain2d) 
           roughness_amp(:,:) = sqrt(roughness_length(:,:)/roughness_kappa)
        endif 
        if(reading_roughness_amp) then 
           roughness_has_been_read=roughness_has_been_read+1
           call read_data('INPUT/roughness_amp.nc','roughness_amp', roughness_amp, Domain%domain2d)
           write (stdoutunit,*) '==>ocean_vert_tidal_mod: Completed read of topographic roughness amplitude on T-grid.'
           call mpp_update_domains(roughness_amp(:,:), Dom%domain2d) 
           roughness_length(:,:) = roughness_kappa*roughness_amp(:,:)*roughness_amp(:,:)
        endif 
        if(roughness_has_been_read > 1) then 
           call mpp_error(FATAL, &
             '==>ocean_vert_tidal_mod: Read in both roughness_amp & roughness_length. Check to be sure what you wish.')
        endif 
        if(roughness_has_been_read==0) then 
           call mpp_error(FATAL, &
             '==>ocean_vert_tidal_mod: To read roughness, reading_roughness_amp or reading_roughness_length must be true.')
        endif 
    else
        write(stdoutunit,'(a)') &
             '==>Note: NOT reading topographic roughness_length for ocean_vert_tidal_mod.'
    endif

    if(read_wave_dissipation) then 
           call read_data('INPUT/wave_dissipation.nc','wave_dissipation', wave_dissipation, Domain%domain2d)
           write (stdoutunit,*) '==>ocean_vert_tidal_mod: Completed read of wave dissipation (W/m^2) on T-grid.'
           call mpp_update_domains(wave_dissipation(:,:), Dom%domain2d) 
    else 
        write(stdoutunit,'(a)') &
             '==>Note: NOT reading wave dissipation for ocean_vert_tidal_mod.'
    endif 

    if(read_leewave_dissipation) then 
        call read_data('INPUT/leewave_dissipation.nc','leewave_dissipation', leewave_dissipation, Domain%domain2d)
        write (stdoutunit,*) '==>ocean_vert_tidal_mod: Completed read of leewave dissipation (W/m^2) on T-grid.'
        call mpp_update_domains(leewave_dissipation(:,:), Dom%domain2d) 
    else 
        write(stdoutunit,'(a)') &
             '==>Note: NOT reading leewave dissipation for ocean_vert_tidal_mod.'
    endif

    ! read tidal speed (m/s) from a tide model, such as the
    ! Global Inverse Solution TPX06.0 created by OSU.
    if(read_tide_speed) then 
        if(tide_speed_data_on_t_grid) then 
           call read_data('INPUT/tideamp.nc','tideamp', tide_speed_t, Domain%domain2d)
           write (stdoutunit,*) '==>ocean_vert_tidal_mod: Completed read of tide_speed on T-grid.'
           call mpp_update_domains(tide_speed_t(:,:), Dom%domain2d) 
        else 
           call read_data('INPUT/tideamp.nc','tideamp', tide_speed_u, Domain%domain2d)
           write (stdoutunit,*) '==>ocean_vert_tidal_mod: Completed read of tide_speed on U-grid.'
           call mpp_update_domains(tide_speed_u(:,:), Dom%domain2d) 
        endif 
    else
        write(stdoutunit,'(a)') &
             '==>Note: NOT reading tide_speed for ocean_vert_tidal_mod.'
        call mpp_error(NOTE, &
             '==>ocean_vert_tidal_mod: Setting tide_speed to default value.')
    endif

    ! map tide speed onto U-cell by 4-point average 
    if(tide_speed_data_on_t_grid) then 
        do j=jsc,jec
           do i=isc,iec
              tide_speed_u(i,j) = onefourth*Grd%umask(i,j,1)     &
                   *(tide_speed_t(i,j)   + tide_speed_t(i+1,j)   &
                    +tide_speed_t(i,j+1) + tide_speed_t(i+1,j+1)) 
           enddo
        enddo
        call mpp_update_domains(tide_speed_u(:,:), Dom%domain2d) 
    endif

    ! speed scale for tides rubbing against bottom 
    ! (defined following eq. (3) in Lee etal)
    if(drag_dissipation_use_cdbot) then 
        write(stdoutunit,'(a)') &
        '==>Note from ocean_vert_tidal: using cdbot_array(i,j) for tide drag_dissipation scheme.'
       do j=jsd,jed
          do i=isd,ied
             rescaled_speed_u(i,j) = sqrt(Velocity%cdbot_array(i,j))*von_karman_inv*tide_speed_u(i,j)
             rescaled_speed_t(i,j) = sqrt(Velocity%cdbot_array(i,j))*von_karman_inv*tide_speed_t(i,j)
             tide_speed_mask(i,j)  = 0.0 
             if(rescaled_speed_u(i,j) > speed_min) then 
                tide_speed_mask(i,j) = 1.0 
             endif 
          enddo
       enddo
    else 
        write(stdoutunit,'(a)') &
        '==>Note from ocean_vert_tidal: using constant bottom drag coefficient for tide drag_dissipation scheme.'
       do j=jsd,jed
          do i=isd,ied
             rescaled_speed_u(i,j) = sqrt_cd*von_karman_inv*tide_speed_u(i,j)
             rescaled_speed_t(i,j) = sqrt_cd*von_karman_inv*tide_speed_t(i,j)
             tide_speed_mask(i,j)  = 0.0 
             if(rescaled_speed_u(i,j) > speed_min) then 
                tide_speed_mask(i,j) = 1.0 
             endif 
          enddo
       enddo
    endif


    ! compute efolding depth scale for use in Lee etal scheme. 
    ! efold depth set as rescaled_speed/(radial tide frequency). 
    ! Choose default radial tide frequency as 2pi/12hrs for semi-diurnal tide.  
    ! let this efolding hold whether using Bgrid or Cgrid.  
    wrk1_2d(:,:) = 0.0
    do j=jsc,jec
       do i=isc,iec
          active_cells =   Grd%umask(i,j,1)   + Grd%umask(i-1,j,1) &
                         + Grd%umask(i,j-1,1) + Grd%umask(i-1,j-1,1) + epsln
          temporary    = (rescaled_speed_u(i,j)   + rescaled_speed_u(i-1,j) + &
                          rescaled_speed_u(i,j-1) + rescaled_speed_u(i-1,j-1))/active_cells
          wrk1_2d(i,j)       = Grd%tmask(i,j,1)*temporary*drag_dissipation_tide_period/(2.0*pi)
          efold_depth_r(i,j) = Grd%tmask(i,j,1)/(wrk1_2d(i,j) + epsln) 
       enddo
    enddo

    ! mask out the deep ocean regions for the drag scheme.  
    if(drag_mask_deep .and. .not. use_legacy_methods) then
        do j=jsc,jec
           do i=isc,iec
              if(Grd%tmask(i,j,1) == 1.0) then 
                  temporary = wrk1_2d(i,j)/(epsln+Grd%ht(i,j))
                  if(temporary > drag_mask_deep_ratio) then 
                      tide_speed_mask(i,j) = 1.0 
                  else 
                      tide_speed_mask(i,j) = 0.0
                  endif
              endif
           enddo
        enddo
        call mpp_update_domains(tide_speed_mask(:,:), Dom%domain2d) 
    endif

    ! compute static piece of the energy flux on T-grid for wave diffusivity 
    do j=jsd,jed
       do i=isd,ied
          wave_term(i,j) = Grd%tmask(i,j,1)*p5rho0*roughness_length(i,j)*tide_speed_t(i,j)**2
       enddo
    enddo

    ! diagnostics 

    id_tide_diff_cbt_back = -1
    id_tide_diff_cbt_back = register_static_field ('ocean_model', 'tide_diff_cbt_back',     &
                     Grid%tracer_axes(1:3), 'static background diff_cbt set in tide module',&
                     'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/))
    if (id_tide_diff_cbt_back > 0) then 
        wrk1(:,:,:) = background_diffusivity*Grid%tmask(:,:,:)
        call diagnose_3d(Time, Grd, id_tide_diff_cbt_back, wrk1(:,:,:))
    endif

    id_tide_visc_cbu_back = -1
    id_tide_visc_cbu_back = register_static_field ('ocean_model', 'tide_visc_cbu_back',     &
                     Grid%vel_axes_wu(1:3), 'static background visc_cbu set in tide module',&
                     'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/))
    if (id_tide_visc_cbu_back > 0) then 
        wrk1(:,:,:) = background_viscosity*Grid%umask(:,:,:)
        call diagnose_3d_u(Time, Grd, id_tide_visc_cbu_back, wrk1(:,:,:))
    endif

    ! e-folding depth for drag dissipation scheme 
    id_drag_diss_efold = register_static_field ('ocean_model', 'drag_diss_efold',&
         Grid%tracer_axes(1:2), 'e-folding depth for drag dissipation scheme',   &
         'm', missing_value=missing_value, range=(/-1.0,1e10/))
    call diagnose_2d(Time, Grd, id_drag_diss_efold, wrk1_2d(:,:))

    ! static input of leewave breaking from Nikurashin
    id_leewave_dissipation = register_static_field ('ocean_model', 'leewave_dissipation',&
         Grid%tracer_axes(1:2), 'specified energy flux input from breaking leewaves',    &
         'W/m^2', missing_value=missing_value, range=(/-1e1,1e15/))
    call diagnose_2d(Time, Grd, id_leewave_dissipation, leewave_dissipation(:,:))

    ! tide speed for breaking internal wave dissipation scheme
    id_tide_speed_wave = register_static_field ('ocean_model', 'tide_speed_wave',                     &
         Grid%tracer_axes(1:2), 'tide speed from tide model for breaking internal wave mixing scheme',&
         'm/s', missing_value=missing_value, range=(/-1e1,1e9/))
    call diagnose_2d(Time, Grd, id_tide_speed_wave, tide_speed_t(:,:))

    ! tide speed scale for barotropic bottom drag dissipation scheme 
    id_tide_speed_drag = register_static_field ('ocean_model', 'tide_speed_drag',              &
         Grid%vel_axes_uv(1:2), 'tide speed from tide model for barotropic drag mixing scheme',&  
         'm/s', missing_value=missing_value, range=(/-1e1,1e9/))
    call diagnose_2d_u(Time, Grd, id_tide_speed_drag, rescaled_speed_u(:,:))

    ! static tide speed mask
    id_tide_speed_mask = register_static_field ('ocean_model', 'tide_speed_mask',           &
         Grid%vel_axes_uv(1:2), 'mask based on tide_speed_drag for barotropic drag mixing', &
         'dimensionless', missing_value=missing_value, range=(/-1e1,1e1/))
    call diagnose_2d_u(Time, Grd, id_tide_speed_mask, tide_speed_mask(:,:))

    ! static roughness amplitude 
    id_roughness_amp = register_static_field ('ocean_model', 'roughness_amp',                  &
         Grid%tracer_axes(1:2), 'roughness amplitude for breaking internal wave mixing scheme',&
         'metre', missing_value=missing_value, range=(/-1e1,1e9/))
    call diagnose_2d(Time, Grd, id_roughness_amp, roughness_amp(:,:))

    ! static roughness length 
    id_roughness_length = register_static_field ('ocean_model', 'roughness_length',         &
         Grid%tracer_axes(1:2), 'roughness length for breaking internal wave mixing scheme',&
         'metre', missing_value=missing_value, range=(/-1e1,1e9/))
    call diagnose_2d(TIme, Grd, id_roughness_length, roughness_length(:,:))


    id_roughness_klevel = register_diag_field ('ocean_model', 'roughness_klevel',                    &
         Grid%tracer_axes(1:2), Time%model_time,                                                     &
         'klevel at top of the bottom layer defined by roughness amplitude for internal tide mixing',&
         'dimensionless', missing_value=missing_value, range=(/-1.0,1.e10/))
    id_energy_flux = register_diag_field ('ocean_model', 'energy_flux',        &
         Grid%tracer_axes(1:2), Time%model_time,                               &
         'energy flux out of barotropic tides for use w/ internal tide mixing',&
         'W/m^2', missing_value=missing_value, range=(/-1e9,1e9/))
    id_power_waves = register_diag_field ('ocean_model', 'power_waves',                           &
         Grid%tracer_axes(1:2), Time%model_time, 'power from barotropic tides to internal tides', &
         'Watt', missing_value=missing_value, range=(/-1e15,1e15/))

    id_power_diss_leewave = register_diag_field ('ocean_model', 'power_diss_leewave',                     &
         Grid%tracer_axes(1:3), Time%model_time, 'power dissipation from mixing due to breaking leewaves',&
         'W/m^2', missing_value=missing_value, range=(/-1e15,1e15/))
    id_power_diss_wave = register_diag_field ('ocean_model', 'power_diss_wave',                         &
         Grid%tracer_axes(1:3), Time%model_time, 'power dissipation from internal wave induced mixing', &
         'W/m^2', missing_value=missing_value, range=(/-1e15,1e15/))
    id_power_diss_drag = register_diag_field ('ocean_model', 'power_diss_drag',                 &
         Grid%tracer_axes(1:3), Time%model_time, 'power dissipation from barotropic tide drag', &
         'W/m^2', missing_value=missing_value, range=(/-1e15,1e15/))
    id_power_diss_tides = register_diag_field ('ocean_model', 'power_diss_tides',& 
         Grid%tracer_axes(1:3), Time%model_time,                                 &
         'power dissipation from barotropic tide drag and baroclinic wave drag', &
         'W/m^2', missing_value=missing_value, range=(/-1e15,1e15/),             &
          standard_name='tendency_of_ocean_potential_energy_content_due_to_tides')

    id_mix_efficiency = register_diag_field ('ocean_model', 'mix_efficiency',                                  &
         Grid%tracer_axes(1:3), Time%model_time, 'efficiency of internal wave dissipation going to mix tracer',&
         'dimensionless', missing_value=missing_value, range=(/-1e5,1e5/))
    id_bvfreq_bottom = register_diag_field ('ocean_model', 'bvfreq_bottom',                     &
         Grid%tracer_axes(1:2), Time%model_time, 'absolute Brunt-Vaisala freq at ocean bottom', &
         's^-1', missing_value=missing_value, range=(/-1e1,1e9/))
    id_bvfreq = register_diag_field ('ocean_model', 'bvfreq',                                         &
         Grid%tracer_axes(1:3), Time%model_time, 'absolute Brunt-Vaisala freq at tracer cell bottom', &
         's^-1', missing_value=missing_value, range=(/-1e1,1e9/))
    id_rinumber_drag = register_diag_field ('ocean_model', 'rinumber_drag',         &
         Grid%tracer_axes(1:3), Time%model_time, 'Richardson number from Lee etal', &
         'dimensionless', missing_value=missing_value, range=(/-1e1,1e18/))
    id_diff_cbt_wave = register_diag_field ('ocean_model', 'diff_cbt_wave',                             &
         Grid%tracer_axes(1:3), Time%model_time, 'diffusivity from breaking internal wave dissipation', &
         'm^2/sec', missing_value=missing_value, range=(/-1.0,1e6/))
    id_diff_cbt_leewave = register_diag_field ('ocean_model', 'diff_cbt_leewave',      &
         Grid%tracer_axes(1:3), Time%model_time, 'diffusivity from breaking leewaves', &
         'm^2/sec', missing_value=missing_value, range=(/-1.0,1e6/))
    id_diff_cbt_drag = register_diag_field ('ocean_model', 'diff_cbt_drag',                             &
         Grid%tracer_axes(1:3), Time%model_time, 'diffusivity from drag of barotropic tides on bottom', &
         'm^2/sec', missing_value=missing_value, range=(/-1.0,1e6/))

    id_visc_cbt_wave = register_diag_field ('ocean_model', 'visc_cbt_wave',                           &
         Grid%tracer_axes(1:3), Time%model_time, 'viscosity from breaking internal wave dissipation', &
         'm^2/sec', missing_value=missing_value, range=(/-1.0,1e6/))
    id_visc_cbt_leewave = register_diag_field ('ocean_model', 'visc_cbt_leewave',    &
         Grid%tracer_axes(1:3), Time%model_time, 'viscosity from breaking leewaves', &
         'm^2/sec', missing_value=missing_value, range=(/-1.0,1e6/))
    id_visc_cbt_drag = register_diag_field ('ocean_model', 'visc_cbt_drag',                           &
         Grid%tracer_axes(1:3), Time%model_time, 'viscosity from drag of barotropic tides on bottom', &
         'm^2/sec', missing_value=missing_value, range=(/-1.0,1e6/))

    id_visc_cbu_wave = register_diag_field ('ocean_model', 'visc_cbu_wave',                           &
         Grid%vel_axes_uv(1:3), Time%model_time, 'viscosity from breaking internal wave dissipation', &
         'm^2/sec', missing_value=missing_value, range=(/-1.0,1e6/))
    id_visc_cbu_leewave = register_diag_field ('ocean_model', 'visc_cbu_leewave',     &
         Grid%vel_axes_uv(1:3), Time%model_time, 'viscosity from breaking leewaves ', &
         'm^2/sec', missing_value=missing_value, range=(/-1.0,1e6/))
    id_visc_cbu_drag = register_diag_field ('ocean_model', 'visc_cbu_drag',                           &
         Grid%vel_axes_uv(1:3), Time%model_time, 'viscosity from drag of barotropic tides on bottom', &
         'm^2/sec', missing_value=missing_value, range=(/-1.0,1e6/))
    id_diff_cbt_tides = register_diag_field ('ocean_model', 'diff_cbt_tides',                                      &
         Grid%tracer_axes(1:3), Time%model_time, 'diffusivity from drag of barotropic tides on bottom + wave drag',&
         'm^2/sec', missing_value=missing_value, range=(/-1.0,1e6/),                                               &
         standard_name='ocean_vertical_tracer_diffusivity_due_to_tides')
    id_visc_cbt_tides = register_diag_field ('ocean_model', 'visc_cbt_tides',                                    &
         Grid%tracer_axes(1:3), Time%model_time, 'viscosity from drag of barotropic tides on bottom + wave drag',&
         'm^2/sec', missing_value=missing_value, range=(/-1.0,1e6/),                                             &
         standard_name='ocean_vertical_tracer_diffusivity_due_to_tides')
    id_visc_cbu_tides = register_diag_field ('ocean_model', 'visc_cbu_tides',                                    &
         Grid%vel_axes_uv(1:3), Time%model_time, 'viscosity from drag of barotropic tides on bottom + wave drag',&
         'm^2/sec', missing_value=missing_value, range=(/-1.0,1e6/),                                             &
         standard_name='ocean_vertical_momentum_diffusivity_due_to_tides')


end subroutine ocean_vert_tidal_init
! </SUBROUTINE> NAME="ocean_vert_tidal_init"


!#######################################################################
! <SUBROUTINE NAME="vert_mix_tidal">
!
! <DESCRIPTION>
! This subroutine computes vertical tracer diffusivity and viscosity
! based on one or both of the following dissipation mechanisms:
!
! 1. internal wave breaking as parameterized by Simmons etal.
!
! 2. barotropic tides feeling the bottom drag, as parameterized by 
!    Lee etal.  
!
! </DESCRIPTION>
!
  subroutine vert_mix_tidal(Time, Thickness, T_prog, Dens, diff_cbt, visc_cbu, visc_cbt, &
                      diff_cbt_wave, diff_cbt_leewave, diff_cbt_drag)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_prog_tracer_type),   intent(in)    :: T_prog(:)
  type(ocean_density_type),       intent(in)    :: Dens
  real, dimension(isd:,jsd:,:,:), intent(inout) :: diff_cbt
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbt
  real, dimension(isd:,jsd:,:),   intent(inout) :: diff_cbt_wave
  real, dimension(isd:,jsd:,:),   intent(inout) :: diff_cbt_leewave
  real, dimension(isd:,jsd:,:),   intent(inout) :: diff_cbt_drag

  integer :: i, j, k, kp1
  real    :: tmp 

  if(.not. use_this_module) return 

  if(use_legacy_methods) then

      call compute_bvfreq_legacy(Time, Thickness, T_prog, Dens)
      if(use_wave_dissipation)  then 
          call vert_mix_wave_legacy(Time, Thickness, diff_cbt, visc_cbu, visc_cbt, diff_cbt_wave)
      endif
      if(use_drag_dissipation)  then
          call vert_mix_drag_legacy(Time, Thickness, diff_cbt, visc_cbu, visc_cbt, diff_cbt_drag)
      endif
      diff_cbt_leewave(:,:,:) = 0.0
      diff_leewave(:,:,:)     = 0.0
      
  else 

      call compute_bvfreq(Time, Thickness, Dens)
      if(use_wave_dissipation)  then 
          call vert_mix_wave(Time, Thickness, Dens, diff_cbt, visc_cbu, visc_cbt, diff_cbt_wave, diff_cbt_leewave)
      endif
      if(use_drag_dissipation .and. horz_grid == MOM_BGRID)  then
          call vert_mix_drag_bgrid(Time, Thickness, diff_cbt, visc_cbu, visc_cbt, diff_cbt_drag)
      endif
      if(use_drag_dissipation .and. horz_grid == MOM_CGRID)  then
          call vert_mix_drag_cgrid(Time, Thickness, diff_cbt, visc_cbu, visc_cbt, diff_cbt_drag)
      endif

  endif

  ! add the background diffusivity and viscosity 
  do k=1,nk
     kp1 = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec
          diff_cbt(i,j,k,1) = Grd%tmask(i,j,kp1)*(diff_cbt(i,j,k,1) + background_diffusivity)
          diff_cbt(i,j,k,2) = Grd%tmask(i,j,kp1)*(diff_cbt(i,j,k,2) + background_diffusivity)
          visc_cbt(i,j,k)   = Grd%tmask(i,j,kp1)*(visc_cbt(i,j,k)   + background_viscosity)
          visc_cbu(i,j,k)   = Grd%umask(i,j,kp1)*(visc_cbu(i,j,k)   + background_viscosity)
        enddo
     enddo
  enddo

  ! compute power dissipated by mixing against stratification 
  if(id_power_diss_wave  > 0 .or. id_power_diss_drag    > 0 .or. &
     id_power_diss_tides > 0 .or. id_power_diss_leewave > 0 ) then 
      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               tmp = Thickness%dzt(i,j,k)*rho_N2(i,j,k)*Grd%tmask(i,j,k)
               wrk1(i,j,k) = tmp*diff_wave(i,j,k) 
               wrk2(i,j,k) = tmp*diff_drag(i,j,k) 
               wrk3(i,j,k) = tmp*diff_leewave(i,j,k) 
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_power_diss_wave, wrk1(:,:,:))
      call diagnose_3d(Time, Grd, id_power_diss_drag, wrk2(:,:,:))
      call diagnose_3d(Time, Grd, id_power_diss_leewave, wrk3(:,:,:))
      if (id_power_diss_tides > 0) then 
         call diagnose_3d(Time, Grd, id_power_diss_tides, wrk1(:,:,:)+wrk2(:,:,:))
      endif
  endif

  if (id_diff_cbt_tides > 0) then 
     call diagnose_3d(Time, Grd, id_diff_cbt_tides, diff_cbt_wave(:,:,:)+diff_cbt_drag(:,:,:))
  endif 

  ! recall unit Prandtl number 
  if (id_visc_cbt_tides > 0) then  
     call diagnose_3d(Time, Grd, id_visc_cbt_tides, diff_cbt_wave(:,:,:)+diff_cbt_drag(:,:,:))
  endif 

  if (id_visc_cbu_tides > 0) then
      wrk1=0.0
      do k=1,nk-1
         kp1=k+1
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k)          = Grd%umask(i,j,kp1)*onefourth  &
                   *(diff_cbt_wave(i,j,k)  +diff_cbt_wave(i+1,j,k)  &
                    +diff_cbt_wave(i,j+1,k)+diff_cbt_wave(i+1,j+1,k)&
                    +diff_cbt_drag(i,j,k)  +diff_cbt_drag(i+1,j,k)  &
                    +diff_cbt_drag(i,j+1,k)+diff_cbt_drag(i+1,j+1,k))
            enddo
         enddo
      enddo
      call diagnose_3d_u(Time, Grd, id_visc_cbu_tides, wrk1(:,:,:))
  endif


end subroutine vert_mix_tidal
! </SUBROUTINE> NAME="vert_mix_tidal"



!#######################################################################
! <SUBROUTINE NAME="compute_bvfreq">
!
! <DESCRIPTION>
! This subroutine computes the absolute value of rho*N^2 and abs of 
! N^2, with N^2 the squared Brunt-Vaisala (or buoyancy) frequency. 
!
! </DESCRIPTION>
!
  subroutine compute_bvfreq(Time, Thickness, Dens)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_density_type),     intent(in) :: Dens

  real    :: bottom
  real    :: tmp, rho_N2_prev
  integer :: i, j, k, m, kbot
  integer :: tau
  real, dimension(isd:ied,jsd:jed) :: roughness_klevel 

  tau                = Time%tau
  wrk1(:,:,:)        = 0.0
  wrk2(:,:,:)        = 0.0
  bvfreq_bottom(:,:) = 0.0
  

  ! absolute(rho*N^2) computed from ocean_density module calculation. 
  ! use the value at T-cell centre as this produces a smoother and
  ! better behaved bvfreq near the bottom, than does drhodz_wt. 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           rho_N2(i,j,k) = max(0.0,-grav*Dens%drhodz_zt(i,j,k)*Grd%tmask(i,j,k)) 
        enddo
     enddo
  enddo

  ! smooth rho_N2 in the vertical using a 1-2-1 filter
  if (smooth_rho_N2) then
      do m=1,num_121_passes
         do j=jsd,jed
            do i=isd,ied
               rho_N2_prev = onefourth*rho_N2(i,j,1)
               kbot=Grd%kmt(i,j)
               if (kbot>3) then
                   do k=2,kbot-2
                      tmp           = rho_N2(i,j,k)
                      rho_N2(i,j,k) = rho_N2_prev + onehalf*rho_N2(i,j,k) + onefourth*rho_N2(i,j,k+1)
                      rho_N2_prev   = onefourth*tmp
                   enddo
               endif
            enddo
         enddo
      enddo
  endif

  ! compute absolute value of buoyancy frequency, 
  ! using rho0r as an approximation to 1/rho.  
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           bvfreq(i,j,k) = sqrt(rho0r*rho_N2(i,j,k))
        enddo
     enddo
  enddo

  ! determine k-level at top of the bottom roughness boundary layer  
  wrk1_2d(:,:) = 0.0
  do j=jsd,jed
     do i=isd,ied
        kbot                  = Grd%kmt(i,j)
        roughness_klevel(:,:) = kbot
        if(kbot > 1) then 
            bottom = Thickness%depth_zwt(i,j,kbot)
            kloop: do k=kbot,1,-1
               tmp = Thickness%depth_zwt(i,j,k) + roughness_amp(i,j)
               if(tmp <= bottom) then 
                   wrk1_2d(i,j)          = k
                   roughness_klevel(i,j) = k
                   exit kloop
               endif
            enddo  kloop
        endif
     enddo
  enddo

  
  ! set bvfreq in bottom equal to value at kmt-1  
  do j=jsd,jed
     do i=isd,ied
        bvfreq_bottom(i,j) = 0.0
        if(Grd%kmt(i,j) > 1) then 
            kbot=Grd%kmt(i,j)-1
            bvfreq_bottom(i,j) = bvfreq(i,j,kbot)
        endif
     enddo
  enddo

  ! horizontal laplacian smoothing on the bottom bvfreq
  if(smooth_bvfreq_bottom) then 
     bvfreq_bottom(:,:) = bvfreq_bottom(:,:) + dtime*LAP_T(bvfreq_bottom(:,:),smooth_lap(:,:))
     call mpp_update_domains(bvfreq_bottom(:,:), Dom%domain2d) 
  endif 


  ! diagnostics 
  call diagnose_2d(TIme, Grd, id_roughness_klevel, wrk1_2d(:,:))
  call diagnose_2d(Time, Grd, id_bvfreq_bottom, bvfreq_bottom(:,:))
  call diagnose_3d(Time, Grd, id_bvfreq, bvfreq(:,:,:))


end subroutine compute_bvfreq
! </SUBROUTINE> NAME="compute_bvfreq"



!#######################################################################
! <SUBROUTINE NAME="vert_mix_wave">
!
! <DESCRIPTION>
! This subroutine computes dia-surface tracer diffusivity based on the 
! methods of Simmons et al., which consider dissipation from breaking
! internal gravity waves and their conversion into local dia-surface
! mixing, which is parameterized as diffusion.  
!
! Also compute a prototype parameterization of mixing due to 
! breaking leewaves from Nikurashin. 
!
! We assume a unit Prandtl number.
!
! Note that if umask(i,j,k) is 1.0, then so is 
! tmask(i,j,k), tmask(i+1,j,k), tmask(i,j+1,k), and tmask(i+1,j+1,k).
! So there is no need to compute the "active_cells" when doing the 
! space average to go from t-cell to u-cell to compute visc_cbu.
!
! </DESCRIPTION>
!
  subroutine vert_mix_wave(Time, Thickness, Dens, diff_cbt, visc_cbu, visc_cbt, diff_cbt_wave, diff_cbt_leewave)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  real, dimension(isd:,jsd:,:,:), intent(inout) :: diff_cbt
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbt
  real, dimension(isd:,jsd:,:),   intent(inout) :: diff_cbt_wave
  real, dimension(isd:,jsd:,:),   intent(inout) :: diff_cbt_leewave

  integer :: i, j, k, kbot, kp1, tau
  real    :: deposition, factor, tmp, rho_tmp 

  tau = Time%tau
  diff_wave(:,:,:)    = 0.0 ! diffusivity from wave scheme 
  diff_leewave(:,:,:) = 0.0 ! diffusivity from leewave scheme 
  wrk1(:,:,:)         = 0.0 ! viscosity   from wave scheme 
  wrk2(:,:,:)         = 0.0 ! mix_efficiency / rho_N2 

  ! compute mask for regions that are deemed too shallow for this scheme 
  do j=jsd,jed
     do i=isd,ied
        kbot=Grd%kmt(i,j)
        tmask_deep(i,j) = 0.0
        if(kbot > 1) then 
            if(Thickness%depth_zwt(i,j,kbot) > shelf_depth_cutoff) tmask_deep(i,j) = 1.0
        endif
     enddo
  enddo

  ! energy flux array (W/m2) (Simmons etal equation (1))
  if(fixed_wave_dissipation) then 
      do j=jsd,jed
         do i=isd,ied
            energy_flux(i,j) = min(wave_energy_flux_max, wave_dissipation(i,j)*tmask_deep(i,j)) 
         enddo
      enddo
  else 
      do j=jsd,jed
         do i=isd,ied
            energy_flux(i,j) = min(wave_energy_flux_max, wave_term(i,j)*bvfreq_bottom(i,j)*tmask_deep(i,j)) 
         enddo
      enddo
  endif

  ! compute mixing efficiency function 
  if(mixing_efficiency_n2depend) then
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               rho_tmp               = Dens%rho(i,j,k,tau) + epsln
               mix_efficiency(i,j,k) = mixing_efficiency*Grd%tmask(i,j,k) &
                                       *rho_N2(i,j,k)/(rho_N2(i,j,k) + rho_tmp*omega_earth2)
               wrk2(i,j,k)           = mixing_efficiency*Grd%tmask(i,j,k) &
                                       /(rho_N2(i,j,k) + rho_tmp*omega_earth2)
            enddo
         enddo
      enddo
  else 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               mix_efficiency(i,j,k) = mixing_efficiency*Grd%tmask(i,j,k)
               wrk2(i,j,k)           = mixing_efficiency*Grd%tmask(i,j,k) &
                                       /(rho_N2(i,j,k) + epsln)
            enddo
         enddo
      enddo
  endif


  ! diffusivity calculation (Simmons etal equation (3))
  do j=jsd,jed
     do i=isd,ied

        kbot=Grd%kmt(i,j)
        if(kbot > 1) then

            ! normalization of vertical structure function.
            ! Ensure it integrates to unity on the discrete grid.    
            ! "factor" approx decay_scale_inv/(exp[(H+eta)*decay_scale_inv]-1.0)  
            factor = 0.0
            do k=1,kbot-1
               factor = factor + Thickness%dzt(i,j,k)*exp(decay_scale_inv*Thickness%depth_zwt(i,j,k))
            enddo
            factor = 1.0/factor 

            do k=1,kbot-1
               deposition          = factor*exp(decay_scale_inv*Thickness%depth_zwt(i,j,k))
               tmp                 = Grd%tmask(i,j,k+1)*wrk2(i,j,k)*tidal_diss_efficiency*deposition
               diff_wave(i,j,k)    = tmp*energy_flux(i,j)
               diff_leewave(i,j,k) = tmp*leewave_dissipation(i,j)
               diff_wave(i,j,k)    = min(diff_wave(i,j,k),max_wave_diffusivity) 
               diff_leewave(i,j,k) = min(diff_leewave(i,j,k),max_wave_diffusivity) 
            enddo

        endif

     enddo
  enddo


  ! ensure diffusivity monotonically decreases as move upward in column.
  ! recall that diff_wave(i,j,k) is the diffusivity at the bottom of cell-k, 
  ! where diff_wave(i,j,kbot)=0.0 by definition. This prompts the kbot-2,1,-1
  ! loop limits.   
  if(wave_diffusivity_monotonic) then 
      do j=jsd,jed
         do i=isd,ied
            kbot=Grd%kmt(i,j)
            if(kbot > 1) then 
                do k=kbot-2,1,-1
                   diff_wave(i,j,k)    = min(diff_wave(i,j,k),diff_wave(i,j,k+1))
                   diff_leewave(i,j,k) = min(diff_leewave(i,j,k),diff_leewave(i,j,k+1))
                enddo
            endif
         enddo
      enddo
  endif

  ! add wave induced diffusivity and viscosity to diff_cbt and visc_cbu 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  do k=1,nk-1
     kp1=k+1
     do j=jsc,jec
        do i=isc,iec
           diff_cbt_wave(i,j,k)    = diff_wave(i,j,k) 
           diff_cbt_leewave(i,j,k) = diff_leewave(i,j,k) 
           diff_cbt(i,j,k,1)    = diff_cbt(i,j,k,1) + diff_wave(i,j,k) + diff_leewave(i,j,k)
           diff_cbt(i,j,k,2)    = diff_cbt(i,j,k,2) + diff_wave(i,j,k) + diff_leewave(i,j,k)
           visc_cbt(i,j,k)      = visc_cbt(i,j,k)   + diff_wave(i,j,k) + diff_leewave(i,j,k)
           wrk1(i,j,k)          = Grd%umask(i,j,kp1)*onefourth                   &
                                 *(diff_wave(i,j,k)     +diff_wave(i+1,j,k)      &
                                  +diff_wave(i,j+1,k)   +diff_wave(i+1,j+1,k))
           wrk2(i,j,k)          = Grd%umask(i,j,kp1)*onefourth                   &
                                 *(diff_leewave(i,j,k)  +diff_leewave(i+1,j,k)   &
                                  +diff_leewave(i,j+1,k)+diff_leewave(i+1,j+1,k))
           visc_cbu(i,j,k)      = visc_cbu(i,j,k) + wrk1(i,j,k) + wrk2(i,j,k)
        enddo
     enddo
  enddo


  ! send some diagnostics 

  call diagnose_3d(Time, Grd, id_mix_efficiency, mix_efficiency(:,:,:))
  call diagnose_2d(Time, Grd, id_energy_flux, energy_flux(:,:))
  call diagnose_2d(Time, Grd, id_power_waves, Grd%dat(:,:)*energy_flux(:,:))
  call diagnose_3d(Time, Grd, id_diff_cbt_wave, diff_wave(:,:,:))
  call diagnose_3d(Time, Grd, id_diff_cbt_leewave, diff_leewave(:,:,:))
  call diagnose_3d(Time, Grd, id_visc_cbt_wave, diff_wave(:,:,:))
  call diagnose_3d(Time, Grd, id_visc_cbt_leewave, diff_leewave(:,:,:))
  call diagnose_3d_u(Time, Grd, id_visc_cbu_wave, wrk1(:,:,:))
  call diagnose_3d_u(Time, Grd, id_visc_cbu_leewave, wrk2(:,:,:))


end subroutine vert_mix_wave
! </SUBROUTINE> NAME="vert_mix_wave"


!#######################################################################
! <SUBROUTINE NAME="vert_mix_drag_bgrid">
!
! <DESCRIPTION>
! This subroutine computes dia-surface tracer diffusivity based on the 
! methods of Lee etal., which consider the dissipation from barotropic tides
! rubbing against the ocean bottom. 
!
! We assume B-grid layout for the velocity 
!
! We assume a unit Prandtl number, so compute the viscosity as a four-point
! average of the diffusivity. 
!
! We perform various averages here in order to smooth Richardson number.
!
! 1. compute Richardson number on U-cell by averaging bvfreq from T-cell
! 2. average U-cell Richardson number to then get T-cell diffusivity 
! 3. average T-cell diffusivity to get U-cell viscosity. 
!
! Note that if umask(i,j,k)==1.0, then so is tmask(i,j,k), tmask(i+1,j,k), 
! tmask(i,j+1,k), and tmask(i+1,j+1,k). So there is no need to compute 
! active_cells when averaging from T-cell to U-cell.
!
! </DESCRIPTION>
!
  subroutine vert_mix_drag_bgrid(Time, Thickness, diff_cbt, visc_cbu, visc_cbt, diff_cbt_drag)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  real, dimension(isd:,jsd:,:,:), intent(inout) :: diff_cbt
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbt
  real, dimension(isd:,jsd:,:),   intent(inout) :: diff_cbt_drag

  integer :: i, j, k, kbot, kp1
  real    :: height, bottom
  real    :: bvfreq_u, speedr, active_cells

  wrk1(:,:,:)      =0.0  ! Richardson number on U-cell
  wrk2(:,:,:)      =0.0  ! Richardson number on T-cell
  wrk3(:,:,:)      =0.0  ! viscosity from drag scheme 
  diff_drag(:,:,:) =0.0  ! diffusivity from drag scheme 

  ! Richardson number on U-cell. 
  ! perform a 4-point average of T-cell bvfreq 
  ! and then divide by the U-cell tidal speed term. 
  ! tide_speed_mask is useful to reduce overflows
  ! in later calculation of the diffusivity.
  do j=jsd,jed-1
     do i=isd,ied-1
        kbot=Grd%kmu(i,j)
        if(kbot>1) then 
            bottom = Thickness%depth_zwu(i,j,kbot)
            speedr = tide_speed_mask(i,j)/(epsln+rescaled_speed_u(i,j))
            do k=1,kbot-1
               kp1=k+1 
               height      = bottom-Thickness%depth_zwu(i,j,k)
               bvfreq_u    = onefourth*(bvfreq(i,j,k)+bvfreq(i+1,j,k)+bvfreq(i,j+1,k)+bvfreq(i+1,j+1,k))            
               wrk1(i,j,k) = 2.0*Grd%umask(i,j,kp1)*(bvfreq_u*height*speedr)**2
            enddo
        endif
     enddo
  enddo

  ! Richardson number on bottom of T-cells.
  ! need active_cells for averaging operation.
  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec
           active_cells =   Grd%umask(i,j,k)   + Grd%umask(i-1,j,k) &
                          + Grd%umask(i,j-1,k) + Grd%umask(i-1,j-1,k) + epsln
           wrk2(i,j,k)  =  (wrk1(i,j,k) + wrk1(i-1,j,k) + wrk1(i,j-1,k) + wrk1(i-1,j-1,k))/active_cells
        enddo
     enddo
  enddo


  ! compute drag induced diffusivity 
  ! (Lee etal equations (1), (2), and (3))
  ! Multiply by tide_speed_mask so to zero out 
  ! regions with tiny tide speeds, which are regions 
  ! where we do not wish to have any enhanced mixing
  ! arising from the barotropic tide mixing parameterization
  ! anyhow.  
  do k=1,nk-1
     kp1=k+1
     do j=jsc,jec
        do i=isc,iec
           diff_drag(i,j,k) = Grd%tmask(i,j,kp1)*tide_speed_mask(i,j)*max_drag_diffusivity &
                              *(1.0 + munk_anderson_sigma*wrk2(i,j,k))**(-munk_anderson_p)
        enddo
     enddo
  enddo

  if(drag_dissipation_efold) then 
      do j=jsc,jec
         do i=isc,iec
            kbot=Grd%kmt(i,j)
            if(kbot>1) then 
                bottom = Thickness%depth_zwt(i,j,kbot)
                do k=1,kbot-1
                   kp1=k+1 
                   height           = bottom-Thickness%depth_zwt(i,j,k)
                   diff_drag(i,j,k) = diff_drag(i,j,k)*exp(-height*efold_depth_r(i,j))
                enddo
            endif
         enddo
      enddo
  endif

  call mpp_update_domains(diff_drag(:,:,:), Dom%domain2d) 


  ! add drag induced diffusivity and viscosity to diff_cbt and visc_cbu. 
  ! average t-cell diffusivities to get u-cell viscosity.
  do k=1,nk-1
     kp1=k+1
     do j=jsc,jec
        do i=isc,iec
           diff_cbt_drag(i,j,k) = diff_drag(i,j,k)
           diff_cbt(i,j,k,1)    = diff_cbt(i,j,k,1) + diff_drag(i,j,k)
           diff_cbt(i,j,k,2)    = diff_cbt(i,j,k,2) + diff_drag(i,j,k)
           visc_cbt(i,j,k)      = visc_cbt(i,j,k)   + diff_drag(i,j,k)
           wrk3(i,j,k)          = Grd%umask(i,j,kp1)*onefourth              &
                                 *(diff_drag(i,j,k)  +diff_drag(i+1,j,k)    &
                                  +diff_drag(i,j+1,k)+diff_drag(i+1,j+1,k))
           visc_cbu(i,j,k)      = visc_cbu(i,j,k) + wrk3(i,j,k)
        enddo
     enddo
  enddo


  call diagnose_3d(Time, Grd, id_rinumber_drag, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_diff_cbt_drag, diff_drag(:,:,:))
  call diagnose_3d_u(Time, Grd, id_visc_cbu_drag, wrk3(:,:,:))

end subroutine vert_mix_drag_bgrid
! </SUBROUTINE> NAME="vert_mix_drag_bgrid"


!#######################################################################
! <SUBROUTINE NAME="vert_mix_drag_cgrid">
!
! <DESCRIPTION>
! This subroutine computes dia-surface tracer diffusivity based on the 
! methods of Lee etal., which consider the dissipation from barotropic tides
! rubbing against the ocean bottom. 
!
! We assume a unit Prandtl number, so compute the viscosity as a four-point
! average of the diffusivity. 
!
! We assume C-grid layout for the velocity, which renders slight
! distinctions for the calculation of Richardson number. Otherwise, the 
! calculations are the same as the Bgrid. We introduce this separate
! routine, however, to enable easier bitwise agreement with older 
! model results.  Also, further development of this scheme may lead
! to more distinctions from the Bgrid.  
!
! </DESCRIPTION>
!
  subroutine vert_mix_drag_cgrid(Time, Thickness, diff_cbt, visc_cbu, visc_cbt, diff_cbt_drag)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  real, dimension(isd:,jsd:,:,:), intent(inout) :: diff_cbt
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbt
  real, dimension(isd:,jsd:,:),   intent(inout) :: diff_cbt_drag

  integer :: i, j, k, kbot, kp1
  real    :: height, bottom
  real    :: speedr, active_cells

  wrk1(:,:,:)      = 0.0  ! raw Richardson number on T-cell
  wrk2(:,:,:)      = 0.0  ! smoothed Richardson number on T-cell
  wrk3(:,:,:)      = 0.0  ! visc_cbu from drag scheme 
  diff_drag(:,:,:) = 0.0  ! diffusivity from drag scheme 

  ! Richardson number on T-cell. 
  do j=jsd,jed
     do i=isd,ied
        kbot=Grd%kmt(i,j)
        if(kbot>1) then 
            bottom = Thickness%depth_zwt(i,j,kbot)
            speedr = tide_speed_mask(i,j)/(epsln+rescaled_speed_t(i,j))
            do k=1,kbot-1
               kp1=k+1 
               height      = bottom-Thickness%depth_zwt(i,j,k)
               wrk1(i,j,k) = 2.0*Grd%tmask(i,j,kp1)*(bvfreq(i,j,k)*height*speedr)**2
               wrk2(i,j,k) = wrk1(i,j,k)
            enddo
        endif
     enddo
  enddo

  ! perform 9point average to smooth, and to be more consistent 
  ! with the Bgrid approach.
  ! need active_cells for averaging operation.
  if(smooth_ri_drag_cgrid) then 
     do k=1,nk-1
        do j=jsc,jec
           do i=isc,iec
              active_cells =   Grd%tmask(i-1,j-1,k) + Grd%tmask(i,j-1,k) + Grd%tmask(i+1,j-1,k) &
                             + Grd%tmask(i,j-1,k)   + Grd%tmask(i,j,k)   + Grd%tmask(i+1,j,k)   &
                             + Grd%tmask(i-1,j+1,k) + Grd%tmask(i,j+1,k) + Grd%tmask(i+1,j+1,k) &
                             + epsln
              wrk2(i,j,k)  =  (  wrk1(i-1,j-1,k) + wrk1(i,j-1,k) + wrk1(i+1,j-1,k)   &
                               + wrk1(i,j-1,k)   + wrk1(i,j,k)   + wrk1(i+1,j,k)     &
                               + wrk1(i-1,j+1,k) + wrk1(i,j+1,k) + wrk1(i+1,j+1,k)) / active_cells
           enddo
        enddo
     enddo
  endif 

  ! compute drag induced diffusivity 
  ! (Lee etal equations (1), (2), and (3))
  ! Multiply by tide_speed_mask so to zero out 
  ! regions with tiny tide speeds, which are regions 
  ! where we do not wish to have any enhanced mixing
  ! arising from the barotropic tide mixing parameterization
  ! anyhow.  
  do k=1,nk-1
     kp1=k+1
     do j=jsc,jec
        do i=isc,iec
           diff_drag(i,j,k) = Grd%tmask(i,j,kp1)*tide_speed_mask(i,j)*max_drag_diffusivity &
                              *(1.0 + munk_anderson_sigma*wrk2(i,j,k))**(-munk_anderson_p)
        enddo
     enddo
  enddo

  if(drag_dissipation_efold) then 
      do j=jsc,jec
         do i=isc,iec
            kbot=Grd%kmt(i,j)
            if(kbot>1) then 
                bottom = Thickness%depth_zwt(i,j,kbot)
                do k=1,kbot-1
                   kp1=k+1 
                   height           = bottom-Thickness%depth_zwt(i,j,k)
                   diff_drag(i,j,k) = diff_drag(i,j,k)*exp(-height*efold_depth_r(i,j))
                enddo
            endif
         enddo
      enddo
  endif

  call mpp_update_domains(diff_drag(:,:,:), Dom%domain2d) 


  ! add drag induced diffusivity and viscosity to diff_cbt and visc_cbu. 
  ! average t-cell diffusivities to get u-cell viscosity.
  do k=1,nk-1
     kp1=k+1
     do j=jsc,jec
        do i=isc,iec
           diff_cbt_drag(i,j,k) = diff_drag(i,j,k)
           diff_cbt(i,j,k,1)    = diff_cbt(i,j,k,1) + diff_drag(i,j,k)
           diff_cbt(i,j,k,2)    = diff_cbt(i,j,k,2) + diff_drag(i,j,k)
           visc_cbt(i,j,k)      = visc_cbt(i,j,k)   + diff_drag(i,j,k)
           wrk3(i,j,k)          = Grd%umask(i,j,kp1)*onefourth              &
                                 *(diff_drag(i,j,k)  +diff_drag(i+1,j,k)    &
                                  +diff_drag(i,j+1,k)+diff_drag(i+1,j+1,k))
           visc_cbu(i,j,k)      = visc_cbu(i,j,k) + wrk3(i,j,k)
        enddo
     enddo
  enddo


  call diagnose_3d(Time, Grd, id_rinumber_drag, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_diff_cbt_drag, diff_drag(:,:,:))
  call diagnose_3d_u(Time, Grd, id_visc_cbu_drag, wrk3(:,:,:))

end subroutine vert_mix_drag_cgrid
! </SUBROUTINE> NAME="vert_mix_drag_cgrid"



!#######################################################################
! <SUBROUTINE NAME="compute_bvfreq_legacy">
!
! <DESCRIPTION>
! This subroutine computes the absolute value of rho*N^2 and abs of 
! N^2, with N^2 the squared Brunt-Vaisala (or buoyancy) frequency. 
!
! This routine employs a legacy approach, which is not recommended.
! It remains solely to allow exact reproduction of older results.
!
! </DESCRIPTION>
!
  subroutine compute_bvfreq_legacy(Time, Thickness, T_prog, Dens)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  type(ocean_density_type),     intent(in) :: Dens

  real    :: rho_inv, drhodz
  real    :: tmp, rho_N2_prev, rho_tmp
  integer :: i, j, k, m, kp1, kbot
  integer :: tau

  tau = Time%tau
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0

  ! partial derivatives of density wrt to temperature and salinity 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           drhodT(i,j,k) = Dens%drhodT(i,j,k)
           drhodS(i,j,k) = Dens%drhodS(i,j,k)
        enddo
     enddo
  enddo


  ! vertical derivative of temperature and salinity at bottom of tracer cells 
  do k=1,nk
     kp1=min(k+1,nk)
     do j=jsd,jed
        do i=isd,ied
           tmp         = Grd%tmask(i,j,kp1)/Thickness%dzwt(i,j,k)
           wrk1(i,j,k) = tmp*(T_prog(index_temp)%field(i,j,k,tau)-T_prog(index_temp)%field(i,j,kp1,tau)) 
           wrk2(i,j,k) = tmp*(Dens%rho_salinity(i,j,k,tau)-Dens%rho_salinity(i,j,kp1,tau)) 
        enddo
     enddo
  enddo

  ! absolute(rho*N^2) computed from vertical derivative of "neutral density" 
  do k=1,nk
     kp1 = min(k+1,nk)
     do j=jsd,jed
        do i=isd,ied
           drhodz        = onehalf*( (drhodT(i,j,k)+drhodT(i,j,kp1))*wrk1(i,j,k) &
                                    +(drhodS(i,j,k)+drhodS(i,j,kp1))*wrk2(i,j,k) )
           drhodz        = min(drhodz,-drhodz_min)*Grd%tmask(i,j,kp1) 
           rho_N2(i,j,k) = -grav*drhodz
        enddo
     enddo
  enddo

  ! smooth rho_N2 in the vertical using a 1-2-1 filter
  if (smooth_rho_N2) then
      do m=1,num_121_passes
         do j=jsd,jed
            do i=isd,ied
               rho_N2_prev = onefourth*rho_N2(i,j,1)
               kbot=Grd%kmt(i,j)
               if (kbot>3) then
                   do k=2,kbot-2
                      tmp           = rho_N2(i,j,k)
                      rho_N2(i,j,k) = rho_N2_prev + onehalf*rho_N2(i,j,k) + onefourth*rho_N2(i,j,k+1)
                      rho_N2_prev   = onefourth*tmp
                   enddo
               endif
            enddo
         enddo
      enddo
  endif

  ! compute buoyancy frequency 
  do k=1,nk
     kp1 = min(k+1,nk)
     do j=jsd,jed
        do i=isd,ied
           rho_inv       = 2.0/(epsln + Dens%rho(i,j,k,tau) + Dens%rho(i,j,kp1,tau)) 
           bvfreq(i,j,k) = sqrt(rho_inv*rho_N2(i,j,k))
        enddo
     enddo
  enddo

  ! bvfreq at the bottom.
  ! set kbot=kmt-1 rather than kbot=kmt, since N^2=0 
  ! at bottom of bottom-most tracer cell, by definition.
  do j=jsd,jed
     do i=isd,ied
        bvfreq_bottom(i,j) = 0.0
        if(Grd%kmt(i,j) > 1) then 
            kbot=Grd%kmt(i,j)-1
            bvfreq_bottom(i,j) = bvfreq(i,j,kbot)
        endif
     enddo
  enddo

  ! horizontal laplacian smoothing on the bottom bvfreq to reduce noise 
  if(smooth_bvfreq_bottom) then 
     bvfreq_bottom(:,:) = bvfreq_bottom(:,:) + dtime*LAP_T(bvfreq_bottom(:,:),smooth_lap(:,:))
     call mpp_update_domains(bvfreq_bottom(:,:), Dom%domain2d) 
  endif 

  ! compute mixing efficiency 
  mix_efficiency(:,:,:) = mixing_efficiency
  if(mixing_efficiency_n2depend) then
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               rho_tmp = Dens%rho(i,j,k,tau) + epsln
               mix_efficiency(i,j,k) = mixing_efficiency*rho_N2(i,j,k)/(rho_N2(i,j,k) + rho_tmp*omega_earth2)
            enddo
         enddo
      enddo
  endif

  call diagnose_2d(Time, Grd, id_bvfreq_bottom, bvfreq_bottom(:,:))
  call diagnose_3d(Time, Grd, id_bvfreq, bvfreq(:,:,:))
  call diagnose_3d(Time, Grd, id_mix_efficiency, mix_efficiency(:,:,:))

end subroutine compute_bvfreq_legacy
! </SUBROUTINE> NAME="compute_bvfreq_legacy"


!#######################################################################
! <SUBROUTINE NAME="vert_mix_wave_legacy">
!
! <DESCRIPTION>
!
! Legacy routine maintained only to exactly reproduce older results.
! It is not recommended for new experiments, as it uses some obsolete
! methods. 
!
! This subroutine computes dia-surface tracer diffusivity based on the 
! methods of Simmons etal., which consider the dissipation from breaking
! internal gravity waves and their conversion into local dia-surface
! diffusion.  
!
! We assume a unit Prandtl number, so compute the viscosity as a four-point
! average of the diffusivity. 
!
! Note that if umask(i,j,k) is 1.0, then so is 
! tmask(i,j,k), tmask(i+1,j,k), tmask(i,j+1,k), and tmask(i+1,j+1,k).
! So there is no need to compute the "active_cells" when doing the 
! space average to go from t-cell to u-cell to compute viscosity.
!
! </DESCRIPTION>
!
  subroutine vert_mix_wave_legacy(Time, Thickness, diff_cbt, visc_cbu, visc_cbt, diff_cbt_wave)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  real, dimension(isd:,jsd:,:,:), intent(inout) :: diff_cbt
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbt
  real, dimension(isd:,jsd:,:),   intent(inout) :: diff_cbt_wave

  integer :: i, j, k, kbot, kp1
  real    :: deposition, factor 

  diff_wave(:,:,:) = 0.0 ! diffusivity from wave scheme 
  wrk1(:,:,:)      = 0.0 ! viscosity   from wave scheme 

  ! compute mask for regions that are too shallow for this scheme 
  do j=jsd,jed
     do i=isd,ied
        kbot=Grd%kmt(i,j)
        tmask_deep(i,j) = 0.0
        if(kbot > 1) then 
            if(Thickness%depth_zwt(i,j,kbot) > shelf_depth_cutoff) tmask_deep(i,j) = 1.0
        endif
     enddo
  enddo

  ! compute the wave energy flux array (W/m2) and save for diagnostics
  ! (Simmons etal equation (1))
  if(fixed_wave_dissipation) then 
      do j=jsd,jed
         do i=isd,ied
            energy_flux(i,j) = min(wave_energy_flux_max, wave_dissipation(i,j)*tmask_deep(i,j)) 
         enddo
      enddo
  else 
      do j=jsd,jed
         do i=isd,ied
            energy_flux(i,j) = min(wave_energy_flux_max, wave_term(i,j)*bvfreq_bottom(i,j)*tmask_deep(i,j)) 
         enddo
      enddo
  endif

  ! compute wave induced diffusivity 
  ! (Simmons etal equation (3))
  do j=jsd,jed
     do i=isd,ied
        kbot=Grd%kmt(i,j)
        if(kbot > 1) then

            ! normalization of vertical structure function...ensure it 
            ! integrates to unity on the discrete grid.   
            factor = 0.0
            do k=1,kbot-1
               factor = factor + Thickness%dzt(i,j,k)*exp(decay_scale_inv*Thickness%depth_zwt(i,j,k))
            enddo
            factor = 1.0/factor 

            ! calculate diffusivity 
            do k=1,kbot-1
               deposition       = factor*exp(decay_scale_inv*Thickness%depth_zwt(i,j,k))
               diff_wave(i,j,k) = Grd%tmask(i,j,k+1)*mix_efficiency(i,j,k)*tidal_diss_efficiency &
                                  *energy_flux(i,j)*deposition/(epsln+rho_N2(i,j,k))
               diff_wave(i,j,k) = min(diff_wave(i,j,k),max_wave_diffusivity) 
            enddo

        endif
     enddo
  enddo

  ! ensure diffusivity monotonically decreases as move upward in column.
  ! recall that diff_wave(i,j,k) is the diffusivity at the bottom of cell-k, 
  ! where diff_wave(i,j,kbot)=0.0 by definition. This prompts the kbot-2,1,-1
  ! loop limits.   
  if(wave_diffusivity_monotonic) then 
      do j=jsd,jed
         do i=isd,ied
            kbot=Grd%kmt(i,j)
            if(kbot > 1) then 
                do k=kbot-2,1,-1
                   diff_wave(i,j,k) = min(diff_wave(i,j,k),diff_wave(i,j,k+1))
                enddo
            endif
         enddo
      enddo
  endif

  ! add wave induced diffusivity and viscosity to diff_cbt and visc_cbu 
  do k=1,nk-1
     kp1=k+1
     do j=jsc,jec
        do i=isc,iec
           diff_cbt_wave(i,j,k) = diff_wave(i,j,k) 
           diff_cbt(i,j,k,1)    = diff_cbt(i,j,k,1) + diff_wave(i,j,k)
           diff_cbt(i,j,k,2)    = diff_cbt(i,j,k,2) + diff_wave(i,j,k)
           visc_cbt(i,j,k)      = visc_cbt(i,j,k)   + diff_wave(i,j,k)
           wrk1(i,j,k)          = Grd%umask(i,j,kp1)*onefourth             &
                                 *(diff_wave(i,j,k)  +diff_wave(i+1,j,k)   &
                                  +diff_wave(i,j+1,k)+diff_wave(i+1,j+1,k))
           visc_cbu(i,j,k)      = visc_cbu(i,j,k) + wrk1(i,j,k)
        enddo
     enddo
  enddo


  call diagnose_2d(Time, Grd, id_energy_flux, energy_flux(:,:))
  call diagnose_2d(Time, Grd, id_power_waves, Grd%dat(:,:)*energy_flux(:,:))
  call diagnose_3d(Time, Grd, id_diff_cbt_wave, diff_wave(:,:,:))
  call diagnose_3d(Time, Grd, id_visc_cbt_wave, diff_wave(:,:,:))
  call diagnose_3d_u(Time, Grd, id_visc_cbu_wave, wrk1(:,:,:))

end subroutine vert_mix_wave_legacy
! </SUBROUTINE> NAME="vert_mix_wave_legacy"


!#######################################################################
! <SUBROUTINE NAME="vert_mix_drag_legacy">
!
! <DESCRIPTION>
!
! Legacy routine maintained only to exactly reproduce older results.
! It is not recommended for new experiments, as it uses some obsolete
! methods. 
!
! This subroutine computes dia-surface tracer diffusivity based on the 
! methods of Lee etal., which consider the dissipation from barotropic tides
! rubbing against the ocean bottom. 
!
! We assume a unit Prandtl number, so compute the viscosity as a four-point
! average of the diffusivity. 
!
! We perform various averages here in order to smooth Richardson number.
!
! 1. compute Richardson number on U-cell by averaging bvfreq from T-cell
! 2. average U-cell Richardson number to then get T-cell diffusivity 
! 3. average T-cell diffusivity to get U-cell viscosity. 
!
! Note that if umask(i,j,k)==1.0, then so is tmask(i,j,k), tmask(i+1,j,k), 
! tmask(i,j+1,k), and tmask(i+1,j+1,k). So there is no need to compute 
! active_cells when averaging from T-cell to U-cell.
!
! </DESCRIPTION>
!
  subroutine vert_mix_drag_legacy(Time, Thickness, diff_cbt, visc_cbu, visc_cbt, diff_cbt_drag)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  real, dimension(isd:,jsd:,:,:), intent(inout) :: diff_cbt
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbt
  real, dimension(isd:,jsd:,:),   intent(inout) :: diff_cbt_drag

  integer :: i, j, k, kbot, kp1
  real    :: height, bottom, bvfreq_u, speedr, active_cells

  wrk1(:,:,:)      =0.0  ! Richardson number on U-cell
  wrk2(:,:,:)      =0.0  ! Richardson number on T-cell
  wrk3(:,:,:)      =0.0  ! viscosity from drag scheme 
  diff_drag(:,:,:) =0.0  ! diffusivity from drag scheme 

  ! Richardson number on U-cell. 
  ! perform a 4-point average of T-cell bvfreq 
  ! and then divide by the U-cell tidal speed term. 
  ! tide_speed_mask is useful to reduce overflows
  ! in later calculation of the diffusivity.
  do j=jsd,jed-1
     do i=isd,ied-1
        kbot=Grd%kmu(i,j)
        if(kbot>1) then 
            bottom = Thickness%depth_zwu(i,j,kbot)
            speedr = tide_speed_mask(i,j)/(epsln+rescaled_speed_u(i,j))
            do k=1,kbot-1
               kp1=k+1 
               height      = bottom-Thickness%depth_zwu(i,j,k)
               bvfreq_u    = onefourth*(bvfreq(i,j,k)+bvfreq(i+1,j,k)+bvfreq(i,j+1,k)+bvfreq(i+1,j+1,k))            
               wrk1(i,j,k) = 2.0*Grd%umask(i,j,kp1)*(bvfreq_u*height*speedr)**2
            enddo
        endif
     enddo
  enddo

  ! Richardson number on bottom of T-cells.
  ! need active_cells for averaging operation.
  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec
           active_cells =   Grd%umask(i,j,k)   + Grd%umask(i-1,j,k) &
                          + Grd%umask(i,j-1,k) + Grd%umask(i-1,j-1,k) + epsln
           wrk2(i,j,k)  =  (wrk1(i,j,k) + wrk1(i-1,j,k) + wrk1(i,j-1,k) + wrk1(i-1,j-1,k))/active_cells
        enddo
     enddo
  enddo

  ! compute drag induced diffusivity 
  ! (Lee etal equations (1), (2), and (3))
  ! Multiply by tide_speed_mask so to zero out 
  ! regions with tiny tide speeds, which are regions 
  ! where we do not wish to have any enhanced mixing
  ! arising from the barotropic tide mixing parameterization
  ! anyhow.  
  do k=1,nk-1
     kp1=k+1
     do j=jsc,jec
        do i=isc,iec
           diff_drag(i,j,k) = Grd%tmask(i,j,kp1)*tide_speed_mask(i,j)*max_drag_diffusivity &
                              *(1.0 + munk_anderson_sigma*wrk2(i,j,k))**(-munk_anderson_p)
        enddo
     enddo
  enddo
  call mpp_update_domains(diff_drag(:,:,:), Dom%domain2d) 


  ! add drag induced diffusivity and viscosity to diff_cbt and visc_cbu. 
  ! average t-cell diffusivities to get u-cell viscosity.
  do k=1,nk-1
     kp1=k+1
     do j=jsc,jec
        do i=isc,iec
           diff_cbt_drag(i,j,k) = diff_drag(i,j,k)
           diff_cbt(i,j,k,1)    = diff_cbt(i,j,k,1) + diff_drag(i,j,k)
           diff_cbt(i,j,k,2)    = diff_cbt(i,j,k,2) + diff_drag(i,j,k)
           visc_cbt(i,j,k)      = visc_cbt(i,j,k)   + diff_drag(i,j,k)
           wrk3(i,j,k)          = Grd%umask(i,j,kp1)*onefourth              &
                                 *(diff_drag(i,j,k)  +diff_drag(i+1,j,k)    &
                                  +diff_drag(i,j+1,k)+diff_drag(i+1,j+1,k))
           visc_cbu(i,j,k)      = visc_cbu(i,j,k) + wrk3(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_rinumber_drag, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_diff_cbt_drag, diff_drag(:,:,:))
  call diagnose_3d(Time, Grd, id_visc_cbt_drag, diff_drag(:,:,:))
  call diagnose_3d_u(Time, Grd, id_visc_cbu_drag, wrk3(:,:,:))


end subroutine vert_mix_drag_legacy
! </SUBROUTINE> NAME="vert_mix_drag_legacy"



end module ocean_vert_tidal_mod

module ocean_vert_mix_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> A. Rosati
!</REVIEWER>
!
!<OVERVIEW>
! Time tendency from vertical mixing 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module does the following:
! 
! --computes the vertical mixing coefficients for tracer and velocity,
! --computes thickness weighted and density weighted 
!   time tendency for tracer due to vertical diffusion processes,
! --computes the thickness weighted and density weighted 
!   acceleration for velocity due to vertical friction processes.
!
! Account is taken of the differences between Bgrid and Cgrid
! implementations of MOM.
!
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Kirk Bryan and L.J. Lewis
! A water mass model of the world ocean
! Journal of Geophysical Research (1979) vol 84, pages 2503--2517
! </REFERENCE>
!
! <REFERENCE>
! Elements of MOM (2012)
! S.M. Griffies 
! </REFERENCE>
!
! <REFERENCE>
! Henyey, F.S., J. Wright, and S.M. Flatte, 1986: Energy and
! action flow through the internal wave field: an eikonal approach.
! Journal of Geophysical Research, {\bf 91}, Issue C7, 8487--8496.
! </REFERENCE>
!
! <NOTE> 
! The Bryan-Lewis vertical diffusivity is small in the upper ocean and 
! increases with depth according to an inverse tangent profile.  The default
! values are from roughly 0.05e-5 m^2/sec to roughly 1.0e-4 m^2/sec.
! Latitudinally dependent Bryan-Lewis values are available.  
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_vert_mix_nml">
!
!  <DATA NAME="debug_this_module" TYPE="logical">
! For debugging purposes. 
!  </DATA> 
!
!  <DATA NAME="aidif" TYPE="real">
!  aidif=1 for implicit in time solution of the vertical mixing equation.
!  aidif=0 for explicit in time solution of the vertical mixing equation.
!  semi-implicit method with 0 < aidif < 1 is not fully supported in MOM. 
!  </DATA> 
!
!  <DATA NAME="use_explicit_vert_diffuse" TYPE="logical">
! Must be true to use time-explicit vertical tracer diffusion.
!  </DATA> 
!  <DATA NAME="verbose_init" TYPE="logical">
!  For verbose writes during initialization. 
!  </DATA> 
!  <DATA NAME="vert_mix_scheme" TYPE="character">
!  To determine the vertical mixing scheme: 
!  "const", "kpp", "kpp_mom4p0","kpp_mom4p1", "chen", "pp", or "gotm".
!  </DATA> 
!
!  <DATA NAME="vert_diff_back_via_max" TYPE="logical">
! If .true. then include a static background diffusivity 
! via the max function, as used in mom4p0d.  The alternative
! is via simply adding the background to the diffusivity 
! obtained via other approaches.  This option remains for 
! legacy. Default is vert_diff_back_via_max=.true.
!  </DATA> 
!
!  <DATA NAME="vert_visc_back" TYPE="logical">
! If .true. then include a static depth dependent vertical
! viscosity which is used only if running w/ constant 
! vertical viscosity scheme. Standard application is when 
! have a model with fine vertical resolution, yet no mixed 
! layer scheme.  Wind stress must be spread deeper than the 
! top cell, or the model may go unstable, or at the least it 
! will produce spuriously large vertical shears.  
!  </DATA> 
!
!  <DATA NAME="visc_cbu_back_max" UNITS="m^2/sec" TYPE="real">
!  For use in setting background vertical viscosity. 
!  </DATA> 
!  <DATA NAME="visc_cbu_back_min" UNITS="m^2/sec" TYPE="real">
!  For use in setting background vertical viscosity. 
!  </DATA> 
!  <DATA NAME="visc_cbu_back_zmid" UNITS="m" TYPE="real">
!  Mid-point of tanh function used to define background vertical viscosity. 
!  </DATA> 
!  <DATA NAME="visc_cbu_back_zwid" UNITS="m" TYPE="real">
!  Width of tanh function used to define background vertical viscosity. 
!  </DATA> 
!
!  <DATA NAME="diff_cbt_tanh" TYPE="logical">
!  For enabling tanh background vertical diffusivity profile. 
!  Default diff_cbt_tanh=.false.
!  </DATA> 
!  <DATA NAME="diff_cbt_tanh_max" UNITS="m^2/sec" TYPE="real">
!  For use in setting background vertical diffusivity.
!  Default diff_cbt_tanh_max=1e-3.
!  </DATA> 
!  <DATA NAME="diff_cbt_tanh_min" UNITS="m^2/sec" TYPE="real">
!  For use in setting background vertical diffusivity.
!  Default diff_cbt_tanh_min=2e-5.
!  </DATA> 
!  <DATA NAME="diff_cbt_tanh_zmid" UNITS="m" TYPE="real">
!  Mid-point of tanh function used to define background vertical diffusivity.
!  Default diff_cbt_tanh_zmid=150.0.
!  </DATA> 
!  <DATA NAME="diff_cbt_tanh_zwid" UNITS="m" TYPE="real">
!  Width of tanh function used to define background vertical diffusivity.
!  Default diff_cbt_tanh_zwid=30.0.
!  </DATA> 
!
!  <DATA NAME="hwf_diffusivity" TYPE="logical">
!  3D background diffusivity which gets smaller in equatorial region.
!  Based on the work of Henyey etal (1986).   
!  This scheme should NOT be used if Bryan-Lewis is used.  
!  Default hwf_diffusivity=.false. 
!  </DATA> 
!  <DATA NAME="hwf_diffusivity_3d" TYPE="logical">
!  3D background diffusivity which gets smaller in equatorial region.
!  Based on the work of Henyey etal (1986).   
!  This form has not been used much at GFDL, with preference given to 
!  a simpler two-dimensional (depth independent) form assessed with the
!  default hwf_diffusivity_3d=.false. 
!  </DATA> 
!  <DATA NAME="hwf_depth_transition"  UNITS="m"  TYPE="real">
!  Depth of transition for hwf scheme.  The HWF method actually has
!  no depth dependence.  But we include the atan depth dependency 
!  from Bryan-Lewis, for those cases where we wish to replace 
!  Bryan-Lewis with the HWF scheme.  To get the usual Bryan-Lewis
!  transition, set hwf_depth_transition=2500.0.  However, since 
!  we often use hwf in the presence of tide mixing, we do not wish to 
!  have any depth dependence, in which case the default is 
!  hwf_depth_transition=2500.0e4.  
!  </DATA> 
!  <DATA NAME="hwf_min_diffusivity"  UNITS="m^2/sec"  TYPE="real">
!  Minimum diffusivity for the HWF scheme.  
!  Default hwf_min_diffusivity=2e-6.
!  </DATA> 
!  <DATA NAME="hwf_30_diffusivity"  UNITS="m^2/sec"  TYPE="real">
!  Diffusivity at 30deg latitude for the HWF scheme.  
!  Default hwf_30_diffusivity=2e-5.
!  </DATA> 
!  <DATA NAME="hwf_N0_2Omega"  UNITS="dimensionless"  TYPE="real">
!  Ratio of the typical Buoyancy frequency to 
!  twice the earth's rotation period.
!  Default hwf_N0_2Omega=20.0.
!  </DATA> 
!
!  <DATA NAME="bryan_lewis_diffusivity" TYPE="logical">
!  If .true. then add a Bryan-Lewis background to the 
!  diffusivity.  This background is a time-independent function
!  of depth.  This diffusivity is NOT used when have 
!  use_pp_vert_mix_coeff=.true.
!  This scheme should NOT be used if HWF is used.  
!  </DATA> 
!  <DATA NAME="bryan_lewis_lat_depend" TYPE="logical">
!  If .true. then allow for Bryan-Lewis background to be different 
!  outside of a tropical band than inside the band. 
!  </DATA> 
!  <DATA NAME="bryan_lewis_lat_transition" TYPE="real">
!  North/South latitude where transition from Bryan-Lewis values
!  in the tropic to those in the higher latitudes. 
!  </DATA> 
!  <DATA NAME="afkph_90, dfkph_90, sfkph_90, zfkph_90" UNITS="dimensionless" TYPE="real">
!  Parameters setting the Bryan-Lewis vertical diffusivity profile. 
!  When use bryan_lewis_lat_depend, these are the values used in the pole.
!  </DATA> 
!  <DATA NAME="afkph_00, dfkph_00, sfkph_00, zfkph_00" UNITS="dimensionless" TYPE="real">
!  Parameters setting the Bryan-Lewis vertical diffusivity profile in the tropics. 
!  When use bryan_lewis_lat_depend=.true. , these are the values used in the tropics.  
!  When use bryan_lewis_lat_depend=.false., these are the values used globally. 
!  </DATA> 
!
!  <DATA NAME="use_diff_cbt_table" TYPE="logical">
!  If .true., then read in a table that specifies (i,j,ktop-->kbottom) 
!  and the diffusivity. This method is useful when aiming to mix vertically
!  at points where do cross-land insertion or where may wish to enhance 
!  mixing at river mouths.  
!  </DATA> 
!  <DATA NAME="linear_taper_diff_cbt_table" TYPE="logical">
!  If .true., then linear taper the diff_cbt_table value from 
!  so that it gets smaller with depth. 
!  </DATA> 
!
!  <DATA NAME="vmix_rescale_nonbouss" TYPE="logical">
!  To rescale the vertical mixing coefficients by rho0/rho
!  in order to bring the effects from vertical diffusion 
!  in a non-Boussinesq model more in line with that from a 
!  Boussinesq model.
!  Default vmix_rescale_nonbouss=.false. 
!  </DATA> 
!
!  <DATA NAME="vmix_set_min_dissipation" TYPE="logical">
!  To set a minimum dissipation rate.  This scheme will compute the
!  dissipation from the full diffusivity.  If the resulting dissipation
!  is smaller than a specified dissipation, then the diffusivity will 
!  be locally increased so that the min dissipation is maintained.  
!  Default vmix_set_min_dissipation=.false. 
!  </DATA> 
!  <DATA NAME="vmix_min_diss_const"  UNITS="W/m^3"  TYPE="real">
!  Minimum dissipation rate as a constant.  
!  Default vmix_min_diss_const=1e-7. 
!  </DATA> 
!  <DATA NAME="vmix_min_diss_bvfreq_scale" UNITS="J/m^3"  TYPE="real">
!  Scaling use to set the minimum dissipation rate as determined by the 
!  local stratification.  
!  Default vmix_min_diss_bvfreq_scale=6e-4.
!  </DATA> 
!  <DATA NAME="vmix_min_diss_flux_ri_max"  UNITS="dimensionless"  TYPE="real">
!  Maximum flux Richardson number for computation of diffusivity from dissipation.
!  Default vmix_min_diss_flux_ri_max=0.2. 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,       only: epsln, pi, deg_to_rad
use diag_manager_mod,    only: register_diag_field, register_static_field, send_data
use field_manager_mod,   only: MODEL_OCEAN, parse, find_field_index, get_field_methods, method_type, get_field_info
use fms_mod,             only: open_namelist_file, check_nml_error, close_file, write_version_number
use fms_mod,             only: FATAL, WARNING, NOTE, stdout, stdlog
use mpp_domains_mod,     only: mpp_update_domains
use mpp_mod,             only: input_nml_file, mpp_error
use mpp_mod,             only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE

use ocean_domains_mod,         only: get_local_indices
use ocean_parameters_mod,      only: missing_value, rho0, rho0r, omega_earth, onehalf, onefourth, grav 
use ocean_parameters_mod,      only: VERTMIX_CONST, VERTMIX_PP, VERTMIX_CHEN, VERTMIX_GOTM
use ocean_parameters_mod,      only: VERTMIX_KPP_TEST, VERTMIX_KPP_MOM4P0, VERTMIX_KPP_MOM4P1
use ocean_parameters_mod,      only: DEPTH_BASED
use ocean_parameters_mod,      only: MOM_BGRID, MOM_CGRID 
use ocean_tracer_util_mod,     only: tracer_prog_chksum, tracer_min_max
use ocean_types_mod,           only: ocean_grid_type, ocean_domain_type, ocean_thickness_type
use ocean_types_mod,           only: ocean_time_type, ocean_time_steps_type, ocean_options_type
use ocean_types_mod,           only: ocean_velocity_type, ocean_adv_vel_type, ocean_density_type
use ocean_types_mod,           only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_util_mod,            only: invtri, invtri_bmf, write_timestamp
use ocean_util_mod,            only: write_chksum_3d, diagnose_2d, diagnose_3d, diagnose_3d_u, diagnose_2d_u, diagnose_sum, diagnose_2d_en
use ocean_tracer_util_mod,     only: diagnose_3d_rho
use ocean_vert_const_mod,      only: ocean_vert_const_init, vert_mix_const 
use ocean_vert_chen_mod,       only: ocean_vert_chen_init, vert_mix_chen, ocean_vert_chen_end 
use ocean_vert_chen_mod,       only: ocean_vert_chen_restart
use ocean_vert_gotm_mod,       only: ocean_vert_gotm_init, vert_mix_gotm_bgrid, ocean_vert_gotm_end 
use ocean_vert_gotm_mod,       only: ocean_vert_gotm_restart
use ocean_vert_kpp_test_mod,   only: ocean_vert_kpp_test_init, vert_mix_kpp_test 
use ocean_vert_kpp_mom4p0_mod, only: ocean_vert_kpp_mom4p0_init, vert_mix_kpp_mom4p0 
use ocean_vert_kpp_mom4p1_mod, only: ocean_vert_kpp_mom4p1_init, vert_mix_kpp_mom4p1 
use ocean_vert_pp_mod,         only: ocean_vert_pp_init, vert_mix_pp 
use ocean_vert_tidal_mod,      only: ocean_vert_tidal_init, vert_mix_tidal
use ocean_vert_util_mod,       only: ocean_vert_util_init
use ocean_workspace_mod,       only: wrk1, wrk2, wrk3, wrk4, wrk5
use ocean_workspace_mod,       only: wrk1_v, wrk2_v
use ocean_workspace_mod,       only: wrk1_2d, wrk2_2d
use ocean_workspace_mod,       only: wrk1_v2d, wrk2_v2d


implicit none

private

public ocean_vert_mix_init
public ocean_vert_mix_end
public vert_mix_coeff
public vert_friction_bgrid
public vert_friction_cgrid 
public vert_friction_implicit_bgrid
public vert_friction_implicit_cgrid
public vert_diffuse
public vert_diffuse_implicit
public ocean_vert_mix_restart

private vert_friction_init
private bryan_lewis_init
private diff_cbt_tanh_init
private hwf_init
private diff_cbt_table_init
private invcosh
private vmix_min_dissipation
private watermass_diag_init
private watermass_diag
private vert_diffuse_implicit_diag
private vert_diffuse_watermass_diag


! scheme used to compute the vertical diffusivity and vertical viscosity
character(len=10) :: vert_mix_scheme='const' ! "const", "kpp", "kpp_mom4p0", "kpp_mom4p1", "pp", "chen", or "gotm"
integer           :: MIX_SCHEME=1            ! set according to ocean_parameters_mod 

! aidif=1.0 for implicit vertical mixing and aidif=0.0 for explicit
real :: aidif=1.0        

! for Bgrid or Cgrid 
integer :: horz_grid 

! for time steps 
real :: dtime_t
real :: dtime_u
real :: dtime_tr
real :: dtime_ur

! clock ids
integer :: id_clock_vert_const
integer :: id_clock_vert_pp 
integer :: id_clock_vert_kpp_test
integer :: id_clock_vert_kpp_mom4p0
integer :: id_clock_vert_kpp_mom4p1
integer :: id_clock_vert_chen
integer :: id_clock_vert_gotm
integer :: id_clock_vert_tidal 
integer :: id_clock_watermass_diag

! for number of prognostic tracers 
integer :: num_prog_tracers=0

integer :: index_temp
integer :: index_salt

! for global area normalization
real    :: cellarea_r

! for determining how to incorporate the background diffusivity
logical :: vert_diff_back_via_max=.true.

! internally set for computing watermass diagnostics 
logical :: compute_watermass_diag=.false.

! for Bryan-Lewis background vertical diffusivity profile 
logical :: use_explicit_vert_diffuse=.true.             ! to use time-explicit vertical tracer diffusion 
logical :: bryan_lewis_diffusivity=.false.              ! Bryan-Lewis vertical diffusivity 
logical :: bryan_lewis_lat_depend=.false.               ! for different Bryan-Lewis in tropics 
real    :: bryan_lewis_lat_transition=35.0              ! latitude where transition Bryan-Lewis
real    :: afkph_00=0.55,dfkph_00=1.05                  ! Bryan-Lewis parameters for equator (or global)
real    :: sfkph_00=4.5e-5,zfkph_00=2500.0e2            ! Bryan-Lewis parameters for equator (or global)
real    :: afkph_90=0.55,dfkph_90=1.05                  ! Bryan-Lewis parameters for pole
real    :: sfkph_90=4.5e-5,zfkph_90=2500.0e2            ! Bryan-Lewis parameters for pole
real, dimension(:), allocatable :: diff_bryan_lewis_00  ! Bryan-Lewis diffusivity at equator (or global)
real, dimension(:), allocatable :: diff_bryan_lewis_90  ! Bryan-Lewis diffusivity at pole

! tanh vertical profile for background vertical viscosity
real, dimension(:), allocatable ::  visc_cbu_back
real, dimension(:), allocatable ::  visc_cbt_back
logical :: vert_visc_back=.false.
real    :: visc_cbu_back_max=1e-2
real    :: visc_cbu_back_min=1e-3
real    :: visc_cbu_back_zmid=50.0
real    :: visc_cbu_back_zwid=30.0

! tanh vertical profile for background vertical diffusivity
real, dimension(:), allocatable ::  diff_cbt_tanh_back
logical :: diff_cbt_tanh=.false.
real    :: diff_cbt_tanh_max=1e-3
real    :: diff_cbt_tanh_min=2e-5
real    :: diff_cbt_tanh_zmid=150.0
real    :: diff_cbt_tanh_zwid=30.0

! for OBC 
logical :: have_obc=.false.   

! for earth rotation rate squared
real :: omega_earth2 


! for reading in points where wish to set a static background diffusivity 
logical :: use_diff_cbt_table=.false.                     
logical :: linear_taper_diff_cbt_table=.false.
integer, dimension(:), allocatable   :: itable
integer, dimension(:), allocatable   :: jtable     
integer, dimension(:,:), allocatable :: ktable     
real, dimension(:), allocatable      :: diff_cbt_table

! static background vertical diffusivity 
real, dimension(:,:,:), allocatable :: diff_cbt_back  

! diagnosing the vertical diffusion of local reference potrho
real, dimension(:,:,:), allocatable :: wrk1_vmix
real, dimension(:,:,:), allocatable :: wrk2_vmix
real, dimension(:,:,:), allocatable :: diff_cbt_wave
real, dimension(:,:,:), allocatable :: diff_cbt_leewave
real, dimension(:,:,:), allocatable :: diff_cbt_drag

! for background diffusivity that gets smaller in equatorial region
logical :: hwf_diffusivity      = .false.
logical :: hwf_diffusivity_3d   = .false.
real    :: hwf_depth_transition = 2500.0e4
real    :: hwf_min_diffusivity  = 2.e-6
real    :: hwf_30_diffusivity   = 2.e-5
real    :: hwf_N0_2Omega        = 20.0

! to rescale vmix coefficients to bring 
! non-Bouss more in line with Boussinesq
logical :: vmix_rescale_nonbouss=.false.

! for setting a minimum to the dissipation 
logical :: vmix_set_min_dissipation   =.false.
real    :: vmix_min_diss_const        =1e-7
real    :: vmix_min_diss_bvfreq_scale =6e-4
real    :: vmix_min_diss_flux_ri_max  =0.2
logical :: smooth_rho_N2              = .true.   ! for smoothing the rho_N2 field in vertical with 1-2-1
integer :: num_121_passes             = 1        ! number of 1-2-1 passes 

! for diagnostics 
logical :: used
integer :: id_vfrict_expl_u          =-1
integer :: id_vfrict_expl_v          =-1
integer :: id_vfrict_impl_u          =-1
integer :: id_vfrict_impl_v          =-1
integer :: id_visc_cbu_back          =-1
integer :: id_visc_cbt_back          =-1
integer :: id_visc_cbu               =-1
integer :: id_visc_cbt               =-1
integer :: id_visc_cbu_diabatic      =-1
integer :: id_visc_cbt_diabatic      =-1
integer :: id_power_diss_back        =-1 
integer :: id_power_diss             =-1 

integer :: id_visc_cbu_before_min   =-1
integer :: id_diff_cbt_vmix_min     =-1
integer :: id_diff_cbt_t_before_min =-1
integer :: id_diff_cbt_s_before_min =-1
integer :: id_diff_cbt_t            =-1
integer :: id_diff_cbt_s            =-1
integer :: id_diff_bryan_lewis_00   =-1
integer :: id_diff_bryan_lewis_90   =-1
integer :: id_hwf_diffusivity       =-1
integer :: id_diff_cbt_tanh         =-1
integer :: id_diff_cbt_table        =-1
integer :: id_diff_cbt_back         =-1
integer :: id_vmix_min_diss         =-1

integer :: id_wind_power(2)   =-1
integer :: id_bottom_power(2) =-1  

integer :: id_neut_temp_vdiffuse  =-1
integer :: id_neut_salt_vdiffuse  =-1
integer :: id_wdian_temp_vdiffuse =-1
integer :: id_wdian_salt_vdiffuse =-1

integer :: id_neut_rho_vdiffuse          =-1
integer :: id_wdian_rho_vdiffuse         =-1
integer :: id_tform_rho_vdiffuse         =-1
integer :: id_neut_rho_vdiffuse_on_nrho  =-1
integer :: id_wdian_rho_vdiffuse_on_nrho =-1
integer :: id_tform_rho_vdiffuse_on_nrho =-1

integer :: id_neut_rho_vmix          =-1
integer :: id_wdian_rho_vmix         =-1
integer :: id_tform_rho_vmix         =-1
integer :: id_neut_rho_vmix_on_nrho  =-1
integer :: id_wdian_rho_vmix_on_nrho =-1
integer :: id_tform_rho_vmix_on_nrho =-1

integer :: id_eta_tend_k33_tend      =-1
integer :: id_eta_tend_k33_tend_glob =-1

integer :: id_neut_rho_k33          =-1
integer :: id_wdian_rho_k33         =-1
integer :: id_tform_rho_k33         =-1
integer :: id_neut_rho_k33_on_nrho  =-1
integer :: id_wdian_rho_k33_on_nrho =-1
integer :: id_tform_rho_k33_on_nrho =-1

integer :: id_neut_temp_k33          =-1
integer :: id_wdian_temp_k33         =-1
integer :: id_tform_temp_k33         =-1
integer :: id_neut_temp_k33_on_nrho  =-1
integer :: id_wdian_temp_k33_on_nrho =-1
integer :: id_tform_temp_k33_on_nrho =-1

integer :: id_neut_salt_k33          =-1
integer :: id_wdian_salt_k33         =-1
integer :: id_tform_salt_k33         =-1
integer :: id_neut_salt_k33_on_nrho  =-1
integer :: id_wdian_salt_k33_on_nrho =-1
integer :: id_tform_salt_k33_on_nrho =-1

integer :: id_neut_rho_diff_cbt          =-1
integer :: id_wdian_rho_diff_cbt         =-1
integer :: id_tform_rho_diff_cbt         =-1
integer :: id_neut_rho_diff_cbt_on_nrho  =-1
integer :: id_wdian_rho_diff_cbt_on_nrho =-1
integer :: id_tform_rho_diff_cbt_on_nrho =-1

integer :: id_neut_temp_diff_cbt          =-1
integer :: id_wdian_temp_diff_cbt         =-1
integer :: id_tform_temp_diff_cbt         =-1
integer :: id_neut_temp_diff_cbt_on_nrho  =-1
integer :: id_wdian_temp_diff_cbt_on_nrho =-1
integer :: id_tform_temp_diff_cbt_on_nrho =-1

integer :: id_neut_salt_diff_cbt          =-1
integer :: id_wdian_salt_diff_cbt         =-1
integer :: id_tform_salt_diff_cbt         =-1
integer :: id_neut_salt_diff_cbt_on_nrho  =-1
integer :: id_wdian_salt_diff_cbt_on_nrho =-1
integer :: id_tform_salt_diff_cbt_on_nrho =-1

integer :: id_eta_tend_diff_cbt_flx       =-1
integer :: id_eta_tend_diff_cbt_flx_glob  =-1
integer :: id_eta_tend_diff_cbt_tend      =-1
integer :: id_eta_tend_diff_cbt_tend_glob =-1

integer :: id_eta_tend_diff_wave_tend      =-1
integer :: id_eta_tend_diff_wave_tend_glob =-1

integer :: id_neut_rho_diff_wave          =-1
integer :: id_wdian_rho_diff_wave         =-1
integer :: id_tform_rho_diff_wave         =-1
integer :: id_neut_rho_diff_wave_on_nrho  =-1
integer :: id_wdian_rho_diff_wave_on_nrho =-1
integer :: id_tform_rho_diff_wave_on_nrho =-1

integer :: id_neut_temp_diff_wave          =-1
integer :: id_wdian_temp_diff_wave         =-1
integer :: id_tform_temp_diff_wave         =-1
integer :: id_neut_temp_diff_wave_on_nrho  =-1
integer :: id_wdian_temp_diff_wave_on_nrho =-1
integer :: id_tform_temp_diff_wave_on_nrho =-1

integer :: id_neut_salt_diff_wave          =-1
integer :: id_wdian_salt_diff_wave         =-1
integer :: id_tform_salt_diff_wave         =-1
integer :: id_neut_salt_diff_wave_on_nrho  =-1
integer :: id_wdian_salt_diff_wave_on_nrho =-1
integer :: id_tform_salt_diff_wave_on_nrho =-1

integer :: id_eta_tend_diff_drag_tend      =-1
integer :: id_eta_tend_diff_drag_tend_glob =-1

integer :: id_neut_rho_diff_drag          =-1
integer :: id_wdian_rho_diff_drag         =-1
integer :: id_tform_rho_diff_drag         =-1
integer :: id_neut_rho_diff_drag_on_nrho  =-1
integer :: id_wdian_rho_diff_drag_on_nrho =-1
integer :: id_tform_rho_diff_drag_on_nrho =-1

integer :: id_neut_temp_diff_drag          =-1
integer :: id_wdian_temp_diff_drag         =-1
integer :: id_tform_temp_diff_drag         =-1
integer :: id_neut_temp_diff_drag_on_nrho  =-1
integer :: id_wdian_temp_diff_drag_on_nrho =-1
integer :: id_tform_temp_diff_drag_on_nrho =-1

integer :: id_neut_salt_diff_drag          =-1
integer :: id_wdian_salt_diff_drag         =-1
integer :: id_tform_salt_diff_drag         =-1
integer :: id_neut_salt_diff_drag_on_nrho  =-1
integer :: id_wdian_salt_diff_drag_on_nrho =-1
integer :: id_tform_salt_diff_drag_on_nrho =-1


integer :: id_eta_tend_diff_lee_tend     =-1
integer :: id_eta_tend_diff_lee_tend_glob=-1

integer :: id_neut_rho_diff_lee          =-1
integer :: id_wdian_rho_diff_lee         =-1
integer :: id_tform_rho_diff_lee         =-1
integer :: id_neut_rho_diff_lee_on_nrho  =-1
integer :: id_wdian_rho_diff_lee_on_nrho =-1
integer :: id_tform_rho_diff_lee_on_nrho =-1

integer :: id_neut_temp_diff_lee          =-1
integer :: id_wdian_temp_diff_lee         =-1
integer :: id_tform_temp_diff_lee         =-1
integer :: id_neut_temp_diff_lee_on_nrho  =-1
integer :: id_wdian_temp_diff_lee_on_nrho =-1
integer :: id_tform_temp_diff_lee_on_nrho =-1

integer :: id_neut_salt_diff_lee          =-1
integer :: id_wdian_salt_diff_lee         =-1
integer :: id_tform_salt_diff_lee         =-1
integer :: id_neut_salt_diff_lee_on_nrho  =-1
integer :: id_wdian_salt_diff_lee_on_nrho =-1
integer :: id_tform_salt_diff_lee_on_nrho =-1


integer :: id_eta_tend_diff_back_tend      =-1
integer :: id_eta_tend_diff_back_tend_glob =-1

integer :: id_neut_rho_diff_back          =-1
integer :: id_wdian_rho_diff_back         =-1
integer :: id_tform_rho_diff_back         =-1
integer :: id_neut_rho_diff_back_on_nrho  =-1
integer :: id_wdian_rho_diff_back_on_nrho =-1
integer :: id_tform_rho_diff_back_on_nrho =-1

integer :: id_neut_temp_diff_back          =-1
integer :: id_wdian_temp_diff_back         =-1
integer :: id_tform_temp_diff_back         =-1
integer :: id_neut_temp_diff_back_on_nrho  =-1
integer :: id_wdian_temp_diff_back_on_nrho =-1
integer :: id_tform_temp_diff_back_on_nrho =-1

integer :: id_neut_salt_diff_back          =-1
integer :: id_wdian_salt_diff_back         =-1
integer :: id_tform_salt_diff_back         =-1
integer :: id_neut_salt_diff_back_on_nrho  =-1
integer :: id_wdian_salt_diff_back_on_nrho =-1
integer :: id_tform_salt_diff_back_on_nrho =-1


integer  :: id_neut_rho_sbc              =-1
integer  :: id_neut_rho_sbc_temp         =-1
integer  :: id_neut_rho_sbc_salt         =-1
integer  :: id_neut_rho_sbc_on_nrho      =-1
integer  :: id_neut_rho_sbc_temp_on_nrho =-1
integer  :: id_neut_rho_sbc_salt_on_nrho =-1

integer  :: id_wdian_rho_sbc_temp         =-1
integer  :: id_wdian_rho_sbc_salt         =-1
integer  :: id_wdian_rho_sbc              =-1
integer  :: id_wdian_rho_sbc_temp_on_nrho =-1
integer  :: id_wdian_rho_sbc_salt_on_nrho =-1
integer  :: id_wdian_rho_sbc_on_nrho      =-1

integer  :: id_tform_rho_sbc              =-1
integer  :: id_tform_rho_sbc_temp         =-1
integer  :: id_tform_rho_sbc_salt         =-1
integer  :: id_tform_rho_sbc_on_nrho      =-1
integer  :: id_tform_rho_sbc_temp_on_nrho =-1
integer  :: id_tform_rho_sbc_salt_on_nrho =-1

integer  :: id_neut_rho_bbc_temp          =-1
integer  :: id_wdian_rho_bbc_temp         =-1
integer  :: id_tform_rho_bbc_temp         =-1
integer  :: id_neut_rho_bbc_temp_on_nrho  =-1
integer  :: id_wdian_rho_bbc_temp_on_nrho =-1
integer  :: id_tform_rho_bbc_temp_on_nrho =-1

integer, dimension(:), allocatable :: id_zflux_diff
integer, dimension(:), allocatable :: id_vdiffuse
integer, dimension(:), allocatable :: id_vdiffuse_impl
integer, dimension(:), allocatable :: id_vdiffuse_k33
integer, dimension(:), allocatable :: id_vdiffuse_diff_cbt
integer, dimension(:), allocatable :: id_vdiffuse_sbc
integer, dimension(:), allocatable :: id_vdiffuse_bbc
integer, dimension(:), allocatable :: id_vdiffuse_diss

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

real, dimension(:,:,:), allocatable :: flux_z

character(len=128) :: version = &
     '$Id: ocean_vert_mix.F90,v 20.0 2013/12/14 00:16:48 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

logical :: module_is_initialized   = .false.
logical :: debug_this_module       = .false. 
logical :: verbose_init            = .true.
logical :: quebec_2009_10_bug      = .false. 

namelist /ocean_vert_mix_nml/ debug_this_module, vert_mix_scheme, verbose_init, aidif,                              &
 vert_diff_back_via_max, use_explicit_vert_diffuse, use_diff_cbt_table, linear_taper_diff_cbt_table,                &
 bryan_lewis_diffusivity, bryan_lewis_lat_depend, bryan_lewis_lat_transition,                                       &
 afkph_90, dfkph_90, sfkph_90, zfkph_90, afkph_00, dfkph_00, sfkph_00, zfkph_00,                                    &
 vert_visc_back, visc_cbu_back_max, visc_cbu_back_min, visc_cbu_back_zmid, visc_cbu_back_zwid,                      &
 hwf_diffusivity, hwf_depth_transition, hwf_min_diffusivity, hwf_30_diffusivity, hwf_N0_2Omega, hwf_diffusivity_3d, &
 diff_cbt_tanh, diff_cbt_tanh_max, diff_cbt_tanh_min, diff_cbt_tanh_zmid, diff_cbt_tanh_zwid,                       &
 quebec_2009_10_bug, vmix_rescale_nonbouss,                                                                         &
 vmix_set_min_dissipation, vmix_min_diss_const, vmix_min_diss_bvfreq_scale, vmix_min_diss_flux_ri_max,              &
 smooth_rho_N2, num_121_passes

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_vert_mix_init">
!
! <DESCRIPTION>
! Initialization for the vertical mixing module.
! </DESCRIPTION>
!
subroutine ocean_vert_mix_init (Grid, Domain, Time, Dens, Velocity, Time_steps, Ocean_options, T_prog, T_diag, &
                                vert_coordinate, vert_coordinate_class, hor_grid, obc, debug)

  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_density_type),     intent(in)           :: Dens
  type(ocean_velocity_type),    intent(in)           :: Velocity
  type(ocean_time_steps_type),  intent(inout)        :: Time_steps
  type(ocean_options_type),     intent(inout)        :: Ocean_options
  type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
  type(ocean_diag_tracer_type), intent(in)           :: T_diag(:)
  integer,                      intent(in)           :: vert_coordinate
  integer,                      intent(in)           :: vert_coordinate_class
  integer,                      intent(in)           :: hor_grid 
  logical,                      intent(in)           :: obc
  logical,                      intent(in), optional :: debug 

  integer :: ioun, io_status, ierr
  integer :: n
  integer :: num_methods=0 

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_vert_mix_mod (ocean_vert_mix_init)" module already initialized')
  endif 

  module_is_initialized = .TRUE.

  num_prog_tracers = size(T_prog(:))
  do n=1,num_prog_tracers
     if (trim(T_prog(n)%name) == 'temp') index_temp = n
     if (trim(T_prog(n)%name) == 'salt') index_salt = n
  enddo

  call write_version_number( version, tagname )

#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_vert_mix_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_vert_mix_nml')
#else
  ioun = open_namelist_file()
  read  (ioun, ocean_vert_mix_nml,iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_vert_mix_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_vert_mix_nml)  
  write (stdlogunit, ocean_vert_mix_nml)

  if (PRESENT(debug) .and. .not. debug_this_module) then 
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') &
    '==>Note: running ocean_vert_mix_mod with debug_this_module=.true.'  
  endif 

  Dom => Domain
  Grd => Grid
  Time_steps%aidif = aidif
  dtime_t          = Time_steps%dtime_t
  dtime_u          = Time_steps%dtime_u
  dtime_tr         = 1.0/dtime_t
  dtime_ur         = 1.0/dtime_u
  have_obc         = obc
  omega_earth2     = omega_earth**2
  cellarea_r       = 1.0/(epsln + Grd%tcellsurf)
  horz_grid        = hor_grid 

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  if(0.0 < aidif .and. aidif < 1.0) then 
      call mpp_error(FATAL, &
      '==>Error: ocean_vert_mix_mod: set aidif==0.0 OR aidif==1.0. Other values unsupported.')  
  endif

  if(vmix_rescale_nonbouss) then 
      if(vert_coordinate_class==DEPTH_BASED) then 
          call mpp_error(WARNING, &
               '==>Error: ocean_vert_mix_mod: Rescaling of vmix is not intended for depth-based simulations.')  
      else
          write(stdoutunit,'(a)') &
               '==>Note: ocean_vert_mix_mod: rescaling vertical mixing coefficients by rho0/rho.'  
      endif
  endif

  allocate (diff_cbt_back(isd:ied,jsd:jed,nk))
  diff_cbt_back(:,:,:) = 0.0 

  allocate(flux_z(isd:ied,jsd:jed,nk))
  flux_z = 0.0

  allocate (wrk1_vmix(isd:ied,jsd:jed,nk))
  allocate (wrk2_vmix(isd:ied,jsd:jed,nk))
  wrk1_vmix(:,:,:) = 0.0 
  wrk2_vmix(:,:,:) = 0.0 

  allocate (diff_cbt_wave(isd:ied,jsd:jed,nk))
  allocate (diff_cbt_leewave(isd:ied,jsd:jed,nk))
  allocate (diff_cbt_drag(isd:ied,jsd:jed,nk))
  diff_cbt_wave(:,:,:)    = 0.0 
  diff_cbt_leewave(:,:,:) = 0.0 
  diff_cbt_drag(:,:,:)    = 0.0 

  ! initialize clock ids 
  id_clock_vert_const       = mpp_clock_id('(Ocean vmix: constant)'        ,grain=CLOCK_ROUTINE)
  id_clock_vert_pp          = mpp_clock_id('(Ocean vmix: PPvmix)'          ,grain=CLOCK_ROUTINE)
  id_clock_vert_kpp_test    = mpp_clock_id('(Ocean vmix: KPP_test)'        ,grain=CLOCK_ROUTINE)
  id_clock_vert_kpp_mom4p0  = mpp_clock_id('(Ocean vmix: KPP_mom4p0)'      ,grain=CLOCK_ROUTINE)
  id_clock_vert_kpp_mom4p1  = mpp_clock_id('(Ocean vmix: KPP_mom4p1)'      ,grain=CLOCK_ROUTINE)
  id_clock_vert_chen        = mpp_clock_id('(Ocean vmix: chenvmix)'        ,grain=CLOCK_ROUTINE)
  id_clock_vert_gotm        = mpp_clock_id('(Ocean vmix: gotmvmix)'        ,grain=CLOCK_ROUTINE)
  id_clock_vert_tidal       = mpp_clock_id('(Ocean vmix: tidal)'           ,grain=CLOCK_ROUTINE)
  id_clock_watermass_diag   = mpp_clock_id('(Ocean vmix: watermass_diag)'  ,grain=CLOCK_ROUTINE)

  ! initialize the vertical mixing utility
  call ocean_vert_util_init(Grid, Domain, Time)

  ! read in static tracer diffusivities from a table and add to diff_cbt_back
  call diff_cbt_table_init(Time)

  ! initialize static Bryan-Lewis background diffusivity and add to diff_cbt_back
  call bryan_lewis_init(Time, Ocean_options)

  ! initialize static HWF background diffusivity and add to diff_cbt_back
  call hwf_init(Time, Ocean_options)

  ! initialize static tanh background diffusivity and add to diff_cbt_back
  call diff_cbt_tanh_init(Time, Ocean_options)

  ! initialize module used to compute diffusivity based on tidal dissipation
  call ocean_vert_tidal_init(Grid, Domain, Time, T_prog, Velocity, Ocean_options, dtime_t, vert_mix_scheme, horz_grid)

  ! initialize vertical viscosity 
  call vert_friction_init(Time)

  ! initialize "boundary layer" diffusivity and viscosity schemes 
  if(vert_mix_scheme == 'const') then 
    num_methods=num_methods+1
    write(stdoutunit,'(/a/)') &
     '==>Note from ocean_vert_mix: constant scheme used for vertical diffusivities and viscosities'
    MIX_SCHEME = VERTMIX_CONST
    Ocean_options%vert_mix = &
    'Constant scheme used for vertical diffusivities & viscosity.'
    Ocean_options%vertmix  = VERTMIX_CONST
    call ocean_vert_const_init (Grid, Domain, Time, Time_steps, T_prog)

  elseif(vert_mix_scheme == 'kpp_mom4p0') then
    if(horz_grid == MOM_CGRID) then 
       write(stdoutunit,'(a)') &
       '==>Error: ocean_vert_mix_mod: kpp_mom4p0 is only for B-grid. Use kpp_test for KPP + C-grid.'
       call mpp_error(FATAL, &
       '==>Error: ocean_vert_mix_mod: kpp_mom4p0 is only for B-grid. Use kpp_test for KPP + C-grid.')
    endif 
    num_methods=num_methods+1
    write(stdoutunit,'(/a/)') &
     '==>Note from ocean_vert_mix: KPP_mom4p0 for vert diffusivity, viscosity, nonlocal, and barotropic tide drag.'
    MIX_SCHEME = VERTMIX_KPP_MOM4P0
    Ocean_options%vert_mix = 'KPP_mom4p0 used for vertical diffusivities & viscosity & barotropic tide drag.'
    Ocean_options%vertmix  = VERTMIX_KPP_MOM4P0
    call ocean_vert_kpp_mom4p0_init (Grid, Domain, Time, Time_steps, Dens, T_prog, T_diag, vert_coordinate)

  elseif(vert_mix_scheme == 'kpp_mom4p1') then
    if(horz_grid == MOM_CGRID) then 
       write(stdoutunit,'(a)') &
       '==>Error: ocean_vert_mix_mod: kpp_mom4p1 is only for B-grid. Use kpp_test for KPP + C-grid.'
       call mpp_error(FATAL, &
       '==>Error: ocean_vert_mix_mod: kpp_mom4p1 is only for B-grid. Use kpp_test for KPP + C-grid.')
    endif 
    num_methods=num_methods+1
    write(stdoutunit,'(/a/)') &
     '==>Note from ocean_vert_mix: KPP_mom4p1 for vert diffusivity, viscosity, nonlocal, and barotropic tide drag.'
    MIX_SCHEME = VERTMIX_KPP_MOM4P1
    Ocean_options%vert_mix = 'KPP_mom4p1 used for vertical diffusivities & viscosity & barotropic tide drag.'
    Ocean_options%vertmix  = VERTMIX_KPP_MOM4P1
    call ocean_vert_kpp_mom4p1_init (Grid, Domain, Time, Time_steps, Dens, T_prog, T_diag, vert_coordinate, vert_coordinate_class)

  elseif(vert_mix_scheme == 'kpp_test') then
    num_methods=num_methods+1
    write(stdoutunit,'(/a/)') &
     '==>Note from ocean_vert_mix: KPP-test used for vert diffusivity, viscosity, and nonlocal transport.'
    MIX_SCHEME = VERTMIX_KPP_TEST
    Ocean_options%vert_mix = 'KPP-test used for vertical diffusivities & viscosity.'
    Ocean_options%vertmix  = VERTMIX_KPP_TEST
    call ocean_vert_kpp_test_init (Grid, Domain, Time, Time_steps, Dens, T_prog, T_diag, &
                                   vert_coordinate, vert_coordinate_class, horz_grid)

  elseif(vert_mix_scheme == 'pp') then
    num_methods=num_methods+1
    write(stdoutunit,'(/a/)') &
     '==>Note from ocean_vert_mix: Pacanowski-Philander for vert diffusivity and viscosity'
    MIX_SCHEME = VERTMIX_PP
    Ocean_options%vert_mix = 'PP used for vertical diffusivities & viscosity.'
    Ocean_options%vertmix  = VERTMIX_PP
    call ocean_vert_pp_init (Grid, Domain, Time, Time_steps, T_prog, horz_grid)

  elseif(vert_mix_scheme == 'chen') then
    num_methods=num_methods+1
    write(stdoutunit,'(/a/)') &
     '==>Note from ocean_vert_mix: Chen for vertical diffusivities and vertical viscosities'
    Ocean_options%vert_mix = 'Chen used for vertical diffusivities & viscosity.'
    Ocean_options%vertmix  = VERTMIX_CHEN 
    MIX_SCHEME = VERTMIX_CHEN 
    call ocean_vert_chen_init (Grid, Domain, Time, Time_steps, T_prog, horz_grid, debug)

  elseif(vert_mix_scheme == 'gotm') then
    num_methods=num_methods+1
    write(stdoutunit,'(/a/)') &
     '==>Note from ocean_vert_mix: GOTM for vertical diffusivities and vertical viscosities'
    write(stdoutunit,'(a)') ' '
    Ocean_options%vert_mix = 'GOTM used for vertical diffusivities & viscosity.'
    Ocean_options%vertmix  = VERTMIX_GOTM 
    MIX_SCHEME = VERTMIX_GOTM
    call ocean_vert_gotm_init (Grid, Domain, Time, Time_steps, T_prog, have_obc, debug)
  endif 

  if(num_methods==0) then 
    write(stdoutunit,'(a)') &
    '==>Error: ocean_vert_mix_mod: choose vert mixing: "const", "kpp_test", "kpp_mom4p0", "kpp_mom4p1", "pp", "chen", or "gotm".'
    call mpp_error(FATAL, &
    '==>Error: ocean_vert_mix_mod: choose vert mixing: "const", "kpp_test", "kpp_mom4p0", "kpp_mom4p1", "pp", "chen", or "gotm".')
  elseif(num_methods>1) then 
    write(stdoutunit,'(a)')'==>Error: ocean_vert_mix_mod: choose only one vertical mixing scheme.'
    call mpp_error(FATAL,'==>Error: ocean_vert_mix_mod: choose only one vertical mixing scheme.')
  endif 

  if(MIX_SCHEME == VERTMIX_GOTM .and. horz_grid == MOM_CGRID) then 
    write(stdoutunit,'(a)')'==>Warning: ocean_vert_mix_mod: GOTM has yet to be updated for Cgrid. Will use bgrid routine anyhow.'
    write(stdoutunit,'(a)')'                                Vertical diffusivities will be slightly "wrong".' 
  endif 

  ! register thickness weighted vertical diffusion tendency for diagnostic output 
  allocate( id_zflux_diff(num_prog_tracers) )
  allocate( id_vdiffuse(num_prog_tracers) )
  allocate( id_vdiffuse_impl(num_prog_tracers) )
  allocate( id_vdiffuse_k33(num_prog_tracers) )
  allocate( id_vdiffuse_diff_cbt(num_prog_tracers) )
  allocate( id_vdiffuse_sbc(num_prog_tracers) )
  allocate( id_vdiffuse_bbc(num_prog_tracers) )
  allocate( id_vdiffuse_diss(num_prog_tracers) )
  id_zflux_diff       =-1
  id_vdiffuse         =-1
  id_vdiffuse_impl    =-1
  id_vdiffuse_k33     =-1
  id_vdiffuse_diff_cbt=-1
  id_vdiffuse_sbc     =-1
  id_vdiffuse_bbc     =-1
  id_vdiffuse_diss    =-1

  do n=1,num_prog_tracers

     if (trim(T_prog(n)%name) == 'temp') then
         id_zflux_diff(n) = register_diag_field ('ocean_model',          &
              trim(T_prog(n)%name)//'_zflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'z-diffusive heat flux', 'Watts/m^2',     &
              missing_value=missing_value, range=(/-1.e16,1.e16/))
         id_vdiffuse(n) = register_diag_field ('ocean_model',                 &
              trim(T_prog(n)%name)//'_vdiffuse', Grd%tracer_axes(1:3),        &
              Time%model_time, 'explicit vert diffusion of heat', 'Watts/m^2',&
              missing_value=missing_value, range=(/-1.e16,1.e16/)) 
         id_vdiffuse_impl(n) = register_diag_field ('ocean_model',            &
              trim(T_prog(n)%name)//'_vdiffuse_impl', Grd%tracer_axes(1:3),   &
              Time%model_time, 'implicit vert diffusion of heat', 'Watts/m^2',&
              missing_value=missing_value, range=(/-1.e16,1.e16/)) 
         id_vdiffuse_k33(n) = register_diag_field ('ocean_model',                                      &
              trim(T_prog(n)%name)//'_vdiffuse_k33', Grd%tracer_axes(1:3),                             &
              Time%model_time, 'vert diffusion of heat due to K33 from neutral diffusion', 'Watts/m^2',&
              missing_value=missing_value, range=(/-1.e16,1.e16/)) 
         id_vdiffuse_diff_cbt(n) = register_diag_field ('ocean_model',               &
              trim(T_prog(n)%name)//'_vdiffuse_diff_cbt', Grd%tracer_axes(1:3),      &
              Time%model_time, 'vert diffusion of heat due to diff_cbt', 'Watts/m^2',&
              missing_value=missing_value, range=(/-1.e16,1.e16/)) 
         id_vdiffuse_sbc(n) = register_diag_field ('ocean_model',                        &
              trim(T_prog(n)%name)//'_vdiffuse_sbc', Grd%tracer_axes(1:3),               &
              Time%model_time, 'vert diffusion of heat due to surface flux', 'Watts/m^2',&
              missing_value=missing_value, range=(/-1.e16,1.e16/)) 
         id_vdiffuse_bbc(n) = register_diag_field ('ocean_model',                       &
              trim(T_prog(n)%name)//'_vdiffuse_bbc', Grd%tracer_axes(1:3),              &
              Time%model_time, 'vert diffusion of heat due to bottom flux', 'Watts/m^2',&
              missing_value=missing_value, range=(/-1.e16,1.e16/)) 
         id_vdiffuse_diss(n) = register_diag_field ('ocean_model',                               &
              trim(T_prog(n)%name)//'_vdiffuse_diss', Grd%tracer_axes(1:3),                      &
              Time%model_time, 'dissipation of squared temp via vert diffusion', '(Watts/m^2)^2',&
              missing_value=missing_value, range=(/-1.e20,1.e20/)) 
     else 
         id_zflux_diff(n) = register_diag_field ('ocean_model',&
         trim(T_prog(n)%name)//'_zflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'k-diffusive flux of '//trim(T_prog(n)%longname), &
              trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1e20,1e20/))
         id_vdiffuse(n) = register_diag_field ('ocean_model',&
         trim(T_prog(n)%name)//'_vdiffuse', Grd%tracer_axes(1:3), &
              Time%model_time, 'explicit vert diffusion of '//trim(T_prog(n)%longname), &
              trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1e20,1e20/))         
         id_vdiffuse_impl(n) = register_diag_field ('ocean_model',&
         trim(T_prog(n)%name)//'_vdiffuse_impl', Grd%tracer_axes(1:3), &
              Time%model_time, 'implicit vert diffusion of '//trim(T_prog(n)%longname), &
              trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1e20,1e20/))         
         id_vdiffuse_k33(n) = register_diag_field ('ocean_model',                                &
         trim(T_prog(n)%name)//'_vdiffuse_k33', Grd%tracer_axes(1:3), Time%model_time,           &
              'vert diffusion due to K33 from neutral diffusion for '//trim(T_prog(n)%longname), &
              trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1e20,1e20/))         
         id_vdiffuse_diff_cbt(n) = register_diag_field ('ocean_model',                        &
         trim(T_prog(n)%name)//'_vdiffuse_diff_cbt', Grd%tracer_axes(1:3), Time%model_time,   &
              'vert diffusion due to diff_cbt for '//trim(T_prog(n)%longname),                &
              trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1e20,1e20/))         
         id_vdiffuse_sbc(n) = register_diag_field ('ocean_model',                     &
         trim(T_prog(n)%name)//'_vdiffuse_sbc', Grd%tracer_axes(1:3), Time%model_time,&
              'vert diffusion due to surface flux for '//trim(T_prog(n)%longname),    &
              trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1e20,1e20/))         
         id_vdiffuse_bbc(n) = register_diag_field ('ocean_model',                     &
         trim(T_prog(n)%name)//'_vdiffuse_bbc', Grd%tracer_axes(1:3), Time%model_time,&
              'vert diffusion due to bottom flux for '//trim(T_prog(n)%longname),     &
              trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1e20,1e20/))         
         id_vdiffuse_diss(n) = register_diag_field ('ocean_model',                              &
         trim(T_prog(n)%name)//'_vdiffuse_diss', Grd%tracer_axes(1:3), Time%model_time,         &
              'dissipation of squared tracer via vert diffusion for '//trim(T_prog(n)%longname),&
              '[kg/(m^2*sec)]^2', missing_value=missing_value, range=(/-1e20,1e20/))         
    endif  

  enddo
  
  ! register and send static background diffusivity 
  id_diff_cbt_back = -1
  id_diff_cbt_back = register_static_field ('ocean_model', 'diff_cbt_back', Grid%tracer_axes(1:3), &
       'static background vertical diffusivity diff_cbt', 'm^2/s',missing_value=missing_value,     &
       range=(/-10.0,1e6/),                                                                        &
       standard_name='ocean_vertical_tracer_diffusivity_due_to_background')
  call diagnose_3d(Time, Grd, id_diff_cbt_back, diff_cbt_back(:,:,:))

  ! power dissipated by background mixing acting against stratification 
    id_power_diss_back = register_diag_field ('ocean_model', 'power_diss_back',                                 &
    Grid%tracer_axes(1:3), Time%model_time, 'power dissipation from prescribed background vertical diffusivity',&
    'W/m^2', missing_value=missing_value, range=(/-1e15,1e15/),                                                 &
    standard_name='tendency_of_ocean_potential_energy_content_due_to_background')

  ! power dissipated by total vertical mixing acting against temp and salinity stratification 
    id_power_diss = register_diag_field ('ocean_model', 'power_diss',                                 &
    Grid%tracer_axes(1:3), Time%model_time, 'power dissipation from vertical heat and salt diffusion',&
    'W/m^2', missing_value=missing_value, range=(/-1e15,1e15/),                                       &
    standard_name='tendency_of_ocean_potential_energy_content')

  ! register vertical diffusivity
  id_diff_cbt_vmix_min = -1
  id_diff_cbt_vmix_min = register_diag_field ('ocean_model', 'diff_cbt_vmix_min',&
                     Grid%tracer_axes(1:3), Time%model_time,                     &
                     'floored diffusivity based on minimum dissipation',         &
                     'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/))

  id_diff_cbt_t_before_min = -1
  id_diff_cbt_t_before_min = register_diag_field ('ocean_model', 'diff_cbt_t_before_min',&
                     Grid%tracer_axes(1:3), Time%model_time,                             &
                     'vert diff_cbt(temp) before apply min dissipation step',            &
                     'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/))

  id_diff_cbt_s_before_min = -1
  id_diff_cbt_s_before_min = register_diag_field ('ocean_model', 'diff_cbt_s_before_min',&
                     Grid%tracer_axes(1:3), Time%model_time,                             &
                     'vert diff_cbt(salt) before apply min dissipation step',            &
                     'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/))

  id_diff_cbt_t = -1
  id_diff_cbt_t    = register_diag_field ('ocean_model', 'diff_cbt_t', Grid%tracer_axes(1:3), &
                     Time%model_time, 'total vert diff_cbt(temp) (w/o neutral included)',     &
                     'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/),                &
                     standard_name='ocean_vertical_heat_diffusivity')

  id_diff_cbt_s = -1
  id_diff_cbt_s    = register_diag_field ('ocean_model', 'diff_cbt_s', Grid%tracer_axes(1:3), &
                     Time%model_time, 'total vert diff_cbt(salt) (w/o neutral included)',     &
                     'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/),                &
                     standard_name='ocean_vertical_salt_diffusivity')

  id_vmix_min_diss = -1
  id_vmix_min_diss = register_diag_field ('ocean_model', 'vmix_min_diss',&
                     Grid%tracer_axes(1:3), Time%model_time,             &
                     'floored dissipation from vmix_min scheme',         &
                     'W/m^3',missing_value=missing_value, range=(/-10.0,1e12/))

  ! register non-Bouss steric sea level tendency arising from vertical diffusion 
  id_eta_tend_diff_cbt_flx = register_diag_field ('ocean_model',               &
          'eta_tend_diff_cbt_flx', Grid%tracer_axes(1:2), Time%model_time,     &
          'non-Bouss steric sea level tendency from diff_cbt diffusive fluxes',&
          'm/s',missing_value=missing_value, range=(/-1e20,1e20/))
  id_eta_tend_diff_cbt_flx_glob = register_diag_field ('ocean_model',                      &
          'eta_tend_diff_cbt_flx_glob', Time%model_time,                                   &
          'global mean non-Bouss steric sea level tendency from diff_cbt diffusive fluxes',&
          'm/s',missing_value=missing_value, range=(/-1e20,1e20/))

  call watermass_diag_init(Time, Dens)

end subroutine ocean_vert_mix_init
! </SUBROUTINE> NAME="ocean_vert_mix_init"


!#######################################################################
! <SUBROUTINE NAME="diff_cbt_table_init">
!
! <DESCRIPTION>
! Read in static diffusivities that have been entered to the diff_cbt_table.
! </DESCRIPTION>
!
subroutine diff_cbt_table_init(Time)
  type(ocean_time_type), intent(in) :: Time

  type(method_type), allocatable, dimension(:) :: diff_cbt_table_methods

  character(len=32) :: fld_type, fld_name
  real    :: taper_table, thickness_table  
  integer :: ntable, ntable_points, model, parse_ok
  integer :: itbl, jtbl, ntbl
  integer :: i,j,k

  integer :: stdoutunit 
  stdoutunit=stdout() 

  ntable = find_field_index(MODEL_OCEAN,'diff_cbt_enhance')
  if (ntable < 1) then 
    call mpp_error(NOTE, &
    '==>Warning: ocean_vert_mix_init NO table for enhanced diff_cbt. No diffusivities read.')  
  endif
  if(ntable > 1 .and. .not. use_diff_cbt_table) then 
    call mpp_error(NOTE, &
    '==>Warning: ocean_vert_mix_init found diff_cbt_table, yet use_diff_cbt_table=.false.')  
  endif 
  if(ntable > 1 .and. use_diff_cbt_table) then 
    write(stdoutunit,*) &
    '==>Note: ocean_vert_mix_init will read background diffusivities from diff_cbt_table.'  
    if(linear_taper_diff_cbt_table) then 
       write(stdoutunit,*) &
       '==>Note: ocean_vert_mix_init will linearly taper diff_cbt_table values to zero with depth.'
    endif
  endif 


  ! read in static diffusivities specified from the field table "diff_cbt_table"
  wrk2(:,:,:) = 0.0 
  if(use_diff_cbt_table .and. ntable > 1) then 

      call get_field_info(ntable,fld_type,fld_name,model,ntable_points)
      allocate(diff_cbt_table_methods(ntable_points))

      allocate (itable(ntable_points))
      itable(:) = 0
      allocate (jtable(ntable_points))
      jtable(:) = 0
      allocate (ktable(ntable_points,2))
      ktable(:,:) = 0
      allocate (diff_cbt_table(ntable_points))
      diff_cbt_table(:) = 0.0

      call get_field_methods(ntable,diff_cbt_table_methods)
      do ntbl=1,ntable_points
         parse_ok = parse(diff_cbt_table_methods(ntbl)%method_control,'itable',itable(ntbl))
         if (parse_ok == 0) call mpp_error(FATAL,&
         '==>Error ocean_vert_mix_init: diff_cbt_table entry "itable" error')
         parse_ok = parse(diff_cbt_table_methods(ntbl)%method_control,'jtable',jtable(ntbl))
         if (parse_ok == 0) call mpp_error(FATAL,&
         '==>Error ocean_vert_mix_init: diff_cbt_table entry "jtable" error')
         parse_ok = parse(diff_cbt_table_methods(ntbl)%method_control,'ktable_1',ktable(ntbl,1))
         if (parse_ok == 0) call mpp_error(FATAL,&
         '==>Error ocean_vert_mix_init: diff_cbt_table entry "ktable_1" error')
         parse_ok = parse(diff_cbt_table_methods(ntbl)%method_control,'ktable_2',ktable(ntbl,2))
         if (parse_ok == 0) call mpp_error(FATAL,&
         '==>Error ocean_vert_mix_init: diff_cbt_table entry "ktable_2" error')
         parse_ok = parse(diff_cbt_table_methods(ntbl)%method_control,'diff_cbt_table',diff_cbt_table(ntbl))
         if (parse_ok == 0) call mpp_error(FATAL,&
         '==>Error ocean_vert_mix_init: diff_cbt_table entry "diff_cbt_table" error')
         if(ktable(ntbl,2) < ktable(ntbl,1) ) then 
             call mpp_error(FATAL,&
             '==>Error ocean_vert_mix_init: diff_cbt_table entry needs ktable_2 >= ktable_1')
         endif
      enddo

      do ntbl=1,ntable_points

         if(on_comp_domain(ntbl)) then

             itbl=itable(ntbl)-Dom%ioff
             jtbl=jtable(ntbl)-Dom%joff

             if(linear_taper_diff_cbt_table) then  
                 thickness_table = Grd%zw(ktable(ntbl,2)) - Grd%zw(ktable(ntbl,1))
             endif

             do k=ktable(ntbl,1),ktable(ntbl,2)

                if(Grd%tmask(itbl,jtbl,k) == 0.0) then 
                    write(*,'(a,i4,a,i4,a,i4,a)')                                      &
                    '==>Warning ocean_vert_mix_init: ignored nonzero diff_cbt_table(', &
                    itable(ntbl),',',jtable(ntbl),',',k,') set over land. Is your table correct?'
                endif

                if(linear_taper_diff_cbt_table .and. thickness_table > 0.0) then              
                    taper_table = (Grd%zw(ktable(ntbl,2)) -Grd%zw(k))/thickness_table 
                else 
                    taper_table = 1.0
                endif

                wrk2(itbl,jtbl,k)          = taper_table*Grd%tmask(itbl,jtbl,k)*diff_cbt_table(ntbl)  
                diff_cbt_back(itbl,jtbl,k) = diff_cbt_back(itbl,jtbl,k) + wrk2(itbl,jtbl,k)

             enddo

         endif

      enddo

  else 

      ntable_points=1 
      allocate (itable(ntable_points))
      itable(:) = 0
      allocate (jtable(ntable_points))
      itable(:) = 0
      allocate (ktable(ntable_points,2))
      ktable(:,:) = 0
      allocate (diff_cbt_table(ntable_points))
      diff_cbt_table(:) = 0.0

  endif

  k=1
  do j=jsc,jec
     do i=isc,iec
        if ((dtime_t*diff_cbt_back(i,j,k))/Grd%dzt(k)**2 >= 0.5  .and. aidif /= 1.0) then
            write(stdoutunit,'(a,a,i3,a,i3,a,i3,a,f12.3)')                &
                 '==> Warning: vertical diffusive criteria exceeded w/',&
                 'diff_cbt_table(',i+Dom%ioff,',',j+Dom%joff,',',k,') =',diff_cbt_back(i,j,k)
        endif
     enddo
  enddo

  id_diff_cbt_table = -1
  id_diff_cbt_table = register_static_field ('ocean_model', 'diff_cbt_table', &
                      Grd%tracer_axes(1:3), 'diff_cbt from table', 'm^2/s',   &
                      missing_value=missing_value, range=(/-10.0,1e6/))
  call diagnose_3d(Time, Grd, id_diff_cbt_table, wrk2(:,:,:))

end subroutine diff_cbt_table_init
! </SUBROUTINE> NAME="diff_cbt_table_init"



!#######################################################################
! <SUBROUTINE NAME="bryan_lewis_init">
!
! <DESCRIPTION>
! Initialize the Bryan-Lewis static background diffusivities.
! </DESCRIPTION>
!
subroutine bryan_lewis_init(Time, Ocean_options)
  type(ocean_time_type),    intent(in)    :: Time
  type(ocean_options_type), intent(inout) :: Ocean_options

  integer :: i,j,k
  real    :: weight 

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(bryan_lewis_diffusivity) then
     Ocean_options%bryan_lewis_mix = 'Used Bryan-Lewis background vertical diffusivity.'
  else
     Ocean_options%bryan_lewis_mix = 'Did NOT use Bryan-Lewis background vertical diffusivity.'
     return 
  endif 

  if(hwf_diffusivity) then
    call mpp_error(WARNING,&
    '==>Warning from ocean_vert_mix_mod: Bryan-Lewis and HWF should NOT both be enabled, but they are.')
  endif 

  allocate (diff_bryan_lewis_00(nk))
  diff_bryan_lewis_00(:) = 0.0 
  allocate (diff_bryan_lewis_90(nk))
  diff_bryan_lewis_90(:) = 0.0 

  write(stdoutunit,*) &
  '=>Note: USING Bryan Lewis vertical background tracer diffusivity.' 
  write(stdoutunit,*) &
  '        If using diffusivity from vert_tidal_mod, then should turn off Bryan Lewis.' 
  if(bryan_lewis_lat_depend) then 
      write(stdoutunit,*) &
      '=>Note: USING bryan_lewis_lat_depend.  This will modify the background latitudinally, '
      write(stdoutunit,*)'        with a transition latitude at bryan_lewis_lat_transition = ',&
                                bryan_lewis_lat_transition
  endif

  do k=1,nk

     diff_bryan_lewis_90(k) = 1.e-4*(afkph_90 + (dfkph_90/pi)*(atan(sfkph_90*(100.0*Grd%zw(k) - zfkph_90))))
     diff_bryan_lewis_00(k) = 1.e-4*(afkph_00 + (dfkph_00/pi)*(atan(sfkph_00*(100.0*Grd%zw(k) - zfkph_00))))

     if(diff_bryan_lewis_00(k) <= 0.0) then
         call mpp_error(FATAL,&
         '==>Error in ocean_vert_mix_mod: Bryan-Lewis parameters give diffusivities < 0')
     endif
     if ((dtime_t*diff_bryan_lewis_00(k))/Grd%dzt(k)**2 >= 0.5  .and. aidif /= 1.0) then
         write (stdoutunit,'(a,a,i3)')                                                 &
         '==> Warning: vertical diffusive criteria exceeded for "diff_bryan_lewis".',&
              ' use a smaller "dtts" and/or "diff_bryan_lewis" at level k=',k
     endif
  enddo

  if(verbose_init) then 
      write(stdoutunit,'(/a)') ' Bryan-Lewis diffusivities'
      do k=1,nk
         write(stdoutunit,'(a,i4,a,e14.7)') &
         'diff_bryan_lewis_00(m^2/sec)(',k,') = ',diff_bryan_lewis_00(k)
      enddo
      if(bryan_lewis_lat_depend) then 
          write(stdoutunit,'(/a)') ' Bryan-Lewis diffusivities at pole'
          do k=1,nk
             write(stdoutunit,'(a,i4,a,e14.7)') &
             'diff_bryan_lewis_90(m^2/sec)(',k,') = ',diff_bryan_lewis_90(k)
          enddo
      endif
  endif


  ! add Bryan-Lewis to background diffusivity 
  if(bryan_lewis_lat_depend) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               weight = (1+tanh((abs(Grd%yt(i,j))-bryan_lewis_lat_transition)/5.0))/2.0
               diff_cbt_back(i,j,k) = diff_cbt_back(i,j,k) +  &
               Grd%tmask(i,j,k)*((1-weight)*diff_bryan_lewis_00(k) + weight*diff_bryan_lewis_90(k))
            enddo
         enddo
      enddo
  else 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               diff_cbt_back(i,j,k) = diff_cbt_back(i,j,k) +  Grd%tmask(i,j,k)*diff_bryan_lewis_00(k)
            enddo
         enddo
      enddo
  endif


  ! register and send static diffusivities 

  id_diff_bryan_lewis_00 = -1
  id_diff_bryan_lewis_00 = register_static_field ('ocean_model', 'diff_bryan_lewis_00', &
                           Grd%tracer_axes(1:3), 'Bryan-Lewis vertical diffusivity 00', &
                           'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/))
  if (id_diff_bryan_lewis_00 > 0) then 
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,:) = diff_bryan_lewis_00(:)*Grd%tmask(i,j,:)
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_diff_bryan_lewis_00, wrk1(:,:,:))
  endif

  id_diff_bryan_lewis_90 = -1
  id_diff_bryan_lewis_90 = register_static_field ('ocean_model', 'diff_bryan_lewis_90', &
                           Grd%tracer_axes(1:3), 'Bryan-Lewis vertical diffusivity 90', &
                           'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/))
  if (id_diff_bryan_lewis_90 > 0) then 
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,:) = diff_bryan_lewis_90(:)*Grd%tmask(i,j,:)
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_diff_bryan_lewis_90, wrk1(:,:,:))
  endif


end subroutine bryan_lewis_init
! </SUBROUTINE> NAME="bryan_lewis_init"


!#######################################################################
! <SUBROUTINE NAME="hwf_init">
!
! <DESCRIPTION>
! Initialize the HWF static background diffusivity.
!
! Two forms are available:
!
! 1/ Depth dependent form, meant to emulate the Bryan-Lewis approach.  
!    This form is not generally used at GFDL.
!
! 2/ Depth independent form motivated by use in the CM2G isopycnal
!    ocean climate model at GFDL.   
! 
! </DESCRIPTION>
!
subroutine hwf_init(Time, Ocean_options)
  type(ocean_time_type),    intent(in)    :: Time
  type(ocean_options_type), intent(inout) :: Ocean_options

  integer :: i,j,k
  real    :: K_sfc, abs_sin, I_x30
  real    :: I_trans, atan_fn_sfc, I_atan_fn, atan_fn_depth

  if(hwf_diffusivity) then
     Ocean_options%hwf_mix = 'Used HWF background vertical diffusivity.'
  else
     Ocean_options%hwf_mix = 'Did NOT use HWF background vertical diffusivity.'
     return 
  endif 

  if(bryan_lewis_diffusivity) then
    call mpp_error(WARNING,&
    '==>Warning from ocean_vert_mix_mod: Bryan-Lewis and HWF should NOT both be enabled, but they are.')
  endif 
  
  wrk1(:,:,:) = 0.0

  if(hwf_diffusivity_3d) then 

      ! form including depth dependence 
      I_trans     = 1.0/222.22222222
      atan_fn_sfc = atan(hwf_depth_transition*I_trans)
      I_atan_fn   = 1.0/(2.0*atan(1.0) + atan_fn_sfc)
      I_x30       = 1.0/(sin(30.0*deg_to_rad))/invcosh(hwf_N0_2Omega/sin(30.0*deg_to_rad))

      do k=1,nk
         atan_fn_depth = atan((hwf_depth_transition-Grd%zw(k))*I_trans)
         do j=jsc,jec
            do i=isc,iec
               abs_sin  = abs(sin(Grd%yt(i,j)*deg_to_rad))
               K_sfc    = &
                    max(hwf_min_diffusivity,hwf_30_diffusivity*((abs_sin*invcosh(hwf_N0_2Omega/max(1.e-10,abs_sin)))*I_x30))
               wrk1(i,j,k) = K_sfc + (hwf_30_diffusivity - K_sfc)*(atan_fn_sfc - atan_fn_depth)*I_atan_fn
               wrk1(i,j,k) = wrk1(i,j,k)*Grd%tmask(i,j,k)
            enddo
         enddo
      enddo

  else

      ! simplified form without depth dependence, 
      ! used in GFDL CM2G isopycnal ocean climate model 
      I_x30 = (1.0/(sin(30.0*deg_to_rad)))/invcosh(hwf_N0_2Omega/sin(30.0*deg_to_rad))
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               abs_sin     = abs(sin(Grd%yt(i,j)*deg_to_rad))
               K_sfc       = hwf_30_diffusivity*I_x30*abs_sin*invcosh(hwf_N0_2Omega/(epsln+abs_sin))
               K_sfc       = max(hwf_min_diffusivity,K_sfc)
               wrk1(i,j,k) = K_sfc*Grd%tmask(i,j,k)
            enddo
         enddo
      enddo

  endif
  
  ! add to background diffusivity 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           diff_cbt_back(i,j,k) = diff_cbt_back(i,j,k) + Grd%tmask(i,j,k)*wrk1(i,j,k)
        enddo
     enddo
  enddo


  ! register and send static HWF diffusivity
  id_hwf_diffusivity = -1
  id_hwf_diffusivity = register_static_field ('ocean_model', 'hwf_diffusivity',     &
                       Grd%tracer_axes(1:3), 'HWF background vertical diffusivity', &
                       'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/))
  call diagnose_3d(Time, Grd, id_hwf_diffusivity, wrk1(:,:,:))

end subroutine hwf_init
! </SUBROUTINE> NAME="hwf_init"


!#######################################################################
! <SUBROUTINE NAME="diff_cbt_tanh_init">
!
! <DESCRIPTION>
! Initialize the tanh background diffusivity.
! </DESCRIPTION>
!
subroutine diff_cbt_tanh_init(Time, Ocean_options)
  type(ocean_time_type),    intent(in)    :: Time
  type(ocean_options_type), intent(inout) :: Ocean_options

  integer :: i,j,k

  if(diff_cbt_tanh) then
     Ocean_options%tanh_diff_cbt = 'Used tanh background vertical diffusivity.'
  else
     Ocean_options%tanh_diff_cbt = 'Did NOT use tanh background vertical diffusivity.'
     return 
  endif 

  if(bryan_lewis_diffusivity .or. hwf_diffusivity) then
    call mpp_error(WARNING,&
    '==>Warning from ocean_vert_mix_mod: Enabling tanh diff_cbt on top of either Bryan-Lewis and HWF. Not recommended.')
  endif 
  
  ! static background vertical diffusivity 
  allocate (diff_cbt_tanh_back(nk))
  diff_cbt_tanh_back(:) = 0.0 

  ! compute tanh diffusivity and add to background diffusivity 
  do k=1,nk
     diff_cbt_tanh_back(k) = 0.5*diff_cbt_tanh_max &
          *(1.0-tanh((Grd%zt(k)-diff_cbt_tanh_zmid)/diff_cbt_tanh_zwid))
     diff_cbt_tanh_back(k) = max(diff_cbt_tanh_min,diff_cbt_tanh_back(k))
  enddo

  wrk1(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k)          = diff_cbt_tanh_back(k)
           diff_cbt_back(i,j,k) = diff_cbt_back(i,j,k) +  Grd%tmask(i,j,k)*wrk1(i,j,k)
        enddo
     enddo
  enddo

  ! register and send static tanh diffusivity
  id_diff_cbt_tanh = -1
  id_diff_cbt_tanh = register_static_field ('ocean_model', 'diff_cbt_tanh',          &
                       Grd%tracer_axes(1:3), 'Tanh background vertical diffusivity', &
                       'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/))
  call diagnose_3d(Time, Grd, id_diff_cbt_tanh, wrk1(:,:,:))

end subroutine diff_cbt_tanh_init
! </SUBROUTINE> NAME="diff_cbt_tanh_init"


!#######################################################################
! <SUBROUTINE NAME="vert_friction_init">
!
! <DESCRIPTION>
! Initialize vertical friction portion of ocean_vert_mix_mod
! </DESCRIPTION>
!
subroutine vert_friction_init(Time)
  type(ocean_time_type), intent(in) :: Time

  integer :: i,j,k
  integer :: stdoutunit 
  stdoutunit=stdout() 

  ! static background vertical viscosity 
  allocate (visc_cbu_back(nk))
  allocate (visc_cbt_back(nk))
  visc_cbu_back(:) = 0.0 
  visc_cbt_back(:) = 0.0 

  if(vert_visc_back) then 
      do k=1,nk
         visc_cbu_back(k) = 0.5*visc_cbu_back_max &
         *(1.0-tanh((Grd%zt(k)-visc_cbu_back_zmid)/visc_cbu_back_zwid))
         visc_cbu_back(k) = max(visc_cbu_back_min,visc_cbu_back(k))
         visc_cbt_back(k) = visc_cbu_back(k)
      enddo
  endif

  id_visc_cbu_back = -1
  id_visc_cbu_back = register_static_field ('ocean_model', 'visc_cbu_back',   &
                     Grd%vel_axes_wu(1:3), 'static background visc_cbu',      &
                     'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/),&
                     standard_name='ocean_vertical_momentum_diffusivity_due_to_background')

  if (id_visc_cbu_back > 0) then 
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,:) = visc_cbu_back(:)*Grd%umask(i,j,:)
         enddo
      enddo
      call diagnose_3d_u(Time, Grd, id_visc_cbu_back, wrk1(isc:iec,jsc:jec,:))
  endif

  id_visc_cbt_back = -1
  id_visc_cbt_back = register_static_field ('ocean_model', 'visc_cbt_back',     &
                          Grd%tracer_axes_wt(1:3), 'static background visc_cbt',&
                          'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/))
  if(id_visc_cbt_back > 0) then 
     wrk1(:,:,:) = 0.0
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              wrk1(i,j,k) = visc_cbt_back(k)*Grd%tmask(i,j,k)
           enddo
        enddo
     enddo
     call diagnose_3d(Time, Grd, id_visc_cbt_back, wrk1(:,:,:))
  endif

  if(aidif < 1.0) then 
      do k=1,nk
         if (dtime_u*visc_cbu_back(k)/Grd%dzt(k)**2 >= 0.5) then
             write (stdoutunit,'(a,a,i3)')                                               &
                  '==> Warning: vertical stability criteria exceeded for visc_cbu_back.',&
                  ' Use aidif=1.0, smaller dtuv, or smaller visc_cbu_back at level k=',k
         endif
      enddo
  endif


  ! register acceleration due to vertical friction
  id_vfrict_expl_u = register_diag_field('ocean_model','vfrict_expl_u',Grd%vel_axes_u(1:3),&
                     Time%model_time, 'explicit vert friction on u-velocity',              &
                     '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1e20,1e20/))
  id_vfrict_expl_v = register_diag_field('ocean_model','vfrict_expl_v',Grd%vel_axes_v(1:3),&
                     Time%model_time, 'explicit vert friction on v-velocity',              &
                     '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1e20,1e20/))

  id_vfrict_impl_u = register_diag_field ('ocean_model', 'vfrict_impl_u', Grd%vel_axes_u(1:3),&
                     Time%model_time, 'implicit vertical u-mixing', '(kg/m^3)*(m/s^2)',       &
                     missing_value=missing_value, range=(/-1e20,1e20/))
  id_vfrict_impl_v = register_diag_field ('ocean_model', 'vfrict_impl_v', Grd%vel_axes_v(1:3),&
                     Time%model_time, 'implicit vertical v-mixing', '(kg/m^3)*(m/s^2)',       &
                     missing_value=missing_value, range=(/-1e20,1e20/))

  ! register vertical viscosity 
  id_visc_cbu_diabatic = -1
  id_visc_cbu_diabatic = register_diag_field ('ocean_model', 'visc_cbu_diabatic',&
                 Grd%vel_axes_wu(1:3),                                           &
                 Time%model_time, 'vertical viscosity (w/o form drag)', 'm^2/s', &
                 missing_value=missing_value, range=(/-10.0,1e6/))

  id_visc_cbt_diabatic = -1
  id_visc_cbt_diabatic = register_diag_field ('ocean_model', 'visc_cbt_diabatic',&
                 Grd%tracer_axes_wt(1:3),                                        &
                 Time%model_time, 'vertical viscosity (w/o form drag) on T-grid',&
                 'm^2/s',  missing_value=missing_value, range=(/-10.0,1e6/))

  id_visc_cbu_before_min = -1
  id_visc_cbu_before_min = register_diag_field ('ocean_model', 'visc_cbu_before_min',   &
      Grd%vel_axes_wu(1:3),                                                             &
      Time%model_time, 'vertical viscosity before being boosted from a min dissipation',&
      'm^2/s',missing_value=missing_value, range=(/-10.0,1e6/))

  id_visc_cbu  = -1
  id_visc_cbu  = register_diag_field ('ocean_model', 'visc_cbu', Grd%vel_axes_wu(1:3), &
                 Time%model_time, 'total vertical viscosity', 'm^2/s',                 &
                 missing_value=missing_value, range=(/-10.0,1e6/),                     &
                 standard_name='ocean_vertical_momentum_diffusivity')

  id_visc_cbt  = -1
  id_visc_cbt  = register_diag_field ('ocean_model', 'visc_cbt', Grd%tracer_axes_wt(1:3),&
                   Time%model_time, 'total vertical viscosity on T-grid', 'm^2/s',       &
                   missing_value=missing_value, range=(/-10.0,1e6/))

  ! register power from bottom drag and winds 
  id_bottom_power = -1
  id_bottom_power(1) = register_diag_field ('ocean_model', 'bottom_power_u', Grd%vel_axes_u(1:2), &
                       Time%model_time, 'Power dissipation to bottom drag in i-direction',        &
                       'Watt',missing_value=missing_value, range=(/-1e15,1e15/))
  id_bottom_power(2) = register_diag_field ('ocean_model', 'bottom_power_v', Grd%vel_axes_v(1:2), &
                       Time%model_time, 'Power dissipation to bottom drag in j-direction',        &
                       'Watt',missing_value=missing_value, range=(/-1e15,1e15/))

  id_wind_power    = -1
  id_wind_power(1) = register_diag_field ('ocean_model', 'wind_power_u', Grd%vel_axes_u(1:2), &
                     Time%model_time, 'Power from wind stress in i-direction',                &
                     'Watt',missing_value=missing_value, range=(/-1e15,1e15/))
  id_wind_power(2) = register_diag_field ('ocean_model', 'wind_power_v', Grd%vel_axes_v(1:2), &
                     Time%model_time, 'Power from wind stress in j-direction',                &
                     'Watt',missing_value=missing_value, range=(/-1e15,1e15/))


end subroutine vert_friction_init
! </SUBROUTINE> NAME="vert_friction_init"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag_init">
!
! <DESCRIPTION>
! Initialization of watermass diagnostic output files. 
! Also determine the logical compute_watermass_diag.  
! </DESCRIPTION>
!
subroutine watermass_diag_init(Time, Dens)

  type(ocean_time_type),    intent(in) :: Time
  type(ocean_density_type), intent(in) :: Dens

  integer :: stdoutunit 
  stdoutunit=stdout() 
  compute_watermass_diag = .false. 

  !-----------------------------------------------------
  ! vertical diffusion on temp, salt, and rho 

  id_neut_temp_vdiffuse = register_diag_field ('ocean_model', 'neut_temp_vdiffuse',    &
          Grd%tracer_axes(1:3), Time%model_time, 'drhodT * vertical diffusion of temp',&
          '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_temp_vdiffuse > 0) compute_watermass_diag = .true. 

  id_neut_salt_vdiffuse = register_diag_field ('ocean_model', 'neut_salt_vdiffuse',    &
          Grd%tracer_axes(1:3), Time%model_time, 'drhodS * vertical diffusion of salt',&
          '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_salt_vdiffuse > 0) compute_watermass_diag = .true. 

  id_neut_rho_vdiffuse = register_diag_field ('ocean_model', 'neut_rho_vdiffuse',  &
          Grd%tracer_axes(1:3), Time%model_time,                                   &
          'vertical diffusion of locally ref potrho (including stf, btf, and K33)',&
          '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_vdiffuse > 0) compute_watermass_diag = .true. 

  id_wdian_temp_vdiffuse = register_diag_field ('ocean_model', 'wdian_temp_vdiffuse',             &
          Grd%tracer_axes(1:3), Time%model_time,                                                  &
          'dianeutral mass transport due to vert diffusion of temp (including stf, btf, and K33)',&
          'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_temp_vdiffuse > 0) compute_watermass_diag = .true. 

  id_wdian_salt_vdiffuse = register_diag_field ('ocean_model', 'wdian_salt_vdiffuse',             &
          Grd%tracer_axes(1:3), Time%model_time,                                                  &
          'dianeutral mass transport due to vert diffusion of salt (including stf, btf, and K33)',&
          'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_salt_vdiffuse > 0) compute_watermass_diag = .true. 

  id_wdian_rho_vdiffuse = register_diag_field ('ocean_model', 'wdian_rho_vdiffuse',       &
          Grd%tracer_axes(1:3), Time%model_time,                                          &
          'dianeutral mass transport due to vert diffusion (including stf, btf, and K33)',&
          'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_vdiffuse > 0) compute_watermass_diag = .true. 

  id_tform_rho_vdiffuse = register_diag_field ('ocean_model', 'tform_rho_vdiffuse',                               & 
          Grd%tracer_axes(1:3), Time%model_time,                                                                  &
          'watermass transform due to vert diffusion (including stf, btf, and K33) on levels (pre-layer binning)',&
          'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_vdiffuse > 0) compute_watermass_diag = .true. 

  id_neut_rho_vdiffuse_on_nrho = register_diag_field ('ocean_model',                                     &
          'neut_rho_vdiffuse_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
          'vertical diffusion (including stf, btf, and K33) of loc ref potrho binned to neutral density',&
          '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_vdiffuse_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_vdiffuse_on_nrho = register_diag_field ('ocean_model',                                        &
   'wdian_rho_vdiffuse_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                                 &
   'dianeutral mass transport due to vert diffusion (including stf, btf, and K33) binned to neutral density',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_vdiffuse_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_vdiffuse_on_nrho = register_diag_field ('ocean_model',                                         &
   'tform_rho_vdiffuse_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                                  &
   'watermass transform due to vert diffusion (including stf, btf, and K33) binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_vdiffuse_on_nrho > 0) compute_watermass_diag = .true. 


  !-----------------------------------------------------
  ! diff_cbt diffusion on rho 

  id_neut_rho_diff_cbt = register_diag_field ('ocean_model', &
   'neut_rho_diff_cbt',                                      &
    Grd%tracer_axes(1:3), Time%model_time,                  &
   'vertical diffusion of locally ref potrho from diff_cbt',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_diff_cbt > 0) compute_watermass_diag = .true. 

  id_wdian_rho_diff_cbt = register_diag_field ('ocean_model', 'wdian_rho_diff_cbt',&
    Grd%tracer_axes(1:3), Time%model_time,                                         &
   'dianeutral mass transport from diff_cbt',                                      &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_diff_cbt > 0) compute_watermass_diag = .true. 

  id_tform_rho_diff_cbt = register_diag_field ('ocean_model', 'tform_rho_diff_cbt',&
   Grd%tracer_axes(1:3), Time%model_time,                                          &
   'watermass transform due from diff_cbt on levels (pre-layer binning)',          &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_diff_cbt > 0) compute_watermass_diag = .true. 

  id_neut_rho_diff_cbt_on_nrho = register_diag_field ('ocean_model',              &
   'neut_rho_diff_cbt_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,       &
   'vertical diffusion of loc ref potrho from diff_cbt binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_diff_cbt_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_diff_cbt_on_nrho = register_diag_field ('ocean_model',         &
    'wdian_rho_diff_cbt_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time, &
    'dianeutral mass transport due diff_cbt binned to neutral density layers',&
    'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_diff_cbt_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_diff_cbt_on_nrho = register_diag_field ('ocean_model',        &
    'tform_rho_diff_cbt_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,&
    'watermass transform from diff_cbt as binned to neutral density layers', &
    'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_diff_cbt_on_nrho > 0) compute_watermass_diag = .true. 


  !-----------------------------------------------------
  ! temp contribution to diff_cbt diffusion on rho 

  id_neut_temp_diff_cbt = register_diag_field ('ocean_model',            &
   'neut_temp_diff_cbt',                                                 &
   Grd%tracer_axes(1:3), Time%model_time,                                &
   'temp related vertical diffusion of locally ref potrho from diff_cbt',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_temp_diff_cbt > 0) compute_watermass_diag = .true. 

  id_wdian_temp_diff_cbt = register_diag_field ('ocean_model', 'wdian_temp_diff_cbt',&
   Grd%tracer_axes(1:3), Time%model_time,                                            &
   'temp related dianeutral mass transport due from diff_cbt',                       &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_temp_diff_cbt > 0) compute_watermass_diag = .true. 

  id_tform_temp_diff_cbt = register_diag_field ('ocean_model', 'tform_temp_diff_cbt',&
   Grd%tracer_axes(1:3), Time%model_time,                                            &
   'temp related watermass transform due to diff_cbt on levels (pre-layer binning)', &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_temp_diff_cbt > 0) compute_watermass_diag = .true. 

  id_neut_temp_diff_cbt_on_nrho = register_diag_field ('ocean_model',                          &
   'neut_temp_diff_cbt_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                   &
   'temp related vertical diffusion of loc ref potrho from diff_cbt binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_temp_diff_cbt_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_temp_diff_cbt_on_nrho = register_diag_field ('ocean_model',                     &
   'wdian_temp_diff_cbt_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,              &
   'temp related dianeutral mass transport from diff_cbt binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_temp_diff_cbt_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_temp_diff_cbt_on_nrho = register_diag_field ('ocean_model',                  &
   'tform_temp_diff_cbt_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,           &
   'temp related watermass transform from diff_cbt as binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_temp_diff_cbt_on_nrho > 0) compute_watermass_diag = .true. 


  !-----------------------------------------------------
  ! salinity contribution to diff_cbt diffusion on rho 

  id_neut_salt_diff_cbt = register_diag_field ('ocean_model',            &
   'neut_salt_diff_cbt',                                                 &
   Grd%tracer_axes(1:3), Time%model_time,                                &
   'salt related vertical diffusion of locally ref potrho from diff_cbt',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_salt_diff_cbt > 0) compute_watermass_diag = .true. 

  id_wdian_salt_diff_cbt = register_diag_field ('ocean_model', 'wdian_salt_diff_cbt',&
   Grd%tracer_axes(1:3), Time%model_time,                                            &
   'salt related dianeutral mass transport due from diff_cbt',                       &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_salt_diff_cbt > 0) compute_watermass_diag = .true. 

  id_tform_salt_diff_cbt = register_diag_field ('ocean_model', 'tform_salt_diff_cbt',&
   Grd%tracer_axes(1:3), Time%model_time,                                            &
   'salt related watermass transform due to diff_cbt on levels (pre-layer binning)', &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_salt_diff_cbt > 0) compute_watermass_diag = .true. 

  id_neut_salt_diff_cbt_on_nrho = register_diag_field ('ocean_model',                          &
   'neut_salt_diff_cbt_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                   &
   'salt related vertical diffusion of loc ref potrho from diff_cbt binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_salt_diff_cbt_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_salt_diff_cbt_on_nrho = register_diag_field ('ocean_model',                     &
   'wdian_salt_diff_cbt_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,              &
   'salt related dianeutral mass transport from diff_cbt binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_salt_diff_cbt_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_salt_diff_cbt_on_nrho = register_diag_field ('ocean_model',                  &
   'tform_salt_diff_cbt_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,           &
   'salt related watermass transform from diff_cbt as binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_salt_diff_cbt_on_nrho > 0) compute_watermass_diag = .true. 


  !-----------------------------------------------------
  ! k33-implicit diffusion on rho 

  id_neut_rho_k33 = register_diag_field ('ocean_model',&
   'neut_rho_k33',                                     &
   Grd%tracer_axes(1:3), Time%model_time,              &
   'tendency of locally ref potrho from K33impl',      &
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_k33 > 0) compute_watermass_diag = .true. 

  id_wdian_rho_k33 = register_diag_field ('ocean_model', 'wdian_rho_k33',&
   Grd%tracer_axes(1:3), Time%model_time,                                &
   'dianeutral mass transport from K33impl',                             &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_k33 > 0) compute_watermass_diag = .true. 

  id_tform_rho_k33 = register_diag_field ('ocean_model', 'tform_rho_k33',&
   Grd%tracer_axes(1:3), Time%model_time,                                &
   'watermass transform from K33impl on levels (pre-layer binning)',     &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_k33 > 0) compute_watermass_diag = .true. 

  id_neut_rho_k33_on_nrho = register_diag_field ('ocean_model',        &
   'neut_rho_k33_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time, &
   'tendency of loc ref potrho from K33impl binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_k33_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_k33_on_nrho = register_diag_field ('ocean_model',             &
   'wdian_rho_k33_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,      &
   'dianeutral mass transport from K33impl binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_k33_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_k33_on_nrho = register_diag_field ('ocean_model',          &
   'tform_rho_k33_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,   &
   'watermass transform from K33impl as binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_k33_on_nrho > 0) compute_watermass_diag = .true. 

  !-----------------------------------------------------
  ! k33-implicit diffusion on rho due to effects from temperature 

  id_neut_temp_k33 = register_diag_field ('ocean_model',      & 
   'neut_temp_k33',                                           &
   Grd%tracer_axes(1:3), Time%model_time,                     &
   'temp related tendency of locally ref potrho from K33impl',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_temp_k33 > 0) compute_watermass_diag = .true. 

  id_wdian_temp_k33 = register_diag_field ('ocean_model', 'wdian_temp_k33',&
   Grd%tracer_axes(1:3), Time%model_time,                                  &
   'temp related dianeutral mass transport due to vert K33impl',           &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_temp_k33 > 0) compute_watermass_diag = .true. 

  id_tform_temp_k33 = register_diag_field ('ocean_model', 'tform_temp_k33',        &
   Grd%tracer_axes(1:3), Time%model_time,                                          &
   'temp related watermass transform due to K33impl on levels (pre-layer binning)',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_temp_k33 > 0) compute_watermass_diag = .true. 

  id_neut_temp_k33_on_nrho = register_diag_field ('ocean_model',                 &
   'neut_temp_k33_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,          &
   'temp related loc ref potrho tendency from K33impl binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_temp_k33_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_temp_k33_on_nrho = register_diag_field ('ocean_model',                            &
    'wdian_temp_k33_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                    &
    'temp related dianeutral mass transport due to K33impl binned to neutral density layers',&
    'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_temp_k33_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_temp_k33_on_nrho = register_diag_field ('ocean_model',                        &
   'tform_temp_k33_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                 &
   'temp related watermass transform due to K33impl as binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_temp_k33_on_nrho > 0) compute_watermass_diag = .true. 


  !-----------------------------------------------------
  ! k33-implicit diffusion on rho due to effects from salinity 

  id_neut_salt_k33 = register_diag_field ('ocean_model',                   & 
   'neut_salt_k33',                                                        &
   Grd%tracer_axes(1:3), Time%model_time,                                  &
   'salt related tendency of locally ref potrho from K33impl (no stf,btf)',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_salt_k33 > 0) compute_watermass_diag = .true. 

  id_wdian_salt_k33 = register_diag_field ('ocean_model', 'wdian_salt_k33',&
   Grd%tracer_axes(1:3), Time%model_time,                                  &
   'salt related dianeutral mass transport due from K33impl',              &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_salt_k33 > 0) compute_watermass_diag = .true. 

  id_tform_salt_k33 = register_diag_field ('ocean_model', 'tform_salt_k33',        &
   Grd%tracer_axes(1:3), Time%model_time,                                          &
   'salt related watermass transform due to K33impl on levels (pre-layer binning)',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_salt_k33 > 0) compute_watermass_diag = .true. 
  
  id_neut_salt_k33_on_nrho = register_diag_field ('ocean_model',                    &
   'neut_salt_k33_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,             &
   'salt related tendency of loc ref potrho from K33impl binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_salt_k33_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_salt_k33_on_nrho = register_diag_field ('ocean_model',                          &
    'wdian_salt_k33_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
    'salt related dianeutral mass transport from K33impl binned to neutral density layers',&
    'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_salt_k33_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_salt_k33_on_nrho = register_diag_field ('ocean_model',                      &
   'tform_salt_k33_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,               &
   'salt related watermass transform from K33impl as binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_salt_k33_on_nrho > 0) compute_watermass_diag = .true. 



  !-----------------------------------------------------
  ! vertical mixing from diff_cbt + k33-implicit diffusion on rho 

  id_neut_rho_vmix = register_diag_field ('ocean_model',                     &
          'neut_rho_vmix',                                                   &
          Grd%tracer_axes(1:3), Time%model_time,                             &
          'vertical diffusion of locally ref potrho from diff_cbt + K33impl',&
          '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_vmix > 0) compute_watermass_diag = .true. 

  id_wdian_rho_vmix = register_diag_field ('ocean_model', 'wdian_rho_vmix',         &
          Grd%tracer_axes(1:3), Time%model_time,                                    &
          'dianeutral mass transport due to vert diffusion from diff_cbt + K33impl',&
          'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_vmix > 0) compute_watermass_diag = .true. 

  id_tform_rho_vmix = register_diag_field ('ocean_model', 'tform_rho_vmix',                               &
          Grd%tracer_axes(1:3), Time%model_time,                                                          &
          'watermass transform due to vert diffuse from diff_cbt + K33impl on levels (pre-layer binning)',&
          'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_vmix > 0) compute_watermass_diag = .true. 

  id_neut_rho_vmix_on_nrho = register_diag_field ('ocean_model',                                        &
          'neut_rho_vmix_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                          &
          'vertical diffusion of loc ref potrho from from diff_cbt + K33impl binned to neutral density',&
          '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_vmix_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_vmix_on_nrho = register_diag_field ('ocean_model',                                              &
    'wdian_rho_vmix_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                                      &
    'dianeutral mass transport due to vert diffusion from diff_cbt + K33impl binned to neutral density layers',&
    'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_vmix_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_vmix_on_nrho = register_diag_field ('ocean_model',                                           &
    'tform_rho_vmix_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                                   &
    'watermass transform due to vert diffusion from diff_cbt + K33impl as binned to neutral density layers',&
    'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_vmix_on_nrho > 0) compute_watermass_diag = .true. 



  !-----------------------------------------------------
  ! diff_cbt_wave diffusion on rho 

  id_neut_rho_diff_wave = register_diag_field ('ocean_model',    &
   'neut_rho_diff_wave',                                         &
    Grd%tracer_axes(1:3), Time%model_time,                       &
   'vertical diffusion of locally ref potrho from diff_cbt_wave',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_diff_wave > 0) compute_watermass_diag = .true. 

  id_wdian_rho_diff_wave = register_diag_field ('ocean_model', 'wdian_rho_diff_wave',&
    Grd%tracer_axes(1:3), Time%model_time,                                           &
   'dianeutral mass transport from diff_cbt_wave',                                   &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_diff_wave > 0) compute_watermass_diag = .true. 

  id_tform_rho_diff_wave = register_diag_field ('ocean_model', 'tform_rho_diff_wave',&
   Grd%tracer_axes(1:3), Time%model_time,                                            &
   'watermass transform due from diff_cbt_wave on levels (pre-layer binning)',       &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_diff_wave > 0) compute_watermass_diag = .true. 

  id_neut_rho_diff_wave_on_nrho = register_diag_field ('ocean_model',                  &
   'neut_rho_diff_wave_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,           &
   'vertical diffusion of loc ref potrho from diff_cbt_wave binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_diff_wave_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_diff_wave_on_nrho = register_diag_field ('ocean_model',             &
    'wdian_rho_diff_wave_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,     &
    'dianeutral mass transport due diff_cbt_wave binned to neutral density layers',&
    'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_diff_wave_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_diff_wave_on_nrho = register_diag_field ('ocean_model',           &
    'tform_rho_diff_wave_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,   &
    'watermass transform from diff_cbt_wave as binned to neutral density layers',&
    'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_diff_wave_on_nrho > 0) compute_watermass_diag = .true. 


  !-----------------------------------------------------
  ! temp contribution to diff_cbt_wave diffusion on rho 

  id_neut_temp_diff_wave = register_diag_field ('ocean_model',                &
   'neut_temp_diff_wave',                                                     &
   Grd%tracer_axes(1:3), Time%model_time,                                     &
   'temp related vertical diffusion of locally ref potrho from diff_cbt_wave',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_temp_diff_wave > 0) compute_watermass_diag = .true. 

  id_wdian_temp_diff_wave = register_diag_field ('ocean_model', 'wdian_temp_diff_wave',&
   Grd%tracer_axes(1:3), Time%model_time,                                              &
   'temp related dianeutral mass transport due from diff_cbt_wave',                    &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_temp_diff_wave > 0) compute_watermass_diag = .true. 

  id_tform_temp_diff_wave = register_diag_field ('ocean_model', 'tform_temp_diff_wave',  &
   Grd%tracer_axes(1:3), Time%model_time,                                                &
   'temp related watermass transform due to diff_cbt_wave on levels (pre-layer binning)',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_temp_diff_wave > 0) compute_watermass_diag = .true. 

  id_neut_temp_diff_wave_on_nrho = register_diag_field ('ocean_model',                              &
   'neut_temp_diff_wave_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
   'temp related vertical diffusion of loc ref potrho from diff_cbt_wave binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_temp_diff_wave_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_temp_diff_wave_on_nrho = register_diag_field ('ocean_model',                         &
   'wdian_temp_diff_wave_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
   'temp related dianeutral mass transport from diff_cbt_wave binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_temp_diff_wave_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_temp_diff_wave_on_nrho = register_diag_field ('ocean_model',                      &
   'tform_temp_diff_wave_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,               &
   'temp related watermass transform from diff_cbt_wave as binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_temp_diff_wave_on_nrho > 0) compute_watermass_diag = .true. 


  !-----------------------------------------------------
  ! salinity contribution to diff_cbt_wave diffusion on rho 

  id_neut_salt_diff_wave = register_diag_field ('ocean_model',                &
   'neut_salt_diff_wave',                                                     &
   Grd%tracer_axes(1:3), Time%model_time,                                     &
   'salt related vertical diffusion of locally ref potrho from diff_cbt_wave',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_salt_diff_wave > 0) compute_watermass_diag = .true. 

  id_wdian_salt_diff_wave = register_diag_field ('ocean_model', 'wdian_salt_diff_wave',&
   Grd%tracer_axes(1:3), Time%model_time,                                              &
   'salt related dianeutral mass transport due from diff_cbt_wave',                    &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_salt_diff_wave > 0) compute_watermass_diag = .true. 

  id_tform_salt_diff_wave = register_diag_field ('ocean_model', 'tform_salt_diff_wave',  &
   Grd%tracer_axes(1:3), Time%model_time,                                                &
   'salt related watermass transform due to diff_cbt_wave on levels (pre-layer binning)',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_salt_diff_wave > 0) compute_watermass_diag = .true. 

  id_neut_salt_diff_wave_on_nrho = register_diag_field ('ocean_model',                              &
   'neut_salt_diff_wave_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
   'salt related vertical diffusion of loc ref potrho from diff_cbt_wave binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_salt_diff_wave_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_salt_diff_wave_on_nrho = register_diag_field ('ocean_model',                         &
   'wdian_salt_diff_wave_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
   'salt related dianeutral mass transport from diff_cbt_wave binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_salt_diff_wave_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_salt_diff_wave_on_nrho = register_diag_field ('ocean_model',                      &
   'tform_salt_diff_wave_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,               &
   'salt related watermass transform from diff_cbt_wave as binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_salt_diff_wave_on_nrho > 0) compute_watermass_diag = .true. 



  !-----------------------------------------------------
  ! diff_cbt_drag diffusion on rho 

  id_neut_rho_diff_drag = register_diag_field ('ocean_model',    &
   'neut_rho_diff_drag',                                         &
    Grd%tracer_axes(1:3), Time%model_time,                       &
   'vertical diffusion of locally ref potrho from diff_cbt_drag',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_diff_drag > 0) compute_watermass_diag = .true. 

  id_wdian_rho_diff_drag = register_diag_field ('ocean_model', 'wdian_rho_diff_drag',&
    Grd%tracer_axes(1:3), Time%model_time,                                           &
   'dianeutral mass transport from diff_cbt_drag',                                   &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_diff_drag > 0) compute_watermass_diag = .true. 

  id_tform_rho_diff_drag = register_diag_field ('ocean_model', 'tform_rho_diff_drag',&
   Grd%tracer_axes(1:3), Time%model_time,                                            &
   'watermass transform due from diff_cbt_drag on levels (pre-layer binning)',       &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_diff_drag > 0) compute_watermass_diag = .true. 

  id_neut_rho_diff_drag_on_nrho = register_diag_field ('ocean_model',                  &
   'neut_rho_diff_drag_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,           &
   'vertical diffusion of loc ref potrho from diff_cbt_drag binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_diff_drag_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_diff_drag_on_nrho = register_diag_field ('ocean_model',             &
    'wdian_rho_diff_drag_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,     &
    'dianeutral mass transport due diff_cbt_drag binned to neutral density layers',&
    'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_diff_drag_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_diff_drag_on_nrho = register_diag_field ('ocean_model',           &
    'tform_rho_diff_drag_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,   &
    'watermass transform from diff_cbt_drag as binned to neutral density layers',&
    'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_diff_drag_on_nrho > 0) compute_watermass_diag = .true. 


  !-----------------------------------------------------
  ! temp contribution to diff_cbt_drag diffusion on rho 

  id_neut_temp_diff_drag = register_diag_field ('ocean_model',                &
   'neut_temp_diff_drag',                                                     &
   Grd%tracer_axes(1:3), Time%model_time,                                     &
   'temp related vertical diffusion of locally ref potrho from diff_cbt_drag',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_temp_diff_drag > 0) compute_watermass_diag = .true. 

  id_wdian_temp_diff_drag = register_diag_field ('ocean_model', 'wdian_temp_diff_drag',&
   Grd%tracer_axes(1:3), Time%model_time,                                              &
   'temp related dianeutral mass transport due from diff_cbt_drag',                    &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_temp_diff_drag > 0) compute_watermass_diag = .true. 

  id_tform_temp_diff_drag = register_diag_field ('ocean_model', 'tform_temp_diff_drag',  &
   Grd%tracer_axes(1:3), Time%model_time,                                                &
   'temp related watermass transform due to diff_cbt_drag on levels (pre-layer binning)',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_temp_diff_drag > 0) compute_watermass_diag = .true. 

  id_neut_temp_diff_drag_on_nrho = register_diag_field ('ocean_model',                              &
   'neut_temp_diff_drag_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
   'temp related vertical diffusion of loc ref potrho from diff_cbt_drag binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_temp_diff_drag_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_temp_diff_drag_on_nrho = register_diag_field ('ocean_model',                         &
   'wdian_temp_diff_drag_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
   'temp related dianeutral mass transport from diff_cbt_drag binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_temp_diff_drag_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_temp_diff_drag_on_nrho = register_diag_field ('ocean_model',                      &
   'tform_temp_diff_drag_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,               &
   'temp related watermass transform from diff_cbt_drag as binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_temp_diff_drag_on_nrho > 0) compute_watermass_diag = .true. 


  !-----------------------------------------------------
  ! salinity contribution to diff_cbt_drag diffusion on rho 

  id_neut_salt_diff_drag = register_diag_field ('ocean_model',                &
   'neut_salt_diff_drag',                                                     &
   Grd%tracer_axes(1:3), Time%model_time,                                     &
   'salt related vertical diffusion of locally ref potrho from diff_cbt_drag',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_salt_diff_drag > 0) compute_watermass_diag = .true. 

  id_wdian_salt_diff_drag = register_diag_field ('ocean_model', 'wdian_salt_diff_drag',&
   Grd%tracer_axes(1:3), Time%model_time,                                              &
   'salt related dianeutral mass transport due from diff_cbt_drag',                    &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_salt_diff_drag > 0) compute_watermass_diag = .true. 

  id_tform_salt_diff_drag = register_diag_field ('ocean_model', 'tform_salt_diff_drag',  &
   Grd%tracer_axes(1:3), Time%model_time,                                                &
   'salt related watermass transform due to diff_cbt_drag on levels (pre-layer binning)',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_salt_diff_drag > 0) compute_watermass_diag = .true. 

  id_neut_salt_diff_drag_on_nrho = register_diag_field ('ocean_model',                              &
   'neut_salt_diff_drag_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
   'salt related vertical diffusion of loc ref potrho from diff_cbt_drag binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_salt_diff_drag_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_salt_diff_drag_on_nrho = register_diag_field ('ocean_model',                         &
   'wdian_salt_diff_drag_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
   'salt related dianeutral mass transport from diff_cbt_drag binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_salt_diff_drag_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_salt_diff_drag_on_nrho = register_diag_field ('ocean_model',                      &
   'tform_salt_diff_drag_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,               &
   'salt related watermass transform from diff_cbt_drag as binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_salt_diff_drag_on_nrho > 0) compute_watermass_diag = .true. 



  !-----------------------------------------------------
  ! diff_cbt_lee diffusion on rho 

  id_neut_rho_diff_lee = register_diag_field ('ocean_model',    &
   'neut_rho_diff_lee',                                         &
    Grd%tracer_axes(1:3), Time%model_time,                      &
   'vertical diffusion of locally ref potrho from diff_cbt_lee',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_diff_lee > 0) compute_watermass_diag = .true. 

  id_wdian_rho_diff_lee = register_diag_field ('ocean_model', 'wdian_rho_diff_lee',&
    Grd%tracer_axes(1:3), Time%model_time,                                         &
   'dianeutral mass transport from diff_cbt_lee',                                  &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_diff_lee > 0) compute_watermass_diag = .true. 

  id_tform_rho_diff_lee = register_diag_field ('ocean_model', 'tform_rho_diff_lee',&
   Grd%tracer_axes(1:3), Time%model_time,                                          &
   'watermass transform due from diff_cbt_lee on levels (pre-layer binning)',      &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_diff_lee > 0) compute_watermass_diag = .true. 

  id_neut_rho_diff_lee_on_nrho = register_diag_field ('ocean_model',                  &
   'neut_rho_diff_lee_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,           &
   'vertical diffusion of loc ref potrho from diff_cbt_lee binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_diff_lee_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_diff_lee_on_nrho = register_diag_field ('ocean_model',             &
    'wdian_rho_diff_lee_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,     &
    'dianeutral mass transport due diff_cbt_lee binned to neutral density layers',&
    'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_diff_lee_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_diff_lee_on_nrho = register_diag_field ('ocean_model',           &
    'tform_rho_diff_lee_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,   &
    'watermass transform from diff_cbt_lee as binned to neutral density layers',&
    'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_diff_lee_on_nrho > 0) compute_watermass_diag = .true. 


  !-----------------------------------------------------
  ! temp contribution to diff_cbt_lee diffusion on rho 

  id_neut_temp_diff_lee = register_diag_field ('ocean_model',                &
   'neut_temp_diff_lee',                                                     &
   Grd%tracer_axes(1:3), Time%model_time,                                    &
   'temp related vertical diffusion of locally ref potrho from diff_cbt_lee',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_temp_diff_lee > 0) compute_watermass_diag = .true. 

  id_wdian_temp_diff_lee = register_diag_field ('ocean_model', 'wdian_temp_diff_lee',&
   Grd%tracer_axes(1:3), Time%model_time,                                            &
   'temp related dianeutral mass transport due from diff_cbt_lee',                   &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_temp_diff_lee > 0) compute_watermass_diag = .true. 

  id_tform_temp_diff_lee = register_diag_field ('ocean_model', 'tform_temp_diff_lee',   &
   Grd%tracer_axes(1:3), Time%model_time,                                               &
   'temp related watermass transform due to diff_cbt_lee on levels (pre-layer binning)',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_temp_diff_lee > 0) compute_watermass_diag = .true. 

  id_neut_temp_diff_lee_on_nrho = register_diag_field ('ocean_model',                              &
   'neut_temp_diff_lee_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
   'temp related vertical diffusion of loc ref potrho from diff_cbt_lee binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_temp_diff_lee_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_temp_diff_lee_on_nrho = register_diag_field ('ocean_model',                         &
   'wdian_temp_diff_lee_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
   'temp related dianeutral mass transport from diff_cbt_lee binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_temp_diff_lee_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_temp_diff_lee_on_nrho = register_diag_field ('ocean_model',                      &
   'tform_temp_diff_lee_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,               &
   'temp related watermass transform from diff_cbt_lee as binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_temp_diff_lee_on_nrho > 0) compute_watermass_diag = .true. 


  !-----------------------------------------------------
  ! salinity contribution to diff_cbt_lee diffusion on rho 

  id_neut_salt_diff_lee = register_diag_field ('ocean_model',                &
   'neut_salt_diff_lee',                                                     &
   Grd%tracer_axes(1:3), Time%model_time,                                    &
   'salt related vertical diffusion of locally ref potrho from diff_cbt_lee',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_salt_diff_lee > 0) compute_watermass_diag = .true. 

  id_wdian_salt_diff_lee = register_diag_field ('ocean_model', 'wdian_salt_diff_lee',&
   Grd%tracer_axes(1:3), Time%model_time,                                            &
   'salt related dianeutral mass transport due from diff_cbt_lee',                   &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_salt_diff_lee > 0) compute_watermass_diag = .true. 

  id_tform_salt_diff_lee = register_diag_field ('ocean_model', 'tform_salt_diff_lee',   &
   Grd%tracer_axes(1:3), Time%model_time,                                               &
   'salt related watermass transform due to diff_cbt_lee on levels (pre-layer binning)',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_salt_diff_lee > 0) compute_watermass_diag = .true. 

  id_neut_salt_diff_lee_on_nrho = register_diag_field ('ocean_model',                              &
   'neut_salt_diff_lee_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
   'salt related vertical diffusion of loc ref potrho from diff_cbt_lee binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_salt_diff_lee_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_salt_diff_lee_on_nrho = register_diag_field ('ocean_model',                         &
   'wdian_salt_diff_lee_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
   'salt related dianeutral mass transport from diff_cbt_lee binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_salt_diff_lee_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_salt_diff_lee_on_nrho = register_diag_field ('ocean_model',                      &
   'tform_salt_diff_lee_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,               &
   'salt related watermass transform from diff_cbt_lee as binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_salt_diff_lee_on_nrho > 0) compute_watermass_diag = .true. 



  !-----------------------------------------------------
  ! diff_cbt_back diffusion on rho 

  id_neut_rho_diff_back = register_diag_field ('ocean_model',    &
   'neut_rho_diff_back',                                         &
    Grd%tracer_axes(1:3), Time%model_time,                       &
   'vertical diffusion of locally ref potrho from diff_cbt_back',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_diff_back > 0) compute_watermass_diag = .true. 

  id_wdian_rho_diff_back = register_diag_field ('ocean_model', 'wdian_rho_diff_back',&
    Grd%tracer_axes(1:3), Time%model_time,                                           &
   'dianeutral mass transport from diff_cbt_back',                                   &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_diff_back > 0) compute_watermass_diag = .true. 

  id_tform_rho_diff_back = register_diag_field ('ocean_model', 'tform_rho_diff_back',&
   Grd%tracer_axes(1:3), Time%model_time,                                            &
   'watermass transform due from diff_cbt_back on levels (pre-layer binning)',       &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_diff_back > 0) compute_watermass_diag = .true. 

  id_neut_rho_diff_back_on_nrho = register_diag_field ('ocean_model',                  &
   'neut_rho_diff_back_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,           &
   'vertical diffusion of loc ref potrho from diff_cbt_back binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_diff_back_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_diff_back_on_nrho = register_diag_field ('ocean_model',             &
    'wdian_rho_diff_back_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,     &
    'dianeutral mass transport due diff_cbt_back binned to neutral density layers',&
    'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_diff_back_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_diff_back_on_nrho = register_diag_field ('ocean_model',           &
    'tform_rho_diff_back_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,   &
    'watermass transform from diff_cbt_back as binned to neutral density layers',&
    'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_diff_back_on_nrho > 0) compute_watermass_diag = .true. 


  !-----------------------------------------------------
  ! temp contribution to diff_cbt_back diffusion on rho 

  id_neut_temp_diff_back = register_diag_field ('ocean_model',                &
   'neut_temp_diff_back',                                                     &
   Grd%tracer_axes(1:3), Time%model_time,                                     &
   'temp related vertical diffusion of locally ref potrho from diff_cbt_back',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_temp_diff_back > 0) compute_watermass_diag = .true. 

  id_wdian_temp_diff_back = register_diag_field ('ocean_model', 'wdian_temp_diff_back',&
   Grd%tracer_axes(1:3), Time%model_time,                                              &
   'temp related dianeutral mass transport due from diff_cbt_back',                    &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_temp_diff_back > 0) compute_watermass_diag = .true. 

  id_tform_temp_diff_back = register_diag_field ('ocean_model', 'tform_temp_diff_back',  &
   Grd%tracer_axes(1:3), Time%model_time,                                                &
   'temp related watermass transform due to diff_cbt_back on levels (pre-layer binning)',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_temp_diff_back > 0) compute_watermass_diag = .true. 

  id_neut_temp_diff_back_on_nrho = register_diag_field ('ocean_model',                              &
   'neut_temp_diff_back_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
   'temp related vertical diffusion of loc ref potrho from diff_cbt_back binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_temp_diff_back_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_temp_diff_back_on_nrho = register_diag_field ('ocean_model',                         &
   'wdian_temp_diff_back_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
   'temp related dianeutral mass transport from diff_cbt_back binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_temp_diff_back_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_temp_diff_back_on_nrho = register_diag_field ('ocean_model',                      &
   'tform_temp_diff_back_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,               &
   'temp related watermass transform from diff_cbt_back as binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_temp_diff_back_on_nrho > 0) compute_watermass_diag = .true. 


  !-----------------------------------------------------
  ! salinity contribution to diff_cbt_back diffusion on rho 

  id_neut_salt_diff_back = register_diag_field ('ocean_model',                &
   'neut_salt_diff_back',                                                     &
   Grd%tracer_axes(1:3), Time%model_time,                                     &
   'salt related vertical diffusion of locally ref potrho from diff_cbt_back',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_salt_diff_back > 0) compute_watermass_diag = .true. 

  id_wdian_salt_diff_back = register_diag_field ('ocean_model', 'wdian_salt_diff_back',&
   Grd%tracer_axes(1:3), Time%model_time,                                              &
   'salt related dianeutral mass transport due from diff_cbt_back',                    &
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_salt_diff_back > 0) compute_watermass_diag = .true. 

  id_tform_salt_diff_back = register_diag_field ('ocean_model', 'tform_salt_diff_back',  &
   Grd%tracer_axes(1:3), Time%model_time,                                                &
   'salt related watermass transform due to diff_cbt_back on levels (pre-layer binning)',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_salt_diff_back > 0) compute_watermass_diag = .true. 

  id_neut_salt_diff_back_on_nrho = register_diag_field ('ocean_model',                              &
   'neut_salt_diff_back_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
   'salt related vertical diffusion of loc ref potrho from diff_cbt_back binned to neutral density',&
   '(kg/m^3)/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_salt_diff_back_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_salt_diff_back_on_nrho = register_diag_field ('ocean_model',                         &
   'wdian_salt_diff_back_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
   'salt related dianeutral mass transport from diff_cbt_back binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_salt_diff_back_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_salt_diff_back_on_nrho = register_diag_field ('ocean_model',                      &
   'tform_salt_diff_back_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,               &
   'salt related watermass transform from diff_cbt_back as binned to neutral density layers',&
   'kg/sec',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_salt_diff_back_on_nrho > 0) compute_watermass_diag = .true. 



  !-----------------------------------------------------
  ! surface boundary fluxes on rho 

  id_neut_rho_sbc = register_diag_field ('ocean_model',                & 
        'neut_rho_sbc', Grd%tracer_axes(1:3), Time%model_time,         &
        'update of local ref potrho due to surface temp and salt flux',&
        '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_sbc > 0) compute_watermass_diag = .true. 

  id_wdian_rho_sbc = register_diag_field ('ocean_model',              & 
        'wdian_rho_sbc', Grd%tracer_axes(1:3), Time%model_time,       &
        'dianeutral mass transport due to surface temp and salt flux',&
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_sbc > 0) compute_watermass_diag = .true. 

  id_tform_rho_sbc = register_diag_field ('ocean_model',                                    & 
        'tform_rho_sbc', Grd%tracer_axes(1:3), Time%model_time,                             &
        'net surface heat and salt flux water mass transform on levels (pre-layer binning)',&
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_sbc > 0) compute_watermass_diag = .true. 

  id_neut_rho_sbc_on_nrho = register_diag_field ('ocean_model',                                          & 
        'neut_rho_sbc_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                              &
        'update of local ref potrho due to surface temp and salt flux binned to neutral density classes',&
        '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_sbc_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_sbc_on_nrho = register_diag_field ('ocean_model',                                           & 
        'wdian_rho_sbc_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                               &
        'dianeutral mass transport due to surface temp and salt flux as binned to neutral density classes',&
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_sbc_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_sbc_on_nrho = register_diag_field ('ocean_model',                              &  
        'tform_rho_sbc_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
        'net surface heat and salt flux water mass transformation in neutral density classes',&
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_sbc_on_nrho > 0) compute_watermass_diag = .true. 

  !-----------------------------------------------------
  ! surface boundary fluxes on temp

  id_neut_rho_sbc_temp = register_diag_field ('ocean_model',       & 
        'neut_rho_sbc_temp', Grd%tracer_axes(1:3), Time%model_time,&
        'update of local ref potrho due to surface temp flux',     &
        '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_sbc_temp > 0) compute_watermass_diag = .true. 

  id_neut_rho_sbc_temp_on_nrho = register_diag_field ('ocean_model',                            & 
        'neut_rho_sbc_temp_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                &
        'update of local ref potrho due to surface temp flux binned to neutral density classes',&
        '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_sbc_temp_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_sbc_temp = register_diag_field ('ocean_model',       & 
        'wdian_rho_sbc_temp', Grd%tracer_axes(1:3), Time%model_time,&
        'dianeutral mass transport due to surface temp flux',       &
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_sbc_temp > 0) compute_watermass_diag = .true. 

  id_wdian_rho_sbc_temp_on_nrho = register_diag_field ('ocean_model',                             & 
        'wdian_rho_sbc_temp_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                 &
        'dianeutral mass transport due to surface temp flux as binned to neutral density classes',&
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_sbc_temp_on_nrho > 0) compute_watermass_diag = .true. 


  id_tform_rho_sbc_temp = register_diag_field ('ocean_model',                 &  
        'tform_rho_sbc_temp', Grd%tracer_axes(1:3), Time%model_time,          &
        'net surface heat flux water mass transformation (pre-layer binning)',&
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_sbc_temp > 0) compute_watermass_diag = .true. 

  id_tform_rho_sbc_temp_on_nrho = register_diag_field ('ocean_model',                &  
        'tform_rho_sbc_temp_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,    &
        'net surface heat flux water mass transformation in neutral density classes',&
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_sbc_temp_on_nrho > 0) compute_watermass_diag = .true. 

  !-----------------------------------------------------
  ! surface boundary fluxes on salt 

  id_neut_rho_sbc_salt = register_diag_field ('ocean_model',       & 
        'neut_rho_sbc_salt', Grd%tracer_axes(1:3), Time%model_time,&
        'update of local ref potrho due to surface salt flux',     &
        '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_sbc_salt > 0) compute_watermass_diag = .true. 

  id_neut_rho_sbc_salt_on_nrho = register_diag_field ('ocean_model',                            & 
        'neut_rho_sbc_salt_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                &
        'update of local ref potrho due to surface salt flux binned to neutral density classes',&
        '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_sbc_salt_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_sbc_salt = register_diag_field ('ocean_model',       & 
        'wdian_rho_sbc_salt', Grd%tracer_axes(1:3), Time%model_time,&
        'dianeutral mass transport due to surface salt flux',       &
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_sbc_salt > 0) compute_watermass_diag = .true. 

  id_wdian_rho_sbc_salt_on_nrho = register_diag_field ('ocean_model',                             & 
        'wdian_rho_sbc_salt_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                 &
        'dianeutral mass transport due to surface salt flux as binned to neutral density classes',&
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_sbc_salt_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_sbc_salt = register_diag_field ('ocean_model',                 &  
        'tform_rho_sbc_salt', Dens%neutralrho_axes(1:3), Time%model_time,     &
        'net surface salt flux water mass transformation (pre-layer binning)',&
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_sbc_salt > 0) compute_watermass_diag = .true. 

  id_tform_rho_sbc_salt_on_nrho = register_diag_field ('ocean_model',                &  
        'tform_rho_sbc_salt_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,    &
        'net surface salt flux water mass transformation in neutral density classes',&
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_sbc_salt_on_nrho > 0) compute_watermass_diag = .true. 


  !-----------------------------------------------------
  ! bottom boundary fluxes on temp

  id_neut_rho_bbc_temp = register_diag_field ('ocean_model',            & 
        'neut_rho_bbc_temp', Grd%tracer_axes(1:3), Time%model_time,     &
        'update of local ref potrho due to bottom geothermal heat flux',&
        '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_bbc_temp > 0) compute_watermass_diag = .true. 

  id_wdian_rho_bbc_temp = register_diag_field ('ocean_model',          & 
        'wdian_rho_bbc_temp', Grd%tracer_axes(1:3), Time%model_time,   &
        'dianeutral mass transport due to bottom geothermal heat flux',&
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_bbc_temp > 0) compute_watermass_diag = .true. 

  id_tform_rho_bbc_temp = register_diag_field ('ocean_model',                                  & 
        'tform_rho_bbc_temp', Grd%tracer_axes(1:3), Time%model_time,                           &
        'watermass transform due to bottom geothermal heat flux on levels (pre-layer binning)',&
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_bbc_temp > 0) compute_watermass_diag = .true. 

  id_tform_rho_bbc_temp_on_nrho = register_diag_field ('ocean_model',                      &   
        'tform_rho_bbc_temp_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,          &
        'bottom geothermal heat flux water mass transformation in neutral density classes',&
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_bbc_temp_on_nrho > 0) compute_watermass_diag = .true. 

  id_neut_rho_bbc_temp_on_nrho = register_diag_field ('ocean_model',                                      & 
        'neut_rho_bbc_temp_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                          &
        'update of local ref potrho due to bottom geothermal heat flux binned to neutral density classes',&
        '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_bbc_temp_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_bbc_temp_on_nrho = register_diag_field ('ocean_model',                            & 
        'wdian_rho_bbc_temp_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                &
        'dianeutral mass transport due to bottom heat flux as binned to neutral density classes',&
        'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_bbc_temp_on_nrho > 0) compute_watermass_diag = .true. 


  ! sea level contributions 
  id_eta_tend_diff_cbt_tend = register_diag_field ('ocean_model',                &
          'eta_tend_diff_cbt_tend', Grd%tracer_axes(1:2), Time%model_time,       &
          'non-Bouss steric sea level tendency from diff_cbt diffusion tendency',&
          'm/s',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_eta_tend_diff_cbt_tend > 0) compute_watermass_diag = .true. 

  id_eta_tend_diff_cbt_tend_glob = register_diag_field ('ocean_model',           &
          'eta_tend_diff_cbt_tend_glob', Time%model_time,                        &
          'non-Bouss steric sea level tendency from diff_cbt diffusion tendency',&
          'm/s',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_eta_tend_diff_cbt_tend_glob > 0) compute_watermass_diag = .true. 

  id_eta_tend_k33_tend = register_diag_field ('ocean_model',                &
          'eta_tend_k33_tend', Grd%tracer_axes(1:2), Time%model_time,       &
          'non-Bouss steric sea level tendency from k33 diffusion tendency',&
          'm/s',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_eta_tend_k33_tend > 0) compute_watermass_diag = .true. 

  id_eta_tend_k33_tend_glob = register_diag_field ('ocean_model',           &
          'eta_tend_k33_tend_glob', Time%model_time,                        &
          'non-Bouss steric sea level tendency from k33 diffusion tendency',&
          'm/s',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_eta_tend_k33_tend_glob > 0) compute_watermass_diag = .true. 


  id_eta_tend_diff_wave_tend = register_diag_field ('ocean_model',                    &
          'eta_tend_diff_wave_tend', Grd%tracer_axes(1:2), Time%model_time,           &
          'non-Bouss steric sea level tendency from diff_cbt_wave diffusion tendency',&
          'm/s',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_eta_tend_diff_wave_tend > 0) compute_watermass_diag = .true. 

  id_eta_tend_diff_wave_tend_glob = register_diag_field ('ocean_model',               &
          'eta_tend_diff_wave_tend_glob', Time%model_time,                            &
          'non-Bouss steric sea level tendency from diff_cbt_wave diffusion tendency',&
          'm/s',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_eta_tend_diff_wave_tend_glob > 0) compute_watermass_diag = .true. 


  id_eta_tend_diff_drag_tend = register_diag_field ('ocean_model',                    &
          'eta_tend_diff_drag_tend', Grd%tracer_axes(1:2), Time%model_time,           &
          'non-Bouss steric sea level tendency from diff_cbt_drag diffusion tendency',&
          'm/s',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_eta_tend_diff_drag_tend > 0) compute_watermass_diag = .true. 

  id_eta_tend_diff_drag_tend_glob = register_diag_field ('ocean_model',               &
          'eta_tend_diff_drag_tend_glob', Time%model_time,                            &
          'non-Bouss steric sea level tendency from diff_cbt_drag diffusion tendency',&
          'm/s',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_eta_tend_diff_drag_tend_glob > 0) compute_watermass_diag = .true. 


  id_eta_tend_diff_lee_tend = register_diag_field ('ocean_model',                    &
          'eta_tend_diff_lee_tend', Grd%tracer_axes(1:2), Time%model_time,           &
          'non-Bouss steric sea level tendency from diff_cbt_lee diffusion tendency',&
          'm/s',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_eta_tend_diff_lee_tend > 0) compute_watermass_diag = .true. 

  id_eta_tend_diff_lee_tend_glob = register_diag_field ('ocean_model',               &
          'eta_tend_diff_lee_tend_glob', Time%model_time,                            &
          'non-Bouss steric sea level tendency from diff_cbt_lee diffusion tendency',&
          'm/s',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_eta_tend_diff_lee_tend_glob > 0) compute_watermass_diag = .true. 


  id_eta_tend_diff_back_tend = register_diag_field ('ocean_model',                    &
          'eta_tend_diff_back_tend', Grd%tracer_axes(1:2), Time%model_time,           &
          'non-Bouss steric sea level tendency from diff_cbt_back diffusion tendency',&
          'm/s',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_eta_tend_diff_back_tend > 0) compute_watermass_diag = .true. 

  id_eta_tend_diff_back_tend_glob = register_diag_field ('ocean_model',               &
          'eta_tend_diff_back_tend_glob', Time%model_time,                            &
          'non-Bouss steric sea level tendency from diff_cbt_back diffusion tendency',&
          'm/s',missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_eta_tend_diff_back_tend_glob > 0) compute_watermass_diag = .true. 



  if(compute_watermass_diag) then 
    write(stdoutunit,'(/a/)') &
    '==>Note: running ocean_vert_mix_mod w/ compute_watermass_diag=.true.'  
  endif 

end subroutine watermass_diag_init
! </SUBROUTINE> NAME="watermass_diag_init"


!#######################################################################
! <SUBROUTINE NAME="vert_mix_coeff">
!
! <DESCRIPTION>
! This subroutine calls the relevant scheme to compute vertical
! diffusivity and vertical viscosity.  Background values are 
! also incorporated here.
! </DESCRIPTION>
!
subroutine vert_mix_coeff(Time, Thickness, Velocity, T_prog,   &
                          T_diag, Dens, swflx, sw_frac_zt, pme,&
                          river, visc_cbu, visc_cbt, diff_cbt, hblt_depth, do_wave)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_velocity_type),      intent(in)    :: Velocity
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  type(ocean_diag_tracer_type),   intent(in)    :: T_diag(:)
  type(ocean_density_type),       intent(in)    :: Dens
  real, dimension(isd:,jsd:),     intent(in)    :: swflx
  real, dimension(isd:,jsd:,:),   intent(in)    :: sw_frac_zt
  real, dimension(isd:,jsd:),     intent(in)    :: pme
  real, dimension(isd:,jsd:),     intent(in)    :: river
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbt
  real, dimension(isd:,jsd:,:,:), intent(inout) :: diff_cbt
  real, dimension(isd:,jsd:),     intent(inout) :: hblt_depth
  logical, intent(in) :: do_wave

  integer :: i,j,k,kp1,n,tau
  real    :: tmp, rescale 
  real    :: dnu_dtheta_dz, dnu_dsalinity_dz, dtheta_dz, dsalinity_dz  
  real    :: global_mean 

  tau = Time%tau 

  ! initialize diffusivity and viscosity to zero 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           diff_cbt(i,j,k,1) = 0.0
           diff_cbt(i,j,k,2) = 0.0
           visc_cbu(i,j,k)   = 0.0
           visc_cbt(i,j,k)   = 0.0
        enddo
     enddo
  enddo

  if(MIX_SCHEME == VERTMIX_CONST) then 
    call mpp_clock_begin(id_clock_vert_const)
    call vert_mix_const(aidif, Time, T_prog, Dens, visc_cbu, visc_cbt, diff_cbt)
    call mpp_clock_end(id_clock_vert_const)

  elseif(MIX_SCHEME == VERTMIX_KPP_MOM4P0) then 
    call mpp_clock_begin(id_clock_vert_kpp_mom4p0)
    call vert_mix_kpp_mom4p0(aidif, Time, Thickness, Velocity, T_prog, T_diag, Dens, &
                      swflx, sw_frac_zt, pme, river, visc_cbu, diff_cbt, hblt_depth, do_wave)
    ! since this scheme is frozen, we do not compute visc_cbt. 
    ! for vertical reynolds diagnostics, we set  
    visc_cbt = visc_cbu
    call mpp_clock_end(id_clock_vert_kpp_mom4p0)

  elseif(MIX_SCHEME == VERTMIX_KPP_MOM4P1) then 
    call mpp_clock_begin(id_clock_vert_kpp_mom4p1)
    call vert_mix_kpp_mom4p1(aidif, Time, Thickness, Velocity, T_prog, T_diag, Dens, &
                      swflx, sw_frac_zt, pme, river, visc_cbu, diff_cbt, hblt_depth, do_wave)
    ! since this scheme is frozen, we do not compute visc_cbt. 
    ! for vertical reynolds diagnostics, we set  
    visc_cbt = visc_cbu
    call mpp_clock_end(id_clock_vert_kpp_mom4p1)

  elseif(MIX_SCHEME == VERTMIX_KPP_TEST) then 
    call mpp_clock_begin(id_clock_vert_kpp_test)
    call vert_mix_kpp_test(aidif, Time, Thickness, Velocity, T_prog, T_diag, Dens, &
                           swflx, sw_frac_zt, pme, river, visc_cbu, visc_cbt, diff_cbt, hblt_depth, do_wave)
    call mpp_clock_end(id_clock_vert_kpp_test)

  elseif(MIX_SCHEME == VERTMIX_PP) then 
    call mpp_clock_begin(id_clock_vert_pp)
    call vert_mix_pp(Time, Thickness, Velocity, T_prog, Dens, visc_cbu, visc_cbt, diff_cbt)
    call mpp_clock_end(id_clock_vert_pp)

  elseif(MIX_SCHEME == VERTMIX_CHEN) then 
    call mpp_clock_begin(id_clock_vert_chen)
    call vert_mix_chen(Time, Thickness, Velocity, T_prog, Dens, &
                       swflx, pme, visc_cbu, visc_cbt, diff_cbt)
    call mpp_clock_end(id_clock_vert_chen)

  elseif(MIX_SCHEME == VERTMIX_GOTM .and. horz_grid == MOM_BGRID) then 
    call mpp_clock_begin(id_clock_vert_gotm)
    call vert_mix_gotm_bgrid(Time, Thickness, Velocity, T_prog, Dens, visc_cbu, visc_cbt, diff_cbt)
    call mpp_clock_end(id_clock_vert_gotm)

  elseif(MIX_SCHEME == VERTMIX_GOTM .and. horz_grid == MOM_CGRID) then 
    ! a Cgrid wrapper for GOTM needs to be written. until then, will use bgrid version.  
    call mpp_clock_begin(id_clock_vert_gotm)
    call vert_mix_gotm_bgrid(Time, Thickness, Velocity, T_prog, Dens, visc_cbu, visc_cbt, diff_cbt)
    call mpp_clock_end(id_clock_vert_gotm)

  endif 

  ! incorporate static background values 
  if(vert_diff_back_via_max) then 
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               diff_cbt(i,j,k,1) = max(diff_cbt(i,j,k,1)  , diff_cbt_back(i,j,k))
               diff_cbt(i,j,k,2) = max(diff_cbt(i,j,k,2)  , diff_cbt_back(i,j,k))
               visc_cbu(i,j,k)   = max(visc_cbu(i,j,k)    , visc_cbu_back(k))
               visc_cbt(i,j,k)   = max(visc_cbt(i,j,k)    , visc_cbt_back(k))
            enddo
         enddo
      enddo
  else 
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               diff_cbt(i,j,k,1) = diff_cbt(i,j,k,1) + diff_cbt_back(i,j,k)
               diff_cbt(i,j,k,2) = diff_cbt(i,j,k,2) + diff_cbt_back(i,j,k)
               visc_cbu(i,j,k)   = visc_cbu(i,j,k)   + visc_cbu_back(k)
               visc_cbt(i,j,k)   = visc_cbt(i,j,k)   + visc_cbt_back(k)
            enddo
         enddo
      enddo
  endif

  ! incorporate diffusivity from tidal dissipation 
  call mpp_clock_begin(id_clock_vert_tidal)
  call vert_mix_tidal(Time, Thickness, T_prog, Dens, diff_cbt, visc_cbu, visc_cbt, &
                      diff_cbt_wave, diff_cbt_leewave, diff_cbt_drag)
  call mpp_clock_end(id_clock_vert_tidal)

  if(vmix_rescale_nonbouss) then 
      wrk1(:,:,:) = 0.0
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               tmp = onehalf*(Dens%rho(i,j,k,tau)+Dens%rho(i,j,k+1,tau)) + epsln
               wrk1(i,j,k) = Grd%tmask(i,j,k+1)*rho0/tmp
            enddo
         enddo
      enddo
      ! use same rescaling for visc_cbu as for diff_cbt to 
      ! avoid problems near land-sea boundaries
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               diff_cbt(i,j,k,1) = diff_cbt(i,j,k,1)*wrk1(i,j,k)
               diff_cbt(i,j,k,2) = diff_cbt(i,j,k,2)*wrk1(i,j,k)
               visc_cbu(i,j,k)   = visc_cbu(i,j,k)*wrk1(i,j,k)
               visc_cbt(i,j,k)   = visc_cbt(i,j,k)*wrk1(i,j,k)
            enddo
         enddo
      enddo
  endif

  if(vmix_set_min_dissipation) then 
      call vmix_min_dissipation(Time, Dens, diff_cbt, visc_cbu, visc_cbt)
  endif 


  !--------------diagnostics--------------------

  call diagnose_3d(Time, Grd, id_diff_cbt_t, diff_cbt(:,:,:,1))
  call diagnose_3d(Time, Grd, id_diff_cbt_s, diff_cbt(:,:,:,2))
  call diagnose_3d_u(Time, Grd, id_visc_cbu_diabatic, visc_cbu(:,:,:))
  call diagnose_3d(Time, Grd, id_visc_cbt_diabatic, visc_cbt(:,:,:))


  ! compute power dissipation by full temp and salt diffusivity
  ! working against temperature and salinity stratification 
  if(id_power_diss > 0) then 
      wrk1(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               tmp =  Dens%drhodT(i,j,k)*Dens%dTdz_zt(i,j,k)*diff_cbt(i,j,k,1) &
                     +Dens%drhodS(i,j,k)*Dens%dSdz_zt(i,j,k)*diff_cbt(i,j,k,2)
               wrk1(i,j,k) = -tmp*grav*Thickness%dzt(i,j,k)*Grd%tmask(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_power_diss, wrk1(:,:,:))
  endif

  ! compute power dissipation by background mixing working against stratification 
  if(id_power_diss_back > 0) then 
      wrk1(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               tmp = grav*Thickness%dzt(i,j,k)*abs(Dens%drhodz_zt(i,j,k))*Grd%tmask(i,j,k)
               wrk1(i,j,k) = tmp*diff_cbt_back(i,j,k) 
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_power_diss_back, wrk1(:,:,:))
  endif

  ! compute contribution to steric sea level evolution
  ! from vertical diffusion in the ocean interior.  
  if(id_eta_tend_diff_cbt_flx  > 0 .or. id_eta_tend_diff_cbt_flx_glob  > 0) then 

      wrk1(:,:,:) = 0.0   !  nu_theta    = alpha/rho  = -drhodT/rho**2
      wrk2(:,:,:) = 0.0   !  nu_salinity = -beta/rho  = -drhodS/rho**2
      wrk3(:,:,:) = 0.0   !  partial_z(nu_theta)    in 1/(deg C * rho * meter)
      wrk4(:,:,:) = 0.0   !  partial_z(nu_salinity) in 1/(ppt   * rho * meter)
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = -Grd%tmask(i,j,k)*Dens%drhodT(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2) 
               wrk2(i,j,k) = -Grd%tmask(i,j,k)*Dens%drhodS(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2) 
            enddo
         enddo
      enddo

      ! compute partial_z(nu_theta) and partial_z(nu_salinity) at bottom of T-cell 
      do k=1,nk-1
         kp1=k+1
         do j=jsc,jec
            do i=isc,iec
               wrk3(i,j,k) = Grd%tmask(i,j,kp1)*(wrk1(i,j,k)-wrk1(i,j,kp1))/Thickness%dzwt(i,j,k)
               wrk4(i,j,k) = Grd%tmask(i,j,kp1)*(wrk2(i,j,k)-wrk2(i,j,kp1))/Thickness%dzwt(i,j,k)
            enddo
         enddo
      enddo

      ! compute vertical derivatives of theta and S at bottom of T-cell
      wrk1(:,:,:) = 0.0   !  dTdz_wt  (deg C/m)
      wrk2(:,:,:) = 0.0   !  dSdz_wt  (ppt  /m)
      do k=1,nk-1
         kp1=k+1
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = Grd%tmask(i,j,kp1)*                                           & 
               (T_prog(index_temp)%field(i,j,k,tau)-T_prog(index_temp)%field(i,j,kp1,tau)) &
               /Thickness%dzwt(i,j,k)
               wrk2(i,j,k) = Grd%tmask(i,j,kp1)*                                           & 
               (T_prog(index_salt)%field(i,j,k,tau)-T_prog(index_salt)%field(i,j,kp1,tau)) &
               /Thickness%dzwt(i,j,k)
            enddo
         enddo
      enddo

      ! sea level tendency from vertical diffusion 
      wrk1_2d(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               tmp = -Thickness%rho_dzt(i,j,k,tau)                 &
                     *(diff_cbt(i,j,k,1)*wrk1(i,j,k)*wrk3(i,j,k) + &
                       diff_cbt(i,j,k,2)*wrk2(i,j,k)*wrk4(i,j,k))
               wrk1_2d(i,j) = wrk1_2d(i,j) + tmp
            enddo
         enddo
      enddo

      call diagnose_2d(Time, Grd, id_eta_tend_diff_cbt_flx, wrk1_2d(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_diff_cbt_flx_glob, wrk1_2d, cellarea_r)
      
  endif  ! id_eta_tend_diff_cbt_flx > 0 .or. id_eta_tend_diff_cbt_flx_glob > 0 



end subroutine vert_mix_coeff
! </SUBROUTINE> NAME="vert_mix_coeff"


!#######################################################################
! <SUBROUTINE NAME="vert_diffuse">
!
! <DESCRIPTION>
! This subroutine computes the thickness weighted and density 
! weighted time tendency for tracer associated with vertical diffusion
! based on explicit time update of the vertical diffusion equation. 
!
! MOM only supports aidif==0.0 or aidif==1.0.
! MOM does not support cases with 0.0 < aidif < 1.0.
!
! The watermass diagnostics have not been ported to this subroutine
! since aidif=0 is rarely used, even in idealized studies.   
!
! This routine is generally never used, since nearly all the vertical
! physical parameterizations allow for large vertical mixing coefficients,
! thus requiring implicit vertical time stepping.
!
! </DESCRIPTION>
!
subroutine vert_diffuse (Time, Thickness, ntracer, Tracer, diff_cbt, diag_flag) 

  type(ocean_time_type),          intent(in)     :: Time 
  type(ocean_thickness_type),     intent(in)     :: Thickness
  integer,                        intent(in)     :: ntracer     
  type(ocean_prog_tracer_type),   intent(inout)  :: Tracer     
  real, dimension(isd:,jsd:,:,:), intent(in)     :: diff_cbt   
  logical, optional,              intent(in)     :: diag_flag          

  real,dimension(isd:ied,jsd:jed)  :: ft1, ft2
  logical                          :: send_diagnostics 
  integer                          :: taum1, tau, nmix
  integer                          :: i, j, k, kp

  if(.not. use_explicit_vert_diffuse) then 
    return 
  endif 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_vert_mix_mod (vert_diffuse): module must be initialized')
  endif 

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  if (Tracer%name == 'salt') then
     nmix = 2
  else
     nmix = 1
  endif

  ! fully implicit
  if (aidif==1.0) then  
    wrk1 = 0.0

  ! fully explicit 
  ! note Tracer%btf assumed to be zero in this calculation; 
  ! must use aidif=1.0 to have nonzero btf.  
  elseif (aidif < 1.0) then

    ft1(isc:iec,jsc:jec) = Tracer%stf(isc:iec,jsc:jec)
    k = 1
    kp = min(k+1,nk)
    ft2(isc:iec,jsc:jec) = rho0*Grd%tmask(isc:iec,jsc:jec,kp)*diff_cbt(isc:iec,jsc:jec,k,nmix)*&
                           (Tracer%field(isc:iec,jsc:jec,k,taum1)-Tracer%field(isc:iec,jsc:jec,kp,taum1))&
                           /Thickness%dzwt(isc:iec,jsc:jec,k)
    wrk1(isc:iec,jsc:jec,k) = ft1(isc:iec,jsc:jec)-ft2(isc:iec,jsc:jec)
    ft1(isc:iec,jsc:jec)    = ft2(isc:iec,jsc:jec)
    do k=2,nk
      kp = min(k+1,nk)
      ft2(isc:iec,jsc:jec) = rho0*Grd%tmask(isc:iec,jsc:jec,kp)*diff_cbt(isc:iec,jsc:jec,k,nmix)*&
                           (Tracer%field(isc:iec,jsc:jec,k,taum1)-Tracer%field(isc:iec,jsc:jec,kp,taum1))&
                           /Thickness%dzwt(isc:iec,jsc:jec,k)
      if (id_zflux_diff(ntracer) > 0) then
          flux_z(isc:iec,jsc:jec,k) = ft2(isc:iec,jsc:jec)
      endif
      wrk1(isc:iec,jsc:jec,k) = ft1(isc:iec,jsc:jec)-ft2(isc:iec,jsc:jec)
      ft1(isc:iec,jsc:jec)    = ft2(isc:iec,jsc:jec)
    enddo
  endif

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           Tracer%th_tendency(i,j,k) = Tracer%th_tendency(i,j,k) + wrk1(i,j,k)
        enddo
     enddo
  enddo

  ! diagnostics
  if(send_diagnostics) then 

      if (id_zflux_diff(ntracer) > 0) then 
         call diagnose_3d(Time, Grd, id_zflux_diff(ntracer), flux_z(:,:,:)*Tracer%conversion)
      endif

      if (id_vdiffuse(ntracer) > 0) then
         call diagnose_3d(Time, Grd, id_vdiffuse(ntracer), wrk1(:,:,:)*Tracer%conversion)
      endif

  endif       ! endif for send_diagnostics


end subroutine vert_diffuse
! </SUBROUTINE> NAME="vert_diffuse"


!#######################################################################
! <SUBROUTINE NAME="vert_diffuse_implicit">
!
! <DESCRIPTION>
! Contributions to thickness weighted and density weighted time
! tendency from time-implicit vertical diffusion.
! </DESCRIPTION>
!
subroutine vert_diffuse_implicit(diff_cbt, index_salt, Time, Thickness, Dens, T_prog)

  real, dimension(isd:,jsd:,:,:), intent(in)    :: diff_cbt
  integer,                        intent(in)    :: index_salt 
  type(ocean_time_type),          intent(in)    :: Time  
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)

  integer :: taup1, tau
  integer :: i, j, k, n, nmix

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(aidif==0.0) return 

  taup1 = Time%taup1
  tau   = Time%tau

  wrk1      = 0.0 
  wrk2      = 0.0
  wrk3      = 0.0
  wrk1_vmix = 0.0
  wrk2_vmix = 0.0
  wrk1_v    = 0.0  

  ! watermass diagnostic calls prior to updating the taup1 temp and salinity values 
  call mpp_clock_begin(id_clock_watermass_diag)
  call watermass_diag(Time, T_prog, Dens, Thickness, diff_cbt)
  call mpp_clock_end(id_clock_watermass_diag)


  do n=1,num_prog_tracers

     ! skip vertical diffusion for "generic" tracers, as we 
     ! compute vertical physical processes for generic tracers 
     ! within the generic tracer module. 
     if(T_prog(n)%type .eq. 'generic') cycle

     if (n==index_salt) then
         nmix=2
     else
         nmix=1
     endif

     ! K33_implicit is already density weighted, but diff_cbt is not.
     ! Need density weighting to have consistent dimensions with stf and btf
     ! sent to invtri. So multiply diff_cbt by rho0.  
     ! save the pre-implicit taup1 value of tracer for diagnostics.  
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              wrk1_vmix(i,j,k) = T_prog(n)%field(i,j,k,taup1)
              wrk2(i,j,k)      = rho0*diff_cbt(i,j,k,nmix) + T_prog(n)%K33_implicit(i,j,k)
           enddo
        enddo
     enddo

     ! invert tridiagonal matrix to update tracer concentration.  
     ! include boundary fluxes from stf and btf in the inversion.  
     ! use rho_dzt(taup1) for mass per area in grid cell. 
     call invtri (T_prog(n)%field(:,:,:,taup1), T_prog(n)%stf, &
                  T_prog(n)%btf, wrk2(:,:,:), dtime_t, Grd%kmt,&
                  Grd%tmask, Thickness%rho_dzt(:,:,:,taup1), Thickness%dzwt, aidif, nk)  

     if (debug_this_module) then
         write(stdoutunit,*) ' ' 
         write(stdoutunit,*) 'tracers after implicit vertical diffusion (taup1) ==>'
         call write_timestamp(Time%model_time)
         call tracer_prog_chksum(Time, T_prog(n), taup1)
         call tracer_min_max(Time, Thickness, T_prog(n))
     endif

     ! diagnose some pieces of implicit vertical diffusion 
     call vert_diffuse_implicit_diag(Time, Thickness, T_prog, diff_cbt, wrk1_vmix, n)

  enddo  ! enddo for num_prog_tracers 

  ! diagnose some pieces of implicit vertical diffusion 
  call vert_diffuse_watermass_diag(Time, Thickness, Dens, T_prog, diff_cbt)

end subroutine vert_diffuse_implicit
! </SUBROUTINE> NAME="vert_diffuse_implicit"


!#######################################################################
! <SUBROUTINE NAME="vert_friction_bgrid">
!
! <DESCRIPTION>
! This subroutine computes the thickness weighted and density weighted
! acceleration (kg/m^3)*(m^2/s^2) associated with vertical friction.
!
! Assumes here that the horizontal grid is B-grid.  
!
! For aidif=1.0, this module does nothing since all vertical friction  
! is in this case handled implicitly in time, and this is computed
! elswewhere. 
!
! MOM only supports aidif==0.0 or aidif==1.0.
! MOM does not support cases with 0.0 < aidif < 1.0.
!
! Note that smf and bmf have units (kg/m^3)*(m^2/s^2)
! So the vertical diffusive fluxes must be in these units
! too.  For this purpose, we multiply the viscosity by 
! rho0. This is an approximation consistent with the 
! Boussinesq approximation.  For non-Boussinesq, we 
! should be using the in situ rho.  But to the level of 
! accuracy that we know the vertical viscosity, and to 
! the extent that the ocean density is close to rho0,
! our use of rho0 for the non-Boussinesq case is of 
! minor consequence for vertical friction calculation.
!
! Note: the form drag contribution to vertical viscosity
! must be handled within aidif=1.0 implicit vertical
! mixing.  It is ignored in this routine, as its 
! contribution would generally cause the model to be 
! unstable. 
! 
! Note: if try to merge this routine with vert_friction_cgrid
! some machines and compilers will change bits by the mere 
! introduction of extra if-test logic into the calculation.
! So we define the separate routines to maintain bit-wise 
! agreement with older results, with bit-wise agreement a 
! useful means to check for errors as the model evolves.
!
! </DESCRIPTION>
!
subroutine vert_friction_bgrid (Time, Thickness, Velocity, visc_cbu, energy_analysis_step)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_velocity_type),      intent(inout) :: Velocity       
  real, dimension(isd:,jsd:,:),   intent(in)    :: visc_cbu       
  logical,                        intent(in)    :: energy_analysis_step

  real,dimension(isd:ied,jsd:jed) :: ft1, ft2
  integer                         :: k, i, j, n
  integer                         :: taum1, tau

  taum1    = Time%taum1
  tau      = Time%tau 
  wrk1_v   = 0.0

  if (aidif < 1.0) then 

      do n=1,2

         do j=jsc,jec
            do i=isc,iec
               ft1(i,j) = Velocity%smf_bgrid(i,j,n)
            enddo
         enddo

         do k=1,nk-1
            do j=jsc,jec
               do i=isc,iec
                   ft2(i,j) =  (1.0-Grd%umask(i,j,k+1))*Velocity%bmf(i,j,n) &
                             + Grd%umask(i,j,k+1)*rho0*visc_cbu(i,j,k)      &
                             *(Velocity%u(i,j,k,n,taum1)-Velocity%u(i,j,k+1,n,taum1))/Thickness%dzwu(i,j,k)
                  wrk1_v(i,j,k,n) = Grd%umask(i,j,k)*(ft1(i,j)-ft2(i,j))
                  ft1(i,j) = ft2(i,j)
               enddo
            enddo
         enddo
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,nk,n) = Grd%umask(i,j,nk)*(ft1(i,j)-Velocity%bmf(i,j,n))
            enddo
         enddo

      enddo

  endif

  if (aidif > 0.0 .and. aidif < 1.0) then ! mixed implicit/explicit
    wrk1_v(isc:iec,jsc:jec,:,:) = wrk1_v(isc:iec,jsc:jec,:,:)*(1.0-aidif)
  endif

  if(energy_analysis_step) then 

      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%wrkv(i,j,k,n) = wrk1_v(i,j,k,n) 
               enddo
            enddo
         enddo
      enddo

  else 

      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%accel(i,j,k,n) = Velocity%accel(i,j,k,n) + wrk1_v(i,j,k,n)  
               enddo
            enddo
         enddo
      enddo

      call diagnose_3d_u(Time, Grd, id_vfrict_expl_u, wrk1_v(:,:,:,1))
      call diagnose_3d_u(Time, Grd, id_vfrict_expl_v, wrk1_v(:,:,:,2))

     if (id_visc_cbu > 0 .and. aidif /= 1.0) then
        call diagnose_3d_u(Time, Grd, id_visc_cbu, visc_cbu(:,:,:))
     endif 


     ! power from wind stress 
     if(id_wind_power(1) > 0 .or. id_wind_power(2) > 0) then 
         wrk1_v2d(:,:,:) = 0.0
         k=1 
         do n=1,2
            do j=jsc,jec
               do i=isc,iec
                  wrk1_v2d(i,j,n) = &
                  Grd%umask(i,j,k)*Grd%dau(i,j)*Velocity%smf_bgrid(i,j,n)*Velocity%u(i,j,k,n,tau)
               enddo
            enddo
         enddo
         call diagnose_2d_u(Time, Grd, id_wind_power(1), wrk1_v2d(:,:,1))
         call diagnose_2d_u(Time, Grd, id_wind_power(2), wrk1_v2d(:,:,2))
     endif

     ! power dissipated to bottom drag 
     if(id_bottom_power(1) > 0 .or. id_bottom_power(2) > 0) then 
         wrk1_v2d(:,:,:) = 0.0
         do n=1,2
            do j=jsc,jec
               do i=isc,iec
                  k=Grd%kmu(i,j)
                  if (k /= 0) then
                      wrk1_v2d(i,j,n) = &
                      -Grd%umask(i,j,k)*Grd%dau(i,j)*Velocity%bmf(i,j,n)*Velocity%u(i,j,k,n,tau)
                  endif
               enddo
            enddo
         enddo
         call diagnose_2d_u(Time, Grd, id_bottom_power(1), wrk1_v2d(:,:,1))
         call diagnose_2d_u(Time, Grd, id_bottom_power(2), wrk1_v2d(:,:,2))
     endif

  endif !  endif for energy analysis step 

end subroutine vert_friction_bgrid
! </SUBROUTINE> NAME="vert_friction_bgrid"



!#######################################################################
! <SUBROUTINE NAME="vert_friction_cgrid">
!
! <DESCRIPTION>
! This subroutine computes the thickness weighted and density weighted
! acceleration (kg/m^3)*(m^2/s^2) associated with vertical friction.
! Assumes horizontal grid is C-grid.
!
! For aidif=1.0, this module does nothing since all vertical friction 
! is instead handled implicitly in time, and this is computed in the 
! routine vert_friction_implicit_cgrid.
!
! MOM only supports aidif==0.0 or aidif==1.0.
! MOM does not support cases with 0.0 < aidif < 1.0.
!
! Note that smf and bmf have units (kg/m^3)*(m^2/s^2)
! So the vertical diffusive fluxes must be in these units
! too.  For this purpose, we multiply the viscosity by 
! rho0. This is an approximation consistent with the 
! Boussinesq approximation.  For non-Boussinesq, we 
! should be using the in situ rho.  But to the level of 
! accuracy that we know the vertical viscosity, and to 
! the extent that the ocean density is close to rho0,
! our use of rho0 for the non-Boussinesq case is of 
! minor consequence for vertical friction calculation.
!
! Note: use visc_cbt for both C-grid velocity components,
! even though the velocity components sit at different 
! sides of the tracer cell.  This choice is for simplicity.
! It also acknowledges that the alternative of introducing 
! distinct visc_cbt_u and visc_cbt_v would presume knowledge
! of subgrid scale features that we really do not have.  
! So again, visc_cbt is used for both u,v C-grid velocity 
! components.  
!
! Note: the form drag contribution to vertical viscosity
! must be handled within aidif=1.0 implicit vertical
! mixing.  It is ignored in this routine, as its 
! contribution would generally cause the model to be 
! unstable. 
! 
! Note: if try to merge this routine with vert_friction_cgrid
! some machines and compilers will change bits by the mere 
! introduction of extra if-test logic into the calculation.
! So we define the separate routines to maintain bit-wise 
! agreement with older results, with bit-wise agreement a 
! useful means to check for errors as the model evolves.
!
! </DESCRIPTION>
!
subroutine vert_friction_cgrid (Time, Thickness, Velocity, visc_cbt, energy_analysis_step)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_velocity_type),      intent(inout) :: Velocity       
  real, dimension(isd:,jsd:,:),   intent(in)    :: visc_cbt
  logical,                        intent(in)    :: energy_analysis_step

  real,dimension(isd:ied,jsd:jed) :: ft1, ft2
  integer                         :: k, i, j, n
  integer                         :: taum1, tau

  taum1  = Time%taum1
  tau    = Time%tau 
  wrk1_v = 0.0

  if (aidif < 1.0) then 

     do n=1,2

        do j=jsc,jec
           do i=isc,iec
              ft1(i,j) = Velocity%smf_cgrid(i,j,n)
           enddo
        enddo

        do k=1,nk-1
           do j=jsc,jec
              do i=isc,iec
                  ft2(i,j) = (1.0-Grd%tmasken(i,j,k+1,n))*Velocity%bmf(i,j,n) &
                            + Grd%tmasken(i,j,k+1,n)*rho0*visc_cbt(i,j,k) &
                            *(Velocity%u(i,j,k,n,taum1)-Velocity%u(i,j,k+1,n,taum1))/(epsln+Thickness%dzten(i,j,k,n))
                 wrk1_v(i,j,k,n) = Grd%tmasken(i,j,k,n)*(ft1(i,j)-ft2(i,j))
                 ft1(i,j) = ft2(i,j)
              enddo
           enddo
        enddo
        do j=jsc,jec
           do i=isc,iec
              wrk1_v(i,j,nk,n) = Grd%tmasken(i,j,nk,n)*(ft1(i,j)-Velocity%bmf(i,j,n))
           enddo
        enddo

     enddo 

  endif

  if (aidif > 0.0 .and. aidif < 1.0) then ! mixed implicit/explicit
    wrk1_v(isc:iec,jsc:jec,:,:) = wrk1_v(isc:iec,jsc:jec,:,:)*(1.0-aidif)
  endif

  if(energy_analysis_step) then 

      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%wrkv(i,j,k,n) = wrk1_v(i,j,k,n) 
               enddo
            enddo
         enddo
      enddo

  else  ! not an analysis step  

      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%accel(i,j,k,n) = Velocity%accel(i,j,k,n) + wrk1_v(i,j,k,n)  
               enddo
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_vfrict_expl_u, wrk1_v(:,:,:,1))
      call diagnose_3d(Time, Grd, id_vfrict_expl_v, wrk1_v(:,:,:,2))

     if (id_visc_cbu > 0 .and. aidif /= 1.0) then
         wrk1(:,:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1(i,j,k)  = onefourth*Grd%umask(i,j,k)* &
                  (visc_cbt(i,j,k) + visc_cbt(i+1,j,k) + visc_cbt(i,j+1,k) + visc_cbt(i+1,j+1,k))
               enddo
           enddo
        enddo
        call diagnose_3d_u(Time, Grd, id_visc_cbu, wrk1(:,:,:))
     endif 

     if (id_visc_cbt > 0 .and. aidif /= 1.0) then
        call diagnose_3d(Time, Grd, id_visc_cbt, visc_cbt(:,:,:))
     endif 

     ! power from wind stress 
     if(id_wind_power(1) > 0 .or. id_wind_power(2) > 0) then 
         wrk1_v2d(:,:,:) = 0.0
         k=1 
         do n=1,2
            do j=jsc,jec
               do i=isc,iec
                  wrk1_v2d(i,j,n) = &
                  Grd%tmask(i,j,k)*Grd%dat(i,j)*Velocity%smf_cgrid(i,j,n)*Velocity%u(i,j,k,n,tau)
               enddo
            enddo
         enddo
         call diagnose_2d_en(Time, Grd, id_wind_power(1), id_wind_power(2), wrk1_v2d(:,:,:))
     endif

     ! power dissipated to bottom drag 
     if(id_bottom_power(1) > 0 .or. id_bottom_power(2) > 0) then 
         wrk1_v2d(:,:,:) = 0.0
         do n=1,2
            do j=jsc,jec
               do i=isc,iec
                  k=Grd%kmt(i,j)
                  if (k /= 0) then
                      wrk1_v2d(i,j,n) = &
                      -Grd%tmasken(i,j,k,n)*Grd%dat(i,j)*Velocity%bmf(i,j,n)*Velocity%u(i,j,k,n,tau)
                  endif
               enddo
            enddo
         enddo
         call diagnose_2d_en(Time, Grd, id_bottom_power(1), id_bottom_power(2), wrk1_v2d(:,:,:))
     endif

  endif !  endif for energy analysis step


end subroutine vert_friction_cgrid
! </SUBROUTINE> NAME="vert_friction_cgrid"


!#######################################################################
! <SUBROUTINE NAME="vert_friction_implicit_bgrid">
!
! <DESCRIPTION>
! Contributions to thickness weighted and density weighted acceleration 
! from implicit vertical friction. 
! 
! Assume that the horizontal grid is B-grid.  
!
! Note that smf and bmf have units N/m^2.  These are the natural units 
! for surface stress.  To include these stresses as boundary terms in the 
! call to invtri, it is necessary to use vertical viscosities with units
! (kg/m^3)*(m2^/s) = N/m^2.  This is achieved by multiplying visc_cbu
! by rho0 when sent to invtri.  For depth-like vertical coordinates, this
! is cancelled exactly by the rho0 in rho_dzu.  For pressure-like 
! vertical coordinates, the rho0*visc_cbu introduces a negligible 
! change in the vertical viscosity that is well within uncertainty
! in this coefficient. 
!
! Include visc_cbu_form_drag to each of the velocity components 
! vertical friction.  
!
! Note: if try to merge this routine with vert_friction_cgrid
! some machines and compilers will change bits by the mere 
! introduction of extra if-test logic into the calculation.
! So we define the separate routines to maintain bit-wise 
! agreement with older results, with bit-wise agreement a 
! useful means to check for errors as the model evolves.
!
! </DESCRIPTION>
!
subroutine vert_friction_implicit_bgrid(visc_cbu, visc_cbu_form_drag, Time, Thickness, Velocity)

  type(ocean_time_type),          intent(in)    :: Time  
  type(ocean_thickness_type),     intent(in)    :: Thickness
  real, dimension(isd:,jsd:,:),   intent(in)    :: visc_cbu
  real, dimension(isd:,jsd:,:,:), intent(in)    :: visc_cbu_form_drag
  type(ocean_velocity_type),      intent(inout) :: Velocity

  integer :: taum1, tau, taup1
  integer :: i, j, k, n

  integer :: stdoutunit 
  stdoutunit=stdout() 

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  Velocity%vfrict_impl = 0.0
  wrk1_v               = 0.0
  wrk1                 = 0.0
  wrk1_2d              = 0.0

  ! acceleration due to time-implicit vertical friction 
  if (aidif /= 0.0) then  

      ! construct updated velocity due to time-explicit contributions, sans barotropic forcing
      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1_v(i,j,k,n) = Thickness%rho_dzur(i,j,k) &
                   *(Thickness%rho_dzu(i,j,k,taum1)*Velocity%u(i,j,k,n,taum1) &
                     + dtime_u*Velocity%accel(i,j,k,n))
               enddo
            enddo
         enddo
      enddo

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               Velocity%vfrict_impl(i,j,k,1) = wrk1_v(i,j,k,1)
               Velocity%vfrict_impl(i,j,k,2) = wrk1_v(i,j,k,2)
               wrk1(i,j,k)                   = Thickness%rho_dzu(i,j,k,taup1)
            enddo
         enddo
      enddo

      ! invert tridiagonal for vertical friction 
      if(Velocity%bmf_implicit) then 

          call invtri_bmf (wrk1_v(:,:,:,1), Velocity%smf_bgrid(:,:,1), Velocity%gamma(:,:),&
               rho0*(visc_cbu(:,:,:)+visc_cbu_form_drag(:,:,:,1)), dtime_u, Grd%kmu,       &
               Grd%umask, wrk1, Thickness%dzwu, aidif, nk)
          call invtri_bmf (wrk1_v(:,:,:,2), Velocity%smf_bgrid(:,:,2), Velocity%gamma(:,:),&
               rho0*(visc_cbu(:,:,:)+visc_cbu_form_drag(:,:,:,2)), dtime_u, Grd%kmu,       &
               Grd%umask, wrk1, Thickness%dzwu, aidif, nk)

      else 

          call invtri (wrk1_v(:,:,:,1), Velocity%smf_bgrid(:,:,1), Velocity%bmf(:,:,1),&
               rho0*(visc_cbu(:,:,:)+visc_cbu_form_drag(:,:,:,1)), dtime_u, Grd%kmu,   &
               Grd%umask, wrk1, Thickness%dzwu, aidif, nk)
          call invtri (wrk1_v(:,:,:,2), Velocity%smf_bgrid(:,:,2), Velocity%bmf(:,:,2),&
               rho0*(visc_cbu(:,:,:)+visc_cbu_form_drag(:,:,:,2)), dtime_u, Grd%kmu,   &
               Grd%umask, wrk1, Thickness%dzwu, aidif, nk)
      endif

      ! compute time-implicit tendency for diagnostics 
      ! update thickness weighted and density weighted 
      ! acceleration due to implicit vertical friction.
      if( .NOT. quebec_2009_10_bug) then !with the riga bug fix
      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%vfrict_impl(i,j,k,n) = &
                   dtime_ur*(wrk1_v(i,j,k,n)-Velocity%vfrict_impl(i,j,k,n))
                  Velocity%accel(i,j,k,n)       = Velocity%accel(i,j,k,n)  &
                  + Thickness%rho_dzu(i,j,k,taup1)*Velocity%vfrict_impl(i,j,k,n)
               enddo
            enddo
         enddo
      enddo
      
      else !without the riga bug fix 

      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%vfrict_impl(i,j,k,n) = &
                   dtime_ur*(wrk1_v(i,j,k,n)-Velocity%vfrict_impl(i,j,k,n))
                  Velocity%accel(i,j,k,n)       = Velocity%accel(i,j,k,n)  &
                       + Thickness%rho_dzu(i,j,k,taup1)*Velocity%vfrict_impl(i,j,k,n)
                  !
                  ! Due to a compiler bug the following line affects answers in the 
                  ! static + production mode (no fltconsistency).
                  ! This line was in the quebec_2009_10 version. To reproduce those
                  ! answers we must keep it.
                  ! Note that if we keep this line we have to modify the 
                  ! send_data statements for this diagnostics variable.
                  !
                  Velocity%vfrict_impl(i,j,k,n) = Thickness%rho_dzu(i,j,k,taup1)*Velocity%vfrict_impl(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

      endif
      if (debug_this_module) then
          write(stdoutunit,*) ' ' 
          write(stdoutunit,*) 'From ocean_vert_mix_mod: chksums after implicit vert frict'
          call write_timestamp(Time%model_time)
          call write_chksum_3d('accel(1)', Velocity%accel(COMP,:,1))
          call write_chksum_3d('accel(2)', Velocity%accel(COMP,:,2))
          call write_chksum_3d('vel(1)', wrk1_v(COMP,:,1))
          call write_chksum_3d('vel(2)', wrk1_v(COMP,:,2))
      endif


  endif  !endif for aidif /= 0

  if( .NOT. quebec_2009_10_bug) then !diagnostics with the bug fix
  ! note that earlier versions of MOM4p1 failed to multiply Velocity%vfrict_impl
  ! by rho_dzu.  we do not perform this multiplication in the lines above, in order
  ! to preserve bitwise reproducing older results.  so we only perform the extra 
  ! multiply here in the case that we save the vfrict_impl for diagnostics.  
  if (id_vfrict_impl_u > 0) call diagnose_3d_u(Time, Grd, id_vfrict_impl_u,                       &
                                    Thickness%rho_dzu(:,:,:,taup1)*Velocity%vfrict_impl(:,:,:,1))
  if (id_vfrict_impl_v > 0) call diagnose_3d_u(Time, Grd, id_vfrict_impl_v,                       &
                                    Thickness%rho_dzu(:,:,:,taup1)*Velocity%vfrict_impl(:,:,:,2))

  else !diagnostics without the bug fix
     call diagnose_3d_u(Time, Grd, id_vfrict_impl_u, Velocity%vfrict_impl(:,:,:,1))
     call diagnose_3d_u(Time, Grd, id_vfrict_impl_v, Velocity%vfrict_impl(:,:,:,2))
  endif

  ! take n=1 as representative of the form drag contribution to viscosity
  if (id_visc_cbu > 0 .and. aidif == 1.0) then
     call diagnose_3d_u(Time, Grd, id_visc_cbu, visc_cbu(:,:,:)+visc_cbu_form_drag(:,:,:,1))
  endif 

end subroutine vert_friction_implicit_bgrid
! </SUBROUTINE> NAME="vert_friction_implicit_bgrid"


!#######################################################################
! <SUBROUTINE NAME="vert_friction_implicit_cgrid">
!
! <DESCRIPTION>
! Contributions to thickness weighted and density weighted acceleration 
! from implicit vertical friction. 
!
! Assume that the horizontal grid is C-grid.  
!
! Note that smf and bmf have units N/m^2.  These are the natural units 
! for surface stress.  To include these stresses as boundary terms in the 
! call to invtri, it is necessary to use vertical viscosities with units
! (kg/m^3)*(m2^/s) = N/m^2.  This is achieved by multiplying viscosity
! by rho0 when sent to invtri.  For depth-like vertical coordinates, this
! rho0 factor is cancelled exactly by the rho0 in rho_dzu.  For pressure-like
! vertical coordinates, the rho0*viscosity introduces a negligible
! change in the vertical viscosity that is well within uncertainty
! in this coefficient; a presumably more accurate approach is rho*viscosity.
!
! Include visc_cbu_form_drag to each of the velocity components 
! vertical friction. Do not worry about averaging visc_cbu_form_drag
! to u,v grid cell faces.
!
! Note: if try to merge this routine with vert_friction_bgrid,
! some machines and compilers will change bits by the mere 
! introduction of extra if-test logic into the calculation.
! So we define the separate routines to maintain bit-wise 
! agreement with older results, with bit-wise agreement a 
! useful means to check for errors as the model evolves.
!
! Note: use visc_cbt for both C-grid velocity components,
! even though the velocity components sit at different 
! sides of the tracer cell.  This choice is for simplicity.
! It also acknowledges that the alternative of introducing 
! distinct visc_cbt_u and visc_cbt_v would presume knowledge
! of subgrid scale features that we do not have. So we choose
! to use visc_cbt for both u,v C-grid velocity components.
!
! </DESCRIPTION>
!
subroutine vert_friction_implicit_cgrid (visc_cbt, visc_cbu_form_drag, Time, Thickness, Adv_vel, Velocity)

  real, dimension(isd:,jsd:,:),   intent(in)    :: visc_cbt
  real, dimension(isd:,jsd:,:,:), intent(in)    :: visc_cbu_form_drag
  type(ocean_time_type),          intent(in)    :: Time  
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_adv_vel_type),       intent(in)    :: Adv_vel
  type(ocean_velocity_type),      intent(inout) :: Velocity

  integer :: taum1, tau, taup1
  integer :: i, j, k, n

  integer :: stdoutunit 
  stdoutunit=stdout() 

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  Velocity%vfrict_impl = 0.0
  wrk1_v               = 0.0
  wrk1                 = 0.0
  wrk1_2d              = 0.0

  ! acceleration due to time-implicit vertical friction 
  if (aidif /= 0.0) then  

      ! construct updated velocity due to time-explicit contributions, sans barotropic forcing
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,1) = Grd%tmasken(i,j,k,1)*                            &
                    (Adv_vel%uhrho_et(i,j,k) + dtime_u*Velocity%accel(i,j,k,1))/  &
                    (Thickness%rho_dzten(i,j,k,1) + epsln)
               wrk1_v(i,j,k,2) = Grd%tmasken(i,j,k,2)*                            &
                    (Adv_vel%vhrho_nt(i,j,k) + dtime_u*Velocity%accel(i,j,k,2))/  &
                    (Thickness%rho_dzten(i,j,k,2) + epsln)
            enddo
         enddo
      enddo

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               Velocity%vfrict_impl(i,j,k,1) = wrk1_v(i,j,k,1)
               Velocity%vfrict_impl(i,j,k,2) = wrk1_v(i,j,k,2)
               wrk1(i,j,k)                   = Thickness%rho_dzt(i,j,k,taup1)
            enddo
         enddo
      enddo

      ! invert tridiagonal for vertical friction. 
      ! vertical grid spacing taken from T-cell. 
      if(Velocity%bmf_implicit) then 

          call invtri_bmf (wrk1_v(:,:,:,1), Velocity%smf_cgrid(:,:,1), Velocity%gamma(:,:),&
               rho0*(visc_cbt(:,:,:)+visc_cbu_form_drag(:,:,:,1)), dtime_u, Grd%kmt,       &
               Grd%tmask, wrk1, Thickness%dzwt, aidif, nk)
          call invtri_bmf (wrk1_v(:,:,:,2), Velocity%smf_cgrid(:,:,2), Velocity%gamma(:,:),&
               rho0*(visc_cbt(:,:,:)+visc_cbu_form_drag(:,:,:,2)), dtime_u, Grd%kmt,       &
               Grd%tmask, wrk1, Thickness%dzwt, aidif, nk)

      else 

          call invtri (wrk1_v(:,:,:,1), Velocity%smf_cgrid(:,:,1), Velocity%bmf(:,:,1),&
               rho0*(visc_cbt(:,:,:)+visc_cbu_form_drag(:,:,:,1)), dtime_u, Grd%kmt,   &
               Grd%tmask, wrk1, Thickness%dzwt, aidif, nk)
          call invtri (wrk1_v(:,:,:,2), Velocity%smf_cgrid(:,:,2), Velocity%bmf(:,:,2),&
               rho0*(visc_cbt(:,:,:)+visc_cbu_form_drag(:,:,:,2)), dtime_u, Grd%kmt,   &
               Grd%tmask, wrk1, Thickness%dzwt, aidif, nk)
      endif

      ! compute time-implicit tendency for diagnostics 
      ! update thickness weighted and density weighted 
      ! acceleration due to implicit vertical friction.
      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%vfrict_impl(i,j,k,n) = &
                   dtime_ur*(wrk1_v(i,j,k,n)-Velocity%vfrict_impl(i,j,k,n))
                  Velocity%accel(i,j,k,n)       = Velocity%accel(i,j,k,n)  &
                  + Grd%tmasken(i,j,k,n)*Thickness%rho_dzten(i,j,k,n)*Velocity%vfrict_impl(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

      if (debug_this_module) then
          write(stdoutunit,*) ' ' 
          write(stdoutunit,*) 'From ocean_vert_mix_mod: chksums after implicit vert frict'
          call write_timestamp(Time%model_time)
          call write_chksum_3d('accel(1)', Velocity%accel(COMP,:,1))
          call write_chksum_3d('accel(2)', Velocity%accel(COMP,:,2))
          call write_chksum_3d('vel(1)', wrk1_v(COMP,:,1))
          call write_chksum_3d('vel(2)', wrk1_v(COMP,:,2))
      endif


  endif  !endif for aidif /= 0

  if (id_vfrict_impl_u > 0) call diagnose_3d(Time, Grd, id_vfrict_impl_u,                       &
                                    Thickness%rho_dzten(:,:,:,1)*Velocity%vfrict_impl(:,:,:,1))
  if (id_vfrict_impl_v > 0) call diagnose_3d(Time, Grd, id_vfrict_impl_v,                       &
                                    Thickness%rho_dzten(:,:,:,2)*Velocity%vfrict_impl(:,:,:,2))

  ! take n=1 as representative of the form drag contribution to viscosity
  if (id_visc_cbu > 0 .and. aidif == 1.0) then
       wrk1(:,:,:) = 0.0
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                wrk1(i,j,k)  = onefourth*Grd%umask(i,j,k)* &
                (visc_cbt(i,j,k) + visc_cbt(i+1,j,k) + visc_cbt(i,j+1,k) + visc_cbt(i+1,j+1,k))
             enddo
         enddo
      enddo
      call diagnose_3d_u(Time, Grd, id_visc_cbu, wrk1(:,:,:))
  endif 
  if (id_visc_cbt > 0 .and. aidif == 1.0) then
     call diagnose_3d(Time, Grd, id_visc_cbt, visc_cbt(:,:,:))
  endif 


end subroutine vert_friction_implicit_cgrid
! </SUBROUTINE> NAME="vert_friction_implicit_cgrid"



!#######################################################################
! <FUNCTION NAME="on_comp_domain">
!
! <DESCRIPTION>
! Determine if the point is in comp-domain for the processor.
! </DESCRIPTION>
!
function on_comp_domain(ntable)

  integer, intent(in) :: ntable
  logical             :: on_comp_domain

  if(isc+Dom%ioff <= itable(ntable) .and. itable(ntable) <= iec+Dom%ioff .and. &
     jsc+Dom%joff <= jtable(ntable) .and. jtable(ntable) <= jec+Dom%joff) then
     on_comp_domain = .true.
  else 
     on_comp_domain = .false.  
  endif

end function on_comp_domain
! </FUNCTION> NAME="on_comp_domain"


!#######################################################################
! <FUNCTION NAME="invcosh">
!
! <DESCRIPTION>
! Inverse cosh function.  Argument must be >=1. 
! </DESCRIPTION>
!
function invcosh(a)

  real, intent(in) :: a 
  real             :: invcosh

  invcosh = log(a + sqrt(a**2-1.0))

end function invcosh
! </FUNCTION> NAME="invcosh"



!#######################################################################
! <SUBROUTINE NAME="vmix_min_dissipation">
! <DESCRIPTION>
!  Impose a floor to the dissipation arising from vertical tracer diffusion. 
! </DESCRIPTION>
subroutine vmix_min_dissipation(Time, Dens, diff_cbt, visc_cbu, visc_cbt)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_density_type),       intent(in)    :: Dens
  real, dimension(isd:,jsd:,:,:), intent(inout) :: diff_cbt
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbt

  real    :: tmp, absbvfreq, mindiss
  real    :: rho_N2_prev
  integer :: i,j,k,kp1,kbot,m

  wrk1(:,:,:)     = 0.0
  wrk2(:,:,:)     = 0.0
  wrk3(:,:,:)     = 0.0
  wrk1_v(:,:,:,:) = 0.0

  ! compute BV frequency in same way as in ocean_vert_tidal.F90. 
  ! absolute(rho*N^2) computed from ocean_density module calculation. 
  ! use the value at T-cell centre as this produces a smoother and
  ! better behaved bvfreq near the bottom, than does drhodz_wt. 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk2(i,j,k) = max(0.0,-grav*Dens%drhodz_zt(i,j,k)*Grd%tmask(i,j,k)) 
        enddo
     enddo
  enddo

  ! smooth rho*N2 in the vertical using a 1-2-1 filter
  if (smooth_rho_N2) then
      do m=1,num_121_passes
         do j=jsd,jed
            do i=isd,ied
               rho_N2_prev = onefourth*wrk2(i,j,1)
               kbot=Grd%kmt(i,j)
               if (kbot>3) then
                   do k=2,kbot-2
                      tmp         = wrk2(i,j,k)
                      wrk2(i,j,k) = rho_N2_prev + onehalf*wrk2(i,j,k) + onefourth*wrk2(i,j,k+1)
                      rho_N2_prev = onefourth*tmp
                   enddo
               endif
            enddo
         enddo
      enddo
  endif

  ! compute floor to diffusivity based on floor to dissipation
  do k=1,nk-1
     kp1=k+1
     do j=jsd,jed
        do i=isd,ied

           ! absolute buoyancy frequency |N| (1/sec)
           absbvfreq = sqrt(rho0r*wrk2(i,j,k))

           ! compute floor to dissipation  (W/m^3) 
           mindiss = Grd%tmask(i,j,k)*(vmix_min_diss_const + vmix_min_diss_bvfreq_scale*absbvfreq)
           wrk3(i,j,k) = mindiss

           ! set floor to diffusivity
           wrk1_v(i,j,k,1) = diff_cbt(i,j,k,1)
           wrk1_v(i,j,k,2) = diff_cbt(i,j,k,2)
           wrk1(i,j,k) = Grd%tmask(i,j,kp1) &
                         *vmix_min_diss_flux_ri_max*mindiss*rho0r/(absbvfreq**2 + omega_earth2) 

        enddo
     enddo
  enddo

  ! increase the vertical viscosity based on boosted diffusivity 
  ! and assumed unit Prandtl number.  
  wrk2(:,:,:) = 0.0
  do k=1,nk-1
     kp1=k+1
     do j=jsc,jec
        do i=isc,iec
           diff_cbt(i,j,k,1) = max(diff_cbt(i,j,k,1),wrk1(i,j,k))
           diff_cbt(i,j,k,2) = max(diff_cbt(i,j,k,2),wrk1(i,j,k))
           visc_cbt(i,j,k)   = max(visc_cbt(i,j,k)  ,wrk1(i,j,k))
           tmp  = Grd%umask(i,j,kp1)*onefourth  &
                 *(wrk1(i,j,k)+wrk1(i+1,j,k)+wrk1(i,j+1,k)+wrk1(i+1,j+1,k))
           wrk2(i,j,k)     = visc_cbu(i,j,k)
           visc_cbu(i,j,k) = max(visc_cbu(i,j,k),tmp)
        enddo
     enddo
  enddo


  ! diagnostics 
  call diagnose_3d(Time, Grd, id_diff_cbt_vmix_min, wrk1(:,:,:))
  call diagnose_3d(Time, Grd, id_vmix_min_diss, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_diff_cbt_t_before_min, wrk1_v(:,:,:,1))
  call diagnose_3d(TIme, Grd, id_diff_cbt_s_before_min, wrk1_v(:,:,:,2))
  call diagnose_3d_u(Time, Grd, id_visc_cbu_before_min, wrk2(:,:,:))

end subroutine vmix_min_dissipation
! </SUBROUTINE> NAME="vmix_min_dissipation"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag">
!
! <DESCRIPTION>
! Diagnose watermass transformation diagnostics for 
! --sbc
! --bbc
! --vmix = diff_cbt + K33-implicit neutral diffusion
! Estimations of the contributions from diff_cbt alone, and 
! from K33-implicit alone, are computed in vert_diffuse_watermass_diag.
!
! The diagnostic for vdiffuse computes all processes in one 
! invtri call, which is actually how the model prognostically
! updates tracers from vertical diffusion.  
! 
! The sum of the sbc + bbc + vmix should equal to the 
! full vdiffuse diagnostic, so that, for example, 
!
! neut_rho_vdiffuse = neut_rho_sbc + neut_rho_bbc + neut_rho_vmix 
!
! Additionally, we should have 
!
! neut_rho_vmix = neut_rho_diff_cbt + neut_rho_k33
!
! with the terms neut_rho_diff_cbt and neut_rho_k33 computed 
! in routine vert_diffuse_watermass_diag.
!
! watermass_diag is called prior to implicit update of the tracer
! fields, so that the initial taup1 value contains only explicit 
! in-time tendencies. The incremental tendencies are diagnosed
! in this routine by various calls to invtri using same methods 
! as for the prognostic calculation. 
! 
! This routine requires the logical compute_watermass_diag=.true.,
! which is determined inside watermass_diag_init.  
!
! </DESCRIPTION>
!
subroutine watermass_diag(Time, T_prog, Dens, Thickness, diff_cbt)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_prog_tracer_type),   intent(in) :: T_prog(:)
  type(ocean_density_type),       intent(in) :: Dens
  type(ocean_thickness_type),     intent(in) :: Thickness
  real, dimension(isd:,jsd:,:,:), intent(in) :: diff_cbt

  integer :: i,j,k,taup1

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_vert_mix (watermass_diag): module needs initialization')
  endif 

  if(.not. compute_watermass_diag) return

  taup1 = Time%taup1


  !-------------------------------------------------------------------------------------
  ! compute temp and salinity tendencies due to bottom boundary fluxes.

  ! for holding the pre-implicit time step version of temp and salt.
  wrk2_v(:,:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk2_v(i,j,k,1) = T_prog(index_temp)%field(i,j,k,taup1)
           wrk2_v(i,j,k,2) = T_prog(index_salt)%field(i,j,k,taup1)
        enddo
     enddo
  enddo

  ! for holding the time-tendencies 
  wrk1_v(:,:,:,:) = 0.0

  ! save temp prior to implicit time update
  wrk2_vmix(:,:,:) = 0.0 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk2_vmix(i,j,k) = wrk2_v(i,j,k,1)
        enddo
     enddo
  enddo

  ! bottom temp flux 
  do j=jsc,jec
     do i=isc,iec
        wrk2_2d(i,j) = T_prog(index_temp)%btf(i,j)
     enddo
  enddo

  ! no diffusivity or surface flux for this diagnostic calculation 
  wrk2(:,:,:)  = 0.0
  wrk1_2d(:,:) = 0.0

  call invtri (wrk2_vmix(:,:,:), wrk1_2d(:,:),     & 
       wrk2_2d(:,:), wrk2(:,:,:), dtime_t, Grd%kmt,&
       Grd%tmask, Thickness%rho_dzt(:,:,:,taup1), Thickness%dzwt, aidif, nk)  

  ! fill the tendency array and update temp
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec  
           wrk1_v(i,j,k,1) = dtime_tr*Thickness%rho_dzt(i,j,k,taup1) &
                             *Dens%drhodT(i,j,k)*(wrk2_vmix(i,j,k)-wrk2_v(i,j,k,1))
           wrk2_v(i,j,k,1) = wrk2_vmix(i,j,k)
        enddo
     enddo
  enddo

  ! fill diagnostic arrays 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec  
           wrk1(i,j,k) = wrk1_v(i,j,k,1)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_bbc_temp, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_bbc_temp, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_bbc_temp, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_bbc_temp_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_bbc_temp_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_bbc_temp_on_nrho, wrk4)

  !-----------------------------------------------------------------------
  ! compute temp and salinity tendencies due to surface boundary fluxes.

  ! we already have the pre-implicit time step version of temp and salt
  ! held in k=1 from the previous step based on the bbc calculation.  

  ! for holding time tendencies of temp and salinity 
  wrk1_v(:,:,:,:) = 0.0

  ! save temp prior to implicit time update
  wrk2_vmix(:,:,:) = 0.0 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk2_vmix(i,j,k) = wrk2_v(i,j,k,1)
        enddo
     enddo
  enddo

  ! surface temp flux 
  do j=jsc,jec
     do i=isc,iec
        wrk1_2d(i,j) = T_prog(index_temp)%stf(i,j)
     enddo
  enddo

  ! no diffusivity or bottom flux for this diagnostic calculation 
  wrk2_2d(:,:) = 0.0
  wrk2(:,:,:)  = 0.0 

  call invtri (wrk2_vmix(:,:,:), wrk1_2d(:,:),     & 
       wrk2_2d(:,:), wrk2(:,:,:), dtime_t, Grd%kmt,&
       Grd%tmask, Thickness%rho_dzt(:,:,:,taup1), Thickness%dzwt, aidif, nk)  

  ! fill the tendency array for temp update, and increment the value for temp
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec  
           wrk1_v(i,j,k,1) = dtime_tr*Thickness%rho_dzt(i,j,k,taup1) &
                             *Dens%drhodT(i,j,k)*(wrk2_vmix(i,j,k)-wrk2_v(i,j,k,1))
           wrk2_v(i,j,k,1) = wrk2_vmix(i,j,k)
        enddo
     enddo
  enddo

  ! save salinity prior to implicit time update
  wrk2_vmix(:,:,:) = 0.0 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk2_vmix(i,j,k) = wrk2_v(i,j,k,2)
        enddo
     enddo
  enddo

  ! surface salt flux 
  do j=jsc,jec
     do i=isc,iec
        wrk1_2d(i,j) = T_prog(index_salt)%stf(i,j)
     enddo
  enddo

  ! no diffusivity or bottom flux for this diagnostic calculation 
  wrk2_2d(:,:) = 0.0
  wrk2(:,:,:)  = 0.0

  call invtri (wrk2_vmix(:,:,:), wrk1_2d(:,:),     & 
       wrk2_2d(:,:), wrk2(:,:,:), dtime_t, Grd%kmt,&
       Grd%tmask, Thickness%rho_dzt(:,:,:,taup1), Thickness%dzwt, aidif, nk)  

  ! fill the tendency array for salinity update, and increment the value for salinity 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec  
           wrk1_v(i,j,k,2) = dtime_tr*Thickness%rho_dzt(i,j,k,taup1) &
                             *Dens%drhodS(i,j,k)*(wrk2_vmix(i,j,k)-wrk2_v(i,j,k,2))
           wrk2_v(i,j,k,2) = wrk2_vmix(i,j,k)
        enddo
     enddo
  enddo

  ! diagnose effects from surface temperature flux 
  if(id_tform_rho_sbc_temp_on_nrho > 0 .or. id_tform_rho_sbc_temp > 0 .or. &
     id_wdian_rho_sbc_temp_on_nrho > 0 .or. id_wdian_rho_sbc_temp > 0 .or. &
     id_neut_rho_sbc_temp_on_nrho  > 0 .or. id_neut_rho_sbc_temp  > 0 ) then 

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec  
               wrk1(i,j,k) = wrk1_v(i,j,k,1)
               wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
               wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_rho_sbc_temp, wrk2(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_rho_sbc_temp, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_rho_sbc_temp, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_rho_sbc_temp_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, id_wdian_rho_sbc_temp_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, id_tform_rho_sbc_temp_on_nrho, wrk4)
  endif


  ! diagnose effects from surface salt flux 
  if(id_tform_rho_sbc_salt_on_nrho > 0 .or. id_tform_rho_sbc_salt > 0 .or. &
     id_wdian_rho_sbc_salt_on_nrho > 0 .or. id_wdian_rho_sbc_salt > 0 .or. &
     id_neut_rho_sbc_salt_on_nrho  > 0 .or. id_neut_rho_sbc_salt  > 0 ) then 

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec  
               wrk1(i,j,k) = wrk1_v(i,j,k,2)
               wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
               wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_rho_sbc_salt, wrk2(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_rho_sbc_salt, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_rho_sbc_salt, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_rho_sbc_salt_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, id_wdian_rho_sbc_salt_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, id_tform_rho_sbc_salt_on_nrho, wrk4)
  endif ! endif for rho_sbc_salt diagnostics 


  ! diagnose effects from surface temp plus salt flux 
  if(id_tform_rho_sbc_on_nrho > 0 .or. id_tform_rho_sbc > 0 .or. &
     id_neut_rho_sbc_on_nrho  > 0 .or. id_neut_rho_sbc  > 0 .or. &
     id_wdian_rho_sbc_on_nrho > 0 .or. id_wdian_rho_sbc > 0 ) then 

      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec  
               wrk1(i,j,k) = wrk1_v(i,j,k,1) + wrk1_v(i,j,k,2)
               wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
               wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_rho_sbc, wrk2(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_rho_sbc, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_rho_sbc, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_rho_sbc_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, id_wdian_rho_sbc_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, id_tform_rho_sbc_on_nrho, wrk4)
  endif ! endif for rho_sbc diagnostics 


  !-------------------------------------------------------------------------------------
  ! effects from vertical diffusion from vmix = diff_cbt + K33

  ! for holding the time tendencies of temp and salt 
  wrk1_v(:,:,:,:) = 0.0

  ! save temp based on increment from previous process, and hold diffusivity in an array
  wrk2_vmix(:,:,:) = 0.0 
  wrk2(:,:,:)      = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk2_vmix(i,j,k) = wrk2_v(i,j,k,1)
           wrk2(i,j,k)      = rho0*diff_cbt(i,j,k,1) + T_prog(index_temp)%K33_implicit(i,j,k)
        enddo
     enddo
  enddo

  ! surface and bottom boundary fluxes are zero for this process 
  wrk1_2d(:,:) = 0.0
  wrk2_2d(:,:) = 0.0

  call invtri (wrk2_vmix(:,:,:), wrk1_2d(:,:),     & 
       wrk2_2d(:,:), wrk2(:,:,:), dtime_t, Grd%kmt,&
       Grd%tmask, Thickness%rho_dzt(:,:,:,taup1), Thickness%dzwt, aidif, nk)  

  ! compute temp tendency
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec  
           wrk1_v(i,j,k,1) = dtime_tr*Thickness%rho_dzt(i,j,k,taup1) &
                             *Dens%drhodT(i,j,k)*(wrk2_vmix(i,j,k)-wrk2_v(i,j,k,1))
           wrk2_v(i,j,k,1) = wrk2_vmix(i,j,k)
        enddo
     enddo
  enddo

  ! save salinity based on increment from previous process, and hold diffusivity in an array
  wrk2_vmix(:,:,:) = 0.0 
  wrk2(:,:,:)      = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk2_vmix(i,j,k) = wrk2_v(i,j,k,2)
           wrk2(i,j,k)      = rho0*diff_cbt(i,j,k,2) + T_prog(index_salt)%K33_implicit(i,j,k)
        enddo
     enddo
  enddo

  ! surface and bottom boundary fluxes are zero for this process 
  wrk1_2d(:,:) = 0.0
  wrk2_2d(:,:) = 0.0

  call invtri (wrk2_vmix(:,:,:), wrk1_2d(:,:),     & 
       wrk2_2d(:,:), wrk2(:,:,:), dtime_t, Grd%kmt,&
       Grd%tmask, Thickness%rho_dzt(:,:,:,taup1), Thickness%dzwt, aidif, nk)  

  ! compute salinity tendency
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec  
           wrk1_v(i,j,k,2) = dtime_tr*Thickness%rho_dzt(i,j,k,taup1) &
                             *Dens%drhodS(i,j,k)*(wrk2_vmix(i,j,k)-wrk2_v(i,j,k,2))
           wrk2_v(i,j,k,2) = wrk2_vmix(i,j,k)
        enddo
     enddo
  enddo

  ! fill diagnostic arrays 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = wrk1_v(i,j,k,1) + wrk1_v(i,j,k,2)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_vmix, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_vmix, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_vmix, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_vmix_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_vmix_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_vmix_on_nrho, wrk4)

 
  !-------------------------------------------------------------------------------------
  ! diagnose effects from full vertical diffusion operator: 
  ! vdiffuse = diff_cbt + K33 + stf + btf 
  if(id_neut_temp_vdiffuse        > 0 .or.  id_wdian_temp_vdiffuse        > 0  .or.  &
     id_neut_salt_vdiffuse        > 0 .or.  id_wdian_salt_vdiffuse        > 0  .or.  &
     id_neut_rho_vdiffuse         > 0 .or.  id_wdian_rho_vdiffuse         > 0  .or.  &
     id_neut_rho_vdiffuse_on_nrho > 0 .or.  id_wdian_rho_vdiffuse_on_nrho > 0  .or.  &
     id_tform_rho_vdiffuse        > 0 .or.  id_tform_rho_vdiffuse_on_nrho > 0) then 

      ! for holding the pre-implicit time step version of temp and salt.
      wrk2_v(:,:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk2_v(i,j,k,1) = T_prog(index_temp)%field(i,j,k,taup1)
               wrk2_v(i,j,k,2) = T_prog(index_salt)%field(i,j,k,taup1)
            enddo
         enddo
      enddo

      ! for holding the time tendencies of temp and salt 
      wrk1_v(:,:,:,:) = 0.0

      ! save temp prior to implicit time update, and hold diffusivity in an array
      wrk2_vmix(:,:,:) = 0.0 
      wrk2(:,:,:)      = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk2_vmix(i,j,k) = wrk2_v(i,j,k,1)
               wrk2(i,j,k)      = rho0*diff_cbt(i,j,k,1) + T_prog(index_temp)%K33_implicit(i,j,k)
            enddo
         enddo
      enddo

      ! surface and bottom boundary fluxes 
      wrk1_2d(:,:) = 0.0
      wrk2_2d(:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec
            wrk1_2d(i,j) = T_prog(index_temp)%stf(i,j)
            wrk2_2d(i,j) = T_prog(index_temp)%btf(i,j)
         enddo
      enddo

      call invtri (wrk2_vmix(:,:,:), wrk1_2d(:,:),     & 
           wrk2_2d(:,:), wrk2(:,:,:), dtime_t, Grd%kmt,&
           Grd%tmask, Thickness%rho_dzt(:,:,:,taup1), Thickness%dzwt, aidif, nk)  

      ! compute temp tendency
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec  
               wrk1_v(i,j,k,1) = dtime_tr*Thickness%rho_dzt(i,j,k,taup1) &
                                 *Dens%drhodT(i,j,k)*(wrk2_vmix(i,j,k)-wrk2_v(i,j,k,1))
            enddo
         enddo
      enddo

      ! save salinity prior to implicit time update, and hold diffusivity in an array
      wrk2_vmix(:,:,:) = 0.0 
      wrk2(:,:,:)      = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk2_vmix(i,j,k) = wrk2_v(i,j,k,2)
               wrk2(i,j,k)      = rho0*diff_cbt(i,j,k,2) + T_prog(index_salt)%K33_implicit(i,j,k)
            enddo
         enddo
      enddo

      ! surface and bottom boundary fluxes 
      wrk1_2d(:,:) = 0.0
      wrk2_2d(:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec
            wrk1_2d(i,j) = T_prog(index_salt)%stf(i,j)
            wrk2_2d(i,j) = T_prog(index_salt)%btf(i,j)
         enddo
      enddo

      call invtri (wrk2_vmix(:,:,:), wrk1_2d(:,:),     & 
           wrk2_2d(:,:), wrk2(:,:,:), dtime_t, Grd%kmt,&
           Grd%tmask, Thickness%rho_dzt(:,:,:,taup1), Thickness%dzwt, aidif, nk)  

      ! compute salinity tendency
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec  
               wrk1_v(i,j,k,2) = dtime_tr*Thickness%rho_dzt(i,j,k,taup1) &
                                 *Dens%drhodS(i,j,k)*(wrk2_vmix(i,j,k)-wrk2_v(i,j,k,2))
            enddo
         enddo
      enddo


      ! fill temp diagnostics 
      if(id_neut_temp_vdiffuse > 0) then 
          wrk1(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   wrk1(i,j,k) = wrk1_v(i,j,k,1)*Dens%rho_dztr_tau(i,j,k)
                enddo
             enddo
          enddo
          call diagnose_3d(Time, Grd, id_neut_temp_vdiffuse, wrk1(:,:,:))
      endif
      if(id_wdian_temp_vdiffuse > 0) then 
          wrk1(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   wrk1(i,j,k) = wrk1_v(i,j,k,1)*Dens%stratification_factor(i,j,k)*Dens%rho_dztr_tau(i,j,k)
                enddo
             enddo
          enddo
          call diagnose_3d(Time, Grd, id_wdian_temp_vdiffuse, wrk1(:,:,:))
      endif


      ! fill salinity diagnostics 
      if(id_neut_salt_vdiffuse > 0) then 
          wrk1(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   wrk1(i,j,k) = wrk1_v(i,j,k,2)*Dens%rho_dztr_tau(i,j,k)
                enddo
             enddo
          enddo
          call diagnose_3d(Time, Grd, id_neut_salt_vdiffuse, wrk1(:,:,:))
      endif
      if(id_wdian_salt_vdiffuse > 0) then 
          wrk1(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   wrk1(i,j,k) = wrk1_v(i,j,k,2)*Dens%stratification_factor(i,j,k)*Dens%rho_dztr_tau(i,j,k)
                enddo
             enddo
          enddo
          call diagnose_3d(Time, Grd, id_wdian_salt_vdiffuse, wrk1(:,:,:))
      endif


      ! fill rho diagnostic arrays 
      wrk1(:,:,:) = 0.0
      wrk2(:,:,:) = 0.0
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = wrk1_v(i,j,k,1) + wrk1_v(i,j,k,2)
               wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
               wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
               wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
            enddo
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_neut_rho_vdiffuse, wrk2(:,:,:))
      call diagnose_3d(Time, Grd, id_wdian_rho_vdiffuse, wrk3(:,:,:))
      call diagnose_3d(Time, Grd, id_tform_rho_vdiffuse, wrk4(:,:,:))
      call diagnose_3d_rho(Time, Dens, id_neut_rho_vdiffuse_on_nrho, wrk2)
      call diagnose_3d_rho(Time, Dens, id_wdian_rho_vdiffuse_on_nrho, wrk3)
      call diagnose_3d_rho(Time, Dens, id_tform_rho_vdiffuse_on_nrho, wrk4)
  endif   ! endif for diff_cbt + K33 + stf + btf


end subroutine watermass_diag
! </SUBROUTINE>  NAME="watermass_diag"


!#######################################################################
! <SUBROUTINE NAME="vert_diffuse_implicit_diag">
!
! <DESCRIPTION>
! Diagnose contributions from time-implicit vertical diffusion.
!
! We perform the diagnostics using the taup1 value of the tracer 
! concentration, obtained after performing the invtri step.  
! The fluxes are diagnosed as if we are performing an explicit-time
! update, but using the taup1 values of the tracer concentrations.  
!
! This diagnostic suffers from time truncation errors relative to the 
! prognostic calculation, since the prognostic time update is performed using
! an invtri calculation.  However, the errors are small.  The alternative 
! method, which is to separate the combined invtri step into individual 
! physical mixing processes is not an option, since this step is not equivalent
! algorithmically to the single invtri mixing step. 
! It is for this reason that we perform the diagnostic step using this 
! "explicit flux computed with time-implicit computed concentration" approach.  
!
! In contrast to the interior mixing, the boundary fluxes can be diagnostically
! split from the invtri step, as these fluxes are determined prior to the 
! vert_diffuse_implicit routine. Hence, the boundary flux contributions can
! be diagnosed either by calling an invtri step passing just the stf and btf 
! terms, or by using the even simpler methods in this subroutine.     
!
! </DESCRIPTION>
!
subroutine vert_diffuse_implicit_diag(Time, Thickness, T_prog, diff_cbt, work, n)

  type(ocean_time_type),          intent(in) :: Time  
  type(ocean_thickness_type),     intent(in) :: Thickness
  type(ocean_prog_tracer_type),   intent(in) :: T_prog(:)
  real, dimension(isd:,jsd:,:,:), intent(in) :: diff_cbt
  real, dimension(isd:,jsd:,:),   intent(in) :: work
  integer,                        intent(in) :: n

  integer :: taup1, tau
  integer :: i, j, k, kp, kp1, nmix

  real :: dTdz, diffusivity 

  tau   = Time%tau
  taup1 = Time%taup1
 
  if (n==index_salt) then
      nmix=2
  else
      nmix=1
  endif

  ! diagnose upgradient vertical diffusive flux
  if (id_zflux_diff(n) > 0) then
     flux_z(:,:,:) = 0.0
     do k=1,nk-1
        kp = k+1
        do j=jsc,jec
           do i=isc,iec
              flux_z(i,j,k) = rho0*Grd%tmask(i,j,kp)*diff_cbt(i,j,k,nmix)                  &
                              *(T_prog(n)%field(i,j,k,taup1)-T_prog(n)%field(i,j,kp,taup1))&
                              /Thickness%dzwt(i,j,k)
          enddo
       enddo
     enddo
     call diagnose_3d(Time, Grd, id_zflux_diff(n), flux_z(:,:,:)*T_prog(n)%conversion)     
  endif


  ! diagnose update due to implicit time stepping
  if (id_vdiffuse_impl(n) > 0) then
      wrk1(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec  
               wrk1(i,j,k) = dtime_tr*Thickness%rho_dzt(i,j,k,taup1) &
                             *(T_prog(n)%field(i,j,k,taup1)-work(i,j,k))
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_vdiffuse_impl(n), wrk1(:,:,:)*T_prog(n)%conversion)
  endif


  ! impacts from surface boundary fluxes 
  if(id_vdiffuse_sbc(n) > 0) then 
      wrk1_2d(:,:) = 0.0
      wrk2_2d(:,:) = 0.0
      wrk1(:,:,:)  = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1_2d(i,j) = T_prog(n)%stf(i,j)
            wrk1(i,j,k)  = Grd%tmask(i,j,1)*(wrk1_2d(i,j)-wrk2_2d(i,j))
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_vdiffuse_sbc(n), wrk1(:,:,:)*T_prog(n)%conversion)
  endif

  ! impacts from bottom boundary fluxes 
  if(id_vdiffuse_bbc(n) > 0) then 
      wrk1_2d(:,:) = 0.0
      wrk2_2d(:,:) = 0.0
      wrk1(:,:,:)  = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk2_2d(i,j) = T_prog(n)%btf(i,j)
            wrk1(i,j,k)  = Grd%tmask(i,j,1)*(wrk1_2d(i,j)-wrk2_2d(i,j))
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_vdiffuse_bbc(n), wrk1(:,:,:)*T_prog(n)%conversion)
  endif


  ! impacts from diff_cbt 
  if(id_vdiffuse_diff_cbt(n) > 0) then 
      wrk1_2d(:,:) = 0.0
      wrk2_2d(:,:) = 0.0
      wrk1(:,:,:)  = 0.0
      do k=1,nk
         kp = min(k+1,nk)
         do j=jsc,jec
            do i=isc,iec
               diffusivity  = rho0*diff_cbt(i,j,k,nmix) 
               wrk2_2d(i,j) = diffusivity*Grd%tmask(i,j,kp)*                    &
                    (T_prog(n)%field(i,j,k,taup1)-T_prog(n)%field(i,j,kp,taup1))&
                    /Thickness%dzwt(i,j,k)
               wrk1(i,j,k)  = wrk1_2d(i,j)-wrk2_2d(i,j)
               wrk1_2d(i,j) = wrk2_2d(i,j)
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_vdiffuse_diff_cbt(n), wrk1(:,:,:)*T_prog(n)%conversion)
  endif

  ! impacts from K33 
  if(id_vdiffuse_k33(n) > 0) then 
      wrk1_2d(:,:) = 0.0
      wrk2_2d(:,:) = 0.0
      wrk1(:,:,:)  = 0.0
      do k=1,nk
         kp = min(k+1,nk)
         do j=jsc,jec
            do i=isc,iec
               diffusivity  = T_prog(n)%K33_implicit(i,j,k) 
               wrk2_2d(i,j) = diffusivity*Grd%tmask(i,j,kp)*                    &
                    (T_prog(n)%field(i,j,k,taup1)-T_prog(n)%field(i,j,kp,taup1))&
                    /Thickness%dzwt(i,j,k)
               wrk1(i,j,k)  = wrk1_2d(i,j)-wrk2_2d(i,j)
               wrk1_2d(i,j) = wrk2_2d(i,j)
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_vdiffuse_k33(n), wrk1(:,:,:)*T_prog(n)%conversion)
  endif


  ! compute dissipation of squared tracer concentration 
  if(id_vdiffuse_diss(n) > 0) then 
      wrk4(:,:,:) = 0.0

      do k=1,1
         kp1 = k+1  
         do j=jsc,jec
            do i=isc,iec
               dTdz = Grd%tmask(i,j,kp1) &
                    *(T_prog(n)%field(i,j,k,taup1)-T_prog(n)%field(i,j,kp1,taup1))/Thickness%dzwt(i,j,k)
               wrk4(i,j,k) = dtime_tr*diff_cbt(i,j,k,nmix)  &
                    *Thickness%rho_dzt(i,j,k,taup1)*Thickness%rho_dzt(i,j,k,tau)*(T_prog(n)%conversion*dTdz)**2
            enddo
         enddo
      enddo
      do k=2,nk-1
         kp1 = k+1
         do j=jsc,jec
            do i=isc,iec
               dTdz = onehalf*Grd%tmask(i,j,k) &
                    *(T_prog(n)%field(i,j,k-1,taup1)-T_prog(n)%field(i,j,k,taup1))/Thickness%dzwt(i,j,k-1)
               wrk4(i,j,k) = dtime_tr*diff_cbt(i,j,k-1,nmix)  &
                    *Thickness%rho_dzt(i,j,k,taup1)*Thickness%rho_dzt(i,j,k,tau)*(T_prog(n)%conversion*dTdz)**2

               dTdz = onehalf*Grd%tmask(i,j,kp1)  &
                    *(T_prog(n)%field(i,j,k,taup1)-T_prog(n)%field(i,j,kp1,taup1))/Thickness%dzwt(i,j,k)
               wrk4(i,j,k) = wrk4(i,j,k) +        &
                    dtime_tr*diff_cbt(i,j,k,nmix) &
                    *Thickness%rho_dzt(i,j,k,taup1)*Thickness%rho_dzt(i,j,k,tau)*(T_prog(n)%conversion*dTdz)**2
            enddo
         enddo
      enddo
      do j=jsc,jec
         do i=isc,iec
            k=Grd%kmt(i,j)
            if(k > 1) then 
                dTdz = Grd%tmask(i,j,k) &
                     *(T_prog(n)%field(i,j,k-1,taup1)-T_prog(n)%field(i,j,k,taup1))/Thickness%dzwt(i,j,k-1)
                wrk4(i,j,k) = dtime_tr*diff_cbt(i,j,k-1,nmix)  &
                     *Thickness%rho_dzt(i,j,k,taup1)*Thickness%rho_dzt(i,j,k,tau)*(T_prog(n)%conversion*dTdz)**2
            endif
         enddo
      enddo

      call diagnose_3d(Time, Grd, id_vdiffuse_diss(n), wrk4(:,:,:))
  endif


end subroutine vert_diffuse_implicit_diag
! </SUBROUTINE> NAME="vert_diffuse_implicit_diag"


!#######################################################################
! <SUBROUTINE NAME="vert_diffuse_watermass_diag">
!
! <DESCRIPTION>
! Diagnose contributions from time-implicit vertical diffusion
! acting on watermasses from diff_cbt and k33.
!
! We perform the diagnostics using the taup1 value of the tracer 
! concentration, obtained after performing the invtri step.  
! The fluxes are diagnosed as if we are performing an explicit-time
! update, but using the taup1 values of the tracer concentrations.  
!
! This diagnostic suffers from time truncation errors relative to the 
! prognostic calculation, since the prognostic time update is performed using
! an invtri calculation.  However, the errors are small.  The alternative 
! method, which is to separate the combined invtri step into individual 
! physical mixing processes is not an option, since this step is not equivalent
! algorithmically to the single invtri mixing step. 
! It is for this reason that we perform the diagnostic step using this 
! "explicit flux computed with time-implicit computed concentration" approach.  
!
! </DESCRIPTION>
!
subroutine vert_diffuse_watermass_diag(Time, Thickness, Dens, T_prog, diff_cbt)

  type(ocean_time_type),          intent(in) :: Time  
  type(ocean_thickness_type),     intent(in) :: Thickness
  type(ocean_density_type),       intent(in) :: Dens
  type(ocean_prog_tracer_type),   intent(in) :: T_prog(:)
  real, dimension(isd:,jsd:,:,:), intent(in) :: diff_cbt

  integer :: taup1, tau
  integer :: i, j, k, kp, nmix

  real,  dimension(isd:ied,jsd:jed) :: eta_tend
  real :: diffusivity 

  if(.not. compute_watermass_diag) return

  tau   = Time%tau
  taup1 = Time%taup1


  !-------------------------------------------------------------------------------------
  ! effects from vertical diffusion using full diffusivity diff_cbt

  ! for holding the time tendencies
  wrk1_v(:,:,:,:) = 0.0

  ! temperature   
  nmix=1
  wrk1_2d(:,:) = 0.0
  wrk2_2d(:,:) = 0.0
  do k=1,nk
     kp = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec
           diffusivity  = rho0*Grd%tmask(i,j,kp)*diff_cbt(i,j,k,nmix) 
           wrk2_2d(i,j) = diffusivity*Grd%tmask(i,j,kp)*                                      &
                (T_prog(index_temp)%field(i,j,k,taup1)-T_prog(index_temp)%field(i,j,kp,taup1))&
                /Thickness%dzwt(i,j,k)
           wrk1_v(i,j,k,1) = Dens%drhodT(i,j,k)*(wrk1_2d(i,j)-wrk2_2d(i,j))
           wrk1_2d(i,j)    = wrk2_2d(i,j)
        enddo
     enddo
  enddo

  ! salinity
  nmix=2
  wrk1_2d(:,:) = 0.0
  wrk2_2d(:,:) = 0.0
  do k=1,nk
     kp = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec
           diffusivity  = rho0*Grd%tmask(i,j,kp)*diff_cbt(i,j,k,nmix) 
           wrk2_2d(i,j) = diffusivity*Grd%tmask(i,j,kp)*                                      &
                (T_prog(index_salt)%field(i,j,k,taup1)-T_prog(index_salt)%field(i,j,kp,taup1))&
                /Thickness%dzwt(i,j,k)
           wrk1_v(i,j,k,2) = Dens%drhodS(i,j,k)*(wrk1_2d(i,j)-wrk2_2d(i,j))
           wrk1_2d(i,j)    = wrk2_2d(i,j)
        enddo
     enddo
  enddo


  ! fill diagnostic arrays for full contribution on rho 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  wrk5(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Grd%tmask(i,j,k)*(wrk1_v(i,j,k,1) + wrk1_v(i,j,k,2))
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
           wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2)   ! for eta_tend
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_diff_cbt, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_diff_cbt, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_diff_cbt, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_diff_cbt_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_diff_cbt_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_diff_cbt_on_nrho, wrk4)

  ! contributions to sea level 
  if(id_eta_tend_diff_cbt_tend > 0 .or. id_eta_tend_diff_cbt_tend_glob > 0) then
      eta_tend(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_diff_cbt_tend, eta_tend(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_diff_cbt_tend_glob, eta_tend, cellarea_r)
  endif


  ! fill diagnostic arrays for temp contribution on rho 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = wrk1_v(i,j,k,1)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_temp_diff_cbt, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_diff_cbt, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_temp_diff_cbt, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_diff_cbt_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_diff_cbt_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_diff_cbt_on_nrho, wrk4)

  ! fill diagnostic arrays for salinity contribution on rho 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = wrk1_v(i,j,k,2)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_salt_diff_cbt, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_diff_cbt, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_salt_diff_cbt, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_diff_cbt_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_diff_cbt_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_diff_cbt_on_nrho, wrk4)


  !-------------------------------------------------------------------------------------
  ! effects from time-implicit piece of K33
  ! this diagnostic is an approximation to the prognostic calculation. 

  ! for holding the time tendencies
  wrk1_v(:,:,:,:) = 0.0

  ! temperature   
  nmix=1
  wrk1_2d(:,:) = 0.0
  wrk2_2d(:,:) = 0.0
  do k=1,nk
     kp = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec
           diffusivity  = T_prog(index_temp)%K33_implicit(i,j,k) 
           wrk2_2d(i,j) = diffusivity*Grd%tmask(i,j,kp)*                                      &
                (T_prog(index_temp)%field(i,j,k,taup1)-T_prog(index_temp)%field(i,j,kp,taup1))&
                /Thickness%dzwt(i,j,k)
           wrk1_v(i,j,k,1) = Dens%drhodT(i,j,k)*(wrk1_2d(i,j)-wrk2_2d(i,j))
           wrk1_2d(i,j)    = wrk2_2d(i,j)
        enddo
     enddo
  enddo

  ! salinity
  nmix=2
  wrk1_2d(:,:) = 0.0
  wrk2_2d(:,:) = 0.0
  do k=1,nk
     kp = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec
           diffusivity  = T_prog(index_salt)%K33_implicit(i,j,k) 
           wrk2_2d(i,j) = diffusivity*Grd%tmask(i,j,kp)*                                      &
                (T_prog(index_salt)%field(i,j,k,taup1)-T_prog(index_salt)%field(i,j,kp,taup1))&
                /Thickness%dzwt(i,j,k)
           wrk1_v(i,j,k,2) = Dens%drhodS(i,j,k)*(wrk1_2d(i,j)-wrk2_2d(i,j))
           wrk1_2d(i,j)    = wrk2_2d(i,j)
        enddo
     enddo
  enddo

  ! fill diagnostic arrays for impacts on rho 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  wrk5(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = wrk1_v(i,j,k,1) + wrk1_v(i,j,k,2)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
           wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2)   ! for eta_tend
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_k33, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_k33, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_k33, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_k33_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_k33_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_k33_on_nrho, wrk4)

  ! contributions to sea level 
  if(id_eta_tend_k33_tend > 0 .or. id_eta_tend_k33_tend_glob > 0) then
      eta_tend(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_k33_tend, eta_tend(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_k33_tend_glob, eta_tend, cellarea_r)
  endif

  ! fill diagnostic arrays for temp contributions
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = wrk1_v(i,j,k,1) 
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_temp_k33, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_k33, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_temp_k33, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_k33_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_k33_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_k33_on_nrho, wrk4)

  ! fill diagnostic arrays for salinity contributions
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = wrk1_v(i,j,k,2) 
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_salt_k33, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_k33, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_salt_k33, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_k33_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_k33_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_k33_on_nrho, wrk4)

  !-------------------------------------------------------------------------------------
  ! effects from vertical diffusion using diff_cbt_wave

  ! for holding the time tendencies
  wrk1_v(:,:,:,:) = 0.0

  ! temperature   
  wrk1_2d(:,:) = 0.0
  wrk2_2d(:,:) = 0.0
  do k=1,nk
     kp = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec
           diffusivity  = rho0*Grd%tmask(i,j,kp)*diff_cbt_wave(i,j,k) 
           wrk2_2d(i,j) = diffusivity*Grd%tmask(i,j,kp)*                                      &
                (T_prog(index_temp)%field(i,j,k,taup1)-T_prog(index_temp)%field(i,j,kp,taup1))&
                /Thickness%dzwt(i,j,k)
           wrk1_v(i,j,k,1) = Dens%drhodT(i,j,k)*(wrk1_2d(i,j)-wrk2_2d(i,j))
           wrk1_2d(i,j)    = wrk2_2d(i,j)
        enddo
     enddo
  enddo

  ! salinity
  wrk1_2d(:,:) = 0.0
  wrk2_2d(:,:) = 0.0
  do k=1,nk
     kp = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec
           diffusivity  = rho0*Grd%tmask(i,j,kp)*diff_cbt_wave(i,j,k) 
           wrk2_2d(i,j) = diffusivity*Grd%tmask(i,j,kp)*                                      &
                (T_prog(index_salt)%field(i,j,k,taup1)-T_prog(index_salt)%field(i,j,kp,taup1))&
                /Thickness%dzwt(i,j,k)
           wrk1_v(i,j,k,2) = Dens%drhodS(i,j,k)*(wrk1_2d(i,j)-wrk2_2d(i,j))
           wrk1_2d(i,j)    = wrk2_2d(i,j)
        enddo
     enddo
  enddo


  ! fill diagnostic arrays for full contribution on rho 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  wrk5(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Grd%tmask(i,j,k)*(wrk1_v(i,j,k,1) + wrk1_v(i,j,k,2))
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
           wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2)   ! for eta_tend
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_diff_wave, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_diff_wave, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_diff_wave, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_diff_wave_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_diff_wave_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_diff_wave_on_nrho, wrk4)

  ! contributions to sea level 
  if(id_eta_tend_diff_wave_tend > 0 .or. id_eta_tend_diff_wave_tend_glob > 0) then
      eta_tend(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(TIme, Grd, id_eta_tend_diff_wave_tend, eta_tend(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_diff_wave_tend_glob, eta_tend, cellarea_r)
  endif


  ! fill diagnostic arrays for temp contribution on rho from diff_cbt_wave 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = wrk1_v(i,j,k,1)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_temp_diff_wave, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_diff_wave, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_temp_diff_wave, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_diff_wave_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_diff_wave_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_diff_wave_on_nrho, wrk4)

  ! fill diagnostic arrays for salinity contribution on rho from diff_cbt_wave 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = wrk1_v(i,j,k,2)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_salt_diff_wave, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_diff_wave, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_salt_diff_wave, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_diff_wave_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_diff_wave_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_diff_wave_on_nrho, wrk4)


  !-------------------------------------------------------------------------------------
  ! effects from vertical diffusion using diff_cbt_drag

  ! for holding the time tendencies
  wrk1_v(:,:,:,:) = 0.0

  ! temperature   
  wrk1_2d(:,:) = 0.0
  wrk2_2d(:,:) = 0.0
  do k=1,nk
     kp = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec
           diffusivity  = rho0*Grd%tmask(i,j,kp)*diff_cbt_drag(i,j,k) 
           wrk2_2d(i,j) = diffusivity*Grd%tmask(i,j,kp)*                                      &
                (T_prog(index_temp)%field(i,j,k,taup1)-T_prog(index_temp)%field(i,j,kp,taup1))&
                /Thickness%dzwt(i,j,k)
           wrk1_v(i,j,k,1) = Dens%drhodT(i,j,k)*(wrk1_2d(i,j)-wrk2_2d(i,j))
           wrk1_2d(i,j)    = wrk2_2d(i,j)
        enddo
     enddo
  enddo

  ! salinity
  wrk1_2d(:,:) = 0.0
  wrk2_2d(:,:) = 0.0
  do k=1,nk
     kp = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec
           diffusivity  = rho0*Grd%tmask(i,j,kp)*diff_cbt_drag(i,j,k) 
           wrk2_2d(i,j) = diffusivity*Grd%tmask(i,j,kp)*                                      &
                (T_prog(index_salt)%field(i,j,k,taup1)-T_prog(index_salt)%field(i,j,kp,taup1))&
                /Thickness%dzwt(i,j,k)
           wrk1_v(i,j,k,2) = Dens%drhodS(i,j,k)*(wrk1_2d(i,j)-wrk2_2d(i,j))
           wrk1_2d(i,j)    = wrk2_2d(i,j)
        enddo
     enddo
  enddo


  ! fill diagnostic arrays for full contribution on rho 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  wrk5(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Grd%tmask(i,j,k)*(wrk1_v(i,j,k,1) + wrk1_v(i,j,k,2))
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
           wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2)   ! for eta_tend
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_diff_drag, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_diff_drag, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_diff_drag, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_diff_drag_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_diff_drag_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_diff_drag_on_nrho, wrk4)

  ! contributions to sea level 
  if(id_eta_tend_diff_drag_tend > 0 .or. id_eta_tend_diff_drag_tend_glob > 0) then
      eta_tend(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_diff_drag_tend, eta_tend(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_diff_drag_tend_glob, eta_tend, cellarea_r)
  endif


  ! fill diagnostic arrays for temp contribution on rho from diff_cbt_drag 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = wrk1_v(i,j,k,1)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_temp_diff_drag, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_diff_drag, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_temp_diff_drag, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_diff_drag_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_diff_drag_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_diff_drag_on_nrho, wrk4)

  ! fill diagnostic arrays for salinity contribution on rho from diff_cbt_drag 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = wrk1_v(i,j,k,2)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_salt_diff_drag, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_diff_drag, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_salt_diff_drag, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_diff_drag_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_diff_drag_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_diff_drag_on_nrho, wrk4)

  !-------------------------------------------------------------------------------------
  ! effects from vertical diffusion using diff_cbt_leewave

  ! for holding the time tendencies
  wrk1_v(:,:,:,:) = 0.0

  ! temperature   
  wrk1_2d(:,:) = 0.0
  wrk2_2d(:,:) = 0.0
  do k=1,nk
     kp = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec
           diffusivity  = rho0*Grd%tmask(i,j,kp)*diff_cbt_leewave(i,j,k) 
           wrk2_2d(i,j) = diffusivity*Grd%tmask(i,j,kp)*                                      &
                (T_prog(index_temp)%field(i,j,k,taup1)-T_prog(index_temp)%field(i,j,kp,taup1))&
                /Thickness%dzwt(i,j,k)
           wrk1_v(i,j,k,1) = Dens%drhodT(i,j,k)*(wrk1_2d(i,j)-wrk2_2d(i,j))
           wrk1_2d(i,j)    = wrk2_2d(i,j)
        enddo
     enddo
  enddo

  ! salinity
  wrk1_2d(:,:) = 0.0
  wrk2_2d(:,:) = 0.0
  do k=1,nk
     kp = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec
           diffusivity  = rho0*Grd%tmask(i,j,kp)*diff_cbt_leewave(i,j,k) 
           wrk2_2d(i,j) = diffusivity*Grd%tmask(i,j,kp)*                                      &
                (T_prog(index_salt)%field(i,j,k,taup1)-T_prog(index_salt)%field(i,j,kp,taup1))&
                /Thickness%dzwt(i,j,k)
           wrk1_v(i,j,k,2) = Dens%drhodS(i,j,k)*(wrk1_2d(i,j)-wrk2_2d(i,j))
           wrk1_2d(i,j)    = wrk2_2d(i,j)
        enddo
     enddo
  enddo


  ! fill diagnostic arrays for full contribution on rho 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  wrk5(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Grd%tmask(i,j,k)*(wrk1_v(i,j,k,1) + wrk1_v(i,j,k,2))
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
           wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2)   ! for eta_tend
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_diff_lee, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_diff_lee, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_diff_lee, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_diff_lee_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_diff_lee_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_diff_lee_on_nrho, wrk4)

  ! contributions to sea level 
  if(id_eta_tend_diff_lee_tend > 0 .or. id_eta_tend_diff_lee_tend_glob > 0) then
      eta_tend(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_diff_lee_tend, eta_tend(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_diff_lee_tend_glob, eta_tend, cellarea_r)
  endif


  ! fill diagnostic arrays for temp contribution on rho from diff_cbt_lee 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = wrk1_v(i,j,k,1)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_temp_diff_lee, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_diff_lee, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_temp_diff_lee, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_diff_lee_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_diff_lee_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_diff_lee_on_nrho, wrk4)

  ! fill diagnostic arrays for salinity contribution on rho from diff_cbt_lee 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = wrk1_v(i,j,k,2)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_salt_diff_lee, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_diff_lee, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_salt_diff_lee, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_diff_lee_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_diff_lee_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_diff_lee_on_nrho, wrk4)


  !-------------------------------------------------------------------------------------
  ! effects from vertical diffusion using diff_cbt_back

  ! for holding the time tendencies
  wrk1_v(:,:,:,:) = 0.0

  ! temperature   
  wrk1_2d(:,:) = 0.0
  wrk2_2d(:,:) = 0.0
  do k=1,nk
     kp = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec
           diffusivity  = rho0*Grd%tmask(i,j,kp)*diff_cbt_back(i,j,k) 
           wrk2_2d(i,j) = diffusivity*Grd%tmask(i,j,kp)*                                      &
                (T_prog(index_temp)%field(i,j,k,taup1)-T_prog(index_temp)%field(i,j,kp,taup1))&
                /Thickness%dzwt(i,j,k)
           wrk1_v(i,j,k,1) = Dens%drhodT(i,j,k)*(wrk1_2d(i,j)-wrk2_2d(i,j))
           wrk1_2d(i,j)    = wrk2_2d(i,j)
        enddo
     enddo
  enddo

  ! salinity
  wrk1_2d(:,:) = 0.0
  wrk2_2d(:,:) = 0.0
  do k=1,nk
     kp = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec
           diffusivity  = rho0*Grd%tmask(i,j,kp)*diff_cbt_back(i,j,k) 
           wrk2_2d(i,j) = diffusivity*Grd%tmask(i,j,kp)*                                      &
                (T_prog(index_salt)%field(i,j,k,taup1)-T_prog(index_salt)%field(i,j,kp,taup1))&
                /Thickness%dzwt(i,j,k)
           wrk1_v(i,j,k,2) = Dens%drhodS(i,j,k)*(wrk1_2d(i,j)-wrk2_2d(i,j))
           wrk1_2d(i,j)    = wrk2_2d(i,j)
        enddo
     enddo
  enddo


  ! fill diagnostic arrays for full contribution on rho 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  wrk5(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Grd%tmask(i,j,k)*(wrk1_v(i,j,k,1) + wrk1_v(i,j,k,2))
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
           wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2)   ! for eta_tend
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_diff_back, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_diff_back, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_diff_back, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_diff_back_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_diff_back_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_diff_back_on_nrho, wrk4)

  ! contributions to sea level 
  if(id_eta_tend_diff_back_tend > 0 .or. id_eta_tend_diff_back_tend_glob > 0) then
      eta_tend(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_diff_back_tend, eta_tend(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_diff_back_tend_glob, eta_tend, cellarea_r)
  endif


  ! fill diagnostic arrays for temp contribution on rho from diff_cbt_back 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = wrk1_v(i,j,k,1)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_temp_diff_back, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_diff_back, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_temp_diff_back, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_diff_back_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_diff_back_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_diff_back_on_nrho, wrk4)

  ! fill diagnostic arrays for salinity contribution on rho from diff_cbt_back 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = wrk1_v(i,j,k,2)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_salt_diff_back, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_diff_back, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_salt_diff_back, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_diff_back_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_diff_back_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_diff_back_on_nrho, wrk4)

end subroutine vert_diffuse_watermass_diag
! </SUBROUTINE> NAME="vert_diffuse_watermass_diag"



!#######################################################################
! <SUBROUTINE NAME="ocean_vert_mix_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_vert_mix_restart(time_stamp)
  character(len=*),           intent(in), optional :: time_stamp

  if(MIX_SCHEME==VERTMIX_CHEN) then
     call ocean_vert_chen_restart(time_stamp)
  endif

  if(MIX_SCHEME==VERTMIX_GOTM) then
     call ocean_vert_gotm_restart(time_stamp)
  end if
end subroutine ocean_vert_mix_restart
! </SUBROUTINE> NAME="ocean_vert_mix_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_vert_mix_end">
!
! <DESCRIPTION>
! Chen Scheme requires output of Krauss mixed layer for
! reproducible results.
!
! GOTM requires fields for advection tendency.
!
! </DESCRIPTION>
!
subroutine ocean_vert_mix_end(Time)
  type(ocean_time_type), intent(in) :: Time

  call ocean_vert_mix_restart()

  if(MIX_SCHEME==VERTMIX_CHEN) then
     call ocean_vert_chen_end(Time)
  endif

  if(MIX_SCHEME==VERTMIX_GOTM) then
     call ocean_vert_gotm_end(Time)
  endif


end subroutine ocean_vert_mix_end
! </SUBROUTINE> NAME="ocean_vert_mix_end"

end module ocean_vert_mix_mod




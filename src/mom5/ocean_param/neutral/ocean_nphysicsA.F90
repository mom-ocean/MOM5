module ocean_nphysicsA_mod
#define COMP isc:iec,jsc:jec
#define COMPXL isc-1:iec,jsc:jec
#define COMPYL isc:iec,jsc-1:jec
! 
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<CONTACT EMAIL="Russell.Fiedler@csiro.au"> Russell Fiedler 
!</CONTACT>
!
!<REVIEWER EMAIL="tim.leslie@gmail.com"> Tim Leslie
!</REVIEWER>
!  
!<OVERVIEW>
! Thickness weighted and density weighted time tendency for tracer 
! from Laplacian neutral diffusion + Laplacian GM skew-diffusion.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the cell thickness weighted and density 
! weighted tracer tendency from small angle Laplacian neutral diffusion
! plus Laplacian GM skew-diffusion.  The algorithms are based on 
! MOM4p0d methods.  The fundamental differences from the ocean_nphysicsB
! methods relate to the handling of fluxes near the domain boundaries.  
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, A. Gnanadesikan, R.C. Pacanowski, V. Larichev, 
! J.K. Dukowicz,  and R.D. Smith
! Isoneutral diffusion in a z-coordinate ocean model
! Journal of Physical Oceanography (1998) vol 28 pages 805-830
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies 
! The Gent-McWilliams Skew-flux 
! Journal of Physical Oceanography (1998) vol 28 pages 831-841
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies 
! Fundamentals of Ocean Climate Models (2004)
! Princeton University Press 
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, Elements of MOM (2012)
! </REFERENCE>
!
! <REFERENCE>
! G. Danabasoglu and J. C. McWilliams
! Sensitivity of the global ocean circulation to 
! parameterizations of mesoscale tracer transports
! Journal of Climate (1995) vol 8 pages 2967--2987 
! </REFERENCE>
!
! <REFERENCE>
! Gerdes, Koberle, and Willebrand
! The influence of numerical advection schemes on the results of ocean
! general circulation models, Climate Dynamics (1991), vol. 5, 
! pages 211--226. 
! </REFERENCE>
!
! <NOTE>
! Numerical implementation of the flux components follows the triad 
! approach documented in the references and implemented in MOM2 and MOM3.  
! The MOM algorithm accounts for partial bottom cells and generalized
! orthogonal horizontal coordinates.
! </NOTE> 
!
! <NOTE>
! neutral_physics_simple=.true. requires aredi_equal_agm=.true.
! neutral_physics_simple=.true. results in down-gradient  
! horizontal flux components. This setting reduces the overall cost 
! of the neutral physics scheme, but it is not used at GFDL
! anymore, since we favor methods whereby treatment of GM and Redi
! in the boundary layers are distinct.  
! </NOTE> 
!
! <NOTE> 
! In steep slope regions, neutral diffusive fluxes are tapered to
! zero with the tanh taper of Danabasoglu and McWilliams (1995) or the 
! quadratic scheme of Gerdes, Koberle, and Willebrand.  However, if
! neutral_physics_simple=.false., the GM skew-diffusive fluxes 
! can remain nonzero if have neutral_linear_gm_taper=.true.
! </NOTE> 
!
! </INFO>
!
!<NAMELIST NAME="ocean_nphysicsA_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For printing starting and ending checksums for restarts
!  </DATA> 
!
!  <DATA NAME="use_gm_skew" TYPE="logical">
!  Must be true to use GM skewsion.  Set to false if wish to 
!  incorporate the "GM-effect" through form drag, as in 
!  ocean_form_drag module. Default use_gm_skew=.true. 
!  </DATA> 
!
!  <DATA NAME="diffusion_all_explicit" TYPE="logical">
!  To compute all contributions from neutral diffusion explicitly in time, including
!  the K33 diagonal piece.  This approach is available only when have small time 
!  steps and/or running with just a single tracer.  It is for testing purposes. 
!  </DATA> 
!
!  <DATA NAME="neutral_physics_simple" TYPE="logical">
!  If .true. then must have aredi_equal_agm=.true..  The horizontal fluxes are then 
!  computed as horizontal downgradient diffusive fluxes regardless the neutral slope.
!  This approach precluds one from being able to have the GM-skew fluxes remain active
!  in the steep sloped regions, thus shutting off their effects to reduce the slopes 
!  of isopycnals in convective and mixed layer regimes.  It is for this reason that
!  neutral_physics_simple=.false. is the recommended default in MOM.  
!  </DATA> 
!
!  <DATA NAME="neutral_physics_limit" TYPE="logical">
!  When tracer falls outside a specified range, revert to horizontal 
!  diffusive fluxes at this cell. This is an ad hoc and incomplete attempt
!  to maintain monotonicity with the neutral physics scheme.  
!  Default neutral_physics_limit=.true.
!  </DATA> 
!  <DATA NAME="tmask_neutral_on" TYPE="logical">
!  If .true. then this logical reduces the neutral fluxes to 
!  horizontal/vertical diffusion next to boundaries.  
!  This approach has been found to reduce spurious 
!  extrema resulting from truncation of triads used to compute 
!  a neutral flux component. Default tmask_neutral_on=.false.   
!  </DATA> 
!
!  <DATA NAME="dm_taper" TYPE="logical">
!  Set to true to use the tanh tapering scheme of Danabasoglu and McWilliams.
!  Default is true. 
!  </DATA> 
!  <DATA NAME="gkw_taper" TYPE="logical">
!  Set to true to use the quadradic tapering scheme of Gerdes, Koberle, and Willebrand.
!  Default is false. 
!  </DATA> 
!
!  <DATA NAME="neutral_linear_gm_taper" TYPE="logical">
!  If .true. then with neutral_physics_simple=.false., will linearly taper GM
!  skew fluxes towards the surface within regions of steep neutral slopes.  
!  This approach leads to a constant horizontal eddy-induced velocity in 
!  the steeply sloping regions and is recommended for realistic simulations. 
!  </DATA> 
!  <DATA NAME="neutral_sine_taper" TYPE="logical">
!  If .true. then with neutral_physics_simple=.false., will apply a sine-taper 
!  to GM and neutral diffusive fluxes in regions where the penetration depth 
!  of eddies is deeper than the grid point. This method is essential for 
!  fine vertical resolution grids. 
!  </DATA> 
!
!  <DATA NAME="turb_blayer_min" TYPE="real">
!  Minimum depth of a surface turbulent boundary layer
!  used in the transition of the neutral physics fluxes
!  to the surface.  Note that in MOM4.0, 
!  turb_blayer_min was always set to zero. 
!  </DATA> 
!
!  <DATA NAME="neutral_blayer_diagnose" TYPE="logical">
!  Diagnose properties of the neutral physics boundary layer, whether have 
!  neutral_linear_gm_taper or neutral_sine_taper true or not.  
!  </DATA> 
!
!  <DATA NAME="neutral_taper_diagonal" UNITS="dimensionless" TYPE="logical">
!  For cases with neutral_physics_simple=.false., then neutral_taper_diagonal=.true.
!  will taper the diagonal pieces of the horizontal flux components when neutral slopes
!  are steep. With neutral_taper_diagonal=.false., then the horizontal flux components will 
!  remain enabled for all slopes, thus producing horizontal downgradient diffusion in 
!  regions of vertical neutral directions.
!  </DATA> 
!
!</NAMELIST>

use constants_mod,           only: epsln, pi
use diag_manager_mod,        only: register_diag_field, register_static_field, need_data
use fms_mod,                 only: FATAL, WARNING, NOTE
use fms_mod,                 only: open_namelist_file, check_nml_error, close_file, write_version_number
use mpp_domains_mod,         only: mpp_update_domains
use mpp_domains_mod,         only: CGRID_NE, WUPDATE, SUPDATE, EUPDATE, NUPDATE
use mpp_domains_mod,         only: cyclic_global_domain, global_data_domain 
use mpp_mod,                 only: input_nml_file, mpp_error, stdout, stdlog
use mpp_mod,                 only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use time_manager_mod,        only: set_time, time_type, increment_time, operator ( + )

use ocean_domains_mod,           only: get_local_indices, set_ocean_domain
use ocean_nphysics_util_mod,     only: ocean_nphysics_coeff_init, ocean_nphysics_coeff_end
use ocean_nphysics_util_mod,     only: neutral_slopes, tracer_derivs 
use ocean_nphysics_util_mod,     only: compute_eady_rate, compute_baroclinicity
use ocean_nphysics_util_mod,     only: compute_rossby_radius, compute_bczone_radius
use ocean_nphysics_util_mod,     only: compute_diffusivity, ocean_nphysics_util_restart
use ocean_nphysics_util_mod,     only: transport_on_nrho_gm, transport_on_rho_gm, transport_on_theta_gm
use ocean_nphysics_util_mod,     only: cabbeling_thermob_tendency
use ocean_nphysics_util_mod,     only: compute_eta_tend_gm90
use ocean_nphysics_util_mod,     only: watermass_diag_init, watermass_diag
use ocean_operators_mod,         only: FAX, FAY, FMX, FMY, BDX_ET, BDY_NT
use ocean_parameters_mod,        only: missing_value, onehalf, onefourth, oneeigth, DEPTH_BASED
use ocean_parameters_mod,        only: rho0r, rho0, grav
use ocean_sigma_transport_mod,   only: tmask_sigma_on, tmask_sigma  
use ocean_types_mod,             only: ocean_grid_type, ocean_domain_type, ocean_density_type
use ocean_types_mod,             only: ocean_prog_tracer_type, ocean_thickness_type
use ocean_types_mod,             only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,             only: tracer_2d_type, tracer_3d_0_nk_type, tracer_3d_1_nk_type 
use ocean_util_mod,              only: write_note, write_line, write_warning
use ocean_util_mod,              only: diagnose_2d, diagnose_3d
use ocean_workspace_mod,         only: wrk1, wrk3, wrk1_v, wrk1_2d

implicit none

public ocean_nphysicsA_init
public ocean_nphysicsA_end
public nphysicsA
public ocean_nphysicsA_restart

private fx_flux
private fy_flux
private fz_flux
private fz_terms
private fx_flux_diag
private fy_flux_diag
private fz_terms_diag 
private fz_flux_diag
private neutral_blayer
private slope_function_gm
private gm_velocity
private nphysics_diagnostics

private 

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_domain_type), save    :: Dom_flux

integer :: num_prog_tracers = 0

! clock ids
integer :: id_clock_neutral_blayer
integer :: id_clock_neutral_blayer_1
integer :: id_clock_neutral_blayer_2
integer :: id_clock_neutral_blayer_3
integer :: id_clock_fz_terms 
integer :: id_clock_fx_flux
integer :: id_clock_fy_flux
integer :: id_clock_fz_flux
integer :: id_clock_fx_flux_diag
integer :: id_clock_fy_flux_diag
integer :: id_clock_fz_flux_diag
integer :: id_clock_fz_terms_diag

! diagnostic manager ids
logical :: used
integer :: id_k33_explicit          =-1
integer :: id_ah_bdy                =-1
integer :: id_ustar                 =-1
integer :: id_vstar                 =-1
integer :: id_wstar                 =-1
integer :: id_tx_trans_gm           =-1
integer :: id_ty_trans_gm           =-1
integer :: id_depth_blayer_base     =-1
integer :: id_eddy_depth            =-1
integer :: id_steep_depth           =-1
integer :: id_slope_blayer_base     =-1
integer :: id_grav_agm_dz_sx_drhodx =-1
integer :: id_grav_agm_dz_sy_drhody =-1
integer :: id_gm_eddy_ke_source     =-1
integer :: id_slopex_drhodx         =-1
integer :: id_slopey_drhody         =-1

integer, dimension(:), allocatable  :: id_neutral_physics
integer, dimension(:), allocatable  :: id_neutral_physics_ndiffuse
integer, dimension(:), allocatable  :: id_neutral_physics_gm
integer, dimension(:), allocatable  :: id_k33_implicit
integer, dimension(:), allocatable  :: id_flux_x       ! for i-directed heat flux from neutral physics 
integer, dimension(:), allocatable  :: id_flux_y       ! for j-directed heat flux from neutral physics 
integer, dimension(:), allocatable  :: id_flux_x_int_z ! for vertically integrated i-directed tracer flux 
integer, dimension(:), allocatable  :: id_flux_y_int_z ! for vertically integrated j-directed tracer flux 

integer, dimension(:), allocatable  :: id_flux_x_gm       ! for i-directed heat flux from GM stirring 
integer, dimension(:), allocatable  :: id_flux_y_gm       ! for j-directed heat flux from GM stirring 
integer, dimension(:), allocatable  :: id_flux_x_gm_int_z ! for vertically integrated i-directed GM tracer flux 
integer, dimension(:), allocatable  :: id_flux_y_gm_int_z ! for vertically integrated j-directed GM tracer flux 

integer, dimension(:), allocatable  :: id_flux_x_ndiffuse       ! for i-directed heat flux from neutral diffusion 
integer, dimension(:), allocatable  :: id_flux_y_ndiffuse       ! for j-directed heat flux from neutral diffusion 
integer, dimension(:), allocatable  :: id_flux_x_ndiffuse_int_z ! for vertically integrated i-directed neutral diffusive flux 
integer, dimension(:), allocatable  :: id_flux_y_ndiffuse_int_z ! for vertically integrated j-directed neutral diffusive flux 

#include <ocean_memory.h>

#ifdef MOM_STATIC_ARRAYS

real, dimension(isd:ied,jsd:jed,nk,0:1) :: delqc   !density weighted (kg/m^3) quarter cell thickness(m)
real, dimension(isd:ied,jsd:jed,0:nk)   :: dzwtr   !(1/dzwt)(m^-1)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtew    !grid distances from T-point to cell faces (m)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtns    !grid distances from T-point to cell faces (m)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtwedyt !horizontal areas (m^2) of quarter cell
real, dimension(isd:ied,jsd:jed,0:1)    :: dxtdtsn !horizontal areas (m^2) of quarter cell

real, dimension(isd:ied,jsd:jed,nk) :: slopex_drhodx         !3D array of slopex * drhodx for diagnostics
real, dimension(isd:ied,jsd:jed,nk) :: slopey_drhody         !3D array of slopey * drhody for diagnostics
real, dimension(isd:ied,jsd:jed,nk) :: grav_agm_dz_sx_drhodx !3D array of grav*agm*dz*slopex*drhodx for diagnostics
real, dimension(isd:ied,jsd:jed,nk) :: grav_agm_dz_sy_drhody !3D array of grav*agm*dz*slopey*drhody for diagnostics

real, dimension(isd:ied,jsd:jed,nk) :: aredi_array       !3D array of redi diffusivities (m^2/sec)     
real, dimension(isd:ied,jsd:jed,nk) :: agm_array         !3D array of gm diffusivities (m^2/sec)        
real, dimension(isd:ied,jsd:jed)    :: ah_array          !2D array of micom horizontal diffusivities (m^2/sec)

real, dimension(isd:ied,jsd:jed)    :: bczone_radius     !for bczone calculation (m) 
real, dimension(isd:ied,jsd:jed)    :: rossby_radius     !first baroclinic Rossby radius (m) 
real, dimension(isd:ied,jsd:jed)    :: rossby_radius_raw !first baroclinic Rossby radius (m) without min/max settings  
real, dimension(isd:ied,jsd:jed,nk) :: eady_termx        !rho_z*(S_x)^2 for computing Eady growth rate 
real, dimension(isd:ied,jsd:jed,nk) :: eady_termy        !rho_z*(S_y)^2 for computing Eady growth rate 
real, dimension(isd:ied,jsd:jed)    :: baroclinic_termx  !intermediate term for computing vert ave baroclinicity
real, dimension(isd:ied,jsd:jed)    :: baroclinic_termy  !intermediate term for computing vert ave baroclinicity 
real, dimension(isd:ied,jsd:jed)    :: grid_length       !grid length scale (m)

real, dimension(isd:ied,jsd:jed,nk)     :: drhodT       !drho/dtheta     (kg/m^3/C)
real, dimension(isd:ied,jsd:jed,nk)     :: drhodS       !drho/dsalinity  (kg/m^3/psu)
real, dimension(isd:ied,jsd:jed,nk,0:1) :: drhodzb      !vertical local ref potrho derivative (kg/m^3/m)
real, dimension(isd:ied,jsd:jed,nk,0:1) :: drhodzh      !vertical local ref potrho derivative (kg/m^3/m)

real, dimension(isd:ied,jsd:jed,nk)     :: K33_implicit !density weighted (kg/m^3) implicit in time 
                                                        !diagonal term in redi tensor (m^2/sec) 
real, dimension(isd:ied,jsd:jed,nk)     :: K33_explicit !density weighted (kg/m^3) explicit in time 
                                                        !diagonal term in redi tensor (m^2/sec) 

real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: tensor_31      !tracer independent portion of mixing tensor
real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: tensor_32      !tracer independent portion of mixing tensor
real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: tensor_31_redi !tracer independent portion of Redi diffusion tensor
real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: tensor_32_redi !tracer independent portion of Redi diffusion tensor

real, dimension(isd:ied,jsd:jed,nk) :: tx_trans_gm       !for diagnosing i-transport due to GM (Sv)
real, dimension(isd:ied,jsd:jed,nk) :: ty_trans_gm       !for diagnosing j-transport due to GM (Sv)

real, dimension(isd:ied,jsd:jed)    :: depth_blayer_base ! depth (m) of boundary layer base for neutral physics
real, dimension(isd:ied,jsd:jed)    :: slope_blayer_base ! abs(slope) at base of neutral boundary layer 
real, dimension(isd:ied,jsd:jed)    :: eddy_depth        ! max of depth(m) mesoscale eddies penetrate & kpp bldepth
real, dimension(isd:ied,jsd:jed)    :: turb_blayer_depth ! depth (m) of surface turbulent boundary layer
real, dimension(isd:ied,jsd:jed)    :: N2_blayer_base    ! squared buoyancy freq (1/s^2) at base of nblayer
integer, dimension(isd:ied,jsd:jed) :: ksurf_blayer      ! k-value at base of surface nblayer 

#else

real, dimension(:,:,:,:),   allocatable :: delqc   !density weighted (kg/m^3) quarter cell thickness(m)
real, dimension(:,:,:),     allocatable :: dzwtr   !(1/dzwt)(m^-1)
real, dimension(:,:,:),     allocatable :: dtew    !grid distances from T-point to cell faces (m)
real, dimension(:,:,:),     allocatable :: dtns    !grid distances from T-point to cell faces (m)
real, dimension(:,:,:),     allocatable :: dtwedyt !horizontal areas (m^2) of quarter cell
real, dimension(:,:,:),     allocatable :: dxtdtsn !horizontal areas (m^2) of quarter cell

real, dimension(:,:,:),     allocatable :: slopex_drhodx         !3D array of slopex * drhodx for diagnostics
real, dimension(:,:,:),     allocatable :: slopey_drhody         !3D array of slopey * drhody for diagnostics
real, dimension(:,:,:),     allocatable :: grav_agm_dz_sx_drhodx !3D array of grav*agm*dz*slopex*drhodx for diagnostics
real, dimension(:,:,:),     allocatable :: grav_agm_dz_sy_drhody !3D array of grav*agm*dz*slopey*drhody for diagnostics

real, dimension(:,:,:),     allocatable :: aredi_array       !3D array of redi diffusivities (m^2/sec) 
real, dimension(:,:,:),     allocatable :: agm_array         !3D array of gm diffusivities (m^2/sec)     
real, dimension(:,:),       allocatable :: ah_array          !2D array of micom horizontal diffusivities (m^2/sec)

real, dimension(:,:),       allocatable :: bczone_radius     !for bzcone calculation (m) 
real, dimension(:,:),       allocatable :: rossby_radius     !first baroclinic Rossby radius (m) 
real, dimension(:,:),       allocatable :: rossby_radius_raw !first baroclinic Rossby radius (m) without min/max settings
real, dimension(:,:,:),     allocatable :: eady_termx        !rho_z*(S_x)^2 for computing Eady growth rate 
real, dimension(:,:,:),     allocatable :: eady_termy        !rho_z*(S_y)^2 for computing Eady growth rate 
real, dimension(:,:),       allocatable :: baroclinic_termx  !intermediate term for computing vert ave baroclinicity 
real, dimension(:,:),       allocatable :: baroclinic_termy  !intermediate term for computing vert ave baroclinicity 
real, dimension(:,:),       allocatable :: grid_length       !grid length scale (m)

real, dimension(:,:,:),     allocatable :: drhodT           !drho/dtheta     (kg/m^3/C)
real, dimension(:,:,:),     allocatable :: drhodS           !drho/dsalinity  (kg/m^3/psu)
real, dimension(:,:,:,:),   allocatable :: drhodzb          !vertical local ref potrho derivative (kg/m^3/m)
real, dimension(:,:,:,:),   allocatable :: drhodzh          !vertical local ref potrho derivative (kg/m^3/m)

real, dimension(:,:,:),     allocatable :: K33_implicit   !density weighted (kg/m^3) implicit in time 
                                                          !diagonal term in redi tensor (m^2/sec) 
real, dimension(:,:,:),     allocatable :: K33_explicit   !density weighted (kg/m^3) explicit in time 

real, dimension(:,:,:,:,:), allocatable :: tensor_31      !tracer independent portion of mixing tensor
real, dimension(:,:,:,:,:), allocatable :: tensor_32      !tracer independent portion of mixing tensor
real, dimension(:,:,:,:,:), allocatable :: tensor_31_redi !tracer independent portion of Redi diffusion tensor
real, dimension(:,:,:,:,:), allocatable :: tensor_32_redi !tracer independent portion of Redi diffusion tensor

real, dimension(:,:,:),  allocatable :: tx_trans_gm       !for diagnosing i-transport due to GM (Sv)
real, dimension(:,:,:),  allocatable :: ty_trans_gm       !for diagnosing j-transport due to GM (Sv)

real, dimension(:,:),    allocatable :: depth_blayer_base ! depth (m) of boundary layer base for neutral physics
real, dimension(:,:),    allocatable :: slope_blayer_base ! abs(slope) at base of neutral boundary layer 
real, dimension(:,:),    allocatable :: eddy_depth        ! max of depth (m) mesoscale eddies penetrate & kpp bldepth
real, dimension(:,:),    allocatable :: turb_blayer_depth ! depth (m) of surface turbulent boundary layer
real, dimension(:,:),    allocatable :: N2_blayer_base    ! squared buoyancy freq (1/s^2) at base of nblayer
integer, dimension(:,:), allocatable :: ksurf_blayer      ! k-value at base of surface nblayer 

#endif


! introduce following derived types so that do not need to know num_prog_tracers at compile time 
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: dTdx    ! tracer partial derivative (tracer/m)
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: dTdy    ! tracer partial derivative (tracer/m)
type(tracer_3d_0_nk_type), dimension(:), allocatable  :: dTdz    ! tracer partial derivative (tracer/m)
type(tracer_3d_1_nk_type), save                       :: dSdx    ! Dens%rho_salinity partial derivative (tracer/m)
type(tracer_3d_1_nk_type), save                       :: dSdy    ! Dens%rho_salinity partial derivative (tracer/m)
type(tracer_3d_0_nk_type), save                       :: dSdz    ! Dens%rho_salinity partial derivative (tracer/m)
type(tracer_2d_type),      dimension(:), allocatable  :: fz1     ! z-flux component for tracers at particular k 
type(tracer_2d_type),      dimension(:), allocatable  :: fz2     ! z-flux component for tracers at particular k 
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: flux_x  ! i-flux component for tracers for all k
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: flux_y  ! j-flux component for tracers for all k


! for diagnosing effects from GM and Redi separately, we add the following fields for neutral diffusion
! and calculation the GM contribution as a difference of the total and the neutral diffusion piece.  
type(tracer_2d_type),      dimension(:), allocatable  :: fz1_redi      ! z-flux redi component for tracers at particular k 
type(tracer_2d_type),      dimension(:), allocatable  :: fz2_redi      ! z-flux redi component for tracers at particular k 
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: flux_x_redi   ! i-flux gm component for tracers for all k
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: flux_y_redi   ! j-flux gm component for tracers for all k
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: tendency_redi ! tendency from Redi (excluding implicit K33 piece)


integer :: index_temp
integer :: index_salt
integer :: neutralrho_nk

character(len=128) :: version=&
     '$Id: ocean_nphysicsA.F90,v 20.0 2013/12/14 00:14:36 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

character(len=*), parameter :: FILENAME=&
     __FILE__

logical :: module_is_initialized = .FALSE.

! time step settings 
real    :: dtime
real    :: two_dtime_inv

! vertical coordinate 
integer :: vert_coordinate_class

! lower and upper depth for vertically averaging ocean properties.
! read into ocean_nphysics_util_nml
real :: agm_closure_upper_depth
real :: agm_closure_lower_depth


!**************nml settings**************

! for using form drag rather than GM-skewsion
! gm_skew=1.0 if use_gm_skew=.true.
! gm_skew=0.0 if use_gm_skew=.false.
logical :: use_gm_skew=.true.
real    :: gm_skew=1.0

! for setting the slope tapering methods 
real    :: smax                        ! set in ocean_neutral_util_nml
real    :: swidth                      ! set in ocean_neutral_util_nml
logical :: dm_taper         = .true.   ! tanh tapering scheme of Danabasoglu and McWilliams
real    :: swidthr                     ! inverse swidth  
real    :: smax_swidthr                ! useful combination of terms 
real    :: dm_taper_const   = 1.0      ! internally set to unity when dm_taper=.true.
logical :: gkw_taper        = .false.  ! quadratic tapering of Gerdes, Koberle, and Willebrand
real    :: gkw_taper_const  = 0.0      ! internally set to unity when gkw_taper=.true.

! neutral_physics_simple=.false. is commonly used at GFDL with MOM4 and later.
! neutral_physics_simple=.true.  was the setting more common in MOM3.  
logical :: neutral_physics_simple=.false.   

! for neutral blayer 
real :: turb_blayer_min = 0.0 ! metres (was always set to zero in MOM4.0)

! for diagnosing properties at the base of the "neutral physics boundary layer" 
logical :: neutral_blayer_diagnose=.false. 

! for maintaining horizontal GM velocity constant with depth within neutral boundary layer 
logical :: neutral_linear_gm_taper=.true.   

! for tapering neutral physics over penetration depth of 
! eddies as determined by slope and Rossby radius
logical :: neutral_sine_taper=.true.        

! to remove diagonal elements to neutral diffusion within the neutral boundar layer 
logical :: neutral_taper_diagonal=.false.   
real    :: taper_diagonal=0.0  ! taper_diagonal=1 if neutral_taper_diagonal=.true. ; =0 otherwise 

! to reduce neutral fluxes to horz/vert diffusion next to model boundaries
logical :: tmask_neutral_on=.false.         

! to compute K33 explicitly in time. 
! diffusion_all_explicit=.false. for realistic simulations.
logical :: diffusion_all_explicit=.false.   

! revert to horizontal diffusion when tracer falls outside specified range 
logical :: neutral_physics_limit=.true.   

! for specifying transport units
! can either be Sv or mks
character(len=32) :: transport_dims ='Sv (10^9 kg/s)' 
real              :: transport_convert=1.0e-9 

! internally set; for diagnosing the effects from GM and Redi separately 
logical :: diagnose_gm_redi = .false.

! for the module as a whole 
logical :: use_this_module   = .false.
logical :: debug_this_module = .false.

!**************end of nml settings**************

namelist /ocean_nphysicsA_nml/ use_this_module, debug_this_module,          &
          neutral_physics_limit, use_gm_skew,                               &
          dm_taper, gkw_taper, tmask_neutral_on,                            &
          neutral_physics_simple, neutral_taper_diagonal,                   &
          neutral_linear_gm_taper, neutral_sine_taper,                      & 
          diffusion_all_explicit, neutral_blayer_diagnose, turb_blayer_min
    

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_nphysicsA_init">
!
! <DESCRIPTION>
! Initialize the neutral physics module by registering fields for 
! diagnostic output and performing some numerical checks to see 
! that namelist settings are appropriate. 
! </DESCRIPTION>
!
subroutine ocean_nphysicsA_init(Grid, Domain, Time, Time_steps, Thickness, Dens, T_prog, &
           ver_coordinate_class, agm_closure_lower_dept, agm_closure_upper_dept,         &
           smx, swidt, cmip_units, debug)

  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_time_steps_type),  intent(in)           :: Time_steps
  type(ocean_thickness_type),   intent(in)           :: Thickness
  type(ocean_density_type),     intent(in)           :: Dens
  type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
  integer,                      intent(in)           :: ver_coordinate_class
  real,                         intent(in)           :: agm_closure_lower_dept
  real,                         intent(in)           :: agm_closure_upper_dept
  real,                         intent(in)           :: smx 
  real,                         intent(in)           :: swidt
  logical,                      intent(in)           :: cmip_units
  logical,                      intent(in), optional :: debug

  logical :: diagnose_gm_redi_input=.false.
  integer :: ioun, io_status, ierr  
  integer :: i, j, n
  
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_nphysicsA_mod (ocean_nphysicsA_init):already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  num_prog_tracers = size(T_prog(:))
  dtime            = Time_steps%dtime_t
  Dom => Domain
  Grd => Grid

  vert_coordinate_class   = ver_coordinate_class
  agm_closure_lower_depth = agm_closure_lower_dept
  agm_closure_upper_depth = agm_closure_upper_dept
  smax                    = smx 
  swidth                  = swidt
  neutralrho_nk           = size(Dens%neutralrho_ref(:))

  ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_nphysicsA_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_nphysicsA_nml')
#else
  ioun =  open_namelist_file()
  read (ioun,ocean_nphysicsA_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_nphysicsA_nml')
  call close_file (ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_nphysicsA_nml)  
  write (stdlogunit,ocean_nphysicsA_nml)

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk
#endif

  if(use_this_module) then 
     call write_note(FILENAME,&
      'USING ocean_nphysicsA_mod.')
  else 
     call write_note(FILENAME,&
     'NOT using ocean_nphysicsA_mod.')
     return
  endif 

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
     call write_note(FILENAME,&
     'running ocean_nphysicsA_mod with debug_this_module=.true.')
  endif 

  write(stdoutunit,'(/1x,a,f10.2)') &
  '==> Note from ocean_nphysicsA_mod: using forward time step of (secs)', dtime 

  if(cmip_units) then
      transport_convert=1.0
      transport_dims   = 'kg/s'
  else
      transport_convert=1.0e-9 
      transport_dims   = 'Sv (10^9 kg/s)'
  endif

  ! Useful constants 
  two_dtime_inv = 0.5/dtime            !for explicit piece of K33 
  swidthr       = 1.0/(swidth + epsln) !for slope taper function when dm_taper used 
  smax_swidthr  = smax*swidthr         !for slope taper function when dm_taper used 

  if(use_gm_skew) then
    gm_skew=1.0
    call write_note(FILENAME,&
    'use_gm_skew=.true. so will use GM-skewsion.')
  else
    gm_skew=0.0
    call write_note(FILENAME,&
    'use_gm_skew=.false. so will NOT use GM-skewsion.' )
    if(neutral_physics_simple) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_nphysicsA_mod: use_gm_skew=.false. incompatible with neutral_physics_simple=.true.')
    endif 
  endif 

  if(neutral_physics_simple) then
     call write_note(FILENAME,&
     'neutral_physics_simple=.true.')
     call write_line('means the horizontal SGS tracer fluxes are downgradient for all neutral slopes.')
     call write_line('Assumes aredi_array = agm_array.')
     call write_line('Analogous to MOM3 approach.')
  endif

  if(neutral_physics_limit) then
     call write_note(FILENAME,&
     'neutral_physics_limit=.true.')
     call write_line('Will revert to horizontal diffusion for points where tracer is outside specified range.')
  endif

  if(dm_taper .and. .not. gkw_taper) then
     call write_note(FILENAME,&
     'dm_taper=.true. Will use the tanh scheme')
     call write_line('of Danabasoglu and McWilliams to taper neutral physics in steep sloped regions')
     dm_taper_const =1.0
     gkw_taper_const=0.0
  endif
  if(gkw_taper .and. .not. dm_taper) then
     call write_note(FILENAME,&
     'gkw_taper=.true. Will use the quadratic scheme')
     call write_line('of Gerdes, Koberle, and Willebrand to taper neutral physics in steep sloped regions')
     dm_taper_const =0.0
     gkw_taper_const=1.0
  endif
  if(gkw_taper .and. dm_taper) then
      dm_taper_const =0.0
      gkw_taper_const=0.0
      call mpp_error(FATAL, &
      '==>Error from ocean_nphysicsA_mod: gkw_taper and dm_taper cannot both be set true--choose only one.')
  endif

  if(.not. neutral_physics_simple) then 
      if(neutral_linear_gm_taper) then
         call write_note(FILENAME,&
         'neutral_linear_gm_taper=.true., so will linearly')
         call write_line('taper GM towards the surface when reaching steep neutral slopes in surface bdy.')
      else
         call write_note(FILENAME,&
         'GM exponentially tapered to zero in neutral bdy layer region.')
      endif
      if(neutral_sine_taper) then
         call write_note(FILENAME,&
         'Running with neutral_sine_taper=.true., and so will use sine-taper')
         call write_line('on fluxes where eddy penetration depth and/or KPP hblt exceeds grid depth.')
      endif
  endif

  if(neutral_blayer_diagnose) then 
     call write_note(FILENAME,&
     'Running with neutral_blayer_diagnose=.true., so will diagnose properties')
     call writE_line('of the neutral physics boundary layer.')
  endif

  if(neutral_linear_gm_taper) then
     call write_note(FILENAME,&
     'Running with a nontrivial GM transport in steep neutral slope regions.')
  endif

  if(diffusion_all_explicit) then
     call write_warning(FILENAME,&
     'Running w/ diffusion_all_explicit=.true., which means compute K33 contribution')
     call write_line('to neutral diffusion explicitly in time.  This method is stable ONLY if taking')
     call write_line('very small time steps and/or running with just a single tracer.')
  endif

  if(neutral_taper_diagonal) then
      taper_diagonal = 1.0 
      call write_note(FILENAME,&
      'Running w/ neutral_taper_diagonal=.true. and so taper_diagonal = 1.0')
      call write_line('Will taper all pieces of horiz neutral flux components in steep neutral slope regions.')
      call write_line('neutral_taper_diagonal=.false. will alternatively keep the diagonal pieces untapered.')
  endif
  if(.not. neutral_taper_diagonal) then
      taper_diagonal = 0.0 
      call write_note(FILENAME,&
      'Running w/ neutral_taper_diagonal=.false. so taper_diagonal = 0.0')
      call write_line('Diagonal pieces of horizontal flux components are untapered regardless the neutral slope.')
  endif


  do n=1,num_prog_tracers
    T_prog(n)%neutral_physics_limit = neutral_physics_limit
  enddo 

  allocate( dTdx(num_prog_tracers) )
  allocate( dTdy(num_prog_tracers) )
  allocate( dTdz(num_prog_tracers) )
  allocate( fz1(num_prog_tracers) )
  allocate( fz2(num_prog_tracers) )
  allocate( flux_x(num_prog_tracers) )
  allocate( flux_y(num_prog_tracers) )


  call set_ocean_domain(Dom_flux,Grid,xhalo=Dom%xhalo,yhalo=Dom%yhalo,name='flux dom neutral',maskmap=Dom%maskmap)

#ifndef MOM_STATIC_ARRAYS
  allocate (dtew(isd:ied,jsd:jed,0:1))
  allocate (dtns(isd:ied,jsd:jed,0:1))
  allocate (dtwedyt(isd:ied,jsd:jed,0:1))
  allocate (dxtdtsn(isd:ied,jsd:jed,0:1))
  allocate (grid_length(isd:ied,jsd:jed))
  allocate (delqc(isd:ied,jsd:jed,nk,0:1))
  allocate (dzwtr(isd:ied,jsd:jed,0:nk))
  allocate (aredi_array(isd:ied,jsd:jed,nk)) 
  allocate (agm_array(isd:ied,jsd:jed,nk))   
  allocate (ah_array(isd:ied,jsd:jed))
  allocate (bczone_radius(isd:ied,jsd:jed))
  allocate (rossby_radius(isd:ied,jsd:jed))
  allocate (rossby_radius_raw(isd:ied,jsd:jed))
  allocate (eady_termx(isd:ied,jsd:jed,nk))
  allocate (eady_termy(isd:ied,jsd:jed,nk))
  allocate (baroclinic_termx(isd:ied,jsd:jed))
  allocate (baroclinic_termy(isd:ied,jsd:jed))
  allocate (tx_trans_gm(isd:ied,jsd:jed,nk))
  allocate (ty_trans_gm(isd:ied,jsd:jed,nk))
  allocate (grav_agm_dz_sx_drhodx(isd:ied,jsd:jed,nk)) 
  allocate (grav_agm_dz_sy_drhody(isd:ied,jsd:jed,nk)) 
  allocate (slopex_drhodx(isd:ied,jsd:jed,nk)) 
  allocate (slopey_drhody(isd:ied,jsd:jed,nk)) 

  allocate (drhodT(isd:ied,jsd:jed,nk))
  allocate (drhodS(isd:ied,jsd:jed,nk))
  allocate (drhodzb(isd:ied,jsd:jed,nk,0:1))
  allocate (drhodzh(isd:ied,jsd:jed,nk,0:1))
  allocate (tensor_31(isd:ied,jsd:jed,nk,0:1,0:1))
  allocate (tensor_32(isd:ied,jsd:jed,nk,0:1,0:1))
  allocate (K33_implicit(isd:ied,jsd:jed,nk)) 
  allocate (K33_explicit(isd:ied,jsd:jed,nk)) 

  do n=1,num_prog_tracers
    allocate ( dTdx(n)%field(isd:ied,jsd:jed,nk) )
    allocate ( dTdy(n)%field(isd:ied,jsd:jed,nk) )
    allocate ( dTdz(n)%field(isd:ied,jsd:jed,0:nk) )
    allocate ( fz1(n)%field(isd:ied,jsd:jed) )
    allocate ( fz2(n)%field(isd:ied,jsd:jed) )
    allocate ( flux_x(n)%field(isd:ied,jsd:jed,nk) )
    allocate ( flux_y(n)%field(isd:ied,jsd:jed,nk) )
  enddo 
  allocate ( dSdx%field(isd:ied,jsd:jed,nk) )
  allocate ( dSdy%field(isd:ied,jsd:jed,nk) )
  allocate ( dSdz%field(isd:ied,jsd:jed,0:nk) )

  allocate(depth_blayer_base(isd:ied,jsd:jed))
  allocate(slope_blayer_base(isd:ied,jsd:jed))
  allocate(eddy_depth(isd:ied,jsd:jed))
  allocate(turb_blayer_depth(isd:ied,jsd:jed))
  allocate(N2_blayer_base(isd:ied,jsd:jed))
  allocate(ksurf_blayer(isd:ied,jsd:jed))

#endif

  do n=1,num_prog_tracers 
    dTdx(n)%field(:,:,:)   = 0.0
    dTdy(n)%field(:,:,:)   = 0.0
    dTdz(n)%field(:,:,:)   = 0.0
    fz1(n)%field(:,:)      = 0.0
    fz2(n)%field(:,:)      = 0.0
    flux_x(n)%field(:,:,:) = 0.0
    flux_y(n)%field(:,:,:) = 0.0

  enddo  
  dSdx%field(:,:,:) = 0.0
  dSdy%field(:,:,:) = 0.0
  dSdz%field(:,:,:) = 0.0

  depth_blayer_base  = 0.0
  slope_blayer_base  = 0.0
  eddy_depth         = 0.0
  turb_blayer_depth  = 0.0
  N2_blayer_base     = 0.0
  ksurf_blayer       = 1
  tx_trans_gm        = 0.0
  ty_trans_gm        = 0.0

  dtew(:,:,0) = Grd%dtw(:,:)
  dtew(:,:,1) = Grd%dte(:,:)
  dtns(:,:,0) = Grd%dts(:,:)
  dtns(:,:,1) = Grd%dtn(:,:)

  dtwedyt(:,:,:) = 0.0
  dtwedyt(:,:,0) = Grd%dte(:,:)*Grd%dyt(:,:)
  do i=isc-1,iec
    dtwedyt(i,:,1) = Grd%dtw(i+1,:)*Grd%dyt(i+1,:)
  enddo

  dxtdtsn(:,:,:) = 0.0
  dxtdtsn(:,:,0) = Grd%dxt(:,:)*Grd%dtn(:,:)
  do j=jsc-1,jec
    dxtdtsn(:,j,1) = Grd%dxt(:,j+1)*Grd%dts(:,j+1)
  enddo

  grid_length(:,:) = 2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(Grd%dxt(:,:)+Grd%dyt(:,:))

  eady_termx(:,:,:)            = 0.0
  eady_termy(:,:,:)            = 0.0
  baroclinic_termx(:,:)        = 0.0
  baroclinic_termy(:,:)        = 0.0
  grav_agm_dz_sx_drhodx(:,:,:) = 0.0
  grav_agm_dz_sy_drhody(:,:,:) = 0.0
  slopex_drhodx(:,:,:)         = 0.0
  slopey_drhody(:,:,:)         = 0.0

  rossby_radius(:,:)           = 0.0
  rossby_radius_raw(:,:)       = 0.0
  bczone_radius(:,:)           = 0.0
  agm_array(:,:,:)             = 0.0
  aredi_array(:,:,:)           = 0.0
  ah_array(:,:)                = 0.0
  call ocean_nphysics_coeff_init(Time, Thickness, rossby_radius, rossby_radius_raw, &
                         bczone_radius, agm_array, aredi_array, ah_array)

  drhodT(:,:,:)             = 0.0
  drhodS(:,:,:)             = 0.0
  drhodzb(:,:,:,:)          = 0.0
  drhodzh(:,:,:,:)          = 0.0
  tensor_31(:,:,:,:,:)      = 0.0
  tensor_32(:,:,:,:,:)      = 0.0
  K33_implicit(:,:,:)       = 0.0           
  K33_explicit(:,:,:)       = 0.0           

  index_temp=-1;index_salt=-1
  do n=1,num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo
  if (index_temp == -1 .or. index_salt == -1) then 
     call mpp_error(FATAL, &
     '==>Error: temp and/or salt not identified in call to ocean_nphysicsA_init')
  endif 

  ! initialize clock ids 
  id_clock_neutral_blayer   = mpp_clock_id('(Ocean neutral: blayer)'       ,grain=CLOCK_ROUTINE)
  id_clock_neutral_blayer_1 = mpp_clock_id('(Ocean neutral: blayer1)'      ,grain=CLOCK_ROUTINE)
  id_clock_neutral_blayer_2 = mpp_clock_id('(Ocean neutral: blayer2)'      ,grain=CLOCK_ROUTINE)
  id_clock_neutral_blayer_3 = mpp_clock_id('(Ocean neutral: blayer3)'      ,grain=CLOCK_ROUTINE)
  id_clock_fz_terms         = mpp_clock_id('(Ocean neutral: fz-terms)'     ,grain=CLOCK_ROUTINE)
  id_clock_fx_flux          = mpp_clock_id('(Ocean neutral: fx-flux)'      ,grain=CLOCK_ROUTINE)
  id_clock_fy_flux          = mpp_clock_id('(Ocean neutral: fy-flux)'      ,grain=CLOCK_ROUTINE)
  id_clock_fz_flux          = mpp_clock_id('(Ocean neutral: fz-flux)'      ,grain=CLOCK_ROUTINE)
  id_clock_fx_flux_diag     = mpp_clock_id('(Ocean neutral: fx-flux-diag)' ,grain=CLOCK_ROUTINE)
  id_clock_fy_flux_diag     = mpp_clock_id('(Ocean neutral: fy-flux-diag)' ,grain=CLOCK_ROUTINE)
  id_clock_fz_flux_diag     = mpp_clock_id('(Ocean neutral: fz-flux-diag)' ,grain=CLOCK_ROUTINE)
  id_clock_fz_terms_diag    = mpp_clock_id('(Ocean neutral: fz-terms-diag)',grain=CLOCK_ROUTINE)


  ! register fields for diagnostic output 

  id_ah_bdy = -1
  id_ah_bdy = register_static_field ('ocean_model', 'ah_bdy',                  &
              Grd%tracer_axes(1:2), 'Horz diffusivity in surface neutral bdy', &
              'm^2/sec', missing_value=missing_value, range=(/-10.0,1.e10/))
  call diagnose_2d(Time, Grd, id_ah_bdy, ah_array(:,:))

  id_k33_explicit = -1
  id_k33_explicit = register_diag_field ('ocean_model', 'k33_explicit', &
                    Grd%tracer_axes_wt(1:3), Time%model_time,           &
                    'K33_explicit tensor element', 'm^2/sec',           &
                    missing_value=missing_value, range=(/-10.0,1.e20/))

  id_tx_trans_gm = -1
  id_tx_trans_gm = register_diag_field ('ocean_model', 'tx_trans_gm',      &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time,           &
                   'T-cell mass i-transport from GM',trim(transport_dims), &
                   missing_value=missing_value, range=(/-1e20,1.e20/))

  id_ty_trans_gm = -1
  id_ty_trans_gm = register_diag_field ('ocean_model', 'ty_trans_gm',     &
                   Grd%tracer_axes_flux_y(1:3), Time%model_time,          &
                   'T-cell mass j-transport from GM',trim(transport_dims),&
                   missing_value=missing_value, range=(/-1e20,1.e20/))

  id_depth_blayer_base = -1
  id_depth_blayer_base = register_diag_field ('ocean_model','depth_blayer_base', &
                         Grd%tracer_axes(1:2), Time%model_time,                  &
                         'Depth of neutral physics surface bdy-layer', 'm',      &
                         missing_value=missing_value, range=(/-1.0,1.e8/))

  id_slope_blayer_base = -1
  id_slope_blayer_base = register_diag_field ('ocean_model','slope_blayer_base',          &
                         Grd%tracer_axes(1:2), Time%model_time,                           &
                         'abs(slope) at neutral physics bdy-layer base', 'dimensionless', &
                         missing_value=missing_value, range=(/-1.0,1.e8/))

  id_eddy_depth = -1
  id_eddy_depth = register_diag_field ('ocean_model','eddy_depth',                 &
                  Grd%tracer_axes(1:2), Time%model_time,                           &
                  'mesoscale eddy penetration depth used in neutral physicsA', 'm',&
                  missing_value=missing_value, range=(/-1.0,1.e8/))

  id_steep_depth = -1
  id_steep_depth = register_diag_field ('ocean_model','steep_depth',     &
                   Grd%tracer_axes(1:2), Time%model_time,                &
                  'Depth at base of steep slope surface bdy layer', 'm', &
                   missing_value=missing_value, range=(/-1.0,1.e8/))

  id_grav_agm_dz_sx_drhodx = -1
  id_grav_agm_dz_sx_drhodx = register_diag_field ('ocean_model', 'grav_agm_dz_sx_drhodx', &
               Grd%tracer_axes(1:3), Time%model_time,                                     &
               'grav*dz*agm*neutral x-slope times drho_dx', 'W/m^2',                      &
               missing_value=missing_value, range=(/-1.e15,1.e15/))

  id_grav_agm_dz_sy_drhody = -1
  id_grav_agm_dz_sy_drhody = register_diag_field ('ocean_model', 'grav_agm_dz_sy_drhody', &
               Grd%tracer_axes(1:3), Time%model_time,                                     &
               'grav*dz*agm*neutral y-slope times drho_dy', 'W/m^2',                      &
               missing_value=missing_value, range=(/-1.e15,1.e15/))

  id_gm_eddy_ke_source = -1
  id_gm_eddy_ke_source = register_diag_field ('ocean_model', 'gm_eddy_ke_source',&
               Grd%tracer_axes(1:3), Time%model_time,                            &
               'rho0*dz*agm*(slope*N)^2 = eddy ke source from GM', 'W/m^2',      &
               missing_value=missing_value, range=(/-1.e1,1.e15/),               &
               standard_name='tendency_of_ocean_eddy_kinetic_energy_content_due_to_bolus_transport')

  id_slopex_drhodx = -1
  id_slopex_drhodx = register_diag_field ('ocean_model', 'slopex_drhodx', &
               Grd%tracer_axes(1:3), Time%model_time,                     &
               'neutral x-slope times drho_dx', 'kg/m^4',                 &
               missing_value=missing_value, range=(/-1.e10,1.e10/))

  id_slopey_drhody = -1
  id_slopey_drhody = register_diag_field ('ocean_model', 'slopey_drhody', &
               Grd%tracer_axes(1:3), Time%model_time,                     &
               'neutral y-slope times drho_dy', 'kg/m^4',                 &
               missing_value=missing_value, range=(/-1.e10,1.e10/))

  id_ustar = -1
  id_ustar = register_diag_field ('ocean_model', 'ugm', &
             Grd%vel_axes_uv(1:3), Time%model_time,     &
             'GM zonal velocity', 'm/sec',              &
             missing_value=missing_value, range=(/-10.0,10.0/))

  id_vstar = -1    
  id_vstar = register_diag_field ('ocean_model', 'vgm', &
             Grd%vel_axes_uv(1:3), Time%model_time,     &
             'GM merid velocity', 'm/sec',              &
             missing_value=missing_value, range=(/-10.0,10.0/))

  id_wstar = -1
  id_wstar = register_diag_field ('ocean_model', 'wgm',   &
             Grd%vel_axes_wt(1:3), Time%model_time,       &
             'GM vert velocity (T-cell bottom)', 'm/sec', &
             missing_value=missing_value, range=(/-10.0,10.0/))

  call watermass_diag_init(Time, Dens, diagnose_gm_redi_input)
  if(diagnose_gm_redi_input) then 
     diagnose_gm_redi = .true.
  endif 


  allocate (id_k33_implicit(num_prog_tracers))
  allocate (id_neutral_physics(num_prog_tracers))
  allocate (id_neutral_physics_ndiffuse(num_prog_tracers))
  allocate (id_neutral_physics_gm(num_prog_tracers))
  id_k33_implicit             = -1
  id_neutral_physics          = -1
  id_neutral_physics_ndiffuse = -1
  id_neutral_physics_gm       = -1
  do n=1,num_prog_tracers
     id_k33_implicit(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_k33_implicit', &
                          Grd%tracer_axes_wt(1:3), Time%model_time,                                  &
                          'K33_implicit tensor element for '//trim(T_prog(n)%name),                  &
                          'm^2/sec', missing_value=missing_value, range=(/-10.0,1.e10/))

     if (T_prog(n)%name == 'temp') then
       id_neutral_physics(n) = register_diag_field ('ocean_model', 'neutral_'//trim(T_prog(n)%name), &
                               Grd%tracer_axes(1:3), Time%model_time,                                &
                               'rho*dzt*cp*explicit neutral tendency (heating)',                     &
                               trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e10,1.e10/))

       id_neutral_physics_ndiffuse(n) = register_diag_field ('ocean_model', 'neutral_diffusion_'//trim(T_prog(n)%name), &
                               Grd%tracer_axes(1:3), Time%model_time,                                                   &
                               'rho*dzt*cp*explicit neutral diffusion tendency (heating)',                              &
                               trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e10,1.e10/))
       if(id_neutral_physics_ndiffuse(n) > 0) diagnose_gm_redi=.true.

       id_neutral_physics_gm(n) = register_diag_field ('ocean_model', 'neutral_gm_'//trim(T_prog(n)%name), &
                               Grd%tracer_axes(1:3), Time%model_time,                                      &
                               'rho*dzt*cp*GM stirring (heating)',                                         &
                               trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e10,1.e10/))
       if(id_neutral_physics_gm(n) > 0) diagnose_gm_redi=.true.

     else 

       id_neutral_physics(n) = register_diag_field ('ocean_model', 'neutral_'//trim(T_prog(n)%name), &
                               Grd%tracer_axes(1:3), Time%model_time,                                &
                               'rho*dzt*explicit neutral tendency for '//trim(T_prog(n)%name),       &
                               trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e10,1.e10/))

       id_neutral_physics_ndiffuse(n) = register_diag_field ('ocean_model', 'neutral_diffusion_'//trim(T_prog(n)%name), &
                               Grd%tracer_axes(1:3), Time%model_time,                                                   &
                               'rho*dzt*explicit neutral diffusion tendency for '//trim(T_prog(n)%name),                &
                               trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e10,1.e10/))
       if(id_neutral_physics_ndiffuse(n) > 0) diagnose_gm_redi=.true.

       id_neutral_physics_gm(n) = register_diag_field ('ocean_model', 'neutral_gm_'//trim(T_prog(n)%name),&
                               Grd%tracer_axes(1:3), Time%model_time,                                     &
                               'rho*dzt*GM stirring tendency for '//trim(T_prog(n)%name),                 &
                               trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e10,1.e10/))
       if(id_neutral_physics_gm(n) > 0) diagnose_gm_redi=.true.

     endif 
  enddo 

  allocate (id_flux_x(num_prog_tracers))
  allocate (id_flux_y(num_prog_tracers))
  allocate (id_flux_x_int_z(num_prog_tracers))
  allocate (id_flux_y_int_z(num_prog_tracers))

  allocate (id_flux_x_gm(num_prog_tracers))
  allocate (id_flux_y_gm(num_prog_tracers))
  allocate (id_flux_x_gm_int_z(num_prog_tracers))
  allocate (id_flux_y_gm_int_z(num_prog_tracers))

  allocate (id_flux_x_ndiffuse(num_prog_tracers))
  allocate (id_flux_y_ndiffuse(num_prog_tracers))
  allocate (id_flux_x_ndiffuse_int_z(num_prog_tracers))
  allocate (id_flux_y_ndiffuse_int_z(num_prog_tracers))

  id_flux_x       = -1
  id_flux_y       = -1
  id_flux_x_int_z = -1
  id_flux_y_int_z = -1

  id_flux_x_gm       = -1
  id_flux_y_gm       = -1
  id_flux_x_gm_int_z = -1
  id_flux_y_gm_int_z = -1

  id_flux_x_ndiffuse       = -1
  id_flux_y_ndiffuse       = -1
  id_flux_x_ndiffuse_int_z = -1
  id_flux_y_ndiffuse_int_z = -1

  do n=1,num_prog_tracers
     if(n == index_temp) then 

         id_flux_x(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_neutral', &
              Grd%tracer_axes_flux_x(1:3), Time%model_time, 'cp*neutral_xflux*dyt*rho_dzt*temp',    &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         id_flux_y(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_neutral', &
              Grd%tracer_axes_flux_y(1:3), Time%model_time, 'cp*neutral_yflux*dxt*rho_dzt*temp',    &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         id_flux_x_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_neutral_int_z', &
              Grd%tracer_axes_flux_x(1:2), Time%model_time, 'z-integral cp*neutral_xflux*dyt*rho_dzt*temp',     &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         id_flux_y_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_neutral_int_z', &
              Grd%tracer_axes_flux_y(1:2), Time%model_time, 'z-integral cp*neutral_yflux*dxt*rho_dzt*temp',     &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))

         id_flux_x_gm(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_gm',&
              Grd%tracer_axes_flux_x(1:3), Time%model_time, 'cp*gm_xflux*dyt*rho_dzt*temp',      &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_x_gm(n) > 0) diagnose_gm_redi=.true.

         id_flux_y_gm(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_gm',&
              Grd%tracer_axes_flux_y(1:3), Time%model_time, 'cp*gm_yflux*dxt*rho_dzt*temp',      &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_y_gm(n) > 0) diagnose_gm_redi=.true.

         id_flux_x_gm_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_gm_int_z',&
              Grd%tracer_axes_flux_x(1:2), Time%model_time, 'z-integral cp*gm_xflux*dyt*rho_dzt*temp',       &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_x_gm_int_z(n) > 0) diagnose_gm_redi=.true.

         id_flux_y_gm_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_gm_int_z',&
              Grd%tracer_axes_flux_y(1:2), Time%model_time, 'z-integral cp*gm_yflux*dxt*rho_dzt*temp',       &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_y_gm_int_z(n) > 0) diagnose_gm_redi=.true.

         id_flux_x_ndiffuse(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_ndiffuse',&
              Grd%tracer_axes_flux_x(1:3), Time%model_time, 'cp*redi_xflux*dyt*rho_dzt*temp',                &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_x_ndiffuse(n) > 0) diagnose_gm_redi=.true.

         id_flux_y_ndiffuse(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_ndiffuse',&
              Grd%tracer_axes_flux_y(1:3), Time%model_time, 'cp*redi_yflux*dxt*rho_dzt*temp',                &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_y_ndiffuse(n) > 0) diagnose_gm_redi=.true.

         id_flux_x_ndiffuse_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_ndiffuse_int_z',&
              Grd%tracer_axes_flux_x(1:2), Time%model_time, 'z-integral cp*redi_xflux*dyt*rho_dzt*temp',                 &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_x_ndiffuse_int_z(n) > 0) diagnose_gm_redi=.true.

         id_flux_y_ndiffuse_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_ndiffuse_int_z',&
              Grd%tracer_axes_flux_y(1:2), Time%model_time, 'z-integral cp*redi_yflux*dxt*rho_dzt*temp',                 &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_y_ndiffuse_int_z(n) > 0) diagnose_gm_redi=.true.

     else

         id_flux_x(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_neutral', &
              Grd%tracer_axes_flux_x(1:3), Time%model_time,                                         &
              'neutral_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name),                         &
              'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
         id_flux_y(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_neutral', &
              Grd%tracer_axes_flux_y(1:3), Time%model_time,                                         &
              'neutral_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name),                         &
              'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
         id_flux_x_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_neutral_int_z', &
              Grd%tracer_axes_flux_x(1:2), Time%model_time,                                                     &
              'z-integral neutral_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name),                          &
              'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
         id_flux_y_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_neutral_int_z', &
              Grd%tracer_axes_flux_y(1:2), Time%model_time,                                                     &
              'z-integral neutral_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name),                          &
              'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

         id_flux_x_gm(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_gm',&
              Grd%tracer_axes_flux_x(1:3), Time%model_time,                                      &
              'gm_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name),                           &
              'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_x_gm(n) > 0) diagnose_gm_redi=.true.

         id_flux_y_gm(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_gm', &
              Grd%tracer_axes_flux_y(1:3), Time%model_time,                                       &
              'gm_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name),                            &
              'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_y_gm(n) > 0) diagnose_gm_redi=.true.

         id_flux_x_gm_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_gm_int_z', &
              Grd%tracer_axes_flux_x(1:2), Time%model_time,                                                   &
              'z-integral gm_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name),                             &
              'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_x_gm_int_z(n) > 0) diagnose_gm_redi=.true.

         id_flux_y_gm_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_gm_int_z', &
              Grd%tracer_axes_flux_y(1:2), Time%model_time,                                                   &
              'z-integral gm_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name),                             &
              'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_y_gm_int_z(n) > 0) diagnose_gm_redi=.true.

         id_flux_x_ndiffuse(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_ndiffuse',&
              Grd%tracer_axes_flux_x(1:3), Time%model_time,                                          &
              'gm_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name),                               &
              'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_x_ndiffuse(n) > 0) diagnose_gm_redi=.true.

         id_flux_y_ndiffuse(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_ndiffuse',&
              Grd%tracer_axes_flux_y(1:3), Time%model_time,                                          &
              'gm_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name),                               &
              'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_y_ndiffuse(n) > 0) diagnose_gm_redi=.true.

         id_flux_x_ndiffuse_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_ndiffuse_int_z',&
              Grd%tracer_axes_flux_x(1:2), Time%model_time,                                                              &
              'z-integral gm_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name),                                        &
              'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_x_ndiffuse_int_z(n) > 0) diagnose_gm_redi=.true.

         id_flux_y_ndiffuse_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_ndiffuse_int_z',&
              Grd%tracer_axes_flux_y(1:2), Time%model_time,                                                              &
              'z-integral gm_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name),                                        &
              'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
         if(id_flux_y_ndiffuse_int_z(n) > 0) diagnose_gm_redi=.true.

     endif

  enddo

  if(diagnose_gm_redi) then 
     if(neutral_physics_simple) then 
        call write_note(FILENAME,&
        'diagnostics from diagnose_gm_redi=.true. are not available w/ neutral_physics_simple=.true.')
     endif 
     call write_note(FILENAME,&
     'ocean_nphysicA has diagnose_gm_redi=.true. to diagnose GM and Redi. It adds memory requirements.')
     allocate( fz1_redi(num_prog_tracers) )
     allocate( fz2_redi(num_prog_tracers) )
     allocate( flux_x_redi(num_prog_tracers) )
     allocate( flux_y_redi(num_prog_tracers) )
     allocate( tendency_redi(num_prog_tracers) )

#ifndef MOM_STATIC_ARRAYS
      allocate (tensor_31_redi(isd:ied,jsd:jed,nk,0:1,0:1))
      allocate (tensor_32_redi(isd:ied,jsd:jed,nk,0:1,0:1))
      do n=1,num_prog_tracers
         allocate ( fz1_redi(n)%field(isd:ied,jsd:jed) )
         allocate ( fz2_redi(n)%field(isd:ied,jsd:jed) )
         allocate ( flux_x_redi(n)%field(isd:ied,jsd:jed,nk) )
         allocate ( flux_y_redi(n)%field(isd:ied,jsd:jed,nk) )
         allocate ( tendency_redi(n)%field(isd:ied,jsd:jed,nk) )
      enddo
#endif 
      tensor_31_redi(:,:,:,:,:) = 0.0
      tensor_32_redi(:,:,:,:,:) = 0.0
      do n=1,num_prog_tracers
         fz1_redi(n)%field(:,:)        = 0.0
         fz2_redi(n)%field(:,:)        = 0.0
         flux_x_redi(n)%field(:,:,:)   = 0.0
         flux_y_redi(n)%field(:,:,:)   = 0.0
         tendency_redi(n)%field(:,:,:) = 0.0
      enddo

  endif


end subroutine ocean_nphysicsA_init
! </SUBROUTINE>  NAME="ocean_nphysicsA_init"


!#######################################################################
! <SUBROUTINE NAME="nphysicsA">
!
! <DESCRIPTION>
! This function computes the thickness weighted and density weighted
! time tendency for tracer from neutral physics.  Full discussion
! and details are provided by Griffies (2004). 
!
! Here is a brief summary.  
!
!---How the neutral diffusive flux components are computed:
!
! The vertical flux component is split into diagonal (3,3) and 
! off-diagonal (3,1) and (3,2) terms. The off-diagonal (3,1) and (3,2) 
! terms are included explicitly in time. The main contribution from the 
! (3,3) term to the time tendency is included implicitly in time 
! along with the usual contribution from diapycnal processes 
! (vertical mixing schemes).  This is the K33_implicit term.
! This approach is necessary with high vertical resolution, as 
! noted by Cox (1987).  However, splitting the vertical flux into 
! an implicit and explicit piece compromises the 
! integrity of the vertical flux component (see Griffies et al. 1998).
! So to minimize the disparity engendered by this split, the portion of 
! K33 that can be stably included explicitly in time is computed along 
! with the (3,1) and (3,2) terms. 
! 
! All other terms in the mixing tensor are included explicitly in time
! using a forward time step as required for temporal stability of 
! numerical diffusive processes.  
!
! The off-diagonal terms in the horizontal flux components, and all terms
! in the vertical flux component, are tapered in regions of steep neutral
! slope according to the requirements of linear stability.  MOM allows for 
! choice of two tapering schemes:
!
! (a) the tanh taper of Danabasoglu and McWilliams (1995)
! (b) the quadratic scheme of Gerdes, Koberle, and Willebrand (1991)
!
! Linear stability is far less stringent on the diagonal (1,1) and (2,2)
! part of the horizontal flux.  Indeed, these terms in practice need
! not be tapered in steep sloped regions. The namelist 
! neutral_taper_diagonal=.false. keeps the diagnonal terms maintained 
! for all neutral slopes. This approach assists in reducing numerical
! noise in regions where the physical system experiences a lot of
! diapycnal mixing anyhow. 
!
!---How the skew diffusive flux components are computed:
!
! The GM skew flux components are purely off-diagonal.  
! They are generally tapered when neutral slope 
! is large (neutral_physics_simple=.false).
! Doing so maintains a nontrivial GM slumping effect even when the 
! neutral slopes are vertical.  The alternative neutral_physics_simple=.true.
! is the approach used in MOM3, whereby GM effects are removed 
! in steep sloped regions.  neutral_physics_simple=.false. is 
! less efficient, but has been seen to yield superior simulations.
!
! </DESCRIPTION>
!
subroutine nphysicsA (Time, Thickness, Dens, rho, T_prog, &
                      surf_blthick, gm_diffusivity,       &
                      baroclinic_rossby_radius)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  real, dimension(isd:,jsd:,:), intent(in)    :: rho
  real, dimension(isd:,jsd:),   intent(in)    :: surf_blthick
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  real, dimension(isd:,jsd:,:), intent(inout) :: gm_diffusivity
  real, dimension(isd:,jsd:),   intent(inout) :: baroclinic_rossby_radius

  real, dimension(isd:ied,jsd:jed) :: tchg

  integer :: i, j, k, n, nn
  integer :: tau, taum1

  if (.not. use_this_module) return

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_nphysicsA (neutral_physics): needs initialization')
  endif 

  if (size(T_prog(:)) /= num_prog_tracers) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_nphysicsA (neutral_physics): inconsistent size of T_prog')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  ! time dependent delqc geometric factor 
  do k=1,nk
    delqc(:,:,k,0) = Grd%fracdz(k,0)*Thickness%rho_dzt(:,:,k,tau)
    delqc(:,:,k,1) = Grd%fracdz(k,1)*Thickness%rho_dzt(:,:,k,tau)
  enddo

  ! time dependent inverse dzwt 
  do k=0,nk
    dzwtr(:,:,k) = 1.0/Thickness%dzwt(:,:,k) 
  enddo

  ! compute density derivatives 
  do k=1,nk
     drhodT(:,:,k) = Dens%drhodT(:,:,k)
     drhodS(:,:,k) = Dens%drhodS(:,:,k)
  enddo

  call tracer_derivs(Time, taum1, dtime, drhodT, drhodS, T_prog, Dens, dzwtr, &
                     dTdx, dTdy, dTdz, dSdx, dSdy, dSdz, drhodzb, drhodzh) 

  call neutral_slopes(Time, dTdx, dTdy, dSdx, dSdy, drhodT, drhodS, drhodzb, tensor_31, tensor_32)

  call cabbeling_thermob_tendency(Time, Thickness, T_prog, Dens,                               &
  dTdx(index_temp)%field(:,:,:), dTdy(index_temp)%field(:,:,:), dTdz(index_temp)%field(:,:,:), &
  dSdx%field(:,:,:), dSdy%field(:,:,:),                                                        &
  drhodzh, dxtdtsn, dtwedyt, dzwtr, delqc, aredi_array)

  call compute_eta_tend_gm90(Time, Thickness, Dens,             &
  dTdx(index_temp)%field(:,:,:), dTdy(index_temp)%field(:,:,:), &
  dSdx%field(:,:,:), dSdy%field(:,:,:),                         &
  drhodzh, dtwedyt, dzwtr, delqc, tensor_31, tensor_32, agm_array)

  call neutral_blayer(Time, surf_blthick)


  if(diagnose_gm_redi) then 
      do n=1,num_prog_tracers  
         fz1_redi(n)%field(:,:) = 0.0
         fz2_redi(n)%field(:,:) = 0.0
         flux_x_redi(n)%field(:,:,:)   = 0.0
         flux_y_redi(n)%field(:,:,:)   = 0.0
         tendency_redi(n)%field(:,:,:) = 0.0
      enddo
      call fz_terms_diag(Time)   ! must be called prior to fz_terms
  endif

  call fz_terms(Time, Thickness, Dens, T_prog, rho)
  do n=1,num_prog_tracers  
    fz1(n)%field(:,:)      = 0.0
    fz2(n)%field(:,:)      = 0.0
    flux_x(n)%field(:,:,:) = 0.0
    flux_y(n)%field(:,:,:) = 0.0
    T_prog(n)%wrk1(:,:,:)  = 0.0
  enddo 

  if(neutral_physics_simple) then 

      ! horizontal flux components are downgradient diffusive 
      do k=1,nk
         do nn=1,num_prog_tracers
            flux_x(nn)%field(:,:,k) = &
                 agm_array(:,:,k)*dTdx(nn)%field(:,:,k)*FMX(Thickness%rho_dzt(:,:,k,tau)*Grd%tmask(:,:,k))
            flux_y(nn)%field(:,:,k) = &
                 agm_array(:,:,k)*dTdy(nn)%field(:,:,k)*FMY(Thickness%rho_dzt(:,:,k,tau)*Grd%tmask(:,:,k))
         enddo
      enddo
      if (Grd%tripolar) then 
          do nn=1,num_prog_tracers
             call mpp_update_domains(flux_x(nn)%field(:,:,:), flux_y(nn)%field(:,:,:), Dom_flux%domain2d, &
                  gridtype=CGRID_NE, complete=T_prog(nn)%complete) 
          enddo
      endif
      do k=1,nk
         call mpp_clock_begin(id_clock_fz_flux)
         call fz_flux(T_prog,k)
         call mpp_clock_end(id_clock_fz_flux)

         do nn=1,num_prog_tracers
            tchg(:,:)  = (BDX_ET(flux_x(nn)%field(:,:,k)) + BDY_NT(flux_y(nn)%field(:,:,k))) 
            T_prog(nn)%wrk1(COMP,k) = Grd%tmask(COMP,k)*(tchg(COMP) + (fz1(nn)%field(COMP)-fz2(nn)%field(COMP)))
            T_prog(nn)%th_tendency(COMP,k) = T_prog(nn)%th_tendency(COMP,k) + T_prog(nn)%wrk1(COMP,k)
            fz1(nn)%field(COMP)       = fz2(nn)%field(COMP)
            flux_x(nn)%field(COMP,k)  = Grd%dyte(COMP)*flux_x(nn)%field(COMP,k)
            flux_y(nn)%field(COMP,k)  = Grd%dxtn(COMP)*flux_y(nn)%field(COMP,k)
         enddo
      enddo

  else  ! neutral_physics_simple=.false.

      do k=1,nk
         call mpp_clock_begin(id_clock_fx_flux)
         call fx_flux(Time, Thickness, T_prog, k)
         call mpp_clock_end(id_clock_fx_flux)

         call mpp_clock_begin(id_clock_fy_flux)
         call fy_flux(Time, Thickness, T_prog, k)
         call mpp_clock_end(id_clock_fy_flux)
      enddo

      if (Grd%tripolar) then 
          do nn=1,num_prog_tracers
             call mpp_update_domains(flux_x(nn)%field(:,:,:), flux_y(nn)%field(:,:,:), &
                  Dom_flux%domain2d, gridtype=CGRID_NE, complete=T_prog(nn)%complete) 
          enddo
      endif

      do k=1,nk
         call mpp_clock_begin(id_clock_fz_flux)
         call fz_flux(T_prog,k)
         call mpp_clock_end(id_clock_fz_flux)

         do nn=1,num_prog_tracers
            do j=jsc,jec
               do i=isc,iec
                  T_prog(nn)%wrk1(i,j,k) =                                                 & 
                       Grd%tmask(i,j,k)                                                    &
                       *(fz1(nn)%field(i,j)-fz2(nn)%field(i,j)                             &
                       +(flux_x(nn)%field(i,j,k)-flux_x(nn)%field(i-1,j,k)                 &
                       + flux_y(nn)%field(i,j,k)-flux_y(nn)%field(i,j-1,k) )*Grd%datr(i,j) &
                       )
                  T_prog(nn)%th_tendency(i,j,k) = T_prog(nn)%th_tendency(i,j,k) + T_prog(nn)%wrk1(i,j,k)
                  fz1(nn)%field(i,j) = fz2(nn)%field(i,j)
               enddo
            enddo
         enddo
      enddo  ! enddo for k=1,nk


      ! some extra diagnostics to split GM and Redi pieces 
      if(diagnose_gm_redi) then 

          do k=1,nk

             call mpp_clock_begin(id_clock_fx_flux_diag)
             call fx_flux_diag(Time, Thickness, T_prog, k)
             call mpp_clock_end(id_clock_fx_flux_diag)

             call mpp_clock_begin(id_clock_fy_flux_diag) 
             call fy_flux_diag(Time, Thickness, T_prog, k)
             call mpp_clock_end(id_clock_fy_flux_diag) 

          enddo
          if (Grd%tripolar) then 
              do nn=1,num_prog_tracers
                 call mpp_update_domains(flux_x_redi(nn)%field(:,:,:), flux_y_redi(nn)%field(:,:,:), &
                      Dom_flux%domain2d, gridtype=CGRID_NE, complete=T_prog(nn)%complete) 
              enddo
          endif

          do k=1,nk
             call mpp_clock_begin(id_clock_fz_flux_diag)
             call fz_flux_diag(T_prog,k)
             call mpp_clock_end(id_clock_fz_flux_diag)
             do nn=1,num_prog_tracers
                do j=jsc,jec
                   do i=isc,iec
                      tendency_redi(nn)%field(i,j,k) =                                                   & 
                           Grd%tmask(i,j,k)                                                              &
                           *(fz1_redi(nn)%field(i,j)-fz2_redi(nn)%field(i,j)                             &
                           +(flux_x_redi(nn)%field(i,j,k)-flux_x_redi(nn)%field(i-1,j,k)                 &
                           + flux_y_redi(nn)%field(i,j,k)-flux_y_redi(nn)%field(i,j-1,k) )*Grd%datr(i,j) &
                           )
                      fz1_redi(nn)%field(i,j) = fz2_redi(nn)%field(i,j)
                   enddo
                enddo
             enddo
          enddo  ! enddo for k=1,nk

      endif  ! endif for diagnose_gm_redi


  endif  ! endif for neutral_physics_simple


  ! compute Eady growth rate and baroclinicity for use in next time step 
  call compute_eady_rate(Time, Thickness, T_prog, Dens, eady_termx, eady_termy)
  call compute_baroclinicity(Time, baroclinic_termx, baroclinic_termy)

  ! compute rossby radius for use in next time step 
  call compute_rossby_radius(Thickness, dTdz, dSdz, Time, drhodT, drhodS, &
                             rossby_radius, rossby_radius_raw)
  baroclinic_rossby_radius(:,:) = rossby_radius_raw(:,:)

  ! compute width of baroclinic zone for use in next time step 
  call compute_bczone_radius(Time, bczone_radius)

  ! update closure-based diffusivity for next time step 
  call compute_diffusivity(Time, ksurf_blayer, Dens%drhodz_zt, rossby_radius, &
                           bczone_radius, agm_array, aredi_array)

  ! gm_diffusivity passed to neutral_physics for use in computing form drag
  gm_diffusivity(:,:,:) = agm_array(:,:,:)

  call nphysics_diagnostics(Time, T_prog, Dens) 

  call gm_velocity(Thickness, Time)


end subroutine nphysicsA
! </SUBROUTINE> NAME="nphysicsA"



!#######################################################################
! <SUBROUTINE NAME="neutral_blayer">
!
! <DESCRIPTION>
! Subroutine computes the boundary layer as determined by 
! 1. steep neutral slopes
! 2. depth within which typical mesoscale eddies are partially outcropped
! 3. depth within which vertical mixing scheme (e.g., kpp) computes a boundary layer
!
! Note: Only consider surface boundary layers here.  
!
! Scheme originally coded for MOM4.0 by Stephen.Griffies
! with help for optimization by Russell.Fiedler@csiro.au.
! 
! </DESCRIPTION>
!
subroutine neutral_blayer(Time, surf_blthick)

  type(ocean_time_type),      intent(in) :: Time
  real, dimension(isd:,jsd:), intent(in) :: surf_blthick

  logical :: within_interior
  logical :: slopeok(isc:iec)
  integer :: i, j, k, ip, jq, kb, kr, kk, kkpkr
  real :: depth, slope, absslope 
  real :: steep_depth(isd:ied,jsd:jed)
  real :: slope_blayer_baseA(isd:ied,jsd:jed)
  real :: slope_blayer_baseB(isd:ied,jsd:jed)
  real :: thick_31(isd:ied,jsd:jed,0:1,0:1)
  real :: thick_32(isd:ied,jsd:jed,0:1,0:1)
  real :: thick_13(isd:ied,jsd:jed,0:1,0:1)
  real :: thick_23(isd:ied,jsd:jed,0:1,0:1)
  real :: slope_31(isd:ied,jsd:jed,0:1,0:1)
  real :: slope_32(isd:ied,jsd:jed,0:1,0:1)
  real :: slope_13(isd:ied,jsd:jed,0:1,0:1)
  real :: slope_23(isd:ied,jsd:jed,0:1,0:1)

  call mpp_clock_begin(id_clock_neutral_blayer)

  thick_31 = Grd%zw(1)
  thick_32 = Grd%zw(1)
  thick_13 = Grd%zw(1)
  thick_23 = Grd%zw(1)
  slope_31 = 0.0
  slope_32 = 0.0
  slope_13 = 0.0
  slope_23 = 0.0

  slope_blayer_baseA = 0.0
  slope_blayer_baseB = 0.0
  steep_depth        = 0.0

  eddy_depth         = 0.0
  slope_blayer_base  = 0.0 
  depth_blayer_base  = 0.0


  if(neutral_sine_taper .or. neutral_blayer_diagnose) then 

  ! Determine depth over which mesoscale eddies feel the ocean 
  ! surface.  This depth is a function of the neutral slope 
  ! and the Rossby radius.  This depth is called "eddy_depth".
  ! The algorithm for computing this depth is taken from 
  ! the appendix to Large etal, 1997 JPO vol 27, 2418-2447. 
  ! 
  ! In addition to considering mesoscale eddy lengths,
  ! include the possibility that the diabatic vertical
  ! mixing (e.g., kpp) produces a mixed layer depth that is 
  ! deeper than the depth that mesoscale eddies feel the ocean 
  ! surface.  Include this surf_blthick in the considerations so 
  ! to determine the depth of this generalized "boundary layer" 
  ! and the neutral slope at the base of the boundary layer. 

  call mpp_clock_begin(id_clock_neutral_blayer_1)

! 31-triads 
      do ip=0,1
         do kr=0,1
            do j=jsc,jec
               slopeok(:)=.false.
kkloop31a:     do kk=1,nk-1
                  do i=isc,iec 

                     absslope = abs(tensor_31(i,j,kk,ip,kr))
                     depth    = max(rossby_radius(i,j)*absslope, max(turb_blayer_min,surf_blthick(i,j)))
                     if(depth > Grd%zw(kk) .and. .not. slopeok(i)) then
                         thick_31(i,j,ip,kr) = Grd%zw(kk)
                         slope_31(i,j,ip,kr) = absslope
                     else 
                         slopeok(i)=.true.
                     endif
                  enddo
                  if(all(slopeok(:))) exit kkloop31a

               enddo kkloop31a
            enddo
         enddo
      enddo

      ! 13-triads 
      do kr=0,1
         do ip=0,1 
            do j=jsc,jec
               slopeok(:)=.false.

kkloop13a:     do kk=1,nk
                  kkpkr = min(kk+kr,nk)
                  do i=isc,iec

                     slope = -Grd%tmask(i+ip,j,kkpkr)                                 &
                          *(drhodT(i+ip,j,kk)*dTdx(index_temp)%field(i,j,kk)         +&
                            drhodS(i+ip,j,kk)*dSdx%field(i,j,kk))                     &
                          /(drhodT(i+ip,j,kk)*dTdz(index_temp)%field(i+ip,j,kk-1+kr) +&
                            drhodS(i+ip,j,kk)*dSdz%field(i+ip,j,kk-1+kr) - epsln)
                           
                     absslope = abs(slope)  

                     depth = max(rossby_radius(i,j)*absslope, max(turb_blayer_min,surf_blthick(i,j)))
                     if(depth > Grd%zw(kk) .and. .not.slopeok(i)) then 
                         thick_13(i,j,ip,kr) = Grd%zw(kk)
                         slope_13(i,j,ip,kr) = absslope
                     else 
                         slopeok(i)=.true.
                     endif

                  enddo
                  if(all(slopeok(:))) exit kkloop13a
               enddo kkloop13a

            enddo
         enddo
      enddo

      ! 32-triads 
      do jq=0,1
         do kr=0,1
            do j=jsc,jec
               slopeok(:)=.false.
kkloop32a:      do kk=1,nk-1
                  do i=isc,iec 

                     absslope = abs(tensor_32(i,j,kk,jq,kr))
                     depth    = max(rossby_radius(i,j)*absslope, max(turb_blayer_min,surf_blthick(i,j)))
                     if(depth > Grd%zw(kk) .and. .not. slopeok(i) ) then
                         thick_32(i,j,jq,kr) = Grd%zw(kk)
                         slope_32(i,j,jq,kr) = absslope
                     else 
                         slopeok(i) = .true.
                     endif
                  enddo
                  if(all(slopeok(:))) exit kkloop32a
               enddo kkloop32a
            enddo
         enddo
      enddo

      ! 23-triads 
      do kr=0,1
         do jq=0,1  
            do j=jsc,jec
               slopeok(:)=.false.

kkloop23a:     do kk=1,nk
                  do i=isc,iec
                     kkpkr = min(kk+kr,nk)

                     slope = -Grd%tmask(i,j+jq,kkpkr)                                 &
                          *(drhodT(i,j+jq,kk)*dTdy(index_temp)%field(i,j,kk)+         &
                            drhodS(i,j+jq,kk)*dSdy%field(i,j,kk))                     & 
                          /(drhodT(i,j+jq,kk)*dTdz(index_temp)%field(i,j+jq,kk-1+kr) +&
                            drhodS(i,j+jq,kk)*dSdz%field(i,j+jq,kk-1+kr) - epsln)  
                     absslope = abs(slope)  

                     depth = max(rossby_radius(i,j)*absslope, max(turb_blayer_min,surf_blthick(i,j)))
                     if(depth > Grd%zw(kk) .and. .not.slopeok(i)) then 
                         thick_23(i,j,jq,kr) = Grd%zw(kk)
                         slope_23(i,j,jq,kr) = absslope
                     else
                         slopeok(i)=.true.
                     endif

                  enddo
                  if(all(slopeok(:))) exit kkloop23a
               enddo kkloop23a
            enddo
         enddo
      enddo

  call mpp_clock_end(id_clock_neutral_blayer_1)


  call mpp_clock_begin(id_clock_neutral_blayer_2)
  ! max of triad depth defines eddy_depth. average slope
  ! defines slope at boundary layer base.  This method 
  ! maximizes eddy_depth, and regulates slope at base.
      do j=jsc,jec
         do i=isc,iec
            eddy_depth(i,j) = max(Grd%zw(1), &
                                  thick_31(i,j,0,0), thick_31(i,j,0,1), &
                                  thick_31(i,j,1,0), thick_31(i,j,1,1), &
                                  thick_32(i,j,0,0), thick_32(i,j,0,1), &
                                  thick_32(i,j,1,0), thick_32(i,j,1,1), &
                                  thick_13(i,j,0,0), thick_13(i,j,0,1), &
                                  thick_13(i,j,1,0), thick_13(i,j,1,1), &
                                  thick_23(i,j,0,0), thick_23(i,j,0,1), &
                                  thick_23(i,j,1,0), thick_23(i,j,1,1))

            ! make sure that none of the slopes are larger than smax 
            do kr=0,1
              do ip=0,1      
                 slope_31(i,j,ip,kr) = min(smax, slope_31(i,j,ip,kr)) 
                 slope_13(i,j,ip,kr) = min(smax, slope_13(i,j,ip,kr)) 
                 slope_32(i,j,ip,kr) = min(smax, slope_32(i,j,ip,kr)) 
                 slope_23(i,j,ip,kr) = min(smax, slope_23(i,j,ip,kr)) 
               enddo
            enddo

            ! average of the 16 slopes defines slope at boundary layer base  
                slope_blayer_baseA(i,j) = &
                                         (slope_31(i,j,0,0) + slope_13(i,j,0,0) + &
                                          slope_31(i,j,1,0) + slope_13(i,j,1,0) + &
                                          slope_31(i,j,0,1) + slope_13(i,j,0,1) + &
                                          slope_31(i,j,1,1) + slope_13(i,j,1,1) + &
                                          slope_32(i,j,0,0) + slope_23(i,j,0,0) + &
                                          slope_32(i,j,0,1) + slope_23(i,j,0,1) + &
                                          slope_32(i,j,1,0) + slope_23(i,j,1,0) + &
                                          slope_32(i,j,1,1) + slope_23(i,j,1,1)) /16.0
         enddo
      enddo

  call mpp_clock_end(id_clock_neutral_blayer_2)

  endif  ! endif for neutral_sine_taper .or. neutral_blayer_diagnose


  
  if(neutral_linear_gm_taper .or. neutral_blayer_diagnose) then 

  ! Determine depth of surface boundary layer as defined 
  ! by depth where neutral slopes are larger than smax.
      
      thick_31 = Grd%zw(1)
      thick_32 = Grd%zw(1)
      thick_13 = Grd%zw(1)
      thick_23 = Grd%zw(1)

      ! 31-triads 
      do ip=0,1
         do kr=0,1
            do j=jsc,jec
               slopeok(:)=.false.
kkloop31b:     do kk=1,nk-1
                  do i=isc,iec 

                     absslope = abs(tensor_31(i,j,kk,ip,kr))
                     if(absslope > smax .and. .not. slopeok(i)) then
                         thick_31(i,j,ip,kr) = Grd%zw(kk)
                         slope_31(i,j,ip,kr) = absslope
                     else 
                         slopeok(i)=.true.
                     endif
                  enddo
                  if(all(slopeok(:))) exit kkloop31b

               enddo kkloop31b
            enddo
         enddo
      enddo

      ! 13-triads 
      do kr=0,1
         do ip=0,1 
            do j=jsc,jec
               slopeok(:)=.false.

kkloop13b:     do kk=1,nk
                  kkpkr = min(kk+kr,nk)
                  do i=isc,iec

                     slope = -Grd%tmask(i+ip,j,kkpkr)&
                          *(drhodT(i+ip,j,kk)*dTdx(index_temp)%field(i,j,kk)         +&
                            drhodS(i+ip,j,kk)*dSdx%field(i,j,kk))                     &
                          /(drhodT(i+ip,j,kk)*dTdz(index_temp)%field(i+ip,j,kk-1+kr) +&
                            drhodS(i+ip,j,kk)*dSdz%field(i+ip,j,kk-1+kr) - epsln)
                           
                     absslope = abs(slope) 
 
                     if(absslope > smax .and. .not.slopeok(i)) then 
                         thick_13(i,j,ip,kr) = Grd%zw(kk)
                         slope_13(i,j,ip,kr) = absslope
                     else 
                         slopeok(i)=.true.
                     endif

                  enddo
                  if(all(slopeok(:))) exit kkloop13b
               enddo kkloop13b

            enddo

         enddo
      enddo

      ! 32-triads 
      do jq=0,1
         do kr=0,1
            do j=jsc,jec
               slopeok(:)=.false.
kkloop32b:     do kk=1,nk-1
                  do i=isc,iec 

                     absslope = abs(tensor_32(i,j,kk,jq,kr))
                     if(absslope > smax .and. .not. slopeok(i)) then
                         thick_32(i,j,jq,kr) = Grd%zw(kk)
                         slope_32(i,j,jq,kr) = absslope
                     else 
                         slopeok(i)=.true.
                     endif
                  enddo
                  if(all(slopeok(:))) exit kkloop32b

               enddo kkloop32b
            enddo
         enddo
      enddo

      ! 23-triads 
      do kr=0,1
         do jq=0,1  
            do j=jsc,jec
               slopeok(:)=.false.

kkloop23b:     do kk=1,nk
                  do i=isc,iec
                     kkpkr = min(kk+kr,nk)

                     slope = -Grd%tmask(i,j+jq,kkpkr)&
                          *(drhodT(i,j+jq,kk)*dTdy(index_temp)%field(i,j,kk)         +&
                            drhodS(i,j+jq,kk)*dSdy%field(i,j,kk))                     & 
                          /(drhodT(i,j+jq,kk)*dTdz(index_temp)%field(i,j+jq,kk-1+kr) +&
                            drhodS(i,j+jq,kk)*dSdz%field(i,j+jq,kk-1+kr) - epsln)  
                     absslope = abs(slope)  

                     if(absslope > smax .and. .not.slopeok(i)) then 
                         thick_23(i,j,jq,kr) = Grd%zw(kk)
                         slope_23(i,j,jq,kr) = absslope
                     else
                         slopeok(i)=.true.
                     endif

                  enddo
                  if(all(slopeok(:))) exit kkloop23b
               enddo kkloop23b

            enddo
         enddo
      enddo

  call mpp_clock_begin(id_clock_neutral_blayer_3)
  ! max of triad depth defines steep_depth. average slope
  ! defines slope at boundary layer base.  This method 
  ! maximizes steep_depth, and regulates slope at base.

      do j=jsc,jec
         do i=isc,iec

            steep_depth(i,j) =  max(Grd%zw(1), &
                                    thick_31(i,j,0,0), thick_31(i,j,0,1), &
                                    thick_31(i,j,1,0), thick_31(i,j,1,1), &
                                    thick_32(i,j,0,0), thick_32(i,j,0,1), &
                                    thick_32(i,j,1,0), thick_32(i,j,1,1), &
                                    thick_13(i,j,0,0), thick_13(i,j,0,1), &
                                    thick_13(i,j,1,0), thick_13(i,j,1,1), &
                                    thick_23(i,j,0,0), thick_23(i,j,0,1), &
                                    thick_23(i,j,1,0), thick_23(i,j,1,1))

            ! make sure that none of the slopes are larger than smax 
            do kr=0,1
              do ip=0,1      
                 slope_31(i,j,ip,kr) = min(smax, slope_31(i,j,ip,kr)) 
                 slope_13(i,j,ip,kr) = min(smax, slope_13(i,j,ip,kr)) 
                 slope_32(i,j,ip,kr) = min(smax, slope_32(i,j,ip,kr)) 
                 slope_23(i,j,ip,kr) = min(smax, slope_23(i,j,ip,kr)) 
               enddo
            enddo

            ! average of the 16 slopes defines slope at boundary layer base  
            slope_blayer_baseB(i,j) = 0.0
            do kr=0,1
              do ip=0,1      
                slope_blayer_baseB(i,j) = slope_blayer_baseB(i,j)                   + &
                                          slope_31(i,j,ip,kr) + slope_13(i,j,ip,kr) + &
                                          slope_32(i,j,ip,kr) + slope_23(i,j,ip,kr)
               enddo
            enddo
            slope_blayer_baseB(i,j) = slope_blayer_baseB(i,j)/16.0

         enddo
      enddo

  call mpp_clock_end(id_clock_neutral_blayer_3)

  ! The maximum of steep_depth and eddy_depth define 
  ! neutral_blayer_depth.  For depths shallower than 
  ! neutral_blayer_depth, horizontal gm velocity is depth  
  ! independent and vertical gm velocity linearly decreases
  ! to zero as the surface is reached.   

      do j=jsc,jec
         do i=isc,iec

            if(steep_depth(i,j) > eddy_depth(i,j)) then 
                slope_blayer_base(i,j) = slope_blayer_baseB(i,j)
                depth_blayer_base(i,j) = steep_depth(i,j)
            else 
                slope_blayer_base(i,j) = slope_blayer_baseA(i,j)
                depth_blayer_base(i,j) = eddy_depth(i,j)
            endif

         enddo
      enddo

      ! save the k-level surface transition layer. 
      ! initialize ksurf_blayer to 1 rather than 0, to avoid 
      ! accessing the zeroth element of an array.   
      ksurf_blayer(:,:) = 1
      wrk1_2d(:,:)      = 0.0
      do j=jsd,jed
         do i=isd,ied
            kb=Grd%kmt(i,j)
            within_interior=.false.
            ksurf_blayer_loop: do k=1,kb
               if(Grd%zw(k) > depth_blayer_base(i,j)) then 
                   ksurf_blayer(i,j) = k
                   wrk1_2d(i,j)      = float(k)
                   exit ksurf_blayer_loop
               endif
            enddo ksurf_blayer_loop
         enddo
      enddo


  endif  ! endif for neutral_linear_gm_taper .or. neutral_blayer_diagnose

  call diagnose_2d(Time, Grd, id_depth_blayer_base, depth_blayer_base(:,:))
  call diagnose_2d(Time, Grd, id_slope_blayer_base, slope_blayer_base(:,:))
  call diagnose_2d(Time, Grd, id_eddy_depth, eddy_depth(:,:))
  call diagnose_2d(Time, Grd, id_steep_depth, steep_depth(:,:))

  ! update domains or zero the depths, depending on the chosen method 
  if(neutral_linear_gm_taper) then 
      call mpp_update_domains(depth_blayer_base, Dom%domain2d) 
      call mpp_update_domains(slope_blayer_base, Dom%domain2d) 
  else 
      depth_blayer_base = 0.0
  endif
  if(neutral_sine_taper) then 
      call mpp_update_domains(eddy_depth, Dom%domain2d) 
  else 
      eddy_depth = 0.0
  endif


  call mpp_clock_end(id_clock_neutral_blayer)

end subroutine neutral_blayer
! </SUBROUTINE> NAME="neutral_blayer"


!#######################################################################
! <SUBROUTINE NAME="fz_terms">
!
! <DESCRIPTION>
! Subroutine computes the tracer independent pieces of the vertical 
! flux component. As a result of this routine, 
! Array tensor_31 = x-diffusivity*slope (m^2/sec) for fz
! Array tensor_32 = y-diffusivity*slope (m^2/sec) for fz 
!
! K33 is the (3,3) term in small angle Redi diffusion tensor.
! It is broken into an explicit in time piece and implicit 
! in time piece.  It is weighted by density for non-Boussinesq
! and rho0 for Boussinesq.  
!
! K33 has units (kg/m^3)*m^2/sec.
!
! Also will compute the squared Eady growth rate, with the maximum
! slope contributing to this growth rate set by smax.
! </DESCRIPTION>
!
subroutine fz_terms(Time, Thickness, Dens, T_prog, rho)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  real, dimension(isd:,jsd:,:), intent(in)    :: rho

  integer :: i, j, k, n, ip, jq, kr, tau

  real :: baroclinic_triad, K33, K33_crit
  real :: sumx(0:1), sumy(0:1)

  real :: slope, absslope, depth_ratio
  real :: taperA, taperB
  real :: taper_slope, taper_slope2, gm_taper_slope 

  real :: aredi_scalar, agm_scalar, aredi_plus_agm 
  real :: absdrhodzb(0:1)
  real :: mindelqc(0:1,0:1), delqc_ijk_1, delqc_ijkp1_0

  call mpp_clock_begin(id_clock_fz_terms)

  tau = Time%tau
  
  eady_termx(:,:,:)     = 0.0
  eady_termy(:,:,:)     = 0.0
  baroclinic_termx(:,:) = 0.0
  baroclinic_termy(:,:) = 0.0
  tx_trans_gm(:,:,:)    = 0.0
  ty_trans_gm(:,:,:)    = 0.0


! Main loop. Short ip,jq,kr loops are explicitly unrolled
! in order to expose independence and allow vectorisation

  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec

           aredi_scalar    = aredi_array(i,j,k)
           agm_scalar      = gm_skew*agm_array(i,j,k)
           aredi_plus_agm  = aredi_scalar + agm_scalar

           absdrhodzb(0) = abs(drhodzb(i,j,k,0))
           absdrhodzb(1) = abs(drhodzb(i,j,k,1))

           delqc_ijk_1      = delqc(i,j,k,1)
           delqc_ijkp1_0    = delqc(i,j,k+1,0)
           
           !mindelqc(ip,kr) = min(delqc(i-1+ip,j,k+kr,1-kr),delqc(i+ip,j,k+kr,1-kr)) 
           ip=0 ; kr=0
           mindelqc(ip,kr)  = min(delqc(i-1,j,k,1),delqc_ijk_1) 
           ip=0 ; kr=1
           mindelqc(ip,kr)  = min(delqc(i-1,j,k+1,0),delqc_ijkp1_0) 
           ip=1 ; kr=0
           mindelqc(ip,kr)  = min(delqc_ijk_1,delqc(i+1,j,k,1)) 
           ip=1 ; kr=1
           mindelqc(ip,kr)  = min(delqc_ijkp1_0,delqc(i+1,j,k+1,0)) 

           tx_trans_gm(i,j,k) = 0.0
           baroclinic_triad   = 0.0

           ip=0   
              sumx(ip) = 0.0

              kr=0
                 slope    = tensor_31(i,j,k,ip,kr)
                 absslope = abs(slope)
                 
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 tx_trans_gm(i,j,k)          = tx_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope       
                 tensor_31(i,j,k,ip,kr)      = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumx(ip)                    = sumx(ip)                 + mindelqc(ip,kr)*taper_slope2
                 eady_termx(i,j,k)           = eady_termx(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad            = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)

              kr=1
                 slope    = tensor_31(i,j,k,ip,kr)
                 absslope = abs(slope)
                           
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 tx_trans_gm(i,j,k)          = tx_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope      
                 tensor_31(i,j,k,ip,kr)      = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumx(ip)                    = sumx(ip)                 + mindelqc(ip,kr)*taper_slope2
                 eady_termx(i,j,k)           = eady_termx(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad            = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)
                 sumx(ip)                    = sumx(ip)*dtew(i,j,ip)

           ip=1   
              sumx(ip) = 0.0

              kr=0
                 slope        = tensor_31(i,j,k,ip,kr)
                 absslope     = abs(slope)
                           
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 tx_trans_gm(i,j,k)          = tx_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope    
                 tensor_31(i,j,k,ip,kr)      = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumx(ip)                    = sumx(ip)                 + mindelqc(ip,kr)*taper_slope2
                 eady_termx(i,j,k)           = eady_termx(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad            = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)

              kr=1
                 slope        = tensor_31(i,j,k,ip,kr)
                 absslope     = abs(slope)
                                                     
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 tx_trans_gm(i,j,k)          = tx_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope   
                 tensor_31(i,j,k,ip,kr)      = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumx(ip)                    = sumx(ip)                 + mindelqc(ip,kr)*taper_slope2
                 eady_termx(i,j,k)           = eady_termx(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad            = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)
                 sumx(ip)                    = sumx(ip)*dtew(i,j,ip)

           if(Thickness%depth_zt(i,j,k) >= agm_closure_upper_depth .and. &
              Thickness%depth_zt(i,j,k) <= agm_closure_lower_depth) then 
             baroclinic_termx(i,j) = baroclinic_termx(i,j) + baroclinic_triad*Thickness%dzwt(i,j,k)
           endif 

           tx_trans_gm(i,j,k)    = 0.25*tx_trans_gm(i,j,k)*Grd%dyu(i,j)

           !mindelqc(jq,kr) = min(delqc(i,j-1+jq,k+kr,1-kr),delqc(i,j+jq,k+kr,1-kr)) 
           jq=0 ; kr=0
           mindelqc(jq,kr)  = min(delqc(i,j-1,k,1),delqc_ijk_1) 
           jq=0 ; kr=1
           mindelqc(jq,kr)  = min(delqc(i,j-1,k+1,0),delqc_ijkp1_0) 
           jq=1 ; kr=0
           mindelqc(jq,kr)  = min(delqc_ijk_1,delqc(i,j+1,k,1)) 
           jq=1 ; kr=1
           mindelqc(jq,kr)  = min(delqc_ijkp1_0,delqc(i,j+1,k+1,0)) 

           ty_trans_gm(i,j,k) = 0.0
           baroclinic_triad   = 0.0 

           jq=0   
              sumy(jq) = 0.0

               kr=0
                 slope        = tensor_32(i,j,k,jq,kr)
                 absslope     = abs(slope)
                           
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 ty_trans_gm(i,j,k)          = ty_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope       
                 tensor_32(i,j,k,jq,kr)      = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumy(jq)                    = sumy(jq)                 + mindelqc(jq,kr)*taper_slope2
                 eady_termy(i,j,k)           = eady_termy(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad            = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)

              kr=1
                 slope        = tensor_32(i,j,k,jq,kr)
                 absslope     = abs(slope)
                           
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 ty_trans_gm(i,j,k)          = ty_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope      
                 tensor_32(i,j,k,jq,kr)      = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumy(jq)                    = sumy(jq)                 + mindelqc(jq,kr)*taper_slope2
                 eady_termy(i,j,k)           = eady_termy(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad            = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)
                 sumy(jq)                    = sumy(jq)*dtns(i,j,jq)

           jq=1   
              sumy(jq) = 0.0

               kr=0
                 slope    = tensor_32(i,j,k,jq,kr)
                 absslope = abs(slope)

                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 ty_trans_gm(i,j,k)          = ty_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope     
                 tensor_32(i,j,k,jq,kr)      = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumy(jq)                    = sumy(jq)                 + mindelqc(jq,kr)*taper_slope2
                 eady_termy(i,j,k)           = eady_termy(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad            = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)

              kr=1
                 slope    = tensor_32(i,j,k,jq,kr)
                 absslope = abs(slope)

                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 ty_trans_gm(i,j,k)          = ty_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope      
                 tensor_32(i,j,k,jq,kr)      = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumy(jq)                    = sumy(jq)                 + mindelqc(jq,kr)*taper_slope2
                 eady_termy(i,j,k)           = eady_termy(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad            = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)

                 sumy(jq) = sumy(jq)*dtns(i,j,jq)

           if(Thickness%depth_zt(i,j,k) >= agm_closure_upper_depth .and. &
              Thickness%depth_zt(i,j,k) <= agm_closure_lower_depth) then
             baroclinic_termy(i,j) = baroclinic_termy(i,j) + baroclinic_triad*Thickness%dzwt(i,j,k)
           endif 

           ty_trans_gm(i,j,k)    = 0.25*ty_trans_gm(i,j,k)*Grd%dxu(i,j)
           
           K33 = aredi_scalar*Grd%tmask(i,j,k+1)*(Grd%dxtr(i,j)*(sumx(0)+sumx(1)) &
                            + Grd%dytr(i,j)*(sumy(0)+sumy(1)))*dzwtr(i,j,k)

           ! determine part of K33 included explicitly in time and that part implicitly. 
           ! explicit calculation is more accurate than implicit, so aim to compute as much 
           ! as stably possible via the explicit method. 
           K33_explicit(i,j,k) = K33
           K33_implicit(i,j,k) = 0.0 
           if(.not. diffusion_all_explicit) then 
               K33_crit = two_dtime_inv*Thickness%rho_dzt(i,j,k,tau)*Thickness%dzt(i,j,k)
               if(K33 >= K33_crit) then 
                   K33_explicit(i,j,k) = K33_crit
                   K33_implicit(i,j,k) = K33-K33_crit
               endif
           endif

        enddo
     enddo
  enddo

  if(tmask_sigma_on) then
     do j=jsc,jec
        do i=isc,iec
           if(tmask_sigma(i,j) > 0.0) then 
              k = Grd%kmt(i,j)-1
              K33_implicit(i,j,k) = K33_implicit(i,j,k)*(1.0-tmask_sigma(i,j))
              K33_explicit(i,j,k) = K33_explicit(i,j,k)*(1.0-tmask_sigma(i,j))
           endif
        enddo
     enddo
 endif
 
 if(tmask_neutral_on) then
     do j=jsc,jec
        do i=isc,iec
           K33_implicit(i,j,1) = 0.0  
           K33_explicit(i,j,1) = 0.0  
           tx_trans_gm(i,j,1)  = 0.0  
           ty_trans_gm(i,j,1)  = 0.0  
           if(Grd%kmt(i,j) > 1) then 
               k = Grd%kmt(i,j)-1
               K33_implicit(i,j,k) = 0.0
               K33_explicit(i,j,k) = 0.0
               tx_trans_gm(i,j,k)  = 0.0  
               ty_trans_gm(i,j,k)  = 0.0  
           endif
        enddo
     enddo
 endif

 do n=1,num_prog_tracers 
   T_prog(n)%K33_implicit(:,:,:) = K33_implicit(:,:,:)
 enddo
 if(neutral_physics_limit) then 
     do n=1,num_prog_tracers 
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                if(T_prog(n)%tmask_limit(i,j,k)==1.0) T_prog(n)%K33_implicit(i,j,k) = 0.0
              enddo
           enddo
        enddo
     enddo
 endif


  do n=1,num_prog_tracers
     if (id_k33_implicit(n) > 0) then
        call diagnose_3d(Time, Grd, id_k33_implicit(n), rho0r*T_prog(n)%K33_implicit(:,:,:))
     endif 
  enddo 


  ! rho factor to get mass transport.
  if(vert_coordinate_class==DEPTH_BASED) then 
      do k=1,nk
         tx_trans_gm(COMP,k) = tx_trans_gm(COMP,k)*rho0
         ty_trans_gm(COMP,k) = ty_trans_gm(COMP,k)*rho0
      enddo
  else
      do k=1,nk
         tx_trans_gm(COMP,k) = tx_trans_gm(COMP,k)*rho(COMP,k)
         ty_trans_gm(COMP,k) = ty_trans_gm(COMP,k)*rho(COMP,k)
      enddo
  endif

  call transport_on_rho_gm (Time, Dens, &
       transport_convert*tx_trans_gm(:,:,:), transport_convert*ty_trans_gm(:,:,:))
  call transport_on_nrho_gm(Time, Dens, &
       transport_convert*tx_trans_gm(:,:,:), transport_convert*ty_trans_gm(:,:,:))
  call transport_on_theta_gm(Time, Dens, T_prog(index_temp), &
       transport_convert*tx_trans_gm(:,:,:), transport_convert*ty_trans_gm(:,:,:))

  ! transport_convert converts kg/s to chosen dimensions (either kg/s or Sv (10^9 kg/s)) 
  if(id_tx_trans_gm > 0) then 
      call mpp_update_domains(tx_trans_gm, Dom%domain2d,flags=EUPDATE)
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               tx_trans_gm(i,j,k) = transport_convert*Grd%dxuer(i,j)*   &
                                   (tx_trans_gm(i+1,j,k)*Grd%due(i,j) + &
                                    tx_trans_gm(i,j,k)  *Grd%duw(i+1,j))
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_tx_trans_gm, tx_trans_gm(:,:,:))
  endif

  ! transport_convert converts kg/s to chosen dimensions (either kg/s or Sv (10^9 kg/s))  
  if(id_ty_trans_gm > 0) then 
      call mpp_update_domains(ty_trans_gm, Dom%domain2d,flags=NUPDATE)
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               ty_trans_gm(i,j,k) = transport_convert*Grd%dyunr(i,j)*           &
                                          (ty_trans_gm(i,j+1,k)*Grd%dun(i,j) +  &
                                           ty_trans_gm(i,j,k)  *Grd%dus(i,j+1)) 
                                        
            enddo
         enddo
      enddo
      call diagnose_3d(Time, Grd, id_ty_trans_gm, ty_trans_gm(:,:,:))
  endif

  call mpp_clock_end(id_clock_fz_terms)

end subroutine fz_terms
! </SUBROUTINE>  NAME="fz_terms"


!#######################################################################
! <SUBROUTINE NAME="fz_terms_diag">
!
! <DESCRIPTION>
! For saving the contributions from GM and Redi separately, it is 
! necessary to compute the tensor_redi component here. 
!
! We do so here, reproducing some lines of code from fz_terms,
! to reduce minimize the need to impinge on the case when NOT
! using this generally expensive (memory and computational)
! diagnostic.  
!
! This routine MUST be called prior to fz_terms, since we use
! tensor_31 and tensor_32 in their raw slope forms here.  
!
! </DESCRIPTION>
!
subroutine fz_terms_diag(Time)

  type(ocean_time_type), intent(in) :: Time

  integer :: i, j, k, ip, jq, kr, tau

  real :: slope, absslope, depth_ratio
  real :: taperA, taperB
  real :: taper_slope
  real :: aredi_scalar

  call mpp_clock_begin(id_clock_fz_terms_diag)

  tau = Time%tau
  
  ! Main loop. Short ip,jq,kr loops are explicitly unrolled
  ! in order to expose independence and allow vectorisation
  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec

           aredi_scalar = aredi_array(i,j,k)

           ip=0   

              kr=0
                 slope    = tensor_31(i,j,k,ip,kr)
                 absslope = abs(slope)
                 
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif
                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif
                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope

                 tensor_31_redi(i,j,k,ip,kr) = aredi_scalar*taper_slope

              kr=1
                 slope    = tensor_31(i,j,k,ip,kr)
                 absslope = abs(slope)
                           
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif
                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif
                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope

                 tensor_31_redi(i,j,k,ip,kr) = aredi_scalar*taper_slope 

           ip=1   

              kr=0
                 slope        = tensor_31(i,j,k,ip,kr)
                 absslope     = abs(slope)
                           
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif
                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif
                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope

                 tensor_31_redi(i,j,k,ip,kr) = aredi_scalar*taper_slope 

              kr=1
                 slope        = tensor_31(i,j,k,ip,kr)
                 absslope     = abs(slope)
                                                     
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif
                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif
                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope

                 tensor_31_redi(i,j,k,ip,kr) = aredi_scalar*taper_slope 

           jq=0   

               kr=0
                 slope        = tensor_32(i,j,k,jq,kr)
                 absslope     = abs(slope)
                           
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif
                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif
                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope

                 tensor_32_redi(i,j,k,jq,kr) = aredi_scalar*taper_slope 

              kr=1
                 slope        = tensor_32(i,j,k,jq,kr)
                 absslope     = abs(slope)
                           
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif
                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif
                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope

                 tensor_32_redi(i,j,k,jq,kr) = aredi_scalar*taper_slope 

           jq=1   

               kr=0
                 slope    = tensor_32(i,j,k,jq,kr)
                 absslope = abs(slope)

                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif
                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif
                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope

                 tensor_32_redi(i,j,k,jq,kr) = aredi_scalar*taper_slope 

              kr=1
                 slope    = tensor_32(i,j,k,jq,kr)
                 absslope = abs(slope)

                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif
                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif
                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope

                 tensor_32_redi(i,j,k,jq,kr) = aredi_scalar*taper_slope 

        enddo
     enddo
  enddo


  call mpp_clock_end(id_clock_fz_terms_diag)

end subroutine fz_terms_diag
! </SUBROUTINE>  NAME="fz_terms_diag"


!#######################################################################
! <SUBROUTINE NAME="fx_flux">
!
! <DESCRIPTION>
! Subroutine computes the i-directed neutral physics tracer flux component.
! Compute this component for all tracers at level k.
!
! fx has physical dimensions (area*diffusivity*density*tracer gradient)
!
! </DESCRIPTION>
!
subroutine fx_flux(Time, Thickness, T_prog, k)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  integer,                      intent(in) :: k

  integer :: nn, i, j, ip, tau
  integer :: kr, kpkr
  real    :: triad_weight(isd:ied,jsd:jed)
  real    :: tensor_11(isd:ied,jsd:jed,0:1,0:1)
  real    :: tensor_13(isd:ied,jsd:jed,0:1,0:1)
  real    :: sumz(isd:ied,jsd:jed,0:1)
  real    :: taperA, taperB
  real    :: slope, absslope, taper_slope, gm_taper_slope
  real    :: depth_ratio
  real    :: drhodx, drhodz 

  tau = Time%tau

  ! initialize arrays to zero 
  tensor_11(:,:,:,:)           = 0.0
  tensor_13(:,:,:,:)           = 0.0
  grav_agm_dz_sx_drhodx(:,:,k) = 0.0
  slopex_drhodx(:,:,k)         = 0.0
  triad_weight(:,:)            = 0.0

  ! tracer-independent part of the calculation 
  do kr=0,1
     kpkr = min(k+kr,nk)
     do ip=0,1 
        do j=jsc,jec
           do i=isc-1,iec

              drhodz   = drhodzh(i+ip,j,k,kr)
              drhodx   = Grd%tmask(i+ip,j,kpkr)                            &
                        *( drhodT(i+ip,j,k)*dTdx(index_temp)%field(i,j,k)  &
                          +drhodS(i+ip,j,k)*dSdx%field(i,j,k))
              slope    = -drhodx/drhodz
              absslope = abs(slope)

              ! taper for steep slope regions 
              if(absslope > smax) then 
                  taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                  taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
              else 
                  taperA = 1.0
              endif

              ! taper for grid depths shallower than eddy_depth
              if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                  depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                  taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
              else 
                  taperB          = 1.0
              endif

              ! taper times slope for use with neutral diffusion
              taper_slope  = taperA*taperB*slope

              ! taper times slope for use with GM skewsion
              if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                  gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
              else 
                  gm_taper_slope = taper_slope
              endif

              ! for diagnosing slope dot grad rho
              triad_weight(i,j)            = triad_weight(i,j)    + Grd%tmask(i+ip,j,kpkr)
              slopex_drhodx(i,j,k)         = slopex_drhodx(i,j,k) + gm_taper_slope*drhodx

              ! fill tensor components 
              tensor_11(i,j,ip,kr)      = aredi_array(i,j,k)*( (1.0-taper_diagonal) + taperA*taperB*taper_diagonal)  &
                                        + ah_array(i,j)*(1.0-taperB)
              tensor_13(i,j,ip,kr)      = aredi_array(i,j,k)*taper_slope-gm_skew*agm_array(i,j,k)*gm_taper_slope

           enddo
        enddo
     enddo
  enddo

  ! tracer-dependent part of the calculation
  do nn=1,num_prog_tracers
     sumz(:,:,:) = 0.0
     do kr=0,1
        do ip=0,1
           do j=jsc,jec
              do i=isc-1,iec
                sumz(i,j,kr) = sumz(i,j,kr) + dtwedyt(i,j,ip)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k) &
                *(tensor_11(i,j,ip,kr)*dTdx(nn)%field(i,j,k) + tensor_13(i,j,ip,kr)*dTdz(nn)%field(i+ip,j,k-1+kr)) &
                * min(delqc(i,j,k,kr),delqc(i+1,j,k,kr))
              enddo
           enddo
        enddo
     enddo
     flux_x(nn)%field(COMPXL,k) = Grd%dxter(COMPXL)*(sumz(COMPXL,0)+sumz(COMPXL,1))
  enddo

  ! apply some masks 
  if(tmask_neutral_on) then
     do nn=1,num_prog_tracers
        do j=jsc,jec
           do i=isc-1,iec
              if(k==Grd%kmt(i,j)) then
                flux_x(nn)%field(i,j,k) = aredi_array(i,j,k)*dTdx(nn)%field(i,j,k)*Grd%dyte(i,j)* &
                                          min(Thickness%rho_dzt(i,j,k,tau),Thickness%rho_dzt(i+1,j,k,tau))
              endif
           enddo
        enddo
     enddo
     if(k==1) then
         do nn=1,num_prog_tracers
            flux_x(nn)%field(:,:,k) = aredi_array(:,:,k)*dTdx(nn)%field(:,:,k)*Grd%dyte(:,:) &
                                      *FMX(Thickness%rho_dzt(:,:,k,tau))
         enddo
     endif
  endif
  if(tmask_sigma_on) then
     do nn=1,num_prog_tracers
        do j=jsc,jec
           do i=isc-1,iec
              if(k==Grd%kmt(i,j)) flux_x(nn)%field(i,j,k) = flux_x(nn)%field(i,j,k)*(1.0-tmask_sigma(i,j))
           enddo
        enddo
     enddo
 endif

 if(neutral_physics_limit) then 
     do nn=1,num_prog_tracers
        do j=jsc,jec
           do i=isc-1,iec
              if(T_prog(nn)%tmask_limit(i,j,k)==1.0) then 
                  flux_x(nn)%field(i,j,k) = aredi_array(i,j,k)*dTdx(nn)%field(i,j,k)*Grd%dyte(i,j) &
                                            *min(Thickness%rho_dzt(i,j,k,tau),Thickness%rho_dzt(i+1,j,k,tau))
              endif
           enddo
        enddo
     enddo
 endif


 ! for diagnosing slopex_drhodx and grav_agm_dz_sx_drhodx 
 do j=jsc,jec
    do i=isc,iec
       if(triad_weight(i,j) > 0.0) then 
           triad_weight(i,j)            = 1.0/(epsln + triad_weight(i,j))
           slopex_drhodx(i,j,k)         = triad_weight(i,j)*slopex_drhodx(i,j,k)
           grav_agm_dz_sx_drhodx(i,j,k) = slopex_drhodx(i,j,k) &
                                          *grav*agm_array(i,j,k)*Thickness%dzt(i,j,k)
       else 
           slopex_drhodx(i,j,k)         = 0.0
           grav_agm_dz_sx_drhodx(i,j,k) = 0.0
       endif
    enddo
 enddo

 if(k==nk) then 
    call diagnose_3d(Time, Grd, id_grav_agm_dz_sx_drhodx, grav_agm_dz_sx_drhodx(:,:,:))
    call diagnose_3d(Time, Grd, id_slopex_drhodx, slopex_drhodx(:,:,:))
 endif

end subroutine fx_flux
! </SUBROUTINE> NAME="fx_flux"


!#######################################################################
! <SUBROUTINE NAME="fx_flux_diag">
!
! <DESCRIPTION>
! Subroutine computes the i-directed neutral physics tracer flux component
! for Redi separately from GM, in order to diagnose GM and Redi 
! fluxes independent of one another. 
!
! Compute this component for all tracers at level k.
!
! fx has physical dimensions (area*diffusivity*density*tracer gradient)
!
! </DESCRIPTION>
!
subroutine fx_flux_diag(Time, Thickness, T_prog, k)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  integer,                      intent(in) :: k

  integer :: nn, i, j, ip, tau
  integer :: kr, kpkr
  real    :: tensor_11(isd:ied,jsd:jed,0:1,0:1)
  real    :: tensor_13_redi(isd:ied,jsd:jed,0:1,0:1)
  real    :: sumz(isd:ied,jsd:jed,0:1)
  real    :: taperA, taperB
  real    :: slope, absslope, taper_slope
  real    :: depth_ratio
  real    :: drhodx, drhodz 

  tau = Time%tau

  ! initialize arrays to zero 
  tensor_13_redi(:,:,:,:) = 0.0


  ! tracer-independent part of the calculation 
  do kr=0,1
     kpkr = min(k+kr,nk)
     do ip=0,1 
        do j=jsc,jec
           do i=isc-1,iec

              drhodz   = drhodzh(i+ip,j,k,kr)
              drhodx   = Grd%tmask(i+ip,j,kpkr)                            &
                        *( drhodT(i+ip,j,k)*dTdx(index_temp)%field(i,j,k)  &
                          +drhodS(i+ip,j,k)*dSdx%field(i,j,k))
              slope    = -drhodx/drhodz
              absslope = abs(slope)

              ! taper for steep slope regions 
              if(absslope > smax) then 
                  taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                  taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
              else 
                  taperA = 1.0
              endif

              ! taper for grid depths shallower than eddy_depth
              if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                  depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                  taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
              else 
                  taperB          = 1.0
              endif

              ! taper times slope for use with neutral diffusion
              taper_slope  = taperA*taperB*slope

              ! fill tensor components 
              tensor_11(i,j,ip,kr)      = aredi_array(i,j,k)*( (1.0-taper_diagonal) + taperA*taperB*taper_diagonal)  &
                                        + ah_array(i,j)*(1.0-taperB)
              tensor_13_redi(i,j,ip,kr) = aredi_array(i,j,k)*taper_slope 

           enddo
        enddo
     enddo
  enddo

  ! tracer dependent portion of calculation 
  do nn=1,num_prog_tracers
     sumz(:,:,:) = 0.0
     do kr=0,1
        do ip=0,1
           do j=jsc,jec
              do i=isc-1,iec
                 sumz(i,j,kr) = sumz(i,j,kr) + dtwedyt(i,j,ip)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k) &
                      *(tensor_11(i,j,ip,kr)*dTdx(nn)%field(i,j,k) + tensor_13_redi(i,j,ip,kr)*dTdz(nn)%field(i+ip,j,k-1+kr)) &
                      * min(delqc(i,j,k,kr),delqc(i+1,j,k,kr))
              enddo
           enddo
        enddo
     enddo
     flux_x_redi(nn)%field(COMPXL,k) = Grd%dxter(COMPXL)*(sumz(COMPXL,0)+sumz(COMPXL,1))
  enddo
 
  ! apply some masks 

  if(tmask_neutral_on) then
      do nn=1,num_prog_tracers
         do j=jsc,jec
            do i=isc-1,iec
               if(k==Grd%kmt(i,j)) then
                   flux_x_redi(nn)%field(i,j,k) = aredi_array(i,j,k)*dTdx(nn)%field(i,j,k)*Grd%dyte(i,j)* &
                        min(Thickness%rho_dzt(i,j,k,tau),Thickness%rho_dzt(i+1,j,k,tau))
               endif
            enddo
         enddo
      enddo
      if(k==1) then
          do nn=1,num_prog_tracers
             flux_x_redi(nn)%field(:,:,k) = aredi_array(:,:,k)*dTdx(nn)%field(:,:,k)*Grd%dyte(:,:) &
                  *FMX(Thickness%rho_dzt(:,:,k,tau))
          enddo
      endif
  endif
  if(tmask_sigma_on) then
      do nn=1,num_prog_tracers
         do j=jsc,jec
            do i=isc-1,iec
               if(k==Grd%kmt(i,j)) then 
                   flux_x_redi(nn)%field(i,j,k) = flux_x_redi(nn)%field(i,j,k)*(1.0-tmask_sigma(i,j))
               endif
            enddo
         enddo
      enddo
  endif
  if(neutral_physics_limit) then 
      do nn=1,num_prog_tracers
         do j=jsc,jec
            do i=isc-1,iec
               if(T_prog(nn)%tmask_limit(i,j,k)==1.0) then 
                   flux_x_redi(nn)%field(i,j,k) = aredi_array(i,j,k)*dTdx(nn)%field(i,j,k)*Grd%dyte(i,j) &
                        *min(Thickness%rho_dzt(i,j,k,tau),Thickness%rho_dzt(i+1,j,k,tau))
               endif
            enddo
         enddo
      enddo
  endif


end subroutine fx_flux_diag
! </SUBROUTINE> NAME="fx_flux_diag"



!#######################################################################
! <SUBROUTINE NAME="fy_flux">
!
! <DESCRIPTION>
! Subroutine computes the j-directed neutral physics tracer flux component.
! Compute this component for all tracers at level k.
!
! fy has physical dimensions (area*diffusivity*density*tracer gradient)
!
! </DESCRIPTION>
!
subroutine fy_flux(Time, Thickness, T_prog, k)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  integer,                      intent(in) :: k

  integer :: nn, i, j, jq, tau
  integer :: kr, kpkr
  real :: triad_weight(isd:ied,jsd:jed)
  real :: tensor_22(isd:ied,jsd:jed,0:1,0:1)
  real :: tensor_23(isd:ied,jsd:jed,0:1,0:1)
  real :: sumz(isd:ied,jsd:jed,0:1)
  real :: taperA, taperB
  real :: slope, absslope, taper_slope, gm_taper_slope
  real :: depth_ratio
  real :: drhody, drhodz 

  tau = Time%tau 

  ! initialize arrays to zero 
  tensor_22(:,:,:,:)           = 0.0
  tensor_23(:,:,:,:)           = 0.0
  grav_agm_dz_sy_drhody(:,:,k) = 0.0
  slopey_drhody(:,:,k)         = 0.0
  triad_weight(:,:)            = 0.0

  ! tracer-independent part of the calculation 
  do kr=0,1
     kpkr = min(k+kr,nk)
     do jq=0,1  
       do j=jsc-1,jec
          do i=isc,iec
 
              drhodz   = drhodzh(i,j+jq,k,kr) 
              drhody   = Grd%tmask(i,j+jq,kpkr)                            &
                        *( drhodT(i,j+jq,k)*dTdy(index_temp)%field(i,j,k)  &
                          +drhodS(i,j+jq,k)*dSdy%field(i,j,k)) 
              slope    = -drhody/drhodz
              absslope = abs(slope)  

              ! taper for steep slope regions 
              if(absslope > smax) then 
                  taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                  taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
              else 
                  taperA = 1.0
              endif

              ! taper for grid depths shallower than eddy_depth
              if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                  depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                  taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
              else 
                  taperB          = 1.0
              endif

              ! taper times slope for use with neutral diffusion
              taper_slope  = taperA*taperB*slope

              ! taper times slope for use with GM skewsion
              if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                  gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
              else 
                  gm_taper_slope = taper_slope
              endif

              ! for diagnosing slope dot grad rho
              triad_weight(i,j)    = triad_weight(i,j)    + Grd%tmask(i,j+jq,kpkr) 
              slopey_drhody(i,j,k) = slopey_drhody(i,j,k) + gm_taper_slope*drhody

              ! fill tensor components 
              tensor_22(i,j,jq,kr)      = aredi_array(i,j,k)*( (1.0-taper_diagonal) + taperA*taperB*taper_diagonal) &
                                        + ah_array(i,j)*(1.0-taperB)
              tensor_23(i,j,jq,kr)      = aredi_array(i,j,k)*taper_slope-gm_skew*agm_array(i,j,k)*gm_taper_slope

           enddo
        enddo
     enddo
  enddo

  ! tracer-dependent part of the calculation
  do nn=1,num_prog_tracers
     sumz(:,:,:) = 0.0
     do kr=0,1
        do jq=0,1
           do j=jsc-1,jec
              do i=isc,iec
                 sumz(i,j,kr) = sumz(i,j,kr) + dxtdtsn(i,j,jq)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k) &
                 *(tensor_22(i,j,jq,kr)*dTdy(nn)%field(i,j,k) + tensor_23(i,j,jq,kr)*dTdz(nn)%field(i,j+jq,k-1+kr)) &
                 * min(delqc(i,j,k,kr),delqc(i,j+1,k,kr))
              enddo
           enddo
        enddo
     enddo
     flux_y(nn)%field(COMPYL,k) = Grd%dytnr(COMPYL)*(sumz(COMPYL,0)+sumz(COMPYL,1))
  enddo


  ! apply some masks 
  if(tmask_neutral_on) then
     do nn=1,num_prog_tracers
        do j=jsc-1,jec
           do i=isc,iec
              if(k==Grd%kmt(i,j)) then
                flux_y(nn)%field(i,j,k) = aredi_array(i,j,k)*dTdy(nn)%field(i,j,k)*Grd%dxtn(i,j) &
                                          *min(Thickness%rho_dzt(i,j,k,tau),Thickness%rho_dzt(i,j+1,k,tau))
              endif
           enddo
        enddo
     enddo
     if(k==1) then
         do nn=1,num_prog_tracers
            flux_y(nn)%field(:,:,k) = &
            aredi_array(:,:,k)*dTdy(nn)%field(:,:,k)*Grd%dxtn(:,:)*FMY(Thickness%rho_dzt(:,:,k,tau))
         enddo
     endif
  endif
  if(tmask_sigma_on) then
      do nn=1,num_prog_tracers
         do j=jsc-1,jec
            do i=isc,iec
               if(k==Grd%kmt(i,j)) flux_y(nn)%field(i,j,k) = flux_y(nn)%field(i,j,k)*(1.0-tmask_sigma(i,j))
            enddo
         enddo
      enddo
  endif

  if(neutral_physics_limit) then 
      do nn=1,num_prog_tracers
         do j=jsc-1,jec
            do i=isc,iec
               if(T_prog(nn)%tmask_limit(i,j,k)==1.0) then 
                   flux_y(nn)%field(i,j,k) = aredi_array(i,j,k)*dTdy(nn)%field(i,j,k)*Grd%dxtn(i,j) &
                                             *min(Thickness%rho_dzt(i,j,k,tau),Thickness%rho_dzt(i,j+1,k,tau))
               endif
            enddo
         enddo
      enddo
  endif

 ! for diagnosing slopey_drhody and grav_agm_dz_sy_drhody
 do j=jsc,jec
    do i=isc,iec
       if(triad_weight(i,j) > 0.0) then 
           triad_weight(i,j)    = 1.0/(epsln + triad_weight(i,j))
           slopey_drhody(i,j,k) = triad_weight(i,j)*slopey_drhody(i,j,k)
           grav_agm_dz_sy_drhody(i,j,k) = slopey_drhody(i,j,k) &
                                          *grav*agm_array(i,j,k)*Thickness%dzt(i,j,k)
       else 
           slopey_drhody(i,j,k)         = 0.0
           grav_agm_dz_sy_drhody(i,j,k) = 0.0
       endif
    enddo
 enddo

 if(k==nk) then 
    call diagnose_3d(Time, Grd, id_grav_agm_dz_sy_drhody, grav_agm_dz_sy_drhody(:,:,:))
    call diagnose_3d(Time, Grd, id_slopey_drhody, slopey_drhody(:,:,:))
     if(id_gm_eddy_ke_source > 0) then 
        call diagnose_3d(Time, Grd, id_gm_eddy_ke_source, grav_agm_dz_sx_drhodx(:,:,:) &
             +grav_agm_dz_sy_drhody(:,:,:))
     endif
 endif

end subroutine fy_flux
! </SUBROUTINE> NAME="fy_flux"



!#######################################################################
! <SUBROUTINE NAME="fy_flux_diag">
!
! <DESCRIPTION>
! Subroutine computes the j-directed neutral physics tracer flux component
! for Redi separately, in order to diagnose GM and Redi contributions
! independent of one another. 
!
! Compute this component for all tracers at level k.
!
! fy has physical dimensions (area*diffusivity*density*tracer gradient)
!
! </DESCRIPTION>
!
subroutine fy_flux_diag(Time, Thickness, T_prog, k)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  integer,                      intent(in) :: k

  integer :: nn, i, j, jq, tau
  integer :: kr, kpkr
  real :: tensor_22(isd:ied,jsd:jed,0:1,0:1)
  real :: tensor_23_redi(isd:ied,jsd:jed,0:1,0:1)
  real :: sumz(isd:ied,jsd:jed,0:1)
  real :: taperA, taperB
  real :: slope, absslope, taper_slope
  real :: depth_ratio
  real :: drhody, drhodz 

  tau = Time%tau 

  ! initialize arrays to zero 
  tensor_22(:,:,:,:)       = 0.0
  tensor_23_redi(:,:,:,:)  = 0.0

  ! tracer-independent part of the calculation 
  do kr=0,1
     kpkr = min(k+kr,nk)
     do jq=0,1  
       do j=jsc-1,jec
          do i=isc,iec
 
              drhodz   = drhodzh(i,j+jq,k,kr) 
              drhody   = Grd%tmask(i,j+jq,kpkr)                            &
                        *( drhodT(i,j+jq,k)*dTdy(index_temp)%field(i,j,k)  &
                          +drhodS(i,j+jq,k)*dSdy%field(i,j,k)) 
              slope    = -drhody/drhodz
              absslope = abs(slope)  

              ! taper for steep slope regions 
              if(absslope > smax) then 
                  taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                  taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
              else 
                  taperA = 1.0
              endif

              ! taper for grid depths shallower than eddy_depth
              if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                  depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                  taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
              else 
                  taperB          = 1.0
              endif

              ! taper times slope for use with neutral diffusion
              taper_slope  = taperA*taperB*slope

              ! fill tensor components 
              tensor_22(i,j,jq,kr)      = aredi_array(i,j,k)*( (1.0-taper_diagonal) + taperA*taperB*taper_diagonal) &
                                        + ah_array(i,j)*(1.0-taperB)
              tensor_23_redi(i,j,jq,kr) = aredi_array(i,j,k)*taper_slope 

           enddo
        enddo
     enddo
  enddo


  ! tracer dependent portion of calculation  
  do nn=1,num_prog_tracers
     sumz(:,:,:) = 0.0
     do kr=0,1
        do jq=0,1
           do j=jsc-1,jec
              do i=isc,iec
                 sumz(i,j,kr) = sumz(i,j,kr) + dxtdtsn(i,j,jq)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)                            &
                      *(tensor_22(i,j,jq,kr)*dTdy(nn)%field(i,j,k) + tensor_23_redi(i,j,jq,kr)*dTdz(nn)%field(i,j+jq,k-1+kr)) &
                      * min(delqc(i,j,k,kr),delqc(i,j+1,k,kr))
              enddo
           enddo
        enddo
     enddo
     flux_y_redi(nn)%field(COMPYL,k) = Grd%dytnr(COMPYL)*(sumz(COMPYL,0)+sumz(COMPYL,1))
  enddo


  ! apply some masks 
  if(tmask_neutral_on) then
      do nn=1,num_prog_tracers
         do j=jsc-1,jec
            do i=isc,iec
               if(k==Grd%kmt(i,j)) then
                   flux_y_redi(nn)%field(i,j,k) = aredi_array(i,j,k)*dTdy(nn)%field(i,j,k)*Grd%dxtn(i,j) &
                        *min(Thickness%rho_dzt(i,j,k,tau),Thickness%rho_dzt(i,j+1,k,tau))
               endif
            enddo
         enddo
      enddo
      if(k==1) then
          do nn=1,num_prog_tracers
             flux_y_redi(nn)%field(:,:,k) = &
                  aredi_array(:,:,k)*dTdy(nn)%field(:,:,k)*Grd%dxtn(:,:)*FMY(Thickness%rho_dzt(:,:,k,tau))
          enddo
      endif
  endif
  if(tmask_sigma_on) then
      do nn=1,num_prog_tracers
         do j=jsc-1,jec
            do i=isc,iec
               if(k==Grd%kmt(i,j)) then 
                   flux_y_redi(nn)%field(i,j,k) = flux_y_redi(nn)%field(i,j,k)*(1.0-tmask_sigma(i,j))
               endif
            enddo
         enddo
      enddo
  endif
  if(neutral_physics_limit) then 
      do nn=1,num_prog_tracers
         do j=jsc-1,jec
            do i=isc,iec
               if(T_prog(nn)%tmask_limit(i,j,k)==1.0) then 
                   flux_y_redi(nn)%field(i,j,k) = aredi_array(i,j,k)*dTdy(nn)%field(i,j,k)*Grd%dxtn(i,j) &
                        *min(Thickness%rho_dzt(i,j,k,tau),Thickness%rho_dzt(i,j+1,k,tau))
               endif
            enddo
         enddo
      enddo
  endif

end subroutine fy_flux_diag
! </SUBROUTINE> NAME="fy_flux_diag"


!#######################################################################
! <SUBROUTINE NAME="fz_flux">
!
! <DESCRIPTION>
! Subroutine computes the vertical neutral physics tracer flux component.
! Compute this component for all tracers at level k.
! Surface and bottom boundary condition fz(k=0)=fz(k=kmt(i,j))=0
!
! fz has physical dimensions (density*diffusivity*tracer gradient)
!
! </DESCRIPTION>
!
subroutine fz_flux(T_prog, k)

  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  integer, intent(in) :: k

  integer :: nn, i, j, ip, jq, kr
  real :: sumx_0, sumx_1, sumy_0, sumy_1
  real :: temparray31(isc:iec,jsc:jec,0:1,0:1)
  real :: temparray32(isc:iec,jsc:jec,0:1,0:1)

  if(tmask_neutral_on .and. k==1) then
     do nn = 1, num_prog_tracers
       fz2(nn)%field(:,:)      = 0.0 
      enddo
     return 
  endif

  if(k==nk) then 

     do nn = 1, num_prog_tracers
        fz2(nn)%field(COMP) = 0.0
     enddo
     return  

  elseif(k < nk) then  

      ! tracer-independent part of the calculation 
      do kr=0,1
         do ip=0,1
            do j=jsc,jec
               do i=isc,iec
                  temparray31(i,j,ip,kr) = tensor_31(i,j,k,ip,kr)*dtew(i,j,ip) &
                       *min(delqc(i-1+ip,j,k+kr,1-kr),delqc(i+ip,j,k+kr,1-kr))
               enddo
            enddo
         enddo
         do jq=0,1
            do j=jsc,jec
               do i=isc,iec
                  temparray32(i,j,jq,kr) = tensor_32(i,j,k,jq,kr)*dtns(i,j,jq) &
                       *min(delqc(i,j-1+jq,k+kr,1-kr),delqc(i,j+jq,k+kr,1-kr))  
               enddo
            enddo
         enddo
      enddo

      ! tracer-dependent part of the calculation  
      do nn=1,num_prog_tracers
         do j=jsc,jec
            do i=isc,iec

               ! compute time explicit portion of the vertical flux
               sumx_0 =  temparray31(i,j,0,0)*dTdx(nn)%field(i-1,j,k) &
                      +  temparray31(i,j,0,1)*dTdx(nn)%field(i-1,j,k+1)
               sumx_1 =  temparray31(i,j,1,0)*dTdx(nn)%field(i,j,k) &
                      +  temparray31(i,j,1,1)*dTdx(nn)%field(i,j,k+1)
               sumy_0 =  temparray32(i,j,0,0)*dTdy(nn)%field(i,j-1,k) &
                      +  temparray32(i,j,0,1)*dTdy(nn)%field(i,j-1,k+1)
               sumy_1 =  temparray32(i,j,1,0)*dTdy(nn)%field(i,j,k) &
                      +  temparray32(i,j,1,1)*dTdy(nn)%field(i,j,k+1)
               fz2(nn)%field(i,j) =   Grd%tmask(i,j,k+1) &
                    *( Grd%dxtr(i,j)*(sumx_0+sumx_1) &
                      +Grd%dytr(i,j)*(sumy_0+sumy_1)) &
                    *dzwtr(i,j,k) &
                    + K33_explicit(i,j,k)*dTdz(nn)%field(i,j,k)

            enddo
         enddo
      enddo

      ! apply some masks       
      if(tmask_neutral_on) then
          do nn=1,num_prog_tracers
             do j=jsc,jec
                do i=isc,iec
                   if(k==(Grd%kmt(i,j)-1))  then 
                       fz2(nn)%field(i,j) = 0.0
                   endif  
                enddo
             enddo
          enddo
      endif
      if(tmask_sigma_on) then
          do nn=1,num_prog_tracers
             do j=jsc,jec
                do i=isc,iec
                   if(tmask_sigma(i,j) > 0.0 .and. k==(Grd%kmt(i,j)-1) ) then
                       fz2(nn)%field(i,j) = fz2(nn)%field(i,j)*(1.0-tmask_sigma(i,j))
                    endif
                enddo
             enddo
          enddo
      endif
      if(neutral_physics_limit) then 
          do nn=1,num_prog_tracers
             do j=jsc,jec
                do i=isc,iec
                   if(T_prog(nn)%tmask_limit(i,j,k)==1.0) then 
                        fz2(nn)%field(i,j) = 0.0  
                   endif 
                enddo
             enddo
          enddo
      endif


  endif  !if-test for k-level 

end subroutine fz_flux
! </SUBROUTINE> NAME="fz_flux"


!#######################################################################
! <SUBROUTINE NAME="fz_flux_diag">
!
! <DESCRIPTION>
! For diagnosing the GM and Redi pieces separately. 
! Compute this component for all tracers at level k.
! Surface and bottom boundary condition fz(k=0)=fz(k=kmt(i,j))=0
!
! fz has physical dimensions (density*diffusivity*tracer gradient)
!
! </DESCRIPTION>
!
subroutine fz_flux_diag(T_prog, k)

  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  integer, intent(in) :: k

  integer :: nn, i, j, ip, jq, kr
  real :: sumx_0, sumx_1, sumy_0, sumy_1

  real :: temparray31_redi(isc:iec,jsc:jec,0:1,0:1)
  real :: temparray32_redi(isc:iec,jsc:jec,0:1,0:1)

 
  if(tmask_neutral_on .and. k==1) then
      do nn=1,num_prog_tracers
         fz2_redi(nn)%field(:,:) = 0.0 
      enddo
      return 
  endif

  if(k==nk) then 

      do nn=1,num_prog_tracers
         fz2_redi(nn)%field(COMP) = 0.0
      enddo
      return  

  elseif(k < nk) then  

      ! tracer-independent part of the calculation 
      do kr=0,1
         do ip=0,1
            do j=jsc,jec
               do i=isc,iec
                  temparray31_redi(i,j,ip,kr) = tensor_31_redi(i,j,k,ip,kr)*dtew(i,j,ip) &
                       *min(delqc(i-1+ip,j,k+kr,1-kr),delqc(i+ip,j,k+kr,1-kr))
               enddo
            enddo
         enddo
         do jq=0,1
            do j=jsc,jec
               do i=isc,iec
                  temparray32_redi(i,j,jq,kr) = tensor_32_redi(i,j,k,jq,kr)*dtns(i,j,jq) &
                       *min(delqc(i,j-1+jq,k+kr,1-kr),delqc(i,j+jq,k+kr,1-kr))  
               enddo
            enddo
         enddo
      enddo

      ! tracer-dependent part of the calculation  
      do nn=1,num_prog_tracers
         do j=jsc,jec
            do i=isc,iec

               ! compute time explicit portion of the vertical flux
               sumx_0 =  temparray31_redi(i,j,0,0)*dTdx(nn)%field(i-1,j,k) &
                      +  temparray31_redi(i,j,0,1)*dTdx(nn)%field(i-1,j,k+1)
               sumx_1 =  temparray31_redi(i,j,1,0)*dTdx(nn)%field(i,j,k) &
                      +  temparray31_redi(i,j,1,1)*dTdx(nn)%field(i,j,k+1)
               sumy_0 =  temparray32_redi(i,j,0,0)*dTdy(nn)%field(i,j-1,k) &
                      +  temparray32_redi(i,j,0,1)*dTdy(nn)%field(i,j-1,k+1)
               sumy_1 =  temparray32_redi(i,j,1,0)*dTdy(nn)%field(i,j,k) &
                      +  temparray32_redi(i,j,1,1)*dTdy(nn)%field(i,j,k+1)
               fz2_redi(nn)%field(i,j) =   Grd%tmask(i,j,k+1) &
                    *( Grd%dxtr(i,j)*(sumx_0+sumx_1)          &
                      +Grd%dytr(i,j)*(sumy_0+sumy_1))         &
                    *dzwtr(i,j,k) &
                    + K33_explicit(i,j,k)*dTdz(nn)%field(i,j,k)

            enddo
         enddo
      enddo


      if(tmask_neutral_on) then
          do nn=1,num_prog_tracers
             do j=jsc,jec
                do i=isc,iec
                   if(k==(Grd%kmt(i,j)-1))  then 
                       fz2_redi(nn)%field(i,j) = 0.0
                   endif
                enddo
             enddo
          enddo
      endif
      if(tmask_sigma_on) then
          do nn=1,num_prog_tracers
             do j=jsc,jec
                do i=isc,iec
                   if(tmask_sigma(i,j) > 0.0 .and. k==(Grd%kmt(i,j)-1) ) then
                       fz2_redi(nn)%field(i,j) = fz2_redi(nn)%field(i,j)*(1.0-tmask_sigma(i,j))
                   endif
                enddo
             enddo
          enddo
      endif
      if(neutral_physics_limit) then 
          do nn=1,num_prog_tracers
             do j=jsc,jec
                do i=isc,iec
                   if(T_prog(nn)%tmask_limit(i,j,k)==1.0) then 
                       fz2_redi(nn)%field(i,j) = 0.0  
                   endif
                enddo
             enddo
          enddo
      endif


  endif  !if-test for k-level 


end subroutine fz_flux_diag
! </SUBROUTINE> NAME="fz_flux_diag"



!#######################################################################
! <SUBROUTINE NAME="gm_velocity">
!
! <DESCRIPTION>
! Subroutine computes GM eddy-induced velocity field for diagnostics.
! Compute ustar and vstar at U-cell point, and wstar at T-cell bottom.   
!
! Do a two-point average rather than more democratic four-point avg
! in order to avoid having to call mpp_update domains on tensor_31 and 
! tensor_32.  The 0.5 factor is due to the two-point average.
!   
! Note that this algorithm is ad hoc.  Researchers interested in this 
! field may wish to test alternatives.  
! </DESCRIPTION>
!
subroutine gm_velocity(Thickness, Time)

  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_time_type), intent(in)      :: Time

  real, dimension(isd:ied,jsd:jed) :: ugm, vgm
  real, dimension(0:1)             :: agmSx, agmSy

  real            :: normalize=0.5
  integer         :: i, j, k, km1, kp1
  integer         :: dtint
  type(time_type) :: next_time

  dtint = nint(dtime)
  next_time =  increment_time(Time%model_time, dtint, 0)

  if(need_data(id_ustar, next_time) .or. &
     need_data(id_vstar, next_time) .or. &
     need_data(id_wstar, next_time)) then 

      wrk1_v(:,:,:,:) = 0.0 
      wrk3(:,:,:)     = 0.0
      ugm(:,:)        = 0.0 
      vgm(:,:)        = 0.0

      do k=1,nk         
         km1 = max(1,k-1) 
         kp1 = min(nk,k+1) 

         do j=jsc,jec
            do i=isc,iec

               agmSx(0) =  agm_array(i,j,k)*(                     &   
                    slope_function_gm(tensor_31(i,j,km1,1,0))   + &
                    slope_function_gm(tensor_31(i,j,km1,1,1)) )   
               agmSx(1) =  agm_array(i,j,k)*(                     &   
                    slope_function_gm(tensor_31(i,j,k,1,0))     + &
                    slope_function_gm(tensor_31(i,j,k,1,1)) )     
               ugm(i,j) = normalize*(agmSx(0)-agmSx(1))/Thickness%dzt(i,j,k)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k) 

               agmSy(0) =  agm_array(i,j,k)*(                     &   
                    slope_function_gm(tensor_32(i,j,km1,1,0))   + &
                    slope_function_gm(tensor_32(i,j,km1,1,1)) )   
               agmSy(1) =  agm_array(i,j,k)*(                     & 
                    slope_function_gm(tensor_32(i,j,k,1,0))     + &
                    slope_function_gm(tensor_32(i,j,k,1,1)) )     
               vgm(i,j) = normalize*(agmSy(0)-agmSy(1))/Thickness%dzt(i,j,k)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k) 

            enddo
         enddo

         call mpp_update_domains (ugm(:,:), Dom%domain2d, flags=EUPDATE) 
         call mpp_update_domains (vgm(:,:), Dom%domain2d, flags=NUPDATE) 
         wrk1_v(:,:,k,1) = FAY(ugm(:,:))
         wrk1_v(:,:,k,2) = FAX(vgm(:,:))
      enddo

      call mpp_update_domains (wrk1_v(:,:,:,1), Dom%domain2d, flags=WUPDATE)
      call mpp_update_domains (wrk1_v(:,:,:,2), Dom%domain2d, flags=SUPDATE)

      k=1
      wrk3(:,:,k) = BDX_ET(wrk1_v(:,:,k,1)*Thickness%dzu(:,:,k)) + &
                    BDY_NT(wrk1_v(:,:,k,2)*Thickness%dzu(:,:,k))
      do k=2,nk
         wrk3(:,:,k) = wrk3(:,:,k-1) + BDX_ET(wrk1_v(:,:,k,1)*Thickness%dzu(:,:,k)) + &
                                       BDY_NT(wrk1_v(:,:,k,2)*Thickness%dzu(:,:,k))
      enddo

      if(need_data(id_ustar, next_time)) then
         call diagnose_3d(Time, Grd, id_ustar, wrk1_v(:,:,:,1))
      endif
      if(need_data(id_vstar, next_time)) then
         call diagnose_3d(Time, Grd, id_vstar, wrk1_v(:,:,:,2))
      endif
      if(need_data(id_wstar, next_time)) then
         call diagnose_3d(Time, Grd, id_wstar, wrk3(:,:,:))
      endif
  else 
      return
  endif


end subroutine gm_velocity
! </SUBROUTINE> NAME="gm_velocity"



!#######################################################################
! <FUNCTION NAME="slope_function_gm">
!
! <DESCRIPTION>
! Function for defining effective slope in diagnostic GM velocity
! calculation. Used only for diagnostic purposes.  
! </DESCRIPTION>
!
function slope_function_gm (slope)

  real, intent(in) :: slope
  real             :: absslope, slope_function_gm

  absslope = abs(slope) 

  if(absslope >= smax) then
      slope_function_gm = smax*sign(1.0,slope)
  else
      slope_function_gm = slope
  endif

end function slope_function_gm
! </FUNCTION> NAME="slope_function_gm"


!#######################################################################
! <SUBROUTINE NAME="nphysics_diagnostics">
! <DESCRIPTION>
!  Send some diagnostics to diagnostics manager. 
! </DESCRIPTION>
subroutine nphysics_diagnostics(Time, T_prog, Dens)

  type(ocean_time_type),        intent(in)  :: Time
  type(ocean_prog_tracer_type), intent(in)  :: T_prog(:)
  type(ocean_density_type),     intent(in)  :: Dens

  integer :: k, nn
  integer :: tau
  real, dimension(isd:ied,jsd:jed) :: tmp_flux

  tau = Time%tau

  if (id_k33_explicit > 0) then 
      call diagnose_3d(Time, Grd, id_k33_explicit, rho0r*K33_explicit(:,:,:))
  endif

  do nn=1,num_prog_tracers
     if(id_neutral_physics(nn) > 0) then 
         call diagnose_3d(Time, Grd, id_neutral_physics(nn), T_prog(nn)%wrk1(:,:,:)*T_prog(nn)%conversion)
     endif

     ! send fluxes to diag_manager 
     ! minus sign is due to a MOM-convention for physics fluxes  
     if(id_flux_x(nn) > 0) then 
        call diagnose_3d(Time, Grd, id_flux_x(nn), -1.0*T_prog(nn)%conversion*flux_x(nn)%field(:,:,:))
     endif
     if(id_flux_y(nn) > 0) then 
        call diagnose_3d(Time, Grd, id_flux_y(nn), -1.0*T_prog(nn)%conversion*flux_y(nn)%field(:,:,:))
     endif
     if(id_flux_x_int_z(nn) > 0) then 
        tmp_flux = 0.0
        do k=1,nk
            tmp_flux(COMP) = tmp_flux(COMP) + flux_x(nn)%field(COMP,k)
        enddo        
        call diagnose_2d(Time, Grd, id_flux_x_int_z(nn), -1.0*T_prog(nn)%conversion*tmp_flux(:,:))
     endif
     if(id_flux_y_int_z(nn) > 0) then 
        tmp_flux = 0.0
        do k=1,nk
            tmp_flux(COMP) = tmp_flux(COMP) + flux_y(nn)%field(COMP,k)
        enddo        
        call diagnose_2d(Time, Grd, id_flux_y_int_z(nn), -1.0*T_prog(nn)%conversion*tmp_flux(:,:))
     endif    
  enddo   ! enddo for nn=1,num_prog_tracers

 
  if(diagnose_gm_redi) then 

      do nn=1,num_prog_tracers

         if(id_neutral_physics_ndiffuse(nn) > 0) then 
            call diagnose_3d(Time, Grd, id_neutral_physics_ndiffuse(nn), &
                 tendency_redi(nn)%field(:,:,:)*T_prog(nn)%conversion)
         endif
         if(id_flux_x_ndiffuse(nn) > 0) then 
            call diagnose_3d(Time, Grd, id_flux_x_ndiffuse(nn), &
                 -1.0*T_prog(nn)%conversion*flux_x_redi(nn)%field(:,:,:))
         endif
         if(id_flux_y_ndiffuse(nn) > 0) then 
            call diagnose_3d(Time, Grd, id_flux_y_ndiffuse(nn), &
                  -1.0*T_prog(nn)%conversion*flux_y_redi(nn)%field(:,:,:))
         endif
         if(id_flux_x_ndiffuse_int_z(nn) > 0) then 
             tmp_flux = 0.0
             do k=1,nk
                tmp_flux(COMP) = tmp_flux(COMP) + flux_x_redi(nn)%field(COMP,k)
             enddo
             call diagnose_2d(Time, Grd, id_flux_x_ndiffuse_int_z(nn), -1.0*T_prog(nn)%conversion*tmp_flux(:,:))
          endif
         if(id_flux_y_ndiffuse_int_z(nn) > 0) then 
             tmp_flux = 0.0
             do k=1,nk
                tmp_flux(COMP) = tmp_flux(COMP) + flux_y_redi(nn)%field(COMP,k)
             enddo
             call diagnose_2d(Time, Grd, id_flux_y_ndiffuse_int_z(nn), -1.0*T_prog(nn)%conversion*tmp_flux(:,:))
         endif


         if(id_neutral_physics_gm(nn) > 0) then 
             wrk1(:,:,:) = 0.0
             do k=1,nk
                wrk1(COMP,k) = T_prog(nn)%wrk1(COMP,k) - tendency_redi(nn)%field(COMP,k)
             enddo
             call diagnose_3d(Time, Grd, id_neutral_physics_gm(nn), &
                  wrk1(:,:,:)*T_prog(nn)%conversion)
         endif

         if(id_flux_x_gm(nn) > 0) then 
             wrk1(:,:,:) = 0.0
             do k=1,nk
                wrk1(COMP,k) = flux_x(nn)%field(COMP,k) - flux_x_redi(nn)%field(COMP,k)
             enddo
             call diagnose_3d(Time, Grd, id_flux_x_gm(nn),          &
                  -1.0*T_prog(nn)%conversion*wrk1(:,:,:))
         endif

         if(id_flux_y_gm(nn) > 0) then 
             wrk1(:,:,:) = 0.0
             do k=1,nk
                wrk1(COMP,k) = flux_y(nn)%field(COMP,k) - flux_y_redi(nn)%field(COMP,k)
             enddo
             call diagnose_2d(Time, Grd, id_flux_y_ndiffuse_int_z(nn), -1.0*T_prog(nn)%conversion*tmp_flux(:,:))
         endif

         if(id_flux_x_gm_int_z(nn) > 0) then 
             wrk1(:,:,:) = 0.0
             do k=1,nk
                wrk1(COMP,k) = flux_x(nn)%field(COMP,k) - flux_x_redi(nn)%field(COMP,k)
             enddo
             tmp_flux = 0.0
             do k=1,nk
                tmp_flux(COMP) = tmp_flux(COMP) + wrk1(COMP,k)
             enddo
             call diagnose_2d(Time, Grd, id_flux_x_gm_int_z(nn), -1.0*T_prog(nn)%conversion*tmp_flux(:,:))
         endif

         if(id_flux_y_gm_int_z(nn) > 0) then 
             wrk1(:,:,:) = 0.0
             do k=1,nk
                wrk1(COMP,k) = flux_y(nn)%field(COMP,k) - flux_y_redi(nn)%field(COMP,k)
             enddo
             tmp_flux = 0.0
             do k=1,nk
                tmp_flux(COMP) = tmp_flux(COMP) + wrk1(COMP,k)
             enddo
             call diagnose_2d(Time, Grd, id_flux_y_gm_int_z(nn), -1.0*T_prog(nn)%conversion*tmp_flux(:,:))
         endif
      enddo   ! enddo for nn=1,num_prog_tracers


      ! more diagnostic output that does not need tracer loop 
      call watermass_diag(Time, T_prog, Dens,  &
      tendency_redi(index_temp)%field(:,:,:), tendency_redi(index_salt)%field(:,:,:))

  endif  ! endif for diagnose_gm_redi logical 



end subroutine nphysics_diagnostics
! </SUBROUTINE> NAME="nphysics_diagnostics"




!#######################################################################
! <SUBROUTINE NAME="ocean_nphysicsA_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_nphysicsA_restart(time_stamp)
  character(len=*), intent(in), optional :: time_stamp

  if(.not. use_this_module) return

  call ocean_nphysics_util_restart(time_stamp)

end subroutine ocean_nphysicsA_restart
! </SUBROUTINE> NAME="ocean_nphysicsA_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_nphysicsA_end">
!
! <DESCRIPTION>
! Write to restart.
! </DESCRIPTION>
!
subroutine ocean_nphysicsA_end(Time)

  type(ocean_time_type), intent(in) :: Time
  
  if(.not. use_this_module) return

  call ocean_nphysics_coeff_end(Time, agm_array, aredi_array, rossby_radius, rossby_radius_raw, bczone_radius)

end subroutine ocean_nphysicsA_end
! </SUBROUTINE> NAME="ocean_nphysicsA_end"


end module ocean_nphysicsA_mod

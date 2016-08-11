module ocean_nphysicsB_mod
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
! plus Laplacian GM skew-diffusion.  The methods here differ from 
! ocean_nphysicsA in the treatment of fluxes in the boundary
! regions. This module is experimental, and should be used with caution. 
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
! R. Ferrari and J.C. McWilliams and Canuto and Dubovikov  
! Parameterization of eddy fluxes near oceanic boundaries
! Journal of Climate (2008). 
! </REFERENCE>
!
! <REFERENCE>
! Large etal (1997), Journal of Physical Oceanography, 
! pages 2418-2447 
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
! <NOTE>
! Revisions made for MOM4p1 in Sept 2005, Jan/Feb 2006,
! and June 2008 by Stephen.Griffies. The June 2008
! revision greatly simplified the boundary layer formulation
! from Ferrari and McWilliams, whereby the quadratic transition
! layer is eliminated, thus removing the need to match vertical
! derivatives of the streamfunction.  The matching conditions
! implied by the transition zone added a tremendous amount 
! of code that was not seen to be critical for the purpose 
! of producing a reasonably smooth streamfunction.  
! </NOTE> 
!
! <NOTE>
! Numerical implementation of the flux components follows the triad 
! approach documented in the references and implemented in MOM2 and MOM3.  
! The MOM algorithm accounts for partial bottom cells and generalized
! orthogonal horizontal coordinates.
! </NOTE> 
!
! <NOTE>
! Note: the option neutral_physics_simple is not supported in this 
! module.  Use nphysicA for that option. 
! </NOTE> 
!
! </INFO>
!
!<NAMELIST NAME="ocean_nphysicsB_nml">
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
!  a neutral flux component.  Default tmask_neutral_on=.false.    
!  </DATA> 
!
!  <DATA NAME="surf_turb_thick_min_k" TYPE="integer">
!  Minimum number of k-levels in surface turbulent boundary
!  layer used in the transition of the neutral physics fluxes
!  to the surface. Default surf_turb_thick_min_k = 2.
!  </DATA> 
!  <DATA NAME="surf_turb_thick_min" TYPE="real">
!  Minimum thickness of surface turbulent boundary layer
!  used in the transition of the neutral physics fluxes
!  to the surface. Default surf_turb_thick_min=20m.
!  </DATA> 
!
!  <DATA NAME="neutral_damping_time" UNITS="days" TYPE="real">
!  The damping time used for determining the effective surface
!  boundary layer thickness from other portions of
!  the model. Default neutral_damping_time=10days.  
!  </DATA> 
!
!  <DATA NAME="nblayer_smooth" TYPE="logical">
!  For smoothing the neutral blayer fields. This is useful
!  when aiming to produce a smooth quasi-stokes streamfunction 
!  within the boundary layers. Default is nblayer_smooth=.true.
!  </DATA> 
!  <DATA NAME="vel_micom_smooth" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM Laplacian mixing 
!  coefficient used in the Laplacian smoothing of neutral blayer fields. 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,           only: epsln, pi
use diag_manager_mod,        only: register_diag_field, register_static_field, need_data
use fms_mod,                 only: FATAL, WARNING, NOTE
use fms_mod,                 only: file_exist
use fms_mod,                 only: open_namelist_file, check_nml_error, close_file, write_version_number
use fms_io_mod,              only: register_restart_field, save_restart, restore_state
use fms_io_mod,              only: restart_file_type
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
use ocean_nphysics_util_mod,     only: compute_eta_tend_gm90
use ocean_nphysics_util_mod,     only: cabbeling_thermob_tendency
use ocean_nphysics_util_mod,     only: watermass_diag_init, watermass_diag 
use ocean_operators_mod,         only: FAX, FAY, FMX, FMY, BDX_ET, BDY_NT, LAP_T
use ocean_parameters_mod,        only: missing_value, onehalf, onefourth, oneeigth, DEPTH_BASED
use ocean_parameters_mod,        only: rho0r, rho0, m3toSv, grav
use ocean_sigma_transport_mod,   only: tmask_sigma_on, tmask_sigma  
use ocean_types_mod,             only: ocean_grid_type, ocean_domain_type, ocean_density_type
use ocean_types_mod,             only: ocean_prog_tracer_type, ocean_thickness_type
use ocean_types_mod,             only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,             only: tracer_2d_type, tracer_3d_0_nk_type, tracer_3d_1_nk_type 
use ocean_util_mod,              only: write_timestamp, write_chksum_2d, write_note
use ocean_util_mod,              only: write_warning, write_line, write_chksum_header
use ocean_util_mod,              only: diagnose_2d, diagnose_3d
use ocean_workspace_mod,         only: wrk1, wrk3, wrk1_v

implicit none

public ocean_nphysicsB_init
public ocean_nphysicsB_end
public nphysicsB
public ocean_nphysicsB_restart

private fx_flux
private fy_flux
private fz_flux
private fz_terms
private fx_flux_diag
private fy_flux_diag
private fz_terms_diag 
private fz_flux_diag
private neutral_blayer
private neutral_chksums
private nphysics_diagnostics

private 

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_domain_type), save    :: Dom_flux

integer :: num_prog_tracers = 0

! clock ids
integer :: id_clock_neutral_blayer
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
integer :: id_tx_trans_gm           =-1
integer :: id_ty_trans_gm           =-1
integer :: id_surf_turb_thick       =-1
integer :: id_surf_trans_thick      =-1
integer :: id_depth_blayer_base     =-1
integer :: id_full_turb_column      =-1
integer :: id_slope_blayer_base     =-1
integer :: id_grav_agm_dz_sx_drhodx =-1
integer :: id_grav_agm_dz_sy_drhody =-1
integer :: id_gm_eddy_ke_source     =-1
integer :: id_slopex_drhodx         =-1
integer :: id_slopey_drhody         =-1

integer, dimension(:), allocatable  :: id_k33_implicit
integer, dimension(:), allocatable  :: id_neutral_physics
integer, dimension(:), allocatable  :: id_neutral_physics_ndiffuse
integer, dimension(:), allocatable  :: id_neutral_physics_gm

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

! for restart
type(restart_file_type), save :: NphysicsB_restart

#include <ocean_memory.h>

#ifdef MOM_STATIC_ARRAYS

real, dimension(isd:ied,jsd:jed,nk,0:1) :: delqc   !density weighted (kg/m^3) quarter cell thickness(m)
real, dimension(isd:ied,jsd:jed,0:nk)   :: dzwtr   !(1/dzwt)(m^-1)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtew    !grid distances from T-point to cell faces (m)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtns    !grid distances from T-point to cell faces (m)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtwedyt !horizontal areas (m^2) of quarter cell
real, dimension(isd:ied,jsd:jed,0:1)    :: dxtdtsn !horizontal areas (m^2) of quarter cell

real, dimension(isd:ied,jsd:jed)    :: smooth_lap  !2D array of micom diffusivities (m^2/sec) for smoothing

real, dimension(isd:ied,jsd:jed)    :: surf_bdlthick ! surface turbulent boundary layer thickness (m)

real, dimension(isd:ied,jsd:jed,nk) :: slopex_drhodx         !3D array of slopex * drhodx for diagnostics
real, dimension(isd:ied,jsd:jed,nk) :: slopey_drhody         !3D array of slopey * drhody for diagnostics
real, dimension(isd:ied,jsd:jed,nk) :: grav_agm_dz_sx_drhodx !3D array of grav*agm*dz*slopex*drhodx for diagnostics
real, dimension(isd:ied,jsd:jed,nk) :: grav_agm_dz_sy_drhody !3D array of grav*agm*dz*slopey*drhody for diagnostics

real, dimension(isd:ied,jsd:jed,nk) :: aredi_array       !3D array of redi diffusivities (m^2/sec)     
real, dimension(isd:ied,jsd:jed,nk) :: agm_array         !3D array of gm diffusivities (m^2/sec)        

real, dimension(isd:ied,jsd:jed)    :: bczone_radius     !bczone calculation (m) 
real, dimension(isd:ied,jsd:jed)    :: rossby_radius     !first baroclinic Rossby radius (m) 
real, dimension(isd:ied,jsd:jed)    :: rossby_radius_raw !first baroclinic Rossby radius (m) without max/min  
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

real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: tensor_31 !portion of mixing tensor
real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: tensor_32 !portion of mixing tensor
real, dimension(isd:ied,jsd:jed,0:1,0:1)    :: tensor_11 !portion of mixing tensor
real, dimension(isd:ied,jsd:jed,0:1,0:1)    :: tensor_13 !portion of mixing tensor
real, dimension(isd:ied,jsd:jed,0:1,0:1)    :: tensor_22 !portion of mixing tensor
real, dimension(isd:ied,jsd:jed,0:1,0:1)    :: tensor_23 !portion of mixing tensor
real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: tensor_31_redi !tracer independent portion of Redi diffusion tensor
real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: tensor_32_redi !tracer independent portion of Redi diffusion tensor

real, dimension(isd:ied,jsd:jed,nk) :: tx_trans_gm       !for diagnosing i-transport due to GM (Sv)
real, dimension(isd:ied,jsd:jed,nk) :: ty_trans_gm       !for diagnosing j-transport due to GM (Sv)

integer, dimension(isd:ied,jsd:jed) :: ksurf_blayer      ! k-value at base of surface nblayer
real, dimension(isd:ied,jsd:jed)    :: surf_turb_thick   ! thickness (m) of surface turbulent boundary layer
real, dimension(isd:ied,jsd:jed)    :: surf_trans_thick  ! thickness (m) of surface transition zone 
real, dimension(isd:ied,jsd:jed)    :: full_turb_column  ! =1 if column fully turbulent, =0 if not. 
real, dimension(isd:ied,jsd:jed)    :: depth_blayer_base ! depth at interface between surf_trans zone and interior 
real, dimension(isd:ied,jsd:jed)    :: slope_blayer_base ! abs(slope) at base of neutral boundary layer 

#else

real, dimension(:,:,:,:),   allocatable :: delqc   !density weighted (kg/m^3) quarter cell thickness(m)
real, dimension(:,:,:),     allocatable :: dzwtr   !(1/dzwt)(m^-1)
real, dimension(:,:,:),     allocatable :: dtew    !grid distances from T-point to cell faces (m)
real, dimension(:,:,:),     allocatable :: dtns    !grid distances from T-point to cell faces (m)
real, dimension(:,:,:),     allocatable :: dtwedyt !horizontal areas (m^2) of quarter cell
real, dimension(:,:,:),     allocatable :: dxtdtsn !horizontal areas (m^2) of quarter cell

real, dimension(:,:),       allocatable :: smooth_lap   !2D array of micom diffusivities (m^2/sec) for smoothing
real, dimension(:,:),       allocatable :: surf_bdlthick !surface turbulent boundary layer thickness (m)

real, dimension(:,:,:),     allocatable :: slopex_drhodx         !3D array of slopex * drhodx for diagnostics
real, dimension(:,:,:),     allocatable :: slopey_drhody         !3D array of slopey * drhody for diagnostics
real, dimension(:,:,:),     allocatable :: grav_agm_dz_sx_drhodx !3D array of grav*agm*dz*slopex*drhodx for diagnostics
real, dimension(:,:,:),     allocatable :: grav_agm_dz_sy_drhody !3D array of grav*agm*dz*slopey*drhody for diagnostics

real, dimension(:,:,:),     allocatable :: aredi_array       !3D array of redi diffusivities (m^2/sec) 
real, dimension(:,:,:),     allocatable :: agm_array         !3D array of gm diffusivities (m^2/sec)     

real, dimension(:,:),       allocatable :: bczone_radius     !bczone calculation (m)
real, dimension(:,:),       allocatable :: rossby_radius     !first baroclinic Rossby radius (m) 
real, dimension(:,:),       allocatable :: rossby_radius_raw !first baroclinic Rossby radius (m) without max/min 
real, dimension(:,:,:),     allocatable :: eady_termx        !rho_z*(S_x)^2 for computing Eady growth rate 
real, dimension(:,:,:),     allocatable :: eady_termy        !rho_z*(S_y)^2 for computing Eady growth rate 
real, dimension(:,:),       allocatable :: baroclinic_termx  !intermediate term for computing vert ave baroclinicity 
real, dimension(:,:),       allocatable :: baroclinic_termy  !intermediate term for computing vert ave baroclinicity 
real, dimension(:,:),       allocatable :: grid_length       !grid length scale (m)

real, dimension(:,:,:),     allocatable :: drhodT           !drho/dtheta     (kg/m^3/C)
real, dimension(:,:,:),     allocatable :: drhodS           !drho/dsalinity  (kg/m^3/psu)
real, dimension(:,:,:,:),   allocatable :: drhodzb          !vertical local ref potrho derivative (kg/m^3/m)
real, dimension(:,:,:,:),   allocatable :: drhodzh          !vertical local ref potrho derivative (kg/m^3/m)

real, dimension(:,:,:,:,:), allocatable :: tensor_31 !portion of mixing tensor
real, dimension(:,:,:,:,:), allocatable :: tensor_32 !portion of mixing tensor
real, dimension(:,:,:,:), allocatable   :: tensor_11 !portion of mixing tensor
real, dimension(:,:,:,:), allocatable   :: tensor_13 !portion of mixing tensor
real, dimension(:,:,:,:), allocatable   :: tensor_22 !portion of mixing tensor
real, dimension(:,:,:,:), allocatable   :: tensor_23 !portion of mixing tensor
real, dimension(:,:,:,:,:), allocatable :: tensor_31_redi !tracer independent portion of Redi diffusion tensor
real, dimension(:,:,:,:,:), allocatable :: tensor_32_redi !tracer independent portion of Redi diffusion tensor

real, dimension(:,:,:),     allocatable :: K33_implicit   !density weighted (kg/m^3) implicit in time 
                                                          !diagonal term in redi tensor (m^2/sec) 
real, dimension(:,:,:),     allocatable :: K33_explicit   !density weighted (kg/m^3) explicit in time 

integer, dimension(:,:), allocatable :: ksurf_blayer      ! k-value at base of surface nblayer 
real, dimension(:,:), allocatable    :: surf_turb_thick   ! thickness (m) of surface turbulent boundary layer
real, dimension(:,:), allocatable    :: surf_trans_thick  ! thickness (m) of surface transition region 
real, dimension(:,:), allocatable    :: full_turb_column  ! =1 if column fully turbulent, =0 if not. 
real, dimension(:,:), allocatable    :: depth_blayer_base ! depth at interface between surf_trans zone and interior 
real, dimension(:,:), allocatable    :: slope_blayer_base ! abs(slope) at base of neutral boundary layer 

real, dimension(:,:,:),  allocatable :: tx_trans_gm       !for diagnosing i-transport due to GM (Sv)
real, dimension(:,:,:),  allocatable :: ty_trans_gm       !for diagnosing j-transport due to GM (Sv)

#endif

! use these derived types so that do not need to know num_prog_tracers at compile time 
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
     '$Id: ocean_nphysicsB.F90,v 20.0 2013/12/14 00:14:38 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

character(len=*), parameter :: FILENAME=&
     __FILE__

logical :: module_is_initialized = .FALSE.

! some useful constants 
real :: gravrho0r
real :: gravrho0r_buoyr
real :: fivedeg
real :: gamma_damp
   
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

! time (days) to damp the evolution of the nonconstant diffusivity,
! as well as for updating the surface turbulent boundary 
! layer thickness that is read into this module. This damping 
! allows for the neutral physics module to use boundary layer thicknesses
! and diffusivities that are evolving on a slow (e.g., many days) time. 
real :: neutral_damping_time=10.0    

! for smoothing the neutral blayer fields
logical :: nblayer_smooth   = .true.
real    :: vel_micom_smooth = 0.2   ! m/sec

! for setting the slope tapering methods 
real    :: smax                ! set in ocean_neutral_util_nml
real    :: swidth              ! set in ocean_neutral_util_nml
logical :: dm_taper=.true.     ! tanh taper scheme of Danabasoglu/McWilliams (only one available here)
real    :: dm_taper_const=1.0  ! internally set to unity when dm_taper=.true.
real    :: swidthr                     ! inverse swidth  
real    :: smax_swidthr                ! useful combination of terms 
logical :: gkw_taper        = .false.  ! quadratic tapering of Gerdes, Koberle, and Willebrand
real    :: gkw_taper_const  = 0.0      ! internally set to unity when gkw_taper=.true.

! for determining the neutral blayer regimes 
integer :: surf_turb_thick_min_k  = 5    ! min klevel in surf_turb_thick_min 
real    :: surf_turb_thick_min    = 50.0 ! metres 

! to compute K33 explicitly in time. 
! diffusion_all_explicit=.false. for realistic simulations.
logical :: diffusion_all_explicit=.false.   

! revert to horizontal diffusion when tracer falls outside specified range 
logical :: neutral_physics_limit=.true.   

! to reduce neutral fluxes to horz/vert diffusion next to model boundaries
logical :: tmask_neutral_on=.false.         

! for specifying transport units
! can either be Sv or mks
character(len=32) :: transport_dims ='Sv (10^9 kg/s)' 
real              :: transport_convert=1.0e-9 

! internally set; for diagnosing the effects from GM and Redi separately 
logical :: diagnose_gm_redi = .false.

! for the module as a whole 
logical :: use_this_module   = .false.
logical :: debug_this_module = .false.

namelist /ocean_nphysicsB_nml/                            &
          use_this_module, debug_this_module,             &
          neutral_physics_limit,  diffusion_all_explicit, &
          tmask_neutral_on, use_gm_skew,                  &
          surf_turb_thick_min, surf_turb_thick_min_k,     &
          nblayer_smooth, vel_micom_smooth,               &
          neutral_damping_time, dm_taper, gkw_taper

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_nphysicsB_init">
!
! <DESCRIPTION>
! Initialize the neutral physics module by registering fields for 
! diagnostic output and performing some numerical checks to see 
! that namelist settings are appropriate. 
! </DESCRIPTION>
!
subroutine ocean_nphysicsB_init(Grid, Domain, Time, Time_steps, Thickness, Dens, T_prog,&
           ver_coordinate_class, agm_closure_lower_dept, agm_closure_upper_dept,        &
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
  integer :: i, j, k, n
  integer :: surf_turb_min_k
  integer :: id_restart
  real    :: ah_array(isd:ied,jsd:jed)
  character(len=64) :: file_name
  
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_nphysicsB_mod (ocean_nphysicsB_init):already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number(version, tagname)

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
  read (input_nml_file, nml=ocean_nphysicsB_nml, iostat=io_status) 
  ierr = check_nml_error(io_status,'ocean_nphysicsB_nml')
#else
  ioun =  open_namelist_file()
  read (ioun,ocean_nphysicsB_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_nphysicsB_nml')
  call close_file (ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_nphysicsB_nml)  
  write (stdlogunit,ocean_nphysicsB_nml)

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk
#endif

  if(use_this_module) then 
     call write_note(FILENAME,&
     'USING ocean_nphysicsB_mod.')
  else 
     call write_note(FILENAME,&
     'NOT using ocean_nphysicsB_mod.')
     return
  endif 

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
     call write_note(FILENAME,&
     'running ocean_nphysicsB_mod with debug_this_module=.true.')
  endif 

  write(stdoutunit,'(/1x,a,f10.2)') &
  '==> Note from ocean_nphysicsB_mod: using forward time step of (secs)', dtime 

  if(cmip_units) then
      transport_convert=1.0
      transport_dims   = 'kg/s'
  else
      transport_convert=1.0e-9 
      transport_dims   = 'Sv (10^9 kg/s)'
  endif

  ! Useful constants 
  two_dtime_inv = 0.5/dtime            !for explicit piece of K33 
  swidthr       = 1.0/(swidth + epsln)                 !for slope taper function with dm_taper
  smax_swidthr  = smax*swidthr                         !for slope taper function with dm_taper
  gamma_damp    = dtime/(86400.0*neutral_damping_time) !for damping blayers and diffusivities 

  if(use_gm_skew) then
    gm_skew=1.0
    call write_note(FILENAME,&
    'use_gm_skew=.true. so will use GM-skewsion.')
  else
    gm_skew=0.0
    call write_note(FILENAME,&
    'use_gm_skew=.false. so will NOT use GM-skewsion.')
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
      '==>Error from ocean_nphysicsB_mod: gkw_taper and dm_taper cannot both be set true--choose only one.')
  endif


  if(neutral_physics_limit) then
     call write_note(FILENAME,&
     'neutral_physics_limit=.true.')
     call write_line('Will revert to horizontal diffusion for points where tracer is outside specified range.')
  endif

  if(dm_taper) then
     call write_note(FILENAME,&
     'dm_taper=.true. Will use the tanh scheme')
     call write_line('of Danabasoglu and McWilliams to taper neutral physics in steep sloped regions.')
     call write_line('This is the only tapering scheme available with ocean_nphysicsB_mod.')
  endif

  if(diffusion_all_explicit) then
     call write_warning(FILENAME,&
     'Running w/ diffusion_all_explicit=.true., which means compute K33 contribution')
     call write_line('to neutral diffusion explicitly in time.  This method is stable only if taking')
     call write_line('very small time steps and/or running with just a single tracer.')
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

  call set_ocean_domain(Dom_flux, Grid, xhalo=Dom%xhalo, yhalo=Dom%yhalo, &
                        name='flux dom neutral',maskmap=Dom%maskmap)

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
  allocate (smooth_lap(isd:ied,jsd:jed))
  allocate (surf_bdlthick(isd:ied,jsd:jed))
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
  allocate (tensor_11(isd:ied,jsd:jed,0:1,0:1))
  allocate (tensor_13(isd:ied,jsd:jed,0:1,0:1))
  allocate (tensor_22(isd:ied,jsd:jed,0:1,0:1))
  allocate (tensor_23(isd:ied,jsd:jed,0:1,0:1))
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

  allocate(ksurf_blayer(isd:ied,jsd:jed))
  allocate(surf_turb_thick(isd:ied,jsd:jed))
  allocate(surf_trans_thick(isd:ied,jsd:jed))
  allocate(depth_blayer_base(isd:ied,jsd:jed))
  allocate(full_turb_column(isd:ied,jsd:jed))
  allocate(slope_blayer_base(isd:ied,jsd:jed))

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

  ! fields related to boundary layer regimes 
  ksurf_blayer(:,:)        = 1
  surf_turb_thick(:,:)     = Grd%tmask(:,:,1)*surf_turb_thick_min
  surf_trans_thick(:,:)    = 0.0
  depth_blayer_base(:,:)   = surf_turb_thick(:,:) + surf_trans_thick(:,:)
  full_turb_column(:,:)    = 0.0
  slope_blayer_base(:,:)   = 0.0

  tx_trans_gm(:,:,:)    = 0.0
  ty_trans_gm(:,:,:)    = 0.0

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

  surf_bdlthick(:,:)           = 0.0
  smooth_lap(:,:)              = vel_micom_smooth*grid_length(:,:)
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
  aredi_array(:,:,:)           = 0.0
  agm_array(:,:,:)             = 0.0
  ah_array(:,:)                = 0.0

  call ocean_nphysics_coeff_init(Time, Thickness, rossby_radius, rossby_radius_raw, &
                           bczone_radius, agm_array, aredi_array, ah_array)

  drhodT(:,:,:)        = 0.0
  drhodS(:,:,:)        = 0.0
  drhodzb(:,:,:,:)     = 0.0
  drhodzh(:,:,:,:)     = 0.0
  tensor_31(:,:,:,:,:) = 0.0
  tensor_32(:,:,:,:,:) = 0.0
  tensor_11(:,:,:,:)   = 0.0
  tensor_13(:,:,:,:)   = 0.0
  tensor_22(:,:,:,:)   = 0.0
  tensor_23(:,:,:,:)   = 0.0
  K33_implicit(:,:,:)  = 0.0           
  K33_explicit(:,:,:)  = 0.0    
  file_name = 'ocean_neutralB.res.nc'     

  ! make mandatory=.false. to facilitate backward compatibility with older restart files. 
  id_restart = register_restart_field(NphysicsB_restart, file_name, "surf_turb_thick", surf_turb_thick, &
       domain=Dom%domain2d, mandatory=.false.)
  id_restart = register_restart_field(NphysicsB_restart, file_name, "surf_trans_thick", surf_trans_thick, &
       domain=Dom%domain2d, mandatory=.false.)
  if(.NOT.file_exist('INPUT/ocean_neutralB.res.nc')) then

      if (.NOT. Time%init) then 
          call mpp_error(FATAL, &
          'Expecting file INPUT/ocean_neutralB.res.nc to exist.&
           &This file was not found and Time%init=.false.')
      endif 

  else
      call restore_state(NphysicsB_restart)
      call mpp_update_domains(surf_turb_thick,Dom%domain2d)
      call mpp_update_domains(surf_trans_thick,Dom%domain2d)

      call write_chksum_header(FILENAME,&
      'initial neutral checksums', Time%model_time)
      call neutral_chksums

  endif

  index_temp=-1;index_salt=-1
  do n=1,num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo
  if (index_temp == -1 .or. index_salt == -1) then 
     call mpp_error(FATAL, &
     '==>Error: temp and/or salt not identified in call to ocean_nphysicsB_init')
  endif 

  ! minimum number of grid cells in surface turbulent boundary layer
  write(stdoutunit,'(1x,a,f10.2)') &
  '==>Note: Surface turb boundary layer for nphysicsB assumed to have min thick (m) =', &
  surf_turb_thick_min 
  kturb_loop : do k=1,nk-1
     surf_turb_min_k = k
     if(Grd%zw(k) <= surf_turb_thick_min .and. Grd%zw(k+1) > surf_turb_thick_min) then
         exit 
     endif
  enddo kturb_loop
  surf_turb_thick_min_k = max(surf_turb_thick_min_k,surf_turb_min_k)
  write(stdoutunit,'(1x,a,i3,a)') &
  '==>Note: Surface turb boundary layer assumed to have at least ',surf_turb_thick_min_k,' k-levels.' 


  ! initialize clock ids 
  id_clock_neutral_blayer = mpp_clock_id('(Ocean neutral: neutral blayer)',grain=CLOCK_ROUTINE)
  id_clock_fz_terms       = mpp_clock_id('(Ocean neutral: fz-terms)'      ,grain=CLOCK_ROUTINE)
  id_clock_fx_flux        = mpp_clock_id('(Ocean neutral: fx-flux)'       ,grain=CLOCK_ROUTINE)
  id_clock_fy_flux        = mpp_clock_id('(Ocean neutral: fy-flux)'       ,grain=CLOCK_ROUTINE)
  id_clock_fz_flux        = mpp_clock_id('(Ocean neutral: fz-flux)'       ,grain=CLOCK_ROUTINE)
  id_clock_fx_flux_diag   = mpp_clock_id('(Ocean neutral: fx-flux-diag)' ,grain=CLOCK_ROUTINE)
  id_clock_fy_flux_diag   = mpp_clock_id('(Ocean neutral: fy-flux-diag)' ,grain=CLOCK_ROUTINE)
  id_clock_fz_flux_diag   = mpp_clock_id('(Ocean neutral: fz-flux-diag)' ,grain=CLOCK_ROUTINE)
  id_clock_fz_terms_diag  = mpp_clock_id('(Ocean neutral: fz-terms-diag)',grain=CLOCK_ROUTINE)


  ! register fields for diagnostic output 

  id_k33_explicit = -1
  id_k33_explicit = register_diag_field ('ocean_model', 'k33_explicit', &
                    Grd%tracer_axes_wt(1:3), Time%model_time,           &
                    'K33_explicit tensor element', 'm^2/sec',           &
                    missing_value=missing_value, range=(/-10.0,1.e20/))

  id_tx_trans_gm = -1
  id_tx_trans_gm = register_diag_field ('ocean_model', 'tx_trans_gm',     &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time,          &
                   'T-cell mass i-transport from GM',trim(transport_dims),&
                   missing_value=missing_value, range=(/-1e20,1.e20/))

  id_ty_trans_gm = -1
  id_ty_trans_gm = register_diag_field ('ocean_model', 'ty_trans_gm',     &
                   Grd%tracer_axes_flux_y(1:3), Time%model_time,          &
                   'T-cell mass j-transport from GM',trim(transport_dims),&
                   missing_value=missing_value, range=(/-1e20,1.e20/))

  id_slope_blayer_base = -1
  id_slope_blayer_base = register_diag_field ('ocean_model','slope_blayer_base',          &
                         Grd%tracer_axes(1:2), Time%model_time,                           &
                         'abs(slope) at neutral physics bdy-layer base', 'dimensionless', &
                         missing_value=missing_value, range=(/-1.0,1.e8/))

  id_surf_turb_thick = -1
  id_surf_turb_thick = register_diag_field ('ocean_model','surf_turb_thick', &
                         Grd%tracer_axes(1:2), Time%model_time,              &
                  'surface turbulent thickness for neutral physics', 'm',    &
                         missing_value=missing_value, range=(/-1.0,1.e8/))

  id_surf_trans_thick = -1
  id_surf_trans_thick = register_diag_field ('ocean_model','surf_trans_thick', &
                  Grd%tracer_axes(1:2), Time%model_time,                       &
                         'Surface transition layer thickness', 'm',            &
                  missing_value=missing_value, range=(/-1.0,1.e8/))

  id_depth_blayer_base = -1
  id_depth_blayer_base = register_diag_field ('ocean_model','depth_blayer_base',          &
                   Grd%tracer_axes(1:2), Time%model_time,                                 &
                         'Depth at bottom of neutral physics surface boundary layer', 'm',&
                   missing_value=missing_value, range=(/-1.0,1.e8/))

  id_full_turb_column = -1
  id_full_turb_column = register_diag_field ('ocean_model','full_turb_column',               &
                        Grd%tracer_axes(1:2), Time%model_time,                               &
                        'regions deemed fully turbulent by neutral physics', 'dimensioness', &
                        missing_value=missing_value, range=(/-1.e8,1.e8/))

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
     call write_note(FILENAME,&
     'ocean_nphysicB has diagnose_gm_redi=.true. to diagnose GM and Redi. It adds memory requirements.')
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



end subroutine ocean_nphysicsB_init
! </SUBROUTINE>  NAME="ocean_nphysicsB_init"


!#######################################################################
! <SUBROUTINE NAME="nphysicsB">
!
! <DESCRIPTION>
! This function computes the thickness weighted and density weighted
! time tendency for tracer from neutral physics.  Full discussion
! and details are provided by Griffies (2004,2005). 
!
! Here is a brief summary of the temporal treatment.  
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
! </DESCRIPTION>
!
subroutine nphysicsB (Time, Thickness, Dens, rho, T_prog, &
                      surf_blthick, gm_diffusivity, baroclinic_rossby_radius)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  real, dimension(isd:,jsd:,:), intent(in)    :: rho
  real, dimension(isd:,jsd:),   intent(in)    :: surf_blthick
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  real, dimension(isd:,jsd:,:), intent(inout) :: gm_diffusivity
  real, dimension(isd:,jsd:),   intent(inout) :: baroclinic_rossby_radius

  integer :: i, j, k, n, nn
  integer :: tau, taum1

  if (.not. use_this_module) return

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_nphysicsB (nphysicsB): needs initialization')
  endif 

  if (size(T_prog(:)) /= num_prog_tracers) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_nphysicsB (nphysicsB): inconsistent size of T_prog')
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

  ! boundary layer thicknesses read from other 
  ! modules are needed here on the data domain
  surf_bdlthick(COMP) = surf_blthick(COMP)
  call mpp_update_domains(surf_bdlthick(:,:), Dom%domain2d) 

  drhodT(:,:,:) = Dens%drhodT(:,:,:)
  drhodS(:,:,:) = Dens%drhodS(:,:,:)

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

  call neutral_blayer(Time, Thickness)

  if(diagnose_gm_redi) then 
      do n=1,num_prog_tracers  
         fz1_redi(n)%field(:,:) = 0.0
         fz2_redi(n)%field(:,:) = 0.0
         flux_x_redi(n)%field(:,:,:)   = 0.0
         flux_y_redi(n)%field(:,:,:)   = 0.0
         tendency_redi(n)%field(:,:,:) = 0.0
      enddo
      call fz_terms_diag(Time, Thickness)   ! must be called prior to fz_terms
  endif

  call fz_terms(Time, Thickness, Dens, T_prog, rho)

  do n=1,num_prog_tracers  
    fz1(n)%field(:,:)      = 0.0
    fz2(n)%field(:,:)      = 0.0
    flux_x(n)%field(:,:,:) = 0.0
    flux_y(n)%field(:,:,:) = 0.0
    T_prog(n)%wrk1(:,:,:)  = 0.0
  enddo 

  
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
         call mpp_update_domains(flux_x(nn)%field(:,:,:), flux_y(nn)%field(:,:,:), Dom_flux%domain2d, &
              gridtype=CGRID_NE, complete=T_prog(nn)%complete) 
      enddo
  endif

  do k=1,nk
     call mpp_clock_begin(id_clock_fz_flux)
     call fz_flux(T_prog,k)
     call mpp_clock_end(id_clock_fz_flux)
     do nn=1,num_prog_tracers
        do j=jsc,jec
           do i=isc,iec
              T_prog(nn)%wrk1(i,j,k) = & 
                   Grd%tmask(i,j,k)  &
                   *(fz1(nn)%field(i,j)-fz2(nn)%field(i,j) &
                   +( flux_x(nn)%field(i,j,k)-flux_x(nn)%field(i-1,j,k) &
                   +flux_y(nn)%field(i,j,k)-flux_y(nn)%field(i,j-1,k) )*Grd%datr(i,j) &
                   )
           enddo
        enddo
        T_prog(nn)%th_tendency(COMP,k) = T_prog(nn)%th_tendency(COMP,k) + T_prog(nn)%wrk1(COMP,k)
        fz1(nn)%field(COMP) = fz2(nn)%field(COMP)
     enddo
  enddo

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
               enddo
            enddo
            fz1_redi(nn)%field(COMP) = fz2_redi(nn)%field(COMP)
         enddo
      enddo  ! enddo for k=1,nk
  endif

  ! compute Eady growth rate for use in next time step 
  call compute_eady_rate(Time, Thickness, T_prog, Dens, eady_termx, eady_termy)

  ! compute baroclinicity for use in next time step 
  call compute_baroclinicity(Time, baroclinic_termx, baroclinic_termy)

  ! compute rossby radius for use in next time step.
  ! rossby radius is needed for diffusivity 
  ! calculation and transition layer calculation.   
  call compute_rossby_radius(Thickness, dTdz, dSdz, Time, drhodT, drhodS, &
                             rossby_radius, rossby_radius_raw)
  baroclinic_rossby_radius(:,:) = rossby_radius_raw(:,:)

  ! compute width of baroclinic zone for use in next time step 
  call compute_bczone_radius(Time, bczone_radius)

  ! update closure-based diffusivity for next time step 
  call compute_diffusivity(Time, ksurf_blayer, Dens%drhodz_zt, rossby_radius, &
                           bczone_radius, agm_array, aredi_array)

  ! gm_diffusivity passed to neutral_physics for computing form drag viscosity 
  gm_diffusivity(:,:,:) = agm_array(:,:,:)

  call nphysics_diagnostics(Time, T_prog, Dens)

end subroutine nphysicsB
! </SUBROUTINE> NAME="nphysicsB"


!#######################################################################
! <SUBROUTINE NAME="neutral_blayer">
!
! <DESCRIPTION>
!
! This subroutine computes the "neutral boundary layers" based on
! the formulation of Ferrari and McWilliams (2006). See full
! details and discussion in Elements of MOM4p1 by Griffies (2009).  
!
! Five vertical regions are identified by Ferrari and McWilliams:
! We simplify these regimes by melding the turbulent and transition 
! regimes into an overall neutral boundary layer regime, within which
! the streamfunction is linearly tapers to zero moving towards the 
! boundary.  We also ignore the bottom regimes, as these are poorly 
! resolved in most models, and the neutral physics fluxes are 
! typically small at the bottom.  
! 
! (1) Surface turbulent region:
!    Depth ("h" in Ferrari and McWilliams notation) dominated by 
!    3d turbulent processes.  This depth is taken from surf_blthick,
!    as set by the KPP scheme or another mixed layer scheme.
!    A minimum is set as surf_turb_thick_min and is specified 
!    as a nml parameter in ocean_nphysicsB_nml. 
!
!    In order to use a low frequency version of the boundary layer
!    thickness, we damp its evolution with a damping time scale
!    neutral_damping_time (days). 
!
!    In the code, "h_surf"= surf_turb_thick  
!
! (2) Surface transition region:
!    Thickness ("D" in Ferrari and McWilliams notation)
!    between the turbulent surface boundary layer and the interior.
!    This transition layer thickness is determined by the product of the 
!    neutral slope and first baroclinic Rossby radius.  This specification 
!    is ad hoc, and more fundamental theories are welcome. 
!
!    In the code, "D_surf"= surf_trans_thick 
!    
!!!!!!
!    Within a "boundary layer" region set by the sum of 
!    surf_turb_thick plus surf_trans_thick, the eddy 
!    induced velocity is assumed to have zero vertical shear, 
!    which means the quasi-Stokes streamfunction is linear with 
!    depth.  The neutral diffusive fluxes are reduced to horizontal 
!    downgradient diffusion, with "horizontal" defined according
!    to surfaces of constant vertical coordinate. 
!!!!!!
!
! (3) Interior region:
!    Where neutral diffusion and GM skew-diffusion are taken 
!    from their unmodified form. 
!
!  Only use the 31 and 32 triads for this computation since the 
!  13 and 23 triads require extra slope calculations, and 
!  so will add lots of computational cost. It is felt that the
!  31 and 32 triads are sufficient for this calculation, in 
!  a similar manner that they are used for the calculation of 
!  the non-constant diffusivities.  
! 
! Scheme coded for MOM4p1 by Stephen.Griffies
! Version: March2006
! Simplified version: June2008
! 
! </DESCRIPTION>
!
subroutine neutral_blayer(Time, Thickness)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_thickness_type), intent(in) :: Thickness

  real, dimension(isd:ied,jsd:jed) :: tmp_ksurf_blayer

  logical :: below_blayer
  logical :: within_interior

  integer :: i, j, k, tau
  integer :: kkp1, kkp2, kkp3
  integer :: kb, kk

  real :: tmp_thickness(4)
  real :: region_thickness(4)
  real :: tmp_surf_turb_thick
  real :: tmp_thick_surf

  call mpp_clock_begin(id_clock_neutral_blayer)

  tau = Time%tau

  ! save the max of abs(slope) in wrk1. this array is
  ! needed to define the surface turb region and transition
  ! region.  
  wrk1(:,:,:) = 0.0
  do k=1,nk
     wrk1(COMP,:) = max( abs(tensor_31(COMP,:,0,0)),abs(tensor_31(COMP,:,0,1)), &
          abs(tensor_31(COMP,:,1,0)),abs(tensor_31(COMP,:,1,1)), & 
          abs(tensor_32(COMP,:,0,0)),abs(tensor_32(COMP,:,0,1)), &
          abs(tensor_32(COMP,:,1,0)),abs(tensor_32(COMP,:,1,1)) )
  enddo

  ! Update thickness of the surface turbulent region for use 
  ! in the neutral_blayer calculation. This thickness is damped relative
  ! to that available from the mixed layer physics, as we are interested
  ! here in the slowly varying form of this thickness. 
  do j=jsd,jed
     do i=isd,ied
        kb=Grd%kmt(i,j)
        if(kb>1) then 
            tmp_surf_turb_thick =  &
              Grd%tmask(i,j,1)*max(surf_turb_thick_min,surf_bdlthick(i,j))
            surf_turb_thick(i,j) = &
              surf_turb_thick(i,j) - gamma_damp*(surf_turb_thick(i,j)-tmp_surf_turb_thick)
                     endif
         enddo
      enddo

  ! spatially smooth surface boundary layer thickness 
  if(nblayer_smooth) then
     call mpp_update_domains(surf_turb_thick(:,:), Dom%domain2d) 
     surf_turb_thick(:,:) = surf_turb_thick(:,:) + dtime*LAP_T(surf_turb_thick(:,:),smooth_lap(:,:))
  endif 


  ! flag for case when full column is turbulent 
  ! (e.g., shallow shelves or deep convection)
  do j=jsc,jec
     do i=isc,iec
        kb=max(1,Grd%kmt(i,j))    
        full_turb_column(i,j) = 0.0
        if(kb>1 .and. surf_turb_thick(i,j) >= Thickness%depth_zwt(i,j,kb)) then 
            full_turb_column(i,j) = 1.0            
        endif
     enddo
  enddo


  ! Determine surface transition region. 
  ! Thise region is determined by the thickness over 
  ! which mesoscale eddies feel the turbulent surface layer.
  ! This thickness is a function of the neutral slope 
  ! times the Rossby radius. The algorithm for computing
  ! this depth is taken from the appendix to Large etal, 
  ! (1997) JPO, pages 2418-2447. It is also described 
  ! in Ferrari and McWilliams (2006) and Griffies (2004).
  ! The thickness of this region is always set > 0. 
  ! 
  ! NOTE: We examine more than one slope*rossby in the vertical
  ! since generally there can be regions of small slopes intermixed
  ! within the large slopes present in the transition layer.
  ! The non-local examination allows for a smoother transition out
  ! of the boundary layer into the interior.  
  !
  ! Also note we take the max from the 8-surrounding triad slopes,
  ! so to maximize the transition layer thicknesses. 
  
  do j=jsc,jec
     do i=isc,iec

        kb             = Grd%kmt(i,j) 
        below_blayer   = .false.
        tmp_thick_surf = 0.0

        ! only compute transition regions if there are enough
        ! vertical cells, and if the column is not fully turbulent. 
        if(kb>surf_turb_thick_min_k+2 .and. full_turb_column(i,j)==0.0) then 


            ! surface region  
            kkloop_surf: do kk=surf_turb_thick_min_k+1,kb-1
               kkp1=min(kk+1,kb)
               kkp2=min(kk+2,kb)
               kkp3=min(kk+3,kb)

               ! to ensure we have captured the boundary layer, 
               ! we check kk, kkp1, kkp2 and kkp3 for good measure. 
               ! this is a bit arbitrary...but useful.  
               region_thickness(1) = rossby_radius(i,j)*wrk1(i,j,kk)
               region_thickness(2) = rossby_radius(i,j)*wrk1(i,j,kkp1)
               region_thickness(3) = rossby_radius(i,j)*wrk1(i,j,kkp2)
               region_thickness(4) = rossby_radius(i,j)*wrk1(i,j,kkp3)
               tmp_thickness(1)    = region_thickness(1) + surf_turb_thick(i,j) 
               tmp_thickness(2)    = region_thickness(2) + surf_turb_thick(i,j) 
               tmp_thickness(3)    = region_thickness(3) + surf_turb_thick(i,j) 
               tmp_thickness(4)    = region_thickness(4) + surf_turb_thick(i,j) 

               if(tmp_thickness(1) <= Thickness%depth_zwt(i,j,kk) .and. &
                    tmp_thickness(2) <= Thickness%depth_zwt(i,j,kk) .and. &
                    tmp_thickness(3) <= Thickness%depth_zwt(i,j,kk) .and. & 
                    tmp_thickness(4) <= Thickness%depth_zwt(i,j,kk) .and. & 
                    .not. below_blayer) then 
                   below_blayer   = .true.
                   tmp_thick_surf = Thickness%depth_zwt(i,j,kk)-surf_turb_thick(i,j)
                   exit kkloop_surf 
               endif
            enddo kkloop_surf

            ! never got below the surface transition layer.
            ! this typically means the slopes are huge or the 
            ! column is very shallow. so column is "fully turbulent".  
            if(.not. below_blayer) then 
                full_turb_column(i,j) = 1.0            
            endif

        endif

        ! column is too small, or is fully turbulent.
        if(kb>1 .and. (full_turb_column(i,j)==1.0 .or. kb<=surf_turb_thick_min_k+2)) then 
            tmp_thick_surf        = 0.0
            full_turb_column(i,j) = 1.0 
        endif

        ! time damping to smooth the field 
        surf_trans_thick(i,j) = surf_trans_thick(i,j) - gamma_damp*(surf_trans_thick(i,j) -tmp_thick_surf) 

     ! end of j,i loops 
     enddo
  enddo


  ! apply spatial smoothing to the boundary field...smoothing may
  ! be needed in space, even with the temporal smoothing used above. 
  if(nblayer_smooth) then
      call mpp_update_domains(surf_trans_thick(:,:), Dom%domain2d) 
      surf_trans_thick(:,:) = surf_trans_thick(:,:) &
           + dtime*LAP_T(surf_trans_thick(:,:),smooth_lap(:,:))
  endif

  ! compute boundary layer depth
  do j=jsc,jec
     do i=isc,iec
        kb=Grd%kmt(i,j)
        if(kb>1) then 
            depth_blayer_base(i,j) = surf_trans_thick(i,j) + surf_turb_thick(i,j)
        endif
     enddo
  enddo
  call mpp_update_domains(depth_blayer_base, Dom%domain2d) 


  ! save the k-level surface transition layer. 
  ! initialize ksurf_blayer to 1 rather than 0, to avoid 
  ! accessing the zeroth element of an array.   
  ksurf_blayer(:,:)     = 1
  tmp_ksurf_blayer(:,:) = 0.0
  do j=jsd,jed
     do i=isd,ied
        kb=Grd%kmt(i,j)
        within_interior=.false.
        ksurf_blayer_loop: do k=surf_turb_thick_min_k,kb
           if(Thickness%depth_zwt(i,j,k) > depth_blayer_base(i,j)) then 
               ksurf_blayer(i,j)     = k
               tmp_ksurf_blayer(i,j) = float(k)
               exit ksurf_blayer_loop
           endif
        enddo ksurf_blayer_loop
     enddo
  enddo

  ! save abs(slope) at base of the surface boundary layer for use in 
  ! computing the streamfunction at this point and then linear 
  ! tapering to the surface. It is important to place a limit 
  ! on slope_blayer_base, since in some regions the blayer_base
  ! could still have infinite slopes...
  do j=jsc,jec
     do i=isc,iec
        k = ksurf_blayer(i,j)
        if(k > 1) then 
            slope_blayer_base(i,j)           = &
                 (abs(tensor_31(i,j,k,0,0))  + & 
                  abs(tensor_31(i,j,k,1,0))  + &
                  abs(tensor_31(i,j,k,0,1))  + & 
                  abs(tensor_31(i,j,k,1,1))  + &
                  abs(tensor_32(i,j,k,0,0))  + &
                  abs(tensor_32(i,j,k,0,1))  + &
                  abs(tensor_32(i,j,k,1,0))  + &
                  abs(tensor_32(i,j,k,1,1)))*oneeigth
            slope_blayer_base(i,j) = min(smax,slope_blayer_base(i,j))
        endif
     enddo
  enddo
  call mpp_update_domains(slope_blayer_base, Dom%domain2d) 


  ! diagnostics 
  call diagnose_2d(Time, Grd, id_slope_blayer_base, slope_blayer_base(:,:))
  call diagnose_2d(Time, Grd, id_surf_turb_thick, surf_turb_thick(:,:))
  call diagnose_2d(Time, Grd, id_surf_trans_thick, surf_trans_thick(:,:))
  call diagnose_2d(Time, Grd, id_depth_blayer_base, depth_blayer_base(:,:))
  call diagnose_2d(Time, Grd, id_full_turb_column, full_turb_column(:,:))

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

                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                     depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
                     taperB      = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB      = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then      
                     gm_taper_slope = (Thickness%depth_zt(i,j,k)/depth_blayer_base(i,j)) &
                                       *slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif


                 tx_trans_gm(i,j,k)     = tx_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope       
                 tensor_31(i,j,k,ip,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumx(ip)               = sumx(ip)                 + mindelqc(ip,kr)*taper_slope2
                 eady_termx(i,j,k)      = eady_termx(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)

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

                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then
                     depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
                     taperB      = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB      = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then      
                     gm_taper_slope = (Thickness%depth_zt(i,j,k)/depth_blayer_base(i,j)) &
                                      *slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 tx_trans_gm(i,j,k)     = tx_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope      
                 tensor_31(i,j,k,ip,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumx(ip)               = sumx(ip)                 + mindelqc(ip,kr)*taper_slope2
                 eady_termx(i,j,k)      = eady_termx(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)
                 sumx(ip)               = sumx(ip)*dtew(i,j,ip)

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

                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                     depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
                     taperB      = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB      = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then      
                     gm_taper_slope = (Thickness%depth_zt(i,j,k)/depth_blayer_base(i,j)) &
                                      *slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 tx_trans_gm(i,j,k)     = tx_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope    
                 tensor_31(i,j,k,ip,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumx(ip)               = sumx(ip)                 + mindelqc(ip,kr)*taper_slope2
                 eady_termx(i,j,k)      = eady_termx(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)

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

                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                     depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
                     taperB      = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB      = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then      
                     gm_taper_slope = (Thickness%depth_zt(i,j,k)/depth_blayer_base(i,j)) &
                                      *slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 tx_trans_gm(i,j,k)     = tx_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope   
                 tensor_31(i,j,k,ip,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumx(ip)               = sumx(ip)                 + mindelqc(ip,kr)*taper_slope2
                 eady_termx(i,j,k)      = eady_termx(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)
                 sumx(ip)               = sumx(ip)*dtew(i,j,ip)

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

                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                     depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
                     taperB      = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB      = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then      
                     gm_taper_slope = (Thickness%depth_zt(i,j,k)/depth_blayer_base(i,j)) &
                                      *slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 ty_trans_gm(i,j,k)     = ty_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope       
                 tensor_32(i,j,k,jq,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumy(jq)               = sumy(jq)                 + mindelqc(jq,kr)*taper_slope2
                 eady_termy(i,j,k)      = eady_termy(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)

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

                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                     depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
                     taperB      = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB      = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then      
                     gm_taper_slope = (Thickness%depth_zt(i,j,k)/depth_blayer_base(i,j)) &
                                      *slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 ty_trans_gm(i,j,k)     = ty_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope      
                 tensor_32(i,j,k,jq,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumy(jq)               = sumy(jq)                 + mindelqc(jq,kr)*taper_slope2
                 eady_termy(i,j,k)      = eady_termy(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)
                 sumy(jq)               = sumy(jq)*dtns(i,j,jq)

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

                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                     depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
                     taperB      = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB      = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then      
                     gm_taper_slope = (Thickness%depth_zt(i,j,k)/depth_blayer_base(i,j)) &
                                      *slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 ty_trans_gm(i,j,k)     = ty_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope     
                 tensor_32(i,j,k,jq,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumy(jq)               = sumy(jq)                 + mindelqc(jq,kr)*taper_slope2
                 eady_termy(i,j,k)      = eady_termy(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)

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

                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                     depth_ratio     = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then      
                     gm_taper_slope = (Thickness%depth_zt(i,j,k)/depth_blayer_base(i,j)) &
                                      *slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 ty_trans_gm(i,j,k)     = ty_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope      
                 tensor_32(i,j,k,jq,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumy(jq)               = sumy(jq)                 + mindelqc(jq,kr)*taper_slope2
                 eady_termy(i,j,k)      = eady_termy(i,j,k)        + absdrhodzb(kr)*absslope*absslope
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)

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
     tx_trans_gm(COMP,:) = tx_trans_gm(COMP,:)*rho0
     ty_trans_gm(COMP,:) = ty_trans_gm(COMP,:)*rho0
  else
     tx_trans_gm(COMP,:) = tx_trans_gm(COMP,:)*rho(COMP,:)
     ty_trans_gm(COMP,:) = ty_trans_gm(COMP,:)*rho(COMP,:)
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
subroutine fz_terms_diag(Time, Thickness)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_thickness_type), intent(in) :: Thickness

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
                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                     depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
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
                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then
                     depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
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
                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                     depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
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
                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                     depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
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
                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                     depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
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
                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                     depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
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
                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                     depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
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
                 ! taper for grid depths shallower than depth_blayer_base
                 if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                     depth_ratio     = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
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

              ! taper for grid depths shallower than depth_blayer_base
              if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                  depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
                  taperB      = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
              else 
                  taperB      = 1.0
              endif

              ! taper times slope for use with neutral diffusion
              taper_slope  = taperA*taperB*slope

              ! taper times slope for use with GM skewsion
              if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then      
                  gm_taper_slope = (Thickness%depth_zt(i,j,k)/depth_blayer_base(i,j)) &
                                   *slope_blayer_base(i,j)*sign(1.0,slope) 
              else 
                  gm_taper_slope = taper_slope
              endif

              ! for diagnosing slope dot grad rho
              triad_weight(i,j)    = triad_weight(i,j)    + Grd%tmask(i+ip,j,kpkr)
              slopex_drhodx(i,j,k) = slopex_drhodx(i,j,k) + gm_taper_slope*drhodx

              ! fill tensor components 
              tensor_11(i,j,ip,kr) = aredi_array(i,j,k)
              tensor_13(i,j,ip,kr) = aredi_array(i,j,k)*taper_slope-gm_skew*agm_array(i,j,k)*gm_taper_slope

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

              ! taper for grid depths shallower than depth_blayer_base
              if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                  depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
                  taperB      = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
              else 
                  taperB      = 1.0
              endif

              ! taper times slope for use with neutral diffusion
              taper_slope  = taperA*taperB*slope

              ! fill tensor components 
              tensor_11(i,j,ip,kr)      = aredi_array(i,j,k)
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

              ! taper for grid depths shallower than depth_blayer_base
              if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                  depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
                  taperB      = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
              else 
                  taperB      = 1.0
              endif

              ! taper times slope for use with neutral diffusion
              taper_slope  = taperA*taperB*slope

              ! taper times slope for use with GM skewsion
              if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then      
                  gm_taper_slope = (Thickness%depth_zt(i,j,k)/depth_blayer_base(i,j)) &
                                   *slope_blayer_base(i,j)*sign(1.0,slope) 
              else 
                  gm_taper_slope = taper_slope
              endif

              ! for diagnosing slope dot grad rho
              triad_weight(i,j)    = triad_weight(i,j)    + Grd%tmask(i,j+jq,kpkr) 
              slopey_drhody(i,j,k) = slopey_drhody(i,j,k) + gm_taper_slope*drhody

              ! fill tensor components 
              tensor_22(i,j,jq,kr) = aredi_array(i,j,k)
              tensor_23(i,j,jq,kr) = aredi_array(i,j,k)*taper_slope-gm_skew*agm_array(i,j,k)*gm_taper_slope

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
              ! taper for grid depths shallower than depth_blayer_base
              if(depth_blayer_base(i,j) >= Thickness%depth_zt(i,j,k)) then                                           
                  depth_ratio = Thickness%depth_zt(i,j,k)/(epsln+depth_blayer_base(i,j))
                  taperB      = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
              else 
                  taperB          = 1.0
              endif

              ! taper times slope for use with neutral diffusion
              taper_slope  = taperA*taperB*slope

              ! fill tensor components 
              tensor_22(i,j,jq,kr)      = aredi_array(i,j,k)
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
            call diagnose_3d(Time, Grd, id_neutral_physics_ndiffuse(nn),&
                  tendency_redi(nn)%field(:,:,:)*T_prog(nn)%conversion)
         endif
         if(id_flux_x_ndiffuse(nn) > 0) then 
            call diagnose_3d(Time, Grd, id_flux_x_ndiffuse(nn),&
                  -1.0*T_prog(nn)%conversion*flux_x_redi(nn)%field(:,:,:))
         endif
         if(id_flux_y_ndiffuse(nn) > 0) then 
            call diagnose_3d(Time, Grd, id_flux_y_ndiffuse(nn),&
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
             wrk1(COMP,:) = T_prog(nn)%wrk1(COMP,:) - tendency_redi(nn)%field(COMP,:)
             call diagnose_3d(Time, Grd, id_neutral_physics_gm(nn), &
                  wrk1(:,:,:)*T_prog(nn)%conversion)
         endif

         if(id_flux_x_gm(nn) > 0) then 
             wrk1(COMP,:) = flux_x(nn)%field(COMP,:) - flux_x_redi(nn)%field(COMP,:)
             call diagnose_3d(Time, Grd, id_flux_x_gm(nn),          &
                  -1.0*T_prog(nn)%conversion*wrk1(:,:,:))
         endif

         if(id_flux_y_gm(nn) > 0) then 
             wrk1(COMP,:) = flux_y(nn)%field(COMP,:) - flux_y_redi(nn)%field(COMP,:)
             call diagnose_3d(Time, Grd, id_flux_y_gm(nn),          &
                  -1.0*T_prog(nn)%conversion*wrk1(:,:,:))
         endif

         if(id_flux_x_gm_int_z(nn) > 0) then 
             wrk1(COMP,:) = flux_x(nn)%field(COMP,:) - flux_x_redi(nn)%field(COMP,:)
             tmp_flux = 0.0
             do k=1,nk
                tmp_flux(COMP) = tmp_flux(COMP) + wrk1(COMP,k)
             enddo
             call diagnose_2d(Time, Grd, id_flux_x_gm_int_z(nn), -1.0*T_prog(nn)%conversion*tmp_flux(:,:))
         endif

         if(id_flux_y_gm_int_z(nn) > 0) then 
             wrk1(COMP,:) = flux_y(nn)%field(COMP,:) - flux_y_redi(nn)%field(COMP,:)
             tmp_flux = 0.0
             do k=1,nk
                tmp_flux(COMP) = tmp_flux(COMP) + wrk1(COMP,k)
             enddo
             call diagnose_2d(Time, Grd, id_flux_y_gm_int_z(nn), -1.0*T_prog(nn)%conversion*tmp_flux(:,:))
         endif
      enddo   ! enddo for nn=1,num_prog_tracers


      ! more diagnostic output that does not need tracer loop 

      if (id_k33_explicit > 0) then 
         call diagnose_3d(Time, Grd, id_k33_explicit, rho0r*K33_explicit(:,:,:))
      endif

      call watermass_diag(Time, T_prog, Dens,  &
      tendency_redi(index_temp)%field(:,:,:), tendency_redi(index_salt)%field(:,:,:))


  endif  ! endif for diagnose_gm_redi logical 


end subroutine nphysics_diagnostics
! </SUBROUTINE> NAME="nphysics_diagnostics"



!#######################################################################
! <SUBROUTINE NAME="neutral_chksums">
!
! <DESCRIPTION>
! Write some checksums.
! </DESCRIPTION>
!
subroutine neutral_chksums

  call write_chksum_2d('surf_turb_thick', surf_turb_thick(COMP)*Grd%tmask(COMP,1))
  call write_chksum_2d('surf_trans_thick', surf_trans_thick(COMP)*Grd%tmask(COMP,1))

end subroutine neutral_chksums
! </SUBROUTINE> NAME="neutral_chksums"


!#######################################################################
! <SUBROUTINE NAME="ocean_nphysicsB_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_nphysicsB_restart(time_stamp)
  character(len=*), intent(in), optional :: time_stamp

  if(.not. use_this_module) return
  call save_restart(NphysicsB_restart, time_stamp)

  call ocean_nphysics_util_restart(time_stamp) 

end subroutine ocean_nphysicsB_restart
! </SUBROUTINE> NAME="ocean_nphysicsB_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_nphysicsB_end">
!
! <DESCRIPTION>
! Write to restart.
! </DESCRIPTION>
!
subroutine ocean_nphysicsB_end(Time)

  type(ocean_time_type), intent(in) :: Time
   
  if(.not. use_this_module) return

  call ocean_nphysics_coeff_end(Time, agm_array, aredi_array, rossby_radius, rossby_radius_raw, bczone_radius)

  call write_chksum_header(FILENAME,&
    'ending neutral checksums', Time%model_time)
  call neutral_chksums

end subroutine ocean_nphysicsB_end
! </SUBROUTINE> NAME="ocean_nphysicsB_end"


end module ocean_nphysicsB_mod

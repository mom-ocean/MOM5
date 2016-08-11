module ocean_overflow_OFP_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> H.-C. Lee 
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
!</REVIEWER>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Z. Liang 
!</REVIEWER>
!
!<OVERVIEW>
! Modeling the physical processes of deep overflow from regional seas. 
!</OVERVIEW>
!
!<DESCRIPTION>
!(1) Physics of the overflow is based on the paper  
!    of Briegleb, Danabasoglu and Large (2010)
!(2) Geographic input data (I and J points) for each region 
!    is given by user at the Field Table.
!(3) Model can compute the overflow processes for 
!    many places. There is no maximum number of overflows. 
!(4) Maximum production lines (or regions) in this code is ten.
!    To get more, code needs to be written.   
!(5) Main acronyms.  
!    src:source, int:interior, ent:entrainment, prd:production
!
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! Briegleb B. P., G. Danabasoglu and W. G. Large (2010), An Overflow Parameterization
!    for the ocean component of the Community Climate System Model. NCAR Technical Note,
!    NCAR Boulder, CO.
! </REFERENCE>
! <REFERENCE>
! Danabasoglu G., W. G. Large and B. P. Briegleb (2010), Climate impacts of parameterized
!    Nordic Sea Overflow, Journal of Geophysical Research, (submitted)  
! </REFERENCE>
! </INFO>
!
!<NAMELIST NAME="ocean_overflow_OFP_nml">
!  
!  <DATA NAME="use_this_module" TYPE="logical">  
!  For using this module.  Default use_this_module=.false.
!  </DATA> 
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!    For debugging.  Default is false.
!  </DATA> 
!
!  <DATA NAME="crit_Fr_geo_ofp" TYPE="real">
!    Critical geostrophic Froude number.
!    Set the minimum Froude number for mixing process
!    between source and entrainment waters
!    Default is 1.0
!  </DATA> 
!
!  <DATA NAME="crit_Fr_geo_ofp" TYPE="real">
!    Maximum overflow speed at the source region
!    Default is 3.0 m/s
!  </DATA> 
!
!  <DATA NAME="frac_exchange_src" TYPE="real">
!    Areal fraction of the overflow exchange at the source region
!    Default is 1.0
!  </DATA> 
!
!  <DATA NAME="max_vol_trans_ofp" TYPE="real">
!    Maximum volume transport of the overflow [m^3/s] 
!    Default is 10.e6
!  </DATA> 
!
!  <DATA NAME="max_ofp_speed" TYPE="real">
!    Maximum overflow speed [m^/s] 
!    Default is 2.0
!  </DATA> 
!
!  <DATA NAME="do_mass_ofp" TYPE="logical">
!    Considering the mass source in the overflow process 
!    Default is .true.
!  </DATA> 
!
!  <DATA NAME="diag_step" TYPE="integer">
!    Diagnostic time step for OFP. 
!    Default is diag_step = -1
!    The diagnostic output is saved in the ascii directory as the ascii format.
!  </DATA> 
!
!  <DATA NAME="do_entrainment_para_ofp" TYPE="logical">
!    Considering the parameterization of entrainment process 
!    in the overflow process 
!    Default is .true.
!  </DATA> 
!
!</NAMELIST>
!

use constants_mod,       only: epsln, deg_to_rad, omega
use diag_manager_mod,    only: register_diag_field, register_static_field, send_data
use field_manager_mod,   only: MODEL_OCEAN, parse, find_field_index
use field_manager_mod,   only: get_field_methods, method_type, get_field_info
use fms_mod,             only: write_version_number, error_mesg, FATAL, NOTE
use fms_mod,             only: open_namelist_file, check_nml_error, close_file, stdout, stdlog
use mpp_domains_mod,     only: mpp_update_domains, CGRID_NE
use mpp_domains_mod,     only: mpp_global_sum, BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
use mpp_domains_mod,     only: mpp_define_domains
use mpp_domains_mod,     only: cyclic_global_domain, global_data_domain 
use mpp_mod,             only: mpp_error, mpp_max, mpp_sum, mpp_pe, mpp_root_pe, mpp_broadcast
use mpp_mod,             only: input_nml_file
use mpp_mod,             only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod,             only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE, CLOCK_ROUTINE

use ocean_density_mod,    only: density
use ocean_domains_mod,    only: get_local_indices, set_ocean_domain
use ocean_parameters_mod, only: missing_value, rho0r, grav, omega_earth
use ocean_parameters_mod, only: TERRAIN_FOLLOWING
use ocean_types_mod,      only: ocean_domain_type, ocean_grid_type, ocean_time_type, ocean_options_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_density_type, ocean_thickness_type
use ocean_types_mod,      only: ocean_velocity_type, ocean_external_mode_type
use ocean_util_mod,       only: write_timestamp, diagnose_3d
use ocean_tracer_util_mod,only: diagnose_3d_rho
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk4, wrk1_v

implicit none

private 

public ocean_overflow_OFP_init
public overflow_OFP
private watermass_diag_init
private watermass_diag

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

logical :: used

integer :: id_OFP_n1_src_temp     =-1
integer :: id_OFP_n1_src_salt     =-1
integer :: id_OFP_n1_src_trans    =-1
integer :: id_OFP_n1_int_temp     =-1
integer :: id_OFP_n1_int_salt     =-1
integer :: id_OFP_n1_ent_temp     =-1
integer :: id_OFP_n1_ent_salt     =-1
integer :: id_OFP_n1_ent_trans    =-1
integer :: id_OFP_n1_prd_temp     =-1
integer :: id_OFP_n1_prd_salt     =-1
integer :: id_OFP_n1_prd_trans    =-1
integer :: id_OFP_n1_prd_depth    =-1

integer :: id_OFP_n2_src_temp     =-1
integer :: id_OFP_n2_src_salt     =-1
integer :: id_OFP_n2_src_trans    =-1
integer :: id_OFP_n2_int_temp     =-1
integer :: id_OFP_n2_int_salt     =-1
integer :: id_OFP_n2_ent_temp     =-1
integer :: id_OFP_n2_ent_salt     =-1
integer :: id_OFP_n2_ent_trans    =-1
integer :: id_OFP_n2_prd_temp     =-1
integer :: id_OFP_n2_prd_salt     =-1
integer :: id_OFP_n2_prd_trans    =-1
integer :: id_OFP_n2_prd_depth    =-1

integer :: id_OFP_init
integer :: id_OFP_main
integer :: id_OFP_update_1
integer :: id_OFP_update_2
integer :: id_OFP_update_3

integer :: id_neut_rho_overofp          =-1
integer :: id_wdian_rho_overofp         =-1
integer :: id_neut_rho_overofp_on_nrho  =-1
integer :: id_wdian_rho_overofp_on_nrho =-1
integer :: id_tform_rho_overofp_on_nrho =-1

integer :: id_neut_temp_overofp          =-1
integer :: id_wdian_temp_overofp         =-1
integer :: id_neut_temp_overofp_on_nrho  =-1
integer :: id_wdian_temp_overofp_on_nrho =-1
integer :: id_tform_temp_overofp_on_nrho =-1

integer :: id_neut_salt_overofp          =-1
integer :: id_wdian_salt_overofp         =-1
integer :: id_neut_salt_overofp_on_nrho  =-1
integer :: id_wdian_salt_overofp_on_nrho =-1
integer :: id_tform_salt_overofp_on_nrho =-1


#include <ocean_memory.h>

!------- for multi-regions, input from field table
 real, dimension(:), allocatable :: phi_ofp             !  latitude of overflow (degree N)
 real, dimension(:), allocatable :: hu_ofp              !  upstream thickness of the source water (m)
 real, dimension(:), allocatable :: hs_ofp              
 real, dimension(:), allocatable :: ws_ofp              !  width of strait (m)
 real, dimension(:), allocatable :: xssb_ofp            !  distance from strait to shelf-slop break (m)
 real, dimension(:), allocatable :: alpha_ofp           !  max. bottom slope near shelf break
 real, dimension(:), allocatable :: cdbot_ofp           !  bottom shelf drag coeff.

 integer                         :: OFP_xhalo           !  halo index in x-direction
 integer                         :: OFP_yhalo           !  halo index in y-direction
 type(ocean_domain_type), save   :: OFP_domain          !  domain for OFP

!----- n_OFP:num_prog_tracers
 real, dimension(:,:), allocatable   :: src_flux_trace    ! total source tracer flux [trace kg/m^3 m^3/s]
 real, dimension(:,:), allocatable   :: ent_flux_trace    ! total entrainment tracer flux [trace kg/m^3 m^3/s]
 real, dimension(:,:), allocatable   :: prd_flux_trace    ! total product tracer flux [trace kg/m^3 m^3/s]

!----- n_OFP:num_prog_tracers
 real, dimension(:,:),   allocatable   :: src_vol_trace   ! all tracers in volume average
 real, dimension(:,:),   allocatable   :: ent_vol_trace   ! all tracers in volume average
 real, dimension(:,:),   allocatable   :: prd_vol_trace   ! all tracers in volume average

!----- n_OFP  *** should allocating and localizing
 real, dimension(:),   allocatable   :: src_volT   
 real, dimension(:),   allocatable   :: src_volS   
 real, dimension(:),   allocatable   :: src_volP   
 real, dimension(:),   allocatable   :: src_vol    
 real, dimension(:),   allocatable   :: src_dat   
 real, dimension(:),   allocatable   :: src_flux_Mass   
 real, dimension(:),   allocatable   :: src_den   

 real, dimension(:),   allocatable   :: int_volT   
 real, dimension(:),   allocatable   :: int_volS   
 real, dimension(:),   allocatable   :: int_volP   
 real, dimension(:),   allocatable   :: int_vol   
 real, dimension(:),   allocatable   :: int_dat   
 real, dimension(:),   allocatable   :: int_den   

 real, dimension(:),   allocatable   :: grav_prim_src_ofp  
 real, dimension(:),   allocatable   :: coriolis_source_ofp  
 real, dimension(:),   allocatable   :: trans_vol_ofp_src_raw  
 real, dimension(:),   allocatable   :: trans_vol_ofp_src  
 real, dimension(:),   allocatable   :: vel_src_ofp  

 real, dimension(:),   allocatable   :: ent_volT   
 real, dimension(:),   allocatable   :: ent_volS   
 real, dimension(:),   allocatable   :: ent_volP   
 real, dimension(:),   allocatable   :: ent_vol    
 real, dimension(:),   allocatable   :: ent_dat   
 real, dimension(:),   allocatable   :: ent_den   

 real, dimension(:),   allocatable   :: src_den_prim  
 real, dimension(:),   allocatable   :: src_den_prim_int  
 real, dimension(:),   allocatable   :: grav_prim_ent_ofp 
 real, dimension(:),   allocatable   :: trans_vol_ofp_ent 
 real, dimension(:),   allocatable   :: trans_vol_ofp_prd 
 real, dimension(:),   allocatable   :: trans_mass_ofp_ent 
 real, dimension(:),   allocatable   :: vel_ssb_ofp 
 real, dimension(:),   allocatable   :: vel_ave_ofp 
 real, dimension(:),   allocatable   :: coef_b_ofp 
 real, dimension(:),   allocatable   :: coef_c_ofp 
 real, dimension(:),   allocatable   :: coef_d_ofp 
 real, dimension(:),   allocatable   :: hssb_ofp
 real, dimension(:),   allocatable   :: Fr_geo_ofp 
 real, dimension(:),   allocatable   :: mix_theta_ofp 
 real, dimension(:),   allocatable   :: volf_per_tvol_src 
 real, dimension(:),   allocatable   :: volf_per_tvol_ent
 real, dimension(:),   allocatable   :: ent_flux_Mass

 integer, dimension(:), allocatable  :: prd_line_num,prd_inject_line_num
 real, dimension(:,:),   allocatable :: prd_amb_P, prd_dat, prd_amb_rho
 real, dimension(:,:),   allocatable :: prd_vol, prd_amb_T, prd_amb_S
 real, dimension(:),   allocatable   :: prd_volT,prd_volS,prd_vol_rho
 real, dimension(:),   allocatable   :: prd_flux_Mass, prd_inj_vol
 real, dimension(:),   allocatable   :: volf_per_tvol_prd
 real, dimension(:),   allocatable   :: depth_prd_OFP


!----- pre-defined informations (I,J points) of src, int and ent
integer, dimension(:), allocatable :: src_ist,src_ied,src_jst,src_jed,src_kst,src_ked
integer, dimension(:), allocatable :: int_ist,int_ied,int_jst,int_jed,int_kst,int_ked
integer, dimension(:), allocatable :: ent_ist,ent_ied,ent_jst,ent_jed,ent_kst,ent_ked

!----- I,J points of prd lines (regions)
integer, dimension(:,:), allocatable :: prd_ist,prd_ied,prd_jst,prd_jed,prd_kst,prd_ked

!----- computational domains for regions
integer, dimension(:),   allocatable :: src_isc,src_iec,src_jsc,src_jec
integer, dimension(:),   allocatable :: ent_isc,ent_iec,ent_jsc,ent_jec
integer, dimension(:),   allocatable :: int_isc,int_iec,int_jsc,int_jec
integer, dimension(:,:), allocatable :: prd_isc,prd_iec,prd_jsc,prd_jec

!----- maximum I,J bounds
integer, dimension(:),   allocatable :: max_i_bound,min_i_bound,max_j_bound,min_j_bound
integer, dimension(:),   allocatable :: max_ist,max_ied,max_jst,max_jed

!----- checking the working computational domain for each regions
logical                              :: working_comp_domains_all
logical, dimension(:), allocatable   :: working_comp_domains
logical, dimension(:), allocatable   :: working_src_comp_domains
logical, dimension(:), allocatable   :: working_int_comp_domains
logical, dimension(:), allocatable   :: working_ent_comp_domains
logical, dimension(:,:), allocatable :: working_prd_comp_domains
logical, dimension(:), allocatable   :: occur_overflow_OFP



!------ variables in new ocean domain with halo region
integer, dimension(:,:)    , allocatable :: kmto      ! Grd%kmt(i,j)
real,    dimension(:,:)    , allocatable :: dato      ! Grd%dat(i,j)
real,    dimension(:,:)    , allocatable :: dxto      ! Grd%dxt(i,j)
real,    dimension(:,:)    , allocatable :: dyto      ! Grd%dyt(i,j)
real,    dimension(:,:)    , allocatable :: hto       ! Grd%dyt(i,j)
real,    dimension(:,:,:)  , allocatable :: tmsko     ! Grd%tmask(i,j,k)
real,    dimension(:,:,:)  , allocatable :: umsko     ! Grd%umask(i,j,k)
real,    dimension(:,:,:)  , allocatable :: dzto      ! Thickness%dzt(i,j,k) 
real,    dimension(:,:,:)  , allocatable :: depth_zto ! Thickness%depth_zt(i,j,k) 
real,    dimension(:,:,:)  , allocatable :: rho_dzto  ! Thickness%rho_dzt(i,j,k,tau)
real,    dimension(:,:,:)  , allocatable :: tempo     ! T_prog(index_temp)%field(i,j,k,tau) 
real,    dimension(:,:,:)  , allocatable :: salto     ! T_prog(index_salt)%field(i,j,k,tau) 
real,    dimension(:,:,:,:), allocatable :: traceo    ! T_prog(n)%field(i,j,k,tau)
real,    dimension(:,:,:)  , allocatable :: presso    ! Dens%pressure_at_depth(i,j,k) 

!------ number of OFP regions
integer :: n_OFP_info=0 
integer :: n_OFP_src=0 
integer :: n_OFP_int=0 
integer :: n_OFP_ent=0 
integer :: n_OFP_prd_line_01=0 
integer :: n_OFP_prd_line_02=0 
integer :: n_OFP_prd_line_03=0 
integer :: n_OFP_prd_line_04=0 
integer :: n_OFP_prd_line_05=0 
integer :: n_OFP_prd_line_06=0 
integer :: n_OFP_prd_line_07=0 
integer :: n_OFP_prd_line_08=0 
integer :: n_OFP_prd_line_09=0 
integer :: n_OFP_prd_line_10=0 
integer :: n_OFP=0 
!-integer :: prd_line_num=10


! internally set for computing watermass diagnostics
logical :: compute_watermass_diag = .false. 

!----------------------------------------------------------------------------------

character(len=128) :: version=&
       '=>Using: ocean_overflow_OFP.f90 ($Id: ocean_overflow_OFP.F90,v 20.0 2013/12/14 00:16:08 fms Exp $)'
character (len=128) :: tagname=&
     '$Name: tikal $'

! number of prognostic tracers
integer :: num_prog_tracers=0

! initialization flag 
logical :: module_is_initialized=.false.

! flag for mpp_global_sum
integer :: global_sum_flag      

! time step 
real :: dtime 
real :: p5dtimer

! set from nml
logical :: use_this_module         = .false.    ! must be set .true. in nml to enable this scheme
logical :: debug_this_module       = .false.    ! for debugging
real    :: crit_Fr_geo_ofp         = 1.0        ! Critical geostrophic Froude number
real    :: frac_exchange_src       = 1.0        ! Areal fraction of the exchange at the source region
real    :: max_vol_trans_ofp       = 10.e6      ! Maximum volume transport of the overflow
real    :: max_ofp_speed           = 2.0        ! Maximum OFP speed of the overflow
logical :: do_mass_ofp             = .true.     ! considering the mass source in the overflow process
integer :: diag_step               = -1         ! diagnostic time step
logical :: do_entrainment_para_ofp = .true.     ! parameterization of entrainment process


namelist /ocean_overflow_OFP_nml/ use_this_module, debug_this_module, crit_Fr_geo_ofp &
                                 ,crit_Fr_geo_ofp, frac_exchange_src,max_vol_trans_ofp &
                                 ,max_ofp_speed, do_mass_ofp, diag_step, do_entrainment_para_ofp

contains



!#######################################################################
! <SUBROUTINE NAME="ocean_overflow_OFP_init">
!
! <DESCRIPTION>
! Initial set up for OFP
! </DESCRIPTION>
!
  subroutine ocean_overflow_OFP_init(Grid, Domain, Time, Dens, T_prog, Ocean_options, &
                                 vert_coordinate_type, dtim, debug)

    type(ocean_grid_type),        intent(in), target   :: Grid
    type(ocean_domain_type),      intent(in), target   :: Domain
    type(ocean_time_type),        intent(in), target   :: Time
    type(ocean_density_type),     intent(in)           :: Dens
    type(ocean_prog_tracer_type), intent(in)           :: T_prog(:)
    type(ocean_options_type),     intent(inout)        :: Ocean_options
    integer,                      intent(in)           :: vert_coordinate_type 
    real,                         intent(in)           :: dtim
    logical,                      intent(in), optional :: debug

    integer :: parse_ok
    integer :: io_status, ioun, ierr
    integer :: i, j, k, m, model, nindx, mindx

!------- index of fields
    integer :: n_info, n_src, n_int, n_ent
    integer :: n_prd_line_01, n_prd_line_02, n_prd_line_03, n_prd_line_04, n_prd_line_05
    integer :: n_prd_line_06, n_prd_line_07, n_prd_line_08, n_prd_line_09, n_prd_line_10

!------- methods of fields
    type(method_type), allocatable, dimension(:) :: OFP_methods_info
    type(method_type), allocatable, dimension(:) :: OFP_methods_src, OFP_methods_int, OFP_methods_ent
    type(method_type), allocatable, dimension(:) :: OFP_methods_prd_line_01
    type(method_type), allocatable, dimension(:) :: OFP_methods_prd_line_02
    type(method_type), allocatable, dimension(:) :: OFP_methods_prd_line_03
    type(method_type), allocatable, dimension(:) :: OFP_methods_prd_line_04
    type(method_type), allocatable, dimension(:) :: OFP_methods_prd_line_05
    type(method_type), allocatable, dimension(:) :: OFP_methods_prd_line_06
    type(method_type), allocatable, dimension(:) :: OFP_methods_prd_line_07
    type(method_type), allocatable, dimension(:) :: OFP_methods_prd_line_08
    type(method_type), allocatable, dimension(:) :: OFP_methods_prd_line_09
    type(method_type), allocatable, dimension(:) :: OFP_methods_prd_line_10

    character(len=32) :: fld_type, fld_name
    integer :: stdoutunit, stdlogunit

!----------- pe number
    integer :: pe
    integer :: root_pe

!------------------------
        id_OFP_init        = mpp_clock_id('(Ocean_OFP_init) '           ,grain=CLOCK_MODULE)
        id_OFP_main        = mpp_clock_id('(Ocean_OFP_main) '           ,grain=CLOCK_MODULE)
        id_OFP_update_1    = mpp_clock_id('(ocean_OFP_update_1) '       ,grain=CLOCK_MODULE)
        id_OFP_update_2    = mpp_clock_id('(ocean_OFP_update_2) '       ,grain=CLOCK_MODULE)
        id_OFP_update_3    = mpp_clock_id('(ocean_OFP_update_3) '       ,grain=CLOCK_MODULE)
!------------------------
!    OFP diagnostics
!------------------------
id_OFP_n1_src_temp     = register_diag_field ('ocean_model', 'OFP_n1_src_temp',          &
                         Time%model_time, 'OFP_1 source region temperature', 'oC degree',&
                         missing_value=missing_value, range=(/-100.,100./))
id_OFP_n1_src_salt     = register_diag_field ('ocean_model', 'OFP_n1_src_salt', &
                         Time%model_time, 'OFP_1 source region salinity', 'psu',&
                         missing_value=missing_value, range=(/0.,100./))
id_OFP_n1_src_trans    = register_diag_field ('ocean_model', 'OFP_n1_src_trans',   &
                         Time%model_time, 'OFP_1 source region transport', 'm^3/s',&
                         missing_value=missing_value, range=(/0.,100./))

id_OFP_n1_int_temp     = register_diag_field ('ocean_model', 'OFP_n1_int_temp',            &
                         Time%model_time, 'OFP_1 interior region temperature', 'oC degree',&
                         missing_value=missing_value, range=(/-100.,100./))
id_OFP_n1_int_salt     = register_diag_field ('ocean_model', 'OFP_n1_int_salt',   &
                         Time%model_time, 'OFP_1 interior region salinity', 'psu',&
                         missing_value=missing_value, range=(/0.,100./))

id_OFP_n1_ent_temp     = register_diag_field ('ocean_model', 'OFP_n1_ent_temp',               &
                         Time%model_time, 'OFP_1 entrainment region temperature', 'oC degree',&
                         missing_value=missing_value, range=(/-100.,100./))
id_OFP_n1_ent_salt     = register_diag_field ('ocean_model', 'OFP_n1_ent_salt',       &
                         Time%model_time, 'OFP_1 entrainment  region salinity', 'psu',&
                         missing_value=missing_value, range=(/0.,100./))
id_OFP_n1_ent_trans    = register_diag_field ('ocean_model', 'OFP_n1_ent_trans',        &
                         Time%model_time, 'OFP_1 entrainment region transport', 'm^3/s',&
                         missing_value=missing_value, range=(/0.,100./))

id_OFP_n1_prd_temp     = register_diag_field ('ocean_model', 'OFP_n1_prd_temp',              &
                         Time%model_time, 'OFP_1 production region temperature', 'oC degree',&
                         missing_value=missing_value, range=(/-100.,100./))
id_OFP_n1_prd_salt     = register_diag_field ('ocean_model', 'OFP_n1_prd_salt',      &
                         Time%model_time, 'OFP_1 production  region salinity', 'psu',&
                         missing_value=missing_value, range=(/0.,100./))
id_OFP_n1_prd_trans    = register_diag_field ('ocean_model', 'OFP_n1_prd_trans',       &
                         Time%model_time, 'OFP_1 production region transport', 'm^3/s',&
                         missing_value=missing_value, range=(/0.,100./))
id_OFP_n1_prd_depth    = register_diag_field ('ocean_model', 'OFP_n1_prd_depth',&
                         Time%model_time, 'OFP_1 production depth', 'meter',    &
                         missing_value=missing_value, range=(/0.,10000./))

id_OFP_n2_src_temp     = register_diag_field ('ocean_model', 'OFP_n2_src_temp',          &
                         Time%model_time, 'OFP_2 source region temperature', 'oC degree',&
                         missing_value=missing_value, range=(/-100.,100./))
id_OFP_n2_src_salt     = register_diag_field ('ocean_model', 'OFP_n2_src_salt', &
                         Time%model_time, 'OFP_2 source region salinity', 'psu',&
                         missing_value=missing_value, range=(/0.,100./))
id_OFP_n2_src_trans    = register_diag_field ('ocean_model', 'OFP_n2_src_trans',   &
                         Time%model_time, 'OFP_2 source region transport', 'm^3/s',&
                         missing_value=missing_value, range=(/0.,100./))

id_OFP_n2_int_temp     = register_diag_field ('ocean_model', 'OFP_n2_int_temp',            &
                         Time%model_time, 'OFP_2 interior region temperature', 'oC degree',&
                         missing_value=missing_value, range=(/-100.,100./))
id_OFP_n2_int_salt     = register_diag_field ('ocean_model', 'OFP_n2_int_salt',   &
                         Time%model_time, 'OFP_2 interior region salinity', 'psu',&
                         missing_value=missing_value, range=(/0.,100./))

id_OFP_n2_ent_temp     = register_diag_field ('ocean_model', 'OFP_n2_ent_temp',               &
                         Time%model_time, 'OFP_2 entrainment region temperature', 'oC degree',&
                         missing_value=missing_value, range=(/-100.,100./))
id_OFP_n2_ent_salt     = register_diag_field ('ocean_model', 'OFP_n2_ent_salt',       &
                         Time%model_time, 'OFP_2 entrainment  region salinity', 'psu',&
                         missing_value=missing_value, range=(/0.,100./))
id_OFP_n2_ent_trans    = register_diag_field ('ocean_model', 'OFP_n2_ent_trans',        &
                         Time%model_time, 'OFP_2 entrainment region transport', 'm^3/s',&
                         missing_value=missing_value, range=(/0.,100./))

id_OFP_n2_prd_temp     = register_diag_field ('ocean_model', 'OFP_n2_prd_temp',              &
                         Time%model_time, 'OFP_2 production region temperature', 'oC degree',&
                         missing_value=missing_value, range=(/-100.,100./))
id_OFP_n2_prd_salt     = register_diag_field ('ocean_model', 'OFP_n2_prd_salt',      &
                         Time%model_time, 'OFP_2 production  region salinity', 'psu',&
                         missing_value=missing_value, range=(/0.,100./))
id_OFP_n2_prd_trans    = register_diag_field ('ocean_model', 'OFP_n2_prd_trans',       &
                         Time%model_time, 'OFP_2 production region transport', 'm^3/s',&
                         missing_value=missing_value, range=(/0.,100./))
id_OFP_n2_prd_depth    = register_diag_field ('ocean_model', 'OFP_n2_prd_depth',&
                         Time%model_time, 'OFP_2 production depth', 'meter',    &
                         missing_value=missing_value, range=(/0.,10000./))

!------------------------

    stdoutunit=stdout();stdlogunit=stdlog() 

    pe=mpp_pe()
    root_pe=mpp_root_pe()

 if ( module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_overflow_mod (ocean_overflow_init): module already initialized')
 endif 

    module_is_initialized = .TRUE.

    call write_version_number(version, tagname)
#ifndef MOM_STATIC_ARRAYS
    call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
    nk = Grid%nk
#endif
    Dom => Domain
    Grd => Grid

         call mpp_clock_begin(id_OFP_init)

    num_prog_tracers = size(T_prog(:))
    dtime = dtim 
    p5dtimer = 0.5/dtime

 if (debug_this_module) then 
    write (stdoutunit,*) 'isd,ied,jsd,jed,isc,iec,jsc,jec',isd,ied,jsd,jed,isc,iec,jsc,jec
    write (stdoutunit,*) 'Dom%ioff,Dom%joff',Dom%ioff,Dom%joff
 endif

    call watermass_diag_init(Time, Dens)

!--
  n_info = find_field_index(MODEL_OCEAN,'overflow_ofp_info')
  n_src  = find_field_index(MODEL_OCEAN,'overflow_ofp_src')
  n_int  = find_field_index(MODEL_OCEAN,'overflow_ofp_int')
  n_ent  = find_field_index(MODEL_OCEAN,'overflow_ofp_ent')

  n_prd_line_01 = find_field_index(MODEL_OCEAN,'overflow_ofp_prd_line_01')
  n_prd_line_02 = find_field_index(MODEL_OCEAN,'overflow_ofp_prd_line_02')
  n_prd_line_03 = find_field_index(MODEL_OCEAN,'overflow_ofp_prd_line_03')
  n_prd_line_04 = find_field_index(MODEL_OCEAN,'overflow_ofp_prd_line_04')
  n_prd_line_05 = find_field_index(MODEL_OCEAN,'overflow_ofp_prd_line_05')
  n_prd_line_06 = find_field_index(MODEL_OCEAN,'overflow_ofp_prd_line_06')
  n_prd_line_07 = find_field_index(MODEL_OCEAN,'overflow_ofp_prd_line_07')
  n_prd_line_08 = find_field_index(MODEL_OCEAN,'overflow_ofp_prd_line_08')
  n_prd_line_09 = find_field_index(MODEL_OCEAN,'overflow_ofp_prd_line_09')
  n_prd_line_10 = find_field_index(MODEL_OCEAN,'overflow_ofp_prd_line_10')

  if (n_src < 1) then 
    write(stdoutunit,'(a)')' '
    write(stdoutunit,'(a)') &
    '==>Warning: ocean_overflow_OFP_init found n_src < 1 for overflow_OFP table.  Will NOT use ocean_overflow_OFP.'  
    write(stdoutunit,'(a)')' '
    Ocean_options%overflow_OFP = 'Did NOT use overflow_OFP.'
    return
  endif 


!----
  call get_field_info(n_info,fld_type,fld_name,model,n_OFP_info)
  call get_field_info(n_src,fld_type,fld_name,model,n_OFP_src)
  call get_field_info(n_int,fld_type,fld_name,model,n_OFP_int)
    if(n_OFP_int /= n_OFP_src) then 
       call mpp_error(NOTE,'==>From ocean_overflow_OFP_mod: n_OFP_int is not equal to n_OFP_src')
    endif    
  call get_field_info(n_ent,fld_type,fld_name,model,n_OFP_ent)
    if(n_OFP_ent /= n_OFP_src) then 
       call mpp_error(NOTE,'==>From ocean_overflow_OFP_mod: n_OFP_ent is not equal to n_OFP_src')
    endif    

  call get_field_info(n_prd_line_01,fld_type,fld_name,model,n_OFP_prd_line_01)
    if(n_OFP_prd_line_01 /= n_OFP_src) then 
      call mpp_error(NOTE,'==>From ocean_overflow_OFP_mod: n_OFP_prd_line_01 is not equal to n_OFP_src')
    endif    
  call get_field_info(n_prd_line_02,fld_type,fld_name,model,n_OFP_prd_line_02)
    if(n_OFP_prd_line_02 /= n_OFP_src) then 
      call mpp_error(NOTE,'==>From ocean_overflow_OFP_mod: n_OFP_prd_line_02 is not equal to n_OFP_src')
    endif    
  call get_field_info(n_prd_line_03,fld_type,fld_name,model,n_OFP_prd_line_03)
    if(n_OFP_prd_line_03 /= n_OFP_src) then 
      call mpp_error(NOTE,'==>From ocean_overflow_OFP_mod: n_OFP_prd_line_03 is not equal to n_OFP_src')
    endif    
  call get_field_info(n_prd_line_04,fld_type,fld_name,model,n_OFP_prd_line_04)
    if(n_OFP_prd_line_04 /= n_OFP_src) then 
      call mpp_error(NOTE,'==>From ocean_overflow_OFP_mod: n_OFP_prd_line_04 is not equal to n_OFP_src')
    endif    
  call get_field_info(n_prd_line_05,fld_type,fld_name,model,n_OFP_prd_line_05)
    if(n_OFP_prd_line_05 /= n_OFP_src) then 
      call mpp_error(NOTE,'==>From ocean_overflow_OFP_mod: n_OFP_prd_line_05 is not equal to n_OFP_src')
    endif    
  call get_field_info(n_prd_line_06,fld_type,fld_name,model,n_OFP_prd_line_06)
    if(n_OFP_prd_line_06 /= n_OFP_src) then 
      call mpp_error(NOTE,'==>From ocean_overflow_OFP_mod: n_OFP_prd_line_06 is not equal to n_OFP_src')
    endif    
  call get_field_info(n_prd_line_07,fld_type,fld_name,model,n_OFP_prd_line_07)
    if(n_OFP_prd_line_07 /= n_OFP_src) then 
      call mpp_error(NOTE,'==>From ocean_overflow_OFP_mod: n_OFP_prd_line_07 is not equal to n_OFP_src')
    endif    
  call get_field_info(n_prd_line_08,fld_type,fld_name,model,n_OFP_prd_line_08)
    if(n_OFP_prd_line_08 /= n_OFP_src) then 
      call mpp_error(NOTE,'==>From ocean_overflow_OFP_mod: n_OFP_prd_line_08 is not equal to n_OFP_src')
    endif    
  call get_field_info(n_prd_line_09,fld_type,fld_name,model,n_OFP_prd_line_09)
    if(n_OFP_prd_line_09 /= n_OFP_src) then 
      call mpp_error(NOTE,'==>From ocean_overflow_OFP_mod: n_OFP_prd_line_09 is not equal to n_OFP_src')
    endif    
  call get_field_info(n_prd_line_10,fld_type,fld_name,model,n_OFP_prd_line_10)
    if(n_OFP_prd_line_10 /= n_OFP_src) then 
      call mpp_error(NOTE,'==>From ocean_overflow_OFP_mod: n_OFP_prd_line_10 is not equal to n_OFP_src')
    endif    
!
!----- 
     n_OFP=n_OFP_src
!----
!
  allocate(OFP_methods_info(n_OFP))
  allocate(OFP_methods_src(n_OFP))
  allocate(OFP_methods_int(n_OFP))
  allocate(OFP_methods_ent(n_OFP))
  allocate(OFP_methods_prd_line_01(n_OFP))
  allocate(OFP_methods_prd_line_02(n_OFP))
  allocate(OFP_methods_prd_line_03(n_OFP))
  allocate(OFP_methods_prd_line_04(n_OFP))
  allocate(OFP_methods_prd_line_05(n_OFP))
  allocate(OFP_methods_prd_line_06(n_OFP))
  allocate(OFP_methods_prd_line_07(n_OFP))
  allocate(OFP_methods_prd_line_08(n_OFP))
  allocate(OFP_methods_prd_line_09(n_OFP))
  allocate(OFP_methods_prd_line_10(n_OFP))
  
  allocate(phi_ofp(n_OFP))
  allocate(hu_ofp(n_OFP))
  allocate(hs_ofp(n_OFP))
  allocate(ws_ofp(n_OFP))
  allocate(xssb_ofp(n_OFP))
  allocate(alpha_ofp(n_OFP))
  allocate(cdbot_ofp(n_OFP))
 
  phi_ofp            = 0.0
  hu_ofp             = 0.0
  hs_ofp             = 0.0
  ws_ofp             = 0.0
  xssb_ofp           = 0.0
  alpha_ofp          = 0.0
  cdbot_ofp          = 0.0
 
  allocate(src_vol_trace(n_OFP,num_prog_tracers))
  allocate(ent_vol_trace(n_OFP,num_prog_tracers))
  allocate(prd_vol_trace(n_OFP,num_prog_tracers))
 
  src_vol_trace      = 0.0
  ent_vol_trace      = 0.0
  prd_vol_trace      = 0.0
 
  allocate(src_ist(n_OFP))
  allocate(src_ied(n_OFP))
  allocate(src_jst(n_OFP))
  allocate(src_jed(n_OFP))
  allocate(src_kst(n_OFP))
  allocate(src_ked(n_OFP))
  allocate(src_isc(n_OFP))
  allocate(src_iec(n_OFP))
  allocate(src_jsc(n_OFP))
  allocate(src_jec(n_OFP))
 
  src_ist   = 0
  src_ied   = 0
  src_jst   = 0
  src_jed   = 0
  src_kst   = 0
  src_ked   = 0
  src_isc   = 0
  src_iec   = 0
  src_jsc   = 0
  src_jec   = 0
 
  allocate(int_ist(n_OFP))
  allocate(int_ied(n_OFP))
  allocate(int_jst(n_OFP))
  allocate(int_jed(n_OFP))
  allocate(int_kst(n_OFP))
  allocate(int_ked(n_OFP))
  allocate(int_isc(n_OFP))
  allocate(int_iec(n_OFP))
  allocate(int_jsc(n_OFP))
  allocate(int_jec(n_OFP))
 
  int_ist   = 0
  int_ied   = 0
  int_jst   = 0
  int_jed   = 0
  int_kst   = 0
  int_ked   = 0
  int_isc   = 0
  int_iec   = 0
  int_jsc   = 0
  int_jec   = 0
 
  allocate(ent_ist(n_OFP))
  allocate(ent_ied(n_OFP))
  allocate(ent_jst(n_OFP))
  allocate(ent_jed(n_OFP))
  allocate(ent_kst(n_OFP))
  allocate(ent_ked(n_OFP))
  allocate(ent_isc(n_OFP))
  allocate(ent_iec(n_OFP))
  allocate(ent_jsc(n_OFP))
  allocate(ent_jec(n_OFP))
 
  ent_ist   = 0
  ent_ied   = 0
  ent_jst   = 0
  ent_jed   = 0
  ent_kst   = 0
  ent_ked   = 0
  ent_isc   = 0
  ent_iec   = 0
  ent_jsc   = 0
  ent_jec   = 0
 
  allocate(prd_ist(n_OFP,10))
  allocate(prd_ied(n_OFP,10))
  allocate(prd_jst(n_OFP,10))
  allocate(prd_jed(n_OFP,10))
  allocate(prd_kst(n_OFP,10))
  allocate(prd_ked(n_OFP,10))
  allocate(prd_isc(n_OFP,10))
  allocate(prd_iec(n_OFP,10))
  allocate(prd_jsc(n_OFP,10))
  allocate(prd_jec(n_OFP,10))
 
  prd_ist   = 0
  prd_ied   = 0
  prd_jst   = 0
  prd_jed   = 0
  prd_kst   = 0
  prd_ked   = 0
  prd_isc   = 0
  prd_iec   = 0
  prd_jsc   = 0
  prd_jec   = 0
 
 
 
  allocate(working_comp_domains(n_OFP))
  allocate(working_src_comp_domains(n_OFP))
  allocate(working_ent_comp_domains(n_OFP))
  allocate(working_int_comp_domains(n_OFP))
  allocate(working_prd_comp_domains(n_OFP,10))
  allocate(occur_overflow_OFP(n_OFP))
 
  working_comp_domains      = .false.
  working_src_comp_domains  = .false.
  working_ent_comp_domains  = .false.
  working_int_comp_domains  = .false.
  working_prd_comp_domains  = .false.
  occur_overflow_OFP        = .false.
 
  allocate(max_i_bound(n_OFP))
  allocate(min_i_bound(n_OFP))
  allocate(max_j_bound(n_OFP))
  allocate(min_j_bound(n_OFP))
 
  max_i_bound   = 0
  min_i_bound   = 0
  max_j_bound   = 0
  min_j_bound   = 0
 
  allocate(max_ist(n_OFP))
  allocate(max_ied(n_OFP))
  allocate(max_jst(n_OFP))
  allocate(max_jed(n_OFP))
 
  max_ist   = 0
  max_ied   = 0
  max_jst   = 0
  max_jed   = 0

!----- n_OFP

  allocate(src_volT(n_OFP))   
  allocate(src_volS(n_OFP))   
  allocate(src_volP(n_OFP))   
  allocate(src_vol(n_OFP))    
  allocate(src_dat(n_OFP))   
  allocate(src_flux_Mass(n_OFP))  
  allocate(src_den(n_OFP))   

  src_volT        =0.0
  src_volS        =0.0
  src_volP        =0.0
  src_vol         =0.0
  src_dat         =0.0
  src_flux_Mass   =0.0
  src_den         =0.0


  allocate(int_volT(n_OFP))   
  allocate(int_volS(n_OFP))   
  allocate(int_volP(n_OFP))   
  allocate(int_vol(n_OFP))  
  allocate(int_dat(n_OFP))   
  allocate(int_den(n_OFP))   

  int_volT    =0.0
  int_volS    =0.0
  int_volP    =0.0
  int_vol     =0.0
  int_dat     =0.0
  int_den      =0.0

  allocate(grav_prim_src_ofp(n_OFP)) 
  allocate(coriolis_source_ofp(n_OFP))  
  allocate(trans_vol_ofp_src_raw(n_OFP)) 
  allocate(trans_vol_ofp_src(n_OFP))  
  allocate(vel_src_ofp(n_OFP))

  grav_prim_src_ofp     =0.0
  coriolis_source_ofp   =0.0
  trans_vol_ofp_src_raw =0.0
  trans_vol_ofp_src     =0.0
  vel_src_ofp           =0.0


  allocate(ent_volT(n_OFP))   
  allocate(ent_volS(n_OFP))   
  allocate(ent_volP(n_OFP))   
  allocate(ent_vol(n_OFP))    
  allocate(ent_dat(n_OFP))   
  allocate(ent_den(n_OFP))   

  ent_volT  =0.0
  ent_volS  =0.0
  ent_volP  =0.0
  ent_vol   =0.0
  ent_dat   =0.0
  ent_den   =0.0

  allocate(src_den_prim(n_OFP))  
  allocate(src_den_prim_int(n_OFP))  
  allocate(grav_prim_ent_ofp(n_OFP)) 
  allocate(trans_vol_ofp_ent(n_OFP)) 
  allocate(trans_vol_ofp_prd(n_OFP)) 
  allocate(trans_mass_ofp_ent(n_OFP)) 
  allocate(vel_ssb_ofp(n_OFP)) 
  allocate(vel_ave_ofp(n_OFP)) 
  allocate(coef_b_ofp(n_OFP)) 
  allocate(coef_c_ofp(n_OFP)) 
  allocate(coef_d_ofp(n_OFP)) 
  allocate(hssb_ofp(n_OFP))
  allocate(Fr_geo_ofp(n_OFP))
  allocate(mix_theta_ofp(n_OFP))
  allocate(volf_per_tvol_src(n_OFP))
  allocate(volf_per_tvol_ent(n_OFP))
  allocate(ent_flux_Mass(n_OFP))
  allocate(depth_prd_OFP(n_OFP))

  src_den_prim       =0.0
  src_den_prim_int   =0.0
  grav_prim_ent_ofp  =0.0
  trans_vol_ofp_ent  =0.0
  trans_vol_ofp_prd  =0.0
  trans_mass_ofp_ent =0.0     
  vel_ssb_ofp        =0.0     
  vel_ave_ofp        =0.0     
  coef_b_ofp         =0.0     
  coef_c_ofp         =0.0     
  coef_d_ofp         =0.0     
  hssb_ofp           =0.0     
  Fr_geo_ofp         =0.0     
  mix_theta_ofp      =0.0     
  volf_per_tvol_src  =0.0     
  volf_per_tvol_ent  =0.0     
  ent_flux_Mass      =0.0     
  depth_prd_OFP      =0.0     

  allocate(prd_line_num(n_OFP))
  allocate(prd_inject_line_num(n_OFP))
 
  prd_line_num        =0
  prd_inject_line_num =0

  allocate(prd_amb_T(n_OFP,10))
  allocate(prd_amb_S(n_OFP,10))
  allocate(prd_amb_P(n_OFP,10))
  allocate(prd_amb_rho(n_OFP,10))
  allocate(prd_vol(n_OFP,10))
  allocate(prd_dat(n_OFP,10))

  prd_amb_T     = 0.0
  prd_amb_S     = 0.0
  prd_amb_P     = 0.0
  prd_amb_rho   = 0.0
  prd_vol       = 0.0
  prd_dat       = 0.0

  allocate(prd_volT(n_OFP)) 
  allocate(prd_volS(n_OFP)) 
  allocate(prd_vol_rho(n_OFP)) 
  allocate(prd_flux_Mass(n_OFP))
  allocate(prd_inj_vol(n_OFP))
  allocate(volf_per_tvol_prd(n_OFP))

  prd_volT         = 0.0 
  prd_volS         = 0.0 
  prd_vol_rho      = 0.0 
  prd_flux_Mass    = 0.0
  prd_inj_vol      = 0.0
  volf_per_tvol_prd= 0.0

!
   allocate (src_flux_trace(n_OFP,num_prog_tracers))
   allocate (ent_flux_trace(n_OFP,num_prog_tracers))
   allocate (prd_flux_trace(n_OFP,num_prog_tracers))
!
   src_flux_trace = 0.0
   ent_flux_trace = 0.0
   prd_flux_trace = 0.0

!-----------------------------------------------------
    ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_overflow_OFP_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_overflow_OFP_nml')
#else
    ioun =  open_namelist_file()
    read (ioun,ocean_overflow_OFP_nml,IOSTAT=io_status)
    ierr = check_nml_error(io_status,'ocean_overflow_OFP_nml')
    call close_file (ioun)
#endif
    write (stdlogunit,ocean_overflow_OFP_nml)
    write (stdoutunit,'(/)') 
    write (stdoutunit,ocean_overflow_OFP_nml)

    if (PRESENT(debug) .and. .not. debug_this_module) then
       debug_this_module = debug
    endif 
    if(debug_this_module) then 
       write(stdoutunit,'(a)') '==>Note: running ocean_overflow_OFP_mod with debug_this_module=.true.'  
    endif 

    if(.not. use_this_module) then 
      call mpp_error(NOTE,&
      '==>From ocean_overflow_OFP_mod: NOT using overflow_OFP scheme.')
      Ocean_options%overflow_OFP = 'Did NOT use overflow_OFP scheme.'
      return
    else 
      if(vert_coordinate_type == TERRAIN_FOLLOWING) then 
          call mpp_error(FATAL, &
          '==>ocean_overflow_OFP_mod: this module is NOT for use with TERRAIN_FOLLOWING vert coodinates.')
      endif
      Ocean_options%overflow_OFP = 'Used the OFP scheme.'
      call mpp_error(NOTE,&
      '==>From ocean_overflow_OFP_mod: USING OFP scheme.')
    endif
    if(debug_this_module) then 
      call mpp_error(NOTE,'==>From ocean_overflow_OFP_mod: USING debug_this_module')
    endif 

!----- reading the data from field_table
  call get_field_methods(n_info,OFP_methods_info)

  do i=1, n_OFP
    parse_ok = parse(OFP_methods_info(i)%method_control,'plat',phi_ofp(i))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "platitude" error')
    parse_ok = parse(OFP_methods_info(i)%method_control,'hh',hu_ofp(i))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "h_ofp" error')
    parse_ok = parse(OFP_methods_info(i)%method_control,'ww',ws_ofp(i))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "w_ofp" error')
    parse_ok = parse(OFP_methods_info(i)%method_control,'wdth',xssb_ofp(i))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "x_ofp" error')
    parse_ok = parse(OFP_methods_info(i)%method_control,'aa',alpha_ofp(i))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "a_ofp" error')
    parse_ok = parse(OFP_methods_info(i)%method_control,'cd',cdbot_ofp(i))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "cd_ofp" error')

   if (debug_this_module) then 
     write(stdoutunit,*) 'n of n_OFP_info, phi_ofp,hu_ofp,ws_ofp,xssb_ofp,alpha_ofp,cdbot_ofp' &
               ,i,phi_ofp(i),hu_ofp(i),ws_ofp(i),xssb_ofp(i),alpha_ofp(i),cdbot_ofp(i)
   endif

  enddo
!
!-------- src
  call get_field_methods(n_src,OFP_methods_src)
!
  do i=1, n_OFP
     parse_ok = parse(OFP_methods_src(i)%method_control,'s_ist',src_ist(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "s_ist" error')
     parse_ok = parse(OFP_methods_src(i)%method_control,'s_ied',src_ied(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "s_ied" error')
     parse_ok = parse(OFP_methods_src(i)%method_control,'s_jst',src_jst(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "s_jst" error')
     parse_ok = parse(OFP_methods_src(i)%method_control,'s_jed',src_jed(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "s_jed" error')
     parse_ok = parse(OFP_methods_src(i)%method_control,'s_kst',src_kst(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "s_kst" error')
     parse_ok = parse(OFP_methods_src(i)%method_control,'s_ked',src_ked(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "s_ked" error')
!
       if (debug_this_module) then 
           write(*,*) 'n of n_OFP, src_ist,src_ied,src_jst,src_jed,src_kst,src_ked', &
                      i,src_ist(i),src_ied(i),src_jst(i),src_jed(i),src_kst(i),src_ked(i)
       endif
  enddo
!
!-------- int
  call get_field_methods(n_int,OFP_methods_int)
!
  do i=1, n_OFP
     parse_ok = parse(OFP_methods_int(i)%method_control,'i_ist',int_ist(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "i_ist" error')
     parse_ok = parse(OFP_methods_int(i)%method_control,'i_ied',int_ied(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "i_ied" error')
     parse_ok = parse(OFP_methods_int(i)%method_control,'i_jst',int_jst(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "i_jst" error')
     parse_ok = parse(OFP_methods_int(i)%method_control,'i_jed',int_jed(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "i_jed" error')
     parse_ok = parse(OFP_methods_int(i)%method_control,'i_kst',int_kst(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "i_kst" error')
     parse_ok = parse(OFP_methods_int(i)%method_control,'i_ked',int_ked(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "i_ked" error')
!
       if (debug_this_module) then 
           write(stdoutunit,*) 'n of n_OFP, int_ist,int_ied,int_jst,int_jed,int_kst,int_ked', &
                      i,int_ist(i),int_ied(i),int_jst(i),int_jed(i),int_kst(i),int_ked(i)
       endif
  enddo
!
!-------- ent
  call get_field_methods(n_ent,OFP_methods_ent)
!
  do i=1, n_OFP
     parse_ok = parse(OFP_methods_ent(i)%method_control,'e_ist',ent_ist(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "e_ist" error')
     parse_ok = parse(OFP_methods_ent(i)%method_control,'e_ied',ent_ied(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "e_ied" error')
     parse_ok = parse(OFP_methods_ent(i)%method_control,'e_jst',ent_jst(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "e_jst" error')
     parse_ok = parse(OFP_methods_ent(i)%method_control,'e_jed',ent_jed(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "e_jed" error')
     parse_ok = parse(OFP_methods_ent(i)%method_control,'e_kst',ent_kst(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "e_kst" error')
     parse_ok = parse(OFP_methods_ent(i)%method_control,'e_ked',ent_ked(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP table entry "e_ked" error')
!
       if (debug_this_module) then 
           write(stdoutunit,*) 'n of n_OFP, ent_ist,ent_ied,ent_jst,ent_jed,ent_kst,ent_ked', &
                      i,ent_ist(i),ent_ied(i),ent_jst(i),ent_jed(i),ent_kst(i),ent_ked(i)
       endif
  enddo
!
!--------- prd line 01
  call get_field_methods(n_prd_line_01,OFP_methods_prd_line_01)
!
  do i=1, n_OFP
     parse_ok = parse(OFP_methods_prd_line_01(i)%method_control,'p_ist',prd_ist(i,1))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_01 table entry "p_ist" error')
     parse_ok = parse(OFP_methods_prd_line_01(i)%method_control,'p_ied',prd_ied(i,1))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_01 table entry "p_ied" error')
     parse_ok = parse(OFP_methods_prd_line_01(i)%method_control,'p_jst',prd_jst(i,1))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_01 table entry "p_jst" error')
     parse_ok = parse(OFP_methods_prd_line_01(i)%method_control,'p_jed',prd_jed(i,1))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_01 table entry "p_jed" error')
     parse_ok = parse(OFP_methods_prd_line_01(i)%method_control,'p_kst',prd_kst(i,1))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_01 table entry "p_kst" error')
     parse_ok = parse(OFP_methods_prd_line_01(i)%method_control,'p_ked',prd_ked(i,1))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_01 table entry "p_ked" error')
!
       if (debug_this_module) then 
           write(*,*) 'n of n_OFP_prd_line_01, prd_ist,prd_ied,prd_jst,prd_jed,prd_kst,prd_ked', &
                      i,prd_ist(i,1),prd_ied(i,1),prd_jst(i,1),prd_jed(i,1),prd_kst(i,1),prd_ked(i,1)
       endif
  enddo
!
!--------- prd line 02
  call get_field_methods(n_prd_line_02,OFP_methods_prd_line_02)
!
  do i=1, n_OFP
     parse_ok = parse(OFP_methods_prd_line_02(i)%method_control,'p_ist',prd_ist(i,2))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_02 table entry "p_ist" error')
     parse_ok = parse(OFP_methods_prd_line_02(i)%method_control,'p_ied',prd_ied(i,2))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_02 table entry "p_ied" error')
     parse_ok = parse(OFP_methods_prd_line_02(i)%method_control,'p_jst',prd_jst(i,2))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_02 table entry "p_jst" error')
     parse_ok = parse(OFP_methods_prd_line_02(i)%method_control,'p_jed',prd_jed(i,2))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_02 table entry "p_jed" error')
     parse_ok = parse(OFP_methods_prd_line_02(i)%method_control,'p_kst',prd_kst(i,2))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_02 table entry "p_kst" error')
     parse_ok = parse(OFP_methods_prd_line_02(i)%method_control,'p_ked',prd_ked(i,2))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_02 table entry "p_ked" error')
!
       if (debug_this_module) then 
           write(stdoutunit,*) 'n of n_OFP_prd_line_02, prd_ist,prd_ied,prd_jst,prd_jed,prd_kst,prd_ked', &
                      i,prd_ist(i,2),prd_ied(i,2),prd_jst(i,2),prd_jed(i,2),prd_kst(i,2),prd_ked(i,2)
       endif
  enddo
!
!--------- prd line 03
  call get_field_methods(n_prd_line_03,OFP_methods_prd_line_03)
!
  do i=1, n_OFP
     parse_ok = parse(OFP_methods_prd_line_03(i)%method_control,'p_ist',prd_ist(i,3))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_03 table entry "p_ist" error')
     parse_ok = parse(OFP_methods_prd_line_03(i)%method_control,'p_ied',prd_ied(i,3))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_03 table entry "p_ied" error')
     parse_ok = parse(OFP_methods_prd_line_03(i)%method_control,'p_jst',prd_jst(i,3))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_03 table entry "p_jst" error')
     parse_ok = parse(OFP_methods_prd_line_03(i)%method_control,'p_jed',prd_jed(i,3))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_03 table entry "p_jed" error')
     parse_ok = parse(OFP_methods_prd_line_03(i)%method_control,'p_kst',prd_kst(i,3))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_03 table entry "p_kst" error')
     parse_ok = parse(OFP_methods_prd_line_03(i)%method_control,'p_ked',prd_ked(i,3))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_03 table entry "p_ked" error')
!
       if (debug_this_module) then 
           write(stdoutunit,*) 'n of n_OFP_prd_line_03, prd_ist,prd_ied,prd_jst,prd_jed,prd_kst,prd_ked', &
                i,prd_ist(i,3),prd_ied(i,3),prd_jst(i,3),prd_jed(i,3),prd_kst(i,3),prd_ked(i,3)
       endif
  enddo
!
!--------- prd line 04
  call get_field_methods(n_prd_line_04,OFP_methods_prd_line_04)
!
  do i=1, n_OFP
     parse_ok = parse(OFP_methods_prd_line_04(i)%method_control,'p_ist',prd_ist(i,4))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_04 table entry "p_ist" error')
     parse_ok = parse(OFP_methods_prd_line_04(i)%method_control,'p_ied',prd_ied(i,4))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_04 table entry "p_ied" error')
     parse_ok = parse(OFP_methods_prd_line_04(i)%method_control,'p_jst',prd_jst(i,4))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_04 table entry "p_jst" error')
     parse_ok = parse(OFP_methods_prd_line_04(i)%method_control,'p_jed',prd_jed(i,4))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_04 table entry "p_jed" error')
     parse_ok = parse(OFP_methods_prd_line_04(i)%method_control,'p_kst',prd_kst(i,4))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_04 table entry "p_kst" error')
     parse_ok = parse(OFP_methods_prd_line_04(i)%method_control,'p_ked',prd_ked(i,4))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_04 table entry "p_ked" error')
!
       if (debug_this_module) then 
            write(*,*) 'n of n_OFP_prd_line_04, prd_ist,prd_ied,prd_jst,prd_jed,prd_kst,prd_ked', &
                 i,prd_ist(i,4),prd_ied(i,4),prd_jst(i,4),prd_jed(i,4),prd_kst(i,4),prd_ked(i,4)
       endif
  enddo
!
!--------- prd line 05
  call get_field_methods(n_prd_line_05,OFP_methods_prd_line_05)
!
  do i=1, n_OFP
     parse_ok = parse(OFP_methods_prd_line_05(i)%method_control,'p_ist',prd_ist(i,5))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_05 table entry "p_ist" error')
     parse_ok = parse(OFP_methods_prd_line_05(i)%method_control,'p_ied',prd_ied(i,5))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_05 table entry "p_ied" error')
     parse_ok = parse(OFP_methods_prd_line_05(i)%method_control,'p_jst',prd_jst(i,5))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_05 table entry "p_jst" error')
     parse_ok = parse(OFP_methods_prd_line_05(i)%method_control,'p_jed',prd_jed(i,5))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_05 table entry "p_jed" error')
     parse_ok = parse(OFP_methods_prd_line_05(i)%method_control,'p_kst',prd_kst(i,5))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_05 table entry "p_kst" error')
     parse_ok = parse(OFP_methods_prd_line_05(i)%method_control,'p_ked',prd_ked(i,5))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_05 table entry "p_ked" error')
!
       if (debug_this_module) then 
            write(stdoutunit,*) 'n of n_OFP_prd_line_05, prd_ist,prd_ied,prd_jst,prd_jed,prd_kst,prd_ked', &
                 i,prd_ist(i,5),prd_ied(i,5),prd_jst(i,5),prd_jed(i,5),prd_kst(i,5),prd_ked(i,5)
       endif
  enddo
!
!--------- prd line 06
  call get_field_methods(n_prd_line_06,OFP_methods_prd_line_06)
!
  do i=1, n_OFP
     parse_ok = parse(OFP_methods_prd_line_06(i)%method_control,'p_ist',prd_ist(i,6))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_06 table entry "p_ist" error')
     parse_ok = parse(OFP_methods_prd_line_06(i)%method_control,'p_ied',prd_ied(i,6))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_06 table entry "p_ied" error')
     parse_ok = parse(OFP_methods_prd_line_06(i)%method_control,'p_jst',prd_jst(i,6))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_06 table entry "p_jst" error')
     parse_ok = parse(OFP_methods_prd_line_06(i)%method_control,'p_jed',prd_jed(i,6))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_06 table entry "p_jed" error')
     parse_ok = parse(OFP_methods_prd_line_06(i)%method_control,'p_kst',prd_kst(i,6))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_06 table entry "p_kst" error')
     parse_ok = parse(OFP_methods_prd_line_06(i)%method_control,'p_ked',prd_ked(i,6))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_06 table entry "p_ked" error')
!
       if (debug_this_module) then 
            write(stdoutunit,*) 'n of n_OFP_prd_line_06, prd_ist,prd_ied,prd_jst,prd_jed,prd_kst,prd_ked', &
                 i,prd_ist(i,6),prd_ied(i,6),prd_jst(i,6),prd_jed(i,6),prd_kst(i,6),prd_ked(i,6)
       endif
  enddo
!
!--------- prd line 07
  call get_field_methods(n_prd_line_07,OFP_methods_prd_line_07)
!
  do i=1, n_OFP
     parse_ok = parse(OFP_methods_prd_line_07(i)%method_control,'p_ist',prd_ist(i,7))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_07 table entry "p_ist" error')
     parse_ok = parse(OFP_methods_prd_line_07(i)%method_control,'p_ied',prd_ied(i,7))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_07 table entry "p_ied" error')
     parse_ok = parse(OFP_methods_prd_line_07(i)%method_control,'p_jst',prd_jst(i,7))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_07 table entry "p_jst" error')
     parse_ok = parse(OFP_methods_prd_line_07(i)%method_control,'p_jed',prd_jed(i,7))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_07 table entry "p_jed" error')
     parse_ok = parse(OFP_methods_prd_line_07(i)%method_control,'p_kst',prd_kst(i,7))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_07 table entry "p_kst" error')
     parse_ok = parse(OFP_methods_prd_line_07(i)%method_control,'p_ked',prd_ked(i,7))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_07 table entry "p_ked" error')
!
       if (debug_this_module) then 
            write(stdoutunit,*) 'n of n_OFP_prd_line_07, prd_ist,prd_ied,prd_jst,prd_jed,prd_kst,prd_ked', &
                 i,prd_ist(i,7),prd_ied(i,7),prd_jst(i,7),prd_jed(i,7),prd_kst(i,7),prd_ked(i,7)
       endif
  enddo
!
!--------- prd line 08
  call get_field_methods(n_prd_line_08,OFP_methods_prd_line_08)
!
  do i=1, n_OFP
     parse_ok = parse(OFP_methods_prd_line_08(i)%method_control,'p_ist',prd_ist(i,8))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_08 table entry "p_ist" error')
     parse_ok = parse(OFP_methods_prd_line_08(i)%method_control,'p_ied',prd_ied(i,8))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_08 table entry "p_ied" error')
     parse_ok = parse(OFP_methods_prd_line_08(i)%method_control,'p_jst',prd_jst(i,8))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_08 table entry "p_jst" error')
     parse_ok = parse(OFP_methods_prd_line_08(i)%method_control,'p_jed',prd_jed(i,8))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_08 table entry "p_jed" error')
     parse_ok = parse(OFP_methods_prd_line_08(i)%method_control,'p_kst',prd_kst(i,8))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_08 table entry "p_kst" error')
     parse_ok = parse(OFP_methods_prd_line_08(i)%method_control,'p_ked',prd_ked(i,8))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_08 table entry "p_ked" error')
!
       if (debug_this_module) then 
            write(stdoutunit,*) 'n of n_OFP_prd_line_08, prd_ist,prd_ied,prd_jst,prd_jed,prd_kst,prd_ked', &
                 i,prd_ist(i,8),prd_ied(i,8),prd_jst(i,8),prd_jed(i,8),prd_kst(i,8),prd_ked(i,8)
       endif
  enddo
!
!--------- prd line 09
  call get_field_methods(n_prd_line_09,OFP_methods_prd_line_09)
!
  do i=1, n_OFP
     parse_ok = parse(OFP_methods_prd_line_09(i)%method_control,'p_ist',prd_ist(i,9))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_09 table entry "p_ist" error')
     parse_ok = parse(OFP_methods_prd_line_09(i)%method_control,'p_ied',prd_ied(i,9))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_09 table entry "p_ied" error')
     parse_ok = parse(OFP_methods_prd_line_09(i)%method_control,'p_jst',prd_jst(i,9))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_09 table entry "p_jst" error')
     parse_ok = parse(OFP_methods_prd_line_09(i)%method_control,'p_jed',prd_jed(i,9))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_09 table entry "p_jed" error')
     parse_ok = parse(OFP_methods_prd_line_09(i)%method_control,'p_kst',prd_kst(i,9))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_09 table entry "p_kst" error')
     parse_ok = parse(OFP_methods_prd_line_09(i)%method_control,'p_ked',prd_ked(i,9))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_09 table entry "p_ked" error')
!
       if (debug_this_module) then 
            write(stdoutunit,*) 'n of n_OFP_prd_line_09, prd_ist,prd_ied,prd_jst,prd_jed,prd_kst,prd_ked', &
                 i,prd_ist(i,9),prd_ied(i,9),prd_jst(i,9),prd_jed(i,9),prd_kst(i,9),prd_ked(i,9)
       endif
  enddo
!
!--------- prd line 10
  call get_field_methods(n_prd_line_10,OFP_methods_prd_line_10)
!
  do i=1, n_OFP
     parse_ok = parse(OFP_methods_prd_line_10(i)%method_control,'p_ist',prd_ist(i,10))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_10 table entry "p_ist" error')
     parse_ok = parse(OFP_methods_prd_line_10(i)%method_control,'p_ied',prd_ied(i,10))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_10 table entry "p_ied" error')
     parse_ok = parse(OFP_methods_prd_line_10(i)%method_control,'p_jst',prd_jst(i,10))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_10 table entry "p_jst" error')
     parse_ok = parse(OFP_methods_prd_line_10(i)%method_control,'p_jed',prd_jed(i,10))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_10 table entry "p_jed" error')
     parse_ok = parse(OFP_methods_prd_line_10(i)%method_control,'p_kst',prd_kst(i,10))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_10 table entry "p_kst" error')
     parse_ok = parse(OFP_methods_prd_line_10(i)%method_control,'p_ked',prd_ked(i,10))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_overflow_OFP_mod: OFP_10 table entry "p_ked" error')
!
       if (debug_this_module) then 
            write(stdoutunit,*) 'n of n_OFP_prd_line_10, prd_ist,prd_ied,prd_jst,prd_jed,prd_kst,prd_ked', &
                 i,prd_ist(i,10),prd_ied(i,10),prd_jst(i,10),prd_jed(i,10),prd_kst(i,10),prd_ked(i,10)
       endif
  enddo
!
!----------- define whether working domain or not
  do nindx=1,n_OFP
      working_src_comp_domains(nindx) = .false.
      if (src_ist(nindx)*src_ied(nindx)*src_jst(nindx)*src_jed(nindx) /= 0) then
        if (((src_ist(nindx) >= isc .and. src_ist(nindx) <= iec) .or.      &
             (src_ied(nindx) >= isc .and. src_ied(nindx) <= iec))   .and.  &
            ((src_jst(nindx) >= jsc .and. src_jst(nindx) <= jec) .or.      &
             (src_jed(nindx) >= jsc .and. src_jed(nindx) <= jec)))   then
                  working_src_comp_domains(nindx) = .true.
        endif
      endif
      working_ent_comp_domains(nindx) = .false.
      if (ent_ist(nindx)*ent_ied(nindx)*ent_jst(nindx)*ent_jed(nindx) /= 0) then
        if (((ent_ist(nindx) >= isc .and. ent_ist(nindx) <= iec) .or.      &
             (ent_ied(nindx) >= isc .and. ent_ied(nindx) <= iec))   .and.  &
            ((ent_jst(nindx) >= jsc .and. ent_jst(nindx) <= jec) .or.      &
             (ent_jed(nindx) >= jsc .and. ent_jed(nindx) <= jec)))   then
                  working_ent_comp_domains(nindx) = .true.
        endif
      endif
      working_int_comp_domains(nindx) = .false.
      if (int_ist(nindx)*int_ied(nindx)*int_jst(nindx)*int_jed(nindx) /= 0) then
        if (((int_ist(nindx) >= isc .and. int_ist(nindx) <= iec) .or.      &
             (int_ied(nindx) >= isc .and. int_ied(nindx) <= iec))   .and.  &
            ((int_jst(nindx) >= jsc .and. int_jst(nindx) <= jec) .or.      &
             (int_jed(nindx) >= jsc .and. int_jed(nindx) <= jec)))   then
                  working_int_comp_domains(nindx) = .true.
        endif
      endif
    do mindx=1,10
      working_prd_comp_domains(nindx,mindx) = .false.
      if (prd_ist(nindx,mindx)*prd_ied(nindx,mindx)*prd_jst(nindx,mindx)*prd_jed(nindx,mindx) /= 0) then
        if (((prd_ist(nindx,mindx) >= isc .and. prd_ist(nindx,mindx) <= iec) .or.      &
             (prd_ied(nindx,mindx) >= isc .and. prd_ied(nindx,mindx) <= iec))   .and.  &
            ((prd_jst(nindx,mindx) >= jsc .and. prd_jst(nindx,mindx) <= jec) .or.      &
             (prd_jed(nindx,mindx) >= jsc .and. prd_jed(nindx,mindx) <= jec)))   then
                  working_prd_comp_domains(nindx,mindx) = .true.
        endif
      endif
    enddo ! mindx=1,10
  enddo !nindx=1,n_OFP
!
!-----------
   do nindx=1,n_OFP      
       working_comp_domains(nindx)=.false.
       if ( working_src_comp_domains(nindx) ) working_comp_domains(nindx)=.true.
       if ( working_ent_comp_domains(nindx) ) working_comp_domains(nindx)=.true.
       if ( working_int_comp_domains(nindx) ) working_comp_domains(nindx)=.true.
         do mindx=1,10
           if ( working_prd_comp_domains(nindx,mindx) ) working_comp_domains(nindx)=.true.
         enddo ! mindx=1,10
   enddo !nindx=1,n_OFP
!
!-----------
    working_comp_domains_all=.false.
    do nindx=1,n_OFP      
       if (working_comp_domains(nindx)) working_comp_domains_all=.true.
    enddo !nindx=1,n_OFP
!         
!---  define the max and min i,j for xhalo and jhalo
!---  no consideration for the arctic cyclonic boundaries

 if ( working_comp_domains_all) then  
   do nindx=1,n_OFP
      max_i_bound(nindx)=src_ied(nindx)
      min_i_bound(nindx)=src_ist(nindx)
      max_j_bound(nindx)=src_jed(nindx)
      min_j_bound(nindx)=src_jst(nindx)
!------- max_i_bound
      if (src_ied(nindx) /= 0 .and. src_ied(nindx) > max_i_bound(nindx)) max_i_bound(nindx)=src_ied(nindx)
      if (int_ied(nindx) /= 0 .and. int_ied(nindx) > max_i_bound(nindx)) max_i_bound(nindx)=int_ied(nindx)
      if (ent_ied(nindx) /= 0 .and. ent_ied(nindx) > max_i_bound(nindx)) max_i_bound(nindx)=ent_ied(nindx)
         do mindx=1,10
            if (prd_ied(nindx,mindx) /= 0 .and. prd_ied(nindx,mindx) > max_i_bound(nindx)) &
                max_i_bound(nindx)=prd_ied(nindx,mindx)
         enddo ! mindx=1,10
!------- min_i_bound
      if (src_ist(nindx) /= 0 .and. src_ist(nindx) < min_i_bound(nindx)) min_i_bound(nindx)=src_ist(nindx)
      if (int_ist(nindx) /= 0 .and. int_ist(nindx) < min_i_bound(nindx)) min_i_bound(nindx)=int_ist(nindx)
      if (ent_ist(nindx) /= 0 .and. ent_ist(nindx) < min_i_bound(nindx)) min_i_bound(nindx)=ent_ist(nindx)
         do mindx=1,10
            if (prd_ist(nindx,mindx) /= 0 .and. prd_ist(nindx,mindx) < min_i_bound(nindx)) &
                min_i_bound(nindx)=prd_ist(nindx,mindx)
         enddo ! mindx=1,10
!------- max_j_bound
      if (src_jed(nindx) /= 0 .and. src_jed(nindx) > max_j_bound(nindx)) max_j_bound(nindx)=src_jed(nindx)
      if (int_jed(nindx) /= 0 .and. int_jed(nindx) > max_j_bound(nindx)) max_j_bound(nindx)=int_jed(nindx)
      if (ent_jed(nindx) /= 0 .and. ent_jed(nindx) > max_j_bound(nindx)) max_j_bound(nindx)=ent_jed(nindx)
         do mindx=1,10
            if (prd_jed(nindx,mindx) /= 0 .and. prd_jed(nindx,mindx) > max_j_bound(nindx)) &
                max_j_bound(nindx)=prd_jed(nindx,mindx)
         enddo ! mindx=1,10
!------- min_j_bound
      if (src_jst(nindx) /= 0 .and. src_jst(nindx) < min_j_bound(nindx)) min_j_bound(nindx)=src_jst(nindx)
      if (int_jst(nindx) /= 0 .and. int_jst(nindx) < min_j_bound(nindx)) min_j_bound(nindx)=int_jst(nindx)
      if (ent_jst(nindx) /= 0 .and. ent_jst(nindx) < min_j_bound(nindx)) min_j_bound(nindx)=ent_jst(nindx)
         do mindx=1,10
            if (prd_jst(nindx,mindx) /= 0 .and. prd_jst(nindx,mindx) < min_j_bound(nindx)) &
                min_j_bound(nindx)=prd_jst(nindx,mindx)
         enddo ! mindx=1,10
   enddo !nindx=1,n_OFP
!
   if (debug_this_module) then
     do nindx=1,n_OFP
       write(*,*) 'pe,nindx,max_i_bound,min_i_bound,max_j_bound,min_j_bound ' &
                ,pe,nindx,max_i_bound(nindx),min_i_bound(nindx),max_j_bound(nindx),min_j_bound(nindx)
     enddo !nindx=1,n_OFP
    endif
!
 endif ! (working_comp_domains) 


!----- define halo region
!----- find the halo range on the maximum region of OFP at the comp. domain

!--- mpp_sum
!-!
        OFP_xhalo=2
        OFP_yhalo=2

!----- define computation index
 do nindx=1,n_OFP
  
   if ( working_src_comp_domains(nindx) ) then

       src_isc(nindx)  = isc    
       if (src_ist(nindx) >= isc ) src_isc(nindx)=src_ist(nindx)

       src_iec(nindx)  = iec    
       if (src_ied(nindx) <= iec ) src_iec(nindx)=src_ied(nindx)

       src_jsc(nindx)  = jsc    
       if (src_jst(nindx) >= jsc ) src_jsc(nindx)=src_jst(nindx)

       src_jec(nindx)  = jec    
       if (src_jed(nindx) <= jec ) src_jec(nindx)=src_jed(nindx)

      if (debug_this_module) then
         write(*,*) "pe, nindx,src_isc(nindx),src_iec(nindx)",pe, nindx,src_isc(nindx),src_iec(nindx)
      endif
   endif

   if ( working_ent_comp_domains(nindx) ) then

       ent_isc(nindx)  = isc    
       if (ent_ist(nindx) >= isc ) ent_isc(nindx)=ent_ist(nindx)

       ent_iec(nindx)  = iec    
       if (ent_ied(nindx) <= iec ) ent_iec(nindx)=ent_ied(nindx)

       ent_jsc(nindx)  = jsc    
       if (ent_jst(nindx) >= jsc ) ent_jsc(nindx)=ent_jst(nindx)

       ent_jec(nindx)  = jec    
       if (ent_jed(nindx) <= jec ) ent_jec(nindx)=ent_jed(nindx)

      if (debug_this_module) then
         write(*,*) "pe, nindx,ent_isc(nindx),ent_iec(nindx)",pe, nindx,ent_isc(nindx),ent_iec(nindx)
      endif
   endif


   if ( working_int_comp_domains(nindx) ) then

       int_isc(nindx)  = isc    
       if (int_ist(nindx) >= isc ) int_isc(nindx)=int_ist(nindx)

       int_iec(nindx)  = iec    
       if (int_ied(nindx) <= iec ) int_iec(nindx)=int_ied(nindx)

       int_jsc(nindx)  = jsc    
       if (int_jst(nindx) >= jsc ) int_jsc(nindx)=int_jst(nindx)

       int_jec(nindx)  = jec    
       if (int_jed(nindx) <= jec ) int_jec(nindx)=int_jed(nindx)

      if (debug_this_module) then
         write(*,*) "pe, nindx,int_isc(nindx),int_iec(nindx)",pe, nindx,int_isc(nindx),int_iec(nindx)
      endif
   endif

   do mindx=1,10

     if ( working_prd_comp_domains(nindx,mindx) ) then

       prd_isc(nindx,mindx)  = isc    
       if (prd_ist(nindx,mindx) >= isc ) prd_isc(nindx,mindx)=prd_ist(nindx,mindx)

       prd_iec(nindx,mindx)  = iec    
       if (prd_ied(nindx,mindx) <= iec ) prd_iec(nindx,mindx)=prd_ied(nindx,mindx)

       prd_jsc(nindx,mindx)  = jsc    
       if (prd_jst(nindx,mindx) >= jsc ) prd_jsc(nindx,mindx)=prd_jst(nindx,mindx)

       prd_jec(nindx,mindx)  = jec    
       if (prd_jed(nindx,mindx) <= jec ) prd_jec(nindx,mindx)=prd_jed(nindx,mindx)

       if (debug_this_module) then
         write(*,*) "pe, nindx,mindx,prd_isc(nindx,mindx)",pe, nindx,mindx,prd_isc(nindx,mindx)
      endif
    endif
   enddo ! mindx=1,10

 enddo !nindx=1,n_OFP


!-------------------------------------------------------------------------------
  if (debug_this_module) then 
    do nindx=1,n_OFP
      if (working_comp_domains(nindx)) then
        write(stdoutunit,*) ' test for working_comp_domains(nindx)---> .true. pe,nindx*************', pe,nindx
        write(stdoutunit,*) 'pe,max_i_bound,min_i_bound,max_j_bound,min_j_bound ' &
                   ,pe,max_i_bound(nindx),min_i_bound(nindx),max_j_bound(nindx),min_j_bound(nindx)
        write(stdoutunit,*) ' pe,isc,iec,jsc,jec',pe,isc,iec,jsc,jec
        write(stdoutunit,*) ' pe,OFP_xhalo,OFP_yhalo',pe,OFP_xhalo,OFP_yhalo
        if (working_src_comp_domains(nindx)) then
          write(stdoutunit,*) ' working_src_comp_domains(nindx)---> .true. at pe,nindx*************', pe,nindx
          write(stdoutunit,*) ' pe,isc,iec,jsc,jec',pe,isc,iec,jsc,jec
          write(stdoutunit,*) ' pe,src_ist(nindx),src_ied(nindx),src_jst(nindx),src_jed(nindx)' &
                      ,pe,src_ist(nindx),src_ied(nindx),src_jst(nindx),src_jed(nindx)
          write(stdoutunit,*) ' pe,src_isc(nindx),src_iec(nindx),src_jsc(nindx),src_jec(nindx)' &
                      ,pe,src_isc(nindx),src_iec(nindx),src_jsc(nindx),src_jec(nindx)
        endif
        if (working_ent_comp_domains(nindx)) then       
          write(stdoutunit,*) ' working_ent_comp_domains(nindx)---> .true. at pe,,nindx*************', pe,nindx
          write(stdoutunit,*) ' pe,isc,iec,jsc,jec',pe,isc,iec,jsc,jec
          write(stdoutunit,*) ' pe,ent_ist(nindx),ent_ied(nindx),ent_jst(nindx),ent_jed(nindx)' &
                      ,pe,ent_ist(nindx),ent_ied(nindx),ent_jst(nindx),ent_jed(nindx)
          write(stdoutunit,*) ' pe,ent_isc(nindx),ent_iec(nindx),ent_jsc(nindx),ent_jec(nindx)' &
                      ,pe,ent_isc(nindx),ent_iec(nindx),ent_jsc(nindx),ent_jec(nindx)
         endif
         do mindx=1,10
          if (working_prd_comp_domains(nindx,mindx)) then
            write(stdoutunit,*) ' working_prd_comp_domains(nindx,mindx)---> .true. at pe,nindx,mindx**', &
                        pe,nindx,mindx
            write(stdoutunit,*) 'pe,isc,iec,jsc,jec',pe,isc,iec,jsc,jec
            write(stdoutunit,*)   &
                 'pe,prd_ist(nindx,mindx),prd_ied(nindx,mindx),prd_jst(nindx,mindx),prd_jed(nindx,mindx) ' &
                 ,pe,prd_ist(nindx,mindx),prd_ied(nindx,mindx),prd_jst(nindx,mindx),prd_jed(nindx,mindx)
            write(stdoutunit,*)   &
                 'pe,prd_isc(nindx,mindx),prd_iec(nindx,mindx),prd_jsc(nindx,mindx),prd_jec(nindx,mindx) ' &
                 ,pe,prd_isc(nindx,mindx),prd_iec(nindx,mindx),prd_jsc(nindx,mindx),prd_jec(nindx,mindx)
           endif 
          enddo ! mindx=1,10
      endif
    enddo ! nindx=1,n_OFP
  endif


!--- define OFP domain 
  call set_ocean_domain(OFP_domain,Grd,xhalo=OFP_xhalo,yhalo=OFP_yhalo,&
                        name='overflow_OFP',maskmap=Dom%maskmap)

!-------------------------------------------------------------------------------
!
!--------- allocation
!
!--- for mpp_sum

   allocate (kmto(isc-OFP_xhalo:iec+OFP_xhalo,jsc-OFP_yhalo:jec+OFP_yhalo))
   allocate (dato(isc-OFP_xhalo:iec+OFP_xhalo,jsc-OFP_yhalo:jec+OFP_yhalo))
   allocate (dxto(isc-OFP_xhalo:iec+OFP_xhalo,jsc-OFP_yhalo:jec+OFP_yhalo))
   allocate (dyto(isc-OFP_xhalo:iec+OFP_xhalo,jsc-OFP_yhalo:jec+OFP_yhalo))
   allocate (hto(isc-OFP_xhalo:iec+OFP_xhalo,jsc-OFP_yhalo:jec+OFP_yhalo))
   allocate (tmsko(isc-OFP_xhalo:iec+OFP_xhalo,jsc-OFP_yhalo:jec+OFP_yhalo,1:nk))
   allocate (umsko(isc-OFP_xhalo:iec+OFP_xhalo,jsc-OFP_yhalo:jec+OFP_yhalo,1:nk))
   allocate (dzto(isc-OFP_xhalo:iec+OFP_xhalo,jsc-OFP_yhalo:jec+OFP_yhalo,1:nk))
   allocate (rho_dzto(isc-OFP_xhalo:iec+OFP_xhalo,jsc-OFP_yhalo:jec+OFP_yhalo,1:nk))
   allocate (tempo(isc-OFP_xhalo:iec+OFP_xhalo,jsc-OFP_yhalo:jec+OFP_yhalo,1:nk))
   allocate (salto(isc-OFP_xhalo:iec+OFP_xhalo,jsc-OFP_yhalo:jec+OFP_yhalo,1:nk))
   allocate (traceo(isc-OFP_xhalo:iec+OFP_xhalo,jsc-OFP_yhalo:jec+OFP_yhalo,1:nk,1:num_prog_tracers))
   allocate (presso(isc-OFP_xhalo:iec+OFP_xhalo,jsc-OFP_yhalo:jec+OFP_yhalo,1:nk))
   allocate (depth_zto(isc-OFP_xhalo:iec+OFP_xhalo,jsc-OFP_yhalo:jec+OFP_yhalo,1:nk))
!

   kmto      = 0
   dato      = 0.0
   dxto      = 0.0
   dyto      = 0.0
   hto       = 0.0
   dzto      = 0.0
   tmsko     = 0.0
   umsko     = 0.0
   rho_dzto  = 0.0
   tempo     = 0.0
   salto     = 0.0
   traceo    = 0.0
   presso    = 0.0
   depth_zto = 0.0
!
    do j=jsc,jec
      do i=isc,iec
         kmto(i,j)= Grd%kmt(i,j)   
         dato(i,j)= Grd%dat(i,j)
         dxto(i,j)= Grd%dxt(i,j)
         dyto(i,j)= Grd%dyt(i,j)
         hto(i,j) = Grd%ht(i,j)
      enddo
    enddo 

  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
         tmsko(i,j,k)=Grd%tmask(i,j,k)
         umsko(i,j,k)=Grd%umask(i,j,k)
      enddo
    enddo 
  enddo 
!

!-- no update domains

   ! integer variables
   call mpp_update_domains (kmto, OFP_domain%domain2d, complete=.true.)
   ! real variables
   call mpp_update_domains (dato, OFP_domain%domain2d, complete=.false.)
   call mpp_update_domains (dxto, OFP_domain%domain2d, complete=.false.)
   call mpp_update_domains (dyto, OFP_domain%domain2d, complete=.false.)
   call mpp_update_domains (hto,  OFP_domain%domain2d, complete=.true.)
   call mpp_update_domains (tmsko,OFP_domain%domain2d, complete=.false.)
   call mpp_update_domains (umsko,OFP_domain%domain2d, complete=.true.)


!-------------------------------------------------------------------------------
   do nindx=1,n_OFP 
      if (working_comp_domains(nindx)) then

          do j=src_jsc(nindx),src_jec(nindx)
             do i=src_isc(nindx),src_iec(nindx)
                if (src_ked(nindx) > kmto(i,j)) then
                    write (stdoutunit,*) ' ked > kmt, at pe,nindx,i,j', pe,nindx,i,j
                    write (stdoutunit,*) ' src_ked(nindx),kmto(i,j)',src_ked(nindx),kmto(i,j)
                endif
             enddo
          enddo

          do j=int_jsc(nindx),int_jec(nindx)
             do i=int_isc(nindx),int_iec(nindx)
                if (int_ked(nindx) > kmto(i,j)) then
                    write (stdoutunit,*) ' ked > kmt, at pe,nindx,i,j', pe,nindx,i,j
                    write (stdoutunit,*) ' int_ked(nindx),kmto(i,j)',int_ked(nindx),kmto(i,j)
                endif
             enddo
          enddo

          do j=ent_jsc(nindx),ent_jec(nindx)
             do i=ent_isc(nindx),ent_iec(nindx)
                if (ent_ked(nindx) > kmto(i,j)) then
                    write (stdoutunit,*) ' ked > kmt, at pe,nindx,i,j', pe,nindx,i,j
                    write (stdoutunit,*) ' ent_ked(nindx),kmto(i,j)',ent_ked(nindx),kmto(i,j)
                endif
             enddo
          enddo

          do m=1,10
             if (prd_jst(nindx,m) /= 0) then
                 do j=prd_jsc(nindx,m),prd_jec(nindx,m)
                    do i=prd_isc(nindx,m),prd_iec(nindx,m)
                       if (prd_ked(nindx,m) > kmto(i,j)) then
                           write (stdoutunit,*) ' ked > kmt, at pe,m, nindx,i,j',pe, m, nindx,i,j
                           write (stdoutunit,*) ' prd_ked(nindx,m),kmto(i,j)',prd_ked(nindx,m),kmto(i,j)
                       endif
                    enddo
                 enddo
             endif !(prd_jst(nindx,m) /= 0)
          enddo


      endif   !(working_comp_domains(nindx))
   enddo      ! nindx=1,n_OFP

   call mpp_clock_end(id_OFP_init)


  end subroutine ocean_overflow_OFP_init
! </SUBROUTINE> NAME="ocean_overflow_OFP_init"
!
!
!#######################################################################
! <SUBROUTINE NAME="overflow_OFP">
!
! <DESCRIPTION>
! Compute overflow process 
!
!
! </DESCRIPTION>
!
subroutine overflow_OFP (Time, Thickness, T_prog, Dens, index_temp, index_salt)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(inout) :: Thickness
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  type(ocean_density_type),       intent(in)    :: Dens
  integer,                        intent(in)    :: index_temp
  integer,                        intent(in)    :: index_salt

  real                    :: prd_P,prd_T,prd_S
  real                    :: sum_tmp

  integer                 :: tau, taum1, taup1
  integer                 :: i, j, k, m, n, nid, ijk_check
  integer                 :: pe, root_pe
  integer                 :: stdoutunit

  character(len=256)      :: errMsg
  
  pe         = mpp_pe()
  root_pe    = mpp_root_pe()
  stdoutunit = stdout()


  if(.not. use_this_module) return 

  call mpp_clock_begin(id_OFP_main)

  if(.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_overflow_mod (overflow): module must be initialized')
  endif 

  tau     = Time%tau
  taum1   = Time%taum1
  taup1   = Time%taup1

  dzto(:,:,:)      =0.0
  rho_dzto(:,:,:)  =0.0
  depth_zto(:,:,:) =0.0
  tempo(:,:,:)     =0.0
  salto(:,:,:)     =0.0
  presso(:,:,:)    =0.0
  traceo(:,:,:,:)  =0.0

  wrk1_v(:,:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec   
        do i=isc,iec
           dzto(i,j,k)      = Thickness%dzt(i,j,k)
           depth_zto(i,j,k) = Thickness%depth_zt(i,j,k)
           rho_dzto(i,j,k)  = Thickness%rho_dzt(i,j,k,tau)
           tempo(i,j,k)     = T_prog(index_temp)%field(i,j,k,tau)
           salto(i,j,k)     = T_prog(index_salt)%field(i,j,k,tau)
           presso(i,j,k)    = Dens%pressure_at_depth(i,j,k)
           wrk1_v(i,j,k,1)  = T_prog(index_temp)%th_tendency(i,j,k)
           wrk1_v(i,j,k,2)  = T_prog(index_salt)%th_tendency(i,j,k)
        enddo
     enddo
  enddo

  do n=1,num_prog_tracers
     do k=1,nk
        do j=jsc,jec   
           do i=isc,iec
              traceo(i,j,k,n)=T_prog(n)%field(i,j,k,tau)
           enddo
        enddo
     enddo
  enddo

  call mpp_clock_begin(id_OFP_update_1)
    call mpp_update_domains (dzto(:,:,:),       OFP_domain%domain2d, complete=.false.)
    call mpp_update_domains (depth_zto(:,:,:),  OFP_domain%domain2d, complete=.false.)
    call mpp_update_domains (rho_dzto(:,:,:),   OFP_domain%domain2d, complete=.false.)
    call mpp_update_domains (tempo(:,:,:),      OFP_domain%domain2d, complete=.false.)
    call mpp_update_domains (salto(:,:,:),      OFP_domain%domain2d, complete=.false.)
    call mpp_update_domains (presso(:,:,:),     OFP_domain%domain2d, complete=.false.)
    do n=1,num_prog_tracers
      call mpp_update_domains (traceo(:,:,:,n), OFP_domain%domain2d, complete=(n==num_prog_tracers))
    enddo
  call mpp_clock_end(id_OFP_update_1)


!-----------------------------------------------------------------------
!---- initial lizing for every PEs
  src_vol_trace(:,:) = 0.0
  ent_vol_trace(:,:) = 0.0
  prd_vol_trace(:,:) = 0.0

  src_volT(:)     =0.0
  src_volS(:)     =0.0
  src_volP(:)     =0.0
  src_vol(:)      =0.0
  src_dat(:)      =0.0
  src_den(:)      =0.0
  src_flux_Mass(:)=0.0

  int_volT(:)     =0.0
  int_volS(:)     =0.0
  int_vol(:)      =0.0
  int_volP(:)     =0.0
  int_dat(:)      =0.0
  int_den(:)      =0.0

  grav_prim_src_ofp(:)     =0.0
  hs_ofp(:)                =0.0
  coriolis_source_ofp(:)   =0.0
  trans_vol_ofp_src_raw(:) =0.0
  trans_vol_ofp_src(:)     =0.0

  ent_volT(:)=0.0
  ent_volS(:)=0.0
  ent_volP(:)=0.0
  ent_vol (:)=0.0
  ent_dat (:)=0.0
  ent_den (:)=0.0

  src_den_prim(:)      =0.0 
  src_den_prim_int(:)  =0.0 
  grav_prim_ent_ofp(:) =0.0
  trans_vol_ofp_ent(:) =0.0
  trans_mass_ofp_ent(:)=0.0 
  vel_ssb_ofp(:)       =0.0
  vel_ave_ofp(:)       =0.0
  coef_b_ofp(:)        =0.0
  coef_c_ofp(:)        =0.0
  coef_d_ofp(:)        =0.0
  hssb_ofp(:)          =0.0
  Fr_geo_ofp(:)        =0.0
  mix_theta_ofp(:)     =0.0
  volf_per_tvol_src(:) =0.0
  volf_per_tvol_ent(:) =0.0
  ent_flux_Mass(:)     =0.0
  prd_line_num(:)      =0
  prd_vol(:,:)         =0.0
  prd_amb_T(:,:)       =0.0
  prd_amb_S(:,:)       =0.0
  prd_amb_P(:,:)       =0.0
  prd_dat(:,:)         =0.0
  prd_amb_rho(:,:)     =0.0
  prd_volT(:)          =0.0
  prd_volS(:)          =0.0
  prd_vol_rho(:)       =0.0
  prd_flux_Mass(:)     =0.0
  prd_inj_vol(:)       =0.0
  volf_per_tvol_prd(:) =0.0 
  trans_vol_ofp_prd(:) =0.0
  depth_prd_OFP(:)     =0.0


!-------main loop

do  nid=1,n_OFP

   if (working_src_comp_domains(nid)) then

       !-----Thickness%dzt(i,j,k):dzto is defined in tau
       do k=src_kst(nid),src_ked(nid)
          do j=src_jsc(nid),src_jec(nid)
             do i=src_isc(nid),src_iec(nid)
                src_volT(nid)=src_volT(nid)+tempo(i,j,k)*dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)
                src_volS(nid)=src_volS(nid)+salto(i,j,k)*dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)
                src_vol(nid) =src_vol(nid) + dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)
             enddo
          enddo
       enddo

       do n=1,num_prog_tracers
          sum_tmp=0.0
          do k=src_kst(nid),src_ked(nid)
             do j=src_jsc(nid),src_jec(nid)
                do i=src_isc(nid),src_iec(nid)
                   sum_tmp=sum_tmp+traceo(i,j,k,n)*dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)
                enddo
             enddo
          enddo
          src_vol_trace(nid,n)=sum_tmp
       enddo

       !---  for the pressure, regionally averaged P over the whole boxes 
       do k=src_kst(nid),src_ked(nid)
          do j=src_jsc(nid),src_jec(nid)
             do i=src_isc(nid),src_iec(nid)
                src_volP(nid) = src_volP(nid) + presso(i,j,k)*dato(i,j)*tmsko(i,j,k)
                src_dat(nid)  = src_dat(nid) + dato(i,j)*tmsko(i,j,k)
             enddo
          enddo
       enddo

   endif !(working_src_comp_domains(nid))

 
   call mpp_sum (src_volT(nid))
   call mpp_sum (src_volS(nid))
   call mpp_sum (src_volP(nid))
   call mpp_sum (src_vol (nid))
   call mpp_sum (src_dat (nid))
   call mpp_sum (src_vol_trace(nid,1:num_prog_tracers), num_prog_tracers)

   !--- truncates at 4 bytes
   src_volT(nid) = real( src_volT(nid), kind=4 )
   src_volS(nid) = real( src_volS(nid), kind=4 )
   src_volP(nid) = real( src_volP(nid), kind=4 )
   src_vol (nid) = real( src_vol (nid), kind=4 )
   src_dat (nid) = real( src_dat (nid), kind=4 )

    do n=1,num_prog_tracers
      src_vol_trace(nid,n)=real( src_vol_trace(nid,n), kind=4 )
    enddo


    if(src_vol(nid) .le. 0.0 ) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_overflow_OFP_mod (overflow_OFP): source region has no sea cell')
    endif 

    src_volT(nid) = src_volT(nid)/src_vol(nid)
    src_volS(nid) = src_volS(nid)/src_vol(nid)
    src_volP(nid) = src_volP(nid)/src_dat(nid)
     
    do n=1,num_prog_tracers
      src_vol_trace(nid,n)=src_vol_trace(nid,n)/src_vol(nid)
    enddo

    !--- define averaged density at source region
    src_den(nid)=density(src_volS(nid),src_volT(nid),src_volP(nid))  

   if (debug_this_module .and. working_src_comp_domains(nid)) then 
      write(stdoutunit,*) 'nid,src_volT(nid),src_volS(nid),src_volP(nid),src_den(nid)' & 
                          ,nid,src_volT(nid),src_volS(nid),src_volP(nid),src_den(nid)
   endif

!--- define INT density
!--- thickness is still before the update, so dzt is at tau.

!  int_volT=0.0
!  int_volS=0.0
!  int_vol=0.0
!  int_volP=0.0
!  int_dat=0.0
!

   if (working_int_comp_domains(nid)) then

       do k=int_kst(nid),int_ked(nid)
          do j=int_jsc(nid),int_jec(nid)
             do i=int_isc(nid),int_iec(nid)
                int_volT(nid)=int_volT(nid)+tempo(i,j,k)*dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)
                int_volS(nid)=int_volS(nid)+salto(i,j,k)*dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)
                int_vol(nid) =int_vol(nid) + dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)
             enddo
          enddo
       enddo

       do k=int_kst(nid),int_ked(nid)
          do j=int_jsc(nid),int_jec(nid)
             do i=int_isc(nid),int_iec(nid)
                int_volP(nid) = int_volP(nid) + presso(i,j,k)*dato(i,j)*tmsko(i,j,k)
                int_dat(nid)  = int_dat(nid) + dato(i,j)*tmsko(i,j,k)
             enddo
          enddo
       enddo

   endif !(working_int_comp_domains(nid))

   call mpp_sum (int_volT(nid))
   call mpp_sum (int_volS(nid))
   call mpp_sum (int_volP(nid))
   call mpp_sum (int_vol (nid))
   call mpp_sum (int_dat (nid))

   !--- truncates at 4 bytes
   int_volT(nid) = real( int_volT(nid), kind=4 )
   int_volS(nid) = real( int_volS(nid), kind=4 )
   int_volP(nid) = real( int_volP(nid), kind=4 )
   int_vol (nid) = real( int_vol (nid), kind=4 )
   int_dat (nid) = real( int_dat (nid), kind=4 )


   if(int_vol(nid) .le. 0.0 ) then 
       call mpp_error(FATAL, &
            '==>Error from ocean_overflow_OFP_mod (overflow_OFP): interior region has no sea cell')
   endif

   int_volT(nid) = int_volT(nid)/int_vol(nid)
   int_volS(nid) = int_volS(nid)/int_vol(nid)
   int_volP(nid)=int_volP(nid)/int_dat(nid)

   int_den(nid)=density(int_volS(nid),int_volT(nid),int_volP(nid))  
   src_den_prim_int(nid)=density(src_volS(nid),src_volT(nid),int_volP(nid))

   if (debug_this_module .and. working_int_comp_domains(nid)) then 
       write(stdoutunit,*) 'nid,int_volT,int_volS,int_volP,int_den' & 
            ,nid,int_volT,int_volS,int_volP,int_den
   endif

   grav_prim_src_ofp(nid)=((src_den_prim_int(nid)-int_den(nid))/1027.0)*grav

   if (grav_prim_src_ofp(nid) .gt. 0.0) then
       occur_overflow_OFP(nid)=.true.
   else
       occur_overflow_OFP(nid)=.false.
   endif

   !------- all
   if (occur_overflow_OFP(nid)) then

       !--- still all
       hs_ofp(nid) = 2.0*hu_ofp(nid)/3.0
!
!--- Maximum geostrophic transport at the channel from the source [m^3/s]
!
      coriolis_source_ofp(nid) = 2.0*omega*sin(deg_to_rad*phi_ofp(nid) )
!
      trans_vol_ofp_src_raw(nid) = 0.5*grav_prim_src_ofp(nid)* hu_ofp(nid)**2  &
                                  / (coriolis_source_ofp(nid)+epsln)    
! 
!--- Maximum trasport
!      
      trans_vol_ofp_src(nid) = min(trans_vol_ofp_src_raw(nid), max_vol_trans_ofp)
!
!--- Maximum geostrophic velocity at the channel from the source [m/s]
!
      vel_src_ofp(nid) = trans_vol_ofp_src(nid) / (hs_ofp(nid)*ws_ofp(nid))
!
!---- Entrainment process based on Price & Baringer (1994)
!
!---- Define densities in entrainment region
!   ent_volT=0.0
!   ent_volS=0.0
!   ent_volP=0.0
!   ent_vol=0.0
!   ent_dat=0.0
!

      if (working_ent_comp_domains(nid)) then

          do k=ent_kst(nid),ent_ked(nid)
             do j=ent_jsc(nid),ent_jec(nid)
                do i=ent_isc(nid),ent_iec(nid)
                   ent_volT(nid)=ent_volT(nid)+tempo(i,j,k)*dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)
                   ent_volS(nid)=ent_volS(nid)+salto(i,j,k)*dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)
                   ent_vol(nid) =ent_vol(nid) + dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)
                enddo
             enddo
          enddo

          do n=1,num_prog_tracers
             sum_tmp=0.0
             do k=ent_kst(nid),ent_ked(nid)
                do j=ent_jsc(nid),ent_jec(nid)
                   do i=ent_isc(nid),ent_iec(nid)
                      sum_tmp=sum_tmp+traceo(i,j,k,n)*dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)
                   enddo
                enddo
             enddo
             ent_vol_trace(nid,n)=sum_tmp
          enddo

          do k=ent_kst(nid),ent_ked(nid)
             do j=ent_jsc(nid),ent_jec(nid)
                do i=ent_isc(nid),ent_iec(nid)
                   ent_volP(nid) = ent_volP(nid)+ presso(i,j,k)*dato(i,j)*tmsko(i,j,k)
                   ent_dat(nid)  = ent_dat(nid) + dato(i,j)*tmsko(i,j,k)
                enddo
             enddo
          enddo

      endif !(working_ent_comp_domains(nid))


  call mpp_sum(ent_volT(nid))
  call mpp_sum(ent_volS(nid))
  call mpp_sum(ent_volP(nid))
  call mpp_sum(ent_vol (nid))
  call mpp_sum(ent_dat (nid))
  call mpp_sum(ent_vol_trace(nid,1:num_prog_tracers), num_prog_tracers)

   !--- truncates at 4 bytes
   ent_volT(nid) = real( ent_volT(nid), kind=4 )
   ent_volS(nid) = real( ent_volS(nid), kind=4 )
   ent_volP(nid) = real( ent_volP(nid), kind=4 )
   ent_vol (nid) = real( ent_vol (nid), kind=4 )
   ent_dat (nid) = real( ent_dat (nid), kind=4 )

    do n=1,num_prog_tracers
      ent_vol_trace(nid,n)=real( ent_vol_trace(nid,n), kind=4 )
    enddo
!---

   if(ent_vol(nid) .le. 0.0 ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_overflow_OFP_mod (overflow_OFP): entrainment region has no sea cell')
   endif 

    !---  volume averaged T and S
    !---  for the pressure, regionally averaged P over the whole boxes 
    ent_volT(nid) = ent_volT(nid)/ent_vol(nid)
    ent_volS(nid) = ent_volS(nid)/ent_vol(nid)
    ent_volP(nid) = ent_volP(nid)/ent_dat(nid)
  
   do n=1,num_prog_tracers
     ent_vol_trace(nid,n)=ent_vol_trace(nid,n)/ent_vol(nid)
   enddo
!
!
!---
   ent_den(nid) = density(ent_volS(nid),ent_volT(nid),ent_volP(nid))
!
!--- define the averaged source density at the depth entrainment region
!
   src_den_prim(nid)=density(src_volS(nid),src_volT(nid),ent_volP(nid))
!  
!--- define reduce gravity
!
   grav_prim_ent_ofp(nid)=((src_den_prim(nid)-ent_den(nid))/1027.0)*grav
!
   if (debug_this_module .and. working_ent_comp_domains(nid)) then 
        write(stdoutunit,*) 'occur_overflow_OFP(nid) == .true. : nid', nid
        write(stdoutunit,*) 'hs_ofp(nid),coriolis_source_ofp(nid),trans_vol_ofp_src_raw(nid)' &
                            ,hs_ofp(nid),coriolis_source_ofp(nid),trans_vol_ofp_src_raw(nid)
        write(stdoutunit,*) 'vel_src_ofp(nid),trans_vol_ofp_src(nid)' &
                            ,vel_src_ofp(nid),trans_vol_ofp_src(nid)
        write(stdoutunit,*) 'ent_volT(nid),src_volS(nid),ent_volP(nid),ent_den(nid),grav_prim_ent_ofp(nid)' &
                            ,ent_volT(nid),src_volS(nid),ent_volP(nid),ent_den(nid),grav_prim_ent_ofp(nid)
   endif
!---- 
   if ((grav_prim_ent_ofp(nid) .le. 0.0) .or. (.not. do_entrainment_para_ofp)) then
!
       trans_vol_ofp_ent(nid)=0.0
       trans_mass_ofp_ent(nid)=0.0
       mix_theta_ofp(nid)=0.0
!
   else  !(grav_prim_ent_ofp(nid) .le. 0.0) : entranment process start
!
!---- compute flow speed at the shelf-slop break [m/s]
!
       vel_ssb_ofp(nid)=grav_prim_ent_ofp(nid)*alpha_ofp(nid)/(coriolis_source_ofp(nid)+epsln)
!
!--- checked for instability
!
      vel_ssb_ofp(nid) = min(max_ofp_speed, vel_ssb_ofp(nid))
!
!---- compute flow speed during spreading
!
      vel_ave_ofp(nid)=0.5*(vel_src_ofp(nid)+vel_ssb_ofp(nid))
!
!--- compute overflow thickness at the shelf/slope break
!
      coef_c_ofp(nid) = -(vel_src_ofp(nid)/vel_ssb_ofp(nid))*hs_ofp(nid)*hs_ofp(nid)
!
      coef_b_ofp(nid) = (1.0-vel_src_ofp(nid)/vel_ssb_ofp(nid))*hs_ofp(nid)  &
              +4.0*(cdbot_ofp(nid)*vel_ave_ofp(nid)*xssb_ofp(nid))  &
                  /(coriolis_source_ofp(nid)*ws_ofp(nid))
!
      coef_d_ofp(nid) = coef_b_ofp(nid)**2 - 4.0*coef_c_ofp(nid)
!
!----  chech coef_d_ofp>=0
    if(coef_d_ofp(nid) < 0.0 ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_overflow_OFP_mod (overflow_OFP): coef_d_ofp is less than 0.0')
    endif 
!----
!
    hssb_ofp(nid) = 0.5*(-coef_b_ofp(nid)+sqrt(coef_d_ofp(nid)))
!

!--- compute geostrophic Froude number
!--- check grav_prim_ent_ofp*hssb_ofp>=0
    if(hssb_ofp(nid) < 0.0 ) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_overflow_OFP_mod (overflow_OFP): hssb_ofp is less than 0.0')
    endif 
!----
!
      Fr_geo_ofp(nid)=vel_ssb_ofp(nid)/sqrt(grav_prim_ent_ofp(nid)*hssb_ofp(nid))
!     
!---- occurrence of mixing 
!---- compute entrainment mixing paramenter
!
     if (Fr_geo_ofp(nid) .ge. crit_Fr_geo_ofp) then     ! mixing occur
         mix_theta_ofp(nid)=1.0-1.0/Fr_geo_ofp(nid)**(2./3.)
     else ! (Fr_geo_ofp .ge. crit_Fr_geo_ofp)      ! no mixing
         mix_theta_ofp(nid)=0.0
     endif !(Fr_geo_ofp .ge. crit_Fr_geo_ofp)
!
!---- compute entrainment transport  [m^3/s], [kg/m^3 m^3/s]
!
           trans_vol_ofp_ent(nid)=trans_vol_ofp_src(nid)*mix_theta_ofp(nid)/(1.0-mix_theta_ofp(nid))
!
!---
     if (debug_this_module) then 
       write(stdoutunit,*) 'entranment process start'
       write(stdoutunit,*) 'vel_ave_ofp(nid),coef_c_ofp(nid),coef_b_ofp(nid),coef_d_ofp(nid),hssb_ofp(nid)' &
                           ,vel_ave_ofp(nid),coef_c_ofp(nid),coef_b_ofp(nid),coef_d_ofp(nid),hssb_ofp(nid)
       write(stdoutunit,*) 'Fr_geo_ofp(nid),mix_theta_ofp(nid),trans_vol_ofp_ent(nid)' &
                           ,Fr_geo_ofp(nid),mix_theta_ofp(nid),trans_vol_ofp_ent(nid)
     endif
!---
!
   endif !(grav_prim_ent_ofp .le. 0.0) : entrainment process
!
!---- compute product transport [m^3/s]
!
   trans_vol_ofp_prd(nid)  = trans_vol_ofp_src(nid)  + trans_vol_ofp_ent(nid)

!                   
!--- compute product Temperature and Salinity [degree C], [psu]
!--- with constraint of the volume conservation
!
   prd_volT(nid)=src_volT(nid)*(1.0-mix_theta_ofp(nid))+ent_volT(nid)*mix_theta_ofp(nid)
   prd_volS(nid)=src_volS(nid)*(1.0-mix_theta_ofp(nid))+ent_volS(nid)*mix_theta_ofp(nid)
!
   do n=1,num_prog_tracers
       prd_vol_trace(nid,n) = src_vol_trace(nid,n)*(1.0-mix_theta_ofp(nid))  &
                            + ent_vol_trace(nid,n)*mix_theta_ofp(nid)
   enddo

!           
!-- Volume fluxes [m^3/s] of source and entrainment and product region, are known.
!-- From these fluxes, mass and tracer fluxes are calculated, which should be conserved.
    if (debug_this_module) then 
      write(stdoutunit,*) 'trans_vol_ofp_prd(nid),prd_volT(nid),prd_volS(nid)' &
                          ,trans_vol_ofp_prd(nid),prd_volT(nid),prd_volS(nid)
      write(stdoutunit,*) 'trans_vol_ofp_prd(nid),prd_volT(nid),prd_volS(nid)' &
                          ,trans_vol_ofp_prd(nid),prd_volT(nid),prd_volS(nid)
      write(stdoutunit,*) 'src_volT(nid),ent_volT(nid),mix_theta_ofp(nid),src_volS(nid),ent_volS(nid)' &
                          ,src_volT(nid),ent_volT(nid),mix_theta_ofp(nid),src_volS(nid),ent_volS(nid)
    endif
!
!--- at the source region, volume flux [m^3/s] / total volume [m^3]
!
   volf_per_tvol_src(nid)=trans_vol_ofp_src(nid) / src_vol(nid)  ![1/s]
!   
!   src_flux_Mass(nid) = 0.0  
!
     if (working_src_comp_domains(nid)) then   
!--- total mass flux from the source region  
   do k=src_kst(nid),src_ked(nid)
    do j=src_jsc(nid),src_jec(nid)
     do i=src_isc(nid),src_iec(nid)
       src_flux_Mass(nid)=src_flux_Mass(nid)     &
                         +volf_per_tvol_src(nid)*rho_dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)  ![kg/m^3 m^3/s]
     enddo
    enddo
   enddo
!
     endif !(working_src_comp_domains(nid))

!--- mpp_sum

     call mpp_sum ( src_flux_Mass(nid) )  

!--- truncates at 4 bytes

      src_flux_Mass(nid) = real ( src_flux_Mass(nid), kind=4 )

 if ( working_src_comp_domains(nid) ) then
   if (do_mass_ofp) then

!--- uniformly decrease mass over the source region [-volf_per_tvol_src*vol(i,j,k)/dat(i,j)]
     do k=src_kst(nid),src_ked(nid)
      do j=src_jsc(nid),src_jec(nid)
       do i=src_isc(nid),src_iec(nid)
           Thickness%mass_source(i,j,k) = Thickness%mass_source(i,j,k)  &
                                         -volf_per_tvol_src(nid)*rho_dzto(i,j,k)*tmsko(i,j,k)     ! [kg/m^3 m/s]
       enddo
      enddo
     enddo
!---
    if (debug_this_module .and. working_src_comp_domains(nid)) then 
      write(stdoutunit,*) 'volf_per_tvol_src(nid),src_flux_Mass(nid)',volf_per_tvol_src(nid),src_flux_Mass(nid)
    endif
!---
!
!--- total trace flux from the source region  
!--- in two level time step, taum1=tau
    do n=1,num_prog_tracers
     do k=src_kst(nid),src_ked(nid)
      do j=src_jsc(nid),src_jec(nid)
       do i=src_isc(nid),src_iec(nid)
            src_flux_trace(nid,n)= src_flux_trace(nid,n)  &
                               +T_prog(n)%field(i,j,k,taum1)  &
                               *volf_per_tvol_src(nid)*rho_dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)   ![trace kg/m^3 m^3/s]
!
            T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k)  &
                                          -T_prog(n)%field(i,j,k,taum1)  &
                                          *volf_per_tvol_src(nid)*rho_dzto(i,j,k)*tmsko(i,j,k)  ![trace kg/m^3 m/s]
!
        enddo
       enddo
      enddo
     enddo !n=1,num_prog_tracers
!
!
   endif !(do_mass_ofp)
 endif ! ( working_src_comp_domains(nid) )
!----
!--- at the entrainment region, volume flux [m^3/s] / total volume [m^3]
!
   volf_per_tvol_ent(nid)=trans_vol_ofp_ent(nid) / ent_vol(nid)  ![1/s]
!   
!   ent_flux_Mass(nid) = 0.0  
!
   if ( working_ent_comp_domains(nid) ) then
!----- total mass flux from the entrainment region  
  do k=ent_kst(nid),ent_ked(nid)
   do j=ent_jsc(nid),ent_jec(nid)
    do i=ent_isc(nid),ent_iec(nid)
          ent_flux_Mass(nid)=ent_flux_Mass(nid)     &
                            + volf_per_tvol_ent(nid)*rho_dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)  ![kg/m^3 m^3/s]
    enddo
   enddo
  enddo
!
   endif !( working_ent_comp_domains(nid) )

!--- mpp_sum
     
     call mpp_sum (ent_flux_Mass(nid)) 

!--- truncates at 4 bytes

     ent_flux_Mass(nid) = real( ent_flux_Mass(nid), kind=4 )

 if ( working_ent_comp_domains(nid) ) then
   if (do_mass_ofp) then
!
!--- uniformly decrease mass over the entrainment region 
    do k=ent_kst(nid),ent_ked(nid)
     do j=ent_jsc(nid),ent_jec(nid)
      do i=ent_isc(nid),ent_iec(nid)
         Thickness%mass_source(i,j,k) = Thickness%mass_source(i,j,k)  &
                                       -volf_per_tvol_ent(nid)*rho_dzto(i,j,k)*tmsko(i,j,k)     ! [kg/m^3 m/s]
      enddo
     enddo
    enddo
!
   if (debug_this_module) then 
      write(stdoutunit,*) 'volf_per_tvol_ent(nid),src_flux_Mass(nid)',volf_per_tvol_ent(nid),src_flux_Mass(nid)
   endif
!
!--- total trace flux from the entrainment region  
    do n=1,num_prog_tracers
     do k=ent_kst(nid),ent_ked(nid)
      do j=ent_jsc(nid),ent_jec(nid)
       do i=ent_isc(nid),ent_iec(nid)
         ent_flux_trace(nid,n)= ent_flux_trace(nid,n) &
                           +T_prog(n)%field(i,j,k,taum1)*volf_per_tvol_ent(nid) &
                         *rho_dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)                          ![trace kg/m^3 m^3/s]
!
          T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) &
                                        -T_prog(n)%field(i,j,k,taum1) &
                                        *volf_per_tvol_ent(nid)*rho_dzto(i,j,k)*tmsko(i,j,k)   ![trace kg/m^3 m/s]
!
       enddo
      enddo
     enddo
    enddo !n=1,num_prog_tracers
!----
   endif !(do_mass_ofp)
 endif !( working_ent_comp_domains(nid) )
!
!
!--- At the product region,
!--- to check the regionally averaged density between ambient density of product lines
!--- and the product water density,
!--- and to decide the injection depth and the product line. 
!
     do m=1,10
       ijk_check=prd_kst(nid,m)*prd_ked(nid,m)*prd_jst(nid,m)*prd_jed(nid,m) &
                *prd_ist(nid,m)*prd_ied(nid,m)
       if (ijk_check /=0) prd_line_num(nid)=m
     enddo

!
!--- Calculate ambient averaged density at each line
!--- prd_ked can be kmt, that is, partial cell. So prd_kst (prd_ked-1) is considered.
         
   do m=1,prd_line_num(nid)

       if (working_prd_comp_domains(nid,m)) then

    do k=prd_kst(nid,m),prd_ked(nid,m)
     do j=prd_jsc(nid,m),prd_jec(nid,m)
      do i=prd_isc(nid,m),prd_iec(nid,m)
         prd_vol(nid,m)  = prd_vol(nid,m)    &
                          +dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)
         prd_amb_T(nid,m)= prd_amb_T(nid,m)  &
                          +tempo(i,j,k)*dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)
         prd_amb_S(nid,m)= prd_amb_S(nid,m)  &
                          +salto(i,j,k)*dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)

   if (debug_this_module) then 
      write(stdoutunit,*) 'CP-030411 : m,i,j,k,tempo(i,j,k),salto(i,j,k)',m,i,j,k,tempo(i,j,k),salto(i,j,k)
   endif

     enddo
    enddo
   enddo
!
   do k=prd_kst(nid,m),prd_ked(nid,m)
    do j=prd_jsc(nid,m),prd_jec(nid,m)
     do i=prd_isc(nid,m),prd_iec(nid,m)
       prd_amb_P(nid,m) = prd_amb_P(nid,m) + presso(i,j,k)*dato(i,j)*tmsko(i,j,k)
       prd_dat(nid,m)   = prd_dat(nid,m) + dato(i,j)*tmsko(i,j,k)

       if (debug_this_module) then 
        write(stdoutunit,*) 'PT4.2: m,prd_kst(nid,m), Dens%pressure_at_depth(i,j,prd_kst(nid,m)),Grd%dat(i,j)' &
                            ,m,prd_kst(nid,m), Dens%pressure_at_depth(i,j,prd_kst(nid,m)),Grd%dat(i,j)
       endif

     enddo
    enddo
   enddo

      endif !(working_prd_comp_domains(nid,m))

!--- mpp_sum

     call mpp_sum (prd_amb_T(nid,m))
     call mpp_sum (prd_amb_S(nid,m))
     call mpp_sum (prd_amb_P(nid,m))
     call mpp_sum (prd_dat(nid,m))
     call mpp_sum (prd_vol(nid,m))

!--- truncates at 4 bytes
   prd_amb_T(nid,m) = real( prd_amb_T(nid,m), kind=4 )
   prd_amb_S(nid,m) = real( prd_amb_S(nid,m), kind=4 )
   prd_amb_P(nid,m) = real( prd_amb_P(nid,m), kind=4 )
   prd_vol(nid,m) = real( prd_vol(nid,m), kind=4 )
   prd_dat(nid,m) = real( prd_dat(nid,m), kind=4 )

!
  if(prd_vol(nid,m) .le. 0.0 ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_overflow_OFP_mod (overflow_OFP): interior region has no sea cell')
  endif 
! 
   prd_amb_T(nid,m)= prd_amb_T(nid,m)/(epsln+prd_vol(nid,m))
   prd_amb_S(nid,m)= prd_amb_S(nid,m)/(epsln+prd_vol(nid,m))
   prd_amb_P(nid,m)= prd_amb_P(nid,m)/(epsln+prd_dat(nid,m))
!
   if (debug_this_module) then 
      write(stdoutunit,*) 'CP-072810 : m,prd_amb_T(nid,m),prd_amb_S(nid,m)',m,prd_amb_T(nid,m),prd_amb_S(nid,m)
   endif
! 
!!--- define ambient averaged density at product region
!
    prd_S=prd_amb_S(nid,m)
    prd_T=prd_amb_T(nid,m)
    prd_P=prd_amb_P(nid,m)
    prd_amb_rho(nid,m)=density(prd_S,prd_T,prd_P)  
!
    if (debug_this_module) then 
       write(stdoutunit,*) ' m, prd_amb_T(nid,m), prd_amb_S(nid,m),prd_amb_P(nid,m),prd_amb_rho(nid,m)' &
                           , m, prd_amb_T(nid,m), prd_amb_S(nid,m),prd_amb_P(nid,m),prd_amb_rho(nid,m)
    endif
!
  enddo !m=1,prd_line_num
!
!
!!--- define the injection line from the check of the density
!!--- it is assumed that the depth at m=1 is the shallowest line
!
   prd_inject_line_num(nid)=1
!
   do m=1,prd_line_num(nid)-1
      prd_vol_rho(nid) = density(prd_volS(nid), prd_volT(nid),prd_amb_P(nid,m)) 
      if ( prd_vol_rho(nid) > prd_amb_rho(nid,m) ) then
         prd_inject_line_num(nid) = m+1  
      endif    
   enddo !m=1,prd_line_num
!
   if (debug_this_module) then 
      write(stdoutunit,*) 'CP-072810: prd_inject_line_num(nid)',prd_inject_line_num(nid)
   endif
!
   if (debug_this_module) then 
      write(stdoutunit,*) 'CP-072810: trans_vol_ofp_prd(nid)',trans_vol_ofp_prd(nid)
   endif
!
!--- total mass flux of prd from src and ent
!
   prd_flux_Mass(nid)=src_flux_Mass(nid)+ent_flux_Mass(nid)                             ![kg/m^3 m^3/s]
!
!--- total trace flux of prd from src and ent
!
   do n=1,num_prog_tracers
      prd_flux_trace(nid,n)= src_flux_trace(nid,n)+ent_flux_trace(nid,n)              ![trace kg/m^3 m^3/s]
   enddo !n=1,num_prog_tracers
!
!   prd_inj_vol(nid)=0.0

      if ( working_prd_comp_domains(nid,prd_inject_line_num(nid)) ) then

   do k=prd_kst(nid,prd_inject_line_num(nid)),prd_ked(nid,prd_inject_line_num(nid))
    do j=prd_jsc(nid,prd_inject_line_num(nid)),prd_jec(nid,prd_inject_line_num(nid))
     do i=prd_isc(nid,prd_inject_line_num(nid)),prd_iec(nid,prd_inject_line_num(nid))
        prd_inj_vol(nid) = prd_inj_vol(nid) + dzto(i,j,k)*dato(i,j)*tmsko(i,j,k)                   ! [m^3]
     enddo
    enddo
   enddo
 
      depth_prd_OFP(nid)=Thickness%depth_zt(prd_isc(nid,prd_inject_line_num(nid)) &
                                           ,prd_jsc(nid,prd_inject_line_num(nid)) &
                                           ,prd_kst(nid,prd_inject_line_num(nid)))
      endif !( working_prd_comp_domains(nid,prd_inject_line_num(nid)) )
!--- mpp_sum

     call mpp_sum ( prd_inj_vol(nid) )      
     call mpp_max ( depth_prd_OFP(nid) )      

!--- truncates at 4 bytes

     prd_inj_vol(nid) = real( prd_inj_vol(nid), kind=4 )
     
!---
  if(prd_inj_vol(nid) .le. 0.0 ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_overflow_OFP_mod (overflow_OFP): prd_inj_vol(nid) .le. 0.0')
  endif 

!
!!--- volume flux [m^3/s] / total volume [m^3]
!
     volf_per_tvol_prd(nid) = trans_vol_ofp_prd(nid) / prd_inj_vol(nid)      ! [1/s]
!
!
!--- uniformly increase mass (volume weighted) over the product injection region
!--- defined as (prd_flux_Mass / prd_inj_vol)*vol(i,j,k)/dat(i,j)
!
  if ( working_prd_comp_domains(nid,prd_inject_line_num(nid)) ) then
    if (do_mass_ofp) then
!
      do k=prd_kst(nid,prd_inject_line_num(nid)),prd_ked(nid,prd_inject_line_num(nid))
       do j=prd_jsc(nid,prd_inject_line_num(nid)),prd_jec(nid,prd_inject_line_num(nid))
        do i=prd_isc(nid,prd_inject_line_num(nid)),prd_iec(nid,prd_inject_line_num(nid))
           Thickness%mass_source(i,j,k) = Thickness%mass_source(i,j,k)  &
                                         +(prd_flux_Mass(nid)/prd_inj_vol(nid)) *dzto(i,j,k)*tmsko(i,j,k) ! [kg/m^3 m/s]
        enddo
       enddo
      enddo
!---------------------------------     
     if (prd_vol_trace(nid,index_temp) /=  prd_volT(nid) ) then 
       write (errMsg, FMT='("==>Error from ocean_overflow_mod (ocean_overflow): prd_vol_trace(nid,index_temp) = ",&
            &F20.10," is not same as prd_volT(nid) = ",F20.10)') prd_vol_trace(nid,index_temp), prd_volT(nid)
       call mpp_error(FATAL, trim(errMsg))
     endif 
     if (prd_vol_trace(nid,index_salt) /=  prd_volS(nid) ) then 
       write (errMsg, FMT='("==>Error from ocean_overflow_mod (ocean_overflow): prd_vol_trace(nid,index_salt) = ",&
            &F20.10," is not same as prd_volS(nid) = ",F20.10)') prd_vol_trace(nid,index_salt), prd_volS(nid)
       call mpp_error(FATAL, trim(errMsg))
     endif 
!----------------------------
     do n=1,num_prog_tracers
      do k=prd_kst(nid,prd_inject_line_num(nid)),prd_ked(nid,prd_inject_line_num(nid))
       do j=prd_jsc(nid,prd_inject_line_num(nid)),prd_jec(nid,prd_inject_line_num(nid))
        do i=prd_isc(nid,prd_inject_line_num(nid)),prd_iec(nid,prd_inject_line_num(nid))
           T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) &
                                         +prd_vol_trace(nid,n)*(prd_flux_Mass(nid)  &
                                         /prd_inj_vol(nid)) * dzto(i,j,k) *tmsko(i,j,k)                ![trace kg/m^3 m/s]
        enddo
       enddo
      enddo
     enddo

!

    else !(do_mass_ofp)
!
     do n=1,num_prog_tracers
      do k=prd_kst(nid,prd_inject_line_num(nid)),prd_ked(nid,prd_inject_line_num(nid))
       do j=prd_jsc(nid,prd_inject_line_num(nid)),prd_jec(nid,prd_inject_line_num(nid))
        do i=prd_isc(nid,prd_inject_line_num(nid)),prd_iec(nid,prd_inject_line_num(nid))
           T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) &
                                         +(prd_vol_trace(nid,n)-T_prog(n)%field(i,j,k,tau))  &
                                         *(prd_flux_Mass(nid)/prd_inj_vol(nid)) *dzto(i,j,k)*tmsko(i,j,k) ![trace kg/m^3 m/s]
        enddo
       enddo
      enddo
     enddo 

    endif !(do_mass_ofp)
  endif !( working_prd_comp_domains(nid,prd_inject_line_num) )

!---
!
      if (debug_this_module) then 
         write(stdoutunit,*) '999: pe, nid, src_volT(nid),src_volS(nid),trans_vol_ofp_src(nid)' &
                         ,pe, nid, src_volT(nid),src_volS(nid),trans_vol_ofp_src(nid)
      if (working_comp_domains_all) then 
        write(stdoutunit,'(2i5,12e15.7)') nid, pe, src_volT(nid),src_volS(nid),trans_vol_ofp_src(nid) &
                     , int_volT(nid),int_volS(nid)                        &
                     , ent_volT(nid),ent_volS(nid),trans_vol_ofp_ent(nid) &
                     , prd_volT(nid),prd_volS(nid),trans_vol_ofp_prd(nid) &
                     , depth_prd_OFP(nid)
      endif
      endif
!------------------------- diagnostics
!--- Averaged temperature, salinity, transport and depth of of each region are save 
!--- as the ascii data. In standard script, the output file(.out) is saved
!--- in the directory of the ascii.
!--- Number in the file name (e.g. *01.out, *02.out ...) means the number of the OFP set
!--- (e.g. 01: Denmark Strait, 02: Faroe Bank Channel ...).
!--- diag_step is defined at the name list
!
     if (diag_step > 0 .and. pe==root_pe) then 
      if(mod(Time%itt,diag_step) == 0) then
        write(stdoutunit,'(24x,a)') "Ocean Overflow Parametrization (OFP) summary:"
        write(stdoutunit,'(2i5,12e15.7)') nid, pe, src_volT(nid),src_volS(nid),trans_vol_ofp_src(nid) &
                   , int_volT(nid),int_volS(nid)                        &
                   , ent_volT(nid),ent_volS(nid),trans_vol_ofp_ent(nid) &
                   , prd_volT(nid),prd_volS(nid),trans_vol_ofp_prd(nid) &
                   , depth_prd_OFP(nid)
      endif  
     endif
!---
!
 endif !(working_comp_domains(nid))  
!-------------------------------------------
       if (debug_this_module) then 
         write(stdoutunit,*) '999: pe, nid, src_volT(nid),src_volS(nid),trans_vol_ofp_src(nid)' &
                             ,pe, nid, src_volT(nid),src_volS(nid),trans_vol_ofp_src(nid)
       endif
!------------------------------------------------
enddo ! nid=1,n_OFP
!------------------------------------------------
        call mpp_clock_begin(id_OFP_update_2)       
    do n = 1, num_prog_tracers
       call mpp_update_domains(T_prog(n)%th_tendency, Dom%domain2d, complete=(n==num_prog_tracers))
    enddo
    
        call mpp_clock_end(id_OFP_update_2)  
     
        call mpp_clock_begin(id_OFP_update_3)
        
        call mpp_update_domains(Thickness%mass_source, Dom%domain2d)
        
        call mpp_clock_end(id_OFP_update_3)

        call watermass_diag(Time, Dens, T_prog, index_temp, index_salt, wrk1_v)                

        call mpp_clock_end(id_OFP_main)
        
!-------------------------------------------------
 if(id_OFP_n1_src_temp > 0)  used = send_data(id_OFP_n1_src_temp,src_volT(1),Time%model_time)        
 if(id_OFP_n1_src_salt > 0)  used = send_data(id_OFP_n1_src_salt,src_volS(1),Time%model_time)         
 if(id_OFP_n1_src_trans > 0) used = send_data(id_OFP_n1_src_trans,trans_vol_ofp_src(1),Time%model_time)         
 if(id_OFP_n1_int_temp > 0)  used = send_data(id_OFP_n1_int_temp,int_volT(1),Time%model_time)         
 if(id_OFP_n1_int_salt > 0)  used = send_data(id_OFP_n1_int_salt,int_volS(1),Time%model_time)         
 if(id_OFP_n1_ent_temp > 0)  used = send_data(id_OFP_n1_ent_temp,ent_volT(1),Time%model_time)         
 if(id_OFP_n1_ent_salt > 0)  used = send_data(id_OFP_n1_ent_salt,ent_volS(1),Time%model_time)         
 if(id_OFP_n1_ent_trans > 0) used = send_data(id_OFP_n1_ent_trans,trans_vol_ofp_ent(1),Time%model_time)         
 if(id_OFP_n1_prd_temp > 0)  used = send_data(id_OFP_n1_prd_temp,prd_volT(1),Time%model_time)         
 if(id_OFP_n1_prd_salt > 0)  used = send_data(id_OFP_n1_prd_salt,prd_volS(1),Time%model_time)         
 if(id_OFP_n1_prd_trans > 0) used = send_data(id_OFP_n1_prd_trans,trans_vol_ofp_prd(1),Time%model_time)         
 if(id_OFP_n1_prd_depth > 0) used = send_data(id_OFP_n1_prd_depth,depth_prd_OFP(1),Time%model_time) 
        
 if(id_OFP_n2_src_temp > 0)  used = send_data(id_OFP_n2_src_temp,src_volT(2),Time%model_time)         
 if(id_OFP_n2_src_salt > 0)  used = send_data(id_OFP_n2_src_salt,src_volS(2),Time%model_time)         
 if(id_OFP_n2_src_trans > 0) used = send_data(id_OFP_n2_src_trans,trans_vol_ofp_src(2),Time%model_time)         
 if(id_OFP_n2_int_temp > 0)  used = send_data(id_OFP_n2_int_temp,int_volT(2),Time%model_time)         
 if(id_OFP_n2_int_salt > 0)  used = send_data(id_OFP_n2_int_salt,int_volS(2),Time%model_time)         
 if(id_OFP_n2_ent_temp > 0)  used = send_data(id_OFP_n2_ent_temp,ent_volT(2),Time%model_time)         
 if(id_OFP_n2_ent_salt > 0)  used = send_data(id_OFP_n2_ent_salt,ent_volS(2),Time%model_time)         
 if(id_OFP_n2_ent_trans > 0) used = send_data(id_OFP_n2_ent_trans,trans_vol_ofp_ent(2),Time%model_time)         
 if(id_OFP_n2_prd_temp > 0)  used = send_data(id_OFP_n2_prd_temp,prd_volT(2),Time%model_time)         
 if(id_OFP_n2_prd_salt > 0)  used = send_data(id_OFP_n2_prd_salt,prd_volS(2),Time%model_time)         
 if(id_OFP_n2_prd_trans > 0) used = send_data(id_OFP_n2_prd_trans,trans_vol_ofp_prd(2),Time%model_time)         
 if(id_OFP_n2_prd_depth > 0) used = send_data(id_OFP_n2_prd_depth,depth_prd_OFP(2),Time%model_time) 
!-------------------------------------------------
           
end subroutine overflow_OFP
! </SUBROUTINE> NAME="overflow_OFP"


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

  
  id_neut_rho_overofp = register_diag_field ('ocean_model', 'neut_rho_overofp',&
       Grd%tracer_axes(1:3), Time%model_time,                                  &
       'update of locally referenced potrho from overflow_ofp scheme',         &
       '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_overofp > 0) compute_watermass_diag = .true. 

  id_neut_rho_overofp_on_nrho = register_diag_field ('ocean_model',                     &
   'neut_rho_overofp_on_nrho',Dens%neutralrho_axes(1:3), Time%model_time,               &
       'update of locally ref potrho from overflow_ofp as binned to neutral rho layers',&
       '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_overofp_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_overofp = register_diag_field ('ocean_model',  &
   'wdian_rho_overofp', Grd%tracer_axes(1:3), Time%model_time,&
       'dianeutral mass transport due to overflow_ofp',       &
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_overofp > 0) compute_watermass_diag = .true. 

  id_wdian_rho_overofp_on_nrho = register_diag_field ('ocean_model',                   &
   'wdian_rho_overofp_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,            &
       'dianeutral mass transport due to overflow_ofp as binned to neutral rho layers',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_overofp_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_overofp_on_nrho = register_diag_field ('ocean_model',             &
   'tform_rho_overofp_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,      &
       'watermass transform due to overflow_ofp as binned to neutral rho layers',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_overofp_on_nrho > 0) compute_watermass_diag = .true. 


  ! temperature effects  
  id_neut_temp_overofp = register_diag_field ('ocean_model', 'neut_temp_overofp',&
   Grd%tracer_axes(1:3), Time%model_time,                                        &
   'temp related update of locally referenced potrho from overflow_ofp scheme',  &
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_overofp > 0) compute_watermass_diag = .true. 

  id_neut_temp_overofp_on_nrho = register_diag_field ('ocean_model',                          &
   'neut_temp_overofp_on_nrho',Dens%neutralrho_axes(1:3), Time%model_time,                    &
   'temp related update of locally ref potrho from overflow_ofp binned to neutral rho layers',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_overofp_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_temp_overofp = register_diag_field ('ocean_model',   &
   'wdian_temp_overofp', Grd%tracer_axes(1:3), Time%model_time, &
   'temp related dianeutral mass transport due to overflow_ofp',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_overofp > 0) compute_watermass_diag = .true. 

  id_wdian_temp_overofp_on_nrho = register_diag_field ('ocean_model',                        &
   'wdian_temp_overofp_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                 &
   'temp related dianeutral mass transport due to overflow_ofp binned to neutral rho layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_overofp_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_temp_overofp_on_nrho = register_diag_field ('ocean_model',                     &
   'tform_temp_overofp_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,              &
   'temp related watermass transform due to overflow_ofp as binned to neutral rho layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_temp_overofp_on_nrho > 0) compute_watermass_diag = .true. 


  ! salinity effects  
  id_neut_salt_overofp = register_diag_field ('ocean_model', 'neut_salt_overofp',&
   Grd%tracer_axes(1:3), Time%model_time,                                        &
   'salt related update of locally referenced potrho from overflow_ofp scheme',  &
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_overofp > 0) compute_watermass_diag = .true. 

  id_neut_salt_overofp_on_nrho = register_diag_field ('ocean_model',                          &
   'neut_salt_overofp_on_nrho',Dens%neutralrho_axes(1:3), Time%model_time,                    &
   'salt related update of locally ref potrho from overflow_ofp binned to neutral rho layers',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_overofp_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_salt_overofp = register_diag_field ('ocean_model',   &
   'wdian_salt_overofp', Grd%tracer_axes(1:3), Time%model_time, &
   'salt related dianeutral mass transport due to overflow_ofp',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_overofp > 0) compute_watermass_diag = .true. 

  id_wdian_salt_overofp_on_nrho = register_diag_field ('ocean_model',                        &
   'wdian_salt_overofp_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                 &
   'salt related dianeutral mass transport due to overflow_ofp binned to neutral rho layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_overofp_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_salt_overofp_on_nrho = register_diag_field ('ocean_model',                     &
   'tform_salt_overofp_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,              &
   'salt related watermass transform due to overflow_ofp as binned to neutral rho layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_salt_overofp_on_nrho > 0) compute_watermass_diag = .true. 


  if(compute_watermass_diag) then 
    write(stdoutunit,'(/a/)') &
    '==>Note: running ocean_overflow_OFP_mod w/ compute_watermass_diag=.true.'  
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
subroutine watermass_diag(Time, Dens, T_prog, index_temp, index_salt, wrk1v)

  type(ocean_time_type),           intent(in)    :: Time
  type(ocean_density_type),        intent(in)    :: Dens
  type(ocean_prog_tracer_type),    intent(in)    :: T_prog(:)
  integer,                         intent(in)    :: index_temp
  integer,                         intent(in)    :: index_salt
  real, dimension(isd:,jsd:,:,:), intent(in) :: wrk1v

  integer :: i,j,k
  real    :: tmp1, tmp2

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
           tmp1 = T_prog(index_temp)%th_tendency(i,j,k) - wrk1v(i,j,k,1)
           tmp2 = T_prog(index_salt)%th_tendency(i,j,k) - wrk1v(i,j,k,2)
           wrk1(i,j,k) = Grd%tmask(i,j,k)*(Dens%drhodT(i,j,k)*tmp1 + Dens%drhodS(i,j,k)*tmp2)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)  
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_overofp,wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_overofp, wrk3(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_overofp_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_overofp_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_overofp_on_nrho, wrk4)

  ! temperature effects on rho 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec

           tmp1        = T_prog(index_temp)%th_tendency(i,j,k) - wrk1v(i,j,k,1)
           wrk1(i,j,k) = Grd%tmask(i,j,k)*Dens%drhodT(i,j,k)*tmp1
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)  
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_temp_overofp,wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_overofp, wrk3(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_overofp_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_overofp_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_overofp_on_nrho, wrk4)

  ! salinity effects on rho 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tmp2        = T_prog(index_salt)%th_tendency(i,j,k) - wrk1v(i,j,k,2)
           wrk1(i,j,k) = Grd%tmask(i,j,k)*Dens%drhodT(i,j,k)*tmp2
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)  
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_salt_overofp,wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_overofp, wrk3(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_overofp_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_overofp_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_overofp_on_nrho, wrk4)

end subroutine watermass_diag
! </SUBROUTINE>  NAME="watermass_diag"


end module ocean_overflow_OFP_mod
!!

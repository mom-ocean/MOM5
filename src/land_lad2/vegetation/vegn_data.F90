module vegn_data_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : &
     write_version_number, file_exist, check_nml_error, &
     close_file, stdlog, stdout

use land_constants_mod, only : NBANDS, BAND_VIS, BAND_NIR
use land_tile_selectors_mod, only : &
     tile_selector_type, SEL_VEGN, register_tile_selector
use table_printer_mod

implicit none
private

! ==== public interfaces =====================================================
! ---- public constants
integer, public, parameter :: LU_SEL_TAG = 1 ! tag for the land use selectors
integer, public, parameter :: SP_SEL_TAG = 2 ! tag for the species selectors
integer, public, parameter :: NG_SEL_TAG = 3 ! tag for natural grass selector
  ! by "natural" it means non-human-maintained, so secondary vegetation
  ! grassland will be included. 

integer, public, parameter :: NSPECIES = 5, & ! number of species
 SP_C4GRASS   = 0, & ! c4 grass
 SP_C3GRASS   = 1, & ! c3 grass
 SP_TEMPDEC   = 2, & ! temperate deciduous
 SP_TROPICAL  = 3, & ! non-grass tropical
 SP_EVERGR    = 4    ! non-grass evergreen
character(len=12), parameter :: species_name(0:NSPECIES-1) = &
    (/'c4grass  ',  'c3grass  ' ,  'tempdec  ', 'tropical ','evergreen'/)
character(len=32), parameter :: species_longname(0:NSPECIES-1) = &
    (/'c4 grass                 ', 'c3 grass                 ',  'temperate deciduous trees',&
      'tropical trees           ', 'evergreen trees          '/)

integer, public, parameter :: n_dim_vegn_types = 9
integer, public, parameter :: MSPECIES = NSPECIES+n_dim_vegn_types-1
 
integer, public, parameter :: NCMPT = 6, & ! number of carbon compartments
 CMPT_REPRO   = 1, & ! 
 CMPT_SAPWOOD = 2, & ! sapwood compartment
 CMPT_LEAF    = 3, & ! leaf compartment
 CMPT_ROOT    = 4, & ! fine root compartment
 CMPT_VLEAF   = 5, & ! virtual leaves compartment (labile store)
 CMPT_WOOD    = 6    ! structural wood compartment

integer, public, parameter :: & ! physiology types
 PT_C3        = 0, &
 PT_C4        = 1

integer, public, parameter :: & ! phenology type
 PHEN_DECIDIOUS = 0, &
 PHEN_EVERGREEN = 1

integer, public, parameter :: & ! status of leaves
 LEAF_ON      = 0, &  ! leaves are displayed
 LEAF_OFF     = 5     ! leaves are dropped

integer, public, parameter :: & ! land use types
 N_LU_TYPES = 4, & ! number of different land use types
 LU_PAST    = 1, & ! pasture
 LU_CROP    = 2, & ! crops
 LU_NTRL    = 3, & ! natural vegetation
 LU_SCND    = 4    ! secondary vegetation
character(len=4), public, parameter  :: &
     landuse_name (N_LU_TYPES) = (/ 'past','crop','ntrl','scnd'/)
character(len=32), public, parameter :: &
     landuse_longname (N_LU_TYPES) = (/ 'pasture  ', 'crop     ', 'natural  ', 'secondary' /)

integer, public, parameter :: & ! harvesing pools paraneters
 N_HARV_POOLS        = 6, & ! number of harvesting pools
 HARV_POOL_PAST      = 1, & 
 HARV_POOL_CROP      = 2, &
 HARV_POOL_CLEARED   = 3, &
 HARV_POOL_WOOD_FAST = 4, &
 HARV_POOL_WOOD_MED  = 5, &
 HARV_POOL_WOOD_SLOW = 6
character(len=9), public :: HARV_POOL_NAMES(N_HARV_POOLS)
data HARV_POOL_NAMES &
 / 'past', 'crop', 'cleared', 'wood_fast', 'wood_med', 'wood_slow' /

real, public, parameter :: C2B = 2.0  ! carbon to biomass conversion factor

real, public, parameter :: BSEED = 5e-5 ! seed density for supply/demand calculations, kg C/m2 
! ---- public types
public :: spec_data_type

! ---- public data
public :: &
    vegn_to_use,  input_cover_types, &
    mcv_min, mcv_lai, &
    use_bucket, use_mcm_masking, vegn_index_constant, &
    critical_root_density, &
    ! vegetation data, imported from LM3V
    spdata, &
    min_cosz, &
    agf_bs, K1,K2, fsc_liv, fsc_wood, &
    tau_drip_l, tau_drip_s, & ! canopy water and snow residence times, for drip calculations
    GR_factor, tg_c3_thresh, tg_c4_thresh, &
    fsc_pool_spending_time, ssc_pool_spending_time, harvest_spending_time, &
    l_fract, T_transp_min, soil_carbon_depth_scale, &
    cold_month_threshold, scnd_biomass_bins, &
    phen_ev1, phen_ev2, cmc_eps

! ---- public subroutine
public :: read_vegn_data_namelist
! ==== end of public interfaces ==============================================

! ==== constants =============================================================
character(len=*), parameter   :: &
     version     = '$Id: vegn_data.F90,v 20.0 2013/12/13 23:31:06 fms Exp $', &
     tagname     = '$Name: tikal $', &
     module_name = 'vegn_data_mod'
real, parameter :: TWOTHIRDS  = 2.0/3.0


! ==== types ================================================================
type spec_data_type
  real    :: treefall_disturbance_rate;
  logical :: mortality_kills_balive ! if true, then bl, blv, and br are affected by natural mortality
  integer :: pt           ! photosynthetic physiology of species

  real    :: c1 ! unitless, coefficient for living biomass allocation
  real    :: c2 ! 1/m, coefficient for living biomass allocation
  real    :: c3 ! unitless, coefficient for calculation of sapwood biomass 
                ! fraction times sapwood retirement rate

  real    :: alpha(NCMPT) ! decay rates of plant carbon pools, 1/yr
  real    :: beta (NCMPT) ! respiration rates of plant carbon pools
  
  real    :: dfr          ! fine root diameter ? or parameter relating diameter of fine roots to resistance
  ! the following two parameters are used in the Darcy-law calculations of water supply
  real    :: srl  ! specific root length, m/(kg C)
  real    :: root_r       ! radius of the fine roots, m
  real    :: root_perm    ! fine root membrane permeability per unit area, kg/(m3 s)
!!$  real    :: ltrans       ! leaf translocation fraction
!!$  real    :: rtrans       ! fine root translocation fraction

  real    :: specific_leaf_area ! cm2/(g biomass)
  real    :: leaf_size    ! characteristic leaf size
  real    :: leaf_life_span ! months
  
  real    :: alpha_phot   ! photosynthesis efficiency
  real    :: m_cond       ! factor of stomatal conductance
  real    :: Vmax         ! max rubisco rate
  real    :: gamma_resp
  real    :: wet_leaf_dreg ! wet leaf photosynthesis down-regulation
  real    :: leaf_age_onset, leaf_age_tau

  ! radiation parameters for 2 bands, VIS and NIR
  real    :: leaf_refl (NBANDS) ! reflectance of leaf
  real    :: leaf_tran (NBANDS) ! transmittance of leaf
  real    :: leaf_emis          ! emissivity of leaf 
  real    :: scatter   (NBANDS) ! scattering coefficient of leaf (calculated as leaf_tran+leaf_refl)
  real    :: upscatter_dif (NBANDS)

  ! parameters of leaf angle distribution; see also Bonan, NCAR/TN-417+STR (LSM
  ! 1.0 technical description), p.18
  real    :: ksi    ! departure of leaf angles from a random distribution
  real    :: phi1   ! leaf distribution parameter
  real    :: phi2   ! leaf distribution parameter
  real    :: mu_bar ! average inverse diffuse optical depth per unit leaf are

  ! canopy intercepted water parameters
  real    :: cmc_lai ! max amount of liquid water on vegetation, kg/(m2 of leaf)
  real    :: cmc_pow ! power of wet fraction dependance on amount of canopy water
  real    :: csc_lai ! max amount of snow on vegetation, kg/(m2 of leaf)
  real    :: csc_pow ! power of snow-covered fraction dependance on amount of canopy snow
  real    :: fuel_intensity

  ! critical temperature for leaf drop, was internal to phenology
  real    :: tc_crit
  ! critical soil-water-stress index, used in place of fact_crit_phen and 
  ! cnst_crit_phen. It is used if and only if it's value is greater than 0
  real    :: psi_stress_crit_phen
  real    :: fact_crit_phen, cnst_crit_phen ! wilting factor and offset to 
    ! get critical value for leaf drop -- only one is non-zero at any time
  real    :: fact_crit_fire, cnst_crit_fire ! wilting factor and offset to 
    ! get critical value for fire -- only one is non-zero at the time

  real    :: smoke_fraction ! fraction of carbon lost as smoke during fires

  ! data from LM3W, temporarily here
  real    :: dat_height
  real    :: dat_lai
  real    :: dat_root_density
  real    :: dat_root_zeta
  real    :: dat_rs_min
  real    :: dat_snow_crit
end type

! ==== module data ===========================================================
integer :: idata,jdata ! iterators used in data initialization statements

! ---- namelist --------------------------------------------------------------
type(spec_data_type), save :: spdata(0:MSPECIES)

logical :: use_bucket = .false.
logical :: use_mcm_masking = .false.
real    :: mcv_min = 5.   * 4218.
real    :: mcv_lai = 0.15 * 4218.

! ---- remainder are used only for cold start
character(len=16):: vegn_to_use     = 'single-tile'
       ! 'multi-tile' for tiled vegetation
       ! 'single-tile' for geographically varying vegetation with single type per
       !     model grid cell
       ! 'uniform' for global constant vegetation, e.g., to reproduce MCM
integer :: vegn_index_constant   = 1         ! index of global constant vegn,
                                             ! used when vegn_to_use is 'uniform'
real    :: critical_root_density = 0.125

integer, dimension(1:MSPECIES) :: &
 input_cover_types=(/          -1,   -1,   -1,   -1, &
                          1,    2,    3,    4,    5,    6,    7,    8,    9/)
!character(len=4), dimension(n_dim_vegn_types) :: &
!  tile_names=      (/'be  ','bd  ','bn  ','ne  ','nd  ','g   ','d   ','t   ','a   ' /)

!  BE -- broadleaf evergreen trees
!  BD -- broadleaf deciduous trees
!  BN -- broadleaf/needleleaf trees
!  NE -- needleleaf evergreen trees
!  ND -- needleleaf deciduous trees
!  G  -- grassland
!  D  -- desert
!  T  -- tundra
!  A  -- agriculture


!       c4grass       c3grass    temp-decid      tropical     evergreen      BE     BD     BN     NE     ND      G      D      T      A
real :: dat_height(0:MSPECIES)= &
       (/  0.51,         0.51,          6.6,         19.5,          6.6,   19.5,   6.6,   8.8,   6.6,   5.9,  0.51,   1.0,  0.51,   2.9 /)
real :: dat_lai(0:MSPECIES); data dat_lai(NSPECIES:MSPECIES) &
                                                                        /   5.0,   5.0,   5.0,   5.0,   5.0,  5.0,   .001,   5.0,   5.0 /
! dat_root_density and dat_root_zeta were extended to lm3v species by copying 
! appropriate values from LaD table (e.g., grassland for both C3 and C4 grass)
real :: dat_root_density(0:MSPECIES)= &
!       c4grass       c3grass    temp-decid      tropical     evergreen      BE     BD     BN     NE     ND      G      D      T      A
       (/   1.4,          1.4,          4.2,          4.9,          2.9,    4.9,   4.2,   4.3,   2.9,   2.9,   1.4,   1.0,   1.2,  0.15 /)
real :: dat_root_zeta(0:MSPECIES)= &
       (/  0.26,         0.26,         0.29,         0.26,         0.17,   0.26,  0.29,  0.35,  0.17,  0.17,  0.26,   0.1,  0.11,  0.25 /)
real :: dat_rs_min(0:MSPECIES)= &
       (/ 56.6,          56.6,        131.0,         43.6,         69.7,   43.6, 131.0,  87.1,  69.7, 218.0,  56.6, 100.0, 170.0,  56.6 /)
real :: dat_snow_crit(0:MSPECIES)= &
!       c4grass       c3grass    temp-decid      tropical     evergreen      BE     BD     BN     NE     ND      G      D      T      A
    (/  0.0167,       0.0167,        0.0333,          0.2,       0.1333,    0.2, .0333, .0833, .1333, .1333, .0167, .0167, .0167, .0167 /)

! ==== species data imported from LM3V ======================================

!         c4 grass      c3 grass      c3 temperate  c3 tropical   c3 evergreed
real :: treefall_disturbance_rate(0:MSPECIES); data treefall_disturbance_rate(0:NSPECIES-1) &
        / 0.175,        0.185,        0.015,        0.025,        0.015 /
logical :: mortality_kills_balive(0:MSPECIES); data mortality_kills_balive(0:NSPECIES-1) &
        /.false.,      .false.,      .false.,      .false.,      .false./
integer :: pt(0:MSPECIES)= &
!       c4grass       c3grass    temp-decid      tropical     evergreen      BE     BD     BN     NE     ND      G      D      T      A
       (/ PT_C4,        PT_C3,        PT_C3,        PT_C3,        PT_C3,  PT_C3, PT_C3, PT_C3, PT_C3, PT_C3, PT_C4, PT_C4, PT_C3, PT_C3 /)
real :: alpha(0:MSPECIES,NCMPT) ; data ((alpha(idata,jdata), idata=0,MSPECIES),jdata=1,NCMPT) &
        /   0.0,          0.0,          0.0,          0.0,          0.0,     0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0., & ! reproduction
            0.0,          0.0,          0.0,          0.0,          0.0,     0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0., & ! sapwood
            1.0,          1.0,          1.0,          0.8,         0.12,    0.8,   1.0,   1.0,  0.12,   1.0,   1.0,   1.0,   1.0,   1.0, & ! leaf
            0.9,         0.55,          1.0,          0.8,          0.6,    0.8,   1.0,   1.0,   0.6,   1.0,  0.55,   0.9,  0.55,  0.55, & ! root
            0.0,          0.0,          0.0,          0.0,          0.0,     0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0., & ! virtual leaf
            0.0,          0.0,        0.006,        0.012,        0.006,  0.012, 0.006, 0.006, 0.006, 0.006,    0.,    0.,    0.,    0. /  ! structural

! From Foley (Ibis model) 1996 gbc v10 pp. 603-628
real :: beta(0:MSPECIES,NCMPT) ; data ((beta(idata,jdata), idata=0,NSPECIES-1),jdata=1,NCMPT) &
        /   0.0,          0.0,          0.0,          0.0,          0.0,& ! reproduction
            0.0,          0.0,          0.0,          0.0,          0.0,& ! sapwood
            0.0,          0.0,          0.0,          0.0,          0.0,& ! leaf
           1.25,         1.25,         1.25,         1.25,         1.25,& ! root
            0.0,          0.0,          0.0,          0.0,          0.0,& ! virtual leaf
            0.0,          0.0,          0.0,          0.0,          0.0 / ! structural

! root parameters
real :: dfr(0:MSPECIES) ; data dfr &
!       c4grass       c3grass    temp-decid      tropical     evergreen      BE      BD      BN      NE      ND       G       D       T      A
        /   2.2,          2.2,          5.8,          5.8,          5.8,    5.8,    5.8,    5.8,    5.8,    5.8,    2.2,    2.2,    2.2,    2.2 /
real :: srl(0:MSPECIES); data srl & ! specific root length, m/(kg C)
        / 236e3,        236e3,       24.4e3,       24.4e3,       24.4e3, 24.4e3, 24.4e3, 24.4e3, 24.4e3, 24.4e3,  236e3,  236e3,   60e3,   60e3 /
real :: root_r(0:MSPECIES); data root_r & ! radius of fine roots, m
        /1.1e-4,       1.1e-4,       2.9e-4,       2.9e-4,       2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 1.1e-4, 1.1e-4, 2.2e-4, 2.2e-4 /
real :: root_perm(0:MSPECIES); data root_perm & ! fine root membrane permeability per unit membrane area, kg/(m3 s)
        / 1e-5,          1e-5,         1e-5,         1e-5,         1e-5,   1e-5,   1e-5,   1e-5,   1e-5,   1e-5,   1e-5,   1e-5,   1e-5,   1e-5 /
! Specific root length is from Jackson et al., 1997, PNAS  Vol.94, pp.7362--7366, 
! converted to m/(kg C) from m/(g biomass). Biomass/C mass ratio was assumed 
! to be 2. The fine root radius is from the same source.
!
! Root membrane permeability is "high" value from Siqueira et al., 2008, Water 
! Resource Research Vol. 44, W01432, converted to mass units

real :: c1(0:MSPECIES); data c1(0:NSPECIES-1) &
        /   1.358025,     2.222222,     0.4807692,    0.3333333,    0.1948718 /
real :: c2(0:MSPECIES); data c2(0:NSPECIES-1) &
        /   0.4004486,    0.4004486,    0.4004486,    0.3613833,    0.1509976 /
real :: c3(0:MSPECIES); data c3(0:NSPECIES-1) &
        /   0.5555555,    0.5555555,    0.4423077,    1.230769,     0.5897437 /

! leaf radiation parameters
real :: leaf_refl(0:MSPECIES,NBANDS) ; data leaf_refl & ! leaf reflectance
        /   0.11,        0.11,         0.10,         0.10,         0.07,  0.149, 0.130, 0.132, 0.126, 0.143, 0.182, 0.300, 0.139, 0.160, & ! VIS
            0.45,        0.45,         0.45,         0.45,         0.35,  0.149, 0.130, 0.132, 0.126, 0.143, 0.182, 0.300, 0.139, 0.160  / ! NIR
real :: leaf_tran(0:MSPECIES,NBANDS) ; data leaf_tran & ! leaf transmittance
        /   0.07,        0.07,         0.05,         0.05,         0.05,     0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0., & ! VIS
            0.25,        0.25,         0.25,         0.25,         0.10,     0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.  / ! NIR
real :: leaf_emis(0:MSPECIES)= & ! leaf emissivity
       (/   1.00,        1.00,         1.00,         1.00,         1.00,   0.98,   0.96,  0.97,  0.98,  0.96, 0.96,   1.0,  0.96,  0.96  /)
real :: ksi(0:MSPECIES)= & ! leaf inclination index
       (/      0.,          0.,           0.,           0.,          0.,     0.,     0.,    0.,    0.,    0.,    0.,    0.,    0.,   0.  /)
real :: min_cosz = 0.01 ! minimum allowed value of cosz for vegetation radiation
   ! properties calculations.
   ! It probably doesn't make sense to set it any less than the default value, because the angular 
   ! diameter of the sun is about 0.01 radian (0.5 degree), so the spread of the direct radiation 
   ! zenith angles is about this. Besides, the sub-grid variations of land surface slope are 
   ! probably even larger that that. 

! canopy interception parameters
real :: cmc_lai(0:MSPECIES)= & ! maximum canopy water conntent per unit LAI
       (/    0.1,         0.1,          0.1,          0.1,          0.1,    0.1,    0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,  0.1  /)
real :: cmc_pow(0:MSPECIES)= & ! power of the wet canopy fraction relation 
  (/   TWOTHIRDS,   TWOTHIRDS,    TWOTHIRDS,    TWOTHIRDS,    TWOTHIRDS,     1.,     1.,    1.,    1.,    1.,    1.,    1.,    1.,   1.  /)
real :: csc_lai(0:MSPECIES)= & ! maximum canopy snow conntent per unit LAI
       (/    0.1,         0.1,          0.1,          0.1,          0.1,    0.1,    0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,  0.1  /)
real :: csc_pow(0:MSPECIES)= & ! power of the snow-covered fraction relation 
  (/   TWOTHIRDS,   TWOTHIRDS,    TWOTHIRDS,    TWOTHIRDS,    TWOTHIRDS,     1.,     1.,    1.,    1.,    1.,    1.,    1.,    1.,   1.  /)
real :: cmc_eps = 0.01 ! value of w/w_max for transition to linear function; 
                       ! the same value is used for liquid and snow

real :: fuel_intensity(0:MSPECIES) ; data fuel_intensity(0:NSPECIES-1) &
        /    1.0,         1.0,        0.002,        0.002,        0.004 /
real :: leaf_size(0:MSPECIES)= & ! characteristic leaf size
       (/   0.04,        0.04,         0.04,          0.04,        0.04,    0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1  /)
! photosynthesis parameters
real :: Vmax(0:MSPECIES)= & ! max rubisco rate
       (/  35e-6,       70e-6,        70e-6,         70e-6,       70e-6,  70e-6, 70e-6, 70e-6, 70e-6, 70e-6, 35e-6, 35e-6, 70e-6, 70e-6  /)
real :: m_cond(0:MSPECIES)= & ! factor of stomatal conductance
       (/    4.0,         9.0,          9.0,           9.0,         9.0,    9.0,   9.0,   9.0,   9.0,   9.0,   4.0,   4.0,   9.0,   9.0  /)
real :: alpha_phot(0:MSPECIES)= & ! photosynthesis efficiency
       (/   0.05,        0.06,         0.06,          0.06,        0.06,   0.06,  0.06,  0.06,  0.06,  0.06,  0.05,  0.05,  0.06,  0.06  /)
real :: gamma_resp(0:MSPECIES)= &
       (/   0.03,        0.02,         0.02,          0.02,        0.03,   0.02,  0.02,  0.02,  0.03,  0.03,  0.03,  0.03,  0.03,  0.02  /)
!       c4grass       c3grass    temp-decid      tropical     evergreen      BE     BD     BN     NE     ND      G      D      T      A
real :: tc_crit(0:MSPECIES)= &
       (/ 283.16,      278.16,       283.16,        283.16,      263.16,      0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0. /)
real :: psi_stress_crit_phen(0:MSPECIES)= & ! iff > 0, critical soil-water-stress index for leaf drop, overrides water content
       (/    0.0,         0.0,          0.0,           0.0,         0.0,    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0 /)       
real :: cnst_crit_phen(0:MSPECIES)= & ! constant critical value for leaf drop
       (/    0.1,         0.1,          0.1,           0.1,         0.1,    0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1  /)
real :: fact_crit_phen(0:MSPECIES)= & ! factor for wilting to get critical value for leaf drop
       (/    0.0,         0.0,          0.0,           0.0,         0.0,    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0  /)
real :: cnst_crit_fire(0:MSPECIES)= & ! constant critical value for leaf drop
       (/    0.1,         0.1,          0.1,           0.1,         0.1,    0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1  /)  
real :: fact_crit_fire(0:MSPECIES)= & ! factor for wilting to get critical value for fire
       (/    0.0,         0.0,          0.0,           0.0,         0.0,    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0  /)
real :: wet_leaf_dreg(0:MSPECIES) = & ! wet leaf photosynthesis down-regulation: 0.3 means 
        ! photosynthesis of completely wet leaf will be 30% less than that of dry one,
        ! provided everything else is the same 
       (/    0.3,         0.3,          0.3,           0.3,         0.3,    0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3  /)
real :: leaf_age_onset(0:MSPECIES) = & ! onset of Vmax decrease due to leaf aging, days
       (/  100.0,       100.0,        100.0,         100.0,       100.0,  100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0  /)
real :: leaf_age_tau(0:MSPECIES) = &  ! e-folding time of Vmax decrease due to leaf aging, days (0 or less means no aging)
       (/    0.0,         0.0,          0.0,           0.0,         0.0,    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0  /)

real :: soil_carbon_depth_scale = 0.2   ! depth of active soil for carbon decomposition
real :: cold_month_threshold    = 283.0 ! monthly temperature threshold for calculations of number of cold months
real :: smoke_fraction(0:MSPECIES) = & ! fration of carbon lost as smoke
       (/    0.9,         0.9,          0.9,           0.9,         0.9,    0.9,   0.9,   0.9,   0.9,   0.9,   0.9,   0.9,   0.9,   0.9  /)
real :: agf_bs         = 0.8 ! ratio of above ground stem to total stem
real :: K1 = 10.0, K2 = 0.05 ! soil decomposition parameters
real :: fsc_liv        = 0.8
real :: fsc_wood       = 0.2
real :: tau_drip_l     = 21600.0 ! canopy water residence time, for drip calculations
real :: tau_drip_s     = 86400.0 ! canopy snow residence time, for drip calculations
real :: GR_factor = 0.33 ! growth respiration factor     

real :: tg_c3_thresh = 1.5 ! threshold biomass between tree and grass for C3 plants
real :: tg_c4_thresh = 2.0 ! threshold biomass between tree and grass for C4 plants
real :: fsc_pool_spending_time = 1.0 ! time (yrs) during which intermediate pool of 
                  ! fast soil carbon is entirely converted to the fast soil carbon
real :: ssc_pool_spending_time = 1.0 ! time (yrs) during which intermediate pool of
                  ! slow soil carbon is entirely converted to the slow soil carbon
real :: harvest_spending_time(N_HARV_POOLS) = &
     (/1.0, 1.0, 1.0, 1.0, 10.0, 100.0/)
     ! time (yrs) during which intermediate pool of harvested carbon is completely
     ! released to the atmosphere. 
     ! NOTE: a year in the above *_spending_time definitions is exactly 365*86400 seconds
real :: l_fract      = 0.5 ! fraction of the leaves retained after leaf drop
real :: T_transp_min = 0.0 ! lowest temperature at which transporation is enabled
                           ! 0 means no limit, lm3v value is 268.0
! boundaries of wood biomass bins for secondary veg. (kg C/m2); used to decide 
! whether secondary vegetation tiles can be merged or not. MUST BE IN ASCENDING 
! ORDER.
real  :: scnd_biomass_bins(10) &  
     = (/ 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 1000.0 /)
real :: phen_ev1 = 0.5, phen_ev2 = 0.9 ! thresholds for evergreen/decidious 
      ! differentiation (see phenology_type in cohort.F90)

namelist /vegn_data_nml/ &
  vegn_to_use,  input_cover_types, &
  mcv_min, mcv_lai, &
  use_bucket, use_mcm_masking, vegn_index_constant, &
  critical_root_density, &

  dat_height, dat_lai, dat_root_density, dat_root_zeta, dat_rs_min, dat_snow_crit, &
  ! vegetation data, imported from LM3V
  pt, Vmax, m_cond, alpha_phot, gamma_resp, wet_leaf_dreg, &
  leaf_age_onset, leaf_age_tau, &
  treefall_disturbance_rate, mortality_kills_balive, fuel_intensity, &
  alpha, beta, c1,c2,c3, &
  dfr, &
  srl, root_r, root_perm, &
  cmc_lai, cmc_pow, csc_lai, csc_pow, cmc_eps, &
  min_cosz, &
  leaf_refl, leaf_tran, leaf_emis, ksi, &
  leaf_size, &
  soil_carbon_depth_scale, cold_month_threshold, &

  smoke_fraction, agf_bs, K1,K2, fsc_liv, fsc_wood, &
  tau_drip_l, tau_drip_s, GR_factor, tg_c3_thresh, tg_c4_thresh, &
  fsc_pool_spending_time, ssc_pool_spending_time, harvest_spending_time, &
  l_fract, T_transp_min,  tc_crit, psi_stress_crit_phen, &
  cnst_crit_phen, fact_crit_phen, cnst_crit_fire, fact_crit_fire, &
  scnd_biomass_bins, phen_ev1, phen_ev2


contains ! ###################################################################



! ============================================================================
subroutine read_vegn_data_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: i
  
  type(table_printer_type) :: table

  call write_version_number(version, tagname)
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=vegn_data_nml, iostat=io)
  ierr = check_nml_error(io, 'vegn_data_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=vegn_data_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'vegn_data_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif

  unit=stdlog()

  ! reconcile values of fact_crit_phen and cnst_crit_phen
  cnst_crit_phen = max(0.0,min(1.0,cnst_crit_phen))
  fact_crit_phen = max(0.0,fact_crit_phen)
  where (cnst_crit_phen/=0) fact_crit_phen=0.0
  write(unit,*)'reconciled fact_crit_phen and cnst_crit_phen'

  ! do the same for fire
  cnst_crit_fire = max(0.0,min(1.0,cnst_crit_fire))
  fact_crit_fire = max(0.0,fact_crit_fire)
  where (cnst_crit_fire/=0) fact_crit_fire=0.0 
  write(unit,*)'reconciled fact_crit_fire and cnst_crit_fire'

  ! initialize vegetation data structure

  spdata%dat_height = dat_height
  spdata%dat_lai = dat_lai
  spdata%dat_root_density = dat_root_density
  spdata%dat_root_zeta = dat_root_zeta
  spdata%dat_rs_min = dat_rs_min
  spdata%dat_snow_crit = dat_snow_crit

  spdata%treefall_disturbance_rate = treefall_disturbance_rate
  spdata%mortality_kills_balive    = mortality_kills_balive
  spdata%fuel_intensity            = fuel_intensity
 
  spdata%pt         = pt
  spdata%Vmax       = Vmax
  spdata%m_cond     = m_cond
  spdata%alpha_phot = alpha_phot
  spdata%gamma_resp = gamma_resp
  spdata%wet_leaf_dreg = wet_leaf_dreg
  spdata%leaf_age_onset = leaf_age_onset
  spdata%leaf_age_tau = leaf_age_tau
  spdata%dfr        = dfr

  spdata%srl        = srl
  spdata%root_r     = root_r
  spdata%root_perm  = root_perm
  
  spdata%c1 = c1
  spdata%c2 = c2
  spdata%c3 = c3

  spdata%cmc_lai = cmc_lai
  spdata%cmc_pow = cmc_pow
  spdata%csc_lai = csc_lai
  spdata%csc_pow = csc_pow
  
  spdata%leaf_size = leaf_size

  spdata%tc_crit   = tc_crit
  spdata%psi_stress_crit_phen = psi_stress_crit_phen
  spdata%cnst_crit_phen = cnst_crit_phen
  spdata%fact_crit_phen = fact_crit_phen
  spdata%cnst_crit_fire = cnst_crit_fire
  spdata%fact_crit_fire = fact_crit_fire

  spdata%smoke_fraction = smoke_fraction

  do i = 0, MSPECIES
     spdata(i)%alpha     = alpha(i,:)
     spdata(i)%beta      = beta(i,:)
     spdata(i)%leaf_refl = leaf_refl(i,:)
     spdata(i)%leaf_tran = leaf_tran(i,:)
     spdata(i)%leaf_emis = leaf_emis(i)
     spdata(i)%ksi       = ksi(i)
     call init_derived_species_data(spdata(i))
  enddo

  ! register selectors for land use type-specific diagnostics
  do i=1, N_LU_TYPES
     call register_tile_selector(landuse_name(i), long_name=landuse_longname(i),&
          tag = SEL_VEGN, idata1 = LU_SEL_TAG, idata2 = i )
  enddo
  
  ! register selectors for species-specific diagnostics
  do i=0,NSPECIES-1
     call register_tile_selector(species_name(i), long_name=species_longname(i),&
          tag = SEL_VEGN, idata1 = SP_SEL_TAG, idata2 = i )
  enddo

  ! register selector for natural grass
  call register_tile_selector('ntrlgrass', long_name='natural (non-human-maintained) grass',&
          tag = SEL_VEGN, idata1 = NG_SEL_TAG)

  write (unit, nml=vegn_data_nml)

  call init_with_headers(table,species_name)
  call add_row(table,'Treefall dist. rate', spdata(:)%treefall_disturbance_rate)
  call add_row(table,'Mortality kills balive', spdata(:)%mortality_kills_balive)
  call add_row(table,'Phisiology Type', spdata(:)%pt)
  call add_row(table,'C1',            spdata(:)%c1)
  call add_row(table,'C2',            spdata(:)%c2)
  call add_row(table,'C3',            spdata(:)%c3)

  call add_row(table,'alpha_leaf',    spdata(:)%alpha(CMPT_LEAF))
  call add_row(table,'alpha_root',    spdata(:)%alpha(CMPT_ROOT))
  call add_row(table,'alpha_vleaf',   spdata(:)%alpha(CMPT_VLEAF))
  call add_row(table,'alpha_sapwood', spdata(:)%alpha(CMPT_SAPWOOD))
  call add_row(table,'alpha_wood',    spdata(:)%alpha(CMPT_WOOD))
  call add_row(table,'alpha_repro',   spdata(:)%alpha(CMPT_REPRO))

  call add_row(table,'beta_leaf',     spdata(:)%beta(CMPT_LEAF))
  call add_row(table,'beta_root',     spdata(:)%beta(CMPT_ROOT))
  call add_row(table,'beta_vleaf',    spdata(:)%beta(CMPT_VLEAF))
  call add_row(table,'beta_sapwood',  spdata(:)%beta(CMPT_SAPWOOD))
  call add_row(table,'beta_wood',     spdata(:)%beta(CMPT_WOOD))
  call add_row(table,'beta_repro',    spdata(:)%beta(CMPT_REPRO))

  call add_row(table,'dfr',           spdata(:)%dfr)

  call add_row(table,'srl',           spdata(:)%srl)
  call add_row(table,'root_r',        spdata(:)%root_r)
  call add_row(table,'root_perm',     spdata(:)%root_perm)

  call add_row(table,'specific_leaf_area', spdata(:)%specific_leaf_area)
  call add_row(table,'leaf_size',     spdata(:)%leaf_size)
  call add_row(table,'leaf_life_span',spdata(:)%leaf_life_span)

  call add_row(table,'alpha_phot',    spdata(:)%alpha_phot)
  call add_row(table,'m_cond',        spdata(:)%m_cond)
  call add_row(table,'Vmax',          spdata(:)%Vmax)
  call add_row(table,'gamma_resp',    spdata(:)%gamma_resp)
  call add_row(table,'wet_leaf_dreg', spdata(:)%wet_leaf_dreg)
  call add_row(table,'leaf_age_onset',spdata(:)%leaf_age_onset)
  call add_row(table,'leaf_age_tau',  spdata(:)%leaf_age_tau)

  call add_row(table,'leaf_refl_vis', spdata(:)%leaf_refl(BAND_VIS))
  call add_row(table,'leaf_refl_nir', spdata(:)%leaf_refl(BAND_NIR))
  call add_row(table,'leaf_tran_vis', spdata(:)%leaf_tran(BAND_VIS))
  call add_row(table,'leaf_tran_nir', spdata(:)%leaf_tran(BAND_NIR))
  call add_row(table,'leaf_emis',     spdata(:)%leaf_emis)
  call add_row(table,'ksi',           spdata(:)%ksi)
  call add_row(table,'phi1',          spdata(:)%phi1)
  call add_row(table,'phi2',          spdata(:)%phi2)
  call add_row(table,'mu_bar',        spdata(:)%mu_bar)

  call add_row(table,'cmc_lai',       spdata(:)%cmc_lai)
  call add_row(table,'cmc_pow',       spdata(:)%cmc_pow)
  call add_row(table,'csc_lai',       spdata(:)%csc_lai)
  call add_row(table,'csc_pow',       spdata(:)%csc_pow)
  call add_row(table,'fuel_intensity',spdata(:)%fuel_intensity)

  call add_row(table,'tc_crit',       spdata(:)%tc_crit)
  call add_row(table,'psi_stress_crit_phen', spdata(:)%psi_stress_crit_phen)
  call add_row(table,'fact_crit_phen',spdata(:)%fact_crit_phen)
  call add_row(table,'cnst_crit_phen',spdata(:)%cnst_crit_phen)
  call add_row(table,'fact_crit_fire',spdata(:)%fact_crit_fire)
  call add_row(table,'cnst_crit_fire',spdata(:)%cnst_crit_fire)

  call add_row(table,'smoke_fraction',spdata(:)%smoke_fraction)

  call print(table,stdout())
  call print(table,unit)
  
end subroutine 


! ============================================================================
subroutine init_derived_species_data(sp)
   type(spec_data_type), intent(inout) :: sp

   integer :: j
   
   sp%leaf_life_span     = 12.0/sp%alpha(CMPT_LEAF) ! in months
   ! calculate specific leaf area (cm2/g(biomass))
   ! Global Raich et al 94 PNAS pp 13730-13734
   sp%specific_leaf_area = 10.0**(2.4 - 0.46*log10(sp%leaf_life_span));       
   ! convert to (m2/kg(carbon)
   sp%specific_leaf_area = C2B*sp%specific_leaf_area*1000.0/10000.0

! rho_wood is not used anywhere?
!     ! the relationship from Moorcroft, based on Reich
!     ! units kg C/m^3, hence the factor of 0.001 to convert from g/cm^3
!      sp%rho_wood = (0.5 + 0.2*(sp%leaf_life_span-1))*0.001;
!     if (sp%rho_wood > 500.) sp%rho_wood = 0.5*0.001;

   sp%phi1=0.5-0.633*sp%ksi-0.33*sp%ksi**2;
   sp%phi2=0.877*(1.0-2.0*sp%phi1);
   if(sp%ksi /= 0) then
      sp%mu_bar = &
           (1-sp%phi1/sp%phi2*log(1+sp%phi2/sp%phi1))&
           / sp%phi2
   else
      ! in degenerate case of spherical leaf angular distribution the above 
      ! formula for mu_bar gives an undefined value, so we handle it separately 
      sp%mu_bar = 1.0
   endif
   do j = 1,NBANDS
      sp%scatter(j)       = sp%leaf_refl(j)+sp%leaf_tran(j);
      sp%upscatter_dif(j) = 0.5*(sp%scatter(j) + & 
           (sp%leaf_refl(j)-sp%leaf_tran(j))*(1+sp%ksi)**2/4);
   enddo
end subroutine


end module

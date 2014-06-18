module vegn_tile_mod

use fms_mod, only : &
     write_version_number, stdlog, error_mesg, FATAL
use constants_mod, only : &
     tfreeze, hlf

use land_constants_mod, only : NBANDS
use land_io_mod, only : &
     init_cover_field
use land_tile_selectors_mod, only : &
     tile_selector_type, SEL_VEGN

use vegn_data_mod, only : &
     NSPECIES, MSPECIES, NCMPT, C2B, &
     read_vegn_data_namelist, spdata, &
     vegn_to_use,  input_cover_types, &
     mcv_min, mcv_lai, &
     vegn_index_constant, &
     agf_bs, BSEED, LU_NTRL, LU_SCND, N_HARV_POOLS, &
     LU_SEL_TAG, SP_SEL_TAG, NG_SEL_TAG, &
     SP_C3GRASS, SP_C4GRASS, &
     scnd_biomass_bins

use vegn_cohort_mod, only : vegn_cohort_type, vegn_phys_prog_type, &
     height_from_biomass, lai_from_biomass, update_bio_living_fraction, &
     cohort_uptake_profile, cohort_root_properties, update_biomass_pools

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_tile_type

public :: new_vegn_tile, delete_vegn_tile
public :: vegn_tiles_can_be_merged, merge_vegn_tiles
public :: vegn_is_selected
public :: get_vegn_tile_tag
public :: vegn_tile_stock_pe
public :: vegn_tile_carbon ! returns total carbon per tile
public :: vegn_tile_heat ! returns hate content of the vegetation

public :: read_vegn_data_namelist
public :: vegn_cover_cold_start

public :: vegn_uptake_profile
public :: vegn_root_properties
public :: vegn_data_rs_min
public :: vegn_seed_supply
public :: vegn_seed_demand

public :: vegn_tran_priority ! returns transition priority for land use 

public :: vegn_add_bliving
public :: update_derived_vegn_data  ! given state variables, calculate derived values
! =====end of public interfaces ==============================================
interface new_vegn_tile
   module procedure vegn_tile_ctor
   module procedure vegn_tile_copy_ctor
end interface


! ==== module constants ======================================================
character(len=*), parameter   :: &
     version = '$Id: vegn_tile.F90,v 20.0 2013/12/13 23:31:19 fms Exp $', & 
     tagname = '$Name: tikal $', &
     module_name = 'vegn_tile_mod'

! ==== types =================================================================
type :: vegn_tile_type
   integer :: tag ! kind of the tile
   integer :: landuse = LU_NTRL

   integer :: n_cohorts = 0
   type(vegn_cohort_type), pointer :: cohorts(:)=>NULL()

   real :: age=0.0 ! tile age

   ! fields for smoothing out the contribution of the spike-type processes (e.g. 
   ! harvesting) to the soil carbon pools over some period of time
   real :: fsc_pool=0.0, fsc_rate=0.0 ! for fast soil carbon
   real :: ssc_pool=0.0, ssc_rate=0.0 ! for slow soil carbon

   real :: csmoke_pool=0.0 ! carbon lost through fires, kg C/m2 
   real :: csmoke_rate=0.0 ! rate of release of the above to atmosphere, kg C/(m2 yr)

   real :: harv_pool(N_HARV_POOLS) = 0.0 ! pools of harvested carbon, kg C/m2
   real :: harv_rate(N_HARV_POOLS) = 0.0 ! rates of spending (release to the atmosphere), kg C/(m2 yr)

   ! values for the diagnostic of carbon budget and soil carbon acceleration
   real :: ssc_out=0.0
   real :: fsc_out=0.0
   real :: veg_in=0.0, veg_out=0.0

   real :: disturbance_rate(0:1) = 0 ! 1/year
   real :: lambda = 0.0 ! cumulative drought months per year
   real :: fuel   = 0.0 ! fuel over dry months
   real :: litter = 0.0 ! litter flux

   ! monthly accumulated/averaged values
   real :: theta_av_phen = 0.0 ! relative soil_moisture availability not soil moisture
   real :: theta_av_fire = 0.0
   real :: psist_av = 0.0 ! soil water stress index
   real :: tsoil_av = 0.0 ! bulk soil temperature
   real :: tc_av    = 0.0 ! leaf temperature
   real :: precip_av= 0.0 ! precipitation

   ! accumulation counters for long-term averages (monthly and annual). Having
   ! these counters in the tile is a bit stupid, since the values are the same for
   ! each tile, but it simplifies the current code, and they are going away when we
   ! switch to exponential averaging in any case.
   integer :: n_accum = 0 ! number of accumulated values for monthly averages
   integer :: nmn_acm = 0 ! number of accumulated values for annual averages
   ! annual-mean values
   real :: t_ann  = 0.0 ! annual mean T, degK
   real :: t_cold = 0.0 ! average temperature of the coldest month, degK
   real :: p_ann  = 0.0 ! annual mean precip
   real :: ncm    = 0.0 ! number of cold months
   ! annual accumulated values
   real :: t_ann_acm  = 0.0 ! accumulated annual temperature for t_ann
   real :: t_cold_acm = 0.0 ! temperature of the coldest month in current year
   real :: p_ann_acm  = 0.0 ! accumulated annual precipitation for p_ann
   real :: ncm_acm    = 0.0 ! accumulated number of cold months


   ! it's probably possible to get rid of the fields below
   real :: npp=0.0 ! net primary productivity
   real :: nep=0.0 ! net ecosystem productivity
   real :: rh=0.0 ! soil carbon lost to the atmosphere
   real :: total_biomass !
   real :: area_disturbed_by_treefall
   real :: area_disturbed_by_fire
   real :: total_disturbance_rate
end type vegn_tile_type

! ==== module data ===========================================================
real, public :: &
     cpw = 1952.0, & ! specific heat of water vapor at constant pressure
     clw = 4218.0, & ! specific heat of water (liquid)
     csw = 2106.0    ! specific heat of water (ice)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! ============================================================================
function vegn_tile_ctor(tag) result(ptr)
  type(vegn_tile_type), pointer :: ptr ! return value
  integer, intent(in)  :: tag ! kind of tile

  allocate(ptr)
  ptr%tag = tag
end function vegn_tile_ctor

! ============================================================================
function vegn_tile_copy_ctor(vegn) result(ptr)
  type(vegn_tile_type), pointer :: ptr ! return value
  type(vegn_tile_type), intent(in) :: vegn ! return value

  allocate(ptr)
  ! copy all non-pointer members
  ptr=vegn
  ! copy pointer members (cohorts)
  allocate(ptr%cohorts(ptr%n_cohorts))
  ptr%cohorts(:) = vegn%cohorts(1:ptr%n_cohorts)
end function vegn_tile_copy_ctor

! ============================================================================
subroutine delete_vegn_tile(vegn)
  type(vegn_tile_type), pointer :: vegn

  deallocate(vegn%cohorts)
  deallocate(vegn)
end subroutine delete_vegn_tile

! =============================================================================
function vegn_tiles_can_be_merged(vegn1,vegn2) result(response)
  logical :: response
  type(vegn_tile_type), intent(in) :: vegn1,vegn2

  real    :: b1, b2 
  integer :: i, i1, i2

  if (vegn1%landuse /= vegn2%landuse) then
     response = .false. ! different land use types can't be merged
  else if (vegn1%landuse == LU_SCND) then ! secondary vegetation tiles
     ! get tile wood biomasses
     b1 = get_vegn_tile_bwood(vegn1)
     b2 = get_vegn_tile_bwood(vegn2)
     ! find biomass bins where each the tiles belongs to
     i1 = 0 ; i2 = 0
     do i = 1, size(scnd_biomass_bins(:))
        if (b1>scnd_biomass_bins(i)) i1 = i
        if (b2>scnd_biomass_bins(i)) i2 = i
     enddo
     ! tiles can be merged only if biomasses belong to the same bin
     response = (i1 == i2)
  else
     response = .true. ! non-secondary tiles of the same land use type can always be merged
  endif
end function


! ============================================================================
subroutine merge_vegn_tiles(t1,w1,t2,w2)
  type(vegn_tile_type), intent(in) :: t1
  type(vegn_tile_type), intent(inout) :: t2
  real, intent(in) :: w1, w2 ! relative weights
  
  ! ---- local vars
  real :: x1, x2 ! normalized relative weights
  real :: HEAT1, HEAT2 ! heat stored in respective canopies
  type(vegn_cohort_type), pointer :: c1, c2
  
  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1.0 - x1

  ! the following assumes that there is one, and only one, cohort per tile 
  c1 => t1%cohorts(1)
  c2 => t2%cohorts(1)
  ! define macro for merging cohort values
#define __MERGE__(field) c2%field = x1*c1%field + x2*c2%field
  HEAT1 = (clw*c1%prog%Wl + csw*c1%prog%Ws + c1%mcv_dry)*(c1%prog%Tv-tfreeze)
  HEAT2 = (clw*c2%prog%Wl + csw*c2%prog%Ws + c2%mcv_dry)*(c2%prog%Tv-tfreeze)
  __MERGE__(prog%Wl)
  __MERGE__(prog%Ws)

  __MERGE__(bl)      ! biomass of leaves, kg C/m2
  __MERGE__(blv)     ! biomass of virtual leaves (labile store), kg C/m2
  __MERGE__(br)      ! biomass of fine roots, kg C/m2
  __MERGE__(bsw)     ! biomass of sapwood, kg C/m2
  __MERGE__(bwood)   ! biomass of heartwood, kg C/m2
  __MERGE__(bliving) ! leaves, fine roots, and sapwood biomass
  
  __MERGE__(carbon_gain) ! carbon gain during a day, kg C/m2
  __MERGE__(carbon_loss) ! carbon loss during a day, kg C/m2 [diag only]
  __MERGE__(bwood_gain)  ! heartwood gain during a day, kg C/m2

  ! should we do update_derived_vegn_data here? to get mcv_dry, etc
  call update_biomass_pools(c2)

  ! calculate the resulting dry heat capacity
  c2%mcv_dry = max(mcv_min,mcv_lai*c2%lai)
  ! update canopy temperature -- just merge it based on area weights if the heat 
  ! capacities are zero, or merge it based on the heat content if the heat contents
  ! are non-zero
  if(HEAT1==0.and.HEAT2==0) then
     __MERGE__(prog%Tv)
  else
     c2%prog%Tv = (HEAT1*x1+HEAT2*x2) / &
          (clw*c2%prog%Wl + csw*c2%prog%Ws + c2%mcv_dry) + tfreeze
  endif

#undef  __MERGE__
! re-define macro for tile values
#define __MERGE__(field) t2%field = x1*t1%field + x2*t2%field

  __MERGE__(age);
  
  __MERGE__(fsc_pool); __MERGE__(fsc_rate)
  __MERGE__(ssc_pool); __MERGE__(ssc_rate)

  __MERGE__(csmoke_pool)
  __MERGE__(csmoke_rate)

  __MERGE__(harv_pool)
  __MERGE__(harv_rate)

  ! do we need to merge these?
  __MERGE__(ssc_out)
  __MERGE__(fsc_out)
  __MERGE__(veg_in); __MERGE__(veg_out)
  
  ! or these?
  __MERGE__(disturbance_rate)
  __MERGE__(lambda)     ! cumulative drought months per year
  __MERGE__(fuel)       ! fuel over dry months
  __MERGE__(litter)     ! litter flux

  ! monthly accumulated/averaged values
  __MERGE__(theta_av_phen)   ! relative soil_moisture availability not soil moisture
  __MERGE__(theta_av_fire)
  __MERGE__(psist_av)   ! water potential divided by permanent wilting potential
  __MERGE__(tsoil_av)   ! bulk soil temperature
  __MERGE__(tc_av)      ! leaf temperature
  __MERGE__(precip_av)  ! precipitation

  ! annual-mean values
  __MERGE__(t_ann)      ! annual mean T, degK
  __MERGE__(t_cold)     ! average temperature of the coldest month, degK
  __MERGE__(p_ann)      ! annual mean precip
  __MERGE__(ncm)        ! number of cold months

  ! annual accumulated values
  __MERGE__(t_ann_acm)  ! accumulated annual temperature for t_ann
  __MERGE__(t_cold_acm) ! temperature of the coldest month in current year
  __MERGE__(p_ann_acm)  ! accumulated annual precipitation for p_ann
  __MERGE__(ncm_acm)    ! accumulated number of cold months

#undef __MERGE__

end subroutine merge_vegn_tiles


! ============================================================================
! given a vegetation tile with the state variables set up, calculate derived
! parameters to get a consistent state
! NOTE: this subroutine does not call update_biomass_pools, although some 
! of the calculations are the same. The reason is because this function may 
! be used in the situation when the biomasses are not precisely consistent, for
! example when they come from the data override or from initial conditions.
subroutine update_derived_vegn_data(vegn)
  type(vegn_tile_type), intent(inout) :: vegn
  
  ! ---- local vars 
  type(vegn_cohort_type), pointer :: cc ! pointer to the current cohort
  integer :: i  ! cohort index
  integer :: sp ! shorthand for the vegetation species
  
  ! given that the cohort state variables are initialized, fill in
  ! the intermediate variables
  do i = 1,vegn%n_cohorts
    cc=>vegn%cohorts(i)
    
    sp = cc%species
    ! set the physiology type according to species
    cc%pt     = spdata(sp)%pt
    ! calculate total biomass, calculate height
    cc%b      = cc%bliving + cc%bwood
    cc%height = height_from_biomass(cc%b);
    ! update fractions of the living biomass
    call update_bio_living_fraction(cc)
    cc%bs     = cc%bsw + cc%bwood;   
    cc%bstem  = agf_bs*cc%bs;
    cc%babove = cc%bl + agf_bs*cc%bs; 

    if(sp<NSPECIES) then ! LM3V species
       ! calculate the leaf area index based on the biomass of leaves
       cc%lai = lai_from_biomass(cc%bl, sp)
       ! calculate the root density as the total biomass below ground, in
       ! biomass (not carbon!) units
       cc%root_density = (cc%br + (cc%bsw+cc%bwood+cc%blv)*(1-agf_bs))*C2B
    else
       cc%height        = spdata(sp)%dat_height
       cc%lai           = spdata(sp)%dat_lai
       cc%root_density  = spdata(sp)%dat_root_density
    endif
    cc%sai           = 0.035*cc%height
    cc%leaf_size     = spdata(sp)%leaf_size
    cc%root_zeta     = spdata(sp)%dat_root_zeta
    cc%rs_min        = spdata(sp)%dat_rs_min
    cc%leaf_refl     = spdata(sp)%leaf_refl
    cc%leaf_tran     = spdata(sp)%leaf_tran
    cc%leaf_emis     = spdata(sp)%leaf_emis
    cc%snow_crit     = spdata(sp)%dat_snow_crit
  
    ! putting this initialization within the cohort loop is probably incorrect 
    ! in case of multiple-cohort vegetation, however for a single cohort it works
    cc%Wl_max   = spdata(sp)%cmc_lai*cc%lai
    cc%Ws_max   = spdata(sp)%csc_lai*cc%lai
    cc%mcv_dry = max(mcv_min, mcv_lai*cc%lai)
  enddo
    
end subroutine update_derived_vegn_data

! ============================================================================
! returns the profiles of uptake used in the 'LINEAR' uptake option
subroutine vegn_uptake_profile(vegn, dz, uptake_frac_max, vegn_uptake_term)
  type(vegn_tile_type), intent(in)  :: vegn
  real,                 intent(in)  :: dz(:)
  real,                 intent(out) :: uptake_frac_max(:)
  real,                 intent(out) :: vegn_uptake_term(:)

  call cohort_uptake_profile(vegn%cohorts(1), dz, uptake_frac_max, vegn_uptake_term)
end subroutine


! ============================================================================
subroutine vegn_root_properties (vegn, dz, VRL, K_r, r_r)
  type(vegn_tile_type), intent(in)  :: vegn 
  real,                 intent(in)  :: dz(:)
  real, intent(out) :: &
       vrl(:), & ! volumetric fine root length, m/m3
       K_r,    & ! root membrane permeability per unit area, kg/(m3 s)
       r_r       ! radius of fine roots, m

  call cohort_root_properties(vegn%cohorts(1), dz, VRL, K_r, r_r)
end subroutine 


! ============================================================================
function vegn_data_rs_min ( vegn )
  real :: vegn_data_rs_min
  type(vegn_tile_type), intent(in)  :: vegn
  
  vegn_data_rs_min = vegn%cohorts(1)%rs_min
end function


! ============================================================================
function vegn_seed_supply ( vegn )
  real :: vegn_seed_supply
  type(vegn_tile_type), intent(in) :: vegn

  ! ---- local vars 
  real :: vegn_bliving
  integer :: i
  
  vegn_bliving = 0
  do i = 1,vegn%n_cohorts
     vegn_bliving = vegn_bliving + vegn%cohorts(i)%bliving
  enddo
  vegn_seed_supply = MAX (vegn_bliving-BSEED, 0.0)
  
end function 

! ============================================================================
function vegn_seed_demand ( vegn )
  real :: vegn_seed_demand
  type(vegn_tile_type), intent(in) :: vegn

  integer :: i

  vegn_seed_demand = 0
  do i = 1,vegn%n_cohorts
     if(vegn%cohorts(i)%bliving<BSEED.and.vegn%t_ann>253.16.and.vegn%p_ann>1E-6) then
        vegn_seed_demand = vegn_seed_demand + BSEED
     endif
  enddo
end function 

! ============================================================================
subroutine vegn_add_bliving ( vegn, delta )
  type(vegn_tile_type), intent(inout) :: vegn
  real :: delta ! increment of bliving

  vegn%cohorts(1)%bliving = vegn%cohorts(1)%bliving + delta

  if (vegn%cohorts(1)%bliving < 0)then
     call error_mesg('vegn_add_bliving','resulting bliving is less then 0', FATAL)
  endif
  call update_biomass_pools(vegn%cohorts(1))
end subroutine 





! ============================================================================
! given a vegetation patch, destination kind of transition, and "transition 
! intensity" value, this function returns a fraction of tile that will parti-
! cipate in transition.
!
! this function must be contiguous, monotonic, its value must be within
! interval [0,1]
!
! this function is used to determine what part of each tile is to be converted
! to another land use kind; the equation is solved to get "transition intensity" 
! tau for which total area is equal to requested. Tau is, therefore, a dummy
! parameter, and only relative values of the priority functions for tiles 
! participating in transition have any meaning. For most transitions the priority 
! function is just equal to tau: therefore there is no preference, and all tiles
! contribute equally to converted area. For secondary vegetation harvesting, 
! however, priority also depends on wood biomass, and therefore tiles
! with high wood biomass are harvested first.
function vegn_tran_priority(vegn, dst_kind, tau) result(pri)
  real :: pri
  type(vegn_tile_type), intent(in) :: vegn
  integer             , intent(in) :: dst_kind
  real                , intent(in) :: tau

  real :: vegn_bwood
  integer :: i

  if (vegn%landuse==LU_SCND.and.dst_kind==LU_SCND) then ! secondary biomass harvesting
     vegn_bwood = 0
     do i = 1,vegn%n_cohorts
        vegn_bwood = vegn_bwood + vegn%cohorts(i)%bwood
     enddo
     pri = max(min(tau+vegn_bwood,1.0),0.0)
  else
     pri = max(min(tau,1.0),0.0)
  endif
end function 


! ============================================================================
function vegn_cover_cold_start(land_mask, lonb, latb) result (vegn_frac)
! creates and initializes a field of fractional vegn coverage
  logical, intent(in) :: land_mask(:,:)    ! land mask
  real,    intent(in) :: lonb(:,:), latb(:,:)! boundaries of the grid cells
  real,    pointer    :: vegn_frac (:,:,:) ! output: map of vegn fractional coverage

  allocate( vegn_frac(size(land_mask,1),size(land_mask,2),MSPECIES))

  call init_cover_field(vegn_to_use, 'INPUT/cover_type.nc', 'cover','frac', &
       lonb, latb, vegn_index_constant, input_cover_types, vegn_frac)
  
end function 

! =============================================================================
! returns true if tile fits the specified selector
function vegn_is_selected(vegn, sel)
  logical vegn_is_selected
  type(tile_selector_type),  intent(in) :: sel
  type(vegn_tile_type),      intent(in) :: vegn

  select case (sel%idata1)
  case (LU_SEL_TAG)
     vegn_is_selected = (sel%idata2 == vegn%landuse)
  case (SP_SEL_TAG)
     if (.not.associated(vegn%cohorts)) then
        vegn_is_selected = .FALSE.
     else
        vegn_is_selected = (sel%idata2 == vegn%cohorts(1)%species)
     endif
  case (NG_SEL_TAG)
     if (.not.associated(vegn%cohorts)) then
        vegn_is_selected = .FALSE.
     else
        vegn_is_selected = &
             ((vegn%cohorts(1)%species==SP_C4GRASS) .or.&
              (vegn%cohorts(1)%species==SP_C3GRASS)).and.&
             ((vegn%landuse==LU_NTRL).or. &
              (vegn%landuse==LU_SCND))
     endif
  case default
     vegn_is_selected = .FALSE.
  end select  
     
end function


! ============================================================================
! returns tag of the tile
function get_vegn_tile_tag(vegn) result(tag)
  integer :: tag
  type(vegn_tile_type), intent(in) :: vegn
  
  tag = vegn%tag
end function

! ============================================================================
! returns total wood biomass per tile 
function get_vegn_tile_bwood(vegn) result(bwood)
  real :: bwood
  type(vegn_tile_type), intent(in) :: vegn

  ! ---- local vars
  integer :: i

  bwood = 0
  do i = 1,vegn%n_cohorts
     bwood = bwood + vegn%cohorts(i)%bwood
  enddo
end function

! ============================================================================
subroutine vegn_tile_stock_pe (vegn, twd_liq, twd_sol  )
  type(vegn_tile_type),  intent(in)    :: vegn
  real,                  intent(out)   :: twd_liq, twd_sol
  integer n
  
  twd_liq = 0.
  twd_sol = 0.
  do n=1, vegn%n_cohorts
    twd_liq = twd_liq + vegn%cohorts(n)%prog%wl
    twd_sol = twd_sol + vegn%cohorts(n)%prog%ws
!      vegn_HEAT  = (mcv + clw*cohort%prog%Wl+ csw*cohort%prog%Ws)*(cohort%prog%Tv-tfreeze)

    enddo
end subroutine vegn_tile_stock_pe


! ============================================================================
! returns total carbon in the tile, kg C/m2
function vegn_tile_carbon(vegn) result(carbon) ; real carbon
  type(vegn_tile_type), intent(in)  :: vegn

  integer :: i

  carbon = 0
  do i = 1,vegn%n_cohorts
     carbon = carbon + &
          vegn%cohorts(i)%bl + vegn%cohorts(i)%blv + &
          vegn%cohorts(i)%br + vegn%cohorts(i)%bwood + &
          vegn%cohorts(i)%bsw + &
          vegn%cohorts(i)%carbon_gain + vegn%cohorts(i)%bwood_gain
  enddo
  carbon = carbon + &
       sum(vegn%harv_pool) + vegn%fsc_pool + vegn%ssc_pool + vegn%csmoke_pool
end function


! ============================================================================
! returns heat content of the vegetation, J/m2
function vegn_tile_heat (vegn) result(heat) ; real heat
  type(vegn_tile_type), intent(in)  :: vegn

  integer :: i

  heat = 0
  do i = 1, vegn%n_cohorts
     heat = heat + &
          (clw*vegn%cohorts(i)%prog%Wl + &
             csw*vegn%cohorts(i)%prog%Ws + &
             vegn%cohorts(i)%mcv_dry)*(vegn%cohorts(i)%prog%Tv-tfreeze) - &
           hlf*vegn%cohorts(i)%prog%Ws
  enddo
end function

end module vegn_tile_mod

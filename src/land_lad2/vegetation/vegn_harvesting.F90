module vegn_harvesting_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : write_version_number, string, error_mesg, FATAL, NOTE, &
     mpp_pe, write_version_number, file_exist, close_file, &
     check_nml_error, stdlog, mpp_root_pe
use mpp_io_mod, only : axistype, mpp_get_atts, mpp_get_axis_data, &
     mpp_open, mpp_close, MPP_RDONLY, MPP_WRONLY, MPP_ASCII
use vegn_data_mod, only : &
     N_LU_TYPES, LU_PAST, LU_CROP, LU_NTRL, LU_SCND, &
     HARV_POOL_PAST, HARV_POOL_CROP, HARV_POOL_CLEARED, HARV_POOL_WOOD_FAST, &
     HARV_POOL_WOOD_MED, HARV_POOL_WOOD_SLOW, &
     agf_bs, fsc_liv, fsc_wood
use vegn_tile_mod, only : &
     vegn_tile_type
use vegn_cohort_mod, only : &
     vegn_cohort_type, update_biomass_pools

implicit none
private

! ==== public interface ======================================================
public :: vegn_harvesting_init
public :: vegn_harvesting_end

public :: vegn_harvesting

public :: vegn_graze_pasture
public :: vegn_harvest_cropland
public :: vegn_cut_forest
! ==== end of public interface ===============================================

! ==== module constants =====================================================
character(len=*), parameter   :: &
     version = '$Id: vegn_harvesting.F90,v 20.0 2013/12/13 23:31:12 fms Exp $', &
     tagname = '$Name: tikal $', &
     module_name = 'vegn_harvesting_mod'
real, parameter :: ONETHIRD = 1.0/3.0

! ==== module data ==========================================================

! ---- namelist variables ---------------------------------------------------
logical :: do_harvesting       = .TRUE.  ! if true, then harvesting of crops and pastures is done
real :: grazing_intensity      = 0.25    ! fraction of biomass removed each time by grazing
real :: grazing_residue        = 0.1     ! fraction of the grazed biomass transferred into soil pools
real :: frac_wood_wasted_harv  = 0.25    ! fraction of wood wasted while harvesting
real :: frac_wood_wasted_clear = 0.25    ! fraction of wood wasted while clearing land for pastures or crops
logical :: waste_below_ground_wood = .TRUE. ! If true, all the wood below ground (1-agf_bs fraction of bwood 
        ! and bsw) is wasted. Old behavior assumed this to be FALSE.
real :: frac_wood_fast         = ONETHIRD ! fraction of wood consumed fast
real :: frac_wood_med          = ONETHIRD ! fraction of wood consumed with medium speed
real :: frac_wood_slow         = ONETHIRD ! fraction of wood consumed slowly
real :: crop_seed_density      = 0.1     ! biomass of seeds left after crop harvesting, kg/m2
namelist/harvesting_nml/ do_harvesting, grazing_intensity, grazing_residue, &
     frac_wood_wasted_harv, frac_wood_wasted_clear, waste_below_ground_wood, &
     frac_wood_fast, frac_wood_med, frac_wood_slow, &
     crop_seed_density

contains ! ###################################################################

! ============================================================================
subroutine vegn_harvesting_init
  integer :: unit, ierr, io

  call write_version_number(version, tagname)

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=harvesting_nml, iostat=io)
  ierr = check_nml_error(io, 'harvesting_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=harvesting_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'harvesting_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=harvesting_nml)
  endif

  if (frac_wood_fast+frac_wood_med+frac_wood_slow/=1.0) then
     call error_mesg('vegn_harvesting_init', &
          'sum of frac_wood_fast, frac_wood_med, and frac_wood_slow must be 1.0',&
          FATAL)
  endif
end subroutine vegn_harvesting_init


! ============================================================================
subroutine vegn_harvesting_end
end subroutine vegn_harvesting_end


! ============================================================================
! harvest vegetation in a tile
subroutine vegn_harvesting(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  if (.not.do_harvesting) &
       return ! do nothing if no harvesting requested

  select case(vegn%landuse)
  case(LU_PAST)  ! pasture
     call vegn_graze_pasture    (vegn)
  case(LU_CROP)  ! crop
     call vegn_harvest_cropland (vegn)
  end select
end subroutine


! ============================================================================
subroutine vegn_graze_pasture(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars 
  real ::  bdead0, balive0, btotal0; ! initial combined biomass pools
  real ::  bdead1, balive1, btotal1; ! updated combined biomass pools
  type(vegn_cohort_type), pointer :: cc ! shorthand for the current cohort
  integer :: i

  balive0 = 0 ; bdead0 = 0 ;
  balive1 = 0 ; bdead1 = 0 ;

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1,vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! calculate total biomass pools for the patch
     balive0 = balive0 + cc%bl + cc%blv + cc%br
     bdead0  = bdead0  + cc%bwood + cc%bsw
     ! only potential leaves are consumed
     vegn%harv_pool(HARV_POOL_PAST) = vegn%harv_pool(HARV_POOL_PAST) + &
          cc%bliving*cc%Pl*grazing_intensity*(1-grazing_residue) ;
     cc%bliving = cc%bliving - cc%bliving*cc%Pl*grazing_intensity;

     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(cc);
 
     ! calculate new combined vegetation biomass pools
     balive1 = balive1 + cc%bl + cc%blv + cc%br
     bdead1  = bdead1  + cc%bwood + cc%bsw
  enddo
  btotal0 = balive0 + bdead0
  btotal1 = balive1 + bdead1

  ! update intermediate soil carbon pools
  vegn%fsc_pool = vegn%fsc_pool + &
       (fsc_liv*(balive0-balive1)+fsc_wood*(bdead0-bdead1))*grazing_residue;
  vegn%ssc_pool = vegn%ssc_pool + &
       ((1-fsc_liv)*(balive0-balive1)+ (1-fsc_wood)*(bdead0-bdead1))*grazing_residue;
end subroutine vegn_graze_pasture


! ================================================================================
subroutine vegn_harvest_cropland(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc ! pointer to the current cohort
  real :: fraction_harvested;    ! fraction of biomass harvested this time
  real :: bdead, balive, btotal; ! combined biomass pools
  integer :: i
  
  balive = 0 ; bdead = 0
  ! calculate initial combined biomass pools for the patch
  do i = 1, vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! calculate total biomass pools for the patch
     balive = balive + cc%bl + cc%blv + cc%br
     bdead  = bdead  + cc%bwood + cc%bsw
  enddo
  btotal = balive+bdead;

  ! calculate harvested fraction: cut everything down to seed level
  fraction_harvested = MIN(MAX((btotal-crop_seed_density)/btotal,0.0),1.0);

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1, vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! use for harvest only aboveg round living biomass and waste the correspondent below living and wood
     vegn%harv_pool(HARV_POOL_CROP) = vegn%harv_pool(HARV_POOL_CROP) + &
          cc%bliving*(cc%Pl + cc%Psw*agf_bs)*fraction_harvested;
     vegn%fsc_pool = vegn%fsc_pool + fraction_harvested*(fsc_liv*cc%bliving*cc%Pr + &
          fsc_wood*(cc%bwood + cc%bliving*cc%Psw*(1-agf_bs)));
     vegn%ssc_pool = vegn%ssc_pool + fraction_harvested*((1-fsc_liv)*cc%bliving*cc%Pr + &
          (1-fsc_wood)*(cc%bwood + cc%bliving*cc%Psw*(1-agf_bs)));

     cc%bliving = cc%bliving * (1-fraction_harvested);
     cc%bwood   = cc%bwood   * (1-fraction_harvested);
     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(cc);
  enddo
end subroutine vegn_harvest_cropland


! ============================================================================
! for now cutting forest is the same as harvesting cropland --
! we basically cut down everything, leaving only seeds
subroutine vegn_cut_forest(vegn, new_landuse)
  type(vegn_tile_type), intent(inout) :: vegn
  integer, intent(in) :: new_landuse ! new land use type that gets assigned to
                                     ! the tile after the wood harvesting

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc ! pointer to the current cohort
  real :: frac_harvested;        ! fraction of biomass harvested this time
  real :: frac_wood_wasted       ! fraction of wood wasted during transition
  real :: wood_harvested         ! anount of harvested wood, kgC/m2
  real :: bdead, balive, btotal; ! combined biomass pools
  real :: delta
  integer :: i
  
  balive = 0 ; bdead = 0
  ! calculate initial combined biomass pools for the patch
  do i = 1, vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! calculate total biomass pools for the patch
     balive = balive + cc%bl + cc%blv + cc%br
     bdead  = bdead  + cc%bwood + cc%bsw
  enddo
  btotal = balive+bdead;

  ! calculate harvested fraction: cut everything down to seed level
  frac_harvested = MIN(MAX((btotal-crop_seed_density)/btotal,0.0),1.0);

  ! define fraction of wood wasted, based on the transition type
  if (new_landuse==LU_SCND) then
     frac_wood_wasted = frac_wood_wasted_harv
  else
     frac_wood_wasted = frac_wood_wasted_clear
  endif
  ! take into accont that all wood below ground is wasted; also the fraction
  ! of waste calculated above is lost from the above-ground part of the wood
  if (waste_below_ground_wood) then
     frac_wood_wasted = (1-agf_bs) + agf_bs*frac_wood_wasted
  endif

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)

     ! calculate total amount of harvested wood, minus the wasted part
     wood_harvested = (cc%bwood+cc%bsw)*frac_harvested*(1-frac_wood_wasted)

     ! distribute harvested wood between pools
     if (new_landuse==LU_SCND) then
        ! this is harvesting, distribute between 3 different wood pools 
        vegn%harv_pool(HARV_POOL_WOOD_FAST) = vegn%harv_pool(HARV_POOL_WOOD_FAST) &
             + wood_harvested*frac_wood_fast
        vegn%harv_pool(HARV_POOL_WOOD_MED) = vegn%harv_pool(HARV_POOL_WOOD_MED) &
             + wood_harvested*frac_wood_med
        vegn%harv_pool(HARV_POOL_WOOD_SLOW) = vegn%harv_pool(HARV_POOL_WOOD_SLOW) &
             + wood_harvested*frac_wood_slow
     else
        ! this is land clearance: everything goes into "cleared" pool
        vegn%harv_pool(HARV_POOL_CLEARED) = vegn%harv_pool(HARV_POOL_CLEARED) &
             + wood_harvested
     endif

     ! distribute wood and living biomass between fast and slow intermediate 
     ! soil carbon pools according to fractions specified thorough the namelists
     delta = (cc%bwood+cc%bsw)*frac_harvested*frac_wood_wasted;
     if(delta<0) call error_mesg('vegn_cut_forest', &
          'harvested amount of dead biomass ('//string(delta)//' kgC/m2) is below zero', &
          FATAL)
     vegn%ssc_pool = vegn%ssc_pool + delta*(1-fsc_wood);
     vegn%fsc_pool = vegn%fsc_pool + delta*   fsc_wood ;

     delta = balive * frac_harvested;
     if(delta<0) call error_mesg('vegn_cut_forest', &
          'harvested amount of live biomass ('//string(delta)//' kgC/m2) is below zero', &
          FATAL)
     vegn%ssc_pool = vegn%ssc_pool + delta*(1-fsc_liv) ;
     vegn%fsc_pool = vegn%fsc_pool + delta*   fsc_liv  ;

     cc%bliving = cc%bliving*(1-frac_harvested);
     cc%bwood   = cc%bwood*(1-frac_harvested);
     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(cc);
  enddo
end subroutine vegn_cut_forest

end module 

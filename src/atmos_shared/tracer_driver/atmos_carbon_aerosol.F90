Module atmos_carbon_aerosol_mod

! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Shekar Reddy
! </CONTACT>
use mpp_mod, only: input_nml_file 
use fms_mod,                    only : file_exist, close_file, &
                                       write_version_number, &
                                       mpp_pe, mpp_root_pE, &
                                       open_namelist_file,  &
                                       check_nml_error, error_mesg,  &
                                       stdlog, FATAL, NOTE, WARNING
use time_manager_mod,           only : time_type, &
                                       days_in_month, days_in_year, &
                                       set_date, set_time, get_date_julian, &
                                       print_date, get_date, &
                                       operator(>), operator(+), operator(-)
use time_interp_mod,            only:  fraction_of_year, &
                                       time_interp_init
use diag_manager_mod,           only : send_data, register_diag_field, &
                                       diag_manager_init, get_base_time
use tracer_manager_mod,         only : get_tracer_index, &
                                       set_tracer_atts
use field_manager_mod,          only : MODEL_ATMOS
use atmos_tracer_utilities_mod, only : wet_deposition,       &
                                       dry_deposition
use interpolator_mod,           only:  interpolate_type, interpolator_init, &
                                       obtain_interpolator_time_slices,&
                                       unset_interpolator_time_flag, &
                                       interpolator, interpolator_end, &
                                       CONSTANT, INTERP_WEIGHTED_P
use constants_mod,              only : PI, GRAV, RDGAS
implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_carbon_aerosol_driver,   &
        atmos_carbon_aerosol_init, &
        atmos_carbon_aerosol_time_vary, &
        atmos_carbon_aerosol_endts, &
        atmos_carbon_aerosol_end

!-----------------------------------------------------------------------
! tracer number for carbonaceous aerosols
integer :: nbcphobic=0
integer :: nbcphilic=0
integer :: nomphobic=0
integer :: nomphilic=0

!--- identification numbers for  diagnostic fields and axes ----
integer :: id_bcphob_emis, id_bcphil_emis, id_omphob_emis, id_omphil_emis
integer :: id_om_emis_col, id_bc_emis_col
integer :: id_om_emis_colv2, id_bc_emis_colv2
integer :: id_bcphob_sink, id_omphob_sink
integer :: id_emisbb, id_omemisbb_col
integer :: id_bcemisbf, id_bcemisbb, id_bcemissh, id_bcemisff, id_bcemisav
integer :: id_omemisbf, id_omemisbb, id_omemissh, id_omemisff, id_omemisbg, id_omemisoc
integer :: id_bc_tau
!----------------------------------------------------------------------
!--- Interpolate_type variable containing all the information needed to
! interpolate the emission provided in the netcdf input file.
type(interpolate_type),save  ::bcff_aerosol_interp
type(interpolate_type),save  ::bcbb_aerosol_interp
type(interpolate_type),save  ::bcbf_aerosol_interp
type(interpolate_type),save  ::bcsh_aerosol_interp
type(interpolate_type),save  ::bcav_aerosol_interp
type(interpolate_type),save  ::omff_aerosol_interp
type(interpolate_type),save  ::ombb_aerosol_interp
type(interpolate_type),save  ::ombf_aerosol_interp
type(interpolate_type),save  ::omsh_aerosol_interp
type(interpolate_type),save  ::omna_aerosol_interp
type(interpolate_type),save  ::omss_aerosol_interp
! Initial calendar time for model
type(time_type) :: model_init_time

! Difference between model initial time and source timeseries applied
! at model initial time
type(time_type), save :: bcff_offset
type(time_type), save :: bcbb_offset
type(time_type), save :: bcbf_offset
type(time_type), save :: bcsh_offset
type(time_type), save :: bcav_offset
type(time_type), save :: omff_offset
type(time_type), save :: ombb_offset
type(time_type), save :: ombf_offset
type(time_type), save :: omsh_offset
type(time_type), save :: omna_offset
type(time_type), save :: omss_offset

! timeseries which is mapped to model initial time
type(time_type), save :: bcff_entry
type(time_type), save :: bcbb_entry
type(time_type), save :: bcbf_entry
type(time_type), save :: bcsh_entry
type(time_type), save :: bcav_entry
type(time_type), save :: omff_entry
type(time_type), save :: ombb_entry
type(time_type), save :: ombf_entry
type(time_type), save :: omsh_entry
type(time_type), save :: omna_entry
type(time_type), save :: omss_entry

! The model initial time is later than the XXX_dataset_entry time  ?
logical, save    :: bcff_negative_offset
logical, save    :: bcbb_negative_offset
logical, save    :: bcbf_negative_offset
logical, save    :: bcsh_negative_offset
logical, save    :: bcav_negative_offset
logical, save    :: omff_negative_offset
logical, save    :: ombb_negative_offset
logical, save    :: ombf_negative_offset
logical, save    :: omsh_negative_offset
logical, save    :: omna_negative_offset
logical, save    :: omss_negative_offset

integer, save    :: bcff_time_serie_type
integer, save    :: bcbb_time_serie_type
integer, save    :: bcbf_time_serie_type
integer, save    :: bcsh_time_serie_type
integer, save    :: bcav_time_serie_type
integer, save    :: omff_time_serie_type
integer, save    :: ombb_time_serie_type
integer, save    :: ombf_time_serie_type
integer, save    :: omsh_time_serie_type
integer, save    :: omna_time_serie_type
integer, save    :: omss_time_serie_type
!----------------------------------------------------------------------
!-------- namelist  ---------
character(len=80) :: bcff_filename = 'carbon_aerosol_emission.nc'
character(len=80) :: bcbb_filename = 'carbon_aerosol_emission.nc'
character(len=80) :: bcbf_filename = 'carbon_aerosol_emission.nc'
character(len=80) :: bcsh_filename = 'carbon_aerosol_emission.nc'
character(len=80) :: bcav_filename = 'carbon_aerosol_emission.nc'
character(len=80) :: omff_filename = 'carbon_aerosol_emission.nc'
character(len=80) :: ombb_filename = 'carbon_aerosol_emission.nc'
character(len=80) :: ombf_filename = 'carbon_aerosol_emission.nc'
character(len=80) :: omsh_filename = 'carbon_aerosol_emission.nc'
character(len=80) :: omna_filename = 'carbon_aerosol_emission.nc'
character(len=80) :: omss_filename = 'gocart_emission.nc'

integer :: i
character(len=80), save, dimension(1) :: bcff_emission_name = (/' '/)
character(len=80), save, dimension(6) :: bcbb_emission_name = (/(' ',i=1,6)/)
character(len=80), save, dimension(1) :: bcbf_emission_name = (/' '/)
character(len=80), save, dimension(1) :: bcsh_emission_name = (/' '/)
character(len=80), save, dimension(1) :: bcav_emission_name = (/' '/)
character(len=80), save, dimension(1) :: omff_emission_name = (/' '/)
character(len=80), save, dimension(6) :: ombb_emission_name = (/(' ',i=1,6)/)
character(len=80), save, dimension(1) :: ombf_emission_name = (/' '/)
character(len=80), save, dimension(1) :: omsh_emission_name = (/' '/)
character(len=80), save, dimension(1) :: omna_emission_name = (/' '/)
character(len=80), save, dimension(1) :: omss_emission_name = (/' '/)

character(len=80), dimension(1) :: bcff_input_name = (/' '/)
character(len=80), dimension(6) :: bcbb_input_name = (/(' ',i=1,6)/)
character(len=80), dimension(1) :: bcbf_input_name = (/' '/)
character(len=80), dimension(1) :: bcsh_input_name = (/' '/)
character(len=80), dimension(1) :: bcav_input_name = (/' '/)
character(len=80), dimension(1) :: omff_input_name = (/' '/)
character(len=80), dimension(6) :: ombb_input_name = (/(' ',i=1,6)/)
character(len=80), dimension(1) :: ombf_input_name = (/' '/)
character(len=80), dimension(1) :: omsh_input_name = (/' '/)
character(len=80), dimension(1) :: omna_input_name = (/' '/)
character(len=80), dimension(1) :: omss_input_name = (/' '/)
! Default values for carbon_aerosol_nml
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FOSSIL FUEL source can be either:
! for black carbon: 'cooke_and_wilson_1996' 
!                   'cooke_1999' 
!                   'bond_2004'
character(len=80)     :: bcff_source = ' '
character(len=80)     :: bcff_time_dependency_type = 'constant'
integer, dimension(6) :: bcff_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fossil fuel emission 
! for organic matter: 'cooke_1999' 
!                     'bond_2004'
character(len=80)     :: omff_source = ' '
character(len=80)     :: omff_time_dependency_type = 'constant'
integer, dimension(6) :: omff_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BIOMASS BURNING emissions can be either
! for black carbon: 'cooke_and_wilson_1996' 
!                   'bond_2004' 
!                   'GEIA level 1 and 2'
!                   'AEROCOM level 1 to 6'
character(len=80)     :: bcbb_source = ' '
character(len=80)     :: bcbb_time_dependency_type = 'constant'
integer, dimension(6) :: bcbb_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
! open biomass burning emission
! for organic matter: 'cooke_and_wilson_1996' 
!                     'bond_2004' 
!                     'GEIA level 1 and 2'
!                     'AEROCOM level 1 to 6'
character(len=80)     :: ombb_source = ' '
character(len=80)     :: ombb_time_dependency_type = 'constant'
integer, dimension(6) :: ombb_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! biofuel emission        
! for black carbon     'bond_2004' 
character(len=80)     :: bcbf_source = ' '
character(len=80)     :: bcbf_time_dependency_type = 'constant'
integer, dimension(6) :: bcbf_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ship emissions based on V. Eyrig
character(len=80)     :: bcsh_source = ' '
character(len=80)     :: bcsh_time_dependency_type
integer, dimension(6) :: bcsh_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! aircraft emissions based on Steven L. Baughcum
character(len=80)     :: bcav_source = ' '
character(len=80)     :: bcav_time_dependency_type = 'constant'
integer, dimension(6) :: bcav_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
! Emission index from Henricks et al., Simulating the global atmospheric
! black carbon cycle: a revisit to the contribution of aircraft emissions,
! Atmos. Chem. Phys., 4, 252102541, 2004.
real :: bc_aircraft_EI = 4.e-5   ! kg BC /kg fuel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! biofuel emission for organic matter:  'bond_2004' 
character(len=80)     :: ombf_source = ' '
character(len=80)     :: ombf_time_dependency_type = 'constant'
integer, dimension(6) :: ombf_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ship emissions based on V. Eyrig
character(len=80)     :: omsh_source = ' '
character(len=80)     :: omsh_time_dependency_type = 'constant'
integer, dimension(6) :: omsh_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for natural emssion: Gunther et al., 1996 inventory
character(len=80)     :: omna_source = ' '
character(len=80)     :: omna_time_dependency_type = 'constant'
integer, dimension(6) :: omna_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for sea spray emssion: based on ODowd et al., A combined organic-inorganic sea-spray source function, Geophys. Res. Lett., v35, L01801, doi:10.1029/2007GL030331, 2008
character(len=80)     :: omss_source = ' '
character(len=80)     :: omss_time_dependency_type = 'constant'
integer, dimension(6) :: omss_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
real, save :: coef_omss_emis
real :: omss_coef=-999.
logical               :: do_dynamic_bc = .false.
real                  :: bcage = 1.0
real                  :: bcageslow = 25.
logical               :: do_dynamic_om = .false.
real                  :: omage = 0.5
real                  :: omageslow = 20.  !days
real                  :: frac_bc_phobic = 0.8
real                  :: frac_bc_philic = 0.2
real                  :: frac_bcbb_phobic = 0.8
real                  :: frac_bcbb_philic = 0.2
real                  :: frac_om_phobic = 0.5
real                  :: frac_om_philic = 0.5
!!!!!!!!!!!!!!!!!!!!!!!!!!
namelist /carbon_aerosol_nml/ &
 bcff_source, bcff_input_name, bcff_filename, &
  bcff_time_dependency_type, bcff_dataset_entry, &
 bcbb_source, bcbb_input_name, bcbb_filename, &
  bcbb_time_dependency_type, bcbb_dataset_entry, &
 bcbf_source, bcbf_input_name, bcbf_filename, &
  bcbf_time_dependency_type, bcbf_dataset_entry, &
 bcsh_source, bcsh_input_name, bcsh_filename, &
  bcsh_time_dependency_type, bcsh_dataset_entry, &
 bcav_source, bcav_input_name, bcav_filename, &
  bcav_time_dependency_type, bcav_dataset_entry, bc_aircraft_EI, &
 omff_source, omff_input_name, omff_filename, &
  omff_time_dependency_type, omff_dataset_entry, &
 ombb_source, ombb_input_name, ombb_filename, &
  ombb_time_dependency_type, ombb_dataset_entry, &
 ombf_source, ombf_input_name, ombf_filename, &
  ombf_time_dependency_type, ombf_dataset_entry, &
 omsh_source, omsh_input_name, omsh_filename, &
  omsh_time_dependency_type, omsh_dataset_entry, &
 omna_source, omna_input_name, omna_filename, &
  omna_time_dependency_type, omna_dataset_entry, &
 omss_source, omss_input_name, omss_filename, &
  omss_time_dependency_type, omss_dataset_entry, omss_coef, &
 do_dynamic_bc, bcage, bcageslow,            &
 do_dynamic_om, omage, omageslow,            &
 frac_bc_phobic, frac_bc_philic, frac_om_phobic, frac_om_philic, &
 frac_bcbb_philic, frac_bcbb_phobic

character(len=6), parameter :: module_name = 'tracer'

logical :: module_is_initialized = .FALSE.
logical :: used

!---- version number -----
character(len=128) :: version = '$Id: atmos_carbon_aerosol.F90,v 20.0 2013/12/13 23:23:44 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'
!-----------------------------------------------------------------------

type(time_type)                        :: bcff_time
type(time_type)                        :: bcbb_time
type(time_type)                        :: bcbf_time
type(time_type)                        :: bcsh_time
type(time_type)                        :: bcav_time
type(time_type)                        :: omff_time
type(time_type)                        :: ombb_time
type(time_type)                        :: ombf_time
type(time_type)                        :: omsh_time
type(time_type)                        :: omna_time
type(time_type)                        :: omss_time



contains

!#######################################################################

subroutine atmos_carbon_aerosol_driver(lon, lat, ocn_flx_fraction,  &
                               pfull,phalf, &
                               z_half, z_pbl, &
                               t_surf, w10m, &
                               T, pwt, &
                               bcphob, bcphob_dt,  &
                               bcphil, bcphil_dt,  &
                               omphob, omphob_dt, &
                               omphil, omphil_dt, &
                               oh_conc,&
                               diag_time, is, ie, js, je )

!-----------------------------------------------------------------------
   real, intent(in),  dimension(:,:)   :: lon, lat
   real, intent(in),  dimension(:,:)   :: ocn_flx_fraction
   real, intent(in),  dimension(:,:)   :: z_pbl
   real, intent(in),  dimension(:,:,:) :: z_half
   real, intent(in),  dimension(:,:)   :: w10m, t_surf  ! ocean sea surface temperature and 10 meter wind speed
   real, intent(in),  dimension(:,:,:) :: pwt,pfull,phalf,T
   real, intent(in),  dimension(:,:,:) :: bcphob,bcphil
   real, intent(in),  dimension(:,:,:) :: omphob,omphil
   real, intent(in),  dimension(:,:,:) :: oh_conc
   real, intent(out), dimension(:,:,:) :: bcphob_dt,bcphil_dt
   real, intent(out), dimension(:,:,:) :: omphob_dt,omphil_dt
type(time_type), intent(in)            :: diag_time
integer, intent(in)                    :: is, ie, js, je
!-----------------------------------------------------------------------

real  dtr,bltop,z1,z2,del
real, dimension(size(bcphob,3)) :: fa1, fa2
real, dimension(size(bcphob,3),6) :: fbb
integer :: lf, nlevel_fire
real, dimension(6) :: alt_fire_min, alt_fire_max
! Lower altitude of injection from wild fires 
! These values correspond to the AEROCOM input data (cf. Dentener, ACP, 2006)
integer, parameter :: nlevel_fire_AEROCOM = 6
! Modified GEIA dataset proposed by Reddy and Boucher, J. Geophys. Res., 2004
integer, parameter :: nlevel_fire_GEIA    = 2
real, dimension(nlevel_fire_AEROCOM) :: &
      alt_fire_min_AEROCOM=(/0.,100.,500.,1000.,2000.,3000./)
real, dimension(nlevel_fire_GEIA) :: &
      alt_fire_min_GEIA=(/0.,300./)
! Upper altitude of injection from wild fires 
real, dimension(nlevel_fire_AEROCOM) :: &
      alt_fire_max_AEROCOM=(/100.,500.,1000.,2000.,3000.,6000./)
real, dimension(nlevel_fire_GEIA) :: &
      alt_fire_max_GEIA=(/300.,1500./)
real :: ze1 = 100.
real :: ze2 = 300.
integer        :: i, j,l, id,jd,kd,k
real  :: sst, Schm, SchmCO2, AKw

!
!-------------------------------------------------------------------------
real, dimension(size(bcphob,1),size(bcphob,2)) :: bcemisff_l1, bcemisff_l2
real, dimension(size(omphob,1),size(omphob,2)) :: omemisff_l1, omemisff_l2
real, dimension(size(bcphob,1),size(bcphob,2),6) :: bcemisbb
real, dimension(size(omphob,1),size(omphob,2),6) :: omemisbb

real,dimension(size(bcphob,1),size(bcphob,2)) ::  emisob, omemisob_2d
real,dimension(size(bcphob,1),size(bcphob,2),size(bcphob,3)) ::&
   bcphob_emis, bcphil_emis, bcphob_sink, bcemisob, bcemisff, bcemisav
real,dimension(size(bcphob,1),size(bcphob,2),size(bcphob,3)) :: bc_tau
real,dimension(size(omphob,1),size(omphob,2),size(omphob,3)) ::&
   omphob_emis, omphil_emis, omphob_sink, omemisob, omemisff
real, dimension(size(bcphob,1),size(bcphob,2)) ::          &
   bcemisbf, bcemissh
real, dimension(size(omphob,1),size(omphob,2)) ::          &
   omemisbf, omemissh, omemisbg, dmso, omemisocean
real,dimension(size(bcphob,1),size(bcphob,2)) :: bc_emis, om_emis
real,dimension(size(bcphob,1),size(bcphob,2),size(bcphob,3)) ::zzz1,  &
                                                               oh_conc1
!
!-----------------------------------------------------------------------
!
      id=size(bcphob,1); jd=size(bcphob,2); kd=size(bcphob,3)

      dtr= PI/180.

!----------- compute black carbon source ------------

    bcphob_dt = 0.0
    bcphil_dt = 0.0
    omphob_dt = 0.0
    omphil_dt = 0.0
! emission rates (kg/m2/s)
    bcemisff_l1(:,:) = 0.0
    bcemisff_l2(:,:) = 0.0
    bcemisff(:,:,:)  = 0.0
    bcemisbb(:,:,:)  = 0.0
    bcemisob(:,:,:)  = 0.0
    emisob(:,:) = 0.
    omemisob_2d(:,:) = 0.
    bcemisbf(:,:)    = 0.0
    bcemissh(:,:)    = 0.0
    bcemisav(:,:,:)  = 0.0

    bcphob_emis(:,:,:) = 0.0
    bcphil_emis(:,:,:) = 0.0
    bcphob_sink(:,:,:) = 0.0
    omphob_emis(:,:,:) = 0.0
    omphil_emis(:,:,:) = 0.0
    omphob_sink(:,:,:) = 0.0
    bc_tau(:,:,:)      = 0.0

    omemisff_l1(:,:) = 0.0
    omemisff_l2(:,:) = 0.0
    omemisff(:,:,:)  = 0.0
    omemisbb(:,:,:)  = 0.0
    omemisob(:,:,:)  = 0.0
    omemisbf(:,:)    = 0.0
    omemissh(:,:)    = 0.0
    omemisbg(:,:)    = 0.0
    dmso(:,:)        = 0.0
    omemisocean(:,:) = 0.0

    if ( trim(bcff_source) .ne. ' ') then
     call interpolator(bcff_aerosol_interp, bcff_time, bcemisff_l1, &
          trim(bcff_emission_name(1)), is, js)
   endif

   if ( trim(bcbb_source).ne.' ') then
    nlevel_fire = 1
    alt_fire_min(:) = 0.0
    alt_fire_max(:) = 0.0
    select case (trim(bcbb_source))
      case ('cooke_and_wilson_1996')
        call interpolator(bcbb_aerosol_interp, bcbb_time, bcemisbb(:,:,1), &
                         trim(bcbb_emission_name(1)), is, js)
      case ('bond_2004') 
        call interpolator(bcbb_aerosol_interp, bcbb_time, bcemisbb(:,:,1), &
                         trim(bcbb_emission_name(1)), is, js)
      case ('gocart_2007') 
        call interpolator(bcbb_aerosol_interp, bcbb_time, bcemisbb(:,:,1), &
                         trim(bcbb_emission_name(1)), is, js)
      case ('RETRO') 
        call interpolator(bcbb_aerosol_interp, bcbb_time, bcemisbb(:,:,1), &
                         trim(bcbb_emission_name(1)), is, js)
      case ('GEIA')
        nlevel_fire = nlevel_fire_GEIA
        alt_fire_min(1:nlevel_fire_GEIA) = alt_fire_min_GEIA(1:nlevel_fire_GEIA)
        alt_fire_max(1:nlevel_fire_GEIA) = alt_fire_max_GEIA(1:nlevel_fire_GEIA)
        do lf=1, nlevel_fire 
          call interpolator(bcbb_aerosol_interp, bcbb_time, bcemisbb(:,:,lf), &
                         trim(bcbb_emission_name(lf)), is, js)
        enddo
      case ('AEROCOM') 
! Wildfire emissions at 6 levels from 0 to 6 km
! (cf. AEROCOM web site or Dentener et al., ACPD, 2006)
        nlevel_fire = nlevel_fire_AEROCOM
        alt_fire_min(1:nlevel_fire_AEROCOM) = &
                 alt_fire_min_AEROCOM(1:nlevel_fire_AEROCOM)
        alt_fire_max(1:nlevel_fire_AEROCOM) = &
                 alt_fire_max_AEROCOM(1:nlevel_fire_AEROCOM)
        do lf=1, nlevel_fire 
          call interpolator(bcbb_aerosol_interp, bcbb_time, bcemisbb(:,:,lf), &
                        trim(bcbb_emission_name(lf)), is, js)
        enddo
    end select
   endif
   if ( trim(bcbf_source).ne. ' ') then
     call interpolator(bcbf_aerosol_interp, bcbf_time, bcemisbf, &
                       trim(bcbf_emission_name(1)), is, js)
   endif

   if ( trim(bcsh_source).ne. ' ') then
     call interpolator(bcsh_aerosol_interp, bcsh_time, bcemissh, &
                  trim(bcsh_emission_name(1)), is, js)
   endif
   if ( trim(bcav_source).ne. ' ') then
     call interpolator(bcav_aerosol_interp, bcav_time, phalf, bcemisav, &
                         trim(bcav_emission_name(1)), is, js)
   endif
   if ( trim(omff_source).ne. ' ') then
     call interpolator(omff_aerosol_interp, omff_time, omemisff_l1, &
                      trim(omff_emission_name(1)), is, js)
   endif
   if ( trim(ombb_source).ne. ' ') then
    omemisbb(:,:,:) = 0.0
    nlevel_fire = 1
    alt_fire_min(:) = 0.0
    alt_fire_max(:) = 0.0
    select case (trim(ombb_source))
      case ('cooke_and_wilson_1996')
        call interpolator(ombb_aerosol_interp, ombb_time, omemisbb(:,:,1), &
                           trim(ombb_emission_name(1)), is, js)
        omemisbb(:,:,1) = omemisbb(:,:,1)*7.0
      case ('bond_2004') 
        call interpolator(ombb_aerosol_interp, ombb_time, omemisbb(:,:,1), &
                       trim(ombb_emission_name(1)), is, js)
      case ('gocart_2007')
        call interpolator(ombb_aerosol_interp, ombb_time, omemisbb(:,:,1), &
                       trim(ombb_emission_name(1)), is, js)
      case ('RETRO')
        call interpolator(ombb_aerosol_interp, ombb_time, omemisbb(:,:,1), &
                       trim(ombb_emission_name(1)), is, js)
      case ('GEIA')
        nlevel_fire = nlevel_fire_GEIA
        alt_fire_min(1:nlevel_fire_GEIA) = alt_fire_min_GEIA(1:nlevel_fire_GEIA)
        alt_fire_max(1:nlevel_fire_GEIA) = alt_fire_max_GEIA(1:nlevel_fire_GEIA)
! GEIA emission inventory scaled by ATSR fire counts (Reddy and Boucher, 2004)
! There are 2 levels of emission
        do lf=1, nlevel_fire
          call interpolator(ombb_aerosol_interp, ombb_time, omemisbb(:,:,lf), &
                       trim(ombb_emission_name(lf)), is, js)
        enddo
      case ('AEROCOM')
        nlevel_fire = nlevel_fire_AEROCOM
        alt_fire_min(1:nlevel_fire_AEROCOM) = &
              alt_fire_min_AEROCOM(1:nlevel_fire_AEROCOM)
        alt_fire_max(1:nlevel_fire_AEROCOM) = &
              alt_fire_max_AEROCOM(1:nlevel_fire_AEROCOM)
! Wildfire emissions at 6 levels from 0 to 6 km
! (cf. AEROCOM web site or Dentener et al., ACPD, 2006)
        do lf=1, nlevel_fire
          call interpolator(ombb_aerosol_interp, ombb_time, omemisbb(:,:,lf), &
                        trim(ombb_emission_name(lf)), is, js)
        enddo
    end select
   endif
   if ( trim(ombf_source).ne. ' ') then
     call interpolator(ombf_aerosol_interp, ombf_time, omemisbf, &
                       trim(ombf_emission_name(1)), is, js)
   endif
   if ( trim(omsh_source).ne. ' ') then
     call interpolator(omsh_aerosol_interp, omsh_time, omemissh, &
                  trim(omsh_emission_name(1)), is, js)
   endif
   if ( trim(omna_source).ne. ' ') then
     call interpolator(omna_aerosol_interp, omna_time, omemisbg, &
                       trim(omna_emission_name(1)), is, js)
   endif
   if ( trim(omss_source).ne. ' ') then
     call interpolator(omss_aerosol_interp, omss_time, dmso, &
                       trim(omss_emission_name(1)), is, js)
     do j = 1, jd
     do i = 1, id
       SST = t_surf(i,j)-273.15     ! Sea surface temperature [Celsius]
       if (ocn_flx_fraction(i,j).gt.0.) then
!  < Schmidt number (Saltzman et al., 1993) >
         Schm = 2674.0 - 147.12*SST + 3.726*(SST**2) - 0.038*(SST**3)
         Schm = max(1., Schm)
! ---  Liss and Merlivat (1986) -----------
         SchmCO2 = 600.
         if (w10m(i,j) .le. 3.6) then
           AKw = 0.17 * w10m(i,j)
         else if (w10m(i,j) .le. 13.) then
           AKw = 2.85 * w10m(i,j) - 9.65
              else
           AKw = 5.90 * w10m(i,j) - 49.3
         end if
         if (w10m(i,j) .le. 3.6) then
           AKw = AKw * ((SchmCO2/Schm) ** 0.667)
         else
           AKw = AKw * sqrt(SchmCO2/Schm)
         end if
         omemisocean(i,j) = coef_omss_emis*AKw/100./3600. * 1.e-6*ocn_flx_fraction(i,j)
       end if

     enddo
     enddo


   endif
!
    do j = 1, jd
      do i = 1, id

! --- For fosil fuel emissions, calculate the fraction of emission for
! --- each vertical levels
        fa1(:) = 0.
        fa2(:) = 0.
        do l = kd,2,-1
          Z1 = z_half(i,j,l+1)-z_half(i,j,kd+1)
          Z2 = z_half(i,j,l)-z_half(i,j,kd+1)
          if (Z2.ge.0.and.Z1.lt.ze1) then
            if (Z1.gt.0) then
              if (Z2.lt.ze1) then
                fa1(l)=(Z2-Z1)/ze1
              else
                fa1(l)=(ze1-Z1)/ze1
              endif
            else
              if (Z2.le.ze1) then
                fa1(l)=Z2/ze1
              else
                fa1(l)=1.
              endif
            endif
          endif
          
          if (Z2.ge.ze1.and.z1.lt.ze2) then
            if (Z1.gt.Ze1) then
              if (Z2.lt.ze2) then
                fa2(l)=(z2-z1)/(ze2-ze1)
              else
                fa2(l)=(ze2-z1)/(ze2-ze1)
              endif
            else
              if (Z2.le.ze2) then
                fa2(l)=(z2-ze1)/(ze2-ze1)
              else
                fa2(l)=1.
              endif
            endif
          endif
          if (Z1.gt.Ze2) exit
        enddo
!
! Calculate fraction of emission at every levels for open fires
!
        fbb(:,:)=0.
!
! In case of multiple levels, which are fixed
!
        if (nlevel_fire .gt. 1) then
          do l = kd,2,-1
            Z1 = z_half(i,j,l+1)-z_half(i,j,kd+1)
            Z2 = z_half(i,j,l)-z_half(i,j,kd+1)
            do lf=1,nlevel_fire
              del=alt_fire_max(lf)-alt_fire_min(lf)
              if (del.gt.0. .and. &
                  Z1.lt.alt_fire_max(lf).and.Z2.gt.alt_fire_min(lf) ) then
                if (Z1.ge.alt_fire_min(lf)) then
                  if (Z2 .lt. alt_fire_max(lf)) then
                    fbb(l,lf)=(Z2-Z1)/del
                  else
                    fbb(l,lf)=(alt_fire_max(lf)-z1)/del
                  endif
                else
                  if (Z2.le.alt_fire_max(lf)) then
                    fbb(l,lf) = (Z2-alt_fire_min(lf))/del
                  else
                    fbb(l,lf)=1.
                  endif
                endif
              endif
            enddo
          enddo
        else
!
! --- Inject equally through the boundary layer -------
!
          bltop = z_pbl(i,j)
          do l = kd,1,-1
            z1=z_half(i,j,l+1)-z_half(i,j,kd+1)
            z2=z_half(i,j,l)-z_half(i,j,kd+1)
            if (bltop.lt.z1) exit
            if (bltop.ge.z2) fbb(l,1)=(z2-z1)/bltop
            if (bltop.gt.z1.and.bltop.lt.z2) fbb(l,1) = (bltop-z1)/bltop
          enddo
        endif
! Fossil fuel emission
        bcemisff(i,j,:)=fa1(:) * bcemisff_l1(i,j) + fa2(:) * bcemisff_l2(i,j)
        omemisff(i,j,:)=fa1(:) * omemisff_l1(i,j) + fa2(:) * omemisff_l2(i,j)
! Open biomass burning fires emission
        do lf =1, nlevel_fire
          bcemisob(i,j,:)= bcemisob(i,j,:) + fbb(:,lf)*bcemisbb(i,j,lf)
          omemisob(i,j,:)= omemisob(i,j,:) + fbb(:,lf)*omemisbb(i,j,lf)
        enddo
! Biogenic
! Bio-fuel (if not included in fossil fuel inevntory)
! International shipping
        omphob_emis(i,j,kd) =  omemisbf(i,j) + omemissh(i,j) + &
           omemisbg(i,j) + omemisocean(i,j)
        omphob_emis(i,j,:)= omphob_emis(i,j,:) + omemisff(i,j,:) + &
           omemisob(i,j,:)

        do l=1,kd
          emisob(i,j) = emisob(i,j) + bcemisob(i,j,l) + omemisob(i,j,l)
          omemisob_2d(i,j) = omemisob_2d(i,j) + omemisob(i,j,l)
        end do

        omphil_emis(i,j,:) = omphob_emis(i,j,:) * frac_om_philic/pwt(i,j,:)
        omphob_emis(i,j,:) = omphob_emis(i,j,:) * frac_om_phobic/pwt(i,j,:)
!
      enddo
    enddo

!------------------------------------------------------------------------
! if frac_bcbb_phobic .eq. frac_bc_phobic, then frac_bcbb_philic must equal
!  frac_bc_philic, since sum must be 1.
!------------------------------------------------------------------------
    if (frac_bcbb_phobic == frac_bc_phobic) then
      do j = 1, jd
        do i = 1, id
          bcphob_emis(i,j,kd) =  bcemisbf(i,j) + bcemissh(i,j)
          bcphob_emis(i,j,:)= bcphob_emis(i,j,:) + &
                              bc_aircraft_EI * bcemisav(i,j,:)    + &
                              bcemisff(i,j,:) + bcemisob(i,j,:)
          bcphil_emis(i,j,:) = bcphob_emis(i,j,:)*frac_bc_philic/pwt(i,j,:)
          bcphob_emis(i,j,:) = bcphob_emis(i,j,:)*frac_bc_phobic/pwt(i,j,:)
        enddo
      enddo
    else
      do j = 1, jd
        do i = 1, id
          bcphob_emis(i,j,kd) =  bcemisbf(i,j) + bcemissh(i,j)
          bcphob_emis(i,j,:)= bcphob_emis(i,j,:) + &
                              bc_aircraft_EI * bcemisav(i,j,:)    + &
                              bcemisff(i,j,:) 
          bcphil_emis(i,j,:) = (bcphob_emis(i,j,:)*frac_bc_philic +   &
                               bcemisob(i,j,:)*frac_bcbb_philic)/pwt(i,j,:)
          bcphob_emis(i,j,:) = (bcphob_emis(i,j,:)*frac_bc_phobic +   &
                               bcemisob(i,j,:)*frac_bcbb_phobic)/pwt(i,j,:)
        enddo
      enddo
    endif

!------- compute black carbon phobic sink --------------
!
!  BCphob has a half-life time of 1.0days 
!   (corresponds to an e-folding time of 1.44 days)
!
!  sink = 1./(86400.*1.44) = 8.023e-6

    do l = 1,kd
      zzz1(:,:,l)=z_half(:,:,l)-z_half(:,:,l+1)
    enddo
    oh_conc1 = oh_conc * pwt/zzz1*2.079e19

    bcphob_sink(:,:,:) = 0.0
    if (do_dynamic_bc) then
      where (bcphob > 0.0)
        bcphob_sink = (oh_conc1/2.16e11/bcage + 1./86400./bcageslow)*bcphob
      elsewhere
        bcphob_sink = 0.0
      endwhere
      bc_tau = 1./(1./bcageslow + oh_conc1/(2.5e6*bcage))
    else
      where (bcphob > 0.0)
        bcphob_sink = 8.038e-6*bcphob
      elsewhere
        bcphob_sink = 0.0
      endwhere
      bc_tau = 1.44
    endif
!

!------- tendency ------------------

    bcphob_dt = bcphob_emis - bcphob_sink
    bcphil_dt = bcphil_emis + bcphob_sink
!
!------- compute organic carbon sink --------------
!
!  OCphob has a half-life time of 2.0days 
!   (corresponds to an e-folding time of 2.88 days)
!
!  sink = 1./(86400.*1.44) = 8.023e-6
!
    omphob_sink(:,:,:) = 0.0
    if (do_dynamic_om) then
      where (omphob > 0.0)
        omphob_sink = (oh_conc1/2.16e11/omage + 1./86400./omageslow)*omphob
      elsewhere
        omphob_sink = 0.0
      endwhere
    else
      where (omphob >= 0.0)
        omphob_sink = 8.023e-6*omphob
      elsewhere
        omphob_sink = 0.0
      endwhere
    endif

!------- tendency ------------------

      omphob_dt = omphob_emis - omphob_sink
      omphil_dt = omphil_emis + omphob_sink

!-----------------------------------------------------------------


!
! Send registered results to diag manager
!
      if (id_bc_emis_col > 0) then
! column emissions for bc and om

        bc_emis = 0.
        do k=1,kd
          bc_emis(:,:) = bc_emis(:,:) + bcphob_emis(:,:,k) +  &
                                                       bcphil_emis(:,:,k)
        end do

        used = send_data ( id_bc_emis_col, bc_emis, diag_time, &
              is_in=is,js_in=js)
      endif
 
      if (id_om_emis_col > 0) then
! column emissions for bc and om

        om_emis = 0.
        do k=1,kd
          om_emis(:,:) = om_emis(:,:) + omphob_emis(:,:,k) +    &
                                                       omphil_emis(:,:,k)
        end do

        used = send_data ( id_om_emis_col, om_emis, diag_time, &
              is_in=is,js_in=js)
      endif

!
! Send registered results to diag manager
!
      if (id_bc_emis_colv2 > 0) then
! column emissions for bc (corrected cmip units)

        bc_emis = 0.
        do k=1,kd
          bc_emis(:,:) = bc_emis(:,:) + pwt(:,:,k)*&
                               (bcphob_emis(:,:,k) + bcphil_emis(:,:,k))
        end do
        used = send_data ( id_bc_emis_colv2, bc_emis, diag_time, &
              is_in=is,js_in=js)
      endif
 
      if (id_om_emis_colv2 > 0) then
! column emissions for om (corrected cmip units)

        om_emis = 0.
        do k=1,kd
          om_emis(:,:) = om_emis(:,:) + pwt(:,:,k)* &
                               (omphob_emis(:,:,k) + omphil_emis(:,:,k))
        end do
        used = send_data ( id_om_emis_colv2, om_emis, diag_time, &
              is_in=is,js_in=js)
      endif

!-----------------------------------------------------------------

      if (id_bcphob_emis > 0) then
        used = send_data ( id_bcphob_emis, bcphob_emis, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
      if (id_bcphil_emis > 0) then
        used = send_data ( id_bcphil_emis, bcphil_emis, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
      if (id_omphob_emis > 0) then
        used = send_data ( id_omphob_emis, omphob_emis, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
      if (id_omphil_emis > 0) then
        used = send_data ( id_omphil_emis, omphil_emis, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
      if (id_bcphob_sink > 0) then
        used = send_data ( id_bcphob_sink, bcphob_sink, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
      if (id_omphob_sink > 0) then
        used = send_data ( id_omphob_sink, omphob_sink, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
      if (id_bcemisbf > 0) then
        used = send_data ( id_bcemisbf, bcemisbf, diag_time, &
              is_in=is,js_in=js)
      endif
      if (id_emisbb > 0) then
        used = send_data ( id_emisbb, emisob, diag_time, &
              is_in=is,js_in=js)
      endif
      if (id_omemisbb_col > 0) then
        used = send_data ( id_omemisbb_col, omemisob_2d, diag_time, &
              is_in=is,js_in=js)
      endif
      if (id_bcemisbb > 0) then
        used = send_data ( id_bcemisbb, bcemisob, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
      if (id_bcemissh > 0) then
        used = send_data ( id_bcemissh, bcemissh, diag_time, &
              is_in=is,js_in=js)
      endif
      if (id_bcemisff > 0) then
        used = send_data ( id_bcemisff, bcemisff, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
      if (id_bcemisav > 0) then
        used = send_data ( id_bcemisav, bcemisav*bc_aircraft_EI, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
      if (id_omemisbf > 0) then
        used = send_data ( id_omemisbf, omemisbf, diag_time, &
              is_in=is,js_in=js)
      endif
      if (id_omemisbb > 0) then
        used = send_data ( id_omemisbb, omemisob, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
      if (id_omemissh > 0) then
        used = send_data ( id_omemissh, omemissh, diag_time, &
              is_in=is,js_in=js)
      endif
      if (id_omemisff > 0) then
        used = send_data ( id_omemisff, omemisff, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
      if (id_omemisbg > 0) then
        used = send_data ( id_omemisbg, omemisbg, diag_time, &
              is_in=is,js_in=js)
      endif
      if (id_omemisoc > 0) then
        used = send_data ( id_omemisoc, omemisocean, diag_time, &
              is_in=is,js_in=js)
      endif
      if (id_bc_tau > 0) then
        used = send_data ( id_bc_tau, bc_tau, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
!
 end subroutine atmos_carbon_aerosol_driver

!#######################################################################
!<SUBROUTINE NAME ="atmos_carbon_aerosol_init">

!<OVERVIEW>
! Subroutine to initialize the carbon aerosol module.
!</OVERVIEW>
!<DESCRIPTION>
! This subroutine querys the tracer manager to find the indices for the 
! various carbonaceous aerosol tracers. It also registers the emission 
! fields for diagnostic purposes.
!  
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_carbon_aerosol_init (lonb, latb, r, axes, Time, mask)
!</TEMPLATE>
!   <IN NAME="lonb" TYPE="real" DIM="(:,:)">
!     The longitudes for the local domain.
!   </IN>
!   <IN NAME="latb" TYPE="real" DIM="(:,:)">
!     The latitudes for the local domain.
!   </IN>
!   <INOUT NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer fields dimensioned as (nlon,nlat,nlev,ntrace). 
!   </INOUT>
!   <IN NAME="mask" TYPE="real, optional" DIM="(:,:,:)">
!      optional mask (0. or 1.) that designates which grid points
!           are above (=1.) or below (=0.) the ground dimensioned as
!           (nlon,nlat,nlev).
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array dimensioned as
!      (nlon, nlat, nlev, ntime)
!   </IN>

 subroutine atmos_carbon_aerosol_init (lonb, latb, axes, Time, mask)

!-----------------------------------------------------------------------
real, dimension(:,:),    intent(in) :: lonb, latb
integer        , intent(in)                        :: axes(4)
type(time_type), intent(in)                        :: Time
real,            intent(in),    dimension(:,:,:), optional :: mask
character(len=7), parameter :: mod_name = 'tracers'
integer :: n
integer ::  unit, ierr, io, logunit


   if (module_is_initialized) return
!----------------------------------
!namelist files
      if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=carbon_aerosol_nml, iostat=io)
        ierr = check_nml_error(io,'carbon_aerosol_nml')
#else
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=carbon_aerosol_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'carbon_aerosol_nml')
        end do
10      call close_file (unit)
#endif
      endif
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number(version, tagname)
      logunit=stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                          write (logunit, nml=carbon_aerosol_nml)

!--------------------------------------------------------
!------namelist

!----- set initial value of carbon ------------

   n = get_tracer_index(MODEL_ATMOS,'bcphob')
   if (n>0) then
      nbcphobic = n
      call set_tracer_atts(MODEL_ATMOS,'bcphob','hphobic_bc','mmr')
      if (nbcphobic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (*,30) 'Hydrophobic BC',nbcphobic
      if (nbcphobic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (logunit,30) 'Hydrophobic BC',nbcphobic
   endif

   n = get_tracer_index(MODEL_ATMOS,'bcphil')
   if (n>0) then
      nbcphilic=n
      call set_tracer_atts(MODEL_ATMOS,'bcphil','hphilic_bc','mmr')
      if (nbcphilic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (*,30) 'Hydrophilic BC',nbcphilic
      if (nbcphilic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (logunit,30) 'Hydrophilic BC',nbcphilic
   endif

   n = get_tracer_index(MODEL_ATMOS,'omphob')
   if (n>0) then
      nomphobic=n
      call set_tracer_atts(MODEL_ATMOS,'omphob','phobic_om','mmr')
      if (nomphobic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (*,30) 'Hydrophobic OC',nomphobic
      if (nomphobic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (logunit,30) 'Hydrophobic OC',nomphobic
   endif

   n = get_tracer_index(MODEL_ATMOS,'omphil')
   if (n>0) then
      nomphilic=n
      call set_tracer_atts(MODEL_ATMOS,'omphil','philic_om','mmr')
      if (nomphilic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (*,30) 'Hydrophilic OC',nomphilic
      if (nomphilic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (logunit,30) 'Hydrophilic OC',nomphilic
   endif

30        format (A,' was initialized as tracer number ',i2)
!
!   Register Emissions as static fields (monthly)
!
     id_bcphob_emis = register_diag_field ( mod_name,           &
                    'bcphob_emis', axes(1:3),Time,          &
                    'BC phobic emission rate', 'kg/m2/sec' )

     id_bcphil_emis = register_diag_field ( mod_name,           &
                    'bcphil_emis', axes(1:3),Time,          &
                    'BC phylic emission rate', 'kg/m2/sec' )

     id_omphob_emis = register_diag_field ( mod_name,           &
                    'omphob_emis', axes(1:3),Time,          &
                    'OM phobic emission rate', 'kg/m2/sec' )

     id_omphil_emis = register_diag_field ( mod_name,           &
                    'omphil_emis', axes(1:3),Time,          &
                    'OM phylic emission rate', 'kg/m2/sec' )

     id_bc_emis_col = register_diag_field ( mod_name,           &
                    'bc_emis_col', axes(1:2),Time,          &
                    'total BC column emission rate', 'kg/m2/sec' )

     id_om_emis_col = register_diag_field ( mod_name,           &
                    'om_emis_col', axes(1:2),Time,          &
                    'total OM column emission rate', 'kg/m2/sec' )

     id_bc_emis_colv2 = register_diag_field ( mod_name,           &
                    'bc_emis_colv2', axes(1:2),Time,          &
                    'total BC column emission rate', 'kg/m2/sec' )

     id_om_emis_colv2 = register_diag_field ( mod_name,           &
                    'om_emis_colv2', axes(1:2),Time,          &
                    'total OM column emission rate', 'kg/m2/sec' )

     id_bcphob_sink = register_diag_field ( mod_name,           &
                    'bcphob_sink', axes(1:3),Time,          &
                    'BC phobic sink rate', 'kg/m2/sec' )

     id_omphob_sink = register_diag_field ( mod_name,           &
                    'omphob_sink', axes(1:3),Time,          &
                    'OM phobic sink rate', 'kg/m2/sec' )

     id_bcemisbf    = register_diag_field ( mod_name,           &
                    'bcemisbf', axes(1:2),Time,                 &
                    'BC biofuel emission', 'kg/m2/sec' )

     id_emisbb    = register_diag_field ( mod_name,           &
                    'emisbb', axes(1:2),Time,                 &
                    'column BC + OM open biomass burning emission', 'kg/m2/sec' )

     id_omemisbb_col    = register_diag_field ( mod_name,           &
                    'omemisbb_col', axes(1:2),Time,                 &
                    'column OM open biomass burning emission', 'kg/m2/sec' )

     id_bcemisbb    = register_diag_field ( mod_name,           &
                    'bcemisbb', axes(1:3),Time,                 &
                    'BC open biomass burning emission', 'kg/m2/sec' )

     id_bcemissh    = register_diag_field ( mod_name,           &
                    'bcemissh', axes(1:2),Time,                 &
                    'BC shipping emission', 'kg/m2/sec' )

     id_bcemisff    = register_diag_field ( mod_name,           &
                    'bcemisff', axes(1:3),Time,                 &
                    'BC fossil fuel emission', 'kg/m2/sec' )

     id_bcemisav    = register_diag_field ( mod_name,           &
                    'bcemisav', axes(1:3),Time,                 &
                    'BC aircraft emission', 'kg/m2/sec' )

     id_omemisbf    = register_diag_field ( mod_name,           &
                    'omemisbf', axes(1:2),Time,                 &
                    'OM biofuel emission', 'kg/m2/sec' )

     id_omemisbb    = register_diag_field ( mod_name,           &
                    'omemisbb', axes(1:3),Time,                 &
                    'OM open biomass burning emission', 'kg/m2/sec' )

     id_omemissh    = register_diag_field ( mod_name,           &
                    'omemissh', axes(1:2),Time,                 &
                    'OM shipping emission', 'kg/m2/sec' )

     id_omemisff    = register_diag_field ( mod_name,           &
                    'omemisff', axes(1:3),Time,                 &
                    'OM fossil fuel emission', 'kg/m2/sec' )

     id_omemisbg    = register_diag_field ( mod_name,           &
                    'omemisbg', axes(1:2),Time,                 &
                    'OM biogenic emission over land', 'kg/m2/sec' )

     id_omemisoc    = register_diag_field ( mod_name,           &
                    'omemisoc', axes(1:2),Time,                 &
                    'OM biogenic emission over ocean', 'kg/m2/sec' )

     id_bc_tau    = register_diag_field ( mod_name,           &
                    'bc_tau', axes(1:3),Time,                 &
                    'bcphob aging lifetime', 'days' )

!----------------------------------------------------------------------
!    initialize namelist entries
!----------------------------------------------------------------------
        bcff_offset = set_time (0,0)
        bcbb_offset = set_time (0,0)
        bcbf_offset = set_time (0,0)
        bcsh_offset = set_time (0,0)
        bcav_offset = set_time (0,0)
        omff_offset = set_time (0,0)
        ombb_offset = set_time (0,0)
        ombf_offset = set_time (0,0)
        omsh_offset = set_time (0,0)
        omna_offset = set_time (0,0)
        omss_offset = set_time (0,0)

        bcff_entry = set_time (0,0)
        bcbb_entry = set_time (0,0)
        bcbf_entry = set_time (0,0)
        bcsh_entry = set_time (0,0)
        bcav_entry = set_time (0,0)
        omff_entry = set_time (0,0)
        ombb_entry = set_time (0,0)
        ombf_entry = set_time (0,0)
        omsh_entry = set_time (0,0)
        omna_entry = set_time (0,0)
        omss_entry = set_time (0,0)

        bcff_negative_offset = .false.
        bcbb_negative_offset = .false.
        bcbf_negative_offset = .false.
        bcsh_negative_offset = .false.
        bcav_negative_offset = .false.
        omff_negative_offset = .false.
        ombb_negative_offset = .false.
        ombf_negative_offset = .false.
        omsh_negative_offset = .false.
        omna_negative_offset = .false.
        omss_negative_offset = .false.

        bcff_time_serie_type = 1
        bcbb_time_serie_type = 1
        bcbf_time_serie_type = 1
        bcsh_time_serie_type = 1
        bcav_time_serie_type = 1
        omff_time_serie_type = 1
        ombb_time_serie_type = 1
        ombf_time_serie_type = 1
        omsh_time_serie_type = 1
        omna_time_serie_type = 1
        omss_time_serie_type = 1
!----------------------------------------------------------------------
!    define the model base time  (defined in diag_table)
!----------------------------------------------------------------------
        model_init_time = get_base_time()
   if ( trim(bcff_source) .ne. ' ') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(bcff_time_dependency_type) == 'constant' ) then
        bcff_time_serie_type = 1
        bcff_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcff are constant in atmos_carbon_aerosol module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for bcff is selected.
!---------------------------------------------------------------------
      else if (trim(bcff_time_dependency_type) == 'time_varying') then
        bcff_time_serie_type = 3
        if (bcff_dataset_entry(1) == 1 .and. &
            bcff_dataset_entry(2) == 1 .and. &
            bcff_dataset_entry(3) == 1 .and. &
            bcff_dataset_entry(4) == 0 .and. &
            bcff_dataset_entry(5) == 0 .and. &
            bcff_dataset_entry(6) == 0 ) then
          bcff_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcff_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          bcff_entry  = set_date (bcff_dataset_entry(1), &
                                  bcff_dataset_entry(2), &
                                  bcff_dataset_entry(3), &
                                  bcff_dataset_entry(4), &
                                  bcff_dataset_entry(5), &
                                  bcff_dataset_entry(6))
        endif
        call print_date (bcff_entry , str= &
          'Data from bcff timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        bcff_offset = bcff_entry - model_init_time
        if (model_init_time > bcff_entry) then
          bcff_negative_offset = .true.
        else
          bcff_negative_offset = .false.
        endif
      else if (trim(bcff_time_dependency_type) == 'fixed_year') then
        bcff_time_serie_type = 2
        if (bcff_dataset_entry(1) == 1 .and. &
            bcff_dataset_entry(2) == 1 .and. &
            bcff_dataset_entry(3) == 1 .and. &
            bcff_dataset_entry(4) == 0 .and. &
            bcff_dataset_entry(5) == 0 .and. &
            bcff_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set bcff_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcff_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        bcff_entry  = set_date (bcff_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'bcff is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcff correspond to year :', &
                    bcff_dataset_entry(1)
        endif
     endif
     select case (trim(bcff_source))
       case ('cooke_and_wilson_1996')
         if (trim(bcff_input_name(1)) .eq. ' ') then
           bcff_emission_name(1)='bcff_cw96'
         else
           bcff_emission_name(1)=trim(bcff_input_name(1))
         endif
      case ('cooke_1999')
         if (trim(bcff_input_name(1)) .eq. ' ') then
           bcff_emission_name(1)='bcff_cooke99'
         else
           bcff_emission_name(1)=trim(bcff_input_name(1))
         endif
      case ('bond_2004') 
         if (trim(bcff_input_name(1)) .eq. ' ') then
           bcff_emission_name(1)='bcff_bond'
         else
           bcff_emission_name(1)=trim(bcff_input_name(1))
         endif
      case ('gocart_2007') 
         if (trim(bcff_input_name(1)) .eq. ' ') then
           bcff_emission_name(1)='bc_anthro'
         else
           bcff_emission_name(1)=trim(bcff_input_name(1))
         endif
      case ('AEROCOM') 
         if (trim(bcff_input_name(1)) .eq. ' ') then
           bcff_emission_name(1)='BC1ff'
         else
           bcff_emission_name(1)=trim(bcff_input_name(1))
         endif
     end select
     call interpolator_init (bcff_aerosol_interp,             &
                             trim(bcff_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = bcff_emission_name,        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )
   endif
   if ( trim(bcbb_source) .ne. ' ') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(bcbb_time_dependency_type) == 'constant' ) then
        bcbb_time_serie_type = 1
        bcbb_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcbb are constant in atmos_carbon_aerosol module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for bcbb is selected.
!---------------------------------------------------------------------
      else if (trim(bcbb_time_dependency_type) == 'time_varying') then
        bcbb_time_serie_type = 3
        if (bcbb_dataset_entry(1) == 1 .and. &
            bcbb_dataset_entry(2) == 1 .and. &
            bcbb_dataset_entry(3) == 1 .and. &
            bcbb_dataset_entry(4) == 0 .and. &
            bcbb_dataset_entry(5) == 0 .and. &
            bcbb_dataset_entry(6) == 0 ) then
          bcbb_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcbb_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          bcbb_entry  = set_date (bcbb_dataset_entry(1), &
                                  bcbb_dataset_entry(2), &
                                  bcbb_dataset_entry(3), &
                                  bcbb_dataset_entry(4), &
                                  bcbb_dataset_entry(5), &
                                  bcbb_dataset_entry(6))
        endif
        call print_date (bcbb_entry , str= &
          'Data from bcbb timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        bcbb_offset = bcbb_entry - model_init_time
        if (model_init_time > bcbb_entry) then
          bcbb_negative_offset = .true.
        else
          bcbb_negative_offset = .false.
        endif
      else if (trim(bcbb_time_dependency_type) == 'fixed_year') then
        bcbb_time_serie_type = 2
        if (bcbb_dataset_entry(1) == 1 .and. &
            bcbb_dataset_entry(2) == 1 .and. &
            bcbb_dataset_entry(3) == 1 .and. &
            bcbb_dataset_entry(4) == 0 .and. &
            bcbb_dataset_entry(5) == 0 .and. &
            bcbb_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set bcbb_dataset_entry when using fixed_year source', FATAL)
        endif
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcbb_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        bcbb_entry  = set_date (bcbb_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'bcbb is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcbb correspond to year :', &
                    bcbb_dataset_entry(1)
        endif
     endif
     select case (trim(bcbb_source))
       case ('cooke_and_wilson_1996')
         if (trim(bcbb_input_name(1)) .eq. ' ') then
           bcbb_emission_name(1)='bcbb_cw96'
         else
           bcbb_emission_name(1)=trim(bcbb_input_name(1))
         endif
         call interpolator_init (bcbb_aerosol_interp,           &
           trim(bcbb_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
           data_names=bcbb_emission_name(1:1), vert_interp=(/INTERP_WEIGHTED_P/))
      case ('bond_2004') 
         if (trim(bcbb_input_name(1)) .eq. ' ') then
           bcbb_emission_name(1)='bcob_bond'
         else
           bcbb_emission_name(1)=trim(bcbb_input_name(1))
         endif
         call interpolator_init (bcbb_aerosol_interp,           &
           trim(bcbb_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
           data_names=bcbb_emission_name(1:1), vert_interp=(/INTERP_WEIGHTED_P/))
      case ('gocart_2007') 
         if (trim(bcbb_input_name(1)) .eq. ' ') then
           bcbb_emission_name(1)='bc_biobur'
         else
           bcbb_emission_name(1)=trim(bcbb_input_name(1))
         endif
         call interpolator_init (bcbb_aerosol_interp,           &
           trim(bcbb_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
           data_names=bcbb_emission_name(1:1), vert_interp=(/INTERP_WEIGHTED_P/))
       case ('GEIA')
         if (trim(bcbb_input_name(1)) .eq. ' ') then
           bcbb_emission_name(1)='bc_geia1'
           bcbb_emission_name(2)='bc_geia2'
         else
           bcbb_emission_name(1)=trim(bcbb_input_name(1))
           bcbb_emission_name(2)=trim(bcbb_input_name(2))
         endif
         call interpolator_init (bcbb_aerosol_interp,           &
           trim(bcbb_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
           data_names=bcbb_emission_name(1:2),vert_interp=(/INTERP_WEIGHTED_P/))
       case ('AEROCOM')
         if (trim(bcbb_input_name(1)) .eq. ' ') then
           bcbb_emission_name(1)='GFED_BC_l1'
           bcbb_emission_name(2)='GFED_BC_l2'
           bcbb_emission_name(3)='GFED_BC_l3'
           bcbb_emission_name(4)='GFED_BC_l4'
           bcbb_emission_name(5)='GFED_BC_l5'
           bcbb_emission_name(6)='GFED_BC_l6'
         else
           bcbb_emission_name(1)(:)=trim(bcbb_input_name(1)(:))
           bcbb_emission_name(2)(:)=trim(bcbb_input_name(2)(:))
           bcbb_emission_name(3)(:)=trim(bcbb_input_name(3)(:))
           bcbb_emission_name(4)(:)=trim(bcbb_input_name(4)(:))
           bcbb_emission_name(5)(:)=trim(bcbb_input_name(5)(:))
           bcbb_emission_name(6)(:)=trim(bcbb_input_name(6)(:))
         endif
         call interpolator_init (bcbb_aerosol_interp,           &
           trim(bcbb_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
           data_names=bcbb_emission_name(1:6),vert_interp=(/INTERP_WEIGHTED_P/))
     end select
   endif
   if ( trim(bcsh_source) .ne. ' ') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(bcsh_time_dependency_type) == 'constant' ) then
        bcsh_time_serie_type = 1
        bcsh_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcsh are constant in atmos_carbon_aerosol module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for bcsh is selected.
!---------------------------------------------------------------------
      else if (trim(bcsh_time_dependency_type) == 'time_varying') then
        bcsh_time_serie_type = 3
        if (bcsh_dataset_entry(1) == 1 .and. &
            bcsh_dataset_entry(2) == 1 .and. &
            bcsh_dataset_entry(3) == 1 .and. &
            bcsh_dataset_entry(4) == 0 .and. &
            bcsh_dataset_entry(5) == 0 .and. &
            bcsh_dataset_entry(6) == 0 ) then
          bcsh_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcsh_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          bcsh_entry  = set_date (bcsh_dataset_entry(1), &
                                  bcsh_dataset_entry(2), &
                                  bcsh_dataset_entry(3), &
                                  bcsh_dataset_entry(4), &
                                  bcsh_dataset_entry(5), &
                                  bcsh_dataset_entry(6))
        endif
        call print_date (bcsh_entry , str= &
          'Data from bcsh timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        bcsh_offset = bcsh_entry - model_init_time
        if (model_init_time > bcsh_entry) then
          bcsh_negative_offset = .true.
        else
          bcsh_negative_offset = .false.
        endif
      else if (trim(bcsh_time_dependency_type) == 'fixed_year') then
        bcsh_time_serie_type = 2
        if (bcsh_dataset_entry(1) == 1 .and. &
            bcsh_dataset_entry(2) == 1 .and. &
            bcsh_dataset_entry(3) == 1 .and. &
            bcsh_dataset_entry(4) == 0 .and. &
            bcsh_dataset_entry(5) == 0 .and. &
            bcsh_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set bcsh_dataset_entry when using fixed_year source', FATAL)
        endif
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcsh_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        bcsh_entry  = set_date (bcsh_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'bcsh is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcsh correspond to year :', &
                    bcsh_dataset_entry(1)
        endif
     endif
     if (trim(bcsh_input_name(1)) .eq. ' ') then
       bcsh_emission_name(1)='bc_ship'
     else
       bcsh_emission_name(1)=trim(bcsh_input_name(1))
     endif

     call interpolator_init (bcsh_aerosol_interp,             &
                             trim(bcsh_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = bcsh_emission_name,        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )
   endif
   if ( trim(bcav_source) .ne. ' ') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(bcav_time_dependency_type) == 'constant' ) then
        bcav_time_serie_type = 1
        bcav_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcav are constant in atmos_carbon_aerosol module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for bcav is selected.
!---------------------------------------------------------------------
      else if (trim(bcav_time_dependency_type) == 'time_varying') then
        bcav_time_serie_type = 3
        if (bcav_dataset_entry(1) == 1 .and. &
            bcav_dataset_entry(2) == 1 .and. &
            bcav_dataset_entry(3) == 1 .and. &
            bcav_dataset_entry(4) == 0 .and. &
            bcav_dataset_entry(5) == 0 .and. &
            bcav_dataset_entry(6) == 0 ) then
          bcav_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcav_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          bcav_entry  = set_date (bcav_dataset_entry(1), &
                                  bcav_dataset_entry(2), &
                                  bcav_dataset_entry(3), &
                                  bcav_dataset_entry(4), &
                                  bcav_dataset_entry(5), &
                                  bcav_dataset_entry(6))
        endif
        call print_date (bcav_entry , str= &
          'Data from bcav timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        bcav_offset = bcav_entry - model_init_time
        if (model_init_time > bcav_entry) then
          bcav_negative_offset = .true.
        else
          bcav_negative_offset = .false.
        endif
      else if (trim(bcav_time_dependency_type) == 'fixed_year') then
        bcav_time_serie_type = 2
        if (bcav_dataset_entry(1) == 1 .and. &
            bcav_dataset_entry(2) == 1 .and. &
            bcav_dataset_entry(3) == 1 .and. &
            bcav_dataset_entry(4) == 0 .and. &
            bcav_dataset_entry(5) == 0 .and. &
            bcav_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set bcav_dataset_entry when using fixed_year source', FATAL)
        endif
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcav_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        bcav_entry  = set_date (bcav_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'bcav is defined from a single annual cycle - no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcav correspond to year :', bcav_dataset_entry(1)
        endif
     endif
     if (trim(bcav_input_name(1)) .eq. ' ') then
       bcav_emission_name(1)='bc_aircraft'
     else
       bcav_emission_name(1)=trim(bcav_input_name(1))
     endif
     call interpolator_init (bcav_aerosol_interp,             &
                             trim(bcav_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = bcav_emission_name,        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )
   endif
   if ( trim(bcbf_source) .ne. ' ') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(bcbf_time_dependency_type) == 'constant' ) then
        bcbf_time_serie_type = 1
        bcbf_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcbf are constant in atmos_carbon_aerosol module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for bcbf is selected.
!---------------------------------------------------------------------
      else if (trim(bcbf_time_dependency_type) == 'time_varying') then
        bcbf_time_serie_type = 3
        if (bcbf_dataset_entry(1) == 1 .and. &
            bcbf_dataset_entry(2) == 1 .and. &
            bcbf_dataset_entry(3) == 1 .and. &
            bcbf_dataset_entry(4) == 0 .and. &
            bcbf_dataset_entry(5) == 0 .and. &
            bcbf_dataset_entry(6) == 0 ) then
          bcbf_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcbf_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          bcbf_entry  = set_date (bcbf_dataset_entry(1), &
                                  bcbf_dataset_entry(2), &
                                  bcbf_dataset_entry(3), &
                                  bcbf_dataset_entry(4), &
                                  bcbf_dataset_entry(5), &
                                  bcbf_dataset_entry(6))
        endif
        call print_date (bcbf_entry , str= &
          'Data from bcbf timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        bcbf_offset = bcbf_entry - model_init_time
        if (model_init_time > bcbf_entry) then
          bcbf_negative_offset = .true.
        else
          bcbf_negative_offset = .false.
        endif
      else if (trim(bcbf_time_dependency_type) == 'fixed_year') then
        bcbf_time_serie_type = 2
        if (bcbf_dataset_entry(1) == 1 .and. &
            bcbf_dataset_entry(2) == 1 .and. &
            bcbf_dataset_entry(3) == 1 .and. &
            bcbf_dataset_entry(4) == 0 .and. &
            bcbf_dataset_entry(5) == 0 .and. &
            bcbf_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set bcbf_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcbf_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        bcbf_entry  = set_date (bcbf_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'bcbf is defined from a single annual cycle - no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcbf correspond to year :', &
                    bcbf_dataset_entry(1)
        endif
     endif
     if (trim(bcbf_input_name(1)) .eq. ' ') then
       bcbf_emission_name(1)='bcbf_bond'
     else
       bcbf_emission_name(1)=trim(bcbf_input_name(1))
     endif
     call interpolator_init (bcbf_aerosol_interp,             &
                             trim(bcbf_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = bcbf_emission_name,        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )
   endif
   if ( trim(omff_source) .ne. ' ') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(omff_time_dependency_type) == 'constant' ) then
        omff_time_serie_type = 1
        omff_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'omff are constant in atmos_carbon_aerosol module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for omff is selected.
!---------------------------------------------------------------------
      else if (trim(omff_time_dependency_type) == 'time_varying') then
        omff_time_serie_type = 3
        if (omff_dataset_entry(1) == 1 .and. &
            omff_dataset_entry(2) == 1 .and. &
            omff_dataset_entry(3) == 1 .and. &
            omff_dataset_entry(4) == 0 .and. &
            omff_dataset_entry(5) == 0 .and. &
            omff_dataset_entry(6) == 0 ) then
          omff_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to omff_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          omff_entry  = set_date (omff_dataset_entry(1), &
                                  omff_dataset_entry(2), &
                                  omff_dataset_entry(3), &
                                  omff_dataset_entry(4), &
                                  omff_dataset_entry(5), &
                                  omff_dataset_entry(6))
        endif
        call print_date (omff_entry , str= &
          'Data from omff timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        omff_offset = omff_entry - model_init_time
        if (model_init_time > omff_entry) then
          omff_negative_offset = .true.
        else
          omff_negative_offset = .false.
        endif
      else if (trim(omff_time_dependency_type) == 'fixed_year') then
        omff_time_serie_type = 2
        if (omff_dataset_entry(1) == 1 .and. &
            omff_dataset_entry(2) == 1 .and. &
            omff_dataset_entry(3) == 1 .and. &
            omff_dataset_entry(4) == 0 .and. &
            omff_dataset_entry(5) == 0 .and. &
            omff_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set omff_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to omff_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        omff_entry  = set_date (omff_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'omff is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'omff correspond to year :', &
                    omff_dataset_entry(1)
        endif
     endif
     select case (trim(omff_source))
      case ('cooke_1999')
         if (trim(omff_input_name(1)) .eq. ' ') then
           omff_emission_name(1)='omff_cooke99'
         else
           omff_emission_name(1)=trim(omff_input_name(1))
         endif
      case ('bond_2004') 
         if (trim(omff_input_name(1)) .eq. ' ') then
           omff_emission_name(1)='omff_bond'
         else
           omff_emission_name(1)=trim(omff_input_name(1))
         endif
      case ('gocart_2007') 
         if (trim(omff_input_name(1)) .eq. ' ') then
           omff_emission_name(1)='om_anthro'
         else
           omff_emission_name(1)=trim(omff_input_name(1))
         endif
      case ('AEROCOM') 
         if (trim(omff_input_name(1)) .eq. ' ') then
           omff_emission_name(1)='POMff'
         else
           omff_emission_name(1)=trim(omff_input_name(1))
         endif
     end select
     call interpolator_init (omff_aerosol_interp,           &
       trim(omff_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
       data_names=omff_emission_name(1:1),vert_interp=(/INTERP_WEIGHTED_P/))
   endif
   if ( trim(ombb_source) .ne. ' ') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(ombb_time_dependency_type) == 'constant' ) then
        ombb_time_serie_type = 1
        ombb_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'ombb are constant in atmos_carbon_aerosol module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for ombb is selected.
!---------------------------------------------------------------------
      else if (trim(ombb_time_dependency_type) == 'time_varying') then
        ombb_time_serie_type = 3
        if (ombb_dataset_entry(1) == 1 .and. &
            ombb_dataset_entry(2) == 1 .and. &
            ombb_dataset_entry(3) == 1 .and. &
            ombb_dataset_entry(4) == 0 .and. &
            ombb_dataset_entry(5) == 0 .and. &
            ombb_dataset_entry(6) == 0 ) then
          ombb_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to ombb_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          ombb_entry  = set_date (ombb_dataset_entry(1), &
                                  ombb_dataset_entry(2), &
                                  ombb_dataset_entry(3), &
                                  ombb_dataset_entry(4), &
                                  ombb_dataset_entry(5), &
                                  ombb_dataset_entry(6))
        endif
        call print_date (ombb_entry , str= &
          'Data from ombb timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        ombb_offset = ombb_entry - model_init_time
        if (model_init_time > ombb_entry) then
          ombb_negative_offset = .true.
        else
          ombb_negative_offset = .false.
        endif
      else if (trim(ombb_time_dependency_type) == 'fixed_year') then
        ombb_time_serie_type = 2
        if (ombb_dataset_entry(1) == 1 .and. &
            ombb_dataset_entry(2) == 1 .and. &
            ombb_dataset_entry(3) == 1 .and. &
            ombb_dataset_entry(4) == 0 .and. &
            ombb_dataset_entry(5) == 0 .and. &
            ombb_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set ombb_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to ombb_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        ombb_entry  = set_date (ombb_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'ombb is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'ombb correspond to year :', &
                    ombb_dataset_entry(1)
        endif
     endif
     select case (trim(ombb_source))
       case ('cooke_and_wilson_1996')
         if (trim(ombb_input_name(1)) .eq. ' ') then
           ombb_emission_name(1)='ombb_cw96'
         else
           ombb_emission_name(1)=trim(ombb_input_name(1))
         endif
         call interpolator_init (ombb_aerosol_interp,           &
           trim(ombb_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
           data_names=ombb_emission_name(1:1),vert_interp=(/INTERP_WEIGHTED_P/))
      case ('bond_2004') 
         if (trim(ombb_input_name(1)) .eq. ' ') then
           ombb_emission_name(1)='omob_bond'
         else
           ombb_emission_name(1)=trim(ombb_input_name(1))
         endif
         call interpolator_init (ombb_aerosol_interp,           &
           trim(ombb_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
           data_names=ombb_emission_name(1:1),vert_interp=(/INTERP_WEIGHTED_P/))
      case ('gocart_2007') 
         if (trim(ombb_input_name(1)) .eq. ' ') then
           ombb_emission_name(1)='om_biobur'
         else
           ombb_emission_name(1)=trim(ombb_input_name(1))
         endif
         call interpolator_init (ombb_aerosol_interp,           &
           trim(ombb_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
           data_names=ombb_emission_name(1:1),vert_interp=(/INTERP_WEIGHTED_P/))
       case ('GEIA')
         if (trim(ombb_input_name(1)) .eq. ' ') then
           ombb_emission_name(1)='om_geia1'
           ombb_emission_name(2)='om_geia2'
         else
           ombb_emission_name(1)=trim(ombb_input_name(1))
           ombb_emission_name(2)=trim(ombb_input_name(2))
         endif
         call interpolator_init (ombb_aerosol_interp,           &
           trim(ombb_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
           data_names=ombb_emission_name(1:2),vert_interp=(/INTERP_WEIGHTED_P/))
       case ('AEROCOM')
         if (trim(ombb_input_name(1)) .eq. ' ') then
           ombb_emission_name(1)='GFED_OM_l1'
           ombb_emission_name(2)='GFED_OM_l2'
           ombb_emission_name(3)='GFED_OM_l3'
           ombb_emission_name(4)='GFED_OM_l4'
           ombb_emission_name(5)='GFED_OM_l5'
           ombb_emission_name(6)='GFED_OM_l6'
         else
           ombb_emission_name(1)=trim(ombb_input_name(1))
           ombb_emission_name(2)=trim(ombb_input_name(2))
           ombb_emission_name(3)=trim(ombb_input_name(3))
           ombb_emission_name(4)=trim(ombb_input_name(4))
           ombb_emission_name(5)=trim(ombb_input_name(5))
           ombb_emission_name(6)=trim(ombb_input_name(6))
         endif
         call interpolator_init (ombb_aerosol_interp,           &
           trim(ombb_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
           data_names=ombb_emission_name(1:6),vert_interp=(/INTERP_WEIGHTED_P/))
     end select
   endif
   if ( trim(ombf_source) .ne. ' ') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(ombf_time_dependency_type) == 'constant' ) then
        ombf_time_serie_type = 1
        ombf_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'ombf are constant in atmos_carbon_aerosol module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for ombf is selected.
!---------------------------------------------------------------------
      else if (trim(ombf_time_dependency_type) == 'time_varying') then
        ombf_time_serie_type = 3
        if (ombf_dataset_entry(1) == 1 .and. &
            ombf_dataset_entry(2) == 1 .and. &
            ombf_dataset_entry(3) == 1 .and. &
            ombf_dataset_entry(4) == 0 .and. &
            ombf_dataset_entry(5) == 0 .and. &
            ombf_dataset_entry(6) == 0 ) then
          ombf_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to ombf_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          ombf_entry  = set_date (ombf_dataset_entry(1), &
                                  ombf_dataset_entry(2), &
                                  ombf_dataset_entry(3), &
                                  ombf_dataset_entry(4), &
                                  ombf_dataset_entry(5), &
                                  ombf_dataset_entry(6))
        endif
        call print_date (ombf_entry , str= &
          'Data from ombf timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        ombf_offset = ombf_entry - model_init_time
        if (model_init_time > ombf_entry) then
          ombf_negative_offset = .true.
        else
          ombf_negative_offset = .false.
        endif
      else if (trim(ombf_time_dependency_type) == 'fixed_year') then
        ombf_time_serie_type = 2
        if (ombf_dataset_entry(1) == 1 .and. &
            ombf_dataset_entry(2) == 1 .and. &
            ombf_dataset_entry(3) == 1 .and. &
            ombf_dataset_entry(4) == 0 .and. &
            ombf_dataset_entry(5) == 0 .and. &
            ombf_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set ombf_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to ombf_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        ombf_entry  = set_date (ombf_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'ombf is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'ombf correspond to year :', &
                    ombf_dataset_entry(1)
        endif
     endif
     if (trim(ombf_input_name(1)) .eq. ' ') then
       ombf_emission_name(1)='ombf_bond'
     else
       ombf_emission_name(1)=trim(ombf_input_name(1))
     endif
     call interpolator_init (ombf_aerosol_interp,           &
       trim(ombf_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
       data_names=ombf_emission_name(1:1),vert_interp=(/INTERP_WEIGHTED_P/))
   endif
   if ( trim(omsh_source) .ne. ' ') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(omsh_time_dependency_type) == 'constant' ) then
        omsh_time_serie_type = 1
        omsh_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'omsh are constant in atmos_carbon_aerosol module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for omsh is selected.
!---------------------------------------------------------------------
      else if (trim(omsh_time_dependency_type) == 'time_varying') then
        omsh_time_serie_type = 3
        if (omsh_dataset_entry(1) == 1 .and. &
            omsh_dataset_entry(2) == 1 .and. &
            omsh_dataset_entry(3) == 1 .and. &
            omsh_dataset_entry(4) == 0 .and. &
            omsh_dataset_entry(5) == 0 .and. &
            omsh_dataset_entry(6) == 0 ) then
          omsh_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to omsh_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          omsh_entry  = set_date (omsh_dataset_entry(1), &
                                  omsh_dataset_entry(2), &
                                  omsh_dataset_entry(3), &
                                  omsh_dataset_entry(4), &
                                  omsh_dataset_entry(5), &
                                  omsh_dataset_entry(6))
        endif
        call print_date (omsh_entry , str= &
          'Data from omsh timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        omsh_offset = omsh_entry - model_init_time
        if (model_init_time > omsh_entry) then
          omsh_negative_offset = .true.
        else
          omsh_negative_offset = .false.
        endif
      else if (trim(omsh_time_dependency_type) == 'fixed_year') then
        omsh_time_serie_type = 2
        if (omsh_dataset_entry(1) == 1 .and. &
            omsh_dataset_entry(2) == 1 .and. &
            omsh_dataset_entry(3) == 1 .and. &
            omsh_dataset_entry(4) == 0 .and. &
            omsh_dataset_entry(5) == 0 .and. &
            omsh_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set omsh_dataset_entry when using fixed_year source', FATAL)
        endif
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to omsh_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        omsh_entry  = set_date (omsh_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'omsh is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'omsh correspond to year :', &
                    omsh_dataset_entry(1)
        endif
     endif
     if (trim(omsh_input_name(1)) .eq. ' ') then
       omsh_emission_name(1)='om_ship'
     else
       omsh_emission_name(1)=trim(omsh_input_name(1))
     endif
     call interpolator_init (omsh_aerosol_interp,           &
       trim(omsh_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
       data_names=omsh_emission_name(1:1),vert_interp=(/INTERP_WEIGHTED_P/))
   endif
   if ( trim(omna_source) .ne. ' ') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(omna_time_dependency_type) == 'constant' ) then
        omna_time_serie_type = 1
        omna_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'omna are constant in atmos_carbon_aerosol module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for omna is selected.
!---------------------------------------------------------------------
      else if (trim(omna_time_dependency_type) == 'time_varying') then
        omna_time_serie_type = 3
        if (omna_dataset_entry(1) == 1 .and. &
            omna_dataset_entry(2) == 1 .and. &
            omna_dataset_entry(3) == 1 .and. &
            omna_dataset_entry(4) == 0 .and. &
            omna_dataset_entry(5) == 0 .and. &
            omna_dataset_entry(6) == 0 ) then
          omna_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to omna_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          omna_entry  = set_date (omna_dataset_entry(1), &
                                  omna_dataset_entry(2), &
                                  omna_dataset_entry(3), &
                                  omna_dataset_entry(4), &
                                  omna_dataset_entry(5), &
                                  omna_dataset_entry(6))
        endif
        call print_date (omna_entry , str= &
          'Data from omna timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        omna_offset = omna_entry - model_init_time
        if (model_init_time > omna_entry) then
          omna_negative_offset = .true.
        else
          omna_negative_offset = .false.
        endif
      else if (trim(omna_time_dependency_type) == 'fixed_year') then
        omna_time_serie_type = 2
        if (omna_dataset_entry(1) == 1 .and. &
            omna_dataset_entry(2) == 1 .and. &
            omna_dataset_entry(3) == 1 .and. &
            omna_dataset_entry(4) == 0 .and. &
            omna_dataset_entry(5) == 0 .and. &
            omna_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set omna_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to omna_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        omna_entry  = set_date (omna_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'omna is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'omna correspond to year :', &
                    omna_dataset_entry(1)
        endif
     endif
     if (trim(omna_input_name(1)) .eq. ' ') then
       omna_emission_name(1)='omemisnat'
     else
       omna_emission_name(1)=trim(omna_input_name(1))
     endif
     call interpolator_init (omna_aerosol_interp,           &
       trim(omna_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
       data_names=omna_emission_name(1:1),vert_interp=(/INTERP_WEIGHTED_P/))
   endif
   if ( trim(omss_source) .ne. ' ') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(omss_time_dependency_type) == 'constant' ) then
        omss_time_serie_type = 1
        omss_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'omss are constant in atmos_carbon_aerosol module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for omss is selected.
!---------------------------------------------------------------------
      else if (trim(omss_time_dependency_type) == 'time_varying') then
        omss_time_serie_type = 3
        if (omss_dataset_entry(1) == 1 .and. &
            omss_dataset_entry(2) == 1 .and. &
            omss_dataset_entry(3) == 1 .and. &
            omss_dataset_entry(4) == 0 .and. &
            omss_dataset_entry(5) == 0 .and. &
            omss_dataset_entry(6) == 0 ) then
          omss_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to omss_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          omss_entry  = set_date (omss_dataset_entry(1), &
                                  omss_dataset_entry(2), &
                                  omss_dataset_entry(3), &
                                  omss_dataset_entry(4), &
                                  omss_dataset_entry(5), &
                                  omss_dataset_entry(6))
        endif
        call print_date (omss_entry , str= &
          'Data from omss timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        omss_offset = omss_entry - model_init_time
        if (model_init_time > omss_entry) then
          omss_negative_offset = .true.
        else
          omss_negative_offset = .false.
        endif
      else if (trim(omss_time_dependency_type) == 'fixed_year') then
        omss_time_serie_type = 2
        if (omss_dataset_entry(1) == 1 .and. &
            omss_dataset_entry(2) == 1 .and. &
            omss_dataset_entry(3) == 1 .and. &
            omss_dataset_entry(4) == 0 .and. &
            omss_dataset_entry(5) == 0 .and. &
            omss_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set omss_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to omss_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        omss_entry  = set_date (omss_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'omss is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'omss correspond to year :', &
                    omss_dataset_entry(1)
        endif
     endif
     if (trim(omss_input_name(1)) .eq. ' ') then
       omss_emission_name(1)='DMSo'
     else
       omss_emission_name(1)=trim(omss_input_name(1))
     endif
     call interpolator_init (omss_aerosol_interp,           &
       trim(omss_filename), lonb, latb, data_out_of_bounds=(/CONSTANT/), &
       data_names=omss_emission_name(1:1),vert_interp=(/INTERP_WEIGHTED_P/))
     if (omss_coef .le. -990) then
       coef_omss_emis = 1.
     else
       coef_omss_emis = omss_coef
     endif
   endif


   call write_version_number(version, tagname)
   module_is_initialized = .TRUE.

!-----------------------------------------------------------------------

end subroutine atmos_carbon_aerosol_init



!######################################################################

subroutine atmos_carbon_aerosol_time_vary (model_time)


type(time_type), intent(in) :: model_time


    integer ::  yr, dum, mo_yr, mo, dy, hr, mn, sc, dayspmn


    if ( trim(bcff_source) .ne. ' ') then
!--------------------------------------------------------------------
!    define the time in the bcff data set from which data is to be 
!    taken. if bcff is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(bcff_time_serie_type .eq. 3) then
       if (bcff_negative_offset) then
         bcff_time = model_time - bcff_offset
       else
         bcff_time = model_time + bcff_offset
       endif
     else 
       if(bcff_time_serie_type .eq. 2 ) then
         call get_date (bcff_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(bcff_entry)
           if (dayspmn /= 29) then
             bcff_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             bcff_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           bcff_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         bcff_time = model_time
       endif
     endif
     call obtain_interpolator_time_slices   &
                                   (bcff_aerosol_interp, bcff_time)
   endif

   if ( trim(bcbb_source).ne.' ') then

!--------------------------------------------------------------------
!    define the time in the bcbb data set from which data is to be 
!    taken. if bcbb is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(bcbb_time_serie_type .eq. 3) then
       if (bcbb_negative_offset) then
         bcbb_time = model_time - bcbb_offset
       else
         bcbb_time = model_time + bcbb_offset
       endif
     else 
       if(bcbb_time_serie_type .eq. 2 ) then
         call get_date (bcbb_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(bcbb_entry)
           if (dayspmn /= 29) then
             bcbb_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             bcbb_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           bcbb_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         bcbb_time = model_time
       endif
     endif
     call obtain_interpolator_time_slices   &
                       (bcbb_aerosol_interp, bcbb_time)
   endif

   if ( trim(bcbf_source).ne. ' ') then
!--------------------------------------------------------------------
!    define the time in the bcbf data set from which data is to be 
!    taken. if bcbf is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(bcbf_time_serie_type .eq. 3) then
       if (bcbf_negative_offset) then
         bcbf_time = model_time - bcbf_offset
       else
         bcbf_time = model_time + bcbf_offset
       endif
     else 
       if(bcbf_time_serie_type .eq. 2 ) then
         call get_date (bcbf_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(bcbf_entry)
           if (dayspmn /= 29) then
             bcbf_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             bcbf_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           bcbf_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         bcbf_time = model_time
       endif
     endif
     call obtain_interpolator_time_slices   &
                 (bcbf_aerosol_interp, bcbf_time)
   endif

   if ( trim(bcsh_source).ne. ' ') then
!--------------------------------------------------------------------
!    define the time in the bcsh data set from which data is to be 
!    taken. if bcsh is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(bcsh_time_serie_type .eq. 3) then
       if (bcsh_negative_offset) then
         bcsh_time = model_time - bcsh_offset
       else
         bcsh_time = model_time + bcsh_offset
       endif
     else 
       if(bcsh_time_serie_type .eq. 2 ) then
         call get_date (bcsh_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(bcsh_entry)
           if (dayspmn /= 29) then
             bcsh_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             bcsh_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           bcsh_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         bcsh_time = model_time
       endif
     endif
     call obtain_interpolator_time_slices   &
                    (bcsh_aerosol_interp, bcsh_time)
   endif

   if ( trim(bcav_source).ne. ' ') then
!--------------------------------------------------------------------
!    define the time in the bcav data set from which data is to be 
!    taken. if bcav is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(bcav_time_serie_type .eq. 3) then
       if (bcav_negative_offset) then
         bcav_time = model_time - bcav_offset
       else
         bcav_time = model_time + bcav_offset
       endif
     else 
       if(bcav_time_serie_type .eq. 2 ) then
         call get_date (bcav_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(bcav_entry)
           if (dayspmn /= 29) then
             bcav_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             bcav_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           bcav_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         bcav_time = model_time
       endif
     endif
     call obtain_interpolator_time_slices   &
                     (bcav_aerosol_interp, bcav_time)
   endif

   if ( trim(omff_source).ne. ' ') then
!--------------------------------------------------------------------
!    define the time in the omff data set from which data is to be 
!    taken. if omff is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(omff_time_serie_type .eq. 3) then
       if (omff_negative_offset) then
         omff_time = model_time - omff_offset
       else
         omff_time = model_time + omff_offset
       endif
     else 
       if(omff_time_serie_type .eq. 2 ) then
         call get_date (omff_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(omff_entry)
           if (dayspmn /= 29) then
             omff_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             omff_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           omff_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         omff_time = model_time
       endif
     endif
     call obtain_interpolator_time_slices   &
                     (omff_aerosol_interp, omff_time)
   endif

   if ( trim(ombb_source).ne. ' ') then
!--------------------------------------------------------------------
!    define the time in the ombb data set from which data is to be 
!    taken. if ombb is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(ombb_time_serie_type .eq. 3) then
       if (ombb_negative_offset) then
         ombb_time = model_time - ombb_offset
       else
         ombb_time = model_time + ombb_offset
       endif
     else 
       if(ombb_time_serie_type .eq. 2 ) then
         call get_date (ombb_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(ombb_entry)
           if (dayspmn /= 29) then
             ombb_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             ombb_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           ombb_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         ombb_time = model_time
       endif
     endif
     call obtain_interpolator_time_slices   &
                  (ombb_aerosol_interp, ombb_time)
   endif

   if ( trim(ombf_source).ne. ' ') then
!--------------------------------------------------------------------
!    define the time in the ombf data set from which data is to be 
!    taken. if ombf is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(ombf_time_serie_type .eq. 3) then
       if (ombf_negative_offset) then
         ombf_time = model_time - ombf_offset
       else
         ombf_time = model_time + ombf_offset
       endif
     else 
       if(ombf_time_serie_type .eq. 2 ) then
         call get_date (ombf_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(ombf_entry)
           if (dayspmn /= 29) then
             ombf_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             ombf_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           ombf_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         ombf_time = model_time
       endif
     endif
     call obtain_interpolator_time_slices   &
                 (ombf_aerosol_interp, ombf_time)
   endif

   if ( trim(omsh_source).ne. ' ') then
!--------------------------------------------------------------------
!    define the time in the omsh data set from which data is to be 
!    taken. if omsh is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(omsh_time_serie_type .eq. 3) then
       if (omsh_negative_offset) then
         omsh_time = model_time - omsh_offset
       else
         omsh_time = model_time + omsh_offset
       endif
     else 
       if(omsh_time_serie_type .eq. 2 ) then
         call get_date (omsh_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(omsh_entry)
           if (dayspmn /= 29) then
             omsh_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             omsh_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           omsh_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         omsh_time = model_time
       endif
     endif
     call obtain_interpolator_time_slices   &
                       (omsh_aerosol_interp, omsh_time)
   endif

   if ( trim(omna_source).ne. ' ') then
!--------------------------------------------------------------------
!    define the time in the omna data set from which data is to be 
!    taken. if omna is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(omna_time_serie_type .eq. 3) then
       if (omna_negative_offset) then
         omna_time = model_time - omna_offset
       else
         omna_time = model_time + omna_offset
       endif
     else 
       if(omna_time_serie_type .eq. 2 ) then
         call get_date (omna_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(omna_entry)
           if (dayspmn /= 29) then
             omna_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             omna_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           omna_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         omna_time = model_time
       endif
     endif
     call obtain_interpolator_time_slices   &
                                     (omna_aerosol_interp, omna_time)
   endif
   if ( trim(omss_source).ne. ' ') then
!--------------------------------------------------------------------
!    define the time in the omss data set from which data is to be 
!    taken. if omss is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(omss_time_serie_type .eq. 3) then
       if (omss_negative_offset) then
         omss_time = model_time - omss_offset
       else
         omss_time = model_time + omss_offset
       endif
     else 
       if(omss_time_serie_type .eq. 2 ) then
         call get_date (omss_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(omss_entry)
           if (dayspmn /= 29) then
             omss_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             omss_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           omss_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         omss_time = model_time
       endif
     endif
     call obtain_interpolator_time_slices   &
                                      (omss_aerosol_interp, omss_time)
  endif


end subroutine atmos_carbon_aerosol_time_vary 


!######################################################################

subroutine atmos_carbon_aerosol_endts 


   if ( trim(bcff_source) .ne. ' ') then
     call unset_interpolator_time_flag (bcff_aerosol_interp)
   endif

   if ( trim(bcbb_source).ne.' ') then
     call unset_interpolator_time_flag (bcbb_aerosol_interp)
   endif

   if ( trim(bcbf_source).ne. ' ') then
     call unset_interpolator_time_flag (bcbf_aerosol_interp)
   endif

   if ( trim(bcsh_source).ne. ' ') then
     call unset_interpolator_time_flag (bcsh_aerosol_interp)
   endif

   if ( trim(bcav_source).ne. ' ') then
     call unset_interpolator_time_flag (bcav_aerosol_interp)
   endif

   if ( trim(omff_source).ne. ' ') then
     call unset_interpolator_time_flag (omff_aerosol_interp)
   endif

   if ( trim(ombb_source).ne. ' ') then
     call unset_interpolator_time_flag (ombb_aerosol_interp)
   endif

   if ( trim(ombf_source).ne. ' ') then
     call unset_interpolator_time_flag (ombf_aerosol_interp)
   endif

   if ( trim(omsh_source).ne. ' ') then
     call unset_interpolator_time_flag (omsh_aerosol_interp)
   endif

   if ( trim(omna_source).ne. ' ') then
     call unset_interpolator_time_flag (omna_aerosol_interp)
   endif

   if ( trim(omss_source).ne. ' ') then
     call unset_interpolator_time_flag (omss_aerosol_interp)
  endif


end subroutine atmos_carbon_aerosol_endts 


!######################################################################


!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME="atmos_carbon_aerosol_end">
!<OVERVIEW>
!  The destructor routine for the carbon module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits. 
! </DESCRIPTION>
!<TEMPLATE>
! call atmos_carbon_aero_end
!</TEMPLATE>
 subroutine atmos_carbon_aerosol_end
      call interpolator_end ( bcff_aerosol_interp)
      call interpolator_end ( bcbb_aerosol_interp)
      call interpolator_end ( bcbf_aerosol_interp)
      call interpolator_end ( bcsh_aerosol_interp)
      call interpolator_end ( bcav_aerosol_interp)
      call interpolator_end ( omff_aerosol_interp)
      call interpolator_end ( ombb_aerosol_interp)
      call interpolator_end ( ombf_aerosol_interp)
      call interpolator_end ( omsh_aerosol_interp)
      call interpolator_end ( omna_aerosol_interp)
      call interpolator_end ( omss_aerosol_interp)
      module_is_initialized = .FALSE.

 end subroutine atmos_carbon_aerosol_end
!</SUBROUTINE>
!#######################################################################
end module atmos_carbon_aerosol_mod

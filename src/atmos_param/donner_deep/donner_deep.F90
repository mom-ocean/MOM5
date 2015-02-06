                       module donner_deep_mod

use donner_types_mod,       only: donner_initialized_type, &
                                  donner_save_type, donner_rad_type, &
                                  donner_nml_type, donner_param_type, &
                                  donner_budgets_type, &
                                  donner_column_diag_type, &
                                  MAXMAG, MAXVAL, MINMAG, MINVAL, &
                                  DET_MASS_FLUX, MASS_FLUX,  &
                                  CELL_UPWARD_MASS_FLUX, TEMP_FORCING, &
                                  MOIST_FORCING, PRECIP,  FREEZING, &
                                  RADON_TEND, &
                                  donner_conv_type, donner_cape_type, &
                                  donner_cem_type
use  conv_utilities_k_mod,  only: sd_init_k, sd_end_k, ac_init_k,  &
                                  ac_end_k, uw_params_init_k, &
                                  exn_init_k, exn_end_k, findt_init_k, &
                                  findt_end_k, &
                                  adicloud, sounding, uw_params
use  conv_plumes_k_mod,     only: cp_init_k, cp_end_k, ct_init_k,  &
                                  ct_end_k, cplume, ctend
use fms_donner_mod,         only: fms_donner_process_nml,   &
                                  fms_donner_process_tracers, &
                                  fms_donner_activate_diagnostics, &
                                  fms_donner_col_diag, &
                                  fms_donner_column_control, &
                                  fms_donner_read_restart, &
                                  fms_donner_write_restart, &
                                  fms_donner_deep_netcdf, &
                                  fms_get_pe_number, &
                                  fms_sat_vapor_pres, &
                                  fms_error_mesg, fms_constants, &
                                  fms_close_col_diag_units, &
                                  fms_deallocate_variables, &
                                  fms_donner_process_monitors
use nonfms_donner_mod,      only: nonfms_donner_process_nml,   &
                                  nonfms_donner_process_tracers, &
                                  nonfms_donner_activate_diag, &
                                  nonfms_donner_col_diag, &
                                  nonfms_donner_column_control, &
                                  nonfms_donner_read_restart, &
                                  nonfms_donner_write_restart, &
                                  nonfms_donner_deep_netcdf, &
                                  nonfms_get_pe_number, &
                                  nonfms_sat_vapor_pres, &
                                  nonfms_error_mesg, nonfms_constants,&
                                  nonfms_deallocate_variables, &
                                  nonfms_donner_process_monitors, &
                                  nonfms_close_col_diag_units

implicit none
private

!--------------------------------------------------------------------
!        donner_deep_mod diagnoses the location and computes the 
!        effects of deep convection on the model atmosphere
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------


character(len=128)  :: version =  '$Id: donner_deep.F90,v 20.0 2013/12/13 23:17:16 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!--------------------------------------------------------------------
!---interfaces------

public   &
        donner_deep_init, donner_deep, donner_deep_end,  &
        donner_deep_restart, donner_deep_time_vary, donner_deep_endts

private   & 
        deallocate_variables


!---------------------------------------------------------------------
!---namelist----



!--------------------------------------------------------------------
!--- public data ----------




!--------------------------------------------------------------------
!----private data-----------



!---------------------------------------------------------------------
!  parameters stored in the donner_param derived type variable to facili-
!  tate passage to kernel subroutines:
!

integer,                   &
  parameter                  &
             ::  KPAR=7         
                        ! number of members in cumulus ensemble
integer,                   &
  parameter                  &
             ::  NLEV_HIRES=100       
                        ! number of levels in cloud model
real,                       &
  parameter                 &
             ::  PDEEP_CV = 500.e02 
                        ! minimum pressure difference between level of 
                        ! free convection and level of zero buoyancy 
                        ! needed for deep convection to occur [ Pa ].
real,                       &
  parameter                 &
             ::  MAX_ENTRAINMENT_CONSTANT_GATE = 0.0915
                        ! entrainment constant based on gate data for 
                        ! most entraining ensemble member
real,                       &
  parameter                 &
             ::  MAX_ENTRAINMENT_CONSTANT_KEP  = 0.0915
                        ! entrainment constant based on kep data for most
                        ! entraining ensemble member
real,                       &
  parameter,                 &
  dimension(KPAR)              &
             ::  ENSEMBLE_ENTRAIN_FACTORS_KEP  = (/ 1.0, 1.22, 1.56,  &
                                                    2.05, 2.6, 3.21,   &
                                                    7.84 /)
                        ! ratio of entrainment constant between ensemble
                        ! member 1 and ensemble member i for kep-based
                        ! ensemble
real,                       &
  parameter                 &
             ::  CLD_BASE_VERT_VEL = 0.5                          
                        ! vertical velocity assumed present at cloud base
                        ! [ m / sec ]
real,                       &
  parameter                 &
             ::  PSTOP = 40.0e02   
                        ! lowest possible pressure to which a cloud may 
                        ! extend in the cloud model [ Pa ]
real,                       &
  parameter                 &
             ::  PARCEL_DP = -1.0e02
                        ! pressure increment used for parcel calculations
                        ! [ Pa ]
real,                       &
  parameter                 &
             ::  UPPER_LIMIT_FOR_LFC = 500.e02 
                        ! lowest pressure allowed for level of free conv-
                        ! ection [ Pa ]
real,                       &
  parameter                 &
             ::  DP_OF_CLOUD_MODEL = -10.e02
                        ! pressure thickness (Pa) of the layers in the
                        ! donner parameterization's cloud model.
real,                       &
  parameter                 &
             ::  CLOUD_BASE_RADIUS = 1000.
                        ! radius assumed for cloud ensemble member #1 at
                        ! cloud base [ m ]
real,                       &
  parameter                 &
             ::  WDET = .1   
                        ! vertical velocity at which detrainment from the
                        ! clouds begins [ m/s ]
real,                       &
  parameter                 &
             ::  RBOUND = 0.01    
                        ! value of cumulus radius at which cloud effect-
                        ! ively disappears and cloud model calculation 
                        ! stops [ m ]
real,                       &
  parameter                 &
             ::  WBOUND = 0.01  
                        ! value of cumulus vertical velocity at which 
                        ! cloud model calculation stops [ m / sec ]
real,                       &
  parameter                 &
             ::  FREEZE_FRACTION = 0.52
                        ! fraction of liquid in cloud updraft which may 
                        ! be frozen. (Leary and Houze (JAS,1980)) 
                        ! [ dimensionless ]
real,                       &
  parameter                 &
             ::  VIRT_MASS_CO = 0.5
                        ! virtual mass coefficient [ dimensionless ]
real,                       &
  parameter                 &
             ::  PDEEP_MC = 200.e02 
                        ! pressure thickness [ Pa ] required for meso-
                        ! scale circulation. It refers to the least
                        ! penetrative ensemble member. For this check 
                        ! to function properly, the entrainment coeffic-
                        ! ient in cloud_model for kou=1 must be the 
                        ! largest entrainment coefficient.
real,                       &
  parameter                 &
             ::  TR_INSERT_TIME = 0.0
                        ! fractional point (based on mass increase) 
                        ! during a timestep at which an entraining parcel
                        ! takes on internally-generated tracer 
                        ! [ dimensionless, value between 0.0 and 1.0 ]
real,                       &
  parameter                 &
             ::  AUTOCONV_RATE = 1.0e-03
                        ! rate of autoconversion of cloud to rainwater 
                        ! [ sec**(-1) ]
real,                       &
  parameter                 &
             ::  AUTOCONV_THRESHOLD =  0.5    
                        ! threshold of cloud water at which autoconver-
                        ! sion of cloud to rainwater begins  [ g / m**3 ]
real,                       &
  parameter                 &
             ::  TFRE = 258.  
                        ! temperature at which cloud liquid begins to 
                        ! freeze [ deg K ]
real,                       &
  parameter                 &
             ::  DFRE = 10.   
                        ! range of temperature between the onset and 
                        ! completion of freezing  [ deg K ]
real,                       &
  parameter                 &
             ::  UPPER_LIMIT_FOR_LCL = 500.0E02
                        ! lowest pressure allowable for lifting condens-
                        ! ation level; deep convection will not be pres-
                        ! ent if lcl not reached before this pressure 
                        ! [ Pa ]
integer,                   &
  parameter         &
             ::  ISTART = 1    
                        ! index of level in cape grid from which the 
                        ! parcel originates for the cape calculations
real,                       &
  parameter                 &
             ::  TMIN = 154.       
                        ! cape calculations are terminated when parcel 
                        ! temperature goes below TMIN [ deg K ]
real,                       &
  parameter                 &
             ::  MESO_LIFETIME = 64800.
                        ! assumed lifetime of mesoscale circulation 
                        ! (from Leary and Louze, 1980) [ sec ]
real,                       &
  parameter                 &
             ::  MESO_REF_OMEGA = -0.463
                        ! assumed reference omega for mesoscale updraft 
                        ! (from Leary and Louze, 1980) [ Pa / sec ]
real,                       &
  parameter                 &
             ::  TPRIME_MESO_UPDRFT = 1.0    
                        ! assumed temperature excess of mesoscale updraft
                        ! over its environment [ deg K ]
real,                       &
  parameter                 &
             ::  MESO_SEP = 200.0E+02
                        ! pressure separation between base of mesoscale
                        ! updraft and top of mesoscale downdraft [ Pa ]
real,                       &
  parameter                 &
             ::  REF_PRESS = 1.0E05
                        ! reference pressure used in calculation of exner
                        ! fumction [ Pa ]
real,                       &
  parameter                 &
             ::  R_CONV_LAND  = 10.0    
                        ! assumed convective cloud droplet radius over 
                        ! land [ microns ]   
real,                       &
  parameter                 &
             ::  R_CONV_OCEAN = 16.0  
                        ! assumed convective cloud droplet radius over 
                        ! ocean [ microns ]   
real,                       &
  parameter                 &
             ::  N_LAND = 600*1.0e6 
                        ! assumed droplet number conc over land (m**-3)
real,                       &
  parameter                 &
             ::  N_OCEAN = 150*1.0e6 
                        ! assumed droplet number conc over ocean (m**-3)
real,                       &
  parameter                 &
             ::  DELZ_LAND = 500.0   
                        ! assumed cloud depth over land (m) 
real,                       &
  parameter                 &
             ::  DELZ_OCEAN = 1500.0   
                        ! assumed cloud depth over ocean (m)
real,                       &
  parameter                 &
             ::  CELL_LIQUID_EFF_DIAM_DEF = 15.0    
                        ! default cell liquid eff diameter [ microns ]
real,                       &
  parameter                 &
             ::  CELL_ICE_GENEFF_DIAM_DEF = 18.6   
                        ! default cell ice generalized effective diameter
                        ! [ microns ]
integer,                   &
  parameter         &
             ::  ANVIL_LEVELS = 6  
                        ! number of levels assumed to be in anvil clouds
real,                       &
  parameter,                &
  dimension(ANVIL_LEVELS)   &
             ::  DGEICE  = (/ 38.5, 30.72, 28.28, 25.62, 24.8, 13.3 /)
                        ! generalized effective size of hexagonal ice 
                        ! crystals, defined as in Fu (1996, J. Clim.) 
                        ! values from Table 2 of McFarquhar et al. 
                        ! (1999, JGR) are averaged over all grid boxes 
                        ! for which D_ge is defined for all altitudes 
                        ! between 9.9 and 13.2 km. index 1 at bottom of 
                        ! anvil
real,                       &
  parameter,                &
  dimension(ANVIL_LEVELS)   &
             ::  RELHT  =  (/0.0, 0.3, 0.45, 0.64, 0.76, 1.0/)
                        ! distance from anvil base, normalized by total 
                        ! anvil thickness. from Table 2 of McFarquhar et
                        ! al. (1999, JGR) for grid boxes with data 
                        ! between 9.9 and 13.2 km. index 1 at anvil 
                        ! bottom

integer,                      &
  parameter                   &
             ::  N_WATER_BUDGET = 9
                        ! number of terms in vapor budget

integer,                      &
  parameter                   &
             ::  N_ENTHALPY_BUDGET =  19 
                        ! number of terms in enthalpy budget

integer,               &
  parameter            &
             ::  N_PRECIP_PATHS = 5
                        ! number of paths precip may take from 
                        ! condensing until it reaches the ground
                        ! (liquid; liquid which freezes; liquid which
                        ! freezes and then remelts; ice; ice which 
                        ! melts)

integer,               &
  parameter            &
             ::  N_PRECIP_TYPES = 3
                        ! number of precip types (cell, cell condensate
                        ! tranmsferred to mesoscale circulation,
                        ! mesoscale condensation and deposition)





!---------------------------------------------------------------------
!    derived type variables present for duration of job:
!    (see donner_types.h for documentation of their contents)
!

type(donner_param_type),       save :: Param
type(donner_column_diag_type), save :: Col_diag
type(donner_nml_type),         save :: Nml
type(donner_save_type),        save :: Don_save
type(donner_initialized_type), save :: Initialized
 
type(uw_params),               save :: Uw_p

logical                             :: calc_conv_on_this_step

!-----------------------------------------------------------------------
!   miscellaneous variables
!
!     module_is_initialized       module has been initialized ?
!

logical :: module_is_initialized = .false. 


logical :: running_in_fms = .true.

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                   PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#####################################################################

subroutine donner_deep_init (lonb, latb, pref, axes, secs, days, &
                             tracers_in_donner, do_conservation_checks,&
                             using_unified_closure, using_fms_code)

!---------------------------------------------------------------------
!    donner_deep_init is the constructor for donner_deep_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
real,            dimension(:,:), intent(in)   :: lonb, latb
real,            dimension(:),   intent(in)   :: pref
integer,         dimension(4),   intent(in)   :: axes
integer,                         intent(in)   :: secs, days
logical,         dimension(:),   intent(in)   :: tracers_in_donner
logical,                         intent(in)   :: do_conservation_checks
logical,                         intent(in)   :: using_unified_closure
logical,                         intent(in), optional :: &
                                                 using_fms_code

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      lonb         array of model longitudes on cell corners     
!                   [ radians ]
!      latb         array of model latitudes on cell corners   
!                   [ radians ]
!      pref         array of reference pressures at full levels (plus 
!                   surface value at nlev+1), based on 1013.25 hPa pstar
!                   [ Pa ]
!      axes         data axes for diagnostics
!      Time         current time [ time_type ]
!      tracers_in_donner 
!                   logical array indicating which of the activated 
!                   tracers are to be transported by donner_deep_mod
!
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables:

      integer                             :: idf, jdf, nlev, ntracers
      character(len=200)                  :: ermesg
      integer                             :: erflag
      integer                             :: me, root_pe
  
!-------------------------------------------------------------------
!  local variables:
!
!     idf                    number of columns in the x dimension on the
!                            processors domain
!     jdf                    number of columns in the y dimension on the
!                            processors domain
!     nlev                   number of model layers 
!     ntracers               number of tracers to be transported by
!                            the donner deep convection parameterization
!                         
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!    initialize error message and error flag.
!-------------------------------------------------------------------
      ermesg = '  '
      erflag = 0

!---------------------------------------------------------------------
!    define variable to indicated whether this module is being executed
!    within the FMS infrastructure. by default it is.
!---------------------------------------------------------------------
      if (present (using_fms_code)) then
        if ( .not. using_fms_code) then
          running_in_fms = .false.
        endif
      endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    1. READ NAMELIST AND WRITE IT TO LOG FILE.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      if (running_in_fms) then
        call fms_donner_process_nml (Nml, kpar)
      else 

!---------------------------------------------------------------------
!    for the nonfms case, appropriate code to read the namelist should 
!    be included in nonfms_process_nml.  the current routine (without 
!    such code) modifies nml values with source assignment statements,
!    as needed.
!---------------------------------------------------------------------
        call nonfms_donner_process_nml (Nml, kpar)
      endif 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    2. DO CONSISTENCY / VALIDITY TESTS ON NML AND PARAMETER VARIABLES.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      if (Nml%do_donner_plume) then
        Nml%do_hires_cape_for_closure = .true.
      endif
 
      if (.not. Nml%do_donner_cape .and. &
           Nml%rmuz_for_cape /= 0.0) then
        erflag = 1
        ermesg =  'donner_deep_init: &
            &  a non-zero rmuz_for_cape is allowed only when &
            & do_donner_cape is true'
      endif

      if (Nml%do_donner_cape .and.     &
           .not. Nml%do_hires_cape_for_closure) then
        erflag = 1
        ermesg =  'donner_deep_init: &
            & do_hires_cape_for_closure must be .true. when &
            & do_donner_cape is true'
      endif

      if (.not. (Nml%do_donner_cape) .and.  &
           .not. (Nml%do_hires_cape_for_closure)) then

        if (Nml%rmuz_for_closure /= 0.0) then
          erflag = 1
          ermesg =  'donner_deep_init: &
          &  a non-zero rmuz_for_closure is currently implemented &
          & only for the hi-res cape closure calculation'
        endif
 
        if (trim(Nml%entrainment_scheme_for_closure) /= 'none') then
          erflag = 1
          ermesg =  'donner_deep_init: &
          &  entrainment in the closure calculation is currently &
           & implemented only for the hi-res cape closure calculation'
        endif
              
        if (Nml%modify_closure_plume_condensate ) then
          erflag = 1
          ermesg =  'donner_deep_init: &
          &  condensate modification in the closure calculation is &
          & currently implemented only for the hi-res cape closure &
          & calculation'
        endif
              
        if (Nml%closure_plume_condensate /= -999. ) then
          erflag = 1
          ermesg =  'donner_deep_init: &
          &  condensate modification in the closure calculation is &
          & currently implemented only for the hi-res cape closure &
          & calculation'
        endif
      endif
              
      if (Nml%do_donner_cape .and. Nml%gama /= 0.0) then
        erflag = 1
        ermesg =  'donner_deep_init: &
            & gama must be 0.0 if do_donner_cape is .true.; code for &
            & gama /=  0.0 not yet implemented'
      endif
      if (Nml%deep_closure /= 0 .and. Nml%deep_closure /= 1 ) then
        erflag = 1
        ermesg =  'donner_deep_init: &
              & deep_closure must be 0 or 1; code for &
              & cu_clo_miz not yet implemented'
      endif
      if (Nml%do_rh_trig .and. Nml%do_donner_closure) then
        erflag = 1
        ermesg =  'donner_deep_init: &
             & do_rh_trig must be .false. for donner full  &
               &parameterization; its use not yet implemented'
      endif

!---------------------------------------------------------------------
!    check for a valid value of donner_deep_freq. 
!---------------------------------------------------------------------
      if (Nml%donner_deep_freq > 86400) then
        erflag = 1
        ermesg = 'donner_deep_init: &
         & donner convection must be called at least once per day'
      else if (Nml%donner_deep_freq <= 0) then
        erflag = 1
        ermesg = 'donner_deep_init: &
          & a positive value must be assigned to donner_deep_freq'
      endif

!---------------------------------------------------------------------
!    check for valid value of entrainment_constant_source.
!---------------------------------------------------------------------
      if (trim(Nml%entrainment_constant_source) == 'gate' .or. &
          trim(Nml%entrainment_constant_source) == 'kep' ) then
      else
        erflag = 1
        ermesg = 'donner_deep_init: &
         & invalid string for nml variable entrainment_constant_source'
      endif

!---------------------------------------------------------------------
!    test that PSTOP is smaller than UPPER_LIMIT_FOR_LFC.
!---------------------------------------------------------------------
      if (pstop > upper_limit_for_lfc) then
        erflag = 1
        ermesg = 'donner_deep_init: &
           & pstop must be above the upper limit of &
                                &the level of free convection'
      endif

!---------------------------------------------------------------------
!    test that cell_liquid_size_type has been validly specified, and if
!    it is specified as 'input', an appropriate input value has been
!    supplied.
!---------------------------------------------------------------------
      if (trim(Nml%cell_liquid_size_type) == 'input') then
        Initialized%do_input_cell_liquid_size = .true.
        Initialized%do_bower_cell_liquid_size = .false.
        if (Nml%cell_liquid_eff_diam_input < 0.0) then
          erflag = 1
          ermesg = 'donner_deep_init: &
            & cell liquid size must be input, but no value supplied'
        endif
      else if (trim(Nml%cell_liquid_size_type) == 'bower') then
        Initialized%do_input_cell_liquid_size = .false.
        Initialized%do_bower_cell_liquid_size = .true.
      else
        erflag = 1
        ermesg = 'donner_deep_init: &
           & cell_liquid_size_type must be either input or bower'
      endif

!---------------------------------------------------------------------
!    test that cell_ice_size_type has been validly specified, and if
!    specified as 'input', that cell_ice_geneff_diam_input has also 
!    been appropriately defined.
!---------------------------------------------------------------------
      if (trim(Nml%cell_ice_size_type) == 'input') then
        Initialized%do_input_cell_ice_size = .true.
        Initialized%do_default_cell_ice_size = .false.
        if (Nml%cell_ice_geneff_diam_input <= 0.0) then
          erflag = 1
          ermesg  =  'donner_deep_init: must define a '// &
                     'nonnegative generalized effective'// &
                     'diameter for ice when cell_ice_size_type is input'
        endif
      else if (trim(Nml%cell_ice_size_type) == 'default') then
        Initialized%do_input_cell_ice_size = .false.
        Initialized%do_default_cell_ice_size = .true.
      else
        erflag = 1
        ermesg =  'donner_deep_init: cell_ice_size_type must ' //  &
                  'be input or default'
      endif

!---------------------------------------------------------------------
!    check for consistency between entrainment used in closure 
!    calculation. define logical indicating whether entrainment
!    coefficient is to be constant or ht-dependent.
!---------------------------------------------------------------------
      if (trim(Nml%entrainment_scheme_for_closure) == 'none' .and. &
                               Nml%rmuz_for_closure /= 0.0) then 
        erflag = 1
        ermesg = 'donner_deep_init: do not specify a non-zero ' // &
                  'rmuz_for_closure when no entrainment is desired'
      endif
      if (trim(Nml%entrainment_scheme_for_closure) ==    &
                                                 'ht-dependent' .and. &
                               Nml%rmuz_for_closure == 0.0) then 
        erflag = 1
        ermesg = 'donner_deep_init: must specify rmuz_for_closure ' // &
                  'when ht-dependent entrainment is desired'
      endif
      if (trim(Nml%entrainment_scheme_for_closure) ==    &
                                                 'ht-dependent' ) then
        Initialized%use_constant_rmuz_for_closure = .false.
      else
        Initialized%use_constant_rmuz_for_closure = .true.
      endif

!---------------------------------------------------------------------
!    check that if the closure plume condensate is to be modified that
!    a value is given.
!---------------------------------------------------------------------
      if (Nml%modify_closure_plume_condensate .and. &
          Nml%closure_plume_condensate == -999.) then
        erflag = 1
        ermesg = 'donner_deep_init: must specify ' // &
              'closure_plume_condensate when modification is requested'
      endif
        
!---------------------------------------------------------------------
!    if any errors were encountered, process them.
!---------------------------------------------------------------------
      if (erflag /= 0) then
        if (running_in_fms) then
          call fms_error_mesg (ermesg) 
        else

!---------------------------------------------------------------------
!    appropriate error processing code should be added in subroutine
!    nonfms_error_mesg. currently an error message is printed and a 
!    stop command issued (dangerous on parallel machines!).
!---------------------------------------------------------------------
          call nonfms_error_mesg (ermesg) 
        endif
      endif

!---------------------------------------------------------------------
!    place the logical input argument indicating whether the cloud 
!    base mass flux calculated by uw_conv_mod is also to be used 
!    in defining the closure for donner deep convection in the 
!    donner_initialized_type variable Initialized.
!    place the conservation check flag in the Initialized variable.
!---------------------------------------------------------------------
      Initialized%using_unified_closure = using_unified_closure
      Initialized%do_conservation_checks = do_conservation_checks

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    3. PROCESS TRACERS THAT ARE TO BE TRANSPORTED BY THE DONNER DEEP
!       CONVECTION PARAMETERIZATION.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    determine how many tracers are to be transported by donner_deep 
!    convection. allocate arrays to contain their names and units for use
!    with diagnostics and restarts. define a logical variable indicating
!    if any tracers are to be so transported. obtain the tracer names and
!    units.
!---------------------------------------------------------------------
      ntracers = count(tracers_in_donner)
      allocate ( Don_save%tracername   (ntracers) )
      allocate ( Don_save%tracer_units (ntracers) )
      allocate ( Initialized%wetdep(ntracers) )
      if (ntracers > 0) then
        if (running_in_fms) then
          call fms_donner_process_tracers (Initialized,  &
                                           tracers_in_donner, Don_save)
        else 

!----------------------------------------------------------------------
!    currently tracers are not supported in the nonfms case. if it is
!    desired to transport tracers with donner convection, then the
!    subroutine nonfms_donner_process_tracers should be set up to mimic
!    the functionality of fms_donner_process_tracers for each tracer
!    thus transported (see subroutine fms_donner_process_tracers in 
!    fms_donner.F90).
!----------------------------------------------------------------------
          call nonfms_donner_process_tracers 
        endif 
      else
        Initialized%do_donner_tracer = .false.
      endif
      

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    4. DEFINE PROCESSOR DIMENSIONS AND ALLOCATE SPACE FOR MODULE 
!       VARIABLES.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-------------------------------------------------------------------
!    define the grid dimensions. idf and jdf are the (i,j) dimensions of
!    the domain on this processor, nlev is the number of model layers.
!-------------------------------------------------------------------
      nlev = size(pref(:)) - 1
      idf  = size(lonb,1) - 1
      jdf  = size(latb,2) - 1

!--------------------------------------------------------------------
!    allocate module variables that will be saved across timesteps.
!    these are stored in the derived-type variable Don_save. see 
!    donner_types.h for description of these variables.
!--------------------------------------------------------------------
      allocate ( Don_save%cemetf             (idf, jdf, nlev ) )
      allocate ( Don_save%lag_temp           (idf, jdf, nlev ) )
      allocate ( Don_save%lag_vapor          (idf, jdf, nlev ) )
      allocate ( Don_save%lag_press          (idf, jdf, nlev ) )
      allocate ( Don_save%cememf             (idf, jdf, nlev ) )
      allocate ( Don_save%mass_flux          (idf, jdf, nlev ) )
      allocate ( Don_save%mflux_up           (idf, jdf, nlev ) )
      allocate ( Don_save%cell_up_mass_flux  (idf, jdf, nlev+1 ) )
      allocate ( Don_save%det_mass_flux      (idf, jdf, nlev ) )
      allocate ( Don_save%dql_strat          (idf, jdf, nlev ) )
      allocate ( Don_save%dqi_strat          (idf, jdf, nlev ) )
      allocate ( Don_save%dqa_strat          (idf, jdf, nlev ) )
      allocate ( Don_save%humidity_area      (idf, jdf, nlev ) )
      allocate ( Don_save%humidity_factor    (idf, jdf, nlev ) )
      allocate ( Don_save%tracer_tends       (idf, jdf, nlev, ntracers))
      allocate ( Don_save%parcel_disp        (idf, jdf ) )
      allocate ( Don_save%tprea1             (idf, jdf ) )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    4. INITIALIZE THE NETCDF OUTPUT VARIABLES.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!--------------------------------------------------------------------
!    activate the netcdf diagnostic fields.
!-------------------------------------------------------------------
      if (running_in_fms) then
        call fms_donner_activate_diagnostics (secs, days, axes, &
                  Don_save, Nml, n_water_budget, n_enthalpy_budget, &
                  n_precip_paths, n_precip_types, nlev_hires, kpar)
      else

!---------------------------------------------------------------------
!    subroutine  nonfms_donner_activate_diagnostics should be set up
!    to initialize the procedure needed to output netcdf variable
!    fields in the nonFMS model. by default, it currently does nothing.
!---------------------------------------------------------------------
        call nonfms_donner_activate_diag (secs, days, axes, &
                  Don_save, Nml, n_water_budget, n_enthalpy_budget, &
                  n_precip_paths, n_precip_types, nlev_hires, kpar)
      endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    5. PROCESS THE RESTART FILE.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      if (running_in_fms) then
        call fms_donner_read_restart (Initialized, ntracers,   &
                                      secs, days, Don_save, Nml)
      else

!---------------------------------------------------------------------
!    subroutine nonfms_donner_read_restart should be set up to handle 
!    the reading of the donner_deep.res.nc file in the nonFMS framework.
!    by default, it begins the model run from a coldstart.
!---------------------------------------------------------------------
        call nonfms_donner_read_restart (Initialized, ntracers,   &
                                      secs, days, Don_save, Nml)
      endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    6. INITIALIZE VARIABLES NEEDED FOR COLUMN_DIAGNOSTICS_MOD OUTPUT.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      if (running_in_fms) then
        call fms_donner_col_diag     (lonb, latb, Col_diag, pref)    
      else

!---------------------------------------------------------------------
!    subroutine nonfms_donner_col_diag should be set up to process
!    the column diagnostic output available from donner_deep_mod in 
!    the nonFMS framework. by default, it currently sets variables to
!    disallow that option; if that option is desired, the functionality
!    of fms_donner_col_diag needs to be added to that subroutine.
!---------------------------------------------------------------------
        call nonfms_donner_col_diag     (lonb, latb, Col_diag, pref)   
      endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    7. FILL THE DONNER_PARAM_TYPE VARIABLE WITH VALUES THAT HAVE BEEN 
!       DEFINED HERE.                  
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!----------------------------------------------------------------------
!    define the components of Param that come from constants_mod. see 
!    donner_types.h for their definitions.
!----------------------------------------------------------------------
      if (running_in_fms) then
        call fms_constants (Param)
      else

!---------------------------------------------------------------------
!    subroutine nonfms_constants is designed to obtain model constants
!    and place them in the donner_param_type variable Param. currently
!    values are assigned to them in this subroutine; this process should
!    be modified so that these constants are obtained from their natural
!    location in the nonFMS model.
!---------------------------------------------------------------------
        call nonfms_constants (Param)
      endif


!----------------------------------------------------------------------
!    store the parameters defined in this module into the 
!    donner_parameter_type variables Param. these variables are defined
!    above.
!----------------------------------------------------------------------
      Param%cp_vapor                = 4.0*Param%rvgas
      Param%parcel_dp               = PARCEL_DP
      Param%upper_limit_for_lfc     = UPPER_LIMIT_FOR_LFC
      Param%pstop                   = PSTOP
      Param%cld_base_vert_vel       = CLD_BASE_VERT_VEL
      Param%dp_of_cloud_model       = DP_OF_CLOUD_MODEL
      Param%cloud_base_radius       = CLOUD_BASE_RADIUS
      Param%wdet                    = WDET
      Param%rbound                  = RBOUND
      Param%wbound                  = WBOUND
      Param%freeze_fraction         = FREEZE_FRACTION
      Param%virt_mass_co            = VIRT_MASS_CO
      Param%pdeep_mc                = PDEEP_MC
      Param%tr_insert_time          = TR_INSERT_TIME
      Param%autoconv_rate           = AUTOCONV_RATE
      Param%autoconv_threshold      = AUTOCONV_THRESHOLD
      Param%tfre                    = TFRE
      Param%dfre                    = DFRE
      Param%evap_in_downdrafts      = Nml%EVAP_IN_DOWNDRAFTS
      Param%evap_in_environ         = Nml%EVAP_IN_ENVIRON      
      Param%entrained_into_meso     = Nml%ENTRAINED_INTO_MESO
      Param%d622                    = Param%rdgas/Param%rvgas
      Param%d608                    = Param%rvgas/Param%rdgas - 1.0
      Param%upper_limit_for_lcl     = UPPER_LIMIT_FOR_LCL
      Param%tmin                    = TMIN
      Param%anvil_precip_efficiency = Nml%ANVIL_PRECIP_EFFICIENCY
      Param%meso_lifetime           = MESO_LIFETIME
      Param%meso_ref_omega          = MESO_REF_OMEGA
      Param%tprime_meso_updrft      = TPRIME_MESO_UPDRFT
      Param%meso_sep                = MESO_SEP
      Param%ref_press               = REF_PRESS
      Param%meso_down_evap_fraction = Nml%MESO_DOWN_EVAP_FRACTION
      Param%meso_up_evap_fraction   = Nml%MESO_UP_EVAP_FRACTION
      Param%istart                  = ISTART

      Param%max_entrainment_constant_gate =   &
                                           MAX_ENTRAINMENT_CONSTANT_GATE
      Param%max_entrainment_constant_kep  = MAX_ENTRAINMENT_CONSTANT_KEP
      Param%pdeep_cv                      = PDEEP_CV
      Param%cdeep_cv                      = Nml%CDEEP_CV
      Param%kpar                          = KPAR
      Param%r_conv_land                   = R_CONV_LAND
      Param%r_conv_ocean                  = R_CONV_OCEAN 
      Param%n_land                        = N_LAND
      Param%n_ocean                       = N_OCEAN
      Param%delz_land                     = DELZ_LAND
      Param%delz_ocean                    = DELZ_OCEAN
      Param%cell_liquid_eff_diam_def      = CELL_LIQUID_EFF_DIAM_DEF 
      Param%cell_ice_geneff_diam_def      = CELL_ICE_GENEFF_DIAM_DEF
      Param%anvil_levels                  = ANVIL_LEVELS 

      allocate (Param%arat(kpar))
      allocate (Param%ensemble_entrain_factors_gate(kpar))
      allocate (Param%ensemble_entrain_factors_kep(kpar))
      Param%arat                          = Nml%ARAT
      Param%ensemble_entrain_factors_gate =   &
                                       Nml%ensemble_entrain_factors_gate
      Param%ensemble_entrain_factors_kep  = ENSEMBLE_ENTRAIN_FACTORS_KEP

      allocate (Param%dgeice (ANVIL_LEVELS))
      allocate (Param%relht  (ANVIL_LEVELS))
      Param%dgeice  = DGEICE               
      Param%relht   = RELHT                

!---------------------------------------------------------------------
!    initialize the kernelized modules needed outside of the donner 
!    directory.
!---------------------------------------------------------------------
      if (running_in_fms) then
        call fms_sat_vapor_pres
      else

!--------------------------------------------------------------------
!    this routine is reserved for any initialization involved with the
!    saturation vapor pressure calculation in the nonFMS model. For
!    test purposes, this routine currently uses the FMS routines, but
!    this should be changed to use the nonFMS model's procedures.
!--------------------------------------------------------------------
        call nonfms_sat_vapor_pres
      endif

      if (running_in_fms) then
        call fms_get_pe_number(me, root_pe)
      else

!---------------------------------------------------------------------
!    subroutine nonfms_get_pe_number should be set up to return the
!    current pe's number, if column daignostics are activated in the
!    nonFMS model. By default, the subroutine returns 0 for all pes.
!---------------------------------------------------------------------
        call nonfms_get_pe_number(me, root_pe)
      endif

      call uw_params_init_k (Param%hlv, Param%hls, Param%hlf, &
          Param%cp_air, Param%grav, Param%kappa, Param%rdgas,  &
          Param%ref_press, Param%d622,Param%d608, Param%kelvin -160., & 
          Param%kelvin + 100. , me, root_pe, Uw_p)


      call exn_init_k (Uw_p)
      call findt_init_k (Uw_p)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    9. SET UP CODE TO MONITOR SELECTED OUTPUT VARIABLES.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      if (running_in_fms) then
        call fms_donner_process_monitors (idf, jdf, nlev, ntracers, &
                               axes, secs, days, Initialized, Don_save)
      else

!---------------------------------------------------------------------
!    subroutine nonfms_donner_process_monitors should be set up to 
!    handle the processing of variable monitopr output, if that capa-
!    bility is desired. Currently the subroutine does nothing; see 
!    subroutine fms_donner_process_monitors in fms_donner.F90 for the
!    functionality required to activate this option.
!---------------------------------------------------------------------
        call nonfms_donner_process_monitors 
      endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!   10. END OF SUBROUTINE.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!--------------------------------------------------------------------
!    set flag to indicate that donner_deep_mod has been initialized.
!--------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------



end subroutine donner_deep_init

!###################################################################
         
subroutine donner_deep_time_vary (dt)  
                                  
real, intent(in) :: dt

!--------------------------------------------------------------------
!    decrement the time remaining before the convection calculations. 
!    save the current model physics timestep.
!--------------------------------------------------------------------
      Initialized%conv_alarm  = Initialized%conv_alarm - int(dt)
      Initialized%physics_dt = int(dt)
 
!--------------------------------------------------------------------
!    set a flag to indicate whether the convection calculation is to be 
!    done on this timestep. if this is the first call to donner_deep 
!    (i.e., coldstart), convection cannot be calculated because the
!    lag profiles needed to calculate cape are unavailable, and so
!    a time tendency of cape can not be obtained. otherwise, it is a
!    calculation step or not dependent on whether the convection "alarm"
!    has gone off. 
!---------------------------------------------------------------------
      if (Initialized%coldstart) then
        calc_conv_on_this_step = .false.
      else
        if (Initialized%conv_alarm <= 0) then
          calc_conv_on_this_step = .true.
        else
          calc_conv_on_this_step = .false.
        endif
      endif

!--------------------------------------------------------------------


end subroutine donner_deep_time_vary



!###################################################################
 
subroutine donner_deep_endts

!---------------------------------------------------------------------
!    if this was the first time through the parameterization, set
!    the flag so indicating (coldstart) to be .false.. if this was a 
!    calculation step, set the alarm to define the next time at which 
!    donner convection is to be executed.
!----------------------------------------------------------------------
      if (Initialized%coldstart) Initialized%coldstart = .false.
      if (calc_conv_on_this_step) then
        Initialized%conv_alarm = Initialized%conv_alarm +    &
                                                 Nml%donner_deep_freq
      endif

!--------------------------------------------------------------------


end subroutine donner_deep_endts




!###################################################################

subroutine donner_deep (is, ie, js, je, dt, temp, mixing_ratio, pfull, &
                        phalf, zfull, zhalf, omega, pblht, tkemiz, &
                        qstar, cush, coldT, land, sfc_sh_flux,  &
                        sfc_vapor_flux, tr_flux, tracers, secs, days, &
                        cbmf, cell_cld_frac,  &
                        cell_liq_amt, cell_liq_size, cell_ice_amt,   &
                        cell_ice_size, cell_droplet_number, &
                        meso_cld_frac, meso_liq_amt, &
                        meso_liq_size, meso_ice_amt, meso_ice_size,  &
                        meso_droplet_number, &
                        nsum, precip, delta_temp, delta_vapor, detf, &
                        uceml_inter, mtot, mfluxup, mhalf_3d, &
                        donner_humidity_area,    &
                        donner_humidity_factor, qtrtnd, donner_wetdep,&
                        lheat_precip, vert_motion,        &
                        total_precip, liquid_precip, frozen_precip, &
                        frz_meso, liq_meso, frz_cell, liq_cell, &
                        qlin, qiin, qain,              &      ! optional
                        delta_ql, delta_qi, delta_qa)         ! optional
                        
!-------------------------------------------------------------------
!    donner_deep is the prognostic driver subroutine of donner_deep_mod.
!    it takes as input the temperature (temp), vapor mixing ratio 
!    (mixing_ratio), pressure at full and half-levels (pfull, phalf),
!    vertical velocity at full levels (omega), the large scale cloud 
!    variables (qlin, qiin, qain), the land fraction (land),  the heat 
!    (sfc_sh_flux) , moisture (sfc_vapor_flux) and tracer (tr_flux) 
!    fluxes across the surface that are to be seen by this parameter-
!    ization, the tracers to be transported by the donner convection
!    parameterization (tracers), and the current time (as time_type 
!    variable Time). the routine returns the precipitation (precip),
!    increments to the temperature (delta_temp) and mixing ratio 
!    (delta_vapor), the detrained mass flux (detf), upward cell mass 
!    flux at interface levels  (uceml_inter) and total mass flux at full
!    levels (mtot), two arrays needed to connect the donner convection 
!    and strat cloud parameterizations (donner_humidity_area, 
!    donner_humidity_ratio), increments to the cloudwater (delta_ql), 
!    cloudice (delta_qi) and cloud area (delta_qa) fields and tendencies
!    for those tracers that are to be transported by the donner convect-
!    ion parameterization (qtrtnd). there are an additional eleven arrays
!    defining the donner scheme cloud characteristics needed by the rad-
!    iation package, which are passed in and updated on donner calcul-
!    ation steps.
!-------------------------------------------------------------------

!--------------------------------------------------------------------
integer,                      intent(in)    :: is, ie, js, je
real,                         intent(in)    :: dt
real, dimension(:,:,:),       intent(in)    :: temp, mixing_ratio, &
                                               pfull, phalf, zfull, zhalf, omega
real, dimension(:,:),         intent(in)    :: pblht, tkemiz, qstar,cush
real, dimension(:,:),         intent(in)    :: land
logical, dimension(:,:),      intent(in)    :: coldT
real, dimension(:,:),         intent(in)    :: sfc_sh_flux, &
                                               sfc_vapor_flux
real, dimension(:,:,:),       intent(in)    :: tr_flux 
real, dimension(:,:,:,:),     intent(in)    :: tracers 
integer,                      intent(in)    :: secs, days
real, dimension(:,:),         intent(inout) :: cbmf              
real, dimension(:,:,:),       intent(inout) :: cell_cld_frac,  &
                                               cell_liq_amt,  &
                                               cell_liq_size, &
                                               cell_ice_amt,  &
                                               cell_ice_size, &
                                           cell_droplet_number, &
                                               meso_cld_frac,  &
                                               meso_liq_amt, &
                                               meso_liq_size, &
                                               meso_ice_amt,   &
                                               meso_ice_size, &
                                           meso_droplet_number
integer, dimension(:,:),      intent(inout) :: nsum
real, dimension(:,:),         intent(out)   :: precip, &
                                               lheat_precip, &
                                               vert_motion, &
                                               total_precip
real, dimension(:,:,:),       intent(out)   :: delta_temp, delta_vapor,&
                                               detf, uceml_inter, &
                                               mtot, mfluxup, &
                                               mhalf_3d, &
                                               donner_humidity_area,&
                                               donner_humidity_factor, &
                                               liquid_precip, &
                                               frozen_precip, frz_meso,&
                                            liq_meso, frz_cell, liq_cell
real, dimension(:,:,:,:),     intent(out)   :: qtrtnd 
real, dimension(:,:,:),       intent(out)   :: donner_wetdep
real, dimension(:,:,:),       intent(in),                &
                                   optional :: qlin, qiin, qain
real, dimension(:,:,:),       intent(out),               &
                                   optional :: delta_ql, delta_qi, &
                                               delta_qa

!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!     is, ie         first and last values of i index values of points 
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points 
!                    in this physics window (processor coordinates)
!     dt             physics time step [ sec ]
!     temp           temperature field at model levels [ deg K ]
!     mixing_ratio   vapor mixing ratio field at model levels 
!                    [ kg(h20) / kg(dry air) ]
!     pfull          pressure field on model full levels [ Pa ]
!     phalf          pressure field at half-levels 1:nlev+1  [ Pa ]
!     omega          model omega field at model full levels [ Pa / sec ]
!     qlin           large-scale cloud liquid specific humidity 
!                    [ kg(h2o) / kg (moist air) ]
!     qiin           large-scale cloud ice specific humidity 
!                    [ kg(h2o) / kg (moist air) ]
!     qain           large-scale cloud fraction  
!                    [ fraction ]
!     land           fraction of grid box covered by land
!                    [ fraction ]
!     sfc_sh_flux    sensible heat flux across the surface
!                    [ watts / m**2 ]
!     sfc_vapor_flux water vapor flux across the surface
!                    [ kg(h2o) / (m**2 sec) ]
!     tr_flux        surface flux of tracers transported by
!                    donner_deep_mod [ kg(tracer) / (m**2 sec) ]
!     tracers        tracer mixing ratios
!                    [ kg(tracer) / kg (dry air) ]
!     Time           current time (time_type)
!
!   intent(out) variables:
!
!     precip         precipitation generated by deep convection
!                    [ kg(h2o) / m**2 ]
!     delta_temp     temperature increment due to deep convection 
!                    [ deg K ]
!     delta_vapor    water vapor mixing ratio increment due to deep 
!                    convection [ kg(h2o) / kg (dry air) ]
!     detf           detrained cell mass flux at model levels 
!                    [ (kg / (m**2 sec) ) ]
!     uceml_inter    upward cell mass flux at interface levels 
!                    [ (kg / (m**2 sec) ) ]
!     mtot           mass flux at model full levels, convective plus 
!                    mesoscale, due to donner_deep_mod 
!                    [ (kg / (m**2 sec) ) ]
!     mfluxup        upward mass flux at model full levels, convective 
!                    plus mesoscale, due to donner_deep_mod 
!                    [ (kg / (m**2 sec) ) ]
!     donner_humidity_area
!                    fraction of grid box in which humidity is affected
!                    by the deep convection, defined as 0.0 below cloud
!                    base and above the mesoscale updraft, and as the
!                    sum of the cell and mesoscale cloud areas in 
!                    between. it is used in strat_cloud_mod to determine
!                    the large-scale specific humidity field for the
!                    grid box. DO NOT use for radiation calculation,
!                    since not all of this area includes condensate.
!                    [ fraction ]
!     donner_humidity_ratio
!                    ratio of large-scale specific humidity to specific 
!                    humidity in environment outside convective system
!                    [ dimensionless ]
!     delta_ql       cloud water specific humidity increment due to 
!                    deep convection over the timestep
!                    [ kg (h2o) / kg (moist air) ]
!     delta_qi       cloud ice specific humidity increment due to deep 
!                    convection over the timestep 
!                    [ kg (h2o) / kg (moist air) ]
!     delta_qa       cloud area increment due to deep convection
!                    over the time step [ fraction ]
!     qtrtnd         tracer time tendencies due to deep convection
!                    during the time step
!                    [ kg(tracer) / (kg (dry air) sec) ]
!
!   intent(inout) variables:
!
!     cell_cld_frac  fractional coverage of convective cells in
!                    grid box [ dimensionless ]
!     cell_liq_amt   liquid water content of convective cells
!                    [ kg(h2o) / kg(air) ]
!     cell_liq_size  assumed effective size of cell liquid drops
!                    [ microns ]
!     cell_ice_amt   ice water content of cells
!                    [ kg(h2o) / kg(air) ]
!     cell_ice_size  generalized effective diameter for ice in
!                    convective cells [ microns ]
!     meso_cld_frac  fractional area of mesoscale clouds in grid
!                    box [ dimensionless ]
!     meso_liq_amt   liquid water content in mesoscale clouds
!                    [ kg(h2o) / kg(air) ]
!     meso_liq_size  assumed effective size of mesoscale drops
!                    [ microns ]
!     meso_ice_amt   ice water content of mesoscale elements
!                    [ kg(h2o) / kg(air) ]
!     meso_ice_size  generalized ice effective size for anvil ice
!                    [ microns ]
!     nsum           number of time levels over which the above variables
!                    have so far been summed
!
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!    local variables:

      real,    dimension (size(temp,1), size(temp,2), size(temp,3)) :: &
                       temperature_forcing, moisture_forcing, pmass, &
                       qlin_arg, qiin_arg, qain_arg, delta_ql_arg, & 
                       delta_qi_arg, delta_qa_arg

      real,    dimension (size(temp,1), size(temp,2)) :: parcel_rise, &
                                                         summa

      type(donner_conv_type)            :: Don_conv
      type(donner_budgets_type)         :: Don_budgets
      type(donner_cape_type)            :: Don_cape
      type(donner_rad_type)             :: Don_rad
      type(donner_cem_type)             :: Don_cem
      type(sounding)                    :: sd
      type(adicloud)                    :: ac
      type(cplume)                      :: cp
      type(ctend )                      :: ct
      character(len=128)                :: ermesg
      integer                           :: error
      integer                           :: isize, jsize, nlev_lsm
      integer                           :: ntr, me, root_pe
      logical                           :: cloud_tracers_present
      integer                           :: num_cld_tracers
      integer                           :: k, n   

!--------------------------------------------------------------------
!   local variables:
!
!     temperature_forcing  temperature tendency due to donner convection
!                          [ deg K / sec ]
!     moisture_forcing     vapor mixing ratio tendency due to donner 
!                          convection [ kg(h2o) / (kg(dry air) sec ) ]
!     pmass                mass per unit area within the grid box
!                          [ kg (air) / (m**2) ]
!     parcel_rise          accumulated vertical displacement of a 
!                          near-surface parcel as a result of the lowest
!                          model level omega field [ Pa ]
!     total_precip         total precipitation rate produced by the
!                          donner parameterization [ mm / day ]
!     exit_flag            logical array indicating whether deep conv-
!                          ection exists in a column
!     Don_conv             donner_convection_type derived type variable 
!                          containing diagnostics and intermediate
!                          results describing the nature of the convec-
!                          tion produced by the donner parameterization
!     Don_cape             donner_cape type derived type variable con-
!                          taining diagnostics and intermediate results
!                          related to the cape calculation associated 
!                          with the donner convection parameterization
!     Don_rad              donner_rad_type derived type variable used
!                          to hold those fields needed to connect the
!                          donner deep convection parameterization and
!                          the model radiation package
!     Don_cem              donner_cem_type derived type variable 
!                          containing Donner cumulus ensemble member 
!                          diagnostics
!     ermesg               character string containing any error message
!                          that is returned from a kernel subroutine
!     isize                x-direction size of the current physics window
!     isize, jsize         y-direction size of the current physics window
!     nlev_lsm             number of model layers in large-scale model
!     ntr                  number of tracers to be transported by donner
!                          convection 
!     me                   local pe number
!     calc_conv_on_this_step 
!                          is this a step on which to calculate 
!                          convection ?
!     k                    do-loop index
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    check that the module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized) then
         ermesg = 'donner_deep: &
             &donner_deep_init was not called before subroutine   &
                                                  &donner_deep'
        if (running_in_fms) then
          call fms_error_mesg (ermesg) 
        else

!---------------------------------------------------------------------
!    appropriate error processing code should be added in subroutine
!    nonfms_error_mesg. currently an error message is printed and a 
!    stop command issued (dangerous on parallel machines!).
!---------------------------------------------------------------------
          call nonfms_error_mesg (ermesg) 
        endif
      endif

!----------------------------------------------------------------------
!    determine if the arguments needed when run with the strat_cloud_mod 
!    are present; set cloud_tracers_present appropriately.
!----------------------------------------------------------------------
      num_cld_tracers = count( (/present(qlin), present(qiin),   &
                                 present(qain), present(delta_ql), &
                                 present(delta_qi),present(delta_qa)/) )
      if (num_cld_tracers == 0) then
        cloud_tracers_present = .false.
        qlin_arg = 0.
        qiin_arg = 0.
        qain_arg = 0.
      else if (num_cld_tracers == 6) then
        cloud_tracers_present = .true.
        qlin_arg = qlin 
        qiin_arg = qiin
        qain_arg = qain
      else
        ermesg = 'donner_deep: &
                        &Either none or all of the cloud tracers '// &
                         'and their tendencies must be present'
        if (running_in_fms) then
          call fms_error_mesg (ermesg) 
        else

!---------------------------------------------------------------------
!    appropriate error processing code should be added in subroutine
!    nonfms_error_mesg. currently an error message is printed and a 
!    stop command issued (dangerous on parallel machines!).
!---------------------------------------------------------------------
          call nonfms_error_mesg (ermesg) 
        endif
      endif

!--------------------------------------------------------------------
!    if column diagnostics have been requested for any column, call 
!    donner_column_control to define the components of the 
!    donner_column_diag_type variable for the diagnostic columns in this 
!    window. if column diagnostics have not been requested, the needed
!    variables so indicating have already been set.
!--------------------------------------------------------------------
      if (running_in_fms) then
        if (Col_diag%num_diag_pts > 0) then
          call fms_donner_column_control (is, ie, js, je, secs, days, &
                                          Col_diag)
        endif
      else

!---------------------------------------------------------------------
!    if column diagnostics are desired in the nonFMS model, subroutine
!    nonfms_donner_column_control must be modified to provide the
!    functionality of fms_donner_column_control. by default, column
!    diagnostics are not available with the nonFMS model.
!---------------------------------------------------------------------
        if (Col_diag%num_diag_pts > 0) then
          call nonfms_donner_column_control (is, ie, js, je, secs,  &
                                             days, Col_diag)
        endif
      endif

!-------------------------------------------------------------------
!    define the dimensions for the variables in this physics window.
!    define the pe number of the current pe.
!-------------------------------------------------------------------
      isize     = ie - is + 1
      jsize     = je - js + 1
      nlev_lsm  = size(temp,3)
      ntr       = size(tracers,4) 

      if (running_in_fms) then
        call fms_get_pe_number(me, root_pe)
      else

!---------------------------------------------------------------------
!    subroutine nonfms_get_pe_number should be set up to return the
!    current pe's number, if column daignostics are activated in the
!    nonFMS model. By default, the subroutine returns 0 for all pes.
!---------------------------------------------------------------------
        call nonfms_get_pe_number(me, root_pe)
      endif

      Don_budgets%n_water_budget      = N_WATER_BUDGET
      Don_budgets%n_enthalpy_budget   = N_ENTHALPY_BUDGET
      Don_budgets%n_precip_paths      = N_PRECIP_PATHS     
      Don_budgets%n_precip_types      = N_PRECIP_TYPES     

!-----------------------------------------------------------------------
!    call the kernel subroutine don_d_donner_deep_k to obtain the
!    output fields resulting from the donner deep convection parameter-
!    ization.
!-----------------------------------------------------------------------
      call don_d_donner_deep_k   &
           (is, ie, js, je, isize, jsize, nlev_lsm, NLEV_HIRES, ntr, me,&
            cloud_tracers_present,  cbmf,    &
            dt, Param, Nml, temp, mixing_ratio, pfull,    &
            phalf, zfull, zhalf, omega, pblht, tkemiz, qstar, cush, coldT,&
            qlin_arg, qiin_arg, qain_arg, land, sfc_sh_flux,  &
            sfc_vapor_flux,    &
            tr_flux, tracers, cell_cld_frac, cell_liq_amt,      &
            cell_liq_size, cell_ice_amt, cell_ice_size,   &
            cell_droplet_number, meso_cld_frac,  &
            meso_liq_amt, meso_liq_size, meso_ice_amt, meso_ice_size,  &
            meso_droplet_number, &
            nsum, precip, delta_temp, delta_vapor, detf, uceml_inter,  &
            mtot, mfluxup, donner_humidity_area,  &
            donner_humidity_factor, &
            total_precip, temperature_forcing, moisture_forcing,    &
            parcel_rise, delta_ql_arg, delta_qi_arg, delta_qa_arg,   &
            qtrtnd,         &
            calc_conv_on_this_step, mhalf_3d, ermesg, error, Initialized, Col_diag,   &
            Don_rad, Don_conv, Don_cape, Don_cem, Don_save, &!miz
            sd, Uw_p, ac, cp, ct,  Don_budgets)

!----------------------------------------------------------------------
!    if strat_cloud is active, move the output arguments into the proper
!    locations.
!----------------------------------------------------------------------
      if (cloud_tracers_present) then
        delta_ql = delta_ql_arg
        delta_qi = delta_qi_arg
        delta_qa = delta_qa_arg
      endif

      if (Initialized%do_conservation_checks .or.   &
                                          Nml%do_budget_analysis) then
        lheat_precip = Don_budgets%lheat_precip
        vert_motion = Don_budgets%vert_motion
      else
        lheat_precip = 0.
        vert_motion = 0.
      endif
      liquid_precip = Don_budgets%liq_prcp
      frozen_precip = Don_budgets%frz_prcp

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, process the error message.
!!!  HOW TO DISTINGUISH FATAL, WARNING, NOTE ??
!    FOR NOW, ALL messages considered FATAL.
!----------------------------------------------------------------------
      if (error /= 0) then
        if (running_in_fms) then
          call fms_error_mesg (ermesg) 
        else

!---------------------------------------------------------------------
!    appropriate error processing code should be added in subroutine
!    nonfms_error_mesg. currently an error message is printed and a 
!    stop command issued (dangerous on parallel machines!).
!---------------------------------------------------------------------
          call nonfms_error_mesg (ermesg) 
        endif
      endif

!---------------------------------------------------------------------
!    if this is a calculation step for donner_deep, define a mass
!    weighting factor (mass per unit area) needed for some of the netcdf
!    diagnostics (pmass). call donner_deep_netcdf to send the requested 
!    diagnostic data to the diag_manager for output.
!---------------------------------------------------------------------
      if (calc_conv_on_this_step) then
        do k=1,nlev_lsm
          pmass(:,:,k) = (phalf(:,:,k+1) - phalf(:,:,k))/Param%GRAV   
        end do

!---------------------------------------------------------------------
!   define the column integrated tracer wet deposition associated with 
!   donner convection so that it may be returned to moist_processes.
!---------------------------------------------------------------------

        donner_wetdep = 0.
        do n=1, ntr
          summa = 0.
          do k=1,nlev_lsm
            summa(:,:) = summa(:,:) + &
                                 Don_conv%wetdept(:,:,k,n)*pmass(:,:,k)
          end do
          donner_wetdep(:,:,n) = summa(:,:)
        end do
        if (running_in_fms) then
          call fms_donner_deep_netcdf (is, ie, js, je, Nml, secs, days,&
                                 Param, Initialized, Don_conv,  &
                                 Don_cape, Don_cem,parcel_rise, pmass, &
                                 total_precip,  Don_budgets, &
                                 temperature_forcing, &
                                 moisture_forcing)
        else

!---------------------------------------------------------------------
!    subroutine nonfms_donner_deep_netcdf should be set up to output 
!    the netcdf diagnostic fields. By default, it does nothing.
!---------------------------------------------------------------------
          call nonfms_donner_deep_netcdf 
        endif

!----------------------------------------------------------------------
!    on calculation steps, update the values of the cell and
!    mesoscale cloud variables to be returned to moist_processes_mod. 
!    (on non-calculation steps, the values that were passed in are 
!    simply passed back.)
!----------------------------------------------------------------------
        cell_cld_frac = Don_rad%cell_cloud_frac
        cell_liq_amt  = Don_rad%cell_liquid_amt
        cell_liq_size = Don_rad%cell_liquid_size
        cell_ice_amt  = Don_rad%cell_ice_amt
        cell_ice_size = Don_rad%cell_ice_size
        cell_droplet_number = Don_rad%cell_droplet_number
        meso_cld_frac = Don_rad%meso_cloud_frac
        meso_liq_amt  = Don_rad%meso_liquid_amt
        meso_liq_size = Don_rad%meso_liquid_size
        meso_ice_amt  = Don_rad%meso_ice_amt
        meso_ice_size = Don_rad%meso_ice_size
        meso_droplet_number = Don_rad%meso_droplet_number
        nsum          = Don_rad%nsum

!--------------------------------------------------------------------
!    define the precip fields at each model level for liq and frozen 
!    precip associated with the cell and meso circulations.
!    UNITS KG / KG/ DAY 
!--------------------------------------------------------------------
        if (Initialized%do_conservation_checks .or.   &
                                          Nml%do_budget_analysis) then
          do k=1,nlev_lsm
            frz_meso(:,:,k) = (Don_budgets%precip_budget(:,:,k,2,2) + &
                   Don_budgets%precip_budget(:,:,k,2,3) + &
                   Don_budgets%precip_budget(:,:,k,4,2) + &
                   Don_budgets%precip_budget(:,:,k,4,3))*  &
                                                     Don_conv%a1(:,:)
            liq_meso(:,:,k) = (Don_budgets%precip_budget(:,:,k,1,2) + &
                   Don_budgets%precip_budget(:,:,k,1,3) + &
                   Don_budgets%precip_budget(:,:,k,3,2) + &
                   Don_budgets%precip_budget(:,:,k,3,3) + &
                   Don_budgets%precip_budget(:,:,k,5,2) + &
                   Don_budgets%precip_budget(:,:,k,5,3))*  &
                                                     Don_conv%a1(:,:)
            frz_cell(:,:,k) = (Don_budgets%precip_budget(:,:,k,2,1) + &
                   Don_budgets%precip_budget(:,:,k,4,1))*   &
                                                     Don_conv%a1(:,:)
            liq_cell(:,:,k) = (Don_budgets%precip_budget(:,:,k,1,1) + &
                   Don_budgets%precip_budget(:,:,k,3,1) + &
                   Don_budgets%precip_budget(:,:,k,5,1))*  &
                                                     Don_conv%a1(:,:)
          end do
        endif

!--------------------------------------------------------------------
!    call deallocate_local_variables to deallocate space used by the
!    local derived-type variables.
!--------------------------------------------------------------------
        call don_d_dealloc_loc_vars_k   &
               (Don_conv, Don_cape, Don_rad, Don_cem,Don_budgets, Nml, &
                          Initialized, sd, ac, cp, ct, ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, process the error message.
!----------------------------------------------------------------------
        if (error /= 0) then
          if (running_in_fms) then
            call fms_error_mesg (ermesg) 
          else

!---------------------------------------------------------------------
!    appropriate error processing code should be added in subroutine
!    nonfms_error_mesg. currently an error message is printed and a 
!    stop command issued (dangerous on parallel machines!).
!---------------------------------------------------------------------
            call nonfms_error_mesg (ermesg) 
          endif
        endif
      endif  ! (calc_conv_on_this_step)

!--------------------------------------------------------------------


end subroutine donner_deep


!#######################################################################
! <SUBROUTINE NAME="donner_deep_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine donner_deep_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp
  integer                                :: ntracers

  if (running_in_fms) then
     call fms_donner_write_restart (Initialized, timestamp)
  else 

     !---------------------------------------------------------------------
     !    subroutine nonfms_donner_write_restart should be configured to
     !    write a netcdf restart file in the nonFMS framework (see subroutine
     !    fms_donner_write_restart for the variables which must be included).
     !    by default, the subroutine does nothing.
     !---------------------------------------------------------------------
     if(present(timestamp)) then
        call nonfms_error_mesg('donner_deep_mod: when running_in_fms is false, '// &
             'timestamp should not passed in donner_deep_restart')
     endif
     ntracers = size(Don_save%tracername(:))
     call nonfms_donner_write_restart (ntracers, Don_save, &
          Initialized, Nml)
  endif

end subroutine donner_deep_restart
! </SUBROUTINE> NAME="donner_deep_restart"


!####################################################################

subroutine donner_deep_end

!---------------------------------------------------------------------
!   donner_deep_end is the destructor for donner_deep_mod.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variable

      integer  :: ntracers     ! number of tracers transported by the
                               ! donner deep convection parameterization

!-------------------------------------------------------------------
!    if module has not been initialized, return.
!-------------------------------------------------------------------
      if (.not. module_is_initialized) return

!-------------------------------------------------------------------
!    define the number of tracers that have been transported by the 
!    donner deep convection parameterization.
!-------------------------------------------------------------------
      ntracers = size(Don_save%tracername(:))

!-------------------------------------------------------------------
!    call subroutine to write restart file. NOTE: only the netcdf 
!    restart file is currently supported.
!-------------------------------------------------------------------
      if (running_in_fms) then
        call donner_deep_restart
      else 

!---------------------------------------------------------------------
!    subroutine nonfms_donner_write_restart should be configured to
!    write a netcdf restart file in the nonFMS framework (see subroutine
!    fms_donner_write_restart for the variables which must be included).
!    by default, the subroutine does nothing.
!---------------------------------------------------------------------
        call nonfms_donner_write_restart (ntracers, Don_save, &
                                       Initialized, Nml)
      endif 


!-------------------------------------------------------------------
!    close any column diagnostics units which are open.
!------------------------------------------------------------------
      if (Col_diag%num_diag_pts > 0) then
        if (running_in_fms) then
          call fms_close_col_diag_units 
        else

!--------------------------------------------------------------------
!    subroutine nonfms_close_column_diagnostics_units should be used
!    to close the outpuit units activated for column duiagnostics.
!    by default, column diagnostics are not available in the nonFMS
!    model, and this subroutine does nothing.
!--------------------------------------------------------------------
          call nonfms_close_col_diag_units 
        endif
      endif

!----------------------------------------------------------------------
!    call deallocate_variables to deallocate the module variables.
!----------------------------------------------------------------------
      call deallocate_variables 

      if (running_in_fms) then
        call fms_deallocate_variables (Col_diag)
      else

!---------------------------------------------------------------------
!    subroutine nonfms_deallocate_variables should be used to deallocate
!    local arrays associated with the column diagnostics option,
!    the variable monitoring option, and netcdf output. Currently the 
!    subroutine does nothing, since these features are not available 
!    in the nonFMS model.
!---------------------------------------------------------------------
        call nonfms_deallocate_variables
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.


!---------------------------------------------------------------------

end subroutine donner_deep_end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                   PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      1. ROUTINES CALLED BY DONNER_DEEP_INIT
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
 



!####################################################################


!#####################################################################








!#####################################################################


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      2. ROUTINES CALLED BY DONNER_DEEP
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!######################################################################



!######################################################################


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      3. ROUTINES CALLED BY DONNER_DEEP_END
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





!######################################################################

subroutine deallocate_variables 

!---------------------------------------------------------------------
!    subroutine deallocate_variables deallocates the space used by the
!    module variables.
!---------------------------------------------------------------------
 
      deallocate ( Don_save%cemetf              )
      deallocate ( Don_save%lag_temp            )
      deallocate ( Don_save%lag_vapor           )
      deallocate ( Don_save%lag_press           )
      deallocate ( Don_save%cememf              )
      deallocate ( Don_save%mass_flux           )
      deallocate ( Don_save%mflux_up            )
      deallocate ( Don_save%cell_up_mass_flux   )
      deallocate ( Don_save%det_mass_flux       )
      deallocate ( Don_save%dql_strat           )
      deallocate ( Don_save%dqi_strat           )
      deallocate ( Don_save%dqa_strat           )
      deallocate ( Don_save%humidity_area       )
      deallocate ( Don_save%humidity_factor     )
      deallocate ( Don_save%tracer_tends        )
      deallocate ( Don_save%parcel_disp         )
      deallocate ( Don_save%tprea1              )
      deallocate ( Don_save%tracername          )
      deallocate ( Don_save%tracer_units        )

      deallocate (Param%arat)
      deallocate (Param%ensemble_entrain_factors_gate)
      deallocate (Param%ensemble_entrain_factors_kep )
      deallocate (Param%dgeice)
      deallocate (Param%relht )

      call exn_end_k
      call findt_end_k


!----------------------------------------------------------------------


end subroutine deallocate_variables 




!######################################################################



                     end module donner_deep_mod



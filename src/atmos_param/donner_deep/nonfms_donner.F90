                       module nonfms_donner_mod

use  sat_vapor_pres_k_mod, only: sat_vapor_pres_init_k ! replace with 
                                                     ! non-FMS interface
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

implicit none
private

!--------------------------------------------------------------------
!        donner_deep_mod diagnoses the location and computes the 
!        effects of deep convection on the model atmosphere
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------


character(len=128)  :: version =  '$Id: nonfms_donner.F90,v 19.0 2012/01/06 20:08:40 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!--------------------------------------------------------------------
!---interfaces------

public   &
        nonfms_donner_process_nml,  nonfms_donner_process_tracers, &
        nonfms_donner_process_monitors, &
        nonfms_donner_activate_diag, nonfms_donner_read_restart,&
        nonfms_donner_col_diag, nonfms_donner_write_restart, &
        nonfms_donner_column_control, nonfms_donner_deep_netcdf,    &
        nonfms_sat_vapor_pres, nonfms_get_pe_number, nonfms_error_mesg,&
        nonfms_close_col_diag_units, &
        nonfms_deallocate_variables, nonfms_constants

private   &
!  module subroutines called during initialization:
        process_coldstart


!---------------------------------------------------------------------
!---namelist----

# include "donner_nml.h"


!--------------------------------------------------------------------
!--- public data ----------




!--------------------------------------------------------------------
!----private data-----------








!-----------------------------------------------------------------------
!   miscellaneous variables
!


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                   PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################

subroutine nonfms_donner_process_nml  (Nml, kpar)

!---------------------------------------------------------------------
!    nonfms_donner_process_nml si intended to process the 
!    donner_deep_nml file, using the procedure of the nonFMS model.
!    for now, default values are reset here within the Fortran source
!    and no nml read is done.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
type(donner_nml_type), intent(inout)    :: Nml
integer,               intent(in)       :: kpar

!---------------------------------------------------------------------
!  intent(in) variables:
!
!
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables:

     real, dimension(:), allocatable :: erat_loc, arat_loc
  
!-------------------------------------------------------------------
!  local variables:
!
!                         
!-------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    1. READ NAMELIST AND WRITE IT TO LOG FILE.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!--------------------------------------------------------------------
!    here non-default values are reset (as desired) within the Fortran 
!    source. note that changes to arat / erat are handled at the end of
!    this subroutine.
!--------------------------------------------------------------------
!  THESE SETTINGS MAY BE USED FOR DONNER_LITE:
!    SETTINGS USED IN DATABASE EXPT C48L24_AM3p5-gamma-B6:
       parcel_launch_level = 2
       model_levels_in_sfcbl = 0
       donner_deep_freq = 1800
       allow_mesoscale_circulation = .true.
       do_donner_cape    = .false.
       do_donner_plume   = .false.
       do_donner_closure = .false.
       do_donner_lscloud = .true.
       do_dcape          = .false.
       do_lands          = .false.
       do_freezing_for_cape = .true.
       do_freezing_for_closure = .true.
       gama              = 0.0
       tau               = 28800.
       tke0              = 0.5
       cape0             = 1000.
       lochoice          = 10
       do_capetau_land   = .false.
       use_llift_criteria= .false.
       do_ice            = .true.
       atopevap  = 0.1
       auto_rate = 1.e-3
       auto_th   = 0.5e-3
       frac      = 1.65
       ttend_max = 0.005

       EVAP_IN_DOWNDRAFTS  = 0.00
       EVAP_IN_ENVIRON     = 0.00
       ENTRAINED_INTO_MESO = 1.00

       ANVIL_PRECIP_EFFICIENCY = 0.85
       MESO_DOWN_EVAP_FRACTION = 0.1
       MESO_UP_EVAP_FRACTION   = 0.05

       wmin_ratio      = 0.05
       arat(1:7) =  (/ 1.0, 0.26, 0.35, 0.32, 0.3, 0.54, 0.66 /)
       erat(1:7) =  (/ 1.0, 1.30, 1.80, 2.50, 3.3, 4.50, 10.0 /)
       frc_internal_enthalpy_conserv = .true.
       limit_pztm_to_tropo = .true.
 
!  THESE SETTINGS MAY BE USED FOR DONNER_FULL:
!      parcel_launch_level = 2
!      donner_deep_freq = 1800
!      allow_mesoscale_circulation = .true.
!      do_donner_cape    = .true.
!      do_donner_plume   = .true.
!      do_donner_closure = .true.
!      do_donner_lscloud = .true.
!      do_dcape          = .true.
!      do_freezing_for_cape = .false.
!      do_freezing_for_closure = .false.
!      gama              = 0.0
!      lochoice          = 10
!      use_llift_criteria= .true.
!      EVAP_IN_DOWNDRAFTS  = 0.25
!      EVAP_IN_ENVIRON     = 0.13
!      ENTRAINED_INTO_MESO = 0.62

!      ANVIL_PRECIP_EFFICIENCY = 0.5
!      MESO_DOWN_EVAP_FRACTION = 0.4
!      MESO_UP_EVAP_FRACTION   = 0.1 

!      wmin_ratio      = 0.05
!      arat(1:7) =  (/ 1.0, 0.26, 0.35, 0.32, 0.3, 0.54, 0.3 /)
!      erat(1:7) =  (/ 1.0, 1.30, 1.80, 2.50, 3.3, 4.50, 6.5 /)
!      frc_internal_enthalpy_conserv = .true.
!      limit_pztm_to_tropo = .true.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    8. STORE THE NAMELIST VARIABLES THAT NEED TO BE MADE AVAILABLE 
!       OUTSIDE OF THIS MODULE INTO THE DONNER_NML_TYPE VARIABLE.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      Nml%parcel_launch_level         = parcel_launch_level
      Nml%allow_mesoscale_circulation = allow_mesoscale_circulation
      Nml%do_hires_cape_for_closure   = do_hires_cape_for_closure
      Nml%do_donner_cape              = do_donner_cape    !miz
      Nml%do_donner_plume             = do_donner_plume   !miz
      Nml%do_donner_closure           = do_donner_closure !miz
      Nml%do_dcape                    = do_dcape          !miz
      Nml%do_lands                    = do_lands          !miz
      Nml%tau                         = tau               !miz
      Nml%cape0                       = cape0             !miz
      Nml%rhavg0                      = rhavg0            !miz
      Nml%plev0                       = plev0             !miz
      Nml%do_rh_trig                  = do_rh_trig        !miz
      Nml%do_capetau_land             = do_capetau_land   !miz
      Nml%pblht0                      = pblht0            !miz
      Nml%tke0                        = tke0              !miz
      Nml%lofactor0                   = lofactor0         !miz
      Nml%deephgt0                    = deephgt0          !miz
      Nml%lochoice                    = lochoice          !miz
      Nml%deep_closure                = deep_closure      !miz
      Nml%gama                        = gama              !miz
      Nml%do_ice                      = do_ice            !miz
      Nml%atopevap                    = atopevap          !miz
      Nml%do_donner_lscloud           = do_donner_lscloud !miz
      Nml%auto_rate                   = auto_rate         !miz
      Nml%auto_th                     = auto_th           !miz
      Nml%frac                        = frac              !miz
      Nml%ttend_max                   = ttend_max         !miz
      Nml%mesofactor                  = mesofactor        !miz
      Nml%use_llift_criteria          = use_llift_criteria
      Nml%use_pdeep_cv                = use_pdeep_cv
      Nml%entrainment_constant_source = entrainment_constant_source
      Nml%donner_deep_freq            = donner_deep_freq             
      Nml%model_levels_in_sfcbl       = model_levels_in_sfcbl        
      Nml%cell_liquid_size_type       = cell_liquid_size_type 
      Nml%cell_ice_size_type          = cell_ice_size_type
      Nml%cell_liquid_eff_diam_input  = cell_liquid_eff_diam_input
      Nml%cell_ice_geneff_diam_input  = cell_ice_geneff_diam_input
      Nml%meso_liquid_eff_diam_input  = meso_liquid_eff_diam_input
      Nml%do_average                  = do_average
      Nml%use_memphis_size_limits     = use_memphis_size_limits
      Nml%wmin_ratio                  = wmin_ratio
      Nml%do_freezing_for_cape         = do_freezing_for_cape
      Nml%tfre_for_cape               = tfre_for_cape
      Nml%dfre_for_cape               = dfre_for_cape
      Nml%rmuz_for_cape               = rmuz_for_cape
      Nml%do_freezing_for_closure     = do_freezing_for_closure
      Nml%tfre_for_closure            = tfre_for_closure
      Nml%dfre_for_closure            = dfre_for_closure
      Nml%rmuz_for_closure            = rmuz_for_closure
      Nml%do_budget_analysis          = do_budget_analysis
      Nml%frc_internal_enthalpy_conserv =  &
                                 frc_internal_enthalpy_conserv
      Nml%do_ensemble_diagnostics     = do_ensemble_diagnostics
      Nml%limit_pztm_to_tropo = limit_pztm_to_tropo
      Nml%entrainment_scheme_for_closure =   &
                                        entrainment_scheme_for_closure
      Nml%modify_closure_plume_condensate =   &
                                       modify_closure_plume_condensate
      Nml%closure_plume_condensate = closure_plume_condensate

      Nml%evap_in_downdrafts = evap_in_downdrafts
      Nml%evap_in_environ  = evap_in_environ
      Nml%entrained_into_meso = entrained_into_meso
      Nml%anvil_precip_efficiency = anvil_precip_efficiency
      Nml%meso_down_evap_fraction = meso_down_evap_fraction
      Nml%meso_up_evap_fraction = meso_up_evap_fraction
      Nml%cdeep_cv = cdeep_cv

!---------------------------------------------------------------------
!  if mods are desired for arat / erat when these values are being
!  specified (option = 0), make them to arat_loc / erat_loc. these
!  will be transferred to arat / erat later. if arat / erat come from
!  optional formulae, they will be calculated here used nml-supplied
!  input values.
!---------------------------------------------------------------------
      allocate (arat_loc(kpar))
      allocate (erat_loc(kpar))
 
      if (arat_erat_option == 0) then
        arat_loc = arat
        erat_loc = erat
      else
        call define_arat_erat (arat_erat_option, kpar, eratb, erat0, &
                               erat_min, erat_max,erat_loc,arat_loc)
        print *,'donner_deep_nml: redefined arat and erat using &
                          &arat_erat_option == ', arat_erat_option
        print *,'donner_deep_nml: arat = ',arat_loc
        print *,'donner_deep_nml: erat = ',erat_loc
      endif
      allocate (Nml%arat(kpar))
      allocate (Nml%ensemble_entrain_factors_gate(kpar)) 

      Nml%arat = arat_loc
      Nml%ensemble_entrain_factors_gate = erat_loc

      deallocate (arat_loc, erat_loc)

end subroutine nonfms_donner_process_nml



!#####################################################################

subroutine nonfms_donner_process_tracers 


      return



end subroutine nonfms_donner_process_tracers



!#####################################################################

subroutine nonfms_donner_process_monitors


      return



end subroutine nonfms_donner_process_monitors


!#####################################################################

subroutine nonfms_donner_activate_diag (secs, days, axes, &
                     Don_save, Nml, n_water_budget, n_enthalpy_budget, &
                   n_precip_paths, n_precip_types, nlev_hires, kpar)

integer, intent(in) :: secs, days, n_water_budget,  &
                   n_enthalpy_budget, n_precip_paths, n_precip_types, &
                   nlev_hires, kpar
integer,         dimension(4),   intent(in)   :: axes
type(donner_save_type), intent(inout) :: Don_save
type(donner_nml_type), intent(inout) :: Nml       


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    4. INITIALIZE THE NETCDF OUTPUT VARIABLES.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!--------------------------------------------------------------------
!    activate the netcdf diagnostic fields.
!-------------------------------------------------------------------
!     call register_fields (secs, days, axes, Don_save, Nml)


end subroutine nonfms_donner_activate_diag 


!#####################################################################

subroutine nonfms_donner_read_restart (Initialized, ntracers,   &
                                    secs, days, Don_save, Nml)

type(donner_initialized_type), intent(inout) :: Initialized
type(donner_save_type), intent(inout) :: Don_save
type(donner_nml_type), intent(inout) :: Nml     
integer, intent(in) :: secs, days, ntracers



!--------------------------------------------------------------------
!    if no restart file is present, call subroutine process_coldstart
!    to define the needed variables.
!--------------------------------------------------------------------
      call process_coldstart (secs, days, Initialized, Nml, Don_save)


end subroutine nonfms_donner_read_restart 


!#####################################################################

subroutine nonfms_donner_col_diag (lonb, latb, Col_diag, pref) 

real, dimension(:,:), intent(in) :: lonb, latb
type(donner_column_diag_type), intent(inout) :: Col_diag
real, dimension(:), intent(in) :: pref


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    6. INITIALIZE VARIABLES NEEDED FOR COLUMN_DIAGNOSTICS_MOD OUTPUT.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    define the total number of columns for which diagnostics
!    are desired.
!---------------------------------------------------------------------
      Col_diag%num_diag_pts = num_diag_pts_ij + num_diag_pts_latlon

!---------------------------------------------------------------------
!    initialize the value of the k index associated with diagnostics
!    cutoff.
!---------------------------------------------------------------------
      Col_diag%kstart = -99



end subroutine nonfms_donner_col_diag 




!#####################################################################

subroutine nonfms_donner_write_restart (ntracers, Don_save,  &
                                        Initialized, Nml)

integer, intent(in) :: ntracers
type(donner_initialized_type), intent(inout) :: Initialized
type(donner_save_type), intent(inout) :: Don_save
type(donner_nml_type), intent(inout) :: Nml     

!-------------------------------------------------------------------
!    call subroutine to write restart file. 
!-------------------------------------------------------------------

      return


end subroutine nonfms_donner_write_restart 


!#####################################################################

subroutine nonfms_get_pe_number (me, root_pe)
 
!--------------------------------------------------------------------
!    define pe number (needed for column diagnostics and as dummy arg-
!    ument for donner_lite diagnostics). For now, column
!    diagnostics are unavailable outside of FMS, so the value is set 
!    to 0 for all pes. 
!--------------------------------------------------------------------

integer, intent(out) :: me, root_pe
   
      me = 0
      root_pe = 0

end subroutine nonfms_get_pe_number




!#####################################################################

subroutine nonfms_error_mesg (ermesg)   
                             
character(len=*), intent(in) :: ermesg     

 
!    call error_mesg ('donner_deep_mod', ermesg, FATAL)
!!  NOTE POTENTIAL HANG HERE : USE APPROPRIATE ERROR EXIT ON NONFMS
!!  SYSTEM RATHER THAN 'STOP'
    print *, 'STOPPING DUE TO ERROR:', ermesg
    stop 

 

 
end subroutine nonfms_error_mesg




!#####################################################################

subroutine nonfms_close_col_diag_units 

      return

end subroutine nonfms_close_col_diag_units 



!#####################################################################

subroutine nonfms_deallocate_variables

      return 




end subroutine nonfms_deallocate_variables

!######################################################################

subroutine nonfms_sat_vapor_pres
 
!---------------------------------------------------------------------
!    should contain needed calls to initialize nonfms saturation
!    vapor pressure calculation. currently uses fms interface to allow
!    testing.
!---------------------------------------------------------------------

integer, parameter :: TCMIN = -160
integer, parameter :: TCMAX = 100
integer, parameter :: ESRES = 10
real,    parameter :: HLV = 2.500e6   
real,    parameter :: ES0 = 1.0 
real,    parameter :: RVGAS = 461.50 
integer, parameter :: NSIZE = (TCMAX-TCMIN)*esres + 1
integer, parameter :: NLIM = NSIZE - 1
real, parameter :: TFREEZE = 273.16
logical, parameter :: use_exact_qs_input = .true.
logical, parameter :: do_simple = .false.
logical, parameter :: construct_table_wrt_liq = .false.
logical, parameter :: construct_table_wrt_liq_and_ice = .false.

      real  :: teps, tmin, dtinv
      character(len=128) :: err_msg
!     logical :: dum = .false.

      call sat_vapor_pres_init_k (NSIZE, REAL(TCMIN), REAL(TCMAX), &
                             TFREEZE, HLV, RVGAS, ES0, err_msg,  &
                             use_exact_qs_input, do_simple,  &
                             construct_table_wrt_liq, &
                             construct_table_wrt_liq_and_ice, &
                             teps, tmin, dtinv)
 

end subroutine nonfms_sat_vapor_pres




!######################################################################

subroutine nonfms_constants (Param)

type(donner_param_type), intent(inout)  :: Param


!----------------------------------------------------------------------
!    define the components of Param that come from the fms module
!    constants_mod. see donner_types.h for their definitions.
!----------------------------------------------------------------------
      Param%dens_h2o        = 1000.    
      Param%rdgas           = 287.04
      Param%kappa           = 2. / 7.
      Param%grav            = 9.80 
      Param%cp_air          = Param%rdgas/ Param%kappa
      Param%pie             = 3.14159265358979323846
      Param%rvgas           = 461.5
      Param%seconds_per_day = 86400.          
      Param%hlv             = 2.500e+06
      Param%hlf             = 3.34e+05
      Param%hls             = Param%hlv + Param%hlf
      Param%kelvin          = 273.15 

!----------------------------------------------------------------------


end subroutine nonfms_constants


!######################################################################


 subroutine nonfms_donner_column_control (is, ie, js, je, secs,  &
                                          days, Col_diag)               

!---------------------------------------------------------------------
!    subroutine fms_donner_column_control returns the number, location
!    (processor and window indices) and output units associated with 
!    any diagnostic columns requested within the current physics window.
!---------------------------------------------------------------------

integer,                       intent(in)   :: is, ie, js, je
integer,                       intent(in) :: secs, days
type(donner_column_diag_type), intent(inout) :: Col_diag
      
      return     
 
end subroutine nonfms_donner_column_control



subroutine nonfms_donner_deep_netcdf 



      return



end subroutine nonfms_donner_deep_netcdf 


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

subroutine process_coldstart (secs, days, Initialized, Nml, Don_save)

!-----------------------------------------------------------------------
!    subroutine process_coldstart provides initialization that is needed
!    when the job is a donner_deep coldstart, or if the user-supplied 
!    restart file is not usable for a restart with the current code 
!    version.
!-----------------------------------------------------------------------

integer, intent(in) :: secs, days
type(donner_initialized_type), intent(inout) :: Initialized
type(donner_save_type), intent(inout) :: Don_save
type(donner_nml_type), intent(inout) :: Nml     

!---------------------------------------------------------------------
!   intent(in) variables:
!
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:


!---------------------------------------------------------------------
!    set the coldstart flag to .true.. set the time until the first cal-
!    culation call to donner_deep_mod, donner_deep calculation calls will
!    be every donner_deep_freq seconds after the start of the day.
!---------------------------------------------------------------------
      Initialized%coldstart = .true.
      if (secs == 0) then    ! i.e., 00Z
        Initialized%conv_alarm = Nml%donner_deep_freq
      else 
        Initialized%conv_alarm = Nml%donner_deep_freq -   &
                                 MOD (secs, Nml%donner_deep_freq)
      endif

!----------------------------------------------------------------------
!    initialize the variables which must be returned from donner_deep_mod
!    on the first step when coldstarting.
!----------------------------------------------------------------------
      Don_save%cemetf            = 0.
      Don_save%cememf            = 0.
      Don_save%tracer_tends      = 0.
      Don_save%mass_flux         = 0.
      Don_save%mflux_up          = 0.
      Don_save%cell_up_mass_flux = 0.
      Don_save%det_mass_flux     = 0.
      Don_save%dql_strat         = 0.
      Don_save%dqi_strat         = 0.
      Don_save%dqa_strat         = 0.
      Don_save%humidity_area     = 0.
      Don_save%humidity_factor   = 0.
      Don_save%tprea1            = 0.
      Don_save%parcel_disp       = 0.

!----------------------------------------------------------------------


end subroutine process_coldstart






!#####################################################################





                     end module nonfms_donner_mod


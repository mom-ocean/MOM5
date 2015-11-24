             module fms_donner_mod

use time_manager_mod,       only: time_type, set_time, &
                                  set_date, get_time,   &
                                  get_calendar_type, &
                                  operator(-), &
                                  operator(>=), operator (<)
use diag_manager_mod,       only: register_diag_field, send_data, &
                                  diag_axis_init
use field_manager_mod,      only: MODEL_ATMOS, field_manager_init, &
                                  fm_query_method, get_field_info, &
                                  parse
use tracer_manager_mod,     only: get_tracer_names,get_number_tracers, &
                                  get_tracer_indices, &
!++lwh
                                  query_method
use atmos_tracer_utilities_mod, only : get_wetdep_param
use  sat_vapor_pres_mod,only : sat_vapor_pres_init
!--lwh
use fms_mod,                only: mpp_pe, mpp_root_pe,  &
                                  file_exist,  check_nml_error,  &
                                  error_mesg, FATAL, WARNING, NOTE,  &
                                  close_file, open_namelist_file,    &
                                  stdout, stdlog, write_version_number,  &
                                  field_size, &
                                  read_data, write_data, lowercase
use fms_io_mod,             only: register_restart_field, restart_file_type, &
                                  save_restart, restore_state, get_mosaic_tile_file
use mpp_mod,                only: input_nml_file
use mpp_io_mod,             only: mpp_open, mpp_close, fieldtype,  &
                                  mpp_read_meta, mpp_get_info, &
                                  mpp_get_fields, mpp_read, &
                                  MPP_NETCDF, MPP_SINGLE,   &
                                  MPP_SEQUENTIAL, MPP_RDONLY, MPP_NATIVE, &
                                  mpp_get_field_name
use constants_mod,          only: DENS_H2O, RDGAS, GRAV, CP_AIR,  &
                                  pie=>PI, KAPPA, RVGAS, &
                                  SECONDS_PER_DAY, HLV, HLF, HLS, KELVIN
use column_diagnostics_mod, only: initialize_diagnostic_columns, &
                                  column_diagnostics_header, &
                                  close_column_diagnostics_units
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


character(len=128)  :: version =  '$Id: fms_donner.F90,v 20.0 2013/12/13 23:17:30 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!--------------------------------------------------------------------
!---interfaces------

public   &
        fms_donner_process_nml,                             &
        fms_donner_process_tracers, &
        fms_donner_activate_diagnostics, fms_donner_read_restart, &
        fms_donner_col_diag, fms_donner_write_restart, &
        fms_donner_column_control, &
        fms_sat_vapor_pres, &
        fms_get_pe_number,  fms_error_mesg,  fms_constants, & 
        fms_close_col_diag_units, &
        fms_deallocate_variables, &
        fms_donner_deep_netcdf, fms_donner_process_monitors

private   &
!  module subroutines called by donner_deep_init:
        register_fields, read_restart_nc,  &
        process_coldstart,&
!  module subroutines called by donner_deep:
        donner_deep_netcdf, donner_column_control,     &
!  module subroutines called from donner_deep_end:
        write_restart


!---------------------------------------------------------------------
!---namelist----

# include "donner_nml.h"

!--------------------------------------------------------------------
!--- public data ----------




!--------------------------------------------------------------------
!----private data-----------

!--- for restart file
type(restart_file_type), pointer, save :: Don_restart => NULL()
type(restart_file_type), pointer, save :: Til_restart => NULL()
logical                                :: in_different_file = .false.
!---------------------------------------------------------------------
!  parameters stored in the donner_param derived type variable to facili-
!  tate passage to kernel subroutines:
!


!--------------------------------------------------------------------
!   list of native mode restart versions usable by this module:
!
!   NOTE: none of the earlier versions of restart files can be used to
!         initiate an experiment with this code version due to a change 
!         in the calculation algorithm. experiments begun with this code
!         must be coldstarted, or use a native mode restart file gener-
!         ated by an experiment using this code version (restart version
!         #8), or a netcdf restart file.
!          
!   version 8 has the lag temp, vapor and pressure fields needed to cal-
!             culate the lag time value of cape. tempbl and ratpbl
!             removed. 
!
!   version 9 is reserved for the native mode restart file version cor-
!             responding to the current netcdf restart file. it is up to 
!             the user to generate the code needed to read and write this
!             version, if needed, using the subroutines read_restart and 
!             write_restart that are provided as starting points, since 
!             only netcdf restarts are currently supported.
!
!   version 10 contains donner_humidity_factor rather than 
!             donner_humidity_ratio, a change necessitated by the intro-
!             duction of the uw_conv shallow convection scheme.

integer, dimension(3)  :: restart_versions = (/ 8, 9, 10 /)


!--------------------------------------------------------------------
!   variables associated with netcdf diagnostic output from this module:
!
!   id_xxxx         indices associated with each potential netcdf 
!                   diagnostic field:
!   missing value   value used by netcdf routines if data not present
!   mod_name        module name associated with these diagnostics; used
!                   to connect these diagnostics to the diag_table
!

integer    :: id_leff
integer    :: id_cemetf_deep, id_ceefc_deep, id_cecon_deep, &
              id_cemfc_deep, id_cememf_deep, id_cememf_mod_deep, &
              id_cual_deep, id_fre_deep, id_elt_deep, &
              id_cmus_deep, id_ecds_deep, id_eces_deep, &
              id_emds_deep, id_emes_deep, id_qmes_deep,&
              id_wmps_deep, id_wmms_deep, id_tmes_deep,&
              id_dmeml_deep, id_uceml_deep, id_umeml_deep, &
              id_xice_deep,  id_dgeice_deep, id_dgeliq_deep,  &
              id_xliq_deep,    &
              id_cuqi_deep, id_cuql_deep, &
              id_plcl_deep, id_plfc_deep, id_plzb_deep, &
              id_xcape_deep, id_coin_deep,  &
              id_dcape_deep, id_qint_deep, id_a1_deep, &
              id_amax_deep, id_amos_deep, &
              id_tprea1_deep, id_ampta1_deep, &
              id_omint_deep, id_rcoa1_deep, id_detmfl_deep
integer                  :: id_pfull_cem, id_phalf_cem, &
                            id_zfull_cem, id_zhalf_cem, &
                            id_temp_cem, id_mixing_ratio_cem
integer, dimension(:), allocatable :: id_cpre_cem, id_pb_cem, id_ptma_cem, &
                            id_h1_cem, id_qlw_cem, id_cfi_cem, &
                            id_wv_cem, id_rcl_cem
integer                  :: id_a1_cem, id_cual_cem, id_tfrc_cem, &
                            id_mpre_cem

integer, dimension(:), allocatable :: id_qtren1, id_qtmes1, &
                                      id_wtp1, id_qtceme, &
                                      id_total_wet_dep, &
                                      id_meso_wet_dep, id_cell_wet_dep
integer, dimension(:), allocatable :: id_qtren1_col, id_qtmes1_col, &
                                      id_wtp1_col, id_qtceme_col, &
                                      id_total_wet_dep_col, &
                                      id_meso_wet_dep_col,   &
                                      id_cell_wet_dep_col
integer, dimension(:), allocatable :: id_extremes, id_hits

integer, dimension(:), allocatable    :: id_water_budget, &
                                         id_ci_water_budget        
integer, dimension(:), allocatable    :: id_enthalpy_budget,   &
                                         id_ci_enthalpy_budget
integer, dimension (:,:), allocatable ::            &
                                         id_precip_budget, &
                                         id_ci_precip_budget

integer   :: id_ci_prcp_heat_liq_cell, id_ci_prcp_heat_frz_cell, &
             id_ci_prcp_heat_liq_meso, id_ci_prcp_heat_frz_meso, &
             id_ci_prcp_heat_total, id_ci_prcp_total

real              :: missing_value = -999.
character(len=16) :: mod_name = 'donner_deep'
integer           :: donner_axes(5)

!--------------------------------------------------------------------
!   variables for column diagnostics option
!
!   arrays containing information for all requested diagnostic columns
!   (1:num_diag_pts):
!    col_diag_unit         unit numbers for each column's output file 
!    col_diag_lon          each column's longitude 
!                          [ degrees, 0 < lon < 360 ]
!    col_diag_lat          each column's latitude 
!                          [degrees, -90 < lat < 90 ]
!    col_diag_j            each column's j index (processor coordinates)
!    col_diag_i            each column's i index (processor coordinates) 
!
!    Time_col_diagnostics  time in model simulation at which to activate
!                          column diagnostics 
!

integer, dimension(:), allocatable :: col_diag_unit
real   , dimension(:), allocatable :: col_diag_lon, col_diag_lat   
integer, dimension(:), allocatable :: col_diag_j, col_diag_i        
type(time_type)                    :: Time_col_diagnostics  


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

subroutine fms_donner_process_nml  (Nml, kpar)

!---------------------------------------------------------------------
!    fms_donner_process_nml processes the donner_deep_nml file. 
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

      integer                             :: unit, ierr, io, logunit
  
!-------------------------------------------------------------------
!  local variables:
!
!     unit                   unit number for nml file
!     ierr                   error return flag
!     io                     error return code
!                         
!-------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    1. READ NAMELIST AND WRITE IT TO LOG FILE.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=donner_deep_nml, iostat=io)
      ierr = check_nml_error(io,'donner_deep_nml')
#else   
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read (unit, nml=donner_deep_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'donner_deep_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number(version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() )    &
                                 write (logunit, nml=donner_deep_nml)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    8. STORE THE NAMELIST VARIABLES THAT NEED TO BE MADE AVAILABLE 
!       OUTSIDE OF THIS MODULE INTO THE DONNER_NML_TYPE VARIABLE.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      Nml%parcel_launch_level         = parcel_launch_level
      Nml%allow_mesoscale_circulation = allow_mesoscale_circulation
      Nml%do_hires_cape_for_closure =   do_hires_cape_for_closure
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
      Nml%entrainment_scheme_for_closure =    &
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

     allocate (Nml%arat(kpar))
     allocate (Nml%ensemble_entrain_factors_gate(kpar)) 

     if ( arat_erat_option /= 0 ) then
 
       call define_arat_erat (arat_erat_option, kpar, eratb, erat0, &
                              erat_min, erat_max, erat, arat)
       if (mpp_pe() == mpp_root_pe() ) then
         print *,'donner_deep_nml: redefined arat and erat using &
                        &arat_erat_option == ', arat_erat_option
         print *,'donner_deep_nml: arat = ',arat
         print *,'donner_deep_nml: erat = ',erat
       end if
     endif

     Nml%arat = arat
     Nml%ensemble_entrain_factors_gate = erat

end subroutine fms_donner_process_nml


!#####################################################################

subroutine fms_donner_process_tracers (Initialized, tracers_in_donner,&
                                       Don_save)

type(donner_initialized_type),   intent(inout) :: Initialized
logical, dimension(:),          intent(in) :: tracers_in_donner
type(donner_save_type), intent(inout)      :: Don_save



        integer :: nn, n
      logical                             :: flag
      character(len=200)                  :: method_name, method_control
      real                                :: frac_junk, frac_junk2

        Initialized%do_donner_tracer = .true.
        nn = 1
        do n=1,size(tracers_in_donner(:))
          if (tracers_in_donner(n)) then
            call get_tracer_names (MODEL_ATMOS, n,  &
                                   name = Don_save%tracername(nn), &
                                   units = Don_save%tracer_units(nn))
!++lwh
            Initialized%wetdep(nn)%units = Don_save%tracer_units(nn)
            flag = query_method( 'wet_deposition', MODEL_ATMOS, n, &
                                 method_name, method_control )
            call get_wetdep_param( method_name, method_control, &
                                   Initialized%wetdep(nn)%scheme, &
                                   Initialized%wetdep(nn)%Henry_constant, &
                                   Initialized%wetdep(nn)%Henry_variable, &
                                   frac_junk, frac_junk2, &
                                   Initialized%wetdep(nn)%alpha_r, &
                                   Initialized%wetdep(nn)%alpha_s , &
                                   Initialized%wetdep(nn)%Lwetdep, &
                                   Initialized%wetdep(nn)%Lgas, &
                                   Initialized%wetdep(nn)%Laerosol, &
                                   Initialized%wetdep(nn)%Lice, &
                                   frac_in_cloud_donner=Initialized%wetdep(nn)%frac_in_cloud)
            Initialized%wetdep(nn)%scheme = lowercase( Initialized%wetdep(nn)%scheme )
!-lwh
            nn = nn + 1
          endif
        end do

end subroutine fms_donner_process_tracers




!#####################################################################

subroutine fms_donner_activate_diagnostics (secs, days, axes, &
                             Don_save, Nml, n_water_budget, &
                             n_enthalpy_budget, n_precip_paths, &
                             n_precip_types, nlev_hires, kpar)

integer, intent(in) :: secs, days, n_water_budget, &
                             n_enthalpy_budget, n_precip_paths, &
                             n_precip_types, nlev_hires, kpar
                      
integer,         dimension(4),   intent(in)   :: axes
type(donner_save_type), intent(inout) :: Don_save
type(donner_nml_type), intent(inout) :: Nml       

      type(time_type)    :: Time

      Time = set_time (secs, days)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    4. INITIALIZE THE NETCDF OUTPUT VARIABLES.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!--------------------------------------------------------------------
!    activate the netcdf diagnostic fields.
!-------------------------------------------------------------------
      call register_fields (Time, axes, Don_save, Nml, &
                              n_water_budget, &
                             n_enthalpy_budget, n_precip_paths, &
                             n_precip_types, nlev_hires, kpar)


end subroutine fms_donner_activate_diagnostics 


!#####################################################################

subroutine fms_donner_read_restart (Initialized, ntracers,   &
                                    secs, days, Don_save, Nml)

type(donner_initialized_type), intent(inout) :: Initialized
type(donner_save_type), intent(inout) :: Don_save
type(donner_nml_type), intent(inout) :: Nml     
integer, intent(in) :: secs, days, ntracers

      type(time_type) :: Time
integer :: outunit

     Time = set_time (secs, days)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    5. PROCESS THE RESTART FILE.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!--------------------------------------------------------------------
!    if a netcdf restart file is present, call read_restart_nc to read 
!    it.
!--------------------------------------------------------------------
      !--- register restart field to be ready to be written out.
      call fms_donner_register_restart('donner_deep.res.nc', Initialized, ntracers, Don_save, Nml)

      if (file_exist ('INPUT/donner_deep.res.nc') ) then
       Initialized%coldstart= .false.
!        call read_restart_nc (ntracers, Initialized,Nml, Don_save)
       call restore_state(Don_restart)
       if (in_different_file) call restore_state(Til_restart)

!--------------------------------------------------------------------
!    if a native mode restart file is present, call read_restart 
!    to read it.
!--------------------------------------------------------------------
      else if (file_exist ('INPUT/donner_deep.res') ) then
        Initialized%coldstart= .false.
        call error_mesg ( 'fms_donner_mod', 'Native restart capability has been removed.', &
                                         FATAL)
!--------------------------------------------------------------------
!    if no restart file is present, call subroutine process_coldstart
!    to define the needed variables.
!--------------------------------------------------------------------
      else
        call process_coldstart (Time, Initialized, Nml, Don_save)
      endif


end subroutine fms_donner_read_restart 


!#####################################################################

subroutine fms_donner_col_diag (lonb, latb, Col_diag, pref) 

real, dimension(:,:), intent(in) :: lonb, latb
type(donner_column_diag_type), intent(inout) :: Col_diag
real, dimension(:), intent(in) :: pref

    logical, dimension(size(latb,1)-1, size(latb,2)-1) ::    &
                                                  do_column_diagnostics
    integer :: k, n

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

!---------------------------------------------------------------------
!    if any diagnostics are requested, perform various consistency
!    checks.
!---------------------------------------------------------------------
      if (Col_diag%num_diag_pts > 0) then

!---------------------------------------------------------------------
!    check that array dimensions are sufficiently large for the number 
!    of columns requested.
!---------------------------------------------------------------------
        if (Col_diag%num_diag_pts > MAX_PTS) then
          call error_mesg ('donner_deep_mod', 'donner_deep_init: &
         &must reset MAX_PTS or reduce number of diagnostic points', &
                                                           FATAL)  
        endif

!---------------------------------------------------------------------
!    check that the specified time at which diagnostics are to be 
!    activated has been specified.
!---------------------------------------------------------------------
        do n=1,3
          if (diagnostics_start_time(n) == 0) then
            call error_mesg ('donner_deep_mod', 'donner_deep_init:&
             &year, month and/or day invalidly specified for column '//&
                  'diagnostics starting time', FATAL)
          endif
        end do

!---------------------------------------------------------------------
!    define a time_type variable indicating the requested time to begin
!    outputting diagnostics.
!---------------------------------------------------------------------
        Time_col_diagnostics = set_date (diagnostics_start_time(1), &
                                         diagnostics_start_time(2), &   
                                         diagnostics_start_time(3), &   
                                         diagnostics_start_time(4), &   
                                         diagnostics_start_time(5), &   
                                         diagnostics_start_time(6) )    

!---------------------------------------------------------------------
!    allocate space for the arrays used to specify the diagnostics 
!    columns and the output units. initialize the arrays with bogus
!    values.
!---------------------------------------------------------------------
        allocate (col_diag_unit    (Col_diag%num_diag_pts) )
        allocate (col_diag_lon     (Col_diag%num_diag_pts) )
        allocate (col_diag_lat     (Col_diag%num_diag_pts) )
        allocate (col_diag_i       (Col_diag%num_diag_pts) )
        allocate (col_diag_j       (Col_diag%num_diag_pts) )
        col_diag_unit  = -1
        col_diag_lon   = -1.0
        col_diag_lat   = -1.0
        col_diag_i     = -1
        col_diag_j     = -1

!---------------------------------------------------------------------
!    call initialize_diagnostic_columns to determine the locations 
!    (i,j,lat and lon) of any diagnostic columns in this processor's
!    space and to open output files for the diagnostics.
!---------------------------------------------------------------------
        call initialize_diagnostic_columns   &
                     (mod_name, num_diag_pts_latlon, num_diag_pts_ij, &
                      i_coords_gl, j_coords_gl, lat_coords_gl, &
                      lon_coords_gl, lonb(:,:), latb(:,:),  &
                      do_column_diagnostics, &
                      col_diag_lon, col_diag_lat, col_diag_i,  &
                      col_diag_j, col_diag_unit)

!---------------------------------------------------------------------
!    verify that requested pressure cutoff for column diagnostics output
!    is valid. define the model k index which corresponds (kstart).
!---------------------------------------------------------------------
        do k=1,size(pref(:))
          if (pref(k) >= diagnostics_pressure_cutoff) then
            Col_diag%kstart = k
            exit
          endif
        end do

!----------------------------------------------------------------------
!    if the specified pressure is larger than any pressure level in the
!    model grid, write an error message.
!----------------------------------------------------------------------
        if (Col_diag%kstart == -99) then
          call error_mesg ( 'donner_deep_mod', 'donner_deep_init: &
           &diagnostics_pressure_cutoff is higher than pressure at '//&
                                     'any model level', FATAL)
        endif

!----------------------------------------------------------------------
!   if column diagnostics is not requested, define the components of
!   Col_diag that will be needed.
!----------------------------------------------------------------------
      else
        Col_diag%in_diagnostics_window = .false.
        Col_diag%ncols_in_window = 0
      endif

!----------------------------------------------------------------------
!    allocate space for the array elements of the donner_column_diag_type
!    variable Col_diag. These arrays remain for the life of the job and
!    will be defined for each physics window as it is entered.
!----------------------------------------------------------------------
      allocate (Col_diag%i_dc(Col_diag%num_diag_pts))
      allocate (Col_diag%j_dc(Col_diag%num_diag_pts))
      allocate (Col_diag%unit_dc(Col_diag%num_diag_pts))
      allocate (Col_diag%igl_dc(Col_diag%num_diag_pts))
      allocate (Col_diag%jgl_dc(Col_diag%num_diag_pts))


end subroutine fms_donner_col_diag 




!#####################################################################
! <SUBROUTINE NAME="fms_donner_write_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine fms_donner_write_restart (Initialized, timestamp)
  type(donner_initialized_type), intent(in) :: Initialized
  character(len=*), intent(in), optional :: timestamp

!-------------------------------------------------------------------
!    call subroutine to write restart file. NOTE: only the netcdf 
!    restart file is currently supported.
!-------------------------------------------------------------------
      if (mpp_pe() == mpp_root_pe() ) then
        if (.not. (write_reduced_restart_file) ) then
          call error_mesg ('donner_deep_mod', 'write_restart_nc: &
            &Writing FULL netCDF formatted restart file as requested: &
                 &RESTART/donner_deep.res.nc', NOTE)
        else
          if (Initialized%conv_alarm >= Initialized%physics_dt)  then
            call error_mesg ('donner_deep_mod', 'write_restart_nc: &
            &Writing FULL netCDF formatted restart file; it is needed &
             &to allow seamless restart because next step is not a &
             &donner calculation step: RESTART/donner_deep.res.nc', NOTE)
          else
            call error_mesg ('donner_deep_mod', 'write_restart_nc: &
              &Writing REDUCED netCDF formatted restart file as  &
                &requested: RESTART/donner_deep.res.nc', NOTE)
          endif
        endif
      endif
      call save_restart(Don_restart, timestamp)
      if(in_different_file) call save_restart(Til_restart, timestamp)

end subroutine fms_donner_write_restart 


!#####################################################################

subroutine fms_get_pe_number (me, root_pe)

integer, intent(out) :: me, root_pe

    me = mpp_pe()
    root_pe = mpp_root_pe()

end subroutine fms_get_pe_number



!#####################################################################

subroutine fms_close_col_diag_units


      call close_column_diagnostics_units (col_diag_unit)


end subroutine fms_close_col_diag_units 


!#####################################################################

subroutine fms_deallocate_variables (Col_diag)

type(donner_column_diag_type), intent(inout) :: Col_diag

      if (Col_diag%num_diag_pts > 0) then
        deallocate (Col_diag%i_dc    )
        deallocate (Col_diag%j_dc    ) 
        deallocate (Col_diag%unit_dc )
        deallocate (Col_diag%igl_dc  )
        deallocate (Col_diag%jgl_dc  )
      endif

      if (allocated(col_diag_unit)) then
        deallocate (col_diag_unit  )
        deallocate (col_diag_lon   )
        deallocate (col_diag_lat   )
        deallocate (col_diag_i     )
        deallocate (col_diag_j     )
      endif
     
      if (allocated (id_qtren1)) then
        deallocate (id_qtren1)
        deallocate (id_qtmes1)
        deallocate (id_wtp1  )
        deallocate (id_qtceme)
        deallocate (id_total_wet_dep)
        deallocate (id_meso_wet_dep)
        deallocate (id_cell_wet_dep)
        deallocate (id_qtren1_col)
        deallocate (id_qtmes1_col)
        deallocate (id_wtp1_col  )
        deallocate (id_qtceme_col)
        deallocate (id_total_wet_dep_col)
        deallocate (id_meso_wet_dep_col)
        deallocate (id_cell_wet_dep_col)
      endif


      if (allocated (id_extremes)) then
        deallocate (id_extremes)
        deallocate (id_hits)
      endif


end subroutine fms_deallocate_variables



!#####################################################################

subroutine fms_sat_vapor_pres


      call sat_vapor_pres_init


end subroutine fms_sat_vapor_pres


!#####################################################################

subroutine fms_error_mesg (ermesg)

character(len=*), intent(in) :: ermesg


      call error_mesg ('donner_deep_mod', ermesg, FATAL)



end subroutine fms_error_mesg 


!######################################################################

subroutine fms_donner_deep_netcdf (is, ie, js, je, Nml, secs, days, &
                               Param, Initialized, Don_conv, Don_cape,&
                               Don_cem,parcel_rise, pmass, total_precip, &
                               Don_budgets, &
                               temperature_forcing, moisture_forcing)  

!---------------------------------------------------------------------
!    subroutine donner_deep_netcdf sends the fields requested in the
!    diag_table to diag_manager_mod so that they may be appropriately
!    processed for output.
!---------------------------------------------------------------------

integer,                intent(in) :: is, ie, js, je
integer,                intent(in) :: secs, days
type(donner_nml_type), intent(in) :: Nml   
type(donner_param_type), intent(in) :: Param 
type(donner_initialized_type), intent(inout) :: Initialized
type(donner_conv_type), intent(in) :: Don_conv
type(donner_budgets_type), intent(in) :: Don_budgets
type(donner_cape_type), intent(in) :: Don_cape
type(donner_cem_type),  intent(in) :: Don_cem
real, dimension(:,:,:), intent(in) :: pmass, temperature_forcing,&
                                      moisture_forcing
real, dimension(:,:),   intent(in) :: parcel_rise, total_precip

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     is, ie         first and last values of i index values of points 
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points 
!                    in this physics window (processor coordinates)
!     Time           current time (time_type)
!     Don_conv       donner_convection_type derived type variable con-
!                    taining diagnostics describing the nature of the 
!                    convection produced by the donner parameterization
!     Don_cape       donner_cape type derived type variable containing
!                    diagnostics related to the cape calculation assoc-
!                    iated with the donner convection parameterization
!     Don_cem        donner_cem_type derived type variable containing
!                    Donner cumulus ensemble member diagnostics
!     temperature_forcing  
!                    temperature tendency due to donner convection
!                    [ deg K / sec ]
!     moisture_forcing  
!                    vapor mixing ratio tendency due to donner 
!                    convection [ kg(h2o) / (kg(dry air) sec ) ]
!     pmass          mass per unit area within the grid box
!                    [ kg (air) / (m**2) ]
!     parcel_rise    accumulated vertical displacement of a near-surface
!                    parcel as a result of the lowest model level omega 
!                    field [ Pa ]
!     total_precip   total precipitation rate produced by the
!                    donner parameterization [ mm / day ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      type(time_type)  :: Time

      Time = set_time (secs, days)

      call donner_deep_netcdf (is, ie, js, je, Nml, Time,  Param, &
                               Initialized, Don_conv, Don_cape,&
                               Don_cem,parcel_rise, pmass, total_precip, &
                               Don_budgets, &
                               temperature_forcing, moisture_forcing)  


!----------------------------------------------------------------------


end subroutine fms_donner_deep_netcdf



!###################################################################

subroutine fms_donner_process_monitors (idf, jdf, nlev,  &
                              ntracers, axes, secs, days, Initialized,&
                              Don_save)

integer, intent(in)  :: idf, jdf, nlev, ntracers, secs, days
integer, dimension(4), intent(in) :: axes
type(donner_save_type), intent(inout) :: Don_save
type(donner_initialized_type), intent(inout) :: Initialized


   type(time_type)  :: Time

    Time = set_time (secs,days)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    9. SET UP CODE TO MONITOR SELECTED OUTPUT VARIABLES.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

       call process_monitors (idf, jdf, nlev, ntracers, axes, Time, &
                              Initialized, Don_save)


end subroutine fms_donner_process_monitors

!###################################################################

subroutine fms_donner_column_control (is, ie, js, je, secs, days, Col_diag)          

!---------------------------------------------------------------------
!    subroutine fms_donner_column_control returns the number, location
!    (processor and window indices) and output units associated with 
!    any diagnostic columns requested within the current physics window.
!---------------------------------------------------------------------

integer,                       intent(in)   :: is, ie, js, je
integer,                       intent(in) :: secs, days
type(donner_column_diag_type), intent(inout) :: Col_diag


     type(time_type)                             :: Time

      Time = set_time(secs, days)

      call donner_column_control (is, ie, js, je, Time, Col_diag)                


end subroutine fms_donner_column_control 


!####################################################################


subroutine fms_constants (Param)

type(donner_param_type), intent(inout)  :: Param


!----------------------------------------------------------------------
!    define the components of Param that come from constants_mod. see 
!    donner_types.h for their definitions.
!----------------------------------------------------------------------
      Param%dens_h2o        = DENS_H2O
      Param%rdgas           = RDGAS
      Param%grav            = GRAV
      Param%cp_air          = CP_AIR  
      Param%pie             = PIE
      Param%kappa           = KAPPA
      Param%rvgas           = RVGAS
      Param%seconds_per_day = SECONDS_PER_DAY
      Param%hlv             = HLV
      Param%hlf             = HLF
      Param%hls             = HLS
      Param%kelvin          = KELVIN

!----------------------------------------------------------------------

end subroutine fms_constants 

!####################################################################





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
 
subroutine register_fields (Time, axes, Don_save, Nml, &
                              n_water_budget, &
                             n_enthalpy_budget, n_precip_paths, &
                             n_precip_types, nlev_hires, kpar)

!----------------------------------------------------------------------
!    subroutine register_fields registers all of the potential diagnos-
!    tics written by this module with diag_manager_mod.
!----------------------------------------------------------------------

type(time_type),               intent(in)   :: Time
integer, intent(in) :: n_water_budget, &
                             n_enthalpy_budget, n_precip_paths, &
                             n_precip_types, nlev_hires, kpar
integer,         dimension(4), intent(in)   :: axes
type(donner_save_type), intent(inout) :: Don_save
type(donner_nml_type), intent(inout) :: Nml       

!----------------------------------------------------------------------
!   intent(in) variables:
!
!      Time         current time [ time_type ]
!      axes         data axes for diagnostics
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer  :: ntracers     ! number of tracers transported by the
                               ! donner deep convection parameterization
      integer :: nn            ! do-loop index
      integer :: ncem          ! number of cumulus ensemble members in
                               ! the donner deep convection parameter-
                               ! ization
      character(len=2) :: chvers ! character representation of cumulus 
                                 ! ensemble  member number

!  define a variable for telling "register_fields" to put output on
!  "half-levels" (Reference:  Chris Golaz's subroutine "diag_field_init"
!  in /home/cjg/FMS/nalanda/nnew3/m45_am2p14_nnew3/src/atmos_param/
!                                 moist_processes/moist_processes.F90)
 
      integer, dimension(3) :: half = (/1,2,4/)
      integer, dimension(3) :: cldindices = (/1,2,5/)
      integer               :: id_cldmodel
      real                  :: cldvindx(NLEV_HIRES)
      integer               :: k

      ncem = kpar

!---------------------------------------------------------------------  
!    define the axes for the donner cloud model.
!---------------------------------------------------------------------
      donner_axes(1:4) = axes(1:4)     
      if (Nml%do_donner_plume) then
        do k=1, NLEV_HIRES
          cldvindx(k) = real(k)
        end do
        id_cldmodel = diag_axis_init('cldvindx', cldvindx, 'level#', &
                                     'z', 'cld model vertical index', &
                                     set_name=mod_name )
        donner_axes(5) = id_cldmodel
      endif

!----------------------------------------------------------------------
!    define the number of tracers that are to be transported by the 
!    donner deep convection parameterization.
!-------------------------------------------------------------------
      ntracers = size(Don_save%tracername(:))

!---------------------------------------------------------------------
!    register the various diagnostic fields.
!---------------------------------------------------------------------

    if (Nml%do_budget_analysis) then
      allocate (id_water_budget (n_water_budget))
      allocate (id_ci_water_budget (n_water_budget))
      allocate (id_enthalpy_budget (n_enthalpy_budget))
      allocate (id_ci_enthalpy_budget (n_enthalpy_budget))
      allocate (id_precip_budget (n_precip_paths, n_precip_types))
      allocate (id_ci_precip_budget (n_precip_paths, n_precip_types))
      id_water_budget(1)    = register_diag_field    &
            (mod_name, 'vapor_net_tend', axes(1:3),   &
             Time, 'net water vapor tendency', &
             'g(h2o) / kg(air) / day',    &
             missing_value=missing_value)
      
      id_water_budget(2)    = register_diag_field    &
            (mod_name, 'vapor_cell_dynam', axes(1:3),   &
             Time, 'vapor tendency due to cell dynamics', &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_water_budget(3)    = register_diag_field    &
            (mod_name, 'vapor_meso_depo', axes(1:3),   &
             Time, 'vapor tendency from mesoscale deposition', &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_water_budget(4)    = register_diag_field    &
            (mod_name, 'vapor_meso_cd', axes(1:3),   &
             Time, 'vapor tendency from mesoscale condensation',  &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_water_budget(5)    = register_diag_field    &
            (mod_name, 'vapor_cell_evap', axes(1:3),   &
             Time, 'vapor tendency from cell evaporation',  &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_water_budget(6)    = register_diag_field    &
            (mod_name, 'vapor_cell_meso_trans', axes(1:3),   &
             Time, 'vapor tendency from cell to mesoscale transfer',  &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_water_budget(7)    = register_diag_field    &
            (mod_name, 'vapor_meso_evap', axes(1:3),   &
             Time, 'vapor tendency from mesoscale evaporation', &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_water_budget(8)    = register_diag_field    &
            (mod_name, 'vapor_meso_dynam_up', axes(1:3),   &
             Time, 'vapor tendency from mesoscale updrafts',  &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_water_budget(9)    = register_diag_field    &
            (mod_name, 'vapor_meso_dynam_dn',  axes(1:3),   &
             Time, 'vapor tendency from mesoscale downdrafts',  &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_enthalpy_budget(1)    = register_diag_field    &
            (mod_name, 'enth_net_tend', axes(1:3),   &
             Time, 'net temp tendency', 'deg K  /day',    &
             missing_value=missing_value)

      id_enthalpy_budget(2)    = register_diag_field    &
            (mod_name, 'enth_cell_dynam', axes(1:3),   &
             Time, 'temp tendency due to cell dynamics', &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(3)    = register_diag_field    &
            (mod_name, 'enth_meso_depo_liq', axes(1:3), Time, &
             'temp tendency from mesoscale deposition on liquid&
                    & condensate',  'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(4)    = register_diag_field    &
            (mod_name, 'enth_meso_cd_liq', axes(1:3), Time, &
             ' temp tendency from mesoscale liquid condensation', &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(5)    = register_diag_field    &
            (mod_name, 'enth_cell_evap_liq', axes(1:3),   &
             Time, 'temp tendency from evap of liquid condensate', &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(6)    = register_diag_field    &
            (mod_name, 'enth_meso_evap_liq_up', axes(1:3),   &
             Time, 'temp tendency from evaporation of liquid &
              &condensate in mesoscale updrafts',  &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(7)    = register_diag_field    &
            (mod_name, 'enth_meso_evap_liq_dn', axes(1:3),   &
             Time, 'temp tendency from evaporation of liquid &
              &condensate in mesoscale downdrafts',  &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(8)    = register_diag_field    &
            (mod_name, 'enth_meso_depo_ice', axes(1:3),   &
             Time, ' temp tendency from mesoscale deposition on &
              &ice condensate',  'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(9)    = register_diag_field    &
            (mod_name, 'enth_meso_cd_ice', axes(1:3),   &
             Time, 'temp tendency from mesoscale ice condensation', &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(10)    = register_diag_field    &
            (mod_name, 'enth_cell_evap_ice', axes(1:3),   &
             Time, 'temp tendency from evap of ice condensate', &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(11)    = register_diag_field    &
            (mod_name, 'enth_meso_evap_ice_up', axes(1:3),   &
             Time, 'temp tendency from evaporation of ice condensate &
              &in mesoscale updrafts',  'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(12)    = register_diag_field    &
            (mod_name, 'enth_meso_evap_ice_dn', axes(1:3),   &
             Time, 'temp tendency from evaporation of ice &
               &condensate in mesoscale downdrafts',  'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(13)    = register_diag_field    &
            (mod_name, 'enth_meso_freeze', axes(1:3),   &
             Time, 'temp tendency from the freezing of liquid &
              &condensate when it enters the mesoscale circulation',  &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(14)    = register_diag_field    &
            (mod_name, 'enth_cell_freeze', axes(1:3),   &
             Time, 'temp tendency from the freezing of liquid &
                &cell condensate',  'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(15)    = register_diag_field    &
            (mod_name, 'enth_cell_precip_melt', axes(1:3),   &
             Time, 'temp tendency from the melting of cell frozen &
             &liquid and ice that is precipitating out', 'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(16)    = register_diag_field    &
            (mod_name, 'enth_meso_melt', axes(1:3), Time, &
             'temp tendency from melting bogus frozen condensate',  &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(17)    = register_diag_field    &
            (mod_name, 'enth_meso_precip_melt', axes(1:3),   &
             Time, 'temp tendency from the melting of frozen &
               &mesoscale precipitation',  'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(18)    = register_diag_field    &
            (mod_name, 'enth_meso_dynam_up', axes(1:3),   &
             Time, 'temp tendency from mesoscale updraft', &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(19)    = register_diag_field    &
            (mod_name, 'enth_meso_dynam_dn', axes(1:3),   &
             Time, 'temp tendency from mesoscale downdraft', &
             'deg K / day', &
             missing_value=missing_value)

      id_precip_budget(1,1)    = register_diag_field    &
            (mod_name, 'precip_cell_liq', axes(1:3),   &
             Time, 'precip from cell liquid condensate', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(2,1)    = register_diag_field    &
            (mod_name, 'precip_cell_liq_frz', axes(1:3),   &
             Time, 'precip from cell liquid condensate which froze', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(3,1)    = register_diag_field    &
            (mod_name, 'precip_cell_liq_frz_melt', axes(1:3), Time, &
              'precip from cell liquid condensate which froze &
               &and remelted', 'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(4,1)    = register_diag_field    &
            (mod_name, 'precip_cell_ice', axes(1:3),   &
             Time, 'precip from cell ice condensate', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(5,1)    = register_diag_field    &
            (mod_name, 'precip_cell_ice_melt', axes(1:3),   &
             Time, 'precip from cell ice condensate which melted', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(1,2)    = register_diag_field    &
            (mod_name, 'precip_trans_liq', axes(1:3),   &
             Time, 'precip from cell liquid transferred to meso', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(2,2)    = register_diag_field    &
            (mod_name, 'precip_trans_liq_frz', axes(1:3),   &
             Time, 'precip from cell liquid transferred to meso &
              &which froze', 'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(3,2)    = register_diag_field    &
            (mod_name, 'precip_trans_liq_frz_melt', axes(1:3), Time, &
             'precip from cell liquid transferred to meso which &
              &froze and remelted', 'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(4,2)    = register_diag_field    &
            (mod_name, 'precip_trans_ice', axes(1:3),   &
             Time, 'precip from cell ice transferred to meso', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(5,2)    = register_diag_field    &
            (mod_name, 'precip_trans_ice_melt', axes(1:3),   &
             Time, 'precip from cell ice transferred to meso &
              &which melted', 'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(1,3)    = register_diag_field    &
            (mod_name, 'precip_meso_liq', axes(1:3),   &
             Time, 'precip from meso liq condensate', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(2,3)    = register_diag_field    &
            (mod_name, 'precip_meso_liq_frz', axes(1:3),   &
             Time, 'precip from meso liq  condensate which froze', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(3,3)    = register_diag_field    &
            (mod_name, 'precip_meso_liq_frz_melt', axes(1:3), Time, &
            'precip from meso condensate liq which froze and &
             &remelted', 'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(4,3)    = register_diag_field    &
            (mod_name, 'precip_meso_ice', axes(1:3),   &
             Time, 'precip from meso ice condensate', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(5,3)    = register_diag_field    &
            (mod_name, 'precip_meso_ice_melt', axes(1:3),   &
             Time, 'precip from meso ice condensate which melted', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_ci_precip_budget(1,1)    = register_diag_field    &
            (mod_name, 'ci_precip_cell_liq', axes(1:2),   &
             Time, 'col intg precip from cell liquid condensate', &
             'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(2,1)    = register_diag_field    &
            (mod_name, 'ci_precip_cell_liq_frz', axes(1:2),   &
             Time, 'col intg precip from cell liquid condensate &
             &which froze',  'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(3,1)    = register_diag_field    &
            (mod_name, 'ci_precip_cell_liq_frz_melt', axes(1:2), Time, &
             'col intg precip from cell liquid condensate which &
              &froze and remelted',  'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(4,1)    = register_diag_field    &
            (mod_name, 'ci_precip_cell_ice', axes(1:2),   &
             Time, 'col intg precip from cell ice condensate', &
             'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(5,1)    = register_diag_field    &
            (mod_name, 'ci_precip_cell_ice_melt', axes(1:2),   &
             Time, 'col intg precip from cell ice condensate &
             &which melted',  'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(1,2)    = register_diag_field    &
            (mod_name, 'ci_precip_trans_liq', axes(1:2),   &
             Time, 'col intg precip from cell liquid transferred &
             &to meso',  'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(2,2)    = register_diag_field    &
            (mod_name, 'ci_precip_trans_liq_frz', axes(1:2),   &
             Time, 'col intg precip from cell liquid transferred &
              &to meso  which froze', 'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(3,2)    = register_diag_field    &
            (mod_name, 'ci_precip_trans_liq_frz_melt', axes(1:2), &
             Time, 'col intg precip from cell liquid transferred &
             &to meso which froze and remelted', 'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(4,2)    = register_diag_field    &
            (mod_name, 'ci_precip_trans_ice', axes(1:2),   &
             Time, 'col intg precip from cell ice transferred &
              &to meso', 'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(5,2)    = register_diag_field    &
            (mod_name, 'ci_precip_trans_ice_melt', axes(1:2),   &
             Time, 'col intg precip from cell ice transferred to &
             &meso which melted', 'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(1,3)    = register_diag_field    &
            (mod_name, 'ci_precip_meso_liq', axes(1:2),   &
             Time, 'col intg precip from meso liq condensate', &
             'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(2,3)    = register_diag_field    &
            (mod_name, 'ci_precip_meso_liq_frz', axes(1:2),   &
             Time, 'col intg precip from meso liq  condensate &
             &which froze',  'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(3,3)    = register_diag_field    &
            (mod_name, 'ci_precip_meso_liq_frz_melt', axes(1:2), Time, &
             'col intg precip from meso condensate liq which froze &
               &and remelted', 'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(4,3)    = register_diag_field    &
            (mod_name, 'ci_precip_meso_ice', axes(1:2),   &
             Time, 'col intg precip from meso ice condensate', &
             'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(5,3)    = register_diag_field    &
            (mod_name, 'ci_precip_meso_ice_melt', axes(1:2),   &
             Time, 'col intg precip from meso ice condensate &
              &which melted', 'mm / day', &
             missing_value=missing_value)

      id_ci_water_budget(1)    = register_diag_field    &
            (mod_name, 'ci_vapor_net_tend', axes(1:2),   &
             Time, 'col intg net water vapor tendency', 'mm / day',    &
             missing_value=missing_value)
      
      id_ci_water_budget(2)    = register_diag_field    &
            (mod_name, 'ci_vapor_cell_dynam', axes(1:2),   &
             Time, 'col intg vapor tendency due to cell dynamics', &
              'mm / day',    &
             missing_value=missing_value)
      
      id_ci_water_budget(3)    = register_diag_field    &
            (mod_name, 'ci_vapor_meso_depo', axes(1:2),   &
             Time, 'col intg vapor tendency from mesoscale deposition',&
              'mm / day',    &
             missing_value=missing_value)
      
      id_ci_water_budget(4)    = register_diag_field    &
            (mod_name, 'ci_vapor_meso_cd', axes(1:2),   &
             Time, 'col intg vapor tendency from mesoscale &
              &condensation',  'mm / day',    &
             missing_value=missing_value)
      
      id_ci_water_budget(5)    = register_diag_field    &
            (mod_name, 'ci_vapor_cell_evap', axes(1:2),   &
             Time, 'col intg vapor tendency from cell evaporation', &
              'mm / day', missing_value=missing_value)
      
      id_ci_water_budget(6)    = register_diag_field    &
            (mod_name, 'ci_vapor_cell_meso_trans', axes(1:2),   &
             Time, 'col intg vapor tendency from cell to mesoscale &
              &transfer',  'mm / day',    &
             missing_value=missing_value)
      
      id_ci_water_budget(7)    = register_diag_field    &
            (mod_name, 'ci_vapor_meso_evap', axes(1:2),   &
             Time, 'col intg vapor tendency from mesoscale &
              &evaporation', 'mm / day',    &
             missing_value=missing_value)
      
      id_ci_water_budget(8)    = register_diag_field    &
            (mod_name, 'ci_vapor_meso_dynam_up', axes(1:2),   &
             Time, 'col intg vapor tendency from mesoscale updrafts',  &
              'mm / day',    &
             missing_value=missing_value)
      
      id_ci_water_budget(9)    = register_diag_field    &
            (mod_name, 'ci_vapor_meso_dynam_dn',  axes(1:2),   &
             Time, 'col intg vapor tendency from mesoscale downdrafts',&
              'mm / day',    &
             missing_value=missing_value)
      
      id_ci_enthalpy_budget(1)    = register_diag_field    &
            (mod_name, 'ci_enth_net_tend', axes(1:2),   &
             Time, 'col intg net enthalpy tendency', 'J/m**2 / day',   &
             missing_value=missing_value)

      id_ci_enthalpy_budget(2)    = register_diag_field    &
            (mod_name, 'ci_enth_cell_dynam', axes(1:2),   &
             Time, 'col intg enthalpy tendency due to cell dynamics', &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(3)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_depo_liq', axes(1:2),   &
             Time, 'col intg enthalpy tendency from mesoscale &
             &deposition on liquid condensate',  'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(4)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_cd_liq', axes(1:2),   &
             Time, 'col intg enthalpy tendency from mesoscale &
             &liquid condensation', 'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(5)    = register_diag_field    &
            (mod_name, 'ci_enth_cell_evap_liq', axes(1:2),   &
             Time, 'col intg enthalpy tendency from evap of liquid &
             &condensate',  'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(6)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_evap_liq_up', axes(1:2),   &
             Time, 'col intg enthalpy tendency from evaporation of &
             &liquid condensate in mesoscale updrafts',  &
             'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(7)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_evap_liq_dn', axes(1:2),   &
             Time, 'col intg enthalpy tendency from evaporation &
             &of liquid condensate in mesoscale downdrafts',  &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(8)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_depo_ice', axes(1:2),   &
             Time, 'col intg enthalpy tendency from mesoscale &
              &deposition on ice condensate',  &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(9)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_cd_ice', axes(1:2),   &
             Time, 'col intg enthalpy tendency from mesoscale ice &
             &condensation', 'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(10)    = register_diag_field    &
            (mod_name, 'ci_enth_cell_evap_ice', axes(1:2),   &
             Time, 'col intg enthalpy tendency from evap of ice &
              &condensate', 'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(11)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_evap_ice_up', axes(1:2),   &
             Time, 'col intg enthalpy tendency from evaporation of &
             &ice condensate in mesoscale updrafts',  'J/m**2 / day',  &
             missing_value=missing_value)

      id_ci_enthalpy_budget(12)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_evap_ice_dn', axes(1:2),   &
             Time, 'col intg enthalpy tendency from evaporation of &
             &ice condensate in mesoscale downdrafts',  &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(13)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_freeze', axes(1:2),   &
             Time, 'col intg enthalpy tendency from the freezing of &
             &liquid condensate when it enters the mesoscale &
             &circulation',  'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(14)    = register_diag_field    &
            (mod_name, 'ci_enth_cell_freeze', axes(1:2),   &
             Time, 'col intg enthalpy tendency from the freezing of &
             &liquid cell condensate', 'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(15)    = register_diag_field    &
            (mod_name, 'ci_enth_cell_precip_melt', axes(1:2),   &
             Time, 'col intg enthalpy tendency from the melting of &
             &cell frozen liquid and ice that is precipitating out',  &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(16)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_melt', axes(1:2),   &
             Time, 'col intg enthalpy tendency from melting bogus &
              &frozen condensate',  'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(17)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_precip_melt', axes(1:2),   &
             Time, 'col intg enthalpy tendency from the melting of &
              &frozen mesoscale precipitation',  &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(18)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_dynam_up', axes(1:2),   &
             Time, 'col intg enthalpy tendency from mesoscale updraft',&
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(19)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_dynam_dn', axes(1:2),   &
             Time, 'col intg enthalpy tendency from mesoscale &
             &downdraft',  'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_prcp_heat_frz_cell =  register_diag_field & 
            (mod_name, 'ci_prcp_heat_frz_cell', axes(1:2),   &
             Time, 'col intg heat removed by frozen cell precip', &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_prcp_heat_liq_cell =  register_diag_field & 
            (mod_name, 'ci_prcp_heat_liq_cell', axes(1:2),   &
             Time, 'col intg heat removed by liquid cell precip', &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_prcp_heat_frz_meso =  register_diag_field & 
            (mod_name, 'ci_prcp_heat_frz_meso', axes(1:2),   &
             Time, 'col intg heat removed by frozen meso precip', &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_prcp_heat_liq_meso =  register_diag_field & 
            (mod_name, 'ci_prcp_heat_liq_meso', axes(1:2),   &
             Time, 'col intg heat removed by liquid meso precip', &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_prcp_heat_total =  register_diag_field & 
            (mod_name, 'ci_prcp_heat_total', axes(1:2),   &
             Time, 'col intg total heat removed by precip', &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_prcp_total =  register_diag_field & 
            (mod_name, 'ci_prcp_total', axes(1:2),   &
             Time, 'col intg total precip', &
              'mm / day',    &
             missing_value=missing_value)

    endif

      id_leff          = register_diag_field    &
            (mod_name, 'leff_don', axes(1:2),   &
             Time, 'effective latent heat with donner precip ',  &
             'J/kg(h2o)',  missing_value=missing_value)

!    heating rate:
      id_cemetf_deep = register_diag_field    &
            (mod_name, 'cemetf_deep', axes(1:3),   &
             Time, 'heating rate, c + m ', 'K/s',   &
             missing_value=missing_value)

!    cell entropy flux convergence:
      id_ceefc_deep = register_diag_field   &
            (mod_name, 'ceefc_deep', axes(1:3),   &
             Time, 'cell entrpy flx cnvrgnc', 'K/s',   &
             missing_value=missing_value)

!    cell condensation / evaporation:
      id_cecon_deep = register_diag_field      &
            (mod_name, 'cecon_deep', axes(1:3),   &
             Time, 'cell cond/evap ', 'K/s',   &
             missing_value=missing_value)

!    cell moisture flux convergence:
      id_cemfc_deep = register_diag_field       &
            (mod_name, 'cemfc_deep', axes(1:3),   &
             Time, 'cell moist flx cnvgnc', 'kg(h2o)/kg/s',   &
             missing_value=missing_value)

!    moistening rate:
      id_cememf_deep = register_diag_field        &
            (mod_name, 'cememf_deep', axes(1:3),   &
             Time, 'moistening rate, c + m ', 'kg(h2o)/kg/s',   &
             missing_value=missing_value)

!    moistening rate after adjustment for negative vapor mixing ratio:
      id_cememf_mod_deep = register_diag_field       &
            (mod_name, 'cememf_mod_deep', axes(1:3),&
             Time, 'mod cememf due to negative q ', 'kg(h2o)/kg/s',   &
             missing_value=missing_value)

!    cell + mesoscale cloud fraction:
      id_cual_deep = register_diag_field      &
            (mod_name, 'cual_deep', axes(1:3),   &
             Time, 'c + m cld frac ', 'percent',   &
             missing_value=missing_value)

!    heating rate due to freezing:
      id_fre_deep = register_diag_field     &
            (mod_name, 'fre_deep', axes(1:3),   &
             Time, 'freezing ', 'K/sec',   &
             missing_value=missing_value)

!    heating rate due to melting:
      id_elt_deep = register_diag_field         &
            (mod_name, 'elt_deep', axes(1:3),   &
             Time, 'melting', 'K/sec',   &
             missing_value=missing_value)

!    deposition in mesoscale updraft:
      id_cmus_deep = register_diag_field        &
            (mod_name, 'cmus_deep', axes(1:3),   &
             Time, 'meso-up deposition', 'kg(h2o)/kg/sec)',   &
             missing_value=missing_value)

!    evaporation in convective downdraft:
      id_ecds_deep = register_diag_field    &
            (mod_name, 'ecds_deep', axes(1:3),   &
             Time, 'convective dwndrft evap ', 'kg(h2o)/kg/sec', &
             missing_value=missing_value)

!    evaporation / sublimation in convective updraft:
      id_eces_deep = register_diag_field       &
            (mod_name, 'eces_deep', axes(1:3),   &
             Time, 'convective updrft evap/subl ', 'kg(h2o)/kg/sec',  &
             missing_value=missing_value)

!    sublimation in mesoscale downdraft:
      id_emds_deep = register_diag_field     &
            (mod_name, 'emds_deep', axes(1:3),   &
             Time, 'meso-dwn subl ', 'kg(h2o)/kg/sec',   &
             missing_value=missing_value)

!    sublimation in mesoscale updraft:
      id_emes_deep = register_diag_field        &
            (mod_name, 'emes_deep', axes(1:3),   &
             Time, 'meso-up subl ', 'kg(h2o)/kg/sec',   &
             missing_value=missing_value)

!    mesoscale moisture flux convergence:
      id_qmes_deep = register_diag_field     &
             (mod_name, 'qmes_deep', axes(1:3),   &
              Time, 'meso moist flux conv', 'kg(h2o)/kg/sec',   &
              missing_value=missing_value)

!    transfer of vapor from cells to mesoscale:
      id_wmps_deep = register_diag_field      &
             (mod_name, 'wmps_deep', axes(1:3),   &
              Time, 'meso redistrib of vapor from cells',  &
              'kg(h2o)/kg/sec', missing_value=missing_value)

!    deposition of vapor from cells to mesoscale:
      id_wmms_deep = register_diag_field         &
             (mod_name, 'wmms_deep', axes(1:3),   &
              Time, 'meso depo of vapor from cells',    &
              'kg(h2o)/kg/sec',  missing_value=missing_value)

!    mesoscale entropy flux convergesnce:
      id_tmes_deep = register_diag_field         &
            (mod_name, 'tmes_deep', axes(1:3),   &
             Time, 'meso entropy flux conv',  'K/sec',   &
              missing_value=missing_value)
 
!    mass flux in mesoscale downdrafts:    
      id_dmeml_deep = register_diag_field      &
            (mod_name, 'dmeml_deep', axes(1:3), &  
             Time, 'mass flux meso dwndrfts', 'kg/((m**2) s)',   &
             missing_value=missing_value)

!    mass flux in cell updrafts:
      id_uceml_deep = register_diag_field     &
            (mod_name, 'uceml_deep', axes(1:3), &
             Time, 'mass flux cell updrfts', 'kg/((m**2) s)',   &
             missing_value=missing_value)

!    mass flux in mesoscale updrafts:
      id_umeml_deep = register_diag_field       &
            (mod_name, 'umeml_deep', axes(1:3), &
             Time, 'mass flux meso updrfts', 'kg/((m**2) s)',   &
             missing_value=missing_value)

!    mesoscale ice mass mixing ratio:
      id_xice_deep = register_diag_field     &
            (mod_name, 'xice_deep', axes(1:3),  &
             Time, 'meso ice mass mixing ratio ', 'kg(ice)/kg',   &
             missing_value=missing_value)

!    mesoscale liquid mass mixing ratio:
      id_xliq_deep = register_diag_field       &
            (mod_name, 'xliq_deep', axes(1:3),  &
             Time, 'meso liq mass mixing ratio ', 'kg(liq)/kg',   &
             missing_value=missing_value)

!    detrained mass flux:
      id_detmfl_deep = register_diag_field       &
            (mod_name, 'detmfl_deep', axes(1:3),  &
             Time, 'detrained mass flux ', 'kg/((m**2) s)',   &
             missing_value=missing_value)

!---------------------------------------------------------------------
!    if tracers are being transported by donner_deep_mod, allocate diag-
!    nostic indices for each tracer and register their diagnostics.
!---------------------------------------------------------------------
      if (ntracers > 0) then
        allocate (id_qtren1 (ntracers))
        allocate (id_qtmes1 (ntracers))
        allocate (id_wtp1   (ntracers))
        allocate (id_qtceme (ntracers))
        allocate (id_total_wet_dep (ntracers))
        allocate (id_meso_wet_dep (ntracers))
        allocate (id_cell_wet_dep (ntracers))
        allocate (id_qtren1_col (ntracers))
        allocate (id_qtmes1_col (ntracers))
        allocate (id_wtp1_col   (ntracers))
        allocate (id_qtceme_col (ntracers))
        allocate (id_total_wet_dep_col (ntracers))
        allocate (id_meso_wet_dep_col (ntracers))
        allocate (id_cell_wet_dep_col (ntracers))
        do nn=1,ntracers

!    tracer tendency due to cells:
          id_qtren1(nn) = register_diag_field     &
                (mod_name, trim(Don_save%tracername(nn)) // '_qtren1',  &
                 axes(1:3), Time,  &
                 trim(Don_save%tracername(nn)) // ' cell tendency ', &
                 trim(Don_save%tracer_units(nn))//'/s', &
                 missing_value=missing_value)

!    tracer tendency due to mesoscale circulation:
          id_qtmes1(nn) = register_diag_field    &
                (mod_name, trim(Don_save%tracername(nn)) // '_qtmes1', &
                 axes(1:3), Time,   &
                 trim(Don_save%tracername(nn)) //' mesoscale tendency',&
                 trim(Don_save%tracer_units(nn))//'/s', &
                 missing_value=missing_value)

!    tracer tendency due to mesoscale redistribution:
          id_wtp1(nn) = register_diag_field         &
                (mod_name, trim(Don_save%tracername(nn)) // '_wtp1',  &
                 axes(1:3), Time,  &
                 trim(Don_save%tracername(nn)) //' mesoscale redist',&
                 trim(Don_save%tracer_units(nn))//'/s', &
                 missing_value=missing_value)

 !    tracer tendency due to deep convective wet deposition:
          id_total_wet_dep(nn) = register_diag_field         &
              (mod_name, trim(Don_save%tracername(nn)) // '_totwdep',  &
                 axes(1:3), Time,  &
                 trim(Don_save%tracername(nn)) //' deep conv wet depo',&
                 trim(Don_save%tracer_units(nn))//'/s', &
                 missing_value=missing_value)

!    tracer tendency due to wet deposition in mesoscale updrafts:
         id_meso_wet_dep(nn) = register_diag_field         &
                 (mod_name, trim(Don_save%tracername(nn)) // '_mwdep', &
                  axes(1:3), Time,   &
                 trim(Don_save%tracername(nn)) //' mesoscale wet depo',&
                 trim(Don_save%tracer_units(nn))//'/s', &
                 missing_value=missing_value)

!    tracer tendency due to wet deposition in cells:
          id_cell_wet_dep(nn) = register_diag_field         &
                (mod_name, trim(Don_save%tracername(nn)) // '_cwdep', &
                 axes(1:3), Time,  &
                 trim(Don_save%tracername(nn)) //' cell wet depo',&
                 trim(Don_save%tracer_units(nn))//'/s', &
                 missing_value=missing_value)

!    total tracer tendency:
          id_qtceme(nn) = register_diag_field     &
                (mod_name, trim(Don_save%tracername(nn)) // '_qtceme', &
                 axes(1:3), Time,  &
                 trim(Don_save%tracername(nn)) // ' total tendency ',&
                 trim(Don_save%tracer_units(nn))//'/s', &
                 missing_value=missing_value)

!    column-integrated tracer tendency due to cells:
          id_qtren1_col(nn) = register_diag_field      &
                (mod_name,       &
                 trim(Don_save%tracername(nn)) // '_qtren1_col',  &
                 axes(1:2), Time,  & 
                 'column integrated ' //trim(Don_save%tracername(nn)) //&
                 ' cell tendency ', &
                 trim(Don_save%tracer_units(nn)) // '* kg/(m**2 s) ', &
                 missing_value=missing_value)

!    column-integrated tracer tendency due to mesoscale circulation:
          id_qtmes1_col(nn) = register_diag_field    &
                (mod_name,          &
                 trim(Don_save%tracername(nn)) // '_qtmes1_col',  &
                 axes(1:2), Time,   &
                 'column integrated ' //trim(Don_save%tracername(nn)) //&
                 ' mesoscale tendency',&
                trim(Don_save%tracer_units(nn)) // '* kg/(m**2 s) ', &
                 missing_value=missing_value)

!    column-integrated tracer tendency due to mesoscale redistribution:
          id_wtp1_col(nn) = register_diag_field     &
                (mod_name,  &
                 trim(Don_save%tracername(nn)) // '_wtp1_col',   &
                 axes(1:2), Time,  &
                 'column integrated '//trim(Don_save%tracername(nn)) // &
                 ' mesoscale redist',&
               trim(Don_save%tracer_units(nn)) // '* kg/(m**2 s) ', &
                missing_value=missing_value)

!    column-integrated tracer tendency due to deep convective wet 
!    deposition: 
          id_total_wet_dep_col(nn) = register_diag_field     &
                 (mod_name,  &
                  trim(Don_save%tracername(nn)) // '_totwdep_col',   &
                  axes(1:2), Time,  &
                'column integrated '//trim(Don_save%tracername(nn)) // &
                ' deep convective wet depo',&
                 trim(Don_save%tracer_units(nn)) // '* kg/(m**2 s) ', &
                 missing_value=missing_value)

!    column-integrated tracer tendency due to mesocscale updraft  wet 
!    deposition: 
          id_meso_wet_dep_col(nn) = register_diag_field     &
                (mod_name,  &
                  trim(Don_save%tracername(nn)) // '_mwdep_col',   &
                  axes(1:2), Time,  &
                'column integrated '//trim(Don_save%tracername(nn)) // &
                 ' meso updraft wet depo',&
                  trim(Don_save%tracer_units(nn)) // '* kg/(m**2 s) ', &
                 missing_value=missing_value)

!    column-integrated tracer tendency due to wet deposition in cells:
          id_cell_wet_dep_col(nn) = register_diag_field     &
                (mod_name,  &
                 trim(Don_save%tracername(nn)) // '_cwdep_col',   &
                  axes(1:2), Time,  &
                'column integrated '//trim(Don_save%tracername(nn)) // &
                  ' cell wet depo',&
                 trim(Don_save%tracer_units(nn)) // '* kg/(m**2 s) ', &
                  missing_value=missing_value)

!    column-integrated total tracer tendency:
          id_qtceme_col(nn) = register_diag_field     &
                (mod_name,  &
                 trim(Don_save%tracername(nn)) // '_qtceme_col',  &
                 axes(1:2), Time,  &
                 'column integrated ' //trim(Don_save%tracername(nn)) //&
                 ' total tendency ', &
                  trim(Don_save%tracer_units(nn)) // '* kg/(m**2 s) ', &
                 missing_value=missing_value)
        end do
      endif

!    mesoscale ice generalized effective size:
      id_dgeice_deep = register_diag_field    &
            (mod_name, 'dgeice_deep', axes(1:3), &
             Time, 'meso ice gen eff size ', 'micrometers',   &
             missing_value=missing_value)

!    cell ice mixing ratio:
      id_cuqi_deep = register_diag_field         &
            (mod_name, 'cuqi_deep', axes(1:3),  &
             Time, 'cell ice ', 'kg(H2O)/kg',   &
             missing_value=missing_value)

!    cell liquid mixing ratio:
      id_cuql_deep = register_diag_field     &
            (mod_name, 'cuql_deep', axes(1:3),  &
             Time, 'cell liquid ', 'kg(H2O)/kg',   &
             missing_value=missing_value)

!    cell liquid generalized effective size:
      id_dgeliq_deep = register_diag_field    &
            (mod_name, 'dgeliq_deep', axes(1:3), &
             Time, 'cell liq gen eff size ', 'micrometers',   &
             missing_value=missing_value)

!    pressure at lifting condensation level:
      id_plcl_deep = register_diag_field       &
            (mod_name, 'plcl_deep', axes(1:2),   &
             Time, 'pressure at lcl ', 'Pa ',   &
             missing_value=missing_value)

!    pressure at level of free convection:
      id_plfc_deep = register_diag_field     &
            (mod_name, 'plfc_deep', axes(1:2),   &
             Time, 'pressure at lfc ', 'Pa ',   &
             missing_value=missing_value)

!    pressure at level of zero buoyancy:  
      id_plzb_deep = register_diag_field      &
            (mod_name, 'plzb_deep', axes(1:2),   &
             Time, 'pressure at lzb ', 'Pa ',   &
             missing_value=missing_value)

!    convective available potential energy (cape):
      id_xcape_deep = register_diag_field      &
            (mod_name, 'xcape_deep', axes(1:2),  &
             Time, 'cape', 'J/kg',   &
             missing_value=missing_value)

!    convective inhibition:
      id_coin_deep = register_diag_field      &
            (mod_name, 'coin_deep', axes(1:2),   &
             Time, 'convective inhibition ', 'J/kg',   &
             missing_value=missing_value)

!    time tendency of cape:
      id_dcape_deep = register_diag_field      &
            (mod_name, 'dcape_deep', axes(1:2), &
             Time, 'time tendency of cape ', 'J/kg/sec',   &
             missing_value=missing_value)

!    column integrated water vapor:
      id_qint_deep = register_diag_field    &
            (mod_name, 'qint_deep', axes(1:2),   &
             Time, 'column moisture ', 'kg(h2o)/m**2',   &
             missing_value=missing_value)

!    fractional area of cumulus ensemble member:
      id_a1_deep = register_diag_field           &
            (mod_name, 'a1_deep', axes(1:2),   &
             Time, 'fractional area of cu subensemble ', 'percent',   &
             missing_value=missing_value)

!    fractional area of largest cumulus ensemble member:
      id_amax_deep = register_diag_field      &
            (mod_name, 'amax_deep', axes(1:2),   &
             Time, 'fractional area of largest cu subensemble ',  &
             'percent',  missing_value=missing_value)

!    upper limit onfractional area based on moisture constraint:
      id_amos_deep = register_diag_field      &
            (mod_name, 'amos_deep', axes(1:2),   &
             Time, 'uppr lmt on frac area from moisture', 'percent', &
             missing_value=missing_value)

!    area-weighted total precipitation:
      id_tprea1_deep = register_diag_field         &
            (mod_name, 'tprea1_deep', axes(1:2), &
             Time, 'area wtd total precip ', 'mm/day',   &
             missing_value=missing_value)

!    mesoscale cloud fraction:
      id_ampta1_deep = register_diag_field       &
            (mod_name, 'ampta1_deep', axes(1:2), &
             Time, 'meso cld frac', 'percent',   &
             missing_value=missing_value)

!    accumulated low-level vertical displacement:
      id_omint_deep = register_diag_field      &
            (mod_name, 'omint_deep', axes(1:2), &
             Time, 'accumulated low-lvl displ', 'Pa ',   &
             missing_value=missing_value)

!    area-weighted convective precipitation:
      id_rcoa1_deep = register_diag_field     &
            (mod_name, 'rcoa1_deep', axes(1:2),  &
             Time, 'area wtd cnvctv precip ', 'mm/day',   &
             missing_value=missing_value)

!----------------------------------------------------------------------

    if (do_ensemble_diagnostics) then
!
      allocate ( id_cpre_cem(ncem))
      allocate ( id_pb_cem(ncem))
      allocate ( id_ptma_cem(ncem))
      allocate ( id_h1_cem(ncem))
      allocate ( id_qlw_cem(ncem))
      allocate ( id_cfi_cem(ncem))
      allocate ( id_wv_cem(ncem))
      allocate ( id_rcl_cem(ncem))
!  Donner cumulus ensemble member diagnostics
!
!    GCM model pressure field on full levels:
      id_pfull_cem = register_diag_field  &
            (mod_name, 'p_full', axes(1:3), &
             Time, 'GCM model pressure on full levels (lo-res)', 'Pa', &
             missing_value=missing_value)

!    GCM model pressure field on half levels:
      id_phalf_cem = register_diag_field  &
            (mod_name, 'p_half', axes(half), &
             Time, 'GCM model pressure on half levels (lo-res)', 'Pa', &
             missing_value=missing_value)

!    GCM model height field on full levels:
      id_zfull_cem = register_diag_field  &
            (mod_name, 'z_full', axes(1:3), &
             Time, 'GCM model height on full levels (lo-res)', 'm', &
             missing_value=missing_value)

!    GCM model height field on half levels:
      id_zhalf_cem = register_diag_field  &
            (mod_name, 'z_half', axes(half), &
             Time, 'GCM model height on half levels (lo-res)', 'm', &
             missing_value=missing_value)

!    GCM model temperature field on full levels:
      id_temp_cem = register_diag_field  &
            (mod_name, 'temp', axes(1:3), &
             Time, 'GCM model temperature on full levels (lo-res)', 'K', &
             missing_value=missing_value)

!    GCM model mixing ratio field on full levels:
      id_mixing_ratio_cem = register_diag_field  &
            (mod_name, 'mixing_ratio', axes(1:3), &
             Time, 'GCM model mixing ratio on full levels (lo-res)', &
             'kg(h2o)/kg(dry air)', &
             missing_value=missing_value)

      do nn=1,ncem

        if( nn <= 9 )then
          write( chvers, '(i1)' ) nn
        else if( nn <= 99 )then
          write( chvers, '(i2)' ) nn
        else
          print *, 'Error in subroutine register_fields:'
          print *, '  number of specified cumulus ensemble members = ',ncem
          print *, '  is more than current limit of 99.'
!          stop
        call error_mesg ('fms_donner_mod', 'register_fields: & 
         &Error in subroutine register_fields : number of specified &
         &cumulus ensemble members is more than current limit of 99.',&
                                                                  FATAL) 
        endif

!    area-weighted convective precipitation rate:
        id_cpre_cem(nn) = register_diag_field  &
            (mod_name, 'cpre_cem'//TRIM(chvers), axes(1:2), &
             Time, 'area wtd cnvctv precip rate - member '//TRIM(chvers), &
             'mm/day', &
             missing_value=missing_value)

!    pressure at cloud base:
        id_pb_cem(nn) = register_diag_field  &
            (mod_name, 'pb_cem'//TRIM(chvers), axes(1:2), &
             Time, 'pressure at cloud base - member '//TRIM(chvers), &
             'Pa', &
             missing_value=missing_value)

!    pressure at cloud top:
        id_ptma_cem(nn) = register_diag_field  &
            (mod_name, 'ptma_cem'//TRIM(chvers), axes(1:2), &
             Time, 'pressure at cloud top - member '//TRIM(chvers), &
             'Pa', &
             missing_value=missing_value)

!    condensation rate profile on lo-res grid:
        id_h1_cem(nn) = register_diag_field  &
            (mod_name, 'h1_cem'//TRIM(chvers), axes(1:3), &
             Time, 'condensation rate profile - member '//TRIM(chvers), &
             'kg(h2o)/(kg(dry air) sec)', &
             missing_value=missing_value)
        
! IF LOOP HERE:
       if (.not. do_donner_plume) then
!    cloud water profile on lo-res grid:
        id_qlw_cem(nn) = register_diag_field  &
            (mod_name, 'qlw_cem'//TRIM(chvers), axes(1:3), &
             Time, 'cloud water profile - member '//TRIM(chvers), &
            'kg(h2o)/kg(air)', &
             missing_value=missing_value)

!    fraction of condensate that is ice on lo-res grid:
        id_cfi_cem(nn) = register_diag_field  &
            (mod_name, 'cfi_cem'//TRIM(chvers), axes(1:3), &
             Time, 'condensate ice fraction - member '//TRIM(chvers), &
             'fraction', &
             missing_value=missing_value)

!    vertical velocity profile in plume on lo-res grid:
        id_wv_cem(nn) = register_diag_field  &
            (mod_name, 'wv_cem'//TRIM(chvers), axes(1:3), &
             Time, 'plume vertical velocity - member '//TRIM(chvers), &
             'm / s', &
             missing_value=missing_value)

!    cloud radius profile in plume on lo-res grid:
        id_rcl_cem(nn) = register_diag_field  &
            (mod_name, 'rcl_cem'//TRIM(chvers), axes(1:3), &
             Time, 'plume cloud radius - member '//TRIM(chvers), &
             'm', &
             missing_value=missing_value)

        else
!    cloud water profile on hi-res grid:
        id_qlw_cem(nn) = register_diag_field  &
            (mod_name, 'qlw_cem'//TRIM(chvers), donner_axes(cldindices), &
             Time, 'cloud water profile - member '//TRIM(chvers), &
            'kg(h2o)/kg(air)', &
             missing_value=missing_value)

!    fraction of condensate that is ice on hi-res grid:
        id_cfi_cem(nn) = register_diag_field  &
            (mod_name, 'cfi_cem'//TRIM(chvers), donner_axes(cldindices), &
             Time, 'condensate ice fraction - member '//TRIM(chvers), &
             'fraction', &
             missing_value=missing_value)

!    vertical velocity profile in plume on hi-res grid:
        id_wv_cem(nn) = register_diag_field  &
            (mod_name, 'wv_cem'//TRIM(chvers), donner_axes(cldindices), &
             Time, 'plume vertical velocity - member '//TRIM(chvers), &
             'm / s', &
             missing_value=missing_value)

!    cloud radius profile in plume on hi-res grid:
        id_rcl_cem(nn) = register_diag_field  &
            (mod_name, 'rcl_cem'//TRIM(chvers), donner_axes(cldindices), &
             Time, 'plume cloud radius - member '//TRIM(chvers), &
             'm', &
             missing_value=missing_value)

        endif
      enddo

!    area-weighted mesoscale precipitation rate:
        id_mpre_cem = register_diag_field  &
            (mod_name, 'mpre_cem', axes(1:2), &
             Time, 'area wtd mesoscale precip rate ', &
             'mm/day', &
             missing_value=missing_value)

!    fractional area sum:
      id_a1_cem = register_diag_field  &
            (mod_name, 'a1_cem', axes(1:2), &
             Time, 'fractional area sum', 'fraction', &
             missing_value=missing_value)

!    cloud fraction, cells+meso, normalized by a(1,p_b) on lo-res grid:
      id_cual_cem = register_diag_field  &
            (mod_name, 'cual_cem', axes(1:3), &
             Time, 'cloud fraction, cells+meso, normalized by a(1,p_b)', &
             'fraction', &
             missing_value=missing_value)

!    time tendency of temperature due to deep convection on lo-res grid:
      id_tfrc_cem = register_diag_field  &
            (mod_name, 'tfrc_cem', axes(1:3), &
             Time, 'temperature tendency due to deep convection (lo-res)', &
             'K/sec', missing_value=missing_value)

    endif ! (do_ensemble_diagnostics)

end subroutine register_fields 

!#####################################################################

subroutine process_coldstart (Time, Initialized, Nml, Don_save)

!-----------------------------------------------------------------------
!    subroutine process_coldstart provides initialization that is needed
!    when the job is a donner_deep coldstart, or if the user-supplied 
!    restart file is not usable for a restart with the current code 
!    version.
!-----------------------------------------------------------------------

type(time_type), intent(in) :: Time
type(donner_initialized_type), intent(inout) :: Initialized
type(donner_save_type), intent(inout) :: Don_save
type(donner_nml_type), intent(inout) :: Nml     

!---------------------------------------------------------------------
!   intent(in) variables:
!
!        Time      current time [ time_type, secs and days ]
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      integer  :: days, secs   ! components of current time

!---------------------------------------------------------------------
!    set the coldstart flag to .true.. set the time until the first cal-
!    culation call to donner_deep_mod, donner_deep calculation calls will
!    be every donner_deep_freq seconds after the start of the day.
!---------------------------------------------------------------------
      Initialized%coldstart = .true.
      call get_time (Time, secs, days)
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
! register restart field to be written to restart file.
subroutine fms_donner_register_restart(fname, Initialized, ntracers, Don_save, Nml)
  character(len=*),                 intent(in) :: fname
  type(donner_initialized_type), intent(inout) :: Initialized
  integer,                          intent(in) :: ntracers
  type(donner_save_type),        intent(inout) :: Don_save
  type(donner_nml_type),         intent(inout) :: Nml
  character(len=64)                            :: fname2
  integer :: id_restart, n

   call get_mosaic_tile_file(fname, fname2, .false. ) 
   allocate(Don_restart)
   if(trim(fname2) == trim(fname)) then
      Til_restart => Don_restart
      in_different_file = .false.
   else
      in_different_file = .true.
      allocate(Til_restart)
   endif

   id_restart = register_restart_field(Don_restart, fname, 'conv_alarm', Initialized%conv_alarm, no_domain = .true.)
   id_restart = register_restart_field(Don_restart, fname, 'donner_deep_freq', Nml%donner_deep_freq, no_domain = .true.)

   if (.not. (write_reduced_restart_file) .or. &
        Initialized%conv_alarm >  Initialized%physics_dt)  then
      id_restart = register_restart_field(Til_restart, fname, 'cemetf', Don_save%cemetf)
      id_restart = register_restart_field(Til_restart, fname, 'cememf', Don_save%cememf)
      id_restart = register_restart_field(Til_restart, fname, 'mass_flux', Don_save%mass_flux)
      id_restart = register_restart_field(Til_restart, fname, 'cell_up_mass_flux', Don_save%cell_up_mass_flux)
      id_restart = register_restart_field(Til_restart, fname, 'det_mass_flux', Don_save%det_mass_flux)
      id_restart = register_restart_field(Til_restart, fname, 'dql_strat', Don_save%dql_strat)
      id_restart = register_restart_field(Til_restart, fname, 'dqi_strat', Don_save%dqi_strat)
      id_restart = register_restart_field(Til_restart, fname, 'dqa_strat', Don_save%dqa_strat)
      id_restart = register_restart_field(Til_restart, fname, 'tprea1', Don_save%tprea1)
      id_restart = register_restart_field(Til_restart, fname, 'humidity_area', Don_save%humidity_area)
      id_restart = register_restart_field(Til_restart, fname, 'humidity_factor', Don_save%humidity_factor)
      if (Initialized%do_donner_tracer) then
         do n=1,ntracers
            id_restart = register_restart_field(Til_restart, fname, 'tracer_tends_'// trim(Don_save%tracername(n)), &
                 Don_save%tracer_tends(:,:,:,n))
         end do
      endif
   endif
   id_restart = register_restart_field(Til_restart, fname, 'parcel_disp', Don_save%parcel_disp)
   id_restart = register_restart_field(Til_restart, fname, 'lag_temp', Don_save%lag_temp)
   id_restart = register_restart_field(Til_restart, fname, 'lag_vapor', Don_save%lag_vapor)
   id_restart = register_restart_field(Til_restart, fname, 'lag_press', Don_save%lag_press)

end subroutine fms_donner_register_restart


!#####################################################################
! <SUBROUTINE NAME="read_restart_nc">
!  <OVERVIEW>
!    read_restart_nc reads a netcdf restart file containing donner_deep
!    restart information.
!  </OVERVIEW>
!  <DESCRIPTION>
!    read_restart_nc reads a netcdf restart file containing donner_deep
!    restart information.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call read_restart_nc
!  </TEMPLATE>
! </SUBROUTINE>
!


subroutine read_restart_nc (ntracers, Initialized, Nml, Don_save)

!-----------------------------------------------------------------------
!    subroutine read_restart_nc reads a netcdf restart file to obtain 
!    the variables needed upon experiment restart. 
!-----------------------------------------------------------------------

integer, intent(in) :: ntracers
type(donner_initialized_type), intent(inout) :: Initialized
type(donner_save_type), intent(inout) :: Don_save
type(donner_nml_type), intent(inout) :: Nml     

!----------------------------------------------------------------------
!   intent(in) variables:
!
!      ntracers    number of tracers being transported by the
!                  donner deep convection parameterization in this job
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      logical,         dimension(ntracers)  :: success
      integer,         dimension(:), allocatable :: ntindices
      type(fieldtype), dimension(:), allocatable :: tracer_fields

      character(len=64)     :: fname2='INPUT/donner_deep.res.tile1'
      character(len=64)     :: fname='INPUT/donner_deep.res.nc'
      character(len=128)    :: tname
      integer               :: ndim, natt, nvar, ntime
      integer               :: old_freq
      integer               :: n_alltracers, iuic
      logical               :: is_tracer_in_restart_file
      integer, dimension(4) :: siz
      logical               :: field_found, field_found2, &
                               field_found4
      integer               :: it, jn, nn

!---------------------------------------------------------------------
!   local variables:
!
!        success          logical indicating if needed data for tracer n 
!                         was obtained from restart file
!        ntindices        array of all tracer indices
!        tracer_fields    field_type variable containing information on
!                         all restart file variables
!        fname2           restart file name without ".nc" appended, 
!                         needed as argument in call to mpp_open
!        fname            restart file name
!        tname            contains successive variable names from 
!                         restart file
!        ndim             number of dimensions in restart file
!        natt             number of attributes in restart file
!        nvar             number of variables in restart file
!        ntime            number of time levels in restart file
!        old_freq         donner_deep_freq as read from restart file;
!                         value used during previous job
!        n_alltracers     number of tracers registered with 
!                         tracer_manager_mod
!        iuic             unit number assigned to restart file
!        is_tracer_in_restart_file  
!                         should we stop searching the restart file 
!                         for the current tracer name because it has 
!                         been found ?
!        siz              sizes (each dimension) of netcdf variable 
!        field_found      is the requested variable in the restart file ?
!                         if it is not, then this is a reduced restart
!                         file
!        field_found2     is the requested variable in the restart file ?
!                         if it is not, then Don_save%det_mass_flux and
!                         Don_save%cell_up_mass_flux must be initialized
!        it, jn, nn       do-loop indices
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!    output a message indicating entrance into this routine.
!--------------------------------------------------------------------
      if (mpp_pe() == mpp_root_pe() ) then
        call error_mesg ('donner_deep_mod',  'read_restart_nc:&
             &Reading netCDF formatted restart file: &
                                 &INPUT/donner_deep.res.nc', NOTE)
      endif

!-------------------------------------------------------------------
!    read the values of conv_alarm when the restart file was written and
!    the frequency of calculating donner deep convection effects in the
!    job which wrote the file.
!-------------------------------------------------------------------
      call read_data(fname, 'conv_alarm', Initialized%conv_alarm,   &
                                                       no_domain=.true.)
      call read_data(fname, 'donner_deep_freq', old_freq,   &
                                                       no_domain=.true.)
  
!----------------------------------------------------------------------
!    call field_size to determine if variable cemetf is present in the
!    restart file.
!----------------------------------------------------------------------
      call field_size(fname, 'cemetf', siz, field_found=field_found)

!---------------------------------------------------------------------
!    if the frequency of calculating deep convection has changed, 
!    redefine the time remaining until the next calculation.
!---------------------------------------------------------------------
      if (Nml%donner_deep_freq /= old_freq) then
        Initialized%conv_alarm = Initialized%conv_alarm - old_freq +  &
                                 Nml%donner_deep_freq
        if (mpp_pe() == mpp_root_pe()) then
          call error_mesg ('donner_deep_mod', 'read_restart_nc:  &
                   &donner_deep time step has changed', NOTE)
        endif

!----------------------------------------------------------------------
!    if cemetf is not present, then this is a reduced restart file. it 
!    is not safe to change the frequency of calculating donner 
!    effects when reading a reduced restart file, so a fatal error is
!    generated.
!----------------------------------------------------------------------
        if (.not. field_found) then
          call error_mesg ('donner_deep_mod', 'read_restart_nc: &
           & cannot use reduced restart file and change donner_deep_freq&
           & within experiment and guarantee restart reproducibility', &
                                                                  FATAL)
        endif
      endif  !(donner_deep_freq /= old_freq)

!---------------------------------------------------------------------
!    read the restart data that is present in a full restart but absent
!    in a reduced restart.
!---------------------------------------------------------------------
      if (field_found) then
        call read_data (fname, 'cemetf',  Don_save%cemetf)
        call read_data (fname, 'cememf',  Don_save%cememf)            
        call read_data (fname, 'mass_flux', Don_save%mass_flux)
        call read_data (fname, 'dql_strat', Don_save%dql_strat)
        call read_data (fname, 'dqi_strat', Don_save%dqi_strat)
        call read_data (fname, 'dqa_strat', Don_save%dqa_strat)
        call read_data (fname, 'tprea1', Don_save%tprea1)       
        call read_data (fname, 'humidity_area', Don_save%humidity_area) 

!---------------------------------------------------------------------
!  determine if humidity_factor is in file. if it is, read the values 
!  into Don_Save%humidity_factor. if it is not (it is an older file), 
!  it is only required if donner_deep will not be called on the first 
!  step of this job.
!  if that is the case, stop with a fatal error; otherwise, continue on,
!  since humidity_factor will be calculated before it is used.
!---------------------------------------------------------------------
        call field_size(fname, 'humidity_factor', siz,   &
                                              field_found=field_found4)
        if (field_found4) then
          call read_data (fname, 'humidity_factor',  &
                                              Don_save%humidity_factor)
        else if (Initialized%conv_alarm > 0.0) then
          call error_mesg ('donner_deep_mod', &
             'cannot restart with this restart file unless donner_deep &
                &calculated on first step', FATAL)
        endif

!----------------------------------------------------------------------
!    determine if det_mass_flux is present in the file.
!----------------------------------------------------------------------
        call field_size(fname, 'det_mass_flux', siz,    &
                                               field_found=field_found2)

!----------------------------------------------------------------------
!    if it is present, then read det_mass_flux and cell_up_mass_flux.
!----------------------------------------------------------------------
        if (field_found2) then
          call read_data (fname, 'det_mass_flux', Don_save%det_mass_flux)
          call read_data (fname, 'cell_up_mass_flux',    &
                                              Don_save%cell_up_mass_flux)

!----------------------------------------------------------------------
!    if it is not present (an earlier version of this file), set 
!    det_mass_flux and cell_up_mass_flux to default values.
!----------------------------------------------------------------------
        else
          Don_save%det_mass_flux     = 0.0
          Don_save%cell_up_mass_flux = 0.0
        endif

!------------------------------------------------------------------
!    if tracers are to be transported, see if tendencies are available
!    in the restart file.
!------------------------------------------------------------------
        if (Initialized%do_donner_tracer) then

!---------------------------------------------------------------------
!    initialize a logical array indicating whether the data for each
!    tracer is available.
!---------------------------------------------------------------------
          success = .false.

!---------------------------------------------------------------------
!    open the restart file with mpp_open so that the unit number is 
!    available. obtain needed file characteristics by calling 
!    mpp_read_meta and  mpp_get_info. 
!---------------------------------------------------------------------
          call mpp_open(iuic, fname2, &
               action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_SINGLE )
          call mpp_read_meta (iuic)
          call mpp_get_info (iuic, ndim, nvar, natt, ntime)

!---------------------------------------------------------------------
!    obtain information on the file variables by calling mpp_get_fields.
!    it is returned in a field_type variable tracer_fields; the specific
!    information needed is the variable name.
!---------------------------------------------------------------------
          allocate (tracer_fields(nvar))
          if (mpp_pe() == mpp_root_pe()) then
            call mpp_get_fields (iuic, tracer_fields)
          endif

!---------------------------------------------------------------------
!    call get_number_tracers to determine how many tracers are registered
!    with tracer manager. allocate an array to hold their tracer indices.
!    call get_tracer_indices to retrieve the tracer indices. 
!---------------------------------------------------------------------
          call get_number_tracers (MODEL_ATMOS, num_tracers=n_alltracers)
          allocate (ntindices(n_alltracers))
          call get_tracer_indices (MODEL_ATMOS, ind=ntindices)

!----------------------------------------------------------------------
!    loop over the tracers, obtaining their names via a call to
!    get_tracer_names. bypass those tracers known to not be transported
!    by donner convection.
!----------------------------------------------------------------------
          do it=1,n_alltracers
            call get_tracer_names (MODEL_ATMOS, ntindices(it), tname)
            if (tname == "sphum"  ) cycle
            if (tname == "liq_wat") cycle
            if (tname == "ice_wat") cycle
            if (tname == "cld_amt") cycle

!--------------------------------------------------------------------
!    initialize a logical indicating whether this tracer is in the 
!    restart file.
!--------------------------------------------------------------------
            is_tracer_in_restart_file = .FALSE.

!---------------------------------------------------------------------
!    loop over the variables in the restart file to determine if the
!    current tracer's time tendency field is present.
!---------------------------------------------------------------------
            do jn=1,nvar 
              if (lowercase (trim(mpp_get_field_name(tracer_fields(jn)))) ==   &
                  lowercase ('tracer_tends_' // trim(tname)) ) then 

!---------------------------------------------------------------------
!    if tracer tendency is in restart file, write a message. set the 
!    logical flag indicating such to .true..
!---------------------------------------------------------------------
                if (mpp_pe() == mpp_root_pe() )  then
                  print *,'tracer_tends_' // trim(tname), ' found!'
                endif
                is_tracer_in_restart_file = .TRUE.

!---------------------------------------------------------------------
!    loop over the tracers being transported by donner convection in this
!    job to determine if this tracer is one of those being transported.
!    determine the tracer index in tracername array corresponding to 
!    this tracer.
!---------------------------------------------------------------------
                do nn=1,ntracers
                  if (lowercase( 'tracer_tends_' // trim(tname) ) == &
                      'tracer_tends_' // Don_save%tracername(nn) )  then
                  
!---------------------------------------------------------------------
!    if data for this tracer is needed, read data into proper section of
!    array tracer_tends. set the logical flag for this tracer indicating 
!    successful retrieval. exit this loop.
!---------------------------------------------------------------------
                    call read_data (fname,   &
                                  'tracer_tends_' // trim(tname),   &
                                   Don_save%tracer_tends(:,:,:,nn))
                    success(nn) = .true.
                    exit
                  endif 
                end do  ! (nn)
              endif

!---------------------------------------------------------------------
!    if desired tracer has been found, stop searching the restart file
!    variables for this tracer and cycle to begin searching the restart
!    file for the next field_table tracer.
!---------------------------------------------------------------------
              if (is_tracer_in_restart_file) exit
            end do !  (jn)
          end do ! (it)

!---------------------------------------------------------------------
!    initialize the time tendencies to 0.0 for any tracers that are to
!    be transported and whose time tendencies were not found on the 
!    restart file.  enter a message in the output file.
!---------------------------------------------------------------------
          do nn=1,ntracers
            if (success(nn) ) then
            else
              call error_mesg ('donner_deep_mod', 'read_restart_nc: &
                  &did not find tracer restart data for ' //  &
                  trim(Don_save%tracername(nn)) //  &
                  '; am initializing tendency to 0.0', NOTE)
              Don_save%tracer_tends(:,:,:,nn) = 0.0
            endif   
          end do

!----------------------------------------------------------------------
!    deallocate local variables.
!----------------------------------------------------------------------
          deallocate (ntindices)
          deallocate (tracer_fields)
        endif  ! (do_donner_tracer)
      endif  ! (field_found)

!---------------------------------------------------------------------
!    read the restart data that is present in both full and reduced
!    restart files.
!---------------------------------------------------------------------
      call read_data (fname, 'parcel_disp', Don_save%parcel_disp)
      call read_data (fname, 'lag_temp',    Don_save%lag_temp)     
      call read_data (fname, 'lag_vapor',   Don_save%lag_vapor)     
      call read_data (fname, 'lag_press',   Don_save%lag_press)     

!---------------------------------------------------------------------




end subroutine read_restart_nc



!#####################################################################

subroutine process_monitors (idf, jdf, nlev, ntracers, axes, Time, &
                              Initialized, Don_save)

integer,                       intent(in)  :: idf, jdf, nlev, ntracers
integer,         dimension(4), intent(in)  :: axes
type(time_type),               intent(in)  :: Time
type(donner_initialized_type), intent(inout) :: Initialized
type(donner_save_type), intent(inout) :: Don_save

!-------------------------------------------------------------------
!  local variables:

      integer             :: n, nx, nc
      logical             :: flag, success
      integer             :: nfields, model, num_methods
      character(len=200)  :: method_name, field_type, method_control,&
                             field_name, list_name
      character(len=32)   :: path_name = '/atmos_mod/don_deep_monitor/'

!---------------------------------------------------------------------
!    determine if and how many output variables are to be monitored. 
!    set a flag indicating if monitoring is activated.
!---------------------------------------------------------------------
      call field_manager_init (nfields)
      nx = 0
      do n=1,nfields
        call get_field_info (n, field_type, field_name, model, &
                             num_methods)
        if (trim(field_type) == 'don_deep_monitor') then
          nx = nx + 1
        endif
      end do
      if (nx > 0) then
        Initialized%monitor_output = .true.
      else
        Initialized%monitor_output = .false.
      endif

!---------------------------------------------------------------------
!    allocate arrays needed for each monitored variable. 
!---------------------------------------------------------------------
      if (Initialized%monitor_output) then
        allocate (Initialized%Don_monitor(nx))
        allocate (id_extremes(nx))
        allocate (id_hits(nx))

!---------------------------------------------------------------------
!    read the field_table to determine the nature of the monitors
!    requested.
!---------------------------------------------------------------------
        nx = 1
        do n = 1,nfields
          call get_field_info (n, field_type, field_name, model, &
                               num_methods)

!---------------------------------------------------------------------
!    define the list name used by field_manager_mod to point to 
!    monitored variables.
!---------------------------------------------------------------------
          if (trim(field_type) == 'don_deep_monitor') then
            list_name = trim(path_name) // trim(field_name) // '/'

!--------------------------------------------------------------------
!    place name of field in don_monitor_type variable.
!--------------------------------------------------------------------
            Initialized%Don_monitor(nx)%name = trim(field_name)

!--------------------------------------------------------------------
!    map the field name to the list of acceptable field names. store
!    the index of this field name in the don_monitor_type variable.
!    note that any tracer variables need to have 'tr_' as the first
!    three characters in their name to allow proper processing. store
!    the appropriate tracer index for any tracer arrays.
!--------------------------------------------------------------------
            if (trim(field_name(1:3)) == 'tr_') then
              select case (trim(field_name(4:9)))
                case ('rn_ten')
                  Initialized%Don_monitor(nx)%index = RADON_TEND
                  success = .false.
                  do nc=1,ntracers
                    if (trim(Don_save%tracername(nc)) == 'radon') then
                      Initialized%Don_monitor(nx)%tracer_index = nc
                      success = .true.
                      exit
                    endif
                  end do
                  if (.not. success) then
                    call error_mesg ('donner_deep_mod', &
                     'not able to find "radon" tracer index', FATAL)
                  endif
                case default
                  call error_mesg ('donner_deep_mod', &
                 'tracer variable name in field_table don_deep_monitor &
                                             &type is invalid', FATAL)
              end select

!---------------------------------------------------------------------
!    for non-tracer variables, set the tracer index to an arbitrary 
!    value.
!---------------------------------------------------------------------
            else
              Initialized%Don_monitor(nx)%tracer_index = 0
              select case (trim(field_name(1:6)))
                case ('det_ma')
                  Initialized%Don_monitor(nx)%index = DET_MASS_FLUX
                case ('mass_f')
                  Initialized%Don_monitor(nx)%index = MASS_FLUX
                case ('cell_u')
                  Initialized%Don_monitor(nx)%index =   &
                                                  CELL_UPWARD_MASS_FLUX
                case ('temp_f')
                  Initialized%Don_monitor(nx)%index = TEMP_FORCING
                case ('moistu')
                  Initialized%Don_monitor(nx)%index = MOIST_FORCING
                case ('precip')
                  Initialized%Don_monitor(nx)%index = PRECIP
                case ('freeze')
                  Initialized%Don_monitor(nx)%index = FREEZING
                case default
                  call error_mesg ('donner_deep_mod', &
                      'variable name in field_table don_deep_monitor &
                                              &type is invalid', FATAL)
              end select
            endif

!---------------------------------------------------------------------
!    read the units for this variable from the field_table entry.
!    if the units method is missing, set units to be 'missing'.
!---------------------------------------------------------------------
            flag = fm_query_method (trim(list_name) //  'units',    &
                                    method_name, method_control)
            if (flag) then
              Initialized%Don_monitor(nx)%units = trim(method_name)
            else
              Initialized%Don_monitor(nx)%units = 'missing'
            endif

!---------------------------------------------------------------------
!    determine the type of limit being imposed for this variable from 
!    the field_table entry.
!---------------------------------------------------------------------
            flag = fm_query_method (trim(list_name) // 'limit_type',  &
                                    method_name, method_control)

!----------------------------------------------------------------------
!    include the limit_type for this variable in its don_monitor type
!    variable.
!    register diagnostics associated with the monitored output fields
!    (extreme values and number of times threshold was exceeeded).
!----------------------------------------------------------------------
            if ( flag) then
              if (trim(method_name) == 'maxmag') then
                Initialized%Don_monitor(nx)%initial_value = 0.0
                Initialized%Don_monitor(nx)%limit_type =   MAXMAG
                id_extremes(nx) = register_diag_field (mod_name,   &
                  'maxmag_'// trim(Initialized%Don_monitor(nx)%name),  &
                   axes(1:3),  Time,  'maxmag values of ' // &
                            trim(Initialized%Don_monitor(nx)%name),  &
                   Initialized%Don_monitor(nx)%units,   &
                   mask_variant = .true., missing_value=missing_value)
                id_hits(nx) = register_diag_field (mod_name,   &
                  'num_maxmag_'// &
                              trim(Initialized%Don_monitor(nx)%name) , &
                   axes(1:3),  Time,    &
                   '# of times that magnitude of '&
                     // trim(Initialized%Don_monitor(nx)%name) //  &
                 ' > ' // trim(method_control(2:)) // ' ' // &
                    trim(Initialized%Don_monitor(nx)%units) ,  &
                   'number', mask_variant = .true., & 
                   missing_value=missing_value)
              else if (trim(method_name) == 'minmag') then
                Initialized%Don_monitor(nx)%initial_value = 1.0e30
                Initialized%Don_monitor(nx)%limit_type =   MINMAG
                id_extremes(nx) = register_diag_field (mod_name,   &
                  'minmag_'// trim(Initialized%Don_monitor(nx)%name),  &
                   axes(1:3),  Time,  'minmag values of ' // &
                            trim(Initialized%Don_monitor(nx)%name),  &
                   Initialized%Don_monitor(nx)%units,   &
                   mask_variant = .true., missing_value=missing_value)
                id_hits(nx) = register_diag_field (mod_name,   &
                  'num_minmag_'//     &
                             trim(Initialized%Don_monitor(nx)%name) , &
                  axes(1:3),  Time,    &
                   '# of times that magnitude of '&
                     // trim(Initialized%Don_monitor(nx)%name) //  &
                ' < ' // trim(method_control(2:)) // ' ' // &
                    trim(Initialized%Don_monitor(nx)%units) ,  &
                  'number', mask_variant = .true., & 
                  missing_value=missing_value)
              else if (trim(method_name) == 'minval') then
                Initialized%Don_monitor(nx)%initial_value = 1.0e30
                Initialized%Don_monitor(nx)%limit_type =   MINVAL
                id_extremes(nx) = register_diag_field (mod_name,   &
                  'minval_'// trim(Initialized%Don_monitor(nx)%name),  &
                   axes(1:3),  Time,  'minimum values of ' // &
                            trim(Initialized%Don_monitor(nx)%name),  &
                  Initialized%Don_monitor(nx)%units,   &
                  mask_variant = .true., missing_value=missing_value)
                id_hits(nx) = register_diag_field (mod_name,   &
                  'num_minval_'//   &
                             trim(Initialized%Don_monitor(nx)%name) , &
                   axes(1:3),  Time,    &
                   '# of times that value of '&
                     // trim(Initialized%Don_monitor(nx)%name) //  &
                ' < ' // trim(method_control(2:)) // ' ' // &
                    trim(Initialized%Don_monitor(nx)%units) ,  &
                  'number', mask_variant = .true., & 
                  missing_value=missing_value)
              else if (trim(method_name) == 'maxval') then
                Initialized%Don_monitor(nx)%initial_value = -1.0e30
                Initialized%Don_monitor(nx)%limit_type = MAXVAL 
                id_extremes(nx) = register_diag_field (mod_name,   &
                  'maxval_'// trim(Initialized%Don_monitor(nx)%name),  &
                  axes(1:3),  Time,  'maximum values of ' // &
                            trim(Initialized%Don_monitor(nx)%name),  &
                  Initialized%Don_monitor(nx)%units,  &
                  mask_variant = .true., missing_value=missing_value)
                id_hits(nx) = register_diag_field (mod_name,   &
                  'num_maxval_'//    &
                             trim(Initialized%Don_monitor(nx)%name) , &
                  axes(1:3),  Time,    &
                   '# of times that value of '&
                     // trim(Initialized%Don_monitor(nx)%name) //  &
                    ' > ' // trim(method_control(2:)) // ' ' // &
                    trim(Initialized%Don_monitor(nx)%units) ,  &
                  'number', mask_variant = .true., & 
                  missing_value=missing_value)
              else
                call error_mesg ('donner_deep_mod', &
                    'invalid limit_type for monitored variable', FATAL)
              endif

!----------------------------------------------------------------------
!    if limit_type not in field_table, set it to look for maximum
!    magnitude.
!----------------------------------------------------------------------
            else
              Initialized%Don_monitor(nx)%initial_value = 0.0
              Initialized%Don_monitor(nx)%limit_type =   MAXMAG
              id_extremes(nx) = register_diag_field (mod_name,   &
                  'maxmag_'// trim(Initialized%Don_monitor(nx)%name),  &
                   axes(1:3),  Time,  'maxmag values of ' // &
                            trim(Initialized%Don_monitor(nx)%name),  &
                   Initialized%Don_monitor(nx)%units,   &
                   mask_variant = .true., missing_value=missing_value)
              id_hits(nx) = register_diag_field (mod_name,   &
                  'num_maxmag_'// &
                              trim(Initialized%Don_monitor(nx)%name) , &
                   axes(1:3),  Time,    &
                   '# of times that magnitude of '&
                     // trim(Initialized%Don_monitor(nx)%name) //  &
                ' > ' // trim(method_control(2:)) // ' ' // &
                    trim(Initialized%Don_monitor(nx)%units) ,  &
                   'number', mask_variant = .true., & 
                   missing_value=missing_value)
            endif

!----------------------------------------------------------------------
!    obtain the magnitude of the limit being monitored for this 
!    variable from the field_table. 
!----------------------------------------------------------------------
            flag = parse (method_control, 'value',   &
                            Initialized%Don_monitor(nx)%threshold ) > 0

!----------------------------------------------------------------------
!    if no limit_type and / or value has been given, the
!    field will be flagged for magnitudes  > 0.0, i.e., if deep 
!    convection has affected the point.
!----------------------------------------------------------------------
            if ( .not. flag) then
              Initialized%Don_monitor(nx)%threshold = 0.0
            endif

!-------------------------------------------------------------------
!    allocate and initialize arrays to hold the extrema and a count of 
!    times the threshold was exceeded at each point.
!-------------------------------------------------------------------
            allocate (Initialized%Don_monitor(nx)%extrema(idf,jdf,nlev))
            Initialized%Don_monitor(nx)%extrema(:,:,:) =  &
                         Initialized%Don_monitor(nx)%initial_value
            allocate (Initialized%Don_monitor(nx)%hits(idf,jdf,nlev))
            Initialized%Don_monitor(nx)%hits(:,:,:) = 0.0
            nx = nx + 1
          endif
        end do
      endif 

end subroutine process_monitors



!#####################################################################


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      2. ROUTINES CALLED BY DONNER_DEEP
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#######################################################################

subroutine donner_column_control (is, ie, js, je, Time, Col_diag)                

!---------------------------------------------------------------------
!    subroutine donner_column_control returns the number, location
!    (processor and window indices) and output units associated with 
!    any diagnostic columns requested within the current physics window.
!---------------------------------------------------------------------

integer,                       intent(in)   :: is, ie, js, je
type(time_type),               intent(in)   :: Time
type (donner_column_diag_type), intent(inout) :: Col_diag

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     is, ie         first and last values of i index values of points
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points
!                    in this physics window (processor coordinates)
!     Time           current model time [ time_type, days, seconds ]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer  :: isize      !   i-dimension of physics window
      integer  :: jsize      !   j-dimension of physics window
      integer  :: nn, j, i   !   do-loop indices

!--------------------------------------------------------------------
!    define the sizes of the current physics window's horizontal
!    dimensions.
!--------------------------------------------------------------------
      isize = ie - is + 1
      jsize = je - js + 1

!-------------------------------------------------------------------
!    initialize the output variables.
!-------------------------------------------------------------------
      Col_diag%i_dc(:) = -99
      Col_diag%j_dc(:) = -99
      Col_diag%unit_dc(:) = -1
      Col_diag%jgl_dc(:) = -99
      Col_diag%igl_dc(:) = -99
      Col_diag%ncols_in_window = 0

!--------------------------------------------------------------------
!    if any requested diagnostic columns are present within the current
!    physics window, and if it is at or past the time to start output-
!    ting column diagnostics, save the relevant variables describing
!    those diagnostic columns in arrays to be returned to the calling
!    routine. call column_diagnostics_header to write the file header
!    for the diagnostic columns in this window. 
!--------------------------------------------------------------------
      if (Col_diag%num_diag_pts > 0) then
        if (Time >= Time_col_diagnostics) then
          do nn=1,Col_diag%num_diag_pts
            do j=1,jsize      
              if (js + j - 1 == col_diag_j(nn)) then
                do i=1,isize       
                  if (is + i - 1 == col_diag_i(nn)) then
                    Col_diag%ncols_in_window =   &
                                           Col_diag%ncols_in_window + 1
                    Col_diag%i_dc(Col_diag%ncols_in_window) = i
                    Col_diag%j_dc(Col_diag%ncols_in_window) = j
                    Col_diag%igl_dc(COl_diag%ncols_in_window) =  &
                                                          col_diag_i(nn)
                    Col_diag%jgl_dc(Col_diag%ncols_in_window) =   &
                                                           col_diag_j(nn)
                    Col_diag%unit_dc(Col_diag%ncols_in_window) =   &
                                                        col_diag_unit(nn)
                    call column_diagnostics_header &
                            (mod_name, col_diag_unit(nn), Time, nn,  &
                             col_diag_lon, col_diag_lat, col_diag_i,  &
                             col_diag_j)
                  endif
                end do  ! (i loop)
              endif
            end do  ! (j loop) 
          end do  ! (num_diag_pts loop)
        endif  ! (Time >= starting time)
      endif ! (num_diag_pts > 0)

!---------------------------------------------------------------------

end subroutine donner_column_control



!######################################################################

subroutine donner_deep_netcdf (is, ie, js, je, Nml, Time,  Param, &
                               Initialized, Don_conv, Don_cape,&
                               Don_cem,parcel_rise, pmass, total_precip, &
                               Don_budgets, &
                               temperature_forcing, moisture_forcing)  

!---------------------------------------------------------------------
!    subroutine donner_deep_netcdf sends the fields requested in the
!    diag_table to diag_manager_mod so that they may be appropriately
!    processed for output.
!---------------------------------------------------------------------

integer,                intent(in) :: is, ie, js, je
type(time_type),        intent(in) :: Time
type(donner_param_type), intent(in) :: Param
type(donner_initialized_type), intent(inout) :: Initialized
type(donner_nml_type), intent(in) :: Nml   
type(donner_conv_type), intent(in) :: Don_conv
type(donner_budgets_type), intent(in) :: Don_budgets
type(donner_cape_type), intent(in) :: Don_cape
type(donner_cem_type),  intent(in) :: Don_cem
real, dimension(:,:,:), intent(in) :: pmass, temperature_forcing,&
                                      moisture_forcing
real, dimension(:,:),   intent(in) :: parcel_rise, total_precip

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     is, ie         first and last values of i index values of points 
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points 
!                    in this physics window (processor coordinates)
!     Time           current time (time_type)
!     Don_conv       donner_convection_type derived type variable con-
!                    taining diagnostics describing the nature of the 
!                    convection produced by the donner parameterization
!     Don_cape       donner_cape type derived type variable containing
!                    diagnostics related to the cape calculation assoc-
!                    iated with the donner convection parameterization
!     Don_cem        donner_cem_type derived type variable containing
!                    Donner cumulus ensemble member diagnostics
!     temperature_forcing  
!                    temperature tendency due to donner convection
!                    [ deg K / sec ]
!     moisture_forcing  
!                    vapor mixing ratio tendency due to donner 
!                    convection [ kg(h2o) / (kg(dry air) sec ) ]
!     pmass          mass per unit area within the grid box
!                    [ kg (air) / (m**2) ]
!     parcel_rise    accumulated vertical displacement of a near-surface
!                    parcel as a result of the lowest model level omega 
!                    field [ Pa ]
!     total_precip   total precipitation rate produced by the
!                    donner parameterization [ mm / day ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      real, dimension (ie-is+1, je-js+1)  :: tempdiag, tempdiag2, tempdiag3  
                           ! array used to hold various data fields being
                           ! sent to diag_manager_mod
      logical :: used      ! logical indicating data has been received 
                           ! by diag_manager_mod 
      integer :: nlev      ! number of large-scale model layers
      integer :: ntr       ! number of tracers transported by the
                           ! donner deep convection parameterization
      integer :: k, n, nn  ! do-loop indices
      integer :: ncem      ! number of cumulus ensemble members in the
                           ! donner deep convection parameterization

!----------------------------------------------------------------------
!    define the number of model layers (nlev) and number of transported
!    tracers (ntr).
!----------------------------------------------------------------------
      nlev = size (pmass,3)
      ntr  = size (Don_conv%qtren1,4)

!----------------------------------------------------------------------
!    define the number of cumulus ensemble members in the
!    donner deep convection parameterization.
!----------------------------------------------------------------------
      ncem = size (Don_cem%cell_precip,3)

!---------------------------------------------------------------------
!    send the 3D convective output variables to diag_manager_mod.
!!   NOTE: effective with code mod lima_donnermod3_rsh (7-19-05) the
!!         temperature and moisture forcing fields passed to diag_manager
!!         (id_cemetf_deep, id_cememf_deep) are the total convective
!!         forcings calculated by the donner parameterization. Previous
!!         code versions run in models in which strat_cloud_mod was 
!!         activated output the forcing fields less the terms related to 
!!         the flux convergence of the large-scale condensate and the 
!!         mesoscale detrainment.
!---------------------------------------------------------------------

!   total convective temperature forcing:
      used = send_data (id_cemetf_deep, Don_conv%conv_temp_forcing,  &
                        Time, is, js, 1)

!   cell entropy flux convergence:
      used = send_data (id_ceefc_deep, Don_conv%ceefc, Time, is, js, 1)

!   cell condensation / evaporation:
      used = send_data (id_cecon_deep, Don_conv%cecon, Time, is, js, 1)

!   cell moisture flux convergence:
      used = send_data (id_cemfc_deep, Don_conv%cemfc, Time, is, js, 1)

!   total convective moistening forcing:
      used = send_data (id_cememf_deep, Don_conv%conv_moist_forcing,  &
                        Time, is, js, 1)

!   total convective moistening rate after adjustnment for negative 
!   vapor mixing ratio:
      used = send_data (id_cememf_mod_deep, Don_conv%cememf_mod,   &
                        Time, is, js, 1)

!   cell + mesoscale cloud fraction:
      used = send_data (id_cual_deep, Don_conv%cual, Time, is, js, 1)

!   heating rate due to freezing:
      used = send_data (id_fre_deep, Don_conv%fre, Time, is, js, 1)

!   heating rate due to melting:
      used = send_data (id_elt_deep, Don_conv%elt, Time, is, js, 1)

!   deposition in mesoscale updraft:
      used = send_data (id_cmus_deep, Don_conv%cmus, Time, is, js, 1)

!   evaporation in convective downdrafts:
      used = send_data (id_ecds_deep, Don_conv%ecds, Time, is, js, 1)

!   evaporation / sublimation in convective updrafts:
      used = send_data (id_eces_deep, Don_conv%eces, Time, is, js, 1)

!   sublimation in mesoscale downdrafts:
      used = send_data (id_emds_deep, Don_conv%emds, Time, is, js, 1)

!   sublimation in mesoscale updrafts:
      used = send_data (id_emes_deep, Don_conv%emes, Time, is, js, 1)

!   mesoscale moisture flux convergence:
      used = send_data (id_qmes_deep, Don_conv%mrmes, Time, is, js, 1)

!   transfer of vapor from cells to mesoscale:
      used = send_data (id_wmps_deep, Don_conv%wmps, Time, is, js, 1)

!   deposition of vapor from cells to mesoscale:
      used = send_data (id_wmms_deep, Don_conv%wmms, Time, is, js, 1)

!   mesoscale entropy flux convergence:
      used = send_data (id_tmes_deep, Don_conv%tmes, Time, is, js, 1)

!   mass flux in mesoscale downdrafts:
      used = send_data (id_dmeml_deep, Don_conv%dmeml, Time, is, js, 1)

!   mass flux in cell updrafts:
      used = send_data (id_uceml_deep, Don_conv%uceml, Time, is, js, 1)

!   detrained mass flux:
      used = send_data (id_detmfl_deep, Don_conv%detmfl, Time, is, js, 1)

!   mass flux in mesoscale updrafts:
      used = send_data (id_umeml_deep, Don_conv%umeml, Time, is, js, 1)

!   mesoscale ice mixing ratio:
      used = send_data (id_xice_deep, Don_conv%xice, Time, is, js, 1)

!   mesoscale liquid mass mixing ratio
      used = send_data (id_xliq_deep, Don_conv%xliq, Time, is, js, 1)

!   mesoscale ice generalized effective size:
      used = send_data (id_dgeice_deep, Don_conv%dgeice,      &
                        Time, is, js, 1)

!   cell ice mixing ratio:
      used = send_data (id_cuqi_deep, Don_conv%cuqi, Time, is, js, 1)

!   cell liquid mixing ratio:
      used = send_data (id_cuql_deep, Don_conv%cuql, Time, is, js, 1)

!   cell liquid generalized effective size:
      used = send_data (id_dgeliq_deep, Don_conv%cell_liquid_eff_diam, &
                        Time, is, js, 1)

     if (Nml%do_budget_analysis) then
       do n=1,Don_budgets%N_WATER_BUDGET
         if (id_water_budget(n) > 0) then
            used = send_data (id_water_budget(n), &
                              Don_budgets%water_budget(:,:,:,n), &
                              Time, is, js, 1)
         endif
       end do
       do n=1,Don_budgets%N_PRECIP_TYPES
         do nn=1,Don_budgets%N_PRECIP_PATHS
           if (id_precip_budget(nn,n) > 0) then
             used = send_data (id_precip_budget(nn,n), &
                               Don_budgets%precip_budget(:,:,:,nn,n), &
                               Time, is, js, 1)
           endif
         end do
       end do
       do n=1,Don_budgets%N_ENTHALPY_BUDGET
         if (id_enthalpy_budget(n) > 0) then
           used = send_data (id_enthalpy_budget(n),   &
                             Don_budgets%enthalpy_budget(:,:,:,n), &
                             Time, is, js, 1)
         endif
       end do
       do n=1,Don_budgets%N_WATER_BUDGET
         tempdiag(:,:) = 0.
         do k=1,nlev
           tempdiag(:,:) = tempdiag(:,:) + &
                           Don_budgets%water_budget(:,:,k,n)* &
                                                     pmass(:,:,k)/1000.
         end do
         if (id_ci_water_budget(n) > 0) then
           used = send_data (id_ci_water_budget(n), tempdiag, &
                             Time, is, js)
         endif
       end do
       tempdiag3(:,:) = 0.
       do n=1,Don_budgets%N_PRECIP_TYPES
         do nn=1,Don_budgets%N_PRECIP_PATHS
           tempdiag(:,:) = 0.
           do k=1,nlev
             tempdiag(:,:) = tempdiag(:,:) + &
                             Don_budgets%precip_budget(:,:,k,nn,n)* &
                                                           pmass(:,:,k)
           end do
           if (id_ci_precip_budget(nn,n) > 0) then
             used = send_data (id_ci_precip_budget(nn,n), tempdiag, &
                               Time, is, js)
           endif
           tempdiag3(:,:) = tempdiag3(:,:) + tempdiag(:,:)
         end do
       end do
       do n=1,Don_budgets%N_ENTHALPY_BUDGET
         tempdiag(:,:) = 0.
         do k=1,nlev
           tempdiag(:,:) = tempdiag(:,:) +  &
                           Don_budgets%enthalpy_budget(:,:,k,n)* &
                                                    pmass(:,:,k)*CP_AIR
         end do
         if (id_ci_enthalpy_budget(n) > 0) then
           used = send_data (id_ci_enthalpy_budget(n), tempdiag, &
                             Time, is, js)
         endif
       end do
           
        
       tempdiag2(:,:) = 0.
       tempdiag(:,:) = 0.
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) +  &
                         (Don_budgets%precip_budget(:,:,k,2,1) +  &
                          Don_budgets%precip_budget(:,:,k,4,1))* &
                                                 Param%hls*pmass(:,:,k)
       end do
       if (id_ci_prcp_heat_frz_cell > 0) then
         used = send_data (id_ci_prcp_heat_frz_cell, tempdiag, &
                           Time, is, js)
       endif
       tempdiag2 = tempdiag2 + tempdiag
           
       tempdiag(:,:) = 0.
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) +  &
                         (Don_budgets%precip_budget(:,:,k,1,1) +   &
                          Don_budgets%precip_budget(:,:,k,3,1) + &
                          Don_budgets%precip_budget(:,:,k,5,1))* &
                                                 Param%hlv*pmass(:,:,k)
       end do
       if (id_ci_prcp_heat_liq_cell > 0) then
         used = send_data (id_ci_prcp_heat_liq_cell, tempdiag, &
                           Time, is, js)
       endif
       tempdiag2 = tempdiag2 + tempdiag
           
       tempdiag(:,:) = 0.
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) +  &
                         (Don_budgets%precip_budget(:,:,k,2,2) + &
                          Don_budgets%precip_budget(:,:,k,4,2) + &
                          Don_budgets%precip_budget(:,:,k,2,3) + &
                          Don_budgets%precip_budget(:,:,k,4,3))* &
                                                 Param%hls*pmass(:,:,k)
       end do
       if (id_ci_prcp_heat_frz_meso > 0) then
         used = send_data (id_ci_prcp_heat_frz_meso, tempdiag, &
                           Time, is, js)
       endif
       tempdiag2 = tempdiag2 + tempdiag
           
       tempdiag(:,:) = 0.
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) +  &
                         (Don_budgets%precip_budget(:,:,k,1,2) +   &
                          Don_budgets%precip_budget(:,:,k,3,2) +  &
                          Don_budgets%precip_budget(:,:,k,5,2) + &
                          Don_budgets%precip_budget(:,:,k,1,3) +   &
                          Don_budgets%precip_budget(:,:,k,3,3) +  &
                          Don_budgets%precip_budget(:,:,k,5,3))* &
                                                  Param%hlv*pmass(:,:,k)
       end do
       if (id_ci_prcp_heat_liq_meso > 0) then
         used = send_data (id_ci_prcp_heat_liq_meso, tempdiag, &
                           Time, is, js)
       endif
       tempdiag2 = tempdiag2 + tempdiag
       if ( id_ci_prcp_heat_total > 0) then
         used = send_data (id_ci_prcp_heat_total, tempdiag2, &
                           Time, is, js)
       endif
       if (id_ci_prcp_total > 0) then
         used = send_data (id_ci_prcp_total, tempdiag3, &
                           Time, is, js)
       endif
       if ( id_leff > 0) then
         used = send_data(id_leff, tempdiag2/(tempdiag3+1.0e-40), &
                           Time, is, js)
       endif
           
     endif

!--------------------------------------------------------------------
!    send the tracer-related arrays to diag_manager_mod.
!--------------------------------------------------------------------
      do n=1,ntr    

!   tracer tendency due to cells:
        if (id_qtren1(n) > 0) then
        used = send_data (id_qtren1(n), Don_conv%qtren1(:,:,:,n), &
                          Time, is, js, 1)
        endif

!   tracer tendency due to mesoscale:
         if (id_qtmes1(n) > 0) then
        used = send_data (id_qtmes1(n), Don_conv%qtmes1(:,:,:,n),   &
                          Time, is, js, 1)
        endif

!   tracer tendency due to mesoscale redistribution:
        if (id_wtp1(n) > 0) then
        used = send_data (id_wtp1(n), Don_conv%wtp1(:,:,:,n),     &
                          Time, is, js, 1)
        endif

!   tracer tendency due to deep convective wet deposition:
       if (id_total_wet_dep(n) > 0) then
     used = send_data (id_total_wet_dep(n), Don_conv%wetdept(:,:,:,n), &
                            Time, is, js, 1)
        endif
!   tracer tendency due to wet deposition in mesoscale updrafts:
       if ( id_meso_wet_dep(n) > 0) then
     used = send_data (id_meso_wet_dep(n), Don_conv%wetdepm(:,:,:,n), &
                            Time, is, js, 1)
      endif
 
!   tracer tendency due to wet deposition in cells:
      if (id_cell_wet_dep(n) > 0) then
     used = send_data (id_cell_wet_dep(n), Don_conv%wetdepc(:,:,:,n), &
                           Time, is, js, 1)
      endif

!   total tracer tendency:
      if (id_qtceme(n) > 0) then
        used = send_data (id_qtceme(n), Don_conv%qtceme(:,:,:,n), &
                          Time, is, js, 1)
      endif

!---------------------------------------------------------------------
!    define the column-integrated tracer tendency due to convective
!    cells, in units of kg (tracer) / (m**2 sec). send it to 
!    diag_manager_mod.
!---------------------------------------------------------------------
        tempdiag = 0.0
        do k=1,nlev
          tempdiag(:,:) = tempdiag(:,:) + Don_conv%qtren1(:,:,k,n)* &
                          pmass(:,:,k)
        end do
        if (id_qtren1_col(n) > 0) then
        used = send_data (id_qtren1_col(n), tempdiag, Time, is, js)
        endif

!---------------------------------------------------------------------
!    define the column-integrated tracer tendency due to mesoscale circ-
!    ulation, in units of kg (tracer) / (m**2 sec). send it to 
!    diag_manager_mod.
!---------------------------------------------------------------------
        tempdiag = 0.0
        do k=1,nlev
          tempdiag(:,:) = tempdiag(:,:) + Don_conv%qtmes1(:,:,k,n)* &
                          pmass(:,:,k)
        end do
        if (id_qtmes1_col(n) > 0) then
        used = send_data (id_qtmes1_col(n), tempdiag, Time, is, js)
        endif

!---------------------------------------------------------------------
!    define the column-integrated tracer redistribution due to meso-
!    scale circulation, in units of kg (tracer) / (m**2 sec). send it 
!    to diag_manager_mod.
!---------------------------------------------------------------------
        tempdiag = 0.0
        do k=1,nlev
          tempdiag(:,:) = tempdiag(:,:) + Don_conv%wtp1(:,:,k,n)*   &
                          pmass(:,:,k)
        end do
        if (id_wtp1_col(n) > 0) then
        used = send_data (id_wtp1_col(n), tempdiag, Time, is, js)
        endif

!---------------------------------------------------------------------
!    define the column-integrated tracer change due to wet deposition in
!    deep convection (cells and mesoscale) in units of kg (tracer) / 
!    (m**2 sec). send it to diag_manager_mod.
!---------------------------------------------------------------------
        tempdiag = 0.0
        do k=1,nlev
          tempdiag(:,:) = tempdiag(:,:) + Don_conv%wetdept(:,:,k,n)*   &
                          pmass(:,:,k)
        end do
        if (id_total_wet_dep_col(n) > 0) then
        used = send_data (id_total_wet_dep_col(n), tempdiag, Time,  &
                                                                 is, js)
        endif

!---------------------------------------------------------------------
!    define the column-integrated tracer change due to wet deposition in
!    mesoscale updrafts, in units of kg (tracer) / (m**2 sec). send it 
!    to diag_manager_mod.
!---------------------------------------------------------------------
       tempdiag = 0.0
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) + Don_conv%wetdepm(:,:,k,n)*   &
                         pmass(:,:,k)
       end do
       if (id_meso_wet_dep_col(n) > 0) then
       used = send_data (id_meso_wet_dep_col(n), tempdiag, Time,  &
                                                                 is, js)
       endif

!---------------------------------------------------------------------
!    define the column-integrated tracer change due to wet deposition 
!    by convective cells, in units of kg (tracer) / (m**2 sec). send it 
!    to diag_manager_mod.
!---------------------------------------------------------------------
       tempdiag = 0.0
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) + Don_conv%wetdepc(:,:,k,n)*   &
                         pmass(:,:,k)
       end do
        if (id_cell_wet_dep_col(n) > 0) then
       used = send_data (id_cell_wet_dep_col(n), tempdiag, Time,  &
                                                                 is, js)
        endif

!-----------------------------------------------------------------
!    define the column-integrated total tracer tendency, in units of 
!    kg (tracer) / (m**2 sec). send it to diag_manager_mod.
!---------------------------------------------------------------------
        tempdiag = 0.0
        do k=1,nlev
          tempdiag(:,:) = tempdiag(:,:) + Don_conv%qtceme(:,:,k,n)* &
                          pmass(:,:,k)
        end do
         if (id_qtceme_col(n) > 0) then
        used = send_data (id_qtceme_col(n), tempdiag, Time, is, js)
        endif
      end do

!---------------------------------------------------------------------
!    send the 2D convection-related diagnostics to diag_manager_mod.
!---------------------------------------------------------------------

!   pressure at lifting condensation level:
       if (id_plcl_deep > 0) then
      used = send_data (id_plcl_deep, Don_cape%plcl, Time, is, js)
       endif

!   pressure at level of free convection:
       if (id_plfc_deep > 0) then
      used = send_data (id_plfc_deep, Don_cape%plfc, Time, is, js)
       endif

!   pressure at level of zero buoyancy:
       if (id_plzb_deep > 0) then
      used = send_data (id_plzb_deep, Don_cape%plzb, Time, is, js)
       endif

!   convective available potential energy:
      if (id_xcape_deep > 0) then
      used = send_data (id_xcape_deep, Don_cape%xcape_lag, Time, is, js)
       endif

!   convective inhibition:
      if (id_coin_deep > 0) then
      used = send_data (id_coin_deep, Don_cape%coin, Time, is, js)
       endif

!   time tendency of cape:
      if (id_dcape_deep > 0) then
      used = send_data (id_dcape_deep, Don_conv%dcape, Time, is, js)
       endif

!   column integrated water vapor:
      if (id_qint_deep > 0) then
      used = send_data (id_qint_deep, Don_cape%qint_lag, Time, is, js)
       endif

!   fractional area of cumulus ensemble members:
      if (id_a1_deep > 0) then
      used = send_data (id_a1_deep, Don_conv%a1, Time, is, js)
       endif

!   fractional area of largest cumulus ensemble member:
      if (id_amax_deep > 0) then
      used = send_data (id_amax_deep, Don_conv%amax, Time, is, js)
       endif

!   upper limit of fractional area based on moisture constraint:
      if (id_amos_deep > 0) then
      used = send_data (id_amos_deep, Don_conv%amos, Time, is, js)
       endif

!   area-weighted total precipitation:
      if (id_tprea1_deep > 0) then
      used = send_data (id_tprea1_deep, total_precip, Time, is, js)
       endif

!   mesoscale cloud fraction:
       if (id_ampta1_deep > 0) then
      used = send_data (id_ampta1_deep, Don_conv%ampta1, Time, is, js)
       endif

!   accumulated low-level parcel displacement:
       if (id_omint_deep > 0) then
         used = send_data (id_omint_deep, parcel_rise, Time, is, js)
       endif

!   area weighted convective precipitation:
       if (id_rcoa1_deep > 0) then
      used = send_data (id_rcoa1_deep, Don_conv%cell_precip,    &
                        Time, is, js)
       endif

   if (Nml%do_ensemble_diagnostics) then

!---------------------------------------------------------------------
!  Donner cumulus ensemble member diagnostics
!---------------------------------------------------------------------

!    GCM model pressure field on full levels:
      used = send_data (id_pfull_cem, Don_cem%pfull, &
                        Time, is, js, 1)

!    GCM model pressure field on half levels:
      used = send_data (id_phalf_cem, Don_cem%phalf, &
                        Time, is, js, 1)

!    GCM model height field on full levels:
      used = send_data (id_zfull_cem, Don_cem%zfull, &
                        Time, is, js, 1)

!    GCM model height field on half levels:
      used = send_data (id_zhalf_cem, Don_cem%zhalf, &
                        Time, is, js, 1)

!    GCM model temperature field on full levels:
      used = send_data (id_temp_cem, Don_cem%temp, &
                        Time, is, js, 1)

!    GCM model mixing ratio field on full levels:
      used = send_data (id_mixing_ratio_cem, Don_cem%mixing_ratio, &
                        Time, is, js, 1)

      do n=1,ncem     ! ensemble member number

!    area-weighted convective precipitation rate:
        used = send_data (id_cpre_cem(n), Don_cem%cell_precip(:,:,n), &
                          Time, is, js)

!    pressure at cloud base:
        used = send_data (id_pb_cem(n), Don_cem%pb(:,:,n), &
                          Time, is, js)

!    pressure at cloud top:
        used = send_data (id_ptma_cem(n), Don_cem%ptma(:,:,n), &
                          Time, is, js)

!    condensation rate profile on lo-res grid:
        used = send_data (id_h1_cem(n), Don_cem%h1(:,:,:,n), &
                          Time, is, js, 1)

!    cloud water profile on lo- or hi-res grid:
        used = send_data (id_qlw_cem(n), Don_cem%qlw(:,:,:,n), &
                          Time, is, js, 1)

!    fraction of condensate that is ice on lo- or hi-res grid:
        used = send_data (id_cfi_cem(n), Don_cem%cfracice(:,:,:,n), &
                          Time, is, js, 1)

!    plume vertical velocity profile on lo- or hi-res grid:
        used = send_data (id_wv_cem(n), Don_cem%wv(:,:,:,n), &
                          Time, is, js, 1)

!    plume cloud radius profile on lo- or hi-res grid:
        used = send_data (id_rcl_cem(n), Don_cem%rcl(:,:,:,n), &
                          Time, is, js, 1)
      enddo

!    fractional area sum:
      used = send_data (id_a1_cem, Don_cem%a1, &
                        Time, is, js)

!    area-weighted mesoscale precipitation rate:
        used = send_data (id_mpre_cem, Don_cem%meso_precip, &
                          Time, is, js)

!    cloud fraction, cells+meso, normalized by a(1,p_b) on lo-res grid:
      used = send_data (id_cual_cem, Don_cem%cual, &
                        Time, is, js, 1)

!    time tendency of temperature due to deep convection on lo-res grid:
      used = send_data (id_tfrc_cem, Don_cem%temperature_forcing, &
                        Time, is, js, 1)
   endif  ! (do_ensemble_diagnostics)

!----------------------------------------------------------------------
!    send diagnostics associated with the monitored output fields.
!----------------------------------------------------------------------
      if (Initialized%monitor_output) then
        do n=1,size(Initialized%Don_monitor,1)
          if (id_extremes(n) > 0) then
            used = send_data (id_extremes(n),   &
                    Initialized%Don_monitor(n)%extrema(is:ie,js:je,:), &
                    Time, is, js,1, mask =    &
                  Initialized%Don_monitor(n)%extrema(is:ie,js:je,:) /= &
                              Initialized%Don_monitor(n)%initial_value )
          endif
          if (id_hits(n) > 0) then
            used = send_data (id_hits(n),  &
                       Initialized%Don_monitor(n)%hits(is:ie,js:je,:), &
                       Time, is, js,1, mask =   &
                 Initialized%Don_monitor(n)%extrema(is:ie,js:je,:) /= &
                              Initialized%Don_monitor(n)%initial_value )
          endif
        end do
      endif

!----------------------------------------------------------------------


end subroutine donner_deep_netcdf


!######################################################################


!#####################################################################

subroutine write_restart (ntracers, Don_save, Initialized, Nml)          

!--------------------------------------------------------------------
!    subroutine write_restart is a template to be used if a native mode
!    restart file MUST be generated. currently, if a native mode file is
!    requested, a netcdf file will be witten instead, and an informative
!    message provided.
!--------------------------------------------------------------------
 
integer, intent(in) :: ntracers
type(donner_initialized_type), intent(inout) :: Initialized
type(donner_save_type), intent(inout) :: Don_save
type(donner_nml_type), intent(inout) :: Nml     

!----------------------------------------------------------------------
!   intent(in) variables:
!
!     ntracers               number of tracers to be transported by
!                            the donner deep convection parameterization
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

!     integer :: unit          ! unit number for restart file
!     integer :: n             ! do-loop index

!-------------------------------------------------------------------
!    currently code is provided only for writing netcdf restart files.
!    if a non-netcdf restart file has been requested, this routine will 
!    issue a message, and then call the routine to write the netcdf file.
!    if the user is insistent on a native mode restart file, the code to
!    read and write such files (subroutines write_restart and 
!    read_restart_file) must be updated to be compatible with  the cur-
!    rent versions of write_restart_nc and read_restart_nc, and the 
!    code immediately below eliminated. the commented code below repres-
!    ents a starting point for the write_restart routine; it is not 
!    kept up-to-date as far as the variables which must be written.
!-------------------------------------------------------------------
      call error_mesg ('donner_deep_mod', 'write_restart: &
          &writing a netcdf restart despite request for native &
           &format (not currently supported); if you must have native &
           &mode, then you must update the source code and remove &
                                               &this if loop.', NOTE)
!      call write_restart_nc (ntracers, Don_save, Initialized, Nml) 

!-------------------------------------------------------------------
!    open unit for restart file.
!-------------------------------------------------------------------
!      unit = open_restart_file ('RESTART/donner_deep.res', 'write')

!-------------------------------------------------------------------
!    file writing is currently single-threaded. write out restart
!    version, time remaining until next call to donner_deep_mod and
!    the frequency of calculating donner_deep convection.
!-------------------------------------------------------------------
!     if (mpp_pe() == mpp_root_pe()) then
!       write (unit) restart_versions(size(restart_versions(:)))
!       write (unit) Initialized%conv_alarm, donner_deep_freq
!     endif

!-------------------------------------------------------------------
!    write out the donner_deep restart variables.
!    cemetf    - heating rate due to donner_deep
!    cememf    - moistening rate due to donner_deep
!    xcape_lag - cape value which will be used on next step in
!                calculation od dcape/dt
!-------------------------------------------------------------------
!     call write_data (unit, Don_save%cemetf)
!     call write_data (unit, Don_save%cememf)
      
!--------------------------------------------------------------------
!    the following variables are needed when a prognostic cloud scheme
!    is being used. they are always present in the restart file, having
!    been initialized to zero, if prognostic clouds are not active.
!--------------------------------------------------------------------
!     call write_data (unit, Don_save%mass_flux)
!     call write_data (unit, Don_save%dql_strat )
!     call write_data (unit, Don_save%dqi_strat )
!     call write_data (unit, Don_save%dqa_strat )

!----------------------------------------------------------------------
!    
!-------------------------------------------------------------------
!    write out more donner_deep restart variables.
!    qint_lag   - column integrated water vapor mixing ratio
!    parcel_disp  - time-integrated low-level vertical displacement
!    tprea1     - precipitation due to donner_deep_mod
!----------------------------------------------------------------------
!     call write_data (unit, Don_save%parcel_disp)
!     call write_data (unit, Don_save%tprea1)
!     call write_data (unit, Don_save%lag_temp)
!     call write_data (unit, Don_save%lag_vapor)
!     call write_data (unit, Don_save%lag_press)
!     call write_data (unit, Don_save%humidity_area)
!     call write_data (unit, Don_save%humidity_ratio)

!---------------------------------------------------------------------
!    write out the number of tracers that are being transported by
!    donner_deep_mod.
!---------------------------------------------------------------------
!     if (mpp_pe() == mpp_root_pe()) then
!       write (unit) ntracers
!     endif

!----------------------------------------------------------------------
!    if tracers are being transported, write out their names and 
!    current time tendencies.
!----------------------------------------------------------------------
!     if (Initialized%do_donner_tracer) then
!       do n=1,ntracers
!         if (mpp_pe() == mpp_root_pe()) then
!           write (unit) Don_save%tracername(n)         
!         endif
!         call write_data(unit, Don_save%tracer_tends(:,:,:,n))
!       end do
!     endif

!-------------------------------------------------------------------
!    close restart file unit.
!------------------------------------------------------------------
!     call close_file (unit)

!---------------------------------------------------------------------


end subroutine write_restart




!######################################################################





!######################################################################



                     end module fms_donner_mod


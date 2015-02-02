module tropchem_driver_mod
!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Larry W. Horowitz
! </CONTACT>

! <OVERVIEW>
!     This code calculates tracer tendencies due to tropospheric chemistry
! </OVERVIEW>

! <DESCRIPTION>
!
! This code calculates chemical production and loss of tracers due
! to tropospheric chemistry. It also includes dry deposition, upper
! boundary conditions, emissions. Off-line sulfate concentrations are
! read in for use in calculating heterogeneous reaction rates (if SO4
! is not included as a tracer).
!
! This module is only activated if do_tropchem=T in tropchem_driver_nml
!
! </DESCRIPTION>


!-----------------------------------------------------------------------

use                    mpp_mod, only : input_nml_file 
use                    fms_mod, only : file_exist,   &
                                       field_exist, &
                                       write_version_number, &
                                       mpp_pe,  &
                                       mpp_root_pe, &
                                       lowercase,   &
                                       open_namelist_file, &
                                       close_file,   &
                                       stdlog, &
                                       check_nml_error, &
                                       error_mesg, &
                                       FATAL, &
                                       WARNING, &
                                       NOTE
use           time_manager_mod, only : time_type, &
                                       get_date, &
                                       set_date, &
                                       set_time, &
                                       days_in_year, &
                                       real_to_time_type, &
                                       time_type_to_real, &
                                       operator(+), operator(-)
use           diag_manager_mod, only : send_data,            &
                                       register_diag_field,  &
                                       register_static_field, &
                                       get_base_time
use         tracer_manager_mod, only : get_tracer_index,     &
                                       get_tracer_names,     &
                                       query_method,         &
                                       check_if_prognostic,  &
                                       NO_TRACER
use          field_manager_mod, only : MODEL_ATMOS,          &
                                       parse
use atmos_tracer_utilities_mod, only : dry_deposition
use              constants_mod, only : grav, rdgas, WTMAIR, WTMH2O, AVOGNO, &
                                       PI, DEG_TO_RAD, SECONDS_PER_DAY
use                    mpp_mod, only : mpp_clock_id,         &
                                       mpp_clock_begin,      &
                                       mpp_clock_end
use           interpolator_mod, only : interpolate_type,     &
                                        interpolate_type_eq, &
                                       interpolator_init,    &
                                    obtain_interpolator_time_slices, &
                                    unset_interpolator_time_flag, &
                                       interpolator_end,     &
                                       interpolator,         &
                                       query_interpolator,   &
                                       init_clim_diag,       &
                                       CONSTANT,             &
                                       INTERP_WEIGHTED_P  
use            time_interp_mod, only : time_interp_init, time_interp
use              mo_chemdr_mod, only : chemdr
use             mo_chemini_mod, only : chemini
use             M_TRACNAME_MOD, only : tracnam         
use                MO_GRID_MOD, only : pcnstm1 
use              CHEM_MODS_MOD, only : phtcnt, gascnt
use               MOZ_HOOK_MOD, only : moz_hook_init
use   strat_chem_utilities_mod, only : strat_chem_utilities_init, &
                                       strat_chem_dcly_dt, &
                                       strat_chem_dcly_dt_time_vary, &
                                       strat_chem_dcly_dt_endts, &
                                       strat_chem_get_aerosol, &
                                       psc_type, &
                                       strat_chem_get_h2so4, &
                                       strat_chem_get_psc, &
                                       strat_chem_destroy_psc, &
                                       strat_chem_psc_sediment, &
                                       strat_chem_get_extra_h2o
use           mo_chem_utls_mod, only : get_spc_ndx
use          atmos_sulfate_mod, only : atmos_sulfate_init, &
                                       atmos_sulfate_time_vary, &
                                       atmos_DMS_emission
use       esfsw_parameters_mod, only: Solar_spect, esfsw_parameters_init 
use astronomy_mod,         only : diurnal_solar, universal_time
use horiz_interp_mod, only: horiz_interp_type, horiz_interp_init, &
                            horiz_interp_new, horiz_interp
use fms_io_mod, only: read_data


implicit none

private

!-----------------------------------------------------------------------
!     ... interfaces
!-----------------------------------------------------------------------
public  tropchem_driver, tropchem_driver_init,  &
        tropchem_driver_time_vary, tropchem_driver_endts

!-----------------------------------------------------------------------
!     ...  declare type that will store the field infomation for the 
!          emission file
!-----------------------------------------------------------------------
type,public :: field_init_type
   character(len=64), pointer :: field_names(:)
end type field_init_type


!-----------------------------------------------------------------------
!     ... namelist
!-----------------------------------------------------------------------
integer, parameter :: maxinv = 100
real               :: relaxed_dt = SECONDS_PER_DAY*10.,     & ! relaxation timescale (sec) for the upper boundary values
                      relaxed_dt_lbc = SECONDS_PER_DAY*10., & ! relaxation timescale (sec) for the lower boundary values
                      ub_pres = 100.e2,               & ! pressure (Pa) above which to apply chemical upper boundary conditions
                      lb_pres = 950.e2                  ! pressure (Pa) below which to apply chemical lower boundary conditions
character(len=64)  :: file_sulfate = 'sulfate.nc',    & ! NetCDF file for sulfate concentrations
                      file_conc = 'conc_all.nc',      & ! NetCDF file for tracer concentrations (initial and fixed)
                      file_emis_1 = 'emissions.',     & ! NetCDF file name (beginning) for emissions
                      file_emis_2 = '.nc',            & ! NetCDF file name (end) for emissions
                      file_emis3d_1 = 'emissions3d.', & ! NetCDF file name (beginning) for 3-D emissions
                      file_emis3d_2 = '.nc',          & ! NetCDF file name (end) for 3-D emissions
                      file_ub = 'ub_vals.nc'            ! NetCDF file for chemical upper boundary conditions
character(len=64)  :: file_dry = 'depvel.nc',         & ! NetCDF file for dry deposition velocities
                      file_aircraft = 'aircraft.nc',  & ! NetCDF file for aircraft emissions
                      file_jval_lut = 'jvals.v5',     & ! ascii file for photolysis rate lookup table
                      file_jval_lut_min = ''            ! ascii file for photolysis rate LUT (for solar min)
character(len=10), dimension(maxinv) :: inv_list =''    ! list of invariant (fixed) tracers
real               :: lght_no_prd_factor = 1.           ! lightning NOx scale factor
real               :: strat_chem_age_factor = 1.        ! scale factor for age of air
real               :: strat_chem_dclydt_factor = 1.     ! scale factor for dcly/dt
logical            :: do_tropchem = .false.             ! Do tropospheric chemistry?
logical            :: use_tdep_jvals = .false.          ! Use explicit temperature dependence for photolysis rates
real               :: o3_column_top = 10.               ! O3 column above model top (DU)
real               :: jno_scale_factor = 1.             ! scale factor for NO photolysis rate (jNO)
logical            :: repartition_water_tracers = .false. ! Allow PSC scheme to act on total water (vapor+condensed)
logical            :: allow_negative_cosz = .false.     ! Allow negative values for cosine of solar zenith angle
logical            :: allow_psc_settling_type1 = .false.! Allow Type-I (NAT) PSCs to settle
logical            :: allow_psc_settling_type2 = .false.! Allow Type-II (ice) PSCs to settle
logical            :: force_cly_conservation = .false.  ! Force chemical conservation of Cly
logical            :: rescale_cly_components = .false.  ! Rescale individual Cly components to total Cly VMR
logical            :: set_min_h2o_strat = .false.       ! Don't allow total water concentration in the stratosphere to fall below 2*CH4_trop
character(len=64)  :: ch4_filename = 'ch4_gblannualdata'! Methane timeseries filename
real               :: ch4_scale_factor = 1.             ! Methane scale factor to convert to VMR (mol/mol)
character(len=64)  :: cfc_lbc_filename = 'chemlbf'      ! Input file for CFC lower boundary conditions
logical            :: time_varying_cfc_lbc = .true.     ! Allow time variation of CFC lower boundary conditions
integer, dimension(6) :: cfc_lbc_dataset_entry = (/ 1, 1, 1, 0, 0, 0 /) ! Entry date for CFC lower boundary condition file
real               :: Tdaily_clim = 297.                ! climatological T for use in MEGAN gamma_age calc
real               :: Pdaily_clim = 420.                ! climatological PPFD for MEGAN light correction 
!++amf/van
integer, parameter :: nveg=5, npft=17, nmos=12          ! number of vegetation types, pfts, and months
!--amf/van
integer            :: verbose = 3                       ! level of diagnostic output
logical            :: retain_cm3_bugs = .false.         ! retain bugs present in code used in CM3
logical            :: do_fastjx_photo = .false.         ! use fastjx routine ?
character(len=32)   :: clouds_in_fastjx = 'lsc_only'    ! nature of clouds seen in fastjx calculation; may currently be 'none' or 'lsc_only' (default)
logical            :: check_convergence = .false.       ! if T, non-converged chem tendencies will not be used
 
namelist /tropchem_driver_nml/    &
                               relaxed_dt, &
                               relaxed_dt_lbc, &
                               ub_pres, &
                               lb_pres, &
                               file_sulfate, &
                               file_conc, &
                               file_emis_1, &
                               file_emis_2, &
                               file_emis3d_1, &
                               file_emis3d_2, &
                               file_ub, &
                               file_dry, &
                               inv_list, & 
                               file_aircraft,&
                               lght_no_prd_factor, &
                               strat_chem_age_factor, &
                               strat_chem_dclydt_factor, &
                               do_tropchem, &
                               use_tdep_jvals, &
                               file_jval_lut, &
                               file_jval_lut_min, &
                               o3_column_top, &
                               jno_scale_factor, &
                               repartition_water_tracers, &
                               allow_negative_cosz, &
                               allow_psc_settling_type1, &
                               allow_psc_settling_type2, &
                               force_cly_conservation, &
                               rescale_cly_components, &
                               set_min_h2o_strat, &
                               ch4_filename, &
                               ch4_scale_factor, &
                               cfc_lbc_filename, &
                               time_varying_cfc_lbc, &
                               cfc_lbc_dataset_entry, &
                               Tdaily_clim, &
                               Pdaily_clim, &
                               verbose,   &
                               retain_cm3_bugs, &
                               do_fastjx_photo, &
                               clouds_in_fastjx, &
                               check_convergence
                              

character(len=7), parameter :: module_name = 'tracers'
real, parameter :: g_to_kg    = 1.e-3,    & !conversion factor (kg/g)
                   m2_to_cm2  = 1.e4,     & !conversion factor (cm2/m2)
                   twopi      = 2.*PI
real, parameter :: emis_cons = WTMAIR * g_to_kg * m2_to_cm2 / AVOGNO
logical, dimension(pcnstm1) :: has_emis = .false., &      ! does tracer have surface emissions?
                               has_emis3d = .false., &    ! does tracer have 3-D emissions?
                               has_xactive_emis = .false., & ! does tracer have interactive emissions?
                               diurnal_emis = .false., &   ! diurnally varying emissions?
                               diurnal_emis3d = .false.    ! diurnally varying 3-D emissions?

type(interpolate_type),dimension(pcnstm1), save :: inter_emis, &
                                                   inter_emis3d, &
                                                   inter_aircraft_emis
type(interpolate_type), save :: airc_default
type(field_init_type),dimension(pcnstm1) :: emis_field_names, &
                                            emis3d_field_names
logical, dimension(pcnstm1) :: has_ubc = .false., &
                               has_lbc = .false., &
                               fixed_lbc_time = .false.
type(time_type), dimension(pcnstm1) :: lbc_entry
logical, dimension(pcnstm1) :: has_airc = .false.
character(len=64),dimension(pcnstm1) :: ub_names, airc_names
real, parameter :: small = 1.e-50
integer :: sphum_ndx=0, cl_ndx=0, clo_ndx=0, hcl_ndx=0, hocl_ndx=0, clono2_ndx=0, &
           cl2o2_ndx=0, cl2_ndx=0, clno2_ndx=0, br_ndx=0, bro_ndx=0, hbr_ndx=0, &
           hobr_ndx=0, brono2_ndx=0, brcl_ndx=0, &
           hno3_ndx=0, o3_ndx=0, &
           no_ndx=0, no2_ndx=0, no3_ndx=0, n_ndx=0, n2o5_ndx=0, ho2no2_ndx=0, &
           pan_ndx=0, onit_ndx=0, mpan_ndx=0, isopno3_ndx=0, onitr_ndx=0, &
           extinct_ndx=0, noy_ndx=0, cly_ndx=0, bry_ndx=0, ch4_ndx=0, &
           dms_ndx=0
logical :: do_interactive_h2o = .false.         ! Include chemical sources/sinks of water vapor?
real, parameter :: solarflux_min = 1.09082, &   ! solar minimum flux (band 18) [W/m2]
                   solarflux_max = 1.14694      ! solar maximum flux (band 18) [W/m2]

!-----------------------------------------------------------------------
!     ... identification numbers for diagnostic fields
!-----------------------------------------------------------------------
integer :: id_sul, id_temp, id_dclydt, id_dbrydt, id_dclydt_chem, &
           id_psc_sat, id_psc_nat, id_psc_ice, id_volc_aer, &
           id_imp_slv_nonconv, id_srf_o3, id_coszen, id_h2o_chem
integer :: inqa, inql, inqi !index of the three water species(nqa, nql, nqi)
integer :: age_ndx ! index of age tracer
logical :: module_is_initialized=.false.
logical :: use_lsc_in_fastjx

integer, dimension(pcnstm1) :: indices, id_prod, id_loss, id_chem_tend, &
                               id_emis, id_emis3d, id_xactive_emis, &
                               id_ub, id_lb, id_airc
integer :: id_so2_emis_cmip, id_nh3_emis_cmip
integer :: id_co_emis_cmip, id_no_emis_cmip
integer :: id_co_emis_cmip2, id_no_emis_cmip2
integer :: id_so2_emis_cmip2, id_nh3_emis_cmip2
integer :: id_glaiage, id_gtemp, id_glight, id_tsfc, id_fsds, id_ctas, id_cfsds
integer :: isop_oldmonth = 0
logical :: newmonth            
logical :: has_ts_avg = .true.   ! currently reading in from monthly mean files.
integer, dimension(phtcnt)  :: id_jval
integer, dimension(gascnt)  :: id_rate_const
integer :: id_prodox, id_lossox ! for production and loss of ox.(jmao,1/1/2011)

type(interpolate_type), save :: conc       ! used to read in the concentration of OH and CH4
type(interpolate_type), save :: sulfate    ! used to read in the data for sulate
type(interpolate_type), save :: ub_default ! used for the upper bound data
type(interpolate_type),dimension(pcnstm1), save :: ub

type :: lb_type
   real, dimension(:), pointer :: gas_value
   type(time_type), dimension(:), pointer :: gas_time
end type lb_type
type(lb_type), dimension(pcnstm1) :: lb

type(interpolate_type), save :: drydep_data_default
integer :: clock_id,ndiag

real, allocatable, dimension(:,:,:) :: ecisop, pctpft   !emission capacities, % pft
real, allocatable, dimension(:,:,:,:) :: mlai ! monthly lai for each pft (could eventually tie to LM3)
real, allocatable, dimension(:,:) :: emisop_month !isop emissions with monthly lai & age gammas applied
real, allocatable, dimension(:,:) :: diag_gamma_lai_age   ! for combined lai and age gammas
real, allocatable, dimension(:,:) :: diag_gamma_light, diag_gamma_temp  ! for gamma light / T
real, allocatable, dimension(:,:) :: diag_climtas, diag_climfsds ! climatological tas and fsds

!if set up to read from restart file, then these would be 2D arrays.. for now, use mm clim. values
!real, allocatable, dimension(:,:) :: ts_avg   ! surface air T, averaged over some period (K) 
!real, allocatable, dimension(:,:) :: fsds_avg  ! avg shortwave down in visible (W/m2) 

!monthly mean sfc air T and sw down at surface, from C. Wiedinmyer 2/18/09 - Sheffield inputs (Princeton)
! CW provided 1948-2000; current input files take 1980-2000 average
real, allocatable, dimension(:,:,:) :: ts_avg  ! climatological monthly mean surface air T
real, allocatable, dimension(:,:,:) :: fsds_avg  ! climat. montly mean total shortwave down (W/m2)

type (horiz_interp_type), save :: Interp


!---- version number ---------------------------------------------------
character(len=128), parameter :: version     = '$Id: tropchem_driver.F90,v 20.0 2013/12/13 23:25:19 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: tikal $'
!-----------------------------------------------------------------------

contains


!#######################################################################

! <SUBROUTINE NAME="tropchem_driver">
!   <OVERVIEW>
!     Tropospheric chemistry driver.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine calculates the sources and sinks of tracers
!     due to tropospheric chemistry. It is called from atmos_tracer_driver.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call tropchem_driver (lon, lat, land, pwt, r, chem_dt,           &
!                           Time, phalf, pfull, t, is, ie, js, je, dt, &
!                           z_half, z_full, q, tsurf, albedo, coszen,  &
!                           area, w10m, flux_sw_down_vis_dir, flux_sw_down_vis_dif, &
!                           half_day, &
!                           Time_next, rdiag,  kbot)
!   </TEMPLATE>
!   <IN NAME="lon" TYPE="real" DIM="(:,:)">
!     The longitudes for the local domain.
!   </IN>
!   <IN NAME="lat" TYPE="real" DIM="(:,:)">
!     The latitudes for the local domain.
!   </IN>
!   <IN NAME="land" TYPE="real" DIM="(:,:)">
!     The latitudes for the local domain.
!   </IN>
!   <IN NAME="pwt" TYPE="real" DIM="(:,:,:)">
!     Pressure weighting (air mass) for each layer (kg/m2)
!   </IN>
!   <IN NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer mixing ratios (tropchem tracers in VMR)
!   </IN>
!   <IN NAME="Time, Time_next" TYPE="time_type">
!     Model time
!   </IN>
!   <IN NAME="phalf" TYPE="real" DIM="(:,:,:)">
!     Pressure on the model half levels (Pa)
!   </IN>
!   <IN NAME="pfull" TYPE="real" DIM="(:,:,:)">
!     Pressure on the model full levels (Pa)
!   </IN>
!   <IN NAME="t" TYPE="real" DIM="(:,:,:)">
!     Temperature.
!   </IN>
!   <IN NAME="is, js" TYPE="integer">
!     Local domain start indices
!   </IN>
!   <IN NAME="ie, je" TYPE="integer">
!     Local domain end indices
!   </IN>
!   <IN NAME="dt" TYPE="real">
!     Model physics timestep (s)
!   </IN>
!   <IN NAME="z_half" TYPE="real" DIM="(:,:,:)">
!     Height at model half levels (m)
!   </IN>
!   <IN NAME="z_full" TYPE="real" DIM="(:,:,:)">
!     Height at model full levels (m)
!   </IN>
!   <IN NAME="q" TYPE="real" DIM="(:,:,:)">
!     Specific humidity (kg/kg)
!   </IN>
!   <IN NAME="tsurf" TYPE="real" DIM="(:,:)">
!     Surface temperature (K)
!   </IN>
!   <IN NAME="albedo" TYPE="real" DIM="(:,:)">
!     Surface albedo
!   </IN>
!   <IN NAME="coszen" TYPE="real" DIM="(:,:)">
!     Cosine of the solar zenith angle
!   </IN>
!   <IN NAME="area" TYPE="real" DIM="(:,:)">
!     Grid box area (m^2)
!   </IN>
!   <IN NAME="w10m" TYPE="real" DIM="(:,:)">
!     Windspeed at 10m (m/s)
!   </IN>
!   <IN NAME="flux_sw_down_vis_dir" TYPE="real" DIM="(:,:)">
!     Surface downward visible radiation (W/m2)
!   </IN>
!   <IN NAME="flux_sw_down_vis_dif" TYPE="real" DIM="(:,:)">
!     Surface downward visible radiation (W/m2)
!   </IN>
!   <IN NAME="half_day" TYPE="real" DIM="(:,:)">
!     Half-day length  (dimensionless; 0 to pi)
!   </IN>
!   <OUT NAME="chem_dt" TYPE="real" DIM="(:,:,:,:)">
!     Tracer tendencies from tropospheric chemistry (VMR/s)
!   </OUT>
!   <INOUT NAME="rdiag" TYPE="real" DIM="(:,:,:,:)">
!     Diagnostic tracer mixing ratios (tropchem tracers in VMR),
!     updated on output
!   </INOUT>
!   <IN NAME="kbot" TYPE="integer, optional" DIM="(:,:)">
!     Integer array describing which model layer intercepts the surface.
!   </IN>

subroutine tropchem_driver( lon, lat, land, ocn_flx_fraction, pwt, r, chem_dt,                 &
                            Time, phalf, pfull, t, is, ie, js, je, dt,       &
                            z_half, z_full, q, tsurf, albedo, coszen, rrsun, &
                            area, w10m, flux_sw_down_vis_dir, flux_sw_down_vis_dif, &
                            half_day, &
                            Time_next, rdiag,  kbot )

!-----------------------------------------------------------------------
   real, intent(in),    dimension(:,:)            :: lon, lat
   real, intent(in),    dimension(:,:)            :: land    ! land fraction
   real, intent(in),    dimension(:,:)            :: ocn_flx_fraction ! grid box fraction over which DMS flux from ocean occurs
   real, intent(in),    dimension(:,:,:)          :: pwt
   real, intent(in),    dimension(:,:,:,:)        :: r
   real, intent(out),   dimension(:,:,:,:)        :: chem_dt
   type(time_type), intent(in)                    :: Time, Time_next     
   integer, intent(in)                            :: is, ie, js, je
   real, intent(in),    dimension(:,:,:)          :: phalf,pfull,t
   real, intent(in)                               :: dt      ! timestep (s)
   real, intent(in),    dimension(:,:,:)          :: z_half  ! height in meters at half levels
   real, intent(in),    dimension(:,:,:)          :: z_full  ! height in meters at full levels
   real, intent(in),    dimension(:,:,:)          :: q       ! specific humidity at current time step (kg/kg)
   real, intent(in),    dimension(:,:)            :: tsurf   ! surface temperature (K)
   real, intent(in),    dimension(:,:)            :: albedo  ! surface albedo
   real, intent(in),    dimension(:,:)            :: coszen  ! cosine of solar zenith angle
   real, intent(in)                               :: rrsun   ! earth-sun distance factor (r_avg/r)^2
   real, intent(in),    dimension(:,:)            :: area    ! grid box area (m^2)
   real, intent(in),    dimension(:,:)            :: w10m    ! wind speed at 10m (m/s)
   real, intent(in), dimension(:,:)               :: flux_sw_down_vis_dir !W/m2 direct visible sfc flux
   real, intent(in), dimension(:,:)               :: flux_sw_down_vis_dif !W/m2 diffuse visible sfc flux
   real, intent(in), dimension(:,:)               :: half_day! half-day length (0 to pi)   
   real, intent(inout), dimension(:,:,:,:)        :: rdiag   ! diagnostic tracer concentrations
   integer, intent(in),  dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
   real, dimension(size(r,1),size(r,2),size(r,3)) :: sulfate_data
!   real, dimension(size(r,1),size(r,2),size(r,3)) :: ub_temp,rno
   real, dimension(size(r,1),size(r,2),size(r,3),maxinv) :: inv_data
   real, dimension(size(r,1),size(r,2)) :: emis
   real, dimension(size(r,1),size(r,2), pcnstm1) :: emisz
   real, dimension(size(r,1),size(r,2),size(r,3)) :: emis3d, xactive_emis
   real, dimension(size(r,1),size(r,2),size(r,3)) :: age, cly0, cly, cly_ratio, &
                                                     bry, dclydt, dbrydt, noy, &
                                                     extinct, strat_aerosol
   real, dimension(size(r,1),size(r,2),size(r,3),3) :: psc_vmr_save, dpsc_vmr
   real, dimension(size(r,1),size(r,2)) :: tsfcair, pwtsfc, flux_sw_down_vis
   integer :: i,j,k,n,kb,id,jd,kd,ninv,ntp, index1, index2
!  integer :: nno,nno2
   integer :: inv_index
   integer :: plonl
   logical :: used
   real :: scale_factor, frac
   real,  dimension(size(r,1),size(r,3)) :: pdel, h2so4, h2o_temp, qlocal, cloud_water
   real, dimension(size(r,1),size(r,2),size(r,3),pcnstm1)  :: r_temp, r_in, emis_source, r_ub, airc_emis
   real, dimension(size(r,1),size(r,2),size(r,3)) :: tend_tmp, extra_h2o
   real, dimension(pcnstm1) :: r_lb
   real, dimension(size(land,1), size(land,2)) :: oro ! 0 and 1 rep. of land
   real, dimension(size(r,1),size(r,2)) :: coszen_local, fracday_local
   real :: rrsun_local
   real, dimension(size(r,1),size(r,2),size(r,3),pcnstm1) :: prod, loss
!  add new arrays for ox budget (jmao,1/1/2011)
   real, dimension(size(r,1),size(r,2),size(r,3)):: prodox, lossox
   real, dimension(size(r,1),size(r,2),size(r,3),phtcnt) :: jvals
   real, dimension(size(r,1),size(r,2),size(r,3),gascnt) :: rate_constants
   real, dimension(size(r,1),size(r,2),size(r,3)) :: imp_slv_nonconv
   real :: solar_phase
   type(psc_type) :: psc
   type(time_type) :: lbc_Time
!-----------------------------------------------------------------------

!<ERROR MSG="tropchem_driver_init must be called first." STATUS="FATAL">
!   Tropchem_driver_init needs to be called before tracer_driver.
!</ERROR>
   if (.not. module_is_initialized)  &
      call error_mesg ('Tropchem_driver','tropchem_driver_init must be called first.', FATAL)

   ntp = size(r,4)
   plonl = size(r,1)

   where(land(:,:) >= 0.5)
      oro(:,:) = 1.
   elsewhere
      oro(:,:) = 0.
   endwhere

   id=size(r,1); jd=size(r,2); kd=size(r,3)
   
   ninv=0
   do n = 1, size(inv_list)
      if(inv_list(n) /= '') then
         ninv = ninv + 1
      else
         exit
      end if
   end do
 
   emis_source(:,:,:,:) = 0.0
   airc_emis(:,:,:,:) = 0.0

   tsfcair(:,:) = t(:,:,kd)
   pwtsfc(:,:) = t(:,:,kd)
   do n = 1, pcnstm1
!-----------------------------------------------------------------------
!     ... read in the surface emissions, using interpolator
!-----------------------------------------------------------------------
      if (has_emis(n)) then
         call read_2D_emis_data( inter_emis(n), emis, Time, Time_next, &
                                 emis_field_names(n)%field_names, &
                                 diurnal_emis(n), coszen, half_day, lon, &
                                 is, js, id_emis(n) )
         if (tracnam(n) == 'NO') then
           emisz(:,:,n) = emis(:,:)
           if (id_no_emis_cmip > 0) then
             used = send_data(id_no_emis_cmip,emis*1.0e04*0.030/AVOGNO, &
                                  Time_next, &
                                                  is_in=is,js_in=js)
           endif
         endif
         if (tracnam(n) == 'CO') then
           emisz(:,:,n) = emis(:,:)
           if (id_co_emis_cmip > 0) then
             used = send_data(id_co_emis_cmip,emis*1.0e04*0.028/AVOGNO,Time_next, &
                                                  is_in=is,js_in=js)
           endif
         endif
         if (tracnam(n) == 'SO2') then
           emisz(:,:,n) = emis(:,:)
           if (id_so2_emis_cmip > 0) then
             used = send_data(id_so2_emis_cmip,emis*1.0e04*0.064/AVOGNO,Time_next, &
                                                  is_in=is,js_in=js)
           endif
         endif
         if (tracnam(n) == 'NH3') then
           emisz(:,:,n) = emis(:,:)
           if (id_nh3_emis_cmip > 0) then
             used = send_data(id_nh3_emis_cmip,emis*1.0e04*0.017/AVOGNO,Time_next, &
                                                  is_in=is,js_in=js)
           endif
         endif

         if (present(kbot)) then
            do j=1,jd
               do i=1,id
                  kb=kbot(i,j)
                  emis_source(i,j,kb,n) = emis(i,j)/pwt(i,j,kb) * emis_cons
               end do
            end do
         else
            emis_source(:,:,kd,n) = emis(:,:)/pwt(:,:,kd) * emis_cons
         end if
      end if

!-----------------------------------------------------------------------
!     ... read in the 3-D emissions, using interpolator
!-----------------------------------------------------------------------
      if (has_emis3d(n)) then
         call read_3D_emis_data( inter_emis3d(n), emis3d, Time, Time_next,phalf, &
                                 emis3d_field_names(n)%field_names, &
                                 diurnal_emis3d(n), coszen, half_day, lon, &
                                 is, js, id_emis3d(n) )
      
         emis_source(:,:,:,n) = emis_source(:,:,:,n) &
                              + emis3d(:,:,:)/pwt(:,:,:) * emis_cons
         if (tracnam(n) == 'SO2') then
           do k=1, size(emis3d,3)
           emisz(:,:,n) = emisz(:,:,n) + emis3d(:,:,k)
           end do
         endif
         if (tracnam(n) == 'NO') then
           do k=1, size(emis3d,3)
           emisz(:,:,n) = emisz(:,:,n) + emis3d(:,:,k)
           end do
         endif
         if (tracnam(n) == 'CO') then
           do k=1, size(emis3d,3)
           emisz(:,:,n) = emisz(:,:,n) + emis3d(:,:,k)
           end do
         endif
         if (tracnam(n) == 'NH3') then
           do k=1, size(emis3d,3)
           emisz(:,:,n) = emisz(:,:,n) + emis3d(:,:,k)
           end do
         endif
      end if

!-----------------------------------------------------------------------
!     ... calculate interactive emissions
!-----------------------------------------------------------------------
!     if (has_xactive_emis(n)) then
      if ( has_xactive_emis(n) .or. id_xactive_emis(n)>0 ) then
         select case (trim(tracnam(n)))
         case ('ISOP')

            flux_sw_down_vis = flux_sw_down_vis_dir+flux_sw_down_vis_dif

            call calc_xactive_isop ( n, Time, Time_next, lon, lat, oro, pwtsfc, is, js, &
                 area, land, tsfcair, flux_sw_down_vis, &
                 coszen, emis, id_gamma_lai_age=id_glaiage, &
                 id_gamma_temp=id_gtemp, id_gamma_light=id_glight, &
                 id_tsfcair=id_tsfc, id_fsdvd=id_fsds, &
                 id_climtas=id_ctas, id_climfsds=id_cfsds, id_emis_diag=id_xactive_emis(n) )
            if (has_xactive_emis(n)) then
               if (present(kbot)) then
                  do j=1,jd
                     do i=1,id
                        kb=kbot(i,j)

                        emis_source(i,j,kb,n) = emis_source(i,j,kb,n) &
                                              + emis(i,j)/pwt(i,j,kb)*emis_cons

                     end do
                  end do
               else
                  emis_source(:,:,kd,n) = emis_source(:,:,kd,n) &
                                        + emis(:,:)/pwt(:,:,kd) * emis_cons
               end if
            end if
         case ('DMS')
            call calc_xactive_emis( n, Time, Time_next,lon, lat, pwt, is, ie, js, je, &
                 area, land, ocn_flx_fraction,tsurf, w10m, xactive_emis, &
                 kbot=kbot, id_emis_diag=id_xactive_emis(n) )
            if (has_xactive_emis(n)) then
               emis_source(:,:,:,n) = emis_source(:,:,:,n) + xactive_emis(:,:,:)
            end if
         case default
            if (has_xactive_emis(n)) then
            call error_mesg ('tropchem_driver','Interactive emissions not defined for species: '//trim(tracnam(n)), FATAL)
            end if
         end select
      end if

!-----------------------------------------------------------------------
!     ... read in the aircraft emissions
!-----------------------------------------------------------------------
      if(has_airc(n)) then
         call interpolator( inter_aircraft_emis(n), Time, phalf, &
                            airc_emis(:,:,:,n), trim(airc_names(n)),is,js)
         if(id_airc(n) > 0)&
              used = send_data(id_airc(n),airc_emis(:,:,:,n),Time_next, is_in=is, js_in=js)
    
         if (tracnam(n) == 'CO') then
           do k=1, size(emis3d,3)
           emisz(:,:,n) = emisz(:,:,n) + airc_emis(:,:,k,n)
           end do
         endif
         if (tracnam(n) == 'NO') then
           do k=1, size(emis3d,3)
           emisz(:,:,n) = emisz(:,:,n) + airc_emis(:,:,k,n)
           end do
         endif
         if (tracnam(n) == 'SO2') then
           do k=1, size(emis3d,3)
           emisz(:,:,n) = emisz(:,:,n) + airc_emis(:,:,k,n)
           end do
         endif
         if (tracnam(n) == 'NH3') then
           do k=1, size(emis3d,3)
           emisz(:,:,n) = emisz(:,:,n) + airc_emis(:,:,k,n)
           end do
         endif
         airc_emis(:,:,:,n) = airc_emis(:,:,:,n)/pwt(:,:,:)*emis_cons
!     end if
      end if
         if (tracnam(n) == 'NO') then
           if (id_no_emis_cmip2 > 0) then
             used = send_data(id_no_emis_cmip2,emisz(:,:,n)*1.0e04*0.030/AVOGNO,Time_next, &
                                                 is_in=is,js_in=js)
           endif
         endif
         if (tracnam(n) == 'CO') then
           if (id_co_emis_cmip2 > 0) then
             used = send_data(id_co_emis_cmip2,emisz(:,:,n)*1.0e04*0.028/AVOGNO,Time_next, &
                                                 is_in=is,js_in=js)
           endif
         endif
         if (tracnam(n) == 'SO2') then
           if (id_so2_emis_cmip2 > 0) then
             used = send_data(id_so2_emis_cmip2,emisz(:,:,n)*1.0e04*0.064/AVOGNO,Time_next, &
                                                 is_in=is,js_in=js)
           endif
         endif
         if (tracnam(n) == 'NH3') then
           if (id_nh3_emis_cmip2 > 0) then
             used = send_data(id_nh3_emis_cmip2,emisz(:,:,n)*1.0e04*0.017/AVOGNO,Time_next, &
                                                  is_in=is,js_in=js)
           endif
         endif
   end do

!-----------------------------------------------------------------------
!     ... read in the concentrations of "invariant" (i.e., prescribed)
!         species
!-----------------------------------------------------------------------
   do n = 1,ninv
      call interpolator( conc, Time, phalf, inv_data(:,:,:,n), &
                         trim(inv_list(n)), is, js)
      inv_index = get_tracer_index( MODEL_ATMOS, trim(inv_list(n)) ) - ntp
      rdiag(:,:,:,inv_index) = inv_data(:,:,:,n)
   end do
      
!-----------------------------------------------------------------------
!     ... read in the sulfate aerosol concentrations
!-----------------------------------------------------------------------
   call interpolator(sulfate, Time, phalf, sulfate_data, 'sulfate', is,js)
   used = send_data(id_sul, sulfate_data, Time_next, is_in=is, js_in=js)

!  call mpp_clock_begin(clock_id)

   chem_dt(:,:,:,:) =0.
    
!-----------------------------------------------------------------------
!     ... assign concentrations of prognostic (r) and diagnostic (rdiag)
!         species to r_temp
!-----------------------------------------------------------------------
   do n = 1,pcnstm1
      if(indices(n) <= ntp) then
         r_temp(:,:,:,n) = r(:,:,:,indices(n))
      else
         r_temp(:,:,:,n) = rdiag(:,:,:,indices(n)-ntp)
      end if
   end do

!-----------------------------------------------------------------------
!     ... convert to H2O VMR
!-----------------------------------------------------------------------
   if (sphum_ndx > 0) then
      r_temp(:,:,:,sphum_ndx) = r_temp(:,:,:,sphum_ndx) * WTMAIR / WTMH2O
   end if

!-----------------------------------------------------------------------
!     ... convert volcanic aerosol extinction into aerosol surface area
!-----------------------------------------------------------------------
   if (extinct_ndx > 0 .and. extinct_ndx <= ntp) then
      extinct(:,:,:) = r(:,:,:,extinct_ndx)
   else if (extinct_ndx > ntp) then
      extinct(:,:,:) = rdiag(:,:,:,extinct_ndx-ntp)
   else
      extinct(:,:,:) = 0.
   end if
   call strat_chem_get_aerosol( extinct, strat_aerosol )

!-----------------------------------------------------------------------
!     ... get age of air
!-----------------------------------------------------------------------
   if(age_ndx > 0 .and. age_ndx <= ntp) then
      age(:,:,:) = r(:,:,:,age_ndx)
   else
      age(:,:,:) = 0.
   end if

!-----------------------------------------------------------------------
!     ... Chemical families
!-----------------------------------------------------------------------
   cly0(:,:,:) = 0.
   if (cl_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,cl_ndx)
   end if
   if (clo_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,clo_ndx)
   end if
   if (hcl_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,hcl_ndx)
   end if
   if (hocl_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,hocl_ndx)
   end if
   if (clono2_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,clono2_ndx)
   end if
   if (cl2o2_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,cl2o2_ndx)*2
   end if
   if (cl2_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,cl2_ndx)*2
   end if
   if (clno2_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,clno2_ndx)
   end if
   if (brcl_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,brcl_ndx)
   end if

!-----------------------------------------------------------------------
!     ... cosine of solar zenith angle
!-----------------------------------------------------------------------
   if (allow_negative_cosz) then
      call diurnal_solar( lat, lon, Time, coszen_local, fracday_local, &
                          rrsun_local, dt_time=real_to_time_type(dt), &
                          allow_negative_cosz=.true. )
   else
      coszen_local(:,:) = coszen(:,:)
      rrsun_local = rrsun
   end if

   r_temp(:,:,:,:) = MAX(r_temp(:,:,:,:),small)
  
   do j = 1,jd
      do k = 1,kd
         pdel(:,k) = phalf(:,j,k+1) - phalf(:,j,k)
      end do
      qlocal(:,:) = q(:,j,:)
      
!-----------------------------------------------------------------------
!     ... get stratospheric h2so4
!-----------------------------------------------------------------------
      call strat_chem_get_h2so4( pfull(:,j,:), age(:,j,:), h2so4 )

!-----------------------------------------------------------------------
!     ... compute PSC amounts
!-----------------------------------------------------------------------
      if (sphum_ndx>0) then
         h2o_temp(:,:) = r_temp(:,j,:,sphum_ndx)
      else
         h2o_temp(:,:) = qlocal(:,:) * WTMAIR/WTMH2O
      end if
      cloud_water(:,:) = MAX(r(:,j,:,inql)+r(:,j,:,inqi),0.)
      if (repartition_water_tracers) then
         h2o_temp(:,:) = h2o_temp(:,:) + cloud_water(:,:) * WTMAIR/WTMH2O
      end if
      if (set_min_h2o_strat) then
         call strat_chem_get_extra_h2o( h2o_temp, age(:,j,:), r_temp(:,j,:,ch4_ndx), Time, extra_h2o(:,j,:) )
         h2o_temp(:,:) = h2o_temp(:,:) + extra_h2o(:,j,:)
      end if

      call strat_chem_get_psc( t(:,j,:), pfull(:,j,:), &
                               r_temp(:,j,:,hno3_ndx), h2o_temp(:,:), &
                               h2so4, strat_aerosol(:,j,:), psc, psc_vmr_out=psc_vmr_save(:,j,:,:) )

      if (repartition_water_tracers) then
         cloud_water(:,:) = MAX(0.,cloud_water(:,:) - psc_vmr_save(:,j,:,3)*WTMH2O/WTMAIR) ! reduce cloud_water by amount of type-II PSC
         h2o_temp(:,:) = h2o_temp(:,:) - cloud_water(:,:) * WTMAIR/WTMH2O                  ! remaining water is present as vapor
      end if
      if (sphum_ndx>0) then
         r_temp(:,j,:,sphum_ndx) = h2o_temp(:,:)
      end if
      qlocal(:,:) = h2o_temp(:,:) * WTMH2O/WTMAIR
      r_in(:,j,:,:) = r_temp(:,j,:,:)

!-----------------------------------------------------------------------
!     ... get solar cycle phase (use radiation band #18)
!-----------------------------------------------------------------------
      solar_phase = Solar_spect%solflxbandref(Solar_spect%nbands)
      solar_phase = (solar_phase-solarflux_min)/(solarflux_max-solarflux_min)

!-----------------------------------------------------------------------
!     ... call chemistry driver
!-----------------------------------------------------------------------
      call chemdr(r_temp(:,j,:,:),             & ! species volume mixing ratios (VMR)
                  r(:,j,:,:),                  &
                  phalf(:,j,:),                & ! pressure at boundaries (Pa)
                  pwt(:,j,:) ,                 & ! column air density (Kg/m2)  
                  do_fastjx_photo,             & ! true = use fastjx photo
                  use_lsc_in_fastjx,           & ! use lsc clouds in fastjx
                  j,                           & ! j
!                 0,                           & ! time step index
                  Time_next,                   & ! time
                  lat(:,j),                    & ! latitude
                  lon(:,j),                    & ! longitude
                  dt,                          & ! timestep in seconds
                  phalf(:,j,SIZE(phalf,3)),    & ! surface press ( pascals )  
                  phalf(:,j,1),                & ! model top pressure (pascals)
                  pfull(:,j,:),                & ! midpoint press ( pascals )
                  pdel,                        & ! delta press across midpoints
!                 oro(:,j),                    & ! surface orography flag
!                 tsurf(:,j),                  & ! surface temperature
                  z_full(:,j,:),               & ! height at midpoints ( m )
                  z_half(:,j,:),               & ! height at interfaces ( m )
                  MAX(r(:,j,:,inqa),0.),       & ! cloud fraction
                  cloud_water(:,:),            & ! total cloud water (kg/kg)
                  t(:,j,:),                    & ! temperature
                  inv_data(:,j,:,:),           & ! invariant species
                  qlocal(:,:),                 & ! specific humidity ( kg/kg )
                  albedo(:,j),                 & ! surface albedo
                  coszen_local(:,j),           & ! cosine of solar zenith angle
                  rrsun_local,                 & ! earth-sun distance factor
                  prod(:,j,:,:),               & ! chemical production rate
                  loss(:,j,:,:),               & ! chemical loss rate
                  jvals(:,j,:,:),              & ! photolysis rates (s^-1)
                  rate_constants(:,j,:,:),     & ! kinetic rxn rate constants (cm^3 molec^-1 s^-1 for 2nd order)
                  sulfate_data(:,j,:),         & ! sulfate aerosol
                  psc,                         & ! polar stratospheric clouds (PSCs)
                  do_interactive_h2o,          & ! include h2o sources/sinks?
                  solar_phase,                 & ! solar cycle phase (1=max, 0=min)
                  imp_slv_nonconv(:,j,:),      & ! flag for non-convergence of implicit solver
                  plonl,                       & ! number of longitudes
                  prodox(:,j,:),               & ! production of ox(jmao,1/1/2011)
                  lossox(:,j,:),               & ! loss of ox(jmao,1/1/2011)
                  retain_cm3_bugs, check_convergence)

      call strat_chem_destroy_psc( psc )

   end do

   r_temp(:,:,:,:) = MAX( r_temp(:,:,:,:), small )
   if (allow_psc_settling_type1 .or. allow_psc_settling_type2) then
      call strat_chem_psc_sediment( psc_vmr_save, pfull, dt, dpsc_vmr )
      if (.not. allow_psc_settling_type1) dpsc_vmr(:,:,:,2) = 0.
      if (.not. allow_psc_settling_type2) dpsc_vmr(:,:,:,3) = 0.
   end if

!-----------------------------------------------------------------------
!     ... output diagnostics
!-----------------------------------------------------------------------
   do n = 1,pcnstm1
      if(id_prod(n)>0) then
         used = send_data(id_prod(n),prod(:,:,:,n),Time_next,is_in=is,js_in=js)
      end if
      if(id_loss(n)>0) then
         used = send_data(id_loss(n),loss(:,:,:,n),Time_next,is_in=is,js_in=js)
      end if
      
      if (n == sphum_ndx) then
         scale_factor = WTMAIR/WTMH2O
! add PSC ice back to H2O
!         if (repartition_water_tracers) then
!           r_temp(:,:,:,n) = r_temp(:,:,:,n) + &
!                             MAX( 0.,psc_vmr_save(:,:,:,3) - MAX(r(:,:,:,inql)+r(:,:,:,inqi),0.)*scale_factor ) + &
!                             dpsc_vmr(:,:,:,3)
!        else
!           r_temp(:,:,:,n) = r_temp(:,:,:,n) + psc_vmr_save(:,:,:,3)+dpsc_vmr(:,:,:,3)
!        end if
!        if (set_min_h2o_strat) then
!           r_temp(:,:,:,n) = MAX( r_temp(:,:,:,n) - extra_h2o(:,:,:), small )
!        end if
      else
         scale_factor = 1.
      end if

!     if (n == hno3_ndx) then
!        r_temp(:,:,:,n) = r_temp(:,:,:,n) + psc_vmr_save(:,:,:,2)+dpsc_vmr(:,:,:,2) ! add PSC NAT back to gas-phase HNO3
!     end if

!-----------------------------------------------------------------------
!     ... compute tendency
!-----------------------------------------------------------------------
      tend_tmp(:,:,:) = ( r_temp(:,:,:,n) - r_in(:,:,:,n) )/dt
      if(indices(n) <= ntp) then
!-----------------------------------------------------------------------
!     ... prognostic species
!-----------------------------------------------------------------------
!        tend_tmp(:,:,:) = ( r_temp(:,:,:,n) - MAX(r(:,:,:,indices(n))*scale_factor,small) )/dt
         chem_dt(:,:,:,indices(n)) = airc_emis(:,:,:,n) + emis_source(:,:,:,n) + tend_tmp(:,:,:)
      else
!-----------------------------------------------------------------------
!     ... diagnostic species
!-----------------------------------------------------------------------
!        tend_tmp(:,:,:) = ( r_temp(:,:,:,n) - MAX(rdiag(:,:,:,indices(n)-ntp)*scale_factor,small) )/dt
         rdiag(:,:,:,indices(n)-ntp) = r_temp(:,:,:,n)
      end if
!-----------------------------------------------------------------------
!     ... output diagnostic tendency
!-----------------------------------------------------------------------
      if(id_chem_tend(n)>0) then
         used = send_data( id_chem_tend(n), tend_tmp(:,:,:), Time_next, is_in=is,js_in=js)
      end if
     
!-----------------------------------------------------------------------
!     ... apply upper boundary condition
!-----------------------------------------------------------------------
      if(has_ubc(n)) then
         call interpolator(ub(n), Time, phalf, r_ub(:,:,:,n), trim(ub_names(n)), is, js)
         if(id_ub(n)>0) then
            used = send_data(id_ub(n), r_ub(:,:,:,n), Time_next, is_in=is, js_in=js)
         end if
         where (pfull(:,:,:) < ub_pres)            
            chem_dt(:,:,:,indices(n)) = (r_ub(:,:,:,n) - r(:,:,:,indices(n))) / relaxed_dt
         endwhere
      end if

!-----------------------------------------------------------------------
!     ... apply lower boundary condition
!-----------------------------------------------------------------------
      if(has_lbc(n)) then
         if (fixed_lbc_time(n)) then
            lbc_Time = lbc_entry(n)
         else
            lbc_Time = Time
         end if
         call time_interp( lbc_Time, lb(n)%gas_time(:), frac, index1, index2 )
         r_lb(n) = lb(n)%gas_value(index1) + frac*( lb(n)%gas_value(index2) - lb(n)%gas_value(index1) )
         if(id_lb(n)>0) then
            used = send_data(id_lb(n), r_lb(n), Time_next)
         end if
         where (pfull(:,:,:) > lb_pres)
            chem_dt(:,:,:,indices(n)) = (r_lb(n) - r(:,:,:,indices(n))) / relaxed_dt_lbc
         endwhere
      end if

   end do
!-----------------------------------------------------------------------
!     ...send ox budget(jmao,1/1/2011)
!-----------------------------------------------------------------------
   if(id_prodox>0) then
      used = send_data(id_prodox, prodox(:,:,:), Time_next, is_in=is, js_in=js)
   end if
   if(id_lossox>0) then
      used = send_data(id_lossox, lossox(:,:,:), Time_next, is_in=is, js_in=js)
   end if
!-----------------------------------------------------------------------
!     ... surface concentration diagnostics
!-----------------------------------------------------------------------
      if ( o3_ndx>0 ) then
         used = send_data(id_srf_o3, r_temp(:,:,size(r_temp,3),o3_ndx), Time_next, is_in=is, js_in=js)
      end if

   
!-----------------------------------------------------------------------
!     ... special case(nox = no + no2)
!-----------------------------------------------------------------------
!  nno = get_tracer_index(MODEL_ATMOS,'no')
!  nno2 = get_tracer_index(MODEL_ATMOS,'no2')
!  if((nno /= 0) .and. (nno2 /= 0)) then
!     rno(:,:,:) = r(:,:,:,nno)/ MAX((r(:,:,:,nno) + r(:,:,:,nno2)),1.e-30)
     
!     call interpolator(ub, Time,phalf,ub_temp,'nox',is,js)

!     where(pfull(:,:,:) < ub_pres)
!        chem_dt(:,:,:,nno) =((rno(:,:,:)*ub_temp(:,:,:))-r(:,:,:,nno)) / relaxed_dt
!        chem_dt(:,:,:,nno2) = (((1.-rno(:,:,:))*ub_temp(:,:,:))-r(:,:,:,nno2)) / &
!             relaxed_dt
!     endwhere
!  end if

!-----------------------------------------------------------------------
!     ... Chemical families (Cly)
!-----------------------------------------------------------------------
   cly(:,:,:) = 0.
   if (cl_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,cl_ndx)
   end if
   if (clo_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,clo_ndx)
   end if
   if (hcl_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,hcl_ndx)
   end if
   if (hocl_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,hocl_ndx)
   end if
   if (clono2_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,clono2_ndx)
   end if
   if (cl2o2_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,cl2o2_ndx)*2
   end if
   if (cl2_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,cl2_ndx)*2
   end if
   if (clno2_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,clno2_ndx)
   end if
   if (brcl_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,brcl_ndx)
   end if

!-----------------------------------------------------------------------
!     ... Cly chemical tendency diagnostic
!-----------------------------------------------------------------------
   if (id_dclydt_chem>0) then
      used = send_data(id_dclydt_chem, (cly(:,:,:)-cly0(:,:,:))/dt, Time_next, is_in=is, js_in=js)
   end if

!-----------------------------------------------------------------------
!     ... Cly conservation
!-----------------------------------------------------------------------
   if (force_cly_conservation .or. rescale_cly_components) then
      if (rescale_cly_components) then
         cly_ratio(:,:,:) = r(:,:,:,cly_ndx) / MAX( cly(:,:,:), small )
         cly(:,:,:) = r(:,:,:,cly_ndx)
      else if (force_cly_conservation) then
         cly_ratio(:,:,:) = cly0(:,:,:) / MAX( cly(:,:,:), small )
         cly(:,:,:) = cly0(:,:,:)
      end if
      if (cl_ndx>0) then
         r_temp(:,:,:,cl_ndx) = r_temp(:,:,:,cl_ndx) * cly_ratio(:,:,:)
      end if
      if (clo_ndx>0) then
         r_temp(:,:,:,clo_ndx) = r_temp(:,:,:,clo_ndx) * cly_ratio(:,:,:)
      end if
      if (hcl_ndx>0) then
         r_temp(:,:,:,hcl_ndx) = r_temp(:,:,:,hcl_ndx) * cly_ratio(:,:,:)
      end if
      if (hocl_ndx>0) then
         r_temp(:,:,:,hocl_ndx) = r_temp(:,:,:,hocl_ndx) * cly_ratio(:,:,:)
      end if
      if (clono2_ndx>0) then
         r_temp(:,:,:,clono2_ndx) = r_temp(:,:,:,clono2_ndx) * cly_ratio(:,:,:)
      end if
      if (cl2o2_ndx>0) then
         r_temp(:,:,:,cl2o2_ndx) = r_temp(:,:,:,cl2o2_ndx) * cly_ratio(:,:,:)
      end if
      if (cl2_ndx>0) then
         r_temp(:,:,:,cl2_ndx) = r_temp(:,:,:,cl2_ndx) * cly_ratio(:,:,:)
      end if
      if (clno2_ndx>0) then
         r_temp(:,:,:,clno2_ndx) = r_temp(:,:,:,clno2_ndx) * cly_ratio(:,:,:)
      end if
      if (brcl_ndx>0) then
         r_temp(:,:,:,brcl_ndx) = r_temp(:,:,:,brcl_ndx) * cly_ratio(:,:,:)
      end if
   end if

!-----------------------------------------------------------------------
!     ... Chemical families (Bry, NOy)
!-----------------------------------------------------------------------
   bry(:,:,:) = 0.
   noy(:,:,:) = 0.
   if (clono2_ndx>0) then
      noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,clono2_ndx)
   end if
   if (br_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,br_ndx)
   end if
   if (bro_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,bro_ndx)
   end if
   if (hbr_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,hbr_ndx)
   end if
   if (hobr_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,hobr_ndx)
   end if
   if (brono2_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,brono2_ndx)
      noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,brono2_ndx)
   end if
   if (brcl_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,brcl_ndx)
   end if
   if (n_ndx>0) then
      noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,n_ndx)
   end if
   if (no_ndx>0) then
      noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,no_ndx)
   end if
   if (no2_ndx>0) then
      noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,no2_ndx)
   end if
   if (no3_ndx>0) then
      noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,no3_ndx)
   end if
   if (hno3_ndx>0) then
      noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,hno3_ndx)
   end if
   if (n2o5_ndx>0) then
      noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,n2o5_ndx)*2
   end if
   if (ho2no2_ndx>0) then
      noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,ho2no2_ndx)
   end if
   if (pan_ndx>0) then
      noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,pan_ndx)
   end if
   if (mpan_ndx>0) then
      noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,mpan_ndx)
   end if
   if (onit_ndx>0) then
      noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,onit_ndx)
   end if
   if (isopno3_ndx>0) then
      noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,isopno3_ndx)
   end if
   if (onitr_ndx>0) then
      noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,onitr_ndx)
   end if

!-----------------------------------------------------------------------
!     ... stratospheric Cly and Bry source
!-----------------------------------------------------------------------
   if(age_ndx > 0 .and. age_ndx <= ntp) then
      call strat_chem_dcly_dt(Time, phalf, is, js, age, cly, bry, dclydt, dbrydt)
      do k = 1,kd
         where( coszen(:,:) > 0. )
            dclydt(:,:,k) = 2*dclydt(:,:,k)
            dbrydt(:,:,k) = 2*dbrydt(:,:,k)
         elsewhere
            dclydt(:,:,k) = 0.
            dbrydt(:,:,k) = 0.
         end where
      end do
   else
      dclydt(:,:,:) = 0.
      dbrydt(:,:,:) = 0.
   end if
   if (cl_ndx>0) then
      chem_dt(:,:,:,indices(cl_ndx)) = chem_dt(:,:,:,indices(cl_ndx)) + dclydt(:,:,:)
      used = send_data(id_dclydt, dclydt, Time_next, is_in=is, js_in=js)
   end if
   if (br_ndx>0) then
      chem_dt(:,:,:,indices(br_ndx)) = chem_dt(:,:,:,indices(br_ndx)) + dbrydt(:,:,:)
      used = send_data(id_dbrydt, dbrydt, Time_next, is_in=is, js_in=js)
   end if
   
!-----------------------------------------------------------------------
!     ... Set diagnostic tracers for chemical families
!-----------------------------------------------------------------------
   if (noy_ndx > ntp) then
      rdiag(:,:,:,noy_ndx-ntp) = noy(:,:,:)
   end if
   if (cly_ndx > ntp) then
      rdiag(:,:,:,cly_ndx-ntp) = cly(:,:,:) + dclydt(:,:,:)*dt
   else if (cly_ndx > 0) then
      chem_dt(:,:,:,cly_ndx) = dclydt(:,:,:)
   end if
   if (bry_ndx > ntp) then
      rdiag(:,:,:,bry_ndx-ntp) = bry(:,:,:) + dbrydt(:,:,:)*dt
   end if

!-----------------------------------------------------------------------
!     ... Photolysis rates
!-----------------------------------------------------------------------
   do n = 1,phtcnt
      if(id_jval(n)>0) then
         used = send_data(id_jval(n),jvals(:,:,:,n),Time_next,is_in=is,js_in=js)
      end if
   end do

!-----------------------------------------------------------------------
!     ... Kinetic reaction rates
!-----------------------------------------------------------------------
   do n = 1,gascnt
      if(id_rate_const(n)>0) then
         used = send_data(id_rate_const(n),rate_constants(:,:,:,n),Time_next,is_in=is,js_in=js)
      end if
   end do

!-----------------------------------------------------------------------
!     ... Output diagnostics
!-----------------------------------------------------------------------
   used = send_data(id_volc_aer, strat_aerosol, Time_next, is_in=is, js_in=js)
   used = send_data(id_psc_sat, psc_vmr_save(:,:,:,1), Time_next, is_in=is, js_in=js)
   used = send_data(id_psc_nat, psc_vmr_save(:,:,:,2), Time_next, is_in=is, js_in=js)
   used = send_data(id_psc_ice, psc_vmr_save(:,:,:,3), Time_next, is_in=is, js_in=js)
   if (id_h2o_chem>0) then
      if (sphum_ndx>0) then
         used = send_data(id_h2o_chem, r_temp(:,:,:,sphum_ndx), Time_next, is_in=is, js_in=js)
      else
         used = send_data(id_h2o_chem, q(:,:,:)*WTMAIR/WTMH2O, Time_next, is_in=is, js_in=js)
      end if
   end if      
   used = send_data(id_coszen, coszen_local(:,:), Time_next, is_in=is, js_in=js)
   used = send_data(id_imp_slv_nonconv,imp_slv_nonconv(:,:,:),Time_next,is_in=is,js_in=js)

!-----------------------------------------------------------------------
!     ... convert H2O VMR tendency to specific humidity tendency
!-----------------------------------------------------------------------
   if (sphum_ndx > 0) then
      n = indices(sphum_ndx)
      chem_dt(:,:,:,n) = chem_dt(:,:,:,n) * WTMH2O / WTMAIR
!     chem_dt(:,:,:,n) = 0.
   end if

!  call mpp_clock_end(clock_id)
   
!-----------------------------------------------------------------------
    
end subroutine tropchem_driver
!</SUBROUTINE>

!#######################################################################

! <FUNCTION NAME="tropchem_driver_init">
!   <OVERVIEW>
!     Initializes the tropospheric chemistry driver.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine initializes the tropospheric chemistry module.
!     It is called from atmos_tracer_driver_init.
!     Data sets are read in for dry deposition, upper boundary conditions,
!     and emissions. Off-line sulfate concentrations are also read in for
!     use in calculating heterogeneous reaction rates (if SO4 is not
!     included as a tracer).
!   </DESCRIPTION>
!   <TEMPLATE>
!     Ltropchem = tropchem_driver_init( r, mask, axes, Time, &
!                                       lonb_mod, latb_mod, phalf, &
!                                       drydep_data )
!   </TEMPLATE>
!   <IN NAME="mask" TYPE="real, optional" DIM="(:,:,:)">
!      optional mask that designates which grid points
!      are above (1) or below (0) the ground
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="lonb_mod" TYPE="real" DIM="(:,:)">
!     The longitude corners for the local domain.
!   </IN>
!   <IN NAME="latb_mod" TYPE="real" DIM="(:,:)">
!     The latitude corners for the local domain.
!   </IN>
!   <IN NAME="phalf" TYPE="real" DIM="(:,:,:)">
!     Pressure on the model half levels (Pa)
!   </IN>
!   <OUT NAME="drydep_data" TYPE="interpolate_type" DIM="(:)">
!     Tracer dry deposition velocities
!   </OUT>
!   <OUT NAME="Ltropchem" TYPE="logical">
!     Do tropospheric chemistry? (Output as function value)
!   </OUT>
!   <INOUT NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer mixing ratios (tropchem tracers in VMR)
!   </INOUT>

function tropchem_driver_init( r, mask, axes, Time, &
                               lonb_mod, latb_mod, phalf, &
                               drydep_data ) result(Ltropchem)

!-----------------------------------------------------------------------
!
!   r    = tracer fields dimensioned as (nlon,nlat,nlev,ntrace)
!   mask = optional mask (0. or 1.) that designates which grid points
!          are above (=1.) or below (=0.) the ground dimensioned as
!          (nlon,nlat,nlev).
!
!-----------------------------------------------------------------------
   real, intent(inout), dimension(:,:,:,:) :: r
   real, intent(in),    dimension(:,:,:), optional :: mask
   type(time_type), intent(in) :: Time
   integer        , intent(in) :: axes(4)
   real, intent(in), dimension(:,:) :: lonb_mod
   real, intent(in), dimension(:,:) :: latb_mod
   real, intent(in),dimension(:,:,:) :: phalf
   type(interpolate_type), intent(out) :: drydep_data(:)

   logical :: Ltropchem
   integer :: flag_file, flag_spec, flag_fixed
   integer :: n,i
   integer :: ierr, io, logunit
   character(len=64) :: nc_file,filename,specname
   character(len=256) :: control=''
   character(len=64) :: name=''
   type(interpolate_type) :: init_conc
   character(len=64),dimension(pcnstm1) :: emis_files = '', &
                                           emis3d_files = '', &
                                           conc_files = '', &
                                           ub_files = '', &
                                           lb_files = '', &
                                           dry_files, &
                                           wet_ind, &
                                           conc_names, &
                                           dry_names, &
                                           airc_files
   logical :: tracer_initialized

   integer :: unit
   character(len=16) ::  fld

   integer :: flb, series_length, year, diy
   real :: input_time
   real :: scale_factor, extra_seconds, fixed_year
   type(time_type) :: Year_t
!
!-----------------------------------------------------------------------
!
   if (module_is_initialized) return

!-----------------------------------------------------------------------
!     ... write version number
!-----------------------------------------------------------------------
   call write_version_number(version, tagname)
    
!-----------------------------------------------------------------------
!     ... read namelist
!-----------------------------------------------------------------------
   if(file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=tropchem_driver_nml, iostat=io)
      ierr = check_nml_error(io,'tropchem_driver_nml')
#else
      unit = open_namelist_file('input.nml')
      ierr=1; do while (ierr /= 0)
      read(unit, nml = tropchem_driver_nml, iostat=io, end=10)
      ierr = check_nml_error (io, 'tropchem_driver_nml')
      end do
10    call close_file(unit)
#endif
   end if
  
   logunit = stdlog()
   if(mpp_pe() == mpp_root_pe()) then       
      write(logunit, nml=tropchem_driver_nml)
      verbose = verbose + 1
   end if

   Ltropchem = do_tropchem
   if (.not. Ltropchem) then
      return
   end if
      
!-------------------------------------------------------------------------
!     ... Make sure input value for clouds_in_fastjx is a valid option.
!-------------------------------------------------------------------------
   if (trim(clouds_in_fastjx) == 'none') then
     use_lsc_in_fastjx = .false.
   else if (trim(clouds_in_fastjx) == 'lsc_only') then
     use_lsc_in_fastjx = .true.
   else
     call error_mesg ('tropchem_driver_init', &
                     ' invalid string for clouds_in_fastjx', FATAL)
   endif

!-----------------------------------------------------------------------
!     ... Setup sulfate input/interpolation
!-----------------------------------------------------------------------
   call interpolator_init( sulfate, trim(file_sulfate), lonb_mod, latb_mod, &
                           data_out_of_bounds=(/CONSTANT/),      &
                           vert_interp=(/INTERP_WEIGHTED_P/) )

!-----------------------------------------------------------------------
!     ... Initialize chemistry driver
!-----------------------------------------------------------------------
   call chemini( file_jval_lut, file_jval_lut_min, use_tdep_jvals, &
                 o3_column_top, jno_scale_factor, verbose,   &
                 retain_cm3_bugs, do_fastjx_photo)
   
!-----------------------------------------------------------------------
!     ... set initial value of indices
!-----------------------------------------------------------------------
   indices(:) = 0
   do i=1,pcnstm1
      n = get_tracer_index(MODEL_ATMOS, tracnam(i))
      if (trim(tracnam(i)) == 'H2O') then
         if (n <= 0) then
            n = get_tracer_index(MODEL_ATMOS, 'sphum')
         end if
         sphum_ndx = i
         do_interactive_h2o = .true.
      end if
      if (n >0) then
         indices(i) = n
         if (indices(i) > 0 .and. mpp_pe() == mpp_root_pe()) then 
            write(*,30) tracnam(i),indices(i)
            write(logunit,30) trim(tracnam(i)),indices(i)
         end if
      else
!<ERROR MSG="Tropospheric chemistry tracer not found in field table" STATUS="WARNING">
!   A tropospheric chemistry tracer was not included in the field table
!</ERROR>
         call error_mesg ('tropchem_driver_init', trim(tracnam(i)) // ' is not found', WARNING)
      end if
   end do
30 format (A,' was initialized as tracer number ',i3)

   cl_ndx     = get_spc_ndx('Cl')
   clo_ndx    = get_spc_ndx('ClO')
   hcl_ndx    = get_spc_ndx('HCl')
   hocl_ndx   = get_spc_ndx('HOCl')
   clono2_ndx = get_spc_ndx('ClONO2')
   cl2o2_ndx  = get_spc_ndx('Cl2O2')
   cl2_ndx    = get_spc_ndx('Cl2')
   clno2_ndx  = get_spc_ndx('ClNO2')
   br_ndx     = get_spc_ndx('Br')
   bro_ndx    = get_spc_ndx('BrO')
   hbr_ndx    = get_spc_ndx('HBr')
   hobr_ndx   = get_spc_ndx('HOBr')
   brono2_ndx = get_spc_ndx('BrONO2')
   brcl_ndx   = get_spc_ndx('BrCl')
   hno3_ndx   = get_spc_ndx('HNO3')
   no_ndx     = get_spc_ndx('NO')
   no2_ndx    = get_spc_ndx('NO2')
   no3_ndx    = get_spc_ndx('NO3')
   n_ndx      = get_spc_ndx('N')
   n2o5_ndx   = get_spc_ndx('N2O5')
   ho2no2_ndx = get_spc_ndx('HO2NO2')
   pan_ndx    = get_spc_ndx('PAN')
   onit_ndx   = get_spc_ndx('ONIT')
   mpan_ndx   = get_spc_ndx('MPAN')
   isopno3_ndx= get_spc_ndx('ISOPNO3')
   onitr_ndx  = get_spc_ndx('ONITR')
   o3_ndx     = get_spc_ndx('O3')
   ch4_ndx    = get_spc_ndx('CH4')
   dms_ndx    = get_spc_ndx('DMS')

   extinct_ndx = get_tracer_index(MODEL_ATMOS, 'Extinction')
   noy_ndx     = get_tracer_index(MODEL_ATMOS, 'NOy')
   cly_ndx     = get_tracer_index(MODEL_ATMOS, 'Cly')
   bry_ndx     = get_tracer_index(MODEL_ATMOS, 'Bry')

!-----------------------------------------------------------------------
!     ... Check Cly settings
!-----------------------------------------------------------------------
   if (rescale_cly_components) then
      if (cly_ndx == NO_TRACER .or. .not. check_if_prognostic(MODEL_ATMOS,cly_ndx)) then
         call error_mesg ('tropchem_driver_init', &
                          'rescale_cly_components=T requires Cly to be registered as a prognostic tracer', FATAL)
      end if
      if (force_cly_conservation) then
         call error_mesg ('tropchem_driver_init', &
                          'rescale_cly_components=T incompatible with force_cly_conservation=T setting', FATAL)
      end if
   end if

!-----------------------------------------------------------------------
!     ... Setup dry deposition
!-----------------------------------------------------------------------
   call tropchem_drydep_init( dry_files, dry_names, &
                              lonb_mod, latb_mod, &
                              drydep_data )

!-----------------------------------------------------------------------
!     ... Setup upper boundary condition data
!-----------------------------------------------------------------------
   call interpolator_init( ub_default, trim(file_ub), lonb_mod, latb_mod, &
                           data_out_of_bounds=(/CONSTANT/),          &
                           vert_interp=(/INTERP_WEIGHTED_P/) )

!-----------------------------------------------------------------------
!     ... Set up concentration input/interpolation
!-----------------------------------------------------------------------
   call interpolator_init( conc, trim(file_conc), lonb_mod, latb_mod, &
                           data_out_of_bounds=(/CONSTANT/),&
                           vert_interp=(/INTERP_WEIGHTED_P/) )

!-----------------------------------------------------------------------
!     ... Set up aircraft emissions interpolation
!-----------------------------------------------------------------------
   call interpolator_init( airc_default, trim(file_aircraft), lonb_mod, latb_mod,&
                           data_out_of_bounds=(/CONSTANT/), &
                           vert_interp=(/INTERP_WEIGHTED_P/))

!-----------------------------------------------------------------------
!     ... Setup emissions input/interpolation
!-----------------------------------------------------------------------
   do i = 1,pcnstm1
      nc_file = trim(file_emis_1)//lowercase(trim(tracnam(i)))//trim(file_emis_2)
      call init_emis_data( inter_emis(i), MODEL_ATMOS, 'emissions', indices(i), nc_file, &
                           lonb_mod, latb_mod, emis_field_names(i), &
                           has_emis(i), diurnal_emis(i), axes, Time )
      if( has_emis(i) ) emis_files(i) = trim(nc_file)
        
!-----------------------------------------------------------------------
!     ... Vertically-distributed emissions
!-----------------------------------------------------------------------
      nc_file = trim(file_emis3d_1)//lowercase(trim(tracnam(i)))//trim(file_emis3d_2)
      call init_emis_data( inter_emis3d(i), MODEL_ATMOS, 'emissions3d', indices(i), nc_file, &
                           lonb_mod, latb_mod, emis3d_field_names(i), &
                           has_emis3d(i), diurnal_emis3d(i), axes, Time )
      if( has_emis3d(i) ) emis3d_files(i) = trim(nc_file)

!-----------------------------------------------------------------------
!     ... Interactive emissions
!-----------------------------------------------------------------------
      call init_xactive_emis( MODEL_ATMOS, 'xactive_emissions', indices(i), tracnam(i), &
                              axes, Time, lonb_mod, latb_mod, phalf, &
                              has_xactive_emis(i), id_xactive_emis(i), mask )

!-----------------------------------------------------------------------
!     ... Upper boundary condition
!-----------------------------------------------------------------------
      if( query_method('upper_bound', MODEL_ATMOS,indices(i),name,control) ) then
         if( trim(name)=='file' ) then
            flag_file = parse(control, 'file',filename)
            flag_spec = parse(control, 'name',specname)

            if( flag_file > 0 .and. trim(filename) /= trim(file_ub) ) then
               ub_files(i) = trim(filename)
               call interpolator_init(ub(i), trim(filename), lonb_mod, latb_mod, &
                       data_out_of_bounds=(/CONSTANT/),          &
                       vert_interp=(/INTERP_WEIGHTED_P/))
            else
               ub_files(i) = trim(file_ub)
               ub(i) = ub_default
            end if
            if(flag_spec > 0) then
               ub_names(i) = trim(specname)
            else
               ub_names(i) = trim(lowercase(tracnam(i)))
            end if

            has_ubc(i) = .true.
              
         end if
      end if

!-----------------------------------------------------------------------
!     ... Lower boundary condition
!-----------------------------------------------------------------------
      lbc_entry(i) = get_base_time()     
      if( query_method('lower_bound', MODEL_ATMOS,indices(i),name,control) ) then
         if( trim(name)=='file' ) then
            flag_file = parse(control, 'file', filename)
            flag_spec = parse(control, 'factor', scale_factor)
            flag_fixed = parse(control, 'fixed_year', fixed_year)
            if( flag_file > 0 ) then
               lb_files(i) = 'INPUT/' // trim(filename)
               if( file_exist(lb_files(i)) ) then
                  flb = open_namelist_file( lb_files(i) )
                  read(flb, FMT='(i12)') series_length
                  allocate( lb(i)%gas_value(series_length), &
                            lb(i)%gas_time(series_length) )
!---------------------------------------------------------------------
!    convert the time stamps of the series to time_type variables.     
!---------------------------------------------------------------------
                  do n = 1,series_length
                     read (flb, FMT = '(2f12.4)') input_time, lb(i)%gas_value(n)
                     year = INT(input_time)
                     Year_t = set_date(year,1,1,0,0,0)
                     diy = days_in_year (Year_t)
                     extra_seconds = (input_time - year)*diy*SECONDS_PER_DAY 
                     lb(i)%gas_time(n) = Year_t + set_time(NINT(extra_seconds), 0)
                  end do
                  if (flag_spec > 0) then
                     lb(i)%gas_value(:) = lb(i)%gas_value(:) * scale_factor
                  end if
                  call close_file( flb )
                  if( flag_fixed > 0 ) then
                     fixed_lbc_time(i) = .true.
                     year = INT(fixed_year)
                     Year_t = set_date(year,1,1,0,0,0)
                     diy = days_in_year (Year_t)
                     extra_seconds = (fixed_year - year)*diy*SECONDS_PER_DAY 
                     lbc_entry(i) = Year_t + set_time(NINT(extra_seconds), 0)
                  end if
               else
                  call error_mesg ('tropchem_driver_init', &
                                   'Failed to find input file '//trim(lb_files(i)), FATAL)
               end if
            else
               call error_mesg ('tropchem_driver_init', 'Tracer '//trim(lowercase(tracnam(i)))// &
                                ' has lower_bound specified without a filename', FATAL)
            end if
            has_lbc(i) = .true.
         end if
      end if

!-----------------------------------------------------------------------
!     ... Initial conditions
!-----------------------------------------------------------------------
      tracer_initialized = .false.
      if ( field_exist('INPUT/atmos_tracers.res.nc', lowercase(tracnam(i))) .or. &
           field_exist('INPUT/fv_tracer.res.nc', lowercase(tracnam(i))) .or. &
           field_exist('INPUT/tracer_'//trim(lowercase(tracnam(i)))//'.res', lowercase(tracnam(i))) ) then
         tracer_initialized = .true.
      end if

      if(.not. tracer_initialized) then
         if( query_method('init_conc',MODEL_ATMOS,indices(i),name,control) ) then
            if( trim(name) == 'file' ) then
               flag_file = parse(control, 'file',filename)
               flag_spec = parse(control, 'name',specname)

               if( flag_file>0 .and. trim(filename) /= trim(file_conc) ) then
                  conc_files(i) = trim(filename)
                  call interpolator_init( init_conc,trim(filename),lonb_mod,latb_mod,&
                                          data_out_of_bounds=(/CONSTANT/), &
                                          vert_interp=(/INTERP_WEIGHTED_P/) )
               else
                  conc_files(i) = trim(file_conc)
                  init_conc = conc
               end if
                  
               if( flag_spec > 0 ) then
                  conc_names(i) = trim(lowercase(specname))
                  specname = lowercase(specname)
               else
                  conc_names(i) = trim(lowercase(tracnam(i)))
                  specname = trim(lowercase(tracnam(i)))
               end if
                  
               call interpolator(init_conc, Time, phalf,r(:,:,:,indices(i)),trim(specname))
            end if
         end if
      end if
             
!-----------------------------------------------------------------------
!     ... Aircraft emissions
!-----------------------------------------------------------------------
      if( query_method('aircraft_emis',MODEL_ATMOS,indices(i),name,control) ) then
         has_airc(i) = .true.
         if( trim(name) == 'file' ) then
            flag_file = parse(control,'file',filename)
            flag_spec = parse(control,'name',specname)

            if( flag_file >0 .and. trim(filename) /= trim(lowercase(file_aircraft)) ) then
               airc_files(i) = trim(filename)
               call interpolator_init( inter_aircraft_emis(i),trim(filename), lonb_mod, latb_mod, &
                                       data_out_of_bounds=(/CONSTANT/), &
                                       vert_interp=(/INTERP_WEIGHTED_P/) )
            else
               airc_files(i) = trim(file_aircraft)
               inter_aircraft_emis(i) = airc_default
            end if
               
            if( flag_spec >0 ) then
               airc_names(i) = trim(specname)
            else
               airc_names(i) = trim(lowercase(tracnam(i)))
            end if
         end if
      end if

             
!-----------------------------------------------------------------------
!     ... Wet deposition
!-----------------------------------------------------------------------
      if( query_method('wet_deposition',MODEL_ATMOS,indices(i),name,control) ) then
         wet_ind(i) = 'This species has wet deposition'
      else
         wet_ind(i) = ''
      end if
         
   end do
   
!-----------------------------------------------------------------------
!     ... Print out settings for tracer
!-----------------------------------------------------------------------
   if( mpp_pe() == mpp_root_pe() ) then
      write(logunit,*) '---------------------------------------------------------------------------------------'
      do i = 1,pcnstm1
         write(logunit,*) 'The tracname index is ',i
         write(logunit,*) 'The tracname is ',tracnam(i)
         if(check_if_prognostic(MODEL_ATMOS,indices(i))) then
            write(logunit,*) 'This is a prognostic tracer.'
         else
            write(logunit,*) 'This is a diagnostic tracer.'
         end if
         if(has_emis(i)) then
            write(logunit,*)'Emissions from file: ',trim(emis_files(i))
         end if
         if(has_emis3d(i)) then
            write(logunit,*)'3-D Emissions from file: ',trim(emis3d_files(i))
         end if
         if(has_ubc(i)) then
            write(logunit,*)'Upper BC from file: ',trim(ub_files(i)), &
                             ', with the name of ',trim(ub_names(i))
         end if
         if(has_lbc(i)) then
            write(logunit,*)'Lower BC from file: ',trim(lb_files(i))
            if (fixed_lbc_time(i)) then
               write(logunit,*) '... with fixed year'
            end if
         end if
         if(conc_files(i) /= '') then
            write(logunit,*)'Concentration from file: ',trim(conc_files(i)), &
                             ', with the name of ',trim(conc_names(i))
         end if
         if(dry_files(i) /= '') then
            write(logunit,*)'Dry deposition velocity from file: ',trim(dry_files(i)), &
                             ' with the name of '//trim(dry_names(i))
         end if
         if(wet_ind(i) /= '') then
            write(logunit,*) wet_ind(i)
         end if
         if(has_airc(i)) then
            write(logunit,*)'Aircraft emissions from file: ',trim(airc_files(i)), &
                             ' with the name of '//trim(airc_names(i))
         end if
         write(logunit,*) '---------------------------------------------------------------------------------------'
      end do
   end if


!-----------------------------------------------------------------------
!     ... Get the index number for the cloud variables
!-----------------------------------------------------------------------
   inqa = get_tracer_index(MODEL_ATMOS,'cld_amt') ! cloud fraction
   inql = get_tracer_index(MODEL_ATMOS,'liq_wat') ! cloud liquid specific humidity
   inqi = get_tracer_index(MODEL_ATMOS,'ice_wat') ! cloud ice water specific humidity
      
   age_ndx = get_tracer_index(MODEL_ATMOS,'age')  ! age tracer

!-----------------------------------------------------------------------
!     ... Call the chemistry hook init routine
!-----------------------------------------------------------------------
   call moz_hook_init( lght_no_prd_factor, Time, axes, verbose )

!-----------------------------------------------------------------------
!     ... Initializations for stratospheric chemistry
!-----------------------------------------------------------------------
!++lwh
   if (set_min_h2o_strat) then
      if (ch4_ndx>0) then
         if (.not. has_lbc(ch4_ndx)) then
            call error_mesg ('Tropchem_driver','set_min_h2o_strat=T, but LBC not set for CH4', FATAL)
         end if
      else
         call error_mesg ('Tropchem_driver','set_min_h2o_strat=T, but CH4 not included in chemistry solver', FATAL)
      end if
   end if
   call strat_chem_utilities_init( lonb_mod, latb_mod, &
                                   strat_chem_age_factor, strat_chem_dclydt_factor, &
                                   set_min_h2o_strat, ch4_filename, ch4_scale_factor, &
                                   fixed_lbc_time(ch4_ndx), lbc_entry(ch4_ndx), &
                                   cfc_lbc_filename, time_varying_cfc_lbc, cfc_lbc_dataset_entry )
!--lwh
   id_dclydt      = register_diag_field( module_name, 'cly_chem_dt', axes(1:3), Time, 'cly_chem_dt', 'VMR/s' )
   id_dclydt_chem = register_diag_field( module_name, 'cly_chem_dt_diag', axes(1:3), Time, 'cly_chem_dt_diag', 'VMR/s' )
   id_dbrydt      = register_diag_field( module_name, 'bry_chem_dt', axes(1:3), Time, 'bry_chem_dt', 'VMR/s' )

   id_volc_aer = register_diag_field( module_name, 'volc_aer_SA', axes(1:3), Time, 'volcanic_aerosol_surface_area', 'cm2/cm3' )
   id_psc_sat  = register_diag_field( module_name, 'psc_sat', axes(1:3), Time, 'psc_sat', 'VMR' )
   id_psc_nat  = register_diag_field( module_name, 'psc_nat', axes(1:3), Time, 'psc_nat', 'VMR' )
   id_psc_ice  = register_diag_field( module_name, 'psc_ice', axes(1:3), Time, 'psc_ice', 'VMR' )
   id_h2o_chem = register_diag_field( module_name, 'h2o_chem', axes(1:3), Time, 'h2o_chem', 'VMR' )

!-----------------------------------------------------------------------
!     ... Initialize additional diagnostics
!-----------------------------------------------------------------------
   id_sul = register_diag_field( module_name, 'sulfate', axes(1:3), Time, 'sulfate', 'VMR' )
   id_coszen = register_diag_field( module_name, 'coszen_tropchem', axes(1:2), Time, &
                                             'cosine_sza_tropchem', 'none' )
   id_imp_slv_nonconv = register_diag_field( module_name, 'imp_slv_nonconv', axes(1:3), Time, &
                                             'tropchem_implicit_solver_not_converged', 'VMR' )
   id_srf_o3 = register_diag_field( module_name, 'o3_srf', axes(1:2), Time, 'o3_srf', 'VMR' )

!-----------------------------------------------------------------------
!     ... Register diagnostic fields for species tendencies
!-----------------------------------------------------------------------
   id_co_emis_cmip =     &
        register_diag_field( module_name, 'co_emis_cmip', axes(1:2), &
                             Time, 'co_emis_cmip', 'kg/m2/s')
   id_no_emis_cmip =     &
        register_diag_field( module_name, 'no_emis_cmip', axes(1:2), &
                            Time, 'no_emis_cmip', 'kg/m2/s')  
   id_so2_emis_cmip =     &
        register_diag_field( module_name, 'so2_emis_cmip', axes(1:2), &
                             Time, 'so2_emis_cmip', 'kg/m2/s')
   id_nh3_emis_cmip =     &
        register_diag_field( module_name, 'nh3_emis_cmip', axes(1:2), &
                            Time, 'nh3_emis_cmip', 'kg/m2/s')  

   id_co_emis_cmip2 =     &
        register_diag_field( module_name, 'co_emis_cmip2', axes(1:2), &
                             Time, 'co_emis_cmip2', 'kg/m2/s')
   id_no_emis_cmip2 =     &
        register_diag_field( module_name, 'no_emis_cmip2', axes(1:2), &
                            Time, 'no_emis_cmip2', 'kg/m2/s')  
   id_so2_emis_cmip2 =     &
        register_diag_field( module_name, 'so2_emis_cmip2', axes(1:2), &
                             Time, 'so2_emis_cmip2', 'kg/m2/s')
   id_nh3_emis_cmip2 =     &
        register_diag_field( module_name, 'nh3_emis_cmip2', axes(1:2), &
                            Time, 'nh3_emis_cmip2', 'kg/m2/s')  
!--for Ox(jmao,1/1/2011)
   id_prodox = register_diag_field( module_name, 'Ox_prod', axes(1:3), &
        Time, 'Ox_prod','VMR/s')
   id_lossox = register_diag_field( module_name, 'Ox_loss', axes(1:3), &
        Time, 'Ox_loss','VMR/s')

   do i=1,pcnstm1
      id_chem_tend(i) = register_diag_field( module_name, trim(tracnam(i))//'_chem_dt', axes(1:3), &
                                             Time, trim(tracnam(i))//'_chem_dt','VMR/s' )
      id_prod(i) = register_diag_field( module_name, trim(tracnam(i))//'_prod', axes(1:3), &
                                        Time, trim(tracnam(i))//'_prod','VMR/s')
      id_loss(i) = register_diag_field( module_name, trim(tracnam(i))//'_loss', axes(1:3), &
                                        Time, trim(tracnam(i))//'_loss','VMR/s')
      if( has_emis(i) ) then
         id_emis(i) = register_diag_field( module_name, trim(tracnam(i))//'_emis', axes(1:2), &
                                           Time, trim(tracnam(i))//'_emis', 'molec/cm2/s')
      else
         id_emis(i) = 0
      end if
      if( has_emis3d(i) ) then
         id_emis3d(i) = register_diag_field( module_name, trim(tracnam(i))//'_emis3d', axes(1:3), &
                                             Time, trim(tracnam(i))//'_emis3d', 'molec/cm2/s')
      else
         id_emis3d(i) = 0
      end if
!     if( has_xactive_emis(i) ) then
         select case (trim(tracnam(i)))
         case ('ISOP')
               id_glaiage = register_diag_field( module_name, 'gamma_lai_age', axes(1:2), &
                    Time, 'gamma_lai_age', 'unitless' )
               id_gtemp = register_diag_field( module_name, 'gamma_temp', axes(1:2), &
                    Time, 'gamma_temp', 'unitless' )
               id_glight = register_diag_field( module_name, 'gamma_light', axes(1:2), &
                    Time, 'gamma_light', 'unitless' )
               id_tsfc = register_diag_field( module_name, 'tsfcair', axes(1:2), &
                    Time, 'tsfcair', 'K' )
               id_fsds = register_diag_field( module_name, 'fsdvd', axes(1:2), &
                    Time, 'fsdvd', 'W/m2' )
               id_ctas = register_diag_field( module_name, 'clim_tas', axes(1:2), &
                    Time, 'clim_tas', 'K' )
               id_cfsds = register_diag_field( module_name, 'clim_fsds', axes(1:2), &
                    Time, 'clim_fsds', 'umol/m2/s PAR' )
         case default
         end select
!     else
!        id_xactive_emis(i) = 0
!     end if
      if( has_ubc(i) ) then
         id_ub(i) = register_diag_field( module_name, trim(tracnam(i))//'_up', axes(1:3), &
                                         Time, trim(tracnam(i))//'_up','VMR' )
      else
         id_ub(i) = 0
      end if
      if( has_lbc(i) ) then
         id_lb(i) = register_diag_field( module_name, trim(tracnam(i))//'_lbc', &
                                         Time, trim(tracnam(i))//'_lbc','VMR' )
      else
         id_lb(i) = 0
      end if
      if( has_airc(i) ) then
         id_airc(i) = register_diag_field( module_name, trim(tracnam(i))//'_airc_emis', axes(1:3), &
                                           Time, trim(tracnam(i))//'_airc_emis','molec/cm2/s' )
      else
         id_airc(i) = 0
      end if
   end do

!-----------------------------------------------------------------------
!     ... Register diagnostic fields for photolysis rates
!-----------------------------------------------------------------------
   do i=1,phtcnt
      write(fld,'(''jval_'',I3.3,8x)') i
      id_jval(i) = register_diag_field( module_name, TRIM(fld), axes(1:3), Time, TRIM(fld),'1/s')
   end do

!-----------------------------------------------------------------------
!     ... Register diagnostic fields for kinetic rate constants
!-----------------------------------------------------------------------
   do i=1,gascnt
      write(fld,'(''k_rxn'',I3.3,8x)') i
      id_rate_const(i) = register_diag_field( module_name, TRIM(fld), axes(1:3), Time, TRIM(fld),'cm3/molec/s')
   end do

!-----------------------------------------------------------------------
!     ... initialize time_interp
!-----------------------------------------------------------------------
   call time_interp_init
      

!-----------------------------------------------------------------------
!     ... initialize esfsw_parameters
!-----------------------------------------------------------------------
   call esfsw_parameters_init

!-----------------------------------------------------------------------
!     ... initialize mpp clock id
!-----------------------------------------------------------------------
!  clock_id = mpp_clock_id('Chemistry')
      
   module_is_initialized = .true.
      
      
!-----------------------------------------------------------------------
      
end function tropchem_driver_init
!</FUNCTION>
 
!#####################################################################

subroutine tropchem_driver_time_vary (Time)

type(time_type), intent(in) :: Time

      integer :: yr, mo,day, hr,min, sec
      integer :: n

      do n=1, size(inter_emis,1)
        if (has_emis(n)) then
          call obtain_interpolator_time_slices (inter_emis(n), Time)
        endif
      end do

      do n=1, size(inter_emis3d,1)
        if (has_emis3d(n)) then
          call obtain_interpolator_time_slices (inter_emis3d(n), Time)
        endif
      end do

      do n=1, size(inter_aircraft_emis,1)
        if (has_airc(n)) then
          call obtain_interpolator_time_slices   &
                                         (inter_aircraft_emis(n), Time)
        endif
      end do

      call obtain_interpolator_time_slices (conc, Time)

      call obtain_interpolator_time_slices (sulfate, Time)

      do n=1, size(ub,1)
        if (has_ubc(n)) then
          call obtain_interpolator_time_slices (ub(n), Time)
        endif
      end do

      call strat_chem_dcly_dt_time_vary (Time)

!----------------------------------------------------------------------
!    determine if this time step starts a new month; if so, then
!    new interactive isoprene emission data is needed, and the 
!    necessary flag is set.
!----------------------------------------------------------------------
      do n=1,pcnstm1
        if ( has_xactive_emis(n) .or. id_xactive_emis(n)>0 ) then
          select case (trim(tracnam(n)))
            case ('ISOP')
              call get_date(Time,yr,mo,day,hr,min,sec)  !model GMT
              newmonth = ( mo /= isop_oldmonth )
              if (newmonth) then
                isop_oldmonth = mo
              endif
            case ('DMS')
              call atmos_sulfate_time_vary (Time)
            case default
          end select
        endif
      end do
             

end subroutine tropchem_driver_time_vary 
      



!#####################################################################

subroutine tropchem_driver_endts


      integer :: n

      do n=1, size(inter_emis,1)
        if (has_emis(n)) then
          call unset_interpolator_time_flag(inter_emis(n))
         endif
      end do

      do n=1, size(inter_emis3d,1)
        if (has_emis3d(n)) then
         call unset_interpolator_time_flag(inter_emis3d(n))
        endif
      end do

      do n=1, size(inter_aircraft_emis,1)
        if (has_airc(n)) then
          call unset_interpolator_time_flag(inter_aircraft_emis(n))
        endif
      end do

      call unset_interpolator_time_flag(conc)           
      call unset_interpolator_time_flag(sulfate)         

      do n=1, size(ub,1)
        if (has_ubc(n)) then
          call unset_interpolator_time_flag(ub(n))
        endif
      end do

      call strat_chem_dcly_dt_endts               



end subroutine tropchem_driver_endts


!######################################################################

subroutine tropchem_driver_end

!-----------------------------------------------------------------------
!     ... initialize mpp clock id
!-----------------------------------------------------------------------
      
   module_is_initialized = .false.
      
      
!-----------------------------------------------------------------------
      
end subroutine tropchem_driver_end

!#######################################################################

! <SUBROUTINE NAME="read_2D_emis_data">
!   <OVERVIEW>
!     Read emissions file
!   </OVERVIEW>
!   <DESCRIPTION>
!     Reads tracer surface emissions from a NetCDF file
!   </DESCRIPTION>
!   <TEMPLATE>
!     call read_2D_emis_data( emis_type, emis, Time, &
!                             field_names, &
!                             Ldiurnal, coszen, half_day, lon, &
!                             is, js, id_emis_diag ) 
!   </TEMPLATE>

subroutine read_2D_emis_data( emis_type, emis, Time, Time_next, &
                              field_names, &
                              Ldiurnal, coszen, half_day, lon, &
                              is, js, id_emis_diag )
    
   type(interpolate_type),intent(inout) :: emis_type
   real, dimension(:,:),intent(out) :: emis
   type(time_type),intent(in) :: Time, Time_next
   character(len=*),dimension(:), intent(in) :: field_names
   logical, intent(in) :: Ldiurnal
   real, dimension(:,:), intent(in) :: coszen, half_day, lon
   integer, intent(in) :: is, js
   integer, intent(in),optional :: id_emis_diag ! id for diagnostic


   integer :: i, j, k
   logical :: used
   real, dimension(size(emis,1),size(emis,2)) :: temp_data
   real :: diurnal_scale_factor, gmt, iso_on, iso_off, dayfrac
   real :: local_angle, factor_tmp

   emis(:,:) = 0.
   temp_data(:,:) = 0.
   do k = 1,size(field_names)
      call interpolator(emis_type,Time,temp_data,field_names(k),is,js)
      emis(:,:) = emis(:,:) + temp_data(:,:)
   end do

   if (Ldiurnal) then
      do j=1,size(emis,2)
      do i=1,size(emis,1)
         if( coszen(i,j) < 0. ) then
            diurnal_scale_factor = 0.
         else
            iso_off = .8 * half_day(i,j)
            iso_on  = -iso_off
            dayfrac = iso_off/PI
            gmt = universal_time(Time)
            local_angle = gmt + lon(i,j) - PI
            if (local_angle >= PI) local_angle = local_angle - twopi
            if (local_angle < -PI) local_angle = local_angle + twopi
            if( local_angle >= iso_off .or. local_angle <= iso_on ) then
               diurnal_scale_factor = 0.
            else
               factor_tmp = local_angle - iso_on
               factor_tmp = factor_tmp / MAX(2.*iso_off,1.e-6)
               diurnal_scale_factor = 2. / dayfrac * (sin(PI*factor_tmp))**2
            end if
         end if
         emis(i,j) = emis(i,j) * diurnal_scale_factor
      end do
      end do
   end if

   if (present(id_emis_diag)) then
      if (id_emis_diag > 0) then
         used = send_data(id_emis_diag,emis,Time_next,is_in=is,js_in=js)
      end if
   end if
end subroutine read_2D_emis_data
!</SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="read_3D_emis_data">
!   <OVERVIEW>
!     Read emissions file
!   </OVERVIEW>
!   <DESCRIPTION>
!     Reads tracer 3-D emissions from a NetCDF file
!   </DESCRIPTION>
!   <TEMPLATE>
!     call read_3D_emis_data( emis_type, emis, Time, phalf, &
!                             field_names, &
!                             Ldiurnal, coszen, half_day, lon, &
!                             is, js, id_emis_diag ) 
!   </TEMPLATE>

subroutine read_3D_emis_data( emis_type, emis, Time, Time_next, phalf, &
                              field_names, &
                              Ldiurnal, coszen, half_day, lon, &
                              is, js, id_emis_diag )
    
   type(interpolate_type),intent(inout) :: emis_type
   real, dimension(:,:,:),intent(in) :: phalf
   real, dimension(:,:,:),intent(out) :: emis
   type(time_type),intent(in) :: Time, Time_next
   character(len=*),dimension(:), intent(in) :: field_names
   logical, intent(in) :: Ldiurnal
   real, dimension(:,:), intent(in) :: coszen, half_day, lon
   integer, intent(in) :: is, js
   integer, intent(in),optional :: id_emis_diag ! id for diagnostic


   integer :: i, j, k
   logical :: used
   real, dimension(size(emis,1),size(emis,2),size(emis,3)) :: temp_data
   real :: diurnal_scale_factor, gmt, iso_on, iso_off, dayfrac
   real :: local_angle, factor_tmp

   emis(:,:,:) = 0.
   temp_data(:,:,:) = 0.
   do k = 1,size(field_names)
      call interpolator(emis_type,Time,phalf,temp_data,field_names(k),is,js)
      emis(:,:,:) = emis(:,:,:) + temp_data(:,:,:)
   end do
   if (Ldiurnal) then
      do j=1,size(emis,2)
      do i=1,size(emis,1)
         if( coszen(i,j) < 0. ) then
            diurnal_scale_factor = 0.
         else
            iso_off = .8 * half_day(i,j)
            iso_on  = -iso_off
            dayfrac = iso_off/PI
            gmt = universal_time(Time)
            local_angle = gmt + lon(i,j) + PI
            if (local_angle >= PI) local_angle = local_angle - twopi
            if (local_angle < -PI) local_angle = local_angle + twopi
            if( local_angle >= iso_off .or. local_angle <= iso_on ) then
               diurnal_scale_factor = 0.
            else
               factor_tmp = local_angle - iso_on
               factor_tmp = factor_tmp / MAX(2.*iso_off,1.e-6)
               diurnal_scale_factor = 2. / dayfrac * (sin(PI*factor_tmp))**2
            end if
         end if
         emis(i,j,:) = emis(i,j,:) * diurnal_scale_factor
      end do
      end do
   end if

   if (present(id_emis_diag)) then
      if (id_emis_diag > 0) then
         used = send_data(id_emis_diag,emis,Time_next,is_in=is,js_in=js)
      end if
   end if
end subroutine read_3D_emis_data
!</SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="calc_xactive_emis">
!   <OVERVIEW>
!     Calculates interactive emissions
!   </OVERVIEW>
!   <DESCRIPTION>
!     Calculates interactive emissions
!   </DESCRIPTION>
!   <TEMPLATE>
!     call calc_xactive_emis( index, emis, Time, is, js, id_emis_diag ) 
!   </TEMPLATE>

subroutine calc_xactive_emis( index, Time, Time_next, lon, lat, pwt, is, ie, js, je, &
                              area, land, ocn_flx_fraction, tsurf, w10m, emis, &
                              kbot, id_emis_diag )
    
   integer,intent(in) :: index
   type(time_type),intent(in) :: Time, Time_next
   real, intent(in), dimension(:,:) :: lon, lat
   real, intent(in), dimension(:,:,:) :: pwt
   integer, intent(in) :: is, ie, js, je
   real, intent(in), dimension(:,:) :: area    ! grid box area (m^2)
   real, intent(in), dimension(:,:) :: land    ! land fraction
   real, intent(in), dimension(:,:) :: ocn_flx_fraction 
   real, intent(in), dimension(:,:) :: tsurf   ! surface temperature (K)
   real, intent(in), dimension(:,:) :: w10m    ! wind speed at 10m (m/s)
   real, dimension(:,:,:),intent(out) :: emis  ! VMR/s
   integer, intent(in), dimension(:,:), optional :: kbot
   integer, intent(in),optional :: id_emis_diag ! id for diagnostic

   logical :: used

   
   if (index == dms_ndx) then
      call atmos_DMS_emission( lon, lat, area, ocn_flx_fraction, tsurf, w10m, pwt, &
                               emis, Time, Time_next, is, ie, js, je, kbot )
   else
      call error_mesg ('calc_xactive_emis', &
                       'Interactive emissions not defined for species: '//trim(tracnam(index)), FATAL)
   end if

   if (present(id_emis_diag)) then
      if (id_emis_diag > 0) then
         used = send_data( id_emis_diag, emis, Time_next, is_in=is, js_in=js)
      end if
   end if
end subroutine calc_xactive_emis
!</SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="init_emis_data">
!   <OVERVIEW>
!     Open emissions file
!   </OVERVIEW>
!   <DESCRIPTION>
!     Opens NetCDF file of tracer surface emissions for reading,
!     and set up interpolation to model grid/time
!   </DESCRIPTION>
!   <TEMPLATE>
!     call init_emis_data( emis_type, model, method_type, index, file_name, &
!                          lonb_mod, latb_mod, field_type, flag, diurnal )
!   </TEMPLATE>

subroutine init_emis_data( emis_type, model, method_type, index, file_name, &
                           lonb_mod, latb_mod, field_type, flag, diurnal, &
                           axes, Time )
    
   type(interpolate_type),intent(inout) :: emis_type
   integer, intent(in) :: model,index
   character(len=*),intent(in) :: method_type
   character(len=*),intent(inout) ::file_name
   real,intent(in),dimension(:,:) :: lonb_mod,latb_mod
   type(field_init_type),intent(out) :: field_type
   logical, intent(out) :: flag, diurnal
   integer        , intent(in)  :: axes(4)
   type(time_type), intent(in)  :: Time
    
   character(len=64) :: name, control
   integer :: nfields
   integer :: flag_name, flag_file, flag_diurnal
   character(len=64) :: emis_name, emis_file, control_diurnal

   flag = .false.
   diurnal = .false.
   control = ''
   if( query_method(trim(method_type),model,index,name,control) ) then
      if( trim(name) == 'file' ) then
         flag = .true.
         flag_file = parse(control, 'file', emis_file)
         flag_name = parse(control, 'name', emis_name)
         flag_diurnal = parse(control, 'diurnal', control_diurnal)
         if(flag_file > 0) then
            file_name = emis_file
         else if (flag_name > 0) then
            select case (trim(method_type))
               case ('emissions3d')
                  file_name  = trim(file_emis3d_1)//trim(emis_name)//trim(file_emis3d_2)
               case default
                  file_name  = trim(file_emis_1)//trim(emis_name)//trim(file_emis_2)
            end select
         end if
         diurnal = (flag_diurnal > 0)

         call interpolator_init( emis_type, trim(file_name), &
                                 lonb_mod, latb_mod,  &
                                 data_out_of_bounds=(/CONSTANT/), &
                                 vert_interp=(/INTERP_WEIGHTED_P/) )
         call query_interpolator(emis_type,nfields=nfields)
         allocate(field_type%field_names(nfields))
         call query_interpolator(emis_type,field_names=field_type%field_names)
      end if
   end if
end subroutine init_emis_data
!</SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="init_xactive_emis">
!   <OVERVIEW>
!     Set up interactive emission calculations
!   </OVERVIEW>
!   <DESCRIPTION>
!     Set up interactive emission calculations
!   </DESCRIPTION>
!   <TEMPLATE>
!     call init_xactive_emis( model, method_type, index, species, &
!                             axes, Time, lonb_mod, latb_mod, phalf, &
!                             flag, mask )
!   </TEMPLATE>

subroutine init_xactive_emis( model, method_type, index, species, &
                              axes, Time, lonb_mod, latb_mod, phalf, &
                              flag, id_xemis, mask )
    
   integer,         intent(in)  :: model, index
   character(len=*),intent(in)  :: method_type, species
   integer        , intent(in)  :: axes(4)
   type(time_type), intent(in)  :: Time
   real,            intent(in), dimension(:,:)   :: lonb_mod,latb_mod
   real,            intent(in), dimension(:,:,:) :: phalf
   real,            intent(in), dimension(:,:,:), optional :: mask
   logical,         intent(out) :: flag
   integer,         intent(out) :: id_xemis
    
   character(len=64) :: name, control
   integer :: nhalf, nfull

   flag = .false.
   control = ''
   nhalf = SIZE(phalf,3)
   nfull = nhalf - 1

   flag = query_method(trim(method_type),model,index,name,control)
   
!  if (flag) then
      select case (trim(species))
         case ('DMS')
            id_xemis = &
               register_diag_field( module_name, trim(species)//'_xactive_emis', axes(1:3), &
                                    Time, trim(species)//'_xactive_emis', 'VMR/s')
            if (flag .or. id_xemis>0) then
               call atmos_sulfate_init( lonb_mod, latb_mod, nfull, axes, Time, mask )
            end if
         case ('ISOP')
            id_xemis = &
               register_diag_field( module_name, trim(species)//'_xactive_emis', axes(1:2), &
                                    Time, trim(species)//'_xactive_emis', 'molecules/cm2/s')
            if (flag .or. id_xemis>0) then
               call isop_xactive_init( lonb_mod, latb_mod, axes )
            end if
         case default
            if (flag) then
               call error_mesg ('init_xactive_emis','Interactive emissions not defined for species: '//trim(species), FATAL)
            else
               id_xemis = 0
            end if
      end select
!  end if

end subroutine init_xactive_emis
!</SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="isop_xactive_init">
!   <OVERVIEW>
!     Initialize interactive isoprene emissions
!   </OVERVIEW>
!   <DESCRIPTION>
!     Initialize interactive isoprene emissions: read in emission capacities,
!     LAI for each pft, percentage pft in each grid cell, climatological (1980-2000)
!     monthly mean surface air temp + downward short wave radiation from 
!     Sheffield et al., J. clim, 2006  (obtained from Christine Wiedinmyer, NCAR,
!     Feb, 2009.  Also allocate and initialize arrays needed for diagnostics, etc.
!     Updated to MEGAN v2.1 amf/van
!   </DESCRIPTION>
!   <TEMPLATE>
!     call isop_xactive_init( lonb_mod, latb_mod, axes )
!   </TEMPLATE>

subroutine isop_xactive_init( lonb, latb, axes )

  real, intent(in), dimension(:,:)   :: lonb,latb
  integer        , intent(in)  :: axes(4)
  integer :: nlon, nlat, i, m
!parameters for input file
!number of vegetation types for emission capacities

!++amf/van updated to megan2.1 nveg = 5 types above
! integer, parameter :: nveg=6, npft=17, nmos=12
!--amf/van
  integer, parameter :: nlonin = 720, nlatin = 360
  integer, parameter :: metlonin = 360, metlatin= 180
  real, dimension(nlonin) :: inlon, lonpft, lonlai
  real, dimension(nlatin) :: inlat, latpft, latlai
  real, dimension(nlonin+1) :: inlone, lonpfte, lonlaie
  real, dimension(nlatin+1) :: inlate, latpfte, latlaie
  real, dimension(metlonin) :: metlon
  real, dimension(metlatin) :: metlat
  real, dimension(metlonin+1) :: metlone
  real, dimension(metlatin+1) :: metlate
  integer, dimension(npft) :: pft
  integer, dimension(nmos) :: mos
  real :: edgew,edges,edgen,edgee,dlon,dlat
  character(len=19), parameter :: ecfile = 'INPUT/megan.ISOP.nc'       ! amf/van
  character(len=25), parameter :: laifile = 'INPUT/mksrf_lai.060929.nc'
  character(len=25), parameter :: pftfile = 'INPUT/mksrf_pft.060929.nc'
  character(len=35), parameter :: tasfile = 'INPUT/tas_monthly_clim_1980-2000.nc'
  character(len=37), parameter :: dswfile = 'INPUT/dswrf_monthly_clim_1980-2000.nc'
  character(len=3) :: vegnames(nveg) = (/ 'ntr', 'btr', 'crp', 'grs', 'shr' /)    !amf/van
  character(len=5) :: pftnames(npft) = (/ 'pft01','pft02','pft03','pft04','pft05','pft06', &
       'pft07','pft08','pft09','pft10','pft11','pft12','pft13','pft14','pft15','pft16','pft17' /)
  character(len=5) :: lainames(npft) = (/ 'lai01','lai02','lai03','lai04','lai05','lai06', &
       'lai07','lai08','lai09','lai10','lai11','lai12','lai13','lai14','lai15','lai16','lai17' /)
  character(len=5) :: tasnames(nmos) = (/ 'tas01','tas02','tas03','tas04','tas05','tas06', &
       'tas07','tas08','tas09','tas10','tas11','tas12' /)
  character(len=5) :: dswnames(nmos) = (/ 'dsw01','dsw02','dsw03','dsw04','dsw05','dsw06', &
       'dsw07','dsw08','dsw09','dsw10','dsw11','dsw12' /)
  integer :: id_ec(nveg), id_pft(npft), id_lai(npft), id_tas(nmos), id_dsw(nmos)  
  real, dimension(nlonin,nlatin) :: datain
  real, dimension(nlonin,nlatin,npft) :: datapft
  real, dimension(nlonin,nlatin,npft,nmos) :: datalai
  real, dimension(metlonin,metlatin,nmos) :: tas, dswrf
  logical :: used

  nlon = size(lonb,1)-1
  nlat = size(latb,2)-1

!allocate dimensions of array for storing emission capacities (model grid)
  allocate(ecisop(nlon,nlat,nveg))
  allocate(pctpft(nlon,nlat,npft))
  allocate(mlai(nlon,nlat,npft,nmos))      ! monthly mean LAI
  allocate(emisop_month(nlon,nlat))        ! emission capacities adjusted w/ monthly gamma_lai and gamma_age
  allocate(diag_gamma_lai_age(nlon,nlat))
  allocate(diag_gamma_light(nlon,nlat))
  allocate(diag_gamma_temp(nlon,nlat))
  allocate(diag_climtas(nlon,nlat))
  allocate(diag_climfsds(nlon,nlat))
  allocate(ts_avg(nlon,nlat,nmos))
  allocate(fsds_avg(nlon,nlat,nmos))

!Initalize arrays
  emisop_month(:,:) = 0.
  diag_gamma_lai_age(:,:) = 0.
  diag_gamma_light(:,:) = 0.
  diag_gamma_temp(:,:) = 0.
  diag_climtas(:,:) = 0.
  diag_climfsds(:,:) = 0.
! always get gamma age / gamma lai at start of run.
! newmonth = .true.     

!  --- check existence of input file containing isoprene emission capacities --------
  if (file_exist(ecfile)) then
      
!set up for input grid
     if(mpp_pe() == mpp_root_pe()) call error_mesg ('isop_xactive_init',  &
          'Reading NetCDF formatted input file: iso_bvoc.nc', NOTE)
      
! amf/van -- new file only has lat lon centers, not boundaries...
!read in lat & lon boundaries from input file and convert to radians
!no.domain = true tells it all fields are in one global file (as opposed to e.g., one per cube face)
     call read_data (ecfile, 'lon', inlon, no_domain=.true.)
     call read_data (ecfile, 'lat', inlat, no_domain=.true.)
     inlon = inlon*DEG_TO_RAD
     inlat = inlat*DEG_TO_RAD
     dlat = inlat(2)-inlat(1)
     dlon = inlon(2)-inlon(1)
     inlone(1:nlonin) = inlon-(dlon/2.)
     inlone(nlonin+1) = inlon(nlonin)+(dlon/2.)
     inlate(1:nlatin) = inlat-(dlat/2.)
     inlate(nlatin+1) = inlat(nlatin)+(dlat/2.)

!     print*, 'lat,lon edges for new megan.isop.nc file: ', inlate, inlone
     
     call horiz_interp_init
     call horiz_interp_new ( Interp, inlone, inlate, lonb, latb )

! loop over vegnames
     do i = 1, nveg 
!register diagnostic field
        id_ec(i) = register_static_field( 'tracers', vegnames(i), axes(1:2), &
             vegnames(i), 'molec/cm2/s')
!read this field
        call read_data (ecfile, vegnames(i), datain, no_domain=.true.)
        call horiz_interp (Interp, datain, ecisop(:,:,i), verbose=verbose)
!send data to diagnostic
        if (id_ec(i) > 0) then
           used = send_data(id_ec(i),ecisop(:,:,i)) 
        endif
     end do
    
  else 
     print*, 'what is ecfile ', ecfile
     call error_mesg ('isop_xactive_init',  &
          'ecfile does not exist.', FATAL)
  endif

!  --- check existence of input file containing % pft --------
  if (file_exist(pftfile)) then
      
!set up for input grid
     if(mpp_pe() == mpp_root_pe()) call error_mesg ('isop_xactive_init',  &
          'Reading NetCDF formatted input file: mksrf_pft.060929.nc', NOTE)
      
!read in lat & lon from input file, get boundaries and convert to radians
     call read_data (pftfile, 'lon', lonpft, no_domain=.true.)
     call read_data (pftfile, 'lat', latpft, no_domain=.true.)
     call read_data (pftfile, 'EDGEW', edgew, no_domain=.true.)
     call read_data (pftfile, 'EDGES', edges, no_domain=.true.)
     call read_data (pftfile, 'EDGEE', edgee, no_domain=.true.)
     call read_data (pftfile, 'EDGEN', edgen, no_domain=.true.)

     lonpfte(1) = edgew
     latpfte(1) = edges
     latpfte(nlatin+1) = edgen
     lonpfte(nlonin+1) = edgee
     
     dlon = 2.*(lonpft(1)-edgew)
     dlat = 2.*(latpft(1)-edges)

     do i = 2, nlatin 
        latpfte(i) = latpfte(i-1)+dlat
     end do

     do i = 2, nlonin
        lonpfte(i) = lonpfte(i-1)+dlon
     end do     

     lonpfte = lonpfte*DEG_TO_RAD
     latpfte = latpfte*DEG_TO_RAD
     
     call horiz_interp_init
     call horiz_interp_new ( Interp, lonpfte, latpfte, lonb, latb )

     call read_data (pftfile, 'pft', pft, no_domain=.true.)

!read pct_pft field 
     call read_data (pftfile, 'PCT_PFT', datapft, no_domain=.true.)

! loop over vegnames
     do i = 1, npft       
!register diagnostic field
        id_pft(i) = register_static_field( 'tracers', pftnames(i), axes(1:2), &
             pftnames(i), 'unitless')
        call horiz_interp (Interp, datapft(:,:,i), pctpft(:,:,i), verbose=verbose)
!send data to diagnostic
        if (id_pft(i) > 0) then
           used = send_data(id_pft(i),pctpft(:,:,i)) 
        endif
     end do
!scale pctpft from percentage to fraction 
     pctpft(:,:,:) = .01 * pctpft(:,:,:)
    
  else 
     call error_mesg ('isop_xactive_init',  &
          'pftfile does not exist.', FATAL)
  endif

!  --- check existence of input file containing monthly lai, for each pft --------
  if (file_exist(laifile)) then
      
!set up for input grid
     if(mpp_pe() == mpp_root_pe()) call error_mesg ('isop_xactive_init',  &
          'Reading NetCDF formatted input file: mksrf_lai.060929.nc', NOTE)
      
!read in lat & lon from input file, get boundaries and convert to radians
     call read_data (laifile, 'lon', lonlai, no_domain=.true.)
     call read_data (laifile, 'lat', latlai, no_domain=.true.)

! get lat/lon edges
     lonlaie(1) = edgew
     latlaie(1) = edges
     latlaie(nlatin+1) = edgen
     lonlaie(nlonin+1) = edgee
     
     dlon = 2.*(lonlai(1)-edgew)
     dlat = 2.*(latlai(1)-edges)

     do i = 2, nlatin 
        latlaie(i) = latlaie(i-1)+dlat
     end do

     do i = 2, nlonin
        lonlaie(i) = lonlaie(i-1)+dlon
     end do     

     lonlaie = lonlaie*DEG_TO_RAD
     latlaie = latlaie*DEG_TO_RAD
     
     call horiz_interp_init
     call horiz_interp_new ( Interp, lonlaie, latlaie, lonb, latb )

! read in pft and time dimensions from lai file
     call read_data (laifile, 'time', mos, no_domain=.true.)

! loop over vegnames
     do i = 1, npft 
        call read_data (laifile,lainames(i),datalai(:,:,i,:), no_domain=.true.)
        do m = 1, nmos
            call horiz_interp (Interp, datalai(:,:,i,m), mlai(:,:,i,m), verbose=verbose)
!store diagnostics for one month only - choose July for now
            if (m .eq. 7) then 
!register diagnostic field
               id_lai(i) = register_static_field( 'tracers', lainames(i), axes(1:2), &
                    lainames(i), 'unitless')
!send data to diagnostic
               if (id_lai(i) > 0) then
                  used = send_data(id_lai(i),mlai(:,:,i,m))
               endif
            endif
         end do
      end do
    
   else 
      call error_mesg ('isop_xactive_init',  &
           'laifile does not exist', FATAL)
   endif

!  --- check existence of input file containing climatological (1980-2000) monthly air surface temp --------
  if (file_exist(tasfile)) then
      
!set up for input grid
     if(mpp_pe() == mpp_root_pe()) call error_mesg ('isop_xactive_init',  &
          'Reading NetCDF formatted input file: tas_monthly_clim_1980-2000.nc', NOTE)
      
!read in lat & lon from input file, get boundaries and convert to radians
     call read_data (tasfile, 'lon', metlon, no_domain=.true.)
     call read_data (tasfile, 'lat', metlat, no_domain=.true.)
     
     dlon = 0.5*(metlon(1)-metlon(2))
     dlat = 0.5*(metlat(2)-metlat(1))

     do i = 1, metlatin 
        metlate(i) = metlat(i)-dlat
     end do

     metlate(metlatin+1) = metlat(metlatin)+dlat

     do i = 1, metlonin
        metlone(i) = metlon(i)-dlon
     end do     
     metlone(metlonin+1) = metlon(metlonin)+dlon

     metlone = metlone*DEG_TO_RAD
     metlate = metlate*DEG_TO_RAD
     
     call horiz_interp_init
     call horiz_interp_new ( Interp, metlone, metlate, lonb, latb )

     call read_data (tasfile, 'time', mos, no_domain=.true.) 
     call read_data (tasfile,'tas_clim', tas(:,:,:), no_domain=.true.)
        
     do m = 1, nmos
        call horiz_interp (Interp, tas(:,:,m), ts_avg(:,:,m), verbose=verbose) 
!register diagnostic field
        id_tas(m) = register_static_field( 'tracers', tasnames(m), axes(1:2), &
             tasnames(m), 'unitless')
!send data to diagnostic
        if (id_tas(m) > 0) then
           used = send_data(id_tas(m),ts_avg(:,:,m))
        endif
     end do
  else 
     call error_mesg ('isop_xactive_init',  &
          tasfile, NOTE)
     call error_mesg ('isop_xactive_init',  &
          'tasfile does not exist', FATAL)

  endif

!  --- check existence of input file containing climatological (1980-2000) monthly surface down SW radiation --------
  if (file_exist(dswfile)) then
      
!set up for input grid
     if(mpp_pe() == mpp_root_pe()) call error_mesg ('isop_xactive_init',  &
          'Reading NetCDF formatted input file: dswrf_monthly_clim_1980-2000.nc', NOTE)
      
!read in lat & lon from input file, get boundaries and convert to radians
     call read_data (dswfile, 'lon', metlon, no_domain=.true.)
     call read_data (dswfile, 'lat', metlat, no_domain=.true.)
     
     dlon = 0.5*(metlon(1)-metlon(2))
     dlat = 0.5*(metlat(2)-metlat(1))

     do i = 1, metlatin 
        metlate(i) = metlat(i)-dlat
     end do

     metlate(metlatin+1) = metlat(metlatin)+dlat

     do i = 1, metlonin
        metlone(i) = metlon(i)-dlon
     end do     
     metlone(metlonin+1) = metlon(metlonin)+dlon

     metlone = metlone*DEG_TO_RAD
     metlate = metlate*DEG_TO_RAD
     
     call horiz_interp_init
     call horiz_interp_new ( Interp, metlone, metlate, lonb, latb )

     call read_data (dswfile, 'time', mos, no_domain=.true.)
     call read_data (dswfile,'dswrf_clim', dswrf(:,:,:), no_domain=.true.)
        
     do m = 1, nmos
        call horiz_interp (Interp, dswrf(:,:,m), fsds_avg(:,:,m), verbose=verbose)
!register diagnostic field
        id_dsw(m) = register_static_field( 'tracers', dswnames(m), axes(1:2), &
             dswnames(m), 'unitless')
!send data to diagnostic
        if (id_dsw(m) > 0) then
           used = send_data(id_dsw(m),fsds_avg(:,:,m))
        endif
     end do
    
  else 
     call error_mesg ('isop_xactive_init',  &
          'dswfile does not exist', FATAL)
  endif

end subroutine isop_xactive_init

!</SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="calc_xactive_isop">
!   <OVERVIEW>
!     Calculates interactive isoprene emissions
!   </OVERVIEW>
!   <DESCRIPTION>
!     Calculates interactive isoprene emissions using algorithms from 
!      PCEEA MEGAN model in Guenther, ACP, 2006.
!      Note - gamma soil moisture is assumed constant (at one)
!   </DESCRIPTION>
!   <TEMPLATE>
!     call calc_xactive_isop( index, Time, lon, lat, oro, pwtsfc, is, js, &
!                              area, land, tsfcair, flux_sw_down_vis, &
!                              coszen, emis,   id_gamma_lai_age, &
!                                 id_gamma_temp, id_gamma_light, &
!                                 id_tsfcair, id_fsdvd, id_climtas, id_climfsds, id_emis_diag ) 
!   </TEMPLATE>
             
subroutine calc_xactive_isop( index, Time, Time_next,lon, lat, oro, pwtsfc, is, js, &
                              area, land, tsfcair, flux_sw_down_vis, &
                              coszen, emis,   id_gamma_lai_age, &
                              id_gamma_temp, id_gamma_light, id_tsfcair, &
                              id_fsdvd, id_climtas, id_climfsds, id_emis_diag ) 
    
   integer,intent(in) :: index
   type(time_type),intent(in) :: Time, Time_next
   real, intent(in), dimension(:,:) :: lon, lat
   real, intent(in), dimension(:,:) :: pwtsfc
   integer, intent(in) :: is, js 
   real, intent(in), dimension(:,:) :: area    ! grid box area (m^2)
   real, intent(in), dimension(:,:) :: land    ! land fraction
   real, intent(in), dimension(:,:) :: oro     ! land = 1; ocean = 0
   real, intent(in), dimension(:,:) :: tsfcair   ! surface temperature (K)
   real, intent(in), dimension(:,:) :: flux_sw_down_vis !W/m2 visible (direct+diffuse) sfc flux
   real, intent(in), dimension(:,:) :: coszen  ! cosine of solar zenith angle
   real, dimension(:,:),intent(out) :: emis  ! 
   integer, intent(in), optional :: id_gamma_lai_age, id_gamma_temp, id_gamma_light
   integer, intent(in), optional :: id_climtas, id_climfsds
   integer, intent(in), optional :: id_emis_diag, id_tsfcair
   integer, intent(in),optional ::  id_fsdvd ! id for diagnostic
   real    :: calday
   type(time_type) :: Year_t
   integer :: yr, mo, day, hr, min, sec
   logical :: used
   integer :: ie, je

   ie = is + size(land,1) -1
   je = js + size(land,2) -1

   call get_date(Time,yr,mo,day,hr,min,sec)  !model GMT

!update gamma age and gamma once per month
   if (newmonth) then 
!     if( mpp_pe() == mpp_root_pe() ) then 
!        print *, 'pts_proc, oldmonth', pts_processed, isop_oldmonth
!        print *, 'time', yr,mo,day,hr,min,sec
!        print*, 'AMF calc_xactive_isop: calling get_monthly_gammas'
!        print*, 'id_gamma_lai_age = ', id_gamma_lai_age
!        print*, 'sum(diag_gamma_lai_age', sum(diag_gamma_lai_age(:,:))
!     end if
      
      call get_monthly_gammas( lon, lat, oro, is, js, mo, &
                               id_gamma_lai_age )
!     if( mpp_pe() == mpp_root_pe() ) then 
!        print *, 'pts_proc, oldmonth AFTER', pts_processed, isop_oldmonth
!        print*, 'AMF calc_xactive_isop: after call to  get_monthly_gammas'
!        print*, 'id_gamma_lai_age = ', id_gamma_lai_age
!        print*, 'AMF sum(diag_gamma_lai_age)', sum(diag_gamma_lai_age(:,:))
!     end if

   end if

!Get Julian date (fraction) = calday
   Year_t = set_date(yr,1,1,0,0,0)
   calday = time_type_to_real( Time-Year_t) / SECONDS_PER_DAY

!  if( mpp_pe() == mpp_root_pe() ) then
!     print*, 'AMF: calday = ', calday
!  end if

! get isoprene emissions for this timestep
   call get_isop_emis( lon, lat, is, js, calday, mo, tsfcair, flux_sw_down_vis, &
                       coszen, area, pwtsfc, emis )

!accumulate isoprene emissions in diagnostic - units should be molec/cm2/s
   if (present(id_emis_diag)) then
      if (id_emis_diag > 0) then
         used = send_data( id_emis_diag, emis, Time_next, is_in=is, js_in=js)
      end if
   end if

! also store sw visible direct at surface and surface air temperature diagnostics
   if (present(id_fsdvd)) then 
      if (id_fsdvd > 0) then 
         used = send_data( id_fsdvd, flux_sw_down_vis, Time_next, is_in=is, js_in=js)
      end if
   end if

   if (present(id_tsfcair)) then 
      if (id_tsfcair > 0) then 
         used = send_data( id_tsfcair, tsfcair, Time_next, is_in=is, js_in=js)
      end if
   end if

   if (present(id_gamma_light)) then 
      if (id_gamma_light > 0) then 
         used = send_data( id_gamma_light, diag_gamma_light(is:ie,js:je), Time_next, is_in=is, js_in=js)
      end if
   end if

   if (present(id_gamma_temp)) then 
      if (id_gamma_temp > 0) then 
         used = send_data( id_gamma_temp, diag_gamma_temp(is:ie,js:je), Time_next, is_in=is, js_in=js)
      end if
   end if

   if (present(id_climtas)) then 
      if (id_climtas > 0) then 
         used = send_data( id_climtas, diag_climtas(is:ie,js:je), Time_next, is_in=is, js_in=js)
      end if
   end if

   if (present(id_climfsds)) then 
      if (id_climfsds > 0) then 
         used = send_data( id_climfsds, diag_climfsds(is:ie,js:je), Time_next, is_in=is, js_in=js)
      end if
   end if

!store combined diagnostic of gamma_lai * gamma_age, summed over each vegetation type
   if (present(id_gamma_lai_age)) then
      if (id_gamma_lai_age > 0) then
!         if( mpp_pe() == mpp_root_pe() ) then 
!            print*, 'AMF calc_xactive_isop: after call to send_data for glaiage'
!            print*, 'id_gamma_lai_age = ', id_gamma_lai_age
!            print*, 'sum(diag_gamma_lai_age', sum(diag_gamma_lai_age(:,:))
!         end if
         used = send_data( id_gamma_lai_age, diag_gamma_lai_age(is:ie,js:je), Time_next, is_in=is, js_in=js)
              
      end if
   end if

end subroutine calc_xactive_isop

!</SUBROUTINE>



!#######################################################################

! <SUBROUTINE NAME="get_isop_emis">
!   <OVERVIEW>
!     Get isop emissions for this time step
!   </OVERVIEW>
!   <DESCRIPTION>
!     amf Feb 2009
!     This subroutine calculates isoprene emissions according to
!     the MEGAN PCEEA model [Guenther et al., ACP 6, 3181-3210, 2006.]
!     as implemented in mozart4_v4.5 by J.-F. Lamarque and G. Pfister
!     (bvoc_emissions mozart routine)
!   </DESCRIPTION>
!   <TEMPLATE>
!     get_isop_emis( lon, lat, is, js, calday, mo, tsfcair, flux_sw_down_vis, &
!                          coszen, area, pwtsfc, emis )
!   </TEMPLATE>
!   <IN NAME="lon" TYPE="real" DIM="(:,:)">
!      longitude value for each i,j grid cell 
!   </IN>
!   <IN NAME="lat" TYPE="real" DIM="(:,:)">
!      latitude value for each i,j grid cell 
!   </IN>
!   <IN NAME="is" TYPE="integer">
!      beginning i for location of this group of 
!      grid cells within global grid, used in diag manager 
!   </IN>
!   <IN NAME="js" TYPE="integer">
!      beginning j for location of this group of 
!      grid cells within global grid, used in diag manager 
!   </IN>
!   <IN NAME="calday" TYPE="real">
!      Fractional Julian day of year (model GMT)
!   </IN>
!   <IN NAME="mo" TYPE="integer">
!     current month
!   </IN>
!   <IN NAME="tsfcair" TYPE="real" DIM="(:,:)">
!      surface air temperature (K)
!   </IN>
!   <IN NAME="flux_sw_down_vis" TYPE="real" DIM="(:,:)">
!      downward visible shortwave radiation (W/m2)
!   </IN>
!   <IN NAME="coszen" TYPE="real" DIM="(:,:)">
!     Cosine of the solar zenith angle
!   </IN>
!   <IN NAME="area" TYPE="real" DIM="(:,:)">
!     Grid cell area (m^2)
!   </IN>
!   <IN NAME="pwtsfc" TYPE="real" DIM="(:,:)">
!     kg/m2 of air in the surface layer
!   </IN>
!   <IN NAME="emis" TYPE="real" DIM="(:,:)">
!     isoprene emissions for this timestep
!     (molec/cm2/s)
!   </IN>

subroutine get_isop_emis( lon, lat, is, js, calday, mo, tsfcair, flux_sw_down_vis, &
                          coszen, area, pwtsfc, emis )
!-------------------------------------------------------------------------------------
!       ... biogenic voc isoprene emissions
!-------------------------------------------------------------------------------------  

      implicit none

      real, intent(in), dimension(:,:) :: lon, lat
      integer, intent(in) :: is, js, mo  
      real, intent(in)    :: calday          
      real, intent(in), dimension(:,:)  :: coszen     
      real, intent(in), dimension(:,:)  :: tsfcair        ! surface temperature
      real, intent(in), dimension(:,:)  :: flux_sw_down_vis !W/m2 direct visible sfc flux
      real, intent(in), dimension(:,:)  :: area      ! grid box area (m^2)
      real, intent(in), dimension(:,:)  :: pwtsfc    ! kg/m2 air in surface layer
      real, intent(out), dimension(:,:) :: emis      ! output emissions for this time step

!-------------------------------------------------------------------------------------
!       ... local variables
!-------------------------------------------------------------------------------------
      real, parameter :: ctm1   =  80.   ! p 3192 Guenther et al ACP 2006
      real, parameter :: ctm2   = 200.   ! same as above
      real, parameter :: const0 = 4.766    ! to convert W/m2 to micromoles/m2/s for PAR (C. Wiedinmyer, 2/18/09)
      real, parameter :: rho_iso = 0.96  ! to account for deposition to canopy
      real, parameter :: gamma_sm = 1.   ! soil moisture - eventually estimate as f(sm, wilting point)

      integer :: i, j, nlon, nlat
      real    :: ppfd, x, Eopt, Ptoa, phi
      real    :: Topt
      real    :: fac_par, fac_tmp
      real    :: Tdaily, Pdaily 
      real    :: t_diff
    
      nlon = size(lon,1)
      nlat = size(lat,2)
      do j = 1,nlat
         do i = 1,nlon

!-------------------------------------------------------------------------------------
!       ... PAR correction
!                Guenther et al ACP 2006 Eqns 11, 12, 13
!           Note - this does not include separation of sunny/shady - could add.
!-------------------------------------------------------------------------------------
!Currently tests to see if we have read in the climatological surface air & downward shortwave
! could eventually change to use values saved during run (e.g., previous week to month values for
! each grid cell)
!Note the factor of 0.5 to convert from total shortwave to PPFD
! vs. AM3 which has visible component available.
         if( has_ts_avg ) then  
            Pdaily =  fsds_avg(i+is-1,j+js-1,mo) * const0 * 0.5 
            Tdaily  = ts_avg(i+is-1,j+js-1,mo)
         else
            Pdaily = Pdaily_clim * const0 * 0.5
            Tdaily = Tdaily_clim
         end if
         ppfd   = flux_sw_down_vis(i,j) * const0

         Ptoa = 3000. + 99.* cos( twopi*(calday - 10.)/365. )    !Guenther et al Eqn 13


! Guenther eqns 11a/b 
         if (coszen(i,j) <= 0.) then 
            fac_par = 0.   
         else
            phi = ppfd / (coszen(i,j)*Ptoa)   ! Eqn 12

!-------------------------------------------------------------------------------------
!  phi can get > 1 and then fac_par gets negative with the above equation
!  set phi=1 if phi> 1 as recommended by Alex (MZ4 code)
!-------------------------------------------------------------------------------------
            phi = min( phi,1. )   
            !Eqn 11b - note typo in MZ4 - 2.49 instead of 2.46
            fac_par = coszen(i,j)*(2.46 * (1. + .0005 *(Pdaily - 400.))*phi - .9*phi*phi) 
         end if

!-------------------------------------------------------------------------------------
!       ... temperature correction  equation 14
!           Topt from equation 8
!           p. 3192 Guenther et al 2006
!-------------------------------------------------------------------------------------
         t_diff  = Tdaily - Tdaily_clim
         Topt    = 313. + 0.6*t_diff                     !Eqn 8
         x       = (tsfcair(i,j) - Topt)/(tsfcair(i,j)*Topt*.00831)
         Eopt    = 1.75 * exp( .08*t_diff )
!Eqn 14 --  gammaT
         fac_tmp = Eopt * (ctm2 * exp( ctm1*x ))/(ctm2 - ctm1*(1. - exp( ctm2*x )))  

!-------------------------------------------------------------------------------------
!     emisop_month contains regridded potential emissions including 
!                  application of gamma LAI and gamma AGE for this month 
!     ... change units from microg/m2/h to mol/cm2/s
!      units of MZ4 input file incorrectly labled as mol/cm2/s... 
!-------------------------------------------------------------------------------------
         emis(i,j) = emisop_month(i+is-1,j+js-1) * fac_par * fac_tmp *  2.46e8
         diag_gamma_temp(i+is-1,j+js-1) = fac_tmp
         diag_gamma_light(i+is-1,j+js-1) = fac_par
         diag_climtas(i+is-1,j+js-1) = Tdaily
         diag_climfsds(i+is-1,j+js-1) = Pdaily

!        if( mpp_pe() == mpp_root_pe() ) then
!           print*, 'AMF: i,j,gt,gl,emis = ', i,j,fac_tmp,fac_par, emis(i,j)
!           print*, 'AMF: calday = ', calday
!        endif
      end do !lon
   end do !lat

!-------------------------------------------------------------------------------------
!        AMF - apply uniform canopy deposition and soil moisture corrections
!              from Guenther et al 2006 - could eventually parameterize these
!                    (neither included in MZ4)
!-------------------------------------------------------------------------------------
   emis(:,:) =  emis(:,:) * rho_iso * gamma_sm 

end subroutine get_isop_emis
!</SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="get_monthly_gammas">
!   <OVERVIEW>
!     Each month, update isoprene emissions to include scaling from 
!     gamma_age and gamma_lai
!   </OVERVIEW>
!   <DESCRIPTION>
!     amf Feb 2009
!     This subroutine calculates gamma_age and gamma_lai and applies them
!     for each LAI and PFT, aggregating to the 6 vegetation types for the
!     isoprene emission capacities, following
!     the MEGAN PCEEA model [Guenther et al., ACP 6, 3181-3210, 2006.]
!     as implemented in mozart4_v4.5 by J.-F. Lamarque and G. Pfister
!     (interp_map_bvoc mozart routine)
!   </DESCRIPTION>
!   <TEMPLATE>
!     get_monthly_gammas( lon, lat, oro, is, js, &
!                                   month, id_gamma_lai_age )  
!   </TEMPLATE>

!removed fsds since not used... 
subroutine get_monthly_gammas( lon, lat, oro, is, js, &
                                   month, id_gamma_lai_age )

      implicit none

      real, intent(in), dimension(:,:) :: lon, lat
      real, intent(in), dimension(:,:) :: oro     ! land = 1; ocean = 0
      integer, intent(in) :: is, js    
      integer, intent(in)      :: month
      integer, intent(in),optional :: id_gamma_lai_age  !id for diagnostic
  
!-------------------------------------------------------------------------------------
!       ... local variables
!-------------------------------------------------------------------------------------
      integer                  :: i, j, n, nlon, nlat, nl, nu
      integer                  :: ie, je
      integer                  :: pft_li(nveg)
      integer                  :: pft_lu(nveg)
!      real                     :: wrk_area  
      real                     :: total_iso
      real                     :: work_iso(nveg)    !amf/van
      real                     :: ts_wrk
      real                     :: total_work_gamma

!++amf/van
      real                     :: work_gamma(nveg)  ! for diagnostic
      real                     :: lai_fac(npft)  ! npft = 17 (mksrf files)
      real                     :: lai_work(npft,nmos)
      logical                  :: age_work(npft)
!--amf/van

      nlon = size(lon,1)
      nlat = size(lat,2)
      ie = is + nlon -1
      je = js + nlat -1

!       if( mpp_pe() == mpp_root_pe() ) then 
!          print*, 'AMF get_monthly_gammas: made it here'
!       end if
      
!++amf/van
      age_work(:) = .true.
      age_work((/2,3,5,6,10/)) = .false.

! hardwired indices for matching pfts to ecs 
! pfts 2-4   fine leaf 
! pfts 5-9   broadleaf
! pfts 10-12 shrubs 
! pfts 13-15 grass
! pfts 16-17 crops

      pft_li(:) = (/ 2,5,10,13,16 /)
      pft_lu(1:nveg-1) = (/ 4,9,12,15 /)
      pft_lu(nveg) = npft
!--amf/van
      do j = js,je 
         do i = is,ie 
            total_iso = 0. 
            total_work_gamma = 0.
            if (has_ts_avg) then
              ts_wrk = ts_avg(i,j,month)
            end if

!-----------------------------------------------------------------
!       ... no emissions for ocean grid point 
!-----------------------------------------------------------------
            if( oro(i-is+1,j-js+1) .eq. 1 ) then
                
!++amv/van
            lai_work(:,:) = mlai(i,j,:,:)
            lai_fac(:) = calc_gamma_lai_age( lai_work, ts_wrk, month, age_work)

            do n = 1, nveg
              nl = pft_li(n)  !beginning index for the pfts that fall within this veg type for n types in ecisop
              nu = pft_lu(n)  !end index
              work_iso(n) = dot_product( lai_fac (nl:nu), pctpft(i,j,nl:nu)) * ecisop(i,j,n)
              work_gamma(n) = dot_product( lai_fac (nl:nu), pctpft(i,j,nl:nu))
            end do
            total_work_gamma = total_work_gamma + sum(work_gamma)
            total_iso = total_iso + sum(work_iso)
        
!--amf/van
              emisop_month(i,j) = total_iso
              diag_gamma_lai_age(i,j) = total_work_gamma
!             if( mpp_pe() == mpp_root_pe() ) then
!                print*, 'AMF: i,j,oro, work_iso, emis ', i,j,oro(i,j), work_iso, emisop_month(i,j)
!             end if
          end if ! land box
        end do !lon
      end do !lat

end subroutine get_monthly_gammas
!</SUBROUTINE


!#######################################################################

! <FUNCTION NAME="calc_gamma_lai_age">
!   <OVERVIEW>
!     Monthly exchange ratio from MEGAN
!   </OVERVIEW>
!   <DESCRIPTION>
!     amf Feb 2009
!     This subroutine calculates gamma_age and gamma_lai according to
!     the MEGAN PCEEA model [Guenther et al., ACP 6, 3181-3210, 2006.]
!     as implemented in mozart4_v4.5 by J.-F. Lamarque and G. Pfister
!     (their fac_lai function)
!   </DESCRIPTION>
!   <TEMPLATE>
!     work_iso  = calc_gamma_lai_age( clai, ts_wrk, month, doage )
!   </TEMPLATE>
!   <IN NAME="clai" TYPE="real" DIM="(:)">
!      12 monthly values for lai for current grid cell
!   </IN>
!   <IN NAME="ts_wrk" TYPE="real">
!     the climatological monthly mean surface T for this i, j, m (K)
!   </IN>
!   <IN NAME="month" TYPE="integer">
!     current month
!   </IN>
!   <IN NAME="doage" TYPE="integer">
!     1 = calculate gamma age; 0 = don't
!   </IN>


function calc_gamma_lai_age(clai, ts_wrk, month, doage )

!AMF - COULD ADD OPTION TO USE AVERAGE TS AND FSDS OVER
! PAST MONTH /WEEKS -- would need to save to restart file  
! NOTE - MZ4 FAC_LAI DOESN"T SEEM TO USE FSDS_AVG, NOR DOES IT IN GUENTHER'S EQNS, 
! SO ELIMINATE FROM THIS FUNCTION

  implicit none

  logical, intent(in) :: doage(:)         !calculate gamma_age?,  amf/van
  integer, intent(in) :: month            !current month index
  real, intent(in) :: ts_wrk              ! currently = climatological sfc T for ea i,j,m
  real, intent(in) :: clai(:,:)           ! monthly lai for this grid cell,   amf/van
  
!-------------------------------------------------------------------------------------
!       ... local variables
!-------------------------------------------------------------------------------------
! ggp: equations 18 and 19 from Guenther et al. include a dependence of ti and tm on temperature
! of preceding timestep (i.e. month here); not considered here. 
! ggp/lamar: instead of using a constant ti and tm, it can be calculated based on information 
! of monthly average temperature
!-------------------------------------------------------------------------------------

!! amf: time step in days btw previous month's LAI and current month's LAI (p3192 guenther)
      integer, parameter :: t(12) = (/ 31,31,28,31,30,31,30,31,31,30,31,30 /)

      integer :: n
      integer :: mnthm1
      real    :: x
      real    :: wrk
      real    :: gamma
      real    :: lai_n, lai_p      ! amf/van
      real    :: Fnew
      real    :: Fgro
      real    :: Fmat
      real    :: Fsen              ! amf/van
      real    :: ti, tm  ! ti = # days btw budbreak and induction of isoprene emissions
                         ! tm = # days btw budbreak + initiation of peak isop emissions rates

!-------------------------------------------------------------------------------------
!       ... function declarations
!-------------------------------------------------------------------------------------
      real    :: calc_gamma_lai_age(npft)    ! amf/van

      if( month > 1 ) then      ! amf/van
         mnthm1 = month - 1
      else
         mnthm1 = 12
      end if

!-------------------------------------------------------------------------------------
!       ... calculations following equations 17&18 in Guenther et al. [2006]
!           -- getting terms needed for gamma age
!-------------------------------------------------------------------------------------
      if( has_ts_avg ) then
         if( ts_wrk <= 303. ) then
            ti = 5. + 0.7*(300. - ts_wrk)                     !eqn 18a (see corrigendum)
         else        
            ti = 2.9                                          ! eqn 18b (corrigendum)
         end if
      else
         ti = 5. + 0.7*(300. - Tdaily_clim)                   ! Tdaily_clim in tropchem_driver_nml
      end if
      tm = 2.3*ti                                             ! Eq 19 

!++amf/van
      calc_gamma_lai_age(1) = 0.
      do n = 2, npft
       if (doage(n)) then
        Fnew = 0.
        Fgro = 0.
        Fmat = 0.
        Fsen = 0.
        lai_n = clai(n, month)
        lai_p = clai(n, mnthm1)
        if( lai_n == lai_p ) then                  !previous LAI = current LAI  - p.392 G06
         Fmat = 0.8
         Fsen = 0.1
         Fgro = 0.1
      else if( lai_p > lai_n ) then              !LAip > LAIc
         Fsen = (lai_p - lai_n) / lai_p
         Fmat = 1. - Fsen
      else if( lai_p < lai_n ) then              !LAIp < LAIc
         Fsen = 0.
         x    = lai_p/lai_n
!--amf/van
         wrk  = 1. - x                                        ! = Fnew
         if( t(month) <= tm ) then
            Fmat =  x                                         ! Eqn 17c
         else
            Fmat = x + (((t(month) - tm)/t(month) ) * wrk)    ! Eqn 17d
         end if
         if( t(month) <= ti ) then
            Fnew = wrk                                        ! Eqn 17a
            Fgro = 1. - (Fnew + Fmat)                         ! Eqn 17e
         else
            Fnew = (ti/t(month)) * wrk                        ! Eqn 17b
            Fgro = 1. - (Fnew + Fmat)                         ! Eqn 17e
         end if
      end if

!-------------------------------------------------------------------------------------
!       ... equations 15 and 16 in Guenther et al. [2006]
!            -- get gamma age
!-------------------------------------------------------------------------------------
      gamma   = .05*Fnew + .6*Fgro + 1.125*Fmat + Fsen   !!  Eq 16
      else
        gamma   = 1.
      end if

!-------------------------------------------------------------------------------------
!     gamma_age ("gamma" below) * gamma_lai where gamma_lai is from Eqn 15 
!-------------------------------------------------------------------------------------
      calc_gamma_lai_age(n) = gamma * .49 * clai(n,month) / sqrt( 1. + 0.2 * clai(n, month)*clai(n, month) )

      end do
  

end function calc_gamma_lai_age

!</FUNCTION>
!############################################################################

! <SUBROUTINE NAME="tropchem_drydep_init">
!   <OVERVIEW>
!     Open dry deposition file
!   </OVERVIEW>
!   <DESCRIPTION>
!     Opens NetCDF file of tracer dry deposition velocities for reading,
!     and set up interpolation to model grid/time
!   </DESCRIPTION>
!   <TEMPLATE>
!     call tropchem_drydep_init( dry_files, dry_names, &
!                                lonb_mod, latb_mod, &
!                                drydep_data )
!   </TEMPLATE>

subroutine tropchem_drydep_init( dry_files, dry_names, &
                                 lonb_mod, latb_mod, &
                                 drydep_data )

!-----------------------------------------------------------------------

   real,                   intent(in),  dimension(:,:) :: lonb_mod, latb_mod
   character(len=64),      intent(out), dimension(:) :: dry_files, dry_names
   type(interpolate_type), intent(out)               :: drydep_data(:)

!-----------------------------------------------------------------------

   integer :: i
   integer :: flag_file, flag_spec
   character(len=64) :: filename,specname
   character(len=64) :: name='', control=''

!-----------------------------------------------------------------------

!---------- Set interpolator type for dry deposition
   call interpolator_init( drydep_data_default, trim(file_dry), lonb_mod, latb_mod, &
                           data_out_of_bounds=(/CONSTANT/), &
                           vert_interp=(/INTERP_WEIGHTED_P/))

   do i = 1,pcnstm1
      dry_files(i) = ''
      dry_names(i) = ''
      if( query_method('dry_deposition',MODEL_ATMOS,indices(i),name,control) )then
         if( trim(name) == 'file' ) then
            flag_file = parse(control, 'file',filename)
            flag_spec = parse(control, 'name',specname)
            if(flag_file > 0 .and. trim(filename) /= trim(file_dry)) then
               dry_files(i) = trim(filename)
               call interpolator_init( drydep_data(indices(i)), trim(filename), lonb_mod, latb_mod,&
                                       data_out_of_bounds=(/CONSTANT/), &
                                       vert_interp=(/INTERP_WEIGHTED_P/))
            else
               dry_files(i) = trim(file_dry)
               drydep_data(indices(i)) = drydep_data_default

            end if
            if(flag_spec >0) then
               dry_names(i) = trim(specname)
            else
               dry_names(i) = trim(lowercase(tracnam(i)))
            end if
         end if
      end if
   end do

end subroutine tropchem_drydep_init
!</SUBROUTINE>

!############################################################################
end module tropchem_driver_mod

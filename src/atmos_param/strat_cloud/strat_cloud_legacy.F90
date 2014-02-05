!FDOC_TAG_GFDL
module strat_cloud_legacy_mod

!-----------------------------------------------------------------------
! NOTE: 10/7/10 -- RSH
! This is equivalent to the version of strat_cloud.F90 which existed
! before mods to include an option for double-moment microphysics (coded
! by m. salzmann) were added. 
!
! This version is equivalent to that which was used in CM3 AR5 experiments.
! Treatment of internal diagnostics and the namelist within the module has
! been changed to be compatible with that used in the revised version of 
! strat_cloud_mod, so that this code may be cleanly accessed from 
! strat_cloud_mod while being in its own module, thus reducing compilation
! time.
!
! Addition of the double moment mods required some code changes which 
! result in model answer changes when the code is optimized (-O2). However,
! when the atmos_phys library is compiled with -O0 and  -fltconsistency,
! the results of an AM3 integration using the revised code (available in 
! strat_cloud.F90) and this legacy version, may be seen to be identical. 
!-----------------------------------------------------------------------
 
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Stephen Klein
! </CONTACT>
!/"/>q
! <OVERVIEW>
!   Code to compute time tendencies of stratiform clouds and
! diagnoses
!   rain and snow flux with prognostic scheme.
!   
! </OVERVIEW>
! <DESCRIPTION>
!
!
!       The prognostic scheme returns the time tendencies of liquid,
!       ice, and saturated volume fraction that are suspended in 
!       stratiform clouds.  The scheme diagnoses the fluxes of rain
!       and snow in saturated and unsaturated areas.
!
!       The prognostic cloud scheme is responsible for determing
!       cloud volume fractions, condensed water tendencies, and
!       the stratiform precipitation rate.  It includes processes
!       for evaporation, condensation, deposition, and sublimation
!       of cloud water, conversion of cloud water to precipitation,
!       evaporation of falling precipitation, the bergeron-findeisan 
!       process, freezing of cloud liquid, accretion of cloud water 
!       by precipitation, and melting of falling precipitation.
!
!       This scheme is based on the experience the author had 
!       at the ECMWF in 1997. The saturated volume fraction formalism 
!       and type of solution follow directly from the scheme of Tiedtke
!       (1993): Monthly Weather Review, Volume 121, pages 3040-3061.
!       The form of most of the microphysics follows Rotstayn , 1997:
!       Quart. J. Roy. Met. Soc. vol 123, pages 1227-1282. The partial
!       precipitation area formulism follows Jakob and Klein, 2000:
!       Quart. J. Roy. Met. Soc. vol 126, pages 2525-2544. 
!
!       The statistical cloud scheme treatment, which is used as
!       a replacement for the Tiedtke cloud fraction scheme, is based
!       on a number of publications: Tompkins, A., 2002: J. Atmos. 
!       Sci., 59, 1917-1942, Klein et al., 2005: J. Geophys. Res., 
!       110, D15S06, doi:10.1029/2004JD005017. 
! </DESCRIPTION>
!

! <DATASET NAME="strat_cloud.res">
!   native format of the restart file
! </DATASET>
! <DATASET NAME="strat_cloud.res.nc">
!   netcdf format of the restart file
! </DATASET>


! <INFO>
!   <REFERENCE>           
!The saturation volume fraction formalism comes from:
!
!Tiedtke, M., 1993: Representation of clouds in large-scale models.
! Mon. Wea. Rev., 121, 3040-3061.
!
! </REFERENCE>
!   <REFERENCE>           
!The form of most of the microphysics follows:
!
!Rotstayn, L., 1997: A physically based scheme for the treatment of
! stratiform clouds and precipitation in large-scale models. I:
! Description and evaluation of microphysical processes. Quart. J.
! Roy. Met. Soc. 123, 1227-1282. 
! </REFERENCE>
!   <COMPILER NAME="">     </COMPILER>
!   <PRECOMP FLAG="">      </PRECOMP>
!   <LOADER FLAG="">       </LOADER>
!   <TESTPROGRAM NAME="">  </TESTPROGRAM>
!   <BUG>                  </BUG>
!   <NOTE> 
!1. qmin should be chosen such that the range of {qmin, max(qa,ql
!,qi)} is resolved by the precision of the numbers used. (default =
! 1.E-10)
!   </NOTE>

!   <NOTE> 
!2. Dmin will be MACHINE DEPENDENT and occur when
!   </NOTE>

!   <NOTE> 
!a. 1. -exp(-Dmin) = 0. instead of Dmin in the limit of very small
! Dmin
!   </NOTE>

!AND

!   <NOTE> 
!b. 1. - exp(-D) < D for all D > Dmin
!   </NOTE>
!   <FUTURE>               </FUTURE>

! </INFO>

  use  sat_vapor_pres_mod, only :  compute_qs
  use beta_dist_mod,        only: beta_dist_init, beta_dist_end, &
                                  incomplete_beta
  use  aer_ccn_act_mod,    only : aer_ccn_act_wpdf, aer_ccn_act_init
  use  aer_in_act_mod,     only : Jhete_dep
  use  mpp_mod,            only : mpp_clock_id, mpp_clock_begin, &
                                  mpp_clock_end, CLOCK_LOOP

  use             fms_mod, only :  file_exist, open_namelist_file,&
       & error_mesg, FATAL, NOTE, mpp_pe, mpp_root_pe, close_file,&
       & read_data, write_data, check_nml_error, write_version_number&
       &, stdlog, open_restart_file, open_ieee32_file, mpp_error
  use  fms_io_mod,         only :  get_restart_io_mode, &
                                   register_restart_field, restart_file_type, &
                                   save_restart, restore_state, get_mosaic_tile_file
  use  constants_mod,      only :  rdgas,rvgas,hlv,hlf,hls, cp_air&
       &,grav,tfreeze,dens_h2o
  use  cloud_rad_mod,      only :  cloud_rad_init
  use  diag_manager_mod,   only :  register_diag_field, send_data
  use  time_manager_mod,   only :  time_type, get_date, get_time
   use cloud_generator_mod, only :  do_cloud_generator, &
                                   cloud_generator_init, &
                                   compute_overlap_weighting
  
  use  rad_utilities_mod,  only : aerosol_type


  use strat_cloud_utilities_mod,  &
                            only: strat_cloud_utilities_init, &
                                  diag_id_type, diag_pt_type, &
                                  strat_nml_type, atmos_state_type, &
                                  particles_type


  implicit none
  private

  integer, private :: sc_loop, sc_pre_loop, sc_post_loop

  public  strat_cloud_legacy, strat_cloud_legacy_init,  &
          strat_cloud_legacy_end
 private  aerosol_effects, ppm2m_sak, kmppm_sak, steepz_sak,  &
          cloud_clear_xfer





  !        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !              GLOBAL STORAGE VARIABLES
  !
  !     radturbten    The sum of radiation and turbulent tendencies
  !                   for each grid box. (K/sec)
  !


  !
  !     ------ data for cloud averaging code ------
  !

  real,    allocatable, dimension (:,:,:) :: qlsum, qisum, cfsum
  integer, allocatable, dimension (:,:)   :: nsum
  !
  !     ------ constants used by the scheme -------
  !


  real, parameter :: d608 = (rvgas-rdgas) / rdgas
  real, parameter :: d622 = rdgas / rvgas
  real, parameter :: d378 = 1. - d622
  !        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !
  !       DECLARE CONSTANTS AND SET DEFAULT VALUES FOR PARAMETERS OF 
  !       THE SCHEME
  !
  !
  !
  !                  PHYSICAL CONSTANTS USED IN THE SCHEME
  !
  !
  !         constant              definition                  unit
  !       ------------   -----------------------------   ---------------
  !
  !           grav       gravitational acceleration      m/(s*s)
  !
  !           hlv        latent heat of vaporization     J/kg condensate
  !
  !           hlf        latent heat of fusion           J/kg condensate
  !
  !           hls        latent heat of sublimation      J/kg condensate
  !
  !           rdgas      gas constant of dry air         J/kg air/K
  !
  !           rvgas      gas constant of water vapor     J/kg air/K
  !
  !           cp_air     specific heat of air at         J/kg air/K
  !                      constant pressure
  !
  !           d622       rdgas divided by rvgas          dimensionless
  !
  !           d378       One minus d622                  dimensionless
  !
  !           tfreeze    Triple point of water           K
  !
  !           dens_h2o   Density of pure liquid          kg/(m*m*m)
  !
  !
  !
  !
  !                          PARAMETERS OF THE SCHEME 
  !
  !
  !         parameter              definition                  unit
  !       ------------   -----------------------------   ---------------
  !
  !
  !       rho_ice        mass density of ice crystals    kg/(m*m*m)
  !
  !       ELI            collection efficiency of        dimensionless
  !                      cloud liquid by falling ice
  !
  !
  !       do_average     Average stratiform cloud properties
  !                      before computing clouds used by radiation?
  !
  !       overlap        value of the overlap parameter
  !                      from cloud rad
  !                      overlap = 1 is maximum-random
  !                      overlap = 2 is random
  !
  
  real,   parameter :: rho_ice        =  100.
  real,   parameter :: ELI            =  0.7
  logical           :: do_average     =  .false.
  integer           :: overlap        =  2

  !
  !-----------------------------------------------------------------------
  !-------------------- diagnostics fields -------------------------------



  !--- for netcdf restart
  type(restart_file_type), pointer, save :: Str_restart => NULL()
  type(restart_file_type), pointer, save :: Til_restart => NULL()
  logical                                :: in_different_file = .false.

 !----------------------------------------------------------

  character(len=5) :: mod_name = 'strat'
  real :: missing_value = -999.


real :: U00, rthresh, var_limit, sea_salt_scale, om_to_oc,  N_land, &
          N_ocean, U_evap, eros_scale, eros_scale_c, eros_scale_t, &
          mc_thresh, diff_thresh, qmin, Dmin, efact, vfact, cfact, &
          iwc_crit,  vfall_const2, vfall_exp2, qthalfwidth, N_min, &
          num_mass_ratio1, num_mass_ratio2
 
   logical :: do_netcdf_restart, u00_profile, use_kk_auto, &
              use_online_aerosol,  use_sub_seasalt, eros_choice, &
              super_choice, tracer_advec, do_old_snowmelt, retain_cm3_bug, do_pdf_clouds, &
              do_liq_num, do_dust_berg, pdf_org, do_ice_nucl_wpdf, debugo
 
   integer :: num_strat_pts,  betaP, nsublevels, kmap, kord, &
              super_ice_opt, isamp, jsamp, ksamp

   character(len=64) :: microphys_scheme

  !       max_strat_pts  maximum number of strat pts
  !                      for instantaneous output
  !
  integer, parameter:: max_strat_pts = 5
  integer,dimension(2,max_strat_pts) :: strat_pts = 0

       
  !        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !
  !       DECLARE VERSION NUMBER OF SCHEME
  !

  Character(len=128) :: Version = '$Id: strat_cloud_legacy.F90,v 20.0 2013/12/13 23:22:13 fms Exp $'
  Character(len=128) :: Tagname = '$Name: tikal $'
   logical            :: module_is_initialized = .false.
  integer, dimension(1) :: restart_versions = (/ 1 /)
  integer               :: vers


  logical  :: cloud_generator_on

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !
  !       The module contains the following subroutines:
  !
  !       strat_cloud_legacy_init  initializes module.
  !
  !       strat_cloud_legacy calculations of the cloud scheme are performed
  !                     here.
  !
  !       strat_cloud_end    closes out module.
  !


CONTAINS





!#######################################################################

subroutine strat_cloud_legacy_init (do_pdf_clouds)

logical,     intent(in) :: do_pdf_clouds

      if (module_is_initialized) return

!-----------------------------------------------------------------------
!    write version number to output file.
!-----------------------------------------------------------------------
      call write_version_number (version, tagname)

!-----------------------------------------------------------------------
!    make sure needed modules have been initialized.
!-----------------------------------------------------------------------
      call strat_cloud_utilities_init
      if (do_pdf_clouds) call beta_dist_init
      call cloud_generator_init
      cloud_generator_on = do_cloud_generator()

!-----------------------------------------------------------------------
!    mark this module initialized.
!-----------------------------------------------------------------------
      module_is_initialized = .true.

end subroutine strat_cloud_legacy_init
 

!#######################################################################


! <SUBROUTINE NAME="strat_cloud_legacy">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!       
!  </DESCRIPTION>
!  <TEMPLATE>
!   call strat_cloud_legacy( 
!              cloud_generator_on,Nml,  &
!              diag_id, diag_pt, n_diag_4d, n_diag_4d_kp1,  &
!              Time,is,ie,js,je,dtcloud,pfull,phalf,radturbten2,&
!              T,qv,ql,qi,qa,omega,Mc,diff_t,LAND,              &
!              ST,SQ,SL,SI,SA,rain3d,snow3d,snowclr3d,surfrain,     &
!              surfsnow,qrat,ahuco,limit_conv_cloud_frac,MASK, &
!              qn, Aerosol, SN)
!
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!         Time
!  </IN>
!  <IN NAME="is" TYPE="integer">
!         Indice of starting point in the longitude direction of the slab being passed to strat_cloud
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!         Indice of ending point in the longitude direction of the slab being passed 
!  </IN>
!  <IN NAME="js" TYPE="integer">
!         Indice of starting point in the latitude direction of the slab being passed
!  </IN>
!  <IN NAME="je" TYPE="integer">
!         Indice of ending point in the latitude direction of the slab being passed 
!  </IN>
!  <IN NAME="dtcloud" TYPE="real">
!         Physics time step (sec)
!  </IN>
!  <IN NAME="pfull" TYPE="real">
!         Pressure on model full levels (Pa)
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!         Pressure on model half levels (Pa)
!  </IN>
!  <IN NAME="radturbten2" TYPE="real">
!         Sum of the tendencies of temperature from turbulence and radiation schemes (K/s)
!  </IN>
!  <IN NAME="T" TYPE="real">
!         Temperature (K)         
!  </IN>
!  <IN NAME="qv" TYPE="real">
!         Water vapor specific humidity (kg vapor/kg air)
!  </IN>
!  <IN NAME="ql" TYPE="real">
!         Grid-box mean liquid water specific humidity (kg liquid/kg air)
!  </IN>
!  <IN NAME="qi" TYPE="real">
!         Grid-box mean ice water specific humidity (kg ice/kg air)
!  </IN>
!  <IN NAME="qa" TYPE="real">
!         Cloud fraction (3d array and a prognostic variable) (fraction)
!  </IN>
!  <IN NAME="qn" TYPE="real">
!         Cloud droplet number (3d array and a prognostic variable) (#/kg air)
!  </IN>
!  <IN NAME="omega" TYPE="real">
!         Vertical pressure velocity (Pa/sec)
!  </IN>
!  <IN NAME="Mc" TYPE="real">
!         Cumulus mass flux (defined positive as upward) (kg air/m2/sec)
!  </IN>
!  <IN NAME="diff_t" TYPE="real">
!         Vertical diffusion coefficient for temperature and tracer from vertical diffusion scheme (m2/sec) 
!  </IN>
!  <IN NAME="LAND" TYPE="real">
!         Fraction of surface that contains land (fraction)
!  </IN>
!  <OUT NAME="ST" TYPE="real">
!         Change in temperature due to strat_cloud (K) 
!  </OUT>
!  <OUT NAME="SQ" TYPE="real">
!         Change in water vapor due to strat_cloud (kg vapor/kg air) 
!  </OUT>
!  <OUT NAME="SL" TYPE="real">
!         Change in cloud liquid due to strat_cloud (kg liquid/kg air)
!  </OUT>
!  <OUT NAME="SI" TYPE="real">
!         Change in cloud ice due to strat_cloud (kg ice/kg air)
!  </OUT>
!  <OUT NAME="SA" TYPE="real">
!         Change in cloud fraction due to strat_cloud (fraction)
!  </OUT>
!  <OUT NAME="SN" TYPE="real">
!         Change in cloud droplet number due to strat_cloud (fraction)
!  </OUT>
!  <OUT NAME="surfrain" TYPE="real">
!         Surface rain fall over time step dtcloud (kg liquid/m2)
!  </OUT>
!  <OUT NAME="surfsnow" TYPE="real">
!         Surface snow fall over time step dtcloud (kg ice/m2)
!  </OUT>
!  <OUT NAME="rain3d" TYPE="real">
!         3D rain fall over time step dtcloud (kg liquid/m2)
!  </OUT>
!  <OUT NAME="snow3d" TYPE="real">
!         3D snow fall over time step dtcloud (kg ice/m2)
!  </OUT>
!  <IN NAME="qrat" TYPE="real">
!         Ratio of large-scale specific humidity to specific humidity in 
!         environment outside convective system (from donner_deep) 
!         
!         Will be equal to 1 for all normal AM2 operations (i.e. donner_deep is not activated)              
!         
!         Note that index 1 is nearest ground
!
!  </IN>
!  <IN NAME="ahuco" TYPE="real">
!         The fraction of the grid box containing either cumulus cells 
!         or the mesoscale circulation (from donner_deep).
!
!         Will be equal to 0 for all normal AM2 operations 
!         (i.e. donner_deep is not activated)              
!         
!         Note that index 1 is nearest ground
!
!  </IN>
!  <IN NAME="MASK" TYPE="real">
!         Optional input real array indicating the point is above the surface
!         if equal to 1.0 and indicating the point is below the surface if 
!         equal to 0.
!
!         Used only in eta vertical coordinate model.
!  </IN>
! </SUBROUTINE>
!


subroutine strat_cloud_legacy( Nml,  &
       diag_id, diag_pt, n_diag_4d, n_diag_4d_kp1,  &
        diag_4d, diag_4d_kp1, diag_3d, &
      Time,is,ie,js,je,dtcloud,pfull,phalf,radturbten2,&
    T,qv,ql,qi,qa,omega,Mc,diff_t,LAND,              &
    ST,SQ,SL,SI,SA,f_snow_berg,rain3d,snow3d,snowclr3d,surfrain,     &
    surfsnow,qrat,ahuco,limit_conv_cloud_frac,MASK, &
    qn, Aerosol, SN)


  !        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !
  !
  !       VARIABLES
  !
  !
  !
  !       ------
  !       INPUT:
  !       ------
  !
  !         variable              definition                  unit
  !       ------------   -----------------------------   ---------------
  !
  !       Time           time type variable 
  !
  !       is,ie          starting and ending i indices 
  !                      for data window
  !
  !       js,je          starting and ending j indices 
  !                      for data window
  !
  !       dtcloud        time between this call and      s
  !                      the next call to strat_cloud
  !
  !       pfull          pressure at full model levels   Pa
  !                      IMPORTANT NOTE: p(j)<p(j+1)
  !
  !       phalf          pressure at half model levels   Pa
  !                      phalf(j)<pfull(j)<phalf(j+1)
  !
  !       T              temperature                     K
  !
  !       qv             specific humidity of water      kg vapor/kg air
  !                      vapor
  !
  !       ql             specific humidity of cloud      kg condensate/
  !                      liquid                          kg air
  !
  !       qi             specific humidity of cloud      kg condensate/
  !                      ice                             kg air
  !
  !       qa             saturated volume fraction       fraction
  !
  !       qn             cloud droplet number            #/kg air
  !
  !       qrat           ratio of large-scale spec       fraction
  !                      humidity to spec humidity
  !                      in environment outside
  !                      convective system (from
  !                      donner_deep) 
  !                      index 1 nearest ground
  !
  !       ahuco          fraction, cell+meso, from       fraction
  !                      donner_deep
  !                      index 1 nearest ground
  !
  !       omega          vertical pressure velocity      Pa/s
  !
  !       Mc             Cumulus mass flux defined       kg/(m*m)/s
  !                      on full levels
  !
  !       diff_t         Vertical diffusion coefficient  (m*m)/s
  !                      for temperature and tracer
  !
  !       LAND           the fraction of the grid box    fraction
  !                      covered by land
  !                               
  !
  !       -------
  !       OUTPUT:
  !       -------
  !
  !         variable              definition                  unit
  !       ------------   -----------------------------   ---------------
  !
  !       ST             temperature change due to       K
  !                      all stratiform processes
  !
  !       SQ             water vapor change due to       kg vapor/kg air
  !                      all stratiform processes
  !
  !       SL             cloud liquid change due to      kg condensate/
  !                      all stratiform processes        kg air
  !
  !       SI             cloud ice change due to         kg condensate/
  !                      all stratiform processes        kg air
  !
  !       SA             saturated volume fraction       fraction
  !                      change due to all stratiform 
  !                      processes
  !
  !       SN             cloud droplet number            #/kg air
  !                      change due to all stratiform 
  !                      processes
  !
  !       surfrain       rain that falls through the     kg condensate/
  !                      bottom of the column over       (m*m)
  !                      the time dtcloud
  !
  !       surfsnow       snow that falls through the     kg condensate/
  !                      bottom of the column over       (m*m)
  !                      the time dtcloud
  !
  !       rain3d         rain that falls through the     kg condensate/
  !                      each of the model layer         (m*m)/sec
  !
  !       snow3d         snow that falls through the     kg condensate/
  !                      each of the model layer         (m*m)/sec
  !
  !
  !       ---------------
  !       optional INPUT:
  !       ---------------
  !
  !         variable              definition                  unit
  !       ------------   -----------------------------   ---------------
  !
  !
  !       MASK           real array indicating the 
  !                      point is above the surface
  !                      if equal to 1.0 and 
  !                      indicating the point is below
  !                      the surface if equal to 0.
  !
  !       -------------------
  !       INTERNAL VARIABLES:
  !       -------------------
  !
  !         variable              definition                  unit
  !       ------------   -----------------------------   ---------------
  !
  !       kdim           number of vertical levels
  !
  !       j              model vertical level being 
  !                      processed
  !
  !       ipt,jpt        i and j point indice used only
  !                      in instantaneous diag-
  !                      nostic output
  !
  !       i,unit,nn      temporary integers used in
  !                      instantaneous diagnostic
  !                      output.
  !
  !       inv_dtcloud    1 / dtcloud                     1/sec
  !
  !       airdens        air density                     kg air/(m*m*m)
  !
  !       qs             saturation specific humidity    kg vapor/kg air
  !
  !       dqsdT          T derivative of qs              kg vapor/kg air/K
  !
  !       gamma          (L/cp)*dqsdT                    dimensionless
  !
  !       rain_clr       grid mean flux of rain enter-   kg condensate/
  !                      ing the grid box from above     (m*m)/s
  !                      and entering the unsaturated 
!                      portion of the grid box
!
!       rain_cld       grid mean flux of rain enter-   kg condensate/
!                      ing the grid box from above     (m*m)/s
!                      and entering the saturated 
!                      portion of the grid box
!
!       a_rain_clr     fraction of grid box occupied   fraction
!                      by rain_clr
!
!       a_rain_cld     fraction of grid box occupied   fraction
!                      by rain_cld
!
!       snow_cld       flux of ice entering the        kg condensate/
!                      saturated portion of the        (m*m)/s
!                      grid box from above by means
!                      of gravitational settling 
!
!       snow_clr       flux of ice outside of cloud    kg condensate/
!                      entering the unsaturated        (m*m)/s
!                      portion of the grid box from      
!                      above
!
!       a_snow_clr     area fraction of grid box       fraction
!                      covered by snow flux in
!                      unsaturated air
!
!       a_snow_cld     area fraction of grid box       fraction
!                      covered by the snow flux in
!                      saturated air
!
!       deltpg         pressure thickness of grid box  kg air/(m*m)
!                      divided by gravity
!
!       U              grid box relative humidity      fraction
!
!       U00p           critical relative humidity      fraction
!                      which is a function of pressure 
!
!       dqs_ls         change in saturation specific   kg vapor/kg air
!                      due to large-scale processes,
!                      such as large-scale vertical
!                      motion, compensating convective
!                      mass flux, or radiative cooling
!
!       da_ls          change in saturated volume      fraction
!                      fraction due to large-scale
!                      processes
!
!       C_dt           product of A and dtcloud in     dimensionless in 
!                      in the analytic integration     qa integration
!                      of the qa equation, or C and
!                      dtcloud in the analytic         kg condensate/
!                      integration of the ql and qi    kg air in ql or 
!                      equations.                      qi integration
!
!       D_dt           product of B and dtcloud in     dimensionless in
!                      in the analytic integration     qa, ql, and qi
!                      of the qa equation, or D and    integration
!                      dtcloud in the analytic         
!                      integration of the ql and qi    
!                      equations.                      
!
!       qceq           equilibrium value of cloud      dimensionless or
!                      fraction or cloud condensate    kg condensate /
!                      that the analytic integration   kg air
!                      approaches                      
!
!       qcbar          mean value of cloud fraction    dimensionless or
!                      or cloud condensate over the    kg condensate /
!                      t0 to t0 + dtcloud interval     kg air
!
!       qc0            value of cloud fraction or      dimensionless or
!                      cloud condensate at the         kg condensate /
!                      initial time                    kg air
!        
!       qc1            value of cloud fraction or      dimensionless or
!                      cloud condensate at the final   kg condensate /
!                      time                            kg air
!       
!       D1_dt          first sink in either ql or qi   dimensionless
!                      equation. This is analogous to
!                      D_dt.  In ql equation, this 
!                      sink represents the conversion
!                      of cloud liquid to rain. In the
!                      qi equation it represents the
!                      settling of ice crystals.
!
!       D2_dt          second sink in ql or qi         dimensionless
!                      equation. This is analogous 
!                      to D_dt. In ql equation this
!                      sink represents the conversion
!                      of cloud liquid to ice. In the
!                      qi equation this sink 
!                      represents the melting of 
!                      cloud ice into rain.
!
!       D_eros         Sink in ql, qi and qa equation  dimensionless
!                      due to turbulent erosion of
!                      cloud sides
! 
!       ql_upd         updated value of ql             kg condensate/
!                                                      kg air
!       
!       qi_upd         updated value of qi             kg condensate/
!                                                      kg air
!
!       qa_upd         updated value of qa             fraction
!
!       qa_mean        qa + SA; semi-implicit          fraction
!                      saturated volume fraction
!
!       qa_mean_lst    qa_mean of the level above      fraction
!
!       ql_mean        ql + positive increment         kg condensate/
!                      of ql; i.e. a sort of           kg air
!                      intermediate ql
!
!       qi_mean        ql + positive increment         kg condensate/
!                      of qi; i.e. a sort of           kg air
!                      intermediate qi
!
!       dcond_ls       change in condensate due to     kg condensate/
!                      non-convective condensation.    kg air
!                      After phase determination,
!                      this variable refers only to
!                      liquid condensation.
!
!       dcond_ls_ice   change in ice due to            kg condensate/
!                      non-convective condensation.    kg air
!
!       da_cld2clr     fraction of the area in which   fraction
!                      rain/snow in saturated volume 
!                      above falls into unsaturated 
!                      volume in the current layer.
!
!       da_clr2cld     as in da_cld2clr except for     fraction
!                      the transfer from unsaturated
!                      to saturated volume
!
!       dprec_cld2clr  grid mean flux that is trans-   kg condensate/
!                      ferred from rain/snow in        (m*m)/s
!                      saturated volume to rain/snow 
!                      in unsaturated volume at layer 
!                      interfaces.
!
!       dprec_clr2cld  as in dprec_cld2clr except for  kg condensate/
!                      the transfer from unsaturated   (m*m)/s
!                      to saturated volume.
!
!       N              fixed number of cloud drops     1/(m*m*m)
!                      per unit volume in liquid
!                      clouds
!     
!       rad_liq        mean volume radius of liquid    microns
!                      cloud drops
!
!       A_plus_B       sum of vapor diffusion factor   m*s/kg
!                      and thermal conductivity factor
!                      which is used in various 
!                      microphysical formula for the 
!                      evaporation of rain and snow
!
!       Vfall          fall speed of ice crystals      m/s
!
!       lamda_f        slope factor in the SIZE        1/m
!                      distribution of ice crystals
!
!       U_clr          relative humidity in the clear  fraction
!                      portion of the grid box.
!
!       tmp1,tmp2,tmp3 temporary numbers used at 
!                      several points within the
!                      subroutine
!
!
!                STATISTICAL CLOUD SCHEME VARIABLES
!
!       qag            equilibrium value of cloud      dimensionless 
!                      fraction for statistical 
!                      cloud scheme           
!
!       qcg            equilibrium value of cloud      kg condensate /
!                      condensate that PDF clouds      kg air
!                      wants                          
!
!       qcg_ice        equilibrium value of cloud      kg condensate /
!                      ice condensate that PDF clouds  kg air
!                      wants
!
!       qvg            equilibrium value of water      kg vapor /
!                      vapor in the clear portion      kg air
!                      of the grid box that PDF 
!                      clouds wants                          
!
!       qtbar          total water specific humidity   kg water /
!                      which is equal to the sum of    kg air
!                      liquid water, ice water, and
!                      water vapor
!
!       deltaQ         the width of the total water    kg water /
!                      subgrid distribution (= qtmax   kg air
!                      minus qtmin)
!
!       qtmin          the minimum value to the total  kg water /
!                      sub-grid scale distribution     kg air
!
!       qs_norm        the difference between the      dimensionless
!                      saturation specific humidity    
!                      and qtmin normalized by deltaQ
!
!       icbp           the value of the incomplete     dimensionless
!                      beta function evaluated with
!                      x=qs_norm, p=betaP, and q=betaP
!
!       icbp1          the value of the incomplete     dimensionless
!                      beta function evaluated with                  
!                      x=qs_norm, p=betaP+1, and 
!                      q=betaP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!


!        
!       user Interface variables
!       ------------------------
!

        type(strat_nml_type), intent(in)  :: Nml
       TYPE(diag_id_type), intent(inout) :: diag_id
       TYPE(diag_pt_type), intent(inout) :: diag_pt
       integer, intent(in) :: n_diag_4d, n_diag_4d_kp1
     real, dimension(:,:,:,0:), intent(inout) :: diag_4d, diag_4d_kp1
     real, dimension(:,:,0:),   intent(inout) :: diag_3d
        type(time_type), intent (in)           :: Time
        integer, intent (in)                   :: is,ie,js,je
        real, intent (in)                      :: dtcloud
        real, intent (in),    dimension(:,:,:) :: pfull,phalf
        real, intent (in),    dimension(:,:,:) :: T,qv,ql,qi,qa,omega
        real, intent (in),    dimension(:,:,:) :: Mc, diff_t
        real, intent (in),    dimension(:,:,:) :: qrat,ahuco
        logical, intent(in)                    :: limit_conv_cloud_frac
        real, intent (in),    dimension(:,:)   :: LAND
        real, intent (in),    dimension(:,:,:) :: radturbten2
        real, intent (out),   dimension(:,:,:) :: ST,SQ,SL,SI,SA
        real, intent (out),   dimension(:,:)   :: surfrain,surfsnow
        real, intent (in), optional, dimension(:,:,:) :: MASK
        real, intent (in),  optional, dimension(:,:,:) :: qn
        real, intent (out),   dimension(:,:,:) :: f_snow_berg
        real, intent (out),   dimension(:,:,:) :: rain3d
        real, intent (out),   dimension(:,:,:) :: snow3d
        real, intent (out),   dimension(:,:,:) :: snowclr3d
        type(aerosol_type), intent (in), optional      :: Aerosol  
        real, intent (out), optional, dimension(:,:,:) :: SN

!
!       Internal variables
!       ------------------
!

        integer                                        :: idim,jdim,kdim
        integer                                        :: id,jd,ns
        integer                                        :: j,ipt,jpt
        integer                                        :: i,unit,nn        
        real                                           :: inv_dtcloud, Si0
        real                                           :: icbp, icbp1, pnorm
        real, dimension(size(T,1),size(T,2),size(T,3)) :: airdens
        real, dimension(size(T,1),size(T,2),size(T,3)) :: qs,dqsdT
        real, dimension(size(T,1),size(T,2),size(T,3)) :: gamma
        real, dimension(size(T,1),size(T,2),size(T,3)) :: A_plus_B
        real, dimension(4,size(T,1),size(T,2),size(T,3)) :: qta4,qtqsa4
        real, dimension(size(T,1),size(T,2),size(T,3)) :: delp
        real, dimension(size(T,1),size(T,2))   :: rain_clr,rain_cld
        real, dimension(size(T,1),size(T,2))   :: a_rain_clr,a_rain_cld
        real, dimension(size(T,1),size(T,2))   :: snow_clr,snow_cld
        real, dimension(size(T,1),size(T,2))   :: a_snow_clr,a_snow_cld
        real, dimension(size(T,1),size(T,2))   :: deltpg,U,U00p
        real, dimension(size(T,1),size(T,2))   :: dqs_ls,da_ls
!rab        real, dimension(size(T,1),size(T,2))   :: C_dt, D_dt
        real, dimension(size(T,1),size(T,2))   :: D1_dt,D2_dt,D_eros
        real, dimension(size(T,1),size(T,2))   :: qcg_ice, qvg
!rab        real, dimension(size(T,1),size(T,2))   :: qceq, qcbar, qcg, qag
        real, dimension(size(T,1),size(T,2))   :: qcg, qag
!rab        real, dimension(size(T,1),size(T,2))   :: qagtmp,qcgtmp,qvgtmp
!rab        real, dimension(size(T,1),size(T,2))   :: qc1, qc0
        real, dimension(size(T,1),size(T,2))   :: ql_upd,qi_upd,qa_upd, qn_upd
        real, dimension(size(T,1),size(T,2))   :: ql_mean,qi_mean, qn_mean
        real, dimension(size(T,1),size(T,2))   :: qa_mean,qa_mean_lst
        real, dimension(size(T,1),size(T,2))   :: dcond_ls,dcond_ls_ice
        real, dimension(size(T,1),size(T,2))   :: N,rad_liq, N3D_col
!       real, dimension(size(T,1),size(T,2),size(T,3)) :: N3D, &
        real, dimension(size(T,1),size(T,2),size(T,3)) ::      &
                                                        concen_dust_sub
        real, dimension(size(T,1),size(T,2))   :: Vfall,iwc,lamda_f
        real, dimension(size(T,1),size(T,2))   :: U_clr
        real, dimension(size(T,1),size(T,2))   :: tmp1,tmp2,tmp3,tmp5, delta_cf,drop1,crystal
        real, dimension(size(T,1),size(T,2))   :: sum_freeze, sum_rime, &
                                                  sum_berg
        real, dimension(size(T,1),size(T,2))   :: qtbar,deltaQ
        real, dimension(size(T,1),size(T,2))   :: qtmin,qs_norm          

        logical :: used
        integer            :: k
        real               :: thickness, up_strat, wp2 ! cjg
        real, dimension(size(T,1),size(T,2),size(T,3),4) :: totalmass1
!rab - variables necessary to clean up strat_cloud....
        real               :: freeze_pt, tmp1s, tmp2s, tmp3s, snow_fact
        real               :: qc0s, qc1s, qceqs, qcbars, C_dts, D_dts
        real               :: qagtmps,qcgtmps,qvgtmps
        real               :: qldt_sum

!       real, allocatable, dimension(:,:,:,:) :: diag_4d, diag_4d_kp1              
!       real, allocatable, dimension(:,:,:) :: diag_3d        
        real, dimension(size(T,1),size(T,2),size(T,3)) :: deltpg_3d

        integer :: n8
     call mpp_clock_begin(sc_pre_loop)

    do_netcdf_restart = Nml%do_netcdf_restart
    U00 = Nml%U00 
    U00_profile = Nml%U00_profile 
    rthresh = Nml%rthresh 
    use_kk_auto = Nml%use_kk_auto 
    var_limit = Nml%var_limit 
    use_online_aerosol = Nml%use_online_aerosol 
    sea_salt_scale = Nml%sea_salt_scale 
    om_to_oc = Nml%om_to_oc 
    N_land = Nml%N_land 
    use_sub_seasalt = Nml%use_sub_seasalt 
    N_ocean = Nml%N_ocean 
    U_evap = Nml%U_evap 
    eros_scale = Nml%eros_scale 
    eros_choice = Nml%eros_choice
    eros_scale_t = Nml%eros_scale_t 
    eros_scale_c = Nml%eros_scale_c 
    mc_thresh = Nml%mc_thresh 
    diff_thresh = Nml%diff_thresh 
    super_choice = Nml%super_choice 
    tracer_advec = Nml%tracer_advec 
    qmin = Nml%qmin 
    Dmin = Nml%Dmin 
    num_strat_pts = Nml%num_strat_pts 
    if (num_strat_pts > 0) then
      strat_pts(:,1:) = Nml%strat_pts(:,1:) 
    endif
    efact = Nml%efact 
    vfact = Nml%vfact 
    cfact = Nml%cfact 
    do_old_snowmelt = Nml%do_old_snowmelt 
    retain_cm3_bug  = Nml%retain_cm3_bug  
    do_pdf_clouds = Nml%do_pdf_clouds 
    betaP = Nml%betaP 
    iwc_crit = Nml%iwc_crit 
    vfall_const2 = Nml%vfall_const2
    vfall_exp2 = Nml%vfall_exp2 
    qthalfwidth = Nml%qthalfwidth 
    nsublevels = Nml%nsublevels 
    kmap = Nml%kmap 
    kord = Nml%kord
    do_liq_num = Nml%do_liq_num 
    do_dust_berg = Nml%do_dust_berg 
    N_min = Nml%N_min 
    num_mass_ratio1 = Nml%num_mass_ratio1 
    num_mass_ratio2 = Nml%num_mass_ratio2 
    microphys_scheme = Nml%microphys_scheme
    super_ice_opt = Nml%super_ice_opt 
    pdf_org = Nml%pdf_org 
    do_ice_nucl_wpdf = Nml%do_ice_nucl_wpdf 
    debugo = Nml%debugo 
    isamp = Nml%isamp 
    jsamp = Nml%jsamp 
    ksamp = Nml%ksamp

!-----------------------------------------------------------------------
!       
!
     if (.not.module_is_initialized) then
       call error_mesg('strat_cloud_legacy','strat_cloud_legacy is not initialized',FATAL)
     endif





     if (diag_id%ql_wt > 0) then
       diag_4d(:,:,:,diag_pt%ql_wt) = ql
     endif
!-----------------------------------------------------------------------
!
!       initialize select variables to zero. The variables reset
!       are:
!
!       (1) changes of prognostic variables
!
!       (2) variables dealing with the rain/snow fluxes. 
!
!       (3) qa_mean of the level above the top level.
!        
!       (4) diagnostic output fields

        ST = 0.
        SQ = 0.
        SL = 0.
        SI = 0.
        SA = 0.
        if (present(SN)) SN = 0.
   
        rain_cld   = 0.
        rain_clr   = 0.
        a_rain_cld = 0.
        a_rain_clr = 0.
        snow_cld   = 0.
        snow_clr   = 0.
        a_snow_clr = 0.
        a_snow_cld = 0.
        
        qa_mean_lst= 0.

        dcond_ls      = 0.
        dcond_ls_ice  = 0.
        qcg           = 0.
        qcg_ice       = 0.
        
        
                     
!-----------------------------------------------------------------------
!
!       Determine dimensions of slab

        idim = SIZE(T,1)
        jdim = SIZE(T,2)
        kdim = SIZE(T,3)
        
!-----------------------------------------------------------------------
!
!       compute inverse time step

        inv_dtcloud = 1.0 / dtcloud

!-----------------------------------------------------------------------
!
!       Calculate saturation specific humidity and its temperature 
!       derivative, and thermal conductivity plus vapor diffusivity
!       factor.
!
!       These are calculated according to the formulas:
!
!   (1)  qs   = d622*esat/ [pfull  -  (1.-d622)*esat]
!
!   (2) dqsdT = d622*pfull*(desat/dT)/[pfull-(1.-d622)*esat]**2.
!
!   (3) gamma = (L/cp) * dqsdT
!       
!       where d622 = rdgas/rvgas; esat = saturation vapor pressure;
!       and desat/dT is the temperature derivative of esat.
!       Note that in the calculation of gamma, 
!
!            {             hlv          for T > tfreeze             }
!       L =  { 0.05*(T-tfreeze+20.)*hlv + 0.05*(tfreeze-T)*hls      }
!            {                          for tfreeze-20.< T < tfreeze}
!            {             hls          for T < tfreeze-20.         }
!
!       This linear form is chosen because at tfreeze-20. es = esi, and
!       at tfreeze, es = esl, with linear interpolation in between.
!
!       The conductivity/diffusivity factor, A_plus_B is given by:
!
!   (4) A_plus_B =   { (hlv/Ka/T)*((hlv/rvgas/T)-1.) } + 
!
!                    { (rvgas*T/chi*esat) }
!
!       where Ka is the thermal conductivity of air = 0.024 J/m/s/K
!       and chi is the diffusitivy of water vapor in air which is
!       given by
!
!   (5) chi = 2.21 E-05 (m*m)/s  * (1.E+05)/pfull
!
!       where p is the pressure in Pascals.
!    
!
!       Note that qs, dqsdT, and gamma do not have their proper values
!       until all of the following code has been executed.  That
!       is qs and dqsdT are used to store intermediary results
!       in forming the full solution.

        !calculate water saturated vapor pressure from table
        !and store temporarily in the variable gamma
        !calculate qs and dqsdT
        if (do_pdf_clouds) then 
             call compute_qs( T-((hlv*ql+hls*qi)/cp_air), pfull, qs, &
                             dqsdT=dqsdT, esat=gamma )
        else
             call compute_qs( T, pfull, qs, dqsdT=dqsdT, esat=gamma )
        end if
                     
        !compute A_plus_B
        A_plus_B = ( (hlv/0.024/T) * ((hlv/rvgas/T)-1.) ) +            &
           (rvgas*T*pfull/2.21/gamma)  
         
        !calculate gamma
        if (do_pdf_clouds) then
        gamma = dqsdT *(min(1.,max(0.,0.05*(T-((hlv*ql+hls*qi)/cp_air) &
                                             -tfreeze+20.)))*hlv +     &
                        min(1.,max(0.,0.05*(tfreeze -T+((hlv*ql+hls*qi)&
                                              /cp_air)   )))*hls)/cp_air
        else
        gamma = dqsdT *(min(1.,max(0.,0.05*(T-tfreeze+20.)))*hlv +     &
                        min(1.,max(0.,0.05*(tfreeze -T   )))*hls)/cp_air
        end if             
!-----------------------------------------------------------------------
!
!       Calculate air density

!       airdens = pfull / (rdgas * T * (1. + d608*qv  - ql - qi) )
        airdens = pfull / (rdgas * T * (1.   - ql - qi) )
        where (qrat .gt. 0.) 
             airdens = pfull / (rdgas * T *(1.+(d608*qv/qrat)-ql-qi) )
        end where

!-----------------------------------------------------------------------
!
!       Assign cloud droplet number based on land or ocean point.
        
        N = N_land*LAND + N_ocean*(1.-LAND)

!---------------------------------------------------------------------
!   call aerosol_effects to include impact of aerosols on the cloud
!   droplet number and the bergeron process, if these effects activated.
!---------------------------------------------------------------------
        if (do_liq_num .or. do_dust_berg) then
       call aerosol_effects (is, js, Time, phalf, airdens, T, &
                            diag_4d, diag_id, diag_pt,          &
                            concen_dust_sub, totalmass1, Aerosol, mask)
        endif


!-----------------------------------------------------------------------
!
!       Is a sub-vertical grid scale distribution going to be neededed?  
!       If yes, then do ppm fits
!
        if (do_pdf_clouds) then
        
        !initialize quantities
        do j = 1, kdim
             delp(:,:,j) = phalf(:,:,j+1)-phalf(:,:,j)
        enddo     
        qta4(1,:,:,:) = max(qmin,qv+ql+qi)
        qtqsa4(1,:,:,:) = qta4(1,:,:,:)-qs
        
        if (nsublevels.gt.1) then
            do id = 1,idim
                call ppm2m_sak(qta4(:,id,:,:),delp(id,:,:),kdim,kmap,1,jdim,0,kord)
                call ppm2m_sak(qtqsa4(:,id,:,:),delp(id,:,:),kdim,kmap,1,jdim,0,kord)
            enddo                
        else
            qta4(2,:,:,:) = qta4(1,:,:,:)
            qta4(3,:,:,:) = qta4(1,:,:,:)
            qta4(4,:,:,:) = 0.
            qtqsa4(2,:,:,:) = qtqsa4(1,:,:,:)
            qtqsa4(3,:,:,:) = qtqsa4(1,:,:,:)
            qtqsa4(4,:,:,:) = 0.   
        end if

        end if  !end of do_pdf_clouds section

!rab - statements moved outside of large loop
!      they have no bearing on the overall cloud physics
!      they are merely diagnostic values
!
     if (.not. do_pdf_clouds) then

      if (max(diag_id%qadt_fill,diag_id%qa_fill_col) > 0) then
        where (qa .le. qmin) 
!         qadt_fill = -qa * inv_dtcloud
          diag_4d(:,:,:,diag_pt%qadt_fill) = -qa * inv_dtcloud
        endwhere
      endif

      if (max(diag_id%qidt_fill,diag_id%qi_fill_col) > 0) then
        where (qi.le.qmin .or. qa.le.qmin) 
!         qidt_fill = -qi * inv_dtcloud
          diag_4d(:,:,:,diag_pt%qidt_fill) = -qi * inv_dtcloud
        endwhere
      endif

        if (.not. do_liq_num) then
!a        N3D = 0.
          diag_4d(:,:,:,diag_pt%droplets) = 0.
          if (max(diag_id%qldt_fill,diag_id%ql_fill_col) > 0) then
            where (ql .le. qmin .or. qa .le. qmin) 
!             qldt_fill = -ql * inv_dtcloud
              diag_4d(:,:,:,diag_pt%qldt_fill) = -ql * inv_dtcloud
            endwhere
          endif
        else
!a        N3D = qn*airdens*1.e-6
          diag_4d(:,:,:,diag_pt%droplets) = qn(:,:,:)*airdens(:,:,:)*1.e-6
!         if (diag_id%debug1_3d > 0) debug1 = min(qa,1.)
          if (diag_id%cf_liq_init   > 0) diag_4d(:,:,:,diag_pt%cf_liq_init  ) = min(qa,1.)
          do j=1,kdim
           do k=1,jdim
            do i=1,idim
             if (ql(i,k,j).le.qmin .or. qa(i,k,j).le.qmin .or. qn(i,k,j).le.qmin) then
!a            N3D(i,k,j) = 0.
              diag_4d(i,k,j,diag_pt%droplets) = 0.
!             if (max(diag_id%qldt_fill,diag_id%ql_fill_col) > 0) qldt_fill(i,k,j) = -ql(i,k,j) * inv_dtcloud
              if (max(diag_id%qldt_fill,diag_id%ql_fill_col) > 0) diag_4d(i,k,j,diag_pt%qldt_fill) = -ql(i,k,j) * inv_dtcloud
!             if (max(diag_id%qndt_fill,diag_id%qn_fill_col) > 0) qndt_fill(i,k,j) = -qn(i,k,j) * inv_dtcloud
              if (max(diag_id%qndt_fill,diag_id%qn_fill_col, &
                           diag_id%qldt_fill,diag_id%ql_fill_col) > 0)  &
               diag_4d(i,k,j,diag_pt%qndt_fill) = -qn(i,k,j) * inv_dtcloud
!             if (diag_id%debug1_3d > 0) debug1(i,k,j) = 0.
              if (diag_id%cf_liq_init   > 0) diag_4d(i,k,j,diag_pt%cf_liq_init  ) = 0.
             endif
            enddo
           enddo
          enddo
        endif

     else

      if (max(diag_id%qidt_fill,diag_id%qi_fill_col) > 0) then 
        where (qi .le. qmin) 
!         qidt_fill = -qi * inv_dtcloud
          diag_4d(:,:,:,diag_pt%qidt_fill) = -qi * inv_dtcloud
        endwhere
      endif

        if (.not. do_liq_num) then
!a        N3D = 0.
          diag_4d(:,:,:,diag_pt%droplets) = 0.
          if (max(diag_id%qldt_fill,diag_id%ql_fill_col) > 0) then
            where (ql .le. qmin) 
!             qldt_fill = -ql * inv_dtcloud
              diag_4d(:,:,:,diag_pt%qldt_fill) = -ql * inv_dtcloud
            endwhere
          endif
        else
!a        N3D = qn*airdens*1.e-6
          diag_4d(:,:,:,diag_pt%droplets) = qn(:,:,:)*airdens(:,:,:)*1.e-6
!         if (diag_id%debug1_3d   > 0) debug1 = min(qa,1.)
          if (diag_id%cf_liq_init     > 0) diag_4d(:,:,:,diag_pt%cf_liq_init  ) = min(qa,1.)
          do j=1,kdim
           do k=1,jdim
            do i=1,idim
             if (ql(i,k,j).le.qmin .or. qn(i,k,j).le.qmin) then
!a            N3D(i,k,j) = 0.
              diag_4d(i,k,j,diag_pt%droplets) = 0.
!             if (max(diag_id%qldt_fill,diag_id%ql_fill_col) > 0) qldt_fill(i,k,j) = -ql(i,k,j) * inv_dtcloud
              if (max(diag_id%qldt_fill,diag_id%ql_fill_col) > 0) diag_4d(i,k,j,diag_pt%qldt_fill) = -ql(i,k,j) * inv_dtcloud
! Should this just be diag_id%qndt_fill  + diag_id%qn_fill_col ?
!ORIGINAL:    if (diag_id%qldt_fill  > 0) &
              if (max(diag_id%qndt_fill,diag_id%qn_fill_col, &
                           diag_id%qldt_fill,diag_id%ql_fill_col) > 0)  &
!                           qndt_fill(i,k,j) = -qn(i,k,j) * inv_dtcloud
                            diag_4d(i,k,j,diag_pt%qndt_fill) = -qn(i,k,j) * inv_dtcloud
!             if (diag_id%debug1_3d    > 0) debug1(i,k,j) = 0.
              if (diag_id%cf_liq_init      > 0) diag_4d(i,k,j,diag_pt%cf_liq_init  ) = 0.
             endif
            enddo
           enddo
          enddo
        endif
     endif
!rab - end of statements moved outside large loop

!a   diag_4d(:,:,:, diag_pt%droplets) = N3D
!    diag_4d(:,:,:, diag_pt%droplets_wtd) = N3D*ql
     call mpp_clock_end(sc_pre_loop)
     call mpp_clock_begin(sc_loop)
!-----------------------------------------------------------------------
!
!       Enter the large loop over vertical levels.  Level 1 is the top
!       level of the column and level kdim is the bottom of the model.
!       If MASK is present, each column may not have kdim valid levels.
!
        rain3d = 0.
        snow3d = 0.
        snowclr3d = 0.

        DO j = 1, kdim

        sum_freeze = 0.
        sum_rime = 0.
        sum_berg = 0.
        f_snow_berg(:,:,j) = 0.
       
!-----------------------------------------------------------------------
!
!       Calculate pressure thickness of level and relative humidity 

        !calculate difference in pressure across the grid box divided
        !by gravity
        deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
         
        !calculate GRID box mean relative humidity 
        !       U = min(max(0.,qv(:,:,j)/qs(:,:,j)),1.)
 
        U = 0.
        where (qrat(:,:,j) .gt. 0.)
           U = min(max(0.,(qv(:,:,j)/(qrat(:,:,j)*qs(:,:,j)))),1.)
        end where
        
!-----------------------------------------------------------------------
!
!       Account for the fact that other processes may have created
!       negative tracer or extremely small values of tracer fields.
!       The general reason for the extremely small values of the 
!       tracer fields is due to vertical diffusion, advection of 
!       condensate or cumulus induced subsidence (also a form of 
!       advection) of condensate.
!
!       In this step any values of the prognostic variables which are 
!       less than qmin are reset to zero, while conserving total 
!       moisture.
!
!       Note that this is done slightly different for the Tiedtke
!       cloud fraction than it is for pdf clouds. In the former, 
!       the filling requires that cloud liquid, cloud ice, and 
!       cloud fraction are greater than qmin. For PDF clouds, 
!       cloud fraction need not be considered since it is diagnosed
!       below from the PDF clouds
!
    if (.not. do_pdf_clouds) then

        do k=1,jdim
         do i=1,idim
          qa_upd(i,k) = qa(i,k,j)      
          qi_upd(i,k) = qi(i,k,j)      
          ql_upd(i,k) = ql(i,k,j)      

          if (qa(i,k,j) .le. qmin) then
             SA(i,k,j)   = SA(i,k,j) - qa(i,k,j)
             qa_upd(i,k) = 0.
          end if
!       Correct for qa > RH, which is not permitted under the 
!       assumption that the cloudy air is saturated and the temperature 
!       inside and outside of the cloud are about the same.
          if (qa_upd(i,k) .gt. U(i,k)) then
!            if (max(diag_id%qadt_rhred,diag_id%qa_rhred_col) > 0 ) qadt_rhred(i,k,j) = (qa_upd(i,k)-U(i,k)) * inv_dtcloud
             if (max(diag_id%qadt_rhred,diag_id%qa_rhred_col) > 0 ) diag_4d(i,k,j,diag_pt%qadt_rhred) = (qa_upd(i,k)-U(i,k)) * inv_dtcloud
             SA(i,k,j)   = SA(i,k,j) + U(i,k) - qa_upd(i,k)
             qa_upd(i,k) = U(i,k)      
          end if
        
          if (.not. do_liq_num) then
            if (ql(i,k,j) .le. qmin .or. qa(i,k,j) .le. qmin) then
             SL(i,k,j)   = SL(i,k,j) - ql(i,k,j)
             SQ(i,k,j)   = SQ(i,k,j) + ql(i,k,j)
             ST(i,k,j)   = ST(i,k,j) - hlv*ql(i,k,j)/cp_air
             ql_upd(i,k) = 0.
            end if
          else
            qn_upd(i,k) = qn(i,k,j)
            if (ql(i,k,j) .le. qmin .or. qa(i,k,j) .le. qmin .or. qn(i,k,j) .le. qmin) then
             SL(i,k,j)   = SL(i,k,j) - ql(i,k,j)
             SQ(i,k,j)   = SQ(i,k,j) + ql(i,k,j)
             ST(i,k,j)   = ST(i,k,j) - hlv*ql(i,k,j)/cp_air
             SN(i,k,j)   = SN(i,k,j) - qn(i,k,j)
             ql_upd(i,k) = 0.
             qn_upd(i,k) = 0.
            endif
          endif

          if (qi(i,k,j) .le. qmin .or. qa(i,k,j) .le. qmin) then
             SI(i,k,j)   = SI(i,k,j) - qi(i,k,j)
             SQ(i,k,j)   = SQ(i,k,j) + qi(i,k,j)
             ST(i,k,j)   = ST(i,k,j) - hls*qi(i,k,j)/cp_air
             qi_upd(i,k) = 0.
          endif
         enddo
        enddo
        
    else

        do k=1,jdim
         do i=1,idim
          ql_upd(i,k) = ql(i,k,j)
          qi_upd(i,k) = qi(i,k,j)

          if (.not. do_liq_num) then
           if (ql(i,k,j) .le. qmin) then
             SL(i,k,j)   = SL(i,k,j) - ql(i,k,j)
             SQ(i,k,j)   = SQ(i,k,j) + ql(i,k,j)
             ST(i,k,j)   = ST(i,k,j) - hlv*ql(i,k,j)/cp_air
             ql_upd(i,k) = 0.
           endif
          else
           qn_upd(i,k) = qn(i,k,j)
           if (ql(i,k,j) .le. qmin .or. qn(i,k,j) .le. qmin) then
             SL(i,k,j)   = SL(i,k,j) - ql(i,k,j)
             SQ(i,k,j)   = SQ(i,k,j) + ql(i,k,j)
             ST(i,k,j)   = ST(i,k,j) - hlv*ql(i,k,j)/cp_air
             SN(i,k,j)   = SN(i,k,j) - qn(i,k,j)
             ql_upd(i,k) = 0.
             qn_upd(i,k) = 0.
           endif
          endif

          if (qi(i,k,j) .le. qmin) then
             SI(i,k,j)   = SI(i,k,j) - qi(i,k,j)
             SQ(i,k,j)   = SQ(i,k,j) + qi(i,k,j)
             ST(i,k,j)   = ST(i,k,j) - hls*qi(i,k,j)/cp_air
             qi_upd(i,k) = 0.
          endif
         end do
        end do
        
    end if !for do_pdf_clouds
        
        
!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!                 NON-CONVECTIVE CONDENSATION                          !
!                                                                      !
!                                                                      !
!                                                                      !
!                         METHOD 1                                     !
!                                                                      !
!                                                                      !
!                  TIEDTKE (1993) CLOUD FRACTION                       !
!                                                                      !
!       ANALYTIC INTEGRATION OF SATURATED VOLUME FRACTION EQUATION     !
!                                                                      !
!                                                                      !
!
!       Do non-convective condensation following Tiedtke, pages 3044-5.
!       In this formulation stratiform clouds are only formed/destroyed 
!       when there is upward or downward motion to support/destroy it. 
!
!       The first step is to compute the change in qs due to large-
!       scale processes, dqs_ls.   In Tiedtke, it has contributions from 
!       large-scale uplift, convection induced compensating subsidence,
!       turbulence cooling and radiative cooling.  dqs_ls has the form:
!
!               (((omega+ grav*Mc)/airdens/cp)+radturbten)*dqsdT*dtcloud
!   (6) dqs_ls= --------------------------------------------------------
!                  1.  +   ( qa +  (da_ls/2.) ) * gamma
!
!       Here da_ls is the increase in cloud fraction due to non-
!       convective processes.  Because this increase is also a function
!       of dqs_ls, a quadratic equation must be solved for dqs_ls in
!       the case that da_ls is not equal to zero.
!
!       Note that if the PDF cloud scheme is active the Tiedtke large-
!       scale condensation is bypassed.

    if (.not.do_pdf_clouds) then

     do k=1,jdim
      do i=1,idim
        dqs_ls(i,k) =(((omega(i,k,j)+grav*Mc(i,k,j))/airdens(i,k,j)/cp_air)+&
                 radturbten2(i,k,j))*dtcloud*dqsdT(i,k,j)

        !compute pressure dependent U00 following ECMWF formula if 
        !desired
        U00p(i,k) = U00
        if (u00_profile) then
             if (pfull(i,k,j) .gt. 0.8*phalf(i,k,KDIM+1)) then
                    U00p(i,k) = U00 + (1.-U00)* &
                         (((pfull(i,k,j)-(0.8*phalf(i,k,KDIM+1))) &
                                    /    (0.2*phalf(i,k,KDIM+1)) )**2.)
             endif
        end if       

!ljd
!       modify u00p to account for humidity in convective system
!       See "Tiedtke u00 adjustment" notes, 10/22/02
!ljd
 
        u00p(i,k)=u00p(i,k)+(1.-u00p(i,k))*ahuco(i,k,j)

        if (dqs_ls(i,k).le.0. .and. U(i,k).ge.U00p(i,k) .and. qa_upd(i,k).lt.1.) then
             tmp1s = sqrt( (1.+qa_upd(i,k)*gamma(i,k,j))**2. - (1.-qa_upd(i,k)) * &
                    (1.-qa_upd(i,k))*gamma(i,k,j)*dqs_ls(i,k)/qs(i,k,j)/         &
                    max(1.-U(i,k),qmin) ) - (1.+qa_upd(i,k)*gamma(i,k,j))
             tmp1s = -1. * tmp1s / ((1.-qa_upd(i,k))*(1.-qa_upd(i,k))*gamma(i,k,j)/&
                    qs(i,k,j)/max(1.-U(i,k),qmin)/2.)
             dqs_ls(i,k) = min(tmp1s,dqs_ls(i,k)/(1.+0.5*(1.+qa_upd(i,k))*gamma(i,k,j)))
        else
             dqs_ls(i,k) = dqs_ls(i,k)/(1.+qa_upd(i,k)*gamma(i,k,j))
        endif
      
!       The next step is to compute the change in saturated volume
!       fraction due to non-convective condensation, da_ls.   This 
!       occurs in two conditions:
!
!       (a)  dqs_ls < 0. and U00 < U < 1., where U00 is the threshold
!            relative humidity for non-convective condensation. Note 
!            that if U is greater than or equal to 1., ideally qa = 1,
!            and da_ls = 0.  However this may not be the case for 
!            numerical reasons so this must be assured after analytic 
!            integration of the qa equation.
!
!            For these cases the change in saturated volume fraction is:
!
!   (7)      da_ls = - (1.-qa)*(1.-qa)*dqs_ls/2./qs/(1.-U)
!
!            This formula arises from the assumption that vapor is uni-
!            formly distributed in the range [qv_clr - (qs - qv_clr),qs]
!            where qv_clr is the amount of vapor in the unsaturated 
!            volume and is given from the following equation:
!
!   (8)      qv  =   qa * qs      +   (1.-qa) * qv_clr
!          
!            Implicit in equation (7) is the following assumption:
!            As qsat changes, the distribution of qv+ql+qi 
!            remains constant.  That is as qsat rises, portions where
!            qv+ql+qi > qsat+dqsat remain saturated.  This can only
!            occur if it is assumed that ql+qi evaporate-sublimate or
!            condense-deposit to keep qv = qsat. 
!
!       (b)  dqs_ls > 0.  Ideally some portion of the cloud should
!            evaporate however this is not accounted for at present.
!            

        !compute formula for da_ls
!       where (dqs_ls .le. 0. .and. U .ge. U00p)
!            da_ls = -0.5 * (1.-qa_upd) * (1.-qa_upd) * dqs_ls /       &
!             qs(:,:,j) / max(1.-U,qmin)
!       elsewhere
!            da_ls = 0.
!       end where 
 
        da_ls(i,k) = 0.
        if ((dqs_ls(i,k).le.0. .and. U(i,k).ge.U00p(i,k)) .and. &
              (qa_upd(i,k)+ahuco(i,k,j).le.1.)) then
             da_ls(i,k) = -0.5 * (1.-qa_upd(i,k)-ahuco(i,k,j)) * (1.-qa_upd(i,k)-    &
                ahuco(i,k,j))    * dqs_ls(i,k)/   qs(i,k,j) / max(1.-U(i,k),qmin)
        endif
      enddo
     enddo

!       Turbulent erosion of clouds
!
!       As in Tiedtke (1993) this is calculated using the eros_scale
!       parameter as:
!
!   (9) dql/dt    =  - qa * eros_scale * (qs - qv) * (ql/ ql+qi )
!
!  (10) dqi/dt    =  - qa * eros_scale * (qs - qv) * (qi/ ql+qi )
!
!  (11) dqa/dt    =  - qa * eros_scale * (qs - qv) * (qa/ ql+qi )
!
!       for which the erosion sink term (B in equation 13) is
!
!  (12) B = qa * eros_scale * (qs - qv) / (ql + qi)  
!
!
!       Theory for eros_scale
!
!       If eros_choice equals false, then a single erosion time scale
!       is used in all conditions (eros_scale).  If eros_choice equals
!       true then it is assumed that the timescale for turbulent 
!       evaporation is a function of the conditions in the grid box.  
!       Specifically, if the flow is highly turbulent then the scale is 
!       short, and eros_scale is large.  Likewise if convection is 
!       occurring, then it is assumed that the erosion term is larger 
!       than backround conditions. 
!
!       Here are the typical values for the timescales and the 
!       switches used:
!
!         Mixing type      eros_scale (sec-1)          Indicator
!       ----------------   ------------------     --------------------
!
!       Background            1.e-06              always present
!       Convective layers     5.e-06              Mc > Mc_thresh
!       Turbulent  layers     5.e-05              diff_t > diff_thresh
!

     do k=1,jdim
      do i=1,idim
        !Background erosion scale
        tmp2s = eros_scale

        !Do enhanced erosion in convective or turbulent layers?
        !
        !              IMPORTANT NOTE
        !                
        !Note that convection is considered first, so that if 
        !turbulence and convection occur in the same layer, the
        !erosion rate for turbulence is selected.                
        !
        
        if (eros_choice) then
             !Enhanced erosion in convective layers
             if (Mc(i,k,j) .gt. mc_thresh) tmp2s = eros_scale_c
        
             !Enhanced erosion in turbulent layers
             if (diff_t(i,k,j).gt.diff_thresh  .or.  &
                 diff_t(i,k,min(j+1,KDIM)).gt.diff_thresh) &
                   tmp2s = eros_scale_t
        end if   !for erosion choice

        if (ql_upd(i,k) .gt. qmin .or. qi_upd(i,k) .gt. qmin) then
          D_eros(i,k)=qa_upd(i,k) * tmp2s * dtcloud * qs(i,k,j) *      &
                      (1.-U(i,k)) / (qi_upd(i,k) + ql_upd(i,k))
          if (pfull(i,k,j) .gt. 400.e02) then
             D_eros(i,k)=D_eros(i,k)+efact*D_eros(i,k)*((pfull(i,k,kdim)-  &
                       pfull(i,k,j))/(pfull(i,k,kdim)-400.e02))
          else
             D_eros(i,k)=D_eros(i,k)+efact*D_eros(i,k)
          endif
        else
          D_eros(i,k) = 0.
        endif
      enddo
     enddo
     
!    
!       The next step is to analytically integrate the saturated volume
!       fraction equation.  This follows the Tiedtke approach
!
!       The qa equation is written in the form:
!
!  (13) dqa/dt    =   (1.-qa) * A   -  qa * B 
!
!       Note that over the physics time step, A, B are assumed to be 
!       constants.
!
!       Defining qa(t) = qa0 and qa(t+dtcloud) = qa1, the analytic
!       solution of the above equation is:
!
!  (14) qa1 = qaeq -  (qaeq - qa0) * exp (-(A+B)*dtcloud)
! 
!       where qaeq is the equilibrium cloud fraction that is approached
!       with an time scale of 1/(A+B),
!
!  (15) qaeq  =  A/(A+B)
!
!
!       To diagnose the magnitude of each of the right hand terms of
!       (13) integrated over the time step, define the average cloud
!       fraction in the interval t to t + dtcloud qabar as:
!
!  (16) qabar  = qaeq - [ (qa1-qa0) / ( dtcloud * (A+B) ) ]
! 
!       from which the magnitudes of the A and B terms integrated
!       over the time step are:
!
!       A * (1-qabar)    and    -B * (qabar)
!
!       Additional notes on this analytic integration:
!
!       1.   For large-scale cloud formation or destruction from 
!            the dqs_ls term the contributions to A or B are defined
!            from:
!
!  (19)      A_ls * (1. - qa) = da_ls / dtcloud      if da_ls >= 0.
! 
!  (20)      B_ls * qa        = da_ls / dtcloud      if da_ls < 0.
!
!       2.   Note that to reduce the number of variables, the following
!            equivalency exists:
!
!               Ql or Qi equation              Qa equation
!             --------------------         -------------------
! 
!                     C_dt                        A_dt
!                     D_dt                        B_dt
!                     qceq                        qaeq
!                     qcbar                       qabar
!                     qc1                         qa1
!                     qc0                         qa0
!
!       3.   Qa goes to zero only in the case of ql and qi less than or
!            equal to qmin; see 'cloud destruction code' near the end of 
!            this loop over levels.
!
     do k=1,jdim
      do i=1,idim
        !compute C_dt; This is assigned to the large-scale source term
        !following (18). Reset D_dt.
        C_dts = da_ls(i,k)/max((1.-qa_upd(i,k)),qmin)
        D_dts = D_eros(i,k)
  
        !do analytic integration      
        qc0s   = qa_upd(i,k)
        if ( (C_dts.gt.Dmin) .or. (D_dts.gt.Dmin) ) then
             qceqs  = C_dts  / (C_dts + D_dts)
             qc1s   = qceqs - (qceqs - qc0s) * exp ( -1.*(C_dts+D_dts) )
             qcbars = qceqs - ((qc1s - qc0s)/ (C_dts + D_dts))
        else
             qceqs  = qc0s   
             qc1s   = qc0s   
             qcbars = qc0s  
        endif

        !set total tendency term and update cloud fraction    
        if (limit_conv_cloud_frac) then
!RSH     limit cloud area to be no more than that which is not being
!        taken by convective clouds
          qc1s = MIN(qc1s, 1.0 -ahuco(i,k,j))
        endif
        SA(i,k,j)  = SA(i,k,j) + qc1s - qc0s
        qa_upd(i,k)     = qc1s
        
!       if (max(diag_id%qadt_lsform,diag_id%qa_lsform_col) > 0) qadt_lsform(i,k,j) =  C_dts * (1.-qcbars) * inv_dtcloud 
        if (max(diag_id%qadt_lsform,diag_id%qa_lsform_col) > 0) diag_4d(i,k,j,diag_pt%qadt_lsform) =  C_dts * (1.-qcbars) * inv_dtcloud 
!       if (max(diag_id%qadt_eros,diag_id%qa_eros_col)  > 0) qadt_eros  (i,k,j) =  D_dts *     qcbars  * inv_dtcloud
        if (max(diag_id%qadt_eros,diag_id%qa_eros_col)  > 0) diag_4d(i,k,j,diag_pt%qadt_eros) =  D_dts *     qcbars  * inv_dtcloud
        delta_cf(i,k) = C_dts * (1.-qcbars)

!       The next step is to calculate the change in condensate
!       due to non-convective condensation, dcond_ls. Note that this is
!       not the final change but is used only to apportion condensate
!       change between phases. According to Tiedtke 1993 this takes the
!       form:
!
!  (21) dcond_ls = -1. * (qa +  0.5*da_ls) * dqs_ls
!
!       Here the 0.5*da_ls represents using a midpoint cloud fraction.
!       This is accomplished by using the variable qcbar.

        dcond_ls(i,k) = -1. * qcbars * dqs_ls(i,k)
      enddo
     enddo
            

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!                 NON-CONVECTIVE CONDENSATION                          !
!                                                                      !
!                                                                      !
!                                                                      !
!                         METHOD 2                                     !
!                                                                      !
!                                                                      !
!                STATISTICAL CLOUD FRACTION                            !
!                                                                      !
!                                                                      !

    else                   !for doing PDF cloud scheme
    

        !set Tiedtke erosion term to zero            
        D_eros = 0.
        
        !compute pdf cloud fraction and condensate
        ! 
        !Note that the SYMMETRIC beta distribution is used here.
        !
        !
        ! Initialize grid-box mean values of cloud fraction (qag),
        ! cloud condensate(qcg), and clear sky water vapor (qvg)

        qcg = 0.
        qvg = 0.
        qag = 0.
        
        !Create loop over sub-levels within a grid box
        do ns = 1, nsublevels
        
             !calculate normalized vertical level
             ! 0. = top of gridbox
             ! 1. = bottom of gridbox
        
             pnorm =  (real(ns) - 0.5 )/real(nsublevels)
        
             !First step is to calculating the minimum (qtmin)
             !of the total water distribution and 
             !the width of the qt distribution (deltaQ)
             !
             !For diagnostic variance this is set to (1.-qthalfwidth)*qtbar
             !and 2*qthalfwidth*qtbar, respectively, where qtbar is the
             !mean total water in the grid box.        
             !
             !

             qtbar = qta4(2,:,:,j)+pnorm*( (qta4(3,:,:,j)-qta4(2,:,:,j)) + &
                                       qta4(4,:,:,j)*(1-pnorm) )
             
             qtbar = max(qmin,qtbar)
             deltaQ = 2.*qthalfwidth*qtbar
             qtmin = (1.-qthalfwidth)*qtbar
        
             !From this the variable normalized saturation specific
             !humidity qs_norm is calculated.
             !
             !  qs_norm = (qs(Tl) - qtmin)/(qtmax-qtmin)
             !
             !          = 0.5  - (qtbar - qs(Tl))/deltaQ
             !
             !Note that if qs_norm > 1., the grid box is fully clear.
             !If qs_norm < 0., the grid box is fully cloudy.
        
             qs_norm = qtqsa4(2,:,:,j)+  &
                       pnorm*( (qtqsa4(3,:,:,j)-qtqsa4(2,:,:,j)) + &
                       qtqsa4(4,:,:,j)*(1-pnorm) )
      
             qs_norm = 0.5 - ( qs_norm/deltaQ )
             
             !Calculation of cloud fraction (qagtmp), cloud condensate 
             !(qcgtmp), and water vapor in clear air part of the grid 
             !box (qvgtmp)
             !
             !Formulas (from Tompkins, and personal derivations):
             !
             !  Define icbp  = incomplete_beta(qs_norm,p,q)
             !         icbp1 = incomplete_beta(qs_norm,p+1,q)
             !
             !  qagtmp = 1. - icbp
             !
             !  qcgtmp = aThermo * {  (qtbar-qtmin)*(1.-icbp1) - 
             !                       qs_norm*deltaQ*(1.-icbp ) }
             !
             !
             !  qvgtmp = qtmin + (p/(p+q))*(icbp1/icbp)*deltaQ
             !
             !  
             ! where aThermo = 1./(1.+(L/cp)*dqsdT)
             !
             ! note that in the qvg formula below the factor of 0.5
             ! is equal to (p/(p+q)).
             !

             do jd = 1,jdim
             do id = 1,idim
        
             if (qs_norm(id,jd).le.1.) then
                 
                 icbp = incomplete_beta(max(0.,qs_norm(id,jd)), &
                                      p = betaP    , q = betaP)
                 icbp1= incomplete_beta(max(0.,qs_norm(id,jd)), &
                                      p = betaP + 1, q = betaP)
                 qagtmps = 1.-icbp
                 qcgtmps = (qtbar(id,jd)-qtmin(id,jd))*(1.-icbp1)&
                               - qs_norm(id,jd)*deltaQ(id,jd)*(1.-icbp)    
                 qcgtmps = qcgtmps/(1.+gamma(id,jd,j))
                 qvgtmps = qtmin(id,jd) + &
                               0.5*(icbp1/max(icbp,qmin))*deltaQ(id,jd)
             
                 !bound very very small cloud fractions which may
                 !cause negative cloud condensates due to roundoff 
                 !errors or similar errors in the beta table lookup.
                 if((qagtmps.lt.0.).or.(qcgtmps.le.0.))then
                      qagtmps = 0.
                      qcgtmps = 0.
                      qvgtmps = qtbar(id,jd)
                 end if
                 
             else             
                 qagtmps = 0.
                 qcgtmps = 0.
                 qvgtmps = qtbar(id,jd)             
             end if
    
             !sum vertically
             !
             !note special averaging of clear-sky water vapor
             !this is weighting clear-sky relative humidity by the 
             !clear-sky fraction
         
             qag(id,jd) = qag(id,jd) + qagtmps
             qcg(id,jd) = qcg(id,jd) + qcgtmps
             qvg(id,jd) = qvg(id,jd)+(1.-qagtmps)*min(max(qvgtmps/max(qmin, &
                       (qtbar(id,jd)+((qs_norm(id,jd)-0.5)*deltaQ(id,jd)))),0.),1.)
        
             !compute grid-box average cloud fraction, cloud condensate
             !and water vapor
        
             if (nsublevels.gt.1 .and. ns.eq.nsublevels) then
                  qag(id,jd) = qag(id,jd) / real(nsublevels)
                  qcg(id,jd) = qcg(id,jd) / real(nsublevels)
             
                  !note special averaging of clear-sky water vapor
                  if ((1.-qag(id,jd)).gt.qmin) then
                      qvg(id,jd) =qvg(id,jd)/real(nsublevels)/&
                                  (1.-qag(id,jd))
                      qvg(id,jd) =qvg(id,jd) * qs(id,jd,j)
                  else
                      qvg(id,jd) = qs(id,jd,j)
                  end if
             elseif (nsublevels.eq.1) then
                  ! for nsublevels = 1, qag and qcg already hold their
                  ! final values
                  qvg(id,jd) = qvgtmps
             end if
             enddo
             enddo
             
        enddo !for number of sublevels loop
             
        !do adjustment of cloud fraction
!rab        qc0 = qa(:,:,j)
!rab        qc1 = qag

        !set total tendency term and update cloud fraction    
        SA(:,:,j)  = SA(:,:,j) + qag - qa(:,:,j)
        qa_upd     = qag

!       if (max(diag_id%qadt_lsform,diag_id%qa_lsform_col) > 0) qadt_lsform(:,:,j) =  max(qag-qa(:,:,j),0.) * inv_dtcloud 
        if (max(diag_id%qadt_lsform,diag_id%qa_lsform_col) > 0) diag_4d(:,:,j,diag_pt%qadt_lsform) =  max(qag-qa(:,:,j),0.) * inv_dtcloud 
!       if (max(diag_id%qadt_lsdiss,diag_id%qa_lsdiss_col) > 0) qadt_lsdiss(:,:,j) =  max(qa(:,:,j)-qag,0.) * inv_dtcloud
        if (max(diag_id%qadt_lsdiss,diag_id%qa_lsdiss_col) > 0) diag_4d(:,:,j,diag_pt%qadt_lsdiss) =  max(qa(:,:,j)-qag,0.) * inv_dtcloud

        !define da_ls and delta_cf needed when do_liq_num = .true. (cjg)
        da_ls = max(qag-qa(:,:,j),0.)
        delta_cf = max(qag-qa(:,:,j),0.)

        !compute large-scale condensation / evaporation
        dcond_ls = qcg - (ql_upd + qi_upd)

    end if    !for doing PDF clouds

!       The next step is the apportionment on the non-convective 
!       condensation between liquid and ice phases. Following the
!       suggestion of Rostayn (2000), all condensation at temperatures
!       greater than -40C is in liquid form as ice nuclei are generally 
!       limited in the atmosphere. The droplets may subsequently be 
!       converted to ice by the Bergeron-Findeisan mechanism.  
!
!       One problem with this formulation is that the proper saturation
!       vapor pressure is not used for cold clouds as it should be
!       liquid saturation in the case of first forming liquid, but
!       change to ice saturation as the cloud glaciates.  The current
!       use of ice saturation beneath -20C thus crudely mimics the
!       result that nearly all stratiform clouds are glaciated for
!       temperatures less than -15C.
!
!       In the case of large-scale evaporation (dcond_ls<0.), it is
!       assumed that cloud liquid will evaporate faster than cloud
!       ice because if both are present in the same volume the
!       saturation vapor pressure over the droplet is higher than 
!       that over the ice crystal.
!
!       The fraction of large-scale condensation that is liquid
!       is stored in the temporary variable tmp1.   

        do k=1,jdim
         do i=1,idim
          !assume liquid fractionation 
          tmp1s = 1.

          if (dcond_ls(i,k) .ge. 0.) then

           !For cases of cloud condensation where temperatures are
           !less than -40C create only ice
           if (T(i,k,j) .lt. tfreeze-40.) then
             tmp1s = 0.
           endif

          else

           if (qi_upd(i,k).gt.qmin) then

            if (ql_upd(i,k).gt.qmin) then
             !For cases of cloud evaporation of mixed phase clouds
             !set liquid evaporation to preferentially occur first
             tmp1s = min(-1.*dcond_ls(i,k),ql_upd(i,k))/max(-1.*dcond_ls(i,k),qmin)
            else
             !do evaporation of pure ice cloud
             tmp1s = 0.
            endif

           endif

          endif
          !calculate partitioning among liquid and ice to dcond_ls
          dcond_ls_ice(i,k) = (1.-tmp1s) * dcond_ls(i,k)
          dcond_ls(i,k)     = tmp1s * dcond_ls(i,k)      
         enddo
        enddo

!       The next step is to compute semi-implicit qa,ql,qi which are 
!       used in many of the formulas below.  This gives a somewhat 
!       implicitness to the scheme. In this calculation an estimate 
!       is made of what the cloud fields would be in the absence of 
!       cloud microphysics and cloud erosion.
!
!       In the case of the Tiedtke cloud scheme, the mean cloud 
!       condensate is incremented if large-scale condensation is 
!       occuring. For cloud fraction, the value from the analytic 
!       integration above is used.
!
!       For the statistical cloud scheme these are set equal to the
!       values diagnosed from the beta-distribution apart from the
!       corrections for mixed phase clouds.

        qa_mean = qa_upd
        if (.not. do_pdf_clouds) then
             ql_mean = ql_upd + max(dcond_ls    ,0.)        
             qi_mean = qi_upd + max(dcond_ls_ice,0.)
        else
             ql_mean = max(ql_upd + dcond_ls    ,qmin)
             qi_mean = max(qi_upd + dcond_ls_ice,qmin)  
        end if
        if (do_liq_num) then 
!yim's CCN activation
          do k = 1,jdim
            do i = 1,idim
              if ( da_ls(i,k) > 0.0 ) then
                up_strat = -1.*(((omega(i,k,j)+grav*Mc(i,k,j))/ &
                                           airdens(i,k,j)/grav) + &
                                  radturbten2(i,k,j)*cp_air/grav)
                          
!-->cjg: modification
!               call aer_ccn_act (T(i,k,j), pfull(i,k,j), up_strat, &
!                                 totalmass1(i,k,j,:), drop1(i,k))
                thickness = deltpg(i,k) / airdens(i,k,j)
                wp2 = 2.0/(3.0*0.548**2)* &
                      (0.5*(diff_t(i,k,j) + diff_t(i,k,min(j+1,KDIM)))/&
                                                         thickness )**2
                wp2 = MAX (wp2, var_limit**2)
!rab take care of it when writing diags....
!rab                debug2(i,k,j) = wp2**0.5
!               if (diag_id%debug2_3d > 0) debug2(i,k,j) = wp2
                if (diag_id%subgrid_w_variance > 0 .and. up_strat > 0.0) diag_4d(i,k,j,diag_pt%subgrid_w_variance) = wp2
                call aer_ccn_act_wpdf (T(i,k,j), pfull(i,k,j), &
                                       up_strat, wp2,    &
                                       totalmass1(i,k,j,:), drop1(i,k))
!               if (diag_id%debug3_3d > 0) debug3(i,k,j) = drop1(i,k)
                if (diag_id%potential_droplets > 0 .and. up_strat > 0.0) diag_4d(i,k,j,diag_pt%potential_droplets) = drop1(i,k)
!               if (diag_id%debug2_3d > 0) debug2(i,k,j) = 1.
!<--cjg: end of modification
                qn_mean(i,k) = qn_upd(i,k) + max(delta_cf(i,k),0.)*  &
                               drop1(i,k)*1.e6/airdens(i,k,j)
              else
                drop1(i,k) = 0.                
                qn_mean(i,k) = qn_upd(i,k)
                if (diag_id%subgrid_w_variance > 0) diag_4d(i,k,j,diag_pt%subgrid_w_variance) =  0.0
                if (diag_id%potential_droplets > 0) diag_4d(i,k,j,diag_pt%potential_droplets) = 0.0        
              endif
            end do
          end do        
        endif


    
        !compute diagnostics for cloud fraction
!       if (diag_id%aall > 0) areaall(:,:,j) = qa_mean
        if (diag_id%aall > 0) diag_4d(:,:,j, diag_pt%aall) = qa_mean
!       if (diag_id%aliq > 0 .or. diag_id%rvolume > 0) then
        if (max(diag_id%aliq,diag_id%rvolume) > 0) then
!             where (ql_mean .gt. qmin) arealiq(:,:,j) = qa_mean
              where (ql_mean .gt. qmin) diag_4d(:,:,j,diag_pt%aliq) = qa_mean
        end if
!       if (diag_id%aice > 0 .or. diag_id%vfall > 0) then
        if (max(diag_id%aice,diag_id%vfall) > 0) then
!             where (qi_mean .gt. qmin) areaice(:,:,j) = qa_mean
              where (qi_mean .gt. qmin) diag_4d(:,:,j,diag_pt%aice) = qa_mean
        end if              

  
!-----                                                            -----! 
!                                                                      !
!                                  END OF                              !
!                                                                      !
!                        NON-CONVECTIVE CONDENSATION                   !
!                                                                      !
!                                 SECTION                              !
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!
!
!            TRANSFER OF RAIN FLUXES AND SNOW FLUXES BETWEEN 
!      
!           SATURATED AND UNSATURATED PORTIONS OF THE GRID BOX
!
!                           AT LAYER INTERFACES
!                           
!
!       Do transfer of rain and snow fluxes between clear and cloud 
!       portions of the grid box. This formulism follows the rain/snow
!       parameterization developed by Jakob and Klein (2000).  It 
!       treats the grid mean flux of rain and snow separately according 
!       to the portion that exists in unsaturated air and the portion 
!       that exists in saturated air.  At the top of each level, some 
!       precipitation that was in unsaturated air in the levels above 
!       may enter saturated air and vice versa. These transfers are 
!       calculated here under an assumption for the overlap between 
!       precipitation area and saturated area at the same and adjacent 
!       levels.
!
!       For rain, the area of the grid box in which rain in saturated 
!       air of the level above enters unsaturated air in the current 
!       level is called da_cld2clr and is given by 
!       maximum(a_rain_cld - qa , 0.) in the case of maximum overlap of
!       rain_cld and qa, but equal to a_rain_cld*(1.-qa) in the case of
!       random overlap of cloud and precipitation. The grid mean flux of 
!       rain which is transfered from saturated air to unsaturated air 
!       is given by rain_cld*da_cld2clr/a_rain_cld.
!
!       The area of the grid box where rain in unsaturated air of the
!       level above enters saturated air in the current level is
!       called da_clr2cld and is given by
!       maximum(0.,mininum(a_rain_clr,qa(j)-qa(j-1))) in the case of
!       maximum overlap of cloud and rain_clr, but equal to a_rain_clr*
!       qa for the case of random overlap of cloud and precipitation.  
!       The grid mean flux of rain transfered from unsaturated air to 
!       saturated air is given by rain_clr*da_clr2cld/a_rain_clr.
!
!       NOTE: the overlap assumption used is set by the namelist 
!             variable in cloud_rad
!
!       Overlap values are:      1  = maximum
!                                2  = random
!
!       If cloud_generator is on, the overlap choice is taken from
!       there, by computing a quantity alpha, which is weighting 
!       between maximum and random overlap solutions.
!
!                alpha = 1 ---> maximum,
!                alpha = 0 ---> random
!
!       alpha is stored in tmp3.
!       tmp1 has the maximum overlap solution
!       tmp2 has the random overlap solution
       
!
        if (cloud_generator_on) then
             if (j.gt.1) then
                  tmp3 = compute_overlap_weighting(qa_mean_lst,qa_mean,&
                         pfull(:,:,j-1),pfull(:,:,j))
                  tmp3 = min(1.,max(0.,tmp3))       
             else
                  tmp3 = 0.0
             end if
        end if      
!
!
!       Rain transfers are done first
!

        call cloud_clear_xfer (cloud_generator_on, tmp3, qa_mean, qa_mean_lst, a_rain_clr, a_rain_cld, rain_clr, rain_cld)
!
!       Snow transfers are done second, in a manner exactly like that
!       done for the rain fluxes
!

        call cloud_clear_xfer (cloud_generator_on, tmp3, qa_mean, qa_mean_lst, a_snow_clr, a_snow_cld, snow_clr, snow_cld)

   
               
!-----------------------------------------------------------------------
!
!
!                        MELTING OF CLEAR SKY SNOW FLUX
!
!
!       Melting of falling ice to rain occurs when T > tfreeze. The 
!       amount of melting is limited to the melted amount that would 
!       cool the temperature to tfreeze.
!
!       In the snowmelt bug version, the temperature of melting was 
!       tfreeze + 2. like the original Tiedtke (1993) paper, instead of 
!       tfreeze.

        snow_fact=0.
        if (do_old_snowmelt) snow_fact=2.
     do k=1,jdim
      do i=1,idim
        !compute grid mean change in snow flux to cool the
        !grid box to tfreeze and store in temporary variable tmp1
        tmp1s = cp_air*(T(i,k,j)-tfreeze-snow_fact)*deltpg(i,k)*inv_dtcloud/hlf
        
        ! If snow_clr > tmp1, then the amount of snow melted is
        ! limited to tmp1, otherwise melt snow_clr.  The amount
        ! melted is stored in tmp2
        tmp2s = max(min(snow_clr(i,k),tmp1s),0.)     

        ST(i,k,j) = ST(i,k,j) - hlf*tmp2s*dtcloud/deltpg(i,k)/cp_air                
        rain_clr(i,k)  = rain_clr(i,k) + tmp2s
        
        !raise a_rain_clr to a_snow_clr IF AND only IF melting occurs
        !and a_rain_clr < a_snow_clr
        if (tmp2s .gt. 0. .and. a_snow_clr(i,k) .gt. qmin)  &
             a_rain_clr(i,k) = max(a_rain_clr(i,k),a_snow_clr(i,k))

        ! If all of the snow has melted, then zero out a_snow_clr
        if (snow_clr(i,k).lt.tmp1s .and. a_snow_clr(i,k).gt.qmin) then
             snow_clr(i,k) = 0.
             a_snow_clr(i,k) = 0.
        else
             snow_clr(i,k) = snow_clr(i,k) - tmp2s          
        endif

!       if (max(diag_id%snow_melt,diag_id%snow_melt_col) > 0) snow_melt(i,k,j) = tmp2s/deltpg(i,k)             
        if (max(diag_id%snow_melt,diag_id%snow_melt_col) > 0) diag_4d(i,k,j,diag_pt%snow_melt) = tmp2s/deltpg(i,k)             
             
!-----------------------------------------------------------------------
!
!
!                        MELTING OF CLOUDY SKY SNOW FLUX
!
!
!       Melting of falling ice to rain occurs when T > tfreeze. The 
!       amount of melting is limited to the melted amount that would 
!       cool the temperature to tfreeze.
!

        if (.not.do_old_snowmelt) then

        !compute grid mean change in snow flux to cool the
        !grid box to tfreeze and store in temporary variable tmp1
        !
        !note that tmp1 already has the value of this variable 
        !from the clear-sky melt calculation, so one does not need
        !to repeat the calculation here.
        !
        !However, note that clear-sky snow melt may have already 
        !reduced the temperature of the grid box - this snow melt is in 
        !variable tmp2 from lines above. Thus the amount that one
        !can melt is less.
        
        tmp1s = tmp1s - tmp2s
        
        ! If snow_cld > tmp1, then the amount of snow melted is
        ! limited to tmp1, otherwise melt snow_cld.  The amount
        ! melted is stored in tmp2
        tmp2s = max(min(snow_cld(i,k),tmp1s),0.)     

        ST(i,k,j) = ST(i,k,j) - hlf*tmp2s*dtcloud/deltpg(i,k)/cp_air                
        rain_cld(i,k)  = rain_cld(i,k) + tmp2s
        
        !raise a_rain_cld to a_snow_cld IF AND only IF melting occurs
        !and a_rain_cld < a_snow_cld
        if (tmp2s .gt. 0. .and. a_snow_cld(i,k) .gt. qmin) &
             a_rain_cld(i,k) = max(a_rain_cld(i,k),a_snow_cld(i,k))

        ! If all of the snow has melted, then zero out a_snow_cld
        if (snow_cld(i,k).lt.tmp1s .and. a_snow_cld(i,k).gt.qmin) then
             snow_cld(i,k) = 0.
             a_snow_cld(i,k) = 0.
        else
             snow_cld(i,k) = snow_cld(i,k) - tmp2s          
        endif

!       if (max(diag_id%snow_melt,diag_id%snow_melt_col) > 0) snow_melt(i,k,j) =  snow_melt(i,k,j) + &
!                                               tmp2s/deltpg(i,k)
        if (max(diag_id%snow_melt,diag_id%snow_melt_col) > 0) diag_4d(i,k,j,diag_pt%snow_melt) =  diag_4d(i,k,j,diag_pt%snow_melt) + &
                                                tmp2s/deltpg(i,k)

        end if  !for snowmelt bugfix
       enddo
      enddo
                            
!----------------------------------------------------------------------!
!
!
!              COMPUTE SLOPE FACTOR FOR ICE MICROPHYSICS                  
!
!       [The following microphysics follows that of Rotstayn (1997)]
!       The number concentration of ice crystals of diameter D in the
!       SIZE interval D to D+dD is assumed to be 
!       distributed as in a Marshall Palmer distribution :
!
!  (22) N(D)dD = Nof * Exp( - lamda_f * D)
!
!       The slope factor and intercept are not assumed to be constant,
!       but the slope factor is
!
!  (23) lamda_f = 1.6X10^(3.+0.023(tfreeze-T))
!
!       Integration of (22) over all particle sizes with a constant
!       density of ice crystals , rho_ice, and assumed spherical shape
!       yields a relationship between the intercept parameter and qi
!
!  (24) Nof = airdens*qi_local*(lamda_f^4)/pi*rho_ice
!       
!
!       For the calculation of riming and sublimation of snow, 
!       lamda_f is needed, so it is calculated here.
!
!       Also qi_mean is updated here with the flux of ice that falls 
!       into the cloudy portion of the grid box from above. This permits
!       the Bergeron and riming process to know about the ice that falls
!       into the grid box in the same time step.

        !Increment qi_mean by the ice flux entering the
        !the grid box. To convert ice_flux to units of condensate by
        !multiply by dtcloud and dividing by the mass per unit area 
        !of the grid box. Implicit here is the assumption that the
        !ice entering the cloud will be spread instantaneously over
        !all of the cloudy area.
        qi_mean = qi_mean + snow_cld*dtcloud/deltpg        

        !snow falling into cloud reduces the amount that
        !falls out of cloud: a loss of cloud ice from settling
        !is defined to be positive
!       if (max(diag_id%qidt_fall,diag_id%qi_fall_col) > 0) qidt_fall(:,:,j)= -1.*snow_cld/deltpg
        if (max(diag_id%qidt_fall,diag_id%qi_fall_col) > 0) diag_4d(:,:,j,diag_pt%qidt_fall)= -1.*snow_cld/deltpg
         
        !compute lamda_f
        lamda_f = 1.6 * 10**(3.+0.023*(tfreeze-T(:,:,j)))
        

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!                       LIQUID PHASE MICROPHYSICS                      !
!                                                                      !
!                                 AND                                  !
!                                                                      !
!                  ANALYTIC INTEGRATION OF QL EQUATION                 !
!                                                                      !
!                                                                      !
!                                                                      !
!       Accretion
!
!       The parameterization of collection of cloud liquid drops by
!       rain drops follows the parameterization of Rotstayn (1997).
!       The parameterization starts with the continous-collection 
!       equation of drops by a rain drop of diameter D, falling with
!       speed V(D).   The fall speed of rain drops is taken from
!       Gunn and Kinzer(1949):
!
!  (25) V(D) = 141.4 (m**0.5,s-1) * sqrt(D) * (rho_ref/airdens)**0.5
!
!       where D is the radius of the rain drops in meters. Here
!       rho_ref is a reference air density.  This formula is generally
!       good for 1 mm < D < 4 mm.
!
!       The distribution of rain drops by SIZE follows a Marshall-
!       Palmer distribution:
!
!  (26) N(D) dD = Nor * Exp (- lamda *D)
!
!       where N(D)dD is the number of rain drops with SIZE D in the
!       interval D to D+dD, Nor is the intercept (assumed fixed at
!       8E+04 (1/m*m*m*m).  lamda is the slope intercept parameter
!       and with (21) it can be shown that lamda is a function of
!       the rain rate.
!
!       With these assumptions the local rate of accretion of cloud
!       liquid reduces to:
!
!  (27) dl/dt_local = - CB*Eco*((rain_rate_local/dens_h2o)**(7/9))*
!                        
!                       (rho_ref/airdens)**(1/9)   * ql_local
!
!       where CB is the accretion constant:
! 
!       CB = 65.772565 [(m)**(-7/9)] * [(s)**(-2/9)]
!
!       AND Eco is the collection efficiency of a cloud droplet by a 
!       rain droplet.   A good fit to the Table 8.2 of Rogers and Yau 
!       (1988) for rain drops between SIZE 1mm and 3mm is:
!
!  (28) Eco = rad_liq**2 / (rad_liq**2 + 20.5 microns**2) .
!
!       In generalizing to grid mean conditions (27) becomes:
!
!  (29) dl/dt = - (arain_cld/qa_mean) * CB * Eco * 
!
!                 [(rain_cld/a_rain_cld/dens_h2o)**(7/9)] * ql
!        
!       Note that the very weak dependence on air density is
!       neglected at this point.
!
!       The particle sizes are computed from the following equation
!
!  (30) rad_liq = (3*airdens*ql/(4*pi*liq_dens*qa*N)^(1/3)
!
!       
!       For numerical treatment we write (25) as:
!
!       dl/dt = - D_acc * l 
!
!       and if we do so:
!
!  (31) D_acc   =  (arain_cld/qa_mean) * CB * Eco * 
!
!                 [(rain_cld/a_rain_cld/dens_h2o)**(7/9)] 
!
!       In the work below, D_acc is added to D1_dt, the first processes
!       contributing to the depletion of cloud liquid in the analytic
!       integration.  D1_dt represents the conversion of liquid to rain.

        if (.not. do_liq_num) then
        !compute rad_liq.  The constant below is equal to  
        !1.E+06 * (3/4*pi)^(1/3), where the 1E+06 is
        !the factor to convert meters to microns.
        
          rad_liq= 620350.49 *( (airdens(:,:,j)*ql_mean/ &
                        max(qa_mean,qmin)/N/dens_h2o)**(1./3.))

        !do not let very small cloud fractions contribution to
        !autoconversion or accretion
          where (qa_mean .le. qmin) rad_liq = 0.
        else        
!yim The 1st place droplet number is used
          rad_liq= 620350.49 *( (ql_mean/max(qn_mean,qmin)/dens_h2o)**(1./3.))
          !do not let very small cloud fractions contribution to
          !autoconversion or accretion
          where (qa_mean .le. qmin .or. qn_upd .le.qmin) rad_liq = 0.
        endif
        
!       if (diag_id%rvolume > 0) rvolume(:,:,j) = rad_liq*arealiq(:,:,j)
!       if (diag_id%rvolume > 0) rvolume(:,:,j) = rad_liq*diag_4d(:,:,j,diag_pt%aliq)
!       if (diag_id%rvolume > 0) diag_4d(:,:,j,diag_pt%rvolume) = rad_liq*arealiq(:,:,j)
        if (diag_id%rvolume > 0) diag_4d(:,:,j,diag_pt%rvolume) = rad_liq*diag_4d(:,:,j,diag_pt%aliq)

        !compute accretion D term
        D1_dt =  dtcloud * 65.772565 * (a_rain_cld/max(qa_mean,qmin))* &
                 ( rad_liq*rad_liq / (rad_liq*rad_liq+20.5) ) *        &
                 ((rain_cld/max(a_rain_cld,qmin)/dens_h2o)**(7./9.))
            
!       if (max(diag_id%qldt_accr,diag_id%ql_accr_col)  > 0) qldt_accr(:,:,j) = D1_dt
        if (max(diag_id%qldt_accr,diag_id%ql_accr_col)  > 0) diag_4d(:,:,j,diag_pt%qldt_accr) = D1_dt
    
!       Autoconversion
!
!       The autoconversion parameterization follow that of Manton
!       and Cotton (1977).  This formula has been used in Chen and
!       Cotton (1987) and is used in the CSIRO GCM (Rotstayn 1997)
!       and the LMD GCM (Boucher, Le Treut and Baker 1995).  In this
!       formulation the time rate of change of grid mean liquid is
!
!  (32) dl/dt= -CA * qa * [(ql/qa)^(7/3)] * [(N*dens_h2o)^(-1/3)] 
!
!               * H(rad_liq - rthresh)
!
!       where N is the number of cloud droplets per cubic metre,
!       rthresh is a particle radius threshold needed to for autoconv-
!       ersion to occur, H is the Heaviside function, and CA is
!       a constant which is:
!
!  (33) CA =  0.104 * grav * Ec * (airdens)^(4/3) / mu 
!        
!       where grav is gravitational acceleration, Ec is the collection
!       efficiency, airdens is the density of air, and mu is the 
!       dynamic viscosity of air.   This constant is evaluated 
!       ignoring the temperature dependence of mu and with a fixed 
!       airdens of 1 kg/m3.
!
!       With   Ec = 0.55        (standard choice - see references)
!            grav = 9.81        m/(s*s)
!         airdens = 1.00        kg air/(m*m*m)
!              mu = 1.717  E-05 kg condensate/(m*s) 
!
!              CA = 32681. [(kg air)^(4/3)]/kg liq/m/s
!
!
!       For numerical treatment we write (32) as:
!
!       dl/dt = - D_aut * l 
!
!       and if we do so:
!
!  (34) D_aut   =   CA * [(N*dens_h2o)^(-1/3)] * [(ql/qa)^(4/3)] * &
!                   H(r-rthresh)
!
!       In the work below, D_aut is temporarily stored in the variable
!       tmp1 before being added to D1_dt.  D1_dt represents the 
!       conversion of liquid to rain.
!
!       Following Rotstayn, autoconversion is limited to the amount that
!       would reduce the local liquid cloud condensate to the critical
!       value at which autoconversion begins. This limiter is likely to
!       be invoked frequently and is computed from
!
!  (35) D_dt = log( (rad_liq/rthresh)**3. )
!
!       This limiter is stored in tmp2.
!
!
!
!       -------------------------------------------------
!
!       Khairoutdinov and Kogan (2000) Autoconversion
!
!       Reference: Khairoutdinov, M. and Y. Kogan, 2000: A new cloud 
!                  physics parameterization in a large-eddy simulation
!                  model of marine stratocumulus. Mon. Wea. Rev., 128,
!                  229-243.
!
!       If the namelist parameter use_kk_auto = true, then the 
!       Khairoutdinov and Kogan (KK) autoconversion parameterization
!       is used in place of the Manton and Cotton formula described
!       above.
!
!       In SI units this formula is:
!
!  (32A) dl/dt= -CA * qa * [(ql/qa)^(2.47)] * [(N)^(-1.79)] 
!
!
!       where N is the number of cloud droplets per cubic metre
!       and CA is a constant which is:
!
!  (33A) CA =  7.4188E+13 (kg condensate/kg air)**(-1.47)
!                         (# drops/meters**3)**(1.79)
!                         seconds**(-1) 
!        
!       For numerical treatment we write (32A) as:
!
!       dl/dt = - D_aut * l 
!
!       and if we do so:
!
!  (34A) D_aut   =   CA * [(N)^(-1.79)] * [(ql/qa)^(1.47)]
!

        if (do_liq_num) then
          if (use_kk_auto) then
!*************************yim's version based on Khairoutdinov and Kogan (2000)
!The second place N is used
          tmp1 = dtcloud * 1350. *  &
        (1.e-6*max(qn_mean*airdens(:,:,j),max(qa_mean,qmin)*N_min))**(-1.79)*  &
         (ql_mean)**(1.47)*max(qa_mean,qmin)**(0.32)
!         tmp1 = dtcloud * 1350. *  &
!                (1.e-6*max(qn_mean,N_min)*airdens(:,:,j))**(-1.79)*  &
!                (ql_mean)**(1.47)*max(qa_mean,qmin)**(0.32)
!**************************
          else
!yim fall back to M & C using qn
             !compute autoconversion sink as in (34)
             tmp1 = 32681. * dtcloud * ((max(qn_mean,qmin)*airdens(:,:,j)*dens_h2o)**(-1./3.))*       &
                       (ql_mean**(4./3.))/max(qa_mean,qmin)
  
             !compute limiter as in (35)
             tmp2 =max(3*log(max(rad_liq,qmin)/rthresh),0.)
  
             !limit autoconversion to the limiter
             tmp1 = min(tmp1,tmp2)
          endif
       else
        if ( use_kk_auto ) then
             tmp1 = 0.
             !compute autoconversion sink as in (34A)
             where (ql_mean.gt.qmin)
                  tmp1 = 7.4188E+13 * dtcloud *  (N**(-1.79))*         &
                    ((ql_mean/max(qa_mean,qmin))**(1.47))
             endwhere
        else
             !compute autoconversion sink as in (34)
             tmp1 = 32681. * dtcloud * ((N*dens_h2o)**(-1./3.))*       &
                    ((ql_mean/max(qa_mean,qmin))**(4./3.))
        
             !compute limiter as in (35)
             tmp2 =max(3*log(max(rad_liq,qmin)/rthresh),0.)

             !limit autoconversion to the limiter
             tmp1 = min(tmp1,tmp2)

        endif
        endif


        !add autoconversion to D1_dt
        D1_dt = D1_dt + tmp1

        !auto conversion will change a_rain_cld upto area of cloud
        where (tmp1 .gt. Dmin) a_rain_cld = qa_mean

!       if (max(diag_id%qldt_auto,diag_id%ql_auto_col) > 0) qldt_auto(:,:,j) = tmp1        
        if (max(diag_id%qldt_auto,diag_id%ql_auto_col) > 0) diag_4d(:,:,j,diag_pt%qldt_auto) = tmp1        

        if ( diag_id%aauto > 0 ) then
!            where ( rad_liq .gt. rthresh ) areaautocv(:,:,j) = qa_mean       
             where ( rad_liq .gt. rthresh ) diag_4d(:,:,j,diag_pt%aauto) = qa_mean       
        end if
        

!       Bergeron-Findeisan Process 
!
!       where ice and liquid coexist, the differential saturation
!       vapor pressure between liquid and ice phases encourages
!       the growth of ice crystals at the expense of liquid droplets.
!
!       Rotstayn (2000) derive an equation for the growth of ice by
!       starting with the vapor deposition equation for an ice crystal
!       and write it in terms of ice specific humidity as:
!
!                 {(Ni/airdens)**2/3}*7.8
!  (36) dqi/dt =  ------------------------  X [(esl-esi)/esi] X
!                 [rhoice**1/3]* A_plus_B
!
!                 ((max(qi,Mio*Ni/airdens))**(1/3))*
!
!       Here Ni is the ice crystal number which is taken from the 
!       parameterization of Meyers et al. :
!
!  (37) Ni = 1000 * exp( (12.96* [(esl-esi)/esi]) - 0.639 )
!
!       The use of the maximum operator assumed that there is a 
!       background ice crystal always present on which deposition can 
!       occur.  Mio is an initial ice crystal mass taken to be 10-12kg.
!
!       Figure 9.3 of Rogers and Yau (1998) shows the nearly linear
!       variation of [(esl-esi)/esi] from 0. at 273.16K to 0.5 at 
!       233.16K.  Analytically this is parameterized as (tfreeze-T)/80.
!
!
!       Generalizing (36) to grid mean conditions and writing it in 
!       terms of a loss of cloud liquid yields:
!
!                  (1000*exp((12.96*(tfreeze-T)/80)-0.639)/airdens)**2/3
!  (38) dql/dt = - ----------------------------------------------------- 
!                           [rhoice**1/3]* A_plus_B * ql
!
!       *qa*7.8*((max(qi/qa,Mio*Ni/airdens))**(1/3))*[(tfreeze-T)/80]*ql
!
!       Note that the density of ice is set to 700 kg m3 the value 
!       appropriate for pristine ice crystals.  This value is 
!       necessarily different than the value of ice used in the riming 
!       and sublimation part of the code which is for larger particles 
!       which have much lower densities.
        
        if (do_dust_berg) then
        crystal=0.
        do k = 1,jdim
                do i = 1,idim
                if ( (T(i,k,j) .lt. tfreeze) .and. (ql_mean(i,k) .gt. qmin)      &
                                        .and. (qa_mean(i,k) .gt. qmin))         then
                                Si0=1+0.0125*(tfreeze-T(i,k,j))
                                call Jhete_dep(T(i,k,j),Si0,concen_dust_sub(i,k,j),crystal(i,k))
!                               if (diag_id%debug4_3d > 0) debug4(i,k,j) = 1.                                      
                                if (diag_id%dust_berg_flag > 0) diag_4d(i,k,j,diag_pt%dust_berg_flag) = 1.                                      
                endif
                end do
          end do

!        if (max(diag_id%qndt_cond,diag_id%qn_cond_col) > 0) qndt_cond(:,:,j) = crystal
         if (max(diag_id%qndt_cond,diag_id%qn_cond_col) > 0) diag_4d(:,:,j, diag_pt%qndt_cond) = crystal
         if (max(diag_id%qndt_evap,diag_id%qn_evap_col) > 0) then
!          qndt_evap(:,:,j) = 0.
           diag_4d(:,:,j,diag_pt%qndt_evap) = 0.
           where (T(:,:,j).lt.tfreeze .and. ql_mean.gt.qmin .and. qa_mean.gt.qmin)              
!           qndt_evap(:,:,j) = 1.e-3*exp((12.96*0.0125*(tfreeze-T(:,:,j)))-0.639)
            diag_4d(:,:,j,diag_pt%qndt_evap) = 1.e-3*exp((12.96*0.0125*(tfreeze-T(:,:,j)))-0.639)
           end where
         endif
 
       !do Bergeron process
       D2_dt = 0.0
       where (T(:,:,j) .lt. tfreeze .and. ql_mean .gt. qmin .and. qa_mean .gt. qmin)              
             D2_dt =  dtcloud * qa_mean * ((1.e6*crystal(:,:)/airdens(:,:,j))**(2./ &
                      3.))* 7.8* ((max(qi_mean/qa_mean,1.E-12*1.e6*   &
                      crystal(:,:)     &
                      /airdens(:,:,j)))**(1./3.))*0.0125*              &
                      (tfreeze-T(:,:,j))/((700.**(1./3.))*       &
                      A_plus_B(:,:,j)*ql_mean)
       end where

    else
        !do Bergeron process
        D2_dt = 0.0        
        where (T(:,:,j) .lt. tfreeze .and. ql_mean .gt. qmin .and. qa_mean .gt. qmin)           
             D2_dt =  dtcloud * qa_mean * ((cfact*1000.*exp((12.96*0.0125*   &
                      (tfreeze-T(:,:,j)))-0.639)/airdens(:,:,j))**(2./ &
                      3.))* 7.8* ((max(qi_mean/qa_mean,1.E-12*cfact*1000.*   &
                      exp((12.96*0.0125*(tfreeze-T(:,:,j)))-0.639)     &
                      /airdens(:,:,j)))**(1./3.))*0.0125*              &
                      (tfreeze-T(:,:,j))/((700.**(1./3.))*             &
                      A_plus_B(:,:,j)*ql_mean)
        end where

      endif

        sum_berg = D2_dt
!       if (max(diag_id%qldt_berg,diag_id%ql_berg_col) > 0) qldt_berg(:,:,j) = D2_dt
        if (max(diag_id%qldt_berg,diag_id%ql_berg_col) > 0) diag_4d(:,:,j,diag_pt%qldt_berg) = D2_dt
       
!       Accretion of cloud liquid by ice ('Riming')
!       
!       [The below follows Rotstayn]
!       Accretion of cloud liquid by ice ('Riming') is derived in
!       the same way as the accretion of cloud liquid by rain. That
!       is the continous-collection equation for the growth of an
!       ice crystal is integrated over all ice crystal sizes. This
!       calculation assumes all crystals fall at the mass weighted
!       fall speed, Vfall.  This yields the following equation after
!       accounting for the area of the interaction
!
!  (39) dql/dt = -  ( a_snow_cld / qa/a_snow_cld ) *
!                   ( ELI * lamda_f * snow_cld / 2/ rho_ice ) * ql 
!
!
!       Note that in the version with the snowmelt bug, riming was
!       prevented when temperatures were in excess of freezing.

        !add in accretion of cloud liquid by ice
        tmp1 = 0.0
        if (do_old_snowmelt) then
             where ((a_snow_cld.gt.qmin) .and. (ql_mean.gt.qmin) .and. &
                    (   qa_mean.gt.qmin) .and. (T(:,:,j) .lt. tfreeze) )            
                 tmp1 = dtcloud*0.5*ELI*lamda_f*snow_cld/qa_mean/rho_ice              
             end where
        else
             where ((a_snow_cld.gt.qmin) .and. (ql_mean.gt.qmin) .and. &
                    (   qa_mean.gt.qmin) )            
                 tmp1 = dtcloud*0.5*ELI*lamda_f*snow_cld/qa_mean/rho_ice              
             end where
        end if        
        
        D2_dt = D2_dt + tmp1

        sum_rime = tmp1
!       if (max(diag_id%qldt_rime,diag_id%ql_rime_col) > 0) qldt_rime(:,:,j) = tmp1
        if (max(diag_id%qldt_rime,diag_id%ql_rime_col) > 0) diag_4d(:,:,j,diag_pt%qldt_rime) = tmp1

!       Freezing of cloud liquid to cloud ice occurs when
!       the temperature is less than -40C. At these very cold temper-
!       atures it is assumed that homogenous freezing of cloud liquid
!       droplets will occur.   To accomplish this numerically in one 
!       time step:
!
!  (40) D*dtcloud =  ln( ql / qmin ).
!
!       With this form it is guaranteed that if this is the only
!       process acting that ql = qmin after one integration.
!
               
        !do homogeneous freezing
        where ( T(:,:,j).lt.tfreeze-40..and.(ql_mean.gt.qmin).and.     &
               (qa_mean.gt.qmin))
             D2_dt = log ( ql_mean / qmin )
        end where
        
          do k=1,jdim
           do i=1,idim
             if (T(i,k,j).lt.(tfreeze-40.).and.(ql_mean(i,k).gt.qmin)    &
               .and.(qa_mean(i,k).gt.qmin)) then
               sum_freeze(i,k) = D2_dt(i,k)
               sum_rime(i,k) = 0.
               sum_berg(i,k) = 0.
             endif
           enddo
          enddo

        if (max(diag_id%qldt_freez,diag_id%qldt_rime,diag_id%qldt_berg,diag_id%ql_freez_col,  &
                diag_id%ql_rime_col,diag_id%ql_berg_col) > 0) then
          do k=1,jdim
           do i=1,idim
             if (T(i,k,j).lt.(tfreeze-40.).and.(ql_mean(i,k).gt.qmin)    &
               .and.(qa_mean(i,k).gt.qmin)) then
               if (max(diag_id%qldt_freez,diag_id%ql_freez_col) > 0) diag_4d(i,k,j,diag_pt%qldt_freez) = D2_dt(i,k)
               if (max(diag_id%qldt_rime,diag_id%ql_rime_col) > 0) diag_4d(i,k,j,diag_pt%qldt_rime ) = 0.
               if (max(diag_id%qldt_berg,diag_id%ql_berg_col) > 0) diag_4d(i,k,j,diag_pt%qldt_berg ) = 0.     
             endif
           enddo
          enddo
        end if
  
! For BC aerosol in-cloud scavenging:
          do k=1,jdim
           do i=1,idim
             qldt_sum = sum_berg(i,k) + sum_freeze(i,k) + sum_rime(i,k)
             if (qldt_sum > 0.) f_snow_berg(i,k,j) = sum_berg(i,k)/qldt_sum
           enddo
          enddo

!       Analytic integration of ql equation
!
!
!       The next step is to analytically integrate the cloud liquid
!       condensate equation.  This follows the Tiedtke approach.
!
!       The qc equation is written in the form:
!
!  (41) dqc/dt    =   C   -  qc * D   
!
!       Note that over the physics time step, C and D are assumed to 
!       be constants.
!
!       Defining qc(t) = qc0 and qc(t+dtcloud) = qc1, the analytic
!       solution of the above equation is:
!
!  (42) qc1 = qceq -  (qceq - qc0) * exp (-D*dtcloud)
! 
!       where qceq is the equilibrium cloud condensate that is approached
!       with an time scale of 1/D,
!
!  (43) qceq  =   C / D 
!
!
!       To diagnose the magnitude of each of the right hand terms of
!       (41) integrated over the time step, define the average cloud
!       condensate in the interval t to t + dtcloud qcbar as:
!
!  (44) qcbar  = qceq - [ (qc1-qc0) / ( dtcloud * D ) ]
! 
!       from which the magnitudes of the C and D terms integrated
!       over the time step are:
!
!       C   and   -D * (qcbar)
!   
!
!       Additional notes on this analytic integration:
!
!       1.   Because of finite machine precision it is required that
!            D*dt is greater than a minimum value.  This minimum
!            alue occurs where 1. - exp(-D*dt) = 0. instead of 
!            D*dt.  This value will be machine dependent. See discussion
!            at top of code for Dmin.
!

        !C_dt is set to large-scale condensation. Sink of cloud liquid 
        !is set to the sum of D1 (liquid to rain component), and D2 
        !(liquid to ice component), D_eros (erosion), and large-scale 
        !evaporation (note use of ql mean).        

        do k=1,jdim
        do i=1,idim
        C_dts = max(dcond_ls(i,k),0.)
        D_dts = D1_dt(i,k) + D2_dt(i,k) + D_eros(i,k) +                                &
               (max(-1.*dcond_ls(i,k),0.)/max(ql_mean(i,k),qmin)) 
                             
        !do analytic integration      
        qc0s   = ql_upd(i,k)
        if ( D_dts.gt.Dmin ) then
             qceqs  = C_dts   /  D_dts
             qc1s   = qceqs - (qceqs - qc0s) * exp ( -1.* D_dts )
             qcbars = qceqs - ((qc1s - qc0s)/ D_dts)
        else
             qceqs  = qc0s + C_dts   
             qc1s   = qc0s + C_dts
             qcbars = qc0s + 0.5*C_dts
        endif

        !set total tendency term and update cloud
        !Note that the amount of SL calculated here is stored in tmp1.
        SL(i,k,j)  = SL(i,k,j) + qc1s - qc0s
        ql_upd(i,k)     = qc1s

        !compute the amount each term contributes to the change     
!rab        Dterm  = -D_dt *      qcbar

!       Apportion SL between various processes.  This is necessary to
!       account for how much the temperature changes due to various
!       phase changes.   For example:
!
!       liquid to ice   = (D2/D)*(-Dterm)
!
!       liquid to rain  = (D1/D)*(-Dterm) 
!
!       (no phase change but needed to know how much to increment 
!        rainflux)
!
!       vapor to liquid = - { ((-dcond_ls/ql_mean)+D_eros)/D}*(-Dterm) 
!                        where dcond_ls < 0 
!                                  
!                         but
!
!                        dcond_ls  -(D_eros/D)*(-Dterm)
!                        where dcond_ls > 0
!

        !initialize tmp2 to hold (-Dterm)/D
!rab        tmp2 = -Dterm/max(D_dt,Dmin)

   if (Nml%retain_cm3_bug) then
        tmp2(i,k) = D_dts*qcbars/max(D_dts,Dmin)
   else
        if ( D_dts.gt.Dmin ) then
        tmp2(i,k) = D_dts*qcbars/max(D_dts,Dmin)
        else
        tmp2(i,k) = 0.
        endif
   endif
        
        !do phase changes from large-scale processes and boundary
        !layer condensation/evaporation
 
        ST(i,k,j) = ST(i,k,j) + (hlv*max(dcond_ls(i,k),0.)/cp_air) -          &
             (hlv*(max(-1.*dcond_ls(i,k),0.) /max(ql_mean(i,k),qmin))*tmp2(i,k)/cp_air)
   
        SQ(i,k,j) = SQ(i,k,j) -      max(dcond_ls(i,k),0.)     +            &
                  (max(-1.*dcond_ls(i,k),0.) /max(ql_mean(i,k),qmin))*tmp2(i,k)
            
        !add in liquid to ice and cloud erosion to temperature tendency
        ST(i,k,j) = ST(i,k,j) + (hlf*D2_dt(i,k)-hlv*D_eros(i,k))*tmp2(i,k)/cp_air

        !cloud evaporation adds to water vapor
        SQ(i,k,j) = SQ(i,k,j) + D_eros(i,k)*tmp2(i,k)
             
        !add conversion of liquid to rain to the rainflux
        rain_cld(i,k) = rain_cld(i,k) +D1_dt(i,k)*tmp2(i,k)*deltpg(i,k)*inv_dtcloud
     
        !save liquid converted to ice into tmp3 and increment qi_mean
        tmp3(i,k)    = tmp2(i,k)*D2_dt(i,k)
        qi_mean(i,k) = qi_mean(i,k) + tmp3(i,k)
        enddo
        enddo
     
        
        if (do_liq_num) then
!******************************************************************

!       Analytic integration of qn equation
!
!       The qn equation is written in the form:
!
!  (m1) dqc/dt    =   (1 - qabar) * A * qc^   -  qc * D
!                 =   C - qc * D
!
!       where  qc^ is the large-scale qn calculated from the
!       activation parameterization.   Note that over the physics
!       time step, C and D are assumed to be constants.
!
!       Defining qc(t) = qc0 and qc(t+dtcloud) = qc1, the analytic
!       solution of the above equation is:
!
!  (m2) qc1 = qceq -  (qceq - qc0) * exp (-D*dtcloud)
! 
!       where qceq is the equilibrium cloud droplet number that is approached
!       with an time scale of 1/D,
!
!  (m3) qceq  = C / D
!
!
!       To diagnose the magnitude of each of the right hand terms of
!       (m1) integrated over the time step, define the average cloud
!       condensate in the interval t to t + dtcloud qcbar as:
!
!  (m4) qcbar  = qceq - [ (qc1-qc0) / ( dtcloud * D ) ]
! 
!       from which the magnitudes of the C and D terms integrated
!       over the time step are:
!
!       C and -D * (qcbar)
!

!Calculate C_dt
!       C_dt=max(tmp5,0.)*drop1*1.e6/airdens(:,:,j)
!For replying the review, substract autoconversion
!        D_dt = D1_dt + D2_dt + D_eros 

        !do analytic integration      
!       where ( (D_dt.gt.Dmin) ) 
!            qc0   = qn_upd
!            qceq  = C_dt  / max(D_dt, Dmin)
!            qc1   = qceq - (qceq - qc0) * exp ( -1.* D_dt )
!            qcbar = qceq - ((qc1 - qc0)/ max(D_dt, Dmin))
!       elsewhere
!            qc0   = qn_upd
!            qceq  = qc0 + C_dt   
!            qc1   = qc0 + C_dt
!            qcbar = qc0 + 0.5*C_dt
!       end where

!        if (max(diag_id%qndt_cond,diag_id%qn_cond_col) > 0) qndt_cond(:,:,j) = 0.
         if (max(diag_id%qndt_cond,diag_id%qn_cond_col) > 0) diag_4d(:,:,j,diag_pt%qndt_cond) = 0.

          do k=1,jdim
            do i=1,idim
!Calculate C_dt
              C_dts=max(delta_cf(i,k),0.)*drop1(i,k)*1.e6/airdens(i,k,j)
              D_dts =  num_mass_ratio1*D1_dt(i,k) + (num_mass_ratio2*D2_dt(i,k) + D_eros(i,k))
              qc0s = qn_upd(i,k)
              if (D_dts > Dmin) then
                qceqs = C_dts / D_dts
                qc1s  = qceqs - (qceqs - qc0s)* exp(-1.*D_dts)
                qcbars = qceqs - ((qc1s -qc0s)/D_dts)
              else
                qceqs  = qc0s + C_dts
                qc1s   = qc0s + C_dts
                qcbars = qc0s + 0.5*C_dts
              endif
        !set total tendency term and update cloud
        !Note that the amount of SN calculated here is stored in tmp1.
              SN(i,k,j)  = SN(i,k,j) + qc1s - qc0s
              qn_upd(i,k)     = qc1s


        !compute the amount each term contributes to the change 
!        where ( C_dt .gt. 0 )
!                Cterm  =  C_dt             
!                Dterm  =  D_dt *      qcbar 
!        elsewhere
!                Cterm  =  0.             
!                Dterm  =  D_dt *      qcbar 
!        end where

         if (max(diag_id%qndt_cond,diag_id%qn_cond_col) > 0) then
!         if ( C_dts.gt. 0 ) qndt_cond(i,k,j)  =  C_dts
          if ( C_dts.gt. 0 ) diag_4d(i,k,j,diag_pt%qndt_cond)  =  C_dts
         endif

!        if (max(diag_id%qndt_evap,diag_id%qn_evap_col) > 0) qndt_evap(i,k,j) = D_dts * qcbars !Dterm
         if (max(diag_id%qndt_evap,diag_id%qn_evap_col) > 0) diag_4d(i,k,j,diag_pt%qndt_evap) = D_dts * qcbars !Dterm
            end do
          end do

        endif  ! (do_liq_num)

!****************************************************************************


!
!       diagnostics for cloud liquid tendencies
!       

!       if (max(diag_id%qldt_cond,diag_id%ql_cond_col) > 0) qldt_cond(:,:,j)  = max(dcond_ls,0.) *inv_dtcloud
        if (max(diag_id%qldt_cond,diag_id%ql_cond_col) > 0) diag_4d(:,:,j,diag_pt%qldt_cond)  = max(dcond_ls,0.) *inv_dtcloud
!       if (max(diag_id%qldt_evap,diag_id%ql_evap_col) > 0) qldt_evap(:,:,j)  = (max(0.,-1.*dcond_ls )/max(ql_mean,   &
        if (max(diag_id%qldt_evap,diag_id%ql_evap_col) > 0) diag_4d(:,:,j,diag_pt%qldt_evap)  = (max(0.,-1.*dcond_ls )/max(ql_mean,   &
                                 qmin))           *tmp2*inv_dtcloud
!       if (max(diag_id%qldt_accr,diag_id%ql_accr_col) > 0) qldt_accr(:,:,j)  = qldt_accr (:,:,j)*tmp2*inv_dtcloud
        if (max(diag_id%qldt_accr,diag_id%ql_accr_col) > 0) diag_4d(:,:,j,diag_pt%qldt_accr)  = diag_4d(:,:,j,diag_pt%qldt_accr)*tmp2*inv_dtcloud
!       if (max(diag_id%qldt_auto,diag_id%ql_auto_col) > 0) qldt_auto(:,:,j)  = qldt_auto (:,:,j)*tmp2*inv_dtcloud
        if (max(diag_id%qldt_auto,diag_id%ql_auto_col) > 0) diag_4d(:,:,j,diag_pt%qldt_auto)  = diag_4d(:,:,j,diag_pt%qldt_auto )*tmp2*inv_dtcloud
!       if (max(diag_id%qldt_eros,diag_id%ql_eros_col) > 0) qldt_eros(:,:,j)  = D_eros           *tmp2*inv_dtcloud 
!       if (max(diag_id%qldt_berg,diag_id%ql_berg_col) > 0) qldt_berg(:,:,j)  = qldt_berg (:,:,j)*tmp2*inv_dtcloud
        if (max(diag_id%qldt_berg,diag_id%ql_berg_col) > 0) diag_4d(:,:,j,diag_pt%qldt_berg)  = diag_4d(:,:,j,diag_pt%qldt_berg )*tmp2*inv_dtcloud
!       if (max(diag_id%qldt_rime,diag_id%ql_rime_col) > 0) qldt_rime(:,:,j)  = qldt_rime (:,:,j)*tmp2*inv_dtcloud
        if (max(diag_id%qldt_rime,diag_id%ql_rime_col) > 0) diag_4d(:,:,j,diag_pt%qldt_rime)  = diag_4d(:,:,j,diag_pt%qldt_rime )*tmp2*inv_dtcloud
!       if (max(diag_id%qldt_freez,diag_id%ql_freez_col) > 0) qldt_freez(:,:,j) = qldt_freez(:,:,j)*tmp2*inv_dtcloud
        if (max(diag_id%qldt_freez,diag_id%ql_freez_col) > 0) diag_4d(:,:,j,diag_pt%qldt_freez) = diag_4d(:,:,j,diag_pt%qldt_freez)*tmp2*inv_dtcloud
!       if (max(diag_id%qndt_cond,diag_id%qn_cond_col) > 0) qndt_cond(:,:,j)  = qndt_cond(:,:,j)*inv_dtcloud 
        if (max(diag_id%qndt_cond,diag_id%qn_cond_col) > 0) diag_4d(:,:,j,diag_pt%qndt_cond)  = diag_4d(:,:,j,diag_pt%qndt_cond)*inv_dtcloud 
!       if (max(diag_id%qndt_evap,diag_id%qn_evap_col) > 0) qndt_evap(:,:,j)  = qndt_evap(:,:,j)*inv_dtcloud 
        if (max(diag_id%qndt_evap,diag_id%qn_evap_col) > 0) diag_4d(:,:,j,diag_pt%qndt_evap)  = diag_4d(:,:,j,diag_pt%qndt_evap)*inv_dtcloud 

!       if (max(diag_id%qldt_cond,diag_id%ql_cond_col) > 0) diag_4d(:,:,j,diag_pt%qldt_cond)  = max(dcond_ls,0.) *inv_dtcloud
!       if (max(diag_id%qldt_evap,diag_id%ql_evap_col) > 0) diag_4d(:,:,j,      diag_pt%qldt_evap)  = (max(0.,-1.*dcond_ls )/max(ql_mean,   &
!                                qmin))           *tmp2*inv_dtcloud
!       if (max(diag_id%qldt_accr,diag_id%ql_accr_col) > 0) diag_4d(:,:,j,      diag_pt%qldt_accr)  = qldt_accr (:,:,j)
!       if (max(diag_id%qldt_auto,diag_id%ql_auto_col) > 0) diag_4d(:,:,j,      diag_pt%qldt_auto)  = qldt_auto (:,:,j)
        if (max(diag_id%qldt_eros,diag_id%ql_eros_col) > 0) diag_4d(:,:,j,      diag_pt%qldt_eros)  = D_eros           *tmp2*inv_dtcloud 
!       if (max(diag_id%qldt_berg,diag_id%ql_berg_col) > 0) diag_4d(:,:,j,      diag_pt%qldt_berg)  = qldt_berg (:,:,j)
!       if (max(diag_id%qldt_rime,diag_id%ql_rime_col) > 0) diag_4d(:,:,j,      diag_pt%qldt_rime)  = qldt_rime (:,:,j)
!       if (max(diag_id%qldt_freez,diag_id%ql_freez_col) > 0) diag_4d(:,:,j,diag_pt%qldt_freez) = qldt_freez(:,:,j)
!       if (max(diag_id%qndt_cond,diag_id%qn_cond_col) > 0) diag_4d(:,:,j,      diag_pt%qndt_cond)  = qndt_cond(:,:,j) 
!       if (max(diag_id%qndt_evap,diag_id%qn_evap_col) > 0) diag_4d(:,:,j,      diag_pt%qndt_evap)  = qndt_evap(:,:,j) 
        


!-----                                                            -----! 
!                                                                      !
!                                END OF                                !
!                                                                      !
!                       LIQUID PHASE MICROPHYSICS                      !
!                                                                      !
!                                 AND                                  !
!                                                                      !
!                  ANALYTIC INTEGRATION OF QL EQUATION                 !
!                                                                      !
!                               SECTION                                !
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!                        ICE PHASE MICROPHYSICS                        !
!                                                                      !
!                                 AND                                  !
!                                                                      !
!                  ANALYTIC INTEGRATION OF QI EQUATION                 !
!                                                                      !
!                                                                      !
!                                                                      !
!       Ice settling
!
!       Ice settling is treated as in Heymsfield Donner 1990. 
!       The mass weighted fall speed is parameterized as in equation
!       #49.
!
!       In terms of the analytic integration with respect qi of Tiedtke,
!       the flux in from the top of the grid layer is equated to the
!       source term, and the flux out of the bottom of the layer is 
!       equated to the sink term:
!
!  (47) C_dt =  snow_cld * dtcloud * grav / deltp
!
!  (48) D_dt =  airdens * grav * Vfall * dtcloud / deltp
!
!       All ice crystals are assumed to fall with the same fall speed
!       which is given as in Heymsfield and Donner (1990) as:
!
!  (49) Vfall = 3.29 * ( (airdens*qi_mean/qa_mean)**0.16)
!
!       which is the formula in Heymsfield and Donner.  Note however
!       that because this is uncertain, sensitivity runs will be made
!       with different formulations. Note that when Vfall is computed 
!       the source incremented qi, qi_mean, is used.  This gives some 
!       implicitness to the scheme.

        !compute Vfall
        iwc = airdens(:,:,j)*qi_mean/max(qa_mean,qmin)
        where (iwc >= iwc_crit)
           Vfall = vfact*3.29 * iwc**0.16
        elsewhere
           Vfall = vfact*vfall_const2 * iwc**vfall_exp2
        end where

!       if (diag_id%vfall > 0) vfalldiag(:,:,j) = Vfall(:,:)*areaice(:,:,j)
!       if (diag_id%vfall > 0) vfalldiag(:,:,j) = Vfall(:,:)*diag_4d(:,:,j,diag_pt%aice)
!       if (diag_id%vfall > 0) diag_4d(:,:,j,diag_pt%vfall) = Vfall(:,:)*areaice(:,:,j)
        if (diag_id%vfall > 0) diag_4d(:,:,j,diag_pt%vfall) = Vfall(:,:)*diag_4d(:,:,j,diag_pt%aice)

        !add to ice source the settling ice flux from above
        !also note that tmp3 contains the source
        !of liquid converted to ice from above
        tmp3 = tmp3 + snow_cld*dtcloud/deltpg
        
        !Compute settling of ice. The result is multiplied by 
        !dtcloud/deltp to convert to units of D_dt.  
        !Note that if tracers are not advected then this is done
        !relative to the local vertical motion.
        if (tracer_advec) then
             tmp1 = 0.
        else
             tmp1 = omega(:,:,j)
        end if 
        
        where (qi_mean .gt. qmin .and. qa_mean .gt. qmin)
             D1_dt      = max(0.,((airdens(:,:,j)*Vfall)+(tmp1/grav))* &
                          dtcloud/deltpg )
             a_snow_cld = qa_mean     
        elsewhere
             D1_dt      = 0.
             snow_cld   = 0.
             a_snow_cld = 0.    
        end where 

        
!       Melting of in-cloud ice
!
!       Melting occurs where the temperature is greater than Tfreezing. 
!       This is an instaneous process such that no stratiform ice will
!       remain in the grid at the end of the timestep. 
!       No ice settles out of the grid box (i.e. D1 is set to zero) when
!       the amount of ice to be melted is less than that that would
!       bring the grid box to the freezing point.
!
!       The ice that melts becomes rain.  This is because if an ice
!       crystal of dimension ~100 microns and mass density of 100 kg/m2
!       melts it will become a droplet of SIZE 40 microns which is 
!       clearly a drizzle SIZE drop.  Ice crystals at temperatures near 
!       freezing are assumed to be this large, consistent with the 
!       assumption of particle SIZE temperature dependence.
!

        !compute grid mean change in cloud ice to cool the
        !grid box to 0C and store in temporary variable tmp1
        tmp1 = cp_air*(T(:,:,j)-tfreeze)/hlf

        ! If qi_mean > tmp1, then the amount of ice melted is
        ! limited to tmp1, otherwise melt all qi_mean.  The amount
        ! melted is stored in tmp2
        tmp2  = max(min(qi_mean,tmp1),0.)     
        D2_dt = max(0.,log(max(qi_mean,qmin)/max(qi_mean-tmp2,qmin)))
            
        !melting of ice creates area to a_rain_cld
        where (D2_dt .gt. Dmin) 
             a_rain_cld = qa_mean
        endwhere
                  
        !If all of the ice can melt, then don't permit any ice to fall
        !out of the grid box and set a_snow_cld to zero.
        where (qi_mean .lt. tmp1 .and. qi_mean .gt. qmin) 
             D1_dt = 0.
             snow_cld = 0.
             a_snow_cld = 0.
        end where
         
!       Analytic integration of qi equation
!
!       This repeats the procedure done for the ql equation.
!       See above notes for detail.

        !At this point C_dt already includes the source of cloud ice 
        !falling from above as well as liquid converted to ice. 
        !Therefore add in large_scale deposition.
        !
        !Compute D_dt which has contributions from D1_dt (ice settling)
        !and D2_dt (ice melting), D_eros (cloud erosion), and large-
        !scale sublimation (note use of qi mean).

      do k=1,jdim
       do i=1,idim
        C_dts = tmp3(i,k) + max(dcond_ls_ice(i,k),0.)
        D_dts =  D1_dt(i,k) + D2_dt(i,k) + D_eros(i,k) +                               &
                (max(-1.*dcond_ls_ice(i,k),0.)/max(qi_mean(i,k),qmin))
        
        !do analytic integration      
        qc0s   = qi_upd(i,k)
        if ( D_dts.gt.Dmin ) then
             qceqs  = C_dts / D_dts
             qc1s   = qceqs - (qceqs - qc0s) * exp ( -1.* D_dts )
             qcbars = qceqs - ((qc1s - qc0s)/D_dts)
        else
             qceqs  = qc0s + C_dts   
             qc1s   = qc0s + C_dts
             qcbars = qc0s + 0.5*C_dts
        endif

        !set total tendency term and update cloud
        !Note that the amount of SL calculated here is stored in tmp1.
        SI(i,k,j)  = SI(i,k,j) + qc1s - qc0s
        qi_upd(i,k)     = qc1s

        !compute the amount each term contributes to the change     
!rab        Dterm  = -D_dt *          qcbar 
      
!       Apportion SI between various processes.  This is necessary to
!       account for how much the temperature and water vapor changes 
!       due to various phase changes.   For example:
!
!       ice settling = (D1/D)*(-Dterm)*deltp/grav/dtcloud
!                     
!       vapor to ice =
!           -{ ((-dcond_ls_ice/qi_mean)+D_eros)/ D }*(-Dterm) 
!           where dcond_ls_ice  < 0. 
!
!           but
!       
!           dcond_ls_ice -(D_eros/D)* (-Dterm)
!
!           where dcond_ls_ice > 0.
!       
!       melting of ice = (D2/D)*(-Dterm)*deltp/grav/dtcloud
!

        !initialize tmp2 to hold (-Dterm)/D
!rab        tmp2 = -Dterm/max(D_dt,Dmin)
   if (Nml%retain_cm3_bug) then
        tmp2s = D_dts*qcbars/max(D_dts,Dmin)
   else
        if ( D_dts.gt.Dmin ) then
        tmp2s = D_dts*qcbars/max(D_dts,Dmin)
        else
        tmp2s = 0.
        endif
   endif
        
        !do phase changes from large-scale processes 
        ST(i,k,j) = ST(i,k,j) +  hls*max(dcond_ls_ice(i,k),0.)/cp_air -    &
         hls*(max(-1.*dcond_ls_ice(i,k),0.)/max(qi_mean(i,k),qmin))*tmp2s/cp_air
       
        SQ(i,k,j) = SQ(i,k,j) -      max(dcond_ls_ice(i,k),0.)    +         &
             (max(-1.*dcond_ls_ice(i,k),0.)/max(qi_mean(i,k),qmin))*tmp2s
     
        !cloud erosion changes temperature and vapor
        ST(i,k,j) = ST(i,k,j) - hls*D_eros(i,k)* tmp2s/cp_air
        SQ(i,k,j) = SQ(i,k,j) +     D_eros(i,k)* tmp2s

        !add settling ice flux to snow_cld 
        snow_cld(i,k) = D1_dt(i,k)*tmp2s*deltpg(i,k)*inv_dtcloud
       
        !add melting of ice to temperature tendency
        ST(i,k,j) = ST(i,k,j) - hlf*D2_dt(i,k)*tmp2s/cp_air

        !add melting of ice to the rainflux
        rain_cld(i,k) = rain_cld(i,k) + D2_dt(i,k)*tmp2s*deltpg(i,k)*inv_dtcloud

!
!       diagnostics for cloud ice tendencies
!       
        
!       if (max(diag_id%qidt_dep,diag_id%qi_dep_col) > 0) &
!               qidt_dep (i,k,j) = max(dcond_ls_ice(i,k),0.)*inv_dtcloud
        if (max(diag_id%qidt_dep,diag_id%qi_dep_col) > 0) &
                diag_4d(i,k,j,diag_pt%qidt_dep ) = max(dcond_ls_ice(i,k),0.)*inv_dtcloud
!       if (max(diag_id%qidt_subl,diag_id%qi_subl_col) > 0) &
!               qidt_subl(i,k,j) = (max(0.,-1.*dcond_ls_ice(i,k))/ &
!                       max(qi_mean(i,k),qmin))*tmp2s*inv_dtcloud
        if (max(diag_id%qidt_subl,diag_id%qi_subl_col) > 0) &
                diag_4d(i,k,j,diag_pt%qidt_subl) = (max(0.,-1.*dcond_ls_ice(i,k))/ &
                        max(qi_mean(i,k),qmin))*tmp2s*inv_dtcloud
!       if (max(diag_id%qidt_melt,diag_id%qi_melt_col) > 0) &
!               qidt_melt(i,k,j) = D2_dt(i,k) *tmp2s*inv_dtcloud
        if (max(diag_id%qidt_melt,diag_id%qi_melt_col) > 0) &
                diag_4d(i,k,j,diag_pt%qidt_melt) = D2_dt(i,k) *tmp2s*inv_dtcloud
!       if (max(diag_id%qidt_eros,diag_id%qi_eros_col) > 0) &
!               qidt_eros(i,k,j) = D_eros(i,k)*tmp2s*inv_dtcloud       
        if (max(diag_id%qidt_eros,diag_id%qi_eros_col) > 0) &
                diag_4d(i,k,j,diag_pt%qidt_eros) = D_eros(i,k)*tmp2s*inv_dtcloud       

        if (max(diag_id%lsf_strat,diag_id%lcf_strat,diag_id%mfls_strat) > 0) then
            tmp1s = max(dcond_ls_ice(i,k),0.)*inv_dtcloud
!RSH BUGFIX 6/25/10:
!           if ( (qldt_cond(i,k,j) + tmp1s) .gt. 0. ) lsf_strat(:,:,j) = 1.
!           if ( (qldt_cond(i,k,j) + tmp1s) .gt. 0. ) diag_4d(:,:,j, diag_pt%lsf_strat) = 1.
!           if ( (qldt_cond(i,k,j) + tmp1s) .gt. 0. ) lsf_strat(i,k,j) = 1.
!           if ( (qldt_cond(i,k,j) + tmp1s) .gt. 0. ) diag_4d(i,k,j, diag_pt%lsf_strat) = 1.
!           if ( (diag_4d(i,k,j,diag_pt%qldt_cond) + tmp1s) .gt. 0. ) lsf_strat(i,k,j) = 1.
            if ( (diag_4d(i,k,j,diag_pt%qldt_cond) + tmp1s) .gt. 0. ) diag_4d(i,k,j, diag_pt%lsf_strat) = 1.
        endif
       enddo
      enddo

        
!-----                                                            -----! 
!                                                                      !
!                                END OF                                !
!                                                                      !
!                        ICE PHASE MICROPHYSICS                        !
!                                                                      !
!                                 AND                                  !
!                                                                      !
!                  ANALYTIC INTEGRATION OF QI EQUATION                 !
!                                                                      !
!                                 SECTION                              !
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!
!
!                       RAIN EVAPORATION
!                           
!
!       Rain evaporation is derived by integration of the growth 
!       equation of a droplet over the assumed Marshall-Palmer 
!       distribution of rain drops (equation #22).  This leads to the 
!       following formula:
!
!  (50) dqv/dt_local =  56788.636 * {rain_rate/dens_h2o}^(11/18) *(1-U)/
!
!                      ( SQRT(airdens)* A_plus_B)
!
!       Numerically this equation integrated by use of time-centered
!       values for qs and qv.   This leads to the solution:
!
!  (51) qv_clr(t+1)-qv_clr(t) = K3 *[qs(t)-qv_clr(t)]/
!                               {1.+0.5*K3*(1+gamma)}
!       where 
!
!       K3= 56788.636 * dtcloud * {rain_rate_local/dens_h2o}^(11/18) /
!           ( SQRT(airdens)* A_plus_B * qs)
!
!       and gamma is given by (3). Note that in (51), it is made 
!       explicit that it is the vapor concentration in the unsaturated 
!       part of the grid box that is used in the rain evaporation 
!       formula.
!
!       Now there are several limiters to this formula. First,
!       you cannot evaporate more than is available in a time step.
!       The amount available for evaporation locally is 
!       (rain_clr/a_rain_clr)*(grav*dtcloud/deltp).   Second, to
!       avoid supersaturating the box or exceeding the critical
!       relative humidity above which rain does not evaporate, 
!       the amount of evaporation is limited.
!
!       Finally rain evaporation occurs only if the relative humidity
!       in the unsaturated portion of the grid box, U_clr, is less
!       then a threshold, U_evap.   U_evap, will not necessarily be
!       one.   For example, stratiform precipitation in convective
!       regions rarely saturates subcloud air because of the presence
!       of downdrafts.  If the convection scheme does not have down-
!       drafts then it doesn't make sense to allow the sub-cloud layer
!       to saturate. U_clr may be solved from (8) as:
!
!  (52) U_clr = ( U - qa ) / (1. - qa)
!
!       Some variables are temporarily stored in tmp1.
!
!       Note that for pdf clouds the relative humidity in the clear part
!       of the grid box can be calculated exactly from the beta distr-
!       ibution. 

      do k=1,jdim
       do i=1,idim
        !compute U_clr
        if (.not. do_pdf_clouds) then 
             U_clr(i,k) =  (U(i,k)-qa_mean(i,k))/max((1.-qa_mean(i,k)),qmin)
        else
             U_clr(i,k) = qvg(i,k)/qs(i,k,j)
        end if
        
        !keep U_clr > 0. and U_clr < 1.
        U_clr(i,k) = min(max(U_clr(i,k),0.),1.)
        
        !compute K3
        tmp1s = 56788.636 * dtcloud * ((rain_clr(i,k)/max(a_rain_clr(i,k),qmin)/  &
             dens_h2o)**(11./18.))/SQRT(airdens(i,k,j))/A_plus_B(i,k,j)&
             /qs(i,k,j)

        !compute local change in vapor mixing ratio due to 
        !rain evaporation
        tmp1s = tmp1s*qs(i,k,j)*(1.-U_clr(i,k))/(1.+0.5*tmp1s*(1.+gamma(i,k,j)))

        !limit change in qv to the amount that would raise the relative
        !humidity to U_evap in the clear portion of the grid box
        tmp1s = min(tmp1s,((1.-qa_mean(i,k))/max(a_rain_clr(i,k),qmin))*qs(i,k,j)* &
               max(0.,U_evap-U_clr(i,k))/(1.+(U_evap*(1.-qa_mean(i,k))+qa_mean(i,k))* &
               gamma(i,k,j)) )
        
        !do limiter by amount available
        tmp1s= tmp1s*a_rain_clr(i,k)*deltpg(i,k)*inv_dtcloud
        tmp2s= max(min(rain_clr(i,k),tmp1s),0.)
    
        SQ(i,k,j) = SQ(i,k,j) +     tmp2s*dtcloud/deltpg(i,k)
        ST(i,k,j) = ST(i,k,j) - hlv*tmp2s*dtcloud/deltpg(i,k)/cp_air
        
        !if all of the rain evaporates set things to zero.    
        if (tmp1s.gt.rain_clr(i,k).and.a_rain_clr(i,k).gt.qmin) then
             rain_clr(i,k) = 0.
             a_rain_clr(i,k) = 0.
        else
             rain_clr(i,k) = rain_clr(i,k) - tmp2s   
        endif
        
!       if (max(diag_id%rain_evap,diag_id%rain_evap_col) > 0) rain_evap(i,k,j) = tmp2s/deltpg(i,k)
        if (max(diag_id%rain_evap,diag_id%rain_evap_col) > 0) diag_4d(i,k,j,diag_pt%rain_evap) = tmp2s/deltpg(i,k)
!rab       enddo 
!rab      enddo 

!-----------------------------------------------------------------------
!
!
!
!                              SNOW SUBLIMATION
!                           
!
!       Sublimation of cloud ice
!
!       [The following follows Rotstayn (1997)]
!       Given the assumptions of the Marshall-Palmer distribution of
!       ice crystals (18), the crystal growth equation as a function
!       of the humidity of the air and the diffusivity of water vapor
!       and thermal conductivity of air is integrated over all crystal
!       sizes.   This yields:
!
!  (53) dqi/dt_local = - a_snow_clr* K3 * (qs - qv_clr)
!
!       where the leading factor of a_snow_clr is the portion of the
!       grid box undergoing sublimation. K3 is given by
!
!  (54) K3 = (4/(pi*rho_air*qs*rho_ice*A_plus_B))*
!            ((snow_clr/a_snow_clr/3.29)**1.16 ) *
!           [ 0.65*lamda_f^2 + 
!             198.92227 * (airdens)^0.5 * 
!             ((snow_clr/a_snow_clr)**(1/14.5)) * lamda_f^(3/2) ]
!
!       Note that equation (53) is identical to equation (30) of 
!       Rotstayn.
!
!       Numerically this is integrated as in rain evaporation.


!rab      do k=1,jdim
!rab       do i=1,idim
        !compute K3
        tmp1s = dtcloud * (4./3.14159/rho_ice/airdens(i,k,j)/           &
               A_plus_B(i,k,j)/qs(i,k,j))*((snow_clr(i,k)/max(a_snow_clr(i,k),   &
               qmin)/3.29)**(1./1.16))*(0.65*lamda_f(i,k)*lamda_f(i,k) +         &
               198.92227*lamda_f(i,k)*SQRT(airdens(i,k,j)*lamda_f(i,k))*         &
               ( (snow_clr(i,k)/max(a_snow_clr(i,k),qmin))**(1./14.5) )  )

        !compute local change in vapor mixing ratio due to 
        !snow sublimation
        tmp1s = tmp1s*qs(i,k,j)*(1.-U_clr(i,k))/(1.+0.5*tmp1s*(1.+gamma(i,k,j)))

        !limit change in qv to the amount that would raise the relative
        !humidity to U_evap in the clear portion of the grid box
        tmp1s = min(tmp1s,((1.-qa_mean(i,k))/max(a_snow_clr(i,k),qmin))*qs(i,k,j)* &
               max(0.,U_evap-U_clr(i,k))/(1.+(U_evap*(1.-qa_mean(i,k))+qa_mean(i,k))* &
               gamma(i,k,j)) )
        
        !do limiter by amount available
        tmp1s= tmp1s*a_snow_clr(i,k)*deltpg(i,k)*inv_dtcloud
        tmp2s= max(min(snow_clr(i,k),tmp1s),0.)
    
        SQ(i,k,j) = SQ(i,k,j) +     tmp2s*dtcloud/deltpg(i,k)
        ST(i,k,j) = ST(i,k,j) - hls*tmp2s*dtcloud/deltpg(i,k)/cp_air
        
        !if all of the snow sublimates set things to zero.    
        if (tmp1s.gt.snow_clr(i,k).and.a_snow_clr(i,k).gt.qmin) then
             snow_clr(i,k) = 0.
             a_snow_clr(i,k) = 0.
        else
             snow_clr(i,k) = snow_clr(i,k) - tmp2s     
        endif
         
!       if (max(diag_id%snow_subl,diag_id%snow_subl_col) > 0) snow_subl(i,k,j) = tmp2s/deltpg(i,k)
        if (max(diag_id%qdt_snow_sublim,diag_id%q_snow_sublim_col) > 0) diag_4d(i,k,j,diag_pt%qdt_snow_sublim) = tmp2s/deltpg(i,k)
       enddo 
      enddo 

!-----------------------------------------------------------------------
!
!       Adjustment to try to prevent negative water vapor. Adjustment
!       will not remove more condensate than available. Cloud amount
!       is not adjusted. This is left for the remainder of the
!       desctruction code. (cjg)

!       if (max(diag_id%qadt_destr,diag_id%qa_destr_col) > 0) qadt_destr(:,:,j) = qadt_destr(:,:,j) + SA(:,:,j)*inv_dtcloud
        if (max(diag_id%qadt_destr,diag_id%qa_destr_col) > 0) diag_4d(:,:,j,diag_pt%qadt_destr) = diag_4d(:,:,j,diag_pt%qadt_destr) + SA(:,:,j)*inv_dtcloud
!       if (max(diag_id%qldt_destr,diag_id%ql_destr_col) > 0) qldt_destr(:,:,j) = qldt_destr(:,:,j) + SL(:,:,j)*inv_dtcloud
        if (max(diag_id%qldt_destr,diag_id%ql_destr_col) > 0) diag_4d(:,:,j,diag_pt%qldt_destr) = diag_4d(:,:,j,diag_pt%qldt_destr) + SL(:,:,j)*inv_dtcloud
!       if (max(diag_id%qidt_destr,diag_id%qi_destr_col) > 0) qidt_destr(:,:,j) = qidt_destr(:,:,j) + SI(:,:,j)*inv_dtcloud
        if (max(diag_id%qidt_destr,diag_id%qi_destr_col) > 0) diag_4d(:,:,j,diag_pt%qidt_destr) = diag_4d(:,:,j,diag_pt%qidt_destr) + SI(:,:,j)*inv_dtcloud
!       if (max(diag_id%qndt_destr,diag_id%qn_destr_col) > 0) qndt_destr(:,:,j) = qndt_destr(:,:,j) + SN(:,:,j)*inv_dtcloud
        if (max(diag_id%qndt_destr,diag_id%qn_destr_col) > 0) diag_4d(:,:,j,diag_pt%qndt_destr) = diag_4d(:,:,j,diag_pt%qndt_destr) + SN(:,:,j)*inv_dtcloud

        do k=1,jdim
         do i=1,idim
          tmp1s = qv(i,k,j) + SQ(i,k,j)
          tmp2s = 0.0
          tmp3s = 0.0
          if ( tmp1s.lt.0.0 ) then
           if (T(i,k,j).le.tfreeze-40.) then
            tmp2s = min( -tmp1s, ql_upd(i,k) )        ! liquid to evaporate
            tmp3s = min( -tmp1s-tmp2s, qi_upd(i,k) )   ! ice to sublimate
           else
            tmp3s = min( -tmp1s, qi_upd(i,k) )        ! ice to sublimate
            tmp2s = min( -tmp1s-tmp3s, ql_upd(i,k) )   ! liquid to evaporate
           end if
           ql_upd(i,k) = ql_upd(i,k) - tmp2s
           qi_upd(i,k) = qi_upd(i,k) - tmp3s
           SL(i,k,j) = SL(i,k,j) - tmp2s
           SI(i,k,j) = SI(i,k,j) - tmp3s
           SQ(i,k,j) = SQ(i,k,j) + tmp2s + tmp3s
           ST(i,k,j) = ST(i,k,j) - hlv*tmp2s/cp_air - hls*tmp3s/cp_air
          end if
         enddo
        enddo

!-----------------------------------------------------------------------
!       Cloud Destruction occurs where both ql and qi are .le. qmin, 
!       or if qa is .le. qmin. In this case, ql, qi, and qa are set to 
!       zero conserving moisture and energy.

        if (.not.do_liq_num) then
          do k=1,jdim
           do i=1,idim
            if ((ql_upd(i,k) .le. qmin .and. qi_upd(i,k) .le. qmin)               &
               .or. (qa_upd(i,k) .le. qmin)) then
             SL(i,k,j) = SL(i,k,j) - ql_upd(i,k)
             SI(i,k,j) = SI(i,k,j) - qi_upd(i,k)
             SQ(i,k,j) = SQ(i,k,j) + ql_upd(i,k) + qi_upd(i,k)
             ST(i,k,j) = ST(i,k,j) - (hlv*ql_upd(i,k) + hls*qi_upd(i,k))/cp_air
             SA(i,k,j) = SA(i,k,j) - qa_upd(i,k)
             ql_upd(i,k) = 0.0
             qi_upd(i,k) = 0.0
             qa_upd(i,k) = 0.0
            endif
           enddo
          enddo
        else
          do k=1,jdim
           do i=1,idim
            if ((ql_upd(i,k) .le. qmin .and. qi_upd(i,k) .le. qmin)               &
               .or. (qa_upd(i,k) .le. qmin)) then
             SL(i,k,j) = SL(i,k,j) - ql_upd(i,k)
             SI(i,k,j) = SI(i,k,j) - qi_upd(i,k)
             SQ(i,k,j) = SQ(i,k,j) + ql_upd(i,k) + qi_upd(i,k)
             ST(i,k,j) = ST(i,k,j) - (hlv*ql_upd(i,k) + hls*qi_upd(i,k))/cp_air
             SA(i,k,j) = SA(i,k,j) - qa_upd(i,k)
             SN(i,k,j) = SN(i,k,j) - qn_upd(i,k)
             ql_upd(i,k) = 0.0
             qi_upd(i,k) = 0.0
             qa_upd(i,k) = 0.0
             qn_upd(i,k) = 0.0
            endif
           enddo
          enddo
        endif  

!       if (max(diag_id%qadt_destr,diag_id%qa_destr_col) > 0) qadt_destr(:,:,j) = qadt_destr(:,:,j) - SA(:,:,j)*inv_dtcloud
        if (max(diag_id%qadt_destr,diag_id%qa_destr_col) > 0) diag_4d(:,:,j,diag_pt%qadt_destr) = diag_4d(:,:,j,diag_pt%qadt_destr) - SA(:,:,j)*inv_dtcloud
!       if (max(diag_id%qldt_destr,diag_id%ql_destr_col) > 0) qldt_destr(:,:,j) = qldt_destr(:,:,j) - SL(:,:,j)*inv_dtcloud
        if (max(diag_id%qldt_destr,diag_id%ql_destr_col) > 0) diag_4d(:,:,j,diag_pt%qldt_destr) = diag_4d(:,:,j,diag_pt%qldt_destr) - SL(:,:,j)*inv_dtcloud
!       if (max(diag_id%qidt_destr,diag_id%qi_destr_col) > 0) qidt_destr(:,:,j) = qidt_destr(:,:,j) - SI(:,:,j)*inv_dtcloud
        if (max(diag_id%qidt_destr,diag_id%qi_destr_col) > 0) diag_4d(:,:,j,diag_pt%qidt_destr) = diag_4d(:,:,j,diag_pt%qidt_destr) - SI(:,:,j)*inv_dtcloud
!       if (max(diag_id%qndt_destr,diag_id%qn_destr_col) > 0) qndt_destr(:,:,j) = qndt_destr(:,:,j) - SN(:,:,j)*inv_dtcloud
        if (max(diag_id%qndt_destr,diag_id%qn_destr_col) > 0) diag_4d(:,:,j,diag_pt%qndt_destr) = diag_4d(:,:,j,diag_pt%qndt_destr) - SN(:,:,j)*inv_dtcloud

!-----------------------------------------------------------------------
!
!       Adjustment.  Due to numerical errors in detrainment or advection
!       sometimes the current state of the grid box may be super-
!       saturated. Under the assumption that the temperature is constant
!       in the grid box and that q <= qs, the excess vapor is condensed. 
!       
!       What happens to the condensed water vapor is determined by the
!       namelist parameter super_choice.
!
!       If super_choice = .false. (default), then the condensed water is
!       is added to the precipitation fluxes.  If super_choice = .true.,
!       then the condensed water is added to the cloud condensate field.
!       Note that in this case the cloud fraction is raised to one.
!
!       The phase partitioning depends on super_choice; if super_choice
!       is false then at T < -20C, snow is produced.  If super_choice
!       is true, then at T < -40C, ice is produced.  The latter choice 
!       is consistent with that done in the large-scale condensation
!       section above.        
!
!       If pdf clouds are operating then this section is bypassed - 
!       as statistical clouds should remove supersaturation according
!       to the beta distribution used.
             
        if (.not.do_pdf_clouds) then
               
!       Old code

!       !estimate current qs
!       tmp2 = qs(:,:,j)+dqsdT(:,:,j)*ST(:,:,j)
!
!       !compute excess over saturation
!       tmp1 = max(0.,qv(:,:,j)+SQ(:,:,j)-tmp2)/(1.+gamma(:,:,j))
 
!       New more accurate version

!RSH 9/18/09: should be within the super_choice if block:
        ! updated temperature in tmp1
        tmp1 = T(:,:,j) + ST(:,:,j)

        ! updated qs in tmp2, updated gamma in tmp3, updated dqsdT in tmp5
        call compute_qs(tmp1, pfull(:,:,j), tmp2, dqsdT=tmp5)
 
        tmp3 = tmp5 *(min(1.,max(0.,0.05*(tmp1-tfreeze+20.)))*hlv +     &
                      min(1.,max(0.,0.05*(tfreeze -tmp1   )))*hls)/cp_air


        !compute excess over saturation
        tmp1 = max(0.,qv(:,:,j)+SQ(:,:,j)-tmp2)/(1.+tmp3)

!rab - save off tmp1 for diagnostic ice/liq adjustment fields in liq_adj array
!       if (max(diag_id%ice_adj,diag_id%ice_adj_col,diag_id%liq_adj,diag_id%liq_adj_col) .gt. 0) &
!              liq_adj(:,:,j) = tmp1*inv_dtcloud
        if (max(diag_id%ice_adj,diag_id%ice_adj_col,diag_id%liq_adj,diag_id%liq_adj_col) .gt. 0) &
               diag_4d(:,:,j,diag_pt%liq_adj) = tmp1*inv_dtcloud

        !change vapor content
        SQ(:,:,j)=SQ(:,:,j)-tmp1

        if (super_choice) then
        
             ! Put supersaturation into cloud

             !cloud fraction source diagnostic
             if (max(diag_id%qadt_super,diag_id%qa_super_col) > 0) then
               where (tmp1 .gt. 0.)
!                qadt_super(:,:,j)  = (1.-qa_upd) * inv_dtcloud
                 diag_4d(:,:,j,diag_pt%qadt_super)  = (1.-qa_upd) * inv_dtcloud
               endwhere
             endif

!yim 11/7/07
             if (do_liq_num) then
               do k=1,jdim
                 do i=1,idim
                   if (T(i,k,j) > tfreeze - 40. .and. &
                        tmp1(i,k) > 0.0) then
                     qn_upd(i,k) = qn_upd(i,k) + drop1(i,k)*1.0e6/  &
                                    airdens(i,k,j)*(1. - qa_upd(i,k))
                     SN(i,k,j) = SN(i,k,j) + drop1(i,k)*1.e6/  &
                                   airdens(i,k,j)*(1.-qa_upd(i,k))
                     if (max(diag_id%qndt_super,diag_id%qn_super_col) > 0)  &
!                       qndt_super(i,k,j) = qndt_super(i,k,j) +   &
                        diag_4d(i,k,j,diag_pt%qndt_super) = diag_4d(i,k,j,diag_pt%qndt_super) +   &
                                                drop1(i,k)*1.e6 / &
                            airdens(i,k,j)*(1.-qa_upd(i,k))*inv_dtcloud
                     endif
                  end do
               end do
!                    if (max(diag_id%qndt_super,diag_id%qn_super_col) > 0)  &
!               diag_4d(:,:,j,diag_pt%qndt_super) = qndt_super(:,:,j)
             endif

             !add in excess to cloud condensate, change cloud area and 
             !increment temperature
             do k=1,jdim
              do i=1,idim
               if(tmp1(i,k).gt.0) then
                if (T(i,k,j) .le. tfreeze-40.)then
                  qi_upd(i,k)   = qi_upd(i,k) + tmp1(i,k)
                  SI(i,k,j) = SI(i,k,j) + tmp1(i,k)
                  ST(i,k,j) = ST(i,k,j) + hls*tmp1(i,k)/cp_air
                else   ! where (T(i,k,j) .gt. tfreeze-40.)
                  ql_upd(i,k)   = ql_upd(i,k) + tmp1(i,k)
                  SL(i,k,j) = SL(i,k,j) + tmp1(i,k)
                  ST(i,k,j) = ST(i,k,j) + hlv*tmp1(i,k)/cp_air        
                endif
                if (limit_conv_cloud_frac) then
                  tmp2s = ahuco(i,k,j)
                else
                  tmp2s = 0.
                endif
                SA(i,k,j) = SA(i,k,j) + (1.-qa_upd(i,k)-tmp2s)
                qa_upd(i,k) = 1. - tmp2s
               endif
              enddo
             enddo

        else

             !Put supersaturation into precip

             !add in excess to precipitation fluxes, change their area 
             !and increment temperature
             do k=1,jdim
              do i=1,idim
               if(tmp1(i,k).gt.0) then
                if (T(i,k,j) .le. tfreeze-20.) then
                  snow_cld(i,k) = snow_cld(i,k) + qa_mean(i,k) *tmp1(i,k)*deltpg(i,k)*inv_dtcloud
                  snow_clr(i,k) = snow_clr(i,k) + (1.-qa_mean(i,k))*tmp1(i,k)*deltpg(i,k)*      &
                                                             inv_dtcloud
                  a_snow_cld(i,k) = qa_mean(i,k)
                  a_snow_clr(i,k) = 1.-qa_mean(i,k)
                  ST(i,k,j)  = ST(i,k,j) + hls*tmp1(i,k)/cp_air
                else   ! where (T(i,k,j) .gt. tfreeze-20.)
                  rain_cld(i,k) = rain_cld(i,k) + qa_mean(i,k) *tmp1(i,k)*deltpg(i,k)*inv_dtcloud
                  rain_clr(i,k) = rain_clr(i,k) + (1.-qa_mean(i,k))*tmp1(i,k)*deltpg(i,k)*      &
                                                             inv_dtcloud
                  a_rain_cld(i,k) = qa_mean(i,k)
                  a_rain_clr(i,k) = 1.-qa_mean(i,k)
                  ST(i,k,j)  = ST(i,k,j) + hlv*tmp1(i,k)/cp_air
                endif
               endif
              enddo
             enddo

        end if !super choice
        
        end if !for do_pdf_clouds
                      
!-----------------------------------------------------------------------
!       Final clean up to remove numerical noise
!

!       if (max(diag_id%qadt_destr,diag_id%qa_destr_col) > 0) qadt_destr(:,:,j) = qadt_destr(:,:,j) + SA(:,:,j)*inv_dtcloud
        if (max(diag_id%qadt_destr,diag_id%qa_destr_col) > 0) diag_4d(:,:,j,diag_pt%qadt_destr) = diag_4d(:,:,j,diag_pt%qadt_destr) + SA(:,:,j)*inv_dtcloud
!       if (max(diag_id%qldt_destr,diag_id%ql_destr_col) > 0) qldt_destr(:,:,j) = qldt_destr(:,:,j) + SL(:,:,j)*inv_dtcloud
        if (max(diag_id%qldt_destr,diag_id%ql_destr_col) > 0) diag_4d(:,:,j,diag_pt%qldt_destr) = diag_4d(:,:,j,diag_pt%qldt_destr) + SL(:,:,j)*inv_dtcloud
!       if (max(diag_id%qidt_destr,diag_id%qi_destr_col) > 0) qidt_destr(:,:,j) = qidt_destr(:,:,j) + SI(:,:,j)*inv_dtcloud
        if (max(diag_id%qidt_destr,diag_id%qi_destr_col) > 0) diag_4d(:,:,j,diag_pt%qidt_destr) = diag_4d(:,:,j,diag_pt%qidt_destr) + SI(:,:,j)*inv_dtcloud
!       if (max(diag_id%qndt_destr,diag_id%qn_destr_col) > 0) qndt_destr(:,:,j) = qndt_destr(:,:,j) + SN(:,:,j)*inv_dtcloud
        if (max(diag_id%qndt_destr,diag_id%qn_destr_col) > 0) diag_4d(:,:,j,diag_pt%qndt_destr) = diag_4d(:,:,j,diag_pt%qndt_destr) + SN(:,:,j)*inv_dtcloud

        do k=1,jdim
         do i=1,idim
          ql_upd(i,k) = ql(i,k,j) + SL(i,k,j)
          if ( abs(ql_upd(i,k)) .le. qmin  &
                .and. qv(i,k,j)+SQ(i,k,j)+ql_upd(i,k) > 0.0 ) then
            SL(i,k,j) = -ql(i,k,j)
            SQ(i,k,j) = SQ(i,k,j) + ql_upd(i,k)
            ST(i,k,j) = ST(i,k,j) - hlv*ql_upd(i,k)/cp_air
            ql_upd(i,k) = 0.0
          endif

          qi_upd(i,k) = qi(i,k,j) + SI(i,k,j)
          if ( abs(qi_upd(i,k)) .le. qmin  &
                .and. qv(i,k,j)+SQ(i,k,j)+qi_upd(i,k) > 0.0 ) then
            SI(i,k,j) = -qi(i,k,j)
            SQ(i,k,j) = SQ(i,k,j) + qi_upd(i,k)
            ST(i,k,j) = ST(i,k,j) - hls*qi_upd(i,k)/cp_air
            qi_upd(i,k) = 0.0
          endif

          qa_upd(i,k) = qa(i,k,j) + SA(i,k,j)
          if ( abs(qa_upd(i,k)) .le. qmin ) then
            SA(i,k,j) = -qa(i,k,j)
            qa_upd(i,k) = 0.0
          endif
         enddo
        enddo

!       if (max(diag_id%qadt_destr,diag_id%qa_destr_col) > 0) qadt_destr(:,:,j) = qadt_destr(:,:,j) - SA(:,:,j)*inv_dtcloud
        if (max(diag_id%qadt_destr,diag_id%qa_destr_col) > 0) diag_4d(:,:,j,diag_pt%qadt_destr) = diag_4d(:,:,j,diag_pt%qadt_destr) - SA(:,:,j)*inv_dtcloud
!       if (max(diag_id%qldt_destr,diag_id%ql_destr_col) > 0) qldt_destr(:,:,j) = qldt_destr(:,:,j) - SL(:,:,j)*inv_dtcloud
        if (max(diag_id%qldt_destr,diag_id%ql_destr_col) > 0) diag_4d(:,:,j,diag_pt%qldt_destr) = diag_4d(:,:,j,diag_pt%qldt_destr) - SL(:,:,j)*inv_dtcloud
!       if (max(diag_id%qidt_destr,diag_id%qi_destr_col) > 0) qidt_destr(:,:,j) = qidt_destr(:,:,j) - SI(:,:,j)*inv_dtcloud
        if (max(diag_id%qidt_destr,diag_id%qi_destr_col) > 0) diag_4d(:,:,j,diag_pt%qidt_destr) = diag_4d(:,:,j,diag_pt%qidt_destr) - SI(:,:,j)*inv_dtcloud
!       if (max(diag_id%qndt_destr,diag_id%qn_destr_col) > 0) qndt_destr(:,:,j) = qndt_destr(:,:,j) - SN(:,:,j)*inv_dtcloud
        if (max(diag_id%qndt_destr,diag_id%qn_destr_col) > 0) diag_4d(:,:,j,diag_pt%qndt_destr) = diag_4d(:,:,j,diag_pt%qndt_destr) - SN(:,:,j)*inv_dtcloud

!-----------------------------------------------------------------------
!
!       Save qa_mean of current level into qa_mean_lst.   This is used
!       in transferring rain and snow fluxes between levels.

        qa_mean_lst = qa_mean
        
!-----------------------------------------------------------------------
!
!       add the ice falling out from cloud to qidt_fall
        
!       if (max(diag_id%qidt_fall,diag_id%qi_fall_col) > 0) qidt_fall(:,:,j) = qidt_fall(:,:,j) +    & 
!                                        (snow_cld/deltpg)
        if (max(diag_id%qidt_fall,diag_id%qi_fall_col) > 0) diag_4d(:,:,j,diag_pt%qidt_fall) = diag_4d(:,:,j,diag_pt%qidt_fall) +    & 
                                         (snow_cld/deltpg)
        
!-----------------------------------------------------------------------
!
!       save profiles of rain and snow 
!
        rain3d(:,:,j+1) = rain_clr(:,:) + rain_cld(:,:)
        snow3d(:,:,j+1) = snow_clr(:,:) + snow_cld(:,:)
        snowclr3d(:,:,j+1) = snow_clr(:,:)

!-----------------------------------------------------------------------
!
!       Save rain and snow diagnostics

!       if (diag_id%rain_clr   > 0) rain_clr_diag(:,:,j+1)     = rain_clr
!       if (diag_id%rain_cld   > 0) rain_cld_diag(:,:,j+1)     = rain_cld
!       if (diag_id%a_rain_clr > 0) a_rain_clr_diag(:,:,j+1)   = a_rain_clr
!       if (diag_id%a_rain_cld > 0) a_rain_cld_diag(:,:,j+1)   = a_rain_cld
!       if (diag_id%snow_clr   > 0) snow_clr_diag(:,:,j+1)     = snow_clr
!       if (diag_id%snow_cld   > 0) snow_cld_diag(:,:,j+1)     = snow_cld
!       if (diag_id%a_snow_clr > 0) a_snow_clr_diag(:,:,j+1)   = a_snow_clr
!       if (diag_id%a_snow_cld > 0) a_snow_cld_diag(:,:,j+1)   = a_snow_cld
!       if (diag_id%a_precip_clr > 0) a_precip_clr_diag(:,:,j+1) = max(a_rain_clr,a_snow_clr)
!       if (diag_id%a_precip_cld > 0) a_precip_cld_diag(:,:,j+1) = max(a_rain_cld,a_snow_cld)
        if (diag_id%rain_clr   > 0) diag_4d_kp1(:,:,j+1,diag_pt%rain_clr)     = rain_clr
        if (diag_id%rain_cld   > 0) diag_4d_kp1(:,:,j+1,diag_pt%rain_cld)     = rain_cld
        if (diag_id%a_rain_clr > 0) diag_4d_kp1(:,:,j+1,diag_pt%a_rain_clr)   = a_rain_clr
        if (diag_id%a_rain_cld > 0) diag_4d_kp1(:,:,j+1,diag_pt%a_rain_cld)   = a_rain_cld
        if (diag_id%snow_clr   > 0) diag_4d_kp1(:,:,j+1,diag_pt%snow_clr)     = snow_clr
        if (diag_id%snow_cld   > 0) diag_4d_kp1(:,:,j+1,diag_pt%snow_cld)     = snow_cld
        if (diag_id%a_snow_clr > 0) diag_4d_kp1(:,:,j+1,diag_pt%a_snow_clr)   = a_snow_clr
        if (diag_id%a_snow_cld > 0) diag_4d_kp1(:,:,j+1,diag_pt%a_snow_cld)   = a_snow_cld
        if (diag_id%a_precip_clr > 0) diag_4d_kp1(:,:,j+1,diag_pt%a_precip_clr) = max(a_rain_clr,a_snow_clr)
        if (diag_id%a_precip_cld > 0) diag_4d_kp1(:,:,j+1,diag_pt%a_precip_cld) = max(a_rain_cld,a_snow_cld)

!-----------------------------------------------------------------------
!
!
!       Put rain and ice fluxes into surfrain and surfsnow if the
!       grid point is at the bottom of a column.   If MASK is not
!       present then this code is executed only if j .eq. kdim.
!       IF MASK is present some grid points may be beneath ground. 
!       If a given grid point is at the bottom of the column then
!       the surface values of rain and snow must be created.
!       Also if the MASK is present then the code forces all tenden-
!       cies below ground to be zero. Note that MASK = 1. equals above
!       ground point, MASK = 0. equals below ground point.

        if (present(MASK)) then

             !zero out all tendencies below ground
             ST(:,:,j)=MASK(:,:,j)*ST(:,:,j)
             SQ(:,:,j)=MASK(:,:,j)*SQ(:,:,j)
             SL(:,:,j)=MASK(:,:,j)*SL(:,:,j)
             SI(:,:,j)=MASK(:,:,j)*SI(:,:,j)
             SA(:,:,j)=MASK(:,:,j)*SA(:,:,j)
             if (do_liq_num) then
               SN(:,:,j)=MASK(:,:,j)*SN(:,:,j)
             endif
 
             if (j .lt. kdim) then
                  
                  !bottom of true points in columns which contain some
                  !dummy points
                  where(MASK(:,:,j) .eq. 1. .and. MASK(:,:,j+1) .eq. 0.)
                       surfrain = dtcloud*(rain_clr+rain_cld)
                       surfsnow = dtcloud*(snow_clr+snow_cld)
                       rain_clr = 0.
                       rain_cld = 0.
                       snow_clr = 0.
                       snow_cld = 0.
                       a_rain_clr = 0.
                       a_rain_cld = 0.
                       a_snow_clr = 0.
                       a_snow_cld = 0.
                  end where

             else

                  !bottom of column for those columns which contain no
                  !dummy points
                  where(MASK(:,:,j) .eq. 1.)
                       surfrain = dtcloud*(rain_clr+rain_cld)
                       surfsnow = dtcloud*(snow_clr+snow_cld)
                       rain_clr = 0.
                       rain_cld = 0.
                       snow_clr = 0.
                       snow_cld = 0.
                       a_rain_clr = 0.
                       a_rain_cld = 0.
                       a_snow_clr = 0.
                       a_snow_cld = 0.                  
                  end where

             end if

        else

             !do code if we are at bottom of column
             if (j .eq. kdim) then
                  surfrain = dtcloud*(rain_clr+rain_cld)
                  surfsnow = dtcloud*(snow_clr+snow_cld)
             end if

        end if 
                  

!-----------------------------------------------------------------------
!
!       END LOOP OVER VERTICAL LEVELS
!

        enddo
     call mpp_clock_end(sc_loop)
     call mpp_clock_begin(sc_post_loop)




!-----------------------------------------------------------------------
!
!       INSTANTANEOUS OUTPUT DIAGNOSTICS
!
     
        if (num_strat_pts > 0) then
         do nn=1,num_strat_pts
          if (strat_pts(1,nn) >= is .and. strat_pts(1,nn) <= ie .and.  &
             strat_pts(2,nn) >= js .and. strat_pts(2,nn) <= je) then
                ipt=strat_pts(1,nn); jpt=strat_pts(2,nn)
                i=ipt-is+1; j=jpt-js+1
                unit = open_ieee32_file ('strat.data', action='append')
                write (unit) ipt,jpt,     ql(i,j,:)+SL(i,j,:)
                write (unit) ipt,jpt,     qi(i,j,:)+SI(i,j,:)
                write (unit) ipt,jpt,     qa(i,j,:)+SA(i,j,:)
                write (unit) ipt,jpt,      T(i,j,:)+ST(i,j,:) 
                write (unit) ipt,jpt,     qv(i,j,:)+SQ(i,j,:)
                write (unit) ipt,jpt,     pfull(i,j,:)
                call close_file(unit)
          endif
         enddo
        endif


!-----------------------------------------------------------------------
!
!       DIAGNOSTICS
!
!rab - perform the assignments for ice/liq adjustments diagnostics
        if (max(diag_id%ice_adj,diag_id%ice_adj_col,diag_id%liq_adj,diag_id%liq_adj_col) .gt. 0) then
          freeze_pt=40.
          if (.not. super_choice) freeze_pt=20.
          where (T .le. tfreeze-freeze_pt)
!          ice_adj = liq_adj
           diag_4d(:,:,:,diag_pt%ice_adj) = diag_4d(:,:,:,diag_pt%liq_adj)
!          liq_adj = 0.
           diag_4d(:,:,:,diag_pt%liq_adj) = 0.
          endwhere
        endif

        if ( diag_id%droplets_wtd > 0 ) then
!a   diag_4d(:,:,:, diag_pt%droplets_wtd) = N3D*ql
     diag_4d(:,:,:, diag_pt%droplets_wtd) =   &
                           diag_4d(:,:,:,diag_pt%droplets)*ql(:,:,:)
        endif

        if ( diag_id%subgrid_w_variance > 0 ) then
!         debug2=debug2**0.5
!         diag_4d(:,:,:,diag_pt%debug2_3d) = debug2
          diag_4d(:,:,:,diag_pt%subgrid_w_variance) = diag_4d(:,:,:,diag_pt%subgrid_w_variance)**0.5
        endif

      
!RSH BUGFIX 6/25/10
!  j is undefined here
        if ( max(diag_id%lcf_strat,diag_id%mfls_strat) > 0 ) then
!         where (omega(:,:,j)+grav*Mc(:,:,j) .lt. 0 .and. lsf_strat(:,:,j) .eq.1)
!            lcf_strat(:,:,j) = 1.
!            diag_4d(:,:,j,diag_pt%lcf_strat) = lcf_strat(:,:,j)
!            mfls_strat(:,:,j) = omega(:,:,j)+grav*Mc(:,:,j)
!            diag_4d(:,:,j,diag_pt%mfls_strat) = mfls_strat(:,:,j)
!         where (omega(:,:,:)+grav*Mc(:,:,:) .lt. 0 .and. lsf_strat(:,:,:) .eq.1)
          where (omega(:,:,:)+grav*Mc(:,:,:) .lt. 0 .and. diag_4d(:,:,:,diag_pt%lsf_strat) .eq.1)
!            lcf_strat(:,:,:) = 1.
!            diag_4d(:,:,:,diag_pt%lcf_strat) = lcf_strat(:,:,:)
             diag_4d(:,:,:,diag_pt%lcf_strat) = 1.0               
!            mfls_strat(:,:,:) = omega(:,:,:)+grav*Mc(:,:,:)
!            diag_4d(:,:,:,diag_pt%mfls_strat) = mfls_strat(:,:,:)
             diag_4d(:,:,:,diag_pt%mfls_strat) = omega(:,:,:)+grav*Mc(:,:,:)
          end where
        end if

       do k = 1, kdim
         do j=1,jdim
           do i=1,idim
             deltpg_3d(i,j,k) = (phalf(i,j,k+1)-phalf(i,j,k))/grav
             if (present(MASK)) then
               deltpg_3d(i,j,k)=deltpg_3d(i,j,k)*MASK(i,j,k)
             endif
           end do
         end do
       end do

        !-------write out column integrated diagnostics------!
!GENERIC FORM:
          do n8=1, n_diag_4d
             do j = kdim, 1, -1
!                 diag_3d(:,:,diag_pt%qldt_evap) = diag_3d(:,:,diag_pt%qldt_evap) &
!                                         + diag_4d(:,:,j,diag_pt%qldt_evap)*deltpg_3d(:,:,j) 
                  diag_3d(:,:,n8) = diag_3d(:,:,n8) &
                                  + diag_4d(:,:,j,n8)*deltpg_3d(:,:,j) 
             enddo
             enddo

!SPECIAL CASES:
!yim: in-cloud droplet column burden

        if (diag_id%droplets_col > 0) then
          if (present (qn)) then
            if (do_liq_num ) then
              N3D_col(:,:) = 0.
              do k = 1, kdim
                do j=1,jdim
                  do i=1,idim
!                   deltpg(i,j) = (phalf(i,j,k+1)-phalf(i,j,k))/grav
!                   if (present(MASK)) then
!                     deltpg(i,j)=deltpg(i,j)*MASK(i,j,k)
!                   endif
                    if (ql(i,j,k) > qmin .and. &
                        qa(i,j,k) > qmin .and. &
                        qn(i,j,k) > qmin ) then      
                      if (qa(i,j,k) > 0.05) then
                       N3D_col(i,j) = N3D_col(i,j) + qn(i,j,k)*  &
!                                     airdens(i,j,k)*deltpg(i,j)*  &
!RSH 12/22/11 fix as per email from yim 11/3/11:
!                                     airdens(i,j,k)*deltpg_3d(i,j,k)*  &
!                                     1.e-6/min(qa(i,j,k),1.)
!RSH 12/22/11
!NOTE still differs from new strat_cloud code in that is appplied to 
! input fields rather than output fields
                                                     deltpg_3d(i,j,k)*  &
                                      1.e-4/min(qa(i,j,k),1.)
                       endif
                    endif
                  end do
                end do
              end do
              diag_3d(:,:,diag_pt%droplets_col) = N3D_col
            endif
          endif
        endif




     call mpp_clock_end(sc_post_loop)
        
!-----------------------------------------------------------------------
!
!
!       end of subroutine



end subroutine strat_cloud_legacy


!#####################################################################

subroutine strat_cloud_legacy_end (do_pdf_clouds)

logical,     intent(in) :: do_pdf_clouds

      if (do_pdf_clouds) call beta_dist_end

      module_is_initialized = .false.

end subroutine strat_cloud_legacy_end
 

!#####################################################################

subroutine aerosol_effects (is, js, Time, phalf, airdens, T, diag_4d, &
                            diag_id, diag_pt, &
                            concen_dust_sub, totalmass1, Aerosol, mask)

type(diag_id_type), intent(in) :: diag_id
type(diag_pt_type), intent(in) :: diag_pt
integer, intent (in)                   :: is,js
type(time_type), intent (in)           :: Time
real, dimension(:,:,:), intent(in )   :: phalf, airdens, T 
real, dimension(:,:,:,0:), intent(inout )   :: diag_4d           
real, dimension(:,:,:), intent(out)        :: concen_dust_sub
real, dimension(:,:,:,:), intent(out)            :: totalmass1
type(aerosol_type), intent (in), optional      :: Aerosol  
real, intent (in), optional, dimension(:,:,:) :: mask

      real, dimension(size(T,1),size(T,2),size(T,3)) :: pthickness
      real, dimension(size(T,1),size(T,2),size(T,3)) :: concen, &
                                    concen_all_sub, &
                                    concen_ss_sub, concen_ss_sup,&
                                    concen_om, concen_na, concen_an
      integer  :: i,j,k,  na , s
      integer  :: idim, jdim, kdim
      logical :: used

      idim = size(T,1)
      jdim = size(T,2)
      kdim = size(T,3)

      concen_dust_sub(:,:,:) = 0.
      totalmass1(:,:,:,:) = 0.

      if (diag_id%sulfate > 0) then
        concen_an(:,:,:) = 0.
        concen_na(:,:,:) = 0.
        concen(:,:,:) = 0.
      endif

      if (use_online_aerosol) then
        concen_ss_sub(:,:,:) = 0.
        concen_ss_sup(:,:,:) = 0.
        concen_all_sub(:,:,:) = 0.
      endif

      do k = 1,kdim
        do j = 1,jdim
          do i = 1,idim
            if (phalf(i,j,k) < 1.0) then
              pthickness(i,j,k) = (phalf(i,j,k+1) - phalf(i,j,k))/&
                                               grav/airdens(i,j,k)
            else
              pthickness(i,j,k) = log(phalf(i,j,k+1)/ &
                            phalf(i,j,k))*8.314*T(i,j,k)/(9.8*0.02888)
            end if
          end do
        end do
      end do

     if (present (Aerosol)) then
       if (do_liq_num) then
         if (use_online_aerosol) then
           do na = 1,size(Aerosol%aerosol,4)               
             if (trim(Aerosol%aerosol_names(na)) == 'so4' .or. &
                 trim(Aerosol%aerosol_names(na)) == 'so4_anthro' .or.&
                 trim(Aerosol%aerosol_names(na)) == 'so4_natural')  &
                                                                 then
               do k = 1,kdim
                 do j = 1,jdim
                   do i = 1,idim
                     totalmass1(i,j,k,1) = totalmass1(i,j,k,1) + &
                                           Aerosol%aerosol(i,j,k,na)
                   end do
                 end do
               end do
             else if(trim(Aerosol%aerosol_names(na)) == 'omphilic' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'omphobic') &
                                                                 then
               do k = 1,kdim
                 do j = 1,jdim
                   do i = 1,idim
                     totalmass1(i,j,k,4) = totalmass1(i,j,k,4) +  &
                                           Aerosol%aerosol(i,j,k,na)
                   end do
                 end do
               end do
             else if(trim(Aerosol%aerosol_names(na)) == 'seasalt1' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'seasalt2') &
                                                                   then
               do k = 1,kdim
                 do j = 1,jdim
                   do i = 1,idim
                     concen_ss_sub(i,j,k) = concen_ss_sub(i,j,k) +  &
                                            Aerosol%aerosol(i,j,k,na)
                   end do
                 end do
               end do
             else if(trim(Aerosol%aerosol_names(na)) == 'seasalt3' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'seasalt4' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'seasalt5')  &
                                                                  then
               do k = 1,kdim
                 do j = 1,jdim
                   do i = 1,idim
                     concen_ss_sup(i,j,k) = concen_ss_sup(i,j,k) +  &
                                            Aerosol%aerosol(i,j,k,na)
                   end do
                 end do
               end do
             else if(trim(Aerosol%aerosol_names(na)) == 'bcphilic' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'bcphobic' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'dust1' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'dust2' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'dust3')  &
                                                                  then
               do k = 1,kdim
                 do j = 1,jdim
                   do i = 1,idim
                     concen_all_sub(i,j,k) = concen_all_sub(i,j,k) +  &
                                             Aerosol%aerosol(i,j,k,na)
                   end do
                 end do
               end do
             endif
             if (do_dust_berg) then
               if (trim(Aerosol%aerosol_names(na)) == 'dust1' .or. &
                   trim(Aerosol%aerosol_names(na)) == 'dust2' .or. &
                   trim( Aerosol%aerosol_names(na)) == 'dust3') then
                 do k = 1,kdim
                   do j = 1,jdim
                     do i = 1,idim
                       concen_dust_sub(i,j,k) =    &
                                           concen_dust_sub(i,j,k) +   &
                                              Aerosol%aerosol(i,j,k,na)
                     end do
                   end do
                 end do
               endif
             endif
           end do
           do k = 1,kdim
             do j = 1,jdim
               do i = 1,idim
                 totalmass1(i,j,k,3) = concen_ss_sub(i,j,k)
                 totalmass1(i,j,k,2) = concen_all_sub(i,j,k) + &
                                       totalmass1(i,j,k,4) + &
                                       concen_ss_sub(i,j,k)
               end do
             end do
           end do
           if (use_sub_seasalt) then
           else
             do k = 1,kdim
               do j = 1,jdim
                 do i = 1,idim
                   totalmass1(i,j,k,3) = concen_ss_sub(i,j,k) +  &
                                                  concen_ss_sup(i,j,k)
                 end do
               end do
             end do
           endif

           if (diag_id%sulfate > 0) then
             do k = 1,kdim
               do j = 1,jdim
                 do i = 1,idim
                   concen(i,j,k) = 0.7273*totalmass1(i,j,k,1)/  &
                                       pthickness(i,j,k)*1.0e9
                 end do
               end do
             end do
           endif
         
           do k = 1,kdim
             do j = 1,jdim
               do i = 1,idim
                 concen_ss_sub(i,j,k) = concen_ss_sub(i,j,k)/  &
                                              pthickness(i,j,k)*1.0e9
                 concen_ss_sup(i,j,k) = concen_ss_sup(i,j,k)/  &
                                              pthickness(i,j,k)*1.0e9
               end do
             end do
           end do

         else  ! (use_online_aerosol)
           if (do_dust_berg) then
!     YMice submicron dust (NO. 14 to NO. 18)
             do s = 14,18
               do k = 1,kdim
                 do j = 1,jdim
                   do i = 1,idim
                     concen_dust_sub(i,j,k) = concen_dust_sub(i,j,k)+ &
                                              Aerosol%aerosol(i,j,k,s)
                   end do
                 end do
               end do
             end do
           endif

           if (diag_id%sulfate > 0) then
             do k = 1,kdim
               do j = 1,jdim
                 do i = 1,idim
!     anthro. and natural sulfate concentration (ug so4/m3)
                   concen_an(i,j,k) = 0.7273*Aerosol%aerosol(i,j,k,1)/&
                                                pthickness(i,j,k)*1.0e9
                   concen_na(i,j,k) = 0.7273*Aerosol%aerosol(i,j,k,2)/&
                                                pthickness(i,j,k)*1.0e9
                   concen(i,j,k) = concen_an(i,j,k) + concen_na(i,j,k)
                 end do
               end do
             end do
           endif

           do k = 1,kdim
             do j = 1,jdim
               do i = 1,idim
!offline
! NO. 1 natural Sulfate; NO. 2 anthro. sulfate; NO. 3 Sea Salt; NO. 4 Or        ganics
                 totalmass1(i,j,k,1) = Aerosol%aerosol(i,j,k,2)
                 totalmass1(i,j,k,2) = Aerosol%aerosol(i,j,k,1)
                 totalmass1(i,j,k,3) = sea_salt_scale*  &
                                       Aerosol%aerosol(i,j,k,5)
                 totalmass1(i,j,k,4) = om_to_oc*  &
                                       Aerosol%aerosol(i,j,k,3)
               end do
             end do
           end do
         endif ! (use_online_aerosol)

         do na = 1, 4
           do k = 1,kdim
             do j = 1,jdim
               do i = 1,idim
                 totalmass1(i,j,k,na) = totalmass1(i,j,k,na)/  &
                                        pthickness(i,j,k)*1.0e9*1.0e-12
               end do
             end do
           end do
         end do
         if (do_dust_berg) then
           do k = 1,kdim
             do j = 1,jdim
               do i = 1,idim
! submicron dust concentration (ug/m3) (NO. 2 to NO. 4)
                 concen_dust_sub(i,j,k) = concen_dust_sub(i,j,k)/ &
                                              pthickness(i,j,k)*1.0e9 
               end do
             end do
           end do
         endif


         if (diag_id%sulfate > 0) then
           diag_4d(:,:,:,diag_pt%sulfate) = concen (:,:,:)
         endif
         if (use_online_aerosol) then
           if (diag_id%seasalt_sub > 0) then
             diag_4d(:,:,:,diag_pt%seasalt_sub) = concen_ss_sub
           endif
           if (diag_id%seasalt_sup > 0) then
             diag_4d(:,:,:,diag_pt%seasalt_sup) = concen_ss_sup
           endif
         endif




         if (diag_id%om > 0) then
           do k = 1,kdim
             do j = 1,jdim
               do i = 1,idim
                 concen_om(i,j,k) = totalmass1(i,j,k,2)*1.0e12
                 diag_4d(i,j,k,diag_pt%om) = concen_om(i,j,k)
               end do
             end do
           end do
!          used = send_data (diag_id%om, diag_4d(:,:,:,diag_pt%om), Time, is, js, 1,&
!                            rmask=mask )
         endif
       endif  ! (do_liq_num)
     endif ! (Present(Aerosol))

!----------------------------------------------------------------------


end subroutine aerosol_effects


!#######################################################################
!#######################################################################




!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  ppm2m_sak --- Piecewise parabolic method for fields
!
! !INTERFACE:
 subroutine ppm2m_sak(a4, delp, km, kmap, i1, i2, iv, kord)

 implicit none

! !INPUT PARAMETERS:
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: i1      ! Starting longitude
 integer, intent(in):: i2      ! Finishing longitude
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: kmap    ! partial remap to start
 integer, intent(in):: kord    ! Order (or more accurately method no.):
                               ! 
 real, intent(in):: delp(i1:i2,km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real, intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

! !DESCRIPTION:
!
!   Perform the piecewise parabolic method 
! 
! !REVISION HISTORY: 
!   ??.??.??    Lin        Creation
!   02.04.04    Sawyer     Newest release from FVGCM
! 
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! local arrays:
      real   dc(i1:i2,km)
      real   h2(i1:i2,km)
      real delq(i1:i2,km)
      real  df2(i1:i2,km)
      real   d4(i1:i2,km)

! local scalars:
      integer i, k, km1, lmt
      integer it
      real fac
      real a1, a2, c1, c2, c3, d1, d2
      real qmax, qmin, cmax, cmin
      real qm, dq, tmp
      real qmp, pmp
      real lac

      km1 = km - 1
       it = i2 - i1 + 1

      do k=max(2,kmap-2),km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo
 
      do k=max(2,kmap-2),km1
         do i=i1,i2
            c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
            c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            qmax = max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1)) - a4(1,i,k)
            qmin = a4(1,i,k) - min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))
             dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
            df2(i,k) = tmp
         enddo
      enddo

!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------

      do k=max(3,kmap), km1
      do i=i1,i2
        c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
        a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
        a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
        a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                  ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                                delp(i,k-1)*a1*dc(i,k  ) )
      enddo
      enddo

      if(km>8 .and. kord>3) call steepz_sak(i1, i2, km, kmap, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      if ( kmap <= 2 ) then
      do i=i1,i2
         d1 = delp(i,1)
         d2 = delp(i,2)
         qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
         dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
         c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1**2)
         a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
         a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
         dc(i,1) =  a4(1,i,1) - a4(2,i,1)
! No over- and undershoot condition
         cmax = max(a4(1,i,1), a4(1,i,2))
         cmin = min(a4(1,i,1), a4(1,i,2))
         a4(2,i,2) = max(cmin,a4(2,i,2))
         a4(2,i,2) = min(cmax,a4(2,i,2))
      enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
         dc(i,km) = a4(3,i,km) -  a4(1,i,km)
! No over- and under-shoot condition
         cmax = max(a4(1,i,km), a4(1,i,km1))
         cmin = min(a4(1,i,km), a4(1,i,km1))
         a4(2,i,km) = max(cmin,a4(2,i,km))
         a4(2,i,km) = min(cmax,a4(2,i,km))
      enddo

      do k=max(1,kmap),km1
         do i=i1,i2
            a4(3,i,k) = a4(2,i,k+1)
         enddo
      enddo

! Enforce monotonicity of the "slope" within the top layer
      if ( kmap <= 2 ) then
      do i=i1,i2
         if ( a4(2,i,1) * a4(1,i,1) <= 0. ) then 
              a4(2,i,1) = 0.
                dc(i,1) = a4(1,i,1)
         endif
         if ( dc(i,1) * (a4(2,i,2) - a4(1,i,1)) <= 0. ) then
! Setting DC==0 will force piecewise constant distribution after
! calling kmppm_sak
              dc(i,1) = 0.
         endif
      enddo
      endif

! Enforce constraint on the "slope" at the surface

      do i=i1,i2
         if( a4(3,i,km) * a4(1,i,km) <= 0. ) then
!            a4(3,i,km) = 0.
!              dc(i,km) =  -a4(1,i,km)
               dc(i,km) = 0.
         endif
         if( dc(i,km) * (a4(1,i,km) - a4(2,i,km)) <= 0. ) then
             dc(i,km) = 0.
         endif
      enddo
 
!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      if ( kmap <= 2 ) then
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
            call kmppm_sak(dc(i1,k), a4(1,i1,k), it, 0)
      enddo
      endif

      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=max(2,kmap-1), km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2
!           h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))
!    &               / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )
!    &               * delp(i,k)**2
! Method#3
            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      if( kord == 7 ) then
         fac = 1.5           ! original quasi-monotone
      else
         fac = 0.125         ! full monotone
      endif

      do k=max(3,kmap), km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(3,i,k) = min(max(a4(3,i,k), qmin), qmax)
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(2,i,k) = min(max(a4(2,i,k), qmin), qmax)
! Recompute A6
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to ensure positivity when kord=7
         if (iv == 0 .and. kord == 7) then
             call kmppm_sak(dc(i1,k), a4(1,i1,k), it, 2)
         endif
      enddo

      else
 
         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

      do k=max(3,kmap), km-2
      if( kord /= 4) then
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
      endif
         call kmppm_sak(dc(i1,k), a4(1,i1,k), it, lmt)
      enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm_sak(dc(i1,k), a4(1,i1,k), it, 0)
      enddo
!EOC
 end subroutine ppm2m_sak
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  kmppm_sak --- Perform piecewise parabolic method in vertical
!
! !INTERFACE:
 subroutine kmppm_sak(dm, a4, itot, lmt)

 implicit none

! !INPUT PARAMETERS:
      real, intent(in):: dm(*)     ! the linear slope
      integer, intent(in) :: itot      ! Total Longitudes
      integer, intent(in) :: lmt       ! 0: Standard PPM constraint
                                       ! 1: Improved full monotonicity constraint (Lin)
                                       ! 2: Positive definite constraint
                                       ! 3: do nothing (return immediately)
! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout) :: a4(4,*)   ! PPM array
                                           ! AA <-- a4(1,i)
                                           ! AL <-- a4(2,i)
                                           ! AR <-- a4(3,i)
                                           ! A6 <-- a4(4,i)

! !DESCRIPTION:
!
! !REVISION HISTORY: 
!    00.04.24   Lin       Last modification
!    01.03.26   Sawyer    Added ProTeX documentation
!    02.04.04   Sawyer    Incorporated newest FVGCM version
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      real, parameter:: r12 = 1./12.
      real qmp
      real da1, da2, a6da
      real fmin
      integer i

! Developer: S.-J. Lin, NASA-GSFC
! Last modified: Apr 24, 2000

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin 2003)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin < 0. ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

!EOC
 end subroutine kmppm_sak
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  steepz_sak --- Calculate attributes for PPM
!
! !INTERFACE:
 subroutine steepz_sak(i1, i2, km, kmap, a4, df2, dm, dq, dp, d4)

   implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: km                   ! Total levels
      integer, intent(in) :: kmap                 ! 
      integer, intent(in) :: i1                   ! Starting longitude
      integer, intent(in) :: i2                   ! Finishing longitude
      real, intent(in) ::  dp(i1:i2,km)       ! grid size
      real, intent(in) ::  dq(i1:i2,km)       ! backward diff of q
      real, intent(in) ::  d4(i1:i2,km)       ! backward sum:  dp(k)+ dp(k-1) 
      real, intent(in) :: df2(i1:i2,km)       ! first guess mismatch
      real, intent(in) ::  dm(i1:i2,km)       ! monotonic mismatch

! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout) ::  a4(4,i1:i2,km)  ! first guess/steepened

!
! !DESCRIPTION:
!   This is complicated stuff related to the Piecewise Parabolic Method
!   and I need to read the Collela/Woodward paper before documenting
!   thoroughly.
!
! !REVISION HISTORY: 
!   ??.??.??    Lin?       Creation
!   01.03.26    Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, k
      real alfa(i1:i2,km)
      real    f(i1:i2,km)
      real  rat(i1:i2,km)
      real  dg2

! Compute ratio of dq/dp
      do k=max(2,kmap-1),km
         do i=i1,i2
            rat(i,k) = dq(i,k-1) / d4(i,k)
         enddo
      enddo

! Compute F
      do k=max(2,kmap-1),km-1
         do i=i1,i2
            f(i,k) =   (rat(i,k+1) - rat(i,k))                          &
                     / ( dp(i,k-1)+dp(i,k)+dp(i,k+1) )
         enddo
      enddo

      do k=max(3,kmap),km-2
         do i=i1,i2
         if(f(i,k+1)*f(i,k-1)<0. .and. df2(i,k)/=0.) then
            dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2          &
                   + d4(i,k)*d4(i,k+1) )
            alfa(i,k) = max(0., min(0.5, -0.1875*dg2/df2(i,k))) 
         else
            alfa(i,k) = 0.
         endif
         enddo
      enddo

      do k=max(4,kmap+1),km-2
         do i=i1,i2
            a4(2,i,k) = (1.-alfa(i,k-1)-alfa(i,k)) * a4(2,i,k) +         &
                        alfa(i,k-1)*(a4(1,i,k)-dm(i,k))    +             &
                        alfa(i,k)*(a4(1,i,k-1)+dm(i,k-1))
         enddo
      enddo

!EOC
 end subroutine steepz_sak
!-----------------------------------------------------------------------

        subroutine cloud_clear_xfer ( cloud_generator_on,tmp3, qa_mean, qa_mean_lst, a_clr, a_cld, clr, cld)
        logical, intent(in) :: cloud_generator_on
        real, dimension(:,:), intent(in) :: tmp3, qa_mean, qa_mean_lst
        real, dimension(:,:), intent(inout) :: a_clr, a_cld, clr, cld
        real :: cld2clr, clr2cld, prec_cld2clr, prec_clr2cld, tmp1, tmp2
        integer :: k,i,kdim,idim

        idim = size(tmp3,1)
        kdim = size(tmp3,2)

        do k=1,kdim
         do i=1,idim
        !-------------------------------
        !compute cloud to clear transfer
          if (overlap .eq. 1)                                            &
             cld2clr= min(a_cld(i,k),max(0.,a_cld(i,k) - qa_mean(i,k))   )

          if (overlap .eq. 2)                                            &
             cld2clr= min(a_cld(i,k),max(0.,a_cld(i,k)*(1.-qa_mean(i,k))))

          if (cloud_generator_on) then
             tmp1 =      min(a_cld(i,k),max(0.,a_cld(i,k) - qa_mean(i,k))   )
             tmp2 =      min(a_cld(i,k),max(0.,a_cld(i,k)*(1.-qa_mean(i,k))))
             cld2clr=min(a_cld(i,k),max(0.,tmp3(i,k)*tmp1+(1.-tmp3(i,k))*tmp2))
          end if

        !-------------------------------
        !compute clear to cloud transfer
          if (overlap .eq. 1)                                            &
             clr2cld = min(max(qa_mean(i,k)-qa_mean_lst(i,k),0.),a_clr(i,k))
          if (overlap .eq. 2)                                            &
             clr2cld = min(max( a_clr(i,k)*qa_mean(i,k),0.),a_clr(i,k))

          if (cloud_generator_on) then
             tmp1 =       min(max(qa_mean(i,k)-qa_mean_lst(i,k),0.),a_clr(i,k))
             tmp2 =       min(max( a_clr(i,k)*qa_mean(i,k),0.),a_clr(i,k))
             clr2cld=min(a_clr(i,k),max(0.,tmp3(i,k)*tmp1+(1.-tmp3(i,k))*tmp2))
          end if

        !---------------------------------
        !calculate precipitation transfers
          prec_cld2clr = cld(i,k)*(cld2clr/max(a_cld(i,k),qmin))
          prec_clr2cld = clr(i,k)*(clr2cld/max(a_clr(i,k),qmin))

        !----------------
        !add in transfers
          a_clr(i,k) = a_clr(i,k) + cld2clr - clr2cld
          a_cld(i,k) = a_cld(i,k) - cld2clr + clr2cld
          clr(i,k)   = clr(i,k) + prec_cld2clr - prec_clr2cld
          cld(i,k)   = cld(i,k) - prec_cld2clr + prec_clr2cld
         enddo
        enddo

        end subroutine cloud_clear_xfer

!

!#######################################################################
!#######################################################################



end module strat_cloud_legacy_mod

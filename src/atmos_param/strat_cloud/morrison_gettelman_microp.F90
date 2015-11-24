MODULE morrison_gettelman_microp_mod

use sat_vapor_pres_mod,        only : lookup_des, compute_qs,  &
                                      sat_vapor_pres_init
use polysvp_mod,               only : polysvp_l,  polysvp_i,   &
                                      polysvp_init,  polysvp_end
use constants_mod,             only : pi, grav, rvgas, rdgas, tfreeze, &
                                      hlv, hlf, hls, cp_air
use mpp_mod,                   only : input_nml_file
use fms_mod,                   only : mpp_pe, file_exist, error_mesg,  &
                                      open_namelist_file, FATAL, &
                                      stdlog, write_version_number, &
                                      check_nml_error, close_file, &
                                      mpp_root_pe
use strat_cloud_utilities_mod, only : strat_cloud_utilities_init, &
                                      diag_id_type, diag_pt_type, &
                                      strat_nml_type
use  gamma_mg_mod,             only : gamma_mg, gamma_mg_init, gamma_mg_end
use mg_const_mod,              only : mg_const_init, rhow, rhoi,  &
                                      mg_pr, di_mg, ci_mg
use simple_pdf_mod,            only : simple_pdf, simple_pdf_init, &
                                      simple_pdf_end
  
implicit none
private 

!------------------------------------------------------------------------
!--interfaces------------------------------------------------------------

public morrison_gettelman_microp, morrison_gettelman_microp_init,  &
       morrison_gettelman_microp_end

!------------------------------------------------------------------------
!--version number--------------------------------------------------------

Character(len=128) :: Version = '$Id: morrison_gettelman_microp.F90,v 20.0 2013/12/13 23:21:59 fms Exp $'
Character(len=128) :: Tagname = '$Name: tikal $'

!------------------------------------------------------------------------
!--namelist--------------------------------------------------------------

logical           :: do_morrison_gettelman_eros = .true.
integer           :: auto_conv_ice_choice = 1
real              :: auto_conv_time_scale = 180._mg_pr
real              :: auto_conv_m_thresh =  1.e-4_mg_pr 
integer           :: act_choice  = 1     ! M+G super_act
logical           :: super_act = .false. ! super_act is more consistent 
                                         ! with Yi's formula.
real              :: ai = 700._mg_pr     ! a in the fallspeed relationship
                                         ! V=a*D^b for cloud ice 
                                         ! (mg two moment scheme)
real              :: bi = 1._mg_pr       ! b in the fallspeed relationship
                                         ! V=a*D^b for cloud ice  
                                         ! (mg two moment scheme)
real              :: as = 11.72_mg_pr   ! a in the fallspeed relationship
                                         ! V=a*D^b for cloud ice 
                                         ! (mg two moment scheme)
real              :: bs = 0.41_mg_pr     ! b in the fallspeed relationship
                                         ! V=a*D^b for cloud ice 
                                         ! (mg two moment scheme)
real              :: in_cloud_limit = 5.e-3_mg_pr
logical           :: orig_app_test = .false.
logical           :: sub_cld_var = .true.
logical           :: do_berg_snow = .false.
logical           :: do_excess1   = .false.
logical           :: hd_sedi_sens = .false.

real              :: vfact            = 1.0_mg_pr
real              :: vfact_n          = 1.0_mg_pr
real              :: vfact_m          = 1.0_mg_pr
real              :: max_vt_ice       = 1.2_mg_pr
real              :: max_vt_snow       = 1.2_mg_pr
real              :: nucl_thresh = 208.9e3_mg_pr
integer           :: no_rh_adj_opt = 0

real              :: min_diam_ice = 10.e-6_mg_pr
real              :: min_diam_drop = 2.e-6_mg_pr
real              :: max_diam_drop = 50.e-6_mg_pr

real              :: max_diam_ii = 2400.e-6_mg_pr
real              :: min_diam_ii = 1.e-6_mg_pr
logical           :: rain_evap_opt  = .false.
logical           :: subl_snow    = .false.
logical           :: mass_cons = .false.

logical           :: collect_frzreg = .true.
logical           :: do_contact_frz = .false.
integer           :: n_contact_opt = 2
logical           :: do_bigg_frz = .true.
real              :: limit_bigg_t = 0._mg_pr
logical           :: limit_volri = .true.
logical           :: one_ice = .false.
logical           :: scav_by_cloud_ice = .false.
real              :: tmin_fice =  tfreeze - 40._mg_pr ! min temperature for
                                         ! ice deposition/bergeron process
real              :: Dcs  = 200.e-6_mg_pr
real              :: qsmall =  1.e-14_mg_pr ! smallest mixing ratio 
                                            ! considered in microphysics
logical           :: tiedtke_qa_test = .false.
logical           :: qv_on_qi = .false.
integer           :: sat_adj_opt = 1
real              :: autoconv_ice_thr  = 100.e-6_mg_pr
real              :: Eii = 0.1_mg_pr
real              :: size_hom = 25.e-6_mg_pr
logical           :: limit_berg = .false.
logical           :: do_nevap = .true. 
integer           :: limit_droplet_freeze_opt = 1

real              :: berg_lim = 1.e-6_mg_pr 
real              :: rhosn = 100._mg_pr       ! bulk density snow  (from 
                                              ! Reisner et al. (1998))

logical           :: mg_repartition_first = .true.
logical           :: meyers_test = .false.
logical    :: allow_all_cldtop_collection = .false.
logical    :: rho_factor_in_max_vt = .true.
real       :: max_rho_factor_in_vt = 1.0
real       :: lowest_temp_for_sublimation = 180._mg_pr

namelist / morrison_gettelman_microp_nml /   &
                 do_morrison_gettelman_eros, auto_conv_ice_choice,   &
                 auto_conv_time_scale,  auto_conv_m_thresh,              &
                 ai, bi, as , bs, act_choice, super_act, in_cloud_limit, &
                 do_berg_snow, do_excess1, hd_sedi_sens, vfact, vfact_n, &
                 vfact_m, max_vt_ice, max_vt_snow, sub_cld_var,          &
                 nucl_thresh, no_rh_adj_opt, min_diam_ice,     &
                 min_diam_drop, max_diam_drop, rain_evap_opt, mass_cons, &
                 subl_snow, collect_frzreg, one_ice, tmin_fice, Dcs,     &
                 qsmall, tiedtke_qa_test, qv_on_qi, scav_by_cloud_ice,  &
                 max_diam_ii, min_diam_ii, sat_adj_opt, autoconv_ice_thr, &
                 eii, size_hom, limit_berg, do_nevap, berg_lim, rhosn, &
                 do_contact_frz, n_contact_opt, do_bigg_frz, limit_bigg_t,&
                 limit_volri, limit_droplet_freeze_opt,    &
                 mg_repartition_first, meyers_test, &
                 allow_all_cldtop_collection, rho_factor_in_max_vt,&
!                allow_all_cldtop_collection,                      &
                 max_rho_factor_in_vt, &
                 lowest_temp_for_sublimation

!-----------------------------------------------------------------------
!-----internal parameters and constants---------------------------------
 
!-----------------------------------------------------------------------
! radius of contact nuclei aerosol (m) :
REAL(kind=mg_pr), PARAMETER  :: rin = 0.1e-6_mg_pr

!-----------------------------------------------------------------------
! cloud droplet fall speed parameters, V = aD^b;  V is in m/s
REAL(kind=mg_pr), PARAMETER  :: ac = 3.e7_mg_pr
REAL(kind=mg_pr), PARAMETER  :: bc = 2._mg_pr


!-----------------------------------------------------------------------
! rain drop fall speed parameters, V = aD^b;  V is in m/s
REAL(kind=mg_pr), PARAMETER  :: ar = 841.99667_mg_pr
REAL(kind=mg_pr), PARAMETER  :: br = 0.8_mg_pr

!-----------------------------------------------------------------------
! typical air density at 850 mb
REAL(kind=mg_pr), PARAMETER  :: rhosu = 85000._mg_pr/(rdgas*tfreeze)

!-----------------------------------------------------------------------
! immersion freezing parameters, bigg 1953
REAL(kind=mg_pr), PARAMETER  :: bimm = 100._mg_pr
REAL(kind=mg_pr), PARAMETER  :: aimm = 0.66_mg_pr

!------------------------------------------------------------------------
! mass of new crystal due to aerosol freezing and growth (kg)
REAL(kind=mg_pr)             :: mi0

!------------------------------------------------------------------------
! snow mass-diameter relationship
REAL(kind=mg_pr)             :: cs
REAL(kind=mg_pr), PARAMETER  :: ds = 3._mg_pr

!------------------------------------------------------------------------
! collection efficiency, accretion of cloud water by rain
REAL(kind=mg_pr), PARAMETER  :: Ecr = 1.0_mg_pr

!-----------------------------------------------------------------------
! radius for homogeneously frozen droplets
!REAL(kind=mg_pr), PARAMETER  :: homog_frz_radius = 25.e-6

!-----------------------------------------------------------------------
! ventilation constants for rain
REAL(kind=mg_pr), PARAMETER  :: f1r = 0.78_mg_pr
REAL(kind=mg_pr), PARAMETER  :: f2r = 0.32_mg_pr

!-------------------------------------------------------------------------
! ventilation parameters for snow
! hall and prupacher
REAL(kind=mg_pr), PARAMETER  :: f1s = 0.86_mg_pr
REAL(kind=mg_pr), PARAMETER  :: f2s = 0.28_mg_pr

!------------------------------------------------------------------------
! ratio of h2o to dry air molecular weight
REAL(kind=mg_pr), PARAMETER  :: epsqs = 0.62197_mg_pr

!-----------------------------------------------------------------------
!NOTE:
! latent heats should probably be fixed with temperature 
! for energy conservation with the rest of the model
! (this looks like a +/- 3 or 4% effect, but will mess up energy balance)
! xxlv - latent heat of vaporization
! xlf  - latent heat of freezing
! xxls - xxlv + xlf, latent heat of sublimation
REAL(kind=mg_pr), PARAMETER  :: xxlv = hlv       
REAL(kind=mg_pr), PARAMETER  :: xlf = hlf   
REAL(kind=mg_pr), PARAMETER  :: xxls = hls

!------------------------------------------------------------------------
! water vapor gas contstant
REAL(kind=mg_pr), PARAMETER  :: rv= rvgas  

!------------------------------------------------------------------------
! specific heat of dry air at const. press
REAL(kind=mg_pr), PARAMETER  ::  cpp = cp_air         

!------------------------------------------------------------------------
! constants related to subgrid cloud variance
REAL (kind=mg_pr)            :: sfac1, sfac2, sfac3, sfac4, sfac5


!------------------------------------------------------------------------

logical    :: module_is_initialized = .false.

! 1 / relative variance of sub-grid cloud water distribution
! see morrison and gettelman, 2007, J. Climate for details
real              :: qcvar 


CONTAINS



!#########################################################################

SUBROUTINE morrison_gettelman_microp_init (do_pdf_clouds, qcvar_in)

!-----------------------------------------------------------------------
LOGICAL, INTENT(IN ) :: do_pdf_clouds
Real,    INTENT(IN ) :: qcvar_in      

!-----------------------------------------------------------------------
!--local variables------------------------------------------------------

      INTEGER   :: unit, io, ierr, logunit

!-----------------------------------------------------------------------
      if (module_is_initialized) return

      qcvar = qcvar_in

!-----------------------------------------------------------------------
!    process namelist.
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=morrison_gettelman_microp_nml, iostat=io)
      ierr = check_nml_error(io,'morrison_gettelman_microp_nml')
#else
      if ( file_exist('input.nml')) then
        unit = open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=morrison_gettelman_microp_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'morrison_gettelman_microp_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!-----------------------------------------------------------------------
!    write version and namelist to stdlog.
!-----------------------------------------------------------------------
      call write_version_number(version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe()) &
                       write (logunit, nml=morrison_gettelman_microp_nml)

!------------------------------------------------------------------------
!    be sure needed modules have been initilized.
!------------------------------------------------------------------------
      CALL sat_vapor_pres_init
      call mg_const_init
      call strat_cloud_utilities_init
      CALL polysvp_init
      IF (do_pdf_clouds) THEN
        CALL  simple_pdf_init
      END IF
      CALL gamma_mg_init

!-----------------------------------------------------------------------
!    calculate some constants.
!-----------------------------------------------------------------------
      mi0 = 4._mg_pr/3._mg_pr*pi*rhoi*min_diam_ice**3
      IF (sub_cld_var) THEN
        sfac1 = gamma_mg(qcvar+2.47_mg_pr)/(gamma_mg(qcvar)*  &
                                                       qcvar**2.47_mg_pr)
        sfac2 = gamma_mg(qcvar+2._mg_pr)/(gamma_mg(qcvar)*qcvar**2)
        sfac3 = gamma_mg(qcvar+1._mg_pr)/(gamma_mg(qcvar)*qcvar)
        sfac4 = gamma_mg(qcvar+1.15_mg_pr)/(gamma_mg(qcvar)*   &
                                                       qcvar**1.15_mg_pr)
        sfac5 = gamma_mg(qcvar+bc/3._mg_pr)/(gamma_mg(qcvar)*  &
                                                      qcvar**(bc/3._mg_pr))
      ELSE
        sfac1 = 1._mg_pr
        sfac2 = 1._mg_pr 
        sfac3 = 1._mg_pr
        sfac4 = 1._mg_pr
        sfac5 = 1._mg_pr
      END IF
      cs = rhosn*pi/6._mg_pr

!-----------------------------------------------------------------------
!    mark the module as initialized.
!-----------------------------------------------------------------------
      module_is_initialized = .true.

!----------------------------------------------------------------------


END SUBROUTINE morrison_gettelman_microp_init



!########################################################################

SUBROUTINE morrison_gettelman_microp  (      &
                  tiedtke_macrophysics, total_activation, dqa_activation, &
                  ncall, j, idim, jdim, kdim, Nml, deltatin, pfull, &
                  pdel, t_in, t0, qv_in, qv0, qc_in, qi_in, nc_in, ni_in, &
                  cldn_in, dqcdt, dqidt, drop2, crystal1, rbar_dust,   &
                  ndust, delta_cf, qa_upd, qa_upd_0, SA_0, D_eros_l4, &
                  nerosc4,  D_eros_i4, nerosi4, gamma, inv_dtcloud, qa0, & 
                  ST, SQ,  ssat_disposal, tlat, qvlat, qctend, qitend,    &
                  nctend, nitend, SA, rain3d, snow3d, prect, preci,   &
                  qrout, qsout, lsc_rain_size, lsc_snow_size,   &
                  f_snow_berg, n_diag_4d, diag_4d, diag_id, diag_pt,  &
                  nrefuse, debugo0, debugo1, otun)
       
!------------------------------------------------------------------------
logical,                 intent(in)    ::  tiedtke_macrophysics,  &
                                           total_activation, &
                                           dqa_activation
INTEGER,                 INTENT(IN )   ::  ncall, j, idim, jdim, kdim
type(strat_nml_type),    intent(in)    ::  Nml
REAL(kind=mg_pr),        INTENT(IN )   ::  deltatin        ! time step (s)
REAL(kind=mg_pr), DIMENSION(idim,kdim),   INTENT(IN)    & 
                                       ::  pfull ! air pressure (pa)
REAL(kind=mg_pr), DIMENSION(idim,kdim),   INTENT(IN)    &
                                       ::  pdel  ! pressure difference 
                                                 ! across level (pa)
REAL(kind=mg_pr), DIMENSION(idim,kdim),   INTENT(IN)    &
                                       ::  t_in  ! input temperature (K)
REAL(kind=mg_pr), DIMENSION(idim,kdim),   INTENT(IN)    &
                                       ::  qv_in ! input h20 vapor mixing 
                                                 ! ratio (kg/kg)
REAL(kind=mg_pr), DIMENSION(idim,kdim),   INTENT(IN)    &
                                       ::  T0, qv0
REAL(kind=mg_pr), DIMENSION(idim,kdim),   INTENT(IN)     &
                                       ::  qc_in ! cloud water mixing 
                                                 ! ratio (kg/kg)  
REAL(kind=mg_pr), DIMENSION(idim,kdim),   INTENT(IN)   &
                                       ::  qi_in ! cloud ice mixing    
                                                 ! ratio (kg/kg)
REAL(kind=mg_pr), DIMENSION(idim,kdim),   INTENT(IN)   &
                                       ::  nc_in ! cloud droplet number 
                                                 ! conc (1/kg) at the 
                                                 ! beginning of the 
                                                 ! timestep prior to 
                                                 ! activation calc.
REAL(kind=mg_pr), DIMENSION(idim,kdim),   INTENT(IN)   &
                                       ::  ni_in ! cloud ice number conc 
                                                 ! (1/kg)
REAL(kind=mg_pr), DIMENSION(idim,kdim),   INTENT(IN)    &
                                       ::  cldn_in 
                                                 ! cloud fraction
REAL(kind=mg_pr), DIMENSION(idim,kdim),  INTENT(IN)    &
                                       :: dqcdt, dqidt
REAL(kind=mg_pr), DIMENSION(idim,kdim),  INTENT(INOUT)   &
                                       ::  drop2, crystal1
REAL(kind=mg_pr), DIMENSION(idim,kdim),  INTENT(IN)     &
                                       ::  rbar_dust, ndust
REAL(kind=mg_pr), DIMENSION(idim,kdim),  INTENT(INOUT)   &
                                       ::  delta_cf
REAL(kind=mg_pr), DIMENSION(idim,kdim),  INTENT(INOUT)   &
                                       ::  qa_upd
REAL(kind=mg_pr), DIMENSION(idim,kdim),  INTENT(INOUT )   &
                                       ::  qa_upd_0, SA_0
REAL(kind=mg_pr), DIMENSION(idim,kdim),  INTENT(IN)    &
                                       ::  D_eros_l4, D_eros_i4
REAL(kind=mg_pr), DIMENSION(idim,kdim),  INTENT(IN)    &
                                       ::  nerosc4, nerosi4
REAL(kind=mg_pr), DIMENSION(idim,kdim),  INTENT(IN)    &
                                       ::  gamma
REAL(kind=mg_pr),                        INTENT(IN )   &
                                       ::  inv_dtcloud       
REAL(kind=mg_pr), DIMENSION(idim,kdim),  INTENT(IN)    &
                                       ::  qa0
REAL(kind=mg_pr), DIMENSION(idim,kdim),  INTENT(IN)    &
                                       ::  ST, SQ
REAL(kind=mg_pr), DIMENSION(idim,kdim),intent(out)     &
                                       :: ssat_disposal
REAL(kind=mg_pr), DIMENSION(idim,kdim),intent(out)    &
                                       :: tlat   ! latent heating rate 
                                                 !    (K/s)
REAL(kind=mg_pr), DIMENSION(idim,kdim),intent(out)     &
                                       :: qvlat  ! microphysical  tendency
                                                 !  qv (1/s)
REAL(kind=mg_pr), DIMENSION(idim,kdim),  INTENT(OUT)    &
                                       :: qctend, qitend, nctend, nitend
REAL(kind=mg_pr), DIMENSION(idim,kdim),  INTENT(INOUT)    &
                                       :: SA
! snow and rain mass mixing ratio (kg/kg) 
real, dimension(idim, jdim, kdim+1), INTENT(INOUT)      &
                                       :: rain3d, snow3d
REAL(kind=mg_pr), DIMENSION(idim), intent(out)       &
                                       :: prect  ! surface precip rate 
                                                 ! (m/s)
REAL(kind=mg_pr), DIMENSION(idim), intent(out)    &
                                       :: preci  ! cloud ice/snow precip 
                                                 ! rate (m/s)
! diagnostic rain/snow for output to history
! values are in-precip (local) !!!!
REAL(kind=mg_pr), DIMENSION(idim,kdim), INTENT(OUT )      &
                                       :: qrout ! rain mixing ratio (kg/kg)
REAL(kind=mg_pr), DIMENSION(idim,kdim), INTENT(OUT )      &
                                       :: qsout ! snow mixing ratio (kg/kg)
! snow and rain diameter (micrometer):
! (snow_size currently not used) 
REAL, dimension(idim, kdim), INTENT(INOUT)     &
                                       :: lsc_snow_size, lsc_rain_size
REAL(kind=mg_pr), DIMENSION(idim,kdim),  INTENT(OUT)     &
                                       ::  f_snow_berg
INTEGER,                                 intent(in)     &
                                       :: n_diag_4d
REAL, dimension( idim, jdim, kdim, 0:n_diag_4d ), INTENT(INOUT )    &
                                       ::  diag_4d
TYPE(diag_id_type),                       intent(in)     &
                                       :: diag_id
TYPE(diag_pt_type),                       intent(inout)    &
                                       :: diag_pt
INTEGER,                                  INTENT(INOUT)     &
                                       :: nrefuse
LOGICAL,                                  INTENT(IN )    &
                                       :: debugo0, debugo1 
INTEGER,                                  INTENT(IN )     &
                                       :: otun

!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!---local variables-------------------------------------------------------
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: tn !input temperature

      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  D_eros_l, D_eros_i 

      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: dz ! height difference 
                                               !across model vertical level


      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: qv, qtot

      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: qc, qi

      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nc, ni



! sum of source/sink terms for diagnostic prec
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: qnitend ! snow mixing ratio source/sink term
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nstend  ! snow number concentration source/sink term
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: qrtend ! rain mixing ratio source/sink term
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nrtend  ! rain number concentration source/sink term



! new terms for Bergeron process

      REAL(kind=mg_pr) :: dumnnuc ! provisional ice nucleation rate (for calculating bergeron)
      REAL(kind=mg_pr) :: ninew  ! provisional cloud ice number conc (for calculating bergeron)
      REAL(kind=mg_pr) :: qinew ! provisional cloud ice mixing ratio (for calculating bergeron)
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: berg  ! mixing rat tendency due to bergeron process for cloud ice

!cms++
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: qvdep_qi   ! mixing rat tendency due to vapor deposition on ice   
!cms--



      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: esl ! liquid sat vapor pressure (pa)
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: esi ! ice sat vapor pressure (pa)


      REAL(kind=mg_pr), DIMENSION(kdim) ::  nsubi ! evaporation of cloud ice number
      REAL(kind=mg_pr), DIMENSION(kdim) ::  nsubc ! evaporation of droplet number




      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  nerosi! "erosion" of cloud ice number
      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  nerosc! "erosion" of cloud ice number
      REAL(kind=mg_pr), DIMENSION(kdim) ::  nsubs ! evaporation of snow number
      REAL(kind=mg_pr), DIMENSION(kdim) ::  nsubr ! evaporation of rain number


      REAL(kind=mg_pr), DIMENSION(kdim) :: nnuccd   ! ice nucleation rate from various freezing processes
      REAL(kind=mg_pr), DIMENSION(kdim) :: mnuccd   ! ice nucleation rate from various freezing processes
      REAL(kind=mg_pr), DIMENSION(kdim) :: npccn    ! droplet activation rate



      REAL(kind=mg_pr) :: qrtot ! vertically-integrated rain mixing rat source/sink term
      REAL(kind=mg_pr) :: nrtot ! vertically-integrated rain number conc source/sink term
      REAL(kind=mg_pr) :: qstot ! vertically-integrated snow mixing rat source/sink term
      REAL(kind=mg_pr) :: nstot ! vertically-integrated snow number conc source/sink term

      REAL(kind=mg_pr) :: deltat   ! sub-time step (s)  
  

      REAL(kind=mg_pr) :: omsm    ! number near unity for round-off issues
      REAL(kind=mg_pr) :: mincld  ! minimum allowed cloud fraction

      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: q ! water vapor mixing ratio (kg/kg)
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: t ! temperature (K)
      REAL(kind=mg_pr) , DIMENSION(idim,kdim):: rho ! air density (kg m-3)


      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: dv  ! diffusivity of water vapor in air
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: mu  ! viscocity of air
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: sc  ! schmidt number
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: kap ! thermal conductivity of air
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: rhof! air density correction factor for fallspeed
      REAL(kind=mg_pr) :: dap ! effecvtive diffusivity of contact ice nuclei
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: arn ! air density corrected rain fallspeed parameter
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: asn ! air density corrected snow fallspeed parameter
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: acn ! air density corrected cloud droplet fallspeed parameter
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: ain ! air density corrected cloud ice fallspeed parameter


      REAL(kind=mg_pr), DIMENSION(idim,kdim)  :: qcic ! in-cloud cloud liquid mixing ratio
      REAL(kind=mg_pr), DIMENSION(idim,kdim)  :: qiic ! in-cloud cloud ice mixing ratio
      REAL(kind=mg_pr), DIMENSION(idim,kdim)  :: qniic ! in-precip snow mixing ratio
      REAL(kind=mg_pr), DIMENSION(idim,kdim)  :: qric ! in-precip rain mixing ratio
      REAL(kind=mg_pr), DIMENSION(idim,kdim)  :: ncic ! in-cloud droplet number conc
      REAL(kind=mg_pr), DIMENSION(idim,kdim)  :: niic ! in-cloud cloud ice number conc
      REAL(kind=mg_pr), DIMENSION(idim,kdim)  :: nsic ! in-precip snow number conc
      REAL(kind=mg_pr), DIMENSION(idim,kdim)  :: nric ! in-precip rain number conc

        ! hm, add9/5/07, rain rate for reflectivity calculation
      REAL(kind=mg_pr), DIMENSION(idim,kdim)  :: rainrt
      REAL(kind=mg_pr), DIMENSION(idim,kdim)  :: rainrt1


      REAL(kind=mg_pr), DIMENSION(idim,kdim)  :: atotrt
      REAL(kind=mg_pr), DIMENSION(idim,kdim)  :: atotrt1

      REAL(kind=mg_pr), DIMENSION(idim,kdim)  :: asnowrt
      REAL(kind=mg_pr), DIMENSION(idim,kdim)  :: asnowrt1



      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: cwml   ! cloud water mixing ratio
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: cwmi   ! cloud ice mixing ratio

      REAL(kind=mg_pr), DIMENSION(kdim) :: snow2vapor

      REAL(kind=mg_pr) :: uni ! number-weighted cloud ice fallspeed
      REAL(kind=mg_pr), DIMENSION(kdim) ::  umi ! mass-weighted cloud ice fallspeed
      REAL(kind=mg_pr), DIMENSION(kdim) :: uns  ! number-weighted snow fallspeed
      REAL(kind=mg_pr), DIMENSION(kdim) ::  ums ! mass-weighted snow fallspeed
      REAL(kind=mg_pr), DIMENSION(kdim) ::  unr ! number-weighted rain fallspeed
      REAL(kind=mg_pr), DIMENSION(kdim) ::  umr ! mass-weighted rain fallspeed
      REAL(kind=mg_pr) :: unc ! number-weighted cloud droplet fallspeed
      REAL(kind=mg_pr) :: umc ! mass-weighted cloud droplet fallspeed

      REAL(kind=mg_pr), DIMENSION(kdim) :: pracs ! mixing rat tendency due to collection of rain by snow
      REAL(kind=mg_pr), DIMENSION(kdim) :: npracs ! number conc tendency due to collection of rain by snow ! number conc tendency due to collection of rain by snow
      REAL(kind=mg_pr), DIMENSION(kdim) :: mnuccr ! mixing rat tendency due to freezing of rain
      REAL(kind=mg_pr), DIMENSION(kdim) :: nnuccr ! number conc tendency due to freezing of rain
      REAL(kind=mg_pr), DIMENSION(kdim) :: pra ! mixing rat tendnency due to accretion of droplets by rain
      REAL(kind=mg_pr), DIMENSION(kdim) :: npra ! nc tendnency due to accretion of droplets by rain
      REAL(kind=mg_pr), DIMENSION(kdim) :: nragg ! nr tendency due to self-collection of rain
      REAL(kind=mg_pr), DIMENSION(kdim) :: prci ! mixing rat tendency due to autoconversion of cloud ice to snow
      REAL(kind=mg_pr), DIMENSION(kdim) :: nprci ! number conc tendency due to autoconversion of cloud ice to snow
      REAL(kind=mg_pr), DIMENSION(kdim) :: prai ! mixing rat tendency due to accretion of cloud ice by snow
      REAL(kind=mg_pr), DIMENSION(kdim) :: nprai ! number conc tendency due to accretion of cloud ice by snow
      REAL(kind=mg_pr), DIMENSION(kdim) :: mnuccc ! mixing ratio tendency due to freezing of cloud water
      REAL(kind=mg_pr), DIMENSION(kdim) :: nnuccc ! number conc tendency due to freezing of cloud water
      REAL(kind=mg_pr), DIMENSION(kdim) :: nsagg ! ns tendency due to self-aggregation of snow


      REAL(kind=mg_pr), DIMENSION(kdim) :: pre ! rain mixing rat tendency due to evaporation
      REAL(kind=mg_pr), DIMENSION(kdim) :: prds ! snow mixing rat tendency due to sublimation

      REAL(kind=mg_pr), DIMENSION(kdim) ::psacws ! mixing rat tendency due to collection of droplets by snow
      REAL(kind=mg_pr), DIMENSION(kdim) ::npsacws ! number conc tendency due to collection of droplets by snow

! for the one ice class scheme: 
      REAL(kind=mg_pr), DIMENSION(kdim) ::psacws_o ! mixing rat tendency due to collection of droplets by snow
      REAL(kind=mg_pr), DIMENSION(kdim) ::npsacws_o ! number conc tendency due to collection of droplets by snow

      REAL(kind=mg_pr), DIMENSION(kdim) :: prc    ! qc tendency due to autoconversion of cloud droplets
      REAL(kind=mg_pr), DIMENSION(kdim) :: nprc   ! number conc tendency due to autoconversion of cloud droplets
      REAL(kind=mg_pr), DIMENSION(kdim) :: nprc1  ! qr tendency due to autoconversion of cloud droplets

      REAL(kind=mg_pr), DIMENSION(kdim) ::  bergs ! mixing rat tendency due to bergeron process for snow



      REAL(kind=mg_pr), DIMENSION(kdim) ::  lami       ! slope of cloud ice size distr
      REAL(kind=mg_pr), DIMENSION(kdim) ::  n0i        ! intercept of cloud ice size distr


      REAL(kind=mg_pr), DIMENSION(kdim) :: lamr  ! slope of rain size distr
      REAL(kind=mg_pr), DIMENSION(kdim) :: n0r    ! intercept of rain size distr


      REAL(kind=mg_pr), DIMENSION(kdim) :: lams  ! slope of snow size distr
      REAL(kind=mg_pr), DIMENSION(kdim) :: n0s   ! intercept of snow size distr




      REAL(kind=mg_pr) :: lammax  ! maximum allowed slope of size distr
      REAL(kind=mg_pr) :: lammin  ! minimum allowed slope of size distr

      REAL(kind=mg_pr), DIMENSION(kdim) ::  pgam ! spectral width parameter of droplet size distr

      REAL(kind=mg_pr), DIMENSION(kdim) ::  cdist1 ! size distr parameter to calculate droplet freezing



      REAL(kind=mg_pr), DIMENSION(kdim) ::  lamc ! slope of cloud liquid size distr



      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: cldmax ! precip fraction assuming maximum overlap
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: cldm     ! cloud fraction 


      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  t1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  q1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  qc1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  qi1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  nc1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  ni1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  tlat1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  qvlat1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  qctend1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  qitend1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  nctend1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  nitend1

      REAL(kind=mg_pr), DIMENSION(idim) :: prect1
      REAL(kind=mg_pr), DIMENSION(idim) :: preci1

      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: cmei       ! dep/sublimation rate of cloud ice
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: cmel     ! cond/evap rate of cloud liquid
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: cmel_orig, cmei_orig, &
                                                    berg_orig




      REAL(kind=mg_pr) :: qvi ! ice saturation vapor mixing ratio

      REAL(kind=mg_pr) :: qvqvsi ! ICE SATURAION RATIO

      REAL(kind=mg_pr) :: qvl  ! liquid sat mixing ratio  
      REAL(kind=mg_pr) :: qvs ! liquid saturation vapor mixing ratio

      REAL(kind=mg_pr) :: dqsdt ! change of sat vapor mixing ratio with temperature
      REAL(kind=mg_pr) :: dqsidt ! change of ice sat vapor mixing ratio with temperature

      REAL(kind=mg_pr), DIMENSION(idim,kdim)  ::  abi ! correction factor for snow sublimation to account for latent heat
      REAL(kind=mg_pr) :: qclr ! water vapor mixing ratio in clear air

      REAL(kind=mg_pr) :: ab ! correction factor for rain evap to account for latent heat

      REAL(kind=mg_pr) ::  epsi ! 1/ sat relaxation timecale for cloud ice
      REAL(kind=mg_pr) ::  prd ! provisional deposition rate of cloud ice at water sat   
      REAL(kind=mg_pr)  :: epsr ! 1/ sat relaxation timescale for rain
      REAL(kind=mg_pr)  :: epss ! 1/ sat relaxation timescale for snow
      
      real(kind=mg_pr) :: bergtsf   !bergeron timescale to remove all liquid

      REAL(kind=mg_pr) :: esat0, gamma1, qvmax

! hm add 3/19/07, ice nucleation, droplet activation
      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  dum2i ! number conc of ice nuclei available (1/kg)
      REAL(kind=mg_pr) , DIMENSION(idim,kdim)::  dum2l ! number conc of CCN (1/kg)
      REAL(kind=mg_pr) ::   ncmax
      REAL(kind=mg_pr) ::   nimax




! diagnostic rain/snow for output to history
! values are in-precip (local) !!!!

      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nrout ! rain number concentration (1/m3)
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nsout ! snow number concentration (1/m3)



      REAL(kind=mg_pr) :: dc0  ! mean size droplet size distr
      REAL(kind=mg_pr) :: ds0  ! mean size snow size distr (area weighted)
      REAL(kind=mg_pr) :: eci  ! collection efficiency for riming of snow by droplets


      REAL(kind=mg_pr) :: mtime ! factor to account for droplet activation timescale

! variabels to check for RH after rain evap

      REAL(kind=mg_pr) :: esn
      REAL(kind=mg_pr) :: qsn
      REAL(kind=mg_pr) :: qtmp
      REAL(kind=mg_pr) :: ttmp, tc, rhi

      REAL(kind=mg_pr) :: dum, dum1, dum2, dumd 


      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: dumc ! dummy in-cloud qc
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: dumnc ! dummy in-cloud nc
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: dumi ! dummy in-cloud qi
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: dumni ! dummy in-cloud ni
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: dums ! dummy in-cloud snow mixing rat
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: dumns ! dummy in-cloud snow number conc
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: dumr ! dummy in-cloud rain mixing rat
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: dumnr ! dummy in-cloud rain number conc

      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  cldn ! cloud fraction

      REAL(kind=mg_pr) :: qce ! dummy qc for conservation check
      REAL(kind=mg_pr) :: qie ! dummy qi for conservation check
      REAL(kind=mg_pr) :: nce ! dummy nc for conservation check
      REAL(kind=mg_pr) :: nie ! dummy ni for conservation check

      REAL(kind=mg_pr) :: ratio ! parameter for conservation check



! below are parameters for cloud water and cloud ice sedimentation calculations
      REAL(kind=mg_pr), DIMENSION(kdim) :: fr
      REAL(kind=mg_pr), DIMENSION(kdim) :: fnr
      REAL(kind=mg_pr), DIMENSION(kdim) :: fc
      REAL(kind=mg_pr), DIMENSION(kdim) :: fnc
      REAL(kind=mg_pr), DIMENSION(kdim) :: fi
      REAL(kind=mg_pr), DIMENSION(kdim) :: fni
      REAL(kind=mg_pr), DIMENSION(kdim) :: fs
      REAL(kind=mg_pr), DIMENSION(kdim) :: fns
      REAL(kind=mg_pr), DIMENSION(kdim) :: faloutr
      REAL(kind=mg_pr), DIMENSION(kdim) :: faloutnr
      REAL(kind=mg_pr), DIMENSION(kdim) :: faloutc
      REAL(kind=mg_pr), DIMENSION(kdim) :: faloutnc
      REAL(kind=mg_pr), DIMENSION(kdim) :: falouti
      REAL(kind=mg_pr), DIMENSION(kdim) :: faloutni
      REAL(kind=mg_pr), DIMENSION(kdim) :: falouts
      REAL(kind=mg_pr), DIMENSION(kdim) :: faloutns
      REAL(kind=mg_pr) :: faltndr
      REAL(kind=mg_pr) :: faltndnr
      REAL(kind=mg_pr) :: faltndc
      REAL(kind=mg_pr) :: faltndnc
      REAL(kind=mg_pr) :: faltndi
      REAL(kind=mg_pr) :: faltndni
      REAL(kind=mg_pr) :: faltnds
      REAL(kind=mg_pr) :: faltndns
      REAL(kind=mg_pr) :: faltndqie
      REAL(kind=mg_pr) :: faltndqce

          
      REAL(kind=mg_pr) :: rgvm ! max fallspeed for all species

      REAL(kind=mg_pr) :: dumt1, dumt2

      REAL(kind=mg_pr) :: ql_new, qi_new, qn_new, qni_new, qa_new
   
! hm add 3/19/07, new loop variables for sub-step solution
      integer iter, it, ltrue(idim)

      REAL(kind=mg_pr) :: dum3


      REAL(kind=mg_pr) :: qii_new, nii_new,  nii_min,  nii_max 


      logical           ::  do_berg_ex, do_berg1


! v1.4
! new variables for seifert and beheng warm rain scheme
      REAL(kind=mg_pr), DIMENSION(kdim) :: nu
      integer dumii
      REAL(kind=mg_pr) :: dnu(16)

     

      REAL(kind=mg_pr) :: scalef, m1, m2, m3, m4, m5


      INTEGER :: i, k, n ,nstep, kk






!cms++ 6/6/08 temp?
      real(kind=mg_pr), parameter :: d622 = rdgas / rvgas
      real(kind=mg_pr), parameter :: d378 = 1._mg_pr - d622
!cms--

      REAL (kind=mg_pr) :: tmp1, tmp2, tmp3, tmp7, coffi

      REAL (kind=mg_pr) :: qs_t, qs_d , eslt, esit,  dqsi


      REAL (kind=mg_pr) :: q_e




!------------------------------------

!cloud droplets   

      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: pre1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: prds1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: snow2vapor1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: sedi_ice2vapor1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: sedi_liquid2vapor1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: super_saturation_rm1

      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: cmel1     ! cond/evap rate of cloud liquid
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: D_eros_l1   
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: berg1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: qvdep_qi1

      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: prc1

      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: pra1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: mnuccc1

      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: psacws1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: bergs1


!cloud ice
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: psacws_o1

      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: cmei1   
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: D_eros_i1   

      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: prci1
      REAL(kind=mg_pr), DIMENSION(idim,kdim)  :: prai1

!cloud droplet number
      REAL(kind=mg_pr), DIMENSION(kdim) :: nucclim
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nucclim1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: npccn1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nnuccc1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: npsacws1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: npsacws_o1     
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nsubc1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nerosc1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: npra1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nprc11

!cloud ice number
      REAL(kind=mg_pr), DIMENSION(kdim) :: nucclim1i
      REAL(kind=mg_pr), DIMENSION(kdim) :: nucclim2
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nnuccd1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nsubi1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nerosi1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nprci1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nprai1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nucclim1_1
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: nucclim2_1


      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: pracs1
      real(kind=mg_pr), DIMENSION(idim,kdim) :: mnuccr1
      real(kind=mg_pr), DIMENSION(idim,kdim) :: relhum

     
      REAL(kind=mg_pr), DIMENSION(idim,kdim) :: qs


      REAL :: tmp2s, qa_t

      REAL :: NACNT 
    
      REAL(kind=mg_pr), DIMENSION(idim,kdim) ::  sum_freeze, sum_rime, &
                                                 sum_bergs, sum_ice_adj, &
                                                 sum_berg, sum_cond, &
                                                 sum_freeze2
      real :: qldt_sum
 



!---------------------

!some sanity checking (more needed!!)

         IF ( .NOT. Nml%do_pdf_clouds ) THEN
           IF (auto_conv_ice_choice .EQ. 2 ) THEN
             call error_mesg ( 'morrison_gettelman_microp', &
              'ERROR auto_conv_ice_choice =2 not compatible w. Tiedtke&
                  & param. assumption of in-cloud RH=1 ', FATAL)
           END IF
           IF ( qv_on_qi ) THEN
             call error_mesg ( 'morrison_gettelman_microp', &
               'ERROR qv_on_qi = .true. not compatible w. Tiedtke &
                    &param. assumption of in-cloud RH=1 ', FATAL)
           END IF
         END IF

!---------------------------------------------------------------------
         DO k=1,kdim
           DO i= 1,idim
             qv(i,k) = qv_in(i,k)
             qc(i,k) =   qc_in(i,k)   
             qi(i,k) =   qi_in(i,k)  
             nc(i,k) =   nc_in(i,k)   
             ni(i,k) =   ni_in(i,k)
             tn(i,k) = t_in(i,k)
             cldn(i,k) = cldn_in(i,k) 
           END DO
         END DO

!----------------------------------------------------

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!++ag assign variable deltat for sub-stepping...
         deltat = deltatin

! parameters for scheme

         omsm = 0.99999_mg_pr
         mincld = 0.0001_mg_pr
 
! initialize multi-level fields
         do k=1,kdim
           do i=1,idim
             q(i,k) = qv(i,k)
             t(i,k) = tn(i,k)
           end do
         end do
  

! initialize time-varying parameters

         do k=1,kdim
           do i=1,idim
             rho(i,k) = pfull(i,k)/(rdgas*t(i,k))
             dv(i,k) = 8.794E-5_mg_pr*t(i,k)**1.81_mg_pr/pfull(i,k)
             mu(i,k) = 1.496E-6_mg_pr*t(i,k)**1.5_mg_pr/ &
                                            (t(i,k) + 120._mg_pr)
             sc(i,k) = mu(i,k)/(rho(i,k)*dv(i,k))
             kap(i,k) = 1.414e3_mg_pr*1.496e-6_mg_pr*t(i,k)**1.5_mg_pr/  &
                                                     (t(i,k) + 120._mg_pr) 

! air density adjustment for fallspeed parameters
! hm added 11/18/06, add air density correction factor to the
! power of 0.54 following Heymsfield and Bansemer 2006

!            if (rho_factor_in_max_vt) then
               rhof(i,k)=(rhosu/rho(i,k))**0.54_mg_pr
!            else
!              rhof(i,k) = 1.0
!            endif
             rhof(i,k) = MIN (rhof(i,k), max_rho_factor_in_vt)

!            arn(i,k)=ar*(rhosu/rho(i,k))**0.54_mg_pr
!            asn(i,k)=as*(rhosu/rho(i,k))**0.54_mg_pr
!            acn(i,k)=ac*(rhosu/rho(i,k))**0.54_mg_pr
!            ain(i,k)=ai*(rhosu/rho(i,k))**0.54_mg_pr
             arn(i,k)=ar*rhof(i,k)                        
             asn(i,k)=as*rhof(i,k)                        
             acn(i,k)=ac*rhof(i,k)                        
             ain(i,k)=ai*rhof(i,k)                        
 
! keep dz positive (define as layer k-1 - layer k)

             dz(i,k)= pdel(i,k)/(rho(i,k)*grav)
           end do
         end do

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate condensation based on cloud fraction, t,q tendencies
! 

         do k=1,kdim
           do i=1,idim

! store original variables for sub-stepping

             t1(i,k) = t(i,k)
             q1(i,k) = q(i,k)
             qc1(i,k) = qc(i,k)
             qi1(i,k) = qi(i,k)
             nc1(i,k) = nc(i,k)
             ni1(i,k) = ni(i,k)
!cms++
! initialize tendencies to zero

             tlat1(i,k) = 0._mg_pr
             qvlat1(i,k) = 0._mg_pr
             qctend1(i,k) = 0._mg_pr
             qitend1(i,k) = 0._mg_pr
             nctend1(i,k) = 0._mg_pr
             nitend1(i,k) = 0._mg_pr

             tlat(i,k) = 0._mg_pr
             qvlat(i,k) = 0._mg_pr
             qctend(i,k) = 0._mg_pr
             qitend(i,k) = 0._mg_pr
             qnitend(i,k) = 0._mg_pr
             qrtend(i,k) = 0._mg_pr
             nctend(i,k) = 0._mg_pr
             nitend(i,k) = 0._mg_pr
             nrtend(i,k) = 0._mg_pr
             nstend(i,k) = 0._mg_pr
             prect(i) = 0._mg_pr
             preci(i) = 0._mg_pr
             qniic(i,k) = 0._mg_pr
             qric(i,k) = 0._mg_pr
             nsic(i,k) = 0._mg_pr
             nric(i,k) = 0._mg_pr
 
! hm add 9/5/07
             rainrt(i,k) = 0._mg_pr
!cms for calc. rain3d, snow3d

!cms--
! initialize precip output

             qrout(i,k) = 0._mg_pr
             qsout(i,k) = 0._mg_pr
             nrout(i,k) = 0._mg_pr
             nsout(i,k) = 0._mg_pr
 
!  initialize bergeron fraction arrays
             sum_freeze(i,k) = 0._mg_pr
             sum_freeze2(i,k) = 0._mg_pr
             sum_rime  (i,k) = 0._mg_pr
             sum_berg  (i,k) = 0._mg_pr
             sum_ice_adj(i,k) = 0._mg_pr
             sum_bergs (i,k) = 0._mg_pr
             sum_cond  (i,k) = 0._mg_pr

             rainrt1(i,k) = 0._mg_pr

!diag++
!droplets
             pre1(i,k)                 = 0._mg_pr
             prds1(i,k)                = 0._mg_pr
             snow2vapor1(i,k)          = 0._mg_pr
             snow2vapor(k)             = 0._mg_pr
             sedi_ice2vapor1(i,k)      = 0._mg_pr
             sedi_liquid2vapor1(i,k)   = 0._mg_pr
             super_saturation_rm1(i,k) = 0._mg_pr

             cmel1(i,k) = 0._mg_pr
             D_eros_l1(i,k) = 0._mg_pr 
             berg1(i,k) = 0._mg_pr
             qvdep_qi1(i,k) = 0._mg_pr
             prc1(i,k) = 0._mg_pr
             pra1(i,k) = 0._mg_pr
             mnuccc1(i,k) = 0._mg_pr
             psacws1(i,k) = 0._mg_pr
             bergs1(i,k) = 0._mg_pr
!ice
             psacws_o1(i,k) = 0._mg_pr
             cmei1(i,k) = 0._mg_pr
             D_eros_i1(i,k) = 0._mg_pr   
             prci1(i,k) = 0._mg_pr 
             prai1(i,k) = 0._mg_pr 

!droplet number
             nucclim1(i,k) = 0._mg_pr 
             npccn1(i,k) = 0._mg_pr 
             nnuccc1(i,k) = 0._mg_pr 
             npsacws1(i,k) = 0._mg_pr 
             npsacws_o1(i,k) = 0._mg_pr  
             nsubc1(i,k) = 0._mg_pr 
             nerosc1(i,k) = 0._mg_pr 
             npra1(i,k) = 0._mg_pr 
             nprc11(i,k) = 0._mg_pr 

!ice number
             nnuccd1(i,k) = 0._mg_pr 
             nsubi1(i,k) = 0._mg_pr 
             nerosi1(i,k) = 0._mg_pr 
             nprci1(i,k) = 0._mg_pr 
             nprai1(i,k) = 0._mg_pr 
             nucclim1_1(i,k) = 0._mg_pr 
             nucclim2_1(i,k) = 0._mg_pr 
             pracs1(i,k)    = 0._mg_pr
             mnuccr1(i,k)  = 0._mg_pr
           end do 
         end do

! get cloud fraction, check for minimum
         do k=1,kdim
           do i=1,idim
             cldmax(i,k) = mincld
             cldm(i,k) = max(cldn(i,k), mincld)
           end do
         end do

! initialize avg precip rate
         do i=1,idim
           prect1(i)=0._mg_pr
           preci1(i)=0._mg_pr
         end do

         do k=1,kdim
           do i=1,idim
             esl(i,k) = polysvp_l(t(i,k))
             esi(i,k) = polysvp_i(t(i,k))

! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING
             IF (esi(i,k) .GT. esl(i,k)) esi(i,k) = esl(i,k)
           end do
         end do

         IF ( meyers_test ) THEN
! over-write crystal1
! this is Steve Klein's version, Hugh's WRF version is different ...
           do k=1,kdim
             do i=1,idim
               crystal1(i,k)  = 1000._mg_pr*exp((12.96_mg_pr*    &
                                    0.0125_mg_pr*(tfreeze - T(i,k))) - &
                                                             0.639_mg_pr)
             end do
           end do   
         END IF

         do k=1,kdim
           do i=1,idim
             prd = 0._mg_pr
             if (t(i,k) .lt. tfreeze - 5._mg_pr) then
               IF (Nml%do_ice_nucl_wpdf) THEN

!use number calculated outside this code
                 IF (total_activation) THEN
                   dum = crystal1(i,k)*cldm(i,k)/rho(i,k)
                 ELSE IF (dqa_activation) THEN
                   dum =  crystal1(i,k)/rho(i,k)
                 END IF
               ELSE 
! cooper
                 dum = min(0.005_mg_pr*exp(0.304_mg_pr*   &
                           (tfreeze - t(i,k)))*1000._mg_pr, 208.9e3_mg_pr)
                 dum = dum/rho(i,k)
               END IF 

               IF (total_activation) THEN
                 if (delta_cf(i,k) .gt. 0._mg_pr ) then
                   dumnnuc = (dum - ni(i,k)/cldm(i,k))/deltat*cldm(i,k)
                 else
                   dumnnuc = 0._mg_pr
                 end if
               ELSE  IF (dqa_activation) THEN
                 dumnnuc = max(delta_cf(i,k), 0.)*dum/deltat
               END IF
               dumnnuc=max(dumnnuc,0._mg_pr)

  ! get provisional ni and qi after nucleation in order to calculate
  ! Bergeron process below
               ninew = ni(i,k) + dumnnuc*deltat
               if (tiedtke_macrophysics .or. dqa_activation) then
                 qinew = qi(i,k)
               else
                 qinew =qi(i,k) + dumnnuc*deltat*mi0
               endif
             else 
               ninew = ni(i,k)
               qinew = qi(i,k)
             end if

! NOTE: this approach has been replaced by NCAR's, so the following 
!       comments are no longer valid
! for condensation
! for T < tmin_fice, assume all new condensate is ice
! for T > tmin_fice, put all new condensate into liquid temporarily,
! then calculate transfer to ice through Bergeron process
! note: this approach assumes that ice/liquid is mixed throughout the 
! cloudy portion of the grid cell
! END NOTE

! make sure to initialize bergeron process to zero
             berg(i,k) = 0._mg_pr

             qvdep_qi(i,k) = 0._mg_pr

! cmel and cmei are liquid and ice condensate as clacuklated in nc_cond.F90
! and input here
             cmel(i,k) = dqcdt(i,k)
             cmei(i,k) = dqidt(i,k)
             dum2 = dqcdt(i,k) + dqidt(i,k)

! get in-cloud qi and ni after nucleation
!RSH: Note loop may be unnecessary since icldm >= mincld, which is > 0.
             if (cldm(i,k) .gt. 0._mg_pr) then
               qiic(i,k)=qinew/cldm(i,k)
               niic(i,k)=ninew/cldm(i,k)
             else
               qiic(i,k)=0._mg_pr
               niic(i,k)=0._mg_pr
             endif

             IF (.NOT. limit_berg) THEN
               if (dum2 .ge. 0._mg_pr .and.     &
                      (qc(i,k) + dqcdt(i,k)*deltat > qsmall) ) then
                 do_berg1 = .true.
               else
                 do_berg1 = .false.
               end if
             ELSE
               if (dum2 .ge. 0._mg_pr .and. qinew .gt. berg_lim )  then
                 do_berg1 = .true.
               else
                 do_berg1 = .false.
               end if
             END If

             if (do_berg1) then
               if (t(i,k) .lt. tfreeze ) then
                 if (qi(i,k) > qsmall) then
! calculate Bergeron process

                   bergtsf = 0._mg_pr ! bergeron time scale 
                                      ! (fraction of timestep)
                   qvi = 0.622_mg_pr*esi(i,k)/(pfull(i,k) - d378*esi(i,k))
                   qvl = 0.622_mg_pr*esl(i,k)/(pfull(i,k) - d378*esl(i,k))
                   dqsidt =  xxls*qvi/(rv*t(i,k)**2)
                   abi(i,k) = 1._mg_pr + dqsidt*xxls/cpp

! get ice size distribution parameters
!adding this 10/1 0421AM
                   if (qiic(i,k).ge.qsmall) then
                     lami(k) = (gamma_mg(1._mg_pr + di_mg)*ci_mg* &
                                niic(i,k)/qiic(i,k))**(1._mg_pr/di_mg)
                     n0i(k) = niic(i,k)*lami(k)

! check for slope
                     lammax = 1._mg_pr/min_diam_ice
                     lammin = 1._mg_pr/(2._mg_pr*dcs)

! adjust vars
                     if (lami(k) .lt. lammin) then
                       lami(k) = lammin
                       n0i(k) = lami(k)**(di_mg + 1._mg_pr)*qiic(i,k)/  &
                                      (ci_mg*gamma_mg(1._mg_pr + di_mg))
                     else if (lami(k) .gt. lammax) then
                       lami(k) = lammax
                       n0i(k) = lami(k)**(di_mg + 1._mg_pr)*qiic(i,k)/  &
                                        (ci_mg*gamma_mg(1._mg_pr + di_mg))
                     end if
 
                     epsi = 2._mg_pr*pi*n0i(k)*rho(i,k)*Dv(i,k)/  &
                                                        (lami(k)*lami(k))
                     if (qc(i,k) + dqcdt(i,k)*deltat .gt. qsmall) then
                       prd = epsi*(qvl - qvi)/abi(i,k)
                     else
                       prd = 0._mg_pr
                     end if

! multiply by cloud fraction
                     prd = prd*cldm(i,k)

                     berg(i,k) = max(0._mg_pr, prd)
                   end if

                   if (berg(i,k) .gt. 0._mg_pr) then
                     bergtsf = max(0._mg_pr,    &
                                 ((dqcdt(i,k) + qc(i,k)/deltat)/berg(i,k)))
                     if(bergtsf .lt. 1._mg_pr) berg(i,k) = &
                          max(0._mg_pr, dqcdt(i,k) + qc(i,k)/deltat)
                   endif
!TEST IF NEEDED: 9/27 0133Z
                   if (t(i,k) < tfreeze  - 40._mg_pr) then
                     berg(i,k) = 0._mg_pr
                   endif
!cms++
! vapor deposition onto cloud ice 
                   if (qv_on_qi) then
                     dqsi = MAX(qv(i,k) - qvi, 0._mg_pr)
                     qvdep_qi(i,k) = max(min(prd - cmei(i,k) - berg(i,k), &
                                             dqsi/deltat*omsm), 0._mg_pr)
                   end if
!cms--
                 endif
               end if  ! t < 273.15
             endif ! (do_berg1)

! limit cmel,cmei due for roundoff error
             cmel(i,k) = cmel(i,k)*omsm
             cmei(i,k) = cmei(i,k)*omsm

! define activated ice nuclei
             if (t(i,k) <  tfreeze - 5.0_mg_pr) then

! use number calculated outside this code
               IF (Nml%do_ice_nucl_wpdf) THEN
                 IF (total_activation) THEN
                   if (delta_cf(i,k) .gt. 0._mg_pr) then 
                     dum2i(i,k) = crystal1(i,k)*cldm(i,k)/rho(i,k)
                   else
                     dum2i(i,k) = 0._mg_pr
                   end if
                 ELSEIF (dqa_activation) THEN
                   dum2i(i,k) = crystal1(i,k)/rho(i,k)
                 END IF
               ELSE
! cooper
                 dum2i(i,k) =  min(0.005_mg_pr*exp(0.304_mg_pr*  &
                            (tfreeze - t(i,k)))*1000._mg_pr, 208.9e3_mg_pr)
               END IF
             else
               dum2i(i,k) = 0.0_mg_pr
             endif
           end do
         end do



!pdf clouds option -- validity ??
!re-calculate cloud fraction
         IF (Nml%do_pdf_clouds .AND.   &
            (Nml%super_ice_opt .EQ. 1 .OR. Nml%super_ice_opt .EQ. 2)) THEN
           IF (Nml%super_ice_opt .EQ. 1) THEN
             DO k=1,kdim
               DO i=1,idim
                 ttmp = t(i,k) 
                 IF (ttmp .LT. tfreeze - 40._mg_pr .OR.    &
                    (ttmp .LE. tfreeze .AND. qc(i,k) +    &
                         (cmel(i,k) - berg(i,k) )/deltat  &
                                          .LT. 3._mg_pr*Nml%qmin)) THEN 
                   eslt = polysvp_i(ttmp)
                 ELSE
                   eslt = polysvp_l(ttmp)
                 END IF
                 qs_d = pfull(i,k) - d378*eslt
                 qs_d = max(qs_d, eslt)
                 qs(i,k) = d622*eslt/qs_d 
               END DO
             END DO
           END IF

           IF (Nml%super_ice_opt .EQ. 2) THEN
             DO k=1,kdim
               DO i=1,idim
                 ttmp = t(i,k) 
                 IF (ttmp .LT. tfreeze - 40._mg_pr .OR.   &
                    (ttmp .LE. tfreeze .AND. qc(i,k) +   &
                         (cmel(i,k) - berg(i,k) )/deltat   &
                                       .LT. 3._mg_pr*Nml%qmin) ) THEN 
                   eslt = polysvp_i(ttmp)
                   tc = ttmp - tfreeze
!!!                rhi=MIN( max_super_ice, 0.000195*tc**2+0.00266*tc+1.005)
                   rhi = 0.000195_mg_pr*tc**2 + 0.00266_mg_pr*tc +   &
                                                               1.005_mg_pr
                 ELSE
                   eslt = polysvp_l(ttmp)
                   rhi = 1._mg_pr
                 END IF
                 qs_d = pfull(i,k) - d378*eslt
                 qs_d = max(qs_d, eslt)
                 qs(i,k)= rhi*d622*eslt/qs_d 
               END DO
             END DO
           END IF

           qtot = qv_in + qc_in + qi_in 

           IF ( Nml%pdf_org )    &
                 call error_mesg ( 'morrison_gettelman_microp', &
                                        'ERROR 1 simple_pdf ', FATAL)

           CALL  simple_pdf (j, idim, jdim, kdim, Nml%qmin, qa0, &
                             qtot, qs,  gamma, Nml%qthalfwidth,  &
                             Nml%betaP, inv_dtcloud,                 &
!inout
                             SA_0,                            &
!diag
                             n_diag_4d, diag_4d, diag_id, diag_pt,  &
                             SA, qa_upd)


           do k=1,kdim
             do i=1,idim
               cldm(i,k) = max(qa_upd(i,k), mincld)
             end do
           end do 
         END IF  

!! initialize sub-step precip flux variables
         do i=1,idim
!! flux is zero at top interface, so these should stay as 0.
           do k=1,kdim
 
! initialize normal and sub-step precip flux variables
             atotrt1(i,k)=0._mg_pr
             asnowrt1(i,k)=0._mg_pr
           end do 
         end do 

!! initialize final precip flux variables.
         do i=1,idim
!! flux is zero at top interface, so these should stay as 0.
           do k=1,kdim
! initialize normal and sub-step precip flux variables
             atotrt(i,k)=0._mg_pr
             asnowrt(i,k)=0._mg_pr
           end do 
         end do 

         do i=1,idim
           ltrue(i) = 0
           do k=1,kdim
! hm add 3/19/07 skip microphysical calculations if no cloud water

             if (qc(i,k).ge.qsmall .or. qi(i,k).ge.qsmall .or.   &
                 cmel(i,k).ge.qsmall .or. cmei(i,k).ge.qsmall .or.   &
                 qvdep_qi(i,k).ge.qsmall) ltrue(i)=1

!cms also skip if total water amount is negative anywhere within the column
             if (qc(i,k) + qi(i,k) + qv(i,k) .lt. -1.e-9_mg_pr .OR.  &
                 qv(i,k)  .lt. -1.e-9_mg_pr) then
               ltrue(i) = 0
               nrefuse = nrefuse + 1
             end if
           end do
         end do


! assign number of sub-steps to iter
! use 2 sub-steps, following tests described in 
! Morrison and Gettelman, 2007, J. Clim.
         iter = 2

! get sub-step time step
         deltat = deltat/real(iter)

!!!! skip calculations if no cloud water

         do i=1,idim

           if (ltrue(i) .eq. 0) then
             do k=1,kdim
               tlat(i,k)=0._mg_pr
               qvlat(i,k)=0._mg_pr
               qctend(i,k)=0._mg_pr
               qitend(i,k)=0._mg_pr
               qnitend(i,k)=0._mg_pr
               qrtend(i,k)=0._mg_pr
               nctend(i,k)=0._mg_pr
               nitend(i,k)=0._mg_pr
               nrtend(i,k)=0._mg_pr
               nstend(i,k)=0._mg_pr
               prect(i)=0._mg_pr
               preci(i)=0._mg_pr
               qniic(i,k)=0._mg_pr
               qric(i,k)=0._mg_pr
               nsic(i,k)=0._mg_pr
               nric(i,k)=0._mg_pr

! hm add 9/5/07
               rainrt(i,k)=0._mg_pr
             end do
             goto 300
           end if

           cmel_orig(i,:) = cmel(i,:)
           cmei_orig(i,:) = cmei(i,:)
           berg_orig(i,:) = berg(i,:)

!!!!!!!!! begin sub-step!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.........................................................................

           do it=1,iter
             do k=1,kdim

! initialize sub-step microphysical tendencies
               tlat(i,k)=0._mg_pr
               qvlat(i,k)=0._mg_pr
               qctend(i,k)=0._mg_pr
               qitend(i,k)=0._mg_pr
               qnitend(i,k)=0._mg_pr
               qrtend(i,k)=0._mg_pr
               nctend(i,k)=0._mg_pr
               nitend(i,k)=0._mg_pr
               nrtend(i,k)=0._mg_pr
               nstend(i,k)=0._mg_pr

! initialize diagnostic precipitation to zero

               qniic(i,k)=0._mg_pr
               qric(i,k)=0._mg_pr
               nsic(i,k)=0._mg_pr
               nric(i,k)=0._mg_pr

               nerosc(i,k) = nerosc4(i,k)
               nerosi(i,k) = nerosi4(i,k)
               D_eros_l(i,k) = D_eros_l4(i,k)
               D_eros_i(i,k) = D_eros_i4(i,k)

               cmel(i,k) = cmel_orig(i,k)
               cmei(i,k) = cmei_orig(i,k)
               berg(i,k) = berg_orig(i,k)

! hm, add 9/5/07
 
               rainrt(i,k)=0._mg_pr

             end do

! initialize vertically-integrated rain and snow tendencies

             qrtot = 0._mg_pr
             nrtot = 0._mg_pr
             qstot = 0._mg_pr
             nstot = 0._mg_pr

! initialize precip at surface

             prect(i)=0._mg_pr
             preci(i)=0._mg_pr

             do k=1,kdim

! set cwml and cwmi to current qc and qi
               cwml(i,k) = qc(i,k)
               cwmi(i,k) = qi(i,k)

! initialize precip fallspeeds to zero

               ums(k) = 0._mg_pr 
               uns(k) = 0._mg_pr 
               umr(k) = 0._mg_pr 
               unr(k)= 0._mg_pr

! calculate precip fraction based on maximum overlap assumption

               if (k .eq. 1) then
                 cldmax(i,k) = cldm(i,k)
               else

! hm add sep 6, 2006, if rain or snow mix ratio is smaller than
! threshold, then set cldmax to cloud fraction at current level
                 if (qric(i,k-1) .ge. qsmall .or.   &
                        qniic(i,k-1) .ge. qsmall) then
                   cldmax(i,k) = max(cldmax(i,k-1), cldm(i,k))
                 else
                   cldmax(i,k) = cldm(i,k)
                 end if
               end if

! should this be behind the next block (scaling)?
! decrease in number concentration due to sublimation/evap
! divide by cloud fraction to get in-cloud decrease
! don't reduce Nc due to bergeron process ?????
 
               nsubi(k) = 0._mg_pr
               nsubc(k) = 0._mg_pr
 
               if (.not. tiedtke_macrophysics) then
                 if (cmei(i,k) < 0._mg_pr .and.     &
                     qi(i,k) > qsmall .and.    &
                     cldm(i,k) > mincld) then
                   nsubi(k) = cmei(i,k)/qi(i,k)*ni(i,k)/cldm(i,k)
                 else
                   nsubi(k) = 0._mg_pr
                 end if
                 if (cmel(i,k) < 0._mg_pr  .AND.    &
                     qc(i,k) .ge. qsmall  .and.    &
                     cldm(i,k) > mincld)      then
                   nsubc(k) = cmel(i,k)/qc(i,k)*nc(i,k)/cldm(i,k)
                 else
                   nsubc(k) = 0._mg_pr
                 end if
               end if
       
!c.......................................................................
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! obtain in-cloud values of cloud water/ice mixing ratios and number 
! concentrations for microphysical process calculations
! units are kg/kg for mixing ratio, 1/kg for number conc

! limit in-cloud values to 0.005 kg/kg


               qcic(i,k) = min(cwml(i,k)/cldm(i,k), in_cloud_limit)
               qiic(i,k) = min(cwmi(i,k)/cldm(i,k), in_cloud_limit)
               ncic(i,k) = max(nc(i,k)/cldm(i,k), 0._mg_pr)
               niic(i,k) = max(ni(i,k)/cldm(i,k), 0._mg_pr)

               if (qc(i,k) +    &
                      (cmel(i,k) + D_eros_l(i,k) - berg(i,k))*deltat   &
                                                          .lt. qsmall) then
                 qcic(i,k) = 0._mg_pr
                 ncic(i,k) = 0._mg_pr

                 if (qc(i,k) +   &
                      (cmel(i,k) + D_eros_l(i,k) - berg(i,k))*deltat   &
                                                      .lt. 0._mg_pr) then
                   if (cmel(i,k) .lt. 0._mg_pr) then

!++ first only scale cmel, d_eros
                     dum = -cmel(i,k) - D_eros_l(i,k)
                     if (dum .gt. 1.e-30_mg_pr) then
                       dum3 = qc(i,k)/deltat/dum*omsm
                     else
                       dum3 = 0._mg_pr 
                     end if
                     cmel(i,k) = dum3*cmel(i,k)
                     D_eros_l(i,k) = dum3*D_eros_l(i,k)

!--
                     dum = -cmel(i,k) - D_eros_l(i,k) + berg(i,k)
                     if (dum .gt. 1.e-30_mg_pr) then
                       dum3 = qc(i,k)/deltat/dum*omsm
                     else
                       dum3 = 0._mg_pr 
                     end if
                     cmel(i,k) = dum3*cmel(i,k)
                     D_eros_l(i,k) = dum3*D_eros_l(i,k)
                     berg(i,k) = dum3*berg(i,k)
                   else
                     dum = -D_eros_l(i,k) + berg(i,k)
!                    berg(i,k)=qc(i,k)/deltat*omsm
                     if (dum .gt. 1.e-30_mg_pr) then
                       dum3 = (qc(i,k)/deltat + cmel(i,k))/dum*omsm
                     else
                       dum3 = 0._mg_pr 
                     end if 
                     D_eros_l(i,k) = D_eros_l(i,k)*dum3
                     berg(i,k) = berg(i,k)*dum3
                   end if
                 end if
               end if

               if (qi(i,k) + (cmei(i,k) + D_eros_i(i,k) + berg(i,k) +  &
                                  qvdep_qi(i,k))*deltat .lt. qsmall) then
                 qiic(i,k) = 0._mg_pr
                 niic(i,k) = 0._mg_pr
                 if (qi(i,k) + (cmei(i,k) + berg(i,k) + D_eros_i(i,k) +  &
                                 qvdep_qi(i,k))*deltat .lt. 0._mg_pr) then
                   if (cmei(i,k) .lt. 0._mg_pr) then
                     dum = - cmei(i,k) - D_eros_i(i,k) 
                     if (dum .gt. 1.e-30_mg_pr) then
                       dum3 = (qi(i,k)/deltat + berg(i,k) +   &
                                                  qvdep_qi(i,k))/dum*omsm
                     else
                       dum3 = 0._mg_pr 
                     end if
                     cmei(i,k) = dum3*cmei(i,k)
                     D_eros_i(i,k) = dum3*D_eros_i(i,k)
                   else
                     dum = - D_eros_i(i,k)
                     if (dum .gt. 1.e-30_mg_pr) then
                       dum3 = (qi(i,k)/deltat + cmei(i,k) + berg(i,k) + &
                                                   qvdep_qi(i,k))/dum*omsm
                     else
                       dum3 = 0._mg_pr
                     end if
                     D_eros_i(i,k) = dum3*D_eros_i(i,k)
                   end if
                 end if
               end if

! add to cme output

!cms           cmeout1(i,k) = cmeout1(i,k)+cmel(i,k)+cmei(i,k)

            
               if (qiic(i,k) .ge. qsmall .and.    &
                   t(i,k) .lt. tfreeze - 5._mg_pr) then

! if NCAI > 0. then set numice = ncai (as before)
! note: this is gridbox averaged
 
                 if (total_activation) then
                   nnuccd(k) = (dum2i(i,k) - ni(i,k)/cldm(i,k))/deltat*  &
                                                                  cldm(i,k)
                   nnuccd(k) = max(nnuccd(k), 0._mg_pr)
                 else if (dqa_activation ) then
                   nnuccd(k) = max(delta_cf(i,k), 0.)*dum2i(i,k)/deltatin
                 endif
                 nimax = dum2i(i,k)*cldm(i,k)
               else
                 nnuccd(k) = 0._mg_pr
                 nimax = 0._mg_pr
                 mnuccd(k) = 0._mg_pr
               end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! droplet activation
! calculate potential for droplet activation if cloud water is present

               if (qcic(i,k) .ge. qsmall) then
                 IF (total_activation) THEN
                   if (delta_cf(i,k) .gt. 0._mg_pr) then
                     dum2l(i,k) = drop2(i,k)*cldm(i,k)
                   else 
                     dum2l(i,k) = 0._mg_pr
                   end if

! assume aerosols already activated are equal to number of existing 
! droplets for simplicity
! multiply by cloud fraction to obtain grid-average tendency
                   npccn(k) = (dum2l(i,k) - nc(i,k)/cldm(i,k))/deltat*  &
                                                                 cldm(i,k)

! make sure number activated > 0
                   npccn(k) = max(0._mg_pr, npccn(k))

                 ELSE IF (dqa_activation) THEN
! delta_cf:  A_dt * (1.-qabar)   where A_dt = A*dt , A source rate
! Eq. 7 of Yi's 2007 paper
! dum2l has already been multiplied by 1.e6/airdens(i,k)
                   dum2l(i,k) = drop2(i,k) 
                   npccn(k) = max(delta_cf(i,k), 0.)*dum2l(i,k) /deltatin
                 END IF
                 ncmax = dum2l(i,k)*cldm(i,k)
               else
                 npccn(k) = 0._mg_pr
                 ncmax = 0._mg_pr
               end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! get size distribution parameters based on in-cloud cloud water/ice 
! these calculations also ensure consistency between number and 
! mixing ratio
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!......................................................................
! cloud ice

               if (qiic(i,k).ge.qsmall) then

! add upper limit to in-cloud number concentration to prevent
! numerical error
                 niic(i,k) = min(niic(i,k), qiic(i,k)*1.e20_mg_pr)
                 lami(k) = (gamma_mg(1._mg_pr + di_mg)*ci_mg* &
                                   niic(i,k)/qiic(i,k))**(1._mg_pr/di_mg)
                 n0i(k) = niic(i,k)*lami(k)

! check for slope
                 lammax = 1._mg_pr/min_diam_ice
                 lammin = 1._mg_pr/(2._mg_pr*dcs)

! adjust vars
                 if (lami(k) .lt. lammin) then
                   lami(k) = lammin
                   n0i(k) = lami(k)**(di_mg + 1._mg_pr)*qiic(i,k)/  &
                                        (ci_mg*gamma_mg(1._mg_pr + di_mg))
                   niic(i,k) = n0i(k)/lami(k)
                 else if (lami(k) .gt. lammax) then
                   lami(k) = lammax
                   n0i(k) = lami(k)**(di_mg + 1._mg_pr)*qiic(i,k)/    &
                                        (ci_mg*gamma_mg(1._mg_pr + di_mg))
                   niic(i,k) = n0i(k)/lami(k)
                 end if
               else
                 lami(k) = 0._mg_pr
                 n0i(k) = 0._mg_pr
               end if

               if (qcic(i,k).ge.qsmall) then

! add upper limit to in-cloud number concentration to prevent   
! numerical error
                 ncic(i,k) = min(ncic(i,k), qcic(i,k)*1.e20_mg_pr)

! get pgam from fit to observations of martin et al. 1994

!RSH BUGFIX email of 6/8/10
!                pgam(k) = 0.0005714_mg_pr*(ncic(i,k)/1.e6_mg_pr/   &
!                                                 rho(i,k)) + 0.2714_mg_pr
                 pgam(k) = 0.0005714_mg_pr*(ncic(i,k)/1.e6_mg_pr*  &
                                                  rho(i,k)) + 0.2714_mg_pr
                 pgam(k) = 1._mg_pr/(pgam(k)**2) - 1._mg_pr
                 pgam(k) = max(pgam(k), 2._mg_pr)
                 pgam(k) = min(pgam(k), 15._mg_pr)

! calculate lamc
                 lamc(k) = (pi/6._mg_pr*rhow*ncic(i,k)*  &
                             gamma_mg(pgam(k) + 4._mg_pr)/(qcic(i,k)*    &
                       gamma_mg(pgam(k) + 1._mg_pr)))**(1._mg_pr/3._mg_pr)

! lammin, 40 micron diameter max mean size
                 lammin = (pgam(k) + 1._mg_pr)/max_diam_drop
                 lammax = (pgam(k) + 1._mg_pr)/min_diam_drop

                 if (lamc(k) .lt. lammin) then
                   lamc(k) = lammin
                   ncic(i,k) = 6._mg_pr*lamc(k)**3*qcic(i,k)* &
                                     gamma_mg(pgam(k) + 1._mg_pr)/ &
                                    (pi*rhow*gamma_mg(pgam(k) + 4._mg_pr))
                 else if (lamc(k) .gt. lammax) then
                   lamc(k) = lammax
                   ncic(i,k) = 6._mg_pr*lamc(k)**3*qcic(i,k)* &
                                       gamma_mg(pgam(k) + 1._mg_pr)/ &
                                     (pi*rhow*gamma_mg(pgam(k) + 4._mg_pr))
                 end if

! parameter to calculate droplet freezing
                 cdist1(k) = ncic(i,k)/gamma_mg(pgam(k) + 1._mg_pr) 
               else
                 lamc(k) = 0._mg_pr
                 cdist1(k) = 0._mg_pr
               end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin micropysical process calculations 
!.................................................................
! autoconversion of cloud liquid water to rain
! formula from Khrouditnov and Kogan (2000)
! minimum qc of 1 x 10^-8 prevents floating point error

               if (qcic(i,k) .ge. 1.e-8_mg_pr) then

! nprc is increase in rain number conc due to autoconversion
! nprc1 is decrease in cloud droplet conc due to autoconversion

! assume exponential sub-grid distribution of qc, resulting in additional
! factor related to qcvar below

!                prc(k) = gamma_mg(qcvar+2.47_mg_pr)/(gamma_mg(qcvar)*qcvar**2.47_mg_pr)*1350._mg_pr*qcic(i,k)**2.47_mg_pr* &
                 prc(k) = sfac1*1350._mg_pr*qcic(i,k)**2.47_mg_pr* &
                            (ncic(i,k)/1.e6_mg_pr*rho(i,k))**(-1.79_mg_pr)
                 nprc(k) = prc(k)/(4._mg_pr/3._mg_pr*pi*rhow*   &
                                                      (25.e-6_mg_pr)**3)
                 nprc1(k) = prc(k)/(qcic(i,k)/ncic(i,k))
               else
                 prc(k)=0._mg_pr
                 nprc(k)=0._mg_pr
                 nprc1(k)=0._mg_pr
               end if

! add autoconversion to precip from above to get provisional rain mixing 
! ratio and number concentration (qric and nric)
! 0.45 m/s is fallspeed of new rain drop (80 micron diameter)

               dum = 0.45_mg_pr
               dum1 = 0.45_mg_pr

! hm modify 6/12
               if (k .eq. 1) then
                 qric(i,k) = prc(k)*cldm(i,k)*dz(i,k)/cldmax(i,k)/dum
                 nric(i,k) = nprc(k)*cldm(i,k)*dz(i,k)/cldmax(i,k)/dum
               else
                 if (qric(i,k-1) .ge. qsmall) then
                   dum = umr(k-1)
                   dum1 = unr(k-1)
                 end if

! hm add 4/17/06, no autoconversion of rain number if rain/snow falling 
! from above
! this assumes that new drizzle drops formed by autoconversion are 
! rapidly collected by the existing rain/snow particles from above
 
!RSH 2011:
! NCAR allows no autoconversion of rain number if rain/snow falling from 
! above. this assumes that new drizzle drops formed by autoconversion are 
! rapidly collected by the existing rain/snow particles falling from above.
! Marc's code allowed autoconversion to change rain number, so  variable 
! allow_all_cldtop_collection was introduced, which when .true. would turn
! off this effect. By default, it is .false. for GFDL (as in MG) in both 
! this subroutine and in the NCAR subroutine(cldwat2m_micro.F90), in 
! contrast to the original NCAR code.

                 if (allow_all_cldtop_collection) then
                   if (qric(i,k-1) .ge. 1.e-9_mg_pr .or.   &
                                 qniic(i,k-1) .ge. 1.e-9_mg_pr) then
                     nprc(k) = 0._mg_pr
                   end if
                 endif  !  allow_all_cldtop_collection

                 qric(i,k) = (rho(i,k-1)*umr(k-1)*qric(i,k-1)*  &
                                                        cldmax(i,k-1) + &
                             (rho(i,k)*dz(i,k)*  &
                              ((pra(k-1) + prc(k))*cldm(i,k) +    &
                                 (pre(k-1) - pracs(k-1) - mnuccr(k-1))*  &
                                                  cldmax(i,k))))/  &
                                                (dum*rho(i,k)*cldmax(i,k))
                 nric(i,k) = (rho(i,k-1)*unr(k-1)*nric(i,k-1)*   &
                                                         cldmax(i,k-1) + &
                             (rho(i,k)*dz(i,k)*(nprc(k)*cldm(i,k) +   &
                               (nsubr(k-1) - npracs(k-1) - nnuccr(k-1) + &
                                        nragg(k-1))*cldmax(i,k))))/   &
                                             (dum1*rho(i,k)*cldmax(i,k))
               end if

!  cloudice to snow autoconversion. if parameter one_ice is .true., then
!  cloudice and snow are retained in the same variable (as in the R-K 
!  microphysics) 
               IF (.NOT. one_ice) THEN

                 if (t(i,k) .le. tfreeze .and. qiic(i,k) .ge. qsmall) then

                   IF (auto_conv_ice_choice .EQ. 1) THEN  
!.......................................................................
! Autoconversion of cloud ice to snow
! similar to Ferrier (1994)
! note: assumes autoconversion timescale of 180 sec
                     nprci(k) = n0i(k)/(lami(k)*auto_conv_time_scale)*  &
                                                         exp(-lami(k)*dcs)
                     prci(k) = pi*rhoi*n0i(k)/(6._mg_pr*  &
                                                  auto_conv_time_scale)* &
                               (dcs**3/lami(k) + 3._mg_pr*dcs**2/ &
                                  lami(k)**2 + 6._mg_pr*dcs/lami(k)**3 +  &
                                     6._mg_pr/lami(k)**4)*exp(-lami(k)*dcs)

                   ELSE IF (auto_conv_ice_choice .EQ. 2) THEN
!.......................................................................
! AUTOCONVERSION OF CLOUD ICE TO SNOW
! FOLLOWING HARRINGTON ET AL. (1995) WITH MODIFICATION
! HERE IT IS ASSUMED THAT AUTOCONVERSION CAN ONLY OCCUR WHEN THE
! ICE IS GROWING, I.E. IN CONDITIONS OF ICE SUPERSATURATION
                     prci(k) = 0._mg_pr
                     nprci(k) = 0._mg_pr
                     esi(i,k) = polysvp_i(t(i,k))
                     qvi = 0.622_mg_pr*esi(i,k)/    &
                                            (pfull(i,k) - d378*esi(i,k))
                     dqsidt = xxls*qvi/(rv*t(i,k)**2)
                     abi(i,k) = 1._mg_pr + dqsidt*xxls/cpp
                     IF (Q(i,k) - QVI .GT. qsmall) THEN
                       NPRCI(K) = 4._mg_pr/(DCS*RHOI)*(Q(i,k) - QVI)*  &
                                    rho(i,k)*n0i(K)*   &
                                         EXP(-lami(K)*dcs)*dv(i,k)/abi(i,k)
                       NPRCI(K) = MIN(NPRCI(K), niic(i,K)/deltat)
                       PRCI(K) = MIN(PI*RHOI*DCS**3/6._mg_pr*NPRCI(K),   &
                                                         qiic(i,k)/deltat)
                     END IF

                   ELSE IF (auto_conv_ice_choice  .EQ. 3) THEN  
!.......................................................................
! Autoconversion of cloud ice to snow
! similar to Ferrier (1994)
! note: assumes autoconversion timescale of 180 sec
                     IF (lami(k) .LT. 1._mg_pr/autoconv_ice_thr) THEN 
                       nprci(k) = n0i(k)/(lami(k)*auto_conv_time_scale)* &
                                                         exp(-lami(k)*dcs)
                       prci(k) = pi*rhoi*n0i(k)/(6._mg_pr*    &
                                                 auto_conv_time_scale)* &
                                     (dcs**3/lami(k) +     &
                                          3._mg_pr*dcs**2/lami(k)**2 + &
                                              6._mg_pr*dcs/lami(k)**3 +   &
                                    6._mg_pr/lami(k)**4)*exp(-lami(k)*dcs) 
                     END IF

                   ELSE IF (auto_conv_ice_choice .EQ. 4) THEN  
! ferrier, 1994, 4.54    
                     dum = 1._mg_pr/autoconv_ice_thr 
                     IF (lami(k) .LT. dum) THEN 
                       nprci(k) = niic(i,K)/deltat*(1._mg_pr -   &
                                                      (lami(k)/dum )**3 )
                       prci(k) = qiic(i,k)/deltat*(1._mg_pr -   &
                                     (lami(k)/dum)**3)/   &
                                           (1._mg_pr + (lami(k)/dum)**3)
                     END IF

                   ELSE IF (auto_conv_ice_choice .EQ. 5) THEN  
!.......................................................................
! Autoconversion of cloud ice to snow
! similar to Ferrier (1994)
! note: assumes autoconversion timescale of 180 sec
                     IF (qiic(i,k) .GE. auto_conv_m_thresh) THEN 
                       nprci(k) = n0i(k)/(lami(k)*auto_conv_time_scale)*&
                                                        exp(-lami(k)*dcs)
                       prci(k) = pi*rhoi*n0i(k)/(6._mg_pr*  &
                                                  auto_conv_time_scale)* &
                                  (dcs**3/lami(k) +    &
                                      3._mg_pr*dcs**2/lami(k)**2 + &
                                          6._mg_pr*dcs/lami(k)**3 +  &
                                     6._mg_pr/lami(k)**4)*exp(-lami(k)*dcs)
                     END IF
                   END IF
                 else ! t(i,k) .le. tfreeze .and. qiic(i,k) .ge. qsmall
                   prci(k) = 0._mg_pr
                   nprci(k) = 0._mg_pr
                 end if ! t(i,k) .le. tfreeze .and. qiic(i,k) .ge. qsmall
               ELSE ! one_ice
                 prci(k) = 0._mg_pr
                 nprci(k) = 0._mg_pr
               END IF ! one_ice

! add autoconversion on current level to flux from level above to get 
! provisional snow mixing ratio and number concentration (qniic and nsic)

               dum = (asn(i,k)*dcs**bs)
               dum1 = (asn(i,k)*dcs**bs)

               if (k .eq. 1) then
                 qniic(i,k) = prci(k)*cldm(i,k)*dz(i,k)/cldmax(i,k)/dum
                 nsic(i,k) = nprci(k)*cldm(i,k)*dz(i,k)/cldmax(i,k)/dum
               else
                 if (qniic(i,k-1) .ge. qsmall) then
                   dum = ums(k-1)
                   dum1 = uns(k-1)
                 end if

!++ag fixed snow bug (from Morrison nov.27.2007)
!                qniic(i,k) = (rho(i,k-1)*ums(k-1)*qniic(i,k-1)*  &
!                                                     cldmax(i,k-1) + &
!                              (rho(i,k)*dz(i,k)*(prci(k)*cldm(i,k) +  &
!                                  (prai(k-1) + psacws(k-1) + prci(k-1) + &
!                                        bergs(k-1))*cldmax(i,k))))/  &
!                                          (dum*rho(i,k)*cldmax(i,k))

                 qniic(i,k) = (rho(i,k-1)*ums(k-1)*qniic(i,k-1)*  &
                                                       cldmax(i,k-1) + &
                              (rho(i,k)*dz(i,k)*((prci(k) + prai(k-1) + &
                                  psacws(k-1) + bergs(k-1))*cldm(i,k) +  &
                                 (prds(k-1) + pracs(k-1) + mnuccr(k-1))*  &
                                                         cldmax(i,k))))/  &
                                               (dum*rho(i,k)*cldmax(i,k))
!--ag

                 nsic(i,k) = (rho(i,k-1)*uns(k-1)*nsic(i,k-1)*  &
                                                       cldmax(i,k-1) + &
                             (rho(i,k)*dz(i,k)*(nprci(k)*cldm(i,k) +  &
                                (nsubs(k-1) + nsagg(k-1) + nnuccr(k-1))*  &
                                                        cldmax(i,k))))/   &
                                              (dum1*rho(i,k)*cldmax(i,k))
               end if

! if precip mix ratio is zero so should number concentration
               if (qniic(i,k) .lt. qsmall) then
                 qniic(i,k) = 0._mg_pr
                 nsic(i,k) = 0._mg_pr
               end if
               if (qric(i,k) .lt. qsmall) then
                 qric(i,k) = 0._mg_pr
                 nric(i,k) = 0._mg_pr
               end if

! make sure number concentration is a positive number to avoid 
! taking root of negative later
               nric(i,k) = max(nric(i,k), 0._mg_pr)
               nsic(i,k) = max(nsic(i,k), 0._mg_pr)

!.......................................................................
! get size distribution parameters for precip
!......................................................................
! rain
               if (qric(i,k) .ge. qsmall) then
                 lamr(k) = (pi*rhow*nric(i,k)/qric(i,k))**  &
                                                        (1._mg_pr/3._mg_pr)
                 n0r(k) = nric(i,k)*lamr(k)

! check for slope
                 lammax = 1._mg_pr/20.e-6_mg_pr
                 lammin = 1._mg_pr/500.e-6_mg_pr

! adjust vars
                 if (lamr(k) .lt. lammin) then
                   lamr(k) = lammin
                   n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
                   nric(i,k) = n0r(k)/lamr(k)
                 else if (lamr(k).gt.lammax) then
                   lamr(k) = lammax
                   n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
                   nric(i,k) = n0r(k)/lamr(k)
                 end if

! provisional rain number and mass weighted mean fallspeed (m/s)
                 unr(k) = min (           &
                            arn(i,k)*gamma_mg(1._mg_pr + br)/lamr(k)**br, &
                                                      9.1_mg_pr*rhof(i,k))
                 umr(k) = min (       &
                            arn(i,k)*gamma_mg(4._mg_pr + br)/   &
                                                 (6._mg_pr*lamr(k)**br), &
                                                       9.1_mg_pr*rhof(i,k))
               else
                 lamr(k) = 0._mg_pr
                 n0r(k) = 0._mg_pr
                 umr(k) = 0._mg_pr
                 unr(k) = 0._mg_pr
               end if

!......................................................................
! snow

               if (qniic(i,k) .ge. qsmall) then
                 lams(k) = (gamma_mg(1._mg_pr + ds)*cs*nsic(i,k)/ &
                                              qniic(i,k))**(1._mg_pr/ds)
                 n0s(k) = nsic(i,k)*lams(k)

! check for slope
                 lammax = 1._mg_pr/min_diam_ice
                 lammin = 1._mg_pr/2000.e-6_mg_pr

! adjust vars
                 if (lams(k) .lt. lammin) then
                   lams(k) = lammin
                   n0s(k) = lams(k)**(ds + 1._mg_pr)*qniic(i,k)/   &
                                              (cs*gamma_mg(1._mg_pr + ds))
                   nsic(i,k) = n0s(k)/lams(k)
                 else if (lams(k) .gt. lammax) then
                   lams(k) = lammax
                   n0s(k) = lams(k)**(ds + 1._mg_pr)*qniic(i,k)/    &
                                               (cs*gamma_mg(1._mg_pr + ds))
                   nsic(i,k) = n0s(k)/lams(k)
                 end if

! provisional snow number and mass weighted mean fallspeed (m/s)
                 ums(k) = min (      &
                           asn(i,k)*gamma_mg(4._mg_pr + bs)/   &
                                            (6._mg_pr*lams(k)**bs),  &
                                                    max_vt_snow*rhof(i,k))
                 uns(k) = min (       &
                           asn(i,k)*gamma_mg(1._mg_pr + bs)/lams(k)**bs,  &
                                                    max_vt_snow*rhof(i,k))
               else
                 lams(k) = 0._mg_pr
                 n0s(k) = 0._mg_pr
                 ums(k) = 0._mg_pr
                 uns(k) = 0._mg_pr
               end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! heterogeneous freezing of cloud water

               if (qcic(i,k) .ge. qsmall ) then 
                 mnuccc(k) = 0._mg_pr
                 nnuccc(k) = 0._mg_pr

!contact freezing
                 IF (do_contact_frz) THEN 
                   tc = ttmp - tfreeze
                   IF (tc .LE. -3._mg_pr .AND. tc .GE. -40._mg_pr) THEN
                     IF (rbar_dust(i,k) .GT. 0._mg_pr) THEN

! mean free path
                       dum = 7.37_mg_pr*t(i,k)/(288._mg_pr*10._mg_pr*  &
                                                    pfull(i,k))/100._mg_pr

! effective diffusivity based on Brownian collection
                       dap = 4._mg_pr*pi*1.38e-23_mg_pr*t(i,k)*  &
                             (1._mg_pr + dum/rbar_dust(i,k))/ &
                                     (6._mg_pr*pi*rbar_dust(i,k)*mu(i,k))

! number of contact nucleii similar Young as in  Liu et al.
                       IF (n_contact_opt .eq. 1 ) THEN
                         NACNT = ndust(i,k)/rho(i,k)*  &
                                       (tfreeze - 3._mg_pr - T(i,k))**1.3
                       ELSE IF (n_contact_opt .eq. 2) THEN

!similar Meyers et al., 1992, Eq. 2.6
                         NACNT = ndust(i,k)/rho(i,k)*    &
                                        EXP(-2.8_mg_pr + 0.262_mg_pr*   &
                                                       (tfreeze - T(i,k)))
                       END IF

                       MNUCCC(K) = sfac3*PI*PI/3._mg_pr*RHOW*DAP*NACNT*  &
                                   EXP(LOG(CDIST1(K)) +   &
                                     LOG(GAMMA_mg(PGAM(K) + 5._mg_pr)) - &
                                                    4._mg_pr*LOG(LAMC(K)))
                       NNUCCC(K) = 2._mg_pr*PI*DAP*NACNT*CDIST1(K)*      &
                                     GAMMA_mg(PGAM(K) + 2._mg_pr)/LAMC(K) 

                     END IF
                   END IF
                 END IF  ! (do_contact_frz)

                 IF (do_bigg_frz) THEN 

! immersion freezing (Bigg, 1953)
                   if (t(i,k) .lt. tfreeze - 4._mg_pr) then
!                    mnuccc(k) = gamma_mg(qcvar+2._mg_pr)/    &
!                                        (gamma_mg(qcvar)*qcvar**2)* &
                     mnuccc(k) = mnuccc(k) + sfac2*pi*pi/36._mg_pr*rhow* &
                                          cdist1(k)*  &
                                             gamma_mg(7._mg_pr + pgam(k))*&
!RSH BUGFIX email 8/9/10
                                             bimm*(exp(aimm*(   &
                                          tfreeze - t(i,k))) - 1._mg_pr)/ &
                                                      lamc(k)**3/lamc(k)**3
!                    nnuccc(k) = gamma_mg(qcvar+1._mg_pr)/   &
!                                             (gamma_mg(qcvar)*qcvar)* &
                     nnuccc(k) =  nnuccc(k) + sfac3*pi/6._mg_pr*cdist1(k)*&
                                       gamma_mg(pgam(k) + 4._mg_pr)* &
!RSH BUGFIX email 8/9/10
                                                  bimm*(exp(aimm*(  &
                                 tfreeze - t(i,k))) - 1._mg_pr)/lamc(k)**3

                   end if
                 END IF

                 IF (limit_droplet_freeze_opt .EQ. 1) THEN

! hm add 11/17/06
! make sure number of droplets frozen does not exceed available ice 
! nuclei concentration
! this prevents 'runaway' droplet freezing
                   if (nnuccc(k) .gt. nnuccd(k)/cldm(i,k)) then
                     dum = (nnuccd(k)/cldm(i,k))/nnuccc(k)

! scale mixing ratio of droplet freezing with limit
                     mnuccc(k) = mnuccc(k)*dum
                     nnuccc(k) = nnuccd(k)/cldm(i,k)
                   end if
                 ELSE IF (limit_droplet_freeze_opt .EQ. 2) THEN
                   dum1 = nnuccc(k)*deltat
                   dum2 = ndust(i,k)/rho(i,k)/cldm(i,k)
                   if (dum1 .gt. dum2) then
                     dum = dum2/dum1
 
! scale mixing ratio of droplet freezing with limit
                     mnuccc(k) = mnuccc(k)*dum
                     nnuccc(k) = nnuccc(k)*dum
!!!                  nnuccc(k) = nnuccd(k)/cldm(i,k)
                   end if
                 END IF
               else ! (qcic(i,k) .ge. qsmall )  
                 mnuccc(k) = 0._mg_pr
                 nnuccc(k) = 0._mg_pr
               end if ! (qcic(i,k) .ge. qsmall )  

!.......................................................................
! snow self-aggregation from passarelli, 1978, used by reisner, 1998
! this is hard-wired for bs = 0.4 for now
! ignore self-collection of cloud ice

               if (qniic(i,k) .ge. qsmall .and. t(i,k) .le. tfreeze) then
                 nsagg(k) = -1108._mg_pr*asn(i,k)*Eii* &
                                pi**((1._mg_pr - bs)/3._mg_pr)*   &
                                 rhosn**((-2._mg_pr - bs)/3._mg_pr)*  &
                                  rho(i,k)**((2._mg_pr + bs)/3._mg_pr)*  &
                                  qniic(i,k)**((2._mg_pr + bs)/3._mg_pr)* &
                                    (nsic(i,k)*rho(i,k))**  &
                                           ((4._mg_pr-bs)/3._mg_pr)/ &
                                             (4._mg_pr*720._mg_pr*rho(i,k))
               else
                 nsagg(k)=0._mg_pr
               end if

               IF (.NOT. one_ice) THEN

!.......................................................................
! accretion of cloud droplets onto snow/graupel
! here use continuous collection equation with
! simple gravitational collection kernel
! ignore collisions between droplets/cloud ice
! since minimum size ice particle for accretion is 50 - 150 micron

! ignore collision of snow with droplets above freezing

                 if (qniic(i,k) .ge. qsmall .and.     &
                             t(i,k) .le. tfreeze .and. &
                                  qcic(i,k) .ge. qsmall) then

! put in size dependent collection efficiency
! mean diameter of snow is area-weighted, since
! accretion is function of crystal geometric area
! collection efficiency is from stoke's law (Thompson et al. 2004)

                   dc0 = (pgam(k) + 1._mg_pr)/lamc(k)
                   ds0 = 1._mg_pr/lams(k)
                   dum = dc0*dc0*uns(k)*rhow/(9._mg_pr*mu(i,k)*ds0)
                   eci = dum*dum/((dum + 0.4_mg_pr)*(dum + 0.4_mg_pr))
                   eci = max(eci, 0._mg_pr)
                   eci = min(eci, 1._mg_pr)

! no impact of sub-grid distribution of qc since psacws
! is linear in qc

                   psacws(k) = pi/4._mg_pr*asn(i,k)*qcic(i,k)*rho(i,k)* &
                                   n0s(k)*Eci*gamma_mg(bs + 3._mg_pr)/ &
                                                  lams(k)**(bs + 3._mg_pr)
                   npsacws(k) = pi/4._mg_pr*asn(i,k)*ncic(i,k)*rho(i,k)* &
                                    n0s(k)*Eci*gamma_mg(bs + 3._mg_pr)/ &
                                                  lams(k)**(bs + 3._mg_pr)
                 else
                   psacws(k) = 0._mg_pr
                   npsacws(k) = 0._mg_pr
                 end if

                 psacws_o(k) = 0._mg_pr
                 npsacws_o(k) = 0._mg_pr
               ELSE !one_ice 
                 psacws_o(k) = 0._mg_pr
                 npsacws_o(k) =0._mg_pr
                 IF (lami(k) .gt. 1._mg_pr/50.e-6_mg_pr) THEN 

! provisional snow number and mass weighted mean fallspeed (m/s)
                   ums(k) = min (       &
                              asn(i,k)*gamma_mg(4._mg_pr + bs)/    &
                                           (6._mg_pr*lami(k)**bs),   &
                                                              max_vt_snow)
                   uns(k) = min (      &
                              asn(i,k)*gamma_mg(1._mg_pr + bs)/     &
                                                 lami(k)**bs, max_vt_snow)

!PROBLEM k-1 indices 11/20/11
                   qiic(i,k) = (rho(i,k-1)*ums(k-1)*qiic(i,k-1)*     &
                                                         cldmax(i,k-1) + &
                               (rho(i,k)*dz(i,k)*((psacws_o(k-1))*   &
                                                           cldm(i,k))))/ &
                                                (dum*rho(i,k)*cldmax(i,k))
                   if (qiic(i,k) .ge. qsmall .and.    &
                                 t(i,k) .le. tfreeze .and. &
                                        qcic(i,k) .ge. qsmall) then

! put in size dependent collection efficiency
! mean diameter of snow is area-weighted, since
! accretion is function of crystal geometric area
! collection efficiency is from stoke's law (Thompson et al. 2004)
                     dc0 = (pgam(k) + 1._mg_pr)/lamc(k)
                     ds0 = 1._mg_pr/lami(k)
                     dum = dc0*dc0*uns(k)*rhow/(9._mg_pr*mu(i,k)*ds0)
                     eci = dum*dum/((dum + 0.4_mg_pr)*(dum + 0.4_mg_pr))
                     eci = max(eci, 0._mg_pr)
                     eci = min(eci, 1._mg_pr)

! no impact of sub-grid distribution of qc since psacws
! is linear in qc

                     psacws_o(k) = pi/4._mg_pr*asn(i,k)*qcic(i,k)*  &
                                   rho(i,k)*n0i(k)*Eci*   &
                                     gamma_mg(bs + 3._mg_pr)/ &
                                                 lami(k)**(bs + 3._mg_pr)
                     npsacws_o(k) = pi/4._mg_pr*asn(i,k)*ncic(i,k)*  &
                                   rho(i,k)*n0i(k)*Eci*   &
                                      gamma_mg(bs + 3._mg_pr)/ &
                                                 lami(k)**(bs + 3._mg_pr)
                   else
                     psacws_o(k) = 0._mg_pr
                     npsacws_o(k) = 0._mg_pr
                   end if
                   psacws(k) = 0._mg_pr
                   npsacws(k) = 0._mg_pr
                 END IF !size
               END IF !one_ice 

               IF (.NOT. one_ice) THEN 
!.......................................................................
! accretion of rain water by snow
! formula from ikawa and saito, 1991, used by reisner et al., 1998

                 if (qric(i,k) .ge. 1.e-8_mg_pr .and.     &
                         qniic(i,k) .ge. 1.e-8_mg_pr .and. & 
                                       t(i,k).le. tfreeze) then
                   pracs(k) = pi*pi*ecr*(((1.2_mg_pr*umr(k) -    &
                                                 0.95_mg_pr*ums(k))**2 + &
                               0.08_mg_pr*ums(k)*umr(k))**0.5_mg_pr*rhow* &
                                   rho(i,k)*n0r(k)*n0s(k)* &
                                     (5._mg_pr/(lamr(k)**6*lams(k)) + &
                                      2._mg_pr/(lamr(k)**5*lams(k)**2) + &
                                     0.5_mg_pr/(lamr(k)**4*lams(k)**3)))

                   npracs(k) = pi/2._mg_pr*rho(i,k)*ecr*(1.7_mg_pr*  &
                                                   (unr(k) - uns(k))**2 + &
                                 0.3_mg_pr*unr(k)*uns(k))**0.5_mg_pr*  &
                                                        n0r(k)*n0s(k)* &
                                  (1._mg_pr/(lamr(k)**3*lams(k)) + &
                                   1._mg_pr/(lamr(k)**2*lams(k)**2) + &
                                   1._mg_pr/(lamr(k)*lams(k)**3))
                 else
                   pracs(k)=0._mg_pr
                   npracs(k)=0._mg_pr
                 end if
               ELSE !(.NOT. one_ice) 
                 pracs(k)=0._mg_pr
                 npracs(k)=0._mg_pr
               END IF ! (.NOT. one_ice) 

               IF (.NOT. one_ice) THEN

!.......................................................................
! heterogeneous freezing of rain drops
! follows from Bigg (1953)
                 if (t(i,k) .lt. tfreeze -4._mg_pr .and.    &
                        qric(i,k) .ge. qsmall) then
                   mnuccr(k) = 20._mg_pr*pi*pi*rhow*nric(i,k)*bimm* &
!RSH BUGFIX email 8/9/10
                               (exp(aimm*(tfreeze - t(i,k))) - 1._mg_pr)/ &
                                                  lamr(k)**3/lamr(k)**3

                   nnuccr(k) = pi*nric(i,k)*bimm* &
!RSH BUGFIX email 8/9/10
                               (exp(aimm*(tfreeze - t(i,k))) - 1._mg_pr)/ &
                                                  lamr(k)**3
                 else
                   mnuccr(k) = 0._mg_pr
                   nnuccr(k) = 0._mg_pr
                 end if
               ELSE ! (.NOT. one_ice) 
                 mnuccr(k) = 0._mg_pr
                 nnuccr(k) = 0._mg_pr
               END IF ! (.NOT. one_ice) 

!.......................................................................
! accretion of cloud liquid water by rain
! formula from Khrouditnov and Kogan (2000)
! gravitational collection kernel, droplet fall speed neglected

               if (qric(i,k) .ge. qsmall .and.    &
                               qcic(i,k) .ge. qsmall) then

! include sub-grid distribution of cloud water
!                 pra(k) = gamma_mg(qcvar+1.15_mg_pr)/   &
!                          (gamma_mg(qcvar)*qcvar**1.15_mg_pr) * &
                  pra(k) = sfac4*67._mg_pr*  &
                                         (qcic(i,k)*qric(i,k))**1.15_mg_pr
                  npra(k) = pra(k)/(qcic(i,k)/ncic(i,k))
               else
                 pra(k) = 0._mg_pr
                 npra(k) = 0._mg_pr
               end if

!.......................................................................
! Self-collection of rain drops
! from Beheng(1994)

               if (qric(i,k) .ge. qsmall) then
                 nragg(k) = -8._mg_pr*nric(i,k)*qric(i,k)*rho(i,k)
               else
                 nragg(k) = 0._mg_pr
               end if

!.......................................................................
! Accretion of cloud ice by snow
! For this calculation, it is assumed that the Vs >> Vi
! and Ds >> Di for continuous collection

               IF (.NOT. one_ice) THEN
                 if (qniic(i,k) .ge. qsmall .and.   &
                                qiic(i,k) .ge. qsmall  .and.  &
                                         t(i,k) .le. tfreeze) then
                   prai(k) = pi/4._mg_pr*asn(i,k)*qiic(i,k)*rho(i,k)* &
                                n0s(k)*Eii*gamma_mg(bs + 3._mg_pr)/ &
                                                  lams(k)**(bs + 3._mg_pr)
                   nprai(k) = pi/4._mg_pr*asn(i,k)*niic(i,k)*rho(i,k)*  &
                                   n0s(k)*Eii*gamma_mg(bs + 3._mg_pr)/ &
                                                  lams(k)**(bs + 3._mg_pr)
                 else
                   prai(k) = 0._mg_pr
                   nprai(k) = 0._mg_pr
                 end if
               ELSE
                 prai(k) = 0._mg_pr
                 nprai(k) = 0._mg_pr
               END IF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate evaporation/sublimation of rain and snow
! note: evaporation/sublimation occurs only in cloud-free portion of 
! grid cell
! in-cloud condensation/deposition of rain and snow is neglected
! except for transfer of cloud water to snow through bergeron process

! initialize evap/sub tendncies
               pre(k) = 0._mg_pr
               prds(k) = 0._mg_pr

! evaporation of rain
! only calculate if there is some precip fraction > cloud fraction

!RSH 8/7/12: return to using qsmall to avoid model blowups
!              if (qcic(i,k) + qiic(i,k) .lt. 1.e-6_mg_pr .or.   &
!                                       cldmax(i,k) .gt. cldm(i,k)) then
               if (qcic(i,k) + qiic(i,k) .lt. qsmall .or.   &
                                        cldmax(i,k) .gt. cldm(i,k)) then

! set temporary cloud fraction to zero if cloud water + ice is very small
! this will ensure that evaporation/sublimation of precip occurs over
! entire grid cell, since min cloud fraction is specified otherwise

!RSH 8/7/12: return to using qsmall to avoid model blowups
!                if (qcic(i,k) + qiic(i,k) .lt. 1.e-6_mg_pr) then
                 if (qcic(i,k) + qiic(i,k) .lt. qsmall) then
                   dum = 0._mg_pr
                 else
                   dum = cldm(i,k)
                 end if
                 ttmp = t(i,k)

! recalculate saturation vapor pressure for liquid and ice
                 esl(i,k) = polysvp_l(t(i,k))
                 esi(i,k) = polysvp_i(t(i,k))
                 esn = esl(i,k)
                 qsn = min (   &
                         epsqs*esn/(pfull(i,k) - (1._mg_pr - epsqs)*esn), &
                                                                  1._mg_pr)
                 qsn = max (qsn, 0._mg_pr)

! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING
                 IF (esi(i,k) .GT. esl(i,k)) esi(i,k) = esl(i,k)

! calculate q for out-of-cloud region
                 qclr = (q(i,k) - dum*qsn)/(1._mg_pr - dum)
                 if (qric(i,k) .ge. qsmall) then
                   qvs = 0.622_mg_pr*esl(i,k)/(pfull(i,k) - d378*esl(i,k))
                   dqsdt = xxlv*qvs/(rv*t(i,k)**2)
                   ab = 1._mg_pr + dqsdt*xxlv/cpp
                   epsr = 2._mg_pr*pi*n0r(k)*rho(i,k)*Dv(i,k)* &
                         (f1r/(lamr(k)*lamr(k)) +    &
                             f2r*(arn(i,k)*rho(i,k)/mu(i,k))**0.5_mg_pr* &
                                  sc(i,k)**(1._mg_pr/3._mg_pr)*   &
                            gamma_mg(5._mg_pr/2._mg_pr + br/2._mg_pr)/ &
                              (lamr(k)**(5._mg_pr/2._mg_pr + br/2._mg_pr)))
                   pre(k) = epsr*(qclr - qvs)/ab

! only evaporate in out-of-cloud region
! and distribute across cldmax
                   pre(k) = min(pre(k)*(cldmax(i,k) - dum), 0._mg_pr)
                   pre(k) = pre(k)/cldmax(i,k)
                 end if

                 IF (.NOT. one_ice) THEN

! sublimation of snow
                   if (qniic(i,k) .ge. qsmall) then
                     qvi = 0.622_mg_pr*esi(i,k)/    &
                                               (pfull(i,k) - d378*esi(i,k))
                     dqsidt =  xxls*qvi/(rv*t(i,k)**2)
                     abi(i,k) = 1._mg_pr + dqsidt*xxls/cpp
                     dumt1 = 2._mg_pr*pi*n0s(k)*rho(i,k)*Dv(i,k)
                     dumt2 = sc(i,k)**(1._mg_pr/3._mg_pr)*   &
                              gamma_mg(5._mg_pr/2._mg_pr + bs/2._mg_pr)/ &
                               (lams(k)**(5._mg_pr/2._mg_pr + bs/2._mg_pr))
                     epss = 2._mg_pr*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
                                (f1s/(lams(k)*lams(k)) +     &
                              f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_mg_pr* &
                               sc(i,k)**(1._mg_pr/3._mg_pr)*    &
                               gamma_mg(5._mg_pr/2._mg_pr + bs/2._mg_pr)/ &
                              (lams(k)**(5._mg_pr/2._mg_pr + bs/2._mg_pr)))
                     prds(k) = epss*(qclr - qvi)/abi(i,k)
 
! only sublimate in out-of-cloud region and distribute over cldmax
                     prds(k) = min(prds(k)*(cldmax(i,k) - dum), 0._mg_pr)
                     prds(k) = prds(k)/cldmax(i,k)
                   end if
                 ELSE   ! (.NOT. one_ice) 
                   prds(k) = 0._mg_pr
                 END IF ! (.NOT. one_ice) 

! hm add 2/2/07, make sure RH not pushed above 100%
! get updated RH at end of time step based on 
! cloud water/ice condensation/evap

                 qtmp = q(i,k) - (D_eros_l(i,k) + D_eros_i(i,k) +  &
                        cmel(i,k) + cmei(i,k) + qvdep_qi(i,k) +    &
                                    (pre(k) + prds(k))*cldmax(i,k))*deltat

!bug?  2/12/09
!!$     ttmp=t(i,k)+((D_eros_l(i,k)+cmel(i,k)+pre(k)*cldmax(i,k))*xxlv+ &
!!$                (D_eros_i(i,k)+cmei(i,k)+prds(k))*cldmax(i,k)*xxls)*deltat/cpp

                 ttmp = t(i,k) + ((D_eros_l(i,k) + cmel(i,k) +   &
                                               pre(k)*cldmax(i,k))*xxlv + &
                        (D_eros_i(i,k) + cmei(i,k) + qvdep_qi(i,k) + &
                                     prds(k)*cldmax(i,k))*xxls)*deltat/cpp
                 ttmp = MAX(lowest_temp_for_sublimation,   &
                                                    min(ttmp, 323._mg_pr))
                 eslt = polysvp_l(ttmp)
                 esit = polysvp_i(ttmp)
                 esn = eslt
    
                 qsn = min (epsqs*esn/  &
                           (pfull(i,k) - (1._mg_pr - epsqs)*esn), 1._mg_pr)
                 qsn = max(qsn, 0._mg_pr)
      
! modify precip evaporation rate if q > qsat
                 if (qtmp .gt. qsn ) then
                   if (pre(k) + prds(k) .lt. -1.e-20_mg_pr) then
                     dum1 = pre(k)/(pre(k) + prds(k))

! recalculate q and t after cloud water cond but without precip evap
                     qtmp = q(i,k) - (D_eros_l(i,k) + D_eros_i(i,k) +  &
                             cmel(i,k) + cmei(i,k) + qvdep_qi(i,k))*deltat

!bug 2/12/09
!!$     ttmp=t(i,k)+(D_eros_l(i,k)+cmel(i,k)*xxlv+ &
!!$                D_eros_i(i,k)+cmei(i,k)*xxls)*deltat/cpp

                     ttmp = t(i,k) + ((D_eros_l(i,k) + cmel(i,k))*xxlv +&
                             (D_eros_i(i,k) + cmei(i,k) + qvdep_qi(i,k))*&
                                                           xxls)*deltat/cpp
                     eslt = polysvp_l(ttmp)
                     esit = polysvp_i(ttmp)
                     esn = eslt
            
                     qsn = min( epsqs*esn/   &
                           (pfull(i,k) - (1._mg_pr - epsqs)*esn), 1._mg_pr)
                     qsn=max(qsn, 0._mg_pr)
                     dum = (qtmp - qsn)/   &
                                  (1._mg_pr + xxlv**2*qsn/(cpp*rv*ttmp**2))
                     dum = min(dum, 0._mg_pr)

! modify rates if needed, divide by cldmax to get local (in-precip) value
                     pre(k) = dum*dum1/deltat/cldmax(i,k)
                     qsn = min( epsqs*esit/    &
                                (pfull(i,k) - (1._mg_pr - epsqs)*esit),  &
                                                                  1._mg_pr)
                     dum = (qtmp - qsn)/(1._mg_pr + xxls**2*qsn/   &
                                                         (cpp*rv*ttmp**2)) 
                     dum = min(dum, 0._mg_pr)
                     prds(k) = dum*(1._mg_pr - dum1)/deltat/cldmax(i,k)
                   end if
                 end if
               end if

               IF (.NOT. one_ice) THEN

! bergeron process - evaporation of droplets and deposition onto snow
! bergeron process for snow is neglected for now.............
                 if (do_berg_snow ) then
                   if (qniic(i,k) .ge. qsmall .and.     &
                                qcic(i,k) .ge. qsmall .and.     &
                                        t(i,k) .lt. tfreeze) then
                     qvs = 0.622_mg_pr*esl(i,k)/    &
                                          (pfull(i,k) - d378*esl(i,k))
                     qvi = 0.622_mg_pr*esi(i,k)/   &
                                          (pfull(i,k) - d378*esi(i,k))

!8/1/12: place limits to avoid negative values which may occur at low pfull
!        prevents model blowups
                     qvs = MAX(0._mg_pr, MIN(qvs, 1.0_mg_pr))
                     qvi = MAX(0._mg_pr, MIN(qvi, 1.0_mg_pr))
                     dqsidt = xxls*qvi/(rv*t(i,k)**2)
                     abi(i,k) = 1._mg_pr + dqsidt*xxls/cpp
                     epss = 2._mg_pr*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
                            (f1s/(lams(k)*lams(k)) + &
                              f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_mg_pr* &
                                sc(i,k)**(1._mg_pr/3._mg_pr)*    &
                               gamma_mg(5._mg_pr/2._mg_pr + bs/2._mg_pr)/ &
                             (lams(k)**(5._mg_pr/2._mg_pr + bs/2._mg_pr)))
!cms 2009/3/2        bergs(k) = epss*(qvs - qvi)/abi(i,k)
                     bergs(k) = epss*(qvs - qvi)/abi(i,k)
                   else
                     bergs(k) = 0._mg_pr
                   end if
                 else
                   bergs(k) = 0._mg_pr
                 endif
               ELSE
                 bergs(k) = 0._mg_pr
               END IF

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! conservation to ensure no negative values of cloud water/precipitation
! in case microphysical process rates are large

! make sure and use end-of-time step values for cloud water, ice, due
! condensation/deposition

! note: for check on conservation, processes are multiplied by omsm
! to prevent problems due to round off error

! since activation/nucleation processes are fast, need to take into account
! factor mtime = mixing timescale in cloud / model time step
! mixing time can be interpreted as cloud depth divided by sub-grid 
! vertical velocity
! for now mixing timescale is assumed to be 20 min
! could possibly be estimated better from model variables

               IF (total_activation) THEN 
! since act. is assumed to take place at cloud base:
!cms             mtime = deltat/1200._mg_pr
                 mtime = 1._mg_pr
               ELSE IF (dqa_activation) THEN 
! since Yi's formulation assumes activation at lateral cloud boundaries
                 mtime=1._mg_pr
               ENDIF

               qce = (qc(i,k) + (D_eros_l(i,k) + cmel(i,k) -    &
                                                        berg(i,k))*deltat)
               nce = (nc(i,k) + npccn(k)*deltat*mtime)/cldm(i,k)
               qie = (qi(i,k) + (D_eros_i(i,k) + cmei(i,k) +   &
                                        berg(i,k) + qvdep_qi(i,k))*deltat)
               nie = (ni(i,k) + nnuccd(k)*deltat*mtime)/cldm(i,k)

! conservation of qc
               dum = (prc(k) + pra(k) + mnuccc(k) + &
                       psacws(k) + bergs(k) + psacws_o(k))*cldm(i,k)*deltat
               if (dum .gt. qce) then
                 if (dum .gt. 1.e-30_mg_pr) then
                   ratio = qce/deltat/cldm(i,k)/(prc(k) + pra(k) +    &
                                   mnuccc(k) + psacws(k) + psacws_o(k) +  &
                                                            bergs(k))*omsm
                 else 
                   ratio = 0._mg_pr
                 endif
                 prc(k) = prc(k)*ratio
                 pra(k) = pra(k)*ratio
                 mnuccc(k) = mnuccc(k)*ratio
                 psacws(k) = psacws(k)*ratio
                 psacws_o(k) = psacws_o(k)*ratio
                 bergs(k) = bergs(k)*ratio
               end if

! conservation of nc
               dum = (nprc1(k) + npra(k) + nnuccc(k) + &
                      npsacws(k) + npsacws_o(k) - nsubc(k) -   &
                                                     nerosc(i,k))*deltat
               if (dum .gt. nce) then
                 if (dum .gt. 1.e-30_mg_pr) then
                   ratio = nce/deltat/(nprc1(k) + npra(k) + nnuccc(k) + &
                            npsacws(k) + npsacws_o(k) - nsubc(k) -   &
                                                        nerosc(i,k))*omsm
                 else
                   ratio = 0._mg_pr
                 end if
                 nprc1(k) = nprc1(k)*ratio
                 npra(k) = npra(k)*ratio
                 nnuccc(k) = nnuccc(k)*ratio
                 npsacws(k) = npsacws(k)*ratio
                 npsacws_o(k) = npsacws_o(k)*ratio
                 nsubc(k) = nsubc(k)*ratio
                 nerosc(i,k) =  nerosc(i,k)*ratio
               end if

! conservation of qi
               dum = (-mnuccc(k) + prci(k) + &
                        prai(k) - psacws_o(k))*cldm(i,k)*deltat
               if (dum .gt. qie) then
                 if (dum .gt. 1.e-30_mg_pr) then
                   ratio = (qie/deltat/cldm(i,k) + mnuccc(k) +    &
                                psacws_o(k))/(prci(k) + prai(k))*omsm
                 else
                   ratio = 0._mg_pr
                 end if
                 prci(k) = prci(k)*ratio
                 prai(k) = prai(k)*ratio
                 psacws_o(k) = psacws_o(k)*ratio
               end if

! conservation of ni
               dum = (nprci(k) + &
                        nprai(k) - nsubi(k) - nerosi(i,k))*deltat
               if (dum .gt. nie) then
                 if (dum .gt. 1.e-30_mg_pr) then
                   ratio = (nie/deltat)/(nprci(k) + nprai(k)    &
                                           - nsubi(k) - nerosi(i,k))*omsm
                 else
                   ratio = 0._mg_pr
                 end if
                 nprci(k) = nprci(k)*ratio
                 nprai(k) = nprai(k)*ratio
                 nsubi(k) = nsubi(k)*ratio
                 nerosi(i,k) = nerosi(i,k)*ratio
               end if

! for preciptiation conservation, use logic that vertical integral 
! of tendency from current level to top of model (i.e., qrtot) cannot 
! be negative

! conservation of rain mixing rat

               if (((prc(k) + pra(k))*cldm(i,k) +     &
                      (-mnuccr(k) + pre(k) - pracs(k))*cldmax(i,k))*  &
                             rho(i,k)*dz(i,k) + qrtot .lt. 0._mg_pr) then
                 if (-pre(k) + pracs(k) + mnuccr(k) .ge. qsmall) then
                   ratio = (qrtot/(rho(i,k)*dz(i,k)) + (prc(k) +   &
                                                   pra(k))*cldm(i,k))/&
                                    ((-pre(k) + pracs(k) + mnuccr(k))*   &
                                                        cldmax(i,k))*omsm 
                 else 
                   ratio = 0._mg_pr
                 end if
                 pre(k) = pre(k)*ratio
                 pracs(k) = pracs(k)*ratio
                 mnuccr(k) = mnuccr(k)*ratio
               end if

! conservation of nr
! for now neglect evaporation of nr

               IF (rain_evap_opt) THEN

! calculate evaporation of nr
                 if (pre(k) .lt. 0._mg_pr .and.      &
                                       qric(i,k) .ge. qsmall) then
                   nsubr(k) = pre(k)/qric(i,k)*nric(i,k)
                 else
                   nsubr(k) = 0._mg_pr
                 end if
               ELSE
                 nsubr(k)=0._mg_pr
               END IF

               if ((nprc(k)*cldm(i,k) + (-nnuccr(k) + nsubr(k) -    &
                      npracs(k) + nragg(k))*cldmax(i,k))*rho(i,k)*    &
                                     dz(i,k) + nrtot .lt. 0._mg_pr) then
                 if (-nsubr(k) - nragg(k) + npracs(k) +     &
                                              nnuccr(k) .ge. qsmall) then
                   ratio = (nrtot/(rho(i,k)*dz(i,k)) +     &
                                                nprc(k)*cldm(i,k))/   &
                             ((-nsubr(k) - nragg(k) + npracs(k) +   &
                                              nnuccr(k))*cldmax(i,k))*omsm
                 else 
                   ratio = 0._mg_pr
                 end if
                 nsubr(k) = nsubr(k)*ratio
                 npracs(k) = npracs(k)*ratio
                 nnuccr(k) = nnuccr(k)*ratio
                 nragg(k) = nragg(k)*ratio
               end if

! conservation of snow mix ratio
               if (((bergs(k) + psacws(k) + prai(k) +     &
                                                  prci(k))*cldm(i,k) +    &
                         (pracs(k) + mnuccr(k) + prds(k))*cldmax(i,k))*   &
                             rho(i,k)*dz(i,k) + qstot .lt. 0._mg_pr) then
                 if (-prds(k) .ge. qsmall) then
                   ratio = (qstot/(rho(i,k)*dz(i,k)) + (bergs(k) +   &
                             psacws(k) + prai(k) + prci(k))*cldm(i,k) +  &
                             (pracs(k) + mnuccr(k))*cldmax(i,k))/   &
                                            (-prds(k)*cldmax(i,k))*omsm
                 else
                   ratio =0._mg_pr
                 end if
                 prds(k) = prds(k)*ratio
               end if

! conservation of ns


! calculate loss of number due to sublimation
               IF (subl_snow) THEN
                 if (prds(k) .lt. 0._mg_pr .and.     &
                                             qniic(i,k) .ge. qsmall) then
                   nsubs(k) = prds(k)/qniic(i,k)*nsic(i,k)
                 else
                   nsubs(k) = 0._mg_pr
                 end if

! neglect sublimation of ns
               ELSE
                 nsubs(k) = 0._mg_pr
               END IF 

               if ((nprci(k)*cldm(i,k) + (nnuccr(k) + nsubs(k) +   &
                                                nsagg(k))*cldmax(i,k))*&
                            rho (i,k)*dz(i,k) + nstot .lt. 0._mg_pr) then
                 if (-nsubs(k) - nsagg(k) .ge. qsmall) then
                   ratio = (nstot/(rho(i,k)*dz(i,k)) + nprci(k)*   &
                                                           cldm(i,k) +   &
                                           nnuccr(k)*cldmax(i,k))/    &
                                 ((-nsubs(k) - nsagg(k))*cldmax(i,k))*omsm
                 else 
                   ratio =0._mg_pr
                 end if
                 nsubs(k) = nsubs(k)*ratio
                 nsagg(k) = nsagg(k)*ratio
               end if

! get tendencies due to microphysical conversion processes
! note: tendencies are multiplied by appropaiate cloud/precip 
! fraction to get grid-scale values
! note: cmei,cmel are already grid-average values

               qvlat(i,k) = qvlat(i,k) - &
                             (pre(k) + prds(k))*cldmax(i,k) - cmel(i,k) - &
                              cmei(i,k) - D_eros_l(i,k) - D_eros_i(i,k) - &
                                                             qvdep_qi(i,k)
               tlat(i,k) = tlat(i,k) + ((pre(k)*cldmax(i,k) +    &
                            cmel(i,k) + D_eros_l(i,k))*xxlv +   &
                            (prds(k)*cldmax(i,k) + cmei(i,k) +    &
                             D_eros_i(i,k) + qvdep_qi(i,k))*xxls +    &
                            ((bergs(k) + psacws(k) + psacws_o(k) +   &
                                mnuccc(k))*cldm(i,k) + (mnuccr(k) + &
                                  pracs(k))*cldmax(i,k) + berg(i,k))*xlf)
               qctend(i,k) = qctend(i,k) + &
                              (-pra(k) - prc(k) - mnuccc(k) -    &
                               psacws(k) - psacws_o(k) - bergs(k))*   &
                                                          cldm(i,k) +    &
                              cmel(i,k) - berg(i,k) + D_eros_l(i,k)
               qitend(i,k) = qitend(i,k) +       &
                              (mnuccc(k) - prci(k) - prai(k) +    &
                              psacws_o(k))*cldm(i,k) + cmei(i,k) +   &
                                 berg(i,k) + D_eros_i(i,k) + qvdep_qi(i,k)
               qrtend(i,k) = qrtend(i,k) + &
                              (pra(k) + prc(k))*cldm(i,k) + (pre(k) -   &
                                         pracs(k) - mnuccr(k))*cldmax(i,k)
               qnitend(i,k) = qnitend(i,k) + &
                               (prai(k) + psacws(k) + prci(k) +    &
                                          bergs(k))*cldm(i,k) +    &
                                       (prds(k) + pracs(k) + mnuccr(k))*  &
                                                                cldmax(i,k)

! multiply activation/nucleation by mtime to account for fast timescale

               dumd = nctend(i,k) 
               nctend(i,k) = nctend(i,k) + npccn(k)*mtime +   &
                             (-nnuccc(k) - npsacws(k) - npsacws_o(k) +  &
                                 nsubc(k) + nerosc(i,k) - npra(k) -   &
                                                        nprc1(k))*cldm(i,k)
               nitend(i,k) = nitend(i,k) + nnuccd(k)*mtime +    &
                              (nsubi(k) + nerosi(i,k) - nprci(k) - &
                                                       nprai(k))*cldm(i,k)
               nstend(i,k) = nstend(i,k) + (nsubs(k) +       &
                              nsagg(k) + nnuccr(k))*cldmax(i,k) +  &
                                                        nprci(k)*cldm(i,k)
               nrtend(i,k) = nrtend(i,k) +       &
                              nprc(k)*cldm(i,k) + (nsubr(k) - npracs(k) - &
                                        nnuccr(k) + nragg(k))*cldmax(i,k)

! make sure that nc and ni at advanced time step do not exceed
! maximum (existing N + source terms*dt), which is possible due to
! fast nucleation timescale

! diag++
               IF (diag_id%qndt_nucclim +     &
                                     diag_id%qn_nucclim_col  > 0) THEN
                 nucclim(k) = nctend(i,k)
               END IF
               IF (diag_id%qnidt_nucclim1 +    &
                               diag_id%qni_nucclim1_col > 0 ) THEN
                 nucclim1i(k) = nitend(i,k)
               END IF
! diag--
    
               if (nctend(i,k) .gt. 0._mg_pr .and.    &
                       nc(i,k) + nctend(i,k)*deltat .gt. ncmax) then
                 nctend(i,k) = max(0._mg_pr, (ncmax - nc(i,k))/deltat)
               end if
               if (nitend(i,k) .gt. 0._mg_pr .and.    &
                       ni(i,k) + nitend(i,k)*deltat .gt. nimax) then
                 nitend(i,k) = max(0._mg_pr, (nimax - ni(i,k))/deltat)
               end if

! diag++
               IF (diag_id%qndt_nucclim +     &
                                    diag_id%qn_nucclim_col  > 0) THEN
                 nucclim(k) = nctend(i,k) - nucclim(k)
               END IF
               IF (diag_id%qnidt_nucclim1 +   &
                                diag_id%qni_nucclim1_col > 0) THEN
                 nucclim1i(k) = nitend(i,k) - nucclim1i(k)
               END IF
! diag--

! cms 2009-2-26 also limit volume mean ice radius (optionally)
               IF (limit_volri) THEN 

                 IF (diag_id%qnidt_nucclim2 +    &
                             diag_id%qni_nucclim2_col > 0) THEN
                   nucclim2(k) = nitend(i,k)
                 END IF

                 qii_new = (qi(i,k) + qitend(i,k)*deltat)/cldm(i,k)
                 nii_new = (ni(i,k) + nitend(i,k)*deltat)/cldm(i,k)

!max XXX micron 
                 nii_min = qii_new/rhoi*3._mg_pr/(4._mg_pr*3.14_mg_pr*  &
                                                            max_diam_ii**3)
!min XXX micron
                 nii_max = qii_new/rhoi*3._mg_pr/(4._mg_pr*3.14_mg_pr*  &
                                                            min_diam_ii**3)
                 if ( nii_new .gt. nii_max) then
                   nitend(i,k) = (nii_max - ni(i,k)/cldm(i,k))/deltat*   &
                                                                 cldm(i,k)
                 else if (nii_new .lt. nii_min) then
                   nitend(i,k) = (nii_min - ni(i,k)/cldm(i,k))/deltat*   &
                                                                  cldm(i,k)
                 end if
                 IF (diag_id%qnidt_nucclim2 +     &
                                    diag_id%qni_nucclim2_col > 0) THEN
                   nucclim2(k) = nitend(i,k) - nucclim2(k)
                 END IF
               ELSE
                   nucclim2(k) = 0.
               END IF

! get final values for precipitation q and N, based on
! flux of precip from above, source/sink term, and terminal fallspeed
! see eq. 15-16 in Morrison and Gettelman, 2007, J. Climate

! rain
               if (qric(i,k) .ge. qsmall) then
                 if (k .eq. 1) then
                   qric(i,k) = qrtend(i,k)*dz(i,k)/cldmax(i,k)/umr(k)
                   nric(i,k) = nrtend(i,k)*dz(i,k)/cldmax(i,k)/unr(k)
                 else
                   qric(i,k) = (rho(i,k-1)*umr(k-1)*qric(i,k-1)*   &
                                    cldmax(i,k-1) + (rho(i,k)*dz(i,k)*  &
                                             qrtend(i,k)))/(umr(k)*  &
                                                     rho(i,k)*cldmax(i,k))
                   nric(i,k) = (rho(i,k-1)*unr(k-1)*nric(i,k-1)*    &
                                     cldmax(i,k-1) + (rho(i,k)*dz(i,k)*  &
                                             nrtend(i,k)))/(unr(k)*    &
                                                     rho(i,k)*cldmax(i,k))
                 end if
               else
                 qric(i,k) = 0._mg_pr
                 nric(i,k) = 0._mg_pr
               end if

! snow

               if (qniic(i,k) .ge. qsmall) then
                 if (k .eq. 1) then
                   qniic(i,k) = qnitend(i,k)*dz(i,k)/cldmax(i,k)/ums(k)
                   nsic(i,k) = nstend(i,k)*dz(i,k)/cldmax(i,k)/uns(k)
                 else
                   qniic(i,k) = (rho(i,k-1)*ums(k-1)*qniic(i,k-1)*   &
                                     cldmax(i,k-1) + (rho(i,k)*dz(i,k)*  &
                                          qnitend(i,k)))/(ums(k)*   &
                                                    rho(i,k)*cldmax(i,k))
                   nsic(i,k) = (rho(i,k-1)*uns(k-1)*nsic(i,k-1)*    &
                                      cldmax(i,k-1) + (rho(i,k)*dz(i,k)*  &
                                          nstend(i,k)))/(uns(k)*   &
                                                     rho(i,k)*cldmax(i,k))
                 end if
               else
                 qniic(i,k) = 0._mg_pr
                 nsic(i,k) = 0._mg_pr
               end if

! calculate precipitation flux at surface
! divide by density of water to get units of m/s

               prect(i) = prect(i) + (qrtend(i,k)*rho (i,k)*dz(i,k) +  &
                                 qnitend(i,k)*rho(i,k)*dz(i,k))/rhow
               preci(i) = preci(i) + qnitend(i,k)*rho(i,k)*dz(i,k)/rhow

! hm, add 9/5/07
! convert rain rate from m/s to mm/hr
 
               rainrt(i,k) = qric(i,k)*rho(i,k)*umr(k)/     &
                                              rhow*3600._mg_pr*1000._mg_pr

! vertically-integrated precip source/sink terms (note: grid-averaged)

!              qrtot = max(qrtot + qrtend(i,k)*rho(i,k)*dz(i,k), 0._mg_pr)
!              qstot = max(qstot + qnitend(i,k)*rho(i,k)*dz(i,k),    &
!                                                               0._mg_pr)
!              nrtot = max(nrtot + nrtend(i,k)*rho(i,k)*dz(i,k), 0._mg_pr)
!              nstot = max(nstot + nstend(i,k)*rho(i,k)*dz(i,k), 0._mg_pr)
               qrtot = qrtot + qrtend(i,k)*rho(i,k)*dz(i,k)
               qstot = qstot + qnitend(i,k)*rho(i,k)*dz(i,k)
               nrtot = nrtot + nrtend(i,k)*rho(i,k)*dz(i,k)
               nstot = nstot + nstend(i,k)*rho(i,k)*dz(i,k)

! calculate melting and freezing of precip

! melt snow at +2 C
! NOTE RSH 11/22/11:
! NOTE THAT R-K only starts melting at 0 C -- called a bug to melt at +2 C 
!
               if (t(i,k) + tlat(i,k)/cpp*deltat >     &
                                                 tfreeze + 2._mg_pr) then
                 if (qstot > 0._mg_pr) then

! make sure melting snow doesn't reduce temperature below threshold
                   dum = -xlf/cpp*qstot/(rho(i,k)*dz(i,k))
                   if (t(i,k) + tlat(i,k)/cpp*deltat+dum .lt.    &
                                                 tfreeze + 2._mg_pr) then
                     dum = (t(i,k) + tlat(i,k)/cpp*deltat -   &
                                            (tfreeze + 2._mg_pr))*cpp/xlf
                     dum = dum/(xlf/cpp*qstot/(rho(i,k)*dz(i,k)))
                     dum = max(0._mg_pr, dum)
                     dum = min(1._mg_pr, dum)
                   else
                     dum = 1._mg_pr
                   end if
                   qric(i,k) = qric(i,k) + dum*qniic(i,k)
                   nric(i,k) = nric(i,k) + dum*nsic(i,k)
                   qniic(i,k) = (1._mg_pr - dum)*qniic(i,k)
                   nsic(i,k) = (1._mg_pr - dum)*nsic(i,k)
                   tlat(i,k) = tlat(i,k) - xlf*dum*qstot/    &
                                                       (rho(i,k)*dz(i,k))
                   qrtot = qrtot + dum*qstot
                   qstot = (1._mg_pr - dum)*qstot
                   nrtot = nrtot + dum*nstot
                   nstot = (1._mg_pr - dum)*nstot
                   if (diag_id%snow_melt + diag_id%snow_melt_col > 0) & 
                             diag_4d(i,j,k, diag_pt%snow_melt) =     &
                                    diag_4d(i,j,k, diag_pt%snow_melt) +  &
                                       dum*preci(i)*rhow/(rho(i,k)*dz(i,k))
                   preci(i) = (1._mg_pr - dum)*preci(i)

!cms++
! assume that droplets which would condense due to cooling associated 
! with melting of snow are rapidly collected by (newly formed) rain drops
! (this is optional, controlled by namelist, set to .false. in Marc's final
! parameterization)

                   IF (collect_frzreg) THEN
!only if there is net cooling
                     IF (tlat(i,k) .LT. 0._mg_pr) THEN 
                       ttmp = t(i,k) + tlat(i,k)/cpp*deltat
                       qtmp = q(i,k) + qvlat(i,k)*deltat
                       esn = polysvp_l(ttmp)
                       qvs = 0.622_mg_pr*esn/(pfull(i,k) - d378*esn)
                       dqsdt = xxlv*qvs/(rv*ttmp**2)
                       ab = 1._mg_pr + dqsdt*xxlv/cpp
                       tmp2 = pfull(i,k) - d378*esn
                       tmp2 = max(tmp2,esn)
                       tmp2 = d622*esn/tmp2
                       tmp1 = max(0._mg_pr, (qtmp - tmp2)/ab)
        
! change vapor content and T
                       qvlat(i,k) = qvlat(i,k) -tmp1/deltat 
                       snow2vapor(k) = -tmp1/deltat
                       tlat(i,k)  = tlat(i,k) + hlv*tmp1/deltat

! and rain mixing ratio
                       qric(i,k) = qric(i,k) + tmp1/cldmax(i,k)
                       qrtot = qrtot + tmp1/deltat*pdel(i,k)/grav
                       prect(i) = prect(i) + tmp1/deltat*pdel(i,k)/   &
                                                                grav/rhow 
                     END IF 
                   END IF 
                 end if
               end if

! freeze rain at -5 C
               if (t(i,k) + tlat(i,k)/cpp*deltat <      &
                                                 tfreeze - 5._mg_pr ) then
                 if (qrtot > 0._mg_pr) then

! make sure freezing rain doesn't increase temperature above threshold
                   dum = xlf/cpp*qrtot/(rho(i,k)*dz(i,k))
                   if (t(i,k) + tlat(i,k)/cpp*deltat + dum   &
                                            .gt. tfreeze -5._mg_pr ) then
                     dum = -(t(i,k) + tlat(i,k)/cpp*deltat -    &
                                             (tfreeze - 5._mg_pr) )*cpp/xlf
                     dum = dum/(xlf/cpp*qrtot/(rho(i,k)*dz(i,k)))
                     dum = max(0._mg_pr, dum)
                     dum = min(1._mg_pr, dum)
                   else
                     dum = 1._mg_pr
                   end if
                   qniic(i,k) = qniic(i,k) + dum*qric(i,k)
                   nsic(i,k) = nsic(i,k) + dum*nric(i,k)
                   qric(i,k) = (1._mg_pr - dum)*qric(i,k)
                   nric(i,k) = (1._mg_pr - dum)*nric(i,k)
                   tlat(i,k) = tlat(i,k) + xlf*dum*qrtot/   &
                                                       (rho(i,k)*dz(i,k))
                   qstot = qstot + dum*qrtot
                   qrtot = (1._mg_pr - dum)*qrtot
                   nstot = nstot + dum*nrtot
                   nrtot = (1._mg_pr - dum)*nrtot
                   diag_4d(i,j,k, diag_pt%rain_freeze) =  &
                           diag_4d(i,j,k, diag_pt%rain_freeze) +  &
                                          dum*(prect(i) - preci(i))* &
                                                    rhow/(rho(i,k)*dz(i,k))
                   preci(i) = preci(i) + dum* (prect(i) - preci(i))
                 end if
               end if

! if rain/snow mix ratio is zero so should number concentration
               if (qniic(i,k) .lt. qsmall) then
                 qniic(i,k) = 0._mg_pr
                 nsic(i,k) = 0._mg_pr
               end if
               if (qric(i,k) .lt. qsmall) then
                 qric(i,k) = 0._mg_pr
                 nric(i,k) = 0._mg_pr
               end if

! make sure number concentration is a positive number to avoid 
! taking root of negative
               nric(i,k) = max(nric(i,k), 0._mg_pr)
               nsic(i,k) = max(nsic(i,k), 0._mg_pr)

!.......................................................................
! get size distribution parameters for fallspeed calculations
!......................................................................
! rain
               if (qric(i,k) .ge. qsmall) then
                 lamr(k) = (pi*rhow*nric(i,k)/    &
                                          qric(i,k))**(1._mg_pr/3._mg_pr)
                 n0r(k) = nric(i,k)*lamr(k)

! check for slope
! hm 4/5/07, change lammax and lammin for rain and snow
                 lammax = 1._mg_pr/20.e-6_mg_pr
                 lammin = 1._mg_pr/500.e-6_mg_pr

! adjust vars
                 if (lamr(k) .lt. lammin) then
                   lamr(k) = lammin
                   n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
                   nric(i,k) = n0r(k)/lamr(k)
                 else if (lamr(k) .gt. lammax) then
                   lamr(k) = lammax
                   n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
                   nric(i,k) = n0r(k)/lamr(k)
                 end if

! 'final' values of number and mass weighted mean fallspeed for rain (m/s)
                 unr(k) = min(arn(i,k)*gamma_mg(1._mg_pr + br)/    &
                                          lamr(k)**br, 9.1_mg_pr*rhof(i,k))
                 umr(k) = min(arn(i,k)*gamma_mg(4._mg_pr+br)/    &
                               (6._mg_pr*lamr(k)**br), 9.1_mg_pr*rhof(i,k))
               else
                 lamr(k) = 0._mg_pr
                 n0r(k) = 0._mg_pr
                 umr(k) = 0._mg_pr
                 unr(k) = 0._mg_pr
               end if

!......................................................................
! snow
               if (qniic(i,k) .ge. qsmall) then
                 lams(k) = (gamma_mg(1._mg_pr + ds)*cs*nsic(i,k)/ &
                                            qniic(i,k))**(1._mg_pr/ds)
                 n0s(k) = nsic(i,k)*lams(k)

! check for slope
                 lammax = 1._mg_pr/min_diam_ice
                 lammin = 1._mg_pr/2000.e-6_mg_pr

! adjust vars
                 if (lams(k) .lt. lammin) then
                   lams(k) = lammin
                   n0s(k) = lams(k)**(ds + 1._mg_pr)*qniic(i,k)/   &
                                             (cs*gamma_mg(1._mg_pr + ds))
                   nsic(i,k) = n0s(k)/lams(k)
                 else if (lams(k).gt.lammax) then
                   lams(k) = lammax
                   n0s(k) = lams(k)**(ds + 1._mg_pr)*qniic(i,k)/    &
                                             (cs*gamma_mg(1._mg_pr + ds))
                   nsic(i,k) = n0s(k)/lams(k)
                 end if

! 'final' values of number and mass weighted mean fallspeed for snow (m/s)
                 ums(k) = min(asn(i,k)*gamma_mg(4._mg_pr + bs)/    &
                             (6._mg_pr*lams(k)**bs), max_vt_snow*rhof(i,k))
                 uns(k) = min(asn(i,k)*gamma_mg(1._mg_pr + bs)/     &
                                       lams(k)**bs, max_vt_snow*rhof(i,k))
               else
                 lams(k) = 0._mg_pr
                 n0s(k) = 0._mg_pr
                 ums(k) = 0._mg_pr
                 uns(k) = 0._mg_pr
               end if

! convert rain/snow q and N for output to history, note, 
! output is for gridbox average

               qrout(i,k) = qrout(i,k) + qric(i,k)*cldmax(i,k)
               qsout(i,k) = qsout(i,k) + qniic(i,k)*cldmax(i,k)
               nrout(i,k) = nrout(i,k) + nric(i,k)*rho(i,k)*cldmax(i,k)
               nsout(i,k) = nsout(i,k) + nsic(i,k)*rho(i,k)*cldmax(i,k)

!c........................................................................
! sum over sub-step for average process rates

               tlat1(i,k) = tlat1(i,k) + tlat(i,k)
               qvlat1(i,k) = qvlat1(i,k) + qvlat(i,k)
               qctend1(i,k) = qctend1(i,k) + qctend(i,k)
               qitend1(i,k) = qitend1(i,k) + qitend(i,k)
               nctend1(i,k) = nctend1(i,k) + nctend(i,k)
               nitend1(i,k) = nitend1(i,k) + nitend(i,k)

               t(i,k) = t(i,k) + tlat(i,k)*deltat/cpp
               q(i,k) = q(i,k) + qvlat(i,k)*deltat
               qc(i,k) = qc(i,k) + qctend(i,k)*deltat
               qi(i,k) = qi(i,k) + qitend(i,k)*deltat
               nc(i,k) = nc(i,k) + nctend(i,k)*deltat
               ni(i,k) = ni(i,k) + nitend(i,k)*deltat

! hm add 9/5/07
               rainrt1(i,k) = rainrt1(i,k) + rainrt(i,k)

!CODE MOVED TO HERE:
               if ( k > 1 ) then        
                 IF (.NOT. one_ice) THEN
                   asnowrt(i,k-1) = cldmax(i,k-1)*qniic(i,k-1)*   &
                                                      rho(i,k-1)*ums(k-1)
                   if (k == kdim) then
                     asnowrt(i,kdim) = cldmax(i,kdim)*qniic(i,kdim)*   &
                                                    rho(i,kdim)*ums(kdim)
                   endif
                 ELSE
                   asnowrt(i,k-1)= cldmax(i,k-1)*qiic(i,k-1)*   &
                                                      rho(i,k-1)*ums(k-1)
                   if (k == kdim) then
                     asnowrt(i,kdim) = cldmax(i,kdim)*qiic(i,kdim)*   &
                                                     rho(i,kdim)*ums(kdim)
                   endif
                 END IF
                 atotrt(i,k-1) =  asnowrt(i,k-1) + cldmax(i,k-1)*  &
                                           qric(i,k-1)*rho(i,k-1)*umr(k-1)
                 if (k == kdim) then
                   atotrt(i,kdim) =  asnowrt(i,kdim) + cldmax(i,kdim)*  &
                                         qric(i,kdim)*rho(i,kdim)*umr(kdim)
                 endif
                 atotrt1(i,k-1) = atotrt1(i,k-1) + atotrt(i,k-1)
                 asnowrt1(i,k-1) = asnowrt1(i,k-1) + asnowrt(i,k-1)   
               end if
!END OF MOVED CODE

               if (k == kdim) then
                 atotrt1(i,k) = atotrt1(i,k) + atotrt(i,k)
                 asnowrt1(i,k) = asnowrt1(i,k) + asnowrt(i,k)   
               endif

               pre1(i,k)  = pre1(i,k) - pre(k)*cldmax(i,k)
               prds1(i,k) = prds1(i,k) - prds(k)*cldmax(i,k)
               snow2vapor1(i,k) = snow2vapor1(i,k) + snow2vapor(k)
               cmel1(i,k) = cmel1(i,k) + cmel(i,k)
               D_eros_l1(i,k) = D_eros_l1(i,k) + D_eros_l(i,k)
               berg1(i,k) = berg1(i,k) + berg(i,k)
               qvdep_qi1(i,k) = qvdep_qi1(i,k) + qvdep_qi(i,k)
               prc1(i,k) = prc1(i,k) - prc(k)*cldm(i,k)

               pra1(i,k) = pra1(i,k) - pra(k)*cldm(i,k)
               mnuccc1(i,k) = mnuccc1(i,k) - mnuccc(k)*cldm(i,k)

               psacws1(i,k) = psacws1(i,k) -    &
                                        (psacws(k) + psacws_o(k))*cldm(i,k)
               psacws_o1(i,k) =  psacws_o1(i,k) - psacws_o(k)*cldm(i,k)
               bergs1(i,k) = bergs1(i,k) - bergs(k)*cldm(i,k)
   
               cmei1(i,k) = cmei1(i,k) + cmei(i,k)
               D_eros_i1(i,k) = D_eros_i1(i,k) + D_eros_i(i,k)
               prci1(i,k) = prci1(i,k) - prci(k)*cldm(i,k)
               prai1(i,k) = prai1(i,k) - prai(k)*cldm(i,k)

!droplet number
               npccn1(i,k) = npccn1(i,k) + npccn(k)*mtime
               nnuccc1(i,k) = nnuccc1(i,k) - nnuccc(k)*cldm(i,k)
               npsacws1(i,k) = npsacws1(i,k) - npsacws(k)*cldm(i,k)
               npsacws_o1(i,k) = npsacws_o1(i,k) - npsacws_o(k)*cldm(i,k)
               nsubc1(i,k) = nsubc1(i,k) + nsubc(k)*cldm(i,k)
               nerosc1(i,k) = nerosc1(i,k) + nerosc(i,k)*cldm(i,k)
               npra1(i,k) = npra1(i,k) - npra(k)*cldm(i,k)
               nprc11(i,k) = nprc11(i,k) - nprc1(k)*cldm(i,k)
               nucclim1(i,k) = nucclim1(i,k) + nucclim(k)

!ice number
               nnuccd1(i,k) = nnuccd1(i,k) + nnuccd(k)*mtime
               nsubi1(i,k) = nsubi1(i,k) + nsubi(k)*cldm(i,k)
               nerosi1(i,k) = nerosi1(i,k) + nerosi(i,k)*cldm(i,k)
               nprci1(i,k) = nprci1(i,k) -  nprci(k)*cldm(i,k)
               nprai1(i,k) = nprai1(i,k) - nprai(k)*cldm(i,k)
               nucclim1_1(i,k) = nucclim1_1(i,k) + nucclim1i(k)
               nucclim2_1(i,k) = nucclim2_1(i,k) + nucclim2(k)

               pracs1(i,k) = pracs1(i,k) - pracs(k)*cldmax(i,k)
               mnuccr1(i,k) = mnuccr1(i,k) - mnuccr(k)*cldmax(i,k)
             end do   !  k-loop (do k=1,kdim; starts ~ 1800 lines above)

             prect1(i) = prect1(i) + prect(i)
             preci1(i) = preci1(i) + preci(i)
           end do ! it loop  (do it=1,iter) (sub-step loop)
300        continue  ! skip to end of loop if  no cloud water
         end do  ! i loop (do i=1,idim)

! convert dt from sub-step back to full time step
         deltat=deltat*real(iter)

!c........................................................................

!ADD HERE 6/6/12
         ssat_disposal(:,:) = 0._mg_pr

         do i=1,idim

! skip all calculations if no cloud water
           if (ltrue(i) .eq. 0) then
             if (diag_id%vfall > 0)     &
                               diag_4d(i,j,:,diag_pt%vfall) = 0.0_mg_pr   
             goto 500
           endif
  
! initialize nstep for sedimentation sub-steps
           nstep = 1

! divide precip rate by number of sub-steps to get average over time step
           prect(i) = prect1(i)/real(iter)
           preci(i) = preci1(i)/real(iter)

           diag_4d(i,j,:, diag_pt%snow_melt) =  &
                         diag_4d(i,j,:, diag_pt%snow_melt)/real(iter)
           diag_4d(i,j,:, diag_pt%rain_freeze) =  &
                        diag_4d(i,j,:, diag_pt%rain_freeze)/real(iter)
           do k=1,kdim
             umi(k) = 0._mg_pr 
           end do
 
           do k=1,kdim  

! assign variables back to start-of-timestep values before updating 
! after sub-steps 
             t(i,k) = t1(i,k)
             q(i,k) = q1(i,k)
             qc(i,k) = qc1(i,k)
             qi(i,k) = qi1(i,k)
             nc(i,k) = nc1(i,k)
             ni(i,k) = ni1(i,k)

! divide microphysical tendencies by number of sub-steps to get average 
! over time step
             tlat(i,k) = tlat1(i,k)/real(iter)
             qvlat(i,k) = qvlat1(i,k)/real(iter)
             qctend(i,k) = qctend1(i,k)/real(iter)
             qitend(i,k) = qitend1(i,k)/real(iter)
             nctend(i,k) = nctend1(i,k)/real(iter)
             nitend(i,k) = nitend1(i,k)/real(iter)
   
! hm, add 9/5/07
             rainrt(i,k) = rainrt1(i,k)/real(iter)

             atotrt(i,k) = atotrt1(i,k)/real(iter)
             asnowrt(i,k) = asnowrt1(i,k)/real(iter)

! divide output precip q and N by number of sub-steps to get average 
! over time step
             qrout(i,k) = qrout(i,k)/real(iter)
             qsout(i,k) = qsout(i,k)/real(iter)
             nrout(i,k) = nrout(i,k)/real(iter)
             nsout(i,k) = nsout(i,k)/real(iter)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate sedimentation for cloud water and ice
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! update in-cloud cloud mixing ratio and number concentration 
! with microphsical tendencies to calculate sedimentation, assign 
! to dummy vars
! note: these are in-cloud values***, hence we divide by cloud fraction

             dumc(i,k) = (qc(i,k) + qctend(i,k)*deltat)/cldm(i,k)
             dumi(i,k) = (qi(i,k) + qitend(i,k)*deltat)/cldm(i,k)
             dumnc(i,k) = max((nc(i,k) + nctend(i,k)*deltat)/cldm(i,k), &
                                                                 0._mg_pr)
             dumni(i,k) = max((ni(i,k) + nitend(i,k)*deltat)/cldm(i,k),  &
                                                                 0._mg_pr)

! obtain new slope parameter to avoid possible singularity
             if (dumi(i,k) .ge. qsmall) then

! add upper limit to in-cloud number concentration to prevent 
! numerical error
               dumni(i,k) = min(dumni(i,k), dumi(i,k)*1.e20_mg_pr)

               lami(k) = (gamma_mg(1._mg_pr + di_mg)*ci_mg* &
                                 dumni(i,k)/dumi(i,k))**(1._mg_pr/di_mg)
               lammax = 1._mg_pr/min_diam_ice
               lammin = 1._mg_pr/(2._mg_pr*dcs)
               lami(k) = max(lami(k), lammin)
               lami(k) = min(lami(k), lammax)
             else
               lami(k) = 0._mg_pr
             end if

             if (dumc(i,k) .ge. qsmall) then

! add upper limit to in-cloud number concentration to prevent 
! numerical error
               dumnc(i,k) = min(dumnc(i,k), dumc(i,k)*1.e20_mg_pr)

!RSH BUGFIX email of 6/8/10
!              pgam(k)=0.0005714_mg_pr*(dumnc(i,k)/1.e6_mg_pr/rho(i,k))+  &
!                                                              0.2714_mg_pr
               pgam(k) = 0.0005714_mg_pr*(dumnc(i,k)/1.e6_mg_pr*   &
                                                   rho(i,k))+0.2714_mg_pr
               pgam(k) = 1._mg_pr/(pgam(k)**2) - 1._mg_pr
               pgam(k) = max(pgam(k), 2._mg_pr)
               pgam(k) = min(pgam(k), 15._mg_pr)

               lamc(k) = (pi/6._mg_pr*rhow*dumnc(i,k)*    &
                                         gamma_mg(pgam(k) + 4._mg_pr)/ &
                                                          (dumc(i,k)*    &
                        gamma_mg(pgam(k) + 1._mg_pr)))**(1._mg_pr/3._mg_pr)
               lammin = (pgam(k) + 1._mg_pr)/max_diam_drop
               lammax = (pgam(k) + 1._mg_pr)/min_diam_drop
               lamc(k) = max(lamc(k), lammin)
               lamc(k) = min(lamc(k), lammax)
             else
               lamc(k) = 0._mg_pr
             end if

! calculate number and mass weighted fall velocity for droplets
! include effects of sub-grid distribution of cloud water
             if (dumc(i,k) .ge. qsmall) then

!RSH bugfix email 6/8/10
!              unc= sfac5 * &
               unc =         &
                     acn(i,k)*gamma_mg(1._mg_pr + bc+pgam(k))/ &
                               (lamc(k)**bc*gamma_mg(pgam(k) + 1._mg_pr))
!RSH bugfix email 6/8/10
!              umc =   sfac5 * &
               umc =           &
                          acn(i,k)*gamma_mg(4._mg_pr + bc+pgam(k))/ &
                               (lamc(k)**bc*gamma_mg(pgam(k) + 4._mg_pr))
             else
               umc = 0._mg_pr
               unc = 0._mg_pr
             end if

! calculate number and mass weighted fall velocity for cloud ice
             if (dumi(i,k) .ge. qsmall) then
               IF (.NOT. one_ice) THEN
                 uni = ain(i,k)*gamma_mg(1._mg_pr + bi)/lami(k)**bi
                 umi(k) = ain(i,k)*gamma_mg(4._mg_pr + bi)/    &
                                                   (6._mg_pr*lami(k)**bi)
               ELSE
                 umi(k) = min(asn(i,k)*gamma_mg(4._mg_pr + bs)/   &
                                       (6._mg_pr*lami(k)**bs), max_vt_ice)
                 uni = min(asn(i,k)*gamma_mg(1._mg_pr + bs)/   &
                                                   lami(k)**bs, max_vt_ice)
               END IF
               uni = vfact_n*uni
               umi(k) = vfact_m*umi(k)
               uni = min(uni, max_vt_ice*rhof(i,k))
               umi(k) = min(umi(k), max_vt_ice*rhof(i,k))

               IF (hd_sedi_sens) THEN
                 umi(k) = vfact*3.29_mg_pr*    &
                             ((rho(i,k)*qi(i,k)/cldm(i,k))**0.16_mg_pr)
                 uni = umi(k)
               END IF

               IF (scav_by_cloud_ice) THEN
                 if ( k > 1 ) then  
!RSH: if this activated, should deal with difference between rhoi and rhosn
!RSH  in this additional term in asnowrt :
                   asnowrt(i,k) = asnowrt(i,k) + cldmax(i,k-1)*    &
                                          qiic(i,k-1)*rho(i,k-1)*umi(k-1)
                 end if
               END IF
             else
               umi(k) = 0._mg_pr
               uni = 0._mg_pr
             end if

             if (diag_id%vfall > 0) diag_4d(i,j,k,diag_pt%vfall) = umi(k)

             fi(k) = grav*rho(i,k)*umi(k)
             fni(k) = grav*rho(i,k)*uni
             fc(k) = grav*rho(i,k)*umc
             fnc(k) = grav*rho(i,k)*unc

! calculate number of split time steps to ensure courant stability criteria
! for sedimentation calculations

             rgvm = max(fi(k), fc(k), fni(k), fnc(k))
             nstep = max(int(rgvm*deltat/pdel(i,k) + 1._mg_pr), nstep)

! redefine dummy variables - sedimentation is calculated over grid-scale
! quantities to ensure conservation
             dumc(i,k) = (qc(i,k) + qctend(i,k)*deltat)
             dumi(i,k) = (qi(i,k) + qitend (i,k)*deltat)
             dumnc(i,k) = max((nc(i,k) + nctend(i,k)*deltat), 0._mg_pr)
             dumni(i,k) = max((ni(i,k) + nitend(i,k)*deltat), 0._mg_pr)

             if (dumc(i,k) .lt. qsmall) dumnc(i,k) = 0._mg_pr
             if (dumi(i,k) .lt. qsmall) dumni(i,k) = 0._mg_pr
           end do       !!! vertical loop

           do n = 1,nstep  !! loop over sub-time step to ensure stability
             do k = 1,kdim
               falouti(k) = fi(k)*dumi(i,k)
               faloutni(k) = fni(k)*dumni(i,k)
               faloutc(k) = fc(k)*dumc(i,k)
               faloutnc(k) = fnc(k)*dumnc(i,k)
             end do

! top of model
             k = 1
             faltndi = falouti(k)/pdel(i,k)
             faltndni = faloutni(k)/pdel(i,k)
             faltndc = faloutc(k)/pdel(i,k)
             faltndnc = faloutnc(k)/pdel(i,k)

! add fallout terms to microphysical tendencies
             qitend(i,k) = qitend(i,k) - faltndi/nstep
             nitend(i,k) = nitend(i,k) - faltndni/nstep
             qctend(i,k) = qctend(i,k) - faltndc/nstep
             nctend(i,k) = nctend(i,k) - faltndnc/nstep
             IF (diag_id%qldt_sedi + diag_id%ql_sedi_col > 0) &
                    diag_4d(i,j,k,diag_pt%qldt_sedi) =    &
                          diag_4d(i,j,k,diag_pt%qldt_sedi) - faltndc/nstep
             IF (diag_id%sedi_ice > 0) &
                    diag_4d(i,j,1,diag_pt%sedi_ice) =    &
                                 diag_4d(i,j,1,diag_pt%sedi_ice) +       &
                                             falouti(kdim)/grav/nstep/rhoi
             IF (diag_id%qidt_fall + diag_id%qi_fall_col > 0) &
                    diag_4d(i,j,k,diag_pt%qidt_fall) =     &
                        diag_4d(i,j,k,diag_pt%qidt_fall) - faltndi/nstep
             IF (diag_id%qndt_sedi + diag_id%qn_sedi_col > 0) &
                    diag_4d(i,j,k,diag_pt%qndt_sedi) =    &
                        diag_4d(i,j,k,diag_pt%qndt_sedi) - faltndnc/nstep
             IF (diag_id%qnidt_sedi +  diag_id%qni_sedi_col > 0) &
                     diag_4d(i,j,k,diag_pt%qnidt_sedi) =    &
                       diag_4d(i,j,k,diag_pt%qnidt_sedi) - faltndni/nstep

             dumi(i,k) = dumi(i,k) - faltndi*deltat/nstep
             dumni(i,k) = dumni(i,k) - faltndni*deltat/nstep
             dumc(i,k) = dumc(i,k) - faltndc*deltat/nstep
             dumnc(i,k) = dumnc(i,k) - faltndnc*deltat/nstep

             do k = 2,kdim
 
! for cloud liquid and ice, if cloud fraction increases with height
! then add flux from above to both vapor and cloud water of current level
! this means that flux entering clear portion of cell from above evaporates
! instantly

               dum = cldm(i,k)/cldm(i,k-1)
               dum = min(dum, 1._mg_pr)

               faltndqie = (falouti(k) - falouti(k-1))/pdel(i,k)
               faltndi = (falouti(k) - dum*falouti(k-1))/pdel(i,k)
               faltndni = (faloutni(k) - dum*faloutni(k-1))/pdel(i,k)
               faltndqce = (faloutc(k) - faloutc(k-1))/pdel(i,k)
               faltndc = (faloutc(k) - dum*faloutc(k-1))/pdel(i,k)
               faltndnc = (faloutnc(k) - dum*faloutnc(k-1))/pdel(i,k)

! add fallout terms to eulerian tendencies
               qitend(i,k) = qitend(i,k) - faltndi/nstep
               nitend(i,k) = nitend(i,k) - faltndni/nstep
               qctend(i,k) = qctend(i,k)-  faltndc/nstep
               nctend(i,k) = nctend(i,k) - faltndnc/nstep

               IF (diag_id%qldt_sedi + diag_id%ql_sedi_col > 0) &
                       diag_4d(i,j,k,diag_pt%qldt_sedi) =    &
                         diag_4d(i,j,k,diag_pt%qldt_sedi) - faltndc/nstep
               IF (diag_id%qndt_sedi + diag_id%qn_sedi_col > 0) &
                       diag_4d(i,j,k,diag_pt%qndt_sedi) =    &
                         diag_4d(i,j,k,diag_pt%qndt_sedi) - faltndnc/nstep
               IF (diag_id%qnidt_sedi + diag_id%qni_sedi_col > 0) &
                       diag_4d(i,j,k,diag_pt%qnidt_sedi) =    &
                         diag_4d(i,j,k,diag_pt%qnidt_sedi) - faltndni/nstep

! add terms to to evap/sub of cloud water
               qvlat(i,k) = qvlat(i,k) - (faltndqie - faltndi)/nstep
               if (diag_id%qdt_sedi_ice2vapor > 0)  &
                      sedi_ice2vapor1(i,k)  = sedi_ice2vapor1(i,k) -  &
                                               (faltndqie - faltndi)/nstep
               qvlat(i,k) = qvlat(i,k) - (faltndqce - faltndc)/nstep
               if  (diag_id%qdt_sedi_liquid2vapor> 0)  &
                    sedi_liquid2vapor1(i,k)  = sedi_liquid2vapor1(i,k) -  &
                                              (faltndqce - faltndc)/nstep
               tlat(i,k) = tlat(i,k) + (faltndqie - faltndi)*xxls/nstep
               tlat(i,k) = tlat(i,k) + (faltndqce - faltndc)*xxlv/nstep

               dumi(i,k) = dumi(i,k) - faltndi*deltat/nstep
               dumni(i,k) = dumni(i,k) - faltndni*deltat/nstep
               dumc(i,k) = dumc(i,k) - faltndc*deltat/nstep
               dumnc(i,k) = dumnc(i,k) - faltndnc*deltat/nstep

               Fni(K) = MAX(Fni(K)/pdel(i,K), Fni(K-1)/pdel(i,K-1))*   &
                                                                 pdel(i,K)
               FI(K) = MAX(FI(K)/pdel(i,K), FI(K-1)/pdel(i,K-1))*pdel(i,K)
               fnc(k) = max(fnc(k)/pdel(i,k), fnc(k-1)/pdel(i,k-1))*   &
                                                                  pdel(i,k)
               Fc(K) = MAX(Fc(K)/pdel(i,K), Fc(K-1)/pdel(i,K-1))*pdel(i,K)

               IF (diag_id%qidt_fall + diag_id%qi_fall_col > 0) &
                      diag_4d(i,j,k,diag_pt%qidt_fall) =    &
                          diag_4d(i,j,k,diag_pt%qidt_fall) - faltndi/nstep
             end do   !! k loop

! units below are m/s
! cloud water/ice sedimentation flux at surface is added to precip flux 
! at surface to get total precip (cloud + precip water) rate
             prect(i) = prect(i) + (faloutc(kdim) + falouti(kdim)) &
                                                  /grav/nstep/1000._mg_pr
             preci(i) = preci(i) + (falouti(kdim))/grav/nstep/1000._mg_pr

             IF (diag_id%sedi_sfc > 0) &
                    diag_4d(i,j,1,diag_pt%sedi_sfc) =    &
                                  diag_4d(i,j,1,diag_pt%sedi_sfc)  +  & 
                                             faloutc(kdim)/grav/nstep/rhow
           end do   !! nstep loop

! end sedimentation

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! get new update for variables that includes sedimentation tendency
! note : here dum variables are grid-average, NOT in-cloud

           do k=1,kdim
             dumc(i,k) = max(qc(i,k) + qctend(i,k)*deltat, 0._mg_pr)
             dumi(i,k) = max(qi(i,k) + qitend(i,k)*deltat, 0._mg_pr)
             dumnc(i,k) = max(nc(i,k) + nctend(i,k)*deltat, 0._mg_pr)
             dumni(i,k) = max(ni(i,k) + nitend(i,k)*deltat, 0._mg_pr)

             if (dumc(i,k) .lt. qsmall) dumnc(i,k) = 0._mg_pr 
             if (dumi(i,k) .lt. qsmall) dumni(i,k) = 0._mg_pr

! calculate instantaneous processes (melting, homogeneous freezing)

             if (t(i,k) + tlat(i,k)/cpp*deltat > tfreeze) then
               if (dumi(i,k) > 0._mg_pr) then

! limit so that melting does not push temperature below freezing
                 dum = -dumi(i,k)*xlf/cpp
                 if (t(i,k) + tlat(i,k)/cpp*deltat + dum    &
                                                       .lt.  tfreeze) then
                   dum = (t(i,k) + tlat(i,k)/cpp*deltat-tfreeze)*cpp/xlf
                   dum = dum/dumi(i,k)*xlf/cpp 
                   dum = max(0._mg_pr, dum)
                   dum = min(1._mg_pr, dum)
                 else
                   dum = 1._mg_pr
                 end if

                 qctend(i,k) = qctend(i,k) + dum*dumi(i,k)/deltat

! hm add, 9/15/06, assume melting ice produces droplet
! mean volume radius of 8 micron

                 IF (diag_id%qndt_melt + diag_id%qn_melt_col > 0) &
                       diag_4d(i,j,k,diag_pt%qndt_melt) = nctend(i,k)
                 IF (diag_id%qidt_melt2 + diag_id%qi_melt2_col > 0) &
                        diag_4d(i,j,k,diag_pt%qidt_melt2) = qitend(i,k)
                 IF (diag_id%qnidt_melt + diag_id%qni_melt_col > 0) &
                        diag_4d(i,j,k,diag_pt%qnidt_melt) = nitend(i,k)
                 nctend(i,k) = nctend(i,k) + 3._mg_pr*dum*dumi(i,k)/  &
                                  deltat/(4._mg_pr*pi*5.12e-16_mg_pr*rhow)
                 qitend(i,k) = ((1._mg_pr - dum)*dumi(i,k) - qi(i,k))/  &
                                                                   deltat
                 nitend(i,k) = ((1._mg_pr - dum)*dumni(i,k) - ni(i,k))/   &
                                                                    deltat
                 tlat(i,k) = tlat(i,k) - xlf*dum*dumi(i,k)/deltat

                 IF (diag_id%qndt_melt + diag_id%qn_melt_col > 0) &
                         diag_4d(i,j,k,diag_pt%qndt_melt) =     &
                            nctend(i,k) - diag_4d(i,j,k,diag_pt%qndt_melt)
                 IF (diag_id%qidt_melt2  + diag_id%qi_melt2_col > 0) &
                         diag_4d(i,j,k,diag_pt%qidt_melt2) =    &
                            qitend(i,k) - diag_4d(i,j,k,diag_pt%qidt_melt2)
                 IF (diag_id%qnidt_melt + diag_id%qni_melt_col > 0) &
                        diag_4d(i,j,k,diag_pt%qnidt_melt) =     &
                           nitend(i,k) - diag_4d(i,j,k,diag_pt%qnidt_melt)
               end if
             end if

! homogeneously freeze droplets at -40 C
             if (t(i,k) + tlat(i,k)/cpp*deltat < tmin_fice) then
               if (dumc(i,k) .ge. 0._mg_pr) then
! limit so that freezing does not push temperature above threshold
                 dum = dumc(i,k)*xlf/cpp
                 if (t(i,k) + tlat(i,k)/cpp*deltat + dum .gt.   &
                                                           tmin_fice) then
                   dum = -(t(i,k) + tlat(i,k)/cpp*deltat - tmin_fice)*  &
                                                                   cpp/xlf
                   dum = dum/dumc(i,k)*xlf/cpp
                   dum = max(0._mg_pr, dum)
                   dum = min(1._mg_pr, dum)
                 else
                   dum = 1._mg_pr
                 end if
                 qitend(i,k) = qitend(i,k) + dum*dumc(i,k)/deltat
                 IF (diag_id%qldt_freez  + diag_id%ql_freez_col > 0) &
                        diag_4d(i,j,k,diag_pt%qldt_freez) = qctend(i,k)
                        sum_freeze(i,k) = qctend(i,k)
                 IF (diag_id%qndt_ihom + diag_id%qn_ihom_col > 0) &
                        diag_4d(i,j,k,diag_pt%qndt_ihom) =  nctend(i,k)
                 IF (diag_id%qnidt_ihom + diag_id%qni_ihom_col > 0) &
                        diag_4d(i,j,k,diag_pt%qnidt_ihom) =  nitend(i,k)

! hm add 11/18/06
! assume 25 micron mean volume radius of homogeneously frozen droplets
! consistent with size of detrained ice in stratiform.F90
!
! cms nitend .ne. -nctend
!
! 4/24/12: replace the 1.563 with x**3, here and in nCAR routine.

                 nitend(i,k) = nitend(i,k) + dum*3._mg_pr*dumc(i,k)/  &
                         (4._mg_pr*3.14_mg_pr*1.563e-14_mg_pr*rhoi)/deltat
                 qctend(i,k) = ((1._mg_pr - dum)*dumc(i,k) - qc(i,k))/   &
                                                                     deltat
                 nctend(i,k) = ((1._mg_pr - dum)*dumnc(i,k) - nc(i,k))/  &
                                                                     deltat
                 tlat(i,k) = tlat(i,k) + xlf*dum*dumc(i,k)/deltat

                 IF (diag_id%qldt_freez + diag_id%ql_freez_col > 0) &
                        diag_4d(i,j,k,diag_pt%qldt_freez)  =    &
                          qctend(i,k) - diag_4d(i,j,k,diag_pt%qldt_freez) 
                 sum_freeze(i,k) = -(qctend(i,k) - sum_freeze(i,k))
                 IF (diag_id%qndt_ihom + diag_id%qn_ihom_col > 0) &
                         diag_4d(i,j,k,diag_pt%qndt_ihom) =    &
                            nctend(i,k) - diag_4d(i,j,k,diag_pt%qndt_ihom)
                 IF (diag_id%qnidt_ihom + diag_id%qni_ihom_col > 0) &
                        diag_4d(i,j,k,diag_pt%qnidt_ihom) =   &
                            nitend(i,k) - diag_4d(i,j,k,diag_pt%qnidt_ihom)
               end if
             end if

             ssat_disposal(i,k) = 0._mg_pr

             IF (no_rh_adj_opt .EQ. 0) THEN
               IF (Nml%super_ice_opt .EQ. 0) THEN

!RSH: This is the code section which is intended to mimic the RK treatment.
!RSH  (super_ice_opt = 0)
!--++
!cms remove supersat 

                 rho(i,k) = pfull(i,k)/(rdgas*t_in(i,k))
                 ttmp = t(i,k) + tlat(i,k)/cpp*deltat
                 qtmp = q(i,k) + qvlat(i,k)*deltat
                 eslt = polysvp_l(ttmp)
                 esit = polysvp_i(ttmp)
                 esn = min(esit, eslt)
                 tmp2 = pfull(i,k) - d378*esn
                 tmp2 = max(tmp2, esn)
                 call lookup_des(ttmp, tmp7)
                 tmp7 = d622*pfull(i,k)*tmp7/tmp2/tmp2
                 tmp2 = d622*esn/tmp2

! the following is not consistent with the apportioning below ...
                 tmp3 = tmp7*(min(1., max(0._mg_pr,    &
                         0.05_mg_pr*(ttmp - tfreeze + 20._mg_pr)))*hlv + &
                          min(1., max(0._mg_pr,      &
                               0.05_mg_pr*(tfreeze - ttmp)))*hls)/cp_air

! compute excess over saturation
                 tmp1 = max(0._mg_pr, qtmp - tmp2)/(1._mg_pr + tmp3)

! change vapor content
                 qvlat(i,k) = qvlat(i,k) - tmp1/deltat 
                 if  (diag_id%qdt_super_sat_rm > 0)  &
                         super_saturation_rm1(i,k) =      &
                                  super_saturation_rm1(i,k) - tmp1/deltat
                 if  (diag_id%qdt_super_sat_rm > 0)  &
                         diag_4d(i,j,k,diag_pt%qdt_super_sat_rm) =    &
                                           super_saturation_rm1(i,k)
!CHANGE
! add in excess to cloud condensate, change cloud area and 
! increment temperature
                 if (ttmp .le. tfreeze - 40._mg_pr .and.     &
                                          tmp1 .gt. 0._mg_pr) then
                   IF (Nml%super_choice) THEN 
                     qitend(i,k) = qitend(i,k) + tmp1/deltat
                     ssat_disposal(i,k) = 2._mg_pr
                     IF (diag_id%ice_adj + diag_id%ice_adj_col > 0) &
                            diag_4d(i,j,k,diag_pt%ice_adj) = tmp1/deltat
                     sum_ice_adj(i,k) = tmp1/deltat
                   ELSE
                     ssat_disposal(i,k) = 0._mg_pr
                     prect(i) = prect(i) + tmp1/deltat*pdel(i,k)/grav/rhow
                     preci(i) = preci(i) + tmp1/deltat*pdel(i,k)/grav/rhow
                     qsout(i,k) = qsout(i,k) + tmp1
                     asnowrt(i,k) = asnowrt(i,k) + tmp1/deltat*    &
                                                           pdel(i,k)/grav 
                     atotrt(i,k) = atotrt(i,k) + tmp1/deltat *   &
                                                          pdel(i,k)/grav 
                   END IF
                   tlat(i,k) = tlat(i,k) + hls*tmp1/deltat 
                 end if

                 if (ttmp .gt. tfreeze - 40._mg_pr .and.    &
                                               tmp1 .gt. 0._mg_pr) then   
                   IF (Nml%super_choice) THEN 
                     qctend(i,k) = qctend(i,k) + tmp1/deltat 
                     ssat_disposal(i,k) = 1._mg_pr
                     IF (diag_id%liq_adj + diag_id%liq_adj_col > 0 ) &
                           diag_4d(i,j,k,diag_pt%liq_adj) = tmp1/deltat
                   ELSE
                     ssat_disposal(i,k) = 0._mg_pr
                     prect(i) = prect(i) + tmp1/deltat*pdel(i,k)/   &
                                                               grav/rhow
                     qrout(i,k) = qrout(i,k) + tmp1
                     atotrt(i,k) = atotrt(i,k) +  tmp1/deltat*   &
                                                            pdel(i,k)/grav 
                   END IF
                   tlat(i,k) = tlat(i,k) + hlv*tmp1/deltat             
                 end if
               ELSE IF (Nml%super_ice_opt .GE. 1) THEN 
                 ssat_disposal(i,k) = 0._mg_pr

!RSH: THis is the code section  used with M-G (super_ice_opt >= 1).
                 IF (sat_adj_opt .EQ. 1) THEN
                   ttmp = t(i,k) + tlat(i,k)/cpp*deltat
                   qtmp = q(i,k) + qvlat(i,k)*deltat
                   qs_t = polysvp_l(ttmp)
! calculate denominator in qsat formula
                   qs_d = pfull(i,k) - d378*qs_t
! limit denominator to esat, and thus qs to d622
! this is done to avoid blow up in the upper stratosphere
! where pfull ~ esat  
                   qs_d = max(qs_d, qs_t)
! calculate qs
                   qs_t = d622*qs_t/qs_d

! compute super saturation
                   tmp1 = max(0._mg_pr, (qtmp - qs_t))/(1._mg_pr +    &
                                       hlv*qs_t/(rvgas*ttmp**2)*hlv/cp_air)
 
! change vapor content
                   qvlat(i,k) = qvlat(i,k) - tmp1/deltat

! add in excess to cloud condensate, change cloud area and 
! increment temperature
                   if (ttmp .le. tfreeze - 40._mg_pr .and.   &
                                                 tmp1 .gt. 0._mg_pr) then
                     ssat_disposal(i,k) = 2._mg_pr
                     qitend(i,k) = qitend(i,k) + tmp1/deltat
                     IF (diag_id%ice_adj + diag_id%ice_adj_col > 0) &
                               diag_4d(i,j,k,diag_pt%ice_adj) = tmp1/deltat
                     sum_ice_adj(i,k) = tmp1/deltat
                     if  (diag_id%qdt_super_sat_rm > 0)  &
                              super_saturation_rm1(i,k)  =     &
                                   super_saturation_rm1(i,k) - tmp1/deltat
                     if  (diag_id%qdt_super_sat_rm > 0)  &
                              diag_4d(i,j,k,diag_pt%qdt_super_sat_rm) =  &
                                             super_saturation_rm1(i,k)
                     tlat(i,k) = tlat(i,k) + hls*tmp1/deltat
                   end if 

                   if (ttmp .gt. tfreeze - 40._mg_pr .and.    &
                                                  tmp1 .gt. 0._mg_pr) then
                     qctend(i,k) = qctend(i,k) + tmp1/deltat
                     ssat_disposal(i,k) = 1._mg_pr
                     IF (diag_id%liq_adj + diag_id%liq_adj_col > 0 ) &
                            diag_4d(i,j,k,diag_pt%liq_adj) = tmp1/deltat
                     if  (diag_id%qdt_super_sat_rm > 0)  &
                                super_saturation_rm1(i,k)  =   &
                                    super_saturation_rm1(i,k) - tmp1/deltat
                     if  (diag_id%qdt_super_sat_rm > 0)  &
                              diag_4d(i,j,k,diag_pt%qdt_super_sat_rm) =  &
                                          super_saturation_rm1(i,k)
                     tlat(i,k) = tlat(i,k) + hlv*tmp1/deltat
                   end if   
                 END IF 
               END IF ! (super_ice_opt .EQ. 0)
             END IF ! (no_rh_adj_opt .EQ. 0)
!........................................................................
! do not calculate effective radius here,
! but do some limiting
! update cloud variables after instantaneous processes to get effective 
! radius
! variables are in-cloud to calculate size dist parameters

             dumc(i,k) = max(qc(i,k) + qctend(i,k)*deltat, 0._mg_pr)/   &
                                                                 cldm(i,k)
             dumi(i,k) = max(qi(i,k) + qitend(i,k)*deltat, 0._mg_pr)/  &
                                                                 cldm(i,k)
             dumnc(i,k) = max(nc(i,k) + nctend(i,k)*deltat, 0._mg_pr)/  &
                                                                 cldm(i,k)
             dumni(i,k) = max(ni(i,k) + nitend(i,k)*deltat, 0._mg_pr)/  &
                                                                cldm(i,k)

! limit in-cloud mixing ratio to reasonable value of 5 g kg-1
             dumc(i,k) = min(dumc(i,k), in_cloud_limit)
             dumi(i,k) = min(dumi(i,k), in_cloud_limit)

!...................
! cloud ice effective radius

             if (dumi(i,k) .ge. qsmall) then
               IF (diag_id%qnidt_size_adj +   &
                              diag_id%qni_size_adj_col > 0) THEN
                      diag_4d(i,j,k,diag_pt%qnidt_size_adj ) =  nitend(i,k)
               END IF

! add upper limit to in-cloud number concentration to prevent 
! numerical error
               dumni(i,k) = min(dumni(i,k), dumi(i,k)*1.e20_mg_pr)
               lami(k) = (gamma_mg(1._mg_pr + di_mg)*ci_mg* &
                                dumni(i,k)/dumi(i,k))**(1._mg_pr/di_mg)
               lammax = 1._mg_pr/min_diam_ice
               lammin = 1._mg_pr/(2._mg_pr*dcs)
               if (lami(k) .lt. lammin) then
                 lami(k) = lammin
                 n0i(k) = lami(k)**(di_mg + 1._mg_pr)*dumi(i,k)/   &
                                        (ci_mg*gamma_mg(1._mg_pr + di_mg))
                 niic(i,k) = n0i(k)/lami(k)

! adjust number conc if needed to keep mean size in reasonable range
                 nitend(i,k) = (niic(i,k)*cldm(i,k) - ni(i,k))/deltat
               else if (lami(k) .gt. lammax) then
                 lami(k) = lammax
                 n0i(k) = lami(k)**(di_mg + 1._mg_pr)*dumi(i,k)/   &
                                        (ci_mg*gamma_mg(1._mg_pr + di_mg))
                 niic(i,k) = n0i(k)/lami(k)

! adjust number conc if needed to keep mean size in reasonable range
                 nitend(i,k) = (niic(i,k)*cldm(i,k) - ni(i,k))/deltat
               end if
               IF (diag_id%qnidt_size_adj + diag_id%qni_size_adj_col > 0) &
                    diag_4d(i,j,k,diag_pt%qnidt_size_adj) = nitend(i,k) - &
                                     diag_4d(i,j,k,diag_pt%qnidt_size_adj) 
             end if

!...................
! cloud droplet effective radius

             if (dumc(i,k) .ge. qsmall) then
               IF (diag_id%qndt_size_adj + diag_id%qn_size_adj_col > 0) &
                      diag_4d(i,j,k,diag_pt%qndt_size_adj ) =  nctend(i,k)

! add upper limit to in-cloud number concentration to prevent 
! numerical error
               dumnc(i,k) = min(dumnc(i,k), dumc(i,k)*1.e20_mg_pr)
!RSH BUGFIX email of 6/8/10
!              pgam(k)=0.0005714_mg_pr*(dumnc(i,k)/1.e6_mg_pr/rho(i,k))+ &
!                                                              0.2714_mg_pr
               pgam(k) = 0.0005714_mg_pr*(dumnc(i,k)/1.e6_mg_pr*    &
                                                  rho(i,k)) + 0.2714_mg_pr
               pgam(k)=1._mg_pr/(pgam(k)**2) - 1._mg_pr
               pgam(k) = max(pgam(k), 2._mg_pr)
               pgam(k) = min(pgam(k), 15._mg_pr)
               lamc(k) = (pi/6._mg_pr*rhow*dumnc(i,k)*  &
                                    gamma_mg(pgam(k) + 4._mg_pr)/ &
                                              (dumc(i,k)*    &
                        gamma_mg(pgam(k) + 1._mg_pr)))**(1._mg_pr/3._mg_pr)
!              lammin = (pgam(k)+1._mg_pr)/50.e-6_mg_pr
!              lammax = (pgam(k)+1._mg_pr)/2.e-6_mg_pr
               lammin = (pgam(k) + 1._mg_pr)/max_diam_drop
               lammax = (pgam(k) + 1._mg_pr)/min_diam_drop
               if (lamc(k) .lt. lammin) then
                 lamc(k) = lammin
                 ncic(i,k) = 6._mg_pr*lamc(k)**3*dumc(i,k)* &
                                  gamma_mg(pgam(k) + 1._mg_pr)/ &
                                     (pi*rhow*gamma_mg(pgam(k) + 4._mg_pr))

! adjust number conc if needed to keep mean size in reasonable range
                 nctend(i,k) = (ncic(i,k)*cldm(i,k) - nc(i,k))/deltat
               else if (lamc(k).gt.lammax) then
                 lamc(k) = lammax
                 ncic(i,k) = 6._mg_pr*lamc(k)**3*dumc(i,k)* &
                                  gamma_mg(pgam(k) + 1._mg_pr)/ &

                                     (pi*rhow*gamma_mg(pgam(k) + 4._mg_pr))
! adjust number conc if needed to keep mean size in reasonable range
                 nctend(i,k) = (ncic(i,k)*cldm(i,k) - nc(i,k))/deltat
               end if
               IF (diag_id%qndt_size_adj + diag_id%qn_size_adj_col > 0) &
                      diag_4d(i,j,k,diag_pt%qndt_size_adj ) =    &
                       nctend(i,k) - diag_4d(i,j,k,diag_pt%qndt_size_adj) 
             end if

!...................
! rain drop effective size
             if (qrout(i,k) .gt. 1.e-7_mg_pr .and.    &
                                 nrout(i,k) .gt. 0._mg_pr) then
               lsc_rain_size(i,k) = 3.0_mg_pr*(pi*rhow*nrout(i,k)/    &
                              qrout(i,k))**(-1._mg_pr/3._mg_pr)*1.e6_mg_pr
             else
               lsc_rain_size(i,k) = 100._mg_pr
             endif
           end do   ! (k loop) 
500        CONTINUE


           do k=1,kdim
! if updated q (after microphysics) is zero, then ensure updated n is also zero
             IF (diag_id%qndt_fill2 + diag_id%qn_fill2_col > 0) &
                    diag_4d(i,j,k,diag_pt%qndt_fill2 ) = nctend(i,k)
             IF (diag_id%qnidt_fill2 +  diag_id%qni_fill2_col > 0) &
                    diag_4d(i,j,k,diag_pt%qnidt_fill2 ) = nitend(i,k)
             if (qc(i,k) + qctend(i,k)*deltat .lt. qsmall)    &
                                             nctend(i,k) = -nc(i,k)/deltat
             if (qi(i,k) + qitend(i,k)*deltat .lt. qsmall)     &
                                             nitend(i,k) = -ni(i,k)/deltat
             IF (diag_id%qndt_fill2 + diag_id%qn_fill2_col > 0) &
                      diag_4d(i,j,k,diag_pt%qndt_fill2 ) =    &
                          nctend(i,k) - diag_4d(i,j,k,diag_pt%qndt_fill2) 
             IF (diag_id%qnidt_fill2 + diag_id%qni_fill2_col > 0) &
                     diag_4d(i,j,k,diag_pt%qnidt_fill2 ) =   &
                         nitend(i,k) - diag_4d(i,j,k,diag_pt%qnidt_fill2) 
           end do
         end do ! (do i=1,idim)

         rain3d(:,:,1) = 0._mg_pr
         snow3d(:,:,1) = 0._mg_pr
         DO k=1,kdim
           DO i=1,idim

! STILL NEED TO DEAL WITH SSAT GOING to PRECIP RATHER THAN CLOUD -
             rain3d(i,j,k+1) = MAX((atotrt(i,k) - asnowrt(i,k)), 0._mg_pr)
             snow3d(i,j,k+1) = MAX(asnowrt(i,k), 0._mg_pr)
           END DO
         END DO
 
         if (diag_id%rain_evap + diag_id%rain_evap_col > 0)  &
                diag_4d(:,j,:,diag_pt%rain_evap) = pre1(:,:)/real(iter)

         if (diag_id%qdt_rain_evap > 0)  &
                diag_4d(:,j,:,diag_pt%qdt_rain_evap) = pre1(:,:)/real(iter)

         if (diag_id%qdt_cond > 0) &
                diag_4d(:,j,:,diag_pt%qdt_cond) = -cmel1(:,:)/real(iter)

         if (diag_id%qdt_snow_sublim > 0 .or.  &
                         diag_id%q_snow_sublim_col > 0)  &
                diag_4d(:,j,:,diag_pt%qdt_snow_sublim )  =  &
                                                     prds1(:,:)/real(iter)

         if (diag_id%qdt_deposition > 0)  &
                diag_4d(:,j,:,diag_pt%qdt_deposition)  =   &
                                       -cmei1(:,:)/real(iter)

         if  (diag_id%qdt_eros_l > 0)  &
              diag_4d(:,j,:,diag_pt%qdt_eros_l)  =    &
                                    -D_eros_l1(:,:)/real(iter)

         if  (diag_id%qdt_eros_i > 0)  &
                 diag_4d(:,j,:,diag_pt%qdt_eros_i)  =    &
                                     -D_eros_i1(:,:)/real(iter)

         if  (diag_id%qdt_qv_on_qi > 0)  &
                   diag_4d(:,j,:,diag_pt%qdt_qv_on_qi) =  &
                                    - qvdep_qi1(:,:)/real(iter)

         if  (diag_id%qdt_snow2vapor + diag_id%q_snow2vapor_col > 0)  &
                    diag_4d(:,j,:,diag_pt%qdt_snow2vapor) =   &
                                        snow2vapor1(:,:)/real(iter)

         if  (diag_id%qdt_sedi_ice2vapor > 0)  &
               diag_4d(:,j,:,diag_pt%qdt_sedi_ice2vapor) =   &
                                            sedi_ice2vapor1(:,:)

         if  (diag_id%qdt_sedi_liquid2vapor > 0)  &
                  diag_4d(:,j,:,diag_pt%qdt_sedi_liquid2vapor) =  &
                                         sedi_liquid2vapor1(:,:)

         if  (diag_id%qdt_super_sat_rm > 0)  &
              diag_4d(:,j,:,diag_pt%qdt_super_sat_rm) =    &
                                      super_saturation_rm1(:,:)

!       diagnostics for cloud liquid tendencies
!cloud water  
         if (diag_id%qldt_cond  + diag_id%ql_cond_col > 0) &
              diag_4d(:,j,:,diag_pt%qldt_cond)  =     &
                                    max(cmel1(:,:), 0._mg_pr)/real(iter)

         if (diag_id%qldt_evap  + diag_id%ql_evap_col > 0) &
              diag_4d(:,j,:,diag_pt%qldt_evap)  =     &
                          - max(-1._mg_pr*cmel1(:,:),0._mg_pr)/real(iter) 

         if (diag_id%qldt_eros + diag_id%ql_eros_col > 0) &
              diag_4d(:,j,:,diag_pt%qldt_eros) = D_eros_l1(:,:)/real(iter) 

         if (diag_id%qldt_berg + diag_id%ql_berg_col > 0) &
              diag_4d(:,j,:,diag_pt%qldt_berg) = - berg1(:,:)/real(iter)   

         sum_berg(:,:) =  berg1(:,:)/real(iter)
         if (diag_id%qldt_auto + diag_id%ql_auto_col > 0) &
              diag_4d(:,j,:,diag_pt%qldt_auto) = prc1(:,:)/real(iter)   

         if (diag_id%qldt_freez2 + diag_id%ql_freez2_col > 0) &
              diag_4d(:,j,:,diag_pt%qldt_freez2) = mnuccc1(:,:)/real(iter) 
         sum_freeze2(:,:) = -mnuccc1(:,:)/real(iter)

         if (diag_id%qldt_accr  + diag_id%ql_accr_col > 0) &
              diag_4d(:,j,:,diag_pt% qldt_accr) = pra1(:,:) /real(iter)   

         if (diag_id%qldt_accrs  + diag_id%ql_accrs_col > 0) & 
              diag_4d(:,j,:,diag_pt%qldt_accrs) = psacws1(:,:)/real(iter)
         sum_rime(:,:) = -psacws1(:,:)/real(iter)

         if (diag_id%qldt_bergs + diag_id%ql_bergs_col > 0) &
              diag_4d(:,j,:,diag_pt%qldt_bergs) = bergs1(:,:)/real(iter)
         sum_bergs(:,:) = -bergs1(:,:)/real(iter)

!cloud ice 
         if (diag_id%qidt_dep + diag_id%qi_dep_col > 0)    &
             diag_4d(:,j,:,diag_pt%qidt_dep)  =   &
                                 max(cmei1(:,:), 0._mg_pr) /real(iter)
         sum_cond (:,:) = max(cmei1(:,:), 0._mg_pr) /real(iter)

         if (diag_id%qidt_subl + diag_id%qi_subl_col > 0)  &
                 diag_4d(:,j,:,diag_pt%qidt_subl) =     &
                          -max(-1._mg_pr*cmei1(:,:), 0._mg_pr) /real(iter) 

         if (diag_id%qidt_eros + diag_id%qi_eros_col > 0) &
             diag_4d(:,j,:,diag_pt%qidt_eros) = D_eros_i1(:,:)/real(iter) 

         if (diag_id%qidt_auto + diag_id%qi_auto_col > 0) & 
             diag_4d(:,j,:,diag_pt%qidt_auto) = prci1(:,:)/real(iter) 

         if (diag_id%qidt_accr  + diag_id%qi_accr_col > 0) &
             diag_4d(:,j,:,diag_pt%qidt_accr) = prai1(:,:)/real(iter)

         if (diag_id%qidt_accrs  + diag_id%qi_accrs_col > 0) &
             diag_4d(:,j,:,diag_pt%qidt_accrs) = psacws_o1(:,:)/real(iter)

         if (diag_id%qidt_qvdep + diag_id%qi_qvdep_col > 0)    &
             diag_4d(:,j,:,diag_pt%qidt_qvdep) = qvdep_qi1(:,:)/real(iter)

!cloud droplet number
         if ( diag_id%qndt_nucclim + diag_id%qn_nucclim_col > 0) &
            diag_4d(:,j,:,diag_pt%qndt_nucclim) = nucclim1(:,:)/real(iter)

         if (diag_id%qndt_cond + diag_id%qn_cond_col > 0)    &
             diag_4d(:,j,:,diag_pt%qndt_cond) = npccn1(:,:)/real(iter)

         if (diag_id%qndt_freez + diag_id%qn_freez_col > 0)    &
             diag_4d(:,j,:,diag_pt%qndt_freez) = nnuccc1(:,:)/real(iter)

         if (diag_id%qndt_sacws + diag_id%qn_sacws_col > 0)    &
             diag_4d(:,j,:,diag_pt%qndt_sacws) = npsacws1(:,:)/real(iter)

         if (diag_id%qndt_sacws_o + diag_id%qn_sacws_o_col > 0)    &
             diag_4d(:,j,:,diag_pt%qndt_sacws_o) =     &
                                                 npsacws_o1(:,:)/real(iter)

         if (diag_id%qndt_evap + diag_id%qn_evap_col > 0)    &
             diag_4d(:,j,:,diag_pt%qndt_evap) = nsubc1(:,:)/real(iter)

         if (diag_id%qndt_eros + diag_id%qn_eros_col > 0)    &
             diag_4d(:,j,:,diag_pt%qndt_eros) = nerosc1(:,:)/real(iter)

         if (diag_id%qndt_pra + diag_id%qn_pra_col > 0)    &
             diag_4d(:,j,:,diag_pt%qndt_pra) = npra1(:,:)/real(iter)

         if (diag_id%qndt_auto + diag_id%qn_auto_col > 0)    &
             diag_4d(:,j,:,diag_pt%qndt_auto) = nprc11(:,:)/real(iter)

!cloud droplet number
         if (diag_id%qnidt_nnuccd + diag_id%qni_nnuccd_col > 0)    &
             diag_4d(:,j,:,diag_pt%qnidt_nnuccd) = nnuccd1(:,:)/real(iter)

         if (diag_id%qnidt_nsubi  + diag_id%qni_nsubi_col  > 0)    &
             diag_4d(:,j,:,diag_pt%qnidt_nsubi) = nsubi1(:,:)/real(iter)

         if (diag_id%qnidt_nerosi  + diag_id%qni_nerosi_col  > 0)    &
             diag_4d(:,j,:,diag_pt%qnidt_nerosi) = nerosi1(:,:)/real(iter)

         if (diag_id%qnidt_nprci  + diag_id%qni_nprci_col  > 0)    &
             diag_4d(:,j,:,diag_pt%qnidt_nprci) = nprci1(:,:)/real(iter)

         if (diag_id%qnidt_nprai  + diag_id%qni_nprai_col  > 0)    &
             diag_4d(:,j,:,diag_pt%qnidt_nprai) = nprai1(:,:)/real(iter)

         if (diag_id%qnidt_nucclim1  + diag_id%qni_nucclim1_col  > 0)    &
             diag_4d(:,j,:,diag_pt%qnidt_nucclim1) =    &
                                                 nucclim1_1(:,:)/real(iter)

         if (diag_id%qnidt_nucclim2  + diag_id%qni_nucclim2_col  > 0)    &
             diag_4d(:,j,:,diag_pt%qnidt_nucclim2) =     &
                                                 nucclim2_1(:,:)/real(iter)

         if (diag_id%srfrain_accrs + diag_id%srfrain_accrs_col  > 0)    &
             diag_4d(:,j,:,diag_pt%srfrain_accrs) =   &
                                                 pracs1(:,:)/real(iter)
 
         if (diag_id%srfrain_freez + diag_id%srfrain_freez_col > 0)    &
             diag_4d(:,j,:,diag_pt%srfrain_freez) =    &
                                                 mnuccr1(:,:)/real(iter)
!RSH:
!   calculate fraction of total ice/snow creation that requires
!   ice-forming nuclei
         do k=1,kdim
           do i=1,idim
             qldt_sum = sum_berg(i,k) + sum_rime(i,k) + sum_freeze(i,k) +&
                        MAX(sum_bergs(i,k), 0.0_mg_pr) + sum_cond(i,k) + &
                        sum_ice_adj(i,k) + sum_freeze2(i,k)
             if (qldt_sum > 0.0_mg_pr)  then
               f_snow_berg(i,k) = (sum_berg(i,k) + sum_cond(i,k) +   &
                                   sum_ice_adj(i,k) + sum_freeze(i,k) + &
                                   MAX(sum_bergs(i,k), 0.0))/qldt_sum
             else
               f_snow_berg(i,k) = 0._mg_pr
             endif
           end do
         end do

!-----------------------------------------------------------------------

END SUBROUTINE morrison_gettelman_microp


!#########################################################################

SUBROUTINE morrison_gettelman_microp_end (do_pdf_clouds)

LOGICAL, INTENT(IN ) :: do_pdf_clouds

!------------------------------------------------------------------------
!    call destructor routines for modules used here.
!------------------------------------------------------------------------
      call polysvp_end
      if (do_pdf_clouds) then
        CALL  simple_pdf_end
      endif 
      call gamma_mg_end

!------------------------------------------------------------------------
!    mark the module as uninitialized.
!------------------------------------------------------------------------
      module_is_initialized = .false.

!-----------------------------------------------------------------------


END SUBROUTINE morrison_gettelman_microp_end


!#########################################################################



END MODULE morrison_gettelman_microp_mod

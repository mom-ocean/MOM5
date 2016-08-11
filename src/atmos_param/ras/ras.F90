  MODULE RAS_MOD

!=======================================================================
!  RELAXED ARAKAWA/SCHUBERT (RAS) CUMULUS PARAM SCHEME MODULE          !
!=======================================================================
!  SUBROUTINE RAS_INIT        - INITIALIZE RAS
!  SUBROUTINE RAS             - DRIVER
!  SUBROUTINE RAS_CLOUD       - RAS CU PARAMETERIZATION
!  SUBROUTINE COMP_LCL        - COMPUTE LCL (CLOUD BASE)
!  SUBROUTINE RAS_CEVAP       - EVAPORATION OF CONVECTIVE SCALE PRECIP
!  SUBROUTINE RAS_CLOUD_INDEX - SET ORDER IN WHICH CLOUDS ARE TO BE DONE 
!  SUBROUTINE RAS_CLOUD_EXIST - TEST FOR INSTABILITY IN COLUMN
!  SUBROUTINE RAS_BDGT        - BUDGET CHECK FOR RAS
!  FUNCTION   RAN0            - RANDOM NUMBER GENERATOR
!=======================================================================

 use            mpp_mod, only: mpp_pe,             &
                               mpp_root_pe,        &
                               stdlog
 use Sat_Vapor_Pres_Mod, ONLY: compute_qs
 use      Constants_Mod, ONLY:  HLv, HLs, Cp_Air, Grav, Kappa, rdgas, rvgas
 use   Diag_Manager_Mod, ONLY: register_diag_field, send_data
 use   Time_Manager_Mod, ONLY: time_type
 use            mpp_mod, only: input_nml_file
 use            fms_mod, only: write_version_number, open_namelist_file, &
                               FILE_EXIST, ERROR_MESG, check_nml_error, &
                               CLOSE_FILE, FATAL
 use  field_manager_mod, only: MODEL_ATMOS
 use tracer_manager_mod, only: get_tracer_index,   &
                               get_number_tracers, &
                               get_tracer_names,   &
                               get_tracer_indices, &
                               query_method,       &
                               NO_TRACER
 use  rad_utilities_mod,  only : aerosol_type
 use  aer_ccn_act_mod,    only : aer_ccn_act, aer_ccn_act2
!---------------------------------------------------------------------
 implicit none
 private
 public  :: ras, ras_init, ras_end, ras_bdgt
!---------------------------------------------------------------------

!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 character(len=128) :: version = '$Id: ras.F90,v 19.0 2012/01/06 20:11:58 fms Exp $'
 character(len=128) :: tagname = '$Name: tikal $'
!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 real :: cp_div_grav
 real :: one_plus_kappa, one_minus_kappa
 real :: onebcp, onebg
 real :: rn_pfac

 logical :: module_is_initialized = .false.

 real, parameter :: d622        = rdgas/rvgas
 real, parameter :: d378        = 1.0-d622
 
!---------------------------------------------------------------------
! --- Climatological cloud work function data
!---------------------------------------------------------------------
 
 real                :: actop 
 real, dimension(15) :: ac, ad
 real, dimension(15) :: ph, a

 data ph / 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0,   &
           550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0 /

 data a / 1.6851, 1.1686, 0.7663, 0.5255, 0.4100, 0.3677,            &
          0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664,            &
          0.0553, 0.0445, 0.0633  /

!---------------------------------------------------------------------
! --- NAMELIST
!---------------------------------------------------------------------
!       fracs  - Fraction of PBL mass allowed to be used 
!                by a cloud-type in time DT
!       rasal0 - Base value for cloud type relaxation parameter
!       puplim - Upper limit for cloud tops of deepest clouds
!       aratio - Ratio of critical cloud work function to standard
!                value of cloud work function
!       cufric - Should Cumulus friction (momentum transport) occur?
!      rh_trig - Convection takes place only if the relative humidity
!                of the lowest model level exceeds rh_trig
!      alm_min - Min value for entrainment parameter.
!   Tokioka_on - If true, alm_min computed using Tokioka formulation
!  Tokioka_con - Constant for alm_min computed using Tokioka formulation
! Tokioka_plim - Tokioka applied only to clouds detraining above Tokioka_plim
!   modify_pbl - If true, mass flux in sub cloud layer varies linearly   
!                between value at cloud base and zero at surface, and   
!                tendencies are spread out throughout the entire sub cloud layer.
! prevent_unreasonable - If true, unreasonable states (negatives) for the water 
!                        tracers are prevented. This is achieved by borrowing 
!                        from the sources of tracer within this module. 
!                        Otherwise the tendency is added without correction and
!                        can lead to negative mixing ratios. (Default = .FALSE.)
! --- CLOUD ORDER SPECIFICATION ---
!     ncrnd  - Number of random cloud-types between krmin and krmax to be
!              invoked in a single call to ras
!     iseed  - Integer seed used in generating random numbers
!              -- used only when ncrnd > 0
!     krmin  - Index of the top most level to which random clouds may
!              be invoked
!     krmax  - Index of the bottom most level to which random clouds 
!              may be invoked. krmin should be <= krmax. 
!              If ncrnd is specified as zero, then all cloud-types 
!              below the level krmax will be called sequentially.
!     botop  - A logical variable -- .true. if sequential clouds are 
!              called from bottom to top and .false. if top to bottom.
!              Level krmax will be called sequentially.
! --- PARTION CLOUD LIQUID WATER INTO PRECIP & DETRAINED VAPOR AND LIQUID --- 
!     rn_ptop - rn_frac_top of parcel liquid water converted to precip 
!               for cloud top pressures above rn_ptop
!     rn_pbot - rn_frac_bot of parcel liquid water converted to precip 
!               for cloud top pressures below rn_pbot (linear profile in 
!               between)
!     rn_frac_bot - Fraction of parcel liquid water converted to 
!                   precip for cloud top pressures below rn_pbot
!     rn_frac_top - Fraction of liquid water converted to
!                   precip for cloud top pressures above rn_ptop
! --- EVAPORATION OF CONVECTIVE SCALE PRECIP ---
!     evap_on - Turn on evaporation if true
!     cfrac   - Fraction of grid box assumed to be covered by convection
!     hcevap  - Evap allowed while q <= hcevap * qsat
! --- CRITICAL CLOUD WORK FUNCTION ---
!     ph - array, dim(15), of cloud top pressures
!      a - array, dim(15), of critical cloud work functions 
!          for clouds detraining at level ph
!---------------------------------------------------------------------

 logical :: use_online_aerosol = .false.
 real    :: fracs        = 0.25
 real    :: rasal0       = 0.25
 real    :: puplim       = 20.0E2
 real    :: aratio       = 1.4
 logical :: cufric       = .false.
 real    :: rh_trig      = 0.0      
 real    :: alm_min      = 0.0  
 logical :: Tokioka_on   = .false.
 real    :: Tokioka_con  = 0.05 
 real    :: Tokioka_plim = 500.0E2
 logical :: modify_pbl   = .false.
 logical :: prevent_unreasonable   = .false.
 
! --- cloud order specification ---
 integer :: ncrnd = 0
 integer :: iseed = 123
 integer :: krmin = 2
 integer :: krmax = 2
 logical :: botop = .true.

! --- partion cloud liquid water into precip & detrained water --- 
 real :: rn_ptop = 500.0E2
 real :: rn_pbot = 800.0E2
 real :: rn_frac_bot = 0.8
 real :: rn_frac_top = 1.0

! --- evaporation of convective scale precip ---
 logical :: evap_on = .true.   
 real    :: cfrac   = 0.05
 real    :: hcevap  = 0.80    

 real    :: sea_salt_scale = 0.1
 real    :: om_to_oc = 1.67

 integer :: nqn   ! tracer indices for stratiform clouds

    NAMELIST / ras_nml /                          &
      fracs,   rasal0,  puplim, aratio, cufric,   &
      ncrnd,   iseed,   krmin,   krmax, botop,    &
      rn_ptop, rn_pbot, rn_frac_top, rn_frac_bot, &
      evap_on, cfrac,   hcevap, rh_trig, alm_min, &
      Tokioka_on, Tokioka_con, Tokioka_plim, modify_pbl, prevent_unreasonable, &
      ph, a, sea_salt_scale, om_to_oc, use_online_aerosol

!---------------------------------------------------------------------
! DIAGNOSTICS FIELDS 
!---------------------------------------------------------------------

integer :: id_tdt_revap,  id_qdt_revap,    id_prec_revap,  &
           id_snow_revap, id_prec_conv_3d, id_pcldb, &
           id_det0, &
           id_tdt_conv, id_qdt_conv, id_prec_conv, id_snow_conv, &
           id_q_conv_col, id_t_conv_col, id_mc,&
           id_qldt_conv, id_qidt_conv, id_qadt_conv, id_qndt_conv, &
           id_ql_conv_col, id_qi_conv_col, id_qa_conv_col,         &
           id_mfcb, id_mfct, id_alm, id_cfq_ras
integer, allocatable, dimension(:) :: id_tracer_conv, id_tracer_conv_col

character(len=3) :: mod_name = 'ras'

real :: missing_value = -999.

integer  :: num_ras_tracers = 0
logical  :: do_ras_tracer = .false.

!---------------------------------------------------------------------

 contains

!#####################################################################
!#####################################################################

  SUBROUTINE RAS_INIT( do_strat, do_liq_num, axes, Time, tracers_in_ras )

!=======================================================================
! ***** INITIALIZE RAS
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!---------------------------------------------------------------------

 logical,         intent(in) :: do_strat, do_liq_num
 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: Time
 logical, dimension(:), intent(in), optional :: tracers_in_ras

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

 integer             :: unit, io, ierr
 real                :: actp, facm 
 real, dimension(15) :: au,   tem
 integer, dimension(3) :: half = (/1,2,4/)
 character(len=128)    :: diagname, diaglname, tendunits, name, units
 integer               :: tr
 integer               :: num_tracers
 integer               :: nn, logunit
!=====================================================================

!---------------------------------------------------------------------
! --- Read namelist
!---------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ras_nml, iostat=io)
  ierr = check_nml_error(io,"ras_nml")
#else
  if( FILE_EXIST( 'input.nml' ) ) then
      unit = OPEN_NAMELIST_FILE ()
      ierr = 1 ; do while ( ierr .ne. 0 )
        READ( unit, nml = ras_nml, iostat = io, end = 10 )
        ierr = check_nml_error (io, 'ras_nml')
      end do
10  CALL CLOSE_FILE ( unit )
  end if
#endif

!---------------------------------------------------------------------
! --- Write namelist
!---------------------------------------------------------------------

  call write_version_number(version, tagname)
       logunit = stdlog()
       WRITE( logunit, nml = ras_nml ) 


!---------------------------------------------------------------------
! --- Find the tracer indices 
!---------------------------------------------------------------------
  call get_number_tracers(MODEL_ATMOS,num_prog=num_tracers)

  nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
  if (do_liq_num .and. nqn == NO_TRACER) &
    call error_mesg ('ras_init', &
         'prognostic droplet number scheme requested but tracer not found', FATAL) 


!---------------------------------------------------------------------
! --- Initialize constants
!---------------------------------------------------------------------

  
  cp_div_grav = Cp_Air / Grav

  one_plus_kappa  = 1.0 + Kappa
  one_minus_kappa = 1.0 - Kappa

  onebcp = 1.0 / Cp_Air
  onebg  = 1.0 / Grav
  
  rn_pfac = ( rn_frac_top -  rn_frac_bot)  / ( rn_pbot - rn_ptop ) 
  
! --- Climatological cloud work function

 actp  = 1.7
 facm  = 0.01                               
 actop = actp * facm

  a  = facm * a
  ph = 100.0*ph
   
  ac = 0.0     
  ad = 0.0     
  au = 0.0     
 tem = 0.0     

 tem(2:15) = ph(2:15) -  ph(1:14)
  au(2:15) =  a(1:14) / tem(2:15) 
  ad(2:15) =  a(2:15) / tem(2:15)

  ac(2:15) = ph(2:15) * au(2:15) - ph(1:14) * ad(2:15)
  ad(2:15) = ad(2:15) - au(2:15)

!---------------------------------------------------------------------
! --- initialize quantities for diagnostics output
!---------------------------------------------------------------------

   id_tdt_revap = register_diag_field ( mod_name, &
     'tdt_revap', axes(1:3), Time, &
     'Temperature tendency from RAS prec evap',      'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_revap = register_diag_field ( mod_name, &
     'qdt_revap', axes(1:3), Time, &
     'Spec humidity tendency from RAS prec evap',    'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_revap = register_diag_field ( mod_name, &
     'prec_revap', axes(1:2), Time, &
     'Evap Precip rate from RAS',                   'kg/m2/s' )

   id_snow_revap = register_diag_field ( mod_name, &
     'snow_revap', axes(1:2), Time, &
     'Evap Frozen precip rate from RAS',             'kg/m2/s' )

   id_prec_conv_3d = register_diag_field ( mod_name, &
     'prec_conv_3d', axes(1:3), Time, &
     '3D Precipitation rate from RAS',      'kg/m2/s',  &
                        missing_value=missing_value               )

   id_pcldb = register_diag_field ( mod_name, &
     'pcldb', axes(1:2), Time, &
     'Pressure at cloud base from RAS',              'hPa' )

!---------------------------------------------------------------------

   id_tdt_conv = register_diag_field ( mod_name, &
     'tdt_conv', axes(1:3), Time, &
     'Temperature tendency from RAS',                'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_conv = register_diag_field ( mod_name, &
     'qdt_conv', axes(1:3), Time, &
     'Spec humidity tendency from RAS',              'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_conv = register_diag_field ( mod_name, &
     'prec_conv', axes(1:2), Time, &
    'Precipitation rate from RAS',                  'kg/m2/s' )

   id_snow_conv = register_diag_field ( mod_name, &
     'snow_conv', axes(1:2), Time, &
    'Frozen precip rate from RAS',                  'kg/m2/s' )

   id_q_conv_col = register_diag_field ( mod_name, &
     'q_conv_col', axes(1:2), Time, &
    'Water vapor path tendency from RAS',           'kg/m2/s' )
   
   id_t_conv_col = register_diag_field ( mod_name, &
     't_conv_col', axes(1:2), Time, &
    'Column static energy tendency from RAS',       'W/m2' )


   id_mc = register_diag_field ( mod_name, &
     'mc', axes(half), Time, &
         'Cumulus Mass Flux from RAS',                   'kg/m2/s', &
                        missing_value=missing_value               )

   id_mfcb = register_diag_field ( mod_name, &
     'mfcb', axes(1:3), Time, &
       'Cumulus cloud-base mass flux from RAS', 'kg/m2/s', &
                        missing_value=missing_value) 

   id_mfct = register_diag_field ( mod_name, &
     'mfct', axes(1:3), Time, &
     'Cumulus cloud-top mass flux from RAS', 'kg/m2/s', &
                       missing_value=missing_value) 
       
   id_alm  = register_diag_field ( mod_name, &
     'alm', axes(1:3), Time, &
     'Cumulus entrainment rate from RAS', '1/m', &
                       missing_value=missing_value) 

   id_cfq_ras = register_diag_field ( mod_name,&
     'cfq_ras', axes(1:3), Time, &
     'Cumulus frequency from RAS', 'none', &
                       missing_value=missing_value) 

   id_qndt_conv = register_diag_field ( mod_name, &
     'qndt_conv', axes(1:3), Time, &
     'Cloud droplet tendency from RAS',              '#/kg/s',  &
                        missing_value=missing_value               )

   id_qldt_conv = register_diag_field ( mod_name, &
     'qldt_conv', axes(1:3), Time, &
     'Liquid water tendency from RAS',              'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qidt_conv = register_diag_field ( mod_name, &
     'qidt_conv', axes(1:3), Time, &
     'Ice water tendency from RAS',                 'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qadt_conv = register_diag_field ( mod_name, &
     'qadt_conv', axes(1:3), Time, &
     'Cloud fraction tendency from RAS',            '1/sec',    &
                        missing_value=missing_value               )

   id_ql_conv_col = register_diag_field ( mod_name, &
     'ql_conv_col', axes(1:2), Time, &
    'Liquid water path tendency from RAS',          'kg/m2/s' )
   
   id_qi_conv_col = register_diag_field ( mod_name, &
     'qi_conv_col', axes(1:2), Time, &
    'Ice water path tendency from RAS',             'kg/m2/s' )
   
   id_qa_conv_col = register_diag_field ( mod_name, &
     'qa_conv_col', axes(1:2), Time, &
    'Cloud mass tendency from RAS',                 'kg/m2/s' )
      
   id_det0        = register_diag_field ( mod_name, &
     'ras_det0', axes(1:3), Time, &
    'Detrained mass flux from RAS',                 'kg/m2/s' )
      
!----------------------------------------------------------------------
!    determine how many tracers are to be transported by ras_mod.
!----------------------------------------------------------------------
      if (present(tracers_in_ras)) then
        num_ras_tracers = count(tracers_in_ras)
      else
        num_ras_tracers = 0
      endif
      if (num_ras_tracers > 0) then
        do_ras_tracer = .true.
      else
        do_ras_tracer = .false.
      endif

!---------------------------------------------------------------------
!    allocate the arrays to hold the diagnostics for the ras tracers.
!---------------------------------------------------------------------
     allocate(id_tracer_conv    (num_ras_tracers)) ; id_tracer_conv = 0
     allocate(id_tracer_conv_col(num_ras_tracers)) ; id_tracer_conv_col = 0
     nn = 1
    do tr = 1,num_tracers
      if (tracers_in_ras(tr)) then
        call get_tracer_names(MODEL_ATMOS, tr, name=name, units=units)
 
!----------------------------------------------------------------------
!    for the column tendencies, the name for the diagnostic will be 
!    the name of the tracer followed by 'dt_RAS'. the longname will be 
!    the name of the tracer followed by ' tendency from RAS'. units are
!    the supplied units of the tracer divided by seconds.
!----------------------------------------------------------------------
        diagname = trim(name)//'dt_RAS'
        diaglname = trim(name)//' tendency from RAS'
        tendunits = trim(units)//'/s'
        id_tracer_conv(nn) = register_diag_field ( mod_name, &
                             trim(diagname), axes(1:3), Time, &
                             trim(diaglname), trim(tendunits),  &
                             missing_value=missing_value               )

!----------------------------------------------------------------------
!    for the column integral  tendencies, the name for the diagnostic 
!    will be the name of the tracer followed by 'dt_RAS_col'. the long-
!    name will be the name of the tracer followed by ' path tendency 
!    from RAS'. units are the supplied units of the tracer multiplied
!    by m**2 /kg divided by seconds.
!----------------------------------------------------------------------
      diagname = trim(name)//'dt_RAS_col'
      diaglname = trim(name)//' path tendency from RAS'
      tendunits = trim(units)//'m2/kg/s'
      id_tracer_conv_col(nn) = register_diag_field ( mod_name, &
                               trim(diagname), axes(1:2), Time, &
                               trim(diaglname), trim(tendunits),  &
                               missing_value=missing_value)
       nn = nn + 1
      endif
     end do

  module_is_initialized = .true.

!=====================================================================
  end SUBROUTINE RAS_INIT


!#####################################################################
subroutine ras_end

integer :: log_unit

log_unit = stdlog()
if ( mpp_pe() == mpp_root_pe() ) then
   write (log_unit,'(/,(a))') 'Exiting RAS.'
endif

module_is_initialized = .FALSE.


end subroutine ras_end

!#####################################################################

  SUBROUTINE RAS( is,     js,      Time,      temp0,   qvap0,     &
          uwnd0,  vwnd0,  pres0,   pres0_int, zhalf0,  coldT0,    &
                  dtime,  dtemp0,  dqvap0,    duwnd0,  dvwnd0,    &
                  rain3d, snow3d,                                 &
                  rain0,  snow0,   ras_tracers, qtrras,           &
                  mask,    kbot,  mc0, det0,                      & ! optional
                  ql0, qi0, qa0,                                  & ! optional
                  dl0, di0, da0,                                  & ! optional
                  qn0, dn0, do_strat, Aerosol)                      ! optional

!=======================================================================
! ***** DRIVER FOR RAS
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!     is, js    - Starting indices for window
!     Time      - Time used for diagnostics [time_type]
!     dtime     - Size of time step in seconds
!     pres0     - Pressure
!     pres0_int - Pressure at layer interface
!     zhalf0    - Height at layer interface
!     temp0     - Temperature
!     qvap0     - Water vapor 
!     uwnd0     - U component of wind
!     vwnd0     - V component of wind
!     coldT0    - should the precipitation assume to be frozen?
!     kbot      - OPTIONAL;lowest model level index (integer)
!     mask      - OPTIONAL;used only for diagnostic output
!     ras_tracers- prognostic tracers to move around
!                 note that R0 is assumed to be dimensioned
!                 (nx,ny,nz,nt), where nt is the number of tracers
!---------------------------------------------------------------------
! Arguments (Intent inout)
!     ql0       - OPTIONAL;cloud liquid
!     qi0       - OPTIONAL;cloud ice
!     qa0       - OPTIONAL;cloud/saturated volume fraction
!---------------------------------------------------------------------
  type(time_type), intent(in)                   :: Time
  integer,         intent(in)                   :: is, js
  real,            intent(in)                   :: dtime
  real,            intent(in), dimension(:,:,:) :: pres0, pres0_int, zhalf0
  real,            intent(inout), dimension(:,:,:) :: temp0, qvap0, uwnd0, vwnd0
  logical,         intent(in), dimension(:,:)   :: coldT0
  logical,         intent(in), optional         :: do_strat
  type(aerosol_type), intent (in), optional     :: Aerosol  
  real, intent(in) , dimension(:,:,:), OPTIONAL :: mask
  integer, intent(in), OPTIONAL, dimension(:,:) :: kbot
  real,  intent(inout), OPTIONAL,dimension(:,:,:)   :: ql0, qi0, qa0, qn0
  real,  intent(in), dimension(:,:,:,:)          :: ras_tracers
!---------------------------------------------------------------------
! Arguments (Intent out)
!       rain0  - surface rain
!       snow0  - surface snow
!       rain3d - 3D rain
!       snow3d - 3D snow
!       dtemp0 - Temperature change 
!       dqvap0 - Water vapor change 
!       duwnd0 - U wind      change 
!       dvwnd0 - V wind      change 
!       mc0    - OPTIONAL; cumulus mass flux
!       Dl0    - OPTIONAL; cloud liquid tendency
!       Di0    - OPTIONAL; cloud ice tendency
!       Da0    - OPTIONAL; cloud fraction tendency
!       qtrras - prognostic tracers tendency
!---------------------------------------------------------------------
  real, intent(out), dimension(:,:,:)           :: dtemp0, dqvap0, duwnd0, dvwnd0
  real, intent(out), dimension(:,:)             :: rain0,  snow0
  real, intent(out), dimension(:,:,:,:)         :: qtrras
  real, intent(out), dimension(:,:,:)           :: rain3d,snow3d
  
  real, intent(out), OPTIONAL, dimension(:,:,:) :: mc0
  real, intent(out), OPTIONAL, dimension(:,:,:) :: det0
  real, intent(out), OPTIONAL, dimension(:,:,:) :: dl0, di0, da0, dn0

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

! precipitation flux and evaporation profiles [kg/m2/sec]
 real, dimension(SIZE(temp0,3)) :: flxprec,flxprec_evap       ! sum of all clouds
 real, dimension(SIZE(temp0,3)) :: flxprec_ib,flxprec_ib_evap ! for each cloud

 real, parameter :: p00 = 1000.0E2

 logical :: coldT,  exist    
 real    :: precip, Hl, psfc, dpcu, dtinv
 integer :: ksfc,   klcl, kk

 integer, dimension(SIZE(temp0,3)) :: ic

 real, dimension(SIZE(temp0,3)) :: &
       temp, qvap, uwnd, vwnd, pres, dtemp, dqvap, duwnd, dvwnd, &    
       ql,   qi,   qa,   qn ,  Dl,   Di,    Da,    Dn,    mass,  pi,    theta, &
       cp_by_dp,   dqvap_sat,  qvap_sat,    alpha, beta,  gamma, &
       dtcu, dqcu, ducu, dvcu, Dlcu, Dicu,  Dacu,  Dncu

 real, dimension(SIZE(temp0      ,3),num_ras_tracers) :: tracer, dtracer, dtracercu
 real, dimension(SIZE(temp0,3)+1) ::  pres_int, mc, pi_int, mccu, zhalf
 real, dimension(SIZE(temp0,3)) ::  det                               

 logical, dimension(size(temp0,1),size(temp0,2)) :: rhtrig_mask
 integer, dimension(size(temp0,1),size(temp0,2)) :: kcbase

 real,    dimension(size(temp0,1),size(temp0,2)) :: &
          psfc0, t_parc, q_parc, p_parc, qs_parc   

 real,    dimension(size(temp0,1),size(temp0,2),size(temp0,3)) :: &
          qvap_sat0, dqvap_sat0, pmass

 real,    dimension(size(temp0,1),size(temp0,2),size(temp0,3)+1) :: mask3

 real,    dimension(size(temp0,1),size(temp0,2),size(temp0,3)) :: mfcb0, mfct0, alm0, cfq_ras !miz
 real,    dimension(SIZE(temp0,3)) ::  mfcb, mfct, alm !miz
 real  :: almx

real, dimension(size(temp0,1),size(temp0,2),size(temp0,3)+1) :: mc0_local

 logical :: setras, cloud_tracers_present, do_liq_num
 integer :: i, imax, j, jmax, k, kmax, tr, num_present
 integer :: ncmax, nc, ib, naer, na
 real    :: rasal, frac, zbase
 real    :: dpfac, dtcu_pbl, dqcu_pbl, ducu_pbl, dvcu_pbl, dtracercu_pbl(num_ras_tracers)
 
 real, dimension(size(temp0,1),size(temp0,2),size(temp0,3),3) :: totalmass1
 real, dimension(size(temp0,1),size(temp0,2),size(temp0,3)) :: airdens
 real, dimension(size(temp0,3),3) :: aerosolmass

 real :: thickness

!--- For extra diagnostics

 logical :: used
 real    :: dpevap, precip_ev

 real, dimension(SIZE(temp0,3)) :: &
       cup, dtevap, dqevap, dtemp_ev, dqvap_ev

 real, dimension(size(temp0,1),size(temp0,2)) :: &
       pcldb0, rain_ev0, snow_ev0, tempdiag

 real, dimension(size(temp0,1),size(temp0,2),size(temp0,3)) :: &
       cuprc3d, dtemp_ev0, dqvap_ev0

!=====================================================================

! --- Check to see if ras has been initialized
  if( .not. module_is_initialized ) CALL ERROR_MESG( 'RAS',  &
                                 'ras_init has not been called', FATAL )
  
  num_present = count((/present(ql0), present(qi0), present(qa0), present(Dl0), present(Di0), present(Da0)/))
  if(num_present == 0) then
    cloud_tracers_present = .false.
  else if(num_present == 6) then
    cloud_tracers_present = .true.
  else
    call ERROR_MESG( 'RAS','Either all or none of the cloud tracers and their tendencies'// &
                     ' must be present in the call to subroutine ras',FATAL)
  endif

  if(cloud_tracers_present .and. .not.present(mc0)) then
    call ERROR_MESG( 'RAS','mc0 must be present when cloud tracers are present',FATAL)
  endif


! --- Set dimensions
  imax  = size( temp0, 1 )
  jmax  = size( temp0, 2 )
  kmax  = size( temp0, 3 )

  
  do_liq_num = PRESENT(qn0) ! Test for presence of liquid droplet number array

  if ( do_liq_num ) then
    if (.not.(PRESENT(dn0)) .or. .not.(PRESENT(Aerosol))) &
      call ERROR_MESG( 'RAS','dn0 and Aerosol must be present when liquid droplet number are present',FATAL)
  
     naer = size(Aerosol%aerosol,4)

 if(use_online_aerosol) then
 
    do k = 1,kmax
      do j = 1,jmax
        do i = 1,imax
          if(pres0_int(i,j,k)<1.) then
            thickness=(pres0_int(i,j,k+1)-pres0_int(i,j,k))* &
             8.314*temp0(i,j,k)/(9.8*0.02888*pres0_int(i,j,k))
          else
            thickness=log(pres0_int(i,j,k+1)/ &
             pres0_int(i,j,k))*8.314*temp0(i,j,k)/(9.8*0.02888)
          end if
         do na = 1,naer
            if(Aerosol%aerosol_names(na) == 'so4' .or. &
                Aerosol%aerosol_names(na) == 'so4_anthro' .or. Aerosol%aerosol_names(na) == 'so4_natural') then
                        totalmass1(i,j,k,1)=totalmass1(i,j,k,1)+Aerosol%aerosol(i,j,k,na)
            else if(Aerosol%aerosol_names(na) == 'omphilic' .or. &
            Aerosol%aerosol_names(na) == 'omphobic') then
                        totalmass1(i,j,k,3)=totalmass1(i,j,k,3)+Aerosol%aerosol(i,j,k,na)
                 else if(Aerosol%aerosol_names(na) == 'seasalt1' .or. &
                  Aerosol%aerosol_names(na) == 'seasalt2') then
                        totalmass1(i,j,k,2)=totalmass1(i,j,k,2)+Aerosol%aerosol(i,j,k,na)
                 end if
         end do
  
         do na = 1, 3
              totalmass1(i,j,k,na)=totalmass1(i,j,k,na)/thickness*1.0e9*1.0e-12
          end do
  
         end do
        end do
      end do
  
    else
    
    do k = 1,kmax
      do j = 1,jmax
        do i = 1,imax
          if(pres0_int(i,j,k)<1.) then
            thickness=(pres0_int(i,j,k+1)-pres0_int(i,j,k))* &
            8.314*temp0(i,j,k)/(9.8*0.02888*pres0_int(i,j,k))
          else
            thickness=log(pres0_int(i,j,k+1)/ &
            pres0_int(i,j,k))*8.314*temp0(i,j,k)/(9.8*0.02888)
          end if
          totalmass1(i,j,k,1)=(Aerosol%aerosol(i,j,k,1)+Aerosol%aerosol(i,j,k,2))  &
          /thickness*1.0e9*1.0e-12
          totalmass1(i,j,k,2)=sea_salt_scale*Aerosol%aerosol(i,j,k,5)  &
          /thickness*1.0e9*1.0e-12
          totalmass1(i,j,k,3)=om_to_oc*Aerosol%aerosol(i,j,k,3)  &
          /thickness*1.0e9*1.0e-12
        end do
      end do
    end do
  
end if

    airdens = pres0 / (rdgas * temp0 * (1.   - ql0 - qi0) )
  endif   

 
! --- Initalize

   dtemp0 = 0.0                               
   dqvap0 = 0.0                               
   duwnd0 = 0.0                               
   dvwnd0 = 0.0                                                                        
    rain0 = 0.0   
    snow0 = 0.0

  dtemp_ev0 = 0.0
  dqvap_ev0 = 0.0
   rain_ev0 = 0.0
   snow_ev0 = 0.0
    cuprc3d = 0.0
     rain3d = 0.0
     snow3d = 0.0
      mfcb0 = 0.0 
      mfct0 = 0.0 
      alm0  = 0.0
    cfq_ras = 0.0

  if ( cloud_tracers_present ) then
      Da0 = 0.0
      Dl0 = 0.0
      Di0 = 0.0
      if ( do_liq_num ) Dn0 = 0.0
  end if
  if (present(mc0) ) then
      mc0 = 0.0
      det0 = 0.0
  end if
! Initialize the tracer tendencies
    qtrras = 0.0

  do k=1,kmax
    pmass(:,:,k) = (pres0_int(:,:,k+1)-pres0_int(:,:,k))/GRAV
  end do

    frac  = fracs  / dtime
    rasal = rasal0 / dtime 
       
! --- Compute saturation value of water vapor & its derivative

    call compute_qs (temp0, pres0, qvap_sat0, dqsdT=dqvap_sat0)


! --- Find LCL ---> cloud base

  if (present(kbot ) ) then
     do j = 1,jmax
     do i = 1,imax
        k = kbot(i,j)
          t_parc(i,j) =     temp0(i,j,k)
          q_parc(i,j) =     qvap0(i,j,k)
          p_parc(i,j) =     pres0(i,j,k)
         qs_parc(i,j) = qvap_sat0(i,j,k)
     end do
     end do
  else
          t_parc(:,:) =     temp0(:,:,kmax)
          q_parc(:,:) =     qvap0(:,:,kmax)
          p_parc(:,:) =     pres0(:,:,kmax)
         qs_parc(:,:) = qvap_sat0(:,:,kmax)
  end if

  q_parc  = MAX( q_parc, 1.0E-6   )
  q_parc  = MIN( q_parc, qs_parc )

  CALL COMP_LCL( t_parc, q_parc, p_parc, pres0, kcbase )

! --- set rh trigger
  rhtrig_mask(:,:) = q_parc(:,:) >= rh_trig*qs_parc(:,:) 

! --- Set surface pressure
  if (present(kbot ) ) then
     do j = 1,jmax
     do i = 1,imax
        k = kbot(i,j) + 1
        psfc0(i,j) = pres0_int(i,j,k)
     end do
     end do
  else
        psfc0(:,:) = pres0_int(:,:,kmax+1)
  end if

! --- Save cloud base pressure
  if ( id_pcldb > 0 ) then
     do j = 1,jmax
     do i = 1,imax
     klcl = kcbase(i,j)
            pcldb0(i,j) = pres0_int(i,j,klcl)
     end do
     end do
  end if

  mc0_local = 0.

!---------------------------------------------------------------------

! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
  do j = 1,jmax
  do i = 1,imax
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    flxprec(:) = 0.
    flxprec_evap(:) = 0.

! --- Set order in which clouds are to be done                   
  CALL RAS_CLOUD_INDEX ( ic, ncmax )

! --- rh trigger
  if ( .not.rhtrig_mask(i,j) ) CYCLE

! --- Pack input variables / Initalize output variables 

      temp(:) =      temp0(i,j,:)
      qvap(:) =      qvap0(i,j,:)
      uwnd(:) =      uwnd0(i,j,:)
      vwnd(:) =      vwnd0(i,j,:)
      pres(:) =      pres0(i,j,:)
  pres_int(:) =  pres0_int(i,j,:)
     zhalf(:) =     zhalf0(i,j,:)
  qvap_sat(:) =  qvap_sat0(i,j,:) 
 dqvap_sat(:) = dqvap_sat0(i,j,:) 
      psfc    =      psfc0(i,j)
     coldT    =     coldT0(i,j)
      klcl    =     kcbase(i,j)
      if ( do_liq_num ) &
        aerosolmass(:,:) = totalmass1(i,j,:,:)

     dtemp(:) = 0.0                               
     dqvap(:) = 0.0                               
     duwnd(:) = 0.0                               
     dvwnd(:) = 0.0  
    precip    = 0.0   

          cup(:) = 0.0
     dtemp_ev(:) = 0.0 
     dqvap_ev(:) = 0.0 
    precip_ev    = 0.0 
         mfcb(:) = 0.0
         mfct(:) = 0.0
         alm (:) = 0.0
         almx    = 0.0

  if ( cloud_tracers_present ) then
        qa(:) = qa0(i,j,:)
        ql(:) = ql0(i,j,:)
        qi(:) = qi0(i,j,:)
        Da(:) = 0.0
        Dl(:) = 0.0
        Di(:) = 0.0
        if ( do_liq_num ) then
          qn(:) = qn0(i,j,:)
          Dn(:) = 0.0
        endif     
  end if
  if ( present(mc0) .or. id_mc > 0 ) then
        mc(:) = 0.0
        det(:) = 0.0
  end if

! Get the column of tracer data
    do tr = 1,num_ras_tracers
      dtracer(:,tr) = 0.0
      tracer(:,tr)  = ras_tracers(i,j,:,tr)
    enddo

    if( coldT ) then
     Hl = HLs
  else  
     Hl = HLv
  end if

  if( PRESENT( kbot ) ) then
       ksfc = kbot(i,j) + 1
  else
       ksfc = kmax + 1
  endif

! --- Thickness of layers 
  mass(1:kmax) =  pres_int(2:kmax+1) - pres_int(1:kmax) 
  mass(1:kmax) = MAX( mass(1:kmax), 1.0e-5 )

! --- Sub-cloud layer thickness
   zbase = zhalf(klcl) - zhalf(ksfc)

! --- Compute exner functions 
! --- at layer interfaces
  pi_int(:) = ( pres_int(:) / p00  ) ** Kappa                                   
! --- at full levels
  pi(1:kmax) = ( pi_int(2:kmax+1) * pres_int(2:kmax+1)     &
               - pi_int(1:kmax  ) * pres_int(1:kmax  ) )   &
               / ( mass(1:kmax  ) * one_plus_kappa )
  pi(1:kmax) =  MAX( pi(1:kmax), 1.0e-5 )

! --- Compute potential temperature
  theta(:) = temp(:) / pi(:)                                 

! --- Compute Cp divided by dpres                  
  cp_by_dp(1:kmax) = Cp_Air / mass(1:kmax)

! --- Compute mass of layers              
  mass(:) = mass(:) / Grav

  setras  = .true.

! --- Test for convective instability in column             
 CALL RAS_CLOUD_EXIST( klcl, qvap,   qvap_sat, theta,    &
                       pi,   pi_int, ic,       ncmax,    &
                       Hl,   exist  )
 if ( .not. exist ) CYCLE

!---------------------------------------------------------------------
! Cloud top loop starts
!---------------------------------------------------------------------
  do nc = 1,ncmax
!---------------------------------------------------------------------

     ib  = ic(nc)
 if( ib >= klcl) CYCLE

 if ( setras ) then
! --- Compute some stuff
     alpha(:) =  qvap_sat(:) - dqvap_sat(:) * temp(:)
      beta(:) = dqvap_sat(:) * pi(:)
     gamma(:) = 1.0 / ( ( 1.0 + ( Hl * dqvap_sat(:) / Cp_Air ) ) * pi(:) )
 endif

! --- Do adjustment
  if ( cloud_tracers_present .and. .not.(do_liq_num)) then
  CALL RAS_CLOUD(temp0(i,j,:), pres0(i,j,:),airdens(i,j,:),          &
       klcl,  ib,   rasal, frac, Hl, coldT,                          &
       theta, qvap, uwnd,  vwnd, pres_int, pi_int, pi, psfc,         &
       alpha, beta, gamma, cp_by_dp, zbase, almx,                    &
       dtcu,  dqcu, ducu,  dvcu, dpcu, tracer, dtracercu,            &
       mccu, ql, qi, qa, Dlcu, Dicu, Dacu)
  else if ( cloud_tracers_present .and. do_liq_num) then
  CALL RAS_CLOUD(temp0(i,j,:), pres0(i,j,:),airdens(i,j,:),          &
       klcl,  ib,   rasal, frac, Hl, coldT,                          &
       theta, qvap, uwnd,  vwnd, pres_int, pi_int, pi, psfc,         &
       alpha, beta, gamma, cp_by_dp, zbase, almx,                    &
       dtcu,  dqcu, ducu,  dvcu, dpcu, tracer, dtracercu,            &
       mccu, ql, qi, qa, Dlcu, Dicu, Dacu, aerosolmass, qn, Dncu)
  else if (present(mc0) .or. id_mc > 0) then
  CALL RAS_CLOUD(temp0(i,j,:), pres0(i,j,:),airdens(i,j,:),          &
       klcl,  ib,   rasal, frac, Hl,  coldT,                         &
       theta, qvap, uwnd,  vwnd, pres_int, pi_int, pi, psfc,         &
       alpha, beta, gamma, cp_by_dp, zbase, almx,                    &
       dtcu,  dqcu, ducu,  dvcu, dpcu, tracer, dtracercu, mccu)  
  else
  CALL RAS_CLOUD(temp0(i,j,:), pres0(i,j,:),airdens(i,j,:),          &
       klcl,  ib,   rasal, frac, Hl,  coldT,                         &
       theta, qvap, uwnd,  vwnd, pres_int, pi_int, pi, psfc,         &
       alpha, beta, gamma, cp_by_dp, zbase, almx,                    &
       dtcu,  dqcu, ducu,  dvcu, dpcu, tracer, dtracercu)
  end if

! --- For optional diagnostic output
  cup(ib) = dpcu

! --- Multiply tendencies by size of time step
  dtcu(:) =  dtcu(:) * dtime
  dqcu(:) =  dqcu(:) * dtime
  ducu(:) =  ducu(:) * dtime
  dvcu(:) =  dvcu(:) * dtime
  dpcu    =  dpcu    * dtime

  if ( cloud_tracers_present ) then
  Dacu(:) =  Dacu(:) * dtime
  Dlcu(:) =  Dlcu(:) * dtime
  Dicu(:) =  Dicu(:) * dtime
  if ( do_liq_num ) &
    Dncu(:) =  Dncu(:) * dtime
  end if

    do tr = 1,num_ras_tracers
      dtracercu(:,tr) = dtracercu(:,tr) * dtime
    enddo
  if( modify_pbl ) then   
! TTTTTTTTTTTTTTTTTTTTT      
! --- Adjust mass flux between cloud base and surface       
    if ( present(mc0) .or. id_mc > 0 ) then
      do k = klcl+1,ksfc
        mccu(k) = mccu(klcl) * ( psfc - pres_int(k)    ) / &
                               ( psfc - pres_int(klcl) )
      end do
    end if
! --- Adjust tendencies between cloud base and surface
    dpfac = ( pres_int(klcl+1) - pres_int(klcl) ) / &
             ( psfc             - pres_int(klcl) )   
    dtcu_pbl = dpfac * dtcu(klcl) * pi(klcl)
    dqcu_pbl = dpfac * dqcu(klcl) 
    ducu_pbl = dpfac * ducu(klcl) 
    dvcu_pbl = dpfac * dvcu(klcl) 
    do tr = 1,num_ras_tracers
      dtracercu_pbl(tr) = dpfac * dtracercu(klcl,tr)
    enddo
    do k = klcl,ksfc-1
      dtcu(k) = dtcu_pbl / pi(k) 
      dqcu(k) = dqcu_pbl
      ducu(k) = ducu_pbl
      dvcu(k) = dvcu_pbl
        do tr = 1,num_ras_tracers
          dtracercu(k,tr) = dtracercu_pbl(tr)
        enddo
    end do
! LLLLLLLLLLLLLLLLL
  endif
! --- Evaporate some precip

      dtevap(:) = 0.0
      dqevap(:) = 0.0
      dpevap    = 0.0
!
! -- initialise the precipitation and evaporation flux
      flxprec_ib_evap = 0.0

  if( evap_on .and. ( dpcu > 0.0 ) ) then

  CALL RAS_CEVAP ( ib, temp, qvap, pres, mass, qvap_sat,       &
                   dqvap_sat, psfc, Hl, dtime, ksfc, dpcu,     &
                   dtevap, dqevap,  dpevap, flxprec_ib, flxprec_ib_evap  )

      dtcu(:) =  dtcu(:) + dtevap(:) / pi(:)
      dqcu(:) =  dqcu(:) + dqevap(:)
      dpcu    = MAX(dpcu - dpevap, 0.)

  else
    flxprec_ib(1:ib-1) = 0.0
    flxprec_ib(ib:) = dpcu/dtime
  endif

!---sum up precipitation flux and evaporation from  each cloud type
!
flxprec      = flxprec +flxprec_ib
flxprec_evap = flxprec_evap +flxprec_ib_evap

! --- Update prognostic tracers
!     NOTE: negative values of tracers are not prevented
    do tr = 1,num_ras_tracers
         tracer(:,tr)   = amax1(0.,tracer(:,tr) + dtracercu(:,tr))
    enddo  

if ( prevent_unreasonable ) then
! --- Update cloud liquid, ice, and fraction
!     NOTE: unreasonable states are prevented

  if ( cloud_tracers_present ) then
  
    !cloud fraction---------------
    where ((qa + Dacu) .lt. 0.)
         Dacu = -1.*qa
        qa   = 0.
    elsewhere
         qa   = qa + Dacu
    end where
    where (qa .gt. 1.)
         Dacu = Dacu + (1. - qa)
         qa   = 1.
    end where
    
    !cloud liquid----------------
    where ((ql + Dlcu) .lt. 0.)
         dtcu = dtcu - HLv*(ql+Dlcu)/Cp_Air/pi
        dqcu = dqcu + (ql + Dlcu)
         Dlcu = Dlcu - (ql + Dlcu)
        ql   = 0.
    elsewhere
         ql   = ql + Dlcu
    end where
    
    !cloud ice----------------
    where ((qi + Dicu) .lt. 0.)
         dtcu = dtcu - HLs*(qi+Dicu)/Cp_Air/pi
        dqcu = dqcu + (qi + Dicu)
         Dicu = Dicu - (qi + Dicu)
        qi   = 0.
    elsewhere
         qi   = qi + Dicu
    endwhere
            
  end if
else

  if ( cloud_tracers_present ) then
    ql(:) = ql(:) + Dlcu(:)
    qi(:) = qi(:) + Dicu(:)
    qa(:) = qa(:) + Dacu(:)
    if ( do_liq_num ) &
      qn(:) = qn(:) + Dncu(:)
  end if

endif

! --- Update potential temperature, water vapor, winds
  theta(:) = theta(:) + dtcu(:)
   qvap(:) =  qvap(:) + dqcu(:)
   uwnd(:) =  uwnd(:) + ducu(:)
   vwnd(:) =  vwnd(:) + dvcu(:)

! --- Recover temperature 
  temp(:) = theta(:) * pi(:)

! --- Accumulate precip 
  precip    = precip    + dpcu 
  precip_ev = precip_ev + dpevap

! --- Accumulate tendencies 
     dtemp(:) =    dtemp(:) +   dtcu(:) * pi(:)
     dqvap(:) =    dqvap(:) +   dqcu(:)
     duwnd(:) =    duwnd(:) +   ducu(:) 
     dvwnd(:) =    dvwnd(:) +   dvcu(:)
  dtemp_ev(:) = dtemp_ev(:) + dtevap(:)
  dqvap_ev(:) = dqvap_ev(:) + dqevap(:)

  if ( cloud_tracers_present ) then
    Da(:) = Da(:) + Dacu(:)
    Dl(:) = Dl(:) + Dlcu(:)
    Di(:) = Di(:) + Dicu(:)
    if ( do_liq_num ) &
      Dn(:) = Dn(:) + Dncu(:)
  end if
  mfcb(ib)=mccu(klcl) 
  mfct(ib)=mccu(ib+1) 
  alm (ib)=almx 

  if ( present(mc0) .or. id_mc > 0 ) then
    mc(:) = mc(:) + mccu(:)
    det(ib) = det(ib) + mccu(ib+1)
  end if
  do tr = 1,num_ras_tracers
    dtracer(:,tr) = dtracer(:,tr) + dtracercu(:,tr)
  enddo  

  setras = .false.

 if ( setras ) then
! --- Re-Compute saturation value of water vapor & its derivative
   call compute_qs (temp, pres, qvap_sat, dqsdT=dqvap_sat)
  end if
     
!---------------------------------------------------------------------
  end do ! do nc = 1,ncmax
!---------------------------------------------------------------------
! Cloud top loop ends
!---------------------------------------------------------------------

! --- Unpack variables
      dtemp0(i,j,:) = dtemp(:)
      dqvap0(i,j,:) = dqvap(:) 
      duwnd0(i,j,:) = duwnd(:)
      dvwnd0(i,j,:) = dvwnd(:)
      mfcb0 (i,j,:) = mfcb(:)
      mfct0 (i,j,:) = mfct(:)
      alm0  (i,j,:) = alm (:)

   if ( coldT ) then
      snow0(i,j) = precip
   else
      rain0(i,j) = precip
   end if

!-- unpacking the precipitation flux. The 
    do kk = 1,kmax
      if ( coldT ) then
       snow3d(i,j,kk+1)   = flxprec(kk)      !kg/m2/sec
      else
       rain3d(i,j,kk+1)   = flxprec(kk)      !kg/m2/sec
      end if
    enddo

  if ( cloud_tracers_present ) then
        Da0(i,j,:) = Da(:) 
        Dl0(i,j,:) = Dl(:) 
        Di0(i,j,:) = Di(:)
        if ( do_liq_num ) &
          Dn0(i,j,:) = Dn(:)
  end if
  if ( present(mc0) .or. id_mc > 0 ) then
        mc0_local(i,j,:) = mc(:)
        det0(i,j,:) = det(:)
  end if
  do tr = 1,num_ras_tracers
    qtrras(i,j,:,tr) = dtracer(:,tr)
  enddo    

! --- For extra diagnostics
 
     cuprc3d(i,j,:) =       cup(:) 
   dtemp_ev0(i,j,:) =  dtemp_ev(:)
   dqvap_ev0(i,j,:) =  dqvap_ev(:)

    if ( coldT ) then
       snow_ev0(i,j) = precip_ev
    else
       rain_ev0(i,j) = precip_ev
    end if

! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
  end do ! do i = 1,imax
  end do ! do j = 1,jmax
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

!------- update input values and compute tendency -----------
   dtinv = 1.0/dtime

   temp0=temp0+dtemp0;    qvap0=qvap0+dqvap0
   uwnd0=uwnd0+duwnd0;    vwnd0=vwnd0+dvwnd0

   duwnd0=duwnd0*dtinv; dvwnd0=dvwnd0*dtinv

   dtemp0=dtemp0*dtinv; dqvap0=dqvap0*dtinv
   rain0=rain0*dtinv; snow0=snow0*dtinv

!!------- add on tendency ----------
!      tdt=tdt+ttnd; qdt=qdt+qtnd
!      udt=udt+utnd; vdt=vdt+vtnd

!------- update input values , compute and add on tendency -----------
!-------              in the case of strat                 -----------

   if (cloud_tracers_present) then
     ql0(:,:,:) = ql0(:,:,:)+Dl0(:,:,:)
     qi0(:,:,:) = qi0(:,:,:)+Di0(:,:,:)
     qa0(:,:,:) = qa0(:,:,:)+Da0(:,:,:)
     if ( do_liq_num ) &
       qn0(:,:,:) = qn0(:,:,:)+Dn0(:,:,:)

     Dl0(:,:,:) = Dl0(:,:,:)*dtinv
     Di0(:,:,:) = Di0(:,:,:)*dtinv
     Da0(:,:,:) = Da0(:,:,:)*dtinv
     if ( do_liq_num ) &
       Dn0(:,:,:)=Dn0(:,:,:)*dtinv

   end if
   do tr = 1,num_ras_tracers
     qtrras(:,:,:,tr) = qtrras(:,:,:,tr)*dtinv
   enddo    

   if(present(mc0)) mc0 = mc0_local

!---------------------------------------------------------------------
! --- Extra diagnostics
!---------------------------------------------------------------------
 
   dtemp_ev0(:,:,:) = dtemp_ev0(:,:,:) * dtinv
   dqvap_ev0(:,:,:) = dqvap_ev0(:,:,:) * dtinv
   snow_ev0(:,:)   =  snow_ev0(:,:)   * dtinv
   rain_ev0(:,:)   =  rain_ev0(:,:)   * dtinv

   used = send_data ( id_tdt_revap, dtemp_ev0, Time, is, js, 1 )
   used = send_data ( id_qdt_revap, dqvap_ev0, Time, is, js, 1 )
   used = send_data ( id_prec_revap, rain_ev0+snow_ev0, Time, is, js )
   used = send_data ( id_snow_revap, snow_ev0, Time, is, js )
   used = send_data ( id_prec_conv_3d, cuprc3d, Time, is, js, 1 )
   used = send_data ( id_pcldb, pcldb0, Time, is, js )

!------- diagnostics for dt/dt_ras -------
   used = send_data ( id_tdt_conv, dtemp0, Time, is, js, 1, &
                      rmask=mask )
!------- diagnostics for dq/dt_ras -------
   used = send_data ( id_qdt_conv, dqvap0, Time, is, js, 1, &
                      rmask=mask )
!------- diagnostics for precip_ras -------
   used = send_data ( id_prec_conv, (rain0+snow0), Time, is, js )
!------- diagnostics for snow_ras -------
   used = send_data ( id_snow_conv, snow0, Time, is, js )
!------- diagnostics for cumulus mass flux from ras -------
   if ( id_mc > 0 ) then
        !------- set up mask --------
        mask3(:,:,1:(kmax+1)) = 1.
        if (present(mask)) then
          WHERE (mask(:,:,1:kmax) <= 0.5)
                  mask3(:,:,1:kmax) = 0.
          END WHERE
        endif
        used = send_data ( id_mc, mc0_local, Time, is, js, 1, rmask=mask3 )
   endif
   used = send_data ( id_mfcb, mfcb0(:,:,1:kmax), Time, is, js, 1, rmask=mask )
   used = send_data ( id_mfct, mfct0(:,:,1:kmax), Time, is, js, 1, rmask=mask )
   used = send_data ( id_alm,  alm0 (:,:,1:kmax), Time, is, js, 1, rmask=mask )
   if ( id_cfq_ras > 0 ) then
     where (mfcb0(:,:,:) .gt. 0.)
       cfq_ras(:,:,:) = 1
     end where
     used = send_data ( id_cfq_ras,  cfq_ras (:,:,1:kmax), Time, is, js, 1, rmask=mask )
   endif
    

!------- diagnostics for water vapor path tendency ----------
   if ( id_q_conv_col > 0 ) then
     tempdiag(:,:)=0.
     do k=1,kmax
       tempdiag(:,:) = tempdiag(:,:) + dqvap0(:,:,k)*pmass(:,:,k)*dtinv
     end do
     used = send_data ( id_q_conv_col, tempdiag, Time, is, js )
   end if

!------- diagnostics for dry static energy tendency ---------
   if ( id_t_conv_col > 0 ) then
     tempdiag(:,:)=0.
     do k=1,kmax
       tempdiag(:,:) = tempdiag(:,:) + dtemp0(:,:,k)*Cp_Air*pmass(:,:,k)
     end do
     used = send_data ( id_t_conv_col, tempdiag, Time, is, js )
   end if

   if ( cloud_tracers_present ) then

      !------- diagnostics for dql/dt from RAS -----------------
      used = send_data ( id_qldt_conv, Dl0, Time, is, js, 1, &
                         rmask=mask )
      
      !------- diagnostics for dqi/dt from RAS -----------------
      used = send_data ( id_qidt_conv, Di0, Time, is, js, 1, &
                         rmask=mask )
      
      !------- diagnostics for dqa/dt from RAS -----------------
      used = send_data ( id_qadt_conv, Da0, Time, is, js, 1, &
                         rmask=mask )

      !------- diagnostics for dqn/dt from RAS -----------------
      if (do_liq_num .and. id_qndt_conv > 0 ) &
        used = send_data ( id_qndt_conv, Dn0, Time, is, js, 1, &
                         rmask=mask )

      !------- diagnostics for liquid water path tendency ------
      if ( id_ql_conv_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kmax
          tempdiag(:,:) = tempdiag(:,:) + Dl0(:,:,k)*pmass(:,:,k)
        end do
        used = send_data ( id_ql_conv_col, tempdiag, Time, is, js )
      end if
      
      !------- diagnostics for ice water path tendency ---------
      if ( id_qi_conv_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kmax
          tempdiag(:,:) = tempdiag(:,:) + Di0(:,:,k)*pmass(:,:,k)
        end do
        used = send_data ( id_qi_conv_col, tempdiag, Time, is, js )
      end if
      
      !---- diagnostics for column integrated cloud mass tendency ---
      if ( id_qa_conv_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kmax
          tempdiag(:,:) = tempdiag(:,:) + Da0(:,:,k)*pmass(:,:,k)
        end do
        used = send_data ( id_qa_conv_col, tempdiag, Time, is, js )
      end if
         
      !------- diagnostics for det0 from RAS -----------------
      used = send_data ( id_det0, det0, Time, is, js, 1, &
                         rmask=mask )

   end if !end do strat if

   do tr = 1, num_ras_tracers
      !------- diagnostics for dtracer/dt from RAS -------------
      used = send_data ( id_tracer_conv(tr), qtrras(:,:,:,tr), Time, is, js, 1, &
                         rmask=mask )

      !------- diagnostics for column tracer path tendency -----
      if ( id_tracer_conv_col(tr) > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kmax
          tempdiag(:,:) = tempdiag(:,:) + qtrras(:,:,k,tr)*pmass(:,:,k)
        end do
        used = send_data ( id_tracer_conv_col(tr), tempdiag, Time, is, js )
      end if


   enddo
!=====================================================================



  end SUBROUTINE RAS


!#######################################################################
!#######################################################################

 SUBROUTINE RAS_CLOUD(t, p, dens,                                &
            k, ic, rasal, frac, hl, coldT,                       &
            theta, qvap, uwnd, vwnd, pres_int, pi_int, pi, psfc, &
            alf, bet, gam, cp_by_dp, zbase, almx,                &
            dtcu, dqcu, ducu,  dvcu, dpcu, tracer, dtracercu,    &
            mccu, ql, qi, qa, Dlcu, Dicu, Dacu, aerosolmass, qn, Dncu)   ! optional
!=======================================================================
! RAS Cu Parameterization 
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!     k        : Cloud base index.
!     ic       : Cloud top index.
!     hl       : proper latent heat for the column
!     coldT    : is the precipitation frozen?
!     theta    : Potential temperature
!     qvap     : Water vapor 
!     uwnd     : U component of wind
!     vwnd     : V component of wind
!     pres_int : Pressure       at layer interface
!     pi_int   : Exner function at layer interface
!     pi       : Exner function 
!     psfc     : Surface pressure
!     zbase    : Thickness of sub-cloud layer
!     rasal    : Relaxation parameter for cloud type ic
!     ql       : OPTIONAL, cloud liquid
!     qi       : OPTIONAL, cloud ice
!     qa       : OPTIONAL, cloud fraction or saturated volume fraction
!---------------------------------------------------------------------

 real,    intent(inout):: almx    !miz
 real,    intent(in) :: rasal, frac, zbase
 real,    intent(in) :: hl,    psfc
 integer, intent(in) :: ic,    k
 logical, intent(in) :: coldT
 real, intent(in), dimension(:) :: t,p,dens
 real, intent(in), dimension(:) :: theta, qvap, uwnd, vwnd, pres_int, pi_int
 real, intent(in), dimension(:) :: alf,   bet,  gam,  pi,   cp_by_dp
 real, intent(in), OPTIONAL, dimension(:) :: ql,qi,qa,qn
 real, intent(in), dimension(:,:) :: tracer
 real, intent(in), OPTIONAL, dimension(:,:) :: aerosolmass
!---------------------------------------------------------------------
! Arguments (Intent out)
!     dpcu    : Precip for cloud type ic.
!     dtcu    : Potential temperature change
!     dqcu    : Water vapor change
!     ducu    : Cu friction - u component.
!     dvcu    : Cu friction - v component.
!     mccu    : OPTIONAL ; Cumulus Mass Flux
!     Dacu    : OPTIONAL ; Detrained saturated mass fraction tendency
!     Dlcu    : OPTIONAL ; Detrained cloud liquid tendency
!     Dicu    : OPTIONAL ; Detrained cloud ice tendency
!     dtracercu : OPTIONAL ; Detrained tracer tendency 
!---------------------------------------------------------------------

 real, intent(out)  :: dpcu
 real, intent(out), dimension(:) :: dtcu, dqcu, ducu, dvcu
 real, intent(out), OPTIONAL, dimension(:) :: mccu, Dacu, Dlcu, Dicu, Dncu
 real, intent(out), dimension(:,:) :: dtracercu

!---------------------------------------------------------------------
!    (Intent local)
!---------------------------------------------------------------------

 real, parameter :: rhmax  = 0.9999 

  logical :: lcase1,    lcase2
  real    :: wfn_crit,  ftop, rn_frac
  real    :: wfn,  akm, qs1,  uht, vht, wlq, alm 
  real    :: rasalf
  real    :: wll,  wli, wlN
  integer :: km1,  ic1, l,    iwk
  real    :: xx1,  xx2, xx3,  xx4
  real    :: ssl, dtemp, zzl, hccp, hcc,  dpib, dpit 
  real    :: dut, dvt,   dub, dvb,  wdet, dhic, hkb, hic, sic
 
 real, dimension(size(theta,1)) :: gmh, eta, hol, hst, qol
 real, dimension(SIZE(tracer,2)) :: wlR

 logical :: Ldacu, Lmccu, LRcu, do_liq_num

!thetac in-cloud potential temperature (K)
!qc in-cloud vapor mixing ratio (kg water/kg air)
!qt in-cloud qc + ql (kg water/kg air)
 real, dimension(size(theta,1)) :: thetac, qc, qt
 real     :: tc, te, Nc, up_conv, drop
 real, dimension(3) :: totalmass
 integer :: tr

!=====================================================================

! --- Check for presence of optional arguments
  Ldacu = PRESENT( Dacu )
  Lmccu = PRESENT( mccu ) 
  LRcu  = (num_ras_tracers > 0 ) !.TRUE.
  do_liq_num = PRESENT(qn) 


! Initialize
  dtcu = 0.0
  dqcu = 0.0
  ducu = 0.0
  dvcu = 0.0
  dpcu = 0.0
  if ( Ldacu ) then
  Dacu = 0.0
  Dlcu = 0.0
  Dicu = 0.0
  if ( do_liq_num ) Dncu = 0.0
  end if
  if ( Lmccu ) then
  mccu = 0.0
  end if
  if ( LRcu ) then
  do tr = 1,num_ras_tracers
    dtracercu(:,tr) = 0.0
  enddo
  end if

! Initialize
  wfn=0.0
  akm=0.0
  qs1=0.0
  uht=0.0
  vht=0.0
  wlq=0.0
  alm=0.0
  almx=0.0
  wll=0.0
  wli=0.0
  if ( do_liq_num ) &
    wlN=0.0
  gmh=0.0
  eta=0.0
  hol=0.0
  hst=0.0
  qol=0.0
  km1=0
  ic1=0
  l=0
  rasalf=0.0

  km1 = k  - 1
  ic1 = ic + 1

!=====================================================================
! --- RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL
!=====================================================================

! --- at cloud base
    qs1    = alf(k) + bet(k) * theta(k)
    qol(k) = MIN( qs1 * rhmax, qvap(k) )
    hol(k) = pi_int(k+1) * theta(k) * Cp_Air + qol(k) * Hl
    eta(k) = 0.0
    zzl    = ( pi_int(k+1) - pi_int(k) ) * theta(k) * Cp_Air

! --- between cloud base & cloud top
 if ( ic < km1 ) then
 do l = km1,ic1,-1
    qs1    = alf(l) + bet(l) * theta(l)
    qol(l) = MIN( qs1 * rhmax, qvap(l) )
    ssl    = zzl + pi_int(l+1) * theta(l) * Cp_Air 
    hol(l) = ssl + qol(l) * Hl
    hst(l) = ssl + qs1    * Hl
    dtemp  = ( pi_int(l+1) - pi_int(l) ) * theta(l)
    eta(l) = eta(l+1) + dtemp * cp_div_grav
    zzl    = zzl      + dtemp * Cp_Air    
 end do
 end if

! --- at cloud top
    qs1     = alf(ic) + bet(ic) * theta(ic)
    qol(ic) = MIN( qs1 * rhmax, qvap(ic) ) 
    ssl     = zzl + pi_int(ic1) * theta(ic) * Cp_Air 
    hol(ic) = ssl + qol(ic) * Hl
    hst(ic) = ssl + qs1     * Hl
    dtemp   = ( pi_int(ic1) - pi(ic) ) * theta(ic)
    eta(ic) = eta(ic1) + dtemp * cp_div_grav

!=====================================================================
!     ENTRAINMENT PARAMETER
!=====================================================================

    xx1 = hol(k) - hst(ic)

    xx2 = 0.0
 do l = ic,km1
    xx2 = xx2 + ( hst(ic) - hol(l) ) * ( eta(l) - eta(l+1) )
 end do

   lcase1 = ( xx2    >  0.0      ) .and.  &
            ( xx1    >  0.0      )

   lcase2 = ( xx1    <= 0.0      ) .and. &
            ( hol(k) >  hst(ic1) ) .and. &
            ( ic1    <  k ) 

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 if ( .not.lcase1 .and. .not.lcase2 )  RETURN
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  if ( lcase1 ) then
     alm = xx1 / xx2
  else 
     alm = 0.0
  end if

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  if( Tokioka_on ) alm_min = Tokioka_con / zbase
!  if( alm < alm_min ) RETURN
     xx2 = alm_min
  if( Tokioka_on ) then
     xx1  = 0.5 * ( pres_int(ic) + pres_int(ic1) )
     if(  xx1 <= Tokioka_plim ) xx2 = Tokioka_con / zbase
  endif
  if( alm < xx2 ) then
    almx = 0.0
    RETURN
  endif
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!=====================================================================
!    NORMALIZED MASSFLUX
!=====================================================================

 do l = ic,km1
    eta(l) = 1.0 + alm * eta(l)
 end do
    eta(k) = 1.0

!=====================================================================
!     CLOUD WORKFUNCTION
!=====================================================================

    wfn  = 0.0
    hccp = hol(k)

!=====================================================================
!     IN-CLOUD TEMP AND WATER MIXING RATIO
!=====================================================================

    if ( do_liq_num ) then
    thetac(k) = theta(k)
    zzl  = ( pi_int(k+1) - pi_int(k) ) * theta(k) * Cp_Air
    wlq =  qol(k)
    qt(k) = qol(k)
    endif


 if ( ic1 <= km1 ) then
 do l = km1,ic1,-1
!hcc in-cloud moist static energy (MSE)
    hcc = hccp + ( eta(l) - eta(l+1) ) * hol(l)
    dpib = pi_int(l+1) - pi(l)
    dpit = pi(l)       - pi_int(l)
!environment
    xx1  = ( eta(l+1) * dpib + eta(l) * dpit ) * hst(l)
!in-cloud
    xx2  =       hccp * dpib +    hcc * dpit
    wfn  = wfn + ( xx2 - xx1 ) * gam(l)
    hccp = hcc

    if ( do_liq_num ) then
    thetac(l) = ((hcc-zzl)/eta(l)-alf(l)*Hl)/(pi_int(l+1)*Cp_Air+bet(l)*Hl)
!    qc(l) = alf(l)+bet(l)*thetac(l)
    qc(l) = alf(l)+bet(l)*theta(l)

    ssl  = zzl + pi_int(l+1) * thetac(l) * Cp_Air
    dtemp  = ( pi_int(l+1) - pi_int(l) ) * thetac(l)

    zzl    = zzl      + dtemp * Cp_Air    
    xx1 = eta(l) - eta(l+1)
    wlq = wlq + xx1 *  qol(l)
    qt(l) = wlq/eta(l)
    endif

 end do
 end if

 if ( lcase1 ) then
    wfn = wfn + gam(ic) * ( pi_int(ic1) - pi(ic) ) * &
               ( hccp - hst(ic) * eta(ic1) )
 end if

    if ( do_liq_num ) &
      up_conv=0.7*max(wfn,0.)**0.5

!=====================================================================
!    CRITICAL CLOUD WORK FUNCTION
!=====================================================================

      xx1  = 0.5 * ( pres_int(ic) + pres_int(ic1) )
      xx2  = pres_int(k)
      ftop = 1.0

 if ( lcase2 ) then
   if ( hst(ic1) < hst(ic) ) then
      ftop = ( hst(ic1) - hol(k) ) / ( hst(ic1) - hst(ic) )
   else
      ftop = 0.0
   endif
      xx3 = 0.5 * ( pres_int(ic1) + pres_int(ic1+1) )
      xx1 = xx3 * ( 1.0 - ftop ) + xx1 * ftop
 endif

! $$$  CALL RAS_ACRITN ( xx1, xx2, wfn_crit )

        iwk = xx1 * 0.02E-2 - 0.999999999
 if ( ( iwk > 1 ) .and. ( iwk <= 15 ) ) then
         wfn_crit = ac(iwk) + xx1 * ad(iwk)
 else if ( iwk <= 1 ) then
         wfn_crit = actop
 else if ( iwk > 15 ) then
         wfn_crit = a(15)
 end if
         wfn_crit = aratio * wfn_crit * ( xx2 - xx1 )

  wfn = wfn - wfn_crit

 lcase1 = lcase1 .and. ( wfn > 0.0 )
 lcase2 = lcase2 .and. ( wfn > 0.0 ) .and. ( ftop > 0.0 )

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  if ( .not.lcase1 .and. .not.lcase2 )  then
    almx = 0.0
    RETURN
  endif
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  almx = alm
!=====================================================================

  if ( lcase1 ) then
      dhic = hst(ic) - hol(ic)
  else
      dhic = ( hol(k)  - hol(ic1) ) &
           - ( hol(ic) - hol(ic1) ) * ftop
  end if

  if ( lcase2 ) then
     xx1 = ftop * ( hol(ic) - hol(ic1) ) + hol(ic1)
     xx2 = ftop * ( qol(ic) - qol(ic1) ) + qol(ic1)
     sic = xx1 - xx2 * Hl
     qs1 = xx2 + dhic * ( 1.0 / Hl )
  end if

     hic = hol(ic)
     hkb = hol(k)

!=====================================================================

     wlq =  qol(k) - qs1      * eta(ic)
     uht = uwnd(k) - uwnd(ic) * eta(ic)
     vht = vwnd(k) - vwnd(ic) * eta(ic)

 do l = km1,ic,-1
     xx1 = eta(l) - eta(l+1)
     wlq = wlq + xx1 *  qol(l)
     uht = uht + xx1 * uwnd(l)
     vht = vht + xx1 * vwnd(l)
 end do

!=====================================================================
!         CALCULATE TRACER UPDRAFT PROPERTIES
!            AND CONVECTIVE SOURCE OF TRACER
!=====================================================================


!------------ Prognostic tracers

if ( LRcu ) then

!     wlR = tracer(k,:)
     
! do l = km1,ic,-1
!     xx1 = eta(l) - eta(l+1)
!     wlR = wlR + xx1 * tracer(l,:)
! end do
! 
!!do cloud base level
!     xx1     = 0.5 * cp_by_dp(k) / Cp_Air / onebg
!     dtracercu(k,:) = ( tracer(km1,:) - tracer(k,:) ) * xx1
     
! if ( ic1 <= km1 ) then
! do l = km1,ic1,-1
!     xx1     = 0.5 * cp_by_dp(l) / Cp_Air / onebg
!     dtracercu(l,:) = ( eta(l+1) * ( tracer(l  ,:) - tracer(l+1,:) ) + &
!                   eta(l  ) * ( tracer(l-1,:) - tracer(l  ,:) ) ) * xx1
! end do
! end if
!
! !do cloud top level
!     xx1      = cp_by_dp(ic) / Cp_Air / onebg
!     xx2      = 0.5 * xx1
!     dtracercu(ic,:) = ( eta(ic1) * ( tracer(ic,:) - tracer(ic1,:) ) * xx2 ) + &
!                  ( wlR      - eta(ic) * tracer(ic,:)  ) * xx1
!
 end if

!------------ Cloud liquid, ice, and fraction

 if ( Ldacu ) then

     wll = ql(k)
     wli = qi(k)
     if ( do_liq_num ) &
       wlN = qn(k)
 do tr=1,num_ras_tracers
     wlR(tr) = tracer(k,tr)
 enddo

 do l = km1,ic,-1
     xx1 = eta(l) - eta(l+1)
     wll = wll + xx1 * ql(l)
     wli = wli + xx1 * qi(l)
     do tr=1,num_ras_tracers
       wlR(tr) = wlR(tr) + xx1 * tracer(l,tr)
     enddo
 end do


 if ( do_liq_num ) then
   do l = km1,ic,-1
     xx1 = eta(l) - eta(l+1)
       wlN = wlN + xx1 * qn(l)
   end do
 endif
 
!do cloud base level
     xx1     = 0.5 * cp_by_dp(k) / Cp_Air / onebg
     Dlcu(k) = ( ql(km1) - ql(k) ) * xx1
     Dicu(k) = ( qi(km1) - qi(k) ) * xx1
     Dacu(k) = ( qa(km1) - qa(k) ) * xx1
     if ( do_liq_num ) &
       Dncu(k) = ( qn(km1) - qn(k) ) * xx1
     mccu(k) = eta(k)
     do tr = 1,num_ras_tracers
       dtracercu(k,tr) = ( tracer(km1,tr) - tracer(k,tr) ) * xx1
     enddo

 if ( ic1 <= km1 ) then
   do l = km1,ic1,-1
     xx1     = 0.5 * cp_by_dp(l) / Cp_Air / onebg
     Dlcu(l) = ( eta(l+1) * ( ql(l  ) - ql(l+1) ) + &
                 eta(l  ) * ( ql(l-1) - ql(l  ) ) ) * xx1
     Dicu(l) = ( eta(l+1) * ( qi(l  ) - qi(l+1) ) + &
                 eta(l  ) * ( qi(l-1) - qi(l  ) ) ) * xx1
     Dacu(l) = ( eta(l+1) * ( qa(l  ) - qa(l+1) ) + &
                 eta(l  ) * ( qa(l-1) - qa(l  ) ) ) * xx1
     mccu(l) = eta(l)
     do tr = 1,num_ras_tracers
       dtracercu(l,tr) = ( eta(l+1) * ( tracer(l  ,tr) - tracer(l+1,tr) ) + &
                     eta(l  ) * ( tracer(l-1,tr) - tracer(l  ,tr) ) ) * xx1
     enddo
   end do

   if ( do_liq_num ) then
     do l = km1,ic1,-1
       xx1     = 0.5 * cp_by_dp(l) / Cp_Air / onebg
       Dncu(l) = ( eta(l+1) * ( qn(l  ) - qn(l+1) ) + &
                   eta(l  ) * ( qn(l-1) - qn(l  ) ) ) * xx1
     enddo
   endif
 end if
 

 !do cloud top level
     xx1      = cp_by_dp(ic) / Cp_Air / onebg
     xx2      = 0.5 * xx1
     Dlcu(ic) = ( eta(ic1) * ( ql(ic) - ql(ic1) ) * xx2 ) + &
                ( wll       - eta(ic) * ql(ic)  ) * xx1
     Dicu(ic) = ( eta(ic1) * ( qi(ic) - qi(ic1) ) * xx2 ) + &
                ( wli       - eta(ic) * qi(ic)  ) * xx1
     Dacu(ic) = ( eta(ic1) * ( qa(ic) - qa(ic1) ) * xx2 ) + &
                ( eta(ic)   - eta(ic) * qa(ic)  ) * xx1
     if ( do_liq_num ) &
       Dncu(ic) = ( eta(ic1) * ( qn(ic) - qn(ic1) ) * xx2 ) + &
                  ( wlN       - eta(ic) * qn(ic)  ) * xx1
     do tr = 1,num_ras_tracers
       dtracercu(ic,tr) = ( eta(ic1) * ( tracer(ic,tr) - tracer(ic1,tr) ) * xx2 ) + &
                    ( wlR(tr)      - eta(ic) * tracer(ic,tr)  ) * xx1
     enddo

 end if

 if ( Lmccu .and. .not.Ldacu ) then

 !do cloud base level
     mccu(k) = eta(k)

 if ( ic1 <= km1 ) then
 do l = km1,ic1,-1
   mccu(l) = eta(l)
 end do
 end if

 end if

!=======================================================================
!     CALCULATE GS AND PART OF AKM (THAT REQUIRES ETA)
!=======================================================================

     xx1      = ( theta(km1) - theta(k) ) / ( pi(k) - pi(km1) )
     hol(k)   = xx1 * ( pi_int(k) - pi(km1) ) * pi(k)   * cp_by_dp(k)
     hol(km1) = xx1 * ( pi(k) - pi_int(k)   ) * pi(km1) * cp_by_dp(km1)
     akm      = 0.0

 if ( ic1 <= km1 ) then
 do l = km1,ic1,-1
     xx1      = ( theta(l-1) - theta(l) ) * eta(l) / ( pi(l) - pi(l-1) )
     hol(l)   = xx1 * ( pi_int(l) - pi(l-1) ) * pi(l)   * cp_by_dp(l) + hol(l)  
     hol(l-1) = xx1 * ( pi(l) - pi_int(l)   ) * pi(l-1) * cp_by_dp(l-1)
     akm      = akm - hol(l)  * ( eta(l)   * ( pi(l  ) - pi_int(l) ) +  &
                                  eta(l+1) * ( pi_int(l+1) - pi(l) ) ) / pi(l)
 end do
 end if

!=======================================================================
!  PARTION CLOUD LIQUID WATER INTO PRECIP & DETRAINED WATER 
!=======================================================================

     xx1 = 0.5 * ( pres_int(ic) + pres_int(ic1) )

! $$$ CALL RAS_RNCL( xx1, rn_frac)
       rn_frac = rn_frac_top
   if ( ( xx1 >= rn_ptop ) .and. ( xx1 <= rn_pbot ) ) then
       rn_frac = ( rn_pbot - xx1 ) * rn_pfac +  rn_frac_bot
   end if
   if ( xx1 > rn_pbot ) then
       rn_frac = rn_frac_bot
   end if

       wdet = ( 1.0 - rn_frac ) * wlq
       wlq  =         rn_frac   * wlq

  if ( Ldacu ) then !detrain non-precipitated fraction of condensed vapor
    
      if (coldT) then
         Dicu(ic) = Dicu(ic) + wdet * cp_by_dp(ic)/Cp_Air/onebg
      else
         Dlcu(ic) = Dlcu(ic) + wdet * cp_by_dp(ic)/Cp_Air/onebg
         if ( do_liq_num ) then
!=======================================================================
!     yim's CONVECTIVE NUCLEATION
!=======================================================================

!Use aerosol at cloud base.

!An assumption which treats ss and oc as sulfate
!convert SO4 to AS
!convert OC to OM
           totalmass(1)=aerosolmass(k,1)
           totalmass(2)=aerosolmass(k,2)
           totalmass(3)=aerosolmass(k,3)
 
           call aer_ccn_act(t(k),p(k),up_conv,totalmass,drop)
           zzl=drop*1.0e6/dens(k)

!YM choose not to do above-cloud activation for the sake of computer time
           if ( ic1 <= km1 ) then
             do l = km1,ic1,-1
               totalmass(1)=aerosolmass(l,1)
               totalmass(2)=aerosolmass(l,2)
               totalmass(3)=aerosolmass(l,3)
 
               tc=thetac(l)*pi(l)
               te=theta(l)*pi(l)
               Nc=zzl/eta(l+1)
 
               call aer_ccn_act2(t(l),p(l),up_conv,totalmass,alm,dens(l),  &
                                 Nc,qc(l),qt(l),qol(l),tc,te,drop)
                                 zzl=zzl+drop*1.0e6/dens(l)*(eta(l)-eta(l+1))
             end do
           end if



           Dncu(ic) = Dncu(ic) + zzl*(1.0-rn_frac)* &
           cp_by_dp(ic)/Cp_Air/onebg
         endif ! end of warm_cloud_aerosol interaction
      end if
          wdet = 0.0

  end if

!=======================================================================
!     CALCULATE GH
!=======================================================================

   xx1 = hol(ic)
 if ( lcase2 ) then
   xx1 = xx1 + ( sic - hic + qol(ic) * Hl ) * ( cp_by_dp(ic) / Cp_Air )
 end if

   hol(ic) = xx1 - ( wdet * Hl * cp_by_dp(ic) / Cp_Air )

 if ( lcase1 ) then
   akm = akm - eta(ic1) * ( pi_int(ic1) - pi(ic) ) * xx1 / pi(ic)
 end if

    xx2    = qol(km1) - qol(k)
    gmh(k) = hol(k) + ( xx2 * cp_by_dp(k) * Hl * 0.5 / Cp_Air )
    akm    = akm + gam(km1) * ( pi_int(k) - pi(km1) ) * gmh(k)

 if ( ic1 <= km1 ) then
 do l = km1,ic1,-1
     xx3    = xx2
     xx2    = ( qol(l-1) - qol(l) ) * eta(l)
     xx3    = xx3 + xx2
     gmh(l) = hol(l) + ( xx3 * cp_by_dp(l) * Hl * 0.5 / Cp_Air )
 end do
 end if

 if ( lcase2 ) then
    xx2 = xx2 + 2.0 * ( hkb - dhic - sic - qol(ic) * Hl ) / Hl
 end if

  gmh(ic) = xx1 + cp_by_dp(ic) * onebcp * ( xx2 * Hl * 0.5 + eta(ic) * dhic )

!=======================================================================
!     CALCULATE HC PART OF AKM
!=======================================================================

 if ( ic1 <= km1 ) then
    xx1 = gmh(k)
 do  l = km1,ic1,-1
   xx1 = xx1 + ( eta(l) - eta(l+1) ) * gmh(l)
   xx2 = gam(l-1) * ( pi_int(l) - pi(l-1) )
   if ( lcase2 .and. ( l == ic1 ) ) xx2 = 0.0
   akm = akm + xx1 * ( xx2 + gam(l) * ( pi(l) - pi_int(l) ) )
 end do
 end if

!=======================================================================

 if ( lcase2 ) then

      xx1 = 0.5*( pres_int(ic  ) + pres_int(ic1) )   &
          + 0.5*( pres_int(ic+2) - pres_int(ic ) ) * ( 1.0 - ftop )
      xx2 =       pres_int(ic1 )
      xx3 = 0.5*( pres_int(ic1 ) + pres_int(ic+2) )

 if ( ( xx1 >= xx2 ) .and. ( xx1 < xx3 ) ) then
        ftop = 1.0 - ( xx1 - xx2 ) / ( xx3 - xx2 )
         xx4 = cp_by_dp(ic1) / cp_by_dp(ic)
    hol(ic1) = hol(ic1) + hol(ic) * xx4
    gmh(ic1) = gmh(ic1) + gmh(ic) * xx4
    hol(ic)  = 0.0
    gmh(ic)  = 0.0
 else if ( xx1 < xx2 ) then
     ftop = 1.0
 else
     ftop = 0.0
 end if

 end if

!=======================================================================
!   MASS FLUX
!=======================================================================
 
 if ( ( akm < 0.0 ) .and. ( wlq >= 0.0 ) ) then
!jjs
  rasalf = rasal * ( pres_int(ic+1) - puplim ) /      &
                   (     psfc       - puplim )
  if (puplim > psfc) rasalf=0.
  rasalf = MAX( 0.0, rasalf )
!jjs
     wfn = -ftop * wfn * rasalf / akm
 else
     wfn = 0.0
 end if

!    xx1 = ( pres_int(k+1) - pres_int(k) ) * frac
     xx1 = ( psfc          - pres_int(k) ) * frac
     wfn = MIN( wfn, xx1 )

!=======================================================================
!     FOR SAK CLOUDS
!=======================================================================

 if ( Ldacu ) then
   xx1 = wfn * onebg
   do l = ic,k
     Dacu(l) = Dacu(l) * xx1
     Dlcu(l) = Dlcu(l) * xx1
     Dicu(l) = Dicu(l) * xx1
     mccu(l) = mccu(l) * xx1
   end do
   if ( do_liq_num ) then
     xx1 = wfn * onebg
     do l = ic,k
       Dncu(l) = Dncu(l) * xx1
     enddo
   endif           

 end if


 if ( LRcu ) then
      xx1 = wfn * onebg
     do l = ic,k
       do tr = 1,num_ras_tracers
          dtracercu(l,tr) = dtracercu(l,tr) * xx1
       end do
     end do
 end if

 if ( Lmccu .and. .not.Ldacu ) then
      xx1 = wfn * onebg
     do l = ic,k
          mccu(l) = mccu(l) * xx1
     end do
 end if
 
!=======================================================================
!     PRECIPITATION
!=======================================================================

   dpcu =  wlq * wfn * onebg

!=======================================================================
!     THETA AND Q CHANGE DUE TO CLOUD TYPE IC
!=======================================================================
 
    xx1 = wfn * onebcp
    xx2 = wfn * ( 1.0 / Hl )

 do l = ic,k
   dtcu(l) = xx1 * hol(l) / pi(l)
   dqcu(l) = xx2 * ( gmh(l) - hol(l) )
 end do

!=======================================================================
!  CUMULUS FRICTION 
!=======================================================================
 if (cufric) then

    xx1 = 0.5 * wfn * onebcp

! --- At cloud base
    xx2    = xx1 * cp_by_dp(k)
    dut    = ( uwnd(km1) - uwnd(k) )
    dvt    = ( vwnd(km1) - vwnd(k) )
   ducu(k) = dut * xx2
   dvcu(k) = dvt * xx2

! --- Between cloud base & cloud top
 do l = km1,ic1,-1
    xx2    = xx1 * cp_by_dp(l)
    dub    = dut
    dvb    = dvt
    dut    = ( uwnd(l-1) - uwnd(l) ) * eta(l)
    dvt    = ( vwnd(l-1) - vwnd(l) ) * eta(l)
   ducu(l) = ( dut + dub ) * xx2
   dvcu(l) = ( dvt + dvb ) * xx2
 end do

! --- At cloud top
    xx2      =   xx1 * cp_by_dp(ic)
    ducu(ic) = ( dut + uht + uht ) * xx2
    dvcu(ic) = ( dvt + vht + vht ) * xx2

  end if

!=======================================================================
  end SUBROUTINE RAS_CLOUD

!#####################################################################
!#####################################################################

  SUBROUTINE COMP_LCL( t_parc, q_parc, p_parc, pres, k_lcl )

!=======================================================================
! ***** COMPUTE LCL ( CLOUD BASE )
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!       t_parc   Initial parcel temperature
!       q_parc   Initial parcel mixing ratio
!       p_parc   Initial parcel pressure
!       pres     Pressure in colunm
!---------------------------------------------------------------------
  real, intent(in), dimension(:,:)   :: t_parc, q_parc, p_parc
  real, intent(in), dimension(:,:,:) :: pres

!---------------------------------------------------------------------
! Arguments (Intent out)
!       k_lcl   Index of LCL in column
!---------------------------------------------------------------------
  integer, intent(out), dimension(:,:) :: k_lcl

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

  real, dimension(size(t_parc,1),size(t_parc,2)) :: &
        qsat, rhum, chi, p_lcl

  integer :: k, kmax, k_lcl_min

!=====================================================================

! --- Index of lowest level
  kmax = size( pres, 3 )

! --- Compute relative humidity
  call compute_qs (t_parc, p_parc, qsat, q = q_parc)
  rhum(:,:) = q_parc(:,:) / qsat(:,:)
  rhum(:,:) = MIN( rhum(:,:), 1.0 )

! --- Compute exponent
   chi(:,:) = t_parc(:,:) / ( 1669.0 - 122.0*rhum(:,:) - t_parc(:,:) )

! --- Compute pressure at LCL
    rhum(:,:) =    chi(:,:) * LOG( rhum(:,:) )
   p_lcl(:,:) = p_parc(:,:) * EXP( rhum(:,:) )

! --- Bound p_lcl 
  p_lcl(:,:) = MAX( p_lcl(:,:), pres(:,:,1) )
  p_lcl(:,:) = MIN( p_lcl(:,:), p_parc(:,:) )

! --- Find index of LCL
  do k = 2,kmax
  where ( ( p_lcl(:,:) >= pres(:,:,k-1) ) .and. &
          ( p_lcl(:,:) <= pres(:,:,k) ) )
     k_lcl(:,:) = k
  end where
  end do

! --- Bound k_lcl
  k_lcl_min  = kmax / 2
  k_lcl(:,:) = MAX( k_lcl(:,:), k_lcl_min )

!=====================================================================
  end SUBROUTINE COMP_LCL

!#######################################################################
!#######################################################################

 SUBROUTINE RAS_CEVAP ( type,     temp,      qvap,   pres,   mass,  &
                        qvap_sat, dqvap_sat, psfc,   hl,     dtime, &
                        ksfc,     dpcu,      dtevap, dqevap, dpevap, &
                        flxprec_ib, flxprec_ib_evap )

!=======================================================================
! EVAPORATION OF CONVECTIVE SCALE PRECIP         
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!     type  - cloud type index
!     temp  - Temperature
!     qvap  - Water vapor
!     pres  - Pressure
!     mass  - Mass of layers
!     qvap_sat - Saturation value of qvap
!     dqvap_sat - Temperature derivative of qvap_sat
!     psfc  - surface presure
!     hl    - proper latent heat for the column
!     dtime - size of time step in seconds
!     ksfc  - index of lowest model level
!     dpcu - Precip in mm
!---------------------------------------------------------------------

  integer, intent(in) :: type
  real,    intent(in) :: dtime

  real,    intent(in), dimension(:) :: temp, qvap, pres, mass
  real,    intent(in), dimension(:) :: qvap_sat, dqvap_sat
  real,    intent(in)               :: psfc, hl, dpcu
  integer, intent(in)               :: ksfc

!---------------------------------------------------------------------
! Arguments (Intent out)
!     dtevap - Temperature change due to precip evap
!     dqevap - Water vapor change due to precip evap
!     dpevap - Amount of precip evaporated
!---------------------------------------------------------------------
  real, intent(out), dimension(:) :: dtevap, dqevap
  real, intent(out)               :: dpevap
  real, intent(out), dimension(:) ::  &
                flxprec_ib,    & ! precipation flux profile for cloud ib
                flxprec_ib_evap  ! evaporaton of precip profile for cloud id
                                 ! [kg/m2/s]

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

  real, parameter :: cem  = 0.054
  real, parameter :: ceta = -544.0E-6

  real, dimension(size(temp,1)) ::  temp_new, qvap_new

  real    :: prec, def, evef
  real    :: prec_mmph, pfac, emx
  integer :: itopp1, kmax, k
  real    :: totalprecip

!=======================================================================

  kmax   = size(temp(:))
  itopp1 = type + 1
  flxprec_ib =0.0
  flxprec_ib_evap =0.0

! --- Initalize
  dpevap   = 0.0
  qvap_new = qvap
  temp_new = temp

!=======================================================================

  do k = itopp1,kmax
!-----------------------------------------

! --- Compute whats available for evaporation 
   prec = MAX(dpcu - dpevap, 0.0 )

! --- Compute precipitation efficiency factor
  prec_mmph = prec * 3600.0 / dtime
  pfac      = SQRT( pres(k) / psfc )
  emx       = SQRT( cem * cfrac * prec_mmph * pfac )   
  evef      = 1.0 - EXP( ceta * dtime * emx ) 
  def=0.

! --- Evaporate precip where needed
  if ( ( hcevap*qvap_sat(k) >= qvap(k) ) .and.  &
          ( prec            > 0.0      ) .and.  &
          ( ksfc            > k        ) ) then
            def = ( hcevap*qvap_sat(k) - qvap(k) ) /    &
                  ( 1.0 + (hl * hcevap * dqvap_sat(k) / Cp_Air ) )
            def = evef*def
            def = MIN( def, prec/mass(k) )
    qvap_new(k) = qvap(k) + def
    temp_new(k) = temp(k) - (def * hl/Cp_Air)
         dpevap = dpevap + def * mass(k)
    flxprec_ib_evap(k) = def * mass(k)/dtime
  end if

!-----------------------------------------
  end do   ! itopp1
!
!-------compute precipitation flux at each model layer
! by deducting the evaporation in each layer  from 
! cloud top to surface/lowest model layer
!

  totalprecip = dpcu         ! precipitation at cloud top
!
  flxprec_ib(type) = MAX(totalprecip,0.0)/dtime
  do k = itopp1,kmax
    flxprec_ib(k) = MAX((totalprecip - flxprec_ib_evap(k)*dtime),0.0)/dtime
    totalprecip =   totalprecip - flxprec_ib_evap(k)*dtime

!-----------------------------------------
  end do


! --- Changes to temperature and water vapor from evaporation
  dtevap(:) = temp_new(:) - temp(:) 
  dqevap(:) = qvap_new(:) - qvap(:) 

!=======================================================================
  end SUBROUTINE RAS_CEVAP

!#######################################################################
!#######################################################################

 SUBROUTINE RAS_CLOUD_INDEX ( ic, ncmax )

!=======================================================================
! Set order in which clouds are to be done                   
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent out)
!       ic      Cloud order index
!       ncmax   Max number of clouds
!---------------------------------------------------------------------
 integer, intent(out), dimension(:) :: ic
 integer, intent(out) :: ncmax

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

 integer :: kmax
 integer :: km1, kcr, kfx, nc, irnd

!=======================================================================

 kmax = size( ic(:) )

 km1   = kmax  - 1
 kcr   = MIN( km1, krmax )
 kfx   = km1 - kcr
 ncmax = kfx + ncrnd

! --- Sequential clouds
  if ( kfx > 0) then
     if ( botop ) then
! --- bottom to top
       do nc = 1,kfx
          ic(nc) = kmax - nc
       end do
     else
! --- top to bottom
       do nc = kfx,1,-1
          ic(nc) = kmax - nc
       end do
     endif
  endif

! --- set non-used clouds to zero
  ic((kfx+1):kmax) = 0

! --- Random clouds
 if ( ncrnd > 0 ) then
     do nc = 1,ncrnd
       irnd      = ( RAN0(iseed) - 0.0005 ) * ( kcr - krmin + 1 )
       ic(kfx+nc) = irnd + krmin
     end do
  endif

!=======================================================================
  end SUBROUTINE RAS_CLOUD_INDEX 

!#####################################################################
!#####################################################################


 SUBROUTINE RAS_CLOUD_EXIST( k,      qvap,   qsat, theta,    &
                             pifull, pihalf, nc,   ncmax,    &
                             Hl,     exist  )
!=======================================================================
 implicit none
!=======================================================================

 integer, intent(in)               :: k, ncmax
 real,    intent(in)               :: Hl
 integer, intent(in), dimension(:) :: nc 
 real,    intent(in), dimension(:) :: qvap, qsat, theta, pifull, pihalf 

 logical, intent(out) :: exist

!---------------------------------------------------------------------
 real, parameter :: rhmax  = 0.9999 

 real, dimension(size(theta(:))) :: ssl, hst
 real                         :: zzl, hol_k, qol_k
 integer                      :: km1, ic, ic1, l

!=======================================================================

 ic  = MINVAL( nc(1:ncmax) )

 km1 = k  - 1
 ic1 = ic + 1

! --- at cloud base
    qol_k = MIN( qsat(k) * rhmax, qvap(k) )
    ssl(k) = pihalf(k+1) * theta(k) * Cp_Air 
    hol_k  = ssl(k) +  qol_k  * Hl
    hst(k) = ssl(k) + qsat(k) * Hl
    zzl    = ( pihalf(k+1) - pihalf(k) ) * theta(k) * Cp_Air

! --- between cloud base & cloud top
 do l = km1,ic1,-1
    ssl(l) = zzl + pihalf(l+1) * theta(l) * Cp_Air 
    hst(l) = ssl(l) + qsat(l) * Hl
    zzl    = zzl + ( pihalf(l+1) - pihalf(l) ) * theta(l) * Cp_Air    
 end do

! --- at cloud top
    ssl(ic) = zzl  + pihalf(ic1) * theta(ic) * Cp_Air 
    hst(ic) = ssl(ic) + qsat(ic) * Hl

! --- test
    exist = hol_k > MINVAL( hst(ic:k) )

!=======================================================================
 end SUBROUTINE RAS_CLOUD_EXIST

!#####################################################################
!#####################################################################

  SUBROUTINE RAS_BDGT ( precip, coldT, dtemp, dqvap, duwnd, dvwnd, &
                        pres_int, dql, dqi)

!=====================================================================
! ***** BUDGET CHECK FOR RAS - A DEBUGGING TOOL
!=====================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!       precip   - either rain or snow
!       coldT    - is the precip snow?
!       dtemp    - Temperature change 
!       dqvap    - Water vapor change 
!       duwnd    - U wind change 
!       dvwnd    - V wind change 
!       pres_int - Pressure at layer interface
!       dql      - OPTIONAL, liquid change
!       dqi      - OPTIONAL, ice change
!---------------------------------------------------------------------
  real, intent(in), dimension(:,:,:) :: dtemp, dqvap, duwnd, dvwnd, pres_int
  real, intent(in), dimension(:,:)   :: precip
  logical, intent(in), dimension(:,:):: coldT
  real, intent(in), optional, dimension(:,:,:) :: dql, dqi
!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

 integer :: imax, jmax, kmax, i, j, k
 real    :: sum_dtemp, sum_dqvap, sum_duwnd, sum_dvwnd
 real    :: dqvap_prec, dqvap_dtemp
 real , dimension(size(dtemp,1),size(dtemp,2),size(dtemp,3)) :: mass
 
!=====================================================================

  imax = size ( dtemp, 1 )
  jmax = size ( dtemp, 2 )
  kmax = size ( dtemp, 3 )

  mass(:,:,1:kmax) = pres_int(:,:,2:kmax+1) - pres_int(:,:,1:kmax)
  mass = mass / Grav

  do j = 1,jmax
  do i = 1,imax

    sum_dtemp = 0.                                                          
    sum_dqvap = 0. 
    sum_duwnd = 0. 
    sum_dvwnd = 0. 

  do k = 1,kmax
    sum_dtemp = sum_dtemp + dtemp(i,j,k)*mass(i,j,k)                                   
    sum_dqvap = sum_dqvap + dqvap(i,j,k)*mass(i,j,k)  
    sum_duwnd = sum_duwnd + duwnd(i,j,k)*mass(i,j,k)  
    sum_dvwnd = sum_dvwnd + dvwnd(i,j,k)*mass(i,j,k)
  end do

    dqvap_prec  = sum_dqvap + precip(i,j)
     
    if (present(dql)) then
  do k = 1,kmax
    dqvap_prec = dqvap_prec + (dql(i,j,k)+dqi(i,j,k))*mass(i,j,k)      
  end do
    end if  

     if (coldT(i,j)) then
     dqvap_dtemp = sum_dqvap + Cp_Air*sum_dtemp/HLs
     else
     dqvap_dtemp = sum_dqvap + Cp_Air*sum_dtemp/HLv
     end if

   if ( ( abs( dqvap_prec  ) > 1.0E-4 ) .or.     &
        ( abs( dqvap_dtemp ) > 1.0E-4 ) .or.     &
        ( abs( sum_duwnd   ) > 1.0E-1 ) .or.     &      
        ( abs( sum_dvwnd   ) > 1.0E-1 ) )  then
    print *
    print *, ' RAS BUDGET CHECK AT i,j = ', i,j
    print *, ' dqvap + (Cp_Air/hl)*dtemp = ',  dqvap_dtemp                                                                     
    print *, ' dqvap + precip        = ',  dqvap_prec                                                                     
    print *, ' duwnd                 = ',  sum_duwnd                                                                    
    print *, ' dvwnd                 = ',  sum_dvwnd                                                                     
    print *, 'STOP'
!    STOP
   CALL ERROR_MESG( 'RAS_BDGT', 'RAS budget thresholds exceeded.', FATAL )
   endif

  end do
  end do

!=====================================================================
  end SUBROUTINE RAS_BDGT

!#######################################################################
!#######################################################################

FUNCTION ran0(idum)


!     $Id: ras.F90,v 19.0 2012/01/06 20:11:58 fms Exp $
!     Platform independent random number generator from
!     Numerical Recipies
!     Mark Webb July 1999
      
      REAL :: ran0
      INTEGER, INTENT (INOUT) :: idum
      
      INTEGER :: IA,IM,IQ,IR,k
      REAL :: AM

      IA=16807; IM=2147483647; AM=1.0/IM; IQ=127773; IR=2836
      
      if (idum.eq.0) then
        CALL ERROR_MESG( 'ran0', 'idum=0, ZERO seed not allowed.', FATAL )
      endif

      k=idum/IQ
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      ran0=am*idum

END FUNCTION ran0

!#######################################################################
!#######################################################################
  end MODULE RAS_MOD

MODULE DIAG_CLOUD_MOD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       DIAGNOSTIC CLOUD PREDICTION - Gordon (1992)            
!
!       1999 Feb -> 2000 July
!       Contact persons: Bill Stern (for code structure information)
!                        Tony Gordon (for cloud scheme information)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------------------------------------
!  Calculates cloud fractions diagnostically using relative humidity,
!  omega and stability 
!-------------------------------------------------------------------

use mpp_mod, only: input_nml_file
 use       fms_mod, only: error_mesg, FATAL, NOTE, file_exist,    &
                          check_nml_error, open_namelist_file,       &
                          mpp_pe, mpp_root_pe,  close_file, &
                          read_data, write_data, &
                          write_version_number, stdlog
 use     fms_io_mod, only: register_restart_field, restart_file_type, &
                           save_restart, restore_state
 use  Constants_Mod, only: Cp_Air, rdgas, rvgas, Kappa, HLv
 use time_manager_mod, only:  TIME_TYPE
 use  cloud_zonal_mod, only:  CLOUD_ZONAL_INIT, GETCLD
 use  diag_cloud_rad_mod, only:  CLOUD_TAU_DRIVER, diag_cloud_rad_INIT,&
                                 cloud_pres_thick_for_tau,  &
                                 cloud_opt_prop_tg_lw, &
                                 cloud_opt_prop_tg_sw, &
                                 cloud_optical_depths, &
                                 cloud_optical_depths2
 use  sat_vapor_pres_mod, ONLY: compute_qs
 use  shallow_conv_mod, ONLY: SHALLOW_CONV_INIT,MYLCL

!-----------------------------------------------------------------------
 implicit none
!-----------------------------------------------------------------------

 private


!--------------------- version number ----------------------------------
 character(len=128) :: version = '$Id: diag_cloud.F90,v 19.0 2012/01/06 20:05:06 fms Exp $'
 character(len=128) :: tagname = '$Name: tikal $'
 logical            :: module_is_initialized = .false.
!-----------------------------------------------------------------------

!  parmameter mxband = max number of radiative bands to be considered for some
!              cloud properties
 integer,  parameter :: mxband = 4
!  parmameter nrhc = number of critical relative humidity threshold values
!  as a function of vertical layers
 integer,  parameter :: nrhc = 3

! ****************************************************************************

!  WARNING: The parameter mxband and the namelist variable nband should have
!  the same value. Is there a way to guarantee this or can an error check
!  be added to the code to make sure that this condition is NOT violated?

! ****************************************************************************
!  The following constants are derive from other constants in constants_mod
   real,  parameter :: d622  = rdgas/rvgas, d378  = 1.0-d622
   real,  parameter :: p00 = 1000.0E2
!-----------------------------------------------------------------------

! ****  parameter used and defined in def_hi_mid_low ****

!  real,  parameter :: trshld_camt = 0.25

!      TRSHLD_CAMT - This is a cloud amounht threshold value, used in
!      conjunction with the thick cloud namelist options.  If a thick
!      cloud option is on then the level above the cloud amount max
!      for a particular cloud type is considered to be part of an
!      extended cloud if it is within this fraction of the max cloud amount.
!

!----------------- arrays for cloud predictor averaging code --------------------------

    real,    allocatable, dimension (:,:,:) :: temp_sum,qmix_sum,rhum_sum
    real,    allocatable, dimension (:,:) :: qmix_sum2
    real,    allocatable, dimension (:,:,:) :: omega_sum
    real,    allocatable, dimension (:,:,:) :: lgscldelq_sum,cnvcntq_sum
    real,    allocatable, dimension (:,:)   :: convprc_sum
    integer, allocatable, dimension (:,:)   :: nsum, nsum2

!-----------------------------------------------------------------------
! for netcdf restart
type(restart_file_type), save :: Dia_restart


!---------------------------------------------------------------------


! logical switches for initialization calls: 
! if = .false. -> has not been called
 logical :: do_cpred_init = .false.
 logical :: do_crad_init = .false.

!---------------------------------------------------------------------
! --- NAMELIST (diag_cloud_nml)
!---------------------------------------------------------------------
!     RHC -    critical humidity value (ras = 0.8 - 0.84, mca =0.7)
!              (note:  in vers >= 0.9.3 a function of 3 levels
!                i.e.,"high" , "mid", "low" - but here is more general)
!     PBOUNDS - sets pressure bounds for RHC (dimension = size(rhc) - 1
!     DO_AVERAGE - logical flag for time averaging cloud predictor variables
!     LQUADRA - logical switch for turning on quadratic relation
!             for calculating rhum clouds from rhum,
!             i.e., true for quadratice scheme, false for linear scheme
!     LRHCNV - logical switch for using rhum fields as follows:
!              if true - use rel humidities modified for presence of 
!              convective clouds (rhumcnv), otherwise use original 
!              rel humidities (rhum)
!     LOMEGA - logical switch for turning on omega correction to rhum 
!              clouds - true for omega correction, otherwise false 
!     LCNVCLD - logical switch for turning on calculation of deep convective 
!              clouds - true for deep convective clouds, otherwise false 
!     L_THEQV - logical switch for turning on calculation of shallow convective 
!              clouds - true for shallow convective clouds, otherwise false 
!     LINVERS - logical switch for turning on calculation of marine stratus 
!              clouds - true for marine stratus, otherwise false 
!     LSLINGO - logical variable = true apply Slingo marine stratus 
!                scheme, otherwise = false. 
!     LREGRSC - logical variable = true apply Tim Li marine stratus 
!                scheme, otherwise = false. Slingo & Li schemes may be
!                used in combination, but atleast one scheme must be used. 
!     LTHICK_HIGH - logical variable = true -> allow possibility of raising
!               high cloud tops one sigma level to increase their thickness
!               from 1 to nmax levels; otherwise they remain thin 
!               (1 level)
!     LTHICK_MID - logical variable = true -> allow possibility of raising
!               mid cloud tops one sigma level to increase their thickness
!               from 1 to nmax levels; otherwise they remain thin 
!               (1 level)
!     LTHICK_LOW - logical variable = true -> allow possibility of raising
!               low cloud tops one sigma level to increase their thickness
!               from 1 to nmax levels; otherwise they remain thin 
!               (1 level)
!     HIGH_LEV_CLOUD_INDEX - level above which no clouds are allowed to form 
!              (model level index)
!     NOFOG - logical switch for not allowing rhum clouds (or fog)
!             to occur beneath a certain level (low_lev_cloud_index) -> 
!              nofog = true
!             to allow clouds at the lowest model level -> nofog = false
!     LOW_LEV_CLOUD_INDEX - level below which no clouds are allowed to occur 
!               when nofog = true (model level index)
!     NBAND - max number of radiative bands to be considered for some
!              cloud properties
!     PSHALLOW - top of shallow convective layer (pressure level - n/m**2 )
!     WCUT0 - omega cutoff value for omega cloud depletion factor = 0
!     WCUT1 - omega cutoff value for omega cloud depletion factor = 1
!     t_cold     temperature defining ice-liquid cloud transition
!-----------------------------------------------------------------------

      real        :: t_cold= 263.16
 real :: &
      pshallow=.750E+05, wcut0 = .10, wcut1 = 0.0
 real, dimension(nrhc) :: rhc = (/ 0.8,0.8,0.84 /)
 real, dimension(nrhc-1) :: pbounds = (/ .400E5, .750e5 /)

 integer :: &
      high_lev_cloud_index=3, low_lev_cloud_index=16, nband=4 

 logical :: &
      do_average = .true., &
      lquadra = .true., nofog = .false.,lrhcnv = .false., & 
      lomega = .true.,lcnvcld = .true.,l_theqv = .true., & 
      linvers = .false.,lslingo = .true., lregrsc = .true., &
      lthick_high = .true.,lthick_mid = .true.,lthick_low = .true.

    NAMELIST / diag_cloud_nml /                         &
       rhc,pbounds,do_average,lquadra,lrhcnv,lomega,lcnvcld,l_theqv, & 
       linvers,lslingo,lregrsc,lthick_high,lthick_mid,lthick_low, & 
       high_lev_cloud_index, nofog, low_lev_cloud_index, nband, &
       pshallow, wcut0, wcut1, t_cold

integer :: num_pts, tot_pts

 public diag_cloud_driver, diag_cloud_init, diag_cloud_end
 public diag_cloud_driver2
 public diag_cloud_sum, diag_cloud_avg, diag_cloud_avg2, do_diag_cloud
 public diag_cloud_restart

 contains

!#############################################################################      

 SUBROUTINE DIAG_CLOUD_DRIVER (is,js, &
                    temp,qmix,rhum,omega,lgscldelq,cnvcntq,convprc, &
                    pfull,phalf,psfc,coszen,lat,time, &
                    nclds,cldtop,cldbas,cldamt,r_uv,r_nir,ab_uv,ab_nir, &
                    em_lw, conc_drop, conc_ice, size_drop, size_ice, &
    kbot)

! Arguments (intent in)

 integer, intent(in)   ::  is,js
 type(time_type), intent(in)  :: time
 real, intent(in)  :: lat(:,:)
 real, intent(in), dimension (:,:,:) ::  temp,qmix,rhum,omega
 real, intent(in), dimension (:,:,:) ::  lgscldelq,cnvcntq,pfull, phalf
 real, intent(in), dimension (:,:)   ::  convprc,psfc, coszen
 
 integer, intent(in), OPTIONAL, dimension(:,:) :: kbot


!      INPUT
!      ------

!      IS,JS    starting i,j indices from the full horizontal grid
!      IX, IY   Horizontal dimensions for global storage arrays
!      TEMP     Temperature (Deg K) at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      QMIX     Mixing Ratio at full model levels 
!                   (dimensioned IDIM x JDIM x kx)
!      RHUM     Relative humidity fraction at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      OMEGA  Pressure vertical velocity at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      LGSCLDELQ  Averaged rate of change in mix ratio due to lg scale precip 
!               at full model levels  
!               (dimensioned IDIM x JDIM x kx)
!      CNVCNTQ  Accumulated count of change in mix ratio due to conv precip 
!               at full model levels  
!               (dimensioned IDIM x JDIM x kx)
!      convprc Accumulated conv precip rate summed over all
!               full model levels (mm/day )
!               (dimensioned IDIM x JDIM)
!      PFULL    Pressure at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      PHALF    Pressure at half model levels
!                   (dimensioned IDIM x JDIM x kx+1)
!      PSFC     Surface pressure field
!                   (dimensioned IDIM x JDIM)
!      COSZEN     cosine of the zenith angle
!                   (dimensioned IDIM x JDIM)
!      TIME       time of year (time_type)
!      LAT        latitudes in radians, dimensioned by (1xJDIM)   
!      KBOT      OPTIONAL; lowest model level index array
!                   (dimensioned IDIM x JDIM)
!===================================================================
! Arguments (intent out)

 integer, intent(out), dimension(:,:,:) :: cldtop,cldbas
 integer, intent(out), dimension(:,:)  ::  nclds

   real, intent(out), dimension(:,:,:), optional :: r_uv,r_nir,ab_uv, &
                                                  ab_nir,em_lw, &
                                                  conc_drop, conc_ice, &
                                                  size_drop, size_ice
   real, intent(out), dimension(:,:,:) :: cldamt

!      OUTPUT
!      ------

!       NCLDS   number of (random overlapping) clouds in column and also
!                     the current # for clouds to be operating on
!                   (dimensioned IDIM x JDIM )
!      CLDTOP   index of cloud tops (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      CLDBAS   index of cloud bottoms (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      CLDAMT   cloud amount (fraction) (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      R_UV     fractional amount of ultraviolet radiation
!                     reflected by the clouds (at cloud levels)
!      R_NIRfractional amount of near inrared radiation
!                     reflected by the clouds (at cloud levels)
!      AB_UVfractional amount of ultraviolet radiation
!                     absorbed by the clouds (at cloud levels)
!      AB_NIRfractional amount of near inrared radiation
!                     absorbed by the clouds (at cloud levels)
!      EM_LWemissivity for the clouds (at cloud levels)

!=======================================================================
!  (Intent local)
integer, dimension(size(rhum,1),size(rhum,2)) :: &
         lhight, lhighb, lmidt, lmidb, llowt
integer,  dimension(size(rhum,1),size(rhum,2),size(rhum,3))  :: icld
real, dimension(size(rhum,1),size(rhum,2)) :: qmix_kx
real,  dimension(size(rhum,1),size(rhum,2),size(rhum,3),mxband)  :: tau
real,  dimension(size(rhum,1),size(rhum,2),size(rhum,3))  :: &
                  tempcld,delp_true, delp
integer idim, jdim, kx, lk
logical :: rad_prop, wat_prop
integer :: max_cld
integer :: i,j,k,n

!       TAU        cloud optical depth (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx x MXBAND)
!      TEMPCLD    cloud layer mean temperature (degrees Kelvin)
!                    (at cloud levels)
!      DELP_TRUE  true cloud pressure thickness of distinct cloud layers 
!                    (at cloud levels)
!       QMIX_KX     Lowest level mixing ratio 
!                   (dimensioned IDIM x JDIM)
!       ICLD          marker array of cloud types/heights (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!       LHIGHT        vertical level index upper limit for high cloud tops
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LHIGHB        vertical level index lower limit for high cloud bases
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LMIDT         vertical level index upper limit for mid cloud tops
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LMIDB         vertical level index lower limit for mid cloud bases
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LLOWT         vertical level index upper limit for low cloud tops
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       IERR        Error flag
!       IPRNT       Longitude index for sample printing
!
!  note:  LK - vertical level below which no low cloud bases can exist, is
!         calculated in routines where needed using namelist inputs 
!         low_lev_cloud_index and nofog for the cloud amount prediction code
!         but is passed to the cloud raidiative properties code as an argument.
!===================================================================

!-----------------------------------------------------------------------

      idim = size(rhum,1)
      jdim = size(rhum,2)
      kx = size(rhum,3)
      tau = 0.0

      rad_prop = .false.
      wat_prop = .false.
      if (present (r_uv) .or. present(r_nir) .or. present (ab_uv) &
       .or. present(ab_nir) .or. present (em_lw) ) then
   rad_prop = .true.
   r_uv = 0.
   r_nir = 0.
   ab_uv = 0.
   ab_nir = 0.
   em_lw = 0.
endif
if (present(conc_drop) .or. present(conc_ice) .or.    &
     present(size_drop) .or. present(size_ice) ) then
   wat_prop = .true.
endif
if ( (.not. rad_prop) .and. (.not. wat_prop) ) then
  rad_prop = .true.
endif

  !  define lower limit for low cloud bases
if (nofog) then
  lk = low_lev_cloud_index
else
  lk = kx
endif

!  cldtim drives cloud prediction scheme
      call cldtim ( temp,qmix,rhum,omega,lgscldelq,cnvcntq,convprc, &
                    pfull, phalf,psfc, lat, time, tempcld,delp_true, &
                    cldtop,cldbas,cldamt,                              &
                    lhight,lhighb, lmidt, lmidb, llowt,icld,nclds,kbot)

! lowest level mixing ratio for anomalous absorption in cloud-radiation
      qmix_kx(:,:) = qmix(:,:,kx)

    max_cld  = MAXVAL(nclds(:,:))

    IF (max_cld .gt. 0) then

       call cloud_pres_thick_for_tau (nclds,icld,cldtop,cldbas, &
   &          delp_true,lhight,lhighb, lmidt, lmidb, llowt,lk, delp, &
   &          phalf, psfc )

       call cloud_optical_depths(nclds,icld,cldtop,cldbas,tempcld,delp, &
                          tau,phalf )

!  cloud_tau_driver drives cloud radiative properties scheme
      if (rad_prop) then
      call cloud_tau_driver (            qmix_kx,                   &
                                          tempcld, &
                                               tau, coszen,  &
                                   r_uv=r_uv, r_nir=r_nir, ab_nir=ab_nir, ab_uv=ab_uv, &
                 em_lw=em_lw)

      endif

      endif ! (max_cld > 0)

!-----------------------------------------------------------------------
!  print output for 1-d testing
!            iprnt = 1
!            print *, ' diag_cloud sample output for iprnt,js= ', iprnt,js
!            print *,'cloud layers = ', nclds(iprnt,1)     
!            print *,'cloud amounts = ', cldamt(iprnt,1,:)     
!            print *,'cloud tops = ', cldtop(iprnt,1,:)     
!            print *,'cloud bases = ', cldbas(iprnt,1,:)     
!            print *,'cloud types = ', icld(iprnt,1,:)    
!            print *, ' clouds_rad_tg sample output '
!            print *,'reflectivity uv = ', r_uv(iprnt,1,:)     
!            print *,'reflectivity nir = ', r_nir(iprnt,1,:)     
!            print *,'absorptivity uv = ', ab_uv(iprnt,1,:)     
!            print *,'absorptivity nir = ', ab_nir(iprnt,1,:)     
!            print *,'emissivity = ', em_lw(iprnt,1,:)     
!            print *,'optical depth = ', tau(iprnt,1,:,:)     

!-----------------------------------------------------------------------
   if (present (conc_drop)) then
!-----------------------------------------------------------------------
!!!RSH
!!      NOTE:
!  THE FOLLOWING is here as an INITIAL IMPLEMENTATION to allow compil-
!  ation and model execution, and provide "reasonable ?? " values.
!  Code developed but NOT YET ADDED HERE reflects the current approach.
!  That code is available under the fez release, and will be added to
!  the repository when upgrades to the cloud-radiation modules are com-
!  pleted.
!!!RSH
! obtain drop and ice size and concentration here, consistent with the
! diag_cloud scheme.
!  As a test case, 
!  the following is a simple specification of constant concentration and
!  size in all boxes defined as cloudy, attempting to come close to
!  the prescribd values in microphys_rad.
!  assume ice cld thickness = 2.0 km; then conc_ice=10.0E-03 => 
!    iwp = 20 g/m^2, similar to that prescribed in microphys_rad.
!  assume water cld thickness = 3.5 km; then conc_drop = 20E-03 =>
!    lwp = 70 g / m^2, similar to that prescribed in microphys_rad.
!   use sizes as used in microphys_rad (50 and 20 microns). when done,
!   radiative boundary fluxes are "similar" to non-microphysical results
!   for test case done here, and shows reasonable sensitivity to
!   variations in concentrations.

    conc_ice = 0.
    conc_drop = 0.
    size_ice = 50.
    size_drop = 20.


    IF (max_cld .gt. 0) then
         idim = size(tau,1)
           jdim = size(tau,2)
         do j= 1,jdim
           do i=1,idim
          do n=1,nclds(i,j)
          do k=cldtop(i,j,n), cldbas(i,j,n)
           if (tempcld(i,j,n) < t_cold) then
            conc_ice(i,j,k) = 10.0E-03  ! units : g/m^3
          size_ice(i,j,k) = 50.       ! units : diameter in microns
           else
            conc_drop(i,j,k) = 20.0E-03 ! units : g/m^3
            size_drop(i,j,k) = 20.      ! units : diameter in microns
           endif
           end do
        end do
           end do
           end do

    endif
 endif


end SUBROUTINE DIAG_CLOUD_DRIVER

!---------------------------------------------------------------------

subroutine diag_cloud_driver2 (is, js, press, pflux, lat, time, nclds, &
                               cldtop, cldbas, cldamt, liq_frac, tau, &
                               ice_cloud, kbot) 

!--------------------------------------------------------------------- 
!    diag_cloud_driver2 returns the cloud specification arrays for the 
!    gordon diag cloud scheme. returned are the number of clouds per 
!    column, the cloud top, cloud base and fractional coverage of each 
!    cloud, the amount of the cloud which is liquid, its optical depth 
!    and an indicator as to whether it is ice or liquid (a different 
!    criterion than is used for the liquid fraction determination).
!----------------------------------------------------------------------
 
integer,                     intent(in)             ::  is,js
real,    dimension (:,:,:),  intent(in)             ::  press, pflux 
real,    dimension(:,:),     intent(in)             ::  lat
type(time_type),             intent(in)             ::  time
integer, dimension(:,:),     intent(inout)          ::  nclds
integer, dimension(:,:,:),   intent(out)            ::  cldtop,cldbas
real,    dimension(:,:,:),   intent(out)            ::  cldamt, liq_frac
real,    dimension(:,:,:,:), intent(out)            ::  tau
logical, dimension(:,:,:),   intent(out)            ::  ice_cloud
integer, dimension(:,:),     intent(in), optional   ::  kbot

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,js             starting subdomain i,j indices of data 
!                        in the physics_window being integrated
!      press             pressure at model levels (1:nlev), surface    
!                        pressure is stored at index value nlev+1   
!                        [ (kg /( m s^2) ]
!      pflux             average of pressure at adjacent model levels  
!                        [ (kg /( m s^2) ]
!      lat               latitude of model points  [ radians ]
!      time              time at which radiation calculation is to apply
!                        [ time_type (days, seconds) ]
!
!   intent(inout) variables:
!
!      nclds             total number of clouds in each grid column
!
!   intent(out) variables:
!
!      cldtop            k index of cloud top for each cloud
!      cldbas            k index of cloud base for each cloud
!      cldamt            fractional cloudiness for each cloud
!                        [ dimensionless ]
!      liq_frac          fraction of cloud which is liquid 
!                        [ dimensionless ]
!      tau               cloud optical depth  [ dimensionless ]
!      ice_cloud         logical flag indicating whether cloud is liquid
!                        or ice
!
!    intent(in), optional variables:
!
!      kbot              present when running eta vertical coordinate,
!                        index of lowest model level above ground
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      integer, dimension (size(press,1), size(press,2)) :: &
                                    lhight, lhighb, lmidt, lmidb, llowt

      integer,  dimension (size(press,1), size(press,2),  &
                                          size(press,3)-1)  :: icld

      real,  dimension (size(press,1), size(press,2),   &
                                       size(press,3)-1)  :: &
                         temp, qmix, rhum, omega, lgscldelq, cnvcntq, &
                 delp, tempcld, delp_true, pfull

      real, dimension (size(press,1), size(press,2)) ::  convprc, psfc

      real,  dimension(size(press,1), size(press,2),  &
                       size(press,3))  ::  phalf

      integer     :: kx, lk, ierr, max_cld
      integer     :: i, j, n
      
!---------------------------------------------------------------------
!   local variables:
!
!       lhight     vertical level index upper limit for high cloud tops
!                  a function of lat and lon
!       lhighb     vertical level index lower limit for high cloud bases
!                  a function of lat and lon
!       lmidt      vertical level index upper limit for mid cloud tops
!                  a function of lat and lon
!       lmidb      vertical level index lower limit for mid cloud bases
!                  a function of lat and lon
!       llowt      vertical level index upper limit for low cloud tops
!                  a function of ltt and lon
!       icld       marker array of cloud types/heights (at cloud levels)
!       temp       temperature at full model levels [ deg K ]
!       qmix       mixing ratio at full model levels [ kg H2O / kg air ]
!       rhum       relative humidity fraction at full model levels
!                  [ dimensionless ]
!       omega      pressure vertical velocity at full model levels
!                  [ mb / sec ??????? ]
!       lgscldelq  averaged rate of change in mixing ratio due to large 
!                  scale precip at full model levels  
!       cnvcntq    accumulated count of change in mixing ratio due to 
!                  convective  precip at full model levels  
!       delp       pressure thickness of model layers [ kg / (m s^2) ]
!       tempcld    cloud layer mean temperature, at cloud levels
!                  [ deg K ]
!       delp_true  true cloud pressure thickness of distinct cloud 
!                  layers (at cloud levels) [ kg / (m s^2) ]
!       pfull      pressure at full levels [ kg / (m s^2) ]
!       convprc    accumulated conv precip rate summed over all
!                  full model levels [ mm/day ]
!       psfc       surface pressure field [ kg / (m s^2) ]
!       phalf      pressure at model half levels [ kg / (m s^2) ]


!       kx         number of model layers
!       lk         vertical level below which no low cloud bases can 
!                  exist
!       ierr       error flag
!       max_cld    max number of clouds in any column in the current
!                  physics window
!       i,j,n    do loop indices
!
!--------------------------------------------------------------------- 

!----------------------------------------------------------------------
!    define the number of model layers.
!----------------------------------------------------------------------
      kx   = size(press,3) - 1

!---------------------------------------------------------------------
!    define the needed pressure arrays. 
!---------------------------------------------------------------------
      pfull(:,:,:) = press(:,:,1:kx)
      phalf(:,:,:) = pflux(:,:,:)
      psfc(:,:)    = press(:,:,kx+1)

!--------------------------------------------------------------------
!    call diag_cloud_avg to obtain the appropriate values for the input 
!    arrays needed to define the cloud locations and amounts. these may
!    or may not be time-averaged values.
!---------------------------------------------------------------------
      call diag_cloud_avg (is, js, temp, qmix, rhum, omega, lgscldelq, &
                           cnvcntq, convprc, ierr)

!----------------------------------------------------------------------
!    initialize the output fields produced by this module.
!---------------------------------------------------------------------
      tau = 0.
      liq_frac = 0.
      cldamt = 0.
      cldtop = 0
      cldbas = 0
      ice_cloud = .false.

!----------------------------------------------------------------------
!    if input data was appropriately returned from diag_cloud_avg,
!    proceed with the determination of the cloud field.
!---------------------------------------------------------------------
      if (ierr == 0) then

!---------------------------------------------------------------------
!    define the lowest model level which can be a cloud base. it is 
!    either the lowest model level, or a level determined from namelist
!    input.
!---------------------------------------------------------------------
        if (nofog) then
          lk = low_lev_cloud_index
        else
          lk = kx
        endif

!--------------------------------------------------------------------
!    call cldtim to drive the cloud prediction scheme.
!--------------------------------------------------------------------
        call cldtim (temp, qmix, rhum, omega, lgscldelq, cnvcntq,   &
                     convprc, pfull, phalf, psfc, lat, time, tempcld, &
                     delp_true, cldtop, cldbas, cldamt, lhight, lhighb,&
                     lmidt, lmidb, llowt, icld, nclds, kbot)

!---------------------------------------------------------------------
!    determine the maximum number of clouds in any of the columns in the
!    physics window.
!---------------------------------------------------------------------
        max_cld  = MAXVAL(nclds(:,:))

!--------------------------------------------------------------------
!    if cloud is present anywhere in the window, call 
!    cloud_pres_thick_for_tau to determine the cloud thicknesses and
!    the call cloud_optical_depths2 to determine the optical depths and
!    liquid fraction of each cloud.
!---------------------------------------------------------------------
        if (max_cld > 0) then
          call cloud_pres_thick_for_tau (nclds, icld, cldtop, cldbas, &
                                         delp_true, lhight, lhighb,  &
                                         lmidt, lmidb, llowt, lk, delp,&
                                         phalf, psfc)
          call cloud_optical_depths2 (nclds, icld, cldtop, cldbas,  & 
                                      tempcld, delp, tau, phalf, &
                                      liq_frac)
        endif

!---------------------------------------------------------------------
!    determine whether the cloud temperature will support a liquid or
!    an ice cloud. the parameter t_cold is the cutoff temperature value.
!---------------------------------------------------------------------
        do j= 1,size(press,2)
          do i=1,size(press,1)
            do n=1,nclds(i,j)
              if (tempcld(i,j,n) < t_cold) then
                ice_cloud(i,j,n) = .true.
              endif
            end do
          end do
        end do

!----------------------------------------------------------------------
!    if input data was not acceptably returned from diag_cloud_avg,
!    determine if this represents an error, or is just a manifestation
!    of coldstart behavior. if this is coldstart step, set clouds to
!    zero and continue.
!---------------------------------------------------------------------
      else
        if (num_pts >= tot_pts) then
          call error_mesg ('diag_cloud_mod',  &
             ' no diag cloud data available; ierr /= 0', FATAL)
        else
          num_pts = num_pts + size(press,1)*size(press,2)
          nclds = 0
        endif
      endif ! (ierr=0)


!--------------------------------------------------------------------




end subroutine diag_cloud_driver2



!---------------------------------------------------------------------



!##################################################################      

 SUBROUTINE CLDTIM (temp,qmix,rhum, omega,lgscldelq,cnvcntq,convprc,  &
                    pfull, phalf,psfc, lat, time, tempcld,delp_true, &
                    cldtop,cldbas,cldamt,                              &
                    lhight,lhighb, lmidt, lmidb, llowt,icld,nclds,kbot)

!===================================================================
! Arguments (intent in)

! Namelist variables used in this routine (defined at top of module)
!     HIGH_LEV_CLOUD_INDEX - level above which no clouds are allowed to form 
!              (model level index)
!     NOFOG - logical switch for not allowing rhum clouds (or fog)
!             to occur beneath a certain level (low_lev_cloud_index) -> 
!              nofog = true
!             to allow clouds at the lowest model level -> nofog = false
!     LOW_LEV_CLOUD_INDEX - level below which no clouds are allowed to occur 
!               when nofog = true (model level index)


 type(time_type), intent(in)  :: time
 real, intent(in)  :: lat(:,:)
 real, intent(in), dimension (:,:,:) ::  temp,qmix,rhum,omega
 real, intent(in), dimension (:,:,:) ::  lgscldelq,cnvcntq,pfull, phalf
 real, intent(in), dimension (:,:) ::    convprc,psfc

!
!      INPUT
!      -----
!
!      TEMP     Temperature (Deg K) at full model levels 
!                   (dimensioned IDIM x JDIM x kx)
!      QMIX     Mixing Ratio at full model levels 
!                   (dimensioned IDIM x JDIM x kx)
!      RHUM     Relative humidity fraction at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      OMEGA  Pressure vertical velocity at full model levels 
!                   (dimensioned IDIM x JDIM x kx)
!      LGSCLDELQ  Averaged rate of change in mix ratio due to lg scale precip 
!               at full model levels  
!               (dimensioned IDIM x JDIM x kx)
!      CNVCNTQ  Accumulated count of change in mix ratio due to conv precip 
!               at full model levels  
!               (dimensioned IDIM x JDIM x kx)
!      CONVPRC Accumulated conv precip rate summed over all
!               full model levels (mm/day )
!               (dimensioned IDIM x JDIM)
!      PFULL    Pressure at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      PHALF    Pressure at half model levels
!                   (dimensioned IDIM x JDIM x kx+1)
!      PSFC     Surface pressure field
!                   (dimensioned IDIM x JDIM)
!      TIME       time of year (time_type)
!      LAT        latitudes in radians, dimensioned by (IDIMxJDIM)   
!      KBOT    -  OPTIONAL; lowest model level index array
!                   (dimensioned IDIM x JDIM)
!===================================================================
! Arguments (intent out)

integer, intent(out), dimension(:,:,:) :: cldtop,cldbas,icld
integer, intent(out), dimension(:,:)  :: nclds
integer, intent(out), dimension(:,:)  :: lhight,lhighb, lmidt, lmidb, llowt
   real, intent(out), dimension(:,:,:) :: cldamt
   real, intent(out), dimension(:,:,:) :: tempcld,delp_true


integer, intent(in), OPTIONAL, dimension(:,:) :: kbot


!      OUTPUT
!      ------

!      ICLD          marker array of cloud types (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      CLDTOP     index of cloud tops (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      CLDBAS     index of cloud bottoms (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      CLDAMT   cloud amount (fraction) (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      TEMPCLD    cloud layer mean temperature (degrees Kelvin)
!                    (at cloud levels)
!      DELP_TRUE  true cloud pressure thickness of distinct cloud layers 
!                    (at cloud levels)
!      R_UV    fractional amount of ultraviolet radiation
!                     reflected by the clouds
!      R_NIRfractional amount of near inrared radiation
!                     reflected by the clouds
!      AB_UVfractional amount of ultraviolet radiation
!                     absorbed by the clouds
!      AB_NIRfractional amount of near inrared radiation
!                     absorbed by the clouds
!      EM_LWemissivity for the clouds
!       NCLDS        number of (random overlapping) clouds in column and also
!                        the current # for clouds to be operating on
!       LHIGHT        vertical level index upper limit for high cloud tops
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LHIGHB        vertical level index lower limit for high cloud bases
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LMIDT         vertical level index upper limit for mid cloud tops
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LMIDB         vertical level index lower limit for mid cloud bases
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LLOWT         vertical level index upper limit for low cloud tops
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!
!  note:  vertical level below which no low cloud bases can exist is
!         calculated in routines where needed using namelist inputs 
!         low_lev_cloud_index and nofog

!===================================================================


!=======================================================================
!  (Intent local)
 real , dimension(size(rhum,1),size(rhum,2),size(rhum,3)) :: theta
 real , dimension(size(rhum,1),size(rhum,2),size(rhum,3)) :: rhumcnv, &
      camtcnv,camtrh,camtw,camtsh,camtsc,camt
 real, dimension (size(rhum,1),size(rhum,2)) :: camtsh_mx 
 real, dimension (size(phalf,1),size(phalf,2),size(phalf,3)) :: pnorm
 integer , dimension(size(rhum,1),size(rhum,2),size(rhum,3)) :: icld_k
 integer , dimension(size(rhum,1),size(rhum,2),3) ::  kthclm, kbhclm
 integer, dimension (size(rhum,1),size(rhum,2)) :: kmaxshl,kminshl
! horizontal dimensions
 integer idim,jdim
! bottom level vertical index
 integer kx 
! loop index
 integer k 

!-----------------------------------------------------------------------
!      THETA    Potential temperature 
!               at full model levels (Deg K) )
!               (dimensioned IDIM x JDIM x kx)
!      THETA_E  Equivalent potential temperature 
!               at full model levels (Deg K) )
!               (dimensioned IDIM x JDIM x kx)
!      RHUMCNV  Relative humidity fraction modified for convective clouds
!                   (dimensioned IDIM x JDIM x kx)
!      CAMTCNV  tentative deep convective cloud amounts  
!                   (dimensioned IDIM x JDIM x kx)
!      CAMTRH   tentative rel. humidity cloud amounts  
!                   (dimensioned IDIM x JDIM x kx)
!      CAMTW    tentative cloud amounts omega corrected clouds 
!                   (dimensioned IDIM x JDIM x kx)
!      CAMTSH   tentative shallow convective cloud amounts  
!                   (dimensioned IDIM x JDIM x kx)
!      CAMTSC   tentativemarine stratus cloud amounts  
!                   (dimensioned IDIM x JDIM x kx)
!      CAMTSL   tentative Slingo marine stratus cloud amounts  
!                   (dimensioned IDIM x JDIM x kx)
!      CAMTTL   tentative Tim Li marine stratus cloud amounts  
!                   (dimensioned IDIM x JDIM x kx)
!      CAMTSH_MX maximum value of shallow convective cloud amount within
!                  shallow convective layer. (dimensioned IDIM x JDIM)
!      CAMT     tentative merged cloud amounts  
!                   (dimensioned IDIM x JDIM x kx)
!      PNORM    Normalized pressure at half or full  model levels
!                   (dimensioned IDIM x JDIM x kx+1)
!      ICLD_K     tentative merged cloud marker array  
!                   (dimensioned IDIM x JDIM x kx)
!      KTHCLM,KBHCLM gcm vert coord indices of climo high, mid, low cloud   
!                   tops and bases, (dimensioned IDIM x JDIM x 3)
!      KMINSHL,KMAXSHL indices of vertical levels corresponding to 
!                   final top and base of shallow convective cloud layer
!                   (dimensioned IDIM x JDIM)

!-----------------------------------------------------------------------
! Calculate gcm vetical coord. indices for high, mid and low clouds based
! on observed  mean climatology. They will be a function of latitude and 
! time of year.  Optionally, climatological zonal mean cloud amounts may
! also be returned, presumably for a special "cold" start situation.

! ***** caveats: the current climo cloud setup uses climo average pressure 
! *****          values for locating single layer high and middle clouds
! *****          and the bottom and top of climo average low clouds
! *****          in the future a better approach to locating a vertical
! *****          layer range for each type of cloud will be needed

! *****          camt from getcld is only used to complete the argument list

! need to normalize pressures since table is based on assumption that 
! surface pressure = 101325 Pa


      idim = size(rhum,1)
      jdim = size(rhum,2)
      kx = size(rhum,3)

! for clouds at half levels
!         do k=1,kx+1
!            pnorm(:,:,k)=101325.*phalf(:,:,k)/phalf(:,:,kx+1)
!         enddo
! call getcld (time, lat, phalf, kthclm, kbhclm, camt)

! for clouds at full levels
         do k=1,kx
            pnorm(:,:,k)=101325.*pfull(:,:,k)/phalf(:,:,kx+1)
         enddo
         pnorm(:,:,kx+1) = 101325.
! call getcld (time, lat, pfull, kthclm, kbhclm, camt)
call getcld (time, lat, pnorm, kthclm, kbhclm, camt)

!  re-initialize camt, icld_k, nclds, cldbas,cldtop
      camt(:,:,:) = 0.0
      icld_k(:,:,:) = 0
      cldbas(:,:,:) = 0
      cldtop(:,:,:) = 0
      nclds(:,:) = 0

      lhight(:,:) = high_lev_cloud_index
      lhighb(:,:) = kthclm(:,:,2) - 1
      lmidt(:,:) = kthclm(:,:,2) 
!  if climo mid cloud top level = climo high cloud base level, reset high base 
      where (kthclm(:,:,2) .eq. kbhclm(:,:,1))
         lhighb(:,:) = kbhclm(:,:,1)
         lmidt(:,:) = kbhclm(:,:,1) + 1
      endwhere
      lmidb(:,:) = kthclm(:,:,3) - 1
      llowt(:,:) = kthclm(:,:,3) 
!  if climo low cloud top level = climo mid cloud base level, reset mid base 
      where (kthclm(:,:,3) .eq. kbhclm(:,:,2))
         lmidb(:,:) = kbhclm(:,:,2)
         llowt(:,:) = kbhclm(:,:,2) + 1
      endwhere
         
!      print *, ' climatological cloud tops & bases'
!      print *, ' kthclm = ', kthclm
!      print *, ' kbhclm = ', kbhclm

 
!-----------------------------------------------------------------------
!  calculate potential temperature for use in computing marine stratus and/or 
!  shallow convective clouds

      if (linvers .or. l_theqv) then
        theta(:,:,:) = temp(:,:,:)*(p00/pfull(:,:,:)**Kappa)
      endif

!  calculate deep convective clouds
          if (lcnvcld) then
call cloud_cnv (rhum,cnvcntq,convprc, pfull, phalf, camtcnv, rhumcnv )
          endif

!  calculate rel humidity clouds
          if (lrhcnv .and. lcnvcld) then
call cloud_rhum (rhumcnv, pnorm, camtrh)
          else
call cloud_rhum (rhum, pnorm, camtrh)
          endif

!  calculate omega corrected rel humidity clouds
          if (lomega) then
call cloud_omga (camtrh, omega, llowt, camtw)
          endif


!  calculate marine stratus clouds
! ******* not implemented at this time *******
          if (linvers) then
!!!  call cloud_m_stratus (...,camtsc )
          endif

      kx = size(rhum,3)
!  calculate shallow convective clouds
          if (l_theqv) then
!  need to call shallow_conv_init to set up constants for mylcl
call cloud_shallow_conv (theta,omega,pfull,phalf,temp,qmix,camtrh, &
                    camtsh,camtsh_mx,kminshl,kmaxshl,kbot)
          endif

!  merge stratiform cloud types
call merge_strat_clouds (camtrh,camtw,camtsc,llowt,camt,icld_k)


!  group clouds into high, middle and low
call def_hi_mid_low (phalf,lhight,lhighb,lmidt,lmidb,llowt,camt,icld_k)


!  merge convective cloud types
          if (lcnvcld .or. l_theqv) then
call merge_cnv_clouds (camtcnv,camtsh,camtsh_mx,kminshl,kmaxshl, &
                        lhight,lhighb,lmidt,lmidb,llowt,camt,icld_k)
          endif

!  calculate vertical indices and tot number of distinct cloud layers
!  as a  function of lat and longitude.
call layr_top_base (icld_k,nclds,cldtop,cldbas,icld)

!  calculate vertical, mass weighted cloud amount fraction of distinct cloud
!  layers
call cldamt_mn (nclds,cldtop,cldbas,phalf,camt,cldamt)

!  compute total cloud amount from cloud amount of in distinct cloud layers
! ***** this diagnostice routine is postpooned for now *******
! call total_cld_amt()

!  Define the occurence of anvil cirrus clouds
call anvil_cirrus (lhighb,lmidb,nclds,cldtop,cldbas,cnvcntq,lgscldelq, &
                    pfull,phalf,icld)

!  Compute the mass weighted mean cloud temperature and pressure thickness
!  of each distinct cloud layer
call cld_layr_mn_temp_delp (nclds,cldtop,cldbas,temp,phalf,tempcld,delp_true)

!  Compute the cum mean water vapor mixing ratio of each distinct cloud layer
!  (This is proposed for the future)
!  call cld_lay_cum_qmix()


!-----------------------------------------------------------------------

end subroutine CLDTIM
!=======================================================================

!#############################################################################      
 
subroutine CLOUD_CNV (rhum,cnvcntq,convprc, pfull, phalf, camtcnv, rhumcnv )
                                  


!-------------------------------------------------------------------
!  This subroutine calculates deep convective cloud amounts
!-------------------------------------------------------------------

!===================================================================
! Arguments (intent in)

!-------------------------------------------------------------------
! Namelist variables used in this routine (defined at top of module)
!     HIGH_LEV_CLOUD_INDEX - level above which no clouds are allowed to form 
!              (model level index)
!     NOFOG - logical switch for not allowing rhum clouds (or fog)
!             to occur beneath a certain level (low_lev_cloud_index) -> 
!              nofog = true
!             to allow clouds at the lowest model level -> nofog = false
!     LOW_LEV_CLOUD_INDEX - level below which no clouds are allowed to occur 
!               when nofog = true (model level index)
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  parameters used in cloud_cnv

 real,  parameter :: beta = .125, cprc_thres = 0.14 
 real,  parameter :: cmax = 0.8, ctower = 0.25
 real,  parameter :: cnvdelqmn = 0.0

!      BETA -       convective cloud regression coefficients
!      CPRC_THRES - minimum convective precip amount imposed for
!                   for work array.  This allows the use of the log
!                   function in the convective cloud regression relation.
!      CMAX, CTOWER - factors limiting convective cloud amount
!                           and vertical profile
!      CNVDELQMN - min threshold value of accumulated count of mixing ratio 
!                  change due to convection above which convective clouds are 
!                  allowed
!                  ( It is set = 0.0 as a parameter, because it not anticipated
!                    that it would be changed.)
!-------------------------------------------------------------------

 real, intent(in), dimension (:,:,:) :: rhum, cnvcntq, pfull, phalf
 real, intent(in), dimension (:,:) :: convprc

!-------------------------------------------------------------------
!      RHUM     Relative humidity fraction at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      CNVCNTQ  Accumulated count of change in mix ratio due to conv precip 
!               at full model levels  
!               (dimensioned IDIM x JDIM x kx)
!      convprc Accumulated conv precip rate summed over all
!               full model levels (mm/day )
!               (dimensioned IDIM x JDIM)
!      PFULL    Pressure at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      PHALF    Pressure at half model levels
!                   (dimensioned IDIM x JDIM x kx+1)
!===================================================================
! Arguments (intent out)
 real, intent(out), dimension (:,:,:) :: camtcnv, rhumcnv
!     CAMTCNV - convective cloud amount at full model levels 
!               (dimensioned IDIM x JDIM x kx)
!     RHUMCNV - relative humidity adjusted for convective clouds
!=======================================================================
!  (Intent local)

 real, dimension (size(rhum,1),size(rhum,2),size(rhum,3)) :: delp
 real, dimension (size(rhum,1),size(rhum,2)) :: wrkcprc

!      DELP       pressure thickness of model layers 
!                   (dimensioned IDIM x JDIM x kx)
!      DCPRC - weighted convective precip as a function of vertical layer
!                   (dimensioned IDIM x JDIM x kx)
!      WRKCPRC - convective precip work array
!                   (dimensioned IDIM x JDIM)

 integer i, j, idim, jdim, kx, lk, lkxm1, kcbtop
 real alpha
!-----------------------------------------------------------------------
!  type loop index variable
 integer k 

!===================================================================

! Initialize convective cloud amount array
      camtcnv = 0.0

!  Define regression coefficient alpha base on coeficient beta, defined
!  above as a parameter
      alpha = -beta*log(cprc_thres)

!-----------------------------------------------------------------------
! <><><><><><><><>   set up vertical index range <><><><><><><><>
!-----------------------------------------------------------------------
      idim  = SIZE(rhum,1)
      jdim  = SIZE(rhum,2)
      kx = size(rhum,3)

!  define lower limit of cloud bottoms
      if (nofog) then
        lk = low_lev_cloud_index
      else
        lk = kx
      endif

!  no convective clouds allowed in lowest layer
      lkxm1 = min(lk,kx-1)  

!  define upper limit of cloud tops
      kcbtop = high_lev_cloud_index


! calculate pressure thickness of model layers, for vertcal averaging
      do k=1,kx
        delp (:,:,k) = phalf(:,:,k+1)-phalf(:,:,k)
      end do
!-----------------------------------------------------------------------
! <><><><><><><><>   code to identify convecting levels <><><><><><><><>
!-----------------------------------------------------------------------

! use total convective precip amount passed in as argument
      wrkcprc(:,:) = convprc(:,:)
!-----------------------------------------------------------------------
! <><><><><><><><>   code to compute conv cloud fraction <><><><><><><><>
!-----------------------------------------------------------------------

      do k=kcbtop,lkxm1
! check that accumulated count of convective mix ratio changes is greater
! than a minimum threshold (usually 0 )
        where (cnvcntq (:,:,k) > cnvdelqmn .and. wrkcprc(:,:) > cprc_thres)
          camtcnv(:,:,k) = alpha + beta * log(wrkcprc(:,:) )
        end where
      end do

! Impose contraints on convective cloud

      camtcnv(:,:,:) = min(camtcnv(:,:,:),cmax)


              do k=kcbtop,lkxm1
            do j=1,jdim
            do i=1,idim
      if (camtcnv(i,j,k) .lt. 0.0) then
         print *, ' pe,i,j,k,camtcnv = ', mpp_pe(),i,j,k,  &
                    camtcnv(i,j,k)
         call error_mesg ('cloud_cnv','cloud amount < 0' ,FATAL) 
      endif
            end do
            end do
              end do

      where (camtcnv(:,:,:) > 0.0 .and. pfull(:,:,:) < pshallow)
        camtcnv(:,:,:) = ctower * camtcnv(:,:,:) 
      endwhere

! calculate relative humidity adjusted for convective clouds

       rhumcnv (:,:,:) = (1.0 - camtcnv(:,:,:)) * rhum(:,:,:) 


!-----------------------------------------------------------------------

end subroutine CLOUD_CNV

!#############################################################################      

subroutine CLOUD_RHUM (rhum,pnorm, camtrh)
                                  


!-------------------------------------------------------------------
!  calculates stratiform cloud amounts based on relative humidities
!-------------------------------------------------------------------

!===================================================================
! Arguments (intent in)

! Namelist variables used in this routine (defined at top of module)
!     HIGH_LEV_CLOUD_INDEX - level above which no clouds are allowed to form 
!              (model level index)
!     NOFOG - logical switch for not allowing rhum clouds (or fog)
!             to occur beneath a certain level (low_lev_cloud_index) -> 
!              nofog = true
!             to allow clouds at the lowest model level -> nofog = false
!     LOW_LEV_CLOUD_INDEX - level below which no clouds are allowed to occur 
!               when nofog = true (model level index)
!-------------------------------------------------------------------

real, intent(in), dimension (:,:,:) :: rhum,pnorm

!-------------------------------------------------------------------
!      RHUM     Relative humidity fraction at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      PNORM    Normalized pressure at half or full  model levels
!                   (dimensioned IDIM x JDIM x kx+1)
!===================================================================

! Arguments (intent out)
 real, intent(out), dimension (:,:,:) :: camtrh
!     CAMTRH - Rel Humidity cloud amount at full model levels 
!               (dimensioned IDIM x JDIM x kx)
!=======================================================================
!  (Intent local)
 real, dimension (size(rhum,1),size(rhum,2)) :: rhc_work
! real, dimension (size(rhum,1),size(rhum,2),size(rhum,3)) :: rhum2      
 integer kx, lk, npower
!-----------------------------------------------------------------------
!  type loop index variable
 integer k,kk

!===================================================================

      kx = size(rhum,3)


!  define cloud amt - rel hum relation as linear or quadratic

      if (lquadra) then
        npower = 2
      else
        npower = 1
      endif


!    if (mpp_pe() == 0) then
!      print *, 'npower', npower
!    endif

!  define lower limit rel hum clouds
      if (nofog) then
        lk = low_lev_cloud_index 
      else
        lk = kx
      endif

! Initialize rhum cloud amount array
     camtrh = 0.0
!
!-----------------------------------------------------------------------
! <><><><><><><><>   code to compute rhum cloud fraction <><><><><><><><>
!-----------------------------------------------------------------------

      do k=high_lev_cloud_index,lk
        where (pnorm(:,:,k) .ge. pbounds(nrhc-1))
              rhc_work(:,:) = rhc(nrhc)
        end where
      do kk=nrhc-1,2
        where ((pnorm(:,:,k).lt.pbounds(kk)) .and. &
               (pnorm(:,:,k).ge.pbounds(kk-1)))
              rhc_work(:,:) = rhc(kk)
        end where
      end do
        where (pnorm(:,:,k).lt.pbounds(1)) 
              rhc_work(:,:) = rhc(1)
        end where
!BUGFIX ??
!       where (rhum(:,:,k) == rhc_work(:,:))
!    rhum2(:,:,k) = rhum(:,:,k) - 0.00007
!elsewhere
!    rhum2(:,:,k) = rhum(:,:,k)
!       endwhere
        where (rhum (:,:,k) > rhc_work(:,:))
          camtrh(:,:,k) = min(1.0, (  &
           (rhum (:,:,k) - rhc_work(:,:))/(1-rhc_work(:,:)) )** npower)
        elsewhere
          camtrh(:,:,k) = 0.0
        end where
      end do


end subroutine CLOUD_RHUM

!#############################################################################      

subroutine CLOUD_OMGA (camtrh, omega, llowt, camtw)
                                  


!-------------------------------------------------------------------
!  calculates omega corrected cloud amounts
!  This subroutine reduces rel hum clouds in regions of decending motion
!-------------------------------------------------------------------

!===================================================================
! Arguments (intent in)

! Namelist variables used in this routine (defined at top of module)
!     NOFOG - logical switch for not allowing rhum clouds (or fog)
!             to occur beneath a certain level (low_lev_cloud_index) -> 
!              nofog = true
!             to allow clouds at the lowest model level -> nofog = false
!     LOW_LEV_CLOUD_INDEX - level below which no clouds are allowed to occur 
!               when nofog = true (model level index)
!     WCUT0 - omega cutoff value for omega cloud depletion factor = 0
!     WCUT1 - omega cutoff value for omega cloud depletion factor = 1
!-----------------------------------------------------------------------

 real, intent(in), dimension (:,:,:) :: camtrh, omega
 integer, intent(in), dimension (:,:) :: llowt

!-----------------------------------------------------------------------
!      CAMTRH   tentative rel. humidity cloud amounts 
!                   (dimensioned IDIM x JDIM x kx)
!      OMEGA  Pressure vertical velocity at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      LLOWT    vertical level index upper limit for low cloud tops
!               a function of lon and latitude (dimensioned IDIMxJDIM)
!-----------------------------------------------------------------------

!===================================================================
! Arguments (intent out)
 real, intent(out), dimension (:,:,:) :: camtw
!     CAMTW - omega cloud amount at full model levels 
!               (dimensioned IDIM x JDIM x kx)
!=======================================================================

!  (Intent local)
!-----------------------------------------------------------------------
!  type loop limits 
integer kx, lk, idim, jdim
integer low_top
!-----------------------------------------------------------------------
!  type loop index variables
integer i,j, k 

!===================================================================

      kx = size(camtrh,3)
      jdim = size(camtrh,2)
      idim = size(camtrh,1)

! Initialize omega cloud amount array
     camtw = 0.0

!  define lower limit for calculating omega corrected clouds
      if (nofog) then
        lk = low_lev_cloud_index 
      else
        lk = kx
      endif

!  llowt is used to define an upper bound for calculating omega corrected 
!  clouds, i.e., estimate of highest level of a low cloud top

!-----------------------------------------------------------------------
!===================================================================

!-----------------------------------------------------------------------
! <><><><><><><>   code to compute omega corrected cloud fraction <><><>
!-----------------------------------------------------------------------
     
      do j=1,jdim
      do i=1,idim
        low_top = llowt(i,j)
        do k=low_top,lk
!  For regions of upward motion (or zero) set = rhum clouds
          if (omega (i,j,k) .le. wcut1 ) then
            camtw(i,j,k) = camtrh(i,j,k)
          endif
!  For regions of descending motion reduce clouds linearly to 0 where
!  descent = 3.6mb/hr = .1 n m-2 s-1
          if (omega (i,j,k) .gt. wcut1 .and. omega (i,j,k) .lt. wcut0 ) then
            camtw(i,j,k) = camtrh(i,j,k)*(omega(i,j,k) - wcut0)/(wcut1-wcut0)
          endif
!  For regions of descent >= 3.6mb/hr set to 0 
          if (omega (i,j,k) .ge. wcut0 ) then
            camtw(i,j,k) = 0.0
          endif
        end do
        end do
      end do
!-----------------------------------------------------------------------

end subroutine CLOUD_OMGA


!#############################################################################      
 
subroutine CLOUD_SHALLOW_CONV (theta,omega,pfull,phalf,temp,qmix,camtrh, &
                    camtsh,camtsh_mx,kminshl,kmaxshl,kbot)
                                  

!-------------------------------------------------------------------
!  This subroutine calculates shallow convective cloud amounts
!-------------------------------------------------------------------

!===================================================================
! Arguments (intent in)

! Namelist variables used in this routine (defined at top of module)
!     PSHALLOW - top of shallow convective layer (pressure level - n/m**2 )
!     L_THEQV - logical switch for turning on calculation of shallow convective 
!              clouds - true for shallow convective clouds, otherwise false 
!     WCUT1 - omega cutoff value for omega cloud depletion factor = 1
!     LOMEGA - logical switch for turning on omega correction to rhum 
!              clouds - true for omega correction, otherwise false 
!     LINVERS - logical switch for turning on calculation of marine stratus 
!              clouds - true for marine stratus, otherwise false 
!     HIGH_LEV_CLOUD_INDEX - level above which no clouds are allowed to form 
!              (model level index)
!-------------------------------------------------------------------
!  parameters used in cloud_shallow_conv

 real,  parameter :: twcfac = 0.2, crtcons = 0.0

!      TWCFAC - scaling factor for computing shallow conv cloud amt (= 0.2)
!      CRTCONS - crit value of d(Theta_e)/dP (= 0.0)


 real, intent(in), dimension (:,:,:) :: &
           theta,omega,pfull,phalf,temp,qmix,camtrh

 integer, intent(in), OPTIONAL, dimension(:,:) :: kbot

!     THETA     Potential temperature at full model levels ( Deg K )
!               (dimensioned IDIM x JDIM x kx)
!      OMEGA  Pressure vertical velocity at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      PFULL    Pressure at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      PHALF    Pressure at half model levels
!                   (dimensioned IDIM x JDIM x kx+1)
!      TEMP     Temperature (Deg K) at full model levels 
!                   (dimensioned IDIM x JDIM x kx)
!      QMIX     Mixing Ratio at full model levels 
!                   (dimensioned IDIM x JDIM x kx)
!      CAMTRH   tentative rel. humidity cloud amounts 
!                   (dimensioned IDIM x JDIM x kx)
!      KBOT      OPTIONAL; lowest model level index array
!                   (dimensioned IDIM x JDIM)
!===================================================================
! Arguments (intent out)

 real, intent(out), dimension (:,:,:) :: camtsh
 real, intent(out), dimension (:,:) :: camtsh_mx
 integer, intent(out), dimension (:,:) :: kminshl, kmaxshl

!     CAMTSH - shallow convective cloud amount at full model levels 
!               (dimensioned IDIM x JDIM x kx)
!     CAMTSH_MX  maximum value of shallow convective cloud amount within
!                  shallow convective layer. (dimensioned IDIM x JDIM)
!     KMINSHL,KMAXSHL indices of vertical levels corresponding to 
!                   final top and base of shallow convective cloud layer
!                   (dimensioned IDIM x JDIM)
!=======================================================================
!  (Intent local)

 real,    dimension(size(temp,1),size(temp,2),size(temp,3)) :: &
                   qsat,dcldnrm,theta_e,dtheta_e
 real,    dimension (size(temp,1),size(temp,2)) :: plcl
 integer, dimension (size(temp,1),size(temp,2)) :: ksiglcl,kshallow


!      QSAT - saturation vapor pressure
!                   (dimensioned IDIM x JDIM x KX)
!      PLCL - pressure of the lifting condensation level (LCL)
!                   (dimensioned IDIM x JDIM)
!      KSIGLCL - index of nearest vertical level to the LCL
!                   (dimensioned IDIM x JDIM)
!      DCLDNRM - normalization factor
!                   (dimensioned IDIM x JDIM x KX)
!      THETA_E  Equivalent potential temperature at full model levels ( Deg K )
!               (dimensioned IDIM x JDIM x kx)
!      DTHETA_E - vertical (p) derivative of theta_e
!                   (dimensioned IDIM x JDIM x KX)

! define 2-D work arrays
 real,    dimension(size(temp,1),size(temp,2)) ::  xy1, xy2, xy3

!-----------------------------------------------------------------------
!  type local variables
 integer i, j, k, kmin, kmax, kmin0, kmax0, idim, jdim,lk, lkxm1, kx, &
         kcbtop, kcount, knt_shl

!=======================================================================
!=======================================================================

  idim  = SIZE(Temp,1)
  jdim  = SIZE(Temp,2)
  kx  = SIZE(Temp,3)

 
!  define lower limit of cloud bottom
      if (nofog) then
        lk = low_lev_cloud_index
      else
        lk = kx
      endif

!  no convective clouds allowed in lowest layer
      lkxm1 = min(lk,kx-1) 

!  define upper limit of cloud tops
      kcbtop = high_lev_cloud_index

!  find index of vertical level closest to but not above the 
!  shallow convective lid (pshallow) 
     
      kshallow = 99999
      do k = kx,1,-1
        where (pfull(:,:,k) .ge. pshallow) 
           kshallow(:,:) = k
        end where
      end do

      
! Initialize shallow convective cloud amount and other related arrays
      camtsh = 0.0
      camtsh_mx = 0.0
      dcldnrm = 0.0

! calculate lcl 

! --- saturation mixing ratio 
     call compute_qs (Temp, pfull, qsat)

!  calculate equivalent potential temperature 
     
     theta_e(:,:,:) = theta(:,:,:)*exp(HLv*qmix(:,:,:)/    &
                      (Cp_Air*temp(:,:,:)))
 

!=======================================================================
! --- CALCULATE THE LIFTING CONDENSATION LEVEL, IE CLOUB BASE
!=======================================================================

! The optional argument array kbot allows for use of discontinuous vertical 
! coordinate systems. However, if the kbot index value relates to a level
! which is well above the local surface, the resulting lcl may not be 
! physically correct.

      if( PRESENT( kbot ) ) then
        do j=1,jdim
          do i=1,idim   
             k = kbot(i,j)
             xy1(i,j) =      Temp(i,j,k)
             xy2(i,j) =  MIN(MAX( qmix(i,j,k), 1.0e-6), qsat(i,j,k) )
             xy3(i,j) =     pfull(i,j,k)
          end do
        end do
      else
             xy1(:,:) =      Temp(:,:,kx)
             xy2(:,:) = MIN(MAX( qmix(:,:,kx), 1.0e-6), qsat(:,:,kx) )
             xy3(:,:) =     pfull(:,:,kx)
      end if

 

      CALL MYLCL( xy1, xy2, xy3, phalf, plcl, ksiglcl )


!=======================================================================
! --- Calculate d(theta_e)/dp
!=======================================================================

!  set dtheta_e = 0 at the first level 
      dtheta_e(:,:,1) = 0.0

      do k=2,kx-1
        dtheta_e(:,:,k) = &
       0.5*(theta_e(:,:,k+1)-theta_e(:,:,k))/(pfull(:,:,k+1)-pfull(:,:,k)) + &
       0.5*(theta_e(:,:,k)-theta_e(:,:,k-1))/(pfull(:,:,k)-pfull(:,:,k-1))
      end do

!  set dtheta_e at bottom level = detheta_e at next level up
!  (note:  bottom level value is not used in current scheme) 
      dtheta_e(:,:,kx) = dtheta_e(:,:,kx-1)


!=======================================================================
! --- Shallow Convective Cloud Calculation
!=======================================================================

      kmaxshl(:,:) = min(lkxm1,ksiglcl(:,:))
!  tentatively set kminshl to kshallow
      kminshl(:,:) = kshallow(:,:)
        
      do k = lkxm1,kcbtop,-1
        where (((k.ge.kshallow(:,:)) .and. (k.le.kmaxshl(:,:))) &
         .and. (dtheta_e(:,:,k).ge.crtcons )  &
         .and. (pfull(:,:,k) .le. plcl(:,:) ) )
           dcldnrm(:,:,k) = twcfac
        end where
      end do

!  Apply constraint that the shallow convective layer is confined to the
! lowest contigous layer between plcl and pshallow
     do j=1,jdim
       do i=1,idim
! no shallow conv clouds allowed in regions where lowest full P level
! is above the specified shallow convective lid
                   if (kshallow(i,j) .le. kx) then
         kmax=kmaxshl(i,j) 
         kmin=kshallow(i,j)
         kcount = 0
         knt_shl = 0
         do k = kmax,kmin,-1
             if (pfull(i,j,k) .le. plcl(i,j) ) then
               if (dcldnrm(i,j,k).eq.0.0 .and. knt_shl.ge.1) then
                 dcldnrm(i,j,k-1) = 0.0
                 kcount = kcount + 1
!  kminshl calculation is completed here
                 if (kcount .eq. 1) then
                   kminshl(i,j) = k
                 endif
               else
                 if (dcldnrm(i,j,k) .gt. 0.0) knt_shl = knt_shl + 1
               endif
             endif
         end do  
                   else
         dcldnrm(i,j,:) = 0.0
                   endif
       end do  
      end do  

!  calculate cloud amounts
      do k = lkxm1,kcbtop,-1
        where ((k.ge.kminshl(:,:)) .and. (k .le.kmaxshl(:,:)) &
         .and. dcldnrm(:,:,k).gt.0.0 .and. (omega(:,:,k).gt.wcut1) )
          camtsh(:,:,k) = dcldnrm(:,:,k)*camtrh(:,:,k)
          camtsh_mx (:,:) = max(camtsh_mx (:,:),camtsh (:,:,k) )      
        end where
      end do
!  impose maximum cloud overlap within the contiguous shallow convective 
!  cloud layer. Also check that an rhum cloud exists and that the condition
!  on vertical motion is satisfied
      do k = lkxm1,kcbtop,-1
        where ((k.ge.kminshl(:,:)) .and. (k .le.kmaxshl(:,:)) &
         .and. dcldnrm(:,:,k).gt.0.0 .and. (camtrh(:,:,k).gt.0.0) &
         .and. (omega(:,:,k).gt.wcut1) )
          camtsh(:,:,k) =  camtsh_mx (:,:)
        elsewhere
          camtsh(:,:,k) = 0.0 
        endwhere
      end do
!  redefine kminshl and kmaxshl based on actual shallow conv clouds
        do j=1,jdim
          do i=1,idim
                 if (camtsh_mx(i,j) .gt. 0.0) then

            kmax = kmaxshl(i,j)
            kmax0 = kmaxshl(i,j)
            kmin = kminshl(i,j)
            kmin0 = kminshl(i,j)
! find actual shallow conv cloud base by searching upwards 
            do k=kmax0,kmin0,-1
              if (camtsh(i,j,k) .gt. 0.0) then
                kmax = k
                go to 100
              endif
            end do
100    continue
! find actual shallow conv cloud top by searching downward
            do k=kmin0,kmax
              if (camtsh(i,j,k) .gt. 0.0) then
                kmin = k
                go to 200
              endif
            end do
200    continue
       kmaxshl(i,j) = kmax
       kminshl(i,j) = kmin

                 endif
          end do
        end do


!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

end subroutine CLOUD_SHALLOW_CONV

!#######################################################################

subroutine MERGE_STRAT_CLOUDS (camtrh,camtw,camtsc,llowt,camt,icld_k)



!-------------------------------------------------------------------
!  This subroutine determines a dominant stratiform cloud amount
!-------------------------------------------------------------------

!===================================================================
! Arguments (intent in)

! Namelist variables used in this routine (defined at top of module)
!     LOMEGA - logical switch for turning on omega correction to rhum 
!              clouds - true for omega correction, otherwise false 
!     LINVERS - logical switch for turning on calculation of marine stratus 
!              clouds - true for marine stratus, otherwise false 
!     NOFOG - logical switch for not allowing rhum clouds (or fog)
!             to occur beneath a certain level (low_lev_cloud_index) -> 
!              nofog = true
!             to allow clouds at the lowest model level -> nofog = false
!     LOW_LEV_CLOUD_INDEX - level below which no clouds are allowed to occur 
!               when nofog = true (model level index)
 real, intent(in), dimension (:,:,:) ::  camtrh, camtw, camtsc 
 integer, intent(in), dimension (:,:) :: llowt

!      CAMTRH   tentative cloud amounts from stratiform clouds 
!                   (dimensioned IDIM x JDIM x kx)
!      CAMTW    tentative cloud amounts omega corrected clouds 
!                   (dimensioned IDIM x JDIM x kx)
!      CAMTSC   tentativemarine stratus cloud amounts  
!                   (dimensioned IDIM x JDIM x kx)
!       LLOWT   vertical level index upper limit for low cloud tops
!               a function of lon and latitude (dimensioned IDIMxJDIM)
!===================================================================
! Arguments (intent out)

 real, intent(out), dimension (:,:,:) :: camt
 integer, intent(out), dimension (:,:,:) :: icld_k

!      CAMT     tentative (only stratiform) merged cloud amounts  
!                   (dimensioned IDIM x JDIM x kx)
!      ICLD_K     tentative merged cloud marker array  
!                   (dimensioned IDIM x JDIM x kx)



!=======================================================================
!  (Intent local)
 integer kx, lk, idim, jdim, i, j, k
 integer low_top

!=======================================================================
!=======================================================================

  idim  = SIZE(camtrh,1)
  jdim  = SIZE(camtrh,2)
  kx  = SIZE(camtrh,3)

!  define lower limit for stratiform clouds
      if (nofog) then
        lk = low_lev_cloud_index 
      else
        lk = kx
      endif


!  Initialize output cloud amounts to rel hum cloud amounts and define
!  cloud type marker values

      camt(:,:,:) = camtrh(:,:,:)
      where ( camtrh(:,:,:) .gt. 0.0)
        icld_k(:,:,:) = 1
      endwhere

!  assign omega corrected cloud amounts
        if (lomega) then
      do j=1,jdim
      do i=1,idim
        low_top =llowt(i,j)
        do k=low_top,lk
          if ( camtw(i,j,k).gt. 0.0) then 
            camt(i,j,k) = camtw(i,j,k)
            if ( camtw(i,j,k) .lt. camtrh(i,j,k) ) then
              icld_k(i,j,k) = 2
            endif
          endif
        end do
      end do
      end do
        endif

!  assign marine stratus cloud amounts
        if (linvers) then
      where ( camtsc(:,:,:) .gt. camt(:,:,:) )
        camt(:,:,:) = camtsc(:,:,:)
        icld_k(:,:,:) = 3
      endwhere
        endif

!-----------------------------------------------------------------------

end subroutine MERGE_STRAT_CLOUDS


!#######################################################################

subroutine DEF_HI_MID_LOW (phalf,lhight,lhighb,lmidt,lmidb,llowt,camt,icld_k)


!-------------------------------------------------------------------
!  This subroutine transforms a vertical profile of stratiform cloud
!  amounts into "high", "middle" and "low" cloud layers 
!  if a "thick" cloud option is in effect the current implementation
!  extends the cloud up one layer (with current tentative cloud amounts)
!-------------------------------------------------------------------

!===================================================================
! Arguments (intent in)

! Namelist variables used in this routine (defined at top of module)
!     LTHICK_HIGH - logical variable = true -> allow possibility of raising
!               high cloud tops one sigma level to increase their thickness
!               from 1 to 2 sigma levels; otherwise they remain thin 
!               (1 sigma level)
!     LTHICK_MID - logical variable = true -> allow possibility of raising
!               mid cloud tops one sigma level to increase their thickness
!               from 1 to 2 sigma levels; otherwise they remain thin 
!               (1 sigma level)
!     LTHICK_LOW - logical variable = true -> allow possibility of raising
!               low cloud tops one sigma level to increase their thickness
!               from 1 to 2 sigma levels; otherwise they remain thin 
!               (1 sigma level)
!     NOFOG - logical switch for not allowing rhum clouds (or fog)
!             to occur beneath a certain level (low_lev_cloud_index) -> 
!              nofog = true
!             to allow clouds at the lowest model level -> nofog = false
!     LOW_LEV_CLOUD_INDEX - level below which no clouds are allowed to occur 
!               when nofog = true (model level index)
 real, intent(in), dimension (:,:,:) ::  phalf 
 integer,intent(in), dimension(:,:) :: lhight, lhighb, lmidt, lmidb, llowt

!       PHALF         Pressure at half model levels
!                     (dimensioned IDIM x JDIM x kx+1)
!       LHIGHT        vertical level index upper limit for high cloud tops
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LHIGHB        vertical level index lower limit for high cloud bases
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LMIDT         vertical level index upper limit for mid cloud tops
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LMIDB         vertical level index lower limit for mid cloud bases
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LLOWT         vertical level index upper limit for low cloud tops
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!
! !===================================================================
! Arguments (intent inout)

 real, intent(inout), dimension (:,:,:) :: camt
 integer, intent(inout), dimension (:,:,:) :: icld_k

!      CAMT     tentative merged cloud amounts  
!                   (dimensioned IDIM x JDIM x kx)
!      ICLD_K     tentative merged cloud marker array  
!                   (dimensioned IDIM x JDIM x kx)

!
!=======================================================================
!  (Intent local)

!  parameter used in def_hi_mid_low

 real,  parameter :: trshld_camt = 0.25
 integer,  parameter :: nmax = 2

!      TRSHLD_CAMT - This is a cloud amount threshold value, used in
!      conjunction with the thick cloud namelist options.  If a thick
!      cloud option is on then the level above the cloud amount max
!      for a particular cloud type is considered to be part of an
!      extended cloud if it is within this fraction of the max cloud amount.
!


  real,    dimension(SIZE(camt,1),SIZE(camt,2),SIZE(camt,3)) :: delp,c_work
  integer, dimension(SIZE(camt,1),SIZE(camt,2),SIZE(camt,3)) :: ic_work
  integer, dimension(1) :: kindex

!      DELP       pressure thickness of model layers 
!                   (dimensioned IDIM x JDIM x kx)
!      C_WORK     cloud amount work array 
!                   (dimensioned IDIM x JDIM x kx)
!      IC_WORK    cloud marker work array 
!                   (dimensioned IDIM x JDIM x kx)
!      NMAX       maximum number of model levels allowed in a thick cloud

       
! loop limits and indices
 integer idim,jdim,kx, lk, i, j, k, kbtm, ktop
!  special setting of max levels for special cases of mid & high thick clouds
 integer nmax_mh
! special threshold value for allowing thick middle and high clouds
! when cloud amounts at contiguous levels are the same value (typically 1.0)
 real trshld_mh

  idim  = SIZE(camt,1)
  jdim  = SIZE(camt,2)
  kx  = SIZE(camt,3)

!=======================================================================

 trshld_mh = .0001
 nmax_mh = kx

!=======================================================================


!  define lower limit for stratiform clouds
      if (nofog) then
        lk = low_lev_cloud_index 
      else
        lk = kx
      endif

!-----------------------------------------------------------------------

! caclulate stratiform  l-m-h cloud amounts for thin or thick clouds
! current implementation of thick clouds will extend max stratiform
! cloud layer upward 1 level - assuming the level
! above the max level is within a threshold cloud amount 

! [note: the maxloc function used below will choose the lowest value index
! in case of two identical max values. In the extremely unlikely event of
! two identical max cloud amounts in a particular column section, the
! maxloc selection priority will result in a bias toward higher clouds.]

      c_work(:,:,:) = 0.0
      ic_work(:,:,:) = 0

! calculate pressure thickness of model layers, for vertcal averaging
      do k=1,kx
        delp (:,:,k) = phalf(:,:,k+1)-phalf(:,:,k)
      end do
 

      do j = 1,jdim
        do i = 1,idim

! low clouds
          kbtm = min(lk,kx-1) 
          ktop = llowt(i,j)
          kindex(:) = maxloc (camt(i,j,ktop:kbtm)) 
          k=kindex(1) + ktop - 1
          c_work(i,j,k) = camt(i,j,k)
!  for low clouds add 300 to previous cloud marker type
          ic_work(i,j,k) = icld_k(i,j,k) + 300 

          if (lthick_low) then
!            do k = kbtm,ktop,-1
call THICK_CLOUDS (camt,delp,trshld_camt,i,j,k,ktop,kbtm,nmax,c_work,ic_work)
!            end do   
          endif

! mid clouds
          kbtm = lmidb(i,j) 
          ktop = lmidt(i,j)
          kindex(:) = maxloc (camt(i,j,ktop:kbtm)) 
          k=kindex(1) + ktop - 1
          c_work(i,j,k) = camt(i,j,k)         
!  for mid clouds add 200 to previous cloud marker type
          ic_work(i,j,k) = icld_k(i,j,k) + 200         

          if (lthick_mid) then
call THICK_CLOUDS (camt,delp,trshld_mh,i,j,k,ktop,kbtm,nmax_mh,c_work,ic_work)
! call THICK_CLOUDS (camt,delp,trshld_camt,i,j,k,ktop,kbtm,nmax,c_work,ic_work)
          endif

! high clouds
          kbtm = lhighb(i,j) 
          ktop = lhight(i,j)
          kindex(:) = maxloc (camt(i,j,ktop:kbtm)) 
          k=kindex(1) + ktop - 1
          c_work(i,j,k) = camt(i,j,k)         
!  for high clouds add 100 to previous cloud marker type
          ic_work(i,j,k) = icld_k(i,j,k) + 100         

          if (lthick_high) then
call THICK_CLOUDS (camt,delp,trshld_mh,i,j,k,ktop,kbtm,nmax_mh,c_work,ic_work)
! call THICK_CLOUDS (camt,delp,trshld_camt,i,j,k,ktop,kbtm,nmax,c_work,ic_work)
          endif

        end do
      end do

! Store the resulting values of c_work and ic_work back into camt and icld_k

      where (camt(:,:,:) .eq. 0.0)
        ic_work(:,:,:) = 0
      endwhere
!  Thick clouds are currently not allowed to extend to the lowest level, but
!  fog is allowed there if nofog = f
      do k=1,kx-1
        camt(:,:,k) = c_work(:,:,k)
        icld_k(:,:,k) = ic_work(:,:,k)
      end do
!  at lowest level check for fog and assign low cloud type
      where (icld_k(:,:,kx).gt.0)
       c_work(:,:,kx) = icld_k(:,:,kx) + 300
      endwhere
      icld_k(:,:,kx) = c_work(:,:,kx)


end subroutine DEF_HI_MID_LOW


!#######################################################################

subroutine THICK_CLOUDS (camt,delp,trshld_camt,i,j,kk,ktop,kbtm,nmax, &
                         c_work,ic_work)


!-------------------------------------------------------------------
!  This subroutine returns "thick" cloud amounts to routine def_hi_mid_low
!  if the criteria for those clouds are met.  First the level of the max 
!  cloud amount in a low, middle  or high cloud regime is determined.
!  A thick cloud layer consisting of a maximum of "nmax"  contiguous levels 
!  is determined to exist by looking for a cloud amount value one level
!  above that is within a specified threshold amount (trshld_camt). If
!  a contiguous cloud is found and nmax = 2 the search is over.  If no
!  contiguous cloud is found or nmax > 2 the search continues one level
!  below. If contiguous cloud is found both above and below and nmax > 3
!  the search will continue using the same threshold criteria 2 levels
!  above.  The search will continue to "zigzag" in this manner as long as
!  no maximum cloud thickness criteria has not been violated and as long
!  as there is no break in contiguous cloud amounts as described above.
!  Eventually a pressure thickness value will overide the nmax levels 
!  criteria for determining the maximum thickness of the clouds.  Another
!  factor limiting the thickness of clouds is that the cloud extent must 
!  remain within its respective high, middle, low cloud regime boundaries
!  ( i.e., ktop to kbtm). 
!-------------------------------------------------------------------
!
! !===================================================================
! Arguments (intent in)



 real, intent(in), dimension (:,:,:) :: camt,delp
 real, intent(in) :: trshld_camt
 integer, intent(in) :: i,j,kk,ktop,kbtm,nmax

!      CAMT     tentative merged cloud amounts  
!                   (dimensioned IDIM x JDIM x kx)
!      DELP       pressure thickness of model layers 
!                   (dimensioned IDIM x JDIM x kx)
!      I,J,KK     indicies for longitude, latitude and 
!                 vertical level (at level of max cloud amount)
!      KTOP,KBTM  - vertical index valid range for cloud top, cloud bottom
!      NMAX -       maximum number of model levels allowed in a thick cloud
!      TRSHLD_CAMT - This is a cloud amounht threshold value, used in
!                    conjunction with the thick cloud namelist options.   
!                    If a thick cloud option is on then level(s) adjacent  
!                    to the level of max cloud amount for a particular 
!                    cloud type is considered to be part of an extended cloud 
!                    if it is within this fraction of the max cloud amount.
!

!
!=======================================================================
! Arguments (intent inout)


  real, intent(inout), dimension(:,:,:) :: c_work
  integer,intent(inout), dimension(:,:,:) :: ic_work

!      C_WORK     cloud amount work array 
!                   (dimensioned IDIM x JDIM x kx)
!      IC_WORK    cloud marker work array 
!                   (dimensioned IDIM x JDIM x kx)
!
!=======================================================================
!  (Intent local)
 integer ncount,no_up,no_down,k,c_top,c_btm
! loop index
 integer ki,pass_count
 real delp_sum,thick_cld_amt
!      NCOUNT -       counter for number of model levels in a thick cloud
!      NO_UP -       if = 1 -> no more upward extension of cloud is possible
!      NO_DOWN -      if = 1 -> no more downward extension of cloud is possible
!      K  -          vertical level index of thick cloud as it being extended
!      C_TOP, C_BTM - Top and bottom of calculated thick cloud 
!      DELP_SUM - sum of pressure weights for thcik cloud layer
!      THICK_CLD_AMT - pressure weighted thick cloud amount     
!=======================================================================
!=======================================================================
      ncount = 1
      no_up = 0
      no_down = 0
      k = kk
      c_top = kk
      c_btm = kk

      pass_count = 0

!  no clouds to process
      if (camt(i,j,kk) .eq.0.0) go to 300 

100   continue

      pass_count = pass_count + 1

     
                     if (no_up .eq. 0) then
! First check to extend cloud upward
              k = c_top
              if (camt(i,j,k) .gt.0.0 .and. k.gt.ktop .and. k.le.kbtm) then
                if ((camt(i,j,kk)-camt(i,j,k-1)).le.trshld_camt) then
                  c_top = c_top - 1
                  ic_work(i,j,k-1) = ic_work(i,j,k)
                  ncount = ncount + 1
                else
! no upward extension of cloud possible
                  no_up = 1
                endif
             else
                no_up = 1
             endif
! check for reaching limiting number of levels in a thick cloud
     if (ncount .eq. nmax) go to 200
                     endif

! Next check to extend cloud downward
                     if (no_down .eq. 0) then
              k = c_btm
              if (camt(i,j,k) .gt.0.0 .and. k.ge.ktop .and. k.lt.kbtm) then
                if ((camt(i,j,kk)-camt(i,j,k+1)).le.trshld_camt) then
                  c_btm = c_btm + 1
                  ic_work(i,j,k+1) = ic_work(i,j,k)
                  ncount = ncount + 1
                else
! no downward extension of cloud possible
                  no_down = 1
                endif
             else
                no_down = 1
             endif
! check for reaching limiting number of levels in a thick cloud
      if (ncount .eq. nmax) go to 200
                     endif

!  Further tests for completion of thick cloud calculations
!  Make sure upward or downward extension of the cloud is still allowed
!  and the cloud level index is still within the allowable range
       if ( ( no_up.eq.1 .or. c_top.eq.ktop) .and. &
          ( no_down.eq.1 .or. c_btm.eq.kbtm) ) go to 200

! loop back to do another pass at extending the thick cloud
       go to 100

200   continue

!  calculate press thickness weighted thick clouds
      delp_sum = 0.0
      do ki=c_top,c_btm
        delp_sum = delp_sum + delp(i,j,ki)
      end do
      thick_cld_amt = 0.0
      do ki=c_top,c_btm
        thick_cld_amt = thick_cld_amt + camt(i,j,ki)*delp(i,j,ki)
      end do
      thick_cld_amt = thick_cld_amt/delp_sum

! broadcast thick cloud amount throughout thick cloud layer
      do ki=c_top,c_btm
        c_work(i,j,ki) = thick_cld_amt
      end do

300   continue

          
!===================================================================
end subroutine THICK_CLOUDS

!#######################################################################

subroutine MERGE_CNV_CLOUDS (camtcnv,camtsh,camtsh_mx,kminshl,kmaxshl, &
                           lhight,lhighb,lmidt,lmidb,llowt,camt,icld_k)



!-------------------------------------------------------------------
!  This subroutine adds convective cloud amounts to the already merged 
!  stratiform cloud amounts
!-------------------------------------------------------------------

!===================================================================
! Arguments (intent in)

!-------------------------------------------------------------------
! Namelist variables used in this routine (defined at top of module)
!     LCNVCLD - logical switch for turning on calculation of deep convective 
!              clouds - true for deep convective clouds, otherwise false 
!     L_THEQV - logical switch for turning on calculation of shallow convective 
!              clouds - true for shallow convective clouds, otherwise false 
!     HIGH_LEV_CLOUD_INDEX - level above which no clouds are allowed to form 
!              (model level index)
!     NOFOG - logical switch for not allowing rhum clouds (or fog)
!             to occur beneath a certain level (low_lev_cloud_index) -> 
!              nofog = true
!             to allow clouds at the lowest model level -> nofog = false
!     LOW_LEV_CLOUD_INDEX - level below which no clouds are allowed to occur 
!               when nofog = true (model level index)
!!-------------------------------------------------------------------

 real, intent(in), dimension (:,:,:) :: camtcnv, camtsh
 real, intent(in), dimension (:,:) :: camtsh_mx
 integer, intent(in), dimension (:,:) :: kminshl, kmaxshl
 integer, intent(in), dimension(:,:) :: lhight, lhighb, lmidt, lmidb, llowt

!     CAMTCNV - deep convective cloud amounts at full model levels  
!                   (dimensioned IDIM x JDIM x kx)
!     CAMTSH - shallow convective cloud amount at full model levels 
!               (dimensioned IDIM x JDIM x kx)
!     CAMTSH_MX  maximum value of shallow convective cloud amount within
!                  shallow convective layer. (dimensioned IDIM x JDIM)
!     KMINSHL,KMAXSHL indices of vertical levels corresponding to 
!                   final top and base of shallow convective cloud layer
!                   (dimensioned IDIM x JDIM)
!       LHIGHT        vertical level index upper limit for high cloud tops
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LHIGHB        vertical level index lower limit for high cloud bases
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LMIDT         vertical level index upper limit for mid cloud tops
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LMIDB         vertical level index lower limit for mid cloud bases
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LLOWT         vertical level index upper limit for low cloud tops
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
! !===================================================================
! Arguments (intent inout)

 real, intent(inout), dimension (:,:,:) :: camt
 integer, intent(inout), dimension (:,:,:) :: icld_k

!      CAMT     tentative merged cloud amounts  
!                   (dimensioned IDIM x JDIM x kx)
!      ICLD_K     tentative merged cloud marker array  
!                   (dimensioned IDIM x JDIM x kx)


!=======================================================================
!  (Intent local)
  real,    dimension(SIZE(camt,1),SIZE(camt,2)) ::   cmax_strat
  real,    dimension(SIZE(camt,1),SIZE(camt,2),SIZE(camt,3)) ::   c_work
  integer, dimension(SIZE(camt,1),SIZE(camt,2),SIZE(camt,3)) ::   ic_work

!      CMAX_STRAT   max stratiform cloud amount work array 
!                   (dimensioned IDIM x JDIM)
!      C_WORK     cloud amount work array 
!                   (dimensioned IDIM x JDIM x kx)
!      IC_WORK    cloud marker work array 
!                   (dimensioned IDIM x JDIM x kx)

 integer i,j,k,kx, lk, kcbtop, kbtm , ktop, idim, jdim

!=======================================================================
!=======================================================================

  idim  = SIZE(camt,1)
  jdim  = SIZE(camt,2)
  kx  = SIZE(camt,3)


      if (nofog) then
        lk = low_lev_cloud_index 
      else
        lk = kx
      endif
!  define upper limit of cloud tops
      kcbtop = high_lev_cloud_index

!  initialize work arrays
      c_work(:,:,:) = camt(:,:,:)

!  decide whether or not to replace stratiform low cloud within the
!  contiguous shallow convective layer with shallow convective cloud
      
      if (l_theqv) then
!  find maximum stratiform low cloud amount within shallow conv layer
!  assume shallow convective layer is always below upper bound for low clouds

      do j = 1,jdim
        do i = 1,idim
          kbtm = kmaxshl(i,j)
          ktop = kminshl(i,j)
          cmax_strat(i,j) = 0.0
          do k = kbtm,ktop,-1
            cmax_strat(i,j) = max (cmax_strat(i,j),camt(i,j,k))
          end do
        end do
      end do

!  check to see if shallow convective cloud will be dominant type
!  if so reset cloud amounts to shallow amount and reset cloud marker
      do k = lk,kcbtop,-1
        where ( (k.ge.kminshl(:,:)) .and. (k.le.kmaxshl(:,:)) .and. &
              (camtsh_mx(:,:) .gt. cmax_strat(:,:)) )
          c_work(:,:,k) =  camtsh (:,:,k)
          icld_k(:,:,k) = 304
        endwhere
      end do

      endif

!  randomly overlap stratiform and deep convective cloud amounts.
!  also check to see if the convective cloud is dominant based on
!  a comparison of cloud amounts (pre-random overlap), if so assign
!  cloud marker values based on H-M-L and cloud type (convective = 5)
  
      ic_work(:,:,:) = mod(icld_k(:,:,:),100)
       
        if (lcnvcld) then

      do k = lk,kcbtop,-1

        where ( ic_work(:,:,k) .le. 3)
          c_work(:,:,k) = camtcnv(:,:,k) + (1.0-camtcnv(:,:,k))*camt(:,:,k)
        endwhere

!  low cloud range
          where ( (camtcnv(:,:,k).gt.camt(:,:,k)) .and. &
             (k.ge.llowt(:,:) .and. k.le.lk) .and. (ic_work(:,:,k).le.3) )
            icld_k(:,:,k) = 305
          endwhere
!  mid cloud range
          where ( (camtcnv(:,:,k).gt.camt(:,:,k)) .and. &
             (k.ge.lmidt(:,:) .and. k.le.lmidb(:,:)) )
            icld_k(:,:,k) = 205
          endwhere
!  high cloud range
          where ( (camtcnv(:,:,k).gt.camt(:,:,k)) .and. &
             (k.ge.lhight(:,:) .and. k.le.lhighb(:,:)) )
            icld_k(:,:,k) = 105
          endwhere

      end do

        endif

!  re-define convective cloud amount as the max of deep and shallow convective
!  cloud at levels where shallow convective clouds are currently defined.

        if (lcnvcld .and. l_theqv) then

          where ( (camtcnv(:,:,:).gt.c_work(:,:,:)) .and. &
             (ic_work(:,:,:) .eq. 4 ) )
            c_work(:,:,:) = camtcnv(:,:,:)
            icld_k(:,:,:) = 305
          endwhere

        endif

! Store the resulting values of c_work back into camt

      camt(:,:,:) = c_work(:,:,:)

!-----------------------------------------------------------------------

end subroutine MERGE_CNV_CLOUDS

!#######################################################################

subroutine LAYR_TOP_BASE (icld_k,nclds,cldtop,cldbas,icld)



!-------------------------------------------------------------------
!  This subroutine determines cloud layers, tops and bases
!-------------------------------------------------------------------

!===================================================================
! Arguments (intent in)

! Namelist variables used in this routine (defined at top of module)
!     LTHICK_HIGH - logical variable = true -> allow possibility of raising
!               high cloud tops one sigma level to increase their thickness
!               from 1 to 2 sigma levels; otherwise they remain thin 
!               (1 sigma level)
!     LTHICK_MID - logical variable = true -> allow possibility of raising
!               mid cloud tops one sigma level to increase their thickness
!               from 1 to 2 sigma levels; otherwise they remain thin 
!               (1 sigma level)
!     LTHICK_LOW - logical variable = true -> allow possibility of raising
!               low cloud tops one sigma level to increase their thickness
!               from 1 to 2 sigma levels; otherwise they remain thin 
!               (1 sigma level)
!     HIGH_LEV_CLOUD_INDEX - level above which no clouds are allowed to form 
!              (model level index)
!     NOFOG - logical switch for not allowing rhum clouds (or fog)
!             to occur beneath a certain level (low_lev_cloud_index) -> 
!              nofog = true
!             to allow clouds at the lowest model level -> nofog = false
!     LOW_LEV_CLOUD_INDEX - level below which no clouds are allowed to occur 
!               when nofog = true (model level index)
!-------------------------------------------------------------------

integer, intent(in), dimension (:,:,:) :: icld_k

!      ICLD_K     tentative merged cloud marker array  
!                   (dimensioned IDIM x JDIM x kx)
!===================================================================
! Arguments (intent out)

integer, intent(out), dimension(:,:,:) :: cldtop,cldbas,icld
integer, intent(out), dimension(:,:)  :: nclds

!      NCLDS       number of (random overlapping) clouds in column and also
!                     the current # for clouds to be operating on
!                   (dimensioned IDIM x JDIM )
!      CLDTOP     index of cloud tops (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      CLDBAS    index of cloud bottoms (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      ICLD          marker array of cloud types (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)



!=======================================================================
!  (Intent local)

  integer, dimension(SIZE(icld_k,1),SIZE(icld_k,2),SIZE(icld_k,3)) :: &  
                 ctop_work,cbas_work,ic_work

!      CTOP_WORK    cloud top work array (at model levels)
!                   (dimensioned IDIM x JDIM x kx)
!      CBAS_WORK    cloud bsse work array (at model levels)
!                   (dimensioned IDIM x JDIM x kx)
!      IC_WORK    cloud marker work array 
!                   (dimensioned IDIM x JDIM x kx)
 integer k,kcbtop,i,j,idim, jdim, kx, lk, kp, kpmax, maxcld

!=======================================================================
!=======================================================================

  idim  = SIZE(icld_k,1)
  jdim  = SIZE(icld_k,2)
  kx  = SIZE(icld_k,3)


!  define lower limit of cloud bottom
      if (nofog) then
        lk = low_lev_cloud_index
      else
        lk = kx
      endif
!  define upper limit of cloud tops
      kcbtop = high_lev_cloud_index

!  Initialize output arrays

      nclds(:,:) = 0
      icld(:,:,:) = 0
      cldtop(:,:,:) = 0
      cldbas(:,:,:) = 0

!  Initialize work arrays

      ctop_work(:,:,:) = 0
      cbas_work(:,:,:) = 0
      ic_work(:,:,:) = mod(icld_k(:,:,:),100)


!  locate distinct cloud layers, tops and bases

      do k=kcbtop,lk

!  first locate tops and increment cloud layer counter, nclds

!  check for no clouds above any type of cloud layer 
      where ( ic_work(:,:,k).ge.1 .and. ic_work(:,:,k-1).eq.0)
        nclds(:,:) = nclds(:,:) + 1
        ctop_work(:,:,k) = k
      endwhere
!  check for stratiform cloud layer above convective cloud layer 
!  anvil and super anvil (types 6,7) are stratiform type clouds 
      where ( (ic_work(:,:,k).ge.4 .and. ic_work(:,:,k).le.5) .and. &
             (ic_work(:,:,k-1).ge.1 .and. ic_work(:,:,k-1).le.3) .or. &
             (ic_work(:,:,k-1).ge.6) )
        nclds(:,:) = nclds(:,:) + 1
        ctop_work(:,:,k) = k
      endwhere
!  check for convective cloud layer above stratiform cloud layer 
      where ( (ic_work(:,:,k-1).ge.4 .and. ic_work(:,:,k-1).le.5) .and. &
             (ic_work(:,:,k).ge.1 .and. ic_work(:,:,k).le.3) .or. &
             (ic_work(:,:,k).ge.6) )
        nclds(:,:) = nclds(:,:) + 1
        ctop_work(:,:,k) = k
      endwhere

        if (k.lt.lk) then

!  locate bases

!  check for no clouds above any type of cloud layer 
      where ( ic_work(:,:,k).ge.1 .and. ic_work(:,:,k+1).eq.0)
        cbas_work(:,:,k) = k
      endwhere
!  check for stratiform cloud layer above convective cloud layer 
!  anvil and super anvil (types 6,7) are stratiform type clouds 
      where ( (ic_work(:,:,k+1).ge.4 .and. ic_work(:,:,k+1).le.5) .and. &
             (ic_work(:,:,k).ge.1 .and. ic_work(:,:,k).le.3) .or. &
             (ic_work(:,:,k).ge.6) )
        cbas_work(:,:,k) = k
      endwhere
!  check for convective cloud layer above stratiform cloud layer 
      where ( (ic_work(:,:,k).ge.4 .and. ic_work(:,:,k).le.5) .and. &
             (ic_work(:,:,k+1).ge.1 .and. ic_work(:,:,k+1).le.3) .or. &
             (ic_work(:,:,k+1).ge.6) )
        cbas_work(:,:,k) = k
      endwhere

        endif


          do j = 1,jdim
            do i = 1,idim
        if (k.eq.lk .and. ic_work(i,j,k).ge.1) then

          cbas_work(i,j,k) = lk

        endif
            end do
          end do

      end do

!-----------------------------------------------------------------------

!  compress output cloud tops, bases and cloud markers to cloud levels

      do j = 1,jdim
        do i = 1,idim

          if (nclds(i,j) .ge.1) then

          kpmax = nclds(i,j)
          kp = 1
            do k = kcbtop,lk

              if (ctop_work(i,j,k).ne. 0) then
                cldtop(i,j,kp) = ctop_work(i,j,k)
                icld(i,j,kp) = icld_k(i,j,k)
              endif

              if (cbas_work(i,j,k).ne. 0) then
                cldbas(i,j,kp) = cbas_work(i,j,k)
                kp = kp + 1 
              endif

              if (kp .gt. kpmax) go to 999              

            end do
   
999     continue
     
        endif

        end do
      end do


!-----------------------------------------------------------------------

!  error check for too many cloud layers

! find maximum number of cloud layers
      maxcld  = maxval(nclds(:,:))
                   if (maxcld .gt. kx/2) then
     print *,'pe, NCLDS =', mpp_pe(),nclds
                 call error_mesg ('diag_cloud, layr_top_base',  &
                   'NCLDS too large', FATAL)
                    endif


!-----------------------------------------------------------------------

end subroutine LAYR_TOP_BASE

!#######################################################################


subroutine CLDAMT_MN (nclds,cldtop,cldbas,phalf,camt,cldamt)



!-------------------------------------------------------------------
!  This subroutine computes mass-weighted vertical mean cloud amount 
!  for each distinct cloud layer
!-------------------------------------------------------------------

!===================================================================
! Arguments (intent in)

! Namelist variables used in this routine (defined at top of module)
!     NOFOG - logical switch for not allowing rhum clouds (or fog)
!             to occur beneath a certain level (low_lev_cloud_index) -> 
!              nofog = true
!             to allow clouds at the lowest model level -> nofog = false
!     LOW_LEV_CLOUD_INDEX - level below which no clouds are allowed to occur 
!               when nofog = true (model level index)
!-------------------------------------------------------------------

 integer, intent(out), dimension(:,:,:) :: cldtop,cldbas
 integer, intent(out), dimension(:,:)  :: nclds
 real, intent(in), dimension (:,:,:) ::  phalf

!      NCLDS       number of (random overlapping) clouds in column and also
!                     the current # for clouds to be operating on
!                   (dimensioned IDIM x JDIM )
!      CLDTOP     index of cloud tops (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      CLDBAS    index of cloud bottoms (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      PHALF    Pressure at half model levels
!                   (dimensioned IDIM x JDIM x kx+1)
!===================================================================
! Arguments (intent inout)

 real, intent(inout), dimension (:,:,:) :: camt

!      CAMT     tentative merged cloud amounts  
!                   (dimensioned IDIM x JDIM x kx)
!=======================================================================
 real, intent(out), dimension(:,:,:) :: cldamt

!      CLDAMT   cloud amount (fraction) (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)

!=======================================================================
!  (Intent local)


  real, dimension(SIZE(camt,1),SIZE(camt,2),SIZE(camt,3)) :: delp,c_work,weight

!      DELP       pressure thickness of model layers 
!                   (dimensioned IDIM x JDIM x kx)
!      C_WORK     cloud amount work array (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      WEIGHT     vertical sum of press thickness over model levels for each
!                 cloud level (dimensioned IDIM x JDIM x kx)

! loop limits and indices
 integer idim,jdim,kx, kpr, lk, i, j, k
 integer maxcld

!=======================================================================
!=======================================================================

  idim  = SIZE(camt,1)
  jdim  = SIZE(camt,2)
  kx  = SIZE(camt,3)


!  define lower limitfor stratiform clouds
      if (nofog) then
        lk = low_lev_cloud_index 
      else
        lk = kx
      endif

! find maximum number of cloud layers
      maxcld  = maxval(nclds(:,:))

!-----------------------------------------------------------------------
 ! initialize output and work cloud amount arrays 

      cldamt(:,:,:) = 0.0
      c_work(:,:,:) = 0.0
      weight(:,:,:) = 0.0
!-----------------------------------------------------------------------

! calculate pressure thickness of model layers, for vertcal averaging
      do k=1,kx
        delp (:,:,k) = phalf(:,:,k+1)-phalf(:,:,k)
      end do
!-----------------------------------------------------------------------

! calculate mass weighted vertical mean cloud amount
      do kpr=1,maxcld
        do j=1,jdim
          do i=1,idim 
  
             if (cldtop(i,j,kpr) .lt. cldbas(i,j,kpr)) then

               do k=cldtop(i,j,kpr),cldbas(i,j,kpr) 
                 c_work(i,j,kpr) = c_work(i,j,kpr) + &
                 delp(i,j,k)*camt(i,j,k)
                 weight(i,j,kpr) = weight(i,j,kpr) + delp(i,j,k)
               end do
               cldamt(i,j,kpr) = c_work(i,j,kpr)/weight(i,j,kpr)
! reset cloud amounts at model sigma levels to mass weighted for the layer
               do k=cldtop(i,j,kpr),cldbas(i,j,kpr) 
                 camt(i,j,k) = cldamt(i,j,kpr)
               end do

             else if ( (cldtop(i,j,kpr) .eq. cldbas(i,j,kpr)) .and. &
             (cldtop(i,j,kpr).gt.0 .and. cldtop(i,j,kpr).le.kx)) then
 
               k = cldtop(i,j,kpr)
               cldamt(i,j,kpr) = camt(i,j,k)

             else if ( (cldtop(i,j,kpr) .gt. cldbas(i,j,kpr)) .or. &
             (cldtop(i,j,kpr).lt.0) .or. (cldbas(i,j,kpr).lt.0)) then

     print *,'pe, i,j,kpr, cldtop,cldbas =', mpp_pe(),i,j,kpr, &
         cldtop(i,j,kpr),cldbas(i,j,kpr),cldtop(i,j,kpr),cldbas(i,j,kpr)
                 call error_mesg ('diag_cloud, cldamt_mn',  &
                   'invalid cldtop and/or cldbas', FATAL)

             endif
          end do
        end do
      end do


end subroutine CLDAMT_MN 

!#############################################################################      
 
subroutine ANVIL_CIRRUS (lhighb,lmidb,nclds,cldtop,cldbas,cnvcntq,lgscldelq,  &
                   pfull,phalf,icld)
                                  


!------------------------------------------------------------------------
!  This subroutine calculates anvil and super anvil cirrus cloud amounts
!------------------------------------------------------------------------

!===================================================================
! Arguments (intent in)

!-------------------------------------------------------------------
! Namelist variables used in this routine (defined at top of module)
!     HIGH_LEV_CLOUD_INDEX - level above which no clouds are allowed to form 
!              (model level index)
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  parameters used in anvil_cirrus

 real,  parameter :: cnvdelqmn = 0.0,lgscldelqmn = 1.0e-9 
 real,  parameter :: lcldbas_prc = .55, lcldtop_prc = .225

!      CNVDELQMN - min threshold value of accumulated count of mixing ratio 
!                  change due to convection above which convective clouds are 
!                  allowed
!                  ( It is set = 0.0 as a parameter, because it not anticipated
!                    that it would be changed.)
!      LGSCLDELQMN - min threshold value of averaged mixing ratio rate of
!                  change due to large scale condensation
!      LCLDBAS_PRC - Threshold sigma value for base of an "anvil cirrus 
!                    indicator" layer
!      LCLDBAS_PRC - Threshold sigma value for top of an "anvil cirrus 
!                    indicator" layer
!-------------------------------------------------------------------

 integer, intent(out), dimension(:,:,:) :: cldtop,cldbas
 integer, intent(out), dimension(:,:)  ::  nclds,lhighb,lmidb
 real, intent(in), dimension (:,:,:) ::  lgscldelq,cnvcntq,pfull,phalf

!-------------------------------------------------------------------
!       LHIGHB        vertical level index lower limit for high cloud bases
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       LMIDB         vertical level index lower limit for mid cloud bases
!                     a function of lon and latitude (dimensioned IDIMxJDIM)
!       NCLDS   number of (random overlapping) clouds in column and also
!                     the current # for clouds to be operating on
!                   (dimensioned IDIM x JDIM )
!      CLDTOP   index of cloud tops (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      CLDBAS   index of cloud bottoms (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      LGSCLDELQ  Averaged rate of change in mix ratio due to lg scale precip 
!               at full model levels  
!               (dimensioned IDIM x JDIM x kx)
!      CNVCNTQ  Accumulated count of change in mix ratio due to conv precip 
!               at full model levels  
!               (dimensioned IDIM x JDIM x kx)
!      PFULL    Pressure at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      PHALF    Pressure at half model levels
!                   (dimensioned IDIM x JDIM x kx+1)
!===================================================================
! Arguments (intent inout)
 integer, intent(out), dimension (:,:,:) :: icld
!       ICLD          marker array of cloud types/heights (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!=======================================================================
!  (Intent local)

 integer, dimension (size(pfull,1),size(pfull,2),size(pfull,3)) :: ic_work
!-------------------------------------------------------------------
!      IC_WORK    cloud marker work array 
!                   (dimensioned IDIM x JDIM x kx)
!       LCLDT         Top of anvil cirrus 
!                     (model level index)
!       LCLDB         Base of anvil cirrus 
!                     (model level index)
!       LPRCTOP       Base of region in search for anvil cirrus markers
!                     (model level index)
!       LPRCBAS       Base of region in search for anvil cirrus markers
!                     (model level index)
!       SIG_TOP_BASE  Sigma value of top/base of region in search for 
!                      anvil cirrus markers
!                     

 real sig_top_base
 integer lcldt, lcldb, lprcbas, lprctop
 integer i, j, idim, jdim, kx
!-----------------------------------------------------------------------
!  type loop index variable
 integer k,kl 

!===================================================================

      idim  = SIZE(pfull,1)
      jdim  = SIZE(pfull,2)
      kx = size(pfull,3)

! Initialize cloud type work array
      ic_work (:,:,:) = icld(:,:,:)


!-----------------------------------------------------------------------
! <><><><><><><><>   search for anvil cirrus  <><><><><><><><>
!-----------------------------------------------------------------------
      
      do j=1,jdim
      do i=1,idim
        do k=1,kx
          sig_top_base = pfull(i,j,k)/phalf(i,j,kx+1)
          if ( sig_top_base .gt. lcldtop_prc) then
             lprctop = k-1
             go to 100
          endif
        end do

100     continue

        do k=1,kx
          sig_top_base = pfull(i,j,k)/phalf(i,j,kx+1)
          if ( sig_top_base .gt. lcldbas_prc) then
             lprcbas = k-1
             go to 150
          endif
        end do

150     continue
           
        do k=1,nclds(i,j)
          lcldt = cldtop(i,j,k)
          lcldb = cldbas(i,j,k)


! Must have a high cloud for any anvil cirrus

                    if (icld(i,j,k) .ne. 101) go to 400

! First check whether super-anvil cirrus criteria is satisfied
! This criteria differs from v197, because calculation of convective precip for 
! high cloud region only is unreliable. Instead it is required that there has  
! been a change in mixing ratio due to convection at at least one time step during
! the time interval over which the clouds are being calculated, at ALL levels 
! between the base of the anvil cirrus indicator region (lprcbas- clculated from preset parameter)
! and the higher of either the actual top of the high cloud or the previoulsy computed top of the 
! Anvil Cirrus cloud indicator region.

           if (lcldt .lt. lprctop) then
              lprctop = lcldt
           endif
           do kl=lprcbas,lprctop,-1
             if (cnvcntq(i,j,kl) .lt. 1.0) then
               go to 200
             endif
           end do
           ic_work(i,j,k) = 107
! found super-anvil cirrus, no need to search for anvil cirrus
           go to 400

200        continue

! search for anvil cirrus using criteria that a change in mixing ratio due to 
! either convection or large scale condensation at any level within an actual high 
! cloud takes places during the time interval over which the clouds are being 
! calculated.

          do kl=lcldb,lcldt,-1
            if (cnvcntq(i,j,kl) .gt. 0.0) go to 300
            if (kl.le.lhighb(i,j) .and. lgscldelq(i,j,kl).gt.lgscldelqmn) &
                go to 300
          end do
          go to 400
300       continue
               ic_work(i,j,k) = 106
   
        end do

400   continue

!                    endif

      end do
      end do


!-----------------------------------------------------------------------
! <><><><><><><><>   update cloud marker array  <><><><><><><><>
!-----------------------------------------------------------------------
      
      icld(:,:,:) = ic_work(:,:,:)       


!-----------------------------------------------------------------------

end subroutine ANVIL_CIRRUS



!#######################################################################

subroutine CLD_LAYR_MN_TEMP_DELP &
                        (nclds,cldtop,cldbas,temp,phalf,tempcld,delp_true)



!-------------------------------------------------------------------
!  This subroutine computes the mass-weighted mean cloud temperature
!  and true cloud pressure thickness of distinct cloud layers.
!-------------------------------------------------------------------

!===================================================================
! Arguments (intent in)

! Namelist variables used in this routine (defined at top of module)
!     HIGH_LEV_CLOUD_INDEX - level above which no clouds are allowed to form 
!              (model level index)
!     NOFOG - logical switch for not allowing rhum clouds (or fog)
!             to occur beneath a certain level (low_lev_cloud_index) -> 
!              nofog = true
!             to allow clouds at the lowest model level -> nofog = false
!     LOW_LEV_CLOUD_INDEX - level below which no clouds are allowed to occur 
!               when nofog = true (model level index)
!-------------------------------------------------------------------

! Arguments (intent in)

integer, intent(in), dimension(:,:)  :: nclds
integer, intent(in), dimension(:,:,:) :: cldtop,cldbas
real, intent(in), dimension(:,:,:) :: temp,phalf

!      NCLDS       number of (random overlapping) clouds in column and also
!                     the current # for clouds to be operating on
!                   (dimensioned IDIM x JDIM )
!      CLDTOP     index of cloud tops (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      CLDBAS    index of cloud bottoms (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      TEMP        Temperature (Deg K) at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      PHALF    Pressure at half model levels
!                   (dimensioned IDIM x JDIM x kx+1)
!===================================================================
! Arguments (intent out)

real, intent(out), dimension(:,:,:) :: tempcld,delp_true

!      TEMPCLD    cloud layer mean temperature (degrees Kelvin)
!                    (at cloud levels)
!      DELP_TRUE  true cloud pressure thickness of distinct cloud layers 
!                    (at cloud levels)



!=======================================================================
!  (Intent local)

 integer, dimension (size(cldtop,1),size(cldtop,2)) :: cldt,cldb

!      CLDT       cloud top index work array 
!                   (dimensioned IDIM x JDIM x kx)
!      CLDB       cloud base index work array 
!                   (dimensioned IDIM x JDIM x kx)
 integer i,j,k,kk,kx,idim,jdim, lk, kcbtop, maxcld

!=======================================================================
!=======================================================================

  idim  = SIZE(temp,1)
  jdim  = SIZE(temp,2)
  kx  = SIZE(temp,3)


!  define lower limit of cloud bottom
      if (nofog) then
        lk = low_lev_cloud_index
      else
        lk = kx
      endif
!  define upper limit of cloud tops
      kcbtop = high_lev_cloud_index

!  Initialize output arrays

      tempcld(:,:,:) = 0
      delp_true(:,:,:) = 0

! find maximum number of clouds
      maxcld  = maxval(nclds(:,:))

                   if (maxcld .ge. 1) then

!-----------------------------------------------------------------------
! <><><><><><><><>   calc pressure thickness <><><><><><><><>
!-----------------------------------------------------------------------

      do k = 1,maxcld

! Initialize internal arrays
      cldt = 0
      cldb = 0
         where (nclds(:,:) .ge. k)
           cldt(:,:) = cldtop(:,:,k)      
           cldb(:,:) = cldbas(:,:,k) 
         end where

         do j=1,jdim
         do i=1,idim
           if (k .le. nclds(i,j)) then
             delp_true(i,j,k) = phalf(i,j,cldb(i,j)+1)-phalf(i,j,cldt(i,j))
           endif
         end do
         end do
       end do



!-----------------------------------------------------------------------
! <><><><><><><><>   calc mass-weghted mean cloud temp <><><><><><><><>
!-----------------------------------------------------------------------

      do j=1,jdim
      do i=1,idim
        do kk= 1,nclds(i,j)
          do k=cldtop(i,j,kk),cldbas(i,j,kk)
            tempcld(i,j,kk) = tempcld(i,j,kk) + &
            (phalf(i,j,k+1)-phalf(i,j,k))*temp(i,j,k)/delp_true(i,j,kk)
          end do
        end do
      end do
      end do

                   endif


!-----------------------------------------------------------------------

end subroutine CLD_LAYR_MN_TEMP_DELP

!#######################################################################

  SUBROUTINE DIAG_CLOUD_INIT( ix,iy,kx, ierr )

!=======================================================================
! ***** INITIALIZE Predicted Cloud Scheme
!=======================================================================


!---------------------------------------------------------------------
! Arguments (Intent in)
!  parmameter mxband = max number of radiative bands to be considered for some
!              cloud properties (defined at top of module)
!---------------------------------------------------------------------
 integer, intent(in) :: ix, iy, kx
!      INPUT
!      ------

!      IX, IY, KX   Dimensions for global storage arrays (2- horiz, vert)
!---------------------------------------------------------------------
! Arguments (Intent out)
!---------------------------------------------------------------------
 integer, intent(out) :: ierr

!      OUTPUT
!      ------

!      IERR     Error flag

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------
 integer  unit, io, ierrnml, logunit
 integer  id_restart
 character(len=32) :: fname

!=====================================================================


!---------------------------------------------------------------------
! --- Read namelist
!---------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=diag_cloud_nml, iostat=io)
   ierr = check_nml_error(io,"diag_cloud_nml")
#else
  if( FILE_EXIST( 'input.nml' ) ) then
! -------------------------------------
    unit = open_namelist_file ('input.nml')
    ierrnml = 1
    do while( ierrnml .ne. 0 )
      READ ( unit,  nml = diag_cloud_nml, iostat = io, end = 10 ) 
      ierrnml = check_nml_error(io,'diag_cloud_nml')
    end do
10  call close_file (unit)
! -------------------------------------
  end if
#endif

!---------------------------------------------------------------------
! --- Output namelist
!---------------------------------------------------------------------

  logunit = stdlog()
  if ( mpp_pe() == mpp_root_pe() ) then
    call write_version_number(version, tagname)
    write (logunit, nml=diag_cloud_nml)
  endif     

!---------------------------------------------------------------------
! --- Allocate storage for global cloud quantities
!---------------------------------------------------------------------


  allocate( temp_sum(ix,iy,kx),qmix_sum(ix,iy,kx),rhum_sum(ix,iy,kx) )
  allocate( qmix_sum2(ix,iy) )
  allocate( omega_sum(ix,iy,kx),lgscldelq_sum(ix,iy,kx),cnvcntq_sum(ix,iy,kx) )
  allocate( convprc_sum(ix,iy),nsum(ix,iy), nsum2(ix,iy) )

! need to set up to account for first radiation step without having
! diag cloud info available (radiation called before diag_cloud, and
! diag_cloud being initiated (cold-started) in this job). 

      tot_pts = ix*iy
!---------------------------------------------------------------------
!---------- initialize for cloud averaging -------------------------
!---------------------------------------------------------------------

  fname = 'diag_cloud.res.nc'
  id_restart = register_restart_field(Dia_restart, fname, 'nsum', nsum, no_domain=.true.)
  id_restart = register_restart_field(Dia_restart, fname, 'temp_sum', temp_sum, no_domain=.true.)
  id_restart = register_restart_field(Dia_restart, fname, 'qmix_sum', qmix_sum, no_domain=.true.)
  id_restart = register_restart_field(Dia_restart, fname, 'rhum_sum', rhum_sum, no_domain=.true.)
  id_restart = register_restart_field(Dia_restart, fname, 'omega_sum', omega_sum, no_domain=.true.)
  id_restart = register_restart_field(Dia_restart, fname, 'lgscldelq_sum', lgscldelq_sum, no_domain=.true.)
  id_restart = register_restart_field(Dia_restart, fname, 'cnvcntq_sum', cnvcntq_sum, no_domain=.true.)
  id_restart = register_restart_field(Dia_restart, fname, 'convprc_sum', convprc_sum, no_domain=.true.)

  if( FILE_EXIST( 'INPUT/diag_cloud.res.nc' ) ) then
     if(mpp_pe() == mpp_root_pe() ) call error_mesg ('diag_cloud_mod', &
          'Reading netCDF formatted restart file: INPUT/diag_cloud.res.nc', NOTE)
     call restore_state(Dia_restart)
     nsum2 = nsum
     qmix_sum2(:,:) = qmix_sum(:,:,size(qmix_sum,3))
     ierr = 0
     num_pts = tot_pts
  else if( FILE_EXIST( 'INPUT/diag_cloud.res' ) ) then
      call error_mesg ( 'diag_cloud_mod', 'Native restart capability has been removed.', &
                                         FATAL)
  else

      ierr = 1
      if (mpp_pe() == mpp_root_pe() ) write (logunit,12)
  12  format ('*** WARNING *** No cloud_tg restart file found ***  ' )

      nsum = 0
      nsum2 = 0
      temp_sum = 0.0
      qmix_sum = 0.0
      qmix_sum2 = 0.0
      rhum_sum = 0.0
      omega_sum = 0.0
      lgscldelq_sum = 0.0
      cnvcntq_sum = 0.0
      convprc_sum = 0.0

      num_pts = 0

  end if





!-------------------------------------------------------------------
! initialize zonal cloud routine for climatological zonal mean cloud info
! Passing the number 5 as the argument to cloud_zonal_init initializes the 
! clim zonal clouds allowing seasonal variation

       call cloud_zonal_init (5)

! initialize shallow convection for shallow convective clouds
       call shallow_conv_init( kx )

! Initialize cloud optical and radiative properties scheme 
      if (.not. do_crad_init) then
         call diag_cloud_rad_init (do_crad_init)
      endif
  do_cpred_init = .true.
  module_is_initialized = .true.

 
!=====================================================================
  end SUBROUTINE DIAG_CLOUD_INIT

!#######################################################################

  SUBROUTINE DIAG_CLOUD_END

!=======================================================================
! local
!=======================================================================

  if (mpp_pe() == mpp_root_pe()) then
     call error_mesg ('diag_cloud_mod', 'Writing netCDF formatted restart file: RESTART/diag_cloud.res.nc', NOTE)
  endif
  call diag_cloud_restart

  module_is_initialized = .false.
 
!=====================================================================
  end SUBROUTINE DIAG_CLOUD_END

!#######################################################################
! <SUBROUTINE NAME="diag_cloud_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine diag_cloud_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp

  call save_restart(Dia_restart, timestamp)

end subroutine diag_cloud_restart
! </SUBROUTINE> NAME="diag_cloud_restart"

!#######################################################################

 function do_diag_cloud ( ) result (answer)
   logical :: answer

!  returns logical value for whether diag_cloud has been initialized
!  presumably if initialized then diag_cloud will be used

   answer = do_cpred_init

 end function do_diag_cloud

!#######################################################################

 SUBROUTINE DIAG_CLOUD_SUM (is,js, &
                    temp,qmix,rhum,omega,lgscldelq,cnvcntq,convprc,kbot)

!-----------------------------------------------------------------------
 integer, intent(in)                 :: is,js
 real, intent(in), dimension (:,:,:) ::  temp,qmix,rhum,omega
 real, intent(in), dimension (:,:,:) ::  lgscldelq,cnvcntq
 real, intent(in), dimension (:,:)   ::  convprc

 integer, intent(in), OPTIONAL, dimension(:,:) :: kbot

!      INPUT
!      ------

!      IS,JS    starting i,j indices from the full horizontal grid
!      TEMP     Temperature (Deg K) at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      QMIX     Mixing Ratio at full model levels 
!                   (dimensioned IDIM x JDIM x kx)
!      RHUM     Relative humidity fraction at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      OMEGA  Pressure vertical velocity at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      LGSCLDELQ  Averaged rate of change in mix ratio due to lg scale precip 
!               at full model levels  
!               (dimensioned IDIM x JDIM x kx)
!      CNVCNTQ  Accumulated count of change in mix ratio due to conv precip 
!               at full model levels  
!               (dimensioned IDIM x JDIM x kx)
!      convprc Accumulated conv precip rate summed over all
!               full model levels (mm/day )
!               (dimensioned IDIM x JDIM)
!      KBOT      OPTIONAL; lowest model level index array
!                   (dimensioned IDIM x JDIM)
! ******* kbot will be used to select only those qmix values that are really
! ******* needed (typically this will be the bottom level except for 
! ******* step mountains
!-----------------------------------------------------------------------
   integer :: ie, je

   ie = is + size(rhum,1) - 1
   je = js + size(rhum,2) - 1

!--------- use time-averaged or instantaneous clouds -----------

   if (do_average) then
      nsum(is:ie,js:je)   =  nsum(is:ie,js:je)   +  1
      nsum2(is:ie,js:je)   =  nsum2(is:ie,js:je)   +  1
      temp_sum(is:ie,js:je,:) = temp_sum(is:ie,js:je,:) + temp(:,:,:)
      qmix_sum(is:ie,js:je,:) = qmix_sum(is:ie,js:je,:) + qmix(:,:,:)
      qmix_sum2(is:ie,js:je) = qmix_sum2(is:ie,js:je) + qmix(:,:,size(qmix,3))
      rhum_sum(is:ie,js:je,:) = rhum_sum(is:ie,js:je,:) + rhum(:,:,:)
      omega_sum(is:ie,js:je,:) = omega_sum(is:ie,js:je,:) + omega(:,:,:)
      lgscldelq_sum(is:ie,js:je,:) = lgscldelq_sum(is:ie,js:je,:) &
                                   + lgscldelq(:,:,:)
      cnvcntq_sum(is:ie,js:je,:) = cnvcntq_sum(is:ie,js:je,:) + cnvcntq(:,:,:)
      convprc_sum(is:ie,js:je) = convprc_sum(is:ie,js:je) + convprc(:,:)
   else
      nsum(is:ie,js:je)   =  1
      nsum2(is:ie,js:je)   =  1
      temp_sum(is:ie,js:je,:) = temp(:,:,:)
      qmix_sum(is:ie,js:je,:) = qmix(:,:,:)
      qmix_sum2(is:ie,js:je) = qmix(:,:,size(qmix,3))
      rhum_sum(is:ie,js:je,:) = rhum(:,:,:)
      omega_sum(is:ie,js:je,:) = omega(:,:,:)
      lgscldelq_sum(is:ie,js:je,:) = lgscldelq(:,:,:)
      cnvcntq_sum(is:ie,js:je,:) = cnvcntq(:,:,:)
      convprc_sum(is:ie,js:je) = convprc(:,:)
   endif

!-----------------------------------------------------------------------

 end SUBROUTINE DIAG_CLOUD_SUM

!#######################################################################

 subroutine DIAG_CLOUD_AVG (is, js, temp,qmix,rhum,omega, &
                           lgscldelq,cnvcntq,convprc,      ierr)

!-----------------------------------------------------------------------
   integer, intent(in)                    :: is, js
      real, intent(inout), dimension(:,:,:) :: temp,qmix,rhum,omega
      real, intent(inout), dimension(:,:,:) :: lgscldelq,cnvcntq
      real, intent(inout), dimension(:,:)   :: convprc
   integer, intent(out)                   :: ierr
!-----------------------------------------------------------------------
   integer ::ie, je, num, k
!-----------------------------------------------------------------------

   if (size(rhum,3) .ne. size(rhum_sum,3)) call error_mesg ( &
                               'diag_cloud_avg in diag_cloud_mod',  &
                               'input argument has the wrong size',2)

   ie = is + size(rhum,1) - 1
   je = js + size(rhum,2) - 1
   num = count(nsum(is:ie,js:je) == 0)

   if (num > 0) then

!     ----- no average, return error flag -----

!!!    call error_mesg ('diag_cloud_avg in diag_cloud_mod',  &
!!!                     'dividing by a zero counter', 2)
       ierr = 1

   else

!      ----- compute average -----

       do k = 1, size(rhum,3)
          temp(:,:,k) = temp_sum(is:ie,js:je,k) / float(nsum(is:ie,js:je))
          qmix(:,:,k) = qmix_sum(is:ie,js:je,k) / float(nsum(is:ie,js:je))
          rhum(:,:,k) = rhum_sum(is:ie,js:je,k) / float(nsum(is:ie,js:je))
          omega(:,:,k) = omega_sum(is:ie,js:je,k) / float(nsum(is:ie,js:je))
          lgscldelq(:,:,k) = lgscldelq_sum(is:ie,js:je,k) / &
                             float(nsum(is:ie,js:je))
       enddo
          convprc(:,:) = convprc_sum(is:ie,js:je) / &
                             float(nsum(is:ie,js:je))

! The convective delta qmix count should be a sum, so no average is taken
          cnvcntq(:,:,:) = cnvcntq_sum(is:ie,js:je,:) 

       ierr = 0

   endif

    nsum(is:ie,js:je)   = 0
   temp_sum(is:ie,js:je,:) = 0.0
   qmix_sum(is:ie,js:je,:) = 0.0
   rhum_sum(is:ie,js:je,:) = 0.0
   omega_sum(is:ie,js:je,:) = 0.0
   lgscldelq_sum(is:ie,js:je,:) = 0.0
   cnvcntq_sum(is:ie,js:je,:) = 0.0
   convprc_sum(is:ie,js:je) = 0.0
     
!-----------------------------------------------------------------------

 end SUBROUTINE DIAG_CLOUD_AVG

!#######################################################################

 subroutine DIAG_CLOUD_AVG2 (is, js, qmix, ierr)

!-----------------------------------------------------------------------
   integer, intent(in)                    :: is, js
      real, intent(inout), dimension(:,:) :: qmix
   integer, intent(out)                   :: ierr
!-----------------------------------------------------------------------
   integer ::ie, je, num
!-----------------------------------------------------------------------

!  if (size(qmix,3) .ne. size(qmix_sum2,3)) call error_mesg ( &
!                              'diag_cloud_avg in diag_cloud_mod',  &
!                              'input argument has the wrong size',2)

   ie = is + size(qmix,1) - 1
   je = js + size(qmix,2) - 1
   num = count(nsum2(is:ie,js:je) == 0)

   if (num > 0) then

!     ----- no average, return error flag -----

!!!    call error_mesg ('diag_cloud_avg in diag_cloud_mod',  &
!!!                     'dividing by a zero counter', 2)
       ierr = 1

   else

!      ----- compute average -----

          qmix(:,:) = qmix_sum2(is:ie,js:je) / float(nsum2(is:ie,js:je))

       ierr = 0

   endif

    nsum2(is:ie,js:je)   = 0
   qmix_sum2(is:ie,js:je) = 0.0
     
!-----------------------------------------------------------------------

 end SUBROUTINE DIAG_CLOUD_AVG2

!#######################################################################


end MODULE DIAG_CLOUD_MOD

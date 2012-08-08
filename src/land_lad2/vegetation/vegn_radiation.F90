module vegn_radiation_mod

use fms_mod,            only : write_version_number, error_mesg, FATAL
use constants_mod,      only : stefan

use land_constants_mod, only : NBANDS
use vegn_data_mod,      only : spdata, min_cosz
use vegn_tile_mod,      only : vegn_tile_type
use vegn_cohort_mod,    only : vegn_cohort_type, vegn_data_cover, get_vegn_wet_frac
use snow_mod,           only : snow_radiation

use land_debug_mod,     only : is_watch_point 

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_radiation_init
public :: vegn_radiation
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), private, parameter :: &
   version = '$Id: vegn_radiation.F90,v 17.0 2009/07/21 03:03:28 fms Exp $', &
   tagname = '$Name: siena_201207 $' ,&
   module_name = 'vegn_radiation'
! values for internal vegetation radiation option selector
integer, parameter :: VEGN_RAD_BIGLEAF   = 1 ! "big-leaf" radiation
integer, parameter :: VEGN_RAD_TWOSTREAM = 2 ! two-stream radiation code

! values for internal intercepted snow radiation properties selector -- currently
! works only for VEGN_RAD_TWOSTREAM
integer, parameter :: SNOW_RAD_IGNORE       = 1 ! no influence of intercepted snow
integer, parameter :: SNOW_RAD_PAINT_LEAVES = 2 ! intecepted snow modifies leaf 
   ! reflectance and transmittance

! ==== module variables ======================================================
integer :: vegn_rad_option = -1 ! selector of the current vegetation radiation option
integer :: snow_rad_option = -1 ! selector of the current snow rad properties option


contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


#define __DEBUG__(x) write(*,*) #x , x


! ============================================================================
! initialize vegetation radiation options
subroutine vegn_radiation_init(rad_to_use,snow_rad_to_use)
  character(*), intent(in) :: rad_to_use
  character(*), intent(in) :: snow_rad_to_use

  call write_version_number(version, tagname)

  ! convert symbolic names of radiation options into numeric IDs to
  ! speed up selection during run-time
  if(trim(rad_to_use)=='big-leaf') then
     vegn_rad_option = VEGN_RAD_BIGLEAF
  else if(trim(rad_to_use)=='two-stream') then
     vegn_rad_option = VEGN_RAD_TWOSTREAM
  else
     call error_mesg('vegn_radiation_init',&
          'vegetation radiation option rad_to_use="'//trim(rad_to_use)//'" is invalid, '// &
          'use "big-leaf" or "two-stream"',&
          FATAL)
  endif
  if (trim(snow_rad_to_use)=='ignore') then
     snow_rad_option = SNOW_RAD_IGNORE
  else if (trim(snow_rad_to_use)=='paint-leaves') then
     snow_rad_option = SNOW_RAD_PAINT_LEAVES
  else
     call error_mesg('vegn_radiation_init',&
          'vegetation radiation option snow_rad_to_use="'//trim(snow_rad_to_use)//'" is invalid, '// &
          'use "ignore" or "paint-leaves"',&
          FATAL)
  endif

end subroutine vegn_radiation_init


! ============================================================================
! compute vegetation-only properties needed to do soil-canopy-atmos 
! energy balance. in this version, also derive the aerodynamic soil-canopy
! variables, but probably this should be separated, like radiation is.
subroutine vegn_radiation ( vegn, &
     cosz, snow_depth, snow_refl_dif, snow_emis, &
     vegn_refl_dif, vegn_tran_dif, vegn_refl_dir, vegn_sctr_dir, vegn_tran_dir, &
     vegn_refl_lw, vegn_tran_lw )

  type(vegn_tile_type), intent(inout) :: vegn ! it is only inout because 
                                ! vegn_data_cover modifies cohort; cane we avoid this?
  real, intent(in)  :: cosz     ! cosine of zenith angle of direct (solar) beam
  real, intent(in)  :: snow_depth
  real, intent(in)  :: snow_refl_dif(NBANDS)
  real, intent(in)  :: snow_emis
  real, intent(out) :: &
       vegn_refl_dif(NBANDS), & ! reflectance of canopy for diffuse light
       vegn_tran_dif(NBANDS), & ! transmittance of canopy for diffuse light
       vegn_refl_dir(NBANDS), & ! reflectance of canopy for direct light
       vegn_sctr_dir(NBANDS), & ! part of direct light that is scattered downward
       vegn_tran_dir(NBANDS), & ! part of direct light that passes through the canopy unmolested
       vegn_refl_lw, & ! reflectance of canopy for long-wave (thermal) radiation 
       vegn_tran_lw    ! transmittance of canopy for long-wave (thermal) radiation

  ! ---- local vars 
  real :: vegn_cover, vegn_cover_snow_factor, vegn_lai, &
       vegn_leaf_refl(NBANDS), vegn_leaf_emis, vegn_K

  call vegn_data_cover ( vegn%cohorts(1), snow_depth, vegn_cover, vegn_cover_snow_factor )
  select case(vegn_rad_option)
  case(VEGN_RAD_BIGLEAF)
     call vegn_rad_properties_bigleaf ( vegn%cohorts(1), snow_refl_dif, snow_emis, &
          vegn_leaf_refl, vegn_leaf_emis, vegn_lai, vegn_K )
     vegn_refl_dif = vegn_cover * vegn_leaf_refl
     vegn_refl_dir = vegn_cover * vegn_leaf_refl
     vegn_tran_dif = 1 - vegn_cover 
     vegn_sctr_dir = 0
     vegn_tran_dir = 1 - vegn_cover
  case(VEGN_RAD_TWOSTREAM)
     call vegn_rad_properties_twostream ( vegn%cohorts(1), cosz, &
          vegn_refl_dif, vegn_tran_dif, &
          vegn_refl_dir, vegn_sctr_dir, vegn_tran_dir,&
          vegn_leaf_emis )
!
! ++++ pcm
     vegn_refl_dif = vegn_cover_snow_factor * vegn_refl_dif
     vegn_refl_dir = vegn_cover_snow_factor * vegn_refl_dir
     vegn_sctr_dir = vegn_cover_snow_factor * vegn_sctr_dir
     vegn_tran_dif = vegn_cover_snow_factor * vegn_tran_dif &
                       + (1-vegn_cover_snow_factor)
     vegn_tran_dir = vegn_cover_snow_factor * vegn_tran_dir &
                       + (1-vegn_cover_snow_factor)
! ---- pcm
!
  case default
     call error_mesg('vegn_radiation', &
          'invalid vegetation radiation option', FATAL)
  end select
  vegn_refl_lw       = vegn_cover * (1-vegn_leaf_emis)
  vegn_tran_lw       = 1 - vegn_cover

  ! store the extinction coefficients for use in photosynthesis calculations -- 
  ! currently calculated as if all light were direct
  vegn%cohorts(1)%extinct = &
       (spdata(vegn%cohorts(1)%species)%phi1+spdata(vegn%cohorts(1)%species)%phi2*cosz)&
       / max(cosz,min_cosz)

end subroutine vegn_radiation

! ============================================================================
! compute vegetation-only properties needed to do soil-canopy-atmos 
! energy balance.
subroutine vegn_rad_properties_bigleaf ( cohort, snow_refl, snow_emis, &
     vegn_leaf_refl, vegn_leaf_emis, vegn_lai, vegn_K )
  type(vegn_cohort_type), intent(in) :: cohort
  real, intent(in)  :: snow_refl(NBANDS), snow_emis
  real, intent(out) :: vegn_leaf_refl(NBANDS), vegn_leaf_emis, vegn_lai, vegn_K

  ! ---- local vars 
  real :: a_vs

  if ( cohort%Ws_max > 0 ) then
     a_vs = cohort%prog%Ws / cohort%Ws_max
     ! restrict snow-covered fraction to the interval [0,1]:
     a_vs = min(max(a_vs,0.0), 1.0)
  else
     a_vs = 0
  endif

  ! ---- snow-interception-adjusted radiative properties of vegetation ---------
  vegn_leaf_refl = cohort%leaf_refl + (snow_refl - cohort%leaf_refl)*a_vs
  vegn_leaf_emis = cohort%leaf_emis + (snow_emis - cohort%leaf_emis)*a_vs

  vegn_K = 2.  ! this is a temporary placeholder for now. value does not matter
               ! as long as substrate albedoes for dif/dir are the same.
  vegn_lai      = cohort%lai

end subroutine vegn_rad_properties_bigleaf


! ============================================================================
subroutine vegn_rad_properties_twostream( cohort, cosz, &
          vegn_refl_dif, vegn_tran_dif, &
          vegn_refl_dir, vegn_sctr_dir, vegn_tran_dir, &
          vegn_leaf_emis )
  type(vegn_cohort_type), intent(in) :: cohort
  real, intent(in) :: cosz ! cosine of direct light zenith angle
  real, intent(out), dimension(NBANDS) :: &
       vegn_refl_dif, & ! reflectance for diffuse light
       vegn_tran_dif, & ! transmittance for diffuse light 
       vegn_refl_dir, & ! reflectance for direct light
       vegn_sctr_dir, & ! part of direct light scattered downward (source of 
                        ! diffuse due to direct light scattering)
       vegn_tran_dir    ! transmittance of direct light 
  real, intent(out) :: &
       vegn_leaf_emis   ! emissivity of leaves

  ! ---- local constants
  real, parameter :: albedo_surf = 0.0 ! since we need values for black-background contribution

  ! ---- local vars
  integer :: i
  integer :: sp ! current species, solely to shorten the notation 
  real :: leaf_refl, leaf_tran ! optical properties of partially snow-covered leaves
  real :: snow_refl_dif(NBANDS) ! snow reflectances
  real :: snow_refl_dir(NBANDS), snow_refl_lw, snow_emis ! snow rad. properies (unused)
  real :: fs ! fractional coverage of intercepted snow

  ! get the snow fraction
  select case (snow_rad_option)
  case(SNOW_RAD_PAINT_LEAVES) 
     call get_vegn_wet_frac(cohort, fs=fs)
  case default
     fs = 0
  end select

  ! get the snow radiative properties for current canopy temperature
  call snow_radiation ( cohort%prog%Tv, cosz, snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis )

  sp = cohort%species
  do i = 1, NBANDS
     ! calculate the radiative properties of partially snow-covered leaves, assuming
     ! that the transmittance of snow is zero.
     leaf_refl = (1-fs)*cohort%leaf_refl(i) + fs*snow_refl_dif(i)
     leaf_tran = (1-fs)*cohort%leaf_tran(i)

     call twostream ( max(cosz, min_cosz), &
          spdata(sp)%mu_bar, cohort%lai, albedo_surf, &
          spdata(sp)%phi1, spdata(sp)%phi2, &
          leaf_refl, leaf_tran, &
          vegn_tran_dir(i), vegn_sctr_dir(i), vegn_refl_dir(i), &
          vegn_tran_dif(i), vegn_refl_dif(i)  )
  enddo

  vegn_leaf_emis = cohort%leaf_emis*(1-fs) + snow_emis*fs

end subroutine vegn_rad_properties_twostream


! ============================================================================
subroutine twostream( &
   mu, mu_bar, LAI, albedo_g, phi1, phi2, rl, tl, transm_dir, scatter_dir, albedo_dir, &
   transm_dif, albedo_dif )
   
  real, intent(in)  :: mu         ! cosine of direct light zenith angle
  real, intent(in)  :: mu_bar     ! average inverse diffuse optical depth per unit leaf area
  real, intent(in)  :: LAI        ! leaf area index
  real, intent(in)  :: albedo_g   ! ground surface albedo
  real, intent(in)  :: phi1, phi2 ! coefficients of expression for G_mu
  real, intent(in)  :: rl         ! reflectivity of leaves
  real, intent(in)  :: tl         ! transmittance of leaves
  ! output
  real, intent(out) :: transm_dir ! canopy transmittance for direct beam -- that 
                                  ! is, the part of the beam that passes through 
                                  ! the canopy untouched
  real, intent(out) :: scatter_dir! part of direct beam scattered down, at the 
                                  ! bottom of the canopy 
  real, intent(out) :: albedo_dir ! overall land surface albedo for direct beam 
  real, intent(out) :: transm_dif ! canopy transmittance for diffuse incident light
  real, intent(out) :: albedo_dif ! overall land surface albedo for diffuse incident light

  ! ---- local vars 
  real :: G_mu        ! relative projected leaf area in direction of direct beam
  real :: K           ! optical depth for direct beam per unit LAI
  real :: g1,g2,g3,g4 ! coefficients in the two-stream equation
  real :: kappa       ! eigenvalue of free solution
  real :: a_up, b_up, c_up ! coefficients of upward diffuse light flux
  real :: a_dn, b_dn, c_dn ! coefficients of downward diffuse light flux
  real :: x1,x2       ! intermediate coefficients
  real :: a11,a12,a21,a22, d1,d2 ! coefficients of linear system
  real :: D           ! determinant of the matrix
  real :: A,B         ! coefficients of diffuse light function
  real :: dif_dn_bot, dif_up_bot, dif_up_top
   
  real, parameter :: eta = 6.0; ! this value is suitable only for uniform leaf 
                                ! angular distribution !!!
   
  if(is_watch_point()) then
     write(*,*)'############ twostream input ############'
     __DEBUG__(mu)
     __DEBUG__(mu_bar)
     __DEBUG__(LAI)
     __DEBUG__(albedo_g)
     __DEBUG__(phi1)
     __DEBUG__(phi2)
     __DEBUG__(rl)
     __DEBUG__(tl)
     write(*,*)'############ twostream input ############'
  endif    
  ! calculate coefficients of optical path
  G_mu=phi1+phi2*mu;
  K = G_mu/mu;
 
  ! given optical parameters, calculate coefficients of basic equation
  g1 = (1-(rl+tl)/2+(rl-tl)/eta)/mu_bar;
  g2 = (  (rl+tl)/2+(rl-tl)/eta)/mu_bar;
  g3 = G_mu*((rl+tl)/2+mu*(rl-tl)/eta/G_mu);
  g4 = G_mu*((rl+tl)/2-mu*(rl-tl)/eta/G_mu);
   
  ! calculate eigenvalue of free solution (=exponent coefficient of 
  ! free solution, notes 12)
  kappa = sqrt(g1**2-g2**2);
   
  ! calculate forced term coefficients for diffuse light intensity
  c_up = ( K*g3-g1*g3-g2*g4)/(K*K-g1*g1+g2*g2);
  c_dn = (-K*g4-g1*g4-g2*g3)/(K*K-g1*g1+g2*g2);
  ! calculate intermediate coefficients for solution
  x1 = g1+g2+kappa; x2 = g1+g2-kappa;

  !write(*,*)mu,K,g1,g2,g3,g4,c_up,c_dn
   
  ! calculate coefficients of the matrix
  a11 = x2; 
  a12 = x1; 
  d1  = -c_dn;
  a21 = exp(kappa*LAI)*(x1-albedo_g*x2);
  a22 = exp(-kappa*LAI)*(x2-albedo_g*x1);
  d2  = exp(-K*LAI)*(albedo_g*c_dn + albedo_g*mu - c_up);
  ! solve the equation system 
  D = a11*a22-a12*a21;
  A = (d1*a22-d2*a12)/D;
  B = (a11*d2-a21*d1)/D;
  
  ! calculate coefficients of the diffuse light solution
  a_up = A*x1; b_up=B*x2;
  a_dn = A*x2; b_dn=B*x1;
  
  ! calculate downward diffuse light at the bottom of the canopy
  dif_dn_bot = a_dn*exp(kappa*LAI) + b_dn*exp(-kappa*LAI)+c_dn*exp(-K*LAI);
  dif_up_bot = a_up*exp(kappa*LAI) + b_up*exp(-kappa*LAI)+c_up*exp(-K*LAI);
  
  ! calculate canopy transmittance and scattered part for direct light
  scatter_dir = dif_dn_bot/mu; ! scatter
  transm_dir  = exp(-K*LAI); ! transmittance for unmolested direct light

  ! calculate canopy reflectance for direct light
  dif_up_top = a_up + b_up + c_up;
  albedo_dir = dif_up_top/mu;
     
  ! calculate upward diffuse light at the top of the canopy
  ! no need to recalculate D, since the matrix is the same
  d1 = 1.0;
  d2 = 0.0;
  A = (d1*a22-d2*a12)/D;
  B = (a11*d2-a21*d1)/D;

  ! calculate coefficients of the diffuse light solution
  a_up = A*x1; b_up=B*x2;
  a_dn = A*x2; b_dn=B*x1;
  
  transm_dif = a_dn*exp(kappa*LAI) + b_dn*exp(-kappa*LAI);
  albedo_dif = a_up + b_up;

  if(is_watch_point()) then
     write(*,*)'############ twostream output #############'
     __DEBUG__(transm_dir)
     __DEBUG__(scatter_dir)
     __DEBUG__(albedo_dir)
     __DEBUG__(transm_dif)
     __DEBUG__(albedo_dif)
     write(*,*)'############ end of twostream output #############'
  endif
end subroutine twostream

end module vegn_radiation_mod

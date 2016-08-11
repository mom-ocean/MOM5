MODULE rotstayn_klein_mp_mod

use mpp_mod,             only : input_nml_file
use fms_mod,             only : error_mesg, FATAL, mpp_pe, mpp_root_pe, &
                                open_namelist_file, check_nml_error, &
                                close_file, write_version_number, &
                                file_exist, stdlog
use cloud_generator_mod, only : cloud_generator_init, do_cloud_generator, &
                                compute_overlap_weighting
use  constants_mod,      only : hlv,hlf,hls, rdgas, cp_air, grav, &
                                tfreeze, dens_h2o, rvgas
use  aer_in_act_mod,     only : Jhete_dep
use  sat_vapor_pres_mod, only : compute_qs,  sat_vapor_pres_init 
use polysvp_mod,         only : polysvp_init, polysvp_end, polysvp_l,  &
                                polysvp_i
use strat_cloud_utilities_mod,     &
                         only : strat_cloud_utilities_init, &
                                diag_id_type, diag_pt_type, strat_nml_type

implicit none
private 

!-------------------------------------------------------------------------
!---interfaces------------------------------------------------------------

public   rotstayn_klein_microp, rotstayn_klein_microp_init, &
         rotstayn_klein_microp_end
private  cloud_clear_xfer

!-------------------------------------------------------------------------
!---version number-------------------------------------------------------

Character(len=128) :: Version = '$Id: rotstayn_klein_mp.F90,v 20.0 2013/12/13 23:22:05 fms Exp $'
Character(len=128) :: Tagname = '$Name: tikal $'

!-------------------------------------------------------------------------
!---namelist-------------------------------------------------------------
 
logical :: rk_alt_adj_opt = .false.             ! we should use an altern-
                                                ! ative method to determine
                                                ! supersaturation  using
                                                ! formula used in mg_micro?
logical :: rk_act_only_if_ql_gt_qmin = .false.  ! compute qn source term
                                                ! (activated nuclei)
                                                ! only when ql is > qmin ?
logical :: use_inconsistent_lh = .true.         ! use the legacy latent
                                                ! heat treatment for the 
                                                ! condensation of super-
                                                ! saturated vapor ?
                                                ! (this is non-entropy-
                                                ! conserving and physically
                                                ! unrealistic ??)

namelist / rotstayn_klein_mp_nml /   rk_alt_adj_opt, &
                                     rk_act_only_if_ql_gt_qmin, &
                                     use_inconsistent_lh


!-------------------------------------------------------------------------
!---local module variables and parameters---------------------------------

real, parameter :: d622 = rdgas / rvgas
real, parameter :: d378 = 1. - d622
real, parameter :: rho_ice        =  100.  ! mass density of ice crystals 
                                           ! [ kg/(m*m*m) ]
real, parameter :: ELI            =  0.7   ! collection efficiency of 
                                           ! cloud liquid by falling ice
                                           ! [ dimensionless ]
logical         :: cloud_generator_on      ! stochastic clouds are active?

!----------------------------------------------------------------------- 
logical         :: module_is_initialized = .false.


!----------------------------------------------------------------------


CONTAINS


!########################################################################

SUBROUTINE rotstayn_klein_microp_init

!-----------------------------------------------------------------------
      integer :: unit, ierr, io, logunit

!------------------------------------------------------------------------
!    if module has already been initialized, return.
!------------------------------------------------------------------------
      if (module_is_initialized) return

!-------------------------------------------------------------------------
!    process namelist.
!-------------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=rotstayn_klein_mp_nml, iostat=io)
      ierr = check_nml_error(io,'rotstayn_klein_mp_nml')
#else
      if ( file_exist('input.nml')) then
        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=rotstayn_klein_mp_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'rotstayn_klein_mp_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!-------------------------------------------------------------------------
!    write version and namelist to standard log.
!-------------------------------------------------------------------------
      call write_version_number(version, tagname)
      logunit = stdlog()
      if ( mpp_pe() == mpp_root_pe() ) &
                         write (logunit, nml=rotstayn_klein_mp_nml)

!------------------------------------------------------------------------
!    be sure needed modules have been initialized.
!------------------------------------------------------------------------
      call sat_vapor_pres_init
      call strat_cloud_utilities_init
      call polysvp_init
      call cloud_generator_init
      cloud_generator_on = do_cloud_generator()

!-----------------------------------------------------------------------
!    mark the module as initialized.
!-----------------------------------------------------------------------
      module_is_initialized = .true.


END SUBROUTINE rotstayn_klein_microp_init


!########################################################################

SUBROUTINE rotstayn_klein_microp (&
                     idim, jdim, kdim, Nml, N3D, overlap, dtcloud,  &
                     inv_dtcloud, pfull, deltpg, airdens, mask_present, &
                     mask, esat0, ql, qi, qa, ql_mean, qa_mean, qn_mean, &
                     omega, T, U, qv, qs, D_eros, dcond_ls, dcond_ls_ice, &
                     qvg, gamma, delta_cf, drop1, concen_dust_sub, ql_upd,&
                     qi_upd, qn_upd, qi_mean, qa_upd, ahuco, n_diag_4d, &
                     diag_4d, diag_id, diag_pt, n_diag_4d_kp1,   &
                     diag_4d_kp1, limit_conv_cloud_frac, SA, SN, ST, SQ, &
                     SL, SI, rain3d, snow3d, snowclr3d, surfrain,   &
                     surfsnow, f_snow_berg, otun)                 

!------------------------------------------------------------------------
type(strat_nml_type),                intent(in)   :: Nml
LOGICAL,                             INTENT(IN)   :: mask_present,   &
                                                     limit_conv_cloud_frac
INTEGER,                             INTENT(IN)   :: overlap
INTEGER,                             INTENT(IN)   :: idim, jdim, kdim
INTEGER,                             intent(in)   :: n_diag_4d,  &
                                                     n_diag_4d_kp1
REAL,                                INTENT(IN)   :: dtcloud, inv_dtcloud
REAL, dimension(idim, jdim,kdim),    INTENT(IN)   :: MASK
REAL, dimension(idim,jdim,kdim),     INTENT(IN)   ::   &
                              N3D, pfull, deltpg, airdens, ql, qi, qa, &
                              ql_mean, qn_mean, qa_mean, T, qv, drop1, U, &
                              qvg, esat0, gamma, concen_dust_sub, D_eros, &
                              omega, ahuco
real, dimension(idim,jdim,kdim),     INTENT(INOUT)::   &
                              SN, qi_mean, ST, SQ, SL, SI, SA, qa_upd, &
                              ql_upd, qi_upd, qn_upd, qs, delta_cf, &
                              dcond_ls,dcond_ls_ice    
REAL, dimension(idim,jdim,kdim,0:n_diag_4d),   &
                                     INTENT(INOUT)::  diag_4d
REAL, dimension( idim,jdim,kdim+1,0:n_diag_4d_kp1),   &
                                     INTENT(INOUT)::  diag_4d_kp1
TYPE(diag_id_type),                  intent(in)   :: diag_id
TYPE(diag_pt_type),                  intent(in)   :: diag_pt
real, dimension(idim, jdim, kdim+1), INTENT(OUT)  :: rain3d, snow3d
real, dimension(idim, jdim, kdim+1), INTENT(OUT)  :: snowclr3d
real, dimension(idim,jdim),          INTENT(OUT)  :: surfrain,surfsnow
real, dimension(idim, jdim, kdim  ), INTENT(OUT)  :: f_snow_berg 
INTEGER,                             INTENT(IN)   :: otun
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!---local variables-----------------------------------------------------
!
!       rain_cld       grid mean flux of rain enter-   kg condensate/
!                      ing the grid box from above     (m*m)/s
!                      and entering the saturated 
!                      portion of the grid box
!
!
!       rain_clr       grid mean flux of rain enter-   kg condensate/
!                      ing the grid box from above     (m*m)/s
!                      and entering the unsaturated 
!                      portion of the grid box
!
!       a_rain_clr     fraction of grid box occupied   fraction
!                      by rain_clr
!
!       a_rain_cld     fraction of grid box occupied   fraction
!                      by rain_cld
!
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
!       da_cld2clr     fraction of the area in which   fraction
!                      rain/snow in saturated volume 
!                      above falls into unsaturated 
!                      volume in the current layer.
!
!       da_clr2cld     as in da_cld2clr except for     fraction
!                      the transfer from unsaturated
!                      to saturated volume
!                      non-convective condensation.    kg air
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
!       Vfall          fall speed of ice crystals      m/s
!
!       lamda_f        slope factor in the SIZE        1/m
!                      distribution of ice crystals
!       rad_liq        mean volume radius of liquid    microns
!                      cloud drops
!
!       A_plus_B       sum of vapor diffusion factor   m*s/kg
!                      and thermal conductivity factor
!                      which is used in various 
!                      microphysical formula for the 
!                      evaporation of rain and snow
!
!       U_clr          relative humidity in the clear  fraction
!                      portion of the grid box.
!
!------------------------------------------------------------------------

      real, dimension(idim,jdim,kdim+1)  :: rain_clr, rain_cld, &
                                            a_rain_clr, a_rain_cld, &
                                            snow_clr, snow_cld, &
                                            a_snow_clr, a_snow_cld
      real, dimension(idim,jdim,kdim)    :: A_plus_B, C_dt, D_dt, &
                                            da_cld2clr, da_clr2cld, &
                                            dprec_clr2cld, dprec_cld2clr, &
                                            Vfall, lamda_f, tmp1, tmp2, &
                                            tmp3, tmp8, crystal,          &
                                            rad_liq, D1_dt, D2_dt, qc1, &
                                            qc0, qceq, qcbar, U_clr, &
                                            qs_d_a, tmp2s_a, tmp3s_a, &
                                            tmp5s_a, est_a, sum_freeze, &
                                            sum_cond,  sum_ice_adj, &
                                            sum_rime, sum_berg, tmp5
      real                               :: dum, Si0, qs_t, qs_d, tmp2s, &
                                            tmp3s, tmp5s, est, rhi, tc, &
                                            tcrit, qldt_sum
      integer                            :: i, j, k, km1

 
!------------------------------------------------------------------------  
!    initialize top levels of precip and precip area arrays.
!------------------------------------------------------------------------  
      a_rain_cld(:,:,1) = 0.
      a_rain_clr(:,:,1) = 0.
      rain_cld(:,:,1)   = 0.
      rain_clr(:,:,1)   = 0.
      snow_cld(:,:,1)   = 0.
      snow_clr(:,:,1)   = 0.
      a_snow_clr(:,:,1) = 0.
      a_snow_cld(:,:,1) = 0.

!-----------------------------------------------------------------------
!    compute A_plus_B which is the sum of vapor diffusion factor and 
!    thermal conductivity factor which is used in various microphysical
!    formula for the evaporation of rain and snow [m s / kg ].
!    The conductivity/diffusivity factor, A_plus_B is given by:
!
!    (4) A_plus_B =   { (hlv/Ka/T)*((hlv/rvgas/T)-1.) } + 
!
!                    { (rvgas*T/chi*esat) }
!
!    where Ka is the thermal conductivity of air = 0.024 J/m/s/K
!    and chi is the diffusitivy of water vapor in air which is
!    given by
!
!    (5) chi = 2.21 E-05 (m*m)/s  * (1.E+05)/pfull
!
!    where p is the pressure in Pascals.
!-----------------------------------------------------------------------
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            A_plus_B(i,j,k) = ((hlv/0.024/T(i,j,k))*  &
                              ((hlv/rvgas/T(i,j,k)) - 1.) ) +        &
                              (rvgas*T(i,j,k)*pfull(i,j,k)/2.21/  &
                                                           esat0(i,j,k))  
          end do
        end do
      end do

      sum_freeze = 0.
      sum_rime = 0.
      sum_berg = 0.
      sum_ice_adj = 0.
      sum_cond = 0.

!------------------------------------------------------------------------
!    begin big vertical loop.
!------------------------------------------------------------------------
      do k=1, kdim  
        km1 = MAX(k-1, 1)

!------------------------------------------------------------------------
!    compute weighting factor for overlap when stochastic clouds are
!    activated.
!------------------------------------------------------------------------
        if (cloud_generator_on) then
          if (k .GT. 1) THEN 
            tmp3(:,:,k) = compute_overlap_weighting (   &
                                  qa_mean(:,:,km1),qa_mean(:,:,k),&
                                              pfull(:,:,k-1),pfull(:,:,k))
            tmp3(:,:,k) = min(1., max(0., tmp3(:,:,k)))       
          else
            tmp3(:,:,k) = 0.
          endif 
        endif

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
 
!------------------------------------------------------------------------
!    rain transfers are done first; compute cloud to clear transfer.
!------------------------------------------------------------------------
        call cloud_clear_xfer (k, Nml, cloud_generator_on, overlap, tmp3, &
                               qa_mean, a_rain_clr, a_rain_cld,   &
                                                      rain_clr, rain_cld)
  
!------------------------------------------------------------------------
!    snow transfers are done second, in a manner exactly like that
!    done for the rain fluxes
!------------------------------------------------------------------------
        call cloud_clear_xfer (k, Nml, cloud_generator_on, overlap, tmp3, &
                               qa_mean, a_snow_clr, a_snow_cld,  &
                                                      snow_clr, snow_cld)

!-----------------------------------------------------------------------
!
!
!                   MELTING OF CLEAR SKY SNOW FLUX
!
!
!    Melting of falling ice to rain occurs when T > tfreeze. The amount of
!    melting is limited to the melted amount that would cool the 
!    temperature to tfreeze.
!
!    In the snowmelt bug version, the temperature of melting was 
!    tfreeze + 2. like the original Tiedtke (1993) paper, instead of 
!    tfreeze.
!
!    compute grid mean change in snow flux to cool the grid box to tfreeze 
!    and store in temporary variable tmp1
!-----------------------------------------------------------------------
        if (Nml%do_old_snowmelt) then
          tmp1(:,:,k) = cp_air*(T(:,:,k) - tfreeze - 2.)*   &
                                              deltpg(:,:,k)*inv_dtcloud/hlf
        else
          tmp1(:,:,k) = cp_air*(T(:,:,k) - tfreeze)*  &
                                              deltpg(:,:,k)*inv_dtcloud/hlf
        end if
        
!------------------------------------------------------------------------
!    If snow_clr > tmp1, then the amount of snow melted is limited to tmp1,
!    otherwise melt snow_clr.  The amount melted is stored in tmp2. 
!    update temp and rain to account for the melting and conserve heat and
!    h2o.
!------------------------------------------------------------------------
        tmp2(:,:,k) = max(min(snow_clr(:,:,k), tmp1(:,:,k)), 0.)     
        ST(:,:,k) = ST(:,:,k) - hlf*tmp2(:,:,k)*dtcloud/deltpg(:,:,k)/ &
                                                                    cp_air
        rain_clr(:,:,k) = rain_clr(:,:,k) + tmp2(:,:,k)
 
!------------------------------------------------------------------------
!    change the area of clear-sky rain (a_rain_clr) to be the area of 
!    clear-sky snow (a_snow_clr) IF AND only IF melting occurs and 
!    a_rain_clr < a_snow_clr.
!------------------------------------------------------------------------
        where (tmp2(:,:,k) .gt. 0. .and. a_snow_clr(:,:,k) .gt. Nml%qmin)
          a_rain_clr(:,:,k) = max(a_rain_clr(:,:,k), a_snow_clr(:,:,k))
        end where

!------------------------------------------------------------------------
!    if all of the snow has melted, then zero out a_snow_clr
!------------------------------------------------------------------------
        where (snow_clr(:,:,k).lt.tmp1(:,:,k) .and.   &
                                        a_snow_clr(:,:,k).gt.Nml%qmin)
          snow_clr(:,:,k) = 0.
          a_snow_clr(:,:,k) = 0.
        elsewhere
          snow_clr(:,:,k) = snow_clr(:,:,k) - tmp2(:,:,k)          
        end where

!------------------------------------------------------------------------
!    define snow melt diagnostic.
!------------------------------------------------------------------------
        if (diag_id%snow_melt  + diag_id%snow_melt_col > 0)    &
            diag_4d(:,:,k,diag_pt%snow_melt) = tmp2(:,:,k)/deltpg(:,:,k)  
             
!-----------------------------------------------------------------------
!
!
!                        MELTING OF CLOUDY SKY SNOW FLUX
!
!
!    Melting of falling ice to rain occurs when T > tfreeze. The 
!    amount of melting is limited to the melted amount that would 
!    cool the temperature to tfreeze.
!-------------------------------------------------------------------------
        if (.not.Nml%do_old_snowmelt) then

!----------------------------------------------------------------------
!    compute grid mean change in snow flux to cool the grid box to tfreeze
!    and store in temporary variable tmp1
!
!    note that tmp1 already has the value of this variable from the 
!    clear-sky melt calculation, so one does not need to repeat the 
!    calculation here.
!
!    However, note that clear-sky snow melt may have already reduced the 
!    temperature of the grid box - this snow melt is in variable tmp2 from 
!    lines above. Thus the amount that one can melt is less.
!------------------------------------------------------------------------
          tmp1(:,:,k) = tmp1(:,:,k) - tmp2(:,:,k)
        
!------------------------------------------------------------------------
!    If snow_cld > tmp1, then the amount of snow melted is limited to 
!    tmp1, otherwise melt snow_cld.  The amount melted is stored in tmp2
!    update temp and rain to account for the melting and conserve heat and
!    h2o.
!------------------------------------------------------------------------
          tmp2(:,:,k) = max(min(snow_cld(:,:,k), tmp1(:,:,k)),0.)     
          ST(:,:,k) = ST(:,:,k) - hlf*tmp2(:,:,k)*dtcloud/  &
                                                     deltpg(:,:,k)/cp_air 
          rain_cld(:,:,k) = rain_cld(:,:,k) + tmp2(:,:,k)
 
!-----------------------------------------------------------------------
!    change the area of cloudy-sky rain (a_rain_cld) to be the area of 
!    cloudy-sky snow (a_snow_cld) IF AND only IF melting occurs and 
!    a_rain_cld < a_snow_cld.
!-----------------------------------------------------------------------
          where (tmp2 (:,:,k).gt. 0. .and. a_snow_cld(:,:,k) .gt. Nml%qmin)
            a_rain_cld(:,:,k) = max(a_rain_cld(:,:,k), a_snow_cld(:,:,k))
          end where

!------------------------------------------------------------------------
!    If all of the snow has melted, then zero out a_snow_cld.
!------------------------------------------------------------------------
          where (snow_cld(:,:,k).lt.tmp1(:,:,k) .and.   &
                                         a_snow_cld(:,:,k).gt.Nml%qmin)
            snow_cld(:,:,k)= 0.
            a_snow_cld(:,:,k) = 0.
          elsewhere
            snow_cld(:,:,k) = snow_cld(:,:,k) - tmp2(:,:,k)          
          end where

!------------------------------------------------------------------------
!    define snow melt diagnostic -- add to clear sky melting.
!------------------------------------------------------------------------
          if (diag_id%snow_melt  + diag_id%snow_melt_col > 0)   &
                      diag_4d(:,:,k,diag_pt%snow_melt) =   &
                                diag_4d(:,:,k,diag_pt%snow_melt) + &
                                                tmp2(:,:,k)/deltpg(:,:,k) 
        end if  ! (.not. do_old_snowmelt)  
                            
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

!       Increment qi_mean by the ice flux entering the the grid box. To 
!       convert ice_flux to units of condensate by multiply by dtcloud and
!       dividing by the mass per unit area of the grid box. Implicit here 
!       is the assumption that the ice entering the cloud will be spread 
!       instantaneously over  all of the cloudy area.
!-----------------------------------------------------------------------

!------------------------------------------------------------------------
!    update the qi variable with the amount of snow added to the box. 
!------------------------------------------------------------------------
        qi_mean(:,:,k) = qi_mean(:,:,k) + snow_cld(:,:,k)*  &
                                                    dtcloud/deltpg(:,:,k)  

!-----------------------------------------------------------------------
!    snow falling into cloud reduces the amount that falls out of cloud: 
!    a loss of cloud ice from settling is defined to be negative. define
!    the ice-fall diagnostic.
!-----------------------------------------------------------------------
        if (diag_id%qidt_fall + diag_id%qi_fall_col > 0)  &
                 diag_4d(:,:,k,diag_pt%qidt_fall) =    &
                                            snow_cld(:,:,k)/deltpg(:,:,k)
         
!-----------------------------------------------------------------------
!    compute slope factor lamda_f.
!RSH  -- conditionally calculate this only when needed
!-----------------------------------------------------------------------
        where ( ((a_snow_cld(:,:,k).gt.Nml%qmin) .and.   &
                     (ql_mean(:,:,k).gt.Nml%qmin) .and. &
                (qa_mean(:,:,k).gt.Nml%qmin) ) .or. snow_clr(:,:,k) > 0.0 )
          lamda_f(:,:,k) = 1.6*10**(3.+0.023*(tfreeze - T(:,:,k)))
        elsewhere
          lamda_f(:,:,k) = 0.
        endwhere

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
!------------------------------------------------------------------------

        if (.not. Nml%do_liq_num) then

!-----------------------------------------------------------------------
!    compute rad_liq.  The constant below is equal to   
!    1.E+06 * (3/4*pi)^(1/3), where the 1E+06 is the factor to convert 
!    meters to microns. do not let very small cloud fractions contribute 
!    to autoconversion or accretion.
!-----------------------------------------------------------------------
          where (qa_mean(:,:,k) .gt. Nml%qmin) 
            rad_liq(:,:,k)= 620350.49 *( (airdens(:,:,k)*ql_mean(:,:,k)/ &
                              max(qa_mean(:,:,k), Nml%qmin)/   &
                                          N3D(:,:,k)/dens_h2o)**(1./3.))
          elsewhere
            rad_liq(:,:,k) = 0.0
          end where
        else        
          where (qa_mean(:,:,k) .gt. Nml%qmin .and.   &
                 qn_upd(:,:,k) .gt.Nml%qmin) 
            rad_liq(:,:,k) = 620350.49*( (ql_mean(:,:,k)/  &
                             max(qn_mean(:,:,k),   &
                                            Nml%qmin)/dens_h2o)**(1./3.))
          elsewhere
            rad_liq(:,:,k) = 0.
          endwhere
        endif
        
!------------------------------------------------------------------------
!    save  particle size diagnostic.
!------------------------------------------------------------------------
        if (diag_id%rvolume > 0) diag_4d(:,:,k,diag_pt%rvolume) =   &
                              rad_liq(:,:,k) *diag_4d(:,:,k,diag_pt%aliq)

!------------------------------------------------------------------------
!    compute accretion D term from formula above.
!------------------------------------------------------------------------
        where  (rad_liq(:,:,k) > 0.0) 
          D1_dt(:,:,k) =  dtcloud*65.772565*(a_rain_cld(:,:,k)/  &
                                      max(qa_mean(:,:,k), Nml%qmin))* &
                        (rad_liq(:,:,k)*rad_liq(:,:,k)/  &
                             (rad_liq(:,:,k)*rad_liq(:,:,k) + 20.5) )*   &
                        ((rain_cld(:,:,k)/   &
                               max(a_rain_cld(:,:,k), Nml%qmin)/  &
                                                      dens_h2o)**(7./9.))
        elsewhere
          D1_dt(:,:,k) = 0.
        end where
            
!------------------------------------------------------------------------
!    save  accretion process diagnostics.
!------------------------------------------------------------------------
        if (diag_id%qldt_accr  + diag_id%ql_accr_col > 0)   &
                           diag_4d(:,:,k,diag_pt%qldt_accr) = -D1_dt(:,:,k)
        if (diag_id%qndt_pra  + diag_id%qn_pra_col > 0)    &
                            diag_4d(:,:,k,diag_pt%qndt_pra) = -D1_dt(:,:,k)
    
!------------------------------------------------------------------------
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
!------------------------------------------------------------------------
        if (Nml%do_liq_num) then
          if (Nml%use_kk_auto) then
            where (ql_mean(:,:,k) > 0.0) 
              tmp1(:,:,k) = dtcloud*1350.*(1.e-6*   &
                max(qn_mean(:,:,k)*airdens(:,:,k),  &
                      max(qa_mean(:,:,k), Nml%qmin)*Nml%N_min))**(-1.79)* &
                            (ql_mean(:,:,k))**(1.47)*  &
                                     max(qa_mean(:,:,k),Nml%qmin)**(0.32)
            elsewhere
              tmp1(:,:,k) = 0.
            end where
          else
            where (ql_mean(:,:,k) > 0.0) 
              tmp1(:,:,k) = 32681.*dtcloud*    &
                             ((max(qn_mean(:,:,k), Nml%qmin)*  &
                                    airdens(:,:,k)*dens_h2o)**(-1./3.))*  &
                                       (ql_mean(:,:,k)**(4./3.))/   &
                                             max(qa_mean(:,:,k), Nml%qmin)
!------------------------------------------------------------------------
!    compute limiter as in (35) and  limit autoconversion to the limiter
!------------------------------------------------------------------------
              tmp2(:,:,k) =max(3*log(max(rad_liq(:,:,k), Nml%qmin)/  &
                                                         Nml%rthresh), 0.)
              tmp1(:,:,k) = min(tmp1(:,:,k), tmp2(:,:,k))
            elsewhere
              tmp1(:,:,k) = 0.
            end where
          endif
        else

!------------------------------------------------------------------------
!   for non-predicted cloud number:
!------------------------------------------------------------------------
          if (Nml%use_kk_auto) then

!------------------------------------------------------------------------
!    compute autoconversion sink as in (34A)
!------------------------------------------------------------------------
            where (ql_mean(:,:,k).gt.Nml%qmin)
              tmp1(:,:,k) = 7.4188E+13*dtcloud*(N3D(:,:,k)**(-1.79))*     &
                   ((ql_mean(:,:,k)/max(qa_mean(:,:,k), Nml%qmin))**(1.47))
            elsewhere
              tmp1(:,:,k) = 0.
            endwhere
          else

!------------------------------------------------------------------------
!    compute autoconversion sink as in (34)
!------------------------------------------------------------------------
            where (ql_mean(:,:,k) > 0.0) 
              tmp1(:,:,k) = 32681.*dtcloud*((N3D(:,:,k)*   &
                              dens_h2o)**(-1./3.))*((ql_mean(:,:,k)/   &
                                  max(qa_mean(:,:,k), Nml%qmin))**(4./3.))
        
!------------------------------------------------------------------------
!    compute limiter as in (35) and  limit autoconversion to the limiter
!------------------------------------------------------------------------
              tmp2(:,:,k) = max(3*log(max(rad_liq(:,:,k), Nml%qmin)/  &
                                                       Nml%rthresh), 0.)
              tmp1(:,:,k) = min(tmp1(:,:,k), tmp2(:,:,k))
            elsewhere
              tmp1(:,:,k) = 0.
            end where
          endif
        endif

!------------------------------------------------------------------------
!    add autoconversion to D1_dt
!------------------------------------------------------------------------
        D1_dt(:,:,k) = D1_dt(:,:,k) + tmp1(:,:,k)

!------------------------------------------------------------------------
!    autoconversion will increase a_rain_cld to area of cloud
!------------------------------------------------------------------------
        where (tmp1(:,:,k) .gt. Nml%Dmin) a_rain_cld(:,:,k) =   &
                                                           qa_mean(:,:,k)

!------------------------------------------------------------------------
!    save autoconversion diagnostics.
!------------------------------------------------------------------------
        if  (diag_id%qldt_auto  + diag_id%ql_auto_col > 0)   &
                           diag_4d(:,:,k,diag_pt%qldt_auto) = -tmp1(:,:,k)  
        if  (diag_id%qndt_auto  + diag_id%qn_auto_col > 0)   &
                           diag_4d(:,:,k,diag_pt%qndt_auto) = -tmp1(:,:,k)  
        if ( diag_id%aauto > 0 ) then
          where ( rad_liq(:,:,k) .gt. Nml%rthresh )    &
                            diag_4d(:,:,k,diag_pt%aauto) = qa_mean(:,:,k)
        end if
        
!-----------------------------------------------------------------------
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
!------------------------------------------------------------------------
        if (Nml%do_dust_berg) then
          do j=1,jdim
            do i=1,idim
              if ( (T(i,j,k) .lt. tfreeze) .and.   &
                   (ql_mean(i,j,k) .gt. Nml%qmin) .and.   &
                   (qa_mean(i,j,k) .gt. Nml%qmin))  then
                Si0 = 1. + 0.0125*(tfreeze-T(i,j,k))
                call Jhete_dep (T(i,j,k), Si0, concen_dust_sub(i,j,k), &
                                                            crystal(i,j,k))
                if (diag_id%dust_berg_flag > 0)   &
                             diag_4d(i,k,j,diag_pt%dust_berg_flag) = 1.
              else
                crystal(i,j,k)=0.
              endif
            end do
          end do

!------------------------------------------------------------------------
!    save ice crystal diagnostics.
!------------------------------------------------------------------------
          if  (diag_id%qndt_cond + diag_id%qn_cond_col > 0) &
                        diag_4d(:,:,k,diag_pt%qndt_cond) = crystal(:,:,k)
          
!------------------------------------------------------------------------
!    do Bergeron process
!------------------------------------------------------------------------
          where ( (T(:,:,k) .lt. tfreeze) .and.  &
                  (ql_mean(:,:,k) .gt. Nml%qmin) .and.   &
                  (qa_mean(:,:,k) .gt. Nml%qmin))              
            D2_dt(:,:,k) = dtcloud*qa_mean(:,:,k)*((1.e6*crystal(:,:,k)/  &
                                  airdens(:,:,k))**(2./3.))* 7.8*  &
                      ((max(qi_mean(:,:,k)/qa_mean(:,:,k),   &
                                1.E-12*1.e6*crystal(:,:,k)/   &
                                       airdens(:,:,k)))**(1./3.))*  &
                            0.0125*(tfreeze - T(:,:,k))/((700.**(1./3.))* &
                                           A_plus_B(:,:,k)*ql_mean(:,:,k))
          elsewhere
            D2_dt(:,:,k) = 0.0
          end where

!------------------------------------------------------------------------
!    if dust not included in bergeron process:
!------------------------------------------------------------------------
        else
          where ( (T(:,:,k) .lt. tfreeze) .and.  &
                  (ql_mean(:,:,k) .gt. Nml%qmin) .and.   &
                  (qa_mean(:,:,k) .gt. Nml%qmin))           
            D2_dt(:,:,k) = dtcloud*qa_mean(:,:,k)*((Nml%cfact*1000.*  &
                           exp((12.96*0.0125*(tfreeze - T(:,:,k))) -   &
                               0.639)/airdens(:,:,k))**(2./3.))* 7.8* &
                               ((max(qi_mean(:,:,k)/qa_mean(:,:,k),  &
                                       1.E-12*Nml%cfact*1000.*   &
                            exp((12.96*0.0125*(tfreeze-T(:,:,k)))-0.639)/ &
                                   airdens(:,:,k)))**(1./3.))*0.0125*     &
                                (tfreeze - T(:,:,k))/((700.**(1./3.))*  &
                                          A_plus_B(:,:,k)*ql_mean(:,:,k))
          elsewhere
            D2_dt(:,:,k) = 0.0        
          end where
        endif

        sum_berg(:,:,k) = D2_dt(:,:,k)

!------------------------------------------------------------------------
!    diagnostics for bergeron process
!------------------------------------------------------------------------
        if  (diag_id%qldt_berg  + diag_id%ql_berg_col > 0)   &
                         diag_4d(:,:,k,diag_pt%qldt_berg) = - D2_dt(:,:,k)
        if  (diag_id%qndt_berg  + diag_id%qn_berg_col > 0)    &
                         diag_4d(:,:,k,diag_pt%qndt_berg) = - D2_dt(:,:,k)

!------------------------------------------------------------------------
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
!-----------------------------------------------------------------------

!------------------------------------------------------------------------
!    compute the tendency due to riming.
!------------------------------------------------------------------------
        if (Nml%do_old_snowmelt) then
          where ((a_snow_cld(:,:,k).gt.Nml%qmin) .and.    &
                 (ql_mean(:,:,k).gt.Nml%qmin) .and. &
                 (qa_mean(:,:,k).gt.Nml%qmin) .and.   &
                 (T(:,:,k) .lt. tfreeze) )            
            tmp1(:,:,k) = dtcloud*0.5*ELI*lamda_f(:,:,k)*snow_cld(:,:,k)/ &
                                                   qa_mean(:,:,k)/rho_ice 
          elsewhere
            tmp1(:,:,k) = 0.0
          end where
        else
          where ((a_snow_cld(:,:,k).gt.Nml%qmin) .and. &
                 (ql_mean(:,:,k).gt.Nml%qmin) .and. &
                 (qa_mean(:,:,k).gt.Nml%qmin) )            
            tmp1(:,:,k) = dtcloud*0.5*ELI*lamda_f(:,:,k)*  &
                                    snow_cld(:,:,k)/qa_mean(:,:,k)/rho_ice 
          elsewhere
            tmp1(:,:,k) = 0.0
          end where
        end if        
        D2_dt(:,:,k) = D2_dt(:,:,k) + tmp1(:,:,k)

        sum_rime(:,:,k) = tmp1(:,:,k)

!-----------------------------------------------------------------------
!    save the riming diagnostics.
!-----------------------------------------------------------------------
        if  (diag_id%qldt_rime  + diag_id%ql_rime_col > 0)   &
                            diag_4d(:,:,k,diag_pt%qldt_rime) = -tmp1(:,:,k)

!------------------------------------------------------------------------
!       Freezing of cloud liquid to cloud ice occurs when
!       the temperature is less than -40C. At these very cold temper-
!       atures it is assumed that homogenous freezing of cloud liquid
!       droplets will occur.   To accomplish this numerically in one 
!       time step:
!
!  (40) D*dtcloud =  ln( ql / Nml%qmin ).
!
!       With this form it is guaranteed that if this is the only
!       process acting that ql = qmin after one integration.
!------------------------------------------------------------------------

!-------------------------------------------------------------------------
!    do homogeneous freezing. save freezing diagnostics and modify
!    bergeron and riming diagnostics if homogeneous freezing has occurred.
!-------------------------------------------------------------------------
        do j=1,jdim
          do i=1,idim
            if ((T(i,j,k) .lt. tfreeze - 40.) .and.     &
                (ql_mean(i,j,k) .gt. Nml%qmin).and.     &
                (qa_mean(i,j,k) .gt. Nml%qmin))  then
              D2_dt(i,j,k) = log ( ql_mean(i,j,k)/Nml%qmin )
              sum_freeze(i,j,k) = D2_dt(i,j,k)
              sum_rime(i,j,k) = 0.
              sum_berg(i,j,k) = 0.
              if (diag_id%qldt_freez + diag_id%ql_freez_col > 0)  then
                diag_4d(i,j,k,diag_pt%qldt_freez) = -D2_dt(i,j,k)
              endif
              if  (diag_id%qldt_rime  + diag_id%ql_rime_col > 0) then
                diag_4d(i,j,k,diag_pt%qldt_rime) = 0.
              endif
              if (diag_id%qldt_berg  + diag_id%ql_berg_col > 0) then
                diag_4d(i,j,k,diag_pt%qldt_berg) = 0. 
              endif
            endif      
          end do
        end do
  
!------------------------------------------------------------------------
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
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!    C_dt is set to large-scale condensation. Sink of cloud liquid 
!    is set to the sum of D1 (liquid to rain component), and D2 
!    (liquid to ice component), D_eros (erosion), and large-scale 
!    evaporation (note use of ql mean).        
!------------------------------------------------------------------------
        C_dt(:,:,k) = max(dcond_ls(:,:,k),0.)     
        D_dt(:,:,k) = D1_dt(:,:,k) + D2_dt(:,:,k) + D_eros(:,:,k) + &
                                     (max(-1.*dcond_ls(:,:,k), 0.)/    &
                                          max(ql_mean(:,:,k), Nml%qmin)) 

!------------------------------------------------------------------------
!    do analytic integration      
!------------------------------------------------------------------------
        where ( D_dt(:,:,k).gt.Nml%Dmin ) 
          qc0 (:,:,k) = ql_upd(:,:,k)
          qceq(:,:,k) = C_dt(:,:,k)/D_dt(:,:,k)
          qc1(:,:,k) = qceq(:,:,k) - (qceq(:,:,k) - qc0(:,:,k))*  &
                                                   exp( -1.* D_dt(:,:,k) )
          qcbar(:,:,k) = qceq(:,:,k) - ((qc1(:,:,k) -   &
                                                  qc0(:,:,k))/D_dt(:,:,k))
        elsewhere
          qc0(:,:,k)   = ql_upd(:,:,k)
          qceq(:,:,k)  = qc0(:,:,k) + C_dt(:,:,k)   
          qc1(:,:,k)   = qc0(:,:,k) + C_dt(:,:,k)
          qcbar(:,:,k) = qc0(:,:,k) + 0.5*C_dt(:,:,k)

        end where

!------------------------------------------------------------------------
!    set total tendency term and update cloud
!    Note that the amount of SL calculated here is stored in tmp1.
!------------------------------------------------------------------------
        SL(:,:,k) = SL(:,:,k) + qc1(:,:,k) - qc0(:,:,k)
        tmp1(:,:,k) = qc1(:,:,k) - qc0(:,:,k)        
        ql_upd(:,:,k) = qc1(:,:,k)

!------------------------------------------------------------------------
!    compute the amount each term contributes to the change

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

!       initialize tmp2 to hold (-Dterm)/D
!------------------------------------------------------------------------
     if (Nml%retain_cm3_bug) then
        tmp2(:,:,k) = D_dt(:,:,k)*qcbar(:,:,k)/max(D_dt(:,:,k), Nml%Dmin)
 
     else
          where (D_dt(:,:,k) > Nml%Dmin) 
          tmp2(:,:,k) = D_dt(:,:,k)*qcbar(:,:,k)/max(D_dt(:,:,k), Nml%Dmin)
          elsewhere
          tmp2(:,:,k) = 0.0              
          endwhere
      endif

         if (diag_id%qldt_cond + diag_id%ql_cond_col > 0)  &
                 diag_4d(:,:,k,diag_pt%qldt_cond) =    &
                                       max(dcond_ls(:,:,k),0.) *inv_dtcloud
         if (diag_id%qldt_evap + diag_id%ql_evap_col > 0)   &
              diag_4d(:,:,k,diag_pt%qldt_evap) =    &
                  - (max(0., -1.*dcond_ls(:,:,k) )/  &
                     max(ql_mean(:,:,k), Nml%qmin))*tmp2(:,:,k)*inv_dtcloud
        if (diag_id%qldt_accr  + diag_id%ql_accr_col > 0)   &
            diag_4d(:,:,k,diag_pt%qldt_accr) =    &
                  diag_4d(:,:,k,diag_pt%qldt_accr)*tmp2(:,:,k)*inv_dtcloud
        if (diag_id%qldt_auto  + diag_id%ql_auto_col > 0)  &
            diag_4d(:,:,k,diag_pt%qldt_auto) =    &
                 diag_4d (:,:,k,diag_pt%qldt_auto)*tmp2(:,:,k)*inv_dtcloud
        if (diag_id%qldt_eros + diag_id%ql_eros_col > 0)  &
            diag_4d(:,:,k,diag_pt%qldt_eros) =    &
                                  -  D_eros(:,:,k)*tmp2(:,:,k)*inv_dtcloud 
        if (diag_id%qldt_berg + diag_id%ql_berg_col > 0)   &
           diag_4d(:,:,k,diag_pt%qldt_berg) =    &
                 diag_4d (:,:,k,diag_pt%qldt_berg)*tmp2(:,:,k)*inv_dtcloud
        sum_berg(:,:,k) = sum_berg(:,:,k)*tmp2(:,:,k)*inv_dtcloud
        if (diag_id%qldt_rime  + diag_id%ql_rime_col > 0)   &
           diag_4d(:,:,k,diag_pt%qldt_rime) =   &
                 diag_4d (:,:,k,diag_pt%qldt_rime)*tmp2(:,:,k)*inv_dtcloud
        sum_rime(:,:,k) = sum_rime(:,:,k)*tmp2(:,:,k)*inv_dtcloud
        if (diag_id%qldt_freez + diag_id%ql_freez_col > 0)  &
                 diag_4d(:,:,k,diag_pt%qldt_freez) =   &
                  diag_4d(:,:,k,diag_pt%qldt_freez)*tmp2(:,:,k)*inv_dtcloud
        sum_freeze(:,:,k) = sum_freeze(:,:,k)*tmp2(:,:,k)*inv_dtcloud
!-------------------------------------------------------------------------
!    do phase changes from large-scale processes and boundary
!    layer condensation/evaporation
!-------------------------------------------------------------------------
        ST(:,:,k) = ST(:,:,k) + (hlv*max(dcond_ls(:,:,k), 0.)/cp_air) -   &
                    (hlv*(max(-1.*dcond_ls(:,:,k), 0.)/     &
                       max(ql_mean(:,:,k), Nml%qmin))*tmp2(:,:,k)/cp_air)
        SQ(:,:,k) = SQ (:,:,k) - max(dcond_ls(:,:,k),0.) +            &
                    (max(-1.*dcond_ls(:,:,k), 0.)/    &
                               max(ql_mean(:,:,k), Nml%qmin))*tmp2(:,:,k)

!------------------------------------------------------------------------  
!    if debug is requested, output some fields.
!------------------------------------------------------------------------  
        if (Nml%debugo .and.                        k .eq. Nml%ksamp) then
          write(otun,*) "ccc ST 5 ",   ST(Nml%isamp,Nml%jsamp,Nml%ksamp)
        end if

!------------------------------------------------------------------------
!    add in  temperature tendencies resuklting from liquid to ice  trans-
!    formations and cloud erosion.
!------------------------------------------------------------------------
        ST(:,:,k) = ST(:,:,k) + (hlf*D2_dt(:,:,k) - hlv*D_eros(:,:,k))* &
                                                        tmp2(:,:,k)/cp_air

!------------------------------------------------------------------------  
!    if debug is requested, output some fields.
!------------------------------------------------------------------------  
        if (Nml%debugo .and.                        k .eq. Nml%ksamp) then
          write(otun,*) "ccc D_eros, D2_dt, tmp2 ",   &
                         D_eros(Nml%isamp, Nml%jsamp,Nml%ksamp),   &
                         D2_dt(Nml%isamp,Nml%jsamp,Nml%ksamp),&
                         tmp2(Nml%isamp,Nml%jsamp,Nml%ksamp)
          write(otun,*) "ccc ST 6 ",   ST(Nml%isamp,Nml%jsamp,Nml%ksamp)
        end if

!------------------------------------------------------------------------  
!    cloud evaporation adds to water vapor
!------------------------------------------------------------------------  
        SQ(:,:,k) = SQ(:,:,k) + D_eros(:,:,k)*tmp2(:,:,k)  
             
!------------------------------------------------------------------------  
!    add conversion of liquid to rain to the rainflux
!------------------------------------------------------------------------  
        rain_cld(:,:,k) = rain_cld(:,:,k) + D1_dt(:,:,k)*tmp2(:,:,k)*  &
                                                 deltpg(:,:,k)*inv_dtcloud
     
!------------------------------------------------------------------------  
!    save liquid converted to ice into tmp3 and increment qi_mean
!------------------------------------------------------------------------  
        tmp3(:,:,k)    = tmp2(:,:,k)*D2_dt(:,:,k)
        qi_mean(:,:,k) = qi_mean(:,:,k) + tmp3(:,:,k)
 
        
!******************************************************************
!
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
!------------------------------------------------------------------------
        if (Nml%do_liq_num) then
            
!------------------------------------------------------------------------
!    calculate C_dt
!------------------------------------------------------------------------
          IF (    rk_act_only_if_ql_gt_qmin) THEN
            where (ql_upd(:,:,k) .GT. Nml%qmin ) 
              C_dt(:,:,k) = max (delta_cf(:,:,k), 0.)*drop1(:,:,k)*  &
                                                      1.e6/airdens(:,:,k)
            elsewhere
              C_dt(:,:,k)=0.
            end where
          ELSE
            C_dt(:,:,k)=max(delta_cf(:,:,k), 0.)*drop1(:,:,k)*  &
                                                       1.e6/airdens(:,:,k)
          END IF

!------------------------------------------------------------------------
!    set total tendency term and update cloud
!    Note that the amount of SN calculated here is stored in tmp1.
!------------------------------------------------------------------------
          D_dt(:,:,k) =  Nml%num_mass_ratio1*D1_dt(:,:,k) +   &
                        (Nml%num_mass_ratio2*D2_dt(:,:,k) + D_eros(:,:,k))
          qc0(:,:,k) = qn_upd(:,:,k)
          where (D_dt(:,:,k) > Nml%Dmin) 
            qceq(:,:,k) = C_dt(:,:,k) / D_dt(:,:,k)
            qc1(:,:,k) = qceq(:,:,k) - (qceq(:,:,k) - qc0(:,:,k))*  &
                                                    exp(-1.*D_dt(:,:,k))
            qcbar(:,:,k) = qceq(:,:,k) - ((qc1(:,:,k) -qc0(:,:,k))/ &
                                                            D_dt(:,:,k))
          elsewhere
            qceq(:,:,k) = qc0(:,:,k) + C_dt(:,:,k)

            qc1 (:,:,k) = qc0(:,:,k) + C_dt(:,:,k)
            qcbar(:,:,k) = qc0(:,:,k) + 0.5*C_dt(:,:,k)
          end where

!------------------------------------------------------------------------
!    set total tendency term and update cloud
!    Note that the amount of SN calculated here is stored in tmp1.
!------------------------------------------------------------------------
          SN(:,:,k) = SN(:,:,k) + qc1(:,:,k) - qc0(:,:,k)
          qn_upd(:,:,k)     = qc1(:,:,k)

!------------------------------------------------------------------------
!    compute the amount each term contributes to the change for diagnostics
!    d2_dt stuff not fully implemented ...
!        where ( C_dt .gt. 0 )
!                Cterm  =  C_dt             
!                Dterm  =  D_dt *      qcbar 
!        elsewhere
!                Cterm  =  0.             
!                Dterm  =  D_dt *      qcbar 
!        end where
!------------------------------------------------------------------------
          if (diag_id%qndt_cond + diag_id%qn_cond_col > 0) then
            where ( C_dt(:,:,k) .gt. 0 )
              diag_4d(:,:,k,diag_pt%qndt_cond) = C_dt(:,:,k)
            elsewhere
              diag_4d(:,:,k,diag_pt%qndt_cond) = 0.
            endwhere
          end if  
       if (Nml%retain_cm3_bug) then
          tmp8(:,:,k) = D_dt(:,:,k)*qcbar(:,:,k)/max(D_dt(:,:,k), Nml%Dmin)
       else
          where (D_dt(:,:,k) > Nml%Dmin) 
          tmp8(:,:,k) = D_dt(:,:,k)*qcbar(:,:,k)/max(D_dt(:,:,k), Nml%Dmin)
          elsewhere
          tmp8(:,:,k) = 0.0              
          endwhere
       endif

          if (diag_id%qndt_pra  + diag_id%qn_pra_col > 0) then
            diag_4d(:,:,k,diag_pt%qndt_pra) =    &
                   Nml%num_mass_ratio1*diag_4d(:,:,k,diag_pt%qndt_pra)* &
                                                    tmp8(:,:,k)*inv_dtcloud
          end if
          if (diag_id%qndt_auto  + diag_id%qn_auto_col > 0) then
            diag_4d(:,:,k,diag_pt%qndt_auto) =   &
                   Nml%num_mass_ratio1*diag_4d(:,:,k,diag_pt%qndt_auto)* &
                                                   tmp8(:,:,k) *inv_dtcloud
          end if
          if (diag_id%qndt_berg  + diag_id%qn_berg_col > 0) then
            diag_4d(:,:,k,diag_pt%qndt_berg) =  &
                 Nml%num_mass_ratio2*diag_4d(:,:,k,diag_pt%qndt_berg)*   &
                                                   tmp8(:,:,k) *inv_dtcloud
          end if
          if (diag_id%qndt_eros  + diag_id%qn_eros_col > 0) then
            diag_4d(:,:,k,diag_pt%qndt_eros) = -D_eros(:,:,k)*  &
                                                   tmp8(:,:,k)*inv_dtcloud
          end if
        end if  ! (Nml%do_liq_num)

        if (Nml%do_liq_num) then
          if (diag_id%qndt_cond + diag_id%qn_cond_col > 0)   &
                   diag_4d(:,:,k,diag_pt%qndt_cond) =   &
                               diag_4d(:,:,k,diag_pt%qndt_cond)*inv_dtcloud
          if (diag_id%qndt_evap + diag_id%qn_evap_col > 0)   &
                   diag_4d(:,:,k,diag_pt%qndt_evap) =   &
                              diag_4d(:,:,k,diag_pt%qndt_evap)*inv_dtcloud 
        endif  

!----------------------------------------------------------------------! 
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
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!    compute Vfall
!------------------------------------------------------------------------
        where (qi_mean(:,:,k) > 0.0) 
          Vfall(:,:,k) = Nml%vfact*3.29*((airdens(:,:,k) *qi_mean(:,:,k)/ &
                                     max(qa_mean(:,:,k), Nml%qmin))**0.16)
        elsewhere
          Vfall(:,:,k) = 0.
        end where
        if (diag_id%vfall > 0) diag_4d(:,:,k,diag_pt%vfall) = Vfall(:,:,k)* qa_mean(:,:,k)

!------------------------------------------------------------------------
!    add to ice source the settling ice flux from above
!    also note that tmp3 contains the source of liquid converted to ice 
!    from above
!------------------------------------------------------------------------
        C_dt(:,:,k) = tmp3(:,:,k) + snow_cld(:,:,k)*dtcloud/deltpg(:,:,k)
        
!------------------------------------------------------------------------
!    Compute settling of ice. The result is multiplied by dtcloud/deltp to
!    convert to units of D_dt. Note that if tracers are not advected then 
!    this is done relative to the local vertical motion.
!------------------------------------------------------------------------
        if (Nml%tracer_advec) then
          tmp1(:,:,k) = 0.
        else
          tmp1(:,:,k) = omega(:,:,k)
        end if 
        
        where (qi_mean(:,:,k) .gt. Nml%qmin .and.   &
               qa_mean(:,:,k) .gt. Nml%qmin)
          D1_dt(:,:,k) = max(0., ((airdens(:,:,k)*Vfall(:,:,k)) +   &
                                (tmp1(:,:,k)/grav))*dtcloud/deltpg(:,:,k) )
          a_snow_cld(:,:,k) = qa_mean(:,:,k)     
        elsewhere
          D1_dt(:,:,k)      = 0.
          snow_cld(:,:,k)   = 0.
          a_snow_cld(:,:,k) = 0.    
        end where 

!-------------------------------------------------------------------------
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
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!    compute grid mean change in cloud ice to cool the grid box to 0C 
!    and store in temporary variable tmp1
!------------------------------------------------------------------------
        tmp1(:,:,k) = cp_air*(T(:,:,k)-tfreeze)/hlf

!-------------------------------------------------------------------------
!    If qi_mean > tmp1, then the amount of ice melted is limited to tmp1, 
!    otherwise melt all qi_mean.  The amount melted is stored in tmp2
!------------------------------------------------------------------------
        tmp2(:,:,k)  = max(min(qi_mean(:,:,k), tmp1(:,:,k)),0.)     
        D2_dt(:,:,k) = max(0., log(max(qi_mean(:,:,k),Nml%qmin)/   &
                              max(qi_mean(:,:,k) - tmp2(:,:,k), Nml%qmin)))
            
!-----------------------------------------------------------------------
!    melting of ice adds area to a_rain_cld
!-----------------------------------------------------------------------
        where (D2_dt(:,:,k) .gt. Nml%Dmin) 
          a_rain_cld(:,:,k) = qa_mean(:,:,k)
        endwhere
                  
!------------------------------------------------------------------------
!    If all of the ice can melt, then don't permit any ice to fall out of
!    the grid box and set a_snow_cld to zero.
!-------------------------------------------------------------------------
        where (qi_mean(:,:,k) .lt. tmp1(:,:,k) .and.   &
               qi_mean(:,:,k) .gt. Nml%qmin) 
          D1_dt(:,:,k) = 0.
          snow_cld(:,:,k) = 0.
          a_snow_cld(:,:,k) = 0.
        end where
         
!------------------------------------------------------------------------
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
!-----------------------------------------------------------------------
        C_dt(:,:,k) = C_dt(:,:,k) + max(dcond_ls_ice(:,:,k), 0.)
        sum_cond(:,:,k) =  max(dcond_ls_ice(:,:,k), 0.)*inv_dtcloud
        D_dt(:,:,k) = D1_dt(:,:,k) + D2_dt(:,:,k) + D_eros(:,:,k) +   &
                       (max(-1.*dcond_ls_ice(:,:,k), 0.)/   &
                                max(qi_mean(:,:,k), Nml%qmin))
        
!-----------------------------------------------------------------------
!    do analytic integration      
!-----------------------------------------------------------------------
        where ( D_dt(:,:,k).gt.Nml%Dmin ) 
          qc0(:,:,k) = qi_upd(:,:,k)
          qceq(:,:,k) = C_dt(:,:,k)/D_dt(:,:,k)
          qc1(:,:,k) = qceq(:,:,k) - (qceq(:,:,k) - qc0(:,:,k))*  &
                                                  exp( -1.* D_dt(:,:,k) )
          qcbar(:,:,k) = qceq(:,:,k) - ((qc1(:,:,k) - qc0(:,:,k))/  &
                                                             D_dt(:,:,k))
        elsewhere
          qc0(:,:,k) = qi_upd(:,:,k)
          qceq(:,:,k) = qc0(:,:,k)+ C_dt(:,:,k)   
          qc1(:,:,k) = qc0(:,:,k) + C_dt(:,:,k)
          qcbar(:,:,k) = qc0(:,:,k) + 0.5*C_dt(:,:,k)
        end where

!-----------------------------------------------------------------------
!    set total tendency term and update cloud
!    Note that the amount of SL calculated here is stored in tmp1.
!-----------------------------------------------------------------------
        SI (:,:,k) = SI(:,:,k) + qc1(:,:,k) - qc0(:,:,k)
        tmp1(:,:,k) = qc1(:,:,k) - qc0(:,:,k)        
        qi_upd(:,:,k) = qc1(:,:,k)

!-----------------------------------------------------------------------
!    compute the amount each term contributes to the change     
!    Apportion SI between various processes.  This is necessary to
!    account for how much the temperature and water vapor changes 
!    due to various phase changes.   For example:
!
!    ice settling = (D1/D)*(-Dterm)*deltp/grav/dtcloud
!                     
!    vapor to ice =
!           -{ ((-dcond_ls_ice/qi_mean)+D_eros)/ D }*(-Dterm) 
!           where dcond_ls_ice  < 0. 
!
!    but
!    vapor to ice =
!           dcond_ls_ice -(D_eros/D)* (-Dterm)
!           where dcond_ls_ice > 0.
!       
!    melting of ice = (D2/D)*(-Dterm)*deltp/grav/dtcloud
!-----------------------------------------------------------------------

!------------------------------------------------------------------------
!    initialize tmp2 to hold (-Dterm)/D
!------------------------------------------------------------------------
   if (Nml%retain_cm3_bug) then
        tmp2 (:,:,k) = D_dt(:,:,k)*qcbar(:,:,k)/max(D_dt(:,:,k), Nml%Dmin)
   else
       where (D_dt(:,:,k) > Nml%Dmin)
          tmp2(:,:,k) = D_dt(:,:,k)*qcbar(:,:,k)/max(D_dt(:,:,k), Nml%Dmin)
      elsewhere
          tmp2(:,:,k) = 0.0
      endwhere
   endif

!------------------------------------------------------------------------
!    do phase changes from large-scale processes 
!------------------------------------------------------------------------
        ST(:,:,k) = ST(:,:,k) + hls*max(dcond_ls_ice(:,:,k) ,0.)/cp_air - &
                         hls*(max(-1.*dcond_ls_ice(:,:,k), 0.)/   &
                         max(qi_mean(:,:,k), Nml%qmin))*tmp2(:,:,k)/cp_air
        SQ(:,:,k)  = SQ(:,:,k) - max(dcond_ls_ice(:,:,k) ,0.) +       &
                         (max(-1.*dcond_ls_ice(:,:,k), 0.)/  &
                                max(qi_mean(:,:,k), Nml%qmin))*tmp2(:,:,k) 

!------------------------------------------------------------------------
!    cloud erosion changes temperature and vapor
!------------------------------------------------------------------------
        ST(:,:,k) = ST(:,:,k) - hls*D_eros(:,:,k)*tmp2(:,:,k)/cp_air
        SQ(:,:,k) = SQ (:,:,k) + D_eros(:,:,k)*tmp2(:,:,k) 
 
!-----------------------------------------------------------------------
!    add settling ice flux to snow_cld 
!-----------------------------------------------------------------------
        snow_cld(:,:,k) = D1_dt(:,:,k)*tmp2(:,:,k)*deltpg(:,:,k)*  &
                                                               inv_dtcloud
       
!-----------------------------------------------------------------------
!    add melting of ice to temperature tendency
!-----------------------------------------------------------------------
        ST(:,:,k) = ST(:,:,k) - hlf*D2_dt(:,:,k)*tmp2(:,:,k)/cp_air

!-----------------------------------------------------------------------
!    add melting of ice to the rainflux
!-----------------------------------------------------------------------
        rain_cld(:,:,k) = rain_cld(:,:,k) + D2_dt(:,:,k)*tmp2(:,:,k)*  &
                                                deltpg(:,:,k)*inv_dtcloud

!------------------------------------------------------------------------
!    diagnostics for cloud ice tendencies
!------------------------------------------------------------------------
        if (diag_id%qidt_dep + diag_id%qi_dep_col > 0)      &
                  diag_4d (:,:,k,diag_pt%qidt_dep) =   &
                                  max(dcond_ls_ice(:,:,k), 0.)*inv_dtcloud
        if (diag_id%qidt_subl + diag_id%qi_subl_col > 0)   &
                  diag_4d(:,:,k,diag_pt%qidt_subl) =   &
                        -(max(0., -1.*dcond_ls_ice(:,:,k) )/   &
                           max(qi_mean(:,:,k), Nml%qmin))*tmp2(:,:,k)  &
                                                              *inv_dtcloud
        if (diag_id%qidt_melt + diag_id%qi_melt_col > 0)    &
                  diag_4d(:,:,k,diag_pt%qidt_melt) = -D2_dt(:,:,k)*  &
                                                   tmp2(:,:,k) *inv_dtcloud
        if (diag_id%qidt_eros + diag_id%qi_eros_col > 0)    &
                  diag_4d(:,:,k,diag_pt%qidt_eros) =  - D_eros(:,:,k)*  &
                                                   tmp2(:,:,k) *inv_dtcloud
        
!----------------------------------------------------------------------! 
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
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!    compute U_clr. keep U_clr > 0. and U_clr < 1.
!------------------------------------------------------------------------
        if (.not. Nml%do_pdf_clouds) then 
          where (snow_clr(:,:,k) > 0. .or.   &
                 rain_clr(:,:,k) > 0. )           
            U_clr(:,:,k) = (U(:,:,k) - qa_mean(:,:,k))/  &
                                      max((1. - qa_mean(:,:,k)), Nml%qmin)
          elsewhere
            U_clr(:,:,k) = 0.
          end where
        else
          U_clr(:,:,k) = qvg(:,:,k)/qs(:,:,k)
        end if
        U_clr(:,:,k) = min(max(U_clr(:,:,k), 0.), 1.)
        
!------------------------------------------------------------------------
!    compute K3
!------------------------------------------------------------------------
        where (rain_clr(:,:,k) > 0)
          tmp1(:,:,k) = 56788.636*dtcloud*((rain_clr(:,:,k)/   &
                            max(a_rain_clr(:,:,k), Nml%qmin)/  &
                                           dens_h2o)**(11./18.))/  &
                            SQRT(airdens(:,:,k))/A_plus_B(:,:,k)/qs(:,:,k)

!------------------------------------------------------------------------
!    compute local change in vapor mixing ratio due to rain evaporation
!------------------------------------------------------------------------
          tmp1(:,:,k) = tmp1(:,:,k)*qs(:,:,k)*(1. - U_clr(:,:,k))/  &
                               (1. + 0.5*tmp1(:,:,k)*(1. + gamma(:,:,k)))

!------------------------------------------------------------------------
!    limit change in qv to the amount that would raise the relative
!    humidity to U_evap in the clear portion of the grid box
!------------------------------------------------------------------------
          tmp1(:,:,k) = min(tmp1(:,:,k), ((1.-qa_mean(:,:,k))/  &
                            max(a_rain_clr(:,:,k), Nml%qmin))*qs(:,:,k)* &
                          max(0. ,Nml%U_evap-U_clr(:,:,k))/  &
                            (1. + (Nml%U_evap*(1. - qa_mean(:,:,k)) +  &
                                           qa_mean(:,:,k))*gamma(:,:,k)) )
        
!------------------------------------------------------------------------
!    do limiter by amount available
!------------------------------------------------------------------------
          tmp1(:,:,k) = tmp1(:,:,k)*a_rain_clr(:,:,k)*  &
                                                 deltpg(:,:,k)*inv_dtcloud
          tmp2(:,:,k) = max(min(rain_clr(:,:,k), tmp1(:,:,k)), 0.)
    
          SQ(:,:,k) = SQ(:,:,k) + tmp2(:,:,k)*dtcloud/deltpg(:,:,k)
          ST(:,:,k) = ST(:,:,k) - hlv*tmp2(:,:,k)*dtcloud/   &
                                                      deltpg(:,:,k)/cp_air
        elsewhere
          tmp1(:,:,k) = 0.
          tmp2(:,:,k) = 0.
        end where
                                                        
!------------------------------------------------------------------------
!    if all of the rain evaporates set things to zero.    
!------------------------------------------------------------------------
        where (tmp1(:,:,k) .gt. rain_clr(:,:,k) .and.    &
               a_rain_clr(:,:,k) .gt. Nml%qmin)         
          rain_clr(:,:,k) = 0.
          a_rain_clr(:,:,k) = 0.
        elsewhere
          rain_clr(:,:,k) = rain_clr(:,:,k) - tmp2 (:,:,k)    
        endwhere
        
!------------------------------------------------------------------------
!    save evaporation diagnostics
!------------------------------------------------------------------------
        if (diag_id%rain_evap + diag_id%rain_evap_col > 0)     &
              diag_4d(:,:,k,diag_pt%rain_evap) = tmp2(:,:,k)/deltpg(:,:,k)

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
!-------------------------------------------------------------------------
        where (snow_clr(:,:,k) > 0.0) 

!-------------------------------------------------------------------------
!    compute K3.
!-------------------------------------------------------------------------
          tmp1(:,:,k) = dtcloud*(4./3.14159/rho_ice/airdens(:,:,k)/       &
                        A_plus_B(:,:,k)/qs(:,:,k))*((snow_clr(:,:,k)/  &
                       max(a_snow_clr(:,:,k), Nml%qmin)/3.29)**(1./1.16))*&
                        (0.65*lamda_f(:,:,k)*lamda_f(:,:,k) +         &
                            198.92227*lamda_f(:,:,k)*  &
                                  SQRT(airdens(:,:,k)*lamda_f(:,:,k))*   &
                                          ( (snow_clr(:,:,k)/    &
                           max(a_snow_clr(:,:,k),Nml%qmin))**(1./14.5) )  )

!-------------------------------------------------------------------------
!    compute local change in vapor mixing ratio due to snow sublimation.
!-------------------------------------------------------------------------
          tmp1(:,:,k) = tmp1(:,:,k)*qs(:,:,k)*(1. - U_clr(:,:,k))/  &
                                (1. + 0.5*tmp1(:,:,k)*(1. + gamma(:,:,k)))

!-------------------------------------------------------------------------
!    limit change in qv to the amount that would raise the relative
!    humidity to U_evap in the clear portion of the grid box
!-------------------------------------------------------------------------
          tmp1 (:,:,k) = min(tmp1(:,:,k), ((1. - qa_mean(:,:,k))/  &
                           max(a_snow_clr(:,:,k), Nml%qmin))*qs(:,:,k)* &
                           max(0., Nml%U_evap  -U_clr(:,:,k))/  &
                               (1. + (Nml%U_evap*(1. - qa_mean(:,:,k)) + &
                                           qa_mean(:,:,k))*gamma(:,:,k)) )
        
!-------------------------------------------------------------------------
!    do limiter by amount available
!-------------------------------------------------------------------------
          tmp1(:,:,k) = tmp1(:,:,k)*a_snow_clr(:,:,k)*    &
                                                deltpg(:,:,k)*inv_dtcloud
        elsewhere
          tmp1(:,:,k)= 0.0
        endwhere
        tmp2(:,:,k) = max(min(snow_clr(:,:,k), tmp1(:,:,k)), 0.)
        SQ(:,:,k) = SQ(:,:,k) + tmp2(:,:,k)*dtcloud/deltpg(:,:,k)
        ST(:,:,k) = ST(:,:,k) - hls*tmp2(:,:,k)*dtcloud/deltpg(:,:,k)/  &
                                                                    cp_air

!-------------------------------------------------------------------------
!    if all of the snow sublimates set things to zero.    
!-------------------------------------------------------------------------
        where (tmp1(:,:,k) .gt. snow_clr(:,:,k) .and.   &
               a_snow_clr(:,:,k) .gt. Nml%qmin)         
          snow_clr(:,:,k) = 0.
          a_snow_clr(:,:,k) = 0.
        elsewhere
          snow_clr(:,:,k) = snow_clr(:,:,k) - tmp2(:,:,k)     
        endwhere
         
!-------------------------------------------------------------------------
!    save snow sublimation diagnostics.
!-------------------------------------------------------------------------
        if  (diag_id%qdt_snow_sublim + diag_id%q_snow_sublim_col > 0)     &
              diag_4d(:,:,k,diag_pt%qdt_snow_sublim) = tmp2(:,:,k)/deltpg(:,:,k) 
       
!------------------------------------------------------------------------
!    save diagnostics for the predicted cloud tendencies so that 
!    destruction diagnostics may be computed after adjustment is done.
!------------------------------------------------------------------------
        if (diag_id%qadt_destr + diag_id%qa_destr_col > 0)    &
               diag_4d(:,:,k,diag_pt%qadt_destr) =    &
                 diag_4d(:,:,k,diag_pt%qadt_destr) + SA(:,:,k) *inv_dtcloud
        if   (diag_id%qldt_destr + diag_id%ql_destr_col > 0)    &
              diag_4d(:,:,k,diag_pt%qldt_destr) =   &
                 diag_4d(:,:,k,diag_pt%qldt_destr) + SL(:,:,k) *inv_dtcloud
        if   (diag_id%qidt_destr + diag_id%qi_destr_col > 0)    &
               diag_4d(:,:,k,diag_pt%qidt_destr) =   &
                 diag_4d(:,:,k,diag_pt%qidt_destr) + SI(:,:,k)*inv_dtcloud
        if   (diag_id%qndt_destr + diag_id%qn_destr_col > 0 .and.    &
                    Nml%do_liq_num )      &
                diag_4d(:,:,k,diag_pt%qndt_destr) = &
                 diag_4d(:,:,k,diag_pt%qndt_destr) + SN(:,:,k)*inv_dtcloud
      
!-----------------------------------------------------------------------
!    adjustment to try to prevent negative water vapor. adjustment will not
!    remove more condensate than available. cloud amount is not adjusted. 
!    this is left for the remainder of the destruction code. (cjg)
!    tmp1 = updated vapor
!    tmp2 = liquid which can be evaporated
!    tmp3 = ice which can be sublimated
!-----------------------------------------------------------------------
        tmp1(:,:,k) = qv(:,:,k) + SQ(:,:,k)
        tmp2(:,:,k) = 0.0
        tmp3(:,:,k) = 0.0
        where ( tmp1(:,:,k).lt.0.0 .and. T(:,:,K).le.tfreeze-40.)
          tmp2(:,:,k) = min( -tmp1(:,:,k), ql_upd(:,:,k) ) 
          tmp3 (:,:,k) = min( -tmp1(:,:,k) - tmp2(:,:,k), qi_upd(:,:,k) ) 
          ql_upd(:,:,k) = ql_upd (:,:,k) - tmp2(:,:,k)
          qi_upd(:,:,k) = qi_upd(:,:,k) - tmp3(:,:,k)
          SL(:,:,k) = SL(:,:,k) - tmp2(:,:,k)
          SI(:,:,k) = SI(:,:,k) - tmp3(:,:,k)
          SQ(:,:,k) = SQ(:,:,k) + tmp2(:,:,k) + tmp3(:,:,k)
          ST(:,:,k) = ST(:,:,k) - hlv*tmp2(:,:,k)/cp_air -   &
                                               hls*tmp3(:,:,k)/cp_air
        end where
        where ( tmp1(:,:,k).lt.0.0 .and. T(:,:,k).gt.tfreeze-40.)
          tmp3(:,:,k) = min( -tmp1(:,:,k), qi_upd(:,:,k) ) 
          tmp2(:,:,k) = min( -tmp1(:,:,k) - tmp3(:,:,k), ql_upd(:,:,k) )
          ql_upd(:,:,k) = ql_upd(:,:,k) - tmp2(:,:,k)
          qi_upd(:,:,k) = qi_upd(:,:,k) - tmp3(:,:,k)
          SL(:,:,k) = SL(:,:,k) - tmp2(:,:,k)
          SI(:,:,k) = SI(:,:,k) - tmp3(:,:,k)
          SQ(:,:,k) = SQ(:,:,k) + tmp2(:,:,k) + tmp3(:,:,k)
          ST(:,:,k) = ST(:,:,k) - hlv*tmp2(:,:,k)/cp_air -    &
                                               hls*tmp3(:,:,k)/cp_air
        end where
      end do !  (big k loop) 

!-----------------------------------------------------------------------
!    cloud destruction occurs where both ql and qi are .le. qmin, or if 
!    qa is .le. qmin. In this case, ql, qi, and qa are set to zero 
!    conserving moisture and energy.
!-----------------------------------------------------------------------
      if (.not.Nml%do_liq_num) then
        where ((ql_upd .le. Nml%qmin .and. qi_upd .le. Nml%qmin) .or.   &
               (qa_upd .le. Nml%qmin))
          SL = SL - ql_upd 
          SI = SI - qi_upd 
          SQ = SQ + ql_upd + qi_upd 
          ST = ST - (hlv*ql_upd + hls*qi_upd )/cp_air
          SA = SA  - qa_upd 
          ql_upd = 0.0
          qi_upd = 0.0
          qa_upd = 0.0
        end where
      else
        where ((ql_upd .le. Nml%qmin .and. qi_upd .le. Nml%qmin) .or.  &
               (qa_upd .le. Nml%qmin))
          SL = SL - ql_upd 
          SI = SI - qi_upd 
          SQ = SQ + ql_upd + qi_upd 
          ST = ST - (hlv*ql_upd + hls*qi_upd)/cp_air
          SA = SA - qa_upd 
          SN = SN - qn_upd 
          ql_upd = 0.0
          qi_upd = 0.0
          qa_upd = 0.0
          qn_upd = 0.0
        end where
      endif  

!------------------------------------------------------------------------
!    save destruction diagnostics.
!------------------------------------------------------------------------
      if   (diag_id%qadt_destr + diag_id%qa_destr_col > 0)   &
               diag_4d(:,:,:,diag_pt%qadt_destr) =    &
                     diag_4d(:,:,:,diag_pt%qadt_destr) - SA*inv_dtcloud
      if   (diag_id%qldt_destr + diag_id%ql_destr_col > 0)    &
               diag_4d(:,:,:,diag_pt%qldt_destr) =   &
                     diag_4d(:,:,:,diag_pt%qldt_destr) - SL*inv_dtcloud
      if   (diag_id%qidt_destr + diag_id%qi_destr_col > 0) &
               diag_4d(:,:,:,diag_pt%qidt_destr) =    &
                    diag_4d(:,:,:,diag_pt%qidt_destr) - SI*inv_dtcloud
      if   (diag_id%qndt_destr + diag_id%qn_destr_col > 0)    &
               diag_4d(:,:,:,diag_pt%qndt_destr) =   &
                    diag_4d(:,:,:,diag_pt%qndt_destr) - SN*inv_dtcloud
                                                           
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
!------------------------------------------------------------------------
      if (.not.Nml%do_pdf_clouds) then

!-------------------------------------------------------------------------
!    updated temperature in tmp8, updated qs in tmp2,  
!    updated gamma in tmp3, updated dqsdT in tmp5
!------------------------------------------------------------------------
        tmp8 = T + ST
        IF (rk_alt_adj_opt .and. Nml%super_choice ) THEN
          where ( tmp8 > tfreeze - 40. )      
            est_a = polysvp_l(tmp8, idim, jdim, kdim)   

!------------------------------------------------------------------------
!    calculate denominator in qsat formula
!    limit denominator to esat, and thus qs to d622
!    this is done to avoid blow up in the upper stratosphere
!    where pfull ~ esat  
!------------------------------------------------------------------------
            qs_d_a = pfull - d378*est_a
            qs_d_a = max(qs_d_a, est_a) 
            tmp2s_a = d622*est_a/qs_d_a 
            tmp5s_a = hlv*tmp2s_a/(rvgas*tmp8**2)
            tmp3s_a = tmp5s_a*hlv/cp_air

!-------------------------------------------------------------------------
!     calculate excess water vapor
!-------------------------------------------------------------------------
            tmp1 = max(0., qv + SQ - tmp2s_a)/(1. + tmp3s_a)
          elsewhere
            est_a = polysvp_i  (tmp8, idim, jdim,kdim)   

!------------------------------------------------------------------------
!    calculate denominator in qsat formula
!    limit denominator to esat, and thus qs to d622
!    this is done to avoid blow up in the upper stratosphere
!    where pfull ~ esat  
!------------------------------------------------------------------------
            qs_d_a = pfull - d378*est_a
            qs_d_a = max(qs_d_a, est_a) 
            tmp2s_a = d622*est_a/qs_d_a 
            tmp5s_a = hls*tmp2s_a/(rvgas*tmp8**2)
            tmp3s_a = tmp5s_a*hls/cp_air

!-------------------------------------------------------------------------
!     calculate excess water vapor
!-------------------------------------------------------------------------
            tmp1 = max(0., qv + SQ - tmp2s_a)/(1. + tmp3s_a)
          end where
        ELSE ! rk_alt_adj_if

!------------------------------------------------------------------------
!RSH 10/28/10 -- proposal -- needs additional investigation
!RSH to avoid having non-saturation after adjustment, the value of L used 
!    here must agree with that used below when adjustments made to ST due 
!    to condensation of the supersaturated vapor.
!    Therefore compute_qs should be called with optional 
!    argument esat_over_liq = .true. when t > tcrit  (-40 C when 
!    super_choice is .true.,  -20 C when .false.), and without the optional
!    argument when t is cooler than tcrit, so that in that case es is 
!    calculated wrt to ice. then hlv would be used for temps above -40C 
!    and hls used at t < -40C, consistent with the ST terms below, and the 
!    condition for an increase in droplet number below.
!
!   NOTE: this requires setting sat_vapor_pres_nml variable
!            construct_table_wrt_liq = .true.
!    so that es table wrt liquid at all temps is computed.
!
!   NOTE: The above is essentially what rk_alt_adj_opt is doing, but
!   with different formulae for es.
!------------------------------------------------------------------------
!         if (use_inconsistent_lh) then
            call compute_qs (tmp8, pfull, tmp2, dqsdT=tmp5)
            tmp3 = tmp5*(min(1., max(0.,     &
                                   0.05*(tmp8 - tfreeze + 20.)))*hlv +  &
                              min(1., max(0.,    &
                                      0.05*(tfreeze - tmp8)))*hls)/cp_air

!-------------------------------------------------------------------------
!    compute excess over saturation
!-------------------------------------------------------------------------
            tmp1 = max(0., qv + SQ - tmp2)/(1. + tmp3)
!RSH
!         else
!           if (Nml%super_choice) then
!             tcrit = tfreeze - 40.
!           else
!             tcrit = tfreeze - 20.
!           endif
!           call compute_qs (tmp8, pfull, tmp2, dqsdT=tmp5,  &
!                                                   es_over_liq = .true.)
!           where (tmp8 > tcrit) 
!             tmp3 = tmp5*hlv/cp_air
!-------------------------------------------------------------------------
!    compute excess over saturation
!-------------------------------------------------------------------------
!             tmp1 = max(0., qv + SQ - tmp2)/(1. + tmp3)
!           endwhere
!           call compute_qs (tmp8, pfull, tmp2, dqsdT=tmp5)
!           where (tmp8 <= tcrit) 
!             tmp3 = tmp5*hls/cp_air
!-------------------------------------------------------------------------
!    compute excess over saturation
!-------------------------------------------------------------------------
!             tmp1 = max(0., qv + SQ - tmp2)/(1. + tmp3)
!           end where
!         endif

        END IF !  rk_alt_adj_if

!------------------------------------------------------------------------
!    change vapor content
!------------------------------------------------------------------------
        SQ = SQ - tmp1

!------------------------------------------------------------------------
!    if super_choice is .true., then put supersaturation into cloud,
!------------------------------------------------------------------------
        if (Nml%super_choice) then


!------------------------------------------------------------------------
!    assign additional droplets where supersaturation is predicted. 
!------------------------------------------------------------------------
          if (Nml%do_liq_num) then
            do k=1,kdim
              do j=1,jdim
                do i=1,idim
                  if (T(i,j,k) > tfreeze - 40. .and. &
                      tmp1(i,j,k) > 0.0) then
                    qn_upd(i,j,k) = qn_upd(i,j,k) + drop1(i,j,k)*1.0e6/  &
                                    airdens(i,j,k)*(1. - qa_upd(i,j,k))
                    SN(i,j,k) = SN(i,j,k) + drop1(i,j,k)*1.e6/  &
                                   airdens(i,j,k)*(1.-qa_upd(i,j,k))
                  endif
                end do
              end do
            end do
          endif
          if (diag_id%qndt_super + diag_id%qn_super_col > 0) then
            if (Nml%do_liq_num) then
              where (T .gt. tfreeze - 40. .and. tmp1 .gt. 0.)
                diag_4d(:,:,:,diag_pt%qndt_super) =     &
                         diag_4d(:,:,:,diag_pt%qndt_super) + drop1*1.e6/ &
                                         airdens*(1. - qa_upd)*inv_dtcloud
              endwhere
            endif 
          end if

!------------------------------------------------------------------------
!  THE -40C THRESHOLD IN THE FOLLOWING IS NOT CONSISTENT WITH THE 
!  CALCULATION OF QS ABOVE. IT IS UNPHYSICAL TO ASSUME THE FORMATION OF 
!  SUPERCOOLED DROPLETS BETWEEN QS_ICE AND QS_LIQ! ( ... but similar 
!  assumptions are occasionally made in cloud cover schemes) 
!  RSH:  this may be addressed by proposed mod of 10/28/10.
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!    add excess condensate to appropriate cloud condensate species. 
!    change cloud area and increment temperature.
!------------------------------------------------------------------------
          do k=1,kdim
            do j=1,jdim
              do i=1,idim
                if (tmp1(i,j,k) .gt. 0.) then
                  if (T(i,j,k) .le. tfreeze - 40.)then
                    qi_upd(i,j,k)   = qi_upd(i,j,k) + tmp1(i,j,k)
                    SI(i,j,k) = SI(i,j,k) + tmp1(i,j,k)
                    ST(i,j,k) = ST(i,j,k) + hls*tmp1(i,j,k)/cp_air
                  else   
                    ql_upd(i,j,k)   = ql_upd(i,j,k) + tmp1(i,j,k)
                    SL(i,j,k) = SL(i,j,k) + tmp1(i,j,k)
                    ST(i,j,k) = ST(i,j,k) + hlv*tmp1(i,j,k)/cp_air      
                  endif
                  if (limit_conv_cloud_frac) then
                    tmp2s = ahuco(i,j,k)
                  else
                    tmp2s = 0.
                  endif
!moved from above:
!------------------------------------------------------------------------
!    cloud fraction source diagnostic
!------------------------------------------------------------------------
          if (diag_id%qadt_super + diag_id%qa_super_col > 0) then
              diag_4d(i,j,k,diag_pt%qadt_super) = (1.-qa_upd(i,j,k)-tmp2s)*inv_dtcloud
          end if
                  SA(i,j,k) = SA(i,j,k) + (1.-qa_upd(i,j,k) - tmp2s)  
                  qa_upd(i,j,k)   = 1. - tmp2s        
                endif
              enddo
            enddo
          enddo
                                                               
!------------------------------------------------------------------------
!    save adjustment diagnostics.
!------------------------------------------------------------------------
            where (T .le. tfreeze - 40.)
              sum_ice_adj(:,:,:) = tmp1*inv_dtcloud
            endwhere
          if (diag_id%liq_adj  + diag_id%liq_adj_col +   &
              diag_id%ice_adj + diag_id%ice_adj_col  > 0) then       
            where (T .le. tfreeze - 40.)
              diag_4d(:,:,:,diag_pt%ice_adj) = tmp1*inv_dtcloud
            elsewhere
              diag_4d(:,:,:,diag_pt%liq_adj) = tmp1*inv_dtcloud
            endwhere
          end if

!-----------------------------------------------------------------------
!    put supersaturation into precip. add in excess to precipitation 
!    fluxes, change their area and increment temperature.
!-----------------------------------------------------------------------
        else
          where (T .le. tfreeze - 20. .and. tmp1 .gt. 0.)
            snow_cld = snow_cld + qa_mean*tmp1*deltpg*inv_dtcloud
            snow_clr = snow_clr + (1. - qa_mean)*tmp1*deltpg*inv_dtcloud
            a_snow_cld = qa_mean
            a_snow_clr = 1. - qa_mean
            ST  = ST + hls*tmp1/cp_air
          end where
          where (T .gt. tfreeze - 20. .and. tmp1 .gt. 0.)    
            rain_cld = rain_cld + qa_mean*tmp1*deltpg*inv_dtcloud
            rain_clr = rain_clr + (1. - qa_mean)*tmp1*deltpg*inv_dtcloud
            a_rain_cld = qa_mean
            a_rain_clr = 1. - qa_mean
            ST  = ST + hlv*tmp1/cp_air
          end where

!-------------------------------------------------------------------------
!    save adjustment diagnostics.
!-------------------------------------------------------------------------
            where (T .le. tfreeze - 20.)
              sum_ice_adj(:,:,:) = tmp1*inv_dtcloud
            endwhere
          if (diag_id%liq_adj + diag_id%liq_adj_col +   &
              diag_id%ice_adj + diag_id%ice_adj_col  > 0) then       
            where (T .le. tfreeze - 20.)
              diag_4d(:,:,:,diag_pt%ice_adj) = tmp1*inv_dtcloud
            elsewhere
              diag_4d(:,:,:,diag_pt%liq_adj) = tmp1*inv_dtcloud
            endwhere
          end if
        end if !super choice
      end if !for Nml%do_pdf_clouds

!------------------------------------------------------------------------
!    complete calculation of the adjustment part of the destruction 
!    diagnostics.
!------------------------------------------------------------------------
      if  (diag_id%qadt_destr + diag_id%qa_destr_col > 0)  &
           diag_4d(:,:,:,diag_pt%qadt_destr) =   &
                      diag_4d(:,:,:,diag_pt%qadt_destr) + SA*inv_dtcloud
      if  (diag_id%qldt_destr + diag_id%ql_destr_col > 0)   &
           diag_4d(:,:,:,diag_pt%qldt_destr) =    &
                      diag_4d(:,:,:,diag_pt%qldt_destr) + SL*inv_dtcloud
      if  (diag_id%qidt_destr + diag_id%qi_destr_col > 0)   &
           diag_4d(:,:,:,diag_pt%qidt_destr) =    &
                      diag_4d(:,:,:,diag_pt%qidt_destr) + SI*inv_dtcloud
      if  (diag_id%qndt_destr + diag_id%qn_destr_col > 0 .and.    &
                                                  Nml%do_liq_num )   &
           diag_4d(:,:,:,diag_pt%qndt_destr) =    &
                      diag_4d(:,:,:,diag_pt%qndt_destr) + SN*inv_dtcloud

!-----------------------------------------------------------------------
!    final clean up to remove numerical noise
!----------------------------------------------------------------------
      ql_upd = ql + SL
      where ( abs(ql_upd) .le. Nml%qmin .and. qv + SQ + ql_upd > 0.0 )    
        SL = -ql
        SQ = SQ + ql_upd
        ST = ST - hlv*ql_upd/cp_air
        ql_upd = 0.0
      end where
      qi_upd = qi + SI
      where ( abs(qi_upd) .le. Nml%qmin .and. qv + SQ + qi_upd > 0.0 )   
        SI = -qi
        SQ = SQ + qi_upd
        ST = ST - hls*qi_upd/cp_air
        qi_upd = 0.0
      end where
      qa_upd = qa + SA
      where ( abs(qa_upd) .le. Nml%qmin )         
        SA  = -qa
        qa_upd = 0.0
      end where

!------------------------------------------------------------------------
!    add the noise elimination contribution to the destruction diagnostics.
!------------------------------------------------------------------------
      if (diag_id%qadt_destr + diag_id%qa_destr_col > 0)    &
           diag_4d(:,:,:,diag_pt%qadt_destr) =    &
                   diag_4d(:,:,:,diag_pt%qadt_destr) - SA*inv_dtcloud
           diag_4d(:,:,:,diag_pt%qadt_destr) =    &
                   -diag_4d(:,:,:,diag_pt%qadt_destr) 
      if (diag_id%qldt_destr + diag_id%ql_destr_col > 0)    &
           diag_4d(:,:,:,diag_pt%qldt_destr) =     &
                   diag_4d(:,:,:,diag_pt%qldt_destr) - SL*inv_dtcloud
      if (diag_id%qidt_destr + diag_id%qi_destr_col > 0)    &
           diag_4d(:,:,:,diag_pt%qidt_destr) =    &
                  -( diag_4d(:,:,:,diag_pt%qidt_destr) - SI*inv_dtcloud)
      if (diag_id%qndt_destr + diag_id%qn_destr_col > 0 .and.    &
                                                    Nml%do_liq_num )  &
           diag_4d(:,:,:,diag_pt%qndt_destr) =    &
                   diag_4d(:,:,:,diag_pt%qndt_destr) - SN *inv_dtcloud
       
      if (diag_id%qldt_destr + diag_id%ql_destr_col > 0)    &
           diag_4d(:,:,:,diag_pt%qldt_destr) =     &
                   -diag_4d(:,:,:,diag_pt%qldt_destr) 
      if (diag_id%qndt_destr + diag_id%qn_destr_col > 0 .and.    &
                                                    Nml%do_liq_num )  &
           diag_4d(:,:,:,diag_pt%qndt_destr) =    &
                   -diag_4d(:,:,:,diag_pt%qndt_destr) 

!-----------------------------------------------------------------------
!    add the ice falling out from cloud to  the qidt_fall diagnostic.
!-----------------------------------------------------------------------
        if (diag_id%qidt_fall + diag_id%qi_fall_col > 0)  &
             diag_4d(:,:,:,diag_pt%qidt_fall) =    &
                      diag_4d(:,:,:,diag_pt%qidt_fall) - (snow_cld/deltpg)
        
!-----------------------------------------------------------------------
!    save output fields of profiles of total rain and snow and clear-sky
!    snow.
!-----------------------------------------------------------------------
      rain3d(:,:,1) = 0.
      snow3d(:,:,1) = 0.
      snowclr3d(:,:,1) = 0.
      do k=1,kdim 
        rain3d(:,:,k+1) = rain_clr(:,:,k) + rain_cld(:,:,k)
        snow3d(:,:,k+1) = snow_clr(:,:,k) + snow_cld(:,:,k)
        snowclr3d(:,:,k+1) = snow_clr(:,:,k)

!-----------------------------------------------------------------------
!    save rain and snow diagnostics
!-----------------------------------------------------------------------
        diag_4d_kp1(:,:,k+1,diag_pt%rain_clr) = rain_clr(:,:,k)
        diag_4d_kp1(:,:,k+1,diag_pt%rain_cld) = rain_cld(:,:,k)
        diag_4d_kp1(:,:,k+1,diag_pt%a_rain_clr) = a_rain_clr(:,:,k)
        diag_4d_kp1(:,:,k+1,diag_pt%a_rain_cld) = a_rain_cld(:,:,k)
        diag_4d_kp1(:,:,k+1,diag_pt% snow_clr) = snow_clr(:,:,k)
        diag_4d_kp1(:,:,k+1,diag_pt%snow_cld) = snow_cld(:,:,k)
        diag_4d_kp1(:,:,k+1,diag_pt%a_snow_clr) = a_snow_clr(:,:,k)
        diag_4d_kp1(:,:,k+1,diag_pt%a_snow_cld) = a_snow_cld(:,:,k)
        diag_4d_kp1(:,:,k+1,diag_pt%a_precip_clr) =   &
                                  max(a_rain_clr(:,:,k),a_snow_clr(:,:,k))
        diag_4d_kp1(:,:,k+1,diag_pt%a_precip_cld) =    &
                                  max(a_rain_cld(:,:,k),a_snow_cld(:,:,k))
      end do

!-----------------------------------------------------------------------
!    put rain and ice fluxes into surfrain and surfsnow if the grid point 
!    is at the bottom of a column. If MASK is not present then this code 
!    is executed only if k .eq. kdim.
!    IF MASK is present some grid points may be beneath ground. If a given
!    grid point is at the bottom of the column then the surface values of
!    rain and snow must be created. Also if the MASK is present then the 
!    code forces all tendencies below ground to be zero. Note that
!    MASK = 1. equals above ground point, MASK = 0. equals below ground 
!    point.
!------------------------------------------------------------------------
      if ( mask_present ) then

!------------------------------------------------------------------------
!    zero out all tendencies below ground
!------------------------------------------------------------------------
        ST = MASK*ST
        SQ = MASK*SQ
        SL = MASK*SL
        SI = MASK*SI
        SA = MASK*SA
        rain3d(:,:,1:kdim) = mask*rain3d(:,:,1:kdim)
        snow3d(:,:,1:kdim) = mask*snow3d(:,:,1:kdim)
        rain3d(:,:,kdim+1) = mask(:,:,kdim)*rain3d(:,:,kdim+1)
        snow3d(:,:,kdim+1) = mask(:,:,kdim)*snow3d(:,:,kdim+1)
        if (Nml%do_liq_num) then
          SN = MASK*SN
        endif

!------------------------------------------------------------------------
!    determine the true bottom of columns which contain some below-ground
!    points. grab that value as surfrain or surfsnow and set the other
!    rain and asnow arrays at that level to 0.
!------------------------------------------------------------------------
        DO k=1,kdim-1
          where(MASK(:,:,k) .eq. 1. .and. MASK(:,:,k+1) .eq. 0.)
            surfrain(:,:) = dtcloud*rain3d(:,:,k+1)
            surfsnow(:,:) = dtcloud*snow3d(:,:,k+1)
          end where
        END DO

!-----------------------------------------------------------------------
!    define the precip at the surface in those columns with no below-ground
!    points.
!-----------------------------------------------------------------------
        where(MASK(:,:,kdim) .eq. 1.)
          surfrain(:,:) = dtcloud*rain3d(:,:,kdim+1)
          surfsnow(:,:) = dtcloud*snow3d(:,:,kdim+1)
        end where
      else

!------------------------------------------------------------------------
!    when mask not present, the surface precip is defined as the values
!    at level kdim.
!------------------------------------------------------------------------
        surfrain(:,:) = dtcloud*rain3d(:,:,kdim+1)
        surfsnow(:,:) = dtcloud*snow3d(:,:,kdim+1)
      end if 
                  
!-----------------------------------------------------------------------
!  Used for BC aerosol in-cloud scavenging:
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            qldt_sum = sum_berg(i,j,k) + sum_rime(i,j,k) +   &
                       sum_ice_adj(i,j,k) + sum_cond(i,j,k) +   &
                       sum_freeze(i,j,k)
            if (qldt_sum > 0.)  then
              f_snow_berg(i,j,k) = (sum_berg(i,j,k) + sum_freeze(i,j,k) + &
                          sum_ice_adj(i,j,k) + sum_cond(i,j,k))/qldt_sum 
            else
              f_snow_berg(i,j,k) = 0.
            endif
          end do
        end do
      end do
        



END SUBROUTINE rotstayn_klein_microp


!########################################################################

SUBROUTINE rotstayn_klein_microp_end

!-----------------------------------------------------------------------
!    call destructors for modules used here.
!-----------------------------------------------------------------------
      call polysvp_end

!-----------------------------------------------------------------------
!    mark the module as uninitialized.
!-----------------------------------------------------------------------
      module_is_initialized = .false.

!----------------------------------------------------------------------

END SUBROUTINE rotstayn_klein_microp_end


!########################################################################

subroutine cloud_clear_xfer (k, Nml, cloud_generator_on, overlap,  &
                             tmp3, qa_mean, a_clr, a_cld, clr, cld)

!---------------------------------------------------------------------
integer,                intent(in)    :: k, overlap
type(strat_nml_type),   intent(in)    :: Nml
logical,                intent(in)    :: cloud_generator_on
real, dimension(:,:,:), intent(in)    :: tmp3, qa_mean
real, dimension(:,:,:), intent(inout) :: a_clr, a_cld, clr, cld

!------------------------------------------------------------------------
!---local variables-----------------------------------------------------

      real, dimension (size(tmp3,1), size(tmp3,2)) ::    &
                               da_cld2clr, da_clr2cld, dprec_cld2clr,  &
                                                 dprec_clr2cld, tmp1, tmp2
      integer :: km1

!------------------------------------------------------------------------
      km1 = MAX(k-1,1)
      if (overlap .eq. 1)  then
        da_cld2clr = min(a_cld(:,:,km1), max(0.,   &
                                         a_cld(:,:,km1) - qa_mean(:,:,k)))
      end if        
      if (overlap .eq. 2)  then    
        da_cld2clr = min(a_cld(:,:,km1), max(0.,  &
                                    a_cld(:,:,km1)*(1. - qa_mean(:,:,k))))
      end if
       
      if (cloud_generator_on) then
        tmp1 = min(a_cld(:,:,km1), max(0.,     &
                                      a_cld(:,:,km1) - qa_mean(:,:,k))) 
        tmp2 = min(a_cld(:,:,km1), max(0.,    &
                                     a_cld(:,:,km1)*(1. - qa_mean(:,:,k))))
        da_cld2clr = min(a_cld(:,:,km1), max(0.,   &
                             tmp3(:,:,k)*tmp1 + (1. - tmp3(:,:,k))*tmp2)) 
      end if      

!-------------------------------------------------------------------------
!    compute clear to cloud transfer
!-------------------------------------------------------------------------
      if (overlap .eq. 1)                                            &
        da_clr2cld = min(max(qa_mean(:,:,k) - qa_mean(:,:,km1), 0.), &
                                                           a_clr(:,:,km1))
      if (overlap .eq. 2)                                            &
        da_clr2cld = min(max( a_clr(:,:,km1)*qa_mean(:,:,k), 0.), &
                                                           a_clr(:,:,km1))
      if (cloud_generator_on) then
        tmp1 = min(max(qa_mean(:,:,k) - qa_mean(:,:,km1), 0.),  &
                                                          a_clr(:,:,km1)) 
        tmp2 = min(max( a_clr(:,:,km1)*qa_mean(:,:,k), 0.),   &
                                                           a_clr(:,:,km1))
        da_clr2cld = min(a_clr(:,:,km1), max(0.,   &
                             tmp3(:,:,k)*tmp1 + (1. - tmp3(:,:,k))*tmp2))
      end if      
        
!------------------------------------------------------------------------
!    calculate precipitation transfers
!------------------------------------------------------------------------
      dprec_cld2clr = cld(:,:,km1)*(da_cld2clr/max(a_cld(:,:,km1),  &
                                                                Nml%qmin))
      dprec_clr2cld = clr(:,:,km1)*(da_clr2cld(:,:)/max(a_clr(:,:,km1), &
                                                                 Nml%qmin))
        
!------------------------------------------------------------------------
!    add in transfers
!------------------------------------------------------------------------
      a_clr(:,:,k) = a_clr(:,:,km1) + da_cld2clr - da_clr2cld
      a_cld(:,:,k) = a_cld(:,:,km1) - da_cld2clr + da_clr2cld
      clr(:,:,k) = clr(:,:,km1) + dprec_cld2clr - dprec_clr2cld   
      cld(:,:,k) = cld(:,:,km1) - dprec_cld2clr + dprec_clr2cld

!----------------------------------------------------------------------


end subroutine cloud_clear_xfer


!#######################################################################



END MODULE rotstayn_klein_mp_mod

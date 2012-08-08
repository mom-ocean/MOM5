
!VERSION NUMBER:
!  $Id: donner_meso_k.F90,v 19.0 2012/01/06 20:08:00 fms Exp $

!module donner_meso_inter_mod

!#include "donner_meso_interfaces.h"

!end module donner_meso_inter_mod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine don_m_meso_effects_k    &
         (me, nlev_lsm, nlev_hires, ntr, diag_unit, debug_ijt, Param, Nml,&
          pfull_c, temp_c, mixing_ratio_c, phalf_c, rlsm, emsm, etsm, &
          tracers_c, ensembl_cond, ensmbl_precip, pb, plzb_c, pt_ens, &
          ampta1, ensembl_anvil_cond_liq, ensembl_anvil_cond_liq_frz, &
          ensembl_anvil_cond_ice,  wtp, qtmes, meso_frz_intg_sum,  &
          anvil_precip_melt, meso_cloud_area, cmus_tot, dmeml,  &
          emds_liq, emds_ice, emes_liq, emes_ice,  wmms, wmps, &
          umeml, temptr, tmes, tmes_up, tmes_dn, mrmes, mrmes_up, &
          mrmes_dn, emdi, pmd, pztm, pzm, meso_precip, ermesg, error)

!-------------------------------------------------------------------
!    subroutine don_m_meso_effects_k obtains the mesoscale effects
!    of the composited cloud ensemble on the heat, moisture and tracer 
!    budgets, producing tendency terms which are to be applied to the 
!    large-scale model equations. the scheme employed here is a variation
!    on the procedure of Leary and Houze (JAS, 1980). for more details 
!    on notation, see "Cu Closure A notes," 2/97.
!-------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_nml_type

implicit none

!-------------------------------------------------------------------
integer,                           intent(in)  :: me, nlev_lsm, nlev_hires, &
                                                  ntr, diag_unit
logical,                           intent(in)  :: debug_ijt        
type(donner_param_type),           intent(in)  :: Param
type(donner_nml_type),             intent(in)  :: Nml  
real,   dimension(nlev_lsm),       intent(in)  :: pfull_c, temp_c, &
                                                  mixing_ratio_c
real,   dimension(nlev_lsm+1),     intent(in)  :: phalf_c
real,   dimension(nlev_hires),     intent(in)  :: rlsm, emsm
real,   dimension(nlev_hires,ntr), intent(in)  :: etsm
real,   dimension(nlev_lsm,ntr),   intent(in)  :: tracers_c
logical,                           intent(in)  :: meso_frz_intg_sum
real,                              intent(in)  :: ensembl_cond,   &
                                                  ensmbl_precip, pb, &
                                                  plzb_c, pt_ens,   &
                                                  ampta1, &
                                              ensembl_anvil_cond_liq, &
                                       ensembl_anvil_cond_liq_frz, &
                                                  ensembl_anvil_cond_ice
real,   dimension(nlev_lsm,ntr),   intent(out) :: wtp, qtmes, temptr
real,   dimension(nlev_lsm),       intent(out) :: anvil_precip_melt, &
                                                  meso_cloud_area,    &
                                                  cmus_tot, dmeml, &
                                                  emds_liq, emds_ice, &
                                                        wmms, wmps,   &
                                                  emes_liq, emes_ice, &
                                                  umeml, tmes, mrmes, &
                                                  mrmes_up, mrmes_dn, &
                                                  tmes_up, tmes_dn
real,                              intent(out) :: emdi,           &
                                                  pmd, pztm, pzm, &
                                                  meso_precip
character(len=*),                  intent(out) :: ermesg
integer,                           intent(out) :: error

!----------------------------------------------------------------------
!   intent(in) variables:
!
!       pfull_c      large-scale model pressure full levels [ Pa ]
!       phalf_c      large-scale model pressure half levels [ Pa ]
!       temp_c       large-scale model temperature profile [ deg K ]
!       mixing_ratio_c  
!                    large-scale model mixing ratio profile
!                    [ kg(h2o) / kg(air) ]
!       rlsm         cloud model condensation profile summed over
!                    cloud ensemble
!                    [ kg(h2o) / kg(air) / sec ]
!       emsm         cloud model moisture flux convergence summed over 
!                    the cloud ensemble
!                    [ kg(h2o) / kg(air) / sec ]
!       etsm         cloud model tracer flux convergence summed over
!                    the cloud ensemble 
!                    [ kg(tracer) / kg(air) / sec ]
!       tracers_c    large-scale model tracer mixing ratio profiles
!                    [ kg(tracer) /kg(air) ]
!       ensmbl_cond  total ensemble condensation integral
!                    [ mm / day ]
!       ensmbl_precip   total ensemble precipitation integral
!                    [ mm / day ]
!       ps           surface pressure [ Pa ]
!       pb           cloud-base pressure [ Pa ]
!       plzb_c       level of zero buoyancy [ Pa ]
!       pt_ens       cloud-top pressure [ Pa ]
!       ampta1       fractional area of mesoscale anvil
!                    [ dimensionless ]
!       ensembl_anvil_cond 
!                    condensed water transferred from cells to anvil 
!                    [ mm / day ]
!       debug_ijt    is this a diagnostics column ?
!       diag_unit    output unit number for this diagnostics column
!
!  output variables:
! 
!       meso_cloud_area 
!               fractional mesoscale area, normalized by
!               a(1,p_b) at resolution of GCM
!       meso_precip
!       cmu     water mass condensed in mesoscale updraft
!               (g/kg/day) (normalized by a(1,p_b))
!       cmui    vertical integral of mesoscale-updraft deposition
!               (kg(H2O)/((m**2)*sec) 
!       dmeml   mass flux in mesoscale downdraft (kg/((m**2) s))
!               (normalized by a(1,p_b)) (index 1 at atmosphere bottom)
!               (resolution of GCM)
!       emds    water mass evaporated in mesoscale
!               downdraft (g/kg/day) (normalized by a(1,p_b))
!       emdi    vertical integral of mesoscale-downdraft sublimation
!               (mm/d)
!       emes    water mass evaporated from mesoscale
!               updraft (g/kg/day) (normalized by a(1,p_b))
!       emei    vertical integral of mesoscale-updraft sublimation
!               (kg(h2O)/((m**2)*sec)
!       pmd     pressure at top of mesoscale downdraft (Pa)
!       pztm    pressure at top of mesoscale updraft (Pa)
!       wmms    water vapor removal by condensation of
!               cell vapor source (g/kg/day) (normalized by a(1,p_b))
!       wmps    water vapor redistributed from cell vapor source
!               (g/kg/day) (normalized by a(1,p_b))
!       wtp     tracer redistributed by mesoscale processes
!               (kg/kg/s) (normalized by a(1,p_b))
!       anvil_precip_melt     melting of ice in mesoscale updraft-
!               equivalent (g/kg/day)-which falls as meso sfc precip
!               (normalized by a(1,p_b))
!       tmes    temperature tendency due to mesoscale entropy-flux-
!               convergence (K/day) (normalized by a(1,p_b))
!       mrmes    moisture tendency due to mesoscale moisture-flux
!               convergence (g/kg/day) (normalized by a(1,p_b))
!       qtmes   tracer tendency due to mesoscale tracer-flux
!               convergence (kg/kg/s) (normalized by a(1,p_b))
!       umeml   mass flux in mesoscale updraft (kg/((m**2) s))
!               (normalized by a(1,p_b)) (index 1 at atmosphere bottom)
!               (resolution of GCM)
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      real, dimension (nlev_lsm)     ::  cmu
      real, dimension (nlev_lsm)     ::  out
      real, dimension (nlev_hires)   ::  p_hires
      real                           ::  alp, hfmin, cmui, qtmesum, dp,&
                                         available_condensate, &
                                         available_condensate_liq, &
                                         available_condensate_ice
      real                           ::  emdi_liq, emdi_ice
      integer                        ::  k, kcont, itrop
      real  :: intgl_lo, intgl_hi
      real          :: p2, ptrop

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

      dp = Param%dp_of_cloud_model

!--------------------------------------------------------------------
!    define the pressure at the melting level (p2).
!--------------------------------------------------------------------
      p2 = -10.
      do k=1,nlev_lsm-1
        if ((temp_c(k) >= Param%kelvin) .and.   &
             (temp_c(k+1) <= Param%kelvin))  then
          p2 = phalf_c(k+1)
          exit
        end if
      end do

!---------------------------------------------------------------------
!    if in diagnostics column, output message indicating that sub-
!    routine meso_effects has been entered.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a)') 'in meens: entering meens'
      endif

!--------------------------------------------------------------------
!    define the pressure at the top of the mesoscale updraft (pztm) to 
!    be the pressure at the zero buyancy level, unless the cloud top is
!    above 100 hPa, in which case pztm is set to be one level above the 
!    level of zero buoyancy.  previously pztm was restricted to be  >=
!    100 hPa, cf Ackerman et al (JAS,1988), unless pt_ens <= 10kPa. 
!    result was that stratospheric water vapor was transported too high 
!    in AM2p9 with this pztm, so the constraint was changed to pztm >= 
!    plzb_c + dp
!--------------------------------------------------------------------
      if ((pt_ens + dp) >= 10.e03)  then
        pztm = plzb_c
      else
        pztm = plzb_c + dp
      endif

      if (Nml%limit_pztm_to_tropo) then
        call find_tropopause (nlev_lsm, temp_c, pfull_c, ptrop, itrop)
        pztm = MAX (pztm, ptrop)
      endif

!---------------------------------------------------------------------
!    if in diagnostics column, output the pressure at top of meso-
!    scale circulation (pztm) and the precipitation efficiency 
!    (ensmbl_precip/ensembl_cond).
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a,  e20.12)') 'in meens: pztm = ',pztm 
        write (diag_unit, '(a, e20.12)') 'in meens: gnu= ',   &
                                             ensmbl_precip/ensembl_cond
      endif

!---------------------------------------------------------------------
!    define the pressure at the vertical grid levels of the cloud model
!    grid.
!---------------------------------------------------------------------
      do k=1,nlev_hires        
        p_hires(k) = pb + (k-1)*dp
      end do

!---------------------------------------------------------------------
!    call subroutine meso_updraft to define the needed output fields 
!    associated with the mesoscale updraft.
!---------------------------------------------------------------------
      call don_m_meso_updraft_k   &
           (nlev_lsm, nlev_hires, ntr, diag_unit, debug_ijt, Param, &
            p_hires, rlsm, emsm, etsm, pfull_c,  &
             temp_c, mixing_ratio_c, phalf_c, tracers_c,  &
              pb, pt_ens, ampta1, dp, pztm,  wtp, &
                  qtmes, cmu, wmms, wmps, temptr, tmes_up, mrmes_up,   &
                  meso_cloud_area, umeml,&
                  alp, pzm, hfmin, cmui, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

      if (Nml%frc_internal_enthalpy_conserv) then
!-----------------------------------------------------------------------
!    call don_u_set_column_integral_k to adjust the tmes_up
!    profile below cloud base so that the desired integral value is
!    obtained.
!-----------------------------------------------------------------------
         call don_u_set_column_integral_k    &
              (nlev_lsm, tmes_up   , pb, &
               phalf_c(1), 0.0, phalf_c , intgl_hi,     &
               intgl_lo, out, ermesg, error)

!---------------------------------------------------------------------
!    if column diagnostics are desired, output the integrals and 
!    profiles, both before and after the adjustment to the desired value        .
!---------------------------------------------------------------------
        if (debug_ijt) then
           write (diag_unit, '(a, e20.12)')  &
                   'in set_col_integral: tmes_up column(in)= ',intgl_hi
           write (diag_unit, '(a, e20.12)')  &
                   'in set_col_integral: tmes_up column(out)= ',intgl_lo
           do k=1,nlev_lsm
            if (tmes_up(k)       /= out(k)) then
               write (diag_unit, '(a, i4, 2e20.12)') &
               'in set_col_integral: k,tmesup(in), tmesup(out)= ', k,  &
                     tmes_up(k)      , out(k)
            endif
          end do
        endif
 
!---------------------------------------------------------------------
!    define the adjusted output profile by removing conservation_factor.
!---------------------------------------------------------------------
       tmes_up(:) = out(:)       
    endif

!---------------------------------------------------------------------
!    call subroutine meso_downdraft to define the needed output fields 
!    associated with the mesoscale downdraft.
!---------------------------------------------------------------------
      call don_m_meso_downdraft_k  &
           (nlev_lsm, nlev_hires, diag_unit, debug_ijt, Param, Nml,  &
            p_hires, pfull_c, temp_c, mixing_ratio_c, phalf_c, pb, &
            ampta1, dp, pztm, pzm, alp, hfmin, pmd, tmes_dn,  &
            mrmes_dn, dmeml, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

      if (Nml%frc_internal_enthalpy_conserv) then
!-----------------------------------------------------------------------
!    call don_u_set_column_integral_k to adjust the tmes_dn
!    profile below cloud base so that the desired integral value is
!    obtained.
!-----------------------------------------------------------------------
         call don_u_set_column_integral_k    &
              (nlev_lsm, tmes_dn   , pb, &
               phalf_c(1), 0.0, phalf_c , intgl_hi,     &
               intgl_lo, out, ermesg, error)

!---------------------------------------------------------------------
!    if column diagnostics are desired, output the integrals and 
!    profiles, both before and after the adjustment to the desired value        .
!---------------------------------------------------------------------
        if (debug_ijt) then
           write (diag_unit, '(a, e20.12)')  &
                    'in set_col_integral: tmes_dn column(in)= ',intgl_hi
           write (diag_unit, '(a, e20.12)')  &
                  'in set_col_integral: tmes_dn column(out)= ',intgl_lo
           do k=1,nlev_lsm
            if (tmes_dn(k) /= out(k)) then
               write (diag_unit, '(a, i4, 2e20.12)') &
               'in set_col_integral: k,tmesdn(in), tmesdn(out)= ', k,  &
                     tmes_dn(k)      , out(k)
            endif
          end do
        endif
 
!---------------------------------------------------------------------
!    define the adjusted output profile by removing conservation_factor.
!---------------------------------------------------------------------
       tmes_dn(:) = out(:)       
     endif

!---------------------------------------------------------------------
!    combine the heating and moistening effects from the updraft and
!    downdraft to obtain the total mesoscale effect on the large-scale
!    model temperature and water vapor mixing ratio(?) equations.
!---------------------------------------------------------------------
      tmes = (tmes_up + tmes_dn)*86400.
      tmes_up = tmes_up*86400.
      tmes_dn = tmes_dn*86400.
      mrmes = (mrmes_up + mrmes_dn)*8.64e07
      mrmes_up = mrmes_up*8.64e07
      mrmes_dn = mrmes_dn*8.64e07

!---------------------------------------------------------------------
!    if in a diagnostics column, output the entropy (tmes) and
!    mixing ratio (mrmes) tendencies due to the mesoscale
!    updraft and downdraft.
!---------------------------------------------------------------------
      do k=1,nlev_lsm              
        if (debug_ijt) then
          if (tmes(k) /= 0.0) then
            write (diag_unit, '(a, i4, f19.10, f20.14, 2e20.12)')   &
                    'in meens: jk,pr,tmes,tmes_u, tmes_d,= ', &
                     k, pfull_c(k), tmes(k)/86400., tmes_up(k)/86400., &
                     tmes_dn(k)/86400.
            write (diag_unit, '(a, i4, f19.10, f20.14, 3e20.12)')   &
                    'in meens: jk,pr,mrmes,mrmes_u, mrmes_d= ', &
                     k, pfull_c(k), mrmes(k)/8.64e07,  &
                      mrmes_up(k)/8.64e07, mrmes_dn(k)/8.64e07
          endif
        endif
      end do

!---------------------------------------------------------------------
!    define the column anvil precip (meso_precip) as the precipitation
!    efficiency times the available condensate in the anvil, which is 
!    made up of the deposition in the updraft (cmui) and the condensate
!    transferred from the cells to the anvil (ensembl_anvil_cond). 
!---------------------------------------------------------------------
      available_condensate = cmui + ensembl_anvil_cond_liq + &
                                ensembl_anvil_cond_liq_frz + &
                             ensembl_anvil_cond_ice
! precip from _liq takes hlv with it; precip from _ice takes hls
! with it
      if ( p2 == -10. .or. p2 > pb .or. p2 < pt_ens) then
        if ( .not. meso_frz_intg_sum ) then
!   this implies no melting of precip; cmui and _liq don't freeze.
         available_condensate_liq =  cmui + ensembl_anvil_cond_liq 
         available_condensate_ice =         &
                                ensembl_anvil_cond_liq_frz + &
                             ensembl_anvil_cond_ice
        else
         available_condensate_liq =  0.0                           
         available_condensate_ice =    cmui + ensembl_anvil_cond_liq + & 
                                ensembl_anvil_cond_liq_frz + &
                             ensembl_anvil_cond_ice
        endif
      else
!    all condensate will melt before leaving
      available_condensate_ice = 0.0                               
      available_condensate_liq = cmui + ensembl_anvil_cond_liq + &
                                ensembl_anvil_cond_liq_frz + &
                             ensembl_anvil_cond_ice
      endif

      meso_precip = Param%anvil_precip_efficiency*available_condensate

!---------------------------------------------------------------------
!    if in a diagnostics column, output the total mesoscale-supplied
!    condensate (condensation plus deposition), the cell provided 
!    condensate (ensembl_anvil_cond),  the mesoscale precipitation 
!    (meso_precip) and the cell-scale precipitation (ensmbl_precip).
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 4e20.12)')  &
                     'in meens: cmui,ca (liq,frxliq,ice)=', cmui, &
                 ensembl_anvil_cond_liq, ensembl_anvil_cond_liq_frz,  &
                             ensembl_anvil_cond_ice
        write (diag_unit, '(a, e20.12, a, e20.12)')  &
                     'in meens: rm= ',meso_precip,  'rc= ',ensmbl_precip
      endif

!----------------------------------------------------------------------
!    call subroutine meso_evap to define the amount of condensate that
!    is evaporated in the mesoscale updraft (emes) and mesoscale 
!    downdraft (emds).
!----------------------------------------------------------------------
      call don_m_meso_evap_k  &
           (nlev_lsm, diag_unit, debug_ijt, Param,  &
            available_condensate, available_condensate_liq,  &
            available_condensate_ice, pzm, pztm, phalf_c, emdi_liq, &
            emdi_ice, emds_liq, emds_ice, emes_liq, emes_ice, ermesg, error)
 
      emdi = emdi_liq + emdi_ice

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    call subroutine meso_melt to distribute the melting of precipitat-
!    ing anvil ice within the column (anvil_precip_melt).
!---------------------------------------------------------------------
      call don_m_meso_melt_k   &
           (nlev_lsm, diag_unit, debug_ijt, Param, temp_c, phalf_c, &
            pztm, meso_precip, pb, anvil_precip_melt, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!--------------------------------------------------------------------
!    define cmus_tot   as the profile of total condensate source to the
!    large-scale flow from the mesoscale circulation; the sum of the
!    water mass condensed in the mesoscale updraft plus the vapor
!    transferred from cell to mesoscale and then condensed.
!--------------------------------------------------------------------
      do k=1,nlev_lsm            
        cmus_tot(k) = cmu(k) - wmms(k)
      end do

!---------------------------------------------------------------------
!    if in a diagnostics column, output the profiles of tracer tranfer-
!    red from cells to mesoscale circulation (wtp), mesoscale tracer-
!    flux convergence (qtmes), and cell-scale tracer flux convergence 
!    (qtren). also output the  column integral of the mesoscale 
!    tracer-flux convergence (qtmesum).
!---------------------------------------------------------------------
      if (debug_ijt) then
        qtmesum = 0.
        do k=1,nlev_lsm             
          do kcont=1,ntr           
            write (diag_unit, '(a, 2i4, f19.10, e20.12)')  &
                      'in mulsub: jk, pr,wtp= ',k, kcont,  &
                          pfull_c(k), wtp(k,kcont)
            write (diag_unit, '(a, 2i4, f19.10, e20.12)')  &
                     'in mulsub: jk, pr,qtmes= ', k, kcont,         &
                            pfull_c(k),  qtmes(k,kcont)
            qtmesum = qtmesum + qtmes(k,kcont)*  &
                      (phalf_c(k) - phalf_c(k+1))
            write (diag_unit, '(a, i4, e20.12)')  &
                        'in mulsub: jk,qtmesum= ', k, qtmesum
          end do
        end do
      endif

!--------------------------------------------------------------------


end subroutine don_m_meso_effects_k 


!#######################################################################

subroutine don_m_meso_updraft_k    &
         (nlev_lsm, nlev_hires, ntr, diag_unit, debug_ijt, Param,  &
          p_hires, rlsm, emsm, etsm, pfull_c, temp_c, mixing_ratio_c, &
          phalf_c, tracers_c, pb, pt_ens, ampta1, dp, pztm, wtp, &
          qtmes, cmu, wmms, wmps, temptr, tmes_up, mrmes_up,  &
          meso_cloud_area, umeml, alp, pzm, hfmin, cmui, ermesg, error)

!-------------------------------------------------------------------
!    subroutine meens computes the mesoscale effects of the composited
!    cloud ensemble on the heat, moisture and tracer budgets, producing
!    tendency terms which are to be applied to the large-scale model.
!    scheme employed here is a variation on procedure of Leary and 
!    Houze (JAS, 1980). for more details on notation, see 
!    "Cu Closure A notes," 2/97.
!-------------------------------------------------------------------

use donner_types_mod, only : donner_param_type
use sat_vapor_pres_k_mod, only: compute_mrs_k

implicit none

!-------------------------------------------------------------------
integer,                         intent(in)  :: nlev_lsm, nlev_hires, ntr
integer,                         intent(in)  :: diag_unit
logical,                         intent(in)  :: debug_ijt
type(donner_param_type),         intent(in)  :: Param
real,   dimension(nlev_hires),   intent(in)  :: p_hires, rlsm, emsm
real,   dimension(nlev_hires,ntr),                          &
                                 intent(in)  :: etsm
real,   dimension(nlev_lsm),     intent(in)  :: pfull_c, temp_c,    &
                                                mixing_ratio_c
real,   dimension(nlev_lsm+1),   intent(in)  :: phalf_c
real,   dimension(nlev_lsm,ntr), intent(in)  :: tracers_c
real,                            intent(in)  :: pb, pt_ens, ampta1,   &
                                                dp, pztm
real,   dimension(nlev_lsm,ntr), intent(out) :: wtp, qtmes, temptr
real,   dimension(nlev_lsm),     intent(out) :: cmu, wmms, wmps, &
                                                tmes_up, mrmes_up, &
                                                meso_cloud_area, umeml
real,                            intent(out) :: alp, pzm, hfmin, cmui
character(len=128),              intent(out) :: ermesg
integer,                         intent(out) :: error

!---------------------------------------------------------------------
!   local variables:



      real, dimension (nlev_hires)       :: wmhr, cumh
      real, dimension (nlev_lsm)         :: omv, tempq, owm, tempqa
      real, dimension(nlev_lsm,ntr)      :: otm
      real, dimension(nlev_hires, ntr)   :: wthr
      real, dimension(ntr)               :: q1t


      real      ::  cmfhr, pc1, pc2, omer, pctm, q1, q4, mrsat, &
                    q3, anv, qref, pp, pm, qprip, qprim, eqfp, eqfm, &
                    qmu, hflux, pfmin, owms, wpc, wmc, ta, te, tep, tmu,&
                    qtprip, qtprim, eqtfp, eqtfm, rintsum, rintsum2
      logical   :: do_donner_tracer
      integer   :: ncc, ncztm
      integer   :: kcont, kk
      integer   :: jk, i, jsave, jkm, jkp, k, nbad

!-----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

      if (ntr > 0) then
        do_donner_tracer = .true.
      else
        do_donner_tracer = .false.
      endif

      do i=1,nlev_hires
        if (p_hires(i) < pt_ens) then
          ncc = i
          exit
        endif
      end do
      do i=1,nlev_hires
        if (p_hires(i) < pztm) then
          ncztm = i + 1
          exit
        endif
      end do

      do kcont=1,ntr
        wtp(:,kcont) = 0.
        qtmes(:,kcont) = 0.
        temptr(:,kcont) = tracers_c(:,kcont)
      end do
      tmes_up(:) = 0.
      mrmes_up(:) = 0.
      cmu = 0.
      wmms = 0.
      wmps = 0.
      tempq(:) = mixing_ratio_c(:)
      tempqa(:) = mixing_ratio_c(:)

!----------------------------------------------------------------------
!    initialize the pressure at the base of the mesoscale circulation
!    (pzm).
!----------------------------------------------------------------------
      pzm = 0.

!----------------------------------------------------------------------
!    define the vertical profile of the rate at which water vapor is
!    made available to the mesoscale circulation by the convective 
!    updrafts on the cloud model grid (wmhr). if vapor is being made 
!    available, determine if there is also a vertical flux convergence 
!    of tracer; if so, define the rate at which tracer is being made
!    available to the mesoscale circulation (wthr). define the pressure
!    at the base of the mesoscale circulation (pzm) as the pressure at 
!    the lowest cloud model level where the convective updrafts are 
!    supplying condensate to the mesoscale circulation.
!----------------------------------------------------------------------
      do k=1,nlev_hires
        cmfhr = -rlsm(k) + emsm(k)
        if (cmfhr > 0.) then
          wmhr(k) = -cmfhr
          if (do_donner_tracer) then
            do kcont=1,ntr
              if (etsm(k,kcont) > 0.) then
                wthr(k,kcont) = -etsm(k,kcont)
              else
                wthr(k,kcont) = 0.0               
              endif
            end do
          else
            wthr(k,:) = 0.0               
          endif
          if (pzm == 0.) then
            pzm = p_hires(k)
          endif
        else
          wmhr(k) = 0.0   
          wthr(k,:) = 0.0               
        endif
      end do

!---------------------------------------------------------------------
!    if in diagnostics column, output the profiles of condensation rate
!    (rlsm), water vapor flux convergence (emsm) and water vapor 
!    supplied to the mesoscale (wmhr) on the cloud model grid.
!---------------------------------------------------------------------
      do k=1,nlev_hires
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 2e20.12)')  &
                     'in meens: i,rlhr,emfhr= ',k,rlsm(k),emsm(k)
          write (diag_unit, '(a, i4, e20.12)')  &
                      'in meens: i,wmhr= ',k,wmhr(k)
        endif
      end do

      if (debug_ijt) then
        write (diag_unit, '(a, i4, e20.12)')  &
                        'in meens: ncc+1, pt', ncc+1, pt_ens
        do k=1,ncc+1
          write (diag_unit, '(a, i4, e20.12)')  &
                         'in meens: k,p_hi= ', k, p_hires(k)
        end do
        do k=1,nlev_lsm+1         
          write (diag_unit, '(a, i4, e20.12)')  &
                      'in meens: k,p_lo= ', k, phalf_c(k)
        end do
      endif

!---------------------------------------------------------------------
!    convert the vertical profile of vapor made available to the meso-
!    scale from the updraft to the large-scale model grid (output var-
!    iable is owm). if tracers are being transported by donner conv-
!    ection, convert the vertical profile of tracer made available to 
!    the mesoscale from the updraft to the large-scale model grid 
!    (output variable is otm). 
!---------------------------------------------------------------------
      call don_u_map_hires_c_to_lores_c_k &
           (nlev_lsm, nlev_hires, wmhr, p_hires, pt_ens + dp, phalf_c,&
            owm, rintsum, rintsum2, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)')  &
            'in meens: rintsum(owm) =', rintsum, rintsum2
        call don_u_compare_integrals_k  &
             (rintsum, rintsum2, diag_unit, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return
      endif

      if (do_donner_tracer) then
        do kcont=1,ntr
          call don_u_map_hires_c_to_lores_c_k  &
               (nlev_lsm, nlev_hires, wthr (:,kcont), p_hires,  &
                pt_ens + dp, phalf_c, otm(:,kcont), rintsum,   &
                rintsum2, ermesg, error) 
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return

          if (debug_ijt) then
            write (diag_unit, '(a, 2e20.12)')  &
                 'in meens: rintsum(otm) =', rintsum, rintsum2
            call don_u_compare_integrals_k  &
                 (rintsum, rintsum2, diag_unit, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
            if (error /= 0 ) return
          endif
        end do
      endif

!----------------------------------------------------------------------
!    adjust the value for pressure at base of mesocscale circulation,
!    if necessary.
!----------------------------------------------------------------------
      if (pzm == 0.) pzm = pt_ens
      if (pzm <= pztm - dp) pzm = pztm - dp

!---------------------------------------------------------------------
!    if in diagnostics column, output the pressure at the base of the
!    mesoscale circulation (pzm), and the vertical profile of vapor 
!    supplied to the mesoscale by the updraft on the large-scale model
!    grid (owm).
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, f19.10)') 'in meens: pzm= ',pzm
        do k=1,nlev_lsm
          write (diag_unit, '(a, i4, e20.12)')  &
                                 'in meens: jk,owm= ',k,owm(k)
        end do
      endif

!---------------------------------------------------------------------
!    march up the column, determining the redistribution of the cumulus-
!    updraft-supplied vapor by the mesoscale updraft.
!---------------------------------------------------------------------
      do k=1,nlev_lsm

!---------------------------------------------------------------------
!    if there is  vapor being supplied to the mesoscale by the cumulus
!    updraft at this level, determine the pressure depth over which the
!    mesoscale updraft will distribute that vapor over the lifetime of
!    the mesoscale circulation.
!---------------------------------------------------------------------
        if (owm(k) < 0.) then     

!---------------------------------------------------------------------
!    define the bottom (pc1) and top (pc2) of the current layer. deter-
!    mine the pressure level to which air in this layer will reach when
!    moving at the appropriate mesoscale updraft velocity for the dur-
!    ation of the mesoscale circulation (pctm). this level is limited to
!    be no higher than the top of the mesoscale circulation; if it is 
!    calculated to be higher, redefine the mesoscale updraft velocity 
!    for this layer so that the air in this layer will reach only to
!    the mesoscale circulation top, and no higher.
!---------------------------------------------------------------------
          pc1 = phalf_c(k)
          pc2 = phalf_c(k+1)
          pctm = pc2 + Param%meso_ref_omega*Param%meso_lifetime
          if (pctm <= pztm) then
            omer = (pztm - pc2)/Param%meso_lifetime
            pctm = pc2 + omer*Param%meso_lifetime
          else
            omer = Param%meso_ref_omega
          endif
 
!---------------------------------------------------------------------
!    define the amount of water vapor from this layer (owm(k)* 
!    (pc2 - pc1)*MESO_LIFETIME) which is to be distributed
!    uniformly between pc1 and pctm (q1).
!--------------------------------------------------------------------  
          q1 = owm(k)*(pc2 - pc1)*Param%meso_lifetime/(pc1 - pctm)
          q4 = 0.5*q1

!---------------------------------------------------------------------
!    define the amount of tracer from this layer (otm(k,kcont)* 
!    (pc2 - pc1)*meso_Lifetime) which is to be distributed
!    uniformly between pc1 and pctm (q1t).
!--------------------------------------------------------------------  
          if (do_donner_tracer) then
            do kcont=1,ntr
              q1t(kcont) = otm(k,kcont)*(pc2 - pc1)*Param%meso_lifetime/&
                           (pc1 - pctm)                     
            end do
          endif

!---------------------------------------------------------------------
!    if in diagnostics column, output the topmost pressure reached by
!    the mesoscale updraft from this layer (pctm), the top of the meso-
!    scale circulation (pztm) and the amount of water vapor supplied to
!    each layer between the current vertical level and the top of the 
!    mesoscale updraft originating here (q4).
!---------------------------------------------------------------------
          if (debug_ijt) then
            write (diag_unit, '(a, 3e20.12)')  &
                         'in meens: pctm,pztm,q4= ', pctm, pztm, q4
          endif

!---------------------------------------------------------------------
!    distribute the vapor supplied in the current layer to all layers
!    between the current location and the top of the mesoscale updraft.
!---------------------------------------------------------------------
          do kk=k,nlev_lsm

!--------------------------------------------------------------------
!    exit the loop when above the top of the mesoscale updraft. if still
!    within the mesoscale updraft originating from level k, add the 
!    contribution of water vapor being supplied to the mesoscale circ-
!    ulation at this level (kk) from the current source level (k), 
!    normalized by the anvil fractional area, to the arrays accumulating
!    these moisture sources (tempq, tempqa). these arrays will be used 
!    in the calculation of deposition in the mesoscale updraft.
!--------------------------------------------------------------------
            if (phalf_c(kk) < pctm) exit
            tempq(kk) = tempq(kk) + (q1/ampta1)
            tempqa(kk) = tempqa(kk) + (q4/ampta1)

!--------------------------------------------------------------------
!    add the rate of moisture input to the current layer kk from 
!    the current source layer k to the accumulation array (wmps). if the
!    current model layer extends beyond the top of the mesoscale 
!    updraft, pro-rate the contribution by the ratio of pressure depths.
!--------------------------------------------------------------------
            if (phalf_c(kk+1) <= pctm)  then
              wmps(kk) = wmps(kk) + (q1/Param%meso_lifetime)*  &
                        (phalf_c(kk) - pctm)/  &
                                            (phalf_c(kk) - phalf_c(kk+1))
            else
              wmps(kk) = wmps(kk) + q1/Param%meso_lifetime
            endif

!--------------------------------------------------------------------
!    add the contribution of tracer being supplied to the mesoscale 
!    circulation at this level (kk) from the current source level (k), 
!    normalized by the anvil fractional area, to the array accumulating
!    this tracer source (temptr). this array will be used in the 
!    calculation of tracer deposition in the mesoscale updraft.
!    add the rate of tracer input to the current layer kk from 
!    the current source layer k to the accumulation array (wtp). if the
!    current model layer extends beyond the top of the mesoscale 
!    updraft, pro-rate the contribution by the ratio of pressure depths.
!--------------------------------------------------------------------
            if (do_donner_tracer) then
              do kcont=1,ntr
                temptr(kk,kcont) = temptr(kk,kcont) + (q1t(kcont)/  &
                                   (2.* ampta1))
                if (phalf_c(kk+1) <= pctm) then
                  wtp(kk,kcont) = wtp(kk,kcont) +   &
                                  (q1t(kcont)/Param%meso_lifetime)*  &
                                  (phalf_c(kk)-pctm)/   &
                                              (phalf_c(kk)-phalf_c(kk+1))
                else
                  wtp(kk,kcont) = wtp(kk,kcont) +   &
                                  (q1t(kcont)/Param%meso_lifetime)
                endif
              end do
            endif
          end do

!--------------------------------------------------------------------
!    if in diagnostics column, output the moisture and tracer sources
!    to the mesoscale from the convective scale.
!--------------------------------------------------------------------
          if (debug_ijt) then
            do kk=k,nlev_lsm
              if (phalf_c(kk) < pctm) exit        
              write (diag_unit, '(a, i4, f19.10)') &
                            'in meens: jj,pr= ',kk,pfull_c(kk)
              write (diag_unit, '(a, i4, 3e20.12)')  &
                  'in meens: jj,q1,tempq,wmm= ',kk,q1,tempq(kk),wmms(kk)
              write (diag_unit, '(a, e20.12)')  &
                   'in meens: wmp= ',wmps(kk)
              write (diag_unit, '(a, i4, e20.12)')  &
                   'in meens: jj,tempqa= ',kk,tempqa(kk)
            end do
            write (diag_unit, '(a, i4, 3e20.12)')  &
                  'in meens: jk,q1,tempq,wmm= ',k,q1,tempq(k),wmms(k)
            write (diag_unit, '(a, i4, 2e20.12)')  &
                   'in meens: jk,wmp,owm= ',k,wmps(k),owm(k)
          endif
        endif ! (owm(k) < 0.)

!----------------------------------------------------------------------
!    if in diagnostics column, output the profile of moisture made
!    available to the mesoscale circulation by the cumulus updraft (owm)
!    and the amount deposited in each level (wmps).
!----------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 2e20.12)')  &
                        'in meens: jk,wmp,owm= ',k,wmps(k),owm(k)
        endif

!----------------------------------------------------------------------
!    add the  source level value to the array accumulating the  profile
!    of total updraft source at each level (wmps). the if loop prevents
!    the inclusion of moisture which is available but above the top of 
!    the mesoscale updraft (the level of zero bupoyancy usually).  wmps
!    will only be non-zero at layers within the mesoscale updraft, 
!    but owm may be non-zero in layers above the updraft.
!--------------------------------------------------------------------
        if (wmps(k) /= 0.0) then
          wmps(k) = wmps(k) + owm(k)
          if (do_donner_tracer) then
            wtp(k,:) = wtp(k,:) + otm(k,:)
          endif
        endif
      end do   ! (end of k loop)

!--------------------------------------------------------------------
!    convert various moisture rates from kg(h2o) / kg(air) / sec to
!    g(h2o) / kg(air) / day.
!--------------------------------------------------------------------
      owm(:)  = owm(:)*8.64e07

!---------------------------------------------------------------------
!     calculate the portion of redistributed water vapor that condenses.
!     cycle until lowest level within the region of mesoscale circ-
!     ulation is reached. exit the loop when have marched past top of 
!     the mesoscale circulation.
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (phalf_c(k+1) > pzm) cycle
        if (phalf_c(k) < pztm) exit

!---------------------------------------------------------------------
!    determine if the current level is within the region of the meso-
!    scale circulation (between pzm and pztm).
!---------------------------------------------------------------------
        if ((phalf_c(k+1) <= pzm) .and. (phalf_c(k) >= pztm)) then

!---------------------------------------------------------------------
!    if so, define the top (pc2) of the current layer. deter-
!    mine the pressure level to which air in this layer will reach when
!    moving at the appropriate mesoscale updraft velocity for the dur-
!    ation of the mesoscale circulation (pctm). this level is limited to
!    be no higher than the top of the mesoscale circulation; if it is 
!    calculated to be higher, redefine the mesoscale updraft velocity 
!    for this layer so that the air in this layer will reach only to
!    the mesoscale circulation top, and no higher.
!---------------------------------------------------------------------
          pc2 = phalf_c(k+1)
          pctm = pc2 +Param%meso_ref_omega*Param%meso_lifetime
          if (pctm <= pztm)  then
            omer = (pztm - pc2)/Param%meso_lifetime
          else
            omer = Param%meso_ref_omega
          endif
          pctm = pc2 + omer*Param%meso_lifetime

!---------------------------------------------------------------------
!    define the temperature of the mesoscale updraft at this level.
!    determine its saturation vapor pressure and saturation mixing  
!    ratio. define saturation deficit
!    or excess relative to tempq(k), which is the mixing ratio in the 
!    mesoscale region (environmental mixing ratio plus source from 
!    cumulus updrafts). if there is a moisture excess (and thus conden-
!    sation must occur), define the condensation rate in the mesoscale
!    region, normalized over the mesoscale lifetime and its areal cover-
!    age. if only a portion of the layer is within the mesoscale updraft
!    region, adjust the mesoscale condensation rate appropriately.
!    if tempqa is greater than the saturation specific humidity (ERROR-
!    should be mixing ratio), reset it to the saturation value.
!---------------------------------------------------------------------
          ta = temp_c(k) + Param%tprime_meso_updrft
          call compute_mrs_k (ta, pfull_c(k), Param%d622 , Param%d608 ,&
                              mrsat, nbad)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (nbad /= 0) then
            ermesg = 'subroutine don_m_meso_updraft_k: '// &
                     'temperatures out of range of esat table'
            error = 1
            return
          endif

          q3 = mrsat - tempq(k)
          if (q3 <= 0.) then
            if (phalf_c(k+1) <= pctm)  then
              wmms(k) = (q3*ampta1/Param%meso_lifetime)*    &
                       (phalf_c(k) - pctm)/(phalf_c(k) - phalf_c(k+1))
            else
              wmms(k) = q3*ampta1/Param%meso_lifetime
            endif
          endif
          tempqa(k) = MIN (tempqa(k), mrsat)
        endif
      end do

!---------------------------------------------------------------------
!    determine the large-scale model full level at which parcel contain-
!    ing the water vapor at the base of the mesoscale updraft will reach
!    saturation and begin to condense (jsave).
!---------------------------------------------------------------------
      anv = 0.
      do k=1,nlev_lsm

!---------------------------------------------------------------------
!    determine the water vapor mixing ratio at the base of the mesoscale
!    updraft (qref).
!---------------------------------------------------------------------
        if (pfull_c(k) > pzm) cycle       
        if (anv == 0.) qref = tempqa(k)
        anv = 1.
        if (pfull_c(k) < pztm) exit        

!---------------------------------------------------------------------
!    define the temperature of the mesoscale updraft at this level.
!    determine its saturation vapor pressure and saturation specific
!    humidity. NOTE: should be mixing RATIO. define the level at which
!    mesoscale updraft condensation begins as the current level, in 
!    case the loop will be exited.
!---------------------------------------------------------------------
        te = temp_c(k) + Param%tprime_meso_updrft
        call compute_mrs_k (te, pfull_c(k), Param%d622 , Param%d608 , &
                            mrsat, nbad, mr = qref)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (nbad /= 0) then
          ermesg = 'subroutine don_m_meso_updraft_k: '// &
                   'temperatures out of range of esat table'
          error = 1
          return
        endif

        jsave = k

!---------------------------------------------------------------------
!    if in diagnostics column, output the values of saturation mixing  
!    ratio (mrsat) and mixing ratio in the mesoscale region (tempqa).
!---------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, 2e20.12)')  &
                          'in meens: qs,tempqa= ',mrsat,tempqa(k)
        endif

!---------------------------------------------------------------------
!    if there is a saturation excess at this level then exit, saving the
!    level index as jsave. this is the level at which condensation  in
!    the mesoscale updraft will begin.
!---------------------------------------------------------------------
        if (qref >= mrsat) exit      
      end do

!---------------------------------------------------------------------
!    define the  ???????
!!    What is the 6 ?? how is it related to the 8 below in the omd
!!    definition ???
!---------------------------------------------------------------------
      alp = 6.*Param%meso_ref_omega/((pzm - pztm)**2)

      omv = 0.

!---------------------------------------------------------------------
!    define the forcing terms associated with mesoscale updrafts.
!---------------------------------------------------------------------
      do k=1,nlev_lsm

!-------------------------------------------------------------------
!    if the current level is below the base of the mesoscale updraft,
!    cycle. if the current level is above the top of the mesoscale 
!    updraft, exit the loop.
!-------------------------------------------------------------------
        if (pfull_c(k) .gt. pzm) cycle       
        if (pfull_c(k) .lt. pztm) exit

!--------------------------------------------------------------------
!    define the limits of the current layer, modified from the large-
!    scale model levels when the mesoscale updraft region starts or ends
!    within the layer.
!--------------------------------------------------------------------
        pp = phalf_c(k+1)
        pm = phalf_c(k)
        if (phalf_c(k+1) < pztm) pp = pztm
        if (phalf_c(k) > pzm) pm = pzm

!---------------------------------------------------------------------
!    calculate mesoscale vertical velocity profile.
!---------------------------------------------------------------------
        omv(k) = (pzm + pztm)*((pp**2) - (pm**2))/2.
        omv(k) =  omv(k) - (((pp**3) - (pm**3))/3.)
        omv(k) = omv(k) - pztm*pzm*(pp - pm)
        omv(k) = omv(k)/(phalf_c(k+1) - phalf_c(k))
        omv(k) = omv(k)*alp

!---------------------------------------------------------------------
!    calculate mesoscale entropy-flux convergence. analytic integration
!    used, possible only because mesoscale temperature perturbation is 
!    not function of pressure. see "Vertical Velocity in Mesoscale 
!    Cloud" notes, 11/12/91.
!---------------------------------------------------------------------
        tmes_up(k) = (pzm + pztm)*(Param%rdgas - Param%cp_air)*  &
                     (pp - pm)/Param%cp_air
        tmes_up(k) = tmes_up(k) + ((2.*Param%cp_air - Param%rdgas)*  &
                     ((pp**2) - (pm**2))/(2.*Param%cp_air))
        tmes_up(k) = tmes_up(k) - (Param%rdgas*pztm*pzm/Param%cp_air)* &
                     alog(pp/pm)
        tmes_up(k) = tmes_up(k)/(phalf_c(k+1) - phalf_c(k))
        tmes_up(k) = tmes_up(k)*ampta1*Param%tprime_meso_updrft*alp

!--------------------------------------------------------------------
!    if currently below the level at which condensation in the meso-
!    scale updraft begins, cycle until that level is reached.
!--------------------------------------------------------------------
        if (k < jsave) cycle      

!--------------------------------------------------------------------
!    if into the region where deposition occurs, define the appropriate
!    above and below indices for boundary levels.
!--------------------------------------------------------------------
        if (k == 1) then
          jkm = k
        else
          jkm = k - 1
        endif
        if (k == nlev_lsm) then
          jkp = k
        else
          jkp = k + 1
        endif

!--------------------------------------------------------------------
!    define the temperature of the mesoscale updraft (te). define the
!    associated saturation vapor pressure and specific humidity (ERROR 
!    !!!).
!--------------------------------------------------------------------
        te = temp_c(k) + Param%tprime_meso_updrft
        call compute_mrs_k (te, pfull_c(k), Param%d622 , Param%d608 ,&
                              tempqa(k), nbad)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (nbad /= 0) then
          ermesg = 'subroutine don_m_meso_updraft_k: '// &
                   'temperatures out of range of esat table'
          error = 1
          return
        endif

!--------------------------------------------------------------------
!    if an excess of vapor is present and deposition should occur, 
!    define the mesoscale updraft temperature at the next higher level 
!    (tep). 
!--------------------------------------------------------------------
        if (qref >= tempqa(k)) then
          tep = temp_c(jkp) + Param%tprime_meso_updrft

!--------------------------------------------------------------------
!    if the next higher level is no longer in the mesoscale updraft 
!    layer, define the deposition rate in the mesoscale updraft at 
!    level k as the vapor flux divergence between layer k-1 and layer k.
!--------------------------------------------------------------------
          if (pfull_c(jkp) <= pztm) then
            cmu(k) = -omv(k)*(tempqa(k) - tempqa(jkm))/ &
                     (pfull_c(k) - pfull_c(jkm))

!--------------------------------------------------------------------
!     if level k is the lowest level within the condensation region,
!     determine the saturation specific humidity (ERROR !!!) at the
!     next higher level. define the deposition rate in the mesoscale  
!     updraft at level k as the vapor flux divergence between level k 
!     and level k+1. redefine qref as the amount of vapor remaining
!     in the parcel at the jkp level.
!--------------------------------------------------------------------
          else if (k == jsave) then
            call compute_mrs_k (tep, pfull_c(jkp), Param%d622 , &
                                Param%d608 , tempqa(jkp), nbad)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
            if (nbad /= 0) then
              ermesg = 'subroutine don_m_meso_updraft_k: '// &
                       'temperatures out of range of esat table'
              error = 1
              return
            endif

            cmu(k) = -omv(k)*(tempqa(jkp) - tempqa(k))/  &
                     (pfull_c(jkp) - pfull_c(k))
            qref = tempqa(jkp)

!--------------------------------------------------------------------
!     if level k is within the condensation region, determine the  
!     saturation specific humidity (ERROR !!!) at the next higher level.
!     define the deposition rate in the mesoscale updraft at level k as
!     the vapor flux divergence between level k-1 and level k+1. 
!     redefine qref as the amount of vapor remaining in the parcel at 
!     the jkp level.
!--------------------------------------------------------------------
          else
            call compute_mrs_k (tep, pfull_c(jkp), Param%d622 , &
                                Param%d608 , tempqa(jkp), nbad)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
            if (nbad /= 0) then
              ermesg = 'subroutine don_m_meso_updraft_k: '// &
                       'temperatures out of range of esat table'
              error = 1
              return
            endif

            cmu(k) = -omv(k)*(tempqa(jkp) - tempqa(jkm))/ &
                     (pfull_c(jkp) - pfull_c(jkm))
            qref = tempqa(jkp)
          endif

!---------------------------------------------------------------------
!    make certain that the deposition rate is non-negative.
!---------------------------------------------------------------------
          if (cmu(k) < 0.) cmu(k) = 0.

!---------------------------------------------------------------------
!    if there is insufficient moisture for deposition, set the depo-
!    sition rate to 0.0.
!---------------------------------------------------------------------
        else
          cmu(k) = 0.
        endif

!---------------------------------------------------------------------
!    convert the deposition rate to g(h2o) / kg(air) / day. multiply
!    by the anvil area (ampta1) to obtain a grid-box-mean value of the
!    deposition rate.
!---------------------------------------------------------------------
        cmu(k) = cmu(k)*ampta1*8.64e07

!--------------------------------------------------------------------
!    if in diagnostics column, output the environmental temperature
!    (temp_c) and the mesoscale vertical velocity (omv).
!--------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, i4, f20.14, e20.12)') &
                      'in meens: jk,t,omv= ', k, temp_c(k), omv(k)
        endif
      end do

!---------------------------------------------------------------------
!    calculate the mesoscale moisture-flux and tracer-flux convergence.
!---------------------------------------------------------------------
      do k=1,nlev_lsm 

!---------------------------------------------------------------------
!    if the current level is above the mesoscale updraft, exit the loop.
!    if the next level is still below the base of the mesoscale updraft,
!    cycle to the end of the loop.
!---------------------------------------------------------------------
        if (phalf_c(k) .lt. pztm) exit       
        if (phalf_c(k+1) .gt. pzm) cycle      

!--------------------------------------------------------------------
!    define the appropriate above and below indices for boundary levels.
!--------------------------------------------------------------------
        if (k == 1) then
          jkm = k
        else
          jkm = k - 1
        endif
        if (k == nlev_lsm) then
          jkp = k
        else
          jkp = k + 1
        endif

!---------------------------------------------------------------------
!    define the difference between the environmental vapor mixing ratio 
!    and that in the mesoscale updraft at the two half-levels bracketing
!    the current level.
!---------------------------------------------------------------------
        qprip = (tempqa(jkp) + tempqa(k) -    &
                             mixing_ratio_c(jkp) - mixing_ratio_c(k))/2.
        qprim = (tempqa(k) + tempqa(jkm) -    &
                             mixing_ratio_c(k) - mixing_ratio_c(jkm))/2.

!---------------------------------------------------------------------
!    define the difference between the environmental tracer mixing 
!    ratios and those in the mesoscale updraft at the two half-levels 
!    bracketing the current level.
!---------------------------------------------------------------------
        if (do_donner_tracer) then
          do kcont=1,ntr
            qtprip = (temptr(jkp,kcont) + temptr(k,kcont) - &
                      tracers_c(jkp,kcont) - tracers_c(k,kcont))/2.
            qtprim = (temptr(k,kcont) + temptr(jkm,kcont) -  &
                      tracers_c(k,kcont) - tracers_c(jkm,kcont))/2.
            eqtfp = ampta1*qtprip*alp*(phalf_c(k+1) - pztm)*  &
                    (pzm - phalf_c(k+1))
            eqtfm = ampta1*qtprim*alp*(phalf_c(k) - pztm)*  &
                    (pzm - phalf_c(k))
            if ((phalf_c(k) <= pzm) .and. (phalf_c(k+1) >= pztm)) then
              qtmes(k,kcont) = (eqtfm - eqtfp)/   &
                                              (phalf_c(k+1) - phalf_c(k))
            endif
            if ((pzm <= phalf_c(k)) .and. (pzm >= phalf_c(k+1))) then
              qtmes(k,kcont) = eqtfp/(phalf_c(k) - phalf_c(k+1))
            endif
            if ((pztm >= phalf_c(k+1)) .and. (pztm <= phalf_c(k))) then
              qtmes(k,kcont) = eqtfm/(phalf_c(k+1) - phalf_c(k))
              if ((pzm <= phalf_c(k)) .and. (pzm >= phalf_c(k+1))) then
                qtmes(k,kcont) = 0.
              endif
            endif ! ((pztm >= phalf_c(k+1)) .and. (pztm <= phalf_c(k)))
          end do
        endif

!-------------------------------------------------------------------
!    define the
!-------------------------------------------------------------------
        eqfp = ampta1*qprip*alp*(phalf_c(k+1) - pztm)*   &
                                                    (pzm - phalf_c(k+1))
        eqfm = ampta1*qprim*alp*(phalf_c(k) - pztm)*(pzm - phalf_c(k))
        if ((phalf_c(k) <= pzm) .and. (phalf_c(k+1) >= pztm)) then
          mrmes_up(k) = (eqfm - eqfp)/(phalf_c(k+1) - phalf_c(k))
        endif
        if ((pzm <= phalf_c(k)) .and. (pzm >= phalf_c(k+1))) then
          mrmes_up(k) = eqfp/(phalf_c(k) - phalf_c(k+1))
        endif
        if ((pztm >= phalf_c(k+1)) .and. (pztm <= phalf_c(k))) then
          mrmes_up(k) = eqfm/(phalf_c(k+1) - phalf_c(k))
          if ((pzm <= phalf_c(k)) .and. (pzm >= phalf_c(k+1))) then
            mrmes_up(k) = 0.
          endif
        endif ! ((pztm .ge. phalf_c(k+1)) .and. (pztm .le. phalf_c(k)))

!---------------------------------------------------------------------
!    if in diagnostics column,  output the entropy     (tmes) and
!    specific humidity (?)(mrmes) tendencies due to the mesoscale
!    updraft.
!---------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, i4, f19.10, f20.14, e20.12)')   &
                  'in meens: jk,pr,tmes,qmes= ', k, pfull_c(k),  &
                   tmes_up(k), mrmes_up(k)
        endif
      end do

!---------------------------------------------------------------------
!    calculate the eddy flux of moist static energy in mesoscale
!    updraft (hflux) and identify its minimum (hfmin).
!---------------------------------------------------------------------
      hfmin = 0.
      do jk=1,nlev_lsm
!---------------------------------------------------------------------
!    if the current level is above the mesoscale updraft, exit the loop.
!    if the next level is still below the base of the mesoscale updraft,
!    cycle to the end of the loop.
!---------------------------------------------------------------------
        if (pfull_c(jk) .lt. pztm) exit      
        if (pfull_c(jk) .gt. pzm) cycle      

!--------------------------------------------------------------------
!    define the temperature of the mesoscale updraft (tmu). define the
!    associated saturation vapor pressure and specific humidity (ERROR 
!    !!!).
!--------------------------------------------------------------------
        tmu = temp_c(jk) + Param%TPRIME_MESO_UPDRFT
        call compute_mrs_k (tmu, pfull_c(jk), Param%d622 , &
                                Param%d608 , qmu, nbad)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (nbad /= 0) then
          ermesg = 'subroutine don_m_meso_updraft_k: '// &
                   'temperatures out of range of esat table'
          error = 1
          return
        endif

!---------------------------------------------------------------------
!    define the eddy flux of moist static energy in the mesoscale 
!    updraft (hflux). retain the minimum value in the profile (hfmin)
!    and its pressure level (pfmin).
!---------------------------------------------------------------------
        hflux = omv(jk)*(((Param%cp_air*Param%tprime_meso_updrft ) + &
                                 Param%hlv*(qmu - mixing_ratio_c(jk))))
        if (hflux < hfmin) then
          hfmin = hflux      
          pfmin = pfull_c(jk)
        endif
      end do

!---------------------------------------------------------------------
!    if in a diagnostics column, output the minimum of the eddy moist 
!    static energy flux and its level.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)')  &
                        'in meens: hfmin,pfmin= ', hfmin, pfmin
      endif

!---------------------------------------------------------------------
!    define the mesoscale fractional area (cumh) in the region of the 
!    mesoscale updraft. 
!---------------------------------------------------------------------
      do k=1,nlev_hires
        if ((p_hires(k) <= pzm) .and. (p_hires(k) >= pztm))  then
          cumh(k) = ampta1
        else
          cumh(k) = 0.0 
        endif
      end do

!---------------------------------------------------------------------
!    call map_hi_res_col_to_lo_res_col to map the mesoscale anvil area 
!    from the cloud model to the large-scale 
!    model.
!---------------------------------------------------------------------
      call don_u_map_hires_c_to_lores_c_k  &
           (nlev_lsm, nlev_hires, cumh, p_hires, pztm + dp, phalf_c, &
            meso_cloud_area, rintsum, rintsum2, ermesg, error) 
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)')  &
            'in meens: rintsum(cuml) =', rintsum, rintsum2
        call don_u_compare_integrals_k   &
             (rintsum, rintsum2, diag_unit, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return
      endif

!---------------------------------------------------------------------
!    define the upward mass flux associated with the mesoscale 
!    circulation. 
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        umeml(k) = -omv(k)*ampta1/Param%grav  
        wmms(k)  = wmms(k)*8.64e07
        wmps(k)  = wmps(k)*8.64e07
      end do

!---------------------------------------------------------------------
!    obtain column integrals of deposition rate in the mesoscale (cmui),
!    convective updraft condensation (wmc), cell to mesoscale moisture
!    transfer (wpc), and the moisture made available to the mesoscale
!    by the cumulus updraft (owms). convert to units of mm / day.
!---------------------------------------------------------------------
      cmui = 0.
      wmc  = 0.
      wpc  = 0.
      owms = 0.
      do k=1,nlev_lsm
        wmc  = wmc  + wmms(k)*(phalf_c(k) - phalf_c(k+1))
        owms = owms + owm(k)*(phalf_c(k) - phalf_c(k+1))
        wpc  = wpc  + wmps(k)*(phalf_c(k) - phalf_c(k+1))
        cmui = cmui + cmu(k)*(phalf_c(k) - phalf_c(k+1))
      end do
      wmc  = wmc/(Param%grav*1000.)
      wpc  = wpc/(Param%grav*1000.)
      owms = owms/(Param%grav*1000.)
      cmui = cmui/(Param%grav*1000.)

!---------------------------------------------------------------------
!    if in diagnostics column, output the column-integral moisture 
!    conversion rates (wmc, wpc, owms, cmui). 
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12,a,a,e20.12,a)')  &
               'in meens: wmc=', wmc, ' mm/day', ' wpc=', wpc, 'mm/day'
        write (diag_unit, '(a, e20.12, a, a, e20.12, a)')  &
               'in meens: owms= ', owms, ' mm/day', ' cmui= ',   &
                        cmui, 'mm/day'
      endif

!---------------------------------------------------------------------
!    calculate precipitation resulting from the mesoscale circulation.
!    define the total additional condensate supplied to the column
!    by the mesoscale circulation, the sum of the deposition (wmc) and
!    additional condensation (cmui). 
!---------------------------------------------------------------------
      cmui = cmui - wmc

!--------------------------------------------------------------------


end subroutine don_m_meso_updraft_k



!#####################################################################

subroutine don_m_meso_downdraft_k    &
         (nlev_lsm, nlev_hires, diag_unit, debug_ijt, Param, Nml, &
          p_hires, pfull_c, temp_c, mixing_ratio_c, phalf_c, pb, &
          ampta1, dp, pztm, pzm, alp, hfmin, pmd, tmes_dn, mrmes_dn, &
          dmeml, ermesg, error)

!-------------------------------------------------------------------
!    subroutine meens computes the mesoscale effects of the composited
!    cloud ensemble on the heat, moisture and tracer budgets, producing
!    tendency terms which are to be applied to the large-scale model.
!    scheme employed here is a variation on procedure of Leary and 
!    Houze (JAS, 1980). for more details on notation, see 
!    "Cu Closure A notes," 2/97.
!-------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_nml_type
use sat_vapor_pres_k_mod, only: compute_mrs_k

implicit none

!-------------------------------------------------------------------
integer,                       intent(in)   :: nlev_lsm, nlev_hires, &
                                               diag_unit
logical,                       intent(in)   :: debug_ijt        
type(donner_param_type),       intent(in)   :: Param
type(donner_nml_type),         intent(in)   :: Nml  
real,   dimension(nlev_hires), intent(in)   :: p_hires
real,   dimension(nlev_lsm),   intent(in)   :: pfull_c, temp_c,  &
                                               mixing_ratio_c
real,   dimension(nlev_lsm+1), intent(in)   :: phalf_c
real,                          intent(in)   :: pb, ampta1, dp, pztm, &
                                               pzm, alp, hfmin 
real,                          intent(out)  :: pmd
real,   dimension(nlev_lsm),   intent(out)  :: tmes_dn, mrmes_dn, dmeml
character(len=*),              intent(out)  :: ermesg
integer,                       intent(out)  :: error

!---------------------------------------------------------------------
!   local variables:

      real, dimension(nlev_hires)     :: dmemh
      real, dimension(nlev_lsm)       :: tempt, tempqa
      real, dimension(nlev_lsm+1)     :: emt, emq

!     real :: qlo
      real    ::  es, mrsat, c2, c3, c1, fjk, fjkm, qb, fjkb, qbm, qmd, &
                  qsmd, fjkmd, qmmd, pi, psa, targ, tprimd, tb, qten, tten, &
                  omd, mrsb, wa, wb, tmd, rin, rintsum, rintsum2
      integer :: ncmd
      integer :: jksave, k, nbad

!----------------------------------------------------------------------

      ermesg = ' ' ; error = 0

      tmes_dn = 0.
      mrmes_dn = 0.
      tempt(:) = temp_c(:)
      emt(:) = 0.
      emq(:) = 0.
      tempqa(:) = mixing_ratio_c(:)

!---------------------------------------------------------------------
!    define the top of the mesoscale downdraft (pmd). it is assumed to 
!    be meso_sep Pa below the base of the mesoscale updraft. (no meso-
!    scale motion is assumed between the base of the mesoscale updraft 
!    and the top of the mesoscale downdraft.) make certain it is not 
!    below the surface.
!---------------------------------------------------------------------
      pmd = MIN(pzm + Param%meso_sep, phalf_c(1))
      ncmd = 1
      do k=1,nlev_hires         
        if (p_hires(k) < pmd ) then
          ncmd = k + 1
          exit
        endif
      end do

!---------------------------------------------------------------------
!    calculate mesoscale downdraft speed (omd) at top of mesoscale 
!    downdraft (pmd). follow Leary and Houze (1980,JAS) and set 
!    magnitude to half that in mesoscale updraft; this vertical pressure
!    velocity assumed constant with ht between pzm and cloud base (pb). 
!---------------------------------------------------------------------
      omd = -alp*((pzm-pztm)**2)/8.
      omd = omd/2.

!--------------------------------------------------------------------
!    calculate temperature and specific humidity in mesoscale
!    downdraft. 
!---------------------------------------------------------------------
      do k=1,nlev_lsm

!---------------------------------------------------------------------
!    if the current level is above the top of the mesoscale downdraft, 
!    exit the loop. if the level is below cloud base, cycle to the end
!    of the loop.
!---------------------------------------------------------------------
        if (pfull_c(k) < pmd) exit      
        if (pfull_c(k) > pb) cycle      

!---------------------------------------------------------------------
!    calculate c2, the relative humidity in the mesoscale downdraft,
!    after Table 3 of Leary and Houze (1980, JAS).
!---------------------------------------------------------------------
        c2 = 1. - (.3*(pfull_c(k) - pmd)/(pb - pmd))

!---------------------------------------------------------------------
!    calculate c3, the factor which yields the eddy flux of moist
!    static energy when multiplied by the minimum of moist static
!    energy in the mesoscale updraft. Multiply by 1.3 to take account
!    of convective downdrafts. See Fig. 7 of Leary and Houze
!    (1980,JAS).
!---------------------------------------------------------------------
        c3 = (pfull_c(k) - pmd)/(pb - pmd)
        c3 = 1.3*c3

!---------------------------------------------------------------------
!    see "Moist Static Energy A, 1/26/91" notes.
!---------------------------------------------------------------------
        targ = temp_c(k)
        call compute_mrs_k (targ, pfull_c(k), Param%d622 , &
                                Param%d608 , mrsat, nbad, esat=es)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (nbad /= 0) then
          ermesg = 'subroutine don_m_meso_downdraft_k: '// &
                   'temperatures out of range of esat table'
          error = 1
          return
        endif

        c1 = Param%d622*Param%hlv*es/   &
                                 (pfull_c(k)*Param%rvgas*(temp_c(k)**2))
        tprimd = c3*hfmin/omd
        tprimd = tprimd - Param%hlv*(c2*mrsat - mixing_ratio_c(k))
        tprimd = tprimd/(Param%cp_air + Param%hlv*c1*c2)
        tempt(k) = temp_c(k) + tprimd
        targ = tempt(k)
        call compute_mrs_k (targ, pfull_c(k), Param%d622 , &
                                Param%d608 , tempqa(k), nbad)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (nbad /= 0) then
          ermesg = 'subroutine don_m_meso_downdraft_k: '// &
                   'temperatures out of range of esat table'
          error = 1
          return
        endif

        tempqa(k) = c2*tempqa(k)

!---------------------------------------------------------------------
!    if in diagnostics column, output 
!---------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, 4e20.12)')  &
                    'in meens: tprimd,tempqa,q,qs= ',tprimd,   &
                    tempqa(k), mixing_ratio_c(k), mrsat
          write (diag_unit, '(a, f19.10, 2e20.12)')  &
                    'in meens: pr,rh,factr= ', pfull_c(k), c2, c3
        endif
      end do

!---------------------------------------------------------------------
!    calculate eddy fluxes of potential temperature and specific
!    humidity in mesoscale downdraft.
!---------------------------------------------------------------------
      do k=2,nlev_lsm-1

!---------------------------------------------------------------------
!    if the current level is above the top of the mesoscale downdraft, 
!    exit the loop. if the level is below cloud base, cycle to the end
!    of the loop.
!---------------------------------------------------------------------
        if (phalf_c(k) .lt. pmd) exit
        if (phalf_c(k) .gt. pb) cycle        

!---------------------------------------------------------------------
!    calculate potential temperature and specific humidity (?) fluxes
!    for pressure levels between cloud base and top of mesoscale down-
!    draft.
!---------------------------------------------------------------------
        if ((pfull_c(k-1) <= pb) .and. (pfull_c(k) >= pmd)) then
          fjk = ampta1*omd*((Param%ref_press/pfull_c(k))**     &
                   (Param%rdgas/Param%cp_air))*(tempt(k) - temp_c(k))    
          fjkm = ampta1*omd*((Param%ref_press/pfull_c(k-1))**  &
                   (Param%rdgas/Param%cp_air))*(tempt(k-1) - temp_c(k-1))
          emt(k) = (fjk + fjkm)/2.
          fjk = ampta1*omd*(tempqa(k) - mixing_ratio_c(k))
          fjkm = ampta1*omd*(tempqa(k-1) - mixing_ratio_c(k-1))
          emq(k) = (fjk + fjkm)/2.
        endif

!---------------------------------------------------------------------
!    calculate potential temperature and specific humidity (?) fluxes
!    for pressure levels below cloud base.
!---------------------------------------------------------------------
        if (pfull_c(k-1) >= pb) then
          fjk = ampta1*omd*((Param%ref_press/pfull_c(k))**   &
                 (Param%rdgas/Param%cp_air))*(tempt(k) - temp_c(k))
          call don_u_lo1d_to_hi0d_linear_k  &
               (nlev_lsm, mixing_ratio_c, pfull_c, pb, qb, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return

          if (debug_ijt   ) then
            write (diag_unit, '(a, i4, f19.10, f20.14)')  &
                           'in polat: k,p,x=', k, pb, qb
          endif
          call don_u_lo1d_to_hi0d_linear_k  &
               (nlev_lsm, temp_c, pfull_c, pb, tb, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return

          if (debug_ijt   ) then
            write (diag_unit, '(a, i4, f19.10, f20.14)')  &
                      'in polat: k,p,x=', k, pb, tb
          endif
        call compute_mrs_k (tb, pb, Param%d622 , &
                                Param%d608 , mrsb, nbad, esat=es)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (nbad /= 0) then
            ermesg = 'subroutine don_m_meso_downdraft_k: '// &
                     'temperatures out of range of esat table'
            error = 1
            return
          endif

          tprimd = hfmin/omd
          tprimd = tprimd - Param%hlv*(.7*mrsb - qb)
          c1 = Param%D622  *Param%hlv*es/(pb*Param%rvgas*(tb**2))
          tprimd = tprimd/(Param%cp_air + .7*Param%hlv*c1)
          fjkb = ampta1*omd*((Param%ref_press/pb)**      &
                                    (Param%rdgas/Param%cp_air))*tprimd
          wa = (phalf_c(k) - pfull_c(k))/(pb - pfull_c(k))
          wb = (pb - phalf_c(k))/(pb - pfull_c(k))
          emt(k) = wa*fjkb + wb*fjk
          fjk = ampta1*omd*(tempqa(k) - mixing_ratio_C(k))
          targ = tb + tprimd
          call compute_mrs_k (targ, pb, Param%d622 , &
                                Param%d608 , qbm, nbad)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (nbad /= 0) then
            ermesg = 'subroutine don_m_meso_downdraft_k: '// &
                     'temperatures out of range of esat table'
            error = 1
            return
          endif

          qbm = .7*qbm   
          fjkb = ampta1*omd*(qbm - qb)
          emq(k) = wa*fjkb + wb*fjk
        endif

!---------------------------------------------------------------------
!    calculate potential temperature and specific humidity (?) fluxes
!    for pressure levels at or above the top of the mesoscale downdraft.
!---------------------------------------------------------------------
        if (pfull_c(k) <= pmd) then
          fjkm = ampta1*omd*((Param%ref_press/pfull_c(k-1))**    &
                 (Param%rdgas/Param%cp_air))*(tempt(k-1) - temp_c(k-1))
          call don_u_lo1d_to_hi0d_linear_k  &
               (nlev_lsm, mixing_ratio_c, pfull_c, pmd, qmd, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return

          if (debug_ijt) then
            write (diag_unit, '(a, i4, f19.10, f20.14)')  &
                      'in polat: k,p,x=', k, pmd, qmd
          endif
          call don_u_lo1d_to_hi0d_linear_k  &
               (nlev_lsm, temp_c, pfull_c, pmd, tmd, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return

          if (debug_ijt   ) then
            write (diag_unit, '(a, i4, f19.10, f20.14)')  &
                      'in polat: k,p,x=', k, pmd, tmd
          endif
          call compute_mrs_k (tmd, pmd, Param%d622 , &
                                Param%d608 , qsmd, nbad, esat=es)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (nbad /= 0) then
            ermesg = 'subroutine don_m_meso_downdraft_k: '// &
                     'temperatures out of range of esat table'
            error = 1
            return
          endif

          c1 = Param%d622*Param%hlv*es/(pmd*Param%rvgas*(tmd**2))
          tprimd = -Param%hlv*(qsmd - qmd)/(Param%cp_air + Param%hlv*c1)
          fjkmd = ampta1*omd*((Param%ref_press/pmd)**   &
                                     (Param%rdgas/Param%cp_air))*tprimd
          wa = (pfull_c(k-1) - phalf_c(k))/(pfull_c(k-1) - pmd)
          wb = (phalf_c(k) - pmd)/(pfull_c(k-1) - pmd)
          emt(k) = fjkmd*wa + fjkm*wb
          targ = tmd + tprimd
          call compute_mrs_k (targ, pmd, Param%d622 , &
                                Param%d608 , qmmd, nbad)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (nbad /= 0) then
            ermesg = 'subroutine don_m_meso_downdraft_k: '// &
                     'temperatures out of range of esat table'
            error = 1
            return
          endif

          fjkm = ampta1*omd*(tempqa(k-1) - mixing_ratio_c(k-1))
          fjkmd = ampta1*omd*(qmmd - qmd)
          emq(k) = fjkmd*wa + fjkm*wb
        endif

!---------------------------------------------------------------------
!    if in diagnostics column, output the potential temprature and
!    specific humidity fluxes associated with the mesoscale downdrafts.
!---------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 3e20.12)')  &
                        'in meens: jk,phr,emt,emq= ', k ,phalf_c(k),   &
                         emt(k), emq(k)
        endif

!---------------------------------------------------------------------
!    convert the potential temperature flux to a temperature flux.
!---------------------------------------------------------------------
! RSH : unneeded, causes error
!       emt(k) = ((Param%ref_press/pfull_c(k))**     &
!                                    (Param%rdgas/Param%cp_air))*emt(k)
      end do  ! (end of k loop)

!---------------------------------------------------------------------
!    calculate temperature and specific humidity tendencies due
!    to eddy-flux convergences in mesoscale downdraft.
!---------------------------------------------------------------------
      rin = 0.
      do k=nlev_lsm,1, -1

!---------------------------------------------------------------------
!    define the index of the base of the mesoscale updraft (jksave).
!---------------------------------------------------------------------
        if ((phalf_c(k+1) <= pzm) .and. (phalf_c(k) >= pzm))   &
                                                         jksave = k + 1
        pi = (Param%ref_press/pfull_c(k))**(Param%rdgas/Param%cp_air)
        if ((emt(k+1) /= 0.) .and. (emt(k) == 0.) .and.    &
            (rin == 0.)) then
          tten = -emt(k+1)/(phalf_c(k+1) - phalf_c(1))
          qten = -emq(k+1)/(phalf_c(k+1) - phalf_c(1))
          rin = 1.
        endif
        if (rin == 1.) then
          
          if (.not. Nml%frc_internal_enthalpy_conserv) then
            tmes_dn(k) = tmes_dn(k) + (tten/pi)
          endif
          mrmes_dn(k) = mrmes_dn(k) + qten
        endif
        if ((rin == 0.) .and. (emt(k+1) /= 0.) .and.   &
            (emt(k) /= 0.)) then
          tten = (emt(k+1) - emt(k))/(phalf_c(k+1) - phalf_c(k))
          tten = -tten/pi
          qten = (emq(k+1) - emq(k))/(phalf_c(k+1) - phalf_c(k))
          qten = -qten
          tmes_dn(k) = tmes_dn(k) + tten
          mrmes_dn(k) = mrmes_dn(k) + qten
        endif

!---------------------------------------------------------------------
!    if in diagnostics column,  output the entropy     (tmes) and
!    specific humidity (?)(mrmes) tendencies due to the mesoscale
!    downdraft.
!---------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, i4, f19.10, f20.14, e20.12)')   &
                    'in meens: jk,pr,tmes,qmes= ', k, pfull_c(k),  &
                     tmes_dn(k), mrmes_dn(k)
        endif
      end do

!---------------------------------------------------------------------
!    define the temperature (tten)and moisture (qten) tendencies result-
!    ing from the mesoscale downdraft that are to be applied to the 
!    layers between the top of mesoscale downdraft (where emt is 
!    non-zero, saved as psa), and the base of the mesoscale updraft 
!    given by phalf_c(jksave).
!---------------------------------------------------------------------
      psa = 0.
      do k=1,nlev_lsm
        if ((emt(k) /= 0.) .and. (emt(k+1) == 0.)) then
          tten = emt(k)/(phalf_c(jksave) - phalf_c(k))
          qten = emq(k)/(phalf_c(jksave) - phalf_c(k))
          psa = phalf_c(k)
        endif
      end do

!---------------------------------------------------------------------
!    if in diagnostcs column, output the pressures at the top of the
!    mesoscale downdraft (pmd) and at cloud base (pb).
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 2f19.10)')  &
                                 'in meens: pmd,pb= ', pmd, pb
      endif

!--------------------------------------------------------------------
!    apply these tendencies to the levels between top of mesoscale
!    downdraft and base of mesoscale updraft.
!--------------------------------------------------------------------
      do k=1,nlev_lsm
        if ((pfull_c(k) <= psa) .and.    &
            (pfull_c(k) >= phalf_c(jksave))) then

!---------------------------------------------------------------------
!    if in diagnostcs column, output the pressure bounds of this region
!    (psa, phalf_c(jksave), the tendencies applied (qten, tten), and the 
!    large-scale model entropy     and moisture tendencies 
!    (mrmes, tmes) prior to the addition of these terms. 
!---------------------------------------------------------------------
          if (debug_ijt) then
            write (diag_unit, '(a, 3e20.12)')  &
                   'in meens: po,psa,phr(jksave)= ',  &
                           Param%REF_PRESS, psa, phalf_c(jksave)
            write (diag_unit, '(a, i4, 2e20.12)')  &
                       'in meens: jk,qmes,qten= ', k, mrmes_dn(k), qten
            write (diag_unit, '(a, i4, 2e20.12)')  &
                           'in meens: jk,tmes,tten= ', k, tmes_dn(k), tten
          endif

!---------------------------------------------------------------------
!    update the moisture and entropy tendencies.
!---------------------------------------------------------------------
          mrmes_dn(k) = mrmes_dn(k) + qten
!!! ISN't emt (and therefore tten) already temperature tendency rather 
!   than theta, and so the conversion here is unnecessary ??
          pi=(Param%ref_press/pfull_c(k))**(Param%rdgas/Param%cp_air)
          tmes_dn(k) = tmes_dn(k) + (tten/pi)
        endif
      end do

!---------------------------------------------------------------------
!    define the mass flux of the mesoscale down-
!    draft (dmemh) in the region of the mesoscale downdraft.
!---------------------------------------------------------------------
      do k=1,nlev_hires
        if ((p_hires(k) <= pb) .and. (p_hires(k) >= pmd))  then
          dmemh(k) = -omd*ampta1/Param%grav  
        else
          dmemh(k) = 0.
        endif
      end do

!---------------------------------------------------------------------
!    call map_hi_res_col_to_lo_res_col to map the 
!    mesoscale downdraft flux from the cloud model to the large-scale 
!    model.
!---------------------------------------------------------------------
      call don_u_map_hires_c_to_lores_c_k  &
           (nlev_lsm, nlev_hires, dmemh, p_hires, pmd + dp, phalf_c, &
            dmeml, rintsum, rintsum2, ermesg, error) 
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)')  &
            'in meens: rintsum(dmeml) =', rintsum , rintsum2
        call don_u_compare_integrals_k  &
             (rintsum, rintsum2, diag_unit, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return
      endif

!---------------------------------------------------------------------


end subroutine don_m_meso_downdraft_k



!####################################################################

subroutine don_m_meso_evap_k    &
         (nlev_lsm, diag_unit, debug_ijt, Param, available_condensate,&
          available_condensate_liq, available_condensate_ice, &
          pzm, pztm, phalf_c, emdi_liq, emdi_ice,      &
          emds_liq, emds_ice, emes_liq, emes_ice, ermesg, error)

!---------------------------------------------------------------------
!    subroutine meso_evap calculates the sublimation associated with
!    the mesoscale circulation and partitions both the updraft- and
!    downdraft-induced sublimation within the column.
!---------------------------------------------------------------------

use donner_types_mod, only : donner_param_type

implicit none

!---------------------------------------------------------------------
integer,                     intent(in)  :: nlev_lsm
integer ,                    intent(in)  :: diag_unit
logical ,                    intent(in)  :: debug_ijt 
type(donner_param_type),     intent(in)  :: Param
real,                        intent(in)  :: available_condensate,  &
                                            available_condensate_liq,  &
                                            available_condensate_ice,  &
                                            pzm, pztm
real, dimension(nlev_lsm+1), intent(in)  :: phalf_c
real,                        intent(out) ::       emdi_liq, emdi_ice
real, dimension(nlev_lsm),   intent(out) :: emds_liq, emds_ice
real, dimension(nlev_lsm),   intent(out) :: emes_liq, emes_ice
character(len=128),          intent(out) :: ermesg
integer,                     intent(out) :: error

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       available_condensate      total condensate available in the
!                                 anvil [ mm(h2o) / day ]
!       pzm                       pressure at base of mesoscale updraft
!                                 [ Pa ]
!       pztm                      pressure at top of mesoscale updraft
!                                 [ Pa ]
!       phalf_c                   pressures at large-scale model inter-
!                                 face levels [ Pa ]
!       diag_unit                 output unit number for this 
!                                 diagnostics column
!       debug_ijt                 is this a diagnostics column ?
!
!   intent(out) variables:
!
!       emdi                      vertical integral of mesoscale down-
!                                 draft sublimation [ mm (h2o) / day ]
!       emds                      water mass sublimated in mesoscale
!                                 downdraft [ g(h2o) / (kg(air) day) ] 
!       emes                      water mass sublimated in mesoscale
!                                 updraft [ g(h2o) / (kg(air) day) ] 
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!       emei                      vertical integral of mesoscale updraft
!                                 sublimation [ mm (h2o) / day ]
!       emea                      vertical integral of mesoscale updraft
!                                 sublimation [ g(h2o) / (kg(air) day) ]




     real    ::                   pm, pp, pbot_meso_sub
     real  :: emei_liq, emei_ice, emea_liq, emea_ice
     real  :: emda_liq, emda_ice
     integer :: k


      ermesg = ' ' ; error = 0
      emes_liq = 0.
      emes_ice = 0.

!---------------------------------------------------------------------
!    define the rate of total water sublimation in the mesoscale updraft
!    (emei) using the Leary and Houze coefficient.
!---------------------------------------------------------------------
      emei_liq = Param%meso_up_evap_fraction*available_condensate_liq
      emei_ice = Param%meso_up_evap_fraction*available_condensate_ice

!----------------------------------------------------------------------
!    convert the integral of mesoscale updraft evaporation from mm (h20)
!    per day to g(h2o) / (kg(air) / day (emea). updraft evaporation is 
!    assumed to occur between base of mesoscale updraft (pzm) and the 
!    top of mesoscale updraft (pztm) at a uniform rate. 
!---------------------------------------------------------------------- 
      emea_liq = emei_liq*Param%grav*1000./(pzm - pztm)
      emea_ice = emei_ice*Param%grav*1000./(pzm - pztm)

!--------------------------------------------------------------------
!    if in diagnostics column, output the column integral mesoscale
!    updraft evaporation, in units of both mm / day (emei) and 
!    g(h2o) / kg(air) / day (emea), and the pressures defining the
!    mesoscale updraft region (pzm, pztm).
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)')  &
                         'in meens: emea,emei= ',emea_liq + emea_ice, &
                                                  emei_liq + emei_ice
        write (diag_unit, '(a, 2e20.12)')  &
                         'in meens: LIQemea,emei= ',emea_liq, emei_liq
        write (diag_unit, '(a, 2e20.12)')  &
                         'in meens: ICEemea,emei= ',emea_ice, emei_ice
        write (diag_unit, '(a, 2f19.10)')  &
                         'in meens: pzm, pztm= ',pzm, pztm
      endif

!---------------------------------------------------------------------
!    call map_hi_res_intgl_to_lo_res_col to distribute the integrated 
!    mesoscale updraft sublimation (emea) within the mesoscale updraft
!    region (pzm -> pztm) of the large-scale model column whose layers
!    are defined by interface pressure array phalf_c.
!---------------------------------------------------------------------
      call don_u_map_hires_i_to_lores_c_k &
           (nlev_lsm, emea_liq, pzm, pztm, phalf_c, emes_liq, ermesg, error)
      call don_u_map_hires_i_to_lores_c_k &
           (nlev_lsm, emea_ice, pzm, pztm, phalf_c, emes_ice, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!----------------------------------------------------------------------
!    if in a diagnostics column, output the mesoscale updraft sublim-
!    ation profile (emes) at those levels where it is non-zero.
!----------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_lsm           
          if ((emes_liq(k) + emes_ice(k)) /= 0.0) then
            write (diag_unit, '(a, i4, e20.12)') &
                   'in cm_intgl_to_gcm_col: k,x= ',k,emes_liq(k) + &
                                                       emes_ice(k)
          endif
        end do
      endif

!---------------------------------------------------------------------
!    define the rate of total water sublimation in the mesoscale down-
!    draft (emdi) using the Leary and Houze coefficient.
!---------------------------------------------------------------------
      emdi_liq = Param%meso_down_evap_fraction*available_condensate_liq
      emdi_ice = Param%meso_down_evap_fraction*available_condensate_ice

!----------------------------------------------------------------------
!    convert the integral of mesoscale downdraft sublimation from mm 
!    (h20) per day to g(h2o) / (kg(air) / day (emda). downdraft sublim- 
!    ation is assumed to occur between base of mesoscale updraft (pzm)
!    and pressure pbot_meso_sub at a uniform rate.
!---------------------------------------------------------------------- 
      pbot_meso_sub = phalf_c(1)
      emda_liq = emdi_liq*Param%grav*1000./(pbot_meso_sub - pzm)
      emda_ice = emdi_ice*Param%grav*1000./(pbot_meso_sub - pzm)

!----------------------------------------------------------------------
!    distribute the integrated downdraft sublimation over the approp-
!    riate large-scale model layers.
!----------------------------------------------------------------------
      emds_liq = 0.
      emds_ice = 0.
      do k=1,nlev_lsm            

!---------------------------------------------------------------------
!    if the current level is above the base of the mesoscale updraft, 
!    exit the loop. if the level is below the surface, cycle to the end
!    of the loop.
!---------------------------------------------------------------------
        if (phalf_c(k) < pzm) exit
        if (phalf_c(k+1) > pbot_meso_sub ) cycle
        pm = phalf_c(k)
        pp = phalf_c(k+1)
        if ((phalf_c(k) >= pbot_meso_sub) .and.    &
            (phalf_c(k+1) <= pbot_meso_sub ) )  then
          pm = phalf_c(1) 
        endif
        if ((phalf_c(k) >= pzm) .and. (phalf_c(k+1) <= pzm))  pp = pzm
        emds_liq(k) = emda_liq*(pm - pp)*(pm + pp - 2.*pbot_meso_sub)
        emds_ice(k) = emda_ice*(pm - pp)*(pm + pp - 2.*pbot_meso_sub)

!--------------------------------------------------------------------
!    if in diagnostics column, output the column integral mesoscale
!    downdraft evaporation (emda) and the amount assigned to each layer
!    (emds), in units of g(h2o) / kg(air) / day.
!--------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 2e20.12)')  &
               'in meens: jk,emda,emd= ', k, emda_liq + emda_ice,  &
                                          emds_liq(k) + emds_ice(k)
          write (diag_unit, '(a, i4, 2e20.12)')  &
                  'in meens: jk,LIQemda,emd= ', k, emda_liq, emds_liq(k)
          write (diag_unit, '(a, i4, 2e20.12)')  &
                 'in meens: jk,ICEemda,emd= ', k, emda_ice, emds_ice(k)
        endif

!---------------------------------------------------------------------
!    
!---------------------------------------------------------------------
        emds_liq(k) = emds_liq(k)/((phalf_c(k) - phalf_c(k+1))*    &
                  (pzm - pbot_meso_sub))
        emds_ice(k) = emds_ice(k)/((phalf_c(k) - phalf_c(k+1))*    &
                  (pzm - pbot_meso_sub))
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 2e20.12)')  &
            'in meens: FINALjk,emda,emd= ', k, emda_liq + emda_ice,   &
                                             emds_liq(k) + emds_ice(k)
        endif
      end do

!--------------------------------------------------------------------
!    if in diagnostics column, output the column integral mesoscale
!    downdraft evaporation in units of mm / day (emdi) and the surface
!    pressure (ps).
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12, f19.10)')  &
              'in meens: emdi,ps= ', emdi_liq + emdi_ice, pbot_meso_sub
        write (diag_unit, '(a, e20.12, f19.10)')  &
                    'in meens: LIQemdi,ps= ', emdi_liq, pbot_meso_sub
        write (diag_unit, '(a, e20.12, f19.10)')  &
                    'in meens: ICEemdi,ps= ', emdi_ice, pbot_meso_sub
      endif

!---------------------------------------------------------------------


end subroutine don_m_meso_evap_k



!######################################################################

subroutine don_m_meso_melt_k   &
         (nlev_lsm, diag_unit, debug_ijt, Param, temp_c, phalf_c,  &
          pztm, meso_precip, pb, anvil_precip_melt, ermesg, error)

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

use donner_types_mod, only : donner_param_type

implicit none

!----------------------------------------------------------------------
integer,                     intent(in)  :: nlev_lsm
integer,                     intent(in)  :: diag_unit
logical,                     intent(in)  :: debug_ijt 
type(donner_param_type),     intent(in)  :: Param
real, dimension(nlev_lsm),   intent(in)  :: temp_c
real, dimension(nlev_lsm+1), intent(in)  :: phalf_c
real,                        intent(in)  :: pztm, meso_precip, pb
real, dimension(nlev_lsm),   intent(out) :: anvil_precip_melt
character(len=*),            intent(out) :: ermesg
integer,                     intent(out) :: error

  
      real  ::  p2, rma
      integer :: k

      anvil_precip_melt = 0.
      ermesg = ' ' ; error = 0

!--------------------------------------------------------------------
!    define the pressure at the melting level (p2).
!--------------------------------------------------------------------
      p2 = -10.
      do k=1,nlev_lsm-1
        if (phalf_c(k+1) < pztm ) exit
        if ((temp_c(k) >= Param%kelvin) .and.   &
             (temp_c(k+1) <= Param%kelvin))  then
          p2 = phalf_c(k+1)
          exit
        end if
      end do

!---------------------------------------------------------------------
!    define the rate of melting (mm/day)of anvil precipitation as it 
!    falls between the melting level (p2) and cloud base (pb). if there
!    is a melting level, all anvil precip is assumed to melt.
!---------------------------------------------------------------------
      if (p2 .ne. -10.) then
        rma = meso_precip
      else
        rma = 0.
      endif

!---------------------------------------------------------------------
!    convert the melting rate in mm / day (rm) to g(h2o) / kg(air) /day
!    (rma).
!---------------------------------------------------------------------
      rma = -rma*Param%grav*1000./(pb - p2)

!--------------------------------------------------------------------
!    if there is a melting level, map the melting rate uniformly across
!    the region of melting (melting level to cloud base). if there is
!    no melting level in the cloud, anvil_precip_melt is set to 0.0.
!--------------------------------------------------------------------
      if (pb > p2) then
        if (debug_ijt) then
          write (diag_unit, '(a, e20.12, 2f19.10)')  &
                     'in cm_intgl_to_gcm_col: xav,p1,p2= ',-rma, pb, p2
        endif
        call don_u_map_hires_i_to_lores_c_k  &
             (nlev_lsm, rma, pb, p2, phalf_c, anvil_precip_melt, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return

        if (debug_ijt) then
          do k=1,nlev_lsm               
            if (anvil_precip_melt(k) /= 0.0) then
              write (diag_unit, '(a, i4, e20.12)') &
                'in cm_intgl_to_gcm_col: k,x= ',k,-anvil_precip_melt(k)
            endif
          end do
        endif
      else                           
        anvil_precip_melt = 0.0
      endif

!---------------------------------------------------------------------

end subroutine don_m_meso_melt_k 



!#####################################################################

subroutine don_m_define_anvil_ice_k   &
         (isize, jsize, nlev_lsm, Param, Col_diag, pfull, temp,       &
          exit_flag, Don_conv, ermesg, error)

!----------------------------------------------------------------------
!    subroutine define_anvil_ice obtains the anvil ice profile
!    (Don_conv%xice) and the corresponding effective ice crystal size 
!    profile (Don_conv%dgeice).
!----------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_conv_type, &
                             donner_column_diag_type

implicit none

!-----------------------------------------------------------------------
integer,                               intent(in)    :: isize, jsize,   &
                                                        nlev_lsm
type(donner_param_type),               intent(in)    :: Param
type(donner_column_diag_type),         intent(in)    :: Col_diag
real, dimension(isize,jsize,nlev_lsm), intent(in)    :: pfull, temp
logical, dimension(isize,jsize),       intent(in)    :: exit_flag      
type(donner_conv_type),                intent(inout) :: Don_conv
character(len=*),                      intent(out)   :: ermesg
integer,                               intent(out)   :: error
            
!---------------------------------------------------------------------
!   intent(in) variables:
!
!     is, ie         first and last values of i index values of points 
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points 
!                    in this physics window (processor coordinates)
!     pfull          pressure at model full levels [ Pa ]
!     temp           temperature at model full levels [ deg K ]
!     total_precip   total convective-system precipitation [ mm / day ]
!
!   intent(inout) variables:
!
!     Don_conv
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:
  
     integer :: anvil_top_indx ! vertical index of highest model level
                               ! containing anvil ice 
     integer :: anvil_bot_indx ! vertical index of lowest model level
                               ! containing anvil ice 
     integer :: i, j, k, n     ! do-loop indices


!--------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    if column diagnostics are desired and there is mesoscale cloud 
!    present, output the mesoscale cloud area (ampta1), the mesoscale
!    downdraft sublimation integral (emdi), 
!---------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          if (Don_conv%ampta1(Col_diag%i_dc(n),   &
                              Col_diag%j_dc(n)) /= 0.0) then
            write (Col_diag%unit_dc(n), '(a, 2e20.12)')   &
                     'pre prean: ampta1, emdi, contot, tprei', &
                  Don_conv%ampta1 (Col_diag%i_dc(n),Col_diag%j_dc(n)), &
                  Don_conv%emdi_v    (Col_diag%i_dc(n),Col_diag%j_dc(n))
          endif
        end do
      endif

!--------------------------------------------------------------------
!    determine the vertical ice distribution and the effective ice 
!    size at each level.
!--------------------------------------------------------------------
      do j=1,jsize
        do i=1,isize
          if (.not. exit_flag(i,j)) then
!---------------------------------------------------------------------
!    call subroutine prean to assign the anvil ice content to the 
!    appropriate model layers. if there is no anvil in the column
!    (i.e., no mesoscale area, no mesoscale precip, no mesoscale down-
!    draft sublimation, no total precip), set the ice profile to be 0.0 
!    throughout the column.
!---------------------------------------------------------------------
            if (Don_conv%ampta1(i,j) > 0.0 .and.  &
                Don_conv%emdi_v(i,j) > 0.0 .and.  &
                Don_conv%meso_precip(i,j) > 0.0) then
              call don_m_prean_k   &
                   (i, j, nlev_lsm, Param, Col_diag, &
                    Don_conv%ampta1(i,j), Don_conv%meso_precip(i,j),&
                    Don_conv%emdi_v(i,j), pfull(i,j,:),  &
                    Don_conv%umeml(i,j,:), temp(i,j,:),   &
                    Don_conv%xice(i,j,:), ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
              if (error /= 0 ) return
            else
              Don_conv%xice(i,j,:) = 0.0      
            endif

!!! ANVIL ICE EXTENDS WELL BELOW FREEZING LEVEL  ?????
!!! BEGINS AT BASE OF MESOSCALE UPDRAFT WHICH IS NOT CONSTRAINED TO 
!!! BE AT OR BELOW FREEZING ??
           
!--------------------------------------------------------------------
!    determine the pressure at the top of the anvil (prztm). this will 
!    be the lowest model level pressure at which ice is present. at 
!    levels above this, set the effective ice size to 0.0.
!--------------------------------------------------------------------
            do k=1,nlev_lsm
              anvil_top_indx  = k 
              if ((Don_conv%xice(i,j,k) >= 1.0E-10)) then 
                Don_conv%prztm(i,j) = pfull(i,j,k)
                exit 
              endif
            end do

!--------------------------------------------------------------------
!    determine the pressure at the bottom of the anvil (przm). this 
!    will be the highest model level pressure at which ice is present.
!--------------------------------------------------------------------
            do k=anvil_top_indx+1,nlev_lsm
              if ((Don_conv%xice(i,j,k) < 1.0E-11) ) then
                Don_conv%przm(i,j) = pfull(i,j,k-1)
                anvil_bot_indx = k-1
                exit
              endif
            end do
          else
            Don_conv%xice(i,j,:) = 0.0      
            Don_conv%przm(i,j) = 0.0              
            Don_conv%prztm(i,j) = 0.0            
          endif
        end do
      end do

!---------------------------------------------------------------------
!    if column diagnostics are desired, output the level index, the 
!    pressure and the amount of ice at those levels at which ice is
!    present.
!---------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          do k=1,nlev_lsm
            if (Don_conv%xice(Col_diag%i_dc(n),  &
                              Col_diag%j_dc(n),k) > 0.0) then
              write (Col_diag%unit_dc(n), '(a, i4, e10.3, e20.12)')  &
                    'post prean: pressure, xice', &
                    k, pfull(Col_diag%i_dc(n), Col_diag%j_dc(n),k),    &
                    Don_conv%xice(Col_diag%i_dc(n), Col_diag%j_dc(n),k) 
            endif
          end do
        end do
      endif

!---------------------------------------------------------------------



end subroutine don_m_define_anvil_ice_k

!#####################################################################

subroutine don_m_prean_k     &
         (i, j, nlev_lsm, Param, Col_diag, ampta1_s, meso_precip_s,  &
          emdi_s, pfull_c, umeml_c, temp_c, xice_c, ermesg, error)

!---------------------------------------------------------------------
!    subroutine prean calculates the ice content assigned to the model
!    layers with an upward mesoscale mass flux.
!     Leo Donner
!     GFDL
!     17 May 2001
!---------------------------------------------------------------------
 
use donner_types_mod, only : donner_param_type, donner_column_diag_type

implicit none

!---------------------------------------------------------------------
integer,                       intent(in)  :: i, j, nlev_lsm
type(donner_param_type),       intent(in)  :: Param
type(donner_column_diag_type), intent(in)  :: Col_diag
real,                          intent(in)  :: ampta1_s, meso_precip_s, &
                                              emdi_s
real, dimension(nlev_lsm),     intent(in)  :: pfull_c, umeml_c, temp_c
real, dimension(nlev_lsm),     intent(out) :: xice_c 
character(len=*),              intent(out) :: ermesg
integer,                       intent(out) :: error

!-------------------------------------------------------------------
!   intent(in) variables:
!
!      i, j        horizontal coordinates of current column
!      ampta1_s    fractional area of mesoscale circulation
!                  [ fraction ]
!      emdi_s      vertical integral of  mesoscale-downdraft 
!                  sublimation [ mm / day ]
!      tprei       total convective-system precipitation [ mm / day ]
!      pfull_c     pressure at full model levels [ Pa ]
!      umeml_c     mesoscale updraft mass flux  [ kg / (m**2 sec) ]
!      temp_c      temperature at model full levels [ deg K ]
!
!   intent(out) variables:
!
!      xice_c      anvil ice [ kg(ice) / kg(air) ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:
 
      real      :: rho       !  height-averaged anvil air density 
                             !  [ kg/ (m**3) ]
      real      :: xicet     !  anvil ice work variable 
                             !  [ kg(ice) / kg (air) ]
      integer   :: kou       !  counter
      integer   :: k, kk, n  !  do-loop indices

!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    sum up the air density in the anvil layers.
!---------------------------------------------------------------------
      rho = 0.
      kou = 0
      do k=1,nlev_lsm             
        if (umeml_c(k) /= 0.) then
          kou = kou + 1
!  SHOULD virtual temp be used in defining this density ??
          rho = rho + pfull_c(k)/(Param%rdgas*temp_c(k))
        endif
      end do

!--------------------------------------------------------------------
!    if an anvil exists, determine the ice content of its layers.
!--------------------------------------------------------------------
      if (kou /= 0) then  

!--------------------------------------------------------------------
!    define the mean air density in the anvil region.
!--------------------------------------------------------------------
        rho = rho/kou

!----------------------------------------------------------------------
!    calculate the mesoscale ice content by balancing the fallout at 
!    anvil base with the mesoscale precipitation and the sublimation 
!    in the mesoscale downdraft.
!---------------------------------------------------------------------
!----------------------------------------------------------------------
!    sum up the anvil precipitation and the mesoscale downdraft sub-
!    limation (xicet).
!----------------------------------------------------------------------
        xicet = (meso_precip_s/86400.) + (emdi_s/86400.)
!!!???????????   DON'T KNOW THIS EXPRESSION ????
        xicet=xicet/(3.29*ampta1_s)
        xicet=xicet**.862
        xicet=xicet/rho

!----------------------------------------------------------------------
!    assign anvil ice to all layers with postive mesoscale updraft mass
!    flux.
!---------------------------------------------------------------------
        do k=1,nlev_lsm            
          if (umeml_c(k) > 0.) then
            xice_c(k) = xicet
          else
            xice_c(k) = 0.
          endif
        end do

!---------------------------------------------------------------------
!   if in diagnostics column, output the variables related to the
!   mesoscale ice content.
!---------------------------------------------------------------------
        if (Col_diag%in_diagnostics_window) then
          do n=1,Col_diag%ncols_in_window
            if (j == Col_diag%j_dc(n) .and. i == Col_diag%i_dc(n)) then
              do k=1, nlev_lsm            
                if (xice_c(k) > 0.00) then
                  write (Col_diag%unit_dc(n), '(a, 2e22.12)') &
                        'prean ampu,contot,emdi=', ampta1_s,  emdi_s
                  write (Col_diag%unit_dc(n), '(a, 1e22.12, i5)') &
                              'rho,     kou= ',rho,     kou
                  do kk=1,nlev_lsm           
                    if (xice_c(kk) > 0.00) then
                      write (Col_diag%unit_dc(n), '(a, i5, 2e22.12)') &
                                'k,prf,xice= ',kk,pfull_c(kk),xice_c(kk)
                      write (Col_diag%unit_dc(n), '(a, i5, 2e22.12)') &
                                 'k,prf,trf= ',kk,pfull_c(kk),temp_c(kk)
                      write (Col_diag%unit_dc(n), '(a, i5, 2e22.12)') &
                                   'k,prf,rmuf= ',kk,pfull_c(kk),umeml_c(kk)
                       write (Col_diag%unit_dc(n), '(a, i5, 2e22.12)') &
                     'k, grid box mean meso_precip_s, emdi_s = ', &
                        kk, meso_precip_s, emdi_s
                       write (Col_diag%unit_dc(n), '(a, i5, 2e22.12)') &
  'k,meso_precip_s, cloud area normalized xice in cloud (g / m**3) = ',&
                   kk, meso_precip_s,     &
        xice_c(kk)*ampta1_s*1.0e03*pfull_c(kk)/(temp_c(kk)*Param%rdgas)
                    endif
                  end do
                  exit
                endif
              end do
            endif
          end do
        endif

!---------------------------------------------------------------------
!    if there are no layers with positive mesoscale mass flux, set the
!    ice content at all levels to be 0.0.
!---------------------------------------------------------------------
      else
        xice_c(:) = 0.0
      endif

!--------------------------------------------------------------------


end subroutine don_m_prean_k



!######################################################################


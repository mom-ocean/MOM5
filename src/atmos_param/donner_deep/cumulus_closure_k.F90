
!VERSION NUMBER:
!  $Id: cumulus_closure_k.F90,v 19.0 2012/01/06 20:06:16 fms Exp $


!module cumulus_closure_inter_mod
!
!#include "cumulus_closure_interfaces.h"

!end module cumulus_closure_inter_mod

!######################################################################

subroutine cu_clo_cumulus_closure_k   &
         (nlev_hires, diag_unit, debug_ijt, Param, Initialized, &
          Nml, lofactor, dcape, cape_p, &
          qli0_v, qli1_v, qr_v, qt_v, env_r, ri_v, rl_v, parcel_r,   &
          env_t, parcel_t, a1, ermesg, error)

!---------------------------------------------------------------------
!    subroutine cumulus_closure calculates a_1(p_b) for closing the 
!    cumulus parameterization. see LJD notes, "Cu Closure D," 6/11/97
!---------------------------------------------------------------------
 
use donner_types_mod, only : donner_param_type, donner_nml_type, &
                             donner_initialized_type

implicit none

!---------------------------------------------------------------------
integer,                        intent(in)  :: nlev_hires
integer,                        intent(in)  :: diag_unit
logical,                        intent(in)  :: debug_ijt
type(donner_param_type),        intent(in)  :: Param
type(donner_initialized_type),  intent(in)  :: Initialized
type(donner_nml_type),          intent(in)  :: Nml    
real,                           intent(in)  :: lofactor, dcape
real,    dimension(nlev_hires), intent(in)  :: cape_p, qli0_v, qli1_v, &
                                               qr_v, qt_v, env_r, ri_v, &
                                               rl_v, parcel_r, env_t,   &
                                               parcel_t
real,                           intent(out) :: a1
character(len=*),               intent(out) :: ermesg
integer,                        intent(out) :: error

!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   intent(in) variables:
! 
!        cape_p        pressure on cape grid [ Pa ]
!        qli0_v        normalized component of cumulus condensate 
!                      forcing [ kg(h2o) / (kg(air) sec) ]
!                      defined in "Cu Closure D," p. 4.
!        qli1_v        un-normalized component of cumulus condensate
!                      forcing [ kg(h2o) / (kg(air) sec) ]
!                      defined in "Cu Closure D," p. 4.
!        qr_v          normalized cumulus moisture forcing 
!                      [ kg(h2o) / (kg(air) sec) ]
!                      defined in "Cu Closure D," p. 1.
!        qt_v          normalized cumulus thermal forcing 
!                      [ deg K / sec ]
!                      defined in "Cu Closure D," p. 1.
!        env_r         large-scale water-vapor mixing ratio 
!                      [ kg (h2o) / kg(air) ]
!        ri_v          large-scale ice mixing ratio 
!                      [ kg (h2o) / kg(air) ]
!        rl_v          large-scale liquid mixing ratio 
!                      [ kg (h2o) / kg(air) ]
!        parcel_r      parcel vapor mixing ratio  
!                      [ kg (h2o) / kg(air) ]
!        env_t         large-scale temperature [ deg K ]
!        parcel_t      parcel temperature [ deg K ]
!        dcape         rate of change of convective available potential
!                      energy due to large-scale processes 
!                      [ J / (kg s) ]
!        no_precip     logical array indicating columns in which there
!                      is no precip (and thus no deep convection)
!
!   intent(out) variables:
!
!        a1            fractional area of cumulus  ensemble
!        
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      real, dimension (nlev_hires)  :: rt, tden, tdena,  &
                                       dtpdta, pert_env_t, pert_env_r, &
                                       pert_parcel_t, pert_parcel_r, &
                                       parcel_r_clo, parcel_t_clo

      real     :: tau, cape_c
      logical  :: ctrig
      real     :: tdens, tdensa, ri1, ri2, rild, rile, rilf, ri2b,  &
                  sum2, rilak, rilbk, rilck, rilakm, rilbkm, rilckm, &
                  rila, rilb, rilc, ri2ak, ri2akm, ri2a, sum1, plcl, &
                  plfc, plzb, dumcoin, dumxcape
      integer  :: k     

!--------------------------------------------------------------------
!   local variables:
!
!      rt
!      ta
!      ra
!      tden
!      tdena
!      dtpdta
!      Cape_pert
!      tdens
!      tdensa
!      ri1
!      ri2
!      rild
!      rile
!      rilf
!      ri2d
!      sum2
!      rilak
!      rilbk
!      rilck
!      rilakm
!      rilbkm
!      rilckm
!      rila
!      rilb
!      rilc
!      ri2ak
!      ri2akm
!      ri2a
!      sum1
!      debug_ijt
!      perturbed
!      k       
!
!--------------------------------------------------------------------

      ermesg = ' ' ; error = 0

!--------------------------------------------------------------------
!    initialize the perturbed parcel profiles (pert_parcel_t,    
!    pert_parcel_r) and  the perturbed parcel environmental profiles to 
!    the actual parcel profiles.
!--------------------------------------------------------------------
      do k=1,nlev_hires
        pert_parcel_t(k) = env_t(k)
        pert_parcel_r(k) = env_r(k)
        pert_env_r(k)    = env_r(k)
        pert_env_t(k)    = env_t(k)
      end do

!--------------------------------------------------------------------
!    perturb lowest cape-model level mixing ratio and temperature so 
!    that one may calculate the derivative of parcel density temperature
!    w.r.t. surface large-scale density temperature. here the environ-
!    ment is made 1 deg K cooler and the mixing ratio is reduced to
!    99% of its unperturbed value.
!--------------------------------------------------------------------
      pert_env_r(1) = pert_env_r(1) - 0.01*pert_env_r(1)
      pert_env_r(1) = max(pert_env_r(1), 0.0)
      pert_env_t(1) = env_t(1) - 1.0

!---------------------------------------------------------------------
!    if this is a diagnostics column, output the environmental profiles
!    of temperature (pert_env_t) and vapor mixing ratio (pert_env_r) for 
!    the perturbed parcel, vertical profiles of pressure (cape_p), 
!    cumulus moisture forcing (qr_v), cumulus thermal forcing (qt_v), 
!    environmental moisture (env_r) and temperature (env_t) for the
!    unperturbed parcel, parcel temperature (parcel_t) and moisture 
!    (parcel_r) for the unperturbed parcel, cumulus condensate forcing 
!    (qli0 and qli1), ice condensate (ri_v) and liquid condensate (rl_v).
!---------------------------------------------------------------------
      if (debug_ijt) then 
        do k=1,nlev_hires
          write (diag_unit, '(a, i4, f19.10, f20.14, e20.12)')   &
                    'press, temp, vapor in cape: k, p,t,r = ',  &
                           k, cape_p(k), pert_env_t(k), pert_env_r(k)
        end do
        do k=1,nlev_hires
          if (qr_v(k) /= 0.0 .or. qt_v(k) /= 0.0) then
            write (diag_unit, '(a, i4, f19.10, 3e20.12, f20.14)') &
                  'in cuclo: k,p,qr,qt,r,t  =', k,  &
                     cape_p(k), qr_v(k), qt_v(k), env_r(k), env_t(k)
          endif
        end do
        do k=1,nlev_hires
          write (diag_unit, '(a, i4, f19.10, f20.14, e20.12)') &
                    'in cuclo: k,p,tpc, rpc   =', k,   &
                           cape_p(k), parcel_t(k), parcel_r(k)
        end do
        do k=1,nlev_hires
          if (qli0_v(k) /= 0.0 .or. qli1_v(k) /= 0.0 .or. &
              ri_v(k) /= 0.0 .or. rl_v(k) /= 0.0) then
              write (diag_unit, '(a, i4, f19.10, 4e20.12)')   &
                'in cuclo: k,p,qli0,qli1,ri,rl =', k,  &
                     cape_p(k), qli0_v(k), qli1_v(k), ri_v(k), rl_v(k)
          endif
        end do
      endif

      if (Nml%do_freezing_for_cape .NEQV. Nml%do_freezing_for_closure .or. &   ! kerr
          Nml%tfre_for_cape /= Nml%tfre_for_closure .or. &
          Nml%dfre_for_cape /= Nml%dfre_for_closure .or. &
          .not. (Initialized%use_constant_rmuz_for_closure) .or. &
          Nml%rmuz_for_cape /= Nml%rmuz_for_closure) then
           call don_c_displace_parcel_k   &
               (nlev_hires, diag_unit, debug_ijt, Param,  &
                Nml%do_freezing_for_closure, Nml%tfre_for_closure, &
                Nml%dfre_for_closure, Nml%rmuz_for_closure, &
                Initialized%use_constant_rmuz_for_closure,  &
                Nml%modify_closure_plume_condensate, &
                Nml%closure_plume_condensate, &
                env_t,  &
                env_r, cape_p, .false., plfc, plzb, plcl, dumcoin,  &
                dumxcape, parcel_r_clo,  parcel_t_clo, ermesg, error)
      else
        parcel_r_clo = parcel_r
        parcel_t_clo = parcel_t
      endif

!--------------------------------------------------------------------
!    call subroutine displace_parcel to determine the movement of a 
!    parcel from the lcl through the environment defined by (pert_env_t, 
!    pert_env_r).
!--------------------------------------------------------------------
      if (Nml%do_dcape) then

!--------------------------------------------------------------------
!    don't need to calculate cape when using Zhang closure 
!    (do_dcape is .true). 
!--------------------------------------------------------------------
      call don_c_displace_parcel_k   &
           (nlev_hires, diag_unit, debug_ijt, Param,   &
            Nml%do_freezing_for_closure, Nml%tfre_for_closure, &
            Nml%dfre_for_closure, Nml%rmuz_for_closure, &
            Initialized%use_constant_rmuz_for_closure, &
                Nml%modify_closure_plume_condensate, &
                Nml%closure_plume_condensate, &
            pert_env_t, &
            pert_env_r, cape_p, .false., plfc, plzb, plcl, dumcoin,  &
            dumxcape, pert_parcel_r,  pert_parcel_t, ermesg, error)
      else

!--------------------------------------------------------------------
!    if using cape relaxation closure then need to return cape value.
!--------------------------------------------------------------------
      call don_c_displace_parcel_k   &
           (nlev_hires, diag_unit, debug_ijt, Param,   &
            Nml%do_freezing_for_closure, Nml%tfre_for_closure, &
            Nml%dfre_for_closure, Nml%rmuz_for_closure, &
            Initialized%use_constant_rmuz_for_closure, &
                Nml%modify_closure_plume_condensate, &
                Nml%closure_plume_condensate, &
            pert_env_t, &
            pert_env_r, cape_p, .true., plfc, plzb, plcl, dumcoin,  &
            dumxcape, pert_parcel_r,  pert_parcel_t, ermesg, error)
     endif

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return


!---------------------------------------------------------------------
!    define quantities needed for cape relaxation closure option.
!---------------------------------------------------------------------
      if ( .not. Nml%do_dcape) then
        cape_c = Nml%cape0
        tau    = Nml%tau
        if (Nml%do_lands) then
          cape_c = Nml%cape0 * lofactor
          tau    = Nml%tau   * lofactor
        endif
      endif
      if (Nml%do_rh_trig) then
!  currently do_rh_trig forced to be .false.
!  no rh array available at current tiome in donner full.
!       rhavg=0.; dpsum=0.
!       do k = 1,sd%kmax
!         if (sd%p(k) .gt. Nml%plev0) then
!           rhavg  = rhavg + sd%rh(k)*sd%dp(k)
!           dpsum = dpsum + sd%dp(k)
!         end if
!       end do
!       rhavg = rhavg/dpsum
!       ctrig = rhavg > Nml%rhavg0
!! FOR NOW:
        ctrig= .true.
      else
        ctrig= .true.
      endif

!---------------------------------------------------------------------
!    if in a diagnostics column, output the path of the parcel (T, p
!    coordinates).
!---------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_hires
          write (diag_unit, '(a, i4, f20.14, e20.12)')  &
                    'in cuclo: k,tpca,rpca= ', k,    &
                                 pert_parcel_t(k), pert_parcel_r(k)
        end do
      endif

!---------------------------------------------------------------------
!    calculate the large-scale model profile of total-water mixing 
!    ratio. 
!---------------------------------------------------------------------
      do k=1,nlev_hires
        rt(k) = env_r(k) + ri_v(k) + rl_v(k)
      end do

!----------------------------------------------------------------------
!    calculate profiles of density temperatures, in the parcel (tden) 
!    and in the perturbed parcel (tdena). condensate is not included in
!    this definition of density temperature.
!----------------------------------------------------------------------
      do k=1,nlev_hires
        tden(k)  = parcel_t_clo(k)*(1. + (parcel_r_clo(k)/Param%d622)) 
        tdena(k) = pert_parcel_t(k)*(1. + (pert_parcel_r(k)/Param%d622))
      end do

!---------------------------------------------------------------------
!    define the values of density temperature in the environment at the
!    lowest level of the standard parcel displacement case (tdens) and 
!    for the displacement within the perturbed environment (tdensa).
!---------------------------------------------------------------------
      tdens  = env_t(1)*(1. + (env_r(1)/Param%d622))
      tdensa = pert_env_t(1)*(1. + (pert_env_r(1)/Param%d622))

!----------------------------------------------------------------------
!    evaluate derivative of parcel density temperature w.r.t. cloud-base
!    level environmental density temperature.
!----------------------------------------------------------------------
      do k=1,nlev_hires
        dtpdta(k) = (tdena(k) - tden(k))/(tdensa - tdens)
      end do

!---------------------------------------------------------------------
!    if this is a diagnostics column, output the profiles of unperturbed
!    parcel density temperature (tden) and the perturbed parcel density 
!    temperature (tdena) and the derivative of parcel density temper-
!    ature w.r.t. cloud-base large-scale density temperature (dtpdta).
!------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_hires
          write (diag_unit, '(a, i4, 2f20.14, e20.12)')  &
                   'in cuclo: k,tden(k),tdena(k),dtpdta(k)= ',   &
                          k,tden(k), tdena(k),dtpdta(k)
        end do
      endif

!--------------------------------------------------------------------
!    calculate the I1 and I2 integrals from p. 5 of "Cu Closure D" 
!    notes.
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!    define values at the cloud-base level.
!--------------------------------------------------------------------
      rild = qt_v(1)*(Param%d622 + env_r(1))/(Param%d622*(1. + rt(1)))
      rile = env_t(1)*(1. + rl_v(1) + ri_v(1) - Param%d622)*qr_v(1)
      rile = rile/(Param%d622*((1. + rt(1))**2))
      rilf = -env_t(1)*(Param%d622 + env_r(1))*qli0_v(1)
      rilf = rilf/(Param%d622*((1. + rt(1))**2))
      ri2b = env_t(1)*(Param%d622 + env_r(1))/   &
             (Param%d622*((1. + rt(1))**2))
      ri2b = ri2b*qli1_v(1)
      if (Nml%do_freezing_for_closure .or. &
          NMl%rmuz_for_closure /= 0.0) then
        sum2 = rild + rile + rilf
      else
        sum2 = 0.
      endif


      ri1 = 0.
      ri2 = 0.
      do k=2,nlev_hires
        if (cape_p(k) == 0.) exit       
        rilak = -qt_v(k)*(Param%d622 + env_r(k))/   &
                                     (Param%d622*(1. + rt(k)))
        rilbk = -env_t(k)*  &
                   (1. + rl_v(k) + ri_v(k) - Param%d622)*qr_v(k)
        rilbk = rilbk/(Param%d622*((1. + rt(k))**2))
        rilck = env_t(k)*(Param%d622 + env_r(k))*qli0_v(k)
        rilck = rilck/(Param%d622*((1. + rt(k))**2))
        rilakm = -qt_v(k-1)*(Param%d622 + env_r(k-1))/   &
                                          (Param%d622*(1. + rt(k-1)))
        rilbkm = -env_t(k-1)*  &
                     (1. + rl_v(k-1) + ri_v(k-1) - Param%d622)*qr_v(k-1)
        rilbkm = rilbkm/(Param%d622*((1. + rt(k-1))**2))
        rilckm = env_t(k-1)*(Param%d622 + env_r(k-1))*qli0_v(k-1)
        rilckm  =rilckm/(Param%d622*((1. + rt(k-1))**2))
        rila = .5*(rilak + rilakm)
        rilb = .5*(rilbk + rilbkm)
        rilc = .5*(rilck + rilckm)
        ri2ak = env_t(k)*(Param%d622 + env_r(k))/  &
                                         (Param%d622*((1. + rt(k))**2))
        ri2ak = ri2ak*qli1_v(k)
        ri2akm = env_t(k-1)*(Param%d622 + env_r(k-1))/ &
                                  (Param%d622*((1. + rt(k-1))**2))
        ri2akm = ri2akm*qli1_v(k-1)
        ri2a = .5*(ri2ak + ri2akm)
        sum1 = rila + rilb + rilc
        ri1 = ri1 + (alog(cape_p(k-1)/cape_p(k)))*   &
                                     (sum1 + dtpdta(k)*sum2)
        ri2 = ri2 + (alog(cape_p(k-1)/cape_p(k)))*  &
                                      (ri2a - dtpdta(k)*ri2b)

!----------------------------------------------------------------------
!    if in diagnostics column, output the 
!----------------------------------------------------------------------
        if (debug_ijt) then
          write(diag_unit, '(a, i4, e20.12)')   &
                        'in cuclo: k,dtpdta(k)= ',k,dtpdta(k)
          write (diag_unit,   '(a, 3e20.12)')  &
                           'in cuclo: rila,rilb,rilc= ', rila,rilb,rilc
          write (diag_unit, '(a, 2e20.12)')  &
                         'in cuclo: ri1,ri2= ',ri1,ri2
          write (diag_unit, '(a, 2e20.12)')  &
                       'in cuclo: sum1,sum2= ',sum1,sum2
        endif
      end do

!----------------------------------------------------------------------
!    if in diagnostics column, output the 
!----------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 3e20.12)')  &
                      'in cuclo: rild,rile,rilf= ', rild, rile, rilf
        if (dcape /= 0.0) then
          write (diag_unit, '(a, e20.12)')   &
                    'in cuclo:         dcape=',  dcape     
        endif
      endif

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      if (ri1 >= 0) then
        a1  = 0.
      else
        ri1 = Param%rdgas*ri1
        ri2 = Param%rdgas*ri2
        if (Nml%do_dcape .and. ctrig) then
!   Zhang closure:
          a1  = -(ri2 + dcape)/ri1
        else
          if (dumxcape > cape_c .and. ctrig) then
!   cape relaxation closure:
            a1 = -(ri2 + (dumxcape - cape_c)/tau)/ri1
          else
            a1 = 0.
          endif
        endif
      endif

!--------------------------------------------------------------------


end subroutine cu_clo_cumulus_closure_k



!######################################################################





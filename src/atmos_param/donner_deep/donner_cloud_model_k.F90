
!VERSION NUMBER:
!  $Id: donner_cloud_model_k.F90,v 19.0 2012/01/06 20:06:50 fms Exp $

!module donner_cloud_model_inter_mod

!#include "donner_cloud_model_interfaces.h"

!end module donner_cloud_model_inter_mod



!#####################################################################

subroutine don_cm_cloud_model_k   &
         (nlev_lsm, nlev_hires, ntr, kou, diag_unit, debug_ijt, Param, &
!++lwh
          Col_diag, Initialized, tb, pb, alpp, cld_press, temp_c, &
!--lwh
          mixing_ratio_c, pfull_c, phalf_c, tracers_c, pcsave, &
          exit_flag_c, wv, rcl, dpf, dpftr, qlw, dfr, flux, pt_kou,  &
          dint, cu, cell_precip, apt, cell_melt, pmelt_lsm, summel, &
          efchr, emfhr, cfracice, etfhr, ncc_kou, tcc, ermesg, error)

!--------------------------------------------------------------------
!
!                ONE-DIMENSIONAL CLOUD MODEL 
!                L. DONNER     NCAR     3 OCT 1984
!
!    subroutine cloud_model receives as input cloud base temperature 
!    (tb), pressure (pb), large-scale model vertical profiles of temper-
!    ature (temp_c), mixing ratio (mixing_ratio_c), pressure (pfull_c and
!    phalf_c), and tracer concentrations (tracers_c) for ensemble member
!    kou and produces as output the in-cloud profiles of temperature 
!    (tcc), vertical velocity (wv), cloud radius (rcl), liquid water 
!    (qlw), condensation rate (dpf), wet-deposition rate (dpftr), 
!    freezing rate (dfr), mass flux
!    (flux), tracer concentrations (xclo), the environmental profiles of 
!    temperature (te), mixing ratio (mre) and tracer concentrations 
!    (xtrae), cloud top pressure (pt_kou), and column integrals of pre-
!    cipitation rate (precip), condensation rate (conint) and freezing 
!    rate (dint). if the current column is a column for which diagnostics
!    is desired (debug_ijt = .true.), then output is written to the diag-
!    nostics file (diag_unit).
!--------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_column_diag_type, &
!++lwh
                             donner_initialized_type
!--lwh

implicit none

!--------------------------------------------------------------------
integer,                            intent(in)     :: nlev_lsm,  &
                                                      nlev_hires, ntr, &
                                                      kou, diag_unit
logical,                            intent(in)     :: debug_ijt
type(donner_param_type),            intent(in)     :: Param
type(donner_column_diag_type),      intent(in)     :: Col_diag
!++lwh
type(donner_initialized_type),      intent(in)     :: Initialized
!--lwh
real,                               intent(in)     :: tb, pb, alpp
real,    dimension(nlev_hires),     intent(in)     :: cld_press
real,    dimension(nlev_lsm),       intent(in)     :: temp_c,   &
                                                      mixing_ratio_c, &
                                                      pfull_c
real,    dimension(nlev_lsm+1),     intent(in)     :: phalf_c 
real,    dimension(nlev_lsm,ntr),   intent(in)     :: tracers_c
real,                               intent(inout)  :: pcsave
logical,                            intent(inout)  :: exit_flag_c     
real,    dimension(nlev_hires),     intent(out)    :: wv, rcl, dpf, &
                                                      qlw, dfr, flux, &
                                                      tcc 
real,                               intent(out)    :: pt_kou, dint, cu,&
                                                      cell_precip,  &
                                                             apt
real,    dimension(nlev_lsm),       intent(out)    :: cell_melt
real,    dimension(nlev_hires),     intent(out)    :: efchr, emfhr,   &
                                                      cfracice
real,    dimension(nlev_hires,ntr), intent(out)    :: etfhr, dpftr
integer,                            intent(out)    :: ncc_kou
real,                               intent( in)    :: pmelt_lsm
real,                               intent(out)    :: summel
character(len=*),                   intent(out)    :: ermesg
integer,                            intent(out)    :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     tb             cloud base temperature [ deg K ]
!     pb             cloud base pressure [ Pa ]             
!     temp_c         large-scale model temperature profile (index 1 
!                    nearest the surface) [ deg K ]
!     mixing_ratio_c large-scale model mixing ratio profile (index 1 
!                    nearest the surface) [ kg(h2o)/ kg(air) ]
!     sig            sigma coordinate of large-scale model levels
!                    index 1 nearest the surface) [ dimensionless ]
!     tracers_c      large-scale model tracer concentration profiles 
!                    (index 1 nearest the surface) [ kg/ kg ]  ???
!     kou            current ensemble member index i
!     diag_unit      output unit number for this diagnostics column
!     debug_ijt      is this a diagnostics column ?
!
!   intent(out) variables:
!
!     tcc            in-cloud temperature profile (index 1 at physical 
!                    base of cloud)  [ deg K ]
!     wv             in-cloud vertical velocity profile (index 1 at 
!                    physical base of cloud)  [ m / sec ]
!     rcl            cloud radius profile (index 1 at physical 
!                    base of cloud)  [ m ]
!     dpf            cloud-area-weighted condensation rate profile
!                    (index 1 at physical base of cloud) 
!                    [ (m**2) * kg(h2o) / kg(air) / sec ]
!     dpftr          wet-deposition rate profile,
!                    weighted by ratio of fractional area to
!                    fractional area at base 
!                    (index 1 at physical base of cloud)
!                    [  [units of xclo] /(sec) ]
!     qlw            cloud liquid water content profile (index 1 at
!                    physical base of cloud) [ kg(h2o) / kg(air) ]
!     dfr            cloud-area-weighted freezing rate profile in 
!                    convective updraft (index 1 at physical base of 
!                    cloud) [ (m**2) *g(h2o) / (kg(air) /  day ]
!     flux           upward mass flux profile in cloud (index 1 at 
!                    physical base of cloud) [ kg(air) / sec ]
!     xclo           in-cloud tracer concentration profiles (index 1 at 
!                    physical base of cloud)  [ kg / kg ] ???
!     te             environmental temperature profile on cloud-model 
!                    grid (index 1 at physical base of cloud) [ deg K ]
!     mre            environmental mixing ratio profile on cloud-model 
!                    grid (index 1 at physical base of cloud) 
!                    [ kg(air) / kg(h2o) ]
!     cell_melt      in-cloud melting of condensate associated with 
!                    convective cells. made up of two parts, 1) that 
!                    due to the freezing of liquid carried upwards 
!                    in the cell updraft, 2) that due to the melting of
!                    condensed ice that precipitates out. if meso 
!                    circulation is present, this component is zero; 
!                    melting will be determined in subroutine mesub.
!     xtrae          environmental tracer profiles on cloud-model 
!                    grid (index 1 at physical base of cloud) 
!                    [ kg / kg ]  ?????
!     pt_kou         pressure at cloud top for this ensemble member[ Pa ]
!     precip         total precipitation rate [ kg(h2o) / sec ]
!     conint         total condensation rate [ kg(h2o) / sec ]
!     dint           total freezing rate in convective updraft 
!                    [ kg(h20) / sec ]
!     cu
!     rc
!     lfc_not_reached      level of free convection was reached ?
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      real, dimension (nlev_hires)         ::  pf, te, mre
      real, dimension (nlev_hires,ntr)     ::  xclo, xtrae, pftr 
      real        :: precip, conint,                 pmel
      integer     :: k, kc
      real        :: accond, acpre
      real        :: sumfrea, sumlhr
      logical     :: lfc_not_reached, do_donner_tracer
      integer     ::  cldtop_indx

!--------------------------------------------------------------------
!   local variables:
!
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!--------------------------------------------------------------------
!    if in diagnostic column, output the ensemble member index which
!    is being integrated.
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, i4)')    &
            'PROCESSING CLOUD ENSEMBLE MEMBER # ', kou
      endif

!-------------------------------------------------------------------
!    define  or initialize an array to hold any in-cloud tracer sources.
!-------------------------------------------------------------------
      call don_cm_gen_incloud_profs_k  &
           (nlev_lsm, nlev_hires, ntr, kou, diag_unit, debug_ijt, &
            Col_diag, Param, Initialized, tb, pb, alpp, cld_press,  &
            temp_c, mixing_ratio_c, pfull_c, phalf_c, tracers_c,  &
            pcsave,tcc, wv, rcl, qlw, dfr, flux, pf, pftr, te, mre, &
            xclo, xtrae, dint, accond, acpre, sumfrea, cldtop_indx, &
            do_donner_tracer, lfc_not_reached, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!--------------------------------------------------------------------
!    define the cloud top pressure for the current ensemble member.
!--------------------------------------------------------------------
      pt_kou = pb + cldtop_indx*Param%dp_of_cloud_model

!--------------------------------------------------------------------
!    if in diagnostics column, output the ensemble member number and
!    the level of free convection of the most-entraining (first)
!    ensemble member.
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, i4, f19.10)') &
                        'in cloudm: kou,pcsave= ',kou,pcsave
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      call don_cm_process_condensate_k     &
         (nlev_lsm, nlev_hires, ntr, cldtop_indx, diag_unit, debug_ijt,&
          Param, acpre, accond, pb, pt_kou, pf, pftr, tcc, rcl,  &
          cld_press, phalf_c, conint, dint, pmel, pmelt_lsm, precip, &
          cu, cell_precip, sumlhr, summel, dpf, dpftr, dfr, cell_melt, &
          ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!----------------------------------------------------------------------
!    if the level of free convection was never reached for this ensemble
!    member, set exit_flag_c so that calculations will cease in this 
!    column, and exit the ensemble members loop. if this is a diag-
!    nostics column, print a message.
!----------------------------------------------------------------------
      if (lfc_not_reached) then
        if (debug_ijt) then
          write (diag_unit, '(a)')  &
                           'in mulsub: lfc never reached'     
        endif
        exit_flag_c = .true.
        return
      endif

!----------------------------------------------------------------------
!    if no precipitation was produced by the ensemble member, set 
!    exit_flag_c so that calculations will cease in this column, and exit
!    the ensemble members loop. if this is a diagnostics column, print 
!    a message.
!----------------------------------------------------------------------
      if (precip == 0.) then
        if (debug_ijt) then
          write (diag_unit, '(a)') &
                       'in mulsub: PRECIP=0 AFTER CLOUD MODEL'
        endif
        exit_flag_c   = .true.
        return
      endif

!---------------------------------------------------------------------
!    if condensate is being evaporated at any level within the cloud, 
!    set exit_flag_c so that calculations will cease in this 
!    column, and exit the ensemble members loop. if this is a diag-
!    nostics column, print a message.
!---------------------------------------------------------------------
      do kc=1,nlev_hires-1           
        if (dpf(kc) > 0.) then 
          if (debug_ijt) then
            write (diag_unit, '(a)') 'in mulsub: dpf  .GT. 0.'
          endif
          exit_flag_c   = .true.
          return   
        endif 
      end do
      if (exit_flag_c) return

!--------------------------------------------------------------------
!    define the cloud model vertical index that is just above cloud top
!    (ncc_kou). this index may be used to limit calculations on the cloud
!    model dimensioned arrays.
!--------------------------------------------------------------------
      do k=1,nlev_hires
        if (cld_press(k) < pt_kou) then
          ncc_kou = k 
          exit
        endif
      end do

!--------------------------------------------------------------------
!    call compute_vertical_fluxes to calculate cumulus thermal forcing 
!    and moisture forcing on the cloud-model grid associated with this
!    ensemble member. 
!--------------------------------------------------------------------
      call don_cm_compute_vert_fluxes_k   &
           (nlev_hires, ntr, ncc_kou, kou, diag_unit, debug_ijt, &
            do_donner_tracer, Param, pt_kou, cld_press, rcl, te, mre, &
            wv, tcc, dpf, xclo, xtrae, apt, efchr, emfhr, etfhr, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!----------------------------------------------------------------------
!    define the fraction of condensation that is ice (cfraci) and that
!    which is liquid (cfracl). at temps above tfre, all condensate is
!    liquid;at temps below (tfre-dfre) all condensate is ice; in between
!    the condensate is partitioned linearly with the temperature depart-
!    ure between these two limits.
!----------------------------------------------------------------------
      cfracice = 0.
      do k=1,ncc_kou
        if (qlw(k) > 0.0) then
          if (tcc(k) > Param%tfre) then
            cfracice(k) = 0.
          else if (tcc(k) < (Param%tfre - Param%dfre)) then
            cfracice(k) = 1.
          else
            cfracice(k) = (Param%tfre - tcc(k))/Param%dfre
          endif
        endif
      end do

!---------------------------------------------------------------------


end subroutine don_cm_cloud_model_k



!#####################################################################

subroutine don_cm_gen_incloud_profs_k  &
         (nlev_lsm, nlev_hires, ntr, kou, diag_unit, debug_ijt, &
          Col_diag,  Param, Initialized, tb, pb, alpp, cld_press, &
          temp_c, mixing_ratio_c, pfull_c, phalf_c, tracers_c, pcsave, &
          tcc, wv, rcl, qlwa, dfr, flux, pf, pftr, te, mre, xclo, &
          xtrae, dint, accond, acpre, sumfrea, cldtop_indx,  &
          do_donner_tracer, lfc_not_reached, ermesg, error)

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

!  modified by Leo Donner, GFDL, 5 February 2007
!
use donner_types_mod, only : donner_param_type, &
                             donner_column_diag_type, &
                             donner_initialized_type
use sat_vapor_pres_k_mod, only: compute_mrs_k

implicit none 

!----------------------------------------------------------------------
integer,                            intent(in)    :: nlev_lsm,   &
                                                     nlev_hires,&
                                                     ntr, kou, diag_unit
logical,                            intent(in)    :: debug_ijt
type(donner_column_diag_type),      intent(in)    :: Col_diag
type(donner_param_type),            intent(in)    :: Param      
!++lwh
type(donner_initialized_type),      intent(in)    :: Initialized
!--lwh
real,                               intent(in)    :: tb, pb, alpp
real,    dimension(nlev_hires),     intent(in)    :: cld_press
real,    dimension(nlev_lsm),       intent(in)    :: temp_c,   &
                                                     mixing_ratio_c,  &
                                                     pfull_c
real,    dimension(nlev_lsm+1),     intent(in)    :: phalf_c    
real,    dimension(nlev_lsm,ntr),   intent(in)    :: tracers_c
real,                               intent(inout) :: pcsave         
real,    dimension(nlev_hires),     intent(out)   :: tcc, wv, rcl, qlwa,&
                                                     dfr, flux, pf
real,    dimension(nlev_hires),     intent(out)   :: te, mre
real,    dimension(nlev_hires,ntr), intent(out)   :: xclo, xtrae, pftr
real,                               intent(out)   :: dint, accond, acpre
real,                               intent(out)  :: sumfrea
integer,                            intent(out)   :: cldtop_indx
logical,                            intent(out)   :: do_donner_tracer, &
                                                     lfc_not_reached
character(len=*),                   intent(out)   :: ermesg
integer,                            intent(out)   :: error

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!   intent(in) variables:
!
!     psfc           surface pressure [ Pa ]
!     tb             cloud base temperature [ deg K ]
!     pb             cloud base pressure [ Pa ]             
!     temp_c         large-scale model temperature profile (index 1 
!                    nearest the surface) [ deg K ]
!     mixing_ratio_c large-scale model mixing ratio profile (index 1 
!                    nearest the surface) [ kg(h2o)/ kg(air) ]
!     sig            sigma coordinate of large-scale model levels
!                    index 1 nearest the surface) [ dimensionless ]
!     tracers_c      large-scale model tracer concentration profiles 
!                    (index 1 nearest the surface) [ kg/ kg ]  ???
!     kou            current ensemble member index i
!     diag_unit      output unit number for this diagnostics column
!     debug_ijt      is this a diagnostics column ?
!
!   intent(out) variables:
!
!     tcc            in-cloud temperature profile (index 1 at physical 
!                    base of cloud)  [ deg K ]
!     wv             in-cloud vertical velocity profile (index 1 at 
!                    physical base of cloud)  [ m / sec ]
!     rcl            cloud radius profile (index 1 at physical 
!                    base of cloud)  [ m ]
!     pf             cloud-area-weighted condensation rate profile
!                    (index 1 at physical base of cloud) 
!                    [ (m**2) * kg(h2o) / kg(air) / sec ]
!     pftr           cloud-area-weighted wet-deposition profile
!                    (index 1 at physical base of cloud)
!                    [ (m**2) * [xlco units] /sec ]
!     qlwa           cloud liquid water content profile (index 1 at
!                    physical base of cloud) [ kg(h2o) / kg(air) ]
!     dfr            cloud-area-weighted freezing rate profile in 
!                    convective updraft (index 1 at physical base of 
!                    cloud) [ (m**2) *g(h2o) / (kg(air) /  day ]
!     flux           upward mass flux profile in cloud (index 1 at 
!                    physical base of cloud) [ kg(air) / sec ]
!     xclo           in-cloud tracer concentration profiles (index 1 at 
!                    physical base of cloud)  [ kg / kg ] ???
!     te             environmental temperature profile on cloud-model 
!                    grid (index 1 at physical base of cloud) [ deg K ]
!     mre            environmental mixing ratio profile on cloud-model 
!                    grid (index 1 at physical base of cloud) 
!                    [ kg(air) / kg(h2o) ]
!     cell_melt      in-cloud melting of condensate associated with 
!                    convective cells. made up of two parts, 1) that 
!                    due to the freezing of liquid carried upwards 
!                    in the cell updraft, 2) that due to the melting of
!                    condensed ice that precipitates out. if meso 
!                    circulation is present, this component is zero; 
!                    melting will be determined in subroutine mesub.
!     xtrae          environmental tracer profiles on cloud-model 
!                    grid (index 1 at physical base of cloud) 
!                    [ kg / kg ]  ?????
!     pt             pressure at cloud top [ Pa ]
!     precip         total precipitation rate [ kg(h2o) / sec ]
!     conint         total condensation rate [ kg(h2o) / sec ]
!     dint           total freezing rate in convective updraft 
!                    [ kg(h20) / sec ]
!     rc
!     lfc_not_reached      level of free convection was reached ?
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      real,   dimension(nlev_lsm)         :: sig
      real,   dimension(nlev_hires)       :: dis, rsc
      real,   dimension (nlev_hires,ntr)  ::   clsou
      real     :: qlw, qcw, qrw, dtfr, dtupa, dfrac, rmu,  &
                  rhodt_inv, psfc, actot, dt_micro,          dcw1,  &
                  dqrw3, rbar, rmub, density, densityp, d1, d2, dztr, &
                  entrain, dt_inv
      logical  :: flag
      integer  :: max_cloud_level, ktr
      integer  :: k, nbad
!++lwh
      real, parameter :: g_2_kg = 1.e-3 ! kg/g
      real :: qlw_save, t_avg
      real, dimension( size(xclo,2) ) :: delta_xclo0, delta_xclo1, dwet
      integer :: n
!--lwh

!--------------------------------------------------------------------
!   local variables:
!
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = '  ' ; error = 0

      psfc = phalf_c(1)
      if (ntr /= 0) then
        do_donner_tracer = .true.
      else
        do_donner_tracer = .false.
      endif

!--------------------------------------------------------------------
!    provide appropriate initialization for the output fields. 
!    NOTE: input tracer fields of all zeroes are currently provided if 
!    do_donner_tracer is .false.. one could make the tracer fields
!    optional and define output only when present.
!--------------------------------------------------------------------
      tcc(:)  = 0.
      wv(:)   = 0.
      rcl(:)  = 0.
      qlwa(:) = 0.
      dfr(:)  = 0.
      te(:)   = 0.
      mre(:)  = 0.
      xclo(:,:)  = 0.
      xtrae(:,:) = 0.
      pftr(:,:) = 0.
      flux(:) = 0.
      dint   = 0. 
      lfc_not_reached = .true.

!-------------------------------------------------------------------
!    initialize local variables.
!-------------------------------------------------------------------
      do k=1,nlev_hires
        dis(k)=0.
        pf(k)=0.
      end do

      rsc(:) = 0.
      sig(:) = pfull_c(:)/psfc

!--------------------------------------------------------------------
!    if in diagnostics column, output the cloud base radius, pressure
!    and temperature, and the large-scale model sigma levels.
!--------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_lsm-Col_diag%kstart+1
          write (diag_unit, '(a, i4, e20.12)')  &
                        'in cloudm: k,sig   = ',k, sig(k)    
        end do
      endif

!---------------------------------------------------------------------
!    initialize various integrals.
!---------------------------------------------------------------------

      qcw = 0.
      qlw = 0.
      qrw = 0.

!---------------------------------------------------------------------
!    initialize various integrals.
!---------------------------------------------------------------------
      accond  = 0.
      acpre   = 0.
      dtupa   = 0.
      dfrac   = 0.
      dtfr    = 0.
      sumfrea = 0.

!---------------------------------------------------------------------
!    initialize various conditionals.
!---------------------------------------------------------------------
      cldtop_indx = 0
  
!--------------------------------------------------------------------
!    define the cloud model pressure levels. 
!--------------------------------------------------------------------
      max_cloud_level = nlev_hires + 2
      do k=1,nlev_hires
        if (cld_press(k) < Param%pstop) then
          max_cloud_level = k - 1 
        endif
      end do

!--------------------------------------------------------------------
!    define cloud base (k = 1) values of vertical velocity (wv),
!    temperature (tcc), mixing ratio (rsc), cloud radius (rcl) and
!    tracer concentrations (xclo).
!--------------------------------------------------------------------
      wv(1) = Param%cld_base_vert_vel
      rcl(1) = Param%cloud_base_radius
      tcc(1) = tb
      xclo(1,:) = tracers_c(1,:)
      call compute_mrs_k (tb, pb, Param%D622, Param%D608, rsc(1), nbad)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (nbad /= 0) then
        ermesg = 'subroutine don_cm_gen_incloud_profs_k: '// &
                 'temperatures out of range of esat table'
        error = 1
        return
      endif

!--------------------------------------------------------------------
!    call gcm_to_cm to obtain the large-scale model moisture, temper-
!    ature, and tracer field values at the cloud base pressure to be
!    used as the environmental values in the cloud model.
!--------------------------------------------------------------------
      call don_u_lo1d_to_hi0d_log_k             &
           (nlev_lsm, mixing_ratio_c, sig, psfc, pb, mre(1), ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

      call don_u_lo1d_to_hi0d_log_k             &
           (nlev_lsm, temp_c, sig, psfc, pb, te(1), ermesg, error )
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return


      if (do_donner_tracer) then
        do ktr=1,ntr          
          call don_u_lo1d_to_hi0d_log_k             &
               (nlev_lsm, tracers_c(:,ktr), sig, psfc, pb,   &
                xtrae(1,ktr), ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return
        end do
      endif

!--------------------------------------------------------------------
!    if in diagnostics column, output the cloud base radius, pressure
!    and temperature, and the large-scale model sigma levels.
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12, f19.10, f20.14)')  &
                'in cloudm: RR,PB,TB= ',Param%cloud_base_radius, pb, tb
      endif

!--------------------------------------------------------------------
!    if in diagnostics column, output cloud base values of environ-
!    mental moisture (mre) and temperature (te), and the in-cloud 
!    temperature (tcc).
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12,f20.14)')  &
                       'in cloudm: QE,TE,TCC= ',mre(1), te(1), tcc(1)
      endif
      
!---------------------------------------------------------------------
!    loop over cloud model levels, generating the desired in-cloud and 
!    environmental profiles. under no circumstances extend the cal-
!    culation above max_cloud_level (the level above which cloud is not
!    allowed).
!---------------------------------------------------------------------
      do k=1,max_cloud_level-1

!----------------------------------------------------------------------
!    exit the loop if the pressure is lower than the upper limit allowed
!    for the level of free convection and free convection has not yet 
!    been achieved.
!----------------------------------------------------------------------
        if ((cld_press(k+1) <= Param%upper_limit_for_lfc) .and. &
            (lfc_not_reached)) then
          cldtop_indx = 0
          exit
        endif 

!--------------------------------------------------------------------
!    call gcm_to_cm to obtain the large-scale model moisture, temper-
!    ature, and tracer field values at cloud model pressure level to be
!    used as the environmental values in the cloud model.
!--------------------------------------------------------------------
        call don_u_lo1d_to_hi0d_log_k             &
             (nlev_lsm, mixing_ratio_c, sig, psfc, cld_press(k+1),  &
              mre(k+1), ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return

        call don_u_lo1d_to_hi0d_log_k             &
             (nlev_lsm, temp_c, sig, psfc, cld_press(k+1),   &
              te(k+1), ermesg, error )

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return

        if (do_donner_tracer) then
          do ktr=1, ntr           
            call don_u_lo1d_to_hi0d_log_k             &
                 (nlev_lsm, tracers_c(:,ktr), sig, psfc,   &
                  cld_press(k+1), xtrae(k+1,ktr), ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
            if (error /= 0 ) return
          end do
        endif

!--------------------------------------------------------------------
!    if in diagnostics column, output level k values of vertical
!    velocity, pressure, environmental temperature (te) and moisture 
!    (mre), the in-cloud temperature (tcc) and liquid water (qlw).
!--------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, 3e20.12)')  &
                      'in cloudm: WV,PP,TE= ',wv(k), cld_press(k), te(k)
          write (diag_unit, '(a,f20.14, 2e20.12)')  &
                     'in cloudm: TCC(k),qe(k),QLW= ',tcc(k), wv(k), qlw
        endif

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!++lwh
        qlw_save = qlw
!--lwh

        call don_cm_move_parcel_k    &
             (k, kou, diag_unit, debug_ijt, cld_press(k),  &
              cld_press(k+1), alpp, Param, pcsave, qlw, sumfrea, qrw, &
              qcw, dcw1, dqrw3, wv(k), tcc(k), rcl(k), te(k), mre(k), &
              rcl(k+1), wv(k+1), tcc(k+1), rsc(k+1), te(k+1), mre(k+1),&
              qlwa(k+1), dfr(k+1), rmu, rbar, dfrac, dtfr, dtupa,  &
              lfc_not_reached, flag, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return

        if (flag) exit

!---------------------------------------------------------------------
!    define the in-cloud  density at levels k and k+1 
!    (density, densityp).
!---------------------------------------------------------------------
        density  = cld_press(k)/(Param%rdgas*(tcc(k)*  &
                               (1. + Param%D608*(rsc(k)/(1.0+rsc(k))))))
        densityp = cld_press(k+1)/(Param%rdgas*(tcc(k+1)*  &
                           (1. + Param%D608*(rsc(k+1)/(1.0+rsc(k+1))))))
 
!--------------------------------------------------------------------
!    calculate in-cloud tracer distribution.
!--------------------------------------------------------------------
        if (do_donner_tracer) then
          clsou(k,:) = 0.0
          clsou(k+1,:) = 0.0
          d1 = Param%rdgas*(1.+Param%virt_mass_co)*tcc(k)*te(k)/  &
               (Param%grav*cld_press(k)*(Param%virt_mass_co*te(k) +  &
                                                                tcc(k)))
          d2 = Param%rdgas*(1.+Param%virt_mass_co)*tcc(k+1)*te(k+1)/ &
               (Param%grav*cld_press(k+1)*      &
                                (Param%virt_mass_co*te(k+1) + tcc(k+1)))
          dztr = ((d1 + d2)*(cld_press(k) - cld_press(k+1))/2.)
          rmub = 2.*alpp/rbar
          entrain = dztr*rmub
          dt_micro = dztr/(0.5*(wv(k) + wv(k+1)))
          if (debug_ijt) then
            write (diag_unit, '(a, 3e20.12)')   &
                 'in clotr: d1,d2,dz= ',d1,d2,dztr
            write (diag_unit, '(a, e20.12)')   &
                 'in clotr: ent= ',rmub
          endif

!++lwh
!--------------------------------------------------------------------
!    Call tracer wet deposition
!    
!    Convert dqrw3 from g/m3 to kg/m3
!    from g(h2o) per m**3 to kg(h2o) per kg(air).
!--------------------------------------------------------------------
          t_avg = 0.5*(tcc(k)+tcc(k+1))
          do n = 1,size(xclo,2)
             if (Initialized%wetdep(n)%Lwetdep) then
                call wet_deposition_0D    &
                           (Initialized%wetdep(n)%Henry_constant, &
                            Initialized%wetdep(n)%Henry_variable, &
                            Initialized%wetdep(n)%frac_in_cloud, &
                            Initialized%wetdep(n)%alpha_r, &
                            Initialized%wetdep(n)%alpha_s, &
                            t_avg, cld_press(k), cld_press(k+1), &
                            0.5*(density+densityp), &
                            qlw_save, dqrw3*g_2_kg, 0., &
                            xclo(k,n), &
                            Initialized%wetdep(n)%Lgas, &
                            Initialized%wetdep(n)%Laerosol, &
                            Initialized%wetdep(n)%Lice, &
                            delta_xclo0(n) )
             end if
          end do
!--lwh          
          call don_cm_clotr_k    &
               (ntr, diag_unit, debug_ijt, Param, clsou(k,:),  &
                clsou(k+1,:), xtrae(k,:), xtrae(k+1,:), xclo(k,:), &
                entrain, dt_micro, xclo(k+1,:), ermesg, error)
!++lwh
          do n = 1,size(xclo,2)
             if (Initialized%wetdep(n)%Lwetdep) then
                call wet_deposition_0D   &
                           (Initialized%wetdep(n)%Henry_constant, &
                            Initialized%wetdep(n)%Henry_variable, &
                            Initialized%wetdep(n)%frac_in_cloud, &
                            Initialized%wetdep(n)%alpha_r, &
                            Initialized%wetdep(n)%alpha_s, &
                            t_avg, cld_press(k), cld_press(k+1), &
                            0.5*(density+densityp), &
                            qlw, dqrw3*g_2_kg, 0., &
                            xclo(k+1,n), &
                            Initialized%wetdep(n)%Lgas, &
                            Initialized%wetdep(n)%Laerosol, &
                            Initialized%wetdep(n)%Lice, &
                            delta_xclo1(n) )
                dwet(n) = - 0.5*(delta_xclo0(n)+delta_xclo1(n))
                xclo(k+1,n) = xclo(k+1,n) + dwet(n)
             else
                dwet(n) = 0.
             end if
          end do
!--lwh          
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return
        endif

!--------------------------------------------------------------------
!    define cloud-area-weighted precip removal from the layer (dis) in 
!    units of  m**2 * kg(h2o) per kg(air) per day. define rhodt_inv 
!    (10e-03*g*w/delta p) to be used to convert units of dqrw3 and dcw1
!    from g(h2o) per m**3 to kg(h2o) per kg(air).
!--------------------------------------------------------------------
        rhodt_inv = 1.0e-03*(0.5*(wv(k) + wv(k+1)))*Param%GRAV  /   &
                    Param%dp_of_cloud_model
        dis(k) = dqrw3*(rbar**2)*rhodt_inv
        pf(k)  = dcw1*(rbar**2)*rhodt_inv

!---------------------------------------------------------------------
!     define the removal by wet deposition between levels k and k+1
!---------------------------------------------------------------------
        dt_inv = -.25*(wv(k)+wv(k+1))*(density+densityp)*Param%grav  / &
                 Param%dp_of_cloud_model
        do n = 1,size(xclo,2)
          pftr(k,n) = dwet(n) *(rbar**2)*dt_inv
        end do

!---------------------------------------------------------------------
!    calculate moisture subject to freezing in units of g(h2o) per 
!    kg(air) per day, weighted by the cloud area (dfr). add this layer's
!    contribution to the integrated water mass frozen in the updraft 
!    in units of kg(h2o)/sec (dint). 
!---------------------------------------------------------------------
        if (dfr(k+1) /= 0.) then
          dfr(k+1) = dfr(k+1)*densityp*wv(k+1)*(rcl(k+1)**2)*   &
                                      Param%grav/Param%dp_of_cloud_model 
          dfr(k+1) = -dfr(k+1)*Param%cp_air/Param%hlf   
          dint = dint - dfr(k+1)*Param%dp_of_cloud_model/Param%grav  
          dfr(k+1) = dfr(k+1)*8.64E07
        endif 
 
!---------------------------------------------------------------------
!    add this layer's contribution to the integrals accumulating total
!    condensation (accond), precipitation (acpre) and  total liquid
!    water (actot).
!---------------------------------------------------------------------
        accond = accond + dcw1
        acpre  = acpre  + dqrw3
        actot =  accond + acpre

!--------------------------------------------------------------------
!    if in diagnostics column, output the values of accumulated conden-
!    sation, precipitation and total liquid water at this level.
!--------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 2e20.12)')  &
                        'in cloudm: k,ACCOND,ACPRE= ',k, accond, acpre
          write (diag_unit, '(a, i4, e20.12)')  &
                        'in cloudm:  k,ACTOT= ',k, actot
        endif

!--------------------------------------------------------------------
!    define the mass flux at level k+1 (flux).
!--------------------------------------------------------------------
        if (k == 1) flux(k) = (rcl(k)**2)*wv(k)*density 
        flux(k+1) = (rcl(k+1)**2)*wv(k+1)*densityp

!--------------------------------------------------------------------
!    if in diagnostics column, output the values of cloud radius and
!    entrainment coefficient at this level.
!--------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 2e20.12)')  &
                     'in cloudm: k,rcl,RMU= ',k, rcl(k), rmu
        endif

!---------------------------------------------------------------------
!    if non-negative water is present, save the current k index as the
!    index of the last cloud model level containing cloud (cldtop_indx).
!---------------------------------------------------------------------
        if (qlwa(k+1) >= 0.) then    
          cldtop_indx = k
        endif  ! (qlwa(k+1) >= 0.0)

!---------------------------------------------------------------------
!    end of loop over cloud model levels. upon leaving loop, in-cloud
!    profiles of model variables have been generated.
!---------------------------------------------------------------------
      end do   ! (do 1 loop)


!--------------------------------------------------------------------
!    if in diagnostics column, output the water mass frozen in the
!    updraft in units of cloud area*kg(h2o) / sec (dint) and the 
!    resultant temperature change (sumfrea).
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)')  &
                'in cloudm: DINT IN CLOUDM,sumfrea= ', dint, sumfrea
      endif

!--------------------------------------------------------------------
!    normalize the column integral of water mass frozen in the updraft
!    by the cloud base area of ensemble member #1.
!--------------------------------------------------------------------
      dint = dint/(Param%cloud_base_radius**2)

!--------------------------------------------------------------------
!    if an acceptable cloud was found for this ensemble member, obtain
!    the needed diagnostic / integral output.
!--------------------------------------------------------------------
      if (cldtop_indx /= 0) then

!--------------------------------------------------------------------
!    if in diagnostics column, output the topmost cloud level and
!    the mixing ratios at and one level above cloud top.  
!--------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 2e20.12)')  &
                      'in cloudm: n,rsc(n+2),rsc(n)= ',cldtop_indx,   &
                             rsc(cldtop_indx+1),rsc(cldtop_indx)
        endif
      endif ! (cldtop_indx /= 0)

!--------------------------------------------------------------------
!    if in diagnostics column, output the values of cloud-area-
!    weighted layer-mean condensation and  fallout.
!--------------------------------------------------------------------
        do k = 1, cldtop_indx
          if (debug_ijt) then
            write (diag_unit, '(a, i4, 2e20.12)')   &
                       'in cloudm: K,PF,DIS= ', k, pf(k), dis(k)
          endif
        end do


end subroutine don_cm_gen_incloud_profs_k


!#####################################################################

subroutine don_cm_move_parcel_k    &
              (k, kou, diag_unit, debug_ijt, pbot, ptop, alpp, Param,  &
               pcsave, qlw, sumfrea, qrw, qcw, dcw1, dqrw3, wvbot,  &
               tccbot, rclbot, tebot, mrebot, rcltop, wvtop, tcctop, &
               rsctop, tetop, mretop, qlwatop, dfrtop, rmu, rbar, &
               dfrac, dtfr,  dtupa, lfc_not_reached, flag, ermesg, error)

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

use donner_types_mod, only : donner_param_type
use sat_vapor_pres_k_mod, only: compute_mrs_k

implicit none

!----------------------------------------------------------------------
integer,                 intent(in)    :: k, kou, diag_unit
logical,                 intent(in)    :: debug_ijt
real,                    intent(in)    :: pbot, ptop, alpp
type(donner_param_type), intent(in)    :: Param
real,                    intent(inout) :: pcsave, qlw, sumfrea, qrw, &
                                          qcw, dcw1, dqrw3, wvbot,  &
                                          tccbot, rclbot, tebot, &
                                          mrebot, rcltop, wvtop, &
                                          tcctop, rsctop,&
                                          tetop, mretop, qlwatop,  &
                                          dfrtop, rmu, rbar, dfrac,  &
                                          dtfr, dtupa
logical,                 intent(inout) :: lfc_not_reached
logical,                 intent(out)   :: flag
character(len=*),        intent(out)   :: ermesg
integer,                 intent(out)   :: error

!----------------------------------------------------------------------
      real    ::  dtdp, drdp, dwdp, tcest, west, &
                  rest, qcwest,  qlwest, dtdp2, drdp2, dwdp2, qrwest
      integer :: nbad

!---------------------------------------------------------------------
      flag = .false.

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    define the entrainment coefficient (rmu). 
!    pcsave    = p at which ensemble member 1 becomes buoyant.
!    lfc_not_reached = T below level where vert vel starts to increase 
!    with height 
!    THUS:
!    rmu = 0 when above buoyancy level for member 1, above level where 
!            vert vel starts to increase for this member, AND vert vel
!            is so small that are detraining. Thus rmu = 0. when
!            detraining after having been through a convective tower.
!---------------------------------------------------------------------
      if ((pbot    <= pcsave) .and. (.not. lfc_not_reached) .and. &
          (wvbot <= Param%wdet)) then
        rmu = 0.
      else
        rmu = 2.*alpp/rclbot
      end if 

!--------------------------------------------------------------------
!    call simult to solve for the in-cloud derivatives dT/dp (dtdp), 
!    dw/dp (dwdp) and dr/dp (drdp) for a parcel in motion from cloud 
!    model level k (pressure pp) to pressure level pp + delta p.
!--------------------------------------------------------------------
      call don_cm_simult_k   &
           (diag_unit, debug_ijt, lfc_not_reached, Param, pcsave, rmu, &
            tccbot, rclbot, wvbot, pbot, qlw, tebot, mrebot,   &
            dtdp, drdp, dwdp, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!--------------------------------------------------------------------
!    if in diagnostics column, output cloud base values of environ-
!    mental moisture (mre), liquid water (qlw) and initial cloud
!    radius (cloud_base_radius), and the cloud radius and vertical vel-
!    ocity tendencies returned from subroutine simult.
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 3e20.12)')   &
             'in cloudm: QE,QLW,RR= ',mrebot,QLW,Param%cloud_base_radius
        write (diag_unit, '(a, 2e20.12)')  &
                                'in cloudm: DPD,DWDP= ',drdp, dwdp
      endif

!--------------------------------------------------------------------
!    define estimated values for t, w and cloud radius based on the
!    derivatives returned from subroutine simult. define the average of
!    the original cloud radius and the new estimate (rbar).
!--------------------------------------------------------------------
      tcest = tccbot + dtdp*Param%dp_of_cloud_model
      west  = wvbot  + dwdp*Param%dp_of_cloud_model
      rest  = rclbot + drdp*Param%dp_of_cloud_model
      rbar  = 0.5*(rest + rclbot)

!--------------------------------------------------------------------
!    if in diagnostics column, output the newly estimated values of 
!    cloud radius (rest) and vertical velocity (west).
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)')  &
                         'in cloudm: rest,west= ',rest, west
      endif

!--------------------------------------------------------------------
!    if any of these estimated values are incompatible with the exist-
!    ence of deep convection at this level, exit this vertical loop --
!    cloud does not extend any higher.
!--------------------------------------------------------------------
      if (rbar <= Param%rbound .or. &    
          west <= Param%wbound .or. &     
          rest <= Param%rbound) then
        flag = .true.
        return
      endif

!---------------------------------------------------------------------
!    define estimated values of the liquid water components to be passed
!    to subroutine micro to hold first pass values.
!---------------------------------------------------------------------
      qrwest = qrw
      qcwest = qcw
      qlwest = qlw

!---------------------------------------------------------------------
!    call micro to calculate the microphysical terms which occur during
!    the parcel movement from cld_press(k) to cld_press(k+1).
!---------------------------------------------------------------------
      call don_cm_micro_k   &
           (diag_unit, debug_ijt, Param, tccbot, tcest, pbot, ptop,  &
            tebot, tetop, mrebot, mretop, wvbot, west, rbar, rmu,  &
            qrwest, qcwest, qlwest, dcw1, dqrw3, ermesg, error)  
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!--------------------------------------------------------------------
!    if in diagnostics column, output the newly estimated values of 
!    environmental pressure (cld_press(k+1)), temperature (te) and 
!    mixing ratio (mre), and cloud temperature (tcest), cloud radius 
!    (rest) and liquid water (qlwest)
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, f19.10, f20.14, e20.12)') &
             'in cloudm: P,TEST,QEST= ',ptop     ,TEtop  ,mretop  
        write (diag_unit, '(a, f20.14, 2e20.12)')   &
                   'in cloudm: TCEST,rest,QLWT= ',TCEST,REST,qlwest
      endif

!----------------------------------------------------------------------
!    define the entrainment coefficient (rmu). it is set to 0.0 if the
!    parcel is above the lcl of ensemble member #1, and above the cur-
!    rent ensemble member's cloud base, and if the vertical velocity
!    is below the detrainment threshold; i.e., the current cloud is in
!    its detraining region. if it is not in its detraining zone, set the
!    entrainment coefficient to be the specified entrainment coefficient
!    divided by current cloud radius.
!----------------------------------------------------------------------
      if ((ptop <= pcsave) .and. (.not. lfc_not_reached) .and. &
          (west <= Param%wdet)) then
        rmu = 0.
      else
        rmu = 2.*alpp/rest
      endif 

!----------------------------------------------------------------------
!    call simult to solve for the in-cloud derivatives dT/dp (dtdp2), 
!    dw/dp (dwdp2) and dr/dp (drdp2) using the previously estimated
!    values at p + delta p, including the microphysical contributions.
!----------------------------------------------------------------------
      call don_cm_simult_k   &
           (diag_unit, debug_ijt, lfc_not_reached, Param, pcsave, rmu, &
            tcest, rest, west, ptop, qlwest, tetop, mretop, dtdp2,  &
            drdp2, dwdp2, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!--------------------------------------------------------------------
!    if in diagnostics column, output the lfc flag (lfc_not_reached), and
!    the temperature, cloud radius and vertical velocity tendencies 
!    returned from subroutine simult.
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, l4, 2e20.12)')  &
            'in cloudm: TESTLC,DTDP2,DPD2= ',lfc_not_reached,DTDP2,drdp2
        write (diag_unit, '(a, e20.12)')  'in cloudm: DWDP2= ',DWDP2
      endif

!--------------------------------------------------------------------
!    define new values of cloud radius and vertical velocity at level
!    k+1 using the values at level k and the vertical derivatives
!    returned from subroutine simult. if either of these values is below
!    its acceptable bound, set both values to 0.0 and exit the loop -- 
!    the cloud top has been reached.
!--------------------------------------------------------------------
      rcltop   = rclbot + drdp2*Param%dp_of_cloud_model
      wvtop    = wvbot  + dwdp2*Param%dp_of_cloud_model
      if ((wvtop <= Param%wbound) .or. (rcltop   <= Param%rbound)) then 
        rcltop   = 0.0
        wvtop    = 0.
        flag = .true.
        return
      endif

!--------------------------------------------------------------------
!    average the two derivatives calculated for temperature, vertical
!    velocity and cloud radius. define the cloud temperature at the
!    k+1 level as the value at k + the  effect of the calculated deriv-
!    ative.
!--------------------------------------------------------------------
      dtdp = 0.5*(dtdp + dtdp2)
      dwdp = 0.5*(dwdp2 + dwdp)
      drdp = 0.5*(drdp + drdp2)
      tcctop   = tccbot + dtdp*Param%dp_of_cloud_model

!--------------------------------------------------------------------
!    call subroutine freeze_liquid to define the amount of liquid in 
!    the updraft is frozen in the current layer.
!--------------------------------------------------------------------
      call don_cm_freeze_liquid_k   &
           (k, diag_unit, debug_ijt, Param, tccbot, tcctop, qlwest, &
            dfrac, dtfr, dtupa, dfrtop, sumfrea, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!--------------------------------------------------------------------
!    add the effect of droplet freezing to the cloud temperature.
!--------------------------------------------------------------------
      tcctop   = tcctop   + dfrtop  

!--------------------------------------------------------------------
!    define the values of vertical velocity and cloud radius at level 
!    k+1 using the averaged values of dw/dp and dr/dp. if either of 
!    these values is below its acceptable bound, set both values to 0.0
!    and exit the loop -- the cloud top has been reached.
!--------------------------------------------------------------------
      wvtop    = wvbot  + dwdp*Param%dp_of_cloud_model
      rcltop   = rclbot + drdp*Param%dp_of_cloud_model
      if ((wvtop <= Param%wbound) .or. (rcltop <= Param%rbound)) then 
        rcltop   = 0.0
        wvtop    = 0.
        flag = .true.
        return
      endif
 
!--------------------------------------------------------------------
!    define the pressure at which nsemble member #l becomes buoyant (its
!    level of free convection). all ensemble members must have their 
!    cloud top pressure lower than this value.
!--------------------------------------------------------------------
      if (kou == 1) then
        if ((lfc_not_reached) .and. (wvtop > wvbot))  then
          pcsave = ptop     
        endif
      endif

!--------------------------------------------------------------------
!    check to see if the level of free convection has been reached. if
!    so, set the logical variable indicating such (lfc_not_reached) to  
!    be .false..
!--------------------------------------------------------------------
      if (wvtop   > wvbot)  then
        lfc_not_reached = .false.
      endif

!--------------------------------------------------------------------
!    define the vapor mixing ratio for the new cloud temperature at
!    level k+1.
!--------------------------------------------------------------------
      call compute_mrs_k (tcctop, ptop, Param%D622, Param%D608,  &
                          rsctop, nbad)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (nbad /= 0) then
        ermesg = 'subroutine don_cm_move_parcel_k: '// &
                 'temperatures out of range of esat table'
        error = 1
        return
      endif

!--------------------------------------------------------------------
!    if in diagnostics column, output the values of environmental 
!    temperature (te) and mixing ratio (mre) and in-cloud mixing ratio 
!    (rsc) at level k+1.
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, f20.14, 2e20.12)')  &
                   'in cloudm: TE,QE,RSC= ',tetop  , mretop  , rsctop  
      endif

!--------------------------------------------------------------------
!    define the layer-mean cloud radius (rbar) and entrainment coef-
!    ficient (rmub).  if the cloud radius is less than the lower bound,
!    exit the loop -- the cloud top has been reached.
!--------------------------------------------------------------------
      rbar = 0.5*(rclbot + rcltop  )
      if (rbar <= Param%rbound) then
        flag = .true.
        return   
      endif

!---------------------------------------------------------------------
!    call micro to calculate the microphysical terms which occur during
!    the parcel movement from p to p + delta p using the second itera-
!    tion values for the cloud conditions at level k + 1.
!---------------------------------------------------------------------
      call don_cm_micro_k   &
           (diag_unit, debug_ijt, Param, tccbot, tcctop, pbot, ptop, &
            tebot, tetop, mrebot, mretop, wvbot, wvtop, rbar, rmu,  &
            qrw, qcw, qlw, dcw1, dqrw3, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    save the final values of cloud water, rainwater and total liquid
!    that are found at level k+1 in arrays disc, disd and qlwa.
!---------------------------------------------------------------------
      qlwatop = qlw

!--------------------------------------------------------------------
!    if in diagnostics column, output the values of in-cloud temperature
!    (tcc), vertical velocity (wv), liquid water (qlw), condensation
!    (dcw1) and precipitation (dqrw3) at level k+1.
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, f20.14, 2e20.12)') &
                     'in cloudm: TCC(k+1),WV(k+1),QLW= ',tcctop  ,  &
                      wvtop  , qlw
        write (diag_unit, '(a, 2e20.12)')  &
                        'in cloudm: DCW1,DQRW3= ',dcw1, dqrw3
      endif

!---------------------------------------------------------------------


end subroutine don_cm_move_parcel_k




!######################################################################

subroutine don_cm_lcl_k    &
         (Param, t_init, p_init, mr_init, t_lcl, p_lcl, mr_lcl,  &
          lcl_reached, ermesg, error)

!---------------------------------------------------------------------
!    subroutine don_cm_lcl_k computes the lifting conden-
!    sation level by raising a parcel with temperature t_init and vapor 
!    mixing ratio mr_init adiabatically from pressure p_init until satur-
!    ation is reached at pressure p_lcl, temperature t_lcl and vapor 
!    mixing ratio mr_lcl. a flag lcl_reached is set to .true. to indicate
!    that the lcl was successfully reached.
!---------------------------------------------------------------------

use donner_types_mod, only : donner_param_type
use sat_vapor_pres_k_mod, only: compute_mrs_k

implicit none

!-----------------------------------------------------------------------
type(donner_param_type), intent(in)   :: Param
real,                    intent(in)   :: t_init, p_init, mr_init
real,                    intent(out)  :: t_lcl, p_lcl, mr_lcl
logical,                 intent(out)  :: lcl_reached
character(len=*),        intent(out)  :: ermesg
integer,                 intent(out)  :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      t_init         temperature of parcel at its launch point
!                     [ deg K ]
!      p_init         pressure of the parcel at its launch point [ Pa ]
!      mr_init        mixing ratio of the parcel at its launch point
!                     [ kg(h2o) / kg(air) ]
!
!   intent(out) variables:
!
!      t_lcl          temperature at lifting condensation level 
!                     [ deg K ]
!      p_lcl          pressure at lifting condensation level [ Pa ]
!      mr_lcl         mixing ratio of parcel at lifting condensation
!                     level [ kg(h2o) / kg(air) ]
!      lcl_reached    lcl was reached for this parcel ?
!
!----------------------------------------------------------------------


!---------------------------------------------------------------------
!   local variables:

      real    ::  t_parcel, p_start, p_end, gamma, rs, kappa_moist
      integer ::  max_levels
      integer ::  k, nbad

!---------------------------------------------------------------------
!   local variables:
!
!      t_parcel     temperature of parcel at current level
!      p_start      pressure at start of current iteration [ Pa ]
!      p_end        pressure at end of current iteration [ Pa ]
!      gamma        adiabatic lapse rate for moist air at current 
!                   location [ deg K / Pa ] 
!      es           saturation vapor pressure at temperature t_parcel 
!                   [ Pa ]
!      rs           saturation mixing ratio at t_parcel and p_end 
!                   [ kg(h2o) / kg(air) ]
!      kappa_moist  exponent in expression for moist potential 
!                   temperature [ dimensionless ]
!      max_levels   maximum number of iterations of parcel movement
!                   which must be taken before knowing that the
!                   parcel will not undergo deep convection
!      k            do-loop index
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' '  ; error = 0

!---------------------------------------------------------------------
!    define an approximation for kappa (R/Cp) for moist air 
!    (kappa_moist). 
!---------------------------------------------------------------------
      kappa_moist = (Param%RDGAS/Param%CP_AIR)*(1. +    &
                    (Param%RVGAS/Param%RDGAS -     &
                     Param%CP_VAPOR/Param%CP_AIR)*mr_init)
      
!---------------------------------------------------------------------
!    define the initial values for t_parcel (the parcel temperature) 
!    and p_start and p_end, which define the pressure increment through
!    which the parcel rises during the first iteration.
!---------------------------------------------------------------------
      t_parcel = t_init
      p_start  = p_init
      p_end    = p_start + Param%parcel_dp

!---------------------------------------------------------------------
!    initialize the output variables.
!---------------------------------------------------------------------
      p_lcl  = 0.
      t_lcl  = 0.
      mr_lcl = 0.
      lcl_reached = .false.

!---------------------------------------------------------------------
!    define the maximum number of increments (max_levels) that the 
!    parcel may be raised without reaching condensation before it is 
!    no longer viable as a deep convection parcel.
!---------------------------------------------------------------------
      max_levels = int( (p_init - Param%upper_limit_for_lfc)/  &
                        abs(Param%parcel_dp) ) + 1

!----------------------------------------------------------------------
!    raise the parcel in increments of parcel_dp, following the adia-
!    batic lapse rate for moist air, until either the lcl is reached
!    or the parcel is no longer capable of undergoing deep convection. 
!----------------------------------------------------------------------
      do k=1,max_levels

!---------------------------------------------------------------------
!    exit the loop if the parcel has risen beyond where deep convection
!    is possible.
!---------------------------------------------------------------------
        if (p_end < Param%upper_limit_for_lfc) exit

!---------------------------------------------------------------------
!    define the lapse rate of temp with respect to pressure for moist
!    air (gamma). update the parcel temperature to the value at the
!    top of this layer.
!---------------------------------------------------------------------
        gamma = kappa_moist*t_parcel/p_start
        t_parcel = t_parcel + gamma*Param%parcel_dp 

!---------------------------------------------------------------------
!    determine the saturation mixing ratio for this temperature and
!    pressure. 
!---------------------------------------------------------------------
        call compute_mrs_k (t_parcel, p_end,                  &
                            Param%d622 , Param%d608 , rs, nbad, &
                            mr = mr_init)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (nbad /= 0) then
          ermesg = 'subroutine don_cm_lcl_k: '// &
                   'temperatures out of range of esat table'
          error = 1
          return
        endif

!---------------------------------------------------------------------
!   determine if the parcel is saturated. if it is, define the lcl 
!   values of temperature (t_lcl), mixing ratio (mr_lcl) and pressure 
!   (p_lcl), and return to the calling routine.
!---------------------------------------------------------------------
        if (rs <= mr_init) then
          p_lcl = p_end
          t_lcl = t_parcel
          mr_lcl = mr_init                     
          lcl_reached = .true.
          exit
        endif

!---------------------------------------------------------------------
!    define the starting and ending pressures for the next iteration.
!---------------------------------------------------------------------
        p_start = p_end
        p_end   = p_end + Param%parcel_dp
      end do

!---------------------------------------------------------------------


 end subroutine don_cm_lcl_k  


!#####################################################################


subroutine don_cm_mesub_k     &
         (Nml, pfull_c, nlev_lsm, me, diag_unit, debug_ijt, Param, cu,           &
          ci_liq_cond, ci_ice_cond, pmelt_lsm, cell_precip, &
          dint, plzb_c, pb, pt_kou, temp_c, phalf_c,   &
          ca_liq, ca_ice, ecd, ecd_liq, ecd_ice, ecei_liq,   &
          ece, ece_liq, ece_ice, meso_freeze, meso_melt, ermesg, error)

!----------------------------------------------------------------------
!    subroutine mesub calculates mesoscale heat and moisture sources,
!    using a variation on the Leary and Houze (JAS, 1980) procedure.
!    the defined fields are condensate transferred from cell to anvil
!    (ca), condensate evaporated in convective downdrafts (ecd), conden-
!    sate evaporated in convective updrafts (ece), the condensate 
!    entering the anvil which has not yet been frozen (meso_freeze), and 
!    the amount of condensate which must be melted in the mesoscale down-
!    draft to assure ice conservation (meso_melt). the subroutine is 
!    called separately for each ensemble member. for notation, see 
!    "Cu Closure A notes," 2/97.
!----------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_nml_type

implicit none

!----------------------------------------------------------------------
type(donner_nml_type),         intent(in)    :: Nml
integer,                       intent(in)    :: nlev_lsm, me, diag_unit
logical,                       intent(in)    :: debug_ijt
type(donner_param_type),       intent(in)    :: Param
real,                          intent(in)    :: cu, cell_precip, dint, &
                                                plzb_c, pb, pt_kou
real,   dimension(nlev_lsm),   intent(in)    :: temp_c, pfull_c
real,   dimension(nlev_lsm+1), intent(in)    :: phalf_c
real,                          intent(out)   :: ca_liq, ca_ice
real,                          intent(in)    :: pmelt_lsm, &
                                                ci_liq_cond, ci_ice_cond
real,   dimension(nlev_lsm),   intent(out)   :: ecd, ece, meso_freeze, &
                                                meso_melt, &
                                                ecd_liq, ecd_ice, &
                                                ece_liq, ece_ice
real,                          intent(out)   :: ecei_liq
character(len=*),              intent(out)   :: ermesg
integer,                       intent(out)   :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       cu           column integrated condensation integral
!                    [ mm / day ]
!       cell_precip  column integrated precipitation integral
!                    [ mm / day ]
!       dint???      water mass frozen in convective updraft
!            ??????  plus ice deposited convective updraft
!                    [ kg(h2o) /( (m**2) sec) ]
!                    weighted as cu,cell_precip
!       plzb_c       pressure at level of zero buoyancy [ Pa ]
!       ps           surface pressure [ Pa ]
!       pb           cloud-base pressure [ Pa ]
!       pt_kou       cloud-top pressure [ Pa ]
!       pmelt_lsm    pressure at bottom of layer in which melting 
!                    begins   [ Pa ]
!       phalf_c      large-scale model pressure half-levels (Pa)
!       debug_ijt    is this a diagnostics column ?
!       diag_unit    output unit number for this diagnostics column
!
!   intent(out) variables:
!
!       ca           total condensate transfered from cells to anvil 
!                    by this ensemble member [ mm/day ]
!       ecd          profile of condensate evaporated in convective
!                    downdraft on large-scale model grid 
!                    [ g(h2o) / kg(air) / day ] 
!       ece          profile of condensate evaporated in convective 
!                    updraft on large-scale model grid 
!                    [ g(h2o) / kg(air) / day ] 
!       meso_freeze  profile of condensate which is frozen upon enter-
!                    ing the anvil on the large-scale grid
!                    [ g(h2o) / kg(air) / day ] 
!       meso_melt    profile of condensate which is melted in mesoscale
!                    downdraft on large-scale model grid
!                    [ g(h2o) / kg(air) / day ] 
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
 
      integer ::     k
      real    ::  avail_meso_cd     ! fraction of column integrated
                                    ! condensation available to meso-
                                    ! scale circulation (1. - gnu)
                                    ! [ dimensionless ]
      real    ::  caa               ! amount of condensate which must
                                    ! be frozen when it enters the anvil
                                    ! [ g(h2o) / kg(air) / day ]
      real    ::  dint2             ! amount of condensate which has
                                    ! been frozen in the cumulus updraft
                                    ! before entering the anvil
                                    ! [ g(h2o) / kg(air) / day ]
      real    ::  ecda              ! amount of condensate evaporated 
                                    ! in cumulus downdrafts
                                    ! [ g(h2o) / kg(air) / day ]
      real :: ecda_liq, ecda_ice
      real    ::  ecdi              ! amount of condensate evaporated 
                                    ! in cumulus downdrafts [ mm / day ]
      real :: ecdi_liq, ecdi_ice
      real    ::  ecea              ! amount of condensate evaporated 
                                    ! in cumulus updrafts 
                                    ! [ g(h2o) / kg(air) / day ]
      real :: ecea_liq, ecea_ice
      real    ::  ecei              ! amount of condensate evaporated 
                                    ! in cumulus updrafts [ mm / day ]
      real ::           ecei_ice
      real    ::  elta              ! amount of condensate which must
                                    ! be melted in the mesoscale down-
                                    ! draft to conserve ice mass
                                    ! [ g(h2o) / kg(air) / day ]
      real    ::  gnu               ! fraction of column integrated 
                                    ! condensation which precipitates
                                    ! out [ dimensionless ]
      real    ::  ptt               ! pressure one cloud model delta p 
                                    ! above cloud top [ Pa ]
      real    ::  pzm               ! pressure at base of mesoscale 
                                    ! circulation [ Pa ]
      real    ::  pztm              ! pressure at top of mesoscale cir-
                                    ! culation [ Pa ]
      real    ::  p1                ! lower pressure limit for the layer
                                    ! in which one of the physical
                                    ! processes is occurring [ Pa ]
      real    ::  p2                ! upper pressure limit for the layer
                                    ! in which one of the physical
                                    ! processes is occurring [ Pa ]
      integer :: itrop
      real    :: ptrop

!---------------------------------------------------------------------
!   local variables:
!
!      
      ermesg = '  ' ; error = 0

!---------------------------------------------------------------------
!    define pressure one cloud-model level above cloud top (ptt). 
!    define the pressure at top of mesoscale updraft (pztm, 300 hPa 
!    plus one model-layer pressure thickness above cloud top).
!---------------------------------------------------------------------
      ptt = pt_kou + Param%dp_of_cloud_model
      pztm = ptt - 300.E02

!---------------------------------------------------------------------
!    restrict pztm to >= 100 hPa, cf Ackerman et al (JAS,1988), unless 
!    pt_kou <= 100 hPa. it was found in AM2p9 that the stratospheric 
!    water vapor was excessive with this pztm restriction, so pztm is now
!    set to be no higher than the level of zero buoyancy, or if the
!    cloud top is above the level of zero buoyancy, it is set to one 
!    model layer above the level of zero buoyancy. 
!---------------------------------------------------------------------
      if (pztm < plzb_c) pztm = plzb_c
      if (ptt < plzb_c)  pztm = plzb_c + Param%dp_of_cloud_model

      if (Nml%limit_pztm_to_tropo) then
        call find_tropopause (nlev_lsm, temp_c, pfull_c, ptrop, itrop)
        pztm = MAX (pztm, ptrop)
      endif

!---------------------------------------------------------------------
!    define the base of the mesoscale updraft (pzm), as the layer imm-
!    ediately above cloud top, or, if the top of the mesoscale updraft
!    has been redefined to be at or just above the level of zero 
!    buoyancy, to be one layer below the mesoscale updraft top. 
!---------------------------------------------------------------------
      pzm = ptt
      if (pzm <= pztm) pzm = pztm - Param%dp_of_cloud_model

!---------------------------------------------------------------------
!    if in a diagnostics column, output the convective rain 
!    (cell_precip), convective updraft condensation (cu), and the pres-
!    sure at the level of zero buoyancy (plzb_c).
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)') 'in mesub: rc,cu= ',  &
                                                    cell_precip, cu
        write (diag_unit, '(a,  e20.12)') 'in mesub: plzb = ',plzb_c
      endif

!----------------------------------------------------------------------
!    define the ratio of precipitation to condensation for the current
!    ensemble member (gnu). define the remaining fraction of condens-
!    ation 1 - gnu as the condensate available to the mesoscale circ-
!    ulation (avail_meso_cd). define the mass of this available conden-
!    sate which is evaporated in convective downdrafts (ecdi), the mass
!    evaporated into the cell environment (ecei) and the portion incor-
!    porated into the mesoscale region (ca). this partitioning is 
!    defined by the parameters evap_in_downdraft, evap_in_environ and 
!    entrained_into_meso, taken from the work of Leary and Houze 
!    (JAS, 1980).
!----------------------------------------------------------------------
      gnu = cell_precip/cu
      avail_meso_cd = 1. - gnu
      ecdi  = (Param%evap_in_downdrafts*avail_meso_cd)*cu
      ecdi_liq  = (Param%evap_in_downdrafts*avail_meso_cd)* &
                         (Param%seconds_per_day*ci_liq_cond)
      ecdi_ice  = (Param%evap_in_downdrafts*avail_meso_cd)*     &
                         (Param%seconds_per_day*ci_ice_cond)
      ecei  = (Param%evap_in_environ*avail_meso_cd)*cu
      ecei_liq  = (Param%evap_in_environ*avail_meso_cd)*    &
                         (Param%seconds_per_day*ci_liq_cond)
      ecei_ice  = (Param%evap_in_environ*avail_meso_cd)*     &
                         (Param%seconds_per_day*ci_ice_cond)
      ca_liq    = (Param%entrained_into_meso*avail_meso_cd)*   &
                         (Param%seconds_per_day*ci_liq_cond)
      ca_ice    = (Param%entrained_into_meso*avail_meso_cd)*   &
                         (Param%seconds_per_day*ci_ice_cond)
      if (debug_ijt) then
        write (diag_unit, '(a, 3e20.12)')  &
               'in mesub: cu, h1_liqintg, h1_iceintg= ',  &
            cu, ci_liq_cond*Param%seconds_per_day      ,  &
                 ci_ice_cond*Param%seconds_per_day      
      endif

!---------------------------------------------------------------------
!    if in a diagnostics column, output the ratio of convective rain 
!    to convective updraft condensation (gnu) and the mass entrained
!    into the mesoscale region (ca).
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12)')  'in mesub: gnu= ',gnu
        write (diag_unit, '(a, e20.12)') 'in mesub: ca= ',  &
                                                       ca_liq + ca_ice
        write (diag_unit, '(a, 2e20.12)') 'in mesub: ca_liq,ca_ice= ', &
                                                ca_liq, ca_ice
      endif

!--------------------------------------------------------------------
!    calculate the mass of water which must be frozen as it enters the
!    mesoscale anvil (caa). if no freezing has occurred in the cumulus
!    updraft (i.e., dint2 = 0) then this will be ca, the total mass 
!    available to the anvil. if freezing has occurred, (ie, 
!    dint2 /= 0.), then the amount to be frozen is the total amount 
!    available (ca) plus additional vapor mass deposited on the ice in 
!    the updraft (ecei), less that which has already frozen (dints). 
!    dints and caa are expressed in units of g(h2o) per kg(air) per day.
!--------------------------------------------------------------------
!9/15/07, 1037AM:
      dint2 = avail_meso_cd*(dint                               )*  &
               8.64e07*Param%grav/(pzm - pztm)

      if (dint2 /= 0.)  then 
       caa = ((ca_liq + ecei_liq)*Param%grav*1000./(pzm - pztm)) - dint2
      else
        caa = ca_liq*Param%grav*1000./(pzm - pztm)
      endif

!---------------------------------------------------------------------
!    if in a diagnostics column, output the previously frozen condensate
!    (dint2), the additional amount to be frozen (caa) and the pressure
!    range over which the freezing will occur (pzm, pztm). if 
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a,  e20.12)')  &
                         'in mesub:     dint           =',    dint 
        write (diag_unit, '(a, 2e20.12)')  &
                         'in mesub:     dint2, ecei_liq=',    dint2, &
                                                          ecei_liq
        write (diag_unit, '(a, 3e20.12)')  &
                           'in mesub: caa,pzm,pztm= ',caa,pzm,pztm
      endif

!---------------------------------------------------------------------
!    if there is additional condensate which must be frozen upon enter-
!    ing the anvil, call map_hi_res_intgl_to_lo_res_col to spread this 
!    additional freezing uniformly over the region between anvil base 
!    (pzm) and anvil top (pztm) in the large-scale model. store the out-
!    put in array meso_freeze. if no additional freezing is needed, set 
!    meso_freeze to be 0.0.
!---------------------------------------------------------------------
      if (caa > 0.)  then 
        if (debug_ijt) then
          write (diag_unit, '(a, e20.12, 2f19.10)')  &
                      'in cm_intgl_to_gcm_col: xav,p1,p2= ',caa, pzm, &
                                                  pztm
        endif
        call don_u_map_hires_i_to_lores_c_k   &
             (nlev_lsm, caa, pzm, pztm, phalf_c, meso_freeze, ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return
        if (debug_ijt) then
          do k=1,nlev_lsm       
            if (meso_freeze(k) /= 0.0) then
              write (diag_unit, '(a, i4, e20.12)') &
                    'in cm_intgl_to_gcm_col: k,x= ',k, meso_freeze   (k)
            endif
          end do
        endif
      else
        meso_freeze = 0.
      endif

!---------------------------------------------------------------------
!    define the evaporation which occurs in the convective downdraft.
!    the convective downdraft is assumed to originate one layer above
!    the cloud top (ptt) and extend to the surface (phalf_c(1)). 
!    convert the convective downdraft evaporation to units of
!    g(h20) / kg(air) per day.
!---------------------------------------------------------------------
      ecda = ecdi*Param%grav*1000./(phalf_c(1) - ptt)
      ecda_liq = ecdi_liq*Param%grav*1000./(phalf_c(1) - ptt)
      ecda_ice = ecdi_ice*Param%grav*1000./(phalf_c(1) - ptt)

!---------------------------------------------------------------------
!    if in a diagnostics column, output the convective downdraft evap-
!    oration (ecda) and the large-scale model pressure limits over which
!    this evaporation occurs (phalf_c(1), ptt).
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 3e20.12)')  &
                         'in mesub: ecda,p1,pz0= ',ecda,phalf_c(1),ptt
        write (diag_unit, '(a, 2e20.12)')  &
                         'in mesub: ecda_liq, ecda_ice= ',  &
                                ecda_liq, ecda_ice
      endif

!---------------------------------------------------------------------
!    call map_hi_res_intgl_to_lo_res_col to spread the integrated evap-
!    oration in convective downdrafts uniformly over the region between
!    the surface (phalf_c(1)) and the anvil base (pzm) and the top of 
!    cloud (ptt). output field is ecd.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12, 2f19.10)')  &
             'in cm_intgl_to_gcm_col: xav,p1,p2= ',ecda, phalf_c(1) , &
                                                  ptt 
      endif
      call don_u_map_hires_i_to_lores_c_k   &
           (nlev_lsm, ecda, phalf_c(1), ptt, phalf_c, ecd, ermesg, error)
      call don_u_map_hires_i_to_lores_c_k   &
         (nlev_lsm, ecda_liq, phalf_c(1), ptt, phalf_c, ecd_liq, ermesg, error)
      call don_u_map_hires_i_to_lores_c_k   &
         (nlev_lsm, ecda_ice, phalf_c(1), ptt, phalf_c, ecd_ice, ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return
      if (debug_ijt) then
        do k=1,nlev_lsm       
          if (ecd(k) /= 0.0) then
            write (diag_unit, '(a, i4, e20.12)') &
                'in cm_intgl_to_gcm_col: k,ecd= ',k, ecd   (k)
            write (diag_unit, '(a, i4, 2e20.12)') &
              'in cm_intgl_to_gcm_col: k,ecdliq,ecdice= ',k, &
                       ecd_liq(k), ecd_ice(k)
          endif
        end do
      endif

!---------------------------------------------------------------------
!    be sure that the melting level in the large-scale model (pmelt_lsm)
!    is below the top of the mesoscale circulation (pztm),and above
!    cloud base (pb). if not, no melting will occur; set p2 to be 0.0.
!---------------------------------------------------------------------
      elta = 0.
      if (pmelt_lsm  < pztm                    )  then
        meso_melt = 0.
      if (debug_ijt) then
        write (diag_unit, '(a, 2f19.10)') &
                 ' NO MELTING DONE: melting level above top of &
                    &mesoscale circulation : pmelt_lsm,pztm',     &
                                                pmelt_lsm, pztm      
      endif

!---------------------------------------------------------------------
!    if pmelt_lsm is within the region of the cloud and mesoscale circ-
!    ulation, calculate any melting that must occur in the mesoscale
!    downdraft in order to conserve ice mass; ie, if the amount to be
!    frozen was calculated as more than the available condensate, then
!    the excess must be melted, and is done so in the mesoscale down-
!    draft between the melting level and cloud base.
!---------------------------------------------------------------------
      else if (pmelt_lsm >= pztm .and. pmelt_lsm <= pb) then
        p2 = pmelt_lsm
        p1 = pb
        if (caa <= 0.) then 
          caa = -caa*(pzm - pztm)/(pb - p2)
          elta = caa
        endif
      if (debug_ijt) then
        write (diag_unit, '(a, 3f19.10)') &
                   'MELTING DONE: pmelt_lsm,pb,caa  ',pmelt_lsm, pb, &
                                                            caa  
      endif

!---------------------------------------------------------------------
!    if in diagnostics column, output the melting (elta) and the 
!    pressures defining the layer in which it occurs (pb, p2)
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12, 2f19.10)') &
                           'in mesub: elta,p1,p2= ',elta,p1,p2
      endif

!---------------------------------------------------------------------
!    call map_hi_res_intgl_to_lo_res_col to spread the required melting
!    resulting from excessive freezing over the layer between cloud base
!    and the melting level. output field is meso_melt.
!---------------------------------------------------------------------
      call don_u_map_hires_i_to_lores_c_k   &
           (nlev_lsm, elta, p1, p2, phalf_c, meso_melt, ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return
      if (debug_ijt) then
        do k=1,nlev_lsm       
          if (meso_melt(k) /= 0.0) then
            write (diag_unit, '(a, i4, e20.12)') &
                 'in cm_intgl_to_gcm_col: k,meso_melt= ',k, meso_melt(k)
          endif
        end do
      endif

      else if (pmelt_lsm > pb) then
        meso_melt = 0.
        if (pmelt_lsm == phalf_c(1)) then
          if (debug_ijt) then
             write (diag_unit, '(a)') &
                   'NO MELTING LEVEL PRESENT IN COLUMN'
          endif
        else
! melt below cloud base 
      if (debug_ijt) then
        write (diag_unit, '(a, 2f19.10)') &
            ' NO MELTING DONE: melting level below PB: pmelt_lsm,pb', &
                                                      pmelt_lsm, pb
      endif
      endif
      endif ! (pmelt<pztm or pmelt > pb)

!---------------------------------------------------------------------
!    calculate the evaporation which occurs in the convective 
!    updraft.
!    this is spread between 50 hPa below cloud top and 10 hPa above 
!    cloud top.
!---------------------------------------------------------------------
      p1 = pt_kou + 50.0e02
      p2 = ptt
      ecea = ecei*Param%grav*1000./(p1-p2)
      ecea_liq = ecei_liq*Param%grav*1000./(p1-p2)
      ecea_ice = ecei_ice*Param%grav*1000./(p1-p2)

!---------------------------------------------------------------------
!    if in diagnostics column, output the convective updraft evaporation
!    (ecea, ecei) and the large-scale model pressure layer limits over 
!    which it occurs (p1, p2).
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)')  &
                         'in mesub: ecea,ecei= ',ecea, ecei
        write (diag_unit, '(a, 2e20.12)')  &
                         'in mesub: LIQecea,ecei= ',ecea_liq, ecei_liq
        write (diag_unit, '(a, 2e20.12)')  &
                         'in mesub: ICEecea,ecei= ',ecea_ice, ecei_ice
        write (diag_unit, '(a, e20.12, 2f19.10)')  &
                         'in mesub: ecea,p1,p2= ',ecea, p1, p2
      endif

!---------------------------------------------------------------------
!    call map_hi_res_intgl_to_lo_res_col to spread the integrated evap-
!    oration in convective updrafts uniformly over the designated 
!    region.  output field is ece.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12, 2f19.10)')  &
                      'in cm_intgl_to_gcm_col: xav,p1,p2= ',ecea, p1, &
                                                  p2
      endif
      call don_u_map_hires_i_to_lores_c_k   &
           (nlev_lsm, ecea, p1, p2, phalf_c, ece, ermesg, error)
      call don_u_map_hires_i_to_lores_c_k   &
           (nlev_lsm, ecea_liq, p1, p2, phalf_c, ece_liq, ermesg, error)
      call don_u_map_hires_i_to_lores_c_k   &
           (nlev_lsm, ecea_ice, p1, p2, phalf_c, ece_ice, ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return
      if (debug_ijt) then
        do k=1,nlev_lsm     
          if (ece(k) /= 0.0) then
            write (diag_unit, '(a, i4, e20.12)') &
                           'in cm_intgl_to_gcm_col: k,x= ',k, ece   (k)
          endif
        end do
      endif

!---------------------------------------------------------------------


end subroutine don_cm_mesub_k



!######################################################################

subroutine don_cm_compute_vert_fluxes_k   &
         (nlev_hires, ntr, ncc_kou, kou, diag_unit, debug_ijt, &
          do_donner_tracer, Param, pt_kou, cld_press, rcl, te, mre, wv, &
          tcc, dpf, xclo, xtrae, apt, efchr, emfhr, etfhr, ermesg, error)

!---------------------------------------------------------------------
!    subroutine compute_vertical_fluxes computes vertical flux conver-
!    gence terms and their column integrals for temperature, entropy, 
!    moisture and tracers within the cloud. 
!---------------------------------------------------------------------

use donner_types_mod, only : donner_param_type
use sat_vapor_pres_k_mod, only: compute_qs_k

implicit none 

!---------------------------------------------------------------------
integer,                            intent(in)     :: nlev_hires, ntr,  &
                                                      ncc_kou, kou,  &
                                                      diag_unit
logical,                            intent(in)     :: debug_ijt,   &
                                                      do_donner_tracer
type(donner_param_type),            intent(in)     :: Param
real,                               intent(in)     :: pt_kou
real,    dimension(nlev_hires),     intent(in)     :: cld_press, rcl,  &
                                                      te, mre, wv, tcc, &
                                                      dpf
real,    dimension(nlev_hires,ntr), intent(in)     :: xclo, xtrae      
real,                               intent(out)    :: apt               
real,    dimension(nlev_hires),     intent(out)    :: efchr, emfhr    
real,    dimension(nlev_hires,ntr), intent(out)    :: etfhr            
character(len=*),                   intent(out)    :: ermesg
integer,                            intent(out)    :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      pt_kou      pressure at cloud top [ Pa ]
!      cld_press   pressure at full levels in cloud model [ Pa ]
!      rcl         vertical profile of cloud radius [ m ]
!      te          environmental temperature profile [ deg K ]
!      mre         environmental vapor mixing ratio profile 
!                  [ kg(h2o) / kg (dry air) ]
!      wv          in-cloud vertical velocity profile [ m / sec ]
!      tcc         in-cloud temperature profile [ deg K ] 
!      xclo        in-cloud tracer profiles 
!                  [ kg(tracer) / kg (dry air) ]
!      xtrae       environmental tracer profiles 
!                  [ kg(tracer) / kg (dry air) ]
!      debug_ijt   logical indicating whether diagnostics are desired
!                  for column
!      ncc_kou     cloud model level index of level at or above 
!                  cloud top
!      kou         index of current ensemble member 
!      diag_unit   output unit for diagnostics for this column
!
!   intent(out) variables:
!
!      apt         ratio of cloud area at cloud top to that at cloud 
!                  base [ dimensionless ]
!      efchr       vertical entropy flux convergence [ deg K / sec ]
!      emfhr       vertical moisture flux convergence 
!                  [ kg(h2o) / ( kg(dry air) sec ] 
!      etfhr       vertical tracer flux convergence 
!                  [ kg(tracer) / ( kg(dry air) sec ] 
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      real, dimension(nlev_hires)       ::  dpdz, q_sat, ehf, emf
      real, dimension(nlev_hires,ntr) ::  etf                      
      real, dimension(ntr)      ::  sumetf
      real                        ::  p, pl, ph, dpp, ptt, exf, &
                                      thetf, tv_env, tv_cld,  &
                                      dpdz_cb, sumemf, sumefc, sumthet
      integer                     ::  kcl, kch
      integer                     ::  kc, kcont, nbad


!---------------------------------------------------------------------
!   local variables:
!
!        dpdz            value of dp/dz, evaluated from hydrostatic
!                        equation, with virtual mass coefficient to
!                        account for non-hydrostatic effects. see
!                        Donner, 1986, j. atmos. sci., 43, p 2288.
!                        [ kg(air) / (sec**2 m**2) ]
!        q_sat           saturation specific humidity 
!                        [ kg(h2o) / kg (air) ]
!        ehf             vertical temperature flux 
!                        [ ( kg(air) deg K ) / (m sec**3) ]
!        emf             vertical moisture flux
!                        [ kg (h2o) / (m sec**3) ]
!        etf             vertical tracer flux
!                        [  kg(tracer) / (m sec**3) ]
!        sumetf          column integral of vertical tracer flux 
!                        convergence [  kg(tracer) / (m sec**3) ]
!        p               pressure  at level where subgrid scale temper-
!                        ature flux is calculated [ Pa ]
!        pl              pressure at lower flux interface [ Pa ]
!        ph              pressure at upper flux interface [ Pa ]
!        dpp             pressure depth over which flux is being 
!                        calculated [ Pa ]
!        ptt             pressure one level above cloud top [ Pa ]
!        exf             perturbation vertical subgrid scale temp-
!                        erature flux [ ( kg(air) deg K ) / (m sec**3) ]
!        thetf           vertical flux of potential temperature at
!                        cloud base [ ( kg(air) deg K ) / (m sec**3) ]
!        esat            saturation vapor pressure [ Pa ]
!        tv_env          environmental virtual temperature [ deg K ]
!        tv_clda         in-cloud virtual temperature [ deg K ]
!        dpdz_cb         value of dp/dz over the lowest cloud model
!                        pressure layer [ kg(air) / (sec**2 m**2) ]
!        sumemf          column integral of vertical moisture flux 
!                        convergence [ kg (h2o) / (m sec**3) ]
!        sumefc          column integral of vertical entropy flux 
!                        convergence [ ( kg(air) deg K ) / (m sec**3) ]
!        sumthet         column integral of vertical potential temper-
!                        ature flux convergence 
!                        [ ( kg(air) deg K ) / (m sec**3) ]
!        kcl             k index of lower interface level
!        kch             k index of upper interface level
!        kc, kcont       do-loop indices
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!--------------------------------------------------------------------
!    initialize the output arrays.
!--------------------------------------------------------------------
      do kc=1,nlev_hires
        do kcont=1,ntr
          etfhr(kc,kcont) = 0.
        end do
        efchr(kc) = 0.
        emfhr(kc) = 0.
      end do

!---------------------------------------------------------------------
!    initialize variables which will collect column sums.
!---------------------------------------------------------------------
      sumemf    = 0.
      sumefc    = 0.
      sumthet   = 0.
      sumetf(1:ntr) = 0.

!--------------------------------------------------------------------
!    initialize the ratio of cloud area at cloud top to cloud area at 
!    cloud base. 
!--------------------------------------------------------------------
      apt = 0.

!---------------------------------------------------------------------
!    loop over the cloud model levels from cloud base (kc=1) to the 
!    level just below cloud top (ncc_kou-1), defining variables needed 
!    for the flux convergence calculations.
!---------------------------------------------------------------------
      do kc=1,ncc_kou-1 

!---------------------------------------------------------------------
!    define the environmental and in-cloud virtual temperatures. NOTE: 
!    q_sat is correctly a specific humidity; mre(kc) is incorrectly a 
!    mixing ratio.
!---------------------------------------------------------------------
        call compute_qs_k (tcc(kc), cld_press(kc), Param%D622,  &
                           Param%D608, q_sat(kc), nbad)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (nbad /= 0) then
          ermesg = 'subroutine don_cm_compute_vert_fluxes_k: '// &
                   'temperatures out of range of esat table'
          error = 1
          return
        endif
        tv_env     = te(kc)* (1. + Param%d608*(mre(kc)/(1. + mre(kc))))
        tv_cld     = tcc(kc)*(1. + Param%d608*q_sat(kc))

!---------------------------------------------------------------------
!    compute dp/dz (dpdz) as in equation (B3) in Donner (1986), 
!    j. atm. sci., 43, pp 2277-2288.  here the expression uses the 
!    pressure mid-way between cloud base and level 2. compute vertical 
!    eddy transport of temperature (ehf), moisture (emf) and tracers
!    (etf), using equation (3) from Donner et al (1982), j. atm. sci.,
!    39, pp 2159-2181. multiply by dpdz to put flux in pressure (omega) 
!    units. note (1 - a) is approximated as 1. note that q_sat is a 
!    specific humidity while mre is a mixing ratio.
!---------------------------------------------------------------------
        if (kc == 1) then

!--------------------------------------------------------------------
!    at cloud base level, calculate dpdz valid over the interval between
!    level 1 and level 2 (dpdz_cb). this term is needed to determine 
!    the subgrid scale vertical temperature flux valid in the first 
!    cloud layer.
!--------------------------------------------------------------------
          dpdz_cb = -Param%grav*(0.5*(cld_press(1) + cld_press(2)))* &
                    (Param%virt_mass_co*tv_env + tv_cld)/(Param%rdgas*  &
                    (1. + Param%virt_mass_co)*tv_env*tv_cld)

!----------------------------------------------------------------------
!    when at cloud top, calculate dpdz and the cloud perturbation flux 
!    terms that are valid for the model layer including cloud top. save
!    these values as index ncc_kou of the various arrays.
!--------------------------------------------------------------------
        else if (kc == ncc_kou-1) then
          dpdz(ncc_kou) = -Param%grav*pt_kou*(    &
                          Param%virt_mass_co*tv_env + tv_cld)/   &
                          (Param%rdgas*(1. + Param%virt_mass_co)*  &
                          tv_env*tv_cld)
          ehf(ncc_kou) = ((rcl(ncc_kou-1)/Param%cloud_base_radius)**2)* &
                         wv(ncc_kou-1)*dpdz(ncc_kou)*   &
                         (tcc(ncc_kou-1) - te(ncc_kou-1))
          emf(ncc_kou) = ((rcl(ncc_kou-1)/Param%cloud_base_radius)**2)* &
                         wv(ncc_kou-1)*dpdz(ncc_kou)*    &
                         (q_sat(ncc_kou-1)/(1. - q_sat(ncc_kou-1)) - &
                                                         mre(ncc_kou-1))
          do kcont=1,ntr
            etf(ncc_kou,kcont) = ((rcl(ncc_kou-1)/  &
                                 Param%cloud_base_radius)**2)* &
                                 wv(ncc_kou-1)*dpdz(ncc_kou)*  &
                                 (xclo(ncc_kou-1,kcont) -   &
                                                  xtrae(ncc_kou-1,kcont))
          end do
        endif

!---------------------------------------------------------------------
!    compute dpdz and the vertical eddy flux terms at interior cloud
!    levels.
!---------------------------------------------------------------------
        dpdz(kc) = -Param%grav*cld_press(kc)*  &
                   (Param%virt_mass_co*tv_env + tv_cld)/&
                   (Param%rdgas*(1. + Param%virt_mass_co)*tv_env*tv_cld)
        ehf(kc) = ((rcl(kc)/Param%cloud_base_radius)**2)*wv(kc)*  &
                  dpdz(kc)*(tcc(kc) - te(kc))
        emf(kc) = ((rcl(kc)/Param%cloud_base_radius)**2)*wv(kc)*   &
                  dpdz(kc)*(q_sat(kc)/(1. - q_sat(kc)) - mre(kc))
        do kcont=1,ntr
          etf(kc,kcont) = ((rcl(kc)/Param%cloud_base_radius)**2)*  &
                            wv(kc)*dpdz(kc)*  &
                            (xclo(kc,kcont) - xtrae(kc,kcont))
        end do
      end do

!--------------------------------------------------------------------
!    loop over levels from cloud base to the highest level within
!    the cloud to compute the eddy flux convergence profiles.
!--------------------------------------------------------------------
      do kc=1,ncc_kou-1  

!--------------------------------------------------------------------
!    define the lower interface index (kcl) to be one less than 
!    the current level. at cloud base, define the lower interface index
!    to be the current index. define the lower interface pressure (pl).
!    if the lower interface is at or above cloud top, calculation is
!    complete in this column, so return. define the upper interface 
!    pressure (ph) to be the pressure one level above the current index.
!    define the delta p between the interface pressures (dpp).
!--------------------------------------------------------------------
        kcl = MAX (1, kc-1)
        pl = cld_press(kcl)
        if (pl <= pt_kou) return
        ph = cld_press(kc+1)
        dpp = (ph - pl)/2.0

!---------------------------------------------------------------------
!     if the upper interface pressure is not above cloud top, set the 
!     upper interface index (kch) to be one above the current index. set
!     the current level pressure p to be the average value of the higher
!     and lower level pressures, which will be the current level value,
!     except at the lowest level. 
!---------------------------------------------------------------------
        if (ph >= pt_kou) then
          kch = kc + 1
          p = (pl + ph)/2.0

!---------------------------------------------------------------------
!     if the upper interface pressure is above cloud top, set the 
!     upper interface index (kch) to be the current value (kc), and set
!     the upper interface and current level pressures to be the cloud 
!     top pressure.  define the ratio of cloud-top cloud area to cloud-
!     base cloud area (apt).
!---------------------------------------------------------------------
        else
          kch = kc
          ph = pt_kou
          p = pt_kou
          apt = (rcl(kc)/rcl(1))**2
        endif

!!$$^%S^&A^^^**
!!    SHOULDN'T dpp be defined here, after any adjustnment to ph has 
!     BEEN MADE ???? Problem would be at topmost layer.
!!$$^%S^&A^^^**

!---------------------------------------------------------------------
!    if in a debug column, output values of cloud-level indices and 
!    cloud temperature, along with environmental mixing ratio and
!    temperature at these levels.
!---------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, 3i4)')  &
                             'in mulsub: kc,kcl,kch= ',kc,kcl,kch
          write (diag_unit, '(a, 3e20.14)')  &
                         'in mulsub: TCC = ',TCC(kc),TCC(kcl),TCC(kch)
          write (diag_unit, '(a, 3e20.12)')  &
                              'in mulsub: QE= ',mrE(kc),mrE(kcl),mrE(kch)
          write (diag_unit, '(a, 3e20.12)')   &
                              'in mulsub: TE= ',TE(kc),TE(kcl),TE(kch)
        endif

!----------------------------------------------------------------------
!    define the perturbation vertical subgrid scale temperature flux 
!    associated with the cloud (exf). note the special definitions at
!    cloud bottom and top.
!----------------------------------------------------------------------
        if (kc == 1) then
          exf = Param%rdgas*wv(kc)*dpdz_cb*(tcc(kc) - te(kc))*  &
                ((rcl(kc)/Param%cloud_base_radius)**2)/(Param%cp_air*p)
        else if (kc == kch) then
          exf = Param%rdgas*wv(kc)*dpdz(ncc_kou)*(tcc(kc) - te(kc))*  &
                ((rcl(kc)/Param%cloud_base_radius)**2)/(Param%cp_air*p)
        else             
          exf = Param%rdgas*wv(kc)*dpdz(kc)*(tcc(kc) - te(kc))*   &
                ((rcl(kc)/Param%cloud_base_radius)**2)/(Param%cp_air*p)
        endif

!----------------------------------------------------------------------
!    at the cloud top, define the cumulus vertical-flux convergence of 
!    temperature, moisture and tracers as the flux transport entering 
!    the layer which contains cloud top. this "through the top" trans-
!    port is assumed to be all deposited in this layer which includes 
!    cloud top. set the vertical eddy transports of temperature, moist-
!    ure and tracer at the uppermost interior cloud level to zero.
!----------------------------------------------------------------------
        if (kc == kch) then
          efchr(ncc_kou ) = ehf(ncc_kou)/Param%dp_of_cloud_model
          emfhr(ncc_kou ) = emf(ncc_kou)/Param%dp_of_cloud_model
          ehf(ncc_kou-1) = 0.
          emf(ncc_kou-1) = 0.

!##$%$%$%^&^&!!&&*** 
!!!???? is this correct ? shouldn't the flux across the upper interface
!!! of this layer be accounted for ?  I..e., 

!!!!     efchr(ncc-1) = exf +  (ehf(ncc-2) - ehf(ncc))/(2.*dpp)

!!   As it stands, it appears that all the flux into the ncc-1 layer 
!!   (ehf(ncc-2)) remains. however, some (ehf(ncc)) is exported 
!    into the ncc layer (above cloud top), suggesting a non-conservation
!!!   condition.
!##$%$%$%^&^&!!&&*** 

          efchr(ncc_kou-1) = exf + (ehf(ncc_kou-2))/(2.*dpp)
          emfhr(ncc_kou-1) = emf(ncc_kou-2)/(2.*dpp)
          do kcont=1,ntr
            etfhr(ncc_kou,kcont) = etf(ncc_kou,kcont)/Param%dp_of_cloud_model
            etf(ncc_kou-1,kcont) = 0.
            etfhr(ncc_kou-1,kcont) = etf(ncc_kou-2,kcont)/(2.*dpp)
          end do

!----------------------------------------------------------------------
!    add the appropriately pressure-weighted contributions from the
!    layer containing cloud top to the integrals of the in-cloud vert-
!    ical-flux convergence of entropy, (sumefc) potential temperature
!    (sumthet), moisture (sumemf) and tracers (sumetf).
!----------------------------------------------------------------------
          sumefc = sumefc + efchr(ncc_kou)*Param%dp_of_cloud_model/2.
          ptt = pt_kou + Param%dp_of_cloud_model
          sumthet = sumthet + efchr(ncc_kou)*    &
                                   ((1.0e05/ptt)**Param%kappa)* &
                                               Param%dp_of_cloud_model/2.
          sumemf = sumemf + emfhr(ncc_kou)*Param%dp_of_cloud_model/2.
          do kcont=1,ntr
            sumetf(kcont) = sumetf(kcont) + etfhr(ncc_kou,kcont)*  &
                            Param%dp_of_cloud_model/2.
          end do

!----------------------------------------------------------------------
!    if in diagnostic column, output the entropy and moisture flux 
!    convergence in the layer containing cloud top.
!----------------------------------------------------------------------
          if (debug_ijt) then 
            write (diag_unit, '(a, e20.12)')  &
                      'in mulsub: EFCHR(kc+1)= ',efchr(ncc_kou)
            write (diag_unit, '(a, e20.12)')  &
                          'in mulsub: EMFHRIT=kch ',emfhr(ncc_kou)
          endif

!---------------------------------------------------------------------
!    for all but the topmost layer, define the layer flux convergences
!    of entropy (efchr), moisture (emfhr) and tracers (etfhr). the ver-
!    tical flux convergence of entropy (efchr) is defined as the sum of
!    the in-cloud perturbation vertical subgrid scale temperature flux 
!    (exf) and the eddy flux convergence of temperature across 
!    the layer.
!---------------------------------------------------------------------
        else
          efchr(kc) = exf + (ehf(kcl) - ehf(kch))/(2.*dpp)
          emfhr(kc) = (emf(kcl) - emf(kch))/(2.*dpp)
          do kcont=1,ntr
            etfhr(kc,kcont) = (etf(kcl,kcont) - etf(kch,kcont))/(2.*dpp)
          end do
        endif

!----------------------------------------------------------------------
!    add the appropriately pressure-weighted contributions from the
!    current layer to the integrals of the in-cloud vertical-flux con-
!    vergence of entropy, (sumefc) potential temperature (sumthet), 
!    moisture (sumemf) and tracers (sumetf).
!----------------------------------------------------------------------
        sumefc  = sumefc  + (efchr(kc)*dpp)
        sumthet = sumthet + efchr(kc)*((1.0e05/p)**Param%kappa)*dpp
        sumemf  = sumemf  + (emfhr(kc)*dpp)
        do kcont=1,ntr
          sumetf(kcont) = sumetf(kcont) + etfhr(kc,kcont)*dpp
        end do

!---------------------------------------------------------------------
!    if in diagnostic column, write out the current vertical sum of 
!    the entropy and potential temperature flux convergence. output the
!    vertical flux of theta (thetf) and moisture (emf) at cloud base.
!    output the in-cloud (xclo) and environmental (xtrae) tracer 
!    mixing ratios.  if this is the 4th cumulus subensemble, output 
!    some fluxes, flux convergences, dpdz and other terms related to 
!    the flux calculation.
!---------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, 2e20.12)')  &
                    'in mulsub: SUMTHET,SUMEFC= ',sumthet, sumefc
          do kcont=1,ntr
            write (diag_unit, '(a, e20.12)')'in mulsub: SUMETF=', &
                                           sumetf(kcont)
          end do
          if (kc == 1) then
            thetf = ehf(kcl)*((1.0e05/pl)**Param%kappa)
            write (diag_unit, '(a, 2e20.12)')  &
                             'in mulsub: THETF,EMFF= ',thetf,emf(kcl)
          endif
          do kcont=1,ntr
            write (diag_unit, '(a, 3i4, 2e20.12)')  &
                   'in muls:kou,k,kcont,xclo,xtrae=' &
                    ,kou,kc,kcont,xclo(kc,kcont),xtrae(kc,kcont)
          end do
          if (kou == 4) then
            write (diag_unit, '(a, i4, 2f19.10)') &
                                  'in mulsub: kc,PL,PH= ',kc, pl, ph
            write (diag_unit, '(a, 4e20.12)')  &
                'in mulsub: EHFH,EHFL,EXF,efchr(kc)= ', &
                        ehf(kch), ehf(kcl), exf, efchr(kc)
            write (diag_unit, '(a, 3e20.12)')   &
                          'in mulsub: EMFH,EMFL,EMF= ',  &
                         emf(kch), emf(kcl), emfhr(kc)
            if (do_donner_tracer) then
              write (diag_unit, '(a, 3e20.12)')   &
                       'in mulsub: ETFH,ETFL,ETF=         ', &
                     etf(kch,ntr), etf(kcl,ntr), etfhr(kc,ntr)
            endif
            write (diag_unit, '(a, 3e20.12)')   &
                               'etfh diag: rcl,wv,dpdzh=         ', &
                                    rcl(kch), wv(kch), dpdz(kch)
            write (diag_unit, '(a, 3e20.12)')   &
                                'etfl diag: rcl,wv,dpdzl=         ', &
                                    rcl(kcl), wv(kcl), dpdz(kcl)
            do kcont=1,ntr
              write (diag_unit, '(a, 3e20.12)')  &
                                   'etfh diag: xclo,xtrae= ',         &
                                   xclo(kch,kcont), xtrae(kch,kcont)
              write (diag_unit, '(a, 3e20.12)')  &
                                 'etfl diag: xclo,xtrae= ',         &
                                   xclo(kcl,kcont), xtrae(kcl,kcont)
            end do
            write (diag_unit, '(a, 3e20.12)')   &
                            'in mulsub: WV,RH,QE= ',  &
                                     wv(kch), q_sat(kch), mre(kch)
            write (diag_unit, '(a, 3e20.12)')   &
                            'in mulsub: WV,RL,QE= ',   &
                                     wv(kcl), q_sat(kcl), mre(kcl)
          endif 
        endif
      end do  ! (end of kc loop)

!----------------------------------------------------------------------
!    if this is a debug column, output the moisture and temperature 
!    tendency profiles produced by this ensemble member.
!----------------------------------------------------------------------
      if (debug_ijt) then
        call don_cm_output_member_tends_k  &
             (nlev_hires, ncc_kou, diag_unit, Param, tcc, dpf,  &
              efchr, emfhr, cld_press, ermesg, error) 
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return
      endif

!------------------------------------------------------------------


end subroutine don_cm_compute_vert_fluxes_k



!#####################################################################

subroutine don_cm_output_member_tends_k    &
         (nlev_hires, ncc_kou, diag_unit, Param, tcc, dpf,   &
          efchr, emfhr, cld_press, ermesg, error) 

!----------------------------------------------------------------------
!    subroutine don_cm_output_member_tends_k prints out
!    the temperature (ctfhr) and moisture (cmfhr) tendency profiles 
!    produced by the current ensemble member.
!----------------------------------------------------------------------

use donner_types_mod, only : donner_param_type

implicit none 

!----------------------------------------------------------------------
integer,                     intent(in)  :: nlev_hires, ncc_kou,  &
                                            diag_unit
type(donner_param_type),     intent(in)  :: Param
real, dimension(nlev_hires), intent(in)  :: tcc, dpf, efchr, emfhr, &
                                            cld_press
character(len=*),            intent(out) :: ermesg
integer,                     intent(out) :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     ncc_kou        vertical index on cloud model grid of the level
!                    above cloud top (maximum number of vertical
!                    levels affected by current cloud)
!     diag_unit      unit number for column diagnostics output, if 
!                    diagnostics are requested for the current column
!     tcc            temperature field at cloud model levels
!                    [ degrees K ]
!     dpf            condensation rate profile 
!                    [ kg(h2o) / ( kg(air) sec) ]
!     efchr          profile of entropy tendency due to flux 
!                    convergence [ deg K / sec ]
!     emfhr          profile of moisture tendency due to flux 
!                    convergence [ kg(h2o) / (kg(air) sec) ]
!     cld_press      pressure at cloud model levels [ Pa ]
!     
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      real, dimension(nlev_hires)     :: cmfhr, ctfhr
      real                            :: convrat
      integer                         :: k

!---------------------------------------------------------------------
!   local variables:
!
!      cmfhr     moisture tendency due to flux convergence and 
!                condensation [ kg(h2o) / (kg(air) sec) ]
!      ctfhr     entropy tendency due to flux convergence and 
!                condensation [ deg K / sec ]
!      convrat   latent heat factor used with condensation term in 
!                entropy equation [ deg K ]
!      k         do-loop index
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    define the tendency terms at the levels within the cloud.
!---------------------------------------------------------------------
      do k=1,ncc_kou

!---------------------------------------------------------------------
!    define the appropriate latent heat to be used at the current 
!    cloud level. use the latent heat of vaporization for temperatures
!    at or above the freezing point (Param%tfre) and the latent heat of 
!    sublimation at temperatures lower than the freezing point. 
!---------------------------------------------------------------------
        if (tcc(k) >= Param%tfre) then
          convrat = Param%hlv/Param%cp_air   
        else
          convrat = Param%hls/Param%cp_air   
        endif

!---------------------------------------------------------------------
!    define the entropy tendency as the sum of the vertical flux conver-
!    gence of entropy (efchr) and the latent heat release.
!---------------------------------------------------------------------
        ctfhr(k) = -dpf(k)*convrat + efchr(k)

!---------------------------------------------------------------------
!     define the moisture tendency as the sum of vertical flux conver-
!     gence and total condensation. 
!---------------------------------------------------------------------
        cmfhr(k) = dpf(k) + emfhr(k)
      end do

!---------------------------------------------------------------------
!    set the values at layers above cloud top to 0.0.
!---------------------------------------------------------------------
      cmfhr(ncc_kou+1:nlev_hires) = 0.
      ctfhr(ncc_kou+1:nlev_hires) = 0.

!---------------------------------------------------------------------
!  ctfhr : cloud ensemble entropy tendency due to flux convergence and
!          condensation       
!          [ deg K / sec ]
!---------------------------------------------------------------------
      do k=1,nlev_hires       
        if (ctfhr(k) /= 0.0) then
          write (diag_unit, '(a, i4, 2e20.12)')  &
                    'in mulsub: k, P & cond/efc              =',    &
                     k, cld_press(k),ctfhr(k)
        endif
      end do

!---------------------------------------------------------------------
!  cmfhr : cloud ensemble moisture tendency due to flux convergence and
!          condensation       
!          [ kg(h2o) / (kg(air) sec) ]
!---------------------------------------------------------------------
      do k=1,nlev_hires       
        if (cmfhr(k) /= 0.0) then
          write (diag_unit, '(a, i4, 2e20.12)')  &
                   'in mulsub: k, P & cond/mfc              =',    &
                     k, cld_press(k),cmfhr(k)
        endif
      end do

!---------------------------------------------------------------------


end subroutine don_cm_output_member_tends_k


!#####################################################################


subroutine don_cm_process_condensate_k     &
       (nlev_lsm, nlev_hires, ntr, cldtop_indx, diag_unit, debug_ijt, &
        Param, acpre, accond, pb, pt_kou, pf, pftr, tcc, rcl, &
        cld_press, phalf_c, conint, dint, pmel, pmelt_lsm, precip, cu,&
        cell_precip, sumlhr, summel, dpf, dpftr, dfr, cell_melt, ermesg, error)

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

use donner_types_mod, only : donner_param_type

implicit none

!----------------------------------------------------------------------
integer,                      intent(in)    :: nlev_lsm, nlev_hires, &
                                               ntr, &
                                               cldtop_indx, diag_unit
logical,                      intent(in)    :: debug_ijt  
type(donner_param_type),      intent(in)    :: Param
real,                         intent(in)    :: acpre, accond, pb, pt_kou
real, dimension(nlev_hires),  intent(in)    :: pf, tcc, rcl, cld_press
real, dimension(nlev_hires,ntr), intent(in) :: pftr
real, dimension(nlev_lsm+1),  intent(in)    :: phalf_c
real,                         intent(inout) :: conint, dint,         &
                                               pmel, precip, cu, &
                                               cell_precip, sumlhr, &
                                               summel
real,                          intent(in) :: pmelt_lsm
real, dimension(nlev_hires),  intent(inout) :: dfr
real, dimension(nlev_hires),  intent(out)   :: dpf
real, dimension(nlev_hires,ntr),  intent(out)   :: dpftr
real, dimension(nlev_lsm),    intent(out)   :: cell_melt
character(len=*),             intent(out)   :: ermesg
integer,                      intent(out)   :: error
              

      real    :: dmela
      integer  :: k, kc, n
 
!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

      dpf(1:nlev_hires) = 0.
      dpftr(:,:) = 0.
      precip = 0.
      conint = 0.


      if (cldtop_indx /= 0) then
!--------------------------------------------------------------------
!    loop over the cloud levels to adjust various outputs as needed.
!--------------------------------------------------------------------
        do k=1,cldtop_indx

!--------------------------------------------------------------------
!    convert the layer-mean values of cloud-area-weighted condensation 
!    (pf) and wet deposition (pftr)
!    to values at cloud model interfaces (dpf and dpftr).
!--------------------------------------------------------------------
          if (k == 1) then
            dpf(1) =     pf(1)
            do n=1,ntr
              dpftr(1,n) =     pftr(1,n)
            end do
          else 
            dpf(k) = 0.5*(pf(k) + pf(k-1))
            do n=1,ntr
              dpftr(k,n) = 0.5*(pftr(k,n) + pftr(k-1,n))
            end do
          endif 
          if (k == cldtop_indx) then
            dpf(cldtop_indx+1) = 0.5*pf(cldtop_indx)
            do n=1,ntr
              dpftr(cldtop_indx+1,n) = 0.5*pftr(cldtop_indx,n)
            end do
          endif
        end do
      endif ! (cldtop_indx /= 0)

!--------------------------------------------------------------------
!    if in a diagnostics column, output the freezing (dfr) and conden-
!    sation (dpf) rates and cloud radius (rcl) at each cloud model
!    level, before normalization by the cloud base area.
!----------------------------------------------------------------------
      do k=1,cldtop_indx + 2
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 3e20.12)')  &
                'in mulsub: k,dfr,dpr,rcl= ', k, dfr(k) ,dpf(k), rcl(k)
        endif
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 3e20.12)')  &
                'in mulsub: k,dpf,pf,pfavg= ', k, dpf(k) ,pf(k),  &
! note pf(k-1) when k=1 is garbage !
                                  0.5*(pf(k) + pf(k-1))
        endif
 
!--------------------------------------------------------------------
!    normalize the freezing rate and condensation rate profiles by the
!    cloud base area of ensemble member #1.
!--------------------------------------------------------------------
        dfr(k) = dfr(k)/(Param%cloud_base_radius**2)
        dpf(k) = dpf(k)/(Param%cloud_base_radius**2)
        do n=1,ntr
          dpftr(k,n) = dpftr(k,n)/(Param%cloud_base_radius**2)
        end do
 
!----------------------------------------------------------------------
!    if in a diagnostics column, output the freezing (dfr) and conden-
!    sation (dpf) rates and cloud radius (rcl) at each cloud model
!    level, after normalization by the cloud base area.
!----------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 3e20.12)')  &
               'in mulsub: k,dfr,dpr,rcl= ', k, dfr(k), dpf(k), rcl(k)
          do n=1,ntr
            write (diag_unit, '(a, i4, e20.12)')  &
                'in mulsub: k,dpftr= ', k, dpftr(k,n)
          end do
        endif
      end do

      if (cldtop_indx /= 0) then
        do k=1,cldtop_indx
          conint = conint + pf(k)*Param%dp_of_cloud_model/Param%grav  
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 2e20.12)')  &
                          'in mulsub: kc, PF,coninT= ',k, pf(k)*  &
                                Param%dp_of_cloud_model/Param%grav, &
                                                            conint  
        endif
        end do

!--------------------------------------------------------------------
!    if in diagnostics column, output the ensemble member number (kou), 
!    moisture convergence integral (sub1sum),  and condensation 
!    integral (conint).
!--------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, e20.12, a)')  &
                   'in cloudm: CONINT= ',CONINT,' KG/(M**2)/SEC'
        endif

!--------------------------------------------------------------------
!    normalize the column integral of condensate by the cloud base area
!    of ensemble member #1. 
!--------------------------------------------------------------------
        conint = conint/(Param%cloud_base_radius**2)

!--------------------------------------------------------------------
!    define the precipitation generated by this ensemble member.
!--------------------------------------------------------------------
        precip = conint*acpre/accond
      endif ! (cldtop_indx /= 0)

!--------------------------------------------------------------------
!    define the condensation (cu) and precipitation rates (cell_precip) 
!    in units of mm(h2o) per day. add this ensemble member's contribution
!    to the total precipitation (ensmbl_precip) and condensation
!    (ensmbl_cond), normalized by the ensemble member's cloud base area.
!--------------------------------------------------------------------
      cu  = conint*Param%seconds_per_day
      cell_precip  = precip*Param%seconds_per_day
 
!----------------------------------------------------------------------
!    if in a diagnostics column, output the condensation (conint),
!    precipitataion (precip) and freezing (dint) rates for this ensemble
!    member, in both kg per m**2 per second and in mm per day.
!----------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 3e20.12, a)')  &
                   'in mulsub: CONPRE, CPRE,DINT= ', conint, precip, &
                      dint, ' KG/(M**2)/SEC'
        write (diag_unit, '(a, e20.12, a)')   &
                    'in mulsub: CONDENSATION PRE= ', cu, ' MM/DAY'
        write (diag_unit, '(a, e20.12, a)')  &
                 'in mulsub: CLOUD MODEL PRE= ', cell_precip, ' MM/DAY'
      endif

!---------------------------------------------------------------------
!    compute the contribution to the column integrals of total conden-
!    sation (sumlhr) and frozen condensate from the lowest cloud layer.
!    units are kg(h2o) per square meter per second.
!---------------------------------------------------------------------
      sumlhr = (dpf(1)*((Param%dp_of_cloud_model/2.)/Param%grav))
      if (tcc(1) <= Param%tfre) then
        summel = (Param%dp_of_cloud_model/2.)*dpf(1)/Param%grav  
      else
        summel = 0.
      endif

!---------------------------------------------------------------------
!    if in a diagnostics window, output the condensation rate and
!    accumulated condensate integral at level 1.
!---------------------------------------------------------------------
      if (debug_ijt) then
        kc = 1
        write (diag_unit, '(a, i4, 2e20.12)')  &
                               'in mulsub: kc,DPF, SUMLHR= ',kc, &
               dpf(1)*((Param%dp_of_cloud_model/2.)/Param%grav),   &
                                                sumlhr
      endif

!--------------------------------------------------------------------
!    add the contributions from the remaining cloud layers to the 
!    integrals.
!--------------------------------------------------------------------
      do kc = 2,cldtop_indx+1

!---------------------------------------------------------------------
!    add the density weighted increment of column condensate to the 
!    array accumulating it (sumlhr). if in diagnostic column, write
!    out the condensate at the level and the current vertical sum of 
!    the column condensate.
!---------------------------------------------------------------------
        sumlhr = sumlhr + (dpf(kc)*(Param%dp_of_cloud_model/Param%grav))
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 2e20.12)')  &
                         'in mulsub: kc,DPF,SUMLHR= ',kc,   &
               dpf(kc)*(Param%dp_of_cloud_model/Param%grav),   &
                                                               sumlhr
        endif

!---------------------------------------------------------------------
!    if the temperature is at or below the freezing level, add the 
!    current-level's density-weighted condensation to the array accum-
!    ulating the total frozen condensate in the column (summel).
!---------------------------------------------------------------------
        if (tcc(kc) <= Param%tfre) then
          summel = summel +  Param%dp_of_cloud_model*dpf(kc) /Param%grav 
        endif
      end do

!---------------------------------------------------------------------
!    if in diagnostics column, output the column integral of frozen 
!    condensate (summel). determine the amount of this condensate which
!    precipitates out by multiplying by the factor (cell_precip/cu).  add
!    this amount to the amount frozen in the updraft (dint) to obtain the
!    total amount which must be melted in the column (dints).
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12)')  'summel= ',summel
      endif

!--------------------------------------------------------------------
!    define the cloud model pressure at the level where melting occurs 
!    (pmel). 
!---------------------------------------------------------------------
      pmel = pb             
      do kc=1,cldtop_indx
        if ((tcc(kc) >= Param%kelvin) .and.   &
            (tcc(kc+1) <= Param%kelvin))  pmel = cld_press(kc)
      end do
      if (tcc(cldtop_indx+1) == Param%kelvin) pmel = pt_kou 
          if (debug_ijt) then
            write (diag_unit, '(a, 3e20.12 )')  &
             'pmelt_lsm, pmel from hi-res, pb   = ', pmelt_lsm, pmel, pb
          endif

!---------------------------------------------------------------------
!    if there has been no freezing in the column, then there will be
!    no melting.
!---------------------------------------------------------------------
      dmela = 0.
      if (dint == 0.) then
        cell_melt(1:nlev_lsm) = 0.
      else

!--------------------------------------------------------------------
!    if the melting level is above cloud base, partition the integrated
!    ice melt from the cloud model within the appropriate large-scale 
!    model layers.
!--------------------------------------------------------------------
        if (pb > pmelt_lsm) then

!--------------------------------------------------------------------
!    define the rate of ice melt (dmela) in units of g(h2o) per 
!    kg(air) per day. 
!--------------------------------------------------------------------
          if (cu /= 0.0) then
! melt ice precip plus frozen liq precip
            dmela = -((summel*cell_precip/cu + dint*cell_precip/cu)*  &
                                  Param%grav/(pmelt_lsm - pb))*8.64E07
         endif

!--------------------------------------------------------------------
!    call map_hi_res_intgl_to_lo_res_col to distribute this integrated 
!    column melting over the appropriate pressure interval on the 
!    large-scale model grid. if in a diagnostic column, output the 
!    melting rates as mapped to the large-scale grid.
!--------------------------------------------------------------------
          if (debug_ijt) then
            write (diag_unit, '(a, e20.12, 2f19.10)')  &
              'in cm_intgl_to_gcm_col: dmela,pb,pmelt_lsm= ',  &
                                                 dmela, pb, pmelt_lsm
          endif
          call don_u_map_hires_i_to_lores_c_k  &
            (nlev_lsm, dmela, pb, pmelt_lsm, phalf_c, cell_melt, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
          if (debug_ijt) then
            do k=1,nlev_lsm              
              if (cell_melt(k) /= 0.0) then
                write (diag_unit, '(a, i4, e20.12)') &
                 'in cm_intgl_to_gcm_col: k,cell_melt= ',k,cell_melt(k)
              endif
            end do
          endif

!---------------------------------------------------------------------
!    if the melting level is at or below cloud base, then there is no
!    in-cloud melting of frozen condensate.
!---------------------------------------------------------------------
        else
          cell_melt(1:nlev_lsm) = 0.
        endif  ! ( pb > pmel)
      endif  !  (dint == 0.0)


!---------------------------------------------------------------------
!    change the sign of the integral and convert to units of kg(h2o)
!    per square meter per day (or mm per day). if in diagnostics col-
!    umn, output the column integrated melting (summel), precipitation 
!    (cell_precip) and condensation rates (cu) in units of mm per day.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 3e20.12,a)')   &
             'in mulsub: summel,rc,cu= ',    &
        -(summel*cell_precip/cu)*Param%seconds_per_day,cell_precip,cu, &
                                            'mm/day'
      endif

!---------------------------------------------------------------------


end subroutine don_cm_process_condensate_k



!#####################################################################


subroutine don_cm_simult_k   &
         (diag_unit, debug_ijt, lfc_not_reached, Param, pcsave, rmu,  &
          cloud_temp, cloud_radius, w_vel, cloud_p, liq_wat, env_temp, &
          env_mixing_ratio, dtdp,  drdp, dwdp, ermesg, error)

!--------------------------------------------------------------------
!    subroutine simult returns the vertical derivatives of temperature
!    (dtdp), cloud radius (drdp) and vertical velocity (dwdp) in the
!    cloud.
!    Reference: Donner, JAS, 1986, v43, pp.2277-2288.
!    See LJD "Cloud Model 89" notes on Generalized mu (10/1/89)
!    and dwdz (10/2/89). The value of epsilon is taken as 1.
!    Version where cloud properties independent of cloud area.
!--------------------------------------------------------------------

use donner_types_mod, only : donner_param_type
use sat_vapor_pres_k_mod, only: compute_qs_k

implicit none

!--------------------------------------------------------------------
integer,                 intent(in)    ::  diag_unit
logical,                 intent(in)    ::  debug_ijt, lfc_not_reached
type(donner_param_type), intent(in)    ::  Param
real,                    intent(in)    ::  pcsave, rmu, cloud_temp,  &
                                           cloud_radius, w_vel,   &
                                           cloud_p, liq_wat, env_temp, &
                                           env_mixing_ratio
real,                    intent(out)   ::  dtdp, drdp, dwdp
character(len=*),        intent(out)   ::  ermesg
integer,                 intent(out)   ::  error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!        cloud_temp  cloud temperature [ deg K ]
!        cloud_radius 
!                    cloud radius  [ m ]
!        w_vel       cloud vertical velocity [ m / s ] 
!        cloud_p     pressure at current level [ Pa ]
!        env_temp    environmental temperature [ deg K ]
!        env_mixing_ratio
!                    environmental mixing ratio [ kg / kg ]
!        liq_wat     liquid water [ kg / kg ]
!        pcsave      pressure at which cloud ensemble member 1 becomes
!                    buoyant [ Pa ]
!        rmu         entrainment coefficient [ m**(-1) ]
!        lfc_not_reached   if true, have not yet reached buoyant part of 
!                    cloud (vertical velocity has not yet begun to
!                    increase)
!        debug_ijt   is this a diagnostics column ?
!        diag_unit   output unit number for this diagnostics column
!
!   intent(out) variables:
!
!        dtdp        temperature derivative [ deg K / Pa ]
!        drdp        cloud radius derivative [ m / Pa ]
!        dwdp        vertical velocity derivative [ m / (sec Pa) ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      real  ::  areap, c2, c3, c4, dadp, dw2dz, dwdz, dzdp, dz, es,  &
                entrain_co, htv, htve, lat, rhoaw, rhoawp, &
                rstar, sphum, tcae, teae, test, west, w2test, wtest
      integer :: nbad

!---------------------------------------------------------------------
!   local variables:
!
!         areap             updated cloud area [ m**2 ]
!         c2                term in dt/dz eqn (eqn 5, Donner, JAS, 1993)
!         c3                term in dt/dz eqn (eqn 5, Donner, JAS, 1993)
!         c4                term in dt/dz eqn (eqn 5, Donner, JAS, 1993)
!         dadp              derivative of cloud area with respect to
!                           pressure [ m**2 / Pa ]
!         dw2dz             vertical derivative of vertical velocity
!                           squared with respect to height 
!                           [ m / sec**2 ]
!         dwdz              vertical derivative of vertical velocity
!                           with respect to height [ sec**(-1) ]
!         dzdp              dz/dp when virtual mass coefficient is
!                           used [ m / Pa ]
!         dz                delta z corresponding to dp_of_cloud_model
!                           [ m ]
!         es                saturation vapor pressure [ Pa ]
!         entrain_co        entrainment coefficient [ m**(-1) ]
!         htv               virtual temperature of cloud       [ deg K ]
!         htve              virtual temperature of environment [ deg K ]
!         lat               latent heat relevant at input temperature
!                           cloud_temp [ J/kg ]
!         rhoaw             rho * cloudarea * vertical velocity for
!                           input values [ kg / sec ]
!         rhoawp            rho * cloudarea * vertical velocity for
!                           updated values [ kg / sec ]
!         rstar             gas constant for moist air [ J / (kg degK) ]
!         sphum             saturation specific humidity [ kg / kg ]   
!         tcae              equivalent potential temperature in cloud
!                           [ deg K ]
!         teae              equivalent potential temperature in env-
!                           ironment [ deg K ]
!         test              updated temperature [ deg K ]
!         west              updated vetical velocity [ m / sec ]
!         w2test            updated vertical velocity squared (without
!                           entrainment term) [ m**2 / sec**2 ]
!         wtest             updated vertical velocity  (without
!                           entrainment term) [ m / sec ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!   define appropriate latent heat for the given cloud temperature.
!---------------------------------------------------------------------
      if (cloud_temp < Param%tfre) then
        lat = Param%hls    
      else
        lat = Param%hlv   
      endif
 
!---------------------------------------------------------------------
!    define the specific humidity at the cloud temperature (sphum) and
!    the gas constant for moist air (rstar).
!---------------------------------------------------------------------
      call compute_qs_k (cloud_temp, cloud_p, Param%D622,  &
                           Param%D608, sphum, nbad, esat = es)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (nbad /= 0) then
        ermesg = 'subroutine don_cm_simult_k: '// &
                 'temperatures out of range of esat table'
        error = 1
        return
      endif

      rstar = Param%rdgas*(1. + Param%d608*sphum)

!---------------------------------------------------------------------
!    define the in-cloud (htv) and environmental (htve) virtual 
!    temperature.                                                       
!---------------------------------------------------------------------
      htv  = cloud_temp*(1. + Param%d608*sphum)
      htve = env_temp*(1. + Param%d608*(env_mixing_ratio/   &
                                                (1.0+env_mixing_ratio)))

!---------------------------------------------------------------------
!    if in diagnostics window, output the in-cloud and environmmental
!    virtual temperatures.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)')  &
                      'in simult: htve,htv= ',htve, htv
      endif

!---------------------------------------------------------------------
!    define one of the terms present in the dt/dz eqn (Eqn 5) in 
!    Donner, JAS, 1993, v50, pp.890-906.
!---------------------------------------------------------------------
      c2 = (htv + (lat*sphum/rstar))*Param%grav/(Param%cp_air*cloud_temp)

!---------------------------------------------------------------------
!    if in diagnostics window, output the c2 term in the dt/dz equation
!    and the environmmental temperature, moisture and pressure at this
!    level. 
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12)')  'in simult: c2= ',c2
        write (diag_unit, '(a, f20.14, f19.10, e20.12)')  &
              'in simult: te,p,qe= ',env_temp, cloud_p ,env_mixing_ratio
      endif

!--------------------------------------------------------------------
!    call tae to calculate environmental adiabatic equivalent 
!    temperature (teae).
!--------------------------------------------------------------------
      call don_cm_tae_k    &
           (Param, env_temp, cloud_p, env_mixing_ratio, lat, teae,  &
            ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    if in diagnostics window, output the environmental adiabatic 
!    equivalent temperature.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, f20.14)') 'in simult: teae= ',teae
      endif

!--------------------------------------------------------------------
!    calculate in-cloud adiabatic equivalent temperature (tcae).
!--------------------------------------------------------------------
      tcae = cloud_temp*exp(lat*sphum/(Param%cp_air*cloud_temp))

!----------------------------------------------------------------------
!    define dz/dp in the presence of a virtual mass coefficient (used
!    to roughly account for non-hydrostatic effects). define the cor-
!    resonding dz for the cloud model pressure increment dp.
!    Reference: Donner, JAS, 1986, v43, pp.2277-2288.
!----------------------------------------------------------------------
      dzdp = -Param%rdgas*(1. + Param%virt_mass_co)*htv*htve/    &
                                (Param%grav*cloud_p*  &
                                      (Param%virt_mass_co*htve + htv))
      dz   = Param%dp_of_cloud_model*dzdp

!---------------------------------------------------------------------
!    define the remaining terms (c3, c4) in the dt/dz equation (eqn (5)
!    in the reference cited below).
!    note : the entrainment coefficient (dln(rho*a*w)/dz) is given by 
!    (exp(rmu*dz) - 1.)/(exp(rmu*dz)*dz), with assumption that 
!    rhoaw(z) = rhoaw(0)*exp(rmu*dz)
!    Reference: Donner, JAS, 1993, v50, pp.890-906.
!---------------------------------------------------------------------
      entrain_co = (exp(rmu*dz) - 1.0)/(exp(rmu*dz)*dz)
      c3 = (tcae - teae)/exp(lat*sphum/(Param%cp_air*cloud_temp))*  &
                                                               entrain_co
      c4 = 1. + (Param%d622*lat*(Param%hlv*es/     &
                   (Param%rvgas*(cloud_temp**2)))/(Param%cp_air*cloud_p))

!---------------------------------------------------------------------
!    if in diagnostics window, output the c3 and c4 terms in the dt/dz 
!    equation.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)')  'in simult: c3,c4= ',c3,c4
      endif

!---------------------------------------------------------------------
!    define the dt/dp tendency by combining the component terms 
!    in the dt/dz equation and multiplying by dz/dp. define an 
!    updated estimated value of in-cloud temperature (test).
!---------------------------------------------------------------------
      dtdp = -((c2 + c3)/c4)*dzdp
      test = cloud_temp + dtdp*Param%dp_of_cloud_model

!--------------------------------------------------------------------
!    define the buoyancy and drag terms in the d(w**2)/dz equation. 
!    Reference: Donner, JAS, 1993, v50, pp.890-906.
!--------------------------------------------------------------------
      dw2dz = 2.0*((Param%grav*(htv - htve)/      &
                                  (htve*(1. + Param%virt_mass_co))) -  &
                                                      Param%grav*liq_wat)

!-------------------------------------------------------------------
!    produce an updated vertical velocity squared (wtest) by adding this
!    tendency term to the input value. be sure w**2 is positive; then 
!    define the vertical velocity.
!-------------------------------------------------------------------
      w2test = (w_vel**2) + dw2dz*dz
      if (w2test < 0.) then
        wtest = 0.
      else
        wtest = sqrt(w2test)
      endif

!---------------------------------------------------------------------
!    if in diagnostics window, output the updated vertical velocity.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, l4, e20.12)')  &
                     'in simult: testlc,wtest= ',lfc_not_reached, wtest
      endif

!----------------------------------------------------------------------
!    define the dw/dp derivative from this updated value, ignoring the 
!    entrainment term.
!----------------------------------------------------------------------
      dwdp = (wtest - w_vel)/Param%dp_of_cloud_model

!----------------------------------------------------------------------
!    determine if the necessary conditions for entrainment are present.
!    allow entrainment to change vertical velocity if parcel has 
!       (1) previously achieved initial acceleration 
!           (.not. lfc_not_reached) 
!           and it has ascended above the level where the most entrain-
!           ing parcel (ensemble member #1) initially accelerates or,
!       (2) the parcel is currently accelerating. 
!    the first condition ensures that a shallow, lower-entrainment cloud
!    will not develop below the pressure where the most entraining
!    parcel initially accelerates.
!----------------------------------------------------------------------
      if ( ((.not. lfc_not_reached) .and. (cloud_p <= pcsave)) .or.   &
           (dwdp <= 0.) ) then
        dwdz = -w_vel*entrain_co
        dwdp = (dwdz*dzdp) + dwdp
      endif

!----------------------------------------------------------------------
!    if the parcel has not reached the level of free convection, do
!    not allow its vertical velocity to decrease. BL turbulence or 
!    other sub-grid mechanisms are assumed present to sustain its
!    upward motion.
!----------------------------------------------------------------------
      if ((lfc_not_reached) .and. (dwdp > 0.) )  then
        dwdp = 0.
      endif

!----------------------------------------------------------------------
!    if the parcel is above its level of free convection  but below
!    the pressure level at which the ensemble member #1 initially 
!    accelerates, do not allow its vertical velocity to decrease. this
!    ensures that less entraining clouds will not develop between the
!    ground and the pressure at which the most entraining parcel 
!    initally accelerates. BL turbulence or some other sub-grid 
!    mechanism is assumed to sustain its upward motion.
!----------------------------------------------------------------------
      if ( (.not. lfc_not_reached) .and. (cloud_p > pcsave) .and.    &
           (dwdp > 0.) )  then
        dwdp = 0.
      endif

!---------------------------------------------------------------------
!    if in diagnostics window, output the updated vertical velocity.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a,l4,  2e20.12)')  &
          'in simult: testlc,dwdp,test= ',lfc_not_reached,dwdp,test
      endif

!---------------------------------------------------------------------
!    produce an updated vertical velocity (eqn (6) of reference).
!---------------------------------------------------------------------
      west = w_vel + dwdp*Param%dp_of_cloud_model

!--------------------------------------------------------------------
!    calculate the terms in the cloud area tendency equation.
!    Reference: Donner, JAS, 1993, v50, pp.890-906.
!--------------------------------------------------------------------
      if (west < Param%wdet) then

!--------------------------------------------------------------------
!    if the updated vertical velocity is weaker than the detrainment
!    velocity, set the cloud radius tendency to 0.0,maintaining the
!    input value.
!--------------------------------------------------------------------
        drdp = 0.
      else

!--------------------------------------------------------------------
!    define the value of rho*cloudarea*vertical velocity (rhoaw) for 
!    the input values of these variables. define the same quantity at 
!    a level dz higher (rhoawp), assuming rhoaw(z) = 
!    rhoaw(0)*exp(rmu*dz).
!--------------------------------------------------------------------
        rhoaw = cloud_p*(cloud_radius**2)*w_vel/(Param%rdgas*htv)
        rhoawp = rhoaw*exp((dzdp*rmu)*Param%dp_of_cloud_model)

!---------------------------------------------------------------------
!    if in diagnostics window, output rho*cloudarea*w obtained from
!    input values and at a height dz higher.
!---------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, f19.10, 2e20.12)')  &
                       'in simult: p,fm,fmp= ', cloud_p, rhoaw, rhoawp
        endif

!---------------------------------------------------------------------
!    define the virtual temperature (htv) corresponding to the updated 
!    temperature test.
!----------------------------------------------------------------------
        call compute_qs_k (test, cloud_p, Param%D622,  &
                           Param%D608, sphum, nbad)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (nbad /= 0) then
          ermesg = 'subroutine don_cm_simult_k: '// &
                   'temperatures out of range of esat table'
          error = 1
          return
        endif
        htv = test*(1. + Param%d608*sphum)

!----------------------------------------------------------------------
!    define the updated cloud area (areap). define the area tendency 
!    (dadp) and then the cloud radius tendency (drdp).
!----------------------------------------------------------------------
        areap = rhoawp*Param%rdgas*htv/  &
                               ((cloud_p + Param%dp_of_cloud_model)*west)
        dadp = (areap - (cloud_radius**2))/Param%dp_of_cloud_model
        drdp = dadp/(2.*cloud_radius)
      endif

!--------------------------------------------------------------------



end subroutine don_cm_simult_k




!####################################################################

subroutine don_cm_tae_k    &
         (Param, init_temp, init_pr, parcel_mixing_ratio, latent_heat,  &
          equivalent_temp, ermesg, error)

!--------------------------------------------------------------------
!    subroutine tae determines the saturation temperature of an init-
!    ially non-saturated parcel ( defined by init_temp, 
!    parcel_mixing_ratio, init_pr) by incrementally moving the parcel 
!    upward dry adiabatically until it reaches the temperature at which 
!    saturation would occur for the initial pressure init_pr.  using this
!    saturation temperature (te) and the parcel moisture 
!    parcel_mixing_ratio and pressure init_pr, the equivalent temperature
!    is calculated (the temperature attained if all vapor were condensed 
!    out by moist adiabatic ascent, followed by dry adiabatic descent to
!    the initial pressure level init_pr.
!--------------------------------------------------------------------

use donner_types_mod, only : donner_param_type
use sat_vapor_pres_k_mod, only: compute_mrs_k

implicit none

!--------------------------------------------------------------------
type(donner_param_type), intent(in)    :: Param
real,                    intent(in)    :: init_temp, init_pr,   &
                                          parcel_mixing_ratio,  &
                                          latent_heat
real,                    intent(out)   :: equivalent_temp
character(len=*),        intent(out)   :: ermesg
integer,                 intent(out)   :: error

!--------------------------------------------------------------------
!   intent(in) variables:
!
!        init_temp              temperature   [ deg K ]
!        init_pr                pressure      [ Pa ]
!        parcel_mixing_ratio    mixing ratio [kg(H2O)/kg (dry air) ]
!        latent_heat            applicable latent heat constant  [ J/kg ]
!
!   intent(out) variables:
!
!        equivalent_temp        adiabatic equivalent temperature 
!                               [ deg K ]
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      real      ::   pr, te, mre
      integer   ::   nbad

!--------------------------------------------------------------------
!   local variables:
!
!        pr   pressure at current level  [ Pa ]
!        te   temperature of parcel move dry adiabatically from
!             pressure p to pressure pr
!        es   stauration vapor pressure at temperature te
!        mre  saturation mixing ratio at temperature te and pressure p
!        k    do-loop index
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutinE.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    initialize the output variable equivalent_temp and the pressure at 
!    current level (pr).
!---------------------------------------------------------------------
      equivalent_temp = init_temp
      pr = init_pr

!---------------------------------------------------------------------
!    vertical loop. if the pressure  pr is less than pstop (absolute 
!    lowest pressure for cloud) exit the loop.
!---------------------------------------------------------------------
      do while (pr >= Param%pstop) 

!--------------------------------------------------------------------
!    define the temperature (te) at pressure pr assuming dry adiabatic
!    ascent from initial pressure p. determine the saturation vapor
!    pressure (es) for this temperature. define the saturation mixing
!    ratio for a parcel of this temperature at the initial pressure
!    level (mre).
!--------------------------------------------------------------------
        te = init_temp*((pr/init_pr)**Param%kappa)
        call compute_mrs_k (te, init_pr,                               &
                            Param%d622 , Param%d608 , mre, nbad, &
                            mr = parcel_mixing_ratio)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (nbad /= 0) then
          ermesg = 'subroutine don_cm_tae_k: '// &
                   'temperatures out of range of esat table'
          error = 1
          return
        endif

!--------------------------------------------------------------------
!    determine if saturation at the initial pressure level would occur
!    for this temperature.
!--------------------------------------------------------------------
        if (parcel_mixing_ratio >= mre) then

!--------------------------------------------------------------------
!    if saturation would occur for this temperature (te) at the initial
!    pressure init_pr (i.e., parcel_mixing_ratio >= mre), then use te to
!    calculate the equivalent temperature. exit the loop. otherwise,
!    increment the current pressure and continue within the loop.
!--------------------------------------------------------------------
          equivalent_temp = init_temp*exp(latent_heat*  &
                                   parcel_mixing_ratio/(te*Param%cp_air))
          exit
        endif
        pr = pr + Param%dp_of_cloud_model
      end do

!---------------------------------------------------------------------



end subroutine don_cm_tae_k

!####################################################################

subroutine don_cm_micro_k   &
         (diag_unit, debug_ijt, Param, tc1, tc2, p1, p2, te1, te2,  &
          qe1, qe2, w1, w2, rr, rmu, qrw, qcw, qlw, dcw1, dqrw3, ermesg, error)

!----------------------------------------------------------------------
!    subroutine micro calculates microphysical tendencies (kessler
!    microphysics) of the cloud parcel defined by (tc, w, qrw, qcw, 
!    qlw, rr) in an environment defined by (te, qe, rmu) during the 
!    movement of the parcel from pressure level p1 to level p2.
!--------------------------------------------------------------------

use donner_types_mod, only : donner_param_type
use sat_vapor_pres_k_mod, only: compute_mrs_k

implicit none

!--------------------------------------------------------------------
integer,                  intent(in)    :: diag_unit
logical,                  intent(in)    :: debug_ijt
type(donner_param_type),  intent(in)    :: Param
real,                     intent(in)    :: tc1, tc2, p1, p2, te1, te2, &
                                           qe1, qe2, w1, w2, rr, rmu
real,                     intent(inout) :: qrw, qcw, qlw, dcw1, dqrw3
character(len=*),         intent(out)   :: ermesg
integer,                  intent(out)   :: error

!--------------------------------------------------------------------
!  intent(in) variables:
!
!        tc1       cloud temperature at starting point [ deg K ]
!        tc2       cloud temperature at ending point [ deg K ]
!        p1        pressure at starting point [ Pa ]
!        p2        pressure at ending point [ Pa ]
!        te1       environmental temperature at starting point [ deg K ]
!        te2       environmental temperature at ending point [ deg K ]
!        qe1       environmental mixing ratio at starting point  
!                  [ kg(H2O) / kg(air) ] 
!        qe2       environmental mixing ratio at ending point  
!                  [ kg(H2O) / kg(air) ] 
!        w1        cloud vertical velocity at starting point [ m / sec ]
!        w2        cloud vertical velocity at ending point [ m / sec ]
!        rr        cloud radius [ m ]
!        rmu       entrainment coefficient [ m**(-1) ]
!        debug_ijt is this a diagnostics column ?
!        diag_unit output unit number for this diagnostics column
!
!   intent(inout) variables:
!
!        qrw    rain water              [ g(h2o) / m**3 ]
!        qcw    cloud water             [ g(h2o) / m**3 ]
!        qlw    total liquid water      [ kg(h2o) / kg(air) ]
!        dcw1   condensation increment  [ g(h2o) / m**3 ]
!        dqrw3  precipitation increment [ g(h2o) / m**3 ]
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:
  
      real :: cond, d1, d2, dcw2, dqcw3, dt_micro, dz, ent, pav,  &
              rb, qcwa, qeb, qrwa, red, rho, rs1, rs2, tcb, wav
      integer :: nbad

!--------------------------------------------------------------------
!   local variables:
!
!     cond      condensation in the layer [ kg(h2o) / kg(air) ]
!     d1        dz/dp at the starting level [ m / Pa ]  
!     d2        dz/dp at the ending level [ m / Pa ]  
!     dcw2      change in cloudwater due to autoconversion 
!               [ g (h2o) / m**3 ]
!     dqcw3     change in cloudwater content due to accretion     
!               [ g (h2o) / m**3 ]
!     dt_micro  microphysical time step [ sec ]
!     dz        average delta z in the layer [ m ]
!     ent       ratio of parcel mass after entrainment to mass
!               before entrainment [ dimensionless ]
!     es1       saturation vapor pressure at level p1 [ Pa ]
!     es2       saturation vapor pressure at level p2 [ Pa ]
!     pav       average pressure in the layer [ Pa ]
!     rb        average in-cloud mixing ratio in the layer 
!               [ kg(h2o) / kg(air) ]
!     qcwa      initial value of cloudwater after microphysical term
!               updates; it is then adjusted to avoid negative values
!               [ g(h2o) / m**3 ]
!     qeb       average environmental mixing ratio in the layer 
!               [ kg(h2o) / kg(air) ]
!     qrwa      initial value of rainwater after microphysical term
!               updates; it is then adjusted to avoid negative values
!               [ g(h2o) / m**3 ]
!     red       reducing factor applied to the microphysical loss terms
!               in the cloud and rain equations to avoid creation of
!               negative cloud and rain water [ dimensionless ]
!     rho       atmospheric density [ kg / m**3 ]
!     rs1       saturation mixing ratio in cloud at level p1
!     rs2       saturation mixing ratio in cloud at level p2
!     tcb       mean in-cloud temperature in the layer [ deg K ]
!     teb       mean environmental temperature in the layer [ deg K ]
!     wav       mean vertical velocity in the layer [ m / sec ]
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!--------------------------------------------------------------------
!    define the saturation mixing ratio (rs1) at level p1.
!--------------------------------------------------------------------
      call compute_mrs_k (tc1, p1, Param%d622 , Param%d608 , rs1, nbad)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (nbad /= 0) then
        ermesg = 'subroutine don_cm_micro_k: '// &
                 'temperatures out of range of esat table'
        error = 1
        return
      endif

!--------------------------------------------------------------------
!    define the saturation mixing ratio (rs2) at level p2.
!--------------------------------------------------------------------
      call compute_mrs_k (tc2, p2, Param%d622 , Param%d608 , rs2, nbad)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (nbad /= 0) then
        ermesg = 'subroutine don_cm_micro_k: '// &
                 'temperatures out of range of esat table'
        error = 1
        return
      endif

!--------------------------------------------------------------------
!    if in diagnostic column, output the relevant atmospheric fields 
!    at the two pressure levels.
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)')   &
                         'in micro: qrw,qcw= ',qrw,qcw
        write (diag_unit, '(a, e20.12)')  'in micro: rr= ',rr
        write (diag_unit, '(a, 2e20.12, 2f19.10)') &
                         'in micro: rs1,rs2,p1,p2= ',rs1,rs2,p1,p2
        write (diag_unit, '(a, 2e20.12, 2f20.14)')  &
!                        'in micro: es1,es2,tc1,tc2= ',es1,es2,tc1,tc2
                         'in micro: tc1,tc2= ',tc1,tc2
        write (diag_unit, '(a, 2f20.14)')  &
                          'in micro: te1,te2= ',te1,te2
      endif

!--------------------------------------------------------------------
!    define the layer-mean cloud temperature (tcb) and mixing ratio 
!    (rb), environmental mixing ratio (qeb), vertical velocity (wav), 
!    pressure (pav) and density (rho). define a layer-mean dz, as the  
!    average of the values of dz/dp (when using a virtual mass coeffic-
!    ient) at the two input levels. define the ratio of post-entrainment
!    mass to the pre-entrainment mass (ent).
!--------------------------------------------------------------------
      tcb = 0.5*(tc1 + tc2)
      rb  = 0.5*(rs1 + rs2)
      qeb = 0.5*(qe1 + qe2)
      wav = 0.5*(w1  +  w2)
      pav = 0.5*(p1  +  p2)
      rho = pav/(Param%rdgas*tcb*(1. + Param%d608*(rb/(1.0+rb))))
      d1  = Param%rdgas*(1. + Param%virt_mass_co)*tc1*te1/  &
            (Param%grav*p1*(Param%virt_mass_co*te1 + tc1))
      d2  = Param%rdgas*(1. + Param%virt_mass_co)*tc2*te2/  &
            (Param%grav *p2*(Param%virt_mass_co*te2 + tc2))
      dz  = -0.5*(d1 + d2)*Param%dp_of_cloud_model
      ent = exp(rmu*dz)

!--------------------------------------------------------------------
!    if in diagnostic column, output the entrainment coefficient (rmu)
!    and the environmental mixing ratio (qeb).
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)')  &
                   'in micro: qeb,rmu= ', qeb, rmu
      endif

!----------------------------------------------------------------------
!    calculate the condensation (cond). the condensation is made up of
!    that which would occur due to adiabatic expansion (rs2 - rs1)
!    less the vapor which must be evaporated to saturate the environ-
!    mental air entrained into the cloud ((ent-1.)*(rs2-qeb)).  cond is
!    further modified to allow the condensation to be added to the
!    parcel at any point during the timestep through the variable
!    tr_insert_time. if the condensate is assumed to be added to the 
!    parcel while at mass m (its initial size), tr_insert_time = 0; 
!    if assumed to be added at end of step when parcel mass is ent 
!    times larger than its initial value, tr_insert_time is 1.  force 
!    the condensation to be non-negative.
!----------------------------------------------------------------------
      cond = (rs2 - rs1) + (ent - 1.)*(rs2 - qeb)
      cond = -cond/(1. - Param%tr_insert_time +    &
                         Param%tr_insert_time*ent)
      cond  = max (cond, 0.0)

!--------------------------------------------------------------------
!    if in diagnostic column, output the condensation (cond).
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12)') 'in micro: cond= ', cond
      endif

!--------------------------------------------------------------------
!    convert the condensation from units of kg(h2o)/ kg(air) to units
!    of (g(h2o) / m**3) by multiplying by rho*1000. define the timestep
!    for microphysical processes (dt_micro).
!--------------------------------------------------------------------
      dcw1 = cond*rho*1000.
      dt_micro = dz/wav

!--------------------------------------------------------------------
!    calculate the amount of cloud autoconverted to rain (dcw2). 
!--------------------------------------------------------------------
      if (qcw >= Param%autoconv_threshold)  then
        dcw2 = Param%autoconv_rate*(qcw - Param%autoconv_threshold)* &
                                                                dt_micro 
      else
        dcw2 = 0.0
      endif

!--------------------------------------------------------------------
!    calculate the cloud accretion by rainwater (dqcw3).
!--------------------------------------------------------------------
      if (qcw /= 0.0 .and. qrw /= 0.0) then
        dqcw3 = 5.26e-03*qcw*(qrw**.875)*dt_micro 
      else
        dqcw3 = 0.
      endif

!--------------------------------------------------------------------
!    calculate effect of entrainment on cloud water.
!--------------------------------------------------------------------
      qcw = qcw/ent

!---------------------------------------------------------------------
!    add the microphysical terms to the cloud water equation. if neces-
!    sary adjust the magnitudes of the loss terms so that negative
!    values of cloud water are not created. define the final value of
!    qcw.
!---------------------------------------------------------------------
      qcwa = qcw + dcw1 - dcw2 - dqcw3
      if (qcwa < 0.) then
        red = (qcw + dcw1)/(dcw2 + dqcw3)
        dcw2 = dcw2*red
        dqcw3 = dqcw3*red
      endif
      qcw = qcw + dcw1 - dcw2 - dqcw3

!--------------------------------------------------------------------
!    if in diagnostic column, output the cloud water entrainment  
!    factor (ent), the updated cloud water (qcw) and the condensation 
!    term (dcw1).
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 3e20.12)')  &
                             'in micro: ent,qcw,dcw1= ',ent,qcw,dcw1
      endif

!--------------------------------------------------------------------
!    calculate the flux of rainwater due to fallout.
!--------------------------------------------------------------------
      if (qrw /= 0.0) then
        dqrw3 = (qrw**1.125)*5.1*dt_micro/rr
      else
        dqrw3 = 0.  
      endif
 
!--------------------------------------------------------------------
!     calculate effect of entrainment on rain water
!--------------------------------------------------------------------
      qrw = qrw/ent

!---------------------------------------------------------------------
!    add the microphysical terms to the rain water equation. if neces-
!    sary adjust the magnitudes of the loss term so that negative
!    values of rain water are not created. define the final value of
!    qrw.
!---------------------------------------------------------------------
      qrwa = qrw + dcw2 + dqcw3 - dqrw3
      if (qrwa < 0.) then
        red   = (qrw + dcw2 + dqcw3)/dqrw3
        dqrw3 = red*dqrw3
      endif
      qrw = qrw + dcw2 + dqcw3 - dqrw3

!--------------------------------------------------------------------
!    if in diagnostic column, output the rainwater entrainment effect
!    (ent) and the updated value of rainwater. (qrw).
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)') &
                               'in micro: ent,qrw= ',ent,qrw
      endif

!--------------------------------------------------------------------
!    apply realizability consitions to qcw and qrw. define total liquid
!    water qlw.
!--------------------------------------------------------------------
      qcw = max (qcw, 0.0)
      qrw = max (qrw, 0.0)
      qlw = qrw + qcw

!--------------------------------------------------------------------
!    if in diagnostic column, output the condensation (dcw1), the 
!    rainwater fallout (dqrw3) and the total liquid remaining (qlw).
!--------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 3e20.12)')   &
                      'in micro: exit micro dcw1,dqrw3,qlw= ',   &
                       dcw1,dqrw3,qlw
      endif

!--------------------------------------------------------------------
!    convert the liquid water from g / m**3 to kg / kg.
!--------------------------------------------------------------------
      qlw = 1.0E-03*qlw/rho

!-------------------------------------------------------------------




end subroutine don_cm_micro_k


!######################################################################

subroutine don_cm_freeze_liquid_k    &
         (k, diag_unit, debug_ijt, Param, tbot, ttop, qlwest, dfrac, &
          dtfr, dtupa, dfrtop, sumfrea, ermesg, error)

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

use donner_types_mod, only : donner_param_type

implicit none

!----------------------------------------------------------------------
integer,                  intent(in)    :: k, diag_unit
logical,                  intent(in)    :: debug_ijt 
type(donner_param_type),  intent(in)    :: Param
real,                     intent(in)    :: tbot, ttop, qlwest
real,                     intent(inout) :: dfrac, dtfr, dtupa, sumfrea
real,                     intent(inout) :: dfrtop
character(len=*),         intent(out)   :: ermesg
integer,                  intent(out)   :: error


      real ::  dfraca, dtupb

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!--------------------------------------------------------------------
!    define the amount of liquid in the updraft which may be frozen.
!    multiply by a factor (freeze_fraction) to take account that not all
!    this water will freeze before falling out. use Leary and Houze
!    (JAS,1980) cell_precip/cu ratio to estimate this ratio.
!--------------------------------------------------------------------
      if ((tbot >= Param%tfre) .and. (ttop <= Param%tfre) .and.    &
          (dtfr == 0.)) then
        dtfr = qlwest*Param%hlf/Param%cp_air
        dtfr = Param%freeze_fraction*dtfr
      endif

!---------------------------------------------------------------------
!    define the fraction of liquid which is to be frozen at this level. 
!    freezing is assumed to occur linearly between 258 K and 248 K 
!    (dfre), with all liquid being frozen at 248 K. the fraction is not
!    allowed to decrease from its previous value if the temperature 
!    starts to increase with height within this temperature range. the 
!    fraction of the total which is frozen at the current level is 
!    dtupb; the amount frozen at this level (dfr(k+1) is the difference
!    between this amount and the amount frozen at the previous level 
!    (dtupa). cumulative total of frozen liquid is kept in sumfrea.
!---------------------------------------------------------------------
      if (dtfr > 0.0 .and. dtfr /= dtupa) then 
        dfraca = MIN ((Param%tfre - ttop)/Param%dfre, 1.0)
        dfrac = AMAX1 (dfrac, dfraca)
        dtupb = dtfr*dfrac 
        dfrtop   = dtupb - dtupa
        sumfrea = sumfrea + dfrtop   
        dtupa = dtupb

!--------------------------------------------------------------------
!    if in diagnostics column, output the values of available liquid
!    from freezing (dtfr), the amount frozen at level k+1 (dfr) and the 
!    constants hlf    and cp_air at this level.
!--------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(4(a, e20.12),a, i4)') &
             'in cloudm: DTFR=',DTFR,' dfr=',dfrtop  ,   &
               'LATICE=',Param%hlf, 'CPAIR=',Param%cp_air,'k= ',k
        endif
      endif



end subroutine don_cm_freeze_liquid_k 


!#####################################################################



!######################################################################

subroutine don_cm_clotr_k    &
         (ntr, diag_unit, debug_ijt, Param, sou1, sou2, xe1, xe2, &
          xc1, entrain, dt_micro, xc2, ermesg, error)

!----------------------------------------------------------------------
!    subroutine clotr calculates the in-cloud tracer profiles for 
!    all of the travcers being transported by the donner convection
!    scheme.
!    author :  Leo Donner,  GFDL, 14 Jan 2000
!----------------------------------------------------------------------

use donner_types_mod, only : donner_param_type

implicit none

!----------------------------------------------------------------------
integer,                  intent(in)  ::  ntr
integer,                  intent(in)  ::  diag_unit
logical,                  intent(in)  ::  debug_ijt
type(donner_param_type),  intent(in)  ::  Param
real,   dimension(ntr),   intent(in)  ::  sou1 
real,   dimension(ntr),   intent(in)  ::  sou2 
real,   dimension(ntr),   intent(in)  ::  xe1  
real,   dimension(ntr),   intent(in)  ::  xe2  
real,   dimension(ntr),   intent(in)  ::  xc1  
real,                     intent(in)  ::  entrain, dt_micro
real,   dimension(ntr),   intent(out) ::  xc2  
character(len=*),         intent(out) ::  ermesg
integer,                  intent(out) ::  error

!----------------------------------------------------------------------
!   intent(in) variables:
!
!        sou1          in-cloud source of tracer at bottom of layer 
!                      [ kg / kg(air) / sec ]  
!        sou2          in-cloud source of tracer at top of layer 
!                      [ kg / kg(air) / sec ]  
!        xe1           environmental tracer concentration at bottom of
!                      layer [ kg / kg(air) ]
!        xe2           environmental tracer concentration at top of
!                      layer [ kg / kg(air) ]
!        xc1           in-cloud tracer concentration at bottom of
!                      layer [ kg / kg(air) ]
!        entrain       entrainment factor [ dimensionless ]
!        dt_micro      microphysics time step
!        debug_ijt     is this a diagnostics column ?
!        diag_unit     output unit number for this diagnostics column
!
!   intent(out) variables:
!
!        xc2           in-cloud tracer concentration at top of
!                      layer [ kg / kg(air) ]
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      real    ::  xeb        ! layer-average environmental value of 
                             ! tracer n  [ kg / kg(air) ]
      real    ::  seb        ! layer-average source for tracer n
                             ! [ kg / kg(air) / sec ]
      real    ::  mass_ratio ! ratio of entraining parcel mass at top
                             ! of layer to that at bottom 
                             ! [ dimensionless ]
      real    ::  delta_m    ! fractional amount of mass added to 
                             ! parcel in current model layer 
                             ! [ dimensionless ]
      integer ::  n          ! do-loop index

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!-----------------------------------------------------------------------
!    if in diagnostics column, print out an entry message.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a)') 'in clotr: entering clotr'
      endif

!---------------------------------------------------------------------
!    define ratio of parcel mass at top of layer to that at bottom 
!    (mass_ratio) and the increase in parcel mass across the layer.
!---------------------------------------------------------------------
      mass_ratio = exp(entrain)
      delta_m = mass_ratio - 1.0

!----------------------------------------------------------------------
!    loop over the tracers transported by donner convection.
!----------------------------------------------------------------------
      do n=1,ntr             

!---------------------------------------------------------------------
!    define layer-mean environmental tracer (xeb) and the tracer source 
!    (seb).
!---------------------------------------------------------------------
        xeb = 0.5*(xe1(n) + xe2(n))
        seb = 0.5*(sou1(n) + sou2(n))

!--------------------------------------------------------------------
!    define in-cloud tracer amount at top of layer (xc2) as the sum of 
!    the value at bottom of layer (xc1) plus the amount of environ-
!    mental tracer entrained into the parcel plus the amount produced
!    by the internal source (seb). the time of insertion of the internal
!    source is given by tr_insert_time (in terms of parcel's mass 
!    increase during its traversal of the layer). if source is assumed
!    to be made available at bottom of layer (tr_insert_time = 0), then
!    seb is unmodified; if not available until completion of traversal,
!    then tr_insert_time is 1.0, and amount supplied is reduced by the
!    ratio of parcel mass between layer top and bottom. renormalize the
!    mixing ratio by the total parcel mass at the top of the layer.
!!! BUG :: seb must be multiplied by dt (perhaps as wv/dz ???)
!!!  BUG FIXED 6/2/05
!!!! BUG :: shouldn't the 1 + epm*(   ) be divided into seb  ????
!--------------------------------------------------------------------
        xc2(n) = (xc1(n) + delta_m*xeb + seb*dt_micro*   &
                              (1. + Param%tr_insert_time*delta_m))/mass_ratio 

!--------------------------------------------------------------------
!    assure that the in-cloud tracer mixing ratio is non-negative.
!--------------------------------------------------------------------
        if (xc2(n) < 0.) xc2(n) = 0.

!--------------------------------------------------------------------
!    if in diagnostics column, output the tracer mixing ratio at top
!    of layer (xc2), the layer-mean environmental tracer mixing ratio
!    (xeb) and the layer-mean tracer source (seb).
!--------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, e20.12)')   &
           'in clotr: xc= ',xc2(n)
          write (diag_unit, '(a, e20.12)')   &
               'in clotr: xeb= ',xeb
          write (diag_unit, '(a, e20.12)')   &
               'in clotr: seb= ',seb
        endif
      end do

!--------------------------------------------------------------------


end subroutine don_cm_clotr_k

!#######################################################################


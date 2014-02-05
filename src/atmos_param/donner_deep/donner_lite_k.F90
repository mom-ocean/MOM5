!VERSION NUMBER:
!   $Id: donner_lite_k.F90,v 20.0 2013/12/13 23:17:21 fms Exp $

!######################################################################
!######################################################################


subroutine don_c_def_conv_env_miz    &
         (isize, jsize, nlev_lsm, ntr, dt, Nml, Param, Initialized, &
           Col_diag,  &
          tracers, pblht, tkemiz, qstar, cush, land, coldT,       &!miz
          temp, mixing_ratio, pfull, phalf, zfull, &
          zhalf, lag_cape_temp, lag_cape_vapor, current_displ,   &
          cbmf, Don_cape, Don_conv, sd, Uw_p, ac)

use donner_types_mod,      only : donner_nml_type, donner_param_type, &
                                  donner_column_diag_type,   &
                                  donner_initialized_type, &
                                  donner_cape_type, donner_conv_type
use  conv_utilities_k_mod, only : pack_sd_lsm_k, extend_sd_k,   &
                                  adi_cloud_k, adicloud, sounding, &
                                  uw_params, qt_parcel_k

implicit none

integer,                                    intent(in)    ::   &
                                             isize, jsize, nlev_lsm, ntr
real,                                     intent(in)    :: dt
type(donner_nml_type),                    intent(in)    :: Nml      
type(donner_param_type),                  intent(in)    :: Param
type(donner_initialized_type), intent(in)             :: Initialized
type(donner_column_diag_type), intent(in)    :: Col_diag
real, dimension(isize,jsize,nlev_lsm),    intent(in)    ::    &
                                             temp, mixing_ratio,  &
                                             pfull, zfull,  &
                                             lag_cape_temp,&
                                             lag_cape_vapor
real, dimension(isize,jsize,nlev_lsm,ntr),intent(in)    :: tracers
real, dimension(isize,jsize,nlev_lsm+1),  intent(in)    :: phalf,  &
                                                            zhalf
real, dimension(isize,jsize),             intent(in)    ::   &
                                             current_displ, cbmf, &
                                             pblht, tkemiz, qstar, cush,land
logical, dimension(isize,jsize),          intent(in)    :: coldT
type(donner_cape_type),                   intent(inout) :: Don_cape
type(donner_conv_type),                   intent(inout) :: Don_conv
type(sounding),                           intent(inout) :: sd
type(adicloud),                           intent(inout) :: ac
type(uw_params),                           intent(inout) :: Uw_p

      real, dimension (nlev_lsm) :: mid_cape_temp, mid_cape_vapor
      real, dimension (isize,jsize,nlev_lsm, ntr) :: xgcm_v
      real         :: zsrc, psrc, hlsrc, thcsrc, qctsrc, &
                      lofactor
      integer      :: i, j, k, n

      do n=1,ntr
      do k=1,nlev_lsm
        xgcm_v(:,:,k,n) = tracers(:,:,nlev_lsm-k+1,n)
      end do
      enddo

!--------------------------------------------------------------------
!    if in diagnostics window, write message indicating lag-time cape
!    calculation is being done.
!--------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          write (Col_diag%unit_dc(n), ' (//, a)')  &
         '               CAPE calculation for LAG and MID time profile'
        end do
      endif


      do j=1,jsize
        do i=1,isize
           if (Initialized%using_unified_closure .and.   &
                                          cbmf(i,j) == 0.0) cycle
          if ((current_displ(i,j) .lt. 0) .or.   &
              (.not.Nml%use_llift_criteria)) then
            call pack_sd_lsm_k (Nml%do_lands, land(i,j), coldT(i,j), &
                                dt, pfull(i,j,:), phalf(i,j,:), &
                                zfull(i,j,:), zhalf(i,j,:), &
                                lag_cape_temp(i,j,:),  &
                                lag_cape_vapor(i,j,:),   &
                                xgcm_v(i,j,:,:), sd)
            call extend_sd_k(sd, pblht(i,j), .false., Uw_p)         
            zsrc  =sd%zs (1)
            psrc  =sd%ps (1)
            thcsrc=sd%thc(1)
            qctsrc=sd%qct(1)
            hlsrc =sd%hl (1)
            if (Nml%do_lands) then
               call qt_parcel_k (sd%qs(1), qstar(i,j), pblht(i,j), &
                                 tkemiz(i,j), sd%land, &
                       Nml%gama,Nml%pblht0,Nml%tke0,Nml%lofactor0, &
                                      Nml%lochoice,qctsrc,lofactor)
            end if
            call adi_cloud_k (zsrc, psrc, hlsrc, thcsrc, qctsrc, sd,  &
                              Uw_p, .false., Nml%do_freezing_for_cape, &
                              ac)                
            zsrc  =sd%zs (1)
            Don_cape%xcape_lag(i,j) = ac%cape
            Don_cape%qint_lag (i,j) = sd%qint
!  additional column diagnostics should be added here
          
            mid_cape_temp (:) = temp(i,j,:)
            mid_cape_vapor(:) = mixing_ratio(i,j,:)
!            mid_cape_vapor(:) = mixing_ratio(i,j,:)/  &
!                                                (1.+mixing_ratio(i,j,:))
        
            do k=nlev_lsm - Nml%model_levels_in_sfcbl + 1, nlev_lsm
              mid_cape_temp (k) = lag_cape_temp(i,j,k)
              mid_cape_vapor(k) = lag_cape_vapor(i,j,k)
            end do
        
            call pack_sd_lsm_k (Nml%do_lands, land(i,j), coldT(i,j), &
                                dt, pfull(i,j,:), phalf(i,j,:),   &
                                zfull(i,j,:), zhalf(i,j,:), &
                                mid_cape_temp(:), mid_cape_vapor(:), &
                                xgcm_v(i,j,:,:), sd)
        
            sd%ql(:)=0.; !max(qlin(i,j,:),0.);
            sd%qi(:)=0.; !max(qiin(i,j,:),0.);

            call extend_sd_k (sd, pblht(i,j), .false., Uw_p)
            zsrc  =sd%zs (1)
            psrc  =sd%ps (1)
            thcsrc=sd%thc(1)
            qctsrc=sd%qct(1)
            hlsrc =sd%hl (1)
            if (Nml%do_lands) then
               call qt_parcel_k (sd%qs(1), qstar(i,j), pblht(i,j), &
                                 tkemiz(i,j), sd%land, &
                       Nml%gama,Nml%pblht0,Nml%tke0,Nml%lofactor0, &
                                      Nml%lochoice,qctsrc,lofactor)
            end if

            call adi_cloud_k (zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, &
                              Uw_p, &
                              .false., Nml%do_freezing_for_cape, ac)
            Don_cape%plfc(i,j) = ac%plfc
            Don_cape%plzb(i,j) = ac%plnb
            Don_cape%plcl(i,j) = ac%plcl
!           Don_cape%parcel_r(i,j,:) = ac%qv(:)
            Don_cape%parcel_r(i,j,:) = ac%qv(:)/(1. - ac%qv(:))
            Don_cape%parcel_t(i,j,:) = ac%t (:)
            Don_cape%coin (i,j) = ac%cin
            Don_cape%xcape(i,j) = ac%cape
            Don_cape%qint (i,j) = sd%qint
!           Don_cape%model_r(i,j,:) = sd%qv(:)
            Don_cape%model_r(i,j,:) = sd%qv(:)/(1. - sd%qv(:))
            Don_cape%model_t(i,j,:) = sd%t (:)
            Don_cape%model_p(i,j,:) = sd%p (:)
!           Don_cape%env_r  (i,j,:) = sd%qv(:)
            Don_cape%env_r  (i,j,:) = sd%qv(:)/ (1. - sd%qv(:))
            Don_cape%env_t  (i,j,:) = sd%t (:)
            Don_cape%cape_p (i,j,:) = sd%p (:)
!  additional column diagnostics should be added here
          endif
        end do
      end do

!----------------------------------------------------------------------


end subroutine don_c_def_conv_env_miz


!######################################################################
!######################################################################


subroutine don_d_integ_cu_ensemble_miz             &
        (nlev_lsm, nlev_hires, ntr, me, diag_unit, debug_ijt, Param,   &
         Col_diag, Nml, Initialized, temp_c, mixing_ratio_c, pfull_c, & 
         phalf_c, pblht, tkemiz, qstar, cush, land, coldT, delt, &
         sd, Uw_p, ac, cp, ct, tracers_c, sfc_sh_flux_c,   &
         sfc_vapor_flux_c, sfc_tracer_flux_c, plzb_c, exit_flag_c, &
         ensmbl_precip, ensmbl_cond, ensmbl_anvil_cond_liq,  &
         ensmbl_anvil_cond_liq_frz, ensmbl_anvil_cond_ice, pb, pt_ens, &
         ampta1, amax, emsm, rlsm, cld_press, ensmbl_melt,  &
         ensmbl_melt_meso,  ensmbl_freeze, ensmbl_freeze_meso, &
         ensmbl_wetc, disb, disc_liq, disc_ice, dism_liq,  &
         dism_liq_frz, dism_liq_remelt, dism_ice, dism_ice_melted, &
         disp_liq, disp_ice, disz, disz_remelt, disp_melted, disze1, &
         disze2, disze3, disd, disv, disg_liq, disg_ice, &
         enctf, encmf, enev, ecds_liq, ecds_ice, eces_liq, &
         eces_ice, ensmbl_cloud_area, cuq, cuql_v, &
         detmfl, uceml, qtren, etsm, lmeso, frz_frac, &
         meso_frz_intg_sum, ermesg, error, melting_in_cloud, &
         i, j, Don_cem)

!----------------------------------------------------------------------
!    subroutine integrate_cumulus_ensemble works on a single model 
!    column. all profile arrays used in this subroutine and below have 
!    index 1 nearest the surface. it first determines the lifting con-
!    densation level (if one exists) of a parcel moving from the 
!    specified parcel_launch_level. if an lcl is found, subroutine 
!    donner_cloud_model_cloud_model is called to determine the behavior
!    of each of kpar cloud ensemble members assumed present in the 
!    column (each ensemble member is assumed to have a different en-
!    trainment rate). if all ensemble members produce deep convection, 
!    the ensemble statistics are produced for use in the large-scale 
!    model; otherwise deep convection is not seen in the large-scale 
!    model in this grid column. if the ensemble will support a mesoscale
!    circulation, its impact on the large-scale model fields is also 
!    determined. upon completion, the appropriate output fields needed 
!    by the large-scale model are returned to the calling routine.
!----------------------------------------------------------------------
use donner_types_mod,     only : donner_initialized_type,   &
                                 donner_param_type, donner_nml_type, &
                                 donner_column_diag_type, donner_cem_type
use  conv_utilities_k_mod,only : qsat_k, exn_k, adi_cloud_k, adicloud, &
                                 sounding, uw_params, qt_parcel_k
use  conv_plumes_k_mod,   only : cumulus_plume_k, cumulus_tend_k, &
                                 cplume, ctend, cpnlist, cwetdep_type

implicit none 

!----------------------------------------------------------------------
integer,                           intent(in)    :: nlev_lsm,    &
                                                    nlev_hires, ntr, &
                                                    me, diag_unit
logical,                           intent(in)    :: debug_ijt
type(donner_param_type),           intent(in)    :: Param
type(donner_column_diag_type),     intent(in)    :: Col_diag
type(donner_nml_type),             intent(in)    :: Nml   
type(donner_initialized_type),     intent(in)    :: Initialized
real,    dimension(nlev_lsm),      intent(in)    :: temp_c,   &
                                                    mixing_ratio_c,   &
                                                    pfull_c
real,    dimension(nlev_lsm+1),    intent(in)    :: phalf_c
real,                              intent(in)    :: pblht, tkemiz,  &
                                                    qstar, cush, land, delt
logical,                           intent(in)    :: coldT
type(sounding),                    intent(inout) :: sd
type(uw_params),                    intent(inout) :: Uw_p
type(adicloud),                    intent(inout) :: ac
type(cplume),                      intent(inout) :: cp
type(ctend),                       intent(inout) :: ct
real,    dimension(nlev_lsm,ntr),  intent(in)    :: tracers_c           
real,                              intent(in)    :: sfc_sh_flux_c,   &
                                                    sfc_vapor_flux_c 
real,    dimension(ntr),           intent(in)    :: sfc_tracer_flux_c 
real,                              intent(in)    :: plzb_c
logical,                           intent(inout) :: exit_flag_c  
real,                              intent(out)   ::  &
                                      ensmbl_precip, ensmbl_cond,&
                                      ensmbl_anvil_cond_liq, &
                                      ensmbl_anvil_cond_liq_frz, &
                                      ensmbl_anvil_cond_ice, &
                                      pb, pt_ens, ampta1, amax
real,    dimension(nlev_lsm),      intent(out)   :: emsm, rlsm,  &
                                                    cld_press 
real,    dimension(nlev_lsm),      intent(out)   ::   &
                                      ensmbl_melt, ensmbl_melt_meso,&
                                      ensmbl_freeze,   &
                                      ensmbl_freeze_meso, disb, disd, &
                                      disv, disc_liq, disc_ice,&
                                      dism_liq, dism_ice,  &
                                      dism_ice_melted, dism_liq_frz, &
                                      dism_liq_remelt, disp_liq,  &
                                      disp_ice,disp_melted, &
                                      disz_remelt, disz, disze1,  &
                                      disze2, disze3, enctf, encmf,&
                                      disg_liq, disg_ice, enev,    &
                                      ecds_liq, ecds_ice,&
                                      eces_liq, eces_ice,&
                                      ensmbl_cloud_area, cuq, cuql_v, &
                                      detmfl, uceml
real,    dimension(nlev_lsm,ntr),  intent(out)   :: qtren, ensmbl_wetc
real,    dimension(nlev_lsm,ntr),intent(out)     :: etsm 
logical,                           intent(out)   :: lmeso       
integer,                           intent(in)    :: i, j
type(donner_cem_type),             intent(inout) :: Don_cem
real   ,                           intent(out)   :: frz_frac
logical,                           intent(out)   :: meso_frz_intg_sum
logical ,                          intent(out)    :: melting_in_cloud

character(len=*),                  intent(out)   :: ermesg
integer,                           intent(out)   :: error

!---------------------------------------------------------------------
!   intent(in) variables:
! 
!     nlev_lsm       number of model layers in large-scale model
!     nlev_hires     number of model layers in hi-res cloud model
!                    of the donner deep convection parameterization
!     ntr            number of tracers to be transported by donner
!                    convection
!     me             local pe number
!     diag_unit      unit number for column diagnostics output, if 
!                    diagnostics are requested for the current column
!     debug_ijt      logical indicating whether current column requested
!                    column diagnostics
!     Param          donner_param_type variable containingthe parameters
!                    of the donner deep convection parameterization
!     Col_diag       donner_column_diagtype variable containing the
!                    information defining the columns fro which diagnos-
!                    tics are desired.
!     Nml            donner_nml_type variable containing the donner_nml
!                    variables that are needed outsied of donner_deep_mod
!     temp_c         temperature field at model full levels 
!                    index 1 nearest the surface [ deg K ]
!     mixing_ratio_c        vapor mixing ratio at model full levels 
!                    index 1 nearest the surface
!                    [ kg(h2o) / kg(dry air) ]
!     pfull_c         pressure field at large-scale model full levels 
!                    index 1 nearest the surface [ Pa ]
!     phalf_c        pressure field at large-scale model half-levels 
!                    index 1 nearest the surface [ Pa ]
!     tracers_c      tracer fields that are to be transported by donner
!                    convection.  index 1 nearest the surface 
!                    [ kg (tracer) / kg (dry air) ]
!     sfc_sh_flux_c  sensible heat flux across the surface
!                    [ watts / m**2 ]
!     sfc_vapor_flux_c water vapor flux across the surface
!                    [ kg(h2o) / (m**2 sec) ]
!     sfc_tracer_flux_c  
!                    flux across the surface of tracers transported by
!                    donner_deep_mod [ kg(tracer) / (m**2 sec) ]
!     plzb_c         level of zero buoyancy for a parcel lifted from
!                    the parcel_launch_level.  [ Pa ]
!
!   intent(inout) variables:
!
!     exit_flag_c    logical indicating whether donner convection
!                    is not active (.true.) or is active (.false.) in
!                    current model column 
!
!     cumulus ensemble member fields (see also donner_types.h):
!
!     --- single level ---
!
!     Don_cem_cell_precip 
!                    area weighted convective precipitation rate
!                    [ mm/day ]
!     Don_cem_pb     pressure at cloud base for ensemble (currently,
!                    all ensemble members have same base) [ Pa ]
!     Don_cem_ptma   pressure at cloud top for ensemble [ Pa ]
!
!     --- lo-res multi-level ---
! 
!     Don_cem_h1     condensation rate profile on lo-res grid
!                    for the current ensemble member
!                    [ ( kg(h2o) ) / ( kg( dry air) sec ) ] 
!
!     --- lo-res multi-level ---
!
!     Don_cem_qlw    profile of cloud water for the current ensemble
!                    member [ kg(h2o) / kg(air) ]
!     Don_cem_cfracice
!                    fraction of condensate that is ice [ fraction ]
!     Don_cem_wv     vertical velocity profile [ m / s ]
!     Don_cem_rcl    cloud radius profile [ m ]
!
!   intent(out) variables:
!    
!     ensmbl_precip      sum of precipitation rate over ensemble members,
!                        # 1 to the current, weighted by the area at 
!                        cloud base of each member
!                        [ mm / day ]
!     ensmbl_cond        sum of condensation rate over ensemble members,
!                        # 1 to the current, weighted by the area at 
!                        cloud base of each member
!                        [ mm / day ]
!     ensmbl_anvil_cond  sum of rate of transfer of condensate from cell 
!                        to anvil over ensemble members, # 1 to the c
!                        current, weighted by the area at cloud base of 
!                        each member [ mm / day ]
!     pb                 pressure at cloud base for ensemble (all ensem-
!                        ble members have same base) [ Pa ]
!     pt_ens             pressure at cloud top for the ensemble (top 
!                        pressure of deepest ensemble member) [ Pa ]
!     ampta1             cloudtop anvil area (assumed to be five times
!                        larger than the sum of the cloud top areas of 
!                        the ensemble members, as in Leary and Houze 
!                        (1980).  [ fraction ]
!     amax               maximum allowable area of cloud base that is
!                        allowed; if cloud base area is larger than 
!                        amax, the cloud fractional area somewhere in
!                        the grid box would be greater than one, which 
!                        is non-physical.
!     emsm               vertical profile on the hi-res grid of vertical
!                        moisture flux convergence, summed over ensemble 
!                        members # 1 to the current, each member's cont-
!                        ribution being weighted by its cloud area at 
!                        level k relative to the cloud base area of 
!                        ensemble member #1  
!                        [ kg (h2o) / ( kg(dry air) sec ) ]
!     rlsm               vertical profile on the hi-res grid of conden-
!                        sation rate, summed over ensemble members # 1 to
!                        the current, each member's contribution being 
!                        weighted by its cloud area at level k relative 
!                        to the cloud base area of ensemble member #1
!                        [ ( kg(h2o) ) / ( kg( dry air) sec ) ] 
!     cld_press          pressures at hi-res model levels [ Pa ]
!     ensmbl_melt        vertical profile on the lo-res grid of ice melt,
!                        both from the cells and any mesoscale circul-
!                        ation, summed over ensemble members # 1 to the 
!                        current, each member's contribution being 
!                        weighted by its cloud area at level k relative !
!                        to the cloud base area of ensemble member #1
!                        [ kg(h2o) / kg (dry air) ]
!     ensmbl_freeze      vertical profile on the lo-res grid of freezing,
!                        both from the cells and any mesoscale circul-
!                        ation, summed over ensemble members # 1 to the 
!                        current, each member's contribution being 
!                        weighted by its cloud area at level k relative !
!                        to the cloud base area of ensemble member #1
!                        [ kg(h2o) / kg (dry air) ]
!     disg               vertical profile on the lo-res grid of the      
!                        latent heat term in the temperature equation
!                        associated with the evaporation of condensate
!                        in the convective downdraft and updraft,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ deg K / day ] 
!     enev               vertical profile on the lo-res grid of the      
!                        cloud-area-weighted profile of the potential
!                        cloud water evaporation, summed over ensemble 
!                        members # 1 to the current, each member's con-
!                        tribution being weighted by its cloud area at !
!                        level k relative to the cloud base area of 
!                        ensemble member #1.  this amount of water
!                        must be evaporated if it turns out that there is
!                        no mesoscale circulation generated in the 
!                        column.
!                        [ ( kg(h2o) ) / ( kg(dry air) sec ) ] 
!     enctf              vertical profile on the lo-res grid of the entr-
!                        opy forcing, consisting of the sum of the
!                        vertical entropy flux convergence and the latent
!                        heat release, summed over 
!                        ensemble members # 1 to the current, each mem-
!                        ber's contribution being weighted by its cloud 
!                        area at level k relative to the cloud base area
!                        of ensemble member #1
!                        [ deg K / day ]                        
!     encmf              vertical profile on the lo-res grid of the      
!                        moisture forcing, consisting of the sum of the
!                        vertical moisture flux convergence and the cond-
!                        ensation, summed over ensemble members # 1 to 
!                        the current, each member's contribution being 
!                        weighted by its cloud area at level k relative 
!                        to the cloud base area of ensemble member #1
!                        [ ( kg(h2o) ) / ( kg( dry air) day ) ] 
!     disb               vertical profile on the lo-res grid of the      
!                        temperature flux convergence, summed over 
!                        ensemble members # 1 to the current, each mem-
!                        ber's contribution being weighted by its cloud 
!                        area at level k relative to the cloud base area 
!                        of ensemble member #1.  
!                        [ deg K / day ] 
!     disc               vertical profile on the lo-res grid of the      
!                        latent heat term in the temperature equation, 
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ deg K / day ] 
!     disd               vertical profile on the lo-res grid of the      
!                        vertical moisture flux convergence,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        lo-res grid for the current ensemble member 
!                        [  g(h2o) / ( kg(dry air) day ) ]
!     ecds               vertical profile on the lo-res grid of the      
!                        condensate evaporated in convective downdraft,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ g(h2o) / kg(air) / day ]
!     eces               vertical profile on the lo-res grid of the      
!                        condensate evaporated in convective updraft,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ g(h2o) / kg(air) / day ]
!     ensmbl_cloud_area  total cloud area profile over all ensemble
!                        members on large_scale model grid [ fraction ]
!     cuq                ice water profile on large-scale model grid, 
!                        normalized by ensemble cloud area.
!     cuql_v             liquid water profile on large-scale model grid, 
!                        normalized by ensemble cloud area.
!     uceml              upward mass flux on large_scale model grid     
!                        [ kg (air) / (sec m**2) ]
!     detmfl             detrained mass flux on large-scale model grid
!                        normalized by ensemble cloud area
!                        [ kg (air) / (sec m**2) ]
!     etsm               vertical profile on the hi-res grid of vertical
!                        tracer flux convergence, summed over ensemble 
!                        members # 1 to the current, each member's con-
!                        tribution being weighted by its cloud area at i
!                        level k relative to the cloud base area of 
!                        ensemble member #1 
!                        [ kg (tracer) / ( kg(dry air) sec ) ]
!     qtren              vertical profile on the lo-res grid of the      
!                        vertical tracer flux convergence,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ kg(tracer) / ( kg(dry air) sec ) ]
!     ensmbl_wetc        vertical profile on the lo-res grid of the              
!                        tracer wet deposition,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ kg(tracer) / ( kg(dry air) sec ) ]
!     lmeso              logical variable; if .false., then it has been
!                        determined that a mesoscale circulation cannot
!                        exist in the current column. final value not
!                        determined until all ensemble members have been
!                        integrated. 
!     ermesg             character string containing any error message
!                        that is returned from a kernel subroutine
!
!---------------------------------------------------------------------

!     cmui             normalized vertical integral of mesoscale-updraft
!                      deposition (kg(H2O)/((m**2) sec)
!     cmus(nlev)       normalized mesoscale-updraft deposition
!                      (kg(H2O)/kg/sec)
!     emds(nlev)       normalized mesoscale-downdraft sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     emei             normalized vertical integral of mesoscale-updraft
!                      sublimation (kg(h2O)/((m**2) sec)
!     emes(nlev)       normalized mesoscale-updraft sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     disa(nlev)       normalized thermal forcing, cells+meso (K/sec)
!                      (excludes convergence of surface heat flux)
!                      index 1 at ground. Cumulus thermal forcing defined
!                      as in Fig. 3 of Donner (1993, JAS).
!     disb(nlev)       normalized cell entropy-flux convergence (K/sec)
!                      (excludes convergence of surface flux)
!                      index 1 at ground. Entropy-flux convergence divided
!                      by (p0/p)**(rd/cp).
!     disc(nlev)       normalized cell condensation/deposition
!                      (K/sec)
!                      index 1 at ground.
!     disd(nlev)       normalized cell moisture-flux convergence
!                      (excludes convergence of surface moisture flux)
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     dise(nlev)       normalized moisture forcing, cells+meso (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     dmeml(nlev)      mass flux in mesoscale downdraft (kg/((m**2) s))
!                      (normalized by a(1,p_b)) (index 1 at atmosphere
!                      bottom)
!     elt(nlev)        normalized melting (K/sec)
!                      index 1 at ground.
!     fre(nlev)        normalized freezing (K/sec)
!                      index 1 at ground.
!     pb               pressure at base of cumulus updrafts (Pa)
!     pmd              pressure at top of mesoscale downdraft (Pa)
!     pztm             pressure at top of mesoscale updraft (Pa)
!     mrmes(nlev)       normalized mesoscale moisture-flux convergence
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     qtmes(nlev,ncont)  tracer tendency due to mesoscale tracer-flux
!                        convergence (kg/kg/s) (normalized by a(1,p_b))
!                        index 1 at ground 
!     qtren_v          normalized tracer tendency due to cells...
!                      (lon,lat,vert,tracer index)
!                      Vertical index increases as height increases.
!     sfcq(nlev)       boundary-layer mixing-ratio tendency due to surface
!                      moisture flux (kg(H2O)/kg/sec)
!     sfch(nlev)       boundary-layer heating due to surface heat flux
!                      (K/sec)
!     tmes(nlev)       normalized mesoscale entropy-flux convergence
!                      (K/sec)
!                      Entropy-flux convergence is mesoscale component
!                      of second term in expression for cumulus thermal
!                      forcing in Fig. 3 of Donner (1993, JAS).
!                      index 1 at ground.
!     tpre_v           total normalized precipitation (mm/day)
!     detmfl(nlev)     detrained mass flux from cell updrafts
!                      (normalized by a(1,p_b))
!                      (index 1 near atmosphere bottom)
!                      (kg/((m**2)*s)
!     uceml(nlev)      normalized mass fluxes in cell updrafts
!                      (kg/((m**2)*s) 
!     umeml(nlev)      mass flux in mesoscale updraft (kg/((m**2) s))
!                      (normalized by a(1,p_b)) (index 1 at atmosphere
!                      bottom)
!                      index 1 at ground.
!     wmms(nlev)       normalized mesoscale deposition of water vapor from
!                      cells (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     wmps(nlev)       normalized mesoscale redistribution of water vapor
!                      from cells (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     wtp_v            tracer redistributed by mesoscale processes
!                      (kg/kg/s) (normalized by a(1,p_b))
!                      vertical index increases with increasing height
!                      (lon,lat,vert,tracer index)
!--------------------------------------------------------------------


!!  UNITS
!!    ucemh  [kg /sec / m**2 ]
!!    detmfh [kg /sec / m**2 ]
!!    conint [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    precip [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    q1     [ kg(h2o) / kg(air) / sec ]
!!    h1     [ kg(h2o) / kg(air) / sec ]
!!    cmf    [ g(h2o) / kg(air) /day ]
!!    rlh    [ kg(h2o) / kg(air) / day ]  * [ L / Cp ] = [ deg K / day ]
!!    h1_2   [ deg K / sec ]
!!    efc    [ deg K / day ]
!!    efchr  [ deg K / sec ]
!!    ehfh   [ kg(air) (deg K) / (sec**3 m)
!!    ctf    [ deg K / day ]
!!    disb_v [ deg K / day ]
!!    disc_v [ deg K / day ] 
!!    disn   [ deg K / day ] 
!!    ecd    [ g(h2o) / kg(air) / day ]
!!    ece    [ g(h2o) / kg(air) / day ]
!!    ecds_v [ g(h2o) / kg(air) / day ]
!!    eces_v [ g(h2o) / kg(air) / day ]
!!    pf     [ (m**2 kg(h2o)) / (kg(air) sec) ]
!!    dpf    [ (m**2 kg(h2o)) / (kg(air) sec) ] ==>   
!!                                          [ kg(h2o)) / (kg(air) sec) ]
!!    qlw2   [ kg(h2o)) / (kg(air) sec) ]
!!    qlw    [ kg(h2o)) / kg(air) ]
!!    evap   [ kg(h2o)) / kg(air) ]
!!    evap_rate [ kg(h2o)) / (kg(air) sec) ]




!        cape     convective available potential energy (J/kg)
!        cin      convective inhibtion (J/kg)
!        cpd      specific heat of dry air at constant pressure (J/(kg K))
!        cpv      specific heat of water vapor [J/(kg K)]
!        dcape    local rate of CAPE change by all processes
!                 other than deep convection [J/(kg s)]
!        dqls     local rate of change in column-integrated vapor
!                 by all processes other than deep convection
!                 {kg(H2O)/[(m**2) s]}
!        epsilo   ratio of molecular weights of water vapor to dry air
!        gravm    gravity constant [m/(s**2)]
!        ilon     longitude index
!        jlat     latitude index
!        mcu      frequency (in time steps) of deep cumulus
!        current_displ  integrated low-level displacement (Pa)
!        cape_p   pressure at Cape.F resolution (Pa)
!                 Index 1 at bottom of model.
!        plfc     pressure at level of free convection (Pa)
!        plzb_c   pressure at level of zero buoyancy (Pa)
!        pr       pressure at Skyhi vertical resolution (Pa)
!                 Index 1 nearest ground  
!        q        large-scale vapor mixing ratio at Skyhi vertical resolution
!                 [kg(h2O)/kg]
!                 Index 1 nearest ground 
!        qlsd     column-integrated vapor divided by timestep for cumulus
!                 parameterization {kg(H2O)/[(m**2) s]}
!        r        large-scale vapor mixing ratio at Cape.F resolution
!                 [kg(h2O)/kg]
!                 Index 1 at bottom of model.
!        rpc      parcel vapor mixing ratio from Cape.F [kg(h2O)/kg]
!                 Index 1 at bottom of model.
!        rd       gas constant for dry air (J/(kg K))
!        rlat     latent heat of vaporization (J/kg)
!        rv       gas constant for water vapor (J/(kg K))
!        t        large-scale temperature at Skyhi vertical resolution (K)
!                 Index 1 nearest ground
!        tcape    large-scale temperature at Cape.F resolution (K)
!                 Index 1 at bottom of model.
!        tpc      parcel temperature from from Cape.F (K)
!                 Index 1 at bottom of model.
!

!----------------------------------------------------------------------
!   local variables:

      real,    dimension (nlev_hires)     ::                &
              alp, ucemh, cuql, cuqli, detmfh, tcc

      real,    dimension (nlev_lsm)       ::           &
              q1, cmf, cell_freeze, cell_melt, &
              h1_liq, h1_ice, meso_melt, meso_freeze, h1_2, &
              evap_rate, ecd, ecd_liq, ecd_ice, ece, ece_liq, ece_ice
      real,   dimension (nlev_lsm) :: rcl_miz, dpf_miz, qlw_miz,  &
                                      dfr_miz, flux_miz, efchr_miz, &
                                      emfhr_miz, cfracice_miz, alp_miz,&
                                      cuql_miz, &
                                      cuqli_miz, ucemh_miz, detmfh_miz,&
                                      rlsm_miz, emsm_miz, qvfm_miz,&
                                      qvfm_tot
      real,   dimension (nlev_lsm,ntr) :: etsm_miz, etfhr_miz, dpftr_miz

      real    :: dint_miz, cu_miz, cell_precip_miz,              &
                 apt_miz
      real    :: wt_factor
      integer :: krel, ncc_kou_miz
      real,    dimension (nlev_lsm,ntr)   :: qtr
      real,    dimension (Param%kpar)     :: ptma_miz
      integer, dimension (Param%kpar)     :: ncca

      logical ::   lcl_reached                  
      integer ::   ncc_ens
      integer ::   k,    kou, n
      integer   :: kk
      real    ::   al, dp, mrb, summel, &
                   dmela, ca_liq, ca_ice, &
                   tb, alpp, pcsave, ensmbl_cld_top_area

      real    :: qs, tp, qp, pp, chi, rhtmp, frac0, lofactor !miz
      real    ::   meso_frac, precip_frac
     real    ::            frz_frac_non_precip
     real    ::        bak
     real   ::                meso_frz_frac
     logical   :: meso_frz_intg               
     real :: pmelt_lsm,                 precip_melt_frac
      real :: ecei_liq
      real   :: ci_liq_cond, ci_ice_cond
     real :: local_frz_frac

      real            :: zsrc, psrc, hlsrc, thcsrc, qctsrc
      real            :: rkm, cbmf, wrel, scaleh
      real, dimension (nlev_lsm  ) ::  dpf_warm, dpf_cold
      type(cpnlist)   :: cpn

!----------------------------------------------------------------------
!   local variables:
!
!      ensmbl_cld_top_area  
!                       sum of the cloud top areas over ensemble members 
!                       # 1 to the current, normalized by the cloud base
!                       area of ensemble member # 1 [ dimensionless ]
!
!----------------------------------------------------------------------

      ermesg = ' ' ; error = 0
!---------------------------------------------------------------------
!    if in diagnostics column, output the large-scale model temperature,
!    vapor mixing ratio and full-level pressure profiles (index 1 near-
!    est the surface).
!---------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_lsm-Col_diag%kstart+1
          write (diag_unit, '(a, i4, f20.14, e20.12, f19.10)')&
                 'in mulsub: k,T,Q,P= ',k, temp_c(k),  &
                                       mixing_ratio_c(k),pfull_c(k)
        end do
      endif

!!$!--------------------------------------------------------------------
!!$!    call don_cm_lcl_k to calculate the temperature (tb), a
!!$!    pressure (pb) and mixing ratio (mrb) at the lifting condensation 
!!$!    level for a parcel starting from the parcel_launch_level. if a sat-
!!$!    isfactory lcl is not reached for this parcel, the logical variable 
!!$!    lcl_reached will be set to .false..
!!$!--------------------------------------------------------------------
!!$      call don_cm_lcl_k    &
!!$           (Param, temp_c (Nml%parcel_launch_level),    &
!!$            pfull_c       (Nml%parcel_launch_level),    &
!!$            mixing_ratio_c(Nml%parcel_launch_level),   &
!!$            tb, pb, mrb, lcl_reached, ermesg)     
!!$      if (trim(ermesg) /= ' ') return

!miz
      tp=temp_c(Nml%parcel_launch_level)
      qp=mixing_ratio_c(Nml%parcel_launch_level)/  &
                           (1.+mixing_ratio_c(Nml%parcel_launch_level))
      pp=pfull_c(Nml%parcel_launch_level)
      qs=qsat_k(tp, pp,Uw_p)
      rhtmp=min(qp/qs,1.)
      chi=tp/(1669.0-122.0*rhtmp-tp)
      pb =pp*(rhtmp**chi); !Emanuel's calculation, results nearly identical to RAS
      mrb=mixing_ratio_c(Nml%parcel_launch_level)
      tb =tp/exn_k(pp,Uw_p)*exn_k(pb,Uw_p)
!miz

!--------------------------------------------------------------------
!    if an acceptable lcl was not reached, set exit_flag_c so that the
!    remaining computations for this column are bypassed, and return to
!    calling routine. 
!--------------------------------------------------------------------
      if (pb > 50000.) then
         lcl_reached=.true.
      else
         lcl_reached=.false.
      end if

!--------------------------------------------------------------------
!    if in diagnostics column and an lcl was defined, output the lcl 
!    temperature, pressure and mixing ratio. if an acceptble lcl was 
!    not reached, print a message.
!--------------------------------------------------------------------
      if (debug_ijt) then
        if (lcl_reached) then
          write (diag_unit, '(a, f20.14, f19.10, e20.12)') &
                                 'in mulsub: tb,pb,qb= ',tb, pb, mrb
        else
          write (diag_unit, '(a)') 'in mulsub: lcl not reached'
        endif
      endif

      if (.not. lcl_reached) then
         exit_flag_c = .true.
         return
      endif

!---------------------------------------------------------------------
!    initialize variables which will accumulate scalar sums over all 
!    ensemble members.
!---------------------------------------------------------------------
      ensmbl_precip       = 0.
      ensmbl_cond         = 0.
      ensmbl_anvil_cond_liq   = 0.
      ensmbl_anvil_cond_liq_frz   = 0.
      ensmbl_anvil_cond_ice   = 0.
      ensmbl_cld_top_area = 0.

!---------------------------------------------------------------------
!    initialize the variables which will contain the sum over the 
!    ensemble members of the vertical profiles of various quantities 
!    on the cloud-model grid.
!---------------------------------------------------------------------
      do k=1,nlev_hires
        cuql(k)   = 0.
        cuqli(k)  = 0.
        ucemh(k)  = 0.
        detmfh(k) = 0.
        alp(k)    = 0.
      end do

      do k=1,nlev_lsm
         rlsm(k)   = 0.
         emsm(k)   = 0.
      end do
      do n=1,ntr
      do k=1,nlev_lsm
         etsm(k,n) = 0.
      end do
      end do


!---------------------------------------------------------------------
!    initialize the variables which will contain the sum over the 
!    ensemble members of the vertical profiles of various quantities 
!    on the large-scale model grid.
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        ensmbl_freeze(k)    = 0.
        ensmbl_freeze_meso(k)    = 0.
        ensmbl_melt(k)    = 0.
        ensmbl_melt_meso(k)    = 0.
        disb(k)    = 0.
        disc_liq(k) = 0.
        disc_ice(k) = 0.
        dism_liq(k) = 0.
        dism_liq_frz(k) = 0.
        dism_liq_remelt(k) = 0.
        dism_ice(k) = 0.
        dism_ice_melted(k) = 0.
        disp_liq(k) = 0.
        disp_ice(k) = 0.
        disp_melted(k) = 0.       
        disd(k)    = 0.
        disv(k)    = 0.
        disz(k) = 0.
        disz_remelt(k) = 0.
        disze1(k) = 0.
        disze2(k) = 0.
        disze3(k) = 0.
        ecds_liq(k)    = 0.
        ecds_ice(k)    = 0.
        eces_liq(k)    = 0.
        eces_ice(k)    = 0.
        enctf(k)   = 0.
        encmf(k)   = 0.
        disg_liq(k)    = 0.
        disg_ice(k)    = 0.
        enev(k)    = 0.
      end do
      do n=1,ntr
      do k=1,nlev_lsm
        qtren(k,n) = 0.
        ensmbl_wetc(k,n) = 0.
      end do
      end do

      alp_miz   = 0.
      cuql_miz  = 0.
      cuqli_miz = 0. 
      ucemh_miz = 0.
      detmfh_miz= 0.
      rlsm_miz  = 0.
      emsm_miz  = 0.
      qvfm_miz  = 0.
      qvfm_tot  = 0.
      etsm_miz  = 0.
      etfhr_miz = 0.
      dpftr_miz = 0.
      ptma_miz  = 200000.
      ncca      = 0.

      ensmbl_cloud_area = 0.
      cuq               = 0.
      cuql_v            = 0.
      uceml             = 0.
      detmfl            = 0.

      ece=0.
      ecd=0.
      evap_rate = 0.
      ampta1 = 0.

      if (Nml%allow_mesoscale_circulation) then
        lmeso = .true.
      else
        lmeso = .false.
      endif

!!$      do k=1,nlev_hires
!!$        cld_press(k) = pb + (k-1)*Param%dp_of_cloud_model
!!$      end do

      pcsave = phalf_c(1)

      k=Nml%parcel_launch_level
      zsrc  =sd%zs (k)
      psrc  =sd%ps (k)
      thcsrc=sd%thc(k)
      qctsrc=sd%qct(k)
      hlsrc =sd%hl (k)
      frac0 = Nml%frac
      if (Nml%do_lands) then
         !frac0 = Nml%frac * ( 1.- 0.5 * sd%land) 
         !frac0 = Nml%frac * ( Nml%pblht0 / max(pblht,  Nml%pblht0))
         !frac0 = Nml%frac * ( Nml%tke0   / max(tkemiz, Nml%tke0  ))
         call qt_parcel_k (sd%qs(k), qstar, pblht, tkemiz, sd%land, &
              Nml%gama, Nml%pblht0, Nml%tke0, Nml%lofactor0, Nml%lochoice, qctsrc, lofactor)          
         frac0 = Nml%frac * lofactor
      endif
      call adi_cloud_k (zsrc, psrc, hlsrc, thcsrc, qctsrc, sd,   &
                        Uw_p, &
                        .false., Nml%do_freezing_for_cape, ac)

       meso_frz_intg_sum = .false.
!--------------------------------------------------------------------
!    loop over the KPAR members of the cumulus ensemble.
!--------------------------------------------------------------------
      do kou=1,Param%kpar

        meso_frz_intg     = .false.
        if (trim(Nml%entrainment_constant_source) == 'gate') then
          alpp = Param%max_entrainment_constant_gate/  &
                           Param%ensemble_entrain_factors_gate(kou)
        else if (trim(Nml%entrainment_constant_source) == 'kep') then
          alpp = Param%max_entrainment_constant_kep/  &
                           Param%ensemble_entrain_factors_kep(kou)
        else
          ermesg = 'invalid entrainment_constant_source'
          error = 1
          return
        endif

         if (debug_ijt) then
           write (diag_unit, '(a)')    &
                     'in mulsub: phalf, temp= :'
           do k=1,nlev_lsm
           write (diag_unit, '(i4, 2f19.10)')    &
                       k, phalf_c(k), temp_c(k)
           end do
        endif
 
        pmelt_lsm = 2.0e05
      if( temp_c(1) >  Param%KELVIN ) then
        do k=1,nlev_lsm-1
         if ((temp_c(k) >= Param%KELVIN) .and.    &
            (temp_c(k+1) <= Param%KELVIN)) then
           pmelt_lsm = phalf_c(k+1)
           exit
        endif
      end do
    endif
 
     if (debug_ijt) then
         write (diag_unit, '(a, 2f19.10)')    &
           'before cm_cloud_model call pb,  pmelt_lsm    = ', &
                                   pb, pmelt_lsm
     endif

!!$!test donner_plumes
!!$        call don_cm_cloud_model_k   &
!!$             (nlev_lsm, nlev_hires, ntr, kou, diag_unit, debug_ijt,   &
!!$              Param, Col_diag, tb, pb, alpp, cld_press, temp_c,   &
!!$              mixing_ratio_c, pfull_c, phalf_c, tracers_c, pcsave,  &
!!$              exit_flag_c, rcl, dpf, qlw, dfr, flux, ptma(kou), &
!!$              dint, cu, cell_precip, dints, apt, cell_melt, efchr, &
!!$              emfhr, cfracice, etfhr, ncc_kou, tcc, wv, ermesg) !miz
!!$        if (trim(ermesg) /= ' ') return
!!$        if (exit_flag_c) return
!!$
!!$!---------------------------------------------------------------------
!!$!    define the cloud water from this ensemble member which must be 
!!$!    evaporated if it turns out that there is no mesoscale circulation 
!!$!    associated with the ensemble.
!!$!---------------------------------------------------------------------
!!$        cld_evap(:) = -dpf(:)*(1. - (cell_precip/cu))
!!$        ptt = ptma(kou) + Param%dp_of_cloud_model
!!$        call don_d_def_lores_model_profs_k        &
!!$             (nlev_lsm, nlev_hires, ntr, ncc_kou, diag_unit, debug_ijt, &
!!$              Param, pb, ptt, sfc_vapor_flux_c, sfc_sh_flux_c,  &
!!$              sfc_tracer_flux_c, pfull_c, phalf_c, cld_press, dpf, dfr, &
!!$              cld_evap, qlw, emfhr, efchr, etfhr, cell_freeze,     &
!!$              evap_rate, h1, h1_2, q1, qtr, ermesg)
!!$        if (trim(ermesg) /= ' ') return
!!$
!!$        call don_d_add_to_ensmbl_sum_lores_k    &
!!$             (nlev_lsm, ntr, diag_unit, debug_ijt, lmeso, Param,   &
!!$              Param%arat(kou), dint, cell_freeze, cell_melt, temp_c,   &
!!$              h1_2, ecd, ece, evap_rate, q1, h1, pfull_c, meso_melt, &
!!$              meso_freeze, phalf_c, qtr, ensmbl_melt, ensmbl_freeze, &
!!$              enctf, encmf, enev, disg, disb, disc, ecds, eces, disd, &
!!$              qtren, ermesg)
!!$!test donner_plumes


       tcc = 0.
       dmela = 0.

!begin: testing unified plume
!!! SHOULD ADD SOME COLUMN DIAGNOSTICS WITHIN THIS CODE SEGMENT
        cpn % rle       = 0.1
        cpn % rpen      = 5
        cpn % rmaxfrac  = 5000000000000000.
        cpn % wmin      = 0.0 
        cpn % rbuoy     = 2./3.
        cpn % rdrag     = 3.0
        cpn % frac_drs  = 1.0
        cpn % bigc      = 0.7
        cpn % tcrit     = -45
        cpn % cldhgt_max= 40.e3
        cpn % auto_th0  = Nml%auto_th
        cpn % auto_rate = Nml%auto_rate
        cpn % atopevap  = Nml%atopevap
        cpn % do_ice    = Nml%do_ice
        cpn % do_ppen   = .false.
!!5 miz replaces do_edplume with mixing_assumption
        cpn % mixing_assumption = 1
!       cpn % do_edplume = .false.
!       cpn % do_edplume= .false.
!       cpn % do_micro  = .true.
        cpn % mp_choice = 0
        cpn % do_forcedlifting  = .true.
        cpn % wtwmin_ratio = Nml%wmin_ratio*Nml%wmin_ratio
! Values for cpn for the following variables are not actually used in the donner_lite routine but 
! they should be initialized in order to pass debug tests where NaNs are trapped.
        cpn % rad_crit        = 14
        cpn % deltaqc0        = 0.0005
        cpn % emfrac_max      = 1
        cpn % wrel_min        = 1
        cpn % Nl_land         = 300000000
        cpn % Nl_ocean        = 100000000
        cpn % r_thresh        = 1.2e-05
        cpn % qi_thresh       = 0.0001
        cpn % peff            = 1
        cpn % rh0             = 0.8
        cpn % cfrac           = 0.05
        cpn % hcevap          = 0.8
        cpn % weffect         = 0.5
        cpn % t00             = 295

        cp%maxcldfrac =  cpn%rmaxfrac

        if (ntr>0) then
           allocate(cpn%wetdep(ntr))
           call don_d_copy_wetdep_miz (Initialized%wetdep(:), &
                                       cpn%wetdep(:), &
                                       size(Initialized%wetdep(:)) )
        endif

        rkm = 2.*alpp*frac0; scaleh = 1000.; wrel = 0.5
        if(ac % plcl .gt. sd % pinv)then
           krel    = sd % kinv
        else
           krel    = ac % klcl
        endif
        cbmf  =sd%rho(krel-1)*(Param%cloud_base_radius**2)*wrel
        call cumulus_plume_k     &
               (cpn, sd, ac, cp, rkm, cbmf, wrel, scaleh, Uw_p, error, ermesg)
        call cumulus_tend_k   &
                (cpn, sd, Uw_p, cp, ct, .true.)

        if (ntr>0) then
           deallocate(cpn%wetdep)
        end if

        if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
           exit_flag_c = .true.
           return
        end if

        rcl_miz  (:)=sqrt(cp%ufrc(:))
        dpf_miz  (:)=(ct%qldiv(:)+ct%qidiv(:))/  &
                                          (Param%cloud_base_radius**2)
        dpf_warm (:)=(ct%qldiv(:))/  &
                                          (Param%cloud_base_radius**2)
        dpf_cold (:)=(ct%qidiv(:))/  &
                                          (Param%cloud_base_radius**2)

!BUGFIX 10/27/07
!       qlw_miz  (:)=cp%qlu(:) + cp%qiu(:)
        qlw_miz  (:)=cp%qlu(1:) + cp%qiu(1:)
        dfr_miz  (:)=0. !ct%qidiv (:)/(Param%cloud_base_radius**2)
!BUGFIX 10/27/07
!       flux_miz (:)=cp%umf
        flux_miz (:)=cp%umf(1:)
        efchr_miz(:)=ct%thcten(:)/(Param%cloud_base_radius**2)
        emfhr_miz(:)=ct%qvdiv (:)/(Param%cloud_base_radius**2)
!++++yim
        etfhr_miz   = ct%trten(:,:)/(Param%cloud_base_radius**2)
        dpftr_miz   = ct%trwet(:,:)/(Param%cloud_base_radius**2)

        ptma_miz (kou) = cp%ps(cp%ltop-1)
        dint_miz       = 0.
        cu_miz         = -(ct%conint+ct%freint)/  &
                               (Param%cloud_base_radius**2)*86400.
        cell_precip_miz=(ct%rain+ct%snow)/  &
                               (Param%cloud_base_radius**2)*86400.
        apt_miz        = rcl_miz(cp%ltop-1)/rcl_miz(krel-1)

        if (cu_miz == 0.0 .or. cell_precip_miz == 0.0) then
          exit_flag_c = .true.
          return
        end if

       summel = 0.
        do k=1,nlev_lsm
          dp = phalf_c(k) - phalf_c(k+1)
          summel = summel + (-1.0*dpf_cold(k))*dp/Param%grav
        end do
        if (debug_ijt) then
           write (diag_unit, '(a, f19.10)')    &
          'in mulsub: summel= ', summel
        endif
          
      if (pb > pmelt_lsm) then
          dmela = - ((summel*cell_precip_miz/cu_miz)*  &
                                Param%grav/(pmelt_lsm - pb))*8.64e07
          if (debug_ijt) then
           write (diag_unit, '(a, 3f19.10)')    &
                      'in mulsub: dmela, pmelt_lsm, pb= ', dmela , &
                                                pmelt_lsm, pb
          endif
            
        call don_u_map_hires_i_to_lores_c_k   &
           (nlev_lsm, dmela, pb, pmelt_lsm, phalf_c, cell_melt, ermesg, error)

          if (debug_ijt) then
            do k=1,nlev_lsm
              if (cell_melt(k) /= 0.0) then
                write (diag_unit, '(a, i4,  f19.10)')    &
                      'in mulsub: k, cell_melt= ', k, cell_melt(k)
              endif
            end do
          endif
    else
        cell_melt      = 0.
     endif

        where (qlw_miz(:) == 0.0)
          cfracice_miz   = 0.
        elsewhere
!BUGFIX 10/27/07
!         cfracice_miz(:) = cp%qiu(:)/qlw_miz(:)
          cfracice_miz(:) = cp%qiu(1:)/qlw_miz(:)
        end where
        ncc_kou_miz    = cp%ltop + 1

        cell_freeze= dfr_miz

        if (Nml%do_donner_lscloud) then
          if (cu_miz > 0.0) then
            evap_rate  =-dpf_miz*(1. - (cell_precip_miz/cu_miz))
          else
            evap_rate = 0.0
          endif
        else
           ecd        =ct%qlten/(Param%cloud_base_radius**2)
           ece        =ct%qiten/(Param%cloud_base_radius**2)
           evap_rate  =ct%qaten/(Param%cloud_base_radius**2)
           meso_melt  =ct%tten /(Param%cloud_base_radius**2)
           meso_freeze=ct%qvten/(Param%cloud_base_radius**2) 
        end if

!100 DEDFINE h1_liq, h1_ice
        h1_liq = -dpf_warm
        h1_ice = -dpf_cold
        h1_2       = efchr_miz
        q1         = emfhr_miz
        qtr        = 0.

!end: testing unified plume

    if (Nml%do_ensemble_diagnostics) then
!----------------------------------------------------------------------
!    save "Don_cem" diagnostics for this ensemble member.
!----------------------------------------------------------------------
       Don_cem%cell_precip(i,j,kou) = cell_precip_miz
       Don_cem%pb(i,j,kou) = pb
       Don_cem%ptma(i,j,kou) = ptma_miz(kou)
!  reverse index order
       do k=1,nlev_lsm
         Don_cem%h1(i,j,k,kou) = h1_liq(nlev_lsm-k + 1) +  &
                             h1_ice(nlev_lsm-k + 1)
         Don_cem%qlw(i,j,k,kou) = qlw_miz(nlev_lsm-k+1)
         Don_cem%cfracice(i,j,k,kou) = cfracice_miz(nlev_lsm-k+1)
         Don_cem%wv(i,j,k,kou) = cp%wu(nlev_lsm-k+1)
         Don_cem%rcl(i,j,k,kou) = rcl_miz(nlev_lsm-k+1)
       end do
    endif

!!!  BUGFIX: applies only when allow_mesoscale_circulation is set to
!!!         .false.; in such cases, bug turned off cell melting.
! Will change answers in runs with allow_mesoscale_circulation = .false.
        if (lmeso) then
         if ((pb - ptma_miz(kou)) < Param%pdeep_mc)  then
           lmeso = .false.
         endif
        endif

!----------------------------------------------------------------------
!    if in diagnostics column, output the cloud base (pb) and cloud top
!    (ptma) pressures, and the mesoscale circulation logical variable
!    (lmeso).
!----------------------------------------------------------------------
         if (debug_ijt) then
           write (diag_unit, '(a, 2f19.10,1l4)')    &
                'in mulsub: PB,PT, lmeso= ', pb, ptma_miz(kou), lmeso
        endif

         if (cu_miz /= 0.0) then
          precip_frac = cell_precip_miz/cu_miz
         else
           precip_frac = 0.
        endif
        ci_ice_cond = 0.
        do kk=1,nlev_lsm
          dp = phalf_c(kk) - phalf_c(kk+1)
          ci_ice_cond = ci_ice_cond + h1_ice(kk)*dp
!RSHfix for "s" release:  
!     replace the above line with the following; also comment out line
!     noted below. This fix will eliminate the occurrence of roundoff
!     "snow" falling from donner (~10e-22, + and -) that results from 
!     difference in this calc and that of "summel" above 
!NOTE THAT THE SAME CHANGE NEEDS TO BE MADE IN THE DONNER-FULL CODE.>>>
!         ci_ice_cond = ci_ice_cond + h1_ice(kk)*dp/Param%grav
        end do
       if (pmelt_lsm < pb) then
         melting_in_cloud = .true.
      else
        melting_in_cloud = .false.
     endif
!RSHfix  for "s" release -- comment out this line:
        ci_ice_cond = ci_ice_cond/(Param%grav)
        if (ci_ice_cond /= 0.0) then
          if (melting_in_cloud) then
         precip_melt_frac = summel/ci_ice_cond
           else
               precip_melt_frac = 0.
           endif
       else
         precip_melt_frac = 0.
       endif
       if (debug_ijt) then
         write (diag_unit, '(a, 3e20.12)')  &
            'in mulsub: h1_ice intg, summel, precip_melt_frac', &
                       ci_ice_cond, summel*cell_precip_miz/cu_miz, &
                                                    precip_melt_frac
          endif
       ci_liq_cond = 0.
       do kk=1,nlev_lsm
         dp = phalf_c(kk) - phalf_c(kk+1)
         ci_liq_cond = ci_liq_cond + h1_liq(kk)*dp
       end do
        ci_liq_cond = ci_liq_cond/(Param%grav)

!---------------------------------------------------------------------
!    if this member of the ensemble supports a mesoscale circulation,
!    call mesub to obtain various terms related to moving condensate
!    from the convective tower into the mesoscale anvil for this member.
!---------------------------------------------------------------------
        if (lmeso) then
          call don_cm_mesub_miz     &
               (Nml, pfull_c,nlev_lsm, me, diag_unit, debug_ijt, Param, cu_miz,   &
                ci_liq_cond, ci_ice_cond, pmelt_lsm, &
                cell_precip_miz, dint_miz, plzb_c, pb, ptma_miz(kou), &
                temp_c, phalf_c,     ca_liq, ca_ice, ecd, ecd_liq, &
                ecd_ice, ecei_liq, ece, ece_liq, ece_ice, meso_freeze, &
                meso_melt, ermesg, error)
          if (error /= 0 ) return
        else
          ca_liq = 0.
          ca_ice = 0.
          meso_freeze = 0.
          meso_melt   = 0.
        endif

        if (ci_liq_cond /= 0.0) then
          frz_frac = dint_miz/ci_liq_cond
          frz_frac_non_precip = dint_miz*(1.-precip_frac)/ci_liq_cond
        else
          frz_frac_non_precip = 0.
          frz_frac = 0.
        endif

        if (debug_ijt) then
          write (diag_unit, '(a, 3e20.12)')  &
           'in mulsub pre anvil_cond_frz: h1_liq intg, dint_miz, frz_frac',&
                         ci_liq_cond, dint_miz, frz_frac_non_precip
        endif
        if (debug_ijt) then
          write (diag_unit, '(a, 1e20.12)')  &
                                 'in mulsub : frz_frac', &
                                   frz_frac
        endif

        if (cu_miz /= 0.0) then
          meso_frac = (ca_liq + ca_ice)/cu_miz
        else
          meso_frac = 0.
        endif

        if (lmeso) then
          bak = 0.
          do kk=1,nlev_lsm
            dp = phalf_c(kk) - phalf_c(kk+1)
            bak = bak + meso_freeze(kk)*dp
          end do
          bak = bak/(Param%grav)
          bak = bak/(Param%seconds_per_day*1.0e3)
          if (debug_ijt) then
            write (diag_unit, '(a, 3e20.12)')  &
             'in mulsub: column meso_freeze', bak
          endif
          if (bak > 0.0) meso_frz_intg = .true.
          if (ci_liq_cond /= 0.0) then
            meso_frz_frac = bak/ci_liq_cond
          else
            meso_frz_frac = 0.
          endif
          if (debug_ijt) then
            write (diag_unit, '(a, i4, 2e20.12)')  &
             'in mulsub : kou, frz_frac, meso_frz_frac', &
                                 kou, frz_frac_non_precip, meso_frz_frac
          endif
        else
          meso_frz_intg = .false.
          meso_frz_frac = 0.
       endif

     if (meso_frz_frac == 0. .and.  .not. melting_in_cloud) then
         if (meso_frac < frz_frac_non_precip .and. meso_frac > 0.) then
             do k=1,nlev_lsm
          cell_freeze(k) = cell_freeze(k)*meso_frac/frz_frac_non_precip
             end do
             dint_miz = dint_miz *meso_frac/frz_frac_non_precip
             frz_frac = frz_frac*meso_frac/frz_frac_non_precip
            frz_frac_non_precip = meso_frac
           endif
        endif
 
         if (debug_ijt) then
           write (diag_unit, '(a, i4, 3e20.12)')  &
                'in mulsub : kou, ADJUSTEDfrz_frac, meso_frz_frac, &
                   &precip_melt_frac', &
                         kou, frz_frac_non_precip, meso_frz_frac,  &
                                                      precip_melt_frac
        endif
        
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 3f19.10)')    &
                      'in mulsub: meso_frac, precip_frac,frz_frac:', &
                        kou, meso_frac, precip_frac, frz_frac_non_precip
          write (diag_unit, '(a, i4, 4f19.10)')    &
                      'in mulsub: cu, ca, cell_precip, dint_miz   :', &
                          kou, cu_miz, ca_liq + ca_ice,   &
                             cell_precip_miz, dint_miz*Param%seconds_per_day
        endif

         if (lmeso) then
           summel = 0.
         endif
       if (debug_ijt) then
          write (diag_unit, '(a, 4f19.10)')    &
                      'in mulsub: pmelt_lsm, pb, summel   :', &
                          pmelt_lsm, pb, summel*cell_precip_miz/cu_miz
       endif
        
         do k=1,ncc_kou_miz
!----------------------------------------------------------------------
!    define the factor needed to normalize each ensemble member's con-
!    tribution by the cloud base area of ensemble member #1. wt_factor
!    is the cloud area at level k for ensemble member kou, relative to
!    the cloud area at cloud base (k=1) for ensemble member #1.
!-----------------------------------------------------------------------
        wt_factor = Param%arat(kou)*(rcl_miz(k)/rcl_miz(krel-1))**2
        
!----------------------------------------------------------------------
!    add this ensemble member's appropriately weighted contribution to
!    the ensemble-total cloud area (alp), condensed ice (cuql), condensed
!    liquid (cuqli), cell upward mass flux (ucemh), cell detrained mass 
!    flux (detmfh), condensation rate (rlsm), vertical moisture flux 
!    convergence (emsm) and vertical tracer flux convergence (etsm). the
!    weighting factor area_ratio*(rcl(k)/rcl(1))**2 allows the contrib-
!    utions from each member to be added by normalizing each member's 
!    contribution by the cloud base area of ensemble member #1.
!    NOTE: several of the arrays already have some of the normalizing
!    factors already included and so here need only to be multiplied by 
!    a portion of wt_factor.
!----------------------------------------------------------------------
        alp_miz  (k) = alp_miz  (k) + wt_factor                      
        cuql_miz (k) = cuql_miz (k) + wt_factor*  &
                             (cfracice_miz(k)*qlw_miz(k))
        cuqli_miz(k) = cuqli_miz(k) + wt_factor*  &
                             ((1.0 - cfracice_miz(k))*qlw_miz(k))
        ucemh_miz(k) = ucemh_miz(k) + Param%arat(kou)*flux_miz(k)/ &
                                                    (rcl_miz(krel-1)**2)
        if (k < ncc_kou_miz) then
          if (flux_miz(k+1) < flux_miz(k)) then
            detmfh_miz(k) = detmfh_miz(k) + Param%arat(kou)*   &
                        ((flux_miz(k)-flux_miz(k+1))/(rcl_miz(krel-1)**2))
          endif
        endif
        rlsm_miz(k)   = rlsm_miz(k)   - Param%arat(kou)*dpf_miz (k) 
        emsm_miz(k)   = emsm_miz(k)   + Param%arat(kou)*emfhr_miz(k)
        qvfm_miz(k)   = Param%arat(kou)*(dpf_miz (k) + emfhr_miz(k))
        qvfm_tot(k)   = qvfm_tot(k) + qvfm_miz(k)

        etsm_miz(k,:) = etsm_miz(k,:) + Param%arat(kou)*etfhr_miz(k,:)
     enddo

     do n=1,ntr
     do k=1,ncc_kou_miz
!        etsm_miz(k,n) = etsm_miz(k,n) + (Param%arat(kou)*etfhr_miz(k,n))
        qtren(k,n)    = qtren(k,n) + (Param%arat(kou)*etfhr_miz(k,n))
        ensmbl_wetc(k,n) = ensmbl_wetc(k,n) +  &
                                   (Param%arat(kou)*dpftr_miz(k,n))
     end do
     end do

!--------------------------------------------------------------------
!    if a mesoscale circulation is present, add this member's cloud-
!    base_area-weighted contribution of condensate transferred to the 
!    anvil (ensmbl_anvil_cond) and cloud top cloud fraction 
!    (ensmbl_cld_top_area) to the arrays accumulating the ensemble sums.
!--------------------------------------------------------------------
      if (lmeso) then
        if (meso_frac /= 0.0) then
          local_frz_frac = (frz_frac_non_precip + meso_frz_frac)/  &
                                                              meso_frac
        else
          local_frz_frac = 0.0
        endif
        ensmbl_anvil_cond_liq   = ensmbl_anvil_cond_liq   +   &
                          Param%arat(kou)*ca_liq*(1.-local_frz_frac)
        ensmbl_anvil_cond_liq_frz   = ensmbl_anvil_cond_liq_frz   +   &
                               Param%arat(kou)*ca_liq*local_frz_frac
        ensmbl_anvil_cond_ice   = ensmbl_anvil_cond_ice   +  &
                                            Param%arat(kou)*ca_ice
        ensmbl_cld_top_area = ensmbl_cld_top_area +   &
                                               Param%arat(kou)*apt_miz
      endif

!--------------------------------------------------------------------
!    add this ensemble member's weighted contribution to the total 
!    precipitation (ensmbl_precip) and condensation (ensmbl_cond). 
!--------------------------------------------------------------------
      ensmbl_precip = ensmbl_precip + Param%arat(kou)*cell_precip_miz
      ensmbl_cond   = ensmbl_cond   + Param%arat(kou)*cu_miz

!---------------------------------------------------------------------

!!$        call don_d_add_to_ensmbl_sum_lores_k    &
!!$             (nlev_lsm, ntr, diag_unit, debug_ijt, lmeso, Param,   &
!!$              Param%arat(kou), dint, cell_freeze, cell_melt, temp_c,   &
!!$              h1_2, ecd, ece, evap_rate, q1, h1, pfull_c, meso_melt, &
!!$              meso_freeze, phalf_c, qtr, ensmbl_melt, ensmbl_freeze, &
!!$              enctf, encmf, enev, disg, disb, disc, ecds, eces, disd, &
!!$              qtren, ermesg)
!!$        if (trim(ermesg) /= ' ') return



      do k=1,nlev_lsm       

!---------------------------------------------------------------------
!    define the moisture forcing term (sum of condensation h1 and 
!    vertical flux convergence q1) on the large-scale grid. convert to
!    units of g(h20) per kg(air) per day, requiring multiplication by
!    1.0e3 g(h2o) /kg(h2o) times SECONDS_PER_DAY. add this member's 
!    contribution to the sum over the ensemble (encmf). 
!----------------------------------------------------------------------
        cmf(k) = (-(h1_liq(k) + h1_ice(k)) + q1(k))*  &
                                         (Param%SECONDS_PER_DAY*1.0e03)
        encmf(k) = encmf(k) + Param%arat(kou)*cmf(k)

!----------------------------------------------------------------------
!    define the condensation term in the temperature equation on the 
!    large-scale grid (rlh), using the latent heat of vaporization when 
!    the ambient temperature is above freezing, and the latent heat of 
!    sublimation when ice may be present. add this member's contribution
!    to the sum over the ensemble (disc).
!----------------------------------------------------------------------
        if (.not. melting_in_cloud) then
          disz(k) = disz(k) + Param%arat(kou)*h1_liq(k)*   &
                             Param%seconds_per_day*frz_frac*precip_frac
        else
          disz_remelt(k) = disz_remelt(k) + Param%arat(kou)*h1_liq(k)* &
                             Param%seconds_per_day*frz_frac*precip_frac
        endif
 
        disp_liq(k) = disp_liq(k) + Param%arat(kou)*h1_liq(k)*   &
                      Param%seconds_per_day*(1.0-frz_frac)* precip_frac
        disze1(k) = disze1(k) + Param%arat(kou)*h1_liq(k)*  &
                        Param%seconds_per_day*Param%hls*frz_frac* &  
                                       (1.-precip_frac)/Param%cp_air
        disze2(k) = disze2(k) + Param%arat(kou)*h1_liq(k)*  &
                         Param%seconds_per_day*Param%hlv*  &
                          (1.0-frz_frac)*(1.-precip_frac)/Param%cp_air
        disze3(k) = disze3(k) - Param%arat(kou)*h1_ice(k)*  &
                        Param%seconds_per_day*Param%hls*  &
                                         (1.-precip_frac)/Param%cp_air
        disp_ice(k) = disp_ice(k) + Param%arat(kou)*h1_ice(k)*  &
                         Param%seconds_per_day*(1.0-precip_melt_frac)* &
                                                         precip_frac   
        disp_melted(k) = disp_melted(k) + Param%arat(kou)*h1_ice(k)*  &
                          Param%seconds_per_day* precip_melt_frac * &
                                                           precip_frac
        disc_liq(k) = disc_liq(k) + Param%arat(kou)*h1_liq(k)* &
                        Param%seconds_per_day*Param%hlv/ Param%cp_air
!  if no melting, the frozen liquid stays frozen and carries out hls;
!  if melting, the frozen liquid melts and carries out hlv.

        if (.not. melting_in_cloud) then
!55 should be if meso and cell freezing both 0.0; if cell freezing
!55  non-zero, then all meso entrained liq will have frozen
!         if (meso_frz_intg <= 0. .and. frz_frac == 0.) then
!         if (.not. meso_frz_intg       .and. frz_frac == 0.) then
          if (.not. meso_frz_intg .and. frz_frac_non_precip == 0.) then
            dism_liq(k) = dism_liq(k) + Param%arat(kou)*h1_liq(k)* &
                                      meso_frac*Param%seconds_per_day                          
          else
!! NOT ALL LIQUID FREEZES; only frz_frac + meso_frz_frac
!53  when no melting and freezing, then all condensate is frozen when
!53   it precipitates 
            dism_liq_frz(k) = dism_liq_frz(k) + Param%arat(kou)*  &
                               h1_liq(k)*meso_frac*Param%seconds_per_day
          endif
          dism_ice(k) = dism_ice(k) + Param%arat(kou)*h1_ice(k)*  &
                                     meso_frac*Param%seconds_per_day
        else
          dism_liq_remelt(k) = dism_liq_remelt(k) + Param%arat(kou)*  &
                               h1_liq(k)*meso_frac*Param%seconds_per_day
          dism_ice_melted(k) = dism_ice_melted(k) +   &
                                  Param%arat(kou)*h1_ice(k)*meso_frac* &
                                                  Param%seconds_per_day
        endif
        disc_ice(k) = disc_ice(k) + Param%arat(kou)*h1_ice(k)* &
                   Param%seconds_per_day*Param%hls/ Param%cp_air

         if (debug_ijt) then
              write (diag_unit, '(a, i4, 3e20.12)')  &
                             'in mulsub: precip profiles', &
                                k, disp_liq(k)*Param%hlv/Param%cp_air,&
                        Param%hls/Param%cp_air*disp_ice(k), &
                         Param%hls/Param%cp_air*disz(k)
              write (diag_unit, '(a, i4, 2e20.12)')  &
                            'in mulsub: remelt, melt precip profiles', &
                    k, Param%hlv/Param%cp_air*disz_remelt(k), &
                                  Param%hlv/Param%cp_air*disp_melted(k)
              write (diag_unit, '(a, i4, 3e20.12)')  &
                             'in mulsub: evap   profiles', &
                              k, disze1(k)                       , &
                           disze2(k)                       , &
                         -disze3(k)                       
              write (diag_unit, '(a, i4, 2e20.12)')  &
                              'in mulsub: cd     profiles', &
                                k, disc_liq(k), disc_ice(k)
       endif


!--------------------------------------------------------------------
!    add this member's weighted contribution to the ensemble's temper-
!    ature flux convergence (disb), the ensemble's water vapor flux 
!    convergence (disd) and the ensemble's entropy flux convergence 
!    (enctf). convert the rates to units of per day, and for disd from
!    kg(h2o) per kg(air) to g(h2o) per kg(air).
!--------------------------------------------------------------------
        disb(k) = disb(k) + Param%arat(kou)*(h1_2(k)*  &
                                              Param%SECONDS_PER_DAY)
        disd(k) = disd(k) + Param%arat(kou)*(q1(k)*   &
                                       (Param%SECONDS_PER_DAY*1.0e3))
        disv(k) = disv(k) + Param%arat(kou)*((h1_liq(k) + h1_ice(k))*  &
                             (Param%SECONDS_PER_DAY*1.0e3))
!   change enctf to reflect need for both ice and liq cd in layer of
!   tfre
       enctf(k) = enctf(k) + Param%arat(kou)*    &
                      (h1_2(k)*Param%SECONDS_PER_DAY + &
                        (h1_liq(k)*Param%hlv + h1_ice(k)*Param%hls)*  &
                           Param%seconds_per_day/ Param%cp_air )

!--------------------------------------------------------------------
!    if a mesoscale circulation exists, add this member's contribution
!    to the mesoscale condensate's evaporation associated with convect-
!    ive downdrafts (ecds) and that associated with evaporation into 
!    the environment (eces). if there has been no freezing associated
!    with the mesoscale condensate, define the condensation term for
!    the temperature equation using the latent heat of vaporization
!    (disg). if there has been freezing, then the convective downdraft 
!    heating uses the latent heat of vaporization, whereas the entrain-
!    ment evaporation is of ice and so uses the latent heat of 
!    sublimation.
!--------------------------------------------------------------------
        if (lmeso) then
          ecds_liq(k) = ecds_liq(k) + Param%arat(kou)*ecd_liq(k)
          ecds_ice(k) = ecds_ice(k) + Param%arat(kou)*ecd_ice(k)
          eces_liq(k) = eces_liq(k) + Param%arat(kou)*ece_liq(k)
          eces_ice(k) = eces_ice(k) + Param%arat(kou)*ece_ice(k)
          disg_ice(k) = disg_ice(k) - Param%arat(kou)*((ecd_ice(k) + &
                                         ece_ice(k))*  &
                                Param%hls/(Param%CP_AIR*1000.))
!         if (dint_miz == 0.) then
!           disg_liq(k) = disg_liq(k) - Param%arat(kou)*((ecd_liq(k) + &
!                                      ece_liq(k))*  &
!                               Param%hlv/(Param%CP_AIR*1000.))
!         else
!          if (melting_in_cloud) then
!            disg_liq(k) = disg_liq(k) - Param%arat(kou)*  &
!                       ((ece_liq(k)                          )*  &
!                             Param%HLS/(Param%CP_AIR*1000.))
!          else
!            if (.not. meso_frz_intg       ) then
!               disg_liq(k) = disg_liq(k) - Param%arat(kou)*  &
!                          (((ece_liq(k)*Param%hlv)  )/ &
!                                       (Param%CP_AIR*1000.))
!            else
!               disg_liq(k) = disg_liq(k) - Param%arat(kou)*  &
!                      ((            ece_liq(k)*Param%hls  )/ &
!                                       (Param%CP_AIR*1000.))
!            endif
!          endif

!          disg_liq(k) = disg_liq(k) - Param%arat(kou)*  &
!                              (ecd_liq(k)*Param%HLV/  &
!                      (Param%CP_AIR*1000.))
!        endif
          if (dint_miz /= 0. .and. &
               (melting_in_cloud .or. meso_frz_intg) ) then
            disg_liq(k) = disg_liq(k) - Param%arat(kou)*  &
                   ((ecd_liq(k)*Param%hlv + &
                     ece_liq(k)*Param%hls)/(Param%CP_AIR*1000.))
          else
            disg_liq(k) = disg_liq(k) - Param%arat(kou)*((ecd_liq(k) + &
                                       ece_liq(k))*  &
                                Param%hlv/(Param%CP_AIR*1000.))
          endif
      endif  ! (lmeso)

!---------------------------------------------------------------------
!    add this member's cloud water evaporation rate to the sum over 
!    the ensemble (enev).
!---------------------------------------------------------------------

        if (Nml%do_donner_lscloud) then
           enev(k) = enev(k) + Param%arat(kou)*evap_rate(k)
        else
           enev(k) = enev(k) + Param%arat(kou)*evap_rate(k) !miz for qa detrainment
           ecds_liq(k) = ecds_liq(k) + Param%arat(kou)*ecd_liq(k) !miz: add temporarily for ql detrainment
           ecds_ice(k) = ecds_ice(k) + Param%arat(kou)*ecd_ice(k) !miz: add temporarily for ql detrainment
           eces_liq(k) = eces_liq(k) + Param%arat(kou)*ece_liq(k) !miz: add temporarily for qi detrainment
           eces_ice(k) = eces_ice(k) + Param%arat(kou)*ece_ice(k) !miz: add temporarily for qi detrainment
        end if

!--------------------------------------------------------------------
!    if a mesoscale circulation exists, add the appropriately-weighted
!    anvil freezing and melting terms to the arrays accumulating their 
!    sums over the ensemble (ensmbl_melt, ensmbl_freeze). if in a diag-
!    nostic column, output the anvil (meso_freeze) and ensemble-sum
!    (ensmbl_freeze) freezing profiles.
!--------------------------------------------------------------------
        if (lmeso) then
          ensmbl_melt_meso(k) = ensmbl_melt_meso(k) - Param%arat(kou)*&
                         meso_melt(k)
          ensmbl_freeze_meso(k) = ensmbl_freeze_meso(k) +  &
                                  Param%arat(kou)*meso_freeze(k)
         if (debug_ijt) then
           if (meso_freeze(k) /= 0.0) then
             write (diag_unit, '(a, i4, 2e20.12)')  &
                              'in mulsub: jk,fres,fre= ',   &
                              k, ensmbl_freeze_meso(k), meso_freeze(k)
            endif
         endif

        endif

!--------------------------------------------------------------------
!    add the appropriately-weighted convective cell freezing and 
!    melting terms to the arrays accumulating vertical profiles of 
!    total cloud melting (ensmbl_melt) and freezing (ensmbl_freeze) 
!    over the entire ensemble.  if in diagnostic column, output the 
!    convective cell (cell_freeze) and accumulated (ensmbl_freeze) 
!    freezing profiles.
!--------------------------------------------------------------------
        ensmbl_freeze(k) = ensmbl_freeze(k) +    &
                                        Param%arat(kou)*cell_freeze(k)
        ensmbl_melt(k) = ensmbl_melt(k) - Param%arat(kou)*cell_melt(k)
        if (debug_ijt) then
          if (cell_freeze(k) /= 0.0) then
            write (diag_unit, '(a, i4, 2e20.12)')  &
                     'in mulsub: jk,fres,frea= ',    &
                                    k, ensmbl_freeze(k), cell_freeze(k)
         endif
       endif
      end do

!--------------------------------------------------------------------
!    save the cloud top (ptma) pressures, the total condensation (cuto),
!    total precpitation (preto) and cloud top index (ncca) from this !
!    ensemble member.
!--------------------------------------------------------------------
        if (meso_frz_intg) meso_frz_intg_sum = .true.
        ncca(kou)  = ncc_kou_miz

      end do   ! (kou loop over ensemble members)

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! 31   CONTINUE
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!----------------------------------------------------------------------
!    define ensemble cloud top pressure (pt_ens) to be the cloud top of 
!    the most penetrative ensemble member. this is frequently, but not 
!    always, the ensemble member with the lowest entrainment rate. 
!    cloud base pressure (pb) is the same for all ensemble members. 
!    define the cloud top index(ncc_ens)  as the highest of any ensemble 
!    member.
!----------------------------------------------------------------------
      pt_ens  = MINVAL (ptma_miz)
      ncc_ens = MAXVAL (ncca)

!----------------------------------------------------------------------
!    divide the ensemble mean ice and liquid condensate terms by the 
!    total cloud area to define the average cloud water and cloud ice 
!    concentrations within the cloudy area, as opposed to averaged over 
!    the entire grid box.
!----------------------------------------------------------------------
      do k=1,ncc_ens
        if (alp_miz(k) > 0.) then
          cuql_miz (k) = cuql_miz (k)/alp_miz(k)
          cuqli_miz(k) = cuqli_miz(k)/alp_miz(k)
        endif
      end do

!---------------------------------------------------------------------
!    define the cloudtop anvil area (ampta1), assumed to be mesofactor
!    (default = 5) times larger than the sum of the cloud top areas of 
!    the ensemble members, as in Leary and Houze (1980), 
!---------------------------------------------------------------------
      ampta1 = Nml%mesofactor*ensmbl_cld_top_area

!---------------------------------------------------------------------
!    if there is no precipitation production in this column, set the 
!    inverse of the max cloud area at any layer in the column to be 0.0.
!---------------------------------------------------------------------
      if (ensmbl_precip == 0.0) then
        amax      = 0.0
      else

!---------------------------------------------------------------------
!    if there is precip in the column, determine the maximum convective 
!    cell area at any level in the column (al). the total normalized 
!    cloud area in the column (cell area + mesoscale area) cannot be 
!    greater than 1.0. this constraint imposes a limit on the cloud area
!    at cloud base (amax). this limit will be imposed in subroutine
!    determine_cloud_area. see "a bounds notes" (7/6/97).
!---------------------------------------------------------------------
        al = MAXVAL (alp_miz)
        amax = 1./(al + ampta1)
      endif

!---------------------------------------------------------------------
!    if in diagnostics column, output the total ensemble condensation,
!    (ensmbl_cond), precipitation (ensmbl_precip), and condensate 
!    transferred into the anvil (ensmbl_anvil_cond). also output 
!    surface pressure (phalf_c(1)), ensemble cloud base nd cloud top 
!    pressures (pb, pt_ens), the flag indicating if a mesoscale circul-
!   ation is present in the grid column (lmeso), and the cloud top anvil
!---------------------------------------------------------------------
       if (debug_ijt) then
         write (diag_unit, '(a, e20.12, a, e20.12)')  &
                       'in mulsub: CUTOT=', ensmbl_cond, ' PRETOT=', &
                                      ensmbl_precip
        write (diag_unit, '(a, 4e20.12)') &
               'in mulsub: CATOT, (sum,liq, frzliq, ice)=', &
                 ensmbl_anvil_cond_liq + ensmbl_anvil_cond_liq_frz +   &
                                               ensmbl_anvil_cond_ice, &
                                     ensmbl_anvil_cond_liq, &
                                      ensmbl_anvil_cond_liq_frz, &
                                     ensmbl_anvil_cond_ice
        write (diag_unit, '(a, 3f19.10, 1l4)')  &
                       'in mulsub: ps,pb,pt,lmeso= ',   &
                                       phalf_c(1), pb, pt_ens, lmeso
        write (diag_unit, '(a, e20.12)')  &
                                 'in mulsub: ampt= ',ampta1
     endif



!!$      ptt = pt_ens + Param%dp_of_cloud_model
!!$!--------------------------------------------------------------------
!!$!    call define_ensemble_profiles to produce vertical profiles 
!!$!    representing the ensemble-total cloud area (ensmbl_cloud_area), 
!!$!    cloud liquid (cuql_v), cloud ice (cuq), mass flux(uceml) and
!!$!    detrained mass flux (detmfl).
!!$!--------------------------------------------------------------------
!!$      call don_d_def_ensemble_profs_k    &
!!$           (nlev_lsm, nlev_hires, ncc_ens, diag_unit, debug_ijt, ptt, &
!!$            cld_press, alp, detmfh, ucemh, cuql, cuqli, phalf_c,  &
!!$            ensmbl_cloud_area, cuql_v, cuq, detmfl, uceml, ermesg)
!!$      if (trim(ermesg) /= ' ') return

      
      ensmbl_cloud_area = alp_miz
      cuq               = cuql_miz
      cuql_v            = cuqli_miz
      uceml             = ucemh_miz
      detmfl            = detmfh_miz

      rlsm              =rlsm_miz
      emsm              =emsm_miz
      etsm              =etsm_miz
      cld_press         =sd%p


    end subroutine don_d_integ_cu_ensemble_miz


!######################################################################
!######################################################################


subroutine don_d_determine_cloud_area_miz            &
        (me, nlev_model, ntr, dt, nlev_parcel, diag_unit, debug_ijt,  &
          Param, Initialized,Nml, tracers, pfull, zfull, phalf, zhalf, &
          pblht, tkemiz, qstar, cush, cbmf, land, coldT, sd, Uw_p, &
          ac, max_depletion_rate, cape, dcape, amax, dise_v, disa_v,  &
          pfull_c, temp_c, mixing_ratio_c, env_t, env_r, parcel_t,  &
          parcel_r, cape_p, exit_flag, amos, a1, ermesg, error)

!---------------------------------------------------------------------
!    subroutine determine_cloud_area defines the convective cloud area
!    and so closes the donner_deep parameterization. The arrays 
!    Don_conv%a1 and Don_conv%amos are output by this routine.
!---------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_nml_type, &
                             donner_initialized_type
use conv_utilities_k_mod, only : sounding, adicloud, uw_params !miz

implicit none

!-----------------------------------------------------------------------
!++++yim
integer,                      intent(in)    :: me, nlev_model,nlev_parcel, ntr, diag_unit
real,                         intent(in)    :: dt
logical,                      intent(in)    :: debug_ijt
type(donner_param_type),      intent(in)    :: Param
type(donner_initialized_type),intent(in)    :: Initialized
type(donner_nml_type),        intent(in)    :: Nml      
real,                         intent(in)    :: max_depletion_rate,  &
                                               cape, dcape, amax
real, dimension(nlev_model),    intent(in)    :: dise_v, disa_v, pfull_c, temp_c, mixing_ratio_c 
real, dimension(nlev_model),  intent(in)    :: pfull, zfull
real, dimension(nlev_model+1),  intent(in)    :: phalf, zhalf !miz
real,                         intent(in)    :: pblht, tkemiz,  &
                                               qstar, cush, cbmf, land !miz
logical,                      intent(in)    :: coldT !miz
!++++yim
real, dimension(nlev_model,ntr),    intent(in)    :: tracers
type(sounding),               intent(inout) :: sd !miz
type(uw_params),               intent(inout) :: Uw_p !miz
type(adicloud),               intent(inout) :: ac !miz
real, dimension(nlev_parcel),    intent(in)    :: env_t, env_r, parcel_t, parcel_r, cape_p
logical,                      intent(inout) :: exit_flag
real,                         intent(out)   :: amos, a1
character(len=*),             intent(out)   :: ermesg
integer,                      intent(out)   :: error
 
real, dimension (nlev_model)  :: a1_vk              
real, dimension(nlev_parcel)   :: qli0_v, qli1_v, qt_v, qr_v, rl_v, ri_v
real                        :: qtest, tfint, disbar
integer                     :: k
!----------------------------------------------------------------------
!   local variables:
!
!         a1_vk
!         qli0      normalized component of cumulus condensate forcing
!         qli1      un-normalized component of condensate forcing
!         qt_v      temperature tendency due to deep convection on
!                   cape grid [ deg K / sec ]
!         qr_v      vapor mixing ratio tendency due to deep convection
!                   on cape grid [ kg(h2o) / ( kg(air) sec ]
!         rl_v      large-scale liquid mixing ratio
!         ri_v      large-scale ice mixing ratio 
!         qtest
!         tfint     column integral of moisture time tendency due to
!                   convection  [ mm / sec , or  kg / (m**2 sec ) ]
!         disbar    water vapor time tendency due to deep convection at 
!                   large-scale model interface levels
!                   [ kg(h2o) / ( kg(air) sec ) ]
!         nlev      number of layers in large-scale model
!         k         do-loop index

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

       if (Nml%do_hires_cape_for_closure) then
 
!---------------------------------------------------------------------
!    call map_lo_res_col_to_hi_res_col to interpolate moisture and
!    temperature forcings from large-scale model grid (dise_v, disa_v)
!    to the vertical grid used in the cape calculation (qr_v, qt_v). 
!--------------------------------------------------------------------
        call don_u_lo1d_to_hi1d_k   &
            (nlev_model, nlev_parcel, disa_v, pfull_c, cape_p, qt_v,  &
             ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
       if (error /= 0 ) return

       call don_u_lo1d_to_hi1d_k   &
         (nlev_model, nlev_parcel, dise_v, pfull_c, cape_p, qr_v, &
           ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
       if (error /= 0 ) return

    else ! (do_hires_cape_for_closure)
      qt_v=disa_v!miz
      qr_v=dise_v!miz
    endif ! (do_hires_cape_for_closure)

!--------------------------------------------------------------------
!    if in a diagnostic column, output the temperature and moisture 
!    forcings on both the cape grid (qt_v, qr_v) and the large-scale
!    model grid (disa_v, dise_v).
!--------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_parcel
          if (qr_v(k) /= 0.0 .or. qt_v(k) /= 0.0) then
            write (diag_unit, '(a, i4, e20.12, f20.14)')  &
                      'in cupar: k,qr,qt= ',k, qr_v(k), qt_v(k)
          endif
        end do
        do k=1,nlev_model
          if (dise_v(k) /= 0.0 .or. disa_v(k) /= 0.0) then
            write (diag_unit, '(a, i4, 2e20.12)')  &
                     'in cupar: k,dise,disa= ',k, dise_v(k), disa_v(k)
          endif
        end do
      endif

!--------------------------------------------------------------------
!   define condensate variables on the cape grid (qli0, qli1, rl_v, 
!   ri_v). these variables are not used in the current version of the
!   cumulus closure scheme implemented in subroutine cumulus_closure, 
!   so they are given values of 0.0.
!--------------------------------------------------------------------
      do k=1,nlev_parcel !miz nlev_hires
        qli0_v(k) = 0.
        qli1_v(k) = 0.
        rl_v(k)   = 0.
        ri_v(k)   = 0.
      end do

!--------------------------------------------------------------------
!    call subroutine cumulus_closure to determine cloud base cloud
!    fraction and so close the deep-cumulus parameterization.
!--------------------------------------------------------------------
    if (Nml%deep_closure .eq. 0) then
      call cu_clo_cumulus_closure_miz   &
           (nlev_model, nlev_parcel, ntr, dt, diag_unit, debug_ijt, &
            Initialized, Param, tracers, &
            dcape, pfull, zfull, phalf, zhalf, pblht, tkemiz, qstar, &
            cush, land, &
            coldT, sd, Uw_p, ac, Nml,        &!miz
            cape_p, qli0_v, qli1_v, qr_v, qt_v, env_r, ri_v, &
            rl_v, parcel_r, env_t, parcel_t, a1, ermesg, error)     
     else if (Nml%deep_closure .eq. 1) then
       call cu_clo_cjg    &
            (me, nlev_model, nlev_parcel, ntr, dt, diag_unit, &
             debug_ijt, Initialized, Param, Nml, amax, cape, dcape, cape_p, env_t, &
             env_r, qt_v, qr_v, pfull, zfull, phalf, zhalf, tracers, &
             land, pblht, tkemiz, qstar, coldT, sd, Uw_p, ac, &
             a1, ermesg, error )
     else
        call  cu_clo_miz   &
             (nlev_model, ntr, dt, Initialized, Param, tracers, &
              pfull, zfull, phalf, zhalf, pblht, tkemiz, qstar, cush, cbmf, land,  &
              coldT, sd, Uw_p, ac, Nml, env_r, env_t, a1, ermesg, error)
     endif

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!--------------------------------------------------------------------
!    calculate the vertical integral of normalized moisture forcing 
!    in the column (tfint) in units of kg (h2o) per m**2 per second, or
!    mm (h2o) per second.
!-------------------------------------------------------------------
      tfint = 0.0
      do k=2,nlev_model
        disbar = 0.5*(dise_v(k-1) + dise_v(k))
        tfint = tfint - disbar*(pfull_c(k-1) - pfull_c(k))
      end do
      tfint = tfint/Param%grav

!--------------------------------------------------------------------
!    restrict the cloud-base area fraction produced by subroutine
!    cumulus_closure to be no larger than the cloud base area that 
!    results in total grid box coverage at some higher level (amax). 
!--------------------------------------------------------------------
    if (Nml%deep_closure .eq. 0 .or. Nml%deep_closure .eq. 1) then
      a1 = MIN (amax, a1)
     else
      a1 = MIN (0.25, a1)
    end if

!---------------------------------------------------------------------
!    set the cloud-base area fraction to be 0.0 if there is no net
!    column integral of moisture forcing in the column. this is 
!    referred to as the moisture constraint. see "Moisture Constraint",
!    8/8/97. set the exit_flag to .true., turning off convection in
!    this column, output a message, and return to calling subprogram.
!---------------------------------------------------------------------
      if (tfint == 0.) then      
        a1 = 0.
        exit_flag      = .true.
        if (debug_ijt) then
          write (diag_unit, '(a)')  &
                 'convection turned off in column because of moist&
                  &ure constraint; cloud area being set to 0.0'
        endif
        return
      endif

!---------------------------------------------------------------------
!    if in a diagnostic column, output the column integral of the 
!    moisture forcing (tfint) and the fractional cloud area (a1) after
!    assuring that moisture forcing is present in the column.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12)')  &
                       'in cupar: tfint= ',tfint
        write (diag_unit, '(a, e20.12)')  &
                       'in cupar: a1_v = ',a1
      endif


!---------------------------------------------------------------------
!    restrict cloud fractional area by the moisture constraint. this
!    requirement limits the cloud area so that the moisture tendency 
!    due to the deep convection (tfint - which occurs only within the 
!    cloud fractional area) will not remove more vapor from the column 
!    than is available. here amos is the cloud area over which applic-
!    ation of the convective moisture tendency will result in total
!    vapor depletion in the column.
!---------------------------------------------------------------------
      amos = max_depletion_rate/tfint     
      if (a1 > amos)  then    
        a1 = max(amos, 0.)
      endif 

!---------------------------------------------------------------------
!    for any diagnostic columns in the window in which deep convection
!    was possible, output the column integral of the moisture forcing 
!    (tfint), the max cloud area allowed by the moisture constraint 
!    (amos) and the fractional cloud area after applying the moisture
!    constraint (a1).
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 3e20.12)')  &
                    'in cupar: tfint,amos,a1= ',  &
                                     tfint, amos, a1
      endif

!---------------------------------------------------------------------
!    verify that the current value of a1 will not produce negative
!    value of vapor mixing ratio at any level in the column when the
!    convective moisture tendency is applied. determine the large-scale
!    model mixing ratio for the current value of a1 (qtest). if qtest
!    is negative at any level for this value of a1, reset the value 
!    of a1, so that no negative mixing ratios will be produced.
!--------------------------------------------------------------------
      do k=1,nlev_model
        qtest = mixing_ratio_c(k) + a1*Nml%donner_deep_freq*dise_v(k)
        if (qtest < 0.) then
          a1_vk(k) = -mixing_ratio_c(k)/(dise_v(k)*Nml%donner_deep_freq)
        else
          a1_vk(k) = a1     
        endif
      end do

!--------------------------------------------------------------------
!    define the a1 for the column as the smallest of those defined
!    in the column. 
!--------------------------------------------------------------------
      a1 = MINVAL (a1_vk)

!---------------------------------------------------------------------
!    if in a diagnostic column, output the final value of a1, after 
!    all necessary constraints have been applied.
!---------------------------------------------------------------------
     if (debug_ijt) then
       write (diag_unit, '(a, e20.12)') 'in cupar: a1= ',a1
     endif

!--------------------------------------------------------------------


    end subroutine don_d_determine_cloud_area_miz



!######################################################################


subroutine cu_clo_cumulus_closure_miz   &
         (nlev_model, nlev_parcel, ntr, dt, diag_unit, debug_ijt, &
          Initialized, Param, tracers, &
          dcape, pfull, zfull, phalf, zhalf, pblht, tkemiz, qstar, cush, land,  &
          coldT, sd, Uw_p, ac, Nml, cape_p, &!miz
          qli0_v, qli1_v, qr_v, qt_v, env_r, ri_v, rl_v, parcel_r,   &
          env_t, parcel_t, a1, ermesg, error)

!---------------------------------------------------------------------
!    subroutine cumulus_closure calculates a_1(p_b) for closing the 
!    cumulus parameterization. see LJD notes, "Cu Closure D," 6/11/97
!---------------------------------------------------------------------
 
use donner_types_mod,     only : donner_param_type, donner_nml_type, &
                                 donner_initialized_type
use conv_utilities_k_mod, only : pack_sd_lsm_k, extend_sd_k,  &
                                 adi_cloud_k, sounding, adicloud, &
                                 uw_params, qt_parcel_k

implicit none

!---------------------------------------------------------------------
!++++yim
integer,                        intent(in)  :: nlev_model,nlev_parcel,  ntr
real,                           intent(in)  :: dt
integer,                        intent(in)  :: diag_unit
logical,                        intent(in)  :: debug_ijt
type(donner_param_type),        intent(in)  :: Param
type(donner_initialized_type),  intent(in)  :: Initialized
type(donner_nml_type),          intent(in)  :: Nml
real,                           intent(in)  :: dcape
real,                           intent(in)  :: pblht, tkemiz, qstar,  cush, land
logical,                        intent(in)  :: coldT

real, dimension(nlev_model),      intent(in)    :: pfull, zfull !miz
!++++yim
real, dimension(nlev_model,ntr),      intent(in)    :: tracers

real, dimension(nlev_model+1),    intent(in)    :: phalf, zhalf !miz
type(sounding),                 intent(inout) :: sd           !miz
type(uw_params),                 intent(inout) :: Uw_p         !miz
type(adicloud),                 intent(inout) :: ac

real,    dimension(nlev_parcel), intent(in)  :: cape_p, qli0_v, qli1_v, &
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

      real, dimension (nlev_parcel)  :: rt, tden, tdena,  &
                                     dtpdta, pert_env_t, pert_env_r, &
                                     pert_parcel_t, pert_parcel_r,  &
                                     parcel_r_clo, parcel_t_clo, &
                                     ttt, rrr !miz

      real     :: tdens, tdensa, ri1, ri2, rild, rile, rilf, ri2b,  &
                  sum2, rilak, rilbk, rilck, rilakm, rilbkm, rilckm, &
                  rila, rilb, rilc, ri2ak, ri2akm, ri2a, sum1, plcl, &
                  plfc, plzb, dumcoin, dumxcape
       real     :: zsrc, psrc, hlsrc, thcsrc, qctsrc, cape_c, lofactor,&
                   tau, rhavg, dpsum
      integer  :: k     
      logical  :: ctrig, return_cape

      ermesg = ' ' ; error = 0

!--------------------------------------------------------------------
!    initialize the perturbed parcel profiles (pert_parcel_t,    
!    pert_parcel_r) and  the perturbed parcel environmental profiles to 
!    the actual parcel profiles.
!--------------------------------------------------------------------
      do k=1,nlev_parcel
        pert_parcel_t(k) = env_t(k)
        pert_parcel_r(k) = env_r(k)
        pert_env_r(k)    = env_r(k)
        pert_env_t(k)    = env_t(k)
!       ttt (k)          = env_t(nlev_model-k+1)
!       rrr (k)          = env_r(nlev_model-k+1)/(1.-env_r(nlev_model-k+1))
        ttt (k)          = env_t(nlev_parcel-k+1)
        rrr (k)          = env_r(nlev_parcel-k+1)
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
       do k=1,nlev_parcel
         write (diag_unit, '(a, i4, f19.10, f20.14, e20.12)')   &
                     'press, temp, vapor in cape: k, p,t,r = ',  &
                           k, cape_p(k), pert_env_t(k), pert_env_r(k)
       end do
       do k=1,nlev_parcel
         if (qr_v(k) /= 0.0 .or. qt_v(k) /= 0.0) then
             write (diag_unit, '(a, i4, f19.10, 3e20.12, f20.14)') &
                   'in cuclo: k,p,qr,qt,r,t  =', k,  &
                     cape_p(k), qr_v(k), qt_v(k), env_r(k), env_t(k)
         endif
       end do
       do k=1,nlev_parcel
         write (diag_unit, '(a, i4, f19.10, f20.14, e20.12)') &
                     'in cuclo: k,p,tpc, rpc   =', k,   &
                          cape_p(k), parcel_t(k), parcel_r(k)
       end do
       do k=1,nlev_parcel
         if (qli0_v(k) /= 0.0 .or. qli1_v(k) /= 0.0 .or. &
                ri_v(k) /= 0.0 .or. rl_v(k) /= 0.0) then
           write (diag_unit, '(a, i4, f19.10, 4e20.12)')   &
                  'in cuclo: k,p,qli0,qli1,ri,rl =', k,  &
                      cape_p(k), qli0_v(k), qli1_v(k), ri_v(k), rl_v(k)
         endif
       end do
     endif

     if (Nml%do_hires_cape_for_closure) then
 
      if (Nml%do_donner_cape) then
!  there is an existing parcel profiles; it should be same if these 
!  conditions all met so it need not be recalculated.
      if (Nml%do_freezing_for_cape .NEQV. Nml%do_freezing_for_closure .or. &
          Nml%tfre_for_cape /= Nml%tfre_for_closure .or. &
          Nml%dfre_for_cape /= Nml%dfre_for_closure .or. &
          .not. (Initialized%use_constant_rmuz_for_closure) .or.  &
          Nml%rmuz_for_cape /= Nml%rmuz_for_closure) then
           call don_c_displace_parcel_k   &
               (nlev_parcel, diag_unit, debug_ijt, Param,  &
                Nml%do_freezing_for_closure, Nml%tfre_for_closure, &
                Nml%dfre_for_closure, Nml%rmuz_for_closure,  &
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

  else  ! (do_donner_cape)
!   no existing parcel profiles
!    want to use hires cape calc for closure, but used lores cape calc 
!    for convection
          call don_c_displace_parcel_k   &
                (nlev_parcel, diag_unit, debug_ijt, Param,  &
                 Nml%do_freezing_for_closure, Nml%tfre_for_closure, &
                 Nml%dfre_for_closure, Nml%rmuz_for_closure,   &
                 Initialized%use_constant_rmuz_for_closure,  &
                 Nml%modify_closure_plume_condensate, &
                 Nml%closure_plume_condensate, &
                 env_t,  &
                 env_r, cape_p, .false., plfc, plzb, plcl, dumcoin,  &
                 dumxcape, parcel_r_clo,  parcel_t_clo, ermesg, error)
  endif  ! (do_donner_cape)

!--------------------------------------------------------------------
!    call subroutine displace_parcel to determine the movement of a 
!    parcel from the lcl through the environment defined by 
!    (pert_env_t, pert_env_r). set return_cape to indicate if cape
!    value is to be returned.
!    don't need to calculate cape when using Zhang closure 
!    (do_dcape is .true). 
!    if using cape relaxation closure then need to return cape value.
!--------------------------------------------------------------------
      return_cape = .not. (Nml%do_dcape)
      call don_c_displace_parcel_k   &
             (nlev_parcel, diag_unit, debug_ijt, Param,   &
             Nml%do_freezing_for_closure, Nml%tfre_for_closure, &
             Nml%dfre_for_closure, Nml%rmuz_for_closure, &
             Initialized%use_constant_rmuz_for_closure,   &
             Nml%modify_closure_plume_condensate, &
             Nml%closure_plume_condensate, &
             pert_env_t, pert_env_r, cape_p, return_cape, &
             plfc, plzb, plcl, dumcoin,  &
             dumxcape, pert_parcel_r,  pert_parcel_t, ermesg, error)
      if (return_cape) ac%cape = dumxcape
  
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
         
        if (Nml%do_capetau_land) then
         !cape_c = Nml%cape0 * (1. - sd%land * (1. - Nml%lofactor0))
         if (Nml%lochoice > 3) then
           lofactor = 1.
         else
           error = 1
           ermesg = 'unsupported value of lofactor for do_hires_cape'
           return
         endif
         cape_c = Nml%cape0 * lofactor
         tau    = Nml%tau   * lofactor
       endif
     endif
     if (Nml%do_rh_trig) then
         error = 2
         ermesg = 'do_rh_trig not currently supported for do_hires_cape'
         return
!       rhavg=0.; dpsum=0.
!       do k = 1,sd%kmax
!         if (sd%p(k) .gt. Nml%plev0) then
!           rhavg  = rhavg + sd%rh(k)*sd%dp(k)
!           dpsum = dpsum + sd%dp(k)
!         end if
!       end do
!       rhavg = rhavg/dpsum
!       ctrig = rhavg > Nml%rhavg0
     else
       ctrig= .true.
     endif
     endif

   else  ! (do_hires_cape)
 
!  no hires cape calculation for closure; standard donner_lite path
 
!--------------------------------------------------------------------
!    call subroutine displace_parcel to determine the movement of a 
!    parcel from the lcl through the environment defined by 
!    (pert_env_t, pert_env_r).
!--------------------------------------------------------------------
     call pack_sd_lsm_k      &
               (Nml%do_lands, land, coldT, dt, pfull, phalf, zfull,  &
                zhalf, ttt, rrr, tracers, sd)
      call extend_sd_k (sd, pblht, .false., Uw_p)
      zsrc  =sd%zs (1)
      psrc  =sd%ps (1)
      thcsrc=sd%thc(1)
      qctsrc=sd%qct(1)
      hlsrc =sd%hl (1)
      cape_c = Nml%cape0 
      tau    = Nml%tau
      if (Nml%do_lands) then
        call qt_parcel_k (sd%qs(1), qstar, pblht, tkemiz, sd%land, &
             Nml%gama, Nml%pblht0, Nml%tke0, Nml%lofactor0, Nml%lochoice, qctsrc, lofactor)      
        if (Nml%do_capetau_land) then
          !cape_c = Nml%cape0 * (1. - sd%land * (1. - Nml%lofactor0))
          cape_c = Nml%cape0 * lofactor
          tau    = Nml%tau   * lofactor
        end if
      endif
      if (Nml%do_rh_trig) then
        rhavg=0.; dpsum=0.
        do k = 1,sd%kmax
          if (sd%p(k) .gt. Nml%plev0) then
            rhavg  = rhavg + sd%rh(k)*sd%dp(k)
            dpsum = dpsum + sd%dp(k)
          end if
        end do
        rhavg = rhavg/dpsum
        ctrig= rhavg .gt. Nml%rhavg0
      else
        ctrig=.true.
      end if

      call adi_cloud_k (zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, Uw_p, &
                       .false., Nml%do_freezing_for_closure, ac)
     parcel_r_clo=ac%qv(:)/(1.-ac%qv(:))
     parcel_t_clo=ac%t (:)

     sd%t(1) =sd%t (1)-1.0
     sd%qv(1)=0.99*(sd%qv(1)/(1. - sd%qv(1)))

     call adi_cloud_k (zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, Uw_p, &
                       .false., Nml%do_freezing_for_closure, ac) 
      pert_parcel_r=ac%qv(:)/(1.-ac%qv(:))
      pert_parcel_t=ac%t (:)
!miz

  endif ! (hires cape)
 

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    if in a diagnostics column, output the path of the parcel (T, p
!    coordinates).
!---------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_parcel
          write (diag_unit, '(a, i4, f20.14, e20.12)')  &
                      'in cuclo: k,tpca,rpca= ', k,    &
                                 pert_parcel_t(k), pert_parcel_r(k)
        end do
      endif

!---------------------------------------------------------------------
!    calculate the large-scale model profile of total-water mixing 
!    ratio. 
!---------------------------------------------------------------------
        do k=1,nlev_parcel
        rt(k) = env_r(k) + ri_v(k) + rl_v(k)
      end do

!----------------------------------------------------------------------
!    calculate profiles of density temperatures, in the parcel (tden) 
!    and in the perturbed parcel (tdena). condensate is not included in
!    this definition of density temperature.
!----------------------------------------------------------------------
      do k=1,nlev_parcel
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
      do k=1,nlev_parcel
        dtpdta(k) = (tdena(k) - tden(k))/(tdensa - tdens)
      end do

!---------------------------------------------------------------------
!    if this is a diagnostics column, output the profiles of unperturbed
!    parcel density temperature (tden) and the perturbed parcel density 
!    temperature (tdena) and the derivative of parcel density temper-
!    ature w.r.t. cloud-base large-scale density temperature (dtpdta).
!------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_parcel
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

      if (Nml%model_levels_in_sfcbl == 0) then
        sum2 = rild + rile + rilf
      else
        sum2 = 0.
      endif


      ri1 = 0.
      ri2 = 0.
      do k=2,nlev_parcel
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
          a1  = -(ri2 + dcape)/ri1
        else
          if (ac%cape .gt. cape_c .and. ctrig) then
            a1  = -(ri2 + (ac%cape-cape_c)/tau)/ri1
          else
            a1  = 0.
          end if
        end if
      endif

!--------------------------------------------------------------------


    end subroutine cu_clo_cumulus_closure_miz


!######################################################################
!######################################################################






subroutine don_cm_mesub_miz     &
         (Nml, pfull_c, nlev_lsm, me, diag_unit, debug_ijt, Param, cu, &
          ci_liq_cond, ci_ice_cond, pmelt_lsm, cell_precip, &
          dint, plzb_c, pb, pt_kou, temp_c, phalf_c,   &
          ca_liq, ca_ice, ecd, ecd_liq, ecd_ice, ecei_liq, &
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
real,                          intent(in)    :: pmelt_lsm
real,                          intent(in)    :: ci_liq_cond, &
                                                ci_ice_cond
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
!       dint??       water mass frozen in convective updraft
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
      integer  :: itrop
      real :: ptrop

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
      ecdi_liq  = (Param%evap_in_downdrafts*avail_meso_cd)*  &
                       (Param%seconds_per_day*ci_liq_cond)
      ecdi_ice  = (Param%evap_in_downdrafts*avail_meso_cd)* &
                       (Param%seconds_per_day*ci_ice_cond)
      ecei  = (Param%evap_in_environ*avail_meso_cd)*cu
      ecei_liq  = (Param%evap_in_environ*avail_meso_cd)*  &
                       (Param%seconds_per_day*ci_liq_cond)
      ecei_ice  = (Param%evap_in_environ*avail_meso_cd)*   &
                       (Param%seconds_per_day*ci_ice_cond)
      ca_liq    = (Param%entrained_into_meso*avail_meso_cd)*  &
                       (Param%seconds_per_day*ci_liq_cond)
      ca_ice    = (Param%entrained_into_meso*avail_meso_cd)*  &
                       (Param%seconds_per_day*ci_ice_cond)
     if (debug_ijt) then
        write (diag_unit, '(a, 3e20.12)')  &
             'in mesub: cu, h1_liqintg, h1_iceintg= ', &
             cu, ci_liq_cond*Param%seconds_per_day, &
                 ci_ice_cond*Param%seconds_per_day
      endif

!---------------------------------------------------------------------
!    if in a diagnostics column, output the ratio of convective rain 
!    to convective updraft condensation (gnu) and the mass entrained
!    into the mesoscale region (ca).
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12)')  'in mesub: gnu= ',gnu
      write (diag_unit, '(a, e20.12)') 'in mesub: ca= ',ca_liq + ca_ice 
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
      dint2 = avail_meso_cd*(dint)*8.64e07*Param%grav/(pzm - pztm)
 
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
                           'in cm_intgl_to_gcm_col: k,x= ',k, ecd   (k)
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
                   'MELTING DONE: pmelt_lsm,pb,caa',pmelt_lsm, pb, &
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


end subroutine don_cm_mesub_miz

!#######################################################################
!#######################################################################

subroutine don_m_meso_effects_miz    &
         (me, nlev_lsm, nlev_hires, ntr, diag_unit, debug_ijt, Param, Nml,&
          pfull_c, temp_c, mixing_ratio_c, phalf_c, rlsm, emsm, etsm, &
          tracers_c, ensembl_cond, ensmbl_precip, pb, plzb_c, pt_ens, &
          ampta1, ensembl_anvil_cond_liq, ensembl_anvil_cond_liq_frz, &
          ensembl_anvil_cond_ice,  &
          wtp, qtmes, meso_frz_intg_sum, anvil_precip_melt, &
          meso_cloud_area, cmus_tot, dmeml, emds_liq, emds_ice, &
          emes_liq, emes_ice, wmms, wmps, &
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
real,   dimension(nlev_lsm+1),     intent(in)  :: phalf_c, rlsm,  &
                                                  emsm, etsm !miz
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
                                                  cmus_tot, dmeml,  &
                                                  emds_liq, emds_ice, &
                                                  emes_liq, emes_ice, &
                                                  mrmes_up, mrmes_dn, &
                                                   tmes_up, tmes_dn, &
                                                        wmms, wmps,   &
                                                  umeml, tmes, mrmes
real,                              intent(out) ::  emdi,pmd, pztm, pzm,&
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
      real  :: emdi_liq, emdi_ice
      real  :: intgl_lo, intgl_hi
      integer                        ::  k, kcont, itrop
      real          :: p2, ptrop

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

      emes_liq = 0.
      emes_ice = 0.
      emdi_liq = 0.
      emdi_ice = 0.
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
      call don_m_meso_updraft_miz   &
           (nlev_lsm, nlev_hires, ntr, diag_unit, debug_ijt, Param, &
            pfull_c, rlsm, emsm, etsm, pfull_c,  &
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
      call don_m_meso_downdraft_miz  &
          (nlev_lsm, nlev_hires, diag_unit, debug_ijt, Param, pfull_c,&
           pfull_c, temp_c, mixing_ratio_c, phalf_c, pb, ampta1, dp,  &
           pztm, pzm, alp, hfmin, pmd, tmes_dn, mrmes_dn, dmeml, ermesg, error)
 
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
        if (.not. meso_frz_intg_sum ) then
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
        write (diag_unit, '(a, 5e20.12)')  &
                     'in meens: cmui,ca (sum,liq,frzliq,ice)=',  cmui, &
           ensembl_anvil_cond_liq + ensembl_anvil_cond_liq_frz + &
                              ensembl_anvil_cond_ice, &
                                         ensembl_anvil_cond_liq,  &
                                 ensembl_anvil_cond_liq_frz,  &
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
           (nlev_lsm, diag_unit, debug_ijt, Param,    &
            available_condensate, available_condensate_liq,  &
            available_condensate_ice, pzm, pztm, phalf_c,       &
            emdi_liq, emdi_ice,       &
            emds_liq, emds_ice, &
                  emes_liq, emes_ice, ermesg, error)

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


    end subroutine don_m_meso_effects_miz

!#######################################################################
!#######################################################################



!#######################################################################

subroutine don_m_meso_updraft_miz    &
         (nlev_lsm, nlev_hires, ntr, diag_unit, debug_ijt, Param,  &
          p_hires, rlsm, emsm, etsm, pfull_c, temp_c, mixing_ratio_c, &
          phalf_c, tracers_c, pb, pt_ens, ampta1, dp, pztm, wtp, &
          qtmes, cmu, wmms, wmps, temptr, tmes_up, mrmes_up,   &
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
real,   dimension(nlev_lsm),     intent(in)  :: p_hires, rlsm, emsm !miz
real,   dimension(nlev_lsm,ntr), intent(in)  :: etsm !miz
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



      real, dimension (nlev_lsm)         :: wmhr, cumh !miz
      real, dimension (nlev_lsm)         :: omv, tempq, owm, tempqa
      real, dimension(nlev_lsm,ntr)      :: otm
      real, dimension(nlev_lsm, ntr)     :: wthr !miz
      real, dimension(ntr)               :: q1t


      real      ::  cmfhr, pc1, pc2, omer, pctm, q1, q4, mrsat, &
                    q3, anv, qref, pp, pm, qprip, qprim, eqfp, eqfm, &
                    qmu, hflux, pfmin, owms, wpc, wmc, ta, te, tep, tmu,&
                    qtprip, qtprim, eqtfp, eqtfm
      logical   :: do_donner_tracer
      integer   :: ncc
      integer   :: kcont, kk
      integer   :: jk, jsave, jkm, jkp, k, nbad

!-----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

      if (ntr > 0) then
        do_donner_tracer = .true.
      else
        do_donner_tracer = .false.
      endif
!miz
         do k=1,nlev_hires
           if (p_hires(k) < pt_ens) then
             ncc = k
             exit
           endif
         end do
!!$      do i=1,nlev_hires
!!$        if (p_hires(i) < pztm) then
!!$          ncztm = i + 1
!!$          exit
!!$        endif
!!$      end do

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
      do k=1,nlev_lsm !miz
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
            pzm = pfull_c(k) !miz
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
      do k=1,nlev_lsm  
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

!!$!---------------------------------------------------------------------
!!$!    convert the vertical profile of vapor made available to the meso-
!!$!    scale from the updraft to the large-scale model grid (output var-
!!$!    iable is owm). if tracers are being transported by donner conv-
!!$!    ection, convert the vertical profile of tracer made available to 
!!$!    the mesoscale from the updraft to the large-scale model grid 
!!$!    (output variable is otm). 
!!$!---------------------------------------------------------------------
!!$      call don_u_map_hires_c_to_lores_c_k &
!!$           (nlev_lsm, nlev_hires, wmhr, p_hires, pt_ens + dp, phalf_c,&
!!$            owm, rintsum, rintsum2, ermesg)
!!$      if (trim(ermesg) /= ' ') return

      owm=wmhr
!mizdelete

      if (do_donner_tracer) then
!!$        do kcont=1,ntr
!!$          call don_u_map_hires_c_to_lores_c_k  &
!!$               (nlev_lsm, nlev_hires, wthr (:,kcont), p_hires,  &
!!$                pt_ens + dp, phalf_c, otm(:,kcont), rintsum,   &
!!$                rintsum2, ermesg) 
!!$          if (trim(ermesg) /= ' ') return
!!$!mizdelete
!!$        end do
         otm=wthr
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
!    will only be non-zero at layers within the mesoscale updraft, but 
!    owm may be non-zero in layers above the updraft.
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
         call compute_mrs_k (te, pfull_c(k), Param%d622 , Param%d608 ,&
                             mrsat, nbad , mr = qref)
 
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
            call compute_mrs_k (tep, pfull_c(jkp), Param%d622 ,  &
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
            call compute_mrs_k (tep, pfull_c(jkp), Param%d622 ,  &
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
        call compute_mrs_k (tmu, pfull_c(jk), Param%d622 ,  &
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
      do k=1,nlev_lsm
        if ((p_hires(k) <= pzm) .and. (p_hires(k) >= pztm))  then
          cumh(k) = ampta1
        else
          cumh(k) = 0.0 
        endif
      end do

!!$      call don_u_map_hires_c_to_lores_c_k  &
!!$           (nlev_lsm, nlev_hires, cumh, p_hires, pztm + dp, phalf_c, &
!!$            meso_cloud_area, rintsum, rintsum2, ermesg) 
!!$       if (trim(ermesg) /= ' ') return

       meso_cloud_area=cumh

!mizdelete

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


    end subroutine don_m_meso_updraft_miz



!#####################################################################

subroutine don_m_meso_downdraft_miz    &
         (nlev_lsm, nlev_hires, diag_unit, debug_ijt,  Param, p_hires, &
          pfull_c, temp_c, mixing_ratio_c, phalf_c, pb, ampta1, dp,  &
          pztm, pzm, alp, hfmin, pmd, tmes_dn, mrmes_dn, dmeml, ermesg, error)

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
integer,                       intent(in)   :: nlev_lsm, nlev_hires, &
                                               diag_unit
logical,                       intent(in)   :: debug_ijt        
type(donner_param_type),       intent(in)   :: Param
real,   dimension(nlev_lsm),   intent(in)   :: p_hires !miz
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

      real, dimension(nlev_lsm)       :: dmemh !miz
      real, dimension(nlev_lsm)       :: tempt, tempqa
      real, dimension(nlev_lsm+1)     :: emt, emq

     real    ::  es, mrsat, c2, c3, c1, fjk, fjkm, qb, fjkb, qbm, qmd, &
                  qsmd, fjkmd, qmmd, pi, psa,  &
                        targ, tprimd, tb, qten, tten, omd, mrsb, wa,   &
                  wb, tmd, rin
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
!miz
!!$      ncmd = 1
!!$      do k=1,nlev_hires         
!!$        if (p_hires(k) < pmd ) then
!!$          ncmd = k + 1
!!$          exit
!!$        endif
!!$      end do

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
        call compute_mrs_k (targ, pfull_c(k), Param%d622 ,  &
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
        call compute_mrs_k (targ, pfull_c(k), Param%d622 ,  &
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
        call compute_mrs_k (tb, pb, Param%d622 ,  &
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
          call compute_mrs_k (targ, pb, Param%d622 ,  &
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
          call compute_mrs_k (tmd, pmd, Param%d622 ,  &
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
          call compute_mrs_k (targ, pmd, Param%d622 ,  &
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
          tmes_dn(k) = tmes_dn(k) + (tten/pi)
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
      do k=1,nlev_lsm !miz
        if ((p_hires(k) <= pb) .and. (p_hires(k) >= pmd))  then
          dmemh(k) = -omd*ampta1/Param%grav  
        else
          dmemh(k) = 0.
        endif
      end do

!!$!---------------------------------------------------------------------
!!$!    call map_hi_res_col_to_lo_res_col to map the 
!!$!    mesoscale downdraft flux from the cloud model to the large-scale 
!!$!    model.
!!$!---------------------------------------------------------------------
!!$      call don_u_map_hires_c_to_lores_c_k  &
!!$           (nlev_lsm, nlev_hires, dmemh, p_hires, pmd + dp, phalf_c, &
!!$            dmeml, rintsum, rintsum2, ermesg) 
!!$      if (trim(ermesg) /= ' ') return

      dmeml=dmemh

    end subroutine don_m_meso_downdraft_miz

!######################################################################
!######################################################################


subroutine don_l_lscloud_driver_miz   &
         (isize, jsize, nlev_lsm, cloud_tracers_present, Param,  &
          Col_diag, pfull, temp,  exit_flag,  &
          mixing_ratio, qlin, qiin, qain, phalf, Don_conv, &
          donner_humidity_factor, donner_humidity_area,  &
           dql, dqi, dqa, mhalf_3d, &
          ermesg, error) 

!---------------------------------------------------------------------
!    subroutine don_l_lscloud_driver obtains variables needed by 
!    strat_cloud_mod that are dependent on the donner_deep parameter-
!    ization. specifically, the convective cell plus mesoscale anvil
!    cloud fraction (donner_humidity_area), the ratio of the large-scale 
!    specific humidity to the specific humidity in the environment out-
!    side of the convective system (donner_humidity_ratio), and the 
!    changes in cloud liquid, cloud ice and cloud area due to the con-
!    vective-system vertical mass flux and detrainment from the mesoscale
!    anvil to the large scale (dql, dqi, dqa) are passed out for use in 
!    strat_cloud_mod.
!---------------------------------------------------------------------

use donner_types_mod, only : donner_conv_type, donner_param_type, &
                             donner_column_diag_type

implicit none

!---------------------------------------------------------------------
integer,                    intent(in)    :: isize, jsize, nlev_lsm
logical,                    intent(in)    :: cloud_tracers_present
type(donner_param_type),    intent(in)    :: Param
type(donner_column_diag_type),                    &
                            intent(in)    :: Col_diag
real, dimension(isize,jsize,nlev_lsm),         &
                            intent(in)    :: pfull, temp, mixing_ratio, &
                                             qlin, qiin, qain
logical, dimension(isize, jsize), intent(in) :: exit_flag
real, dimension(isize,jsize,nlev_lsm+1),        &
                            intent(in)    :: phalf 
type(donner_conv_type),     intent(inout) :: Don_conv
real, dimension(isize,jsize,nlev_lsm),           &
                            intent(out)   :: donner_humidity_factor,  &
                                             donner_humidity_area, dql, &
                                             dqi, dqa
real, dimension   (isize, jsize, nlev_lsm+1), &
                            intent(out)   :: mhalf_3d
character(len=*),           intent(out)   :: ermesg
integer,                    intent(out)   :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     pfull          pressure field on model full levels [ Pa ]
!     phalf          pressure field on model half levels [ Pa ]
!     temp           temperature field at model full levels [ deg K ]
!     mixing_ratio   water vapor specific humidity at model full 
!                    levels [ kg(h2o) / kg(air) ]
!     qlin           large-scale cloud liquid specific humidity
!                    [ kg(h2o) / kg(air) ]
!     qiin           large-scale cloud ice specific humidity
!                    [ kg(h2o) / kg(air) ]
!     qain           large-scale cloud fraction [ fraction ]
!
!   intent(inout) variables:
!
!     Don_conv
!
!
!   intent(out) variables:
!
!     donner_humidity_ratio
!                    ratio of large-scale specific humidity to the
!                    specific humidity in the environment outside
!                    of the convective system [ dimensionless ]
!     donner_humidity_area
!                    fractional area of cell plus meso circulation
!                    associated with donner_deep_mod [ fraction ]
!     dql            increment to large-scale cloud liquid field from
!                    donner_deep_mod [ kg(h2o) / kg (air) ]
!     dqi            increment to large-scale cloud ice field from
!                    donner_deep_mod [ kg(h2o) / kg (air) ]
!     dqa            increment to large-scale cloud area field from
!                    donner_deep_mod [ fraction ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:


      real, dimension (isize, jsize,nlev_lsm) ::  mass  
      integer     :: k, n, i, j

      real, dimension   &
            (isize, jsize, nlev_lsm  ) :: dmeso_3d

!---------------------------------------------------------------------
!   local variables:
!
!     dmeso_3d       detrainment rate from convective system 
!                    [ sec**(-1) ]
!     mhalf_3d       mass flux at model half-levels 
!                    [ kg / (m**2 sec) ]
!---------------------------------------------------------------------

      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    call define_donner_mass_flux to define the convective system 
!    detrainment rate (dmeso_3d) and the mass flux at model interface 
!    levels (mhalf_3d) that is associated with deep convection.
!---------------------------------------------------------------------
      call don_l_define_mass_flux_k    &
           (isize, jsize, nlev_lsm, pfull, phalf, Don_conv,   &
            dmeso_3d, mhalf_3d, exit_flag, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

 
!---------------------------------------------------------------------
!    call adjust_tiedtke_inputs to obtain the convective cloud area 
!    (donner_humidity_area) and the ratio of large-scale specific humid-
!    ity to the humidity in the environment of the convective system 
!    (donner_humidity_ratio).
!---------------------------------------------------------------------
      call don_l_adjust_tiedtke_inputs_k    &
           (isize, jsize, nlev_lsm, Param, Col_diag, pfull,temp,   &
            mixing_ratio, phalf, Don_conv, donner_humidity_factor, &
            donner_humidity_area, exit_flag, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return
 
!---------------------------------------------------------------------
!    when strat_cloud is active, call strat_cloud_donner_tend to
!    define increments to cloudice, cloudwater and cloud area associated
!    with deep convective vertical mass flux and detrainment from the
!    mesoscale to the large-scale. 
!---------------------------------------------------------------------
      if (cloud_tracers_present) then
!!$        call don_l_strat_cloud_donner_tend_k   &
!!$             (isize, jsize, nlev_lsm, Param, Col_diag, dmeso_3d,   &
!!$              Don_conv%xliq, Don_conv%xice, qlin, qiin, qain, mhalf_3d, &
!!$              phalf, dql, dqi, dqa, ermesg)
      do k=1,nlev_lsm
        mass(:,:,k) = (phalf(:,:,k+1) - phalf(:,:,k))/Param%grav 
      end do

!---------------------------------------------------------------------
!    define the large scale cloud increments at level 1 to be 0.0.
!---------------------------------------------------------------------
      dql (:,:,1) = 0.
      dqi (:,:,1) = 0.
      dqa (:,:,1) = 0.
    
!---------------------------------------------------------------------
!    define the tendencies of cloud liquid, cloud ice and cloud area
!    due to the vertical mass flux associated with donner_deep con-
!    vection.
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!    add the effects of detrainment from the mesoscale region.
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        dql (:,:,k) = Don_conv%ecds(:,:,k)
        dqi (:,:,k) = Don_conv%eces(:,:,k)
        dqa (:,:,k) = Don_conv%fre (:,:,k)
      end do

      do i=1,isize  
         do j=1,jsize  
            do k=1,nlev_lsm
               if ( (pfull(i,j,k) .le. Don_conv%pzm_v (i,j)) .and. &
                    (pfull(i,j,k) .ge. Don_conv%pztm_v(i,j)) ) then
!                  meso_area_miz(i,j,k)=Don_conv%ampta1(i,j)
!                  meso_updt_miz(i,j,k)=Param%meso_ref_omega
               else
!                  meso_area_miz(i,j,k)=0.
!                  meso_updt_miz(i,j,k)=0.
               end if
            end do
         end do
      end do

 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return
      endif

!---------------------------------------------------------------------


    end subroutine don_l_lscloud_driver_miz

!#####################################################################
!######################################################################

subroutine don_d_remove_normalization_miz   &
      (isize, jsize, nlev_lsm, ntr, exit_flag, Nml, Don_conv,  &
       total_precip, Initialized, temperature_forcing,  &
       moisture_forcing, ermesg, error)

!---------------------------------------------------------------------
!    subroutine remove_normalization removes the normalization by the
!    cloud base fractional area from the various convective diagnostics
!    and output fields so that they are ready fro use in the large-scale
!    model equations.
!---------------------------------------------------------------------

use donner_types_mod, only : donner_conv_type, donner_nml_type, &
                             donner_initialized_type, DET_MASS_FLUX, &
                             MASS_FLUX, CELL_UPWARD_MASS_FLUX, &
                             TEMP_FORCING, MOIST_FORCING, PRECIP, &
                             FREEZING, RADON_TEND

implicit none 

!---------------------------------------------------------------------
integer,                          intent(in)    :: isize, jsize, &
                                                   nlev_lsm, ntr
logical, dimension(isize,jsize),  intent(in)    :: exit_flag
type(donner_nml_type),            intent(in)    :: Nml      
type(donner_conv_type),           intent(inout) :: Don_conv
type(donner_initialized_type),    intent(inout) :: Initialized
real   , dimension(isize,jsize),  intent(inout) :: total_precip
real   , dimension(isize,jsize,nlev_lsm),                 &
                                  intent(inout) :: temperature_forcing, &
                                                   moisture_forcing
character(len=*),                 intent(out)   :: ermesg
integer,                          intent(out)   :: error
!----------------------------------------------------------------------
!   intent(in) variables:
!
!     exit_flag      logical array indicating whether donner convection
!                    is not active (.true.) or is active (.false.) in
!                    each model column 
!
!   intent(inout) variables:
!    
!     Don_conv       donner_convection_type derived type variable 
!                    containing fields produced by the donner_deep
!                    convection mod 
!     total_precip   precipitation generated by deep convection
!                    [ kg / m**2 ]
!     moisture_forcing
!                    time tendency of vapor mixing ratio due to deep 
!                    convection [ kg(h2o) / kg(dry air) / sec ]
!     temperature_forcing
!                    time tendency of temperature due to deep 
!                    convection [ deg K / sec ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer :: i, j, k, n    ! do-loop indices
      real    :: ttend_max
      real, dimension(nlev_lsm) :: variable

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    remove normalization from the cumulus diagnostics and forcing terms
!    by multiplying them by the fractional cloud base area. these values
!    thus become grid-box averages, rather than averages over the cloudy
!    area, and so are appropriate to use in the large-scale model
!    equations. 
!---------------------------------------------------------------------
      do j=1,jsize                          
        do i=1,isize

!---------------------------------------------------------------------
!    if deep convection is present in the column, denormalize the 
!    convective fields.
!--------------------------------------------------------------------
          if (.not. exit_flag(i,j)) then
            if (Initialized%monitor_output) then
              do n=1, size(Initialized%Don_monitor, 1)
                select case (Initialized%Don_monitor(n)%index)
                  case (DET_MASS_FLUX)
                    variable(:) = Don_conv%detmfl(i,j,:)*   &
                                                       Don_conv%a1(i,j)
                    call don_u_process_monitor_k (variable, i, j,  &
                                   nlev_lsm, Initialized%Don_monitor(n))
                  case (MASS_FLUX)
                    variable(:) =   &
                      (Don_conv%umeml(i,j,:) + Don_conv%dmeml(i,j,:) + &
                       Don_conv%uceml(i,j,:))*Don_conv%a1(i,j) 
                    call don_u_process_monitor_k (variable, i, j,  &
                                   nlev_lsm, Initialized%Don_monitor(n))
                     
                  case (CELL_UPWARD_MASS_FLUX)
                    variable(:) = Don_conv%uceml(i,j,:)*Don_conv%a1(i,j)
                    call don_u_process_monitor_k (variable, i, j,  &
                                   nlev_lsm, Initialized%Don_monitor(n))
                  case (TEMP_FORCING)
                    variable(:) =   &
                           temperature_forcing(i,j,:)*Don_conv%a1(i,j) 
                    call don_u_process_monitor_k (variable, i, j,  &
                                   nlev_lsm, Initialized%Don_monitor(n))

                  case (MOIST_FORCING)
                    variable(:) =   &
                               moisture_forcing(i,j,:)*Don_conv%a1(i,j) 
                    call don_u_process_monitor_k (variable, i, j,  &
                                   nlev_lsm, Initialized%Don_monitor(n))
                  case (PRECIP)
                    variable(:) = total_precip(i,j)*Don_conv%a1(i,j) 
                    call don_u_process_monitor_k (variable, i, j,  &
                                   nlev_lsm, Initialized%Don_monitor(n))
                  case (FREEZING)
                    variable(:) = Don_conv%fre(i,j,:)*Don_conv%a1(i,j) 
                    call don_u_process_monitor_k (variable, i, j,  &
                                   nlev_lsm, Initialized%Don_monitor(n))
                end select
              end do
            endif 
                     
!miz
             ttend_max=0;
             do k=1,nlev_lsm 
                ttend_max=max(ttend_max, abs(temperature_forcing(i,j,k))*Don_conv%a1(i,j))
             end do
             if (ttend_max > Nml%ttend_max) then
                Don_conv%a1(i,j)=Don_conv%a1(i,j)*(Nml%ttend_max/ttend_max)
             end if
!miz

            total_precip(i,j) =  total_precip(i,j)*Don_conv%a1(i,j)
            Don_conv%meso_precip(i,j) = Don_conv%meso_precip(i,j)* &
                 Don_conv%a1(i,j)
            Don_conv%ampta1(i,j) = Don_conv%ampta1(i,j)*Don_conv%a1(i,j)
            Don_conv%cell_precip(i,j) =              &
                             Don_conv%cell_precip (i,j)*Don_conv%a1(i,j)
            Don_conv%emdi_v(i,j) = Don_conv%emdi_v(i,j)*Don_conv%a1(i,j)
            do k=1,nlev_lsm                           
              Don_conv%wetdepc(i,j,k,:) = &
                              Don_conv%wetdepc(i,j,k,:)*Don_conv%a1(i,j)
              Don_conv%wetdept(i,j,k,:) = &
                              Don_conv%wetdepc(i,j,k,:)
              temperature_forcing(i,j,k) =   &
                             temperature_forcing(i,j,k)*Don_conv%a1(i,j)
              Don_conv%ceefc(i,j,k) =   &
                                  Don_conv%ceefc(i,j,k)*Don_conv%a1(i,j)
              Don_conv%cecon(i,j,k) =        &
                                  Don_conv%cecon(i,j,k)*Don_conv%a1(i,j)
              Don_conv%cemfc(i,j,k) =      &
                                  Don_conv%cemfc(i,j,k)*Don_conv%a1(i,j)
              moisture_forcing(i,j,k) =      &
                                moisture_forcing(i,j,k)*Don_conv%a1(i,j)
              Don_conv%cual (i,j,k) =       &
                                   Don_conv%cual(i,j,k)*Don_conv%a1(i,j)
              Don_conv%fre(i,j,k) = Don_conv%fre(i,j,k)*Don_conv%a1(i,j)
              Don_conv%elt(i,j,k) = Don_conv%elt(i,j,k)*Don_conv%a1(i,j)
              Don_conv%cmus(i,j,k) =      &
                                   Don_conv%cmus(i,j,k)*Don_conv%a1(i,j)
              Don_conv%ecds(i,j,k) =      &
                                   Don_conv%ecds(i,j,k)*Don_conv%a1(i,j)
              Don_conv%eces(i,j,k) =      &
                                   Don_conv%eces(i,j,k)*Don_conv%a1(i,j)
              Don_conv%emds(i,j,k) =       &
                                   Don_conv%emds(i,j,k)*Don_conv%a1(i,j)
              Don_conv%emes(i,j,k) =       &
                                   Don_conv%emes(i,j,k)*Don_conv%a1(i,j)
              Don_conv%mrmes(i,j,k) =       &
                                  Don_conv%mrmes(i,j,k)*Don_conv%a1(i,j)
              Don_conv%wmps(i,j,k) =       &
                                   Don_conv%wmps(i,j,k)*Don_conv%a1(i,j)
              Don_conv%wmms(i,j,k) =      &
                                   Don_conv%wmms(i,j,k)*Don_conv%a1(i,j)
              Don_conv%tmes(i,j,k) =      &
                                   Don_conv%tmes(i,j,k)*Don_conv%a1(i,j)
              Don_conv%dmeml(i,j,k) =      &
                                  Don_conv%dmeml(i,j,k)*Don_conv%a1(i,j)
              Don_conv%uceml(i,j,k) =      &
                                  Don_conv%uceml(i,j,k)*Don_conv%a1(i,j)
              Don_conv%detmfl(i,j,k) =      &
                                 Don_conv%detmfl(i,j,k)*Don_conv%a1(i,j)
              Don_conv%umeml(i,j,k) =      &
                                  Don_conv%umeml(i,j,k)*Don_conv%a1(i,j)
              Don_conv%qtren1(i,j,k,:) =     &
                               Don_conv%qtren1(i,j,k,:)*Don_conv%a1(i,j)
              Don_conv%qtmes1(i,j,k,:) =     &
                               Don_conv%qtmes1(i,j,k,:)*Don_conv%a1(i,j)
              Don_conv%wtp1(i,j,k,:) =       &
                                 Don_conv%wtp1(i,j,k,:)*Don_conv%a1(i,j)
              Don_conv%qtceme(i,j,k,:) =   &
                     Don_conv%qtmes1(i,j,k,:) + Don_conv%qtren1(i,j,k,:)
            end do
            if (Initialized%monitor_output) then
              do n=1, size(Initialized%Don_monitor, 1)
                select case (Initialized%Don_monitor(n)%index)
                  case (RADON_TEND)
                    variable(:) = Don_conv%qtceme   &
                        (i,j,:,Initialized%Don_monitor(n)%tracer_index) 
                    call don_u_process_monitor_k (variable, i, j,  &
                                   nlev_lsm, Initialized%Don_monitor(n))

                end select
              end do
            endif 
!---------------------------------------------------------------------
!    if deep convection is not present in the column, define the output
!    fields appropriately.
!---------------------------------------------------------------------
          else
            total_precip(i,j) = 0.
            do k=1,nlev_lsm
              temperature_forcing(i,j,k) = 0.
              moisture_forcing(i,j,k) = 0.
            end do
          endif
        end do
      end do

!---------------------------------------------------------------------

    end subroutine don_d_remove_normalization_miz

!#####################################################################


    subroutine don_d_copy_wetdep_miz( wetdep_donner, wetdep_plume, nspecies )
 
   use  conv_plumes_k_mod, only : cwetdep_type
   use donner_types_mod, only : donner_wetdep_type
 
   implicit none
 
   type(donner_wetdep_type), intent(in) :: wetdep_donner(nspecies)
   type(cwetdep_type), intent(inout)    :: wetdep_plume(nspecies)
   integer, intent(in)                  :: nspecies

   integer :: n

   do n = 1,nspecies
      wetdep_plume(n)%scheme = wetdep_donner(n)%scheme
      wetdep_plume(n)%Henry_constant = wetdep_donner(n)%Henry_constant
      wetdep_plume(n)%Henry_variable = wetdep_donner(n)%Henry_variable
      wetdep_plume(n)%frac_in_cloud = wetdep_donner(n)%frac_in_cloud
      wetdep_plume(n)%alpha_r = wetdep_donner(n)%alpha_r
      wetdep_plume(n)%alpha_s = wetdep_donner(n)%alpha_s
      wetdep_plume(n)%Lwetdep = wetdep_donner(n)%Lwetdep
      wetdep_plume(n)%Lgas = wetdep_donner(n)%Lgas
      wetdep_plume(n)%Laerosol = wetdep_donner(n)%Laerosol
      wetdep_plume(n)%Lice = wetdep_donner(n)%Lice
   end do

   end subroutine don_d_copy_wetdep_miz



!####################################################################
!######################################################################


subroutine cu_clo_miz   &
         (nlev_lsm, ntr, dt, Initialized, Param, tracers, &
          pfull, zfull, phalf, zhalf, pblht, tkemiz, qstar, cush, cbmf, land,  &
          coldT, sd, Uw_p, ac, Nml, env_r, env_t, a1, ermesg, error)

use donner_types_mod,     only : donner_param_type, donner_nml_type, &
                                 donner_initialized_type
use conv_utilities_k_mod, only : pack_sd_lsm_k, extend_sd_k,  &
                                 adi_cloud_k, sounding, adicloud, &
                                 uw_params, qt_parcel_k

implicit none

integer,                        intent(in)  :: nlev_lsm, ntr
real,                           intent(in)  :: dt
type(donner_param_type),        intent(in)  :: Param
type(donner_initialized_type),  intent(in)  :: Initialized
type(donner_nml_type),          intent(in)  :: Nml
real,                           intent(in)  :: pblht, tkemiz, qstar, cush, cbmf, land
logical,                        intent(in)  :: coldT
real, dimension(nlev_lsm),      intent(in)    :: pfull, zfull 
real, dimension(nlev_lsm,ntr),  intent(in)    :: tracers
real, dimension(nlev_lsm+1),    intent(in)    :: phalf, zhalf 
type(sounding),                 intent(inout) :: sd           
type(uw_params),                intent(inout) :: Uw_p         
type(adicloud),                 intent(inout) :: ac

real,    dimension(nlev_lsm),   intent(in)  :: env_r, env_t
real,                           intent(out) :: a1
character(len=*),               intent(out) :: ermesg
integer,                        intent(out) :: error

real, dimension (nlev_lsm)  :: ttt, rrr
real    :: zsrc, psrc, hlsrc, thcsrc, qctsrc, cape_c, lofactor
integer :: k
real    :: sigmaw, wcrit, cbmf1, rbuoy, rkfre, wcrit_min, rmaxfrac

   ermesg = ' '; error = 0

   do k=1,nlev_lsm
      ttt (k) = env_t(nlev_lsm-k+1)
      rrr (k) = env_r(nlev_lsm-k+1)/(1.-env_r(nlev_lsm-k+1))
   end do

   call pack_sd_lsm_k      &
        (Nml%do_lands, land, coldT, dt, pfull, phalf, zfull,  &
        zhalf, ttt, rrr, tracers, sd)
   call extend_sd_k (sd, pblht, .false., Uw_p)
   zsrc  =sd%zs (1)
   psrc  =sd%ps (1)
   thcsrc=sd%thc(1)
   qctsrc=sd%qct(1)
   hlsrc =sd%hl (1)
   cape_c=Nml%cape0
   if (Nml%do_lands) then
      call qt_parcel_k (sd%qs(1), qstar, pblht, tkemiz, sd%land, &
           Nml%gama, Nml%pblht0, Nml%tke0, Nml%lofactor0, Nml%lochoice, qctsrc, lofactor)
      cape_c = Nml%cape0 * (1. - sd%land * (1. - Nml%lofactor0))
   endif
   call adi_cloud_k (zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, Uw_p, &
                        .false., Nml%do_freezing_for_closure, ac)

   rbuoy     = 1.0
   rkfre     = 0.05
   wcrit_min = 0.5
   rmaxfrac  = 0.05

   wcrit  = sqrt(2. * ac % cin * rbuoy)
   sigmaw = sqrt(rkfre * tkemiz)
   wcrit  = max(wcrit, wcrit_min*sigmaw)
   cbmf1   = ac % rho0lcl * sigmaw / 2.5066 * exp(-0.5*((wcrit/sigmaw)**2.))
   cbmf1   = min(cbmf, ( sd%ps(0) - ac%plcl ) * 0.25 / sd%delt / Uw_p%GRAV)

   if (cush > Nml%deephgt0) then  
      a1=cbmf/0.5
   else
      a1=0
   end if

!   erfarg=wcrit / (1.4142 * sigmaw)
!   if(erfarg.lt.20.)then
!      ufrc = min(rmaxfrac, 0.5*erfccc(erfarg))
!   else
!      ufrc = 0.
!   endif

 end subroutine cu_clo_miz



!######################################################################
!######################################################################
!
! Alternative closure for Donner lite convection. This subroutine
! closes convective cloud base area using numerical approximations.
! It attempts to find the cloud base area (a1) that will produce
! the desired change in CAPE after the convective tendencies are
! applied.
!
! This subroutine works for both CAPE relaxation (do_dcape = .false.)
! and instanteneous CAPE adjustment (do_dcape = .true.) closures,
! as well as all available options for CAPE calculations routines
! controlled using namelist parameters do_donner_cape and
! do_hires_cape_for_closure
!
! do_donner_cape  do_hires_cape_for_closure  CAPE algorithm
! .false.         .false.                    UW
! .true.          .true.                     High-res Donner routine 
!                                             for trigger and closure
! .false.         .true.                     High-res Donner routine 
!                                             for closure only
! 
! The numerical algorithm for determining a1 can be best described
! as a "guided" bissection search. It is a bissection search except
! for the first three "guided" steps:
!
! Step 1:  If needed, compute environmental CAPE.
! Step 2:  Compute CAPE for an arbitrary small cloud base area 0.001.
!          Compute approximate cloud base area using linear approximation
!          and CAPE values from steps 1 and 2.
! Step 3:  Compute CAPE for linearly interpolated approximate cloud area.
!          Identify two values of cloud base areas that bracket the target
!          CAPE value.
! Step 4+: Continue using standard bissection algorithm.
!
! The algorithm will often converge during step 3, thus requiring
! only two CAPE computations.
!
! Parameters for evaluating convergence are described below.
!

subroutine cu_clo_cjg   &
           (me, nlev_model, nlev_parcel, ntr, dt, diag_unit, debug_ijt, &
            Initialized, &
            Param, Nml, amax, cape, dcape, cape_p, env_t, env_r, qt_v, qr_v, &
            pfull, zfull, phalf, zhalf, tracers, &
            land, pblht, tkemiz, qstar, coldT, sd, Uw_p, ac, &
            a1, ermesg, error )

use donner_types_mod,     only : donner_param_type, donner_nml_type, &
                                 donner_initialized_type
use conv_utilities_k_mod, only : pack_sd_lsm_k, extend_sd_k,  &
                                 adi_cloud_k, qt_parcel_k, sounding, & 
                                 adicloud, uw_params

implicit none

!----------------------------------------------------------------------
!   calling arguments

! Input

integer,                          intent(in) :: me            ! pe number
integer,                          intent(in) :: nlev_model    ! # model levels
integer,                          intent(in) :: nlev_parcel   ! # levels for CAPE
integer,                          intent(in) :: ntr           ! number of tracers
real,                             intent(in) :: dt            ! time step
integer,                          intent(in) :: diag_unit     ! column diagnostics
logical,                          intent(in) :: debug_ijt

type(donner_initialized_type),    intent(in) :: Initialized   ! Donner parameters
type(donner_param_type),          intent(in) :: Param         ! Donner parameters
type(donner_nml_type),            intent(in) :: Nml           ! Donner namelist
real,                             intent(in) :: amax          ! Max allowable cloud base
real,                             intent(in) :: cape          ! Environment CAPE
real,                             intent(in) :: dcape         ! CAPE change for dcape closure
real,  dimension(nlev_parcel),    intent(in) :: cape_p        ! p levels of CAPE profiles
real,  dimension(nlev_parcel),    intent(in) :: env_t, env_r  ! T, r CAPE profiles
real,  dimension(nlev_parcel),    intent(in) :: qt_v,  qr_v   ! Normalized cu tendencies

real, dimension(nlev_model),      intent(in) :: pfull, zfull  ! Full model levels (p and z)
real, dimension(nlev_model+1),    intent(in) :: phalf, zhalf  ! Half model levels (p and z)
real, dimension(nlev_model,ntr),  intent(in) :: tracers       ! Tracer array

! Variables needed for computing CAPE using UW formulation
real,                             intent(in) :: land, pblht, tkemiz, qstar
logical,                          intent(in) :: coldT
type(sounding),                intent(inout) :: sd
type(uw_params),               intent(inout) :: Uw_p
type(adicloud),                intent(inout) :: ac

! Output

real,                            intent(out) :: a1             ! Output cloud base area
character(len=*),                intent(out) :: ermesg         ! Output error status
integer,                         intent(out) :: error

!----------------------------------------------------------------------
! Parameters controlling iterations
!
! - nitermax is the maximum number of iterations allowed (i.e. the 
!   maximum number of CAPE calculations to be performed)
! - accuracy on a1 is controlled by a1_abs_acc. The bissection
!   algorithm will stop once the bracket interval becomes smaller
!   than a1_abs_acc
! - accuracy of the final CAPE value is controlled by cape_rel_acc,
!   which is the relative accuracy desired.
!

integer, parameter           :: nitermax = 10          ! max number of iterations
real, parameter              :: a1_abs_acc = 0.0001    ! absolute a1 accuracy
real, parameter              :: cape_rel_acc = 0.01    ! relative CAPE accuracy 

!----------------------------------------------------------------------
!   local variables

real, dimension(nitermax+1)  :: xa1, xcape

real                         :: x1, x2, xnew
real                         :: f1, f2, fnew
real                         :: tmp

real                         :: cape_target
real                         :: cape_target_acc

integer                      :: k, n, niter, nfirst, exit_flag
real                         :: plfc, plzb, plcl, coin
real, dimension(nlev_parcel) :: parcel_t, parcel_r

real, dimension(nlev_model)  :: tmp_t, tmp_r
real, dimension(nlev_parcel) :: updated_t, updated_r

real                         :: lofactor

    ermesg = ' '; error = 0

! ------------------------------------------------------------------------------
! Initialization

    exit_flag = 0
    xa1   = -999.0
    xcape = -999.0

!   The first iteration has zero cloud base area and corresponds to 
!   the environmental CAPE.

    xa1(1) = 0.0

!   Input calling argument 'cape' contains CAPE value of the environment.
!   We don't need to recompute it except when
!    the conditions in the if statement apply
!   because in that case, CAPE of the environment has been computed using
!   a different algorithm.


    if ( (Nml%do_donner_cape .EQV. .false.  &                        
         .and. Nml%do_hires_cape_for_closure .EQV. .true. ) .or.  &
    Nml%do_freezing_for_cape .NEQV. Nml%do_freezing_for_closure .or. &
           Nml%tfre_for_cape /=     Nml%tfre_for_closure .or. &
           Nml%dfre_for_cape /=     Nml%dfre_for_closure .or. &
           .not. (Initialized%use_constant_rmuz_for_closure) .or.  &
          Nml%rmuz_for_cape /= Nml%rmuz_for_closure) then
      nfirst = 1

    else

!     Environmental CAPE is known
      nfirst = 2
      xcape(1) = cape

!     Target CAPE value and accuracy for closure
      if (Nml%do_dcape) then
        cape_target = xcape(1) - dcape * dt
      else
        cape_target = xcape(1) - (xcape(1)-Nml%cape0)/Nml%tau * dt
      end if
      cape_target_acc = cape_rel_acc * cape_target

!     Next iteration
      xa1(2) = 0.001

    end if

    do n=nfirst,nitermax
      
      xnew = xa1(n)

!     Update profile and compute its CAPE

      updated_t = env_t + xnew * qt_v * dt
      updated_r = env_r + xnew * qr_v * dt

      if ( Nml%do_donner_cape .or. Nml%do_hires_cape_for_closure ) then

!       Using Donner subroutine

        call don_c_displace_parcel_k  &
             ( nlev_parcel, diag_unit, debug_ijt, Param, &
               Nml%do_freezing_for_closure, Nml%tfre_for_closure, Nml%dfre_for_closure, &
               Nml%rmuz_for_closure, &
               Initialized%use_constant_rmuz_for_closure, &
               Nml%modify_closure_plume_condensate, &
               Nml%closure_plume_condensate, &
               updated_t, updated_r, &
               cape_p, .true., plfc, plzb, plcl, coin, fnew, &
               parcel_r, parcel_t, ermesg, error )
  
        if (error /= 0 ) return

      else

!       Using UW code

        do k=1,nlev_model
          tmp_t(k) = updated_t(nlev_model-k+1)
          tmp_r(k) = updated_r(nlev_model-k+1)
        end do
        call pack_sd_lsm_k &
             ( Nml%do_lands, land, coldT, dt, pfull, phalf, zfull, zhalf, &
               tmp_t, tmp_r, tracers, sd )
        call extend_sd_k (sd, pblht, .false., Uw_p) 
        if (Nml%do_lands) then
           call qt_parcel_k (sd%qs(1), qstar, pblht, &
                             tkemiz, sd%land, &
                             Nml%gama, Nml%pblht0, Nml%tke0, Nml%lofactor0, &
                             Nml%lochoice, sd%qct(1), lofactor)
        end if
        call adi_cloud_k (sd%zs(1), sd%ps(1), sd%hl(1), sd%thc(1), sd%qct(1), sd, &
                          Uw_p, .false., Nml%do_freezing_for_cape, ac)
        fnew=ac%cape

      end if

!     Check for convergence and/or select value for next iteration

      xcape(n) = fnew

      if ( n == 1 ) then       ! First iteration

!       Target CAPE value and accuracy for closure
        if (Nml%do_dcape) then
          cape_target = xcape(1) - dcape * dt
        else
          cape_target = xcape(1) - (xcape(1)-Nml%cape0)/Nml%tau * dt
        end if
        cape_target_acc = cape_rel_acc * cape_target

!       Environmental CAPE already below target
        if ( xcape(1) <= cape_target ) then
          exit_flag = 1
          exit
        end if

!       Next iteration
        xa1(2) = 0.001

      else if ( n == 2 ) then  ! Second iteration

!       Re-order 1 and 2 to ensure xcape(1) >= xcape(2)
        if ( xcape(1) < xcape(2) ) then

          tmp = xa1(2)
          xa1(2) = xa1(1)
          xa1(1) = tmp

          tmp = xcape(2)
          xcape(2) = xcape(1)
          xcape(1) = tmp

        end if

!       Linear approximation for next iteration
        xa1(3)  &
        = max( 0.0,  &
               min( xa1(1) + (cape_target-xcape(1))/(xcape(2)-xcape(1))*(xa1(2)-xa1(1)) &
                  , amax ) &
             )

      else if ( n == 3 ) then  ! Third iteration

!       Exit iteration if convergence criteria for CAPE is met
        if ( abs(fnew-cape_target) .lt. cape_target_acc ) then
          exit_flag = 2
          exit
        end if

!       Exit iteration if desired CAPE value cannot be reached without
!       exceeding maximum cloud base area
        if ( xnew >= amax .and. fnew >= cape_target ) then
          exit_flag = 3
          exit
        end if

!       Re-order 2 and 3 so that xcape(1) >= xcape(2) >= xcape(3)
        if ( xcape(2) < xcape(3) ) then

          tmp = xa1(3)
          xa1(3) = xa1(2)
          xa1(2) = tmp

          tmp = xcape(3)
          xcape(3) = xcape(2)
          xcape(2) = tmp

        end if

!       Find where cape_target falls and decide what to do
        if ( cape_target >= xcape(2) ) then
           x1 = xa1(1)
           x2 = xa1(2)
           f1 = xcape(1)
           f2 = xcape(2)
        else if ( cape_target >= xcape(3) ) then
           x1 = xa1(2)
           x2 = xa1(3)
           f1 = xcape(2)
           f2 = xcape(3)
        else
           x1 = xa1(3)
           x2 = 2*amax - x1
           f1 = xcape(3)
           f2 = 0.0
        end if
        
        xa1(n+1) = 0.5*(x1+x2)

      else    ! Beyond third iteration, continue using bissection algorithm

!       Exit iteration if convergence criteria for CAPE is met
        if ( abs(fnew-cape_target) < cape_target_acc ) then
          exit_flag = 4
          exit
        end if

!       Exit iteration if desired CAPE value cannot be reached without
!       exceeding maximum cloud base area
        if ( xnew >= amax .and. fnew >= cape_target ) then
          exit_flag = 5
          exit
        end if

!       Select bracket for next iteration
        if ( fnew < cape_target ) then
          x2 = xnew
          f2 = fnew
        else
          x1 = xnew
          f1 = fnew
        end if

!       If bracket (x1,x2) is smaller than target, choose best point
!       between x1 and x2 and exit
        if ( abs(x1-x2) < a1_abs_acc ) then
          exit_flag = 6
          exit
        end if

!       Continue to next iteration
        xa1(n+1) = 0.5*(x1+x2)

      end if

    end do
    niter = min( n, nitermax )

!   If iteration exited because bracket became smaller than 
!   threshold or because maximum number of iterations was reached,
!   select the best solution between f1 and f2
    if ( exit_flag == 0 .or. exit_flag == 6 ) then

      if ( abs(f1-cape_target) < abs(f2-cape_target) ) then
        xnew = x1
        fnew = f1
      else
        xnew = x2
        fnew = f2
      end if

    end if

!   Copy result into output variable
    a1 = xnew

! ------------------------------------------------------------------------------
! Output debugging information

    if (debug_ijt) then

      write (diag_unit, '(a)' ) &
        ' --- begin cu_clo_cjg debug info ---------------------------------------'
!     write (diag_unit, '(a)')  & 
!       '                      cu_clo_cjg input profile'
!     write (diag_unit, '(a)')  &
!       '   k   press      temp           mixing ratio '
!     write (diag_unit, '(a)')  &
!       '        hPa       deg K    g(h2o) / kg (dry air)  '
!     do k=1,nlev_parcel
!       write (diag_unit, '(i4, 2f10.4, 7x, 1pe13.5)')  &
!            k, 1.0E-02*cape_p(k), env_t(k), 1.0e03*env_r(k)
!     end do

!     write (diag_unit, '(a)')  &
!       '    k, tmp_t(k), 1.0e3*tmp_r(k)'
!     do k=1,nlev_model
!       write (diag_unit, '(i4, 2f10.4)')  &
!            k, tmp_t(k), 1.0e3*tmp_r(k)
!     end do

!     write (diag_unit, '(a)') '  '
!     write (diag_unit, '(a, f20.12)') &
!       'cjg: xcape0a = ',xcape0a
!     write (diag_unit, '(a, f20.12)') &
!       'cjg: xcape0b = ',xcape0b

!     write (diag_unit, '(a)') '  '
!     write (diag_unit, '(a)')  & 
!       '                      cu_clo_cjg updated profile'
!     write (diag_unit, '(a)')  &
!       '   k   press      temp           mixing ratio '
!     write (diag_unit, '(a)')  &
!       '        hPa       deg K    g(h2o) / kg (dry air)  '
!     do k=1,nlev_parcel
!       write (diag_unit, '(i4, 2f10.4, 7x, 1pe13.5)')  &
!            k, 1.0E-02*cape_p(k), updated_t(k), 1.0e03*updated_r(k)
!     end do

!     write (diag_unit, '(a)') '  '
!     write (diag_unit, '(a, f20.12)') &
!       'cjg: xcape1a = ',xcape1a
!     write (diag_unit, '(a, f20.12)') &
!       'cjg: xcape1b = ',xcape1b

      write (diag_unit, '(a)') '  '
      write (diag_unit, '(a, f10.4)') &
        'cjg: amax = ',amax
      write (diag_unit, '(a, f10.4)') &
        'cjg: cape = ',cape
      write (diag_unit, '(a, f10.4)') &
        'cjg: cape_target = ',cape_target
      do n=1,nitermax
        write (diag_unit, '(a, i4, 2f10.4)')  &
             'cjg: ', n, xa1(n), xcape(n)
      end do
      write (diag_unit, '(a)' ) &
        ' --- end cu_clo_cjg debug info -----------------------------------------'

    end if

!RSH   debug_unit = 1000+mpp_pe()
!   debug_unit = 1000+me
!   write (debug_unit, '(a)') '  '
!   write (debug_unit, '(a, f10.4)') &
!     'cjg: amax = ',amax
!   write (debug_unit, '(a, f10.4)') &
!     'cjg: cape = ',cape
!   write (debug_unit, '(a, f10.4)') &
!     'cjg: cape_target = ',cape_target
!   do n=1,nitermax
!     write (debug_unit, '(a, i4, 2f10.4)')  &
!          'cjg: ', n, xa1(n), xcape(n)
!   end do
!   write (debug_unit, '(a, 2f10.4)') &
!     'cjg: final a1, cape = ',xnew,fnew
!   write (debug_unit, '(a, i4)') 'cjg: iterations = ',niter
!   write (debug_unit, '(a, i4)') 'cjg: exit_flag = ',exit_flag
!   do n=nfirst,niter
!     write (debug_unit, '(a)') 'cjg: niter'
!   end do

end subroutine cu_clo_cjg

!######################################################################
!######################################################################




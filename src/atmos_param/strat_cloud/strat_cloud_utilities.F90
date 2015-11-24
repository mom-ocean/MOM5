module strat_cloud_utilities_mod

use  fms_mod,        only : write_version_number

IMPLICIT NONE
PRIVATE

!-------------------------------------------------------------------------
!--interfaces-------------------------------------------------------------

public strat_cloud_utilities_init

PUBLIC diag_id_type, diag_pt_type, strat_nml_type,  &
       atmos_state_type, particles_type, cloud_state_type, &
       precip_state_type, cloud_processes_type, strat_constants_type

!----------------------------------------------------------------------
!----version number----------------------------------------------------

Character(len=128) :: Version = '$Id: strat_cloud_utilities.F90,v 20.0 2013/12/13 23:22:15 fms Exp $'
Character(len=128) :: Tagname = '$Name: tikal $'

logical  :: module_is_initialized = .false.


!########################################################################

TYPE diag_id_type

!  cloud area variables

  integer :: aall, aliq, aice, cf_liq_init, cf_ice_init, aauto
  integer :: SA3d, qadt_lsform, qadt_lsdiss, qadt_rhred, qadt_eros,  &
             qadt_fill, qadt_super, qadt_destr, qadt_limits, qadt_ahuco, &
             SA_imb
  integer :: SA2d, qa_lsform_col, qa_lsdiss_col, qa_rhred_col,  &
             qa_eros_col, qa_fill_col, qa_super_col, qa_destr_col,     &
             qa_limits_col, qa_ahuco_col, SA_imb_col

!  cloud liquid variables

  integer :: SL3d, qldt_cond, qldt_evap, qldt_eros, qldt_berg, qldt_freez,&
             liq_adj, qldt_rime, qldt_accr, qldt_auto, qldt_fill,  &
             qldt_destr, qldt_freez2, qldt_sedi, qldt_accrs, qldt_bergs, &
             qldt_HM_splinter, SL_imb
  integer :: SL2d, ql_cond_col, ql_evap_col, ql_eros_col, ql_berg_col,   &
             ql_freez_col, liq_adj_col, ql_rime_col, ql_accr_col,  &
             ql_auto_col, ql_fill_col, ql_destr_col, ql_freez2_col,  &
             ql_sedi_col, ql_accrs_col, ql_bergs_col, ql_HM_splinter_col, &
             SL_imb_col

!  cloud ice variables

  integer :: SI3d, qidt_dep, qidt_subl, qidt_fall, qidt_eros, qidt_melt, &
             qidt_melt2, qidt_fill, qidt_destr, qidt_qvdep, qidt_auto,  &
             qidt_accr, qidt_accrs, ice_adj,  SI_imb
  integer :: SI2d, qi_dep_col, qi_subl_col, qi_fall_col, qi_eros_col,  &    
             qi_melt_col, qi_melt2_col, qi_fill_col, qi_destr_col,  &
             qi_qvdep_col, qi_auto_col, qi_accr_col, qi_accrs_col,  &
             ice_adj_col, SI_imb_col  

!  cloud droplet variables

  integer :: droplets_col250, gb_droplets_col, potential_droplets, &
             droplets, droplets_wtd, ql_wt, droplets_col, rvolume
  integer :: SN3d, qndt_cond , qndt_evap, qndt_fill, qndt_berg, &
             qndt_destr, qndt_super, qndt_freez, qndt_sacws, qndt_sacws_o, &
             qndt_eros, qndt_pra, qndt_auto, qndt_nucclim, qndt_sedi, &
             qndt_melt, qndt_ihom, qndt_size_adj, qndt_fill2,   &
             qndt_contact_frz, qndt_cleanup, qndt_cleanup2, SN_imb
  integer :: SN2d, qn_cond_col, qn_evap_col, qn_fill_col, qn_berg_col,  &
             qn_destr_col, qn_super_col, qn_freez_col, qn_sacws_col,  &
             qn_sacws_o_col, qn_eros_col, qn_pra_col, qn_auto_col,    &
             qn_nucclim_col, qn_sedi_col, qn_melt_col, qn_ihom_col,  &
             qn_size_adj_col, qn_fill2_col, qn_contact_frz_col,  &
             qn_cleanup_col, qn_cleanup2_col, SN_imb_col

!  cloud ice particle variables

  integer :: nice, nice_col, gb_nice_col, potential_crystals
  integer :: SNi3d, qnidt_fill, qnidt_nnuccd, qnidt_nsubi,  &
             qnidt_nerosi, qnidt_nprci, qnidt_nprai,                 &
             qnidt_nucclim1, qnidt_nucclim2, qnidt_sedi,             &
             qnidt_melt, qnidt_size_adj, qnidt_fill2,                &
             qnidt_super, qnidt_ihom, qnidt_destr,                   &
             qnidt_cleanup, qnidt_cleanup2, qnidt_nsacwi, SNi_imb  
  integer :: SNi2d, qni_fill_col, qni_nnuccd_col, qni_nsubi_col, &
             qni_nerosi_col, qni_nprci_col, qni_nprai_col, &
             qni_nucclim1_col, qni_nucclim2_col, qni_sedi_col, &
             qni_melt_col, qni_size_adj_col, qni_fill2_col, &
             qni_super_col, qni_ihom_col, qni_destr_col, &
             qni_cleanup_col, qni_cleanup2_col,                    &
             qni_nsacwi_col, SNi_imb_col

!  aerosol diagnostics

  integer :: delta_cf, sulfate, seasalt_sub, seasalt_sup, om,   &
             rhcrit, rhcrit_min, rhiin, rhlin, cfin, imass7,     &
             ni_dust, ni_sulf, ni_bc, ndust1, ndust2, ndust3,  &
             ndust4, ndust5, dust_berg_flag, subgrid_w_variance
 
!  rain diagnostics

  integer :: rain3d, qrout, rain_clr, rain_cld, a_rain_clr, a_rain_cld, &
             rain_evap, rain_freeze, srfrain_accrs, srfrain_freez,  &
             srfrain_evap, rain_evap_col, rain_freeze_col,  &
             srfrain_accrs_col, srfrain_freez_col, srfrain_evap_col, &
             rain_mass_conv, rain_imb, rain_imb_col, cld_liq_imb,  &
             cld_liq_imb_col, neg_rain, qrout_col

!  snow diagnostics

  integer :: snow3d, qsout, snow_clr, snow_cld, a_snow_clr, a_snow_cld, &
             snow_melt, snow_melt_col, snow_mass_conv, sedi_ice, snow_imb, &
             snow_imb_col, cld_ice_imb, cld_ice_imb_col, neg_snow, qsout_col
             

!  total precip diagnostics

  integer :: a_precip_cld, a_precip_clr, sedi_sfc

!  temperature diagnostics

  integer ::  ST3d, ST_imb 
  integer ::  ST2d, ST_imb_col

!  vapor diagnostics

  integer :: SQ3d, qdt_liquid_init, qdt_ice_init, qdt_rain_evap,   &
             qdt_cond, qdt_deposition, qdt_eros_l, qdt_eros_i,        &
             qdt_qv_on_qi, qdt_sedi_ice2vapor, qdt_sedi_liquid2vapor,  &
             qdt_super_sat_rm, qdt_destr, qdt_cleanup_liquid,  &
             qdt_cleanup_ice, qdt_snow_sublim, qdt_snow2vapor, SQ_imb  
  integer :: SQ2d, q_liquid_init_col, q_ice_init_col, q_rain_evap_col, &
             q_cond_col, q_deposition_col, q_eros_l_col, q_eros_i_col, &
             q_qv_on_qi_col, q_sedi_ice2vapor_col, q_sedi_liquid2vapor_col,&
             q_super_sat_rm_col, q_destr_col, q_cleanup_liquid_col, &
             q_cleanup_ice_col, q_snow_sublim_col, q_snow2vapor_col,  &
             SQ_imb_col

!   miscellaneous diagnostics

  integer :: f_snow_berg, f_snow_berg_col, &
             lsf_strat, lcf_strat, mfls_strat, &
             dcond, vfall


END TYPE diag_id_type


!#########################################################################

TYPE diag_pt_type

!  cloud area variables

  integer :: aall, aliq, aice, cf_liq_init, cf_ice_init, aauto
  integer :: SA3d, qadt_lsform, qadt_lsdiss, qadt_rhred, qadt_eros,  &
             qadt_fill, qadt_super, qadt_destr, qadt_limits, qadt_ahuco, &
             SA_imb

!  cloud liquid variables

  integer :: SL3d, qldt_cond, qldt_evap, qldt_eros, qldt_berg, qldt_freez,&
             liq_adj, qldt_rime, qldt_accr, qldt_auto, qldt_fill,  &
             qldt_destr, qldt_freez2, qldt_sedi, qldt_accrs, qldt_bergs, &
             qldt_HM_splinter, SL_imb

!  cloud ice variables

  integer :: SI3d, qidt_dep, qidt_subl, qidt_fall, qidt_eros, qidt_melt, &
             qidt_melt2, qidt_fill, qidt_destr, qidt_qvdep, qidt_auto,  &
             qidt_accr, qidt_accrs, ice_adj,  SI_imb

!  cloud droplet variables

  integer :: droplets_col250, gb_droplets_col, potential_droplets, &
             droplets, droplets_wtd, ql_wt, droplets_col, rvolume
  integer :: SN3d, qndt_cond , qndt_evap, qndt_fill, qndt_berg, &
             qndt_destr, qndt_super, qndt_freez, qndt_sacws, qndt_sacws_o, &
             qndt_eros, qndt_pra, qndt_auto, qndt_nucclim, qndt_sedi, &
             qndt_melt, qndt_ihom, qndt_size_adj, qndt_fill2,   &
             qndt_contact_frz, qndt_cleanup, qndt_cleanup2, SN_imb

!  cloud ice particle variables

  integer :: nice, nice_col, gb_nice_col, potential_crystals
  integer :: SNi3d, qnidt_fill, qnidt_nnuccd, qnidt_nsubi,  &
             qnidt_nerosi, qnidt_nprci, qnidt_nprai,                 &
             qnidt_nucclim1, qnidt_nucclim2, qnidt_sedi,             &
             qnidt_melt, qnidt_size_adj, qnidt_fill2,                &
             qnidt_super, qnidt_ihom, qnidt_destr,                   &
             qnidt_cleanup, qnidt_cleanup2, qnidt_nsacwi, SNi_imb  

!  aerosol diagnostics

  integer :: delta_cf, sulfate, seasalt_sub, seasalt_sup, om,   &
             rhcrit, rhcrit_min, rhiin, rhlin, cfin, imass7,     &
             ni_dust, ni_sulf, ni_bc, ndust1, ndust2, ndust3,  &
             ndust4, ndust5, dust_berg_flag, subgrid_w_variance
 
!  rain diagnostics

  integer :: rain3d, qrout, rain_clr, rain_cld, a_rain_clr, a_rain_cld, &
             rain_evap, rain_freeze, srfrain_accrs, srfrain_freez,  &
             srfrain_evap, rain_mass_conv, rain_imb, cld_liq_imb, neg_rain

!  snow diagnostics

  integer :: snow3d, qsout, snow_clr, snow_cld, a_snow_clr, a_snow_cld, &
             snow_melt, snow_mass_conv, sedi_ice, snow_imb, cld_ice_imb, &
             neg_snow

!  total precip diagnostics

  integer :: a_precip_cld, a_precip_clr, sedi_sfc

!  temperature diagnostics

  integer ::  ST3d, ST_imb 

!  vapor diagnostics

  integer :: SQ3d, qdt_liquid_init, qdt_ice_init, qdt_rain_evap,   &
             qdt_cond, qdt_deposition, qdt_eros_l, qdt_eros_i,        &
             qdt_qv_on_qi, qdt_sedi_ice2vapor, qdt_sedi_liquid2vapor,  &
             qdt_super_sat_rm, qdt_destr, qdt_cleanup_liquid,  &
             qdt_cleanup_ice, qdt_snow_sublim, qdt_snow2vapor, SQ_imb  

!   miscellaneous diagnostics

  integer :: f_snow_berg, &
             lsf_strat, lcf_strat, mfls_strat, &
             dcond, vfall

END TYPE diag_pt_type


!#########################################################################

type strat_nml_type

!-------------------------------------------------------------------------
!    see strat_nml.h for a description of these variables.
!-------------------------------------------------------------------------

  real    :: U00, rthresh, var_limit, sea_salt_scale,         &
             om_to_oc,  N_land, N_ocean, U_evap, eros_scale,  &
             eros_scale_c, eros_scale_t, mc_thresh,           &
             diff_thresh, qmin, Dmin, efact, vfact, cfact,    &
             iwc_crit,  vfall_const2, vfall_exp2,             &
             qthalfwidth, N_min, num_mass_ratio1,             &
             num_mass_ratio2, qcvar

  logical :: do_netcdf_restart, u00_profile, use_kk_auto,     &
             use_online_aerosol, use_sub_seasalt,             &
             eros_choice, super_choice, tracer_advec,         &
             do_old_snowmelt, retain_cm3_bug, do_pdf_clouds, do_liq_num,      &
             do_dust_berg, pdf_org, do_ice_nucl_wpdf, debugo, &
             mass_cons, do_hallet_mossop, activate_all_ice_always

  integer :: num_strat_pts, betaP, nsublevels, kmap, kord,    &
             super_ice_opt, isamp, jsamp, ksamp

  character(len=64)                :: microphys_scheme, &
                                      macrophys_scheme, &
                                      aerosol_activation_scheme

  integer, dimension(:,:), pointer :: strat_pts=>NULL()

end type strat_nml_type


!#########################################################################

type atmos_state_type

!       airdens        air density                     kg air/(m*m*m)
!       qs             saturation specific humidity    kg vapor/kg air
!       dqsdT          T derivative of qs              kg vapor/kg air/K
!       gamma          (L/cp)*dqsdT                    dimensionless
!       deltpg         pressure thickness of grid box  kg air/(m*m)
!                      divided by gravity
!       U_ca           grid box relative humidity      fraction

  real, dimension(:,:,:), pointer ::  &
                                        pfull          =>NULL(), &
                                        phalf          =>NULL(), &
                                        zhalf          =>NULL(), &
                                        zfull          =>NULL(), &
                                        radturbten2    =>NULL(), &
                                        T_in           =>NULL(), &
                                        qv_in          =>NULL(), &
                                        omega          =>NULL(), &
                                        Mc             =>NULL(), &
                                        diff_t         =>NULL(), &
                                        qrat           =>NULL(), &
                                        airdens        =>NULL(), &
                                        tn             =>NULL(), & 
                                        qvn            =>NULL(), &
                                        qs             =>NULL(), &
                                        dqsdT          =>NULL(), &
                                        qsi            =>NULL(), &
                                        qsl            =>NULL(), &
                                        rh_crit        =>NULL(), &
                                        rh_crit_min    =>NULL(), &
                                        gamma          =>NULL(), &
                                        esat0          =>NULL(), &
                                        U_ca           =>NULL(), &
                                        ahuco          =>NULL(), &
                                        delp           =>NULL(), &
                                        U01            =>NULL(), &
                                        deltpg         =>NULL()

end type atmos_state_type


!##########################################################################

type  particles_type 

! drop1           number conc                     [1/cm^3]
! drop2           mass concentration              [1/kg]
  real, dimension(:,:,:), pointer ::  &
                                        concen_dust_sub=>NULL(), &
                                        drop1          =>NULL(), &
                                        drop2          =>NULL(), &
                                        crystal1       =>NULL(), &
                                        rbar_dust      =>NULL(), &
                                        ndust          =>NULL(), &
                                        hom            =>NULL()

end type particles_type 


!########################################################################

type cloud_state_type

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
!       ql_mean        ql + positive increment         kg condensate/
!                      of ql; i.e. a sort of           kg air
!                      intermediate ql
!
!       qi_mean        ql + positive increment         kg condensate/
!                      of qi; i.e. a sort of           kg air
!                      intermediate qi
  real, dimension(:,:,:), pointer ::   &
                                        ql_upd         =>NULL(), &
                                        qi_upd         =>NULL(), &
                                        qa_upd         =>NULL(), &
                                        qn_upd         =>NULL(), &
                                        qni_upd        =>NULL(), &
                                        ql_mean        =>NULL(), &
                                        qi_mean        =>NULL(), &
                                        qa_mean        =>NULL(), &
                                        qn_mean        =>NULL(), &
                                        qni_mean       =>NULL(), &
                                        ql_in          =>NULL(), &
                                        qi_in          =>NULL(), &
                                        qa_in          =>NULL(), &
                                        qn_in          =>NULL(), &
                                        qni_in         =>NULL(), &
                                        SL_out         =>NULL(), &
                                        SI_out         =>NULL(), &
                                        SA_out         =>NULL(), &
                                        SN_out         =>NULL(), &
                                        SNi_out        =>NULL(), &
                                        SA_0           =>NULL(), &
                                        qa_upd_0       =>NULL()

end type cloud_state_type
 

!#########################################################################

type precip_state_type

  real, dimension(:,:,:), pointer ::   &
                                        lsc_snow       =>NULL(), &
                                        lsc_rain       =>NULL(), &
                                        lsc_snow_size  =>NULL(), &
                                        lsc_rain_size  =>NULL(), &
!rain and snow mixing ratios from the Morrison scheme (kg/kg) :
                                        qrout3d_mg     =>NULL(), &
                                        qsout3d_mg     =>NULL(), &
                                        rain3d         =>NULL(), &
                                        snow3d         =>NULL(), &
                                        snowclr3d      =>NULL()

  real, dimension(:,:), pointer   ::   &
                                        surfrain       =>NULL(), &
                                        surfsnow       =>NULL()
             
end type precip_state_type


!#########################################################################

type cloud_processes_type

!       da_ls          change in saturated volume      fraction
!                      fraction due to large-scale
!                      processes
!       D_eros         Sink in ql, qi and qa equation  dimensionless
!                      due to turbulent erosion of
!                      cloud sides
!       dcond_ls       change in condensate due to     kg condensate/
!                      non-convective condensation.    kg air
!                      After phase determination,
!                      this variable refers only to
!                      liquid condensation.
!
!       dcond_ls_ice   change in ice due to            kg condensate/
!                      non-convective condensation.    kg air
!       qvg            equilibrium value of water      kg vapor /
!                      vapor in the clear portion      kg air
!                      of the grid box that PDF 
!                      clouds wants                          
!

  real, dimension(:,:,:), pointer ::   &
                                        da_ls          =>NULL(), &
                                        D_eros         =>NULL(), &
                                        qvg            =>NULL(), &
                                        dcond_ls       =>NULL(), &
                                        dcond_ls_ice   =>NULL(), &
                                        dcond_ls_tot   =>NULL(), &
                                        delta_cf       =>NULL(), &
                                        f_snow_berg    =>NULL()

end type cloud_processes_type


!##########################################################################

type strat_constants_type

!       inv_dtcloud   inverse of model timestep [ sec (-1) ]
  real, dimension(:,:,:), pointer :: mask=>NULL()
  real                            :: dtcloud, inv_dtcloud
  integer                         :: overlap
  logical                         :: limit_conv_cloud_frac,      &
                                     mask_present,               &
                                     do_rk_microphys,            &
                                     do_mg_microphys,            &       
                                     do_mg_ncar_microphys,       &
                                     tiedtke_macrophysics,       &
                                     dqa_activation,             &
                                     total_activation,           &
                                     do_predicted_ice_number

end type strat_constants_type


!#########################################################################


CONTAINS

subroutine strat_cloud_utilities_init

      if (module_is_initialized) return

!------------------------------------------------------------------------
!    write version number to output file.
!------------------------------------------------------------------------
      call write_version_number(version, tagname)

!------------------------------------------------------------------------
!    mark this module as initialized.
!------------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------------

end subroutine strat_cloud_utilities_init



!########################################################################


end module strat_cloud_utilities_mod

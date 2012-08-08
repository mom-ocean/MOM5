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

Character(len=128) :: Version = '$Id: strat_cloud_utilities.F90,v 19.0 2012/01/06 20:27:23 fms Exp $'
Character(len=128) :: Tagname = '$Name: siena_201207 $'

logical  :: module_is_initialized = .false.


!########################################################################

TYPE diag_id_type

  integer                                                     &
      qldt_cond, ql_cond_col, qldt_evap, ql_evap_col,         &
      qldt_eros, ql_eros_col, qldt_berg, ql_berg_col,         &
      qldt_freez, ql_freez_col, qldt_rime, ql_rime_col,       &
      qldt_accr, ql_accr_col, qldt_auto, ql_auto_col,         &
      qldt_fill, ql_fill_col, qldt_destr, ql_destr_col,       &
      liq_adj, liq_adj_col, ice_adj, ice_adj_col,             &
      snow_melt, snow_melt_col, snow_subl, snow_subl_col,     &
      qidt_dep, qi_dep_col, qidt_eros, qi_eros_col,           &
      qidt_fall, qi_fall_col, qidt_fill, qi_fill_col,         &
      qidt_subl, qi_subl_col,  qidt_melt, qi_melt_col,        &
      qidt_destr, qi_destr_col, qidt_qvdep,  qi_qvdep_col,    &
      qadt_lsform, qa_lsform_col,                             &
      qadt_lsdiss, qa_lsdiss_col,  qadt_eros, qa_eros_col,    &
      qadt_rhred, qa_rhred_col, qadt_destr, qa_destr_col,     &
      qadt_fill, qa_fill_col,  qadt_super, qa_super_col,      &
      qndt_evap, qn_evap_col, qndt_fill, qn_fill_col,         &
      qndt_destr, qn_destr_col, qndt_super, qn_super_col,     &
      qldt_freez2, ql_freez2_col, qldt_sedi, ql_sedi_col,     &
      qldt_accrs, ql_accrs_col, qldt_bergs, ql_bergs_col,     &
      qidt_auto, qi_auto_col, qidt_accr, qi_accr_col,         &
      qidt_accrs,  qi_accrs_col, qndt_cond , qn_cond_col,     &
      rain_evap, rain_evap_col, debug1_3d, debug2_3d,         &
      debug3_3d, debug4_3d, debug5_3d, tmp5_3d,               &
      qndt_freez, qndt_sacws, qndt_sacws_o,                   &
      qndt_eros, qndt_pra, qndt_auto,                         &
      qndt_berg,                                              &
      qn_freez_col, qn_sacws_col, qn_sacws_o_col,             &
      qn_eros_col, qn_pra_col, qn_auto_col,                   &
      qndt_nucclim, qndt_sedi, qndt_melt, qndt_ihom,          &
      qndt_size_adj, qndt_fill2,                              &  
      qn_nucclim_col, qn_sedi_col, qn_melt_col, qn_ihom_col,  &
      qn_size_adj_col, qn_fill2_col,                          &
      qn_berg_col,                                            &
      rhcrit, rhcrit_min, ni_dust, ni_sulf, ni_bc,            &
      rhiin, rhlin, cfin, imass7,                             &
      ndust1, ndust2, ndust3, ndust4, ndust5,                 &
      qnidt_fill, qnidt_nnuccd, qnidt_nsubi,                  &
      qnidt_nerosi, qnidt_nprci, qnidt_nprai,                 &
      qnidt_nucclim1, qnidt_nucclim2, qnidt_sedi,             &
      qnidt_melt, qnidt_size_adj, qnidt_fill2,                &
      qnidt_super, qnidt_ihom, qnidt_destr,                   &
      qnidt_cleanup 

  integer :: nice, nice_col, gb_nice_col, qrout, qsout
  integer :: rain3d, snow3d
  integer :: droplets_s, droplets_col_s, droplets_col250,     &
             gb_droplets_col, sulfate, seasalt_sub, seasalt_sup, om
  integer :: aliq, aice, aall, autocv, vfall 
  integer :: rain_clr, rain_cld, a_rain_clr, a_rain_cld
  integer :: snow_clr, snow_cld, a_snow_clr, a_snow_cld
  integer :: a_precip_cld, a_precip_clr
  integer :: areaall, arealiq, dcond, areaice, &
             rvolume, areaautocv, vfalldiag
  integer :: droplets, droplets_wtd, ql_wt, droplets_col, &
             lsf_strat, lcf_strat, mfls_strat

END TYPE diag_id_type


!#########################################################################

TYPE diag_pt_type

  integer ::                                                  &
             qldt_cond, qldt_evap, qldt_berg, qldt_freez,     &
             qldt_rime, qldt_accr, qldt_auto, qldt_fill,      &
             qldt_destr, qldt_eros, liq_adj, rain_evap,       &
             qidt_dep, qidt_subl, qidt_fill, qidt_melt,       &
             qidt_fall, qidt_destr, qidt_qvdep, qidt_eros,    &
             ice_adj, snow_subl, snow_melt, qadt_lsform,      &
             qadt_eros, qadt_rhred, qadt_destr, qadt_fill,    &
             qadt_lsdiss, qadt_super, qndt_cond, qndt_evap,   &
             qndt_fill, qndt_destr, qndt_super,               &
             debug1_3d,  debug2_3d,  debug3_3d,  debug4_3d,   &
             debug5_3d, tmp5_3d, qldt_freez2, qldt_sedi,      &
             qldt_accrs, qldt_bergs, qidt_auto, qidt_accr,    &
             qidt_accrs, qndt_freez, qndt_sacws,              &
             qndt_sacws_o, qndt_eros, qndt_pra, qndt_auto,    &
             qndt_nucclim, qndt_sedi, qndt_melt, qndt_ihom,   &
             qndt_size_adj, qndt_fill2, rhcrit, rhcrit_min,   &
             ni_dust, ni_sulf, ni_bc, rhiin, rhlin, cfin,     &
             imass7, ndust1, ndust2, ndust3, ndust4, ndust5,  &
             qnidt_fill, qnidt_nnuccd, qnidt_nsubi,           &
             qnidt_nerosi, qnidt_nprci, qnidt_nprai,          &
             qnidt_nucclim1, qnidt_nucclim2, qnidt_sedi,      &
             qnidt_melt, qnidt_size_adj, qnidt_fill2,         &
             qnidt_super, qnidt_ihom, qnidt_destr,            &
             qnidt_cleanup, qndt_berg
  integer :: areaall, arealiq, dcond, areaice,                &
             rvolume, areaautocv, vfalldiag
  integer :: droplets, droplets_wtd, ql_wt, droplets_col,     &
             lsf_strat, lcf_strat, mfls_strat
  integer :: droplets_s, droplets_col_s, droplets_col250,     &
             gb_droplets_col
  integer :: nice, qrout, qsout, nice_col, gb_nice_col,       &
             rain3d, snow3d, sulfate, seasalt_sub,            &
             seasalt_sup, om, aall, aliq, aice, autocv,       &
             vfall, rain_clr, rain_cld, a_rain_clr,           &
             a_rain_cld, a_precip_clr, a_precip_cld,          &
             snow_clr, snow_cld, a_snow_cld, a_snow_clr

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
             num_mass_ratio2

  logical :: do_netcdf_restart, u00_profile, use_kk_auto,     &
             use_online_aerosol, use_sub_seasalt,             &
             eros_choice, super_choice, tracer_advec,         &
             do_old_snowmelt, do_pdf_clouds, do_liq_num,      &
             do_dust_berg, pdf_org, do_ice_nucl_wpdf, debugo

  integer :: num_strat_pts, betaP, nsublevels, kmap, kord,    &
             super_ice_opt, isamp, jsamp, ksamp

  character(len=64)                :: microphys_scheme

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
                                        tmp5           =>NULL(), &
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
                                     do_mg_microphys       

end type strat_constants_type


!#########################################################################


CONTAINS

subroutine strat_cloud_utilities_init

      if (module_is_initialized) return

!------------------------------------------------------------------------
!    write version number to output file.
!------------------------------------------------------------------------
      call write_version_number (version, tagname)

!------------------------------------------------------------------------
!    mark this module as initialized.
!------------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------------

end subroutine strat_cloud_utilities_init



!########################################################################


end module strat_cloud_utilities_mod

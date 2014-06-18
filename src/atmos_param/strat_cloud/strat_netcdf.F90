!FDOC_TAG_GFDL
module strat_netcdf_mod

use fms_mod,                   only :  write_version_number
use diag_manager_mod,          only :  register_diag_field, send_data
use time_manager_mod,          only :  time_type
use strat_cloud_utilities_mod, only :  strat_cloud_utilities_init, &
                                       diag_id_type, diag_pt_type

implicit none
private

!-----------------------------------------------------------------------
!---interfaces----------------------------------------------------------
public strat_netcdf_init, strat_netcdf, strat_netcdf_end
private diag_field_init

!----------------------------------------------------------------------
!----version number----------------------------------------------------
Character(len=128) :: Version = '$Id: strat_netcdf.F90,v 20.0 2013/12/13 23:22:17 fms Exp $'
Character(len=128) :: Tagname = '$Name: tikal $'

!-----------------------------------------------------------------------
!-------------------- diagnostics variables-----------------------------

character(len=5) :: mod_name = 'strat'
real :: missing_value = -999.


       
logical  :: module_is_initialized = .false.


CONTAINS


!#########################################################################

subroutine strat_netcdf_init (axes, Time, diag_id, diag_pt, n_diag_4d, &
                              n_diag_4d_kp1)

!------------------------------------------------------------------------
integer,             intent(in)    :: axes(4)
type(time_type),     intent(in)    :: Time
type(diag_id_type),  intent(inout) :: diag_id
type(diag_pt_type),  intent(inout) :: diag_pt
integer,             intent(out)   :: n_diag_4d, n_diag_4d_kp1

!------------------------------------------------------------------------
      if (module_is_initialized) return

!-----------------------------------------------------------------------
!    write version info to standard log.
!-----------------------------------------------------------------------
      call write_version_number ( version, tagname )

!------------------------------------------------------------------------
!    make sure needed modules have been initialized.
!------------------------------------------------------------------------
      call strat_cloud_utilities_init

!------------------------------------------------------------------------
!    call diag_field_init to initialize any desired netcdf output fields.
!------------------------------------------------------------------------
      call diag_field_init (axes, Time, diag_id, diag_pt, n_diag_4d, &
                            n_diag_4d_kp1)

!-----------------------------------------------------------------------
      module_is_initialized = .true.


end subroutine strat_netcdf_init
 

!#######################################################################

subroutine strat_netcdf (diag_id, diag_pt, diag_4d, diag_4d_kp1, &
                         diag_3d, Time, is, js, kdim, mask3d)

!------------------------------------------------------------------------
type(diag_id_type),              intent(in) :: diag_id
type(diag_pt_type),              intent(in) :: diag_pt
real,dimension(:,:,:,0:),        intent(in) :: diag_4d, diag_4d_kp1
real,dimension(:,:,0:),          intent(in) :: diag_3d
type(time_type),                 intent(in) :: Time
integer,                         intent(in) :: is, js, kdim
real, dimension(:,:,:),optional, intent(in) :: mask3d


!---local variables------------------------------------------------------
      real, dimension (size(diag_4d,1), size(diag_4d,2), kdim+1) :: mask3
      logical, dimension (size(diag_4d,1), size(diag_4d,2), kdim) :: mask4
      logical :: used

!----------------------------------------------------------------------
!    set up half level mask.
!----------------------------------------------------------------------
      mask3(:,:,1:(kdim+1)) = 1.
      if (present(mask3d)) then
        where (mask3d(:,:,1:kdim) <= 0.5)
             mask3(:,:,2:(kdim+1)) = 0.
        end where
      endif
        
!-----------------------------------------------------------------------
!
!                    3-DIMENSIONAL DIAGNOSTICS
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!    1) variables associated with droplet number and size:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%droplets, diag_4d(:,:,:,diag_pt%droplets),  &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%droplets_wtd, diag_4d(:,:,:,diag_pt%droplets_wtd),&
               Time, is, js, 1, mask=diag_4d(:,:,:,diag_pt%droplets) > 0.0)
      used = send_data    &
              (diag_id%rvolume, diag_4d(:,:,:,diag_pt%rvolume),   &
               Time, is, js, 1, rmask=mask3d)

!-----------------------------------------------------------------------
!    2) variables associated with cloud liquid content:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%ql_wt, diag_4d(:,:,:,diag_pt%ql_wt),   &
               Time, is, js, 1, mask=diag_4d(:,:,:,diag_pt%droplets) > 0.0)

!-----------------------------------------------------------------------
!    3) variables associated with cloud and precipitation processes:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%lsf_strat, diag_4d(:,:,:,diag_pt%lsf_strat),  &
               Time, is, js, 1, rmask=mask3d)
      used = send_data    &
              (diag_id%dcond, diag_4d(:,:,:,diag_pt%dcond),  &
               Time, is, js, 1, rmask=mask3d )
      used = send_data    &
              (diag_id%aauto, diag_4d(:,:,:,diag_pt%aauto),   &
               Time, is, js, 1, rmask=mask3d)
      used = send_data    &
              (diag_id%vfall, diag_4d(:,:,:,diag_pt%vfall),  &
               Time, is, js, 1, rmask=mask3d) 
      used = send_data   &
              (diag_id%delta_cf, diag_4d(:,:,:,diag_pt%delta_cf), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%cf_liq_init  , diag_4d(:,:,:,diag_pt%cf_liq_init  ), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%subgrid_w_variance,   &
               diag_4d(:,:,:,diag_pt%subgrid_w_variance), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%potential_droplets,   &
               diag_4d(:,:,:,diag_pt%potential_droplets), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%potential_crystals,  &
               diag_4d(:,:,:,diag_pt%potential_crystals), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%dust_berg_flag,    &
               diag_4d(:,:,:,diag_pt%dust_berg_flag), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%cf_ice_init, diag_4d(:,:,:,diag_pt%cf_ice_init), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qdt_snow_sublim,    &
               diag_4d(:,:,:,diag_pt%qdt_snow_sublim), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%snow_melt, diag_4d(:,:,:,diag_pt%snow_melt), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%rain_freeze, diag_4d(:,:,:,diag_pt%rain_freeze), &
               Time, is, js, 1, rmask=mask3d)

!-----------------------------------------------------------------------
!    4) variables associated with model convection:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%lcf_strat, diag_4d(:,:,:,diag_pt%lcf_strat),  &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%mfls_strat, diag_4d(:,:,:,diag_pt%mfls_strat), &
               Time, is, js, 1, rmask=mask3d)

!-----------------------------------------------------------------------
!    5) variables associated with ice particle number:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%nice, diag_4d(:,:,:,diag_pt%nice),   &
               Time, is, js, 1, rmask=mask3d)

!-----------------------------------------------------------------------
!    6) variables associated with precipitation and precipitation area:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%qrout, diag_4d(:,:,:,diag_pt%qrout),  &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qsout, diag_4d(:,:,:,diag_pt%qsout),   &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%rain3d, diag_4d_kp1(:,:,:,diag_pt%rain3d),  &
               Time, is, js, 1)
      used = send_data    &
              (diag_id%snow3d, diag_4d_kp1(:,:,:,diag_pt%snow3d),  &
               Time, is, js, 1)
      used = send_data   &
              (diag_id%rain_clr, diag_4d_kp1(:,:,:,diag_pt%rain_clr), &
               Time, is, js, 1, rmask=mask3)
      used = send_data   &
              (diag_id%a_rain_clr, diag_4d_kp1(:,:,:,diag_pt%a_rain_clr), &
               Time, is, js, 1, rmask=mask3)
      used = send_data   &
              (diag_id%rain_cld, diag_4d_kp1(:,:,:,diag_pt%rain_cld), &
               Time, is, js, 1, rmask=mask3)
      used = send_data   &
              (diag_id%a_rain_cld, diag_4d_kp1(:,:,:,diag_pt%a_rain_cld), &
               Time, is, js, 1,rmask=mask3)
      used = send_data   &
              (diag_id%a_precip_clr,   &
                                 diag_4d_kp1(:,:,:,diag_pt%a_precip_clr), &
               Time, is, js, 1, rmask=mask3)
      used = send_data   &
              (diag_id%a_precip_cld,    &
                                 diag_4d_kp1(:,:,:,diag_pt%a_precip_cld), &
               Time, is, js, 1,rmask=mask3)
      used = send_data   &
              (diag_id%snow_clr, diag_4d_kp1(:,:,:, diag_pt%snow_clr), &
               Time, is, js, 1, rmask=mask3)
      used = send_data   &
               (diag_id%a_snow_clr, diag_4d_kp1(:,:,:,diag_pt%a_snow_clr),&
                Time, is, js, 1, rmask=mask3)
      used = send_data   &
              (diag_id%snow_cld, diag_4d_kp1(:,:,:,diag_pt%snow_cld), &
               Time, is, js, 1, rmask=mask3)
      used = send_data   &
              (diag_id%a_snow_cld, diag_4d_kp1(:,:,:,diag_pt%a_snow_cld), &
               Time, is, js, 1, rmask=mask3)

!-----------------------------------------------------------------------
!    7) variables associated with cloud fraction:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%aall, diag_4d(:,:,:,diag_pt%aall),   &
               Time, is, js, 1, rmask=mask3d)
      used = send_data    &
              (diag_id%aliq, diag_4d(:,:,:,diag_pt%aliq),   &
               Time, is, js, 1, rmask=mask3d)
      used = send_data    &
              (diag_id%aice, diag_4d(:,:,:,diag_pt%aice),   &
               Time, is, js, 1, rmask=mask3d)
      mask4 = diag_4d(:,:,:,diag_pt%cfin) .ne. -1.e30
      used = send_data   &
              (diag_id%cfin, diag_4d(:,:,:,diag_pt%cfin), &
               Time, is, js, 1, mask=mask4)

!-----------------------------------------------------------------------
!    8) variables associated with cloud liquid time tendency:
!-----------------------------------------------------------------------
      used = send_data    &
              (diag_id%qldt_cond, diag_4d(:,:,:,diag_pt%qldt_cond), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data    &
              (diag_id%qldt_evap,  diag_4d(:,:,:,diag_pt%qldt_evap), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data    &
              (diag_id%qldt_eros, diag_4d(:,:,:,diag_pt%qldt_eros), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qldt_accr, diag_4d(:,:,:,diag_pt%qldt_accr), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qldt_auto, diag_4d(:,:,:,diag_pt%qldt_auto), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%liq_adj, diag_4d(:,:,:,diag_pt%liq_adj), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qldt_fill, diag_4d(:,:,:,diag_pt%qldt_fill), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data  &
              (diag_id%qldt_berg, diag_4d(:,:,:,diag_pt%qldt_berg), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qldt_freez, diag_4d(:,:,:,diag_pt%qldt_freez), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qldt_rime, diag_4d(:,:,:,diag_pt%qldt_rime), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qldt_destr, diag_4d(:,:,:,diag_pt%qldt_destr), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qldt_freez2, diag_4d(:,:,:,diag_pt%qldt_freez2), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qldt_sedi, diag_4d(:,:,:,diag_pt%qldt_sedi), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qldt_accrs, diag_4d(:,:,:,diag_pt%qldt_accrs), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qldt_bergs, diag_4d(:,:,:,diag_pt%qldt_bergs), &
               Time, is, js, 1, rmask=mask3d)

!-----------------------------------------------------------------------
!    9) variables associated with cloud droplet number time tendency:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%qndt_cond, diag_4d(:,:,:,diag_pt%qndt_cond), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_evap, diag_4d(:,:,:,diag_pt%qndt_evap), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_fill, diag_4d(:,:,:,diag_pt%qndt_fill), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_destr, diag_4d(:,:,:,diag_pt%qndt_destr), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_super, diag_4d(:,:,:,diag_pt%qndt_super), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_berg, diag_4d(:,:,:,diag_pt%qndt_berg), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_freez, diag_4d(:,:,:,diag_pt%qndt_freez), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_sacws, diag_4d(:,:,:,diag_pt%qndt_sacws), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data  &
              (diag_id%qndt_sacws_o, diag_4d(:,:,:,diag_pt%qndt_sacws_o), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_eros, diag_4d(:,:,:,diag_pt%qndt_eros), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_pra, diag_4d(:,:,:,diag_pt%qndt_pra), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_auto, diag_4d(:,:,:,diag_pt%qndt_auto), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_nucclim, diag_4d(:,:,:,diag_pt%qndt_nucclim), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_sedi, diag_4d(:,:,:,diag_pt%qndt_sedi), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_melt, diag_4d(:,:,:,diag_pt%qndt_melt), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_ihom, diag_4d(:,:,:,diag_pt%qndt_ihom), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_size_adj,diag_4d(:,:,:,diag_pt%qndt_size_adj),&
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qndt_fill2, diag_4d(:,:,:,diag_pt%qndt_fill2), &
               Time, is, js, 1, rmask=mask3d)

!-----------------------------------------------------------------------
!    10) variables associated with ice particle number time tendency:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%qnidt_fill, diag_4d(:,:,:,diag_pt%qnidt_fill), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data  &
              (diag_id%qnidt_nnuccd, diag_4d(:,:,:,diag_pt%qnidt_nnuccd), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qnidt_nsubi, diag_4d(:,:,:,diag_pt%qnidt_nsubi), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data  & 
              (diag_id%qnidt_nerosi, diag_4d(:,:,:,diag_pt%qnidt_nerosi), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qnidt_nprci, diag_4d(:,:,:,diag_pt%qnidt_nprci), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qnidt_nprai, diag_4d(:,:,:,diag_pt%qnidt_nprai), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qnidt_nucclim1,   &
                                  diag_4d(:,:,:,diag_pt%qnidt_nucclim1), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data  &
              (diag_id%qnidt_nucclim2,    &
                                  diag_4d(:,:,:,diag_pt%qnidt_nucclim2), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qnidt_sedi, diag_4d(:,:,:,diag_pt%qnidt_sedi), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qnidt_melt, diag_4d(:,:,:,diag_pt%qnidt_melt), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qnidt_size_adj, &
                                diag_4d(:,:,:,diag_pt%qnidt_size_adj), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qnidt_fill2, diag_4d(:,:,:,diag_pt%qnidt_fill2), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qnidt_super, diag_4d(:,:,:,diag_pt%qnidt_super), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qnidt_ihom, diag_4d(:,:,:,diag_pt%qnidt_ihom), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qnidt_destr, diag_4d(:,:,:,diag_pt%qnidt_destr), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data    &
              (diag_id%qnidt_cleanup,     &
                                  diag_4d(:,:,:,diag_pt%qnidt_cleanup), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data    &
              (diag_id%qnidt_cleanup2,     &
                                  diag_4d(:,:,:,diag_pt%qnidt_cleanup2), &
               Time, is, js, 1, rmask=mask3d)

!-----------------------------------------------------------------------
!    11) variables associated with relative humidity:
!-----------------------------------------------------------------------
      used = send_data    &
              (diag_id%rhcrit, diag_4d(:,:,:,diag_pt%rhcrit), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%rhcrit_min, diag_4d(:,:,:,diag_pt%rhcrit_min), &
               Time, is, js, 1, rmask=mask3d)
      mask4 = diag_4d(:,:,:,diag_pt%rhiin) .ne. -1.e30
      used = send_data   &
              (diag_id%rhiin, 100.*diag_4d(:,:,:,diag_pt%rhiin), &
               Time, is, js, 1, mask=mask4)
      mask4 = diag_4d(:,:,:,diag_pt%rhlin) .ne. -1.e30
      used = send_data   &
              (diag_id%rhlin, 100.*diag_4d(:,:,:,diag_pt%rhlin), &
               Time, is, js, 1, mask=mask4)

!-----------------------------------------------------------------------
!    12) variables associated with aerosol nucleation:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%imass7, diag_4d(:,:,:,diag_pt%imass7), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data    &
              (diag_id%ni_dust, diag_4d(:,:,:,diag_pt%ni_dust), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%ni_sulf, diag_4d(:,:,:,diag_pt%ni_sulf), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%ni_bc, diag_4d(:,:,:,diag_pt%ni_bc), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%ndust1, diag_4d(:,:,:,diag_pt%ndust1), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%ndust2, diag_4d(:,:,:,diag_pt%ndust2), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data    &
              (diag_id%ndust3, diag_4d(:,:,:,diag_pt%ndust3), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%ndust4, diag_4d(:,:,:,diag_pt%ndust4), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%ndust5, diag_4d(:,:,:,diag_pt%ndust5), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%sulfate, diag_4d(:,:,:,diag_pt%sulfate),  &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%seasalt_sub, diag_4d(:,:,:,diag_pt%seasalt_sub), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%seasalt_sup, diag_4d(:,:,:,diag_pt%seasalt_sup), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%om, diag_4d(:,:,:,diag_pt%om),    &
                Time, is, js, 1, rmask=mask3d)

!-----------------------------------------------------------------------
!    13) variables associated with water vapor tendency:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%rain_evap, diag_4d(:,:,:,diag_pt%rain_evap), &
               Time, is, js, 1, rmask=mask3d)

!-----------------------------------------------------------------------
!    14) variables associated with cloud ice tendency:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%qidt_dep, diag_4d(:,:,:,diag_pt%qidt_dep), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qidt_subl, diag_4d(:,:,:,diag_pt%qidt_subl), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qidt_eros, diag_4d(:,:,:,diag_pt%qidt_eros), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data  &
              (diag_id%qidt_fall, diag_4d(:,:,:,diag_pt%qidt_fall), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data  &
              (diag_id%qidt_melt, diag_4d(:,:,:,diag_pt%qidt_melt), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data  &
              (diag_id%qidt_melt2, diag_4d(:,:,:,diag_pt%qidt_melt2), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%ice_adj, diag_4d(:,:,:,diag_pt%ice_adj), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qidt_destr, diag_4d(:,:,:,diag_pt%qidt_destr), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qidt_qvdep, diag_4d(:,:,:,diag_pt%qidt_qvdep), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qidt_fill, diag_4d(:,:,:,diag_pt%qidt_fill), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qidt_auto, diag_4d(:,:,:,diag_pt%qidt_auto), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qidt_accr, diag_4d(:,:,:,diag_pt%qidt_accr), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qidt_accrs, diag_4d(:,:,:,diag_pt%qidt_accrs), &
               Time, is, js, 1, rmask=mask3d)

!-----------------------------------------------------------------------
!    15) variables associated with cloud area tendency:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%qadt_lsform, diag_4d(:,:,:,diag_pt%qadt_lsform), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qadt_lsdiss, diag_4d(:,:,:,diag_pt%qadt_lsdiss), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qadt_rhred, diag_4d(:,:,:,diag_pt%qadt_rhred), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qadt_eros, diag_4d(:,:,:,diag_pt%qadt_eros), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qadt_fill, diag_4d(:,:,:,diag_pt%qadt_fill), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qadt_super, diag_4d(:,:,:,diag_pt%qadt_super), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qadt_destr, diag_4d(:,:,:,diag_pt%qadt_destr), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qadt_limits, diag_4d(:,:,:,diag_pt%qadt_limits), &
               Time, is, js, 1, rmask=mask3d)
      used = send_data   &
              (diag_id%qadt_ahuco, diag_4d(:,:,:,diag_pt%qadt_ahuco), &
               Time, is, js, 1, rmask=mask3d)

!------------------------------------------------------------------------
!   16)  variables added by h1g with ncar M-G microphysics
!------------------------------------------------------------------------
      used = send_data ( diag_id%SA3d, diag_4d(:,:,:,diag_pt%SA3d),  &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%ST3d, diag_4d(:,:,:,diag_pt%ST3d),  &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%SQ3d, diag_4d(:,:,:,diag_pt%SQ3d),  &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%SL3d, diag_4d(:,:,:,diag_pt%SL3d),  &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%SI3d, diag_4d(:,:,:,diag_pt%SI3d),  &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%SN3d, diag_4d(:,:,:,diag_pt%SN3d),  &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%SNI3d, diag_4d(:,:,:,diag_pt%SNI3d),  &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qndt_contact_frz,    &
                          diag_4d(:,:,:,diag_pt%qndt_contact_frz), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qndt_cleanup,    &
                          diag_4d(:,:,:,diag_pt%qndt_cleanup), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qndt_cleanup2,   &
                          diag_4d(:,:,:,diag_pt%qndt_cleanup2), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qnidt_nsacwi,   &
                          diag_4d(:,:,:,diag_pt%qnidt_nsacwi), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_liquid_init,    &
                          diag_4d(:,:,:,diag_pt%qdt_liquid_init), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_ice_init,    &
                          diag_4d(:,:,:,diag_pt%qdt_ice_init), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_rain_evap,     &
                          diag_4d(:,:,:,diag_pt%qdt_rain_evap), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_snow_sublim,    &
                          diag_4d(:,:,:,diag_pt%qdt_snow_sublim), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_cond,    &
                          diag_4d(:,:,:,diag_pt%qdt_cond), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_deposition,    &
                          diag_4d(:,:,:,diag_pt%qdt_deposition), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_eros_l,   &
                          diag_4d(:,:,:,diag_pt%qdt_eros_l), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_eros_i,     &
                         diag_4d(:,:,:,diag_pt%qdt_eros_i), &
                         Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_qv_on_qi,     &
                          diag_4d(:,:,:,diag_pt%qdt_qv_on_qi), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_snow2vapor,    &
                          diag_4d(:,:,:,diag_pt%qdt_snow2vapor), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_sedi_ice2vapor,    &
                          diag_4d(:,:,:,diag_pt%qdt_sedi_ice2vapor), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_sedi_liquid2vapor,     &
                          diag_4d(:,:,:,diag_pt%qdt_sedi_liquid2vapor), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_super_sat_rm,    &
                          diag_4d(:,:,:,diag_pt%qdt_super_sat_rm), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_destr,    &
                          diag_4d(:,:,:,diag_pt%qdt_destr), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_cleanup_liquid,    &
                          diag_4d(:,:,:,diag_pt%qdt_cleanup_liquid), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%qdt_cleanup_ice,     &
                          diag_4d(:,:,:,diag_pt%qdt_cleanup_ice), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%srfrain_evap,    &
                          diag_4d(:,:,:,diag_pt%srfrain_evap ), &
                          Time, is, js, 1, rmask=mask3d )
      used = send_data ( diag_id%srfrain_accrs,     &
                          diag_4d(:,:,:,diag_pt%srfrain_accrs ), &
                          Time, is, js, 1, rmask=mask3d )
     used = send_data ( diag_id%srfrain_freez,     &
                         diag_4d(:,:,:,diag_pt%srfrain_freez ), &
                         Time, is, js, 1, rmask=mask3d )
     used = send_data ( diag_id%snow_mass_conv,   &
                        diag_3d(:,:,  diag_pt%snow_mass_conv), &
                        Time, is, js )
     used = send_data ( diag_id%neg_snow,     &
                         diag_3d(:,:,  diag_pt%neg_snow), &
                         Time, is, js )
     used = send_data ( diag_id%rain_mass_conv,    &
                         diag_3d(:,:,  diag_pt%rain_mass_conv), &
                         Time, is, js )
     used = send_data ( diag_id%neg_rain,    &
                         diag_3d(:,:,diag_pt%neg_rain), &
                         Time, is, js )
     used = send_data ( diag_id%qldt_HM_splinter,    &
                         diag_4d(:,:,:,diag_pt%qldt_HM_splinter), &
                         Time, is, js, 1, rmask=mask3d )

!-----------------------------------------------------------------------
!    17) variables associated with budget verification:
!-----------------------------------------------------------------------
     used = send_data ( diag_id%SA_imb, diag_4d(:,:,:,diag_pt%SA_imb), &
                               Time, is, js, 1, rmask=mask3d )
     used = send_data ( diag_id%ST_imb, diag_4d(:,:,:,diag_pt%ST_imb), &
                               Time, is, js, 1, rmask=mask3d )
     used = send_data ( diag_id%SQ_imb, diag_4d(:,:,:,diag_pt%SQ_imb), &
                               Time, is, js, 1, rmask=mask3d )
     used = send_data ( diag_id%SL_imb, diag_4d(:,:,:,diag_pt%SL_imb), &
                               Time, is, js, 1, rmask=mask3d )
     used = send_data ( diag_id%SI_imb, diag_4d(:,:,:,diag_pt%SI_imb), &
                               Time, is, js, 1, rmask=mask3d )
     used = send_data ( diag_id%SN_imb, diag_4d(:,:,:,diag_pt%SN_imb), &
                               Time, is, js, 1, rmask=mask3d )
     used = send_data ( diag_id%SNi_imb, diag_4d(:,:,:,diag_pt%SNi_imb), &
                               Time, is, js, 1, rmask=mask3d )
     used = send_data ( diag_id%rain_imb, diag_4d(:,:,:,diag_pt%rain_imb),&
                               Time, is, js, 1, rmask=mask3d )
     used = send_data ( diag_id%cld_liq_imb,   &
                         diag_4d(:,:,:,diag_pt%cld_liq_imb), &
                         Time, is, js, 1, rmask=mask3d )
     used = send_data ( diag_id%snow_imb,    &
                         diag_4d(:,:,:,diag_pt%snow_imb), &
                         Time, is, js, 1, rmask=mask3d )
     used = send_data ( diag_id%cld_ice_imb,     &
                         diag_4d(:,:,:,diag_pt%cld_ice_imb), &
                         Time, is, js, 1, rmask=mask3d )

!-----------------------------------------------------------------------
!
!                    COLUMN-INTEGRATED DIAGNOSTICS
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!    1) variables associated with droplet number and size:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%droplets_col, diag_3d(:,:,diag_pt%droplets_col), &
               Time, is, js)
      used = send_data   &
              (diag_id%gb_droplets_col,   &
                                  diag_3d(:,:,diag_pt%gb_droplets_col),   &
               Time, is, js)
      used = send_data   &
              (diag_id%droplets_col250,   &
                                 diag_3d(:,:,diag_pt%droplets_col250),   &
               Time, is, js)

!-----------------------------------------------------------------------
!    3) variables associated with cloud and precipitation processes:
!-----------------------------------------------------------------------
      used = send_data  &
              (diag_id%sedi_ice, diag_3d(:,:,diag_pt%sedi_ice), &
               Time, is, js)
      used = send_data   &
              (diag_id%snow_melt_col, diag_3d(:,:, diag_pt%snow_melt), &
               Time, is, js)
      used = send_data   &
              (diag_id%rain_freeze_col, diag_3d(:,:, diag_pt%rain_freeze), &
               Time, is, js)
      used = send_data   &
              (diag_id%q_snow_sublim_col,    &
               diag_3d(:,:, diag_pt%qdt_snow_sublim), &
               Time, is, js)

      used = send_data   &
              (diag_id%sedi_sfc, diag_3d(:,:, diag_pt%sedi_sfc), &
               Time, is, js)

!-----------------------------------------------------------------------
!    5) variables associated with ice particle number:
!-----------------------------------------------------------------------
      used = send_data    &
              (diag_id%nice_col, diag_3d(:,:,diag_pt%nice_col),  &
               Time, is, js)
      used = send_data   &
              (diag_id%gb_nice_col, diag_3d(:,:,diag_pt%gb_nice_col),  &
               Time, is, js)

!-----------------------------------------------------------------------
!    6) variables associated with precipitation and precipitation area:
!-----------------------------------------------------------------------
      used = send_data    &
              (diag_id%qrout_col, diag_3d(:,:,diag_pt%qrout),  &
               Time, is, js)
      used = send_data   &
              (diag_id%qsout_col, diag_3d(:,:,diag_pt%qsout),  &
               Time, is, js)

!-----------------------------------------------------------------------
!    8) variables associated with cloud liquid time tendency:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%ql_cond_col, diag_3d(:,:,diag_pt%qldt_cond), &
               Time, is, js)
      used = send_data   & 
              (diag_id%ql_evap_col, diag_3d(:,:,diag_pt%qldt_evap), &
               Time, is, js)
      used = send_data   &
              (diag_id%ql_eros_col, diag_3d(:,:, diag_pt%qldt_eros), &
               Time, is, js)
      used = send_data   &
              (diag_id%ql_accr_col, diag_3d(:,:, diag_pt%qldt_accr), &
               Time, is, js)
      used = send_data   &
              (diag_id%ql_auto_col, diag_3d(:,:, diag_pt%qldt_auto), &
               Time, is, js)
      used = send_data   &
              (diag_id%ql_berg_col, diag_3d(:,:, diag_pt%qldt_berg), &
               Time, is, js)
      used = send_data   &
              (diag_id%ql_freez_col, diag_3d(:,:, diag_pt%qldt_freez), &
               Time, is, js)
      used = send_data   &
              (diag_id%ql_destr_col, diag_3d(:,:, diag_pt%qldt_destr), &
               Time, is, js)
      used = send_data   &
              (diag_id%ql_rime_col, diag_3d(:,:, diag_pt%qldt_rime), &
               Time, is, js)
      used = send_data   &
              (diag_id%ql_fill_col, diag_3d(:,:, diag_pt%qldt_fill), &
               Time, is, js)
      used = send_data   &
              (diag_id%liq_adj_col, diag_3d(:,:, diag_pt%liq_adj), &
               Time, is, js)
      used = send_data   &
              (diag_id%ql_freez2_col, diag_3d(:,:,diag_pt%qldt_freez2), &
               Time, is, js)
      used = send_data  &
              (diag_id%ql_sedi_col, diag_3d(:,:, diag_pt%qldt_sedi), &
               Time, is, js)
      used = send_data   &
              (diag_id%ql_accrs_col, diag_3d(:,:, diag_pt%qldt_accrs), &
               Time, is, js)
      used = send_data   &
              (diag_id%ql_bergs_col, diag_3d(:,:, diag_pt%qldt_bergs), &
               Time, is, js)

!-----------------------------------------------------------------------
!    9) variables associated with cloud droplet number time tendency:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%qn_cond_col, diag_3d(:,:, diag_pt%qndt_cond), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_evap_col, diag_3d(:,:, diag_pt%qndt_evap), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_fill_col, diag_3d(:,:, diag_pt%qndt_fill), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_berg_col, diag_3d(:,:, diag_pt%qndt_berg), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_destr_col, diag_3d(:,:, diag_pt%qndt_destr), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_super_col, diag_3d(:,:, diag_pt%qndt_super), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_freez_col, diag_3d(:,:, diag_pt%qndt_freez), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_sacws_col, diag_3d(:,:, diag_pt%qndt_sacws), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_sacws_o_col, diag_3d(:,:,diag_pt%qndt_sacws_o), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_eros_col, diag_3d(:,:, diag_pt%qndt_eros), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_pra_col, diag_3d(:,:, diag_pt%qndt_pra), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_auto_col, diag_3d(:,:, diag_pt%qndt_auto), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_nucclim_col, diag_3d(:,:,diag_pt%qndt_nucclim), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_sedi_col, diag_3d(:,:, diag_pt%qndt_sedi), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_melt_col, diag_3d(:,:, diag_pt%qndt_melt), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_ihom_col, diag_3d(:,:, diag_pt%qndt_ihom), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_size_adj_col,   &
                                     diag_3d(:,:,diag_pt%qndt_size_adj), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_fill2_col, diag_3d(:,:, diag_pt%qndt_fill2), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_contact_frz_col,    &
               diag_3d(:,:, diag_pt%qndt_contact_frz), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_cleanup_col,   &
               diag_3d(:,:, diag_pt%qndt_cleanup), &
               Time, is, js)
      used = send_data   &
              (diag_id%qn_cleanup2_col,    &
               diag_3d(:,:, diag_pt%qndt_cleanup2), &
               Time, is, js)

!-----------------------------------------------------------------------
!    13) variables associated with water vapor tendency:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%rain_evap_col, diag_3d(:,:, diag_pt%rain_evap), &
               Time, is, js)

!-----------------------------------------------------------------------
!    14) variables associated with cloud ice time tendency:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%qi_fall_col, diag_3d(:,:, diag_pt%qidt_fall), &
               Time, is, js)
      used = send_data   &
              (diag_id%qi_fill_col, diag_3d(:,:, diag_pt%qidt_fill), &
               Time, is, js)
      used = send_data   &
              (diag_id%qi_dep_col, diag_3d(:,:, diag_pt%qidt_dep),     &
               Time, is, js)
      used = send_data   &
              (diag_id%qi_subl_col, diag_3d(:,:, diag_pt%qidt_subl), &
               Time, is, js)
      used = send_data   &
              (diag_id%qi_eros_col, diag_3d(:,:, diag_pt%qidt_eros), &
               Time, is, js)
      used = send_data   &
              (diag_id%qi_destr_col, diag_3d(:,:, diag_pt%qidt_destr), &
               Time, is, js)
      used = send_data   &
              (diag_id%qi_qvdep_col, diag_3d(:,:, diag_pt%qidt_qvdep), &
               Time, is, js)
      used = send_data   &
              (diag_id%qi_melt_col, diag_3d(:,:, diag_pt%qidt_melt), &
               Time, is, js)
      used = send_data   &
              (diag_id%qi_melt2_col, diag_3d(:,:, diag_pt%qidt_melt2), &
               Time, is, js)
      used = send_data   &
              (diag_id%ice_adj_col,  diag_3d(:,:, diag_pt%ice_adj), &
               Time, is, js)
      used = send_data   &
              (diag_id%qi_auto_col, diag_3d(:,:, diag_pt%qidt_auto),     &
               Time, is, js)
      used = send_data   &
              (diag_id%qi_accr_col, diag_3d(:,:, diag_pt%qidt_accr),     &
               Time, is, js)
      used = send_data   &
              (diag_id%qi_accrs_col, diag_3d(:,:, diag_pt%qidt_accrs), &
               Time, is, js)

!-----------------------------------------------------------------------
!    15) variables associated with cloud area time tendency:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%qa_lsform_col, diag_3d(:,:, diag_pt%qadt_lsform),&
               Time, is, js)
      used = send_data   &
              (diag_id%qa_lsdiss_col, diag_3d(:,:, diag_pt%qadt_lsdiss),&
               Time, is, js)
      used = send_data   &
              (diag_id%qa_rhred_col, diag_3d(:,:, diag_pt%qadt_rhred), &
               Time, is, js)
      used = send_data   &
              (diag_id%qa_eros_col, diag_3d(:,:, diag_pt%qadt_eros), &
               Time, is, js)
      used = send_data   &
              (diag_id%qa_fill_col, diag_3d(:,:, diag_pt%qadt_fill), &
               Time, is, js)
      used = send_data   &
              (diag_id%qa_super_col, diag_3d(:,:, diag_pt%qadt_super), &
               Time, is, js)
      used = send_data   &
              (diag_id%qa_destr_col, diag_3d (:,:, diag_pt%qadt_destr), &
               Time, is, js)
      used = send_data   &
              (diag_id%qa_limits_col, diag_3d (:,:, diag_pt%qadt_limits), &
               Time, is, js)
      used = send_data   &
              (diag_id%qa_ahuco_col, diag_3d (:,:, diag_pt%qadt_ahuco), &
               Time, is, js)

!------------------------------------------------------------------------
!   16)  variables added by h1g with ncar M-G microphysics
!------------------------------------------------------------------------
      used = send_data ( diag_id%SA2d, diag_3d(:,:,diag_pt%SA3d),  &
                          Time, is, js )
      used = send_data ( diag_id%ST2d, diag_3d(:,:,diag_pt%ST3d),  &
                          Time, is, js )
      used = send_data ( diag_id%SQ2d, diag_3d(:,:,diag_pt%SQ3d),  &
                          Time, is, js )
      used = send_data ( diag_id%SL2d, diag_3d(:,:,diag_pt%SL3d),  &
                          Time, is, js )
      used = send_data ( diag_id%SI2d, diag_3d(:,:,diag_pt%SI3d),  &
                          Time, is, js )
      used = send_data ( diag_id%SN2d, diag_3d(:,:,diag_pt%SN3d),  &
                          Time, is, js )
      used = send_data ( diag_id%SNI2d, diag_3d(:,:,diag_pt%SNI3d),  &
                          Time, is, js )
      used = send_data (diag_id%q_liquid_init_col,  &
                  diag_3d(:,:,diag_pt%qdt_liquid_init),  Time, is, js )
      used = send_data (diag_id%q_ice_init_col,  &
                  diag_3d(:,:,diag_pt%qdt_ice_init),  Time, is, js )
      used = send_data (diag_id%q_rain_evap_col,  &
                  diag_3d(:,:,diag_pt%qdt_rain_evap),  Time, is, js )
      used = send_data (diag_id%q_cond_col,  &
                  diag_3d(:,:,diag_pt%qdt_cond),  Time, is, js )
      used = send_data (diag_id%q_deposition_col,  &
                  diag_3d(:,:,diag_pt%qdt_deposition),  Time, is, js )
      used = send_data (diag_id%q_eros_l_col,  &
                  diag_3d(:,:,diag_pt%qdt_eros_l),  Time, is, js )
      used = send_data (diag_id%q_eros_i_col,  &
                  diag_3d(:,:,diag_pt%qdt_eros_i),  Time, is, js )
      used = send_data (diag_id%q_qv_on_qi_col,  &
                  diag_3d(:,:,diag_pt%qdt_qv_on_qi),  Time, is, js )
      used = send_data (diag_id%q_sedi_ice2vapor_col,  &
                  diag_3d(:,:,diag_pt%qdt_sedi_ice2vapor),  Time, is, js )
      used = send_data (diag_id%q_sedi_liquid2vapor_col,  &
                  diag_3d(:,:,diag_pt%qdt_sedi_liquid2vapor), Time, is, js)
      used = send_data (diag_id%q_super_sat_rm_col,  &
                  diag_3d(:,:,diag_pt%qdt_super_sat_rm),  Time, is, js )
      used = send_data (diag_id%q_destr_col,  &
                  diag_3d(:,:,diag_pt%qdt_destr),  Time, is, js )
      used = send_data (diag_id%q_cleanup_liquid_col,  &
                  diag_3d(:,:,diag_pt%qdt_cleanup_liquid),  Time, is, js )
      used = send_data (diag_id%q_cleanup_ice_col,  &
                  diag_3d(:,:,diag_pt%qdt_cleanup_ice),  Time, is, js )
      used = send_data (diag_id%q_snow2vapor_col,  &
                  diag_3d(:,:,diag_pt%qdt_snow2vapor),  Time, is, js )
      used = send_data (diag_id%srfrain_accrs_col,  &
                  diag_3d(:,:,diag_pt%srfrain_accrs), Time, is, js )
      used = send_data (diag_id%srfrain_freez_col,   &
                  diag_3d(:,:,diag_pt%srfrain_freez), Time, is, js )
      used = send_data (diag_id%srfrain_evap_col,   &
                  diag_3d(:,:,diag_pt%srfrain_evap), Time, is, js )
      used = send_data (diag_id%ql_HM_splinter_col,  &
                  diag_3d(:,:,diag_pt%qldt_HM_splinter), Time, is, js )

!-----------------------------------------------------------------------
!    17) variables associated with budget verification:
!-----------------------------------------------------------------------
     used = send_data ( diag_id%SA_imb_col, diag_3d(:,:,diag_pt%SA_imb), &
                               Time, is, js )
     used = send_data ( diag_id%ST_imb_col, diag_3d(:,:,diag_pt%ST_imb), &
                               Time, is, js )
     used = send_data ( diag_id%SQ_imb_col, diag_3d(:,:,diag_pt%SQ_imb), &
                               Time, is, js )
     used = send_data ( diag_id%SL_imb_col, diag_3d(:,:,diag_pt%SL_imb), &
                               Time, is, js )
     used = send_data ( diag_id%SI_imb_col, diag_3d(:,:,diag_pt%SI_imb), &
                               Time, is, js )
     used = send_data ( diag_id%SN_imb_col, diag_3d(:,:,diag_pt%SN_imb), &
                               Time, is, js )
     used = send_data ( diag_id%SNi_imb_col, diag_3d(:,:,diag_pt%SNi_imb),&
                               Time, is, js )
     used = send_data ( diag_id%rain_imb_col,   &
                        diag_3d(:,:,diag_pt%rain_imb),&
                        Time, is, js )
     used = send_data ( diag_id%cld_liq_imb_col,   &
                        diag_3d(:,:,diag_pt%cld_liq_imb),&
                        Time, is, js )
     used = send_data ( diag_id%snow_imb_col,   &
                        diag_3d(:,:,diag_pt%snow_imb),&
                        Time, is, js )
     used = send_data ( diag_id%cld_ice_imb_col,   &
                        diag_3d(:,:,diag_pt%cld_ice_imb),&
                        Time, is, js )

!-----------------------------------------------------------------------
!    17) variables associated with ice particle number time tendency:
!-----------------------------------------------------------------------
      used = send_data   &
              (diag_id%qni_fill_col, diag_3d(:,:,diag_pt%qnidt_fill), &
               Time, is, js)
      used = send_data  &
              (diag_id%qni_nnuccd_col, diag_3d(:,:,diag_pt%qnidt_nnuccd), &
               Time, is, js)
      used = send_data   &
              (diag_id%qni_nsubi_col, diag_3d(:,:,diag_pt%qnidt_nsubi), &
               Time, is, js)
      used = send_data  & 
              (diag_id%qni_nerosi_col, diag_3d(:,:,diag_pt%qnidt_nerosi), &
               Time, is, js)
      used = send_data   &
              (diag_id%qni_nprci_col, diag_3d(:,:,diag_pt%qnidt_nprci), &
               Time, is, js)
      used = send_data   &
              (diag_id%qni_nprai_col, diag_3d(:,:,diag_pt%qnidt_nprai), &
               Time, is, js)
      used = send_data   &
              (diag_id%qni_nucclim1_col,   &
                                  diag_3d(:,:,diag_pt%qnidt_nucclim1), &
               Time, is, js)
      used = send_data  &
              (diag_id%qni_nucclim2_col,    &
                                  diag_3d(:,:,diag_pt%qnidt_nucclim2), &
               Time, is, js)
      used = send_data   &
              (diag_id%qni_sedi_col, diag_3d(:,:,diag_pt%qnidt_sedi), &
               Time, is, js)
      used = send_data   &
              (diag_id%qni_melt_col, diag_3d(:,:,diag_pt%qnidt_melt), &
               Time, is, js)
      used = send_data   &
              (diag_id%qni_size_adj_col, &
                                diag_3d(:,:,diag_pt%qnidt_size_adj), &
               Time, is, js)
      used = send_data   &
              (diag_id%qni_fill2_col, diag_3d(:,:,diag_pt%qnidt_fill2), &
               Time, is, js)
      used = send_data   &
              (diag_id%qni_super_col, diag_3d(:,:,diag_pt%qnidt_super), &
               Time, is, js)
      used = send_data   &
              (diag_id%qni_ihom_col, diag_3d(:,:,diag_pt%qnidt_ihom), &
               Time, is, js)
      used = send_data   &
              (diag_id%qni_destr_col, diag_3d(:,:,diag_pt%qnidt_destr), &
               Time, is, js)
      used = send_data    &
              (diag_id%qni_cleanup_col,     &
                                  diag_3d(:,:,diag_pt%qnidt_cleanup), &
               Time, is, js)
      used = send_data    &
              (diag_id%qni_nsacwi_col,     &
                                  diag_3d(:,:,diag_pt%qnidt_nsacwi), &
               Time, is, js)
      used = send_data    &
              (diag_id%qni_cleanup2_col,     &
                                  diag_3d(:,:,diag_pt%qnidt_cleanup2), &
               Time, is, js)

!------------------------------------------------------------------------


end subroutine strat_netcdf


!########################################################################

subroutine strat_netcdf_end

     module_is_initialized = .false.

end subroutine strat_netcdf_end


!#######################################################################


! <SUBROUTINE NAME="diag_field_init">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!         Initializes netcdf diagnostics.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diag_field_init (axes,Time)
!
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="integer">
!         Integer array containing axes integers.
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!         Time 
!  </IN>
! </SUBROUTINE>
!
!------------------------------------------------------------------------

subroutine diag_field_init (axes, Time, diag_id, diag_pt, n_diag_4d,  &
                            n_diag_4d_kp1)

!-----------------------------------------------------------------------
integer,            intent(in)    :: axes(4)
type(time_type),    intent(in)    :: Time
type(diag_id_type), intent(inout) :: diag_id
type(diag_pt_type), intent(inout) :: diag_pt
integer,            intent(out)   :: n_diag_4d, n_diag_4d_kp1

!-----------------------------------------------------------------------
!---local variables--------------------------------------------------
      integer, dimension(3) :: half = (/1,2,4/)

!-------------------------------------------------------------------------
!    register the available netcdf fields if they have been requested.
!-------------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                        3-DIMENSIONAL DIAGNOSTICS
!
!-----------------------------------------------------------------------

!------------------------------------------------------------------------
!    1) variables associated with droplet number and size:
!------------------------------------------------------------------------
      diag_id%droplets = register_diag_field (mod_name,    &
             'droplets', axes(1:3), Time,   &
             'Droplet number concentration', '/cm3',  &
             missing_value=missing_value )
      diag_id%droplets_wtd = register_diag_field (mod_name,   &
             'droplets_wtd', axes(1:3), Time,  &
             'Droplet number conc*Cld liq', 'kg/(kg*m3)',  &
             mask_variant=.true., missing_value=missing_value)
      diag_id%rvolume = register_diag_field (mod_name,   &
             'rvolume', axes(1:3), Time,   &
             'Cloud liquid mean volume radius', 'microns',       &
             missing_value=missing_value )

!------------------------------------------------------------------------
!    2) variables associated with cloud liquid content:
!------------------------------------------------------------------------
      diag_id%ql_wt = register_diag_field (mod_name,   &
             'ql_wt', axes(1:3), Time,   &
             'Cld liq for weighting droplet conc', 'kg/kg',  &
             mask_variant=.true., missing_value=missing_value)

!------------------------------------------------------------------------
!    3) variables associated with cloud and precipitation processes:
!------------------------------------------------------------------------
      diag_id%lsf_strat = register_diag_field (mod_name, &
             'lsf_strat', axes(1:3), Time, &
             'Condensation/deposition frequency from LS', 'none',   &
             missing_value=missing_value)
      diag_id%dcond = register_diag_field (mod_name,   &
             'dcond', axes(1:3), Time,    &
             'condensation/evaporation', 'kg/kg/s', &
             missing_value=missing_value)
      diag_id%aauto = register_diag_field (mod_name,   &
             'aauto', axes(1:3), Time,   &
             'Cloud fraction where autoconversion is occurring',  &
             'dimensionless',   &
             missing_value=missing_value)
      diag_id%vfall = register_diag_field (mod_name,   &
             'vfall', axes(1:3), Time,   &
             'Ice crystal fall speed', 'meters/second',          &
             missing_value=missing_value)
      diag_id%delta_cf = register_diag_field (mod_name,  &
             'delta_cf', axes(1:3), Time,   &
             'increase in cloud area due to condensation ',    &
             'dimensionless',    &
             missing_value=missing_value)
      diag_id%cf_liq_init   = register_diag_field (mod_name,   &
             'cf_liq_init', axes(1:3), Time,  &
             'adjusted large-scale liquid cloud area before physics', &
             'dimensionless',  &
             missing_value=missing_value)
      diag_id%subgrid_w_variance = register_diag_field (mod_name,  &
             'subgrid_w_variance', axes(1:3), Time,    &
             'standard  dev of subgrid vert vel distribution', 'm/s',   &
             missing_value=missing_value)
      diag_id%potential_droplets = register_diag_field (mod_name,   &
             'potential_droplets', axes(1:3), Time,   &
             'number of droplets which may be activated', '/cm3',   &
             missing_value=missing_value)
      diag_id%potential_crystals = register_diag_field (mod_name,    &
             'potential_crystals', axes(1:3), Time,    &
             'number of crystals which may be activated', '/m3',   &
             missing_value=missing_value)
      diag_id%dust_berg_flag = register_diag_field (mod_name,    &
             'dust_berg_flag', axes(1:3), Time,    &
             'frequency of active dust_berg process', 'none',   &
             missing_value=missing_value)
      diag_id%cf_ice_init = register_diag_field (mod_name,    &
             'cf_ice_init', axes(1:3), Time,   &
             'adjusted large-scale ice cloud area before physics',  &
             'dimensionless',  &
             missing_value=missing_value)
      diag_id%snow_melt = register_diag_field (mod_name, &
             'snow_melt', axes(1:3), Time, &
             'Rain water tendency from snow melting', 'kg/kg/sec',  &
             missing_value=missing_value)
      diag_id%rain_freeze = register_diag_field (mod_name, &
             'rain_freeze', axes(1:3), Time, &
             'Rain loss from rain freezing', 'kg/kg/sec',  &
             missing_value=missing_value)

!------------------------------------------------------------------------
!    4) variables associated with model convection:
!------------------------------------------------------------------------
      diag_id%lcf_strat = register_diag_field (mod_name, &
             'lcf_strat', axes(1:3), Time, & 
             'Convection frequency from LS', 'none',   &
             missing_value=missing_value)
      diag_id%mfls_strat = register_diag_field (mod_name, &
             'mfls_strat', axes(1:3), Time, &
             'Convective mass flux from LS', 'Pascal/s',  &
             missing_value=missing_value)

!------------------------------------------------------------------------
!    5) variables associated with ice particle number:
!------------------------------------------------------------------------
      diag_id%nice = register_diag_field (mod_name,   &
             'nice', axes(1:3), Time,    &
             'ice number concentration', '/cm3',   &
             missing_value=missing_value)

!------------------------------------------------------------------------
!    6) variables associated with precipitation and precipitation area:
!------------------------------------------------------------------------
      diag_id%qrout = register_diag_field (mod_name,   &
             'qrout', axes(1:3), Time,   &
             'rain mixing ratio', 'kg/kg',    &
             missing_value=missing_value)
      diag_id%qsout = register_diag_field (mod_name,    &
             'qsout', axes(1:3), Time,    &
             'snow mixing ratio', 'kg/kg',  &
             missing_value=missing_value)
      diag_id%rain3d = register_diag_field (mod_name,   &
             'rain3d', axes(half), Time,    &
             'rain rate', 'kg/m2/s',   &
             missing_value=missing_value)
      diag_id%snow3d = register_diag_field (mod_name,    &
             'snow3d', axes(half), Time,  &
             'snow rate', 'kg/m2/s',   &
             missing_value=missing_value)
      diag_id%rain_clr = register_diag_field (mod_name, &
             'rain_clr', axes(half), Time, &
             'Clear sky rain rate averaged to grid box mean', 'kg/m2/s', &
             missing_value=missing_value)
      diag_id%rain_cld = register_diag_field (mod_name, &
             'rain_cld', axes(half), Time, &
             'cloudy sky rain rate averaged to grid box mean', 'kg/m2/s', &
             missing_value=missing_value)
      diag_id%a_rain_clr = register_diag_field (mod_name, &
             'a_rain_clr', axes(half), Time, &
             'Clear sky rain fractional coverage', 'fraction',   &
             missing_value=missing_value)
      diag_id%a_rain_cld = register_diag_field (mod_name, &
             'a_rain_cld', axes(half), Time, &
             'cloudy sky rain fractional coverage', 'fraction',   &
             missing_value=missing_value)
      diag_id%a_precip_clr = register_diag_field (mod_name, &
             'a_precip_clr', axes(half), Time, &
             'Clear sky precip fractional coverage', 'fraction',   &
             missing_value=missing_value)
      diag_id%a_precip_cld = register_diag_field (mod_name, &
            'a_precip_cld', axes(half), Time, &
            'cloudy sky precip fractional coverage', 'fraction',   &
            missing_value=missing_value)
      diag_id%snow_clr = register_diag_field (mod_name, &
             'snow_clr', axes(half), Time, &
             'Clear sky snow rate averaged to grid box mean', 'kg/m2/s',  &
             missing_value=missing_value)
      diag_id%snow_cld = register_diag_field (mod_name, &
             'snow_cld', axes(half), Time, &
             'cloudy sky snow rate averaged to grid box mean', 'kg/m2/s', &
             missing_value=missing_value)
      diag_id%a_snow_clr = register_diag_field (mod_name, &
             'a_snow_clr', axes(half), Time, &
             'Clear sky snow fractional coverage', 'fraction',   &
             missing_value=missing_value)
      diag_id%a_snow_cld = register_diag_field (mod_name, &
             'a_snow_cld', axes(half), Time, &
             'cloudy sky snow fractional coverage', 'fraction',   &
             missing_value=missing_value)

!------------------------------------------------------------------------
!    7) variables associated with cloud fraction:
!------------------------------------------------------------------------
      diag_id%aall = register_diag_field (mod_name,   &
             'aall', axes(1:3), Time,    &
             'Cloud fraction for all large-scale clouds at midtimestep',  &
             'dimensionless',   &
             missing_value=missing_value)
      diag_id%aliq = register_diag_field (mod_name,   &
             'aliq', axes(1:3), Time,    &
             'Cloud fraction for large-scale liquid clouds',   &
             'dimensionless',&
             missing_value=missing_value)
      diag_id%aice = register_diag_field (mod_name,   &
             'aice', axes(1:3), Time,   &
             'Cloud fraction for large-scale ice clouds',    &
             'dimensionless',   &
             missing_value=missing_value)
      diag_id%cfin = register_diag_field (mod_name,    &
             'cfin', axes(1:3), Time,    &
             'cloud fraction at T < 268 K', 'dimensionless',   &
             mask_variant = .true., missing_value=missing_value)

!-----------------------------------------------------------------------
!    8) variables associated with cloud liquid time tendency:
!-----------------------------------------------------------------------
      diag_id%qldt_cond = register_diag_field (mod_name, &
             'qldt_cond', axes(1:3), Time, &
             'Liq water specific humidity tendency from LS condensation', &
             'kg/kg/sec',    &
             missing_value=missing_value)
      diag_id%qldt_evap = register_diag_field (mod_name, &
             'qldt_evap', axes(1:3), Time, &
             'Liq water specific humidity tendency from LS evaporation', &
             'kg/kg/sec',   &
             missing_value=missing_value)
      diag_id%qldt_eros = register_diag_field (mod_name, &
             'qldt_eros', axes(1:3), Time, &
             'Liquid water specific humidity tendency from erosion',   &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%qldt_accr = register_diag_field (mod_name, &
             'qldt_accr', axes(1:3), Time, &
             'Liquid water specific humidity tendency from &
                                                    &accretion by rain', &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%qldt_auto = register_diag_field (mod_name, &
             'qldt_auto', axes(1:3), Time, &
             'Liq water specific humidity tendency from   &
                                                &autoconversion to rain',&
             'kg/kg/sec', missing_value=missing_value)
      diag_id%liq_adj = register_diag_field (mod_name, &
             'liq_adj', axes(1:3), Time, &
             'Liquid water specific humidity tendency from   &
                              &the removal of supersaturation',&
             'kg/kg/sec', missing_value=missing_value)
      diag_id%qldt_fill = register_diag_field (mod_name, &
             'qldt_fill', axes(1:3), Time, &
           'Liquid water specific humidity tendency -- &      
                                       &pre-microphysics filler', &
             'kg/kg/sec',  missing_value=missing_value)
      diag_id%qldt_berg = register_diag_field (mod_name, &
             'qldt_berg', axes(1:3), Time, &
             'Liq water specific humidity tendency from Bergeron process',&
             'kg/kg/sec',    &
             missing_value=missing_value)
      diag_id%qldt_freez = register_diag_field (mod_name, &
             'qldt_freez', axes(1:3), Time, &
             'Liq water specific humidity tendency from homogen freezing',&
             'kg/kg/sec',   &
             missing_value=missing_value)
      diag_id%qldt_rime = register_diag_field (mod_name, &
             'qldt_rime', axes(1:3), Time, & 
             'Liquid water specific humidity tendency from riming',    &
             'kg/kg/sec',   &
             missing_value=missing_value)
      diag_id%qldt_destr = register_diag_field (mod_name, &
             'qldt_destr', axes(1:3), Time, &
             'Liq water specific humidity tendency from cld destruction',&
             'kg/kg/sec',   &
             missing_value=missing_value)
      diag_id%qldt_freez2 = register_diag_field (mod_name, &
             'qldt_freez2', axes(1:3), Time, &
             'Liq water spec hum tend from heterogen drop freezing ', &
             'kg/kg/sec',   &
             missing_value=missing_value)
      diag_id%qldt_sedi = register_diag_field (mod_name, &
             'qldt_sedi', axes(1:3), Time, &
             'Liq water spec hum tendency from droplet sedimentation', &
             'kg/kg/sec',   &
             missing_value=missing_value)
      diag_id%sedi_sfc = register_diag_field (mod_name, &
             'sedi_sfc', axes(1:2), Time, &
             'Liquid sedimentation flux at surface', &
             'm/sec',   &
             missing_value=missing_value)
      diag_id%qldt_accrs = register_diag_field (mod_name, &
             'qldt_accrs', axes(1:3), Time, &
             'Liq water spec hum tend from collection of    &
                                             &cloud liquid by snow ', &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%qldt_bergs = register_diag_field (mod_name, &
             'qldt_bergs', axes(1:3), Time, &
             'Liq water spec hum tendency from bergeron process for snow',&
             'kg/kg/sec',   &
             missing_value=missing_value)

!-----------------------------------------------------------------------
!    9) variables associated with cloud droplet number time tendency:
!-----------------------------------------------------------------------
      diag_id%qndt_cond = register_diag_field (mod_name, &
             'qndt_cond', axes(1:3), Time, &
             'Cloud droplet number tendency from LS condensation', &
             '#/kg/sec',   &
             missing_value=missing_value)
      diag_id%qndt_evap = register_diag_field (mod_name, &
             'qndt_evap', axes(1:3), Time, &
             'Cloud droplet number tendency from LS evaporation',   &
             '#/kg/sec',  &
             missing_value=missing_value)
      diag_id%qndt_fill = register_diag_field (mod_name, &
             'qndt_fill', axes(1:3), Time, &
             'Cloud droplet number tendency from filler', '#/kg/sec',   &
             missing_value=missing_value)
      diag_id%qndt_destr = register_diag_field (mod_name, &
             'qndt_destr', axes(1:3), Time, &
             'Cloud droplet number tendency from cloud destruction',   &
             '#/kg/sec', &
             missing_value=missing_value)
      diag_id%qndt_super = register_diag_field (mod_name, &
             'qndt_super', axes(1:3), Time, &
             'Cloud droplet number tendency from supersaturation removal',&
             '#/kg/sec',    &
             missing_value=missing_value)
      diag_id%qndt_berg = register_diag_field (mod_name, &
             'qndt_berg', axes(1:3), Time, &
             'Cloud droplet number tendency from Bergeron', '#/kg/sec',   &
             missing_value=missing_value)
      diag_id%qndt_freez = register_diag_field (mod_name, &
             'qndt_freez', axes(1:3), Time, &
             'Cloud droplet number tendency from heterogeneous freezing', &
             '#/kg/sec',  &
             missing_value=missing_value)
      diag_id%qndt_sacws = register_diag_field (mod_name, &
             'qndt_sacws', axes(1:3), Time, &
             'Cloud droplet number tendency from collection by snow',  &
             '#/kg/sec',&
             missing_value=missing_value)
      diag_id%qndt_sacws_o = register_diag_field (mod_name, &
             'qndt_sacws_o', axes(1:3), Time, &
             'Cld drop number tend from collection by snw   &
                                              &- one ice cat scheme', &
             '#/kg/sec', missing_value=missing_value)
      diag_id%qndt_eros = register_diag_field (mod_name, &
             'qndt_eros', axes(1:3), Time, &
             'Cloud droplet number tendency from erosion', '#/kg/sec',   &
             missing_value=missing_value)
      diag_id%qndt_pra = register_diag_field (mod_name, &
             'qndt_pra', axes(1:3), Time, &
             'Cloud droplet number tendency from collection by rain',   &
             '#/kg/sec',&
             missing_value=missing_value)
      diag_id%qndt_auto = register_diag_field (mod_name, &
             'qndt_auto', axes(1:3), Time, &
             'Cloud droplet number tendency from autoconversion',   &
             '#/kg/sec',   &
              missing_value=missing_value)
      diag_id%qndt_nucclim = register_diag_field (mod_name, &
             'qndt_nucclim', axes(1:3), Time, &
             'Cloud droplet number tendency from nucleation limiter', &
             '#/kg/sec', missing_value=missing_value)     
      diag_id%qndt_sedi = register_diag_field (mod_name, &
             'qndt_sedi', axes(1:3), Time, &
             'Cloud droplet number tendency from sedimentation',   &
             '#/kg/sec', missing_value=missing_value)  
      diag_id%qndt_melt = register_diag_field (mod_name, &
             'qndt_melt', axes(1:3), Time, &
             'Cloud droplet number tendency from melting', '#/kg/sec',   &
             missing_value=missing_value) 
      diag_id%qndt_ihom = register_diag_field (mod_name, &
             'qndt_ihom', axes(1:3), Time, &
             'Cloud drop number tend from homogeneous freezing',   &
             '#/kg/sec', missing_value=missing_value)  
      diag_id%qndt_size_adj = register_diag_field (mod_name, &
             'qndt_size_adj', axes(1:3), Time, &
             'Cloud droplet number tendency from size adjustment',  &
             '#/kg/sec', missing_value=missing_value)  
      diag_id%qndt_fill2 = register_diag_field (mod_name, &
             'qndt_fill2', axes(1:3), Time, &
             'Cloud droplet number tendency from second filler',  &
             '#/kg/sec', missing_value=missing_value)  

  !-----------------------------------------------------------------------
  !    10) variables associated with ice particle number time tendency:
  !-----------------------------------------------------------------------
      diag_id%qnidt_fill = register_diag_field (mod_name, &
             'qnidt_fill', axes(1:3), Time, &
             'Ice particle number tendency from filler', '#/kg/sec',   &
             missing_value=missing_value)  
      diag_id%qnidt_nnuccd = register_diag_field (mod_name, &
             'qnidt_nnuccd', axes(1:3), Time, &
             'Ice particle number tendency from nucleation', '#/kg/sec',  &
             missing_value=missing_value)  
      diag_id%qnidt_nsubi = register_diag_field (mod_name, &
             'qnidt_nsubi', axes(1:3), Time, &
             'Ice particle number tendency from sublimation', '#/kg/sec', &
             missing_value=missing_value)  
      diag_id%qnidt_nerosi = register_diag_field (mod_name, &
             'qnidt_nerosi', axes(1:3), Time, &
             'Ice particle number tendency from erosion', '#/kg/sec',   &
             missing_value=missing_value) 
      diag_id%qnidt_nprci = register_diag_field (mod_name, &
             'qnidt_nprci', axes(1:3), Time, &
             'Ice particle number tendency from autoconversion', &
             '#/kg/sec', missing_value=missing_value) 
      diag_id%qnidt_nprai = register_diag_field (mod_name, &
             'qnidt_nprai', axes(1:3), Time, &
             'Ice particle number tendency from accretion by snow',   &
             '#/kg/sec', missing_value=missing_value) 
      diag_id%qnidt_nucclim1 = register_diag_field (mod_name, &
             'qnidt_nucclim1', axes(1:3), Time, &
             'Ice particle number tendency from first nucleation limiter',&
             '#/kg/sec', missing_value=missing_value) 
      diag_id%qnidt_nucclim2 = register_diag_field (mod_name, &
             'qnidt_nucclim2', axes(1:3), Time, &
             'Ice particle number tendency from second nucleation &
             &limiter','#/kg/sec', missing_value=missing_value) 
      diag_id%qnidt_sedi = register_diag_field (mod_name, &
             'qnidt_sedi', axes(1:3), Time, &
             'Ice particle number tendency from sedimentation', &
             '#/kg/sec', missing_value=missing_value) 
      diag_id%qnidt_melt = register_diag_field (mod_name, &
             'qnidt_melt', axes(1:3), Time, &
             'Ice particle number tendency from melting', '#/kg/sec',  &
             missing_value=missing_value) 
      diag_id%qnidt_size_adj = register_diag_field (mod_name, &
             'qnidt_size_adj', axes(1:3), Time, &
             'Ice particle number tendency from size adjustment', &
             '#/kg/sec', missing_value=missing_value) 
      diag_id%qnidt_fill2 = register_diag_field (mod_name, &
             'qnidt_fill2', axes(1:3), Time, &
             'Ice particle number tendency from second filler', &
             '#/kg/sec', missing_value=missing_value)  
      diag_id%qnidt_super = register_diag_field (mod_name, &
             'qnidt_super', axes(1:3), Time, &
             'Ice particle number tendency from removal of   &
             &supersaturation', '#/kg/sec', missing_value=missing_value)  
      diag_id%qnidt_ihom = register_diag_field (mod_name, &
             'qnidt_ihom', axes(1:3), Time, &
             'Ice particle number tendency from homogeneous freezing', &
             '#/kg/sec', missing_value=missing_value)  
      diag_id%qnidt_destr = register_diag_field (mod_name, &
             'qnidt_destr', axes(1:3), Time, &
             'Ice particle number tendency from cloud destruction', &
             '#/kg/sec', missing_value=missing_value)  
      diag_id%qnidt_cleanup = register_diag_field (mod_name, &
             'qnidt_cleanup', axes(1:3), Time, &
             'Ice particle number tendency from cleanup', '#/kg/sec',   &
             missing_value=missing_value)  
      diag_id%qnidt_cleanup2 = register_diag_field (mod_name, &
             'qnidt_cleanup2', axes(1:3), Time, &
             'Ice particle number tendency from cleanup2', '#/kg/sec',  &
             missing_value=missing_value)  

  !-----------------------------------------------------------------------
  !    11) variables associated with relative humidity:
  !-----------------------------------------------------------------------
      diag_id%rhcrit = register_diag_field (mod_name,    &
             'rhcrit', axes(1:3), Time,   &
             'spectral avgd critical rh used for ice nuclei activation', &
             'percent', missing_value=missing_value)
      diag_id%rhcrit_min = register_diag_field (mod_name,   &
             'rhcrit_min', axes(1:3), Time,   &
             'minimum critical rh for ice nuclei activation in any  &
             &spectral band', 'percent',  missing_value=missing_value)
      diag_id%rhiin = register_diag_field (mod_name,    &
             'rhiin', axes(1:3), Time,   &
             'rh wrt ice calculated for ice nulei activation', 'percent', &
             mask_variant = .true., missing_value=missing_value)
      diag_id%rhlin = register_diag_field (mod_name,   &
             'rhlin', axes(1:3), Time,    &
             'rh wrt liquid calculated for ice nulei activation',   &
             'percent', mask_variant = .true., &
             missing_value=missing_value)

  !-----------------------------------------------------------------------
  !    12) variables associated with aerosol nucleation:
  !-----------------------------------------------------------------------
      diag_id%imass7 = register_diag_field (mod_name,   &
             'imass7', axes(1:3), Time, &
             'sub-micron dust concentration', 'ug dust /m3', &
             missing_value=missing_value)
      diag_id%ni_dust = register_diag_field (mod_name,   &
             'ni_dust', axes(1:3), Time,   &
             'number of activated dust ice nuclei', '# / m3',    &
             missing_value=missing_value)
      diag_id%ni_sulf = register_diag_field (mod_name,    &
             'ni_sulf', axes(1:3), Time,   &
             'number  of activated sulfate/seasalt ice nuclei', '# / m3', &
             missing_value=missing_value)
      diag_id%ni_bc = register_diag_field (mod_name,    &
             'ni_bc', axes(1:3), Time,   &
             'number of activated bc ice nuclei', '# / m3',   &
             missing_value=missing_value)
      diag_id%ndust1 = register_diag_field (mod_name,    &
             'ndust1', axes(1:3), Time,   &
             'number of particles in dust bin1', 'none',   &
             missing_value=missing_value)
      diag_id%ndust2 = register_diag_field (mod_name,    &
             'ndust2', axes(1:3), Time,   &
             'number of particles in dust bin2', 'none',   &
             missing_value=missing_value)
      diag_id%ndust3 = register_diag_field (mod_name,    &
             'ndust3', axes(1:3), Time,    &
             'number of particles in dust bin3', 'none',   &
             missing_value=missing_value)
      diag_id%ndust4 = register_diag_field (mod_name,    &
             'ndust4', axes(1:3), Time,  &
             'number of particles in dust bin4', 'none',   &
             missing_value=missing_value)
      diag_id%ndust5 = register_diag_field (mod_name,   &
             'ndust5', axes(1:3), Time,   &
             'number of particles in dust bin5', 'none',   &
             missing_value=missing_value)
      diag_id%sulfate = register_diag_field (mod_name,   &
             'sulfate', axes(1:3), Time,   &
             'sulfate mass concentration', 'ug so4/m3',   &
             missing_value=missing_value)
      diag_id%seasalt_sub = register_diag_field (mod_name,   &
             'seasalt_sub', axes(1:3), Time,   &
             'sub-micron sea salt mass concentration', 'ug/m3',    &
             missing_value=missing_value)
      diag_id%seasalt_sup = register_diag_field (mod_name,   &
             'seasalt_sup', axes(1:3), Time,   &
             'super-micron sea salt mass concentration', 'ug/m3', &
             missing_value=missing_value)
      diag_id%om = register_diag_field (mod_name,    &
             'OM', axes(1:3), Time,  &
             'OM mass concentration', 'ug/m3',   &
             missing_value=missing_value)

  !-----------------------------------------------------------------------
  !    13) variables associated with water vapor tendency:
  !-----------------------------------------------------------------------
      diag_id%rain_evap = register_diag_field (mod_name, &
             'rain_evap', axes(1:3), Time, &
             'Water vapor tendency from rain evaporation', 'kg/kg/sec',  &
             missing_value=missing_value)

  !-----------------------------------------------------------------------
  !    14) variables associated with cloud ice time tendency:
  !-----------------------------------------------------------------------
      diag_id%qidt_dep = register_diag_field (mod_name, &
             'qidt_dep', axes(1:3), Time, &
             'Ice water specific humidity tendency from LS deposition', &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%qidt_subl = register_diag_field (mod_name, &
             'qidt_subl', axes(1:3), Time, &
             'Ice water specific humidity tendency from LS sublimation', &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%qidt_eros = register_diag_field (mod_name, &
             'qidt_eros', axes(1:3), Time, &
             'Ice water specific humidity tendency from erosion',      &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%qidt_fall = register_diag_field (mod_name, &
             'qidt_fall', axes(1:3), Time, &
             'Ice water specific humidity tendency from ice settling', &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%qidt_melt = register_diag_field (mod_name, &
             'qidt_melt', axes(1:3), Time, &
             'Ice water specific humidity tendency from melting to rain',&
             'kg/kg/sec', missing_value=missing_value)
      diag_id%qidt_melt2 = register_diag_field (mod_name, &
             'qidt_melt2', axes(1:3), Time, &
             'Ice water specific humidity tendency from melting to  &
             &cloud droplets', 'kg/kg/sec', missing_value=missing_value)
      diag_id%ice_adj = register_diag_field (mod_name, &
             'ice_adj', axes(1:3), Time, &
             'Ice water specific humidity tendency from the removal &
             &of supersaturation', 'kg/kg/sec',    &
             missing_value=missing_value)
      diag_id%qidt_destr = register_diag_field (mod_name, &
             'qidt_destr', axes(1:3), Time, &
             'Ice water spec hum tendency from cloud destruction',  &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%qidt_qvdep = register_diag_field (mod_name, &
             'qidt_qvdep', axes(1:3), Time, &
             'Ice water specific humidity tendency from vapor deposition',&
             'kg/kg/sec', missing_value=missing_value)
      diag_id%qidt_fill = register_diag_field (mod_name, &
             'qidt_fill', axes(1:3), Time, &
             'Ice water specific humidity tendency -- pre-microphysics &
             & filler', 'kg/kg/sec', missing_value=missing_value)
      diag_id%qidt_auto = register_diag_field (mod_name, &
             'qidt_auto', axes(1:3), Time, &
             'Ice water specific humidity tendency from autoconversion &
             &to snow', 'kg/kg/sec', missing_value=missing_value)
      diag_id%qidt_accr = register_diag_field (mod_name, &
             'qidt_accr', axes(1:3), Time, &
             'Ice water spec hum tendency from accretion by snow', &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%qidt_accrs = register_diag_field (mod_name, &
             'qidt_accrs', axes(1:3), Time, &
             'Ice wat spec hum tend from selfcollection (1 class scheme)',&
             'kg/kg/sec', missing_value=missing_value)

  !-----------------------------------------------------------------------
  !    15) variables associated with cloud area time tendency:
  !-----------------------------------------------------------------------
      diag_id%qadt_lsform = register_diag_field (mod_name, &
             'qadt_lsform', axes(1:3), Time, & 
             'cloud fraction tendency from LS condensation', '1/sec',  &
             missing_value=missing_value)
      diag_id%qadt_lsdiss = register_diag_field (mod_name, &
             'qadt_lsdiss', axes(1:3), Time, &
             'cloud fraction tendency from LS dissipation', '1/sec',  &
             missing_value=missing_value)
      diag_id%qadt_rhred = register_diag_field (mod_name, &
             'qadt_rhred', axes(1:3), Time, &
             'cloud fraction tendency from RH limiter', '1/sec',   &
             missing_value=missing_value)
      diag_id%qadt_eros = register_diag_field (mod_name, &
             'qadt_eros', axes(1:3), Time, &
             'cloud fraction tendency from erosion', '1/sec',   &
             missing_value=missing_value)
      diag_id%qadt_fill = register_diag_field (mod_name, &
             'qadt_fill', axes(1:3), Time, &
             'cloud fraction tendency from filler', '1/sec',   &
             missing_value=missing_value)
      diag_id%qadt_super = register_diag_field (mod_name, &
             'qadt_super', axes(1:3), Time, &
             'cloud fraction tendency from supersaturation formation', &
             '1/sec', missing_value=missing_value)
      diag_id%qadt_destr = register_diag_field (mod_name, &
             'qadt_destr', axes(1:3), Time, &
             'cloud fraction tendency from cloud destruction', '1/sec',  &
             missing_value=missing_value)
      diag_id%qadt_limits = register_diag_field (mod_name, &
             'qadt_limits', axes(1:3), Time, &
             'cloud fraction tendency from imposing limits', '1/sec',  &
             missing_value=missing_value)
      diag_id%qadt_ahuco = register_diag_field (mod_name, &
             'qadt_ahuco', axes(1:3), Time, &
             'cloud fraction tendency from convective area restriction', &
             '1/sec', missing_value=missing_value)

  !------------------------------------------------------------------------
  !   16)  variables added by h1g with ncar M-G microphysics
  !------------------------------------------------------------------------
      diag_id%SA3d =    register_diag_field ( mod_name, &
             'SA3d',  axes(1:3), Time, &
             'Total cloud area tendency', '1/sec',  &
             missing_value=missing_value)
      diag_id%ST3d =    register_diag_field ( mod_name, &
             'ST3d',  axes(1:3), Time, &
             'Total temperature tendency', &
             'K/sec', missing_value=missing_value)
      diag_id%SQ3d =    register_diag_field ( mod_name, &
             'SQ3d',  axes(1:3), Time, &
             'Total Water vapor specific humidity tendency ', &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%qdt_liquid_init =    register_diag_field ( mod_name, &
             'qdt_liquid_init',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from initial  &
             &adjustment of liquid water', &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%qdt_ice_init =    register_diag_field ( mod_name, &
             'qdt_ice_init',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from initial &
             &adjustment of ice', &
            'kg/kg/sec', missing_value=missing_value)
      diag_id%qdt_rain_evap =    register_diag_field ( mod_name, &
             'qdt_rain_evap',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from rain  &
             &evaporation', 'kg/kg/sec', missing_value=missing_value)
      diag_id%qdt_snow_sublim  =    register_diag_field ( mod_name, &
             'qdt_snow_sublim',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from snow   &
             &sublimation', 'kg/kg/sec', missing_value=missing_value)
      diag_id%qdt_cond  =    register_diag_field ( mod_name, &
             'qdt_cond',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from LS   &
             &condensation', 'kg/kg/sec', missing_value=missing_value)
      diag_id%qdt_deposition  =    register_diag_field ( mod_name, &
             'qdt_deposition',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from LS deposition', &
             'kg/kg/sec', missing_value=missing_value               )
      diag_id%qdt_eros_l  =    register_diag_field ( mod_name, &
             'qdt_eros_l',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from liquid water &
             & erosion', 'kg/kg/sec', missing_value=missing_value)
      diag_id%qdt_eros_i  =    register_diag_field ( mod_name, &
             'qdt_eros_i',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from ice water &
             &erosion', 'kg/kg/sec', missing_value=missing_value)
      diag_id%qdt_qv_on_qi  =    register_diag_field ( mod_name, &
             'qdt_qv_on_qi',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from direct &
             &deposition on ice', 'kg/kg/sec', missing_value=missing_value)
      diag_id%qdt_snow2vapor  =    register_diag_field ( mod_name, &
             'qdt_snow2vapor',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from condensation  &
             &resulting from cooling of the air by melting snow',      &
             'kg/kg/sec', missing_value=missing_value               )
      diag_id%qdt_sedi_ice2vapor =    register_diag_field ( mod_name, &
             'qdt_sedi_ice2vapor',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from the &
             &sublimation of falling cloud ice', 'kg/kg/sec',  &
             missing_value=missing_value)
      diag_id%qdt_sedi_liquid2vapor =    register_diag_field ( mod_name, &
             'qdt_sedi_liquid2vapor',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from the evaporation &
             &of falling cloud liquid', 'kg/kg/sec', &
             missing_value=missing_value)
      diag_id%qdt_super_sat_rm   =    register_diag_field ( mod_name, &
             'qdt_super_sat_rm',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from removing super &
             &saturation', 'kg/kg/sec', missing_value=missing_value)
      diag_id%qdt_destr  =    register_diag_field ( mod_name, &
             'qdt_destr',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from Cloud &
             &Destruction', 'kg/kg/sec', missing_value=missing_value)
      diag_id%qdt_cleanup_liquid  =    register_diag_field ( mod_name, &
             'qdt_cleanup_liquid',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from cleaning up &
             &liquid cloud', 'kg/kg/sec', missing_value=missing_value)
      diag_id%qdt_cleanup_ice  =    register_diag_field ( mod_name, &
             'qdt_cleanup_ice',  axes(1:3), Time, &
             'Water vapor specific humidity tendency from cleaning up &
             &ice cloud', 'kg/kg/sec', missing_value=missing_value)
      diag_id%SL3d = register_diag_field ( mod_name, &
             'SL3d',  axes(1:3), Time, &
             'Total Liquid Water specific humidity tendency ', &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%qldt_HM_splinter = register_diag_field ( mod_name, &
             'qldt_HM_splinter', axes(1:3), Time, &
             'Liquid water specific humidity tendency from H-M &
             &splintering ', 'kg/kg/sec', missing_value=missing_value)
      diag_id%SN3d =    register_diag_field ( mod_name, &
             'SN3d',  axes(1:3), Time, &
             'Total cloud drop number tendency', &
             '#/kg/sec', missing_value=missing_value)
      diag_id%qndt_contact_frz = register_diag_field ( mod_name, &
             'qndt_contact_frz', axes(1:3), Time, &
             'Cloud droplet number tendency from contact freezing of&
             & liquid water', '#/kg/sec', missing_value=missing_value)
      diag_id%qndt_cleanup = register_diag_field ( mod_name, &
             'qndt_cleanup', axes(1:3), Time, &
             'Cloud droplet number tendency from cloud cleanup', &
             '#/kg/sec', missing_value=missing_value)
      diag_id%qndt_cleanup2 = register_diag_field ( mod_name, &
             'qndt_cleanup2', axes(1:3), Time, &
             'Cloud droplet number tendency from cloud cleanup2', &
             '#/kg/sec', missing_value=missing_value)
      diag_id%SNI3d = register_diag_field ( mod_name, &
             'SNI3d',  axes(1:3), Time, &
             'Total ice crystal number tendency', &
             '#/kg/sec', missing_value=missing_value)
      diag_id%qnidt_nsacwi = register_diag_field ( mod_name, &
             'qnidt_nsacwi', axes(1:3), Time, &
             'Cloud ice tendency from HM ice multiplication', &
             '#/kg/sec', missing_value=missing_value)
      diag_id%srfrain_evap =  register_diag_field ( mod_name, &
             'srfrain_evap', axes(1:3), Time, &
             'rain water sink from evaporation', &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%srfrain_accrs=  register_diag_field ( mod_name, &
             'srfrain_accrs', axes(1:3), Time, &
             'rain water sink from collection by snow', &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%srfrain_freez=  register_diag_field ( mod_name, &
             'srfrain_freez', axes(1:3), Time, &
             'rain water sink from freezing', &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%SI3d =    register_diag_field ( mod_name, &
             'SI3d',  axes(1:3), Time, &
             'Total Ice water specific humidity tendency ', &
             'kg/kg/sec', missing_value=missing_value)

  !------------------------------------------------------------------------
  !   17)  variables associated with budget analysis
  !------------------------------------------------------------------------
      diag_id%SA_imb =    register_diag_field ( mod_name, &
             'SA_imb',  axes(1:3), Time, &
             'difference between qa tendency and sum of individ terms', &
             '1/sec', missing_value=missing_value)
      diag_id%ST_imb =    register_diag_field ( mod_name, &
             'ST_imb',  axes(1:3), Time, &
             'difference between temp tendency and sum of individ terms', &
             'K/sec', missing_value=missing_value)
      diag_id%SQ_imb =    register_diag_field ( mod_name, &
             'SQ_imb',  axes(1:3), Time, &
             'difference between qv tendency and sum of individ terms', &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%SL_imb =    register_diag_field ( mod_name, &
             'SL_imb',  axes(1:3), Time, &
             'difference between ql tendency and sum of individ terms', &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%SI_imb =    register_diag_field ( mod_name, &
             'SI_imb',  axes(1:3), Time, &
             'difference between qi tendency and sum of individ terms', &
             'kg/kg/sec', missing_value=missing_value)
      diag_id%SN_imb =    register_diag_field ( mod_name, &
             'SN_imb',  axes(1:3), Time, &
             'difference between qn tendency and sum of individ terms', &
             '#/kg/sec', missing_value=missing_value)
      diag_id%SNi_imb =    register_diag_field ( mod_name, &
             'SNi_imb',  axes(1:3), Time, &
             'difference between qni tendency and sum of individ terms', &
             '#/kg/sec', missing_value=missing_value)
      diag_id%rain_imb =    register_diag_field ( mod_name, &
             'rain_imb',  axes(1:3), Time, &
             'difference between rain rate at sfc and sum of &
             &individ terms', 'kg/kg/sec', missing_value=missing_value)
      diag_id%cld_liq_imb =    register_diag_field ( mod_name, &
             'cld_liq_imb',  axes(1:3), Time, &
             'difference between ql fallout rate at sfc and sum of &
             &individ terms', 'kg/kg/sec', missing_value=missing_value)
      diag_id%snow_imb =    register_diag_field ( mod_name, &
             'snow_imb',  axes(1:3), Time, &
             'difference between snow rate at sfc and sum of &
             &individ terms', 'kg/kg/sec', missing_value=missing_value)
      diag_id%cld_ice_imb =    register_diag_field ( mod_name, &
             'cld_ice_imb',  axes(1:3), Time, &
             'difference between qi fallout rate at sfc and sum of &
             &individ terms', 'kg/kg/sec', missing_value=missing_value)

  !-----------------------------------------------------------------------
  !
  !                    COLUMN-INTEGRATED DIAGNOSTICS
  !
  !-----------------------------------------------------------------------

  !------------------------------------------------------------------------
  !    1) variables associated with droplet number and size:
  !------------------------------------------------------------------------
      diag_id%droplets_col = register_diag_field (mod_name,   &
             'droplets_col', axes(1:2), Time,   &
             'Droplet number in cloud column burden', '/cm2',   &
             missing_value=missing_value)
      diag_id%gb_droplets_col = register_diag_field (mod_name,   &
             'gb_droplets_col', axes(1:2), Time,   &
             'Droplet number grid box column burden', '/cm2',   &
             missing_value=missing_value)
      diag_id%droplets_col250 = register_diag_field (mod_name,   &
             'droplets_col250', axes(1:2), Time,   &
             'Droplet number in cloud column burden for T> 250K', '/cm2', &
             missing_value=missing_value)

  !------------------------------------------------------------------------
  !    3) variables associated with cloud and precipitation processes:
  !------------------------------------------------------------------------
      diag_id%sedi_ice = register_diag_field (mod_name, &
             'sedi_ice', axes(1:2), Time, &
             'ice sedimentation flux to surface', 'm/sec',   &
             missing_value=missing_value)
      diag_id%snow_melt_col = register_diag_field (mod_name, &
             'snow_melt_col', axes(1:2), Time, &
             'Column integrated snow melting', 'kg/m2/sec', &
             missing_value=missing_value)
      diag_id%rain_freeze_col = register_diag_field (mod_name, &
             'rain_freeze_col', axes(1:2), Time, &
             'Column integrated loss of rain from freezing', 'kg/m2/sec', &
             missing_value=missing_value)

  !------------------------------------------------------------------------
  !    5) variables associated with ice particle number:
  !------------------------------------------------------------------------
      diag_id%nice_col = register_diag_field (mod_name,   &
             'nice_col', axes(1:2), Time,   &
             'ice number in cloud column burden', '/cm2',  &
             missing_value=missing_value)
      diag_id%gb_nice_col = register_diag_field (mod_name,   &
             'gb_nice_col', axes(1:2), Time,   &
             'ice number grid box column burden', '/cm2',   &
             missing_value=missing_value)

!------------------------------------------------------------------------
!    6) variables associated with precipitation and precipitation area:
!------------------------------------------------------------------------
      diag_id%qrout_col = register_diag_field (mod_name,   &
             'rain_path', axes(1:2), Time,   &
             'rain water path from stratiform', 'kg/m2',    &
             missing_value=missing_value)
      diag_id%qsout_col = register_diag_field (mod_name,    &
             'snow_path', axes(1:2), Time,    &
             'snow water path from stratiform', 'kg/m2',  &
             missing_value=missing_value)

  !-----------------------------------------------------------------------
  !    8) variables associated with cloud liquid time tendency:
  !-----------------------------------------------------------------------
      diag_id%ql_cond_col = register_diag_field (mod_name, &
             'ql_cond_col', axes(1:2), Time, &
             'Column integrated condensation', 'kg/m2/sec',  &
             missing_value=missing_value)
      diag_id%ql_evap_col = register_diag_field (mod_name, &
             'ql_evap_col', axes(1:2), Time, &
             'Column integrated evaporation', 'kg/m2/sec', &
             missing_value=missing_value)
      diag_id%ql_eros_col = register_diag_field (mod_name, &
             'ql_eros_col', axes(1:2), Time, &
             'Column integrated liquid erosion', 'kg/m2/sec',   &
             missing_value=missing_value)
      diag_id%ql_accr_col = register_diag_field (mod_name, &
             'ql_accr_col', axes(1:2), Time, &
             'Column integrated accretion of cloud water by rain',  &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%ql_auto_col = register_diag_field (mod_name, &
             'ql_auto_col', axes(1:2), Time, &
             'Column integrated autoconversion of cloud water to rain', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%ql_berg_col = register_diag_field (mod_name, &
             'ql_berg_col', axes(1:2), Time, &
             'Column integrated Bergeron process', 'kg/m2/sec',   &
             missing_value=missing_value)
      diag_id%ql_freez_col = register_diag_field (mod_name, &
             'ql_freez_col', axes(1:2), Time, &
             'Column integrated homogeneous freezing', 'kg/m2/sec',  &
             missing_value=missing_value)
      diag_id%ql_destr_col = register_diag_field (mod_name, &
             'ql_destr_col', axes(1:2), Time, &
             'Column integrated liquid destruction', 'kg/m2/sec',   &
             missing_value=missing_value)
      diag_id%ql_rime_col = register_diag_field (mod_name, &
             'ql_rime_col', axes(1:2), Time, &
             'Column integrated riming', 'kg/m2/sec',   &
             missing_value=missing_value)
      diag_id%ql_fill_col = register_diag_field (mod_name, &
             'ql_fill_col', axes(1:2), Time, &
             'Column integrated pre-microphysics liquid filler',  &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%liq_adj_col = register_diag_field (mod_name, &
             'liq_adj_col', axes(1:2), Time, &
             'Column integrated liquid water specific humidity  &
             &tendency from supersaturation removal', 'kg/m2/sec',   &
             missing_value=missing_value)
      diag_id%ql_freez2_col = register_diag_field (mod_name, &
             'ql_freez2_col', axes(1:2), Time, &
             'Column integrated heterogeneous freezing', 'kg/m2/sec',   &
             missing_value=missing_value)
      diag_id%ql_sedi_col = register_diag_field (mod_name, &
             'ql_sedi_col', axes(1:2), Time, &
             'Column integrated sedimentation ', 'kg/m2/sec',   &
             missing_value=missing_value)
      diag_id%ql_accrs_col = register_diag_field (mod_name, &
             'ql_accrs_col', axes(1:2), Time, &
             'Column integrated cloud liquid collection by snow', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%ql_bergs_col = register_diag_field (mod_name, &
             'ql_bergs_col', axes(1:2), Time, &
             'Column integrated bergeron process for snow ', 'kg/m2/sec', &
             missing_value=missing_value)

  !-----------------------------------------------------------------------
  !    9) variables associated with cloud droplet number time tendency:
  !-----------------------------------------------------------------------
      diag_id%qn_cond_col = register_diag_field (mod_name, &
             'qn_cond_col', axes(1:2), Time, &
             'Column integrated drop number condensation', '#/m2/sec', &
             missing_value=missing_value)
      diag_id%qn_evap_col = register_diag_field (mod_name, &
             'qn_evap_col', axes(1:2), Time, &
             'Column integrated drop number evaporation', '#/m2/sec',   &
             missing_value=missing_value)
      diag_id%qn_berg_col = register_diag_field (mod_name, &
             'qn_berg_col', axes(1:2), Time, &
             'Column integrated drop number bergeron tendency',   &
             '#/m2/sec', missing_value=missing_value)
      diag_id%qn_fill_col = register_diag_field (mod_name, &
             'qn_fill_col', axes(1:2), Time, &
             'Column integrated drop number filler', '#/m2/sec',   &
             missing_value=missing_value)
      diag_id%qn_destr_col = register_diag_field (mod_name, &
             'qn_destr_col', axes(1:2), Time, &
             'Column integrated drop number destruction', '#/m2/sec', &
             missing_value=missing_value)
      diag_id%qn_super_col = register_diag_field (mod_name, &
             'qn_super_col', axes(1:2), Time, &
             'Column integrated drop number tendency from &
             &supersaturation removal', '#/m2/sec',&
             missing_value=missing_value)
      diag_id%qn_freez_col = register_diag_field (mod_name, &
             'qn_freez_col', axes(1:2), Time, &
             'Column integrated drop number tendency from &
             &heterogeneous freezing', '#/m2/sec',&
             missing_value=missing_value)
      diag_id%qn_sacws_col = register_diag_field (mod_name, &
             'qn_sacws_col', axes(1:2), Time, &
             'Column integrated drop number tendency from droplet &
             &collection by snow', '#/m2/sec',&
             missing_value=missing_value)
      diag_id%qn_sacws_o_col = register_diag_field (mod_name, &
             'qn_sacws_o_col', axes(1:2), Time, &
             'Column integrated drop number tendency from droplet &
             &collection by snow (one ice cat)', '#/m2/sec',&
             missing_value=missing_value)
      diag_id%qn_eros_col = register_diag_field (mod_name, &
             'qn_eros_col', axes(1:2), Time, &
             'Column integrated drop number tendency from cloud erosion', &
             '#/m2/sec', missing_value=missing_value)
      diag_id%qn_pra_col = register_diag_field (mod_name, &
             'qn_pra_col', axes(1:2), Time, &
             'Column integrated drop number tendency from collection &
             &by rain', '#/m2/sec', missing_value=missing_value)
      diag_id%qn_auto_col = register_diag_field (mod_name, &
             'qn_auto_col', axes(1:2), Time, &
             'Column integrated drop number tendency from autoconversion &
             &to rain', '#/m2/sec', missing_value=missing_value)
      diag_id%qn_nucclim_col = register_diag_field (mod_name, &
             'qn_nucclim_col', axes(1:2), Time, &
             'Column integrated drop number tendency from nucleation &
             &limiter', '#/m2/sec', missing_value=missing_value)
      diag_id%qn_sedi_col = register_diag_field (mod_name, &
             'qn_sedi_col', axes(1:2), Time, &
             'Column integrated drop number tendency from sedimentation &
             &of droplets', '#/m2/sec', missing_value=missing_value)
      diag_id%qn_melt_col = register_diag_field (mod_name, &
             'qn_melt_col', axes(1:2), Time, &
             'Column integrated drop number tendency from melting of ice',&
             '#/m2/sec', missing_value=missing_value)
      diag_id%qn_ihom_col = register_diag_field (mod_name, &
             'qn_ihom_col', axes(1:2), Time, &
             'Column integrated drop number tendency from homogeneous &
             &freezing of droplets', '#/m2/sec',&
             missing_value=missing_value)
      diag_id%qn_size_adj_col = register_diag_field (mod_name, &
             'qn_size_adj_col', axes(1:2), Time, &
             'Column integrated drop number tendency from droplet &
             &size adjustment', '#/m2/sec',&
             missing_value=missing_value)
      diag_id%qn_fill2_col = register_diag_field (mod_name, &
             'qn_fill2_col', axes(1:2), Time, &
             'Column integrated drop number tendency from second filler', &
             '#/m2/sec', missing_value=missing_value)
      diag_id%qn_contact_frz_col = register_diag_field (mod_name, &
             'qn_contact_frz_col', axes(1:2), Time, &
             'Column integrated drop number tendency from contact &
             &freezing', '#/m2/sec', missing_value=missing_value)
      diag_id%qn_cleanup_col = register_diag_field (mod_name, &
             'qn_cleanup_col', axes(1:2), Time, &
             'Column integrated drop number tendency from cleanup', &
             '#/m2/sec', missing_value=missing_value)
      diag_id%qn_cleanup2_col = register_diag_field (mod_name, &
             'qn_cleanup2_col', axes(1:2), Time, &
             'Column integrated drop number tendency from second cleanup',&
             '#/m2/sec', missing_value=missing_value)

!-----------------------------------------------------------------------
!    13) variables associated with water vapor tendency:
!-----------------------------------------------------------------------
      diag_id%rain_evap_col = register_diag_field (mod_name, &
             'rain_evap_col', axes(1:2), Time, &
             'Column integrated rain evaporation', 'kg/m2/sec',   &
             missing_value=missing_value)

!-----------------------------------------------------------------------
!    14) variables associated with cloud ice time tendency:
!-----------------------------------------------------------------------
      diag_id%qi_fall_col = register_diag_field (mod_name, &
             'qi_fall_col', axes(1:2), Time, &
             'Column integrated ice settling', 'kg/m2/sec',   &
             missing_value=missing_value)
      diag_id%qi_fill_col = register_diag_field (mod_name, &
             'qi_fill_col', axes(1:2), Time, &
             'Column integrated pre-microphysics ice filler', 'kg/m2/sec',&
             missing_value=missing_value)
      diag_id%qi_dep_col = register_diag_field (mod_name, &
             'qi_dep_col', axes(1:2), Time, &
             'Column integrated large-scale deposition', 'kg/m2/sec',   &
             missing_value=missing_value)
      diag_id%qi_subl_col = register_diag_field (mod_name, &
             'qi_subl_col', axes(1:2), Time, &
             'Column integrated large-scale sublimation', 'kg/m2/sec', &
             missing_value=missing_value)
      diag_id%qi_eros_col = register_diag_field (mod_name, &
             'qi_eros_col', axes(1:2), Time, &
             'Column integrated ice erosion', 'kg/m2/sec',    &
             missing_value=missing_value)
      diag_id%qi_destr_col = register_diag_field (mod_name, &
             'qi_destr_col', axes(1:2), Time, &
             'Column integrated ice destruction', 'kg/m2/sec',   &
             missing_value=missing_value)
      diag_id%qi_qvdep_col = register_diag_field (mod_name, &
             'qi_qvdep_col', axes(1:2), Time, &
             'Column integrated vapor deposition', 'kg/m2/sec',   &
             missing_value=missing_value)
      diag_id%qi_melt_col = register_diag_field (mod_name, &
             'qi_melt_col', axes(1:2), Time, &
             'Column integrated ice melting to rain', 'kg/m2/sec', &
             missing_value=missing_value)
      diag_id%qi_melt2_col = register_diag_field (mod_name, &
             'qi_melt2_col', axes(1:2), Time, &
             'Column integrated ice melting to cloud droplets',  &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%ice_adj_col = register_diag_field (mod_name, &
             'ice_adj_col', axes(1:2), Time, &
             'Column integrated ice water specific humidity tendency &
             &from the removal of supersaturation', 'kg/m2/sec',   &
             missing_value=missing_value)
      diag_id%qi_auto_col = register_diag_field (mod_name, &
             'qi_auto_col', axes(1:2), Time, &
             'Column integrated autoconversion from cloud ice to snow ', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%qi_accr_col = register_diag_field (mod_name, &
             'qi_accr_col', axes(1:2), Time, &
             'Column integrated accretion of cloud ice by snow ', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%qi_accrs_col = register_diag_field (mod_name, &
             'qi_accrs_col', axes(1:2), Time, &
             'Column integrated self collection (one class scheme)',   &
             'kg/m2/sec', missing_value=missing_value)

!-----------------------------------------------------------------------
!    15) variables associated with cloud area time tendency:
!-----------------------------------------------------------------------
      diag_id%qa_lsform_col = register_diag_field (mod_name, &
             'qa_lsform_col', axes(1:2), Time, &
             'Column integrated cloud fraction tendency from &
             &condensation', 'kg/m2/sec', missing_value=missing_value)
      diag_id%qa_lsdiss_col = register_diag_field (mod_name, &
             'qa_lsdiss_col', axes(1:2), Time, &
             'Column integrated cloud fraction tendency from dissipation',&
             'kg/m2/sec', missing_value=missing_value)
      diag_id%qa_rhred_col = register_diag_field (mod_name, &
             'qa_rhred_col', axes(1:2), Time, & 
             'Column integrated cloud fraction tendency from RH limiter', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%qa_eros_col = register_diag_field (mod_name, &
             'qa_eros_col', axes(1:2), Time, &
             'Column integrated cloud fraction tendency from erosion', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%qa_fill_col = register_diag_field (mod_name, &
             'qa_fill_col', axes(1:2), Time, &
             'Column integrated cloud fraction tendency from filler', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%qa_super_col = register_diag_field (mod_name, &
             'qa_super_col', axes(1:2), Time, &
             'Column integrated cloud fraction tendency from &
             &supersaturation removal', 'kg/m2/sec',    &
             missing_value=missing_value)
      diag_id%qa_destr_col = register_diag_field (mod_name, &
             'qa_destr_col', axes(1:2), Time, &
             'Column integrated cloud fraction tendency from destruction',&
             'kg/m2/sec', missing_value=missing_value)
      diag_id%qa_limits_col = register_diag_field (mod_name, &
             'qa_limits_col', axes(1:2), Time, &
             'Column integrated cloud fraction tendency from limits &
             &imposition', 'kg/m2/sec', missing_value=missing_value)
      diag_id%qa_ahuco_col = register_diag_field (mod_name, &
             'qa_ahuco_col', axes(1:2), Time, &
             'Column integrated cloud fraction tendency from conv &
             &area restriction', 'kg/m2/sec', &
             missing_value=missing_value)

!------------------------------------------------------------------------
!    16) variables added by h1g with ncar M-G microphysics
!------------------------------------------------------------------------
      diag_id%SA2d =    register_diag_field ( mod_name, &
             'SA2d',  axes(1:2), Time, &
             'column integrated total cld fraction tendency', &
             ' kg/ (m2 sec)', missing_value=missing_value)
      diag_id%ST2d =    register_diag_field ( mod_name, &
             'ST2d',  axes(1:2), Time, &
             'column integrated total temperature tendency', &
             'K kg/ (m2 sec)', missing_value=missing_value)
      diag_id%SQ2d =    register_diag_field ( mod_name, &
             'SQ2d',  axes(1:2), Time, &
             'Column integrated total Water vapor specific &
             &humidity tendency ', 'kg/m2/sec', &
             missing_value=missing_value)
      diag_id%SL2d = register_diag_field ( mod_name, &
             'SL2d',  axes(1:2), Time, &
             'Column integrated total Liquid Water specific humidity', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%SI2d =    register_diag_field ( mod_name, &
             'SI2d',  axes(1:2), Time, &
             'column integrated total Ice water specific humidity &
             &tendency ', 'kg/m2/sec', missing_value=missing_value)
      diag_id%SN2d =    register_diag_field ( mod_name, &
             'SN2d',  axes(1:2), Time, &
             'Column integrated total cloud drop number tendency', &
             '#/m2/sec', missing_value=missing_value)
      diag_id%SNI2d = register_diag_field ( mod_name, &
             'SNI2d',  axes(1:2), Time, &
             'Column integrated total ice crystal number tendency', &
             '#/m2/sec', missing_value=missing_value)
      diag_id%srfrain_accrs_col=  register_diag_field ( mod_name, &
             'srfrain_accrs_col', axes(1:2), Time, &
             'Column integrated rain water sink from collection by snow', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%srfrain_freez_col=  register_diag_field ( mod_name, &
             'srfrain_freez_col', axes(1:2), Time, &
             'Column integrated rain water sink from freezing', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%srfrain_evap_col=  register_diag_field ( mod_name, &
             'srfrain_evap_col', axes(1:2), Time, &
             'Column integrated rain water sink from evaporation', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%snow_mass_conv = register_diag_field ( mod_name, &
             'snow_mass_conv', axes(1:2), Time, &
             'Change necessary to surface snowfall to ensure h2o &
             & balance in column', 'kg/m2/sec',   &
             missing_value=missing_value)
      diag_id%neg_snow = register_diag_field ( mod_name, &
             'neg_snow', axes(1:2), Time, &
             'magnitude of negative snow removed in column', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%rain_mass_conv = register_diag_field ( mod_name, &
             'rain_mass_conv', axes(1:2), Time, &
             'Change necessary to surface rainfall to ensure h2o&
             & balance in column', 'kg/m2/sec',  &
             missing_value=missing_value)
      diag_id%neg_rain = register_diag_field ( mod_name, &
             'neg_rain', axes(1:2), Time, &
             'magnitude of negative rain removed in column', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_liquid_init_col  =    register_diag_field ( mod_name, &
             'q_liquid_init_col',  axes(1:2), Time, &
             'Column integrated Water vapor specific humidity tendency &
             &from initial adjustment of liquid cloud', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_ice_init_col  =    register_diag_field ( mod_name, &
             'q_ice_init_col',  axes(1:2), Time, &
             'Column integrated Water vapor specific humidity tendency &
             &from initial adjustment of ice cloud', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_rain_evap_col  =    register_diag_field ( mod_name, &
             'q_rain_evap_col',  axes(1:2), Time, &
             'Column integrated Water vapor specific humidity tendency &
             &from rain evaporation', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_cond_col  =    register_diag_field ( mod_name, &
             'q_cond_col',  axes(1:2), Time, &
             'Column integrated Water vapor specific humidity tendency &
             &from condensation of vapor', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_deposition_col  =    register_diag_field ( mod_name, &
             'q_deposition_col',  axes(1:2), Time, &
             'Column integrated Water vapor specific humidity tendency &
             &from deposition of vapor onto ice',   &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_eros_l_col  =    register_diag_field ( mod_name, &
             'q_eros_l_col',  axes(1:2), Time, &
             'Column integrated Water vapor specific humidity tendency &
             &from erosion of liquid cloud', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_eros_i_col  =    register_diag_field ( mod_name, &
             'q_eros_i_col',  axes(1:2), Time, &
             'Column integrated Water vapor specific humidity tendency &
             &from erosion of ice cloud', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_qv_on_qi_col  =    register_diag_field ( mod_name, &
             'q_qv_on_qi_col',  axes(1:2), Time, &
             'Column integrated Water vapor specific humidity tendency &
             &from direct deposition on ice', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_sedi_ice2vapor_col  =    register_diag_field ( mod_name, &
             'q_sedi_ice2vapor_col',  axes(1:2), Time, &
             'Column integrated Water vapor specific humidity tendency &
             &from the sublimation of falling cloud ice', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_sedi_liquid2vapor_col  =  register_diag_field ( mod_name, &
             'q_sedi_liquid2vapor_col',  axes(1:2), Time, &
             'Column integrated Water vapor specific humidity tendency &
             &from the evaporation of falling cloud liquid', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_super_sat_rm_col  =    register_diag_field ( mod_name, &
             'q_super_sat_rm_col',  axes(1:2), Time, &
             'Column integrated Water vapor specific humidity tendency &
             &from removal of supersaturation', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_destr_col  =    register_diag_field ( mod_name, &
             'q_destr_col',  axes(1:2), Time, &
             'Column integrated Water vapor specific humidity tendency &
             &from cloud destruction', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_cleanup_liquid_col  =    register_diag_field ( mod_name, &
             'q_cleanup_liquid_col',  axes(1:2), Time, &
             'Column integrated Water vapor specific humidity tendency &
             &from cleaning up liquid cloud', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_cleanup_ice_col  =    register_diag_field ( mod_name, &
             'q_cleanup_ice_col',  axes(1:2), Time, &
             'Column integrated Water vapor specific humidity tendency &
             &from cleaning up ice cloud', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_snow_sublim_col  =    register_diag_field ( mod_name, &
             'q_snow_sublim_col',  axes(1:2), Time, &
             'Column integrated Water vapor specific humidity tendency &
             &from sublimation of snow', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%q_snow2vapor_col  =    register_diag_field ( mod_name, &
             'q_snow2vapor_col',  axes(1:2), Time, &
             'Column integrated water vapor specific humidity tendency &
             &from condensation resulting from cooling of the air &
             &by melting snow', 'kg/m2/sec', missing_value=missing_value)
      diag_id%ql_HM_splinter_col = register_diag_field ( mod_name, &
             'ql_HM_splinter_col', axes(1:2), Time, &
             'Column integrated liquid water specific humidity tendency &
             &from H-M splintering ',  'kg/m2/sec', &
             missing_value=missing_value)

!------------------------------------------------------------------------
!   17)  variables associated with budget analysis
!------------------------------------------------------------------------
      diag_id%SA_imb_col =    register_diag_field ( mod_name, &
             'SA_imb_col',  axes(1:2), Time, &
             'difference between  column-integrated qa tendency and &
             &sum of individ terms', 'kg/ (m2 sec)', &
             missing_value=missing_value)
      diag_id%ST_imb_col =    register_diag_field ( mod_name, &
             'ST_imb_col',  axes(1:2), Time, &
             'difference between column-integrated temp tendency and &
             &sum of individ terms', 'K kg/(m2 sec)', &
             missing_value=missing_value)
      diag_id%SQ_imb_col =    register_diag_field ( mod_name, &
             'SQ_imb_col',  axes(1:2), Time, &
             'difference between column-integrated qv tendency and &
             &sum of individ terms', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%SL_imb_col =    register_diag_field ( mod_name, &
             'SL_imb_col',  axes(1:2), Time, &
             'difference between column-integrated ql tendency and &
             &sum of individ terms', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%SI_imb_col =    register_diag_field ( mod_name, &
             'SI_imb_col',  axes(1:2), Time, &
             'difference between column-integrated qi tendency and &
             &sum of individ terms', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%SN_imb_col =    register_diag_field ( mod_name, &
             'SN_imb_col',  axes(1:2), Time, &
             'difference between column-integrated qn tendency and &
             &sum of individ terms', &
             '#/m2/sec', missing_value=missing_value)
      diag_id%SNi_imb_col =    register_diag_field ( mod_name, &
             'SNi_imb_col',  axes(1:2), Time, &
             'difference between column-integrated qni tendency and &
             &sum of individ terms', &
             '#/m2/sec', missing_value=missing_value)
      diag_id%rain_imb_col =    register_diag_field ( mod_name, &
             'rain_imb_col',  axes(1:2), Time, &
             'difference between column-integrated rainfall rate and &
             &sum of individ terms', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%cld_liq_imb_col =    register_diag_field ( mod_name, &
             'cld_liq_imb_col',  axes(1:2), Time, &
             'difference between column-integrated ql fallout rate and &
             &sum of individ terms', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%snow_imb_col =    register_diag_field ( mod_name, &
             'snow_imb_col',  axes(1:2), Time, &
             'difference between column-integrated snowfall rate and &
             &sum of individ terms', &
             'kg/m2/sec', missing_value=missing_value)
      diag_id%cld_ice_imb_col =    register_diag_field ( mod_name, &
             'cld_ice_imb_col',  axes(1:2), Time, &
             'difference between column-integrated qi fallout rate and &
             &sum of individ terms', &
             'kg/m2/sec', missing_value=missing_value)

!-----------------------------------------------------------------------
!    17) variables associated with ice particle number time tendency:
!-----------------------------------------------------------------------
      diag_id%qni_fill_col = register_diag_field (mod_name, &
             'qni_fill_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from filler',&
             '#/m2/sec', missing_value=missing_value)  
      diag_id%qni_nnuccd_col = register_diag_field (mod_name, &
             'qni_nnuccd_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from &
             &nucleation', '#/m2/sec', missing_value=missing_value)  
      diag_id%qni_nsubi_col = register_diag_field (mod_name, &
             'qni_nsubi_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from &
             &sublimation', '#/m2/sec', missing_value=missing_value)  
      diag_id%qni_nerosi_col = register_diag_field (mod_name, &
             'qni_nerosi_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from &
             &erosion', '#/m2/sec', missing_value=missing_value) 
      diag_id%qni_nprci_col = register_diag_field (mod_name, &
             'qni_nprci_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from &
             &autoconversion', '#/m2/sec', missing_value=missing_value) 
      diag_id%qni_nprai_col = register_diag_field (mod_name, &
             'qni_nprai_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from &
             &accretion by snow', '#/m2/sec', missing_value=missing_value) 
      diag_id%qni_nucclim1_col = register_diag_field (mod_name, &
             'qni_nucclim1_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from first &
             &nucleation limiter', '#/m2/sec', missing_value=missing_value)
      diag_id%qni_nucclim2_col = register_diag_field (mod_name, &
             'qni_nucclim2_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from second &
             &nucleation limiter','#/m2/sec', missing_value=missing_value) 
      diag_id%qni_sedi_col = register_diag_field (mod_name, &
             'qni_sedi_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from &
             &sedimentation', '#/m2/sec', missing_value=missing_value) 
      diag_id%qni_melt_col = register_diag_field (mod_name, &
             'qni_melt_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from &
             &melting', '#/m2/sec', missing_value=missing_value) 
      diag_id%qni_size_adj_col = register_diag_field (mod_name, &
             'qni_size_adj_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from size &
             &adjustment', '#/m2/sec', missing_value=missing_value) 
      diag_id%qni_fill2_col = register_diag_field (mod_name, &
             'qni_fill2_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from second &
             &filler', '#/m2/sec', missing_value=missing_value)  
      diag_id%qni_super_col = register_diag_field (mod_name, &
             'qni_super_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from &
             &removal of supersaturation', '#/m2/sec',   &
             missing_value=missing_value)  
      diag_id%qni_ihom_col = register_diag_field (mod_name, &
             'qni_ihom_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from &
             &homogeneous freezing', '#/m2/sec', &
             missing_value=missing_value)  
      diag_id%qni_destr_col = register_diag_field (mod_name, &
             'qni_destr_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from cloud &
             &destruction', '#/m2/sec', missing_value=missing_value)  
      diag_id%qni_cleanup_col = register_diag_field (mod_name, &
             'qni_cleanup_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from &
             &cleanup', '#/m2/sec', missing_value=missing_value)  
      diag_id%qni_cleanup2_col = register_diag_field (mod_name, &
             'qni_cleanup2_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from &
             &cleanup2', '#/m2/sec', missing_value=missing_value)  
      diag_id%qni_nsacwi_col = register_diag_field (mod_name, &
             'qni_nsacwi_col', axes(1:2), Time, &
             'Column integrated ice particle number tendency from HM &
             &ice multiplication', '#/m2/sec', missing_value=missing_value)

!------------------------------------------------------------------------
!    determine the number of activated diagnostic variables on both
!    full and half levels in the vertical. thesde values will be used
!    to allocate arrays containing the activated diagnostics.
!------------------------------------------------------------------------
      n_diag_4d = 1
      n_diag_4d_kp1 = 1

!------------------------------------------------------------------------
!    1) variables associated with droplet number and size:
!------------------------------------------------------------------------
      if (diag_id%droplets  > 0) then
        diag_pt%droplets = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%droplets_wtd  > 0) then
        diag_pt%droplets_wtd = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%rvolume  > 0 ) then
        diag_pt%rvolume = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%droplets_col  > 0) then
        diag_pt%droplets_col = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%droplets_col250  > 0) then
        diag_pt%droplets_col250 = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%gb_droplets_col  > 0) then
        diag_pt%gb_droplets_col = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if

!------------------------------------------------------------------------
!    2) variables associated with cloud liquid content:
!------------------------------------------------------------------------
      if (diag_id%ql_wt  > 0) then
        diag_pt%ql_wt = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if

!------------------------------------------------------------------------
!    3) variables associated with cloud and precipitation processes:
!------------------------------------------------------------------------
      if (diag_id%lsf_strat > 0) then
        diag_pt%lsf_strat = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%dcond  > 0 ) then
        diag_pt%dcond = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%aauto  > 0) then
        diag_pt%aauto = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%vfall  > 0) then
        diag_pt%vfall = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%delta_cf  > 0 ) then
        diag_pt%delta_cf = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%cf_liq_init    > 0 ) then
        diag_pt%cf_liq_init   = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%subgrid_w_variance  > 0 ) then
        diag_pt%subgrid_w_variance = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%potential_droplets  > 0 ) then
        diag_pt%potential_droplets = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%potential_crystals  > 0 ) then
        diag_pt%potential_crystals = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%dust_berg_flag  > 0 ) then
        diag_pt%dust_berg_flag = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%cf_ice_init    > 0 ) then
        diag_pt%cf_ice_init = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%snow_melt  + diag_id%snow_melt_col > 0) then
        diag_pt%snow_melt = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%rain_freeze  + diag_id%rain_freeze_col > 0) then
        diag_pt%rain_freeze = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if

!------------------------------------------------------------------------
!    4) variables associated with model convection:
!------------------------------------------------------------------------
      if (diag_id%lcf_strat > 0) then
        diag_pt%lcf_strat = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%mfls_strat  > 0) then
        diag_pt%mfls_strat = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if

!------------------------------------------------------------------------
!    5) variables associated with ice particle number:
!------------------------------------------------------------------------
      if (diag_id%nice  > 0) then
        diag_pt%nice = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%nice_col  > 0) then
        diag_pt%nice_col = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%gb_nice_col  > 0) then
        diag_pt%gb_nice_col = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if

!------------------------------------------------------------------------
!    6) variables associated with precipitation and precipitation area:
!------------------------------------------------------------------------
      if (diag_id%qrout  + diag_id%qrout_col > 0) then
        diag_pt%qrout = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qsout  + diag_id%qsout_col > 0) then
        diag_pt%qsout = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%rain3d  > 0) then
        diag_pt%rain3d = n_diag_4d_kp1
        n_diag_4d_kp1 = n_diag_4d_kp1 + 1 
      end if
      if (diag_id%snow3d  > 0) then
        diag_pt%snow3d = n_diag_4d_kp1
        n_diag_4d_kp1 = n_diag_4d_kp1 + 1 
      end if
      if (diag_id%rain_clr  > 0) then
        diag_pt%rain_clr = n_diag_4d_kp1
        n_diag_4d_kp1 = n_diag_4d_kp1 + 1 
      end if
      if (diag_id%rain_cld  > 0) then
        diag_pt%rain_cld = n_diag_4d_kp1
        n_diag_4d_kp1 = n_diag_4d_kp1 + 1 
      end if
      if (diag_id%a_rain_clr  > 0) then
        diag_pt%a_rain_clr = n_diag_4d_kp1
        n_diag_4d_kp1 = n_diag_4d_kp1 + 1 
      end if
      if (diag_id%a_rain_cld  > 0) then
        diag_pt%a_rain_cld = n_diag_4d_kp1
        n_diag_4d_kp1 = n_diag_4d_kp1 + 1 
      end if
      if (diag_id%a_precip_clr  > 0) then
        diag_pt%a_precip_clr = n_diag_4d_kp1
        n_diag_4d_kp1 = n_diag_4d_kp1 + 1 
      end if
      if (diag_id%a_precip_cld  > 0) then
        diag_pt%a_precip_cld = n_diag_4d_kp1
        n_diag_4d_kp1 = n_diag_4d_kp1 + 1 
      end if
      if (diag_id%snow_clr  > 0) then
        diag_pt%snow_clr = n_diag_4d_kp1
        n_diag_4d_kp1 = n_diag_4d_kp1 + 1 
      end if
      if (diag_id%snow_cld  > 0) then
        diag_pt%snow_cld = n_diag_4d_kp1
        n_diag_4d_kp1 = n_diag_4d_kp1 + 1 
      end if
      if (diag_id%a_snow_clr  > 0) then
        diag_pt%a_snow_clr = n_diag_4d_kp1
        n_diag_4d_kp1 = n_diag_4d_kp1 + 1 
      end if
      if (diag_id%a_snow_cld  > 0) then
        diag_pt%a_snow_cld = n_diag_4d_kp1
        n_diag_4d_kp1 = n_diag_4d_kp1 + 1 
      end if

!------------------------------------------------------------------------
!    7) variables associated with cloud fraction:
!------------------------------------------------------------------------
      if (diag_id%aall  > 0) then
        diag_pt%aall = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%aliq  > 0) then
        diag_pt%aliq = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%aice  > 0) then
        diag_pt%aice = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%cfin > 0 ) then
        diag_pt%cfin = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if

!-----------------------------------------------------------------------
!    8) variables associated with cloud liquid time tendency:
!-----------------------------------------------------------------------
      if (diag_id%qldt_cond  + diag_id%ql_cond_col > 0) then
        diag_pt%qldt_cond = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qldt_evap  + diag_id%ql_evap_col > 0) then
        diag_pt%qldt_evap = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qldt_eros  + diag_id%ql_eros_col > 0) then
        diag_pt%qldt_eros = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qldt_accr  + diag_id%ql_accr_col > 0) then
        diag_pt%qldt_accr = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qldt_auto  + diag_id%ql_auto_col > 0) then
        diag_pt%qldt_auto = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%liq_adj  + diag_id%liq_adj_col  +   &
          diag_id%ice_adj + diag_id%ice_adj_col  > 0) then
        diag_pt%liq_adj = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
        diag_pt%ice_adj = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qldt_fill  + diag_id%ql_fill_col > 0) then
        diag_pt%qldt_fill = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qldt_berg  + diag_id%ql_berg_col > 0) then
        diag_pt%qldt_berg = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qldt_freez + diag_id%ql_freez_col > 0) then
        diag_pt%qldt_freez = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qldt_rime  + diag_id%ql_rime_col > 0) then
        diag_pt%qldt_rime = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qldt_destr + diag_id%ql_destr_col > 0) then
        diag_pt%qldt_destr = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qldt_freez2 + diag_id%ql_freez2_col > 0) then
        diag_pt%qldt_freez2 = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qldt_sedi + diag_id%ql_sedi_col > 0) then
        diag_pt%qldt_sedi = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%sedi_sfc  > 0) then
        diag_pt%sedi_sfc = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qldt_accrs + diag_id%ql_accrs_col > 0) then
        diag_pt%qldt_accrs = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qldt_bergs + diag_id%ql_bergs_col > 0) then
        diag_pt%qldt_bergs = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if

!-----------------------------------------------------------------------
!    9) variables associated with cloud droplet number time tendency:
!-----------------------------------------------------------------------
      if (diag_id%qndt_cond + diag_id%qn_cond_col > 0) then
        diag_pt%qndt_cond = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_evap  + diag_id%qn_evap_col > 0) then
        diag_pt%qndt_evap = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_fill + diag_id%qn_fill_col +    &
          diag_id%qldt_fill + diag_id%ql_fill_col > 0) then
        diag_pt%qndt_fill = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_destr + diag_id%qn_destr_col > 0) then
        diag_pt%qndt_destr = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_super + diag_id%qn_super_col > 0) then
        diag_pt%qndt_super = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_berg + diag_id%qn_berg_col +    &
                              diag_id%qldt_berg > 0) then
        diag_pt%qndt_berg = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_freez + diag_id%qn_freez_col > 0) then
        diag_pt%qndt_freez = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_sacws + diag_id%qn_sacws_col > 0) then
        diag_pt%qndt_sacws = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_sacws_o + diag_id%qn_sacws_o_col > 0) then
        diag_pt%qndt_sacws_o = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_eros + diag_id%qn_eros_col > 0) then
        diag_pt%qndt_eros = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_pra + diag_id%qn_pra_col > 0) then
        diag_pt%qndt_pra = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_auto + diag_id%qn_auto_col > 0) then
        diag_pt%qndt_auto = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_nucclim + diag_id%qn_nucclim_col > 0) then
        diag_pt%qndt_nucclim = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_sedi + diag_id%qn_sedi_col > 0) then
        diag_pt%qndt_sedi = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_melt + diag_id%qn_melt_col > 0) then
        diag_pt%qndt_melt = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_ihom + diag_id%qn_ihom_col > 0) then
        diag_pt%qndt_ihom = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_size_adj + diag_id%qn_size_adj_col > 0) then
        diag_pt%qndt_size_adj = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_fill2 + diag_id%qn_fill2_col > 0) then
        diag_pt%qndt_fill2 = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if

!-----------------------------------------------------------------------
!    10) variables associated with ice particle number time tendency:
!-----------------------------------------------------------------------
      if (diag_id%qnidt_fill + diag_id%qni_fill_col > 0) then
        diag_pt%qnidt_fill = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_nnuccd + diag_id%qni_nnuccd_col > 0) then
        diag_pt%qnidt_nnuccd = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_nsubi + diag_id%qni_nsubi_col > 0) then
        diag_pt%qnidt_nsubi = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_nerosi + diag_id%qni_nerosi_col > 0) then
        diag_pt%qnidt_nerosi = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_nprci + diag_id%qni_nprci_col > 0) then
        diag_pt%qnidt_nprci = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_nprai + diag_id%qni_nprai_col > 0) then
        diag_pt%qnidt_nprai = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_nucclim1 + diag_id%qni_nucclim1_col > 0) then
        diag_pt%qnidt_nucclim1 = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_nucclim2 + diag_id%qni_nucclim2_col > 0) then
        diag_pt%qnidt_nucclim2 = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_sedi + diag_id%qni_sedi_col > 0) then
        diag_pt%qnidt_sedi = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_melt + diag_id%qni_melt_col > 0) then
        diag_pt%qnidt_melt = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_size_adj + diag_id%qni_size_adj_col > 0) then
        diag_pt%qnidt_size_adj = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_fill2 + diag_id%qni_fill2_col > 0) then
        diag_pt%qnidt_fill2 = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_super + diag_id%qni_super_col > 0) then
        diag_pt%qnidt_super = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_ihom + diag_id%qni_ihom_col > 0) then
        diag_pt%qnidt_ihom = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_destr + diag_id%qni_destr_col > 0) then
        diag_pt%qnidt_destr = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_cleanup + diag_id%qni_cleanup_col > 0) then
        diag_pt%qnidt_cleanup = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_cleanup2 + diag_id%qni_cleanup2_col > 0) then
        diag_pt%qnidt_cleanup2 = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if

!-----------------------------------------------------------------------
!    11) variables associated with relative humidity:
!-----------------------------------------------------------------------
      if (diag_id%rhcrit > 0 ) then
        diag_pt%rhcrit = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%rhcrit_min > 0 ) then
        diag_pt%rhcrit_min = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%rhiin > 0 ) then
        diag_pt%rhiin = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%rhlin > 0 ) then
        diag_pt%rhlin = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if

!-----------------------------------------------------------------------
!    12) variables associated with aerosol nucleation:
!-----------------------------------------------------------------------
      if (diag_id%imass7 > 0 ) then
        diag_pt%imass7 = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%ni_dust > 0 ) then
        diag_pt%ni_dust = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%ni_sulf> 0 ) then
        diag_pt%ni_sulf = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%ni_bc > 0 ) then
        diag_pt%ni_bc = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%ndust1 > 0 ) then
        diag_pt%ndust1 = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%ndust2 > 0 ) then
        diag_pt%ndust2 = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%ndust3 > 0 ) then
        diag_pt%ndust3 = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%ndust4 > 0 ) then
        diag_pt%ndust4 = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%ndust5 > 0 ) then
        diag_pt%ndust5 = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%sulfate > 0) then
        diag_pt%sulfate = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%seasalt_sub  > 0) then
        diag_pt%seasalt_sub = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%seasalt_sup  > 0) then
        diag_pt%seasalt_sup = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%om  > 0) then
        diag_pt%om = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if

!-----------------------------------------------------------------------
!    13) variables associated with water vapor tendency:
!-----------------------------------------------------------------------
      if (diag_id%rain_evap + diag_id%rain_evap_col > 0) then
        diag_pt%rain_evap = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if

!-----------------------------------------------------------------------
!    14) variables associated with cloud ice time tendency:
!-----------------------------------------------------------------------
      if (diag_id%qidt_dep + diag_id%qi_dep_col > 0) then
        diag_pt%qidt_dep = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qidt_subl + diag_id%qi_subl_col > 0) then
        diag_pt%qidt_subl = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qidt_eros + diag_id%qi_eros_col > 0) then
        diag_pt%qidt_eros = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%sedi_ice  > 0) then
        diag_pt%sedi_ice = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qidt_fall + diag_id%qi_fall_col > 0) then
        diag_pt%qidt_fall = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qidt_melt + diag_id%qi_melt_col > 0) then
        diag_pt%qidt_melt = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qidt_melt2 + diag_id%qi_melt2_col > 0) then
        diag_pt%qidt_melt2 = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qidt_destr + diag_id%qi_destr_col > 0) then
        diag_pt%qidt_destr = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qidt_qvdep + diag_id%qi_qvdep_col > 0) then
        diag_pt%qidt_qvdep = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qidt_fill + diag_id%qi_fill_col > 0) then
        diag_pt%qidt_fill = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qidt_auto + diag_id%qi_auto_col > 0) then
        diag_pt%qidt_auto = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qidt_accr + diag_id%qi_accr_col > 0) then
        diag_pt% qidt_accr= n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qidt_accrs + diag_id%qi_accrs_col > 0) then
        diag_pt%qidt_accrs = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if

!-----------------------------------------------------------------------
!    15) variables associated with cloud area time tendency:
!-----------------------------------------------------------------------
      if (diag_id%qadt_lsform + diag_id%qa_lsform_col > 0) then
        diag_pt%qadt_lsform = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qadt_lsdiss + diag_id%qa_lsdiss_col > 0) then
        diag_pt%qadt_lsdiss = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qadt_rhred + diag_id%qa_rhred_col > 0) then
        diag_pt%qadt_rhred = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qadt_eros  + diag_id%qa_eros_col > 0) then
        diag_pt%qadt_eros = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qadt_fill + diag_id%qa_fill_col  > 0) then
        diag_pt%qadt_fill = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qadt_super + diag_id%qa_super_col > 0) then
        diag_pt%qadt_super = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qadt_destr + diag_id%qa_destr_col > 0) then
        diag_pt%qadt_destr = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qadt_limits + diag_id%qa_limits_col > 0) then
        diag_pt%qadt_limits = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qadt_ahuco + diag_id%qa_ahuco_col > 0) then
        diag_pt%qadt_ahuco = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if

!REV#9
!------------------------------------------------------------------------
!   16)  variables added by h1g with ncar M-G microphysics
!------------------------------------------------------------------------
      if (diag_id%SA3d  + diag_id%SA2d  > 0) then
        diag_pt%SA3d = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%ST3d  + diag_id%ST2d  > 0) then
        diag_pt%ST3d = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%SQ3d  + diag_id%SQ2d  > 0) then
        diag_pt%SQ3d = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%SL3d  + diag_id%SL2d  > 0) then
        diag_pt%SL3d = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%SI3d  + diag_id%SI2d  > 0) then
        diag_pt%SI3d = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%SN3d  + diag_id%SN2d  > 0) then
        diag_pt%SN3d = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%SNI3d  + diag_id%SNI2d  > 0) then
        diag_pt%SNI3d = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_contact_frz + diag_id%qn_contact_frz_col > 0) then
        diag_pt%qndt_contact_frz = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_cleanup + diag_id%qn_cleanup_col     > 0) then
        diag_pt%qndt_cleanup     = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qndt_cleanup2 + diag_id%qn_cleanup2_col    > 0) then
        diag_pt%qndt_cleanup2    = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qnidt_nsacwi +  diag_id%qni_nsacwi_col      > 0) then
        diag_pt%qnidt_nsacwi      = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_liquid_init + diag_id%q_liquid_init_col > 0) then
        diag_pt%qdt_liquid_init  = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_ice_init +  diag_id%q_ice_init_col    > 0) then
        diag_pt%qdt_ice_init    = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_rain_evap   + diag_id%q_rain_evap_col   > 0) then
        diag_pt%qdt_rain_evap = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_cond + diag_id%q_cond_col        > 0) then
        diag_pt%qdt_cond        = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_deposition + diag_id%q_deposition_col  > 0) then
        diag_pt%qdt_deposition  = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_eros_l + diag_id%q_eros_l_col      > 0) then
        diag_pt%qdt_eros_l      = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_eros_i  + diag_id%q_eros_i_col      > 0) then
        diag_pt%qdt_eros_i      = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_qv_on_qi + diag_id%q_qv_on_qi_col    > 0) then
        diag_pt%qdt_qv_on_qi    = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_sedi_ice2vapor +   &
                             diag_id%q_sedi_ice2vapor_col > 0) then
        diag_pt%qdt_sedi_ice2vapor = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_sedi_liquid2vapor +   &
                  diag_id%q_sedi_liquid2vapor_col > 0) then
        diag_pt%qdt_sedi_liquid2vapor = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_super_sat_rm + diag_id%q_super_sat_rm_col   > 0) then
        diag_pt%qdt_super_sat_rm   = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_destr + diag_id%q_destr_col          > 0) then
        diag_pt%qdt_destr          = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_cleanup_liquid +    &
                              diag_id%q_cleanup_liquid_col > 0) then
        diag_pt%qdt_cleanup_liquid = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_cleanup_ice +  diag_id%q_cleanup_ice_col    > 0) then
        diag_pt%qdt_cleanup_ice    = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%srfrain_evap +  diag_id%srfrain_evap_col   > 0) then
        diag_pt%srfrain_evap       = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%rain_mass_conv     > 0) then
        diag_pt%rain_mass_conv     = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%neg_rain     > 0) then
        diag_pt%neg_rain     = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%snow_mass_conv     > 0) then
        diag_pt%snow_mass_conv     = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%neg_snow     > 0) then
        diag_pt%neg_snow     = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_snow_sublim + diag_id%q_snow_sublim_col > 0) then
        diag_pt%qdt_snow_sublim = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qdt_snow2vapor       + diag_id%q_snow2vapor_col > 0) then
        diag_pt%qdt_snow2vapor       = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%srfrain_accrs + diag_id%srfrain_accrs_col > 0) then
        diag_pt%srfrain_accrs        = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%srfrain_freez + diag_id%srfrain_freez_col > 0) then
        diag_pt%srfrain_freez        = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%qldt_HM_splinter + diag_id%ql_HM_splinter_col > 0) then
        diag_pt%qldt_HM_splinter      = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
!------------------------------------------------------------------------
!   17)  variables associated with budget analysis
!------------------------------------------------------------------------
      if (diag_id%SA_imb  + diag_id%SA_imb_col > 0) then
        diag_pt%SA_imb      = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%ST_imb  + diag_id%ST_imb_col > 0) then
        diag_pt%ST_imb      = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%SQ_imb  + diag_id%SQ_imb_col > 0) then
        diag_pt%SQ_imb      = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%SL_imb  + diag_id%SL_imb_col > 0) then
        diag_pt%SL_imb      = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%SI_imb  + diag_id%SI_imb_col > 0) then
        diag_pt%SI_imb      = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%SN_imb  + diag_id%SN_imb_col > 0) then
        diag_pt%SN_imb      = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%SNi_imb  + diag_id%SNi_imb_col > 0) then
        diag_pt%SNi_imb      = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%rain_imb  + diag_id%rain_imb_col > 0) then
        diag_pt%rain_imb      = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%cld_liq_imb  + diag_id%cld_liq_imb_col > 0) then
        diag_pt%cld_liq_imb      = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%snow_imb  + diag_id%snow_imb_col > 0) then
        diag_pt%snow_imb      = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if
      if (diag_id%cld_ice_imb  + diag_id%cld_ice_imb_col > 0) then
        diag_pt%cld_ice_imb      = n_diag_4d
        n_diag_4d = n_diag_4d + 1 
      end if

!----------------------------------------------------------------------


end subroutine diag_field_init



!########################################################################




end module strat_netcdf_mod

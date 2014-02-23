                 module cloudrad_diagnostics_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  ds
! </REVIEWER>
! <OVERVIEW>
!    cloudrad_diagnostics_mod generates any desired netcdf output
!    fields involving the cloud fields seen by the radiation package
!    or the cloud radiation interaction variables.
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>

! shared modules:

use mpp_mod,                 only: input_nml_file
use fms_mod,                 only: fms_init, open_namelist_file, &
                                   write_version_number, mpp_pe, &
                                   mpp_root_pe, stdlog, file_exist,  &
                                   check_nml_error, error_mesg,   &
                                   FATAL, NOTE, close_file
use time_manager_mod,        only: time_type, time_manager_init,  &
                                   operator(>)
use diag_manager_mod,        only: register_diag_field, send_data, &
                                   diag_manager_init
use constants_mod,           only: diffac, GRAV, RDGAS

! shared radiation package modules:

use rad_utilities_mod,       only: rad_utilities_init, &
                                   cldrad_properties_type, &
                                   cld_specification_type, &
                                   Lw_control, &
                                   microrad_properties_type, &
                                   microphysics_type, atmos_input_type,&
                                   Cldrad_control

use esfsw_parameters_mod,    only: Solar_spect, esfsw_parameters_init

use microphys_rad_mod,       only: isccp_microphys_sw_driver,   &
                                   isccp_microphys_lw_driver

!  other cloud diagnostics modules

use isccp_clouds_mod,        only: isccp_clouds_init, isccp_clouds_end,&
                                   isccp_output, isccp_cloudtypes,   &
                                   isccp_cloudtypes_stochastic


!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    cloudrad_diagnostics_mod generates any desired netcdf output
!    fields involving the cloud fields seen by the radiation package
!    or the cloud radiation interaction variables.
!
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: cloudrad_diagnostics.F90,v 20.0 2013/12/13 23:19:04 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!-------  interfaces --------

public          &
         cloudrad_diagnostics_init, cloudrad_netcdf, &
         obtain_cloud_tau_and_em, modis_yim, modis_cmip, &
         model_micro_dealloc, cloudrad_diagnostics_end

private          &
!   called from cloudrad_diagnostics_init:
         diag_field_init, &
!   called from cloudrad_netcdf:
         isccp_diag, isccp_diag_stochastic, compute_isccp_clds,  &
!   called from isccp_diag:  
         cloud_optical_properties_diag


!---------------------------------------------------------------------
!-------- namelist  ---------
!
! do_isccp                 should isccp_cloudtypes processing be done?
!
! do_outdated_isccp        should isccp_cloudtypes processing be done,
!                          here, using outdated isccp code? the 
!                          recommended approach is to use isccp 
!                          supplied via the COSP simulator
!
! isccp_actual_radprops    should the GCM's radiative properties be 
!                          used in the isccp_cloudtypes processing?
!                          If false, then use properties diagnosed
!                          locally from cloud_optical_properties_diag.
!       
! isccp_scale_factor       This scale factor is here to remove the
!                          scaling of liquid water and ice water 
!                          paths in the cloud_rad to account for the
!                          plane-parallel homogenous cloud bias.
!
!                          NOTE THAT THIS SCALE FACTOR SHOULD BE
!                          SET IDENTICAL TO THE ONE SET IN THE
!                          NAMELIST TO CLOUD_RAD.f90
!
!                          It is put here because the diagnostics
!                          are on the clouds themselves, not the
!                          radiative fluxes.  The scale factor
!                          only exists to compute radiative transfer
!                          more accurately.    

logical :: do_isccp = .false.
logical :: do_outdated_isccp = .false.
logical :: isccp_actual_radprops = .true.
real    :: isccp_scale_factor = 0.85
logical :: cloud_screen = .false.
real    :: cloud_cover_limit = 0.8
real    :: cod_limit = 2.
real    :: water_ice_ratio =1.

namelist /cloudrad_diagnostics_nml /  do_isccp, isccp_actual_radprops,&
                                      do_outdated_isccp, &
                                      isccp_scale_factor, cloud_screen,&
                                      cloud_cover_limit, cod_limit, &
                                      water_ice_ratio


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

real, parameter  :: taumin = 1.E-06  ! minimum value allowed for 
                                     ! optical depth 
                                     ! [ dimensionless ]

real,  parameter :: mid_btm  = 6.8e4  ! isccp boundaries
real, parameter  :: high_btm = 4.4e4  ! isccp boundaries
      
!----------------------------------------------------------------------
!    minimum and maximum cloud drop and crystal sizes allowable for use 
!    in microphysically-based radiative property parameterizations 
!----------------------------------------------------------------------
real             :: min_cld_drop_rad, max_cld_drop_rad, &
                    min_cld_ice_size, max_cld_ice_size, &
                    mn_drp_diam, mx_drp_diam

!-----------------------------------------------------------------------
!    if true, then donner meso clouds are treated as largescale in the 
!    optical depth diagnostic; if false, then they are included with
!    convective clouds
!-----------------------------------------------------------------------
logical          :: donner_meso_is_largescale

!----------------------------------------------------------------------
!    number of stochastic subcolumns 
!----------------------------------------------------------------------
integer          :: ncol 

integer          :: nswbands, isccpSwBand, isccpLwBand
integer          :: iuv, ivis, inir

!----------------------------------------------------------------------
!    diagnostics variables.     
!----------------------------------------------------------------------
character(len=8)    :: mod_name = 'cloudrad'
real                :: missing_value = -999.

integer :: id_tot_cld_amt, id_cld_amt, &
           id_high_cld_amt, id_mid_cld_amt, id_low_cld_amt,  &
           id_lam_cld_amt
integer :: id_reff_modis, id_reff_modis2, id_reff_modis3
integer :: id_cldtop_reff, id_cldtop_area, id_cldtop_dropnum, &
           id_dropnum_col

! radiative property diagnostics
integer :: id_em_cld_lw, id_em_cld_10u, & 
           id_abs_lsc_cld_10u, id_abs_lsc_cld_lw,  &
           id_abs_cell_cld_10u, id_abs_cell_cld_lw,  &
           id_abs_meso_cld_10u, id_abs_meso_cld_lw,  &
           id_abs_shallow_cld_10u, id_abs_shallow_cld_lw,  &
           id_abs_cld_10u, id_abs_cld_lw,  &
           id_lsc_cld_ext_uv, id_lsc_cld_ext_vis, id_lsc_cld_ext_nir, &
           id_lsc_cld_sct_uv, id_lsc_cld_sct_vis, id_lsc_cld_sct_nir, &
           id_lsc_cld_asymm_uv, id_lsc_cld_asymm_vis,    &
           id_lsc_cld_asymm_nir, &
           id_cell_cld_ext_uv, id_cell_cld_ext_vis,    &
           id_cell_cld_ext_nir, &
           id_cell_cld_sct_uv, id_cell_cld_sct_vis,    &
           id_cell_cld_sct_nir, &
           id_cell_cld_asymm_uv, id_cell_cld_asymm_vis,    &
           id_cell_cld_asymm_nir, &
           id_meso_cld_ext_uv, id_meso_cld_ext_vis,   &
           id_meso_cld_ext_nir, &
           id_meso_cld_sct_uv, id_meso_cld_sct_vis,   &
           id_meso_cld_sct_nir, &
           id_meso_cld_asymm_uv, id_meso_cld_asymm_vis,    &
           id_meso_cld_asymm_nir, &
           id_shallow_cld_ext_uv, id_shallow_cld_ext_vis,    &
           id_shallow_cld_ext_nir, &
           id_shallow_cld_sct_uv, id_shallow_cld_sct_vis,    &
           id_shallow_cld_sct_nir, &
           id_shallow_cld_asymm_uv, id_shallow_cld_asymm_vis,    &
           id_shallow_cld_asymm_nir,    &
           id_ext_cld_uv,   id_sct_cld_uv,  id_asymm_cld_uv, &
           id_ext_cld_vis,  id_sct_cld_vis, id_asymm_cld_vis, &
           id_ext_cld_nir,  id_sct_cld_nir, id_asymm_cld_nir, &
           id_alb_uv_cld, id_alb_nir_cld, id_abs_uv_cld, id_abs_nir_cld, &
           id_strat_opdepth, id_meso_opdepth, id_cell_opdepth,   &
           id_shallow_opdepth, id_largescale_opdepth, id_convect_opdepth,&
           id_total_opdepth
   
! strat cloud microphysical properties diagnostics
integer::  id_strat_area_liq, id_strat_conc_drop, id_strat_size_drop,&
           id_ra_strat_size_drop, id_strat_area_ice, &
           id_strat_conc_ice, id_strat_size_ice, &
           id_strat_droplet_number, id_gb_strat_conc_ice,   &
           id_strat_ice_number, &
!
           id_strat_ice_number_l, &
           id_strat_droplet_number_l, &
           id_strat_size_ice_l,         &
           id_strat_size_drop_l,       &
           id_drop_size_ave_l, &
           id_ice_size_ave_l, &
           id_droplet_number_ave_l, &
           id_ice_number_ave_l, &
!
           id_lsc_cld_col, &
           id_shallow_cld_col, &
           id_meso_cld_col, &
           id_cell_cld_col, &
           id_conv_cld_col,  &
           id_gb_strat_conc_drop, id_lsc_cld_amt,  id_lsc_lwp, &
           id_gb_lsc_lwp, id_lsc_iwp, id_gb_lsc_iwp

! donner meso cloud microphysical properties diagnostics
integer::  id_meso_area_liq, id_meso_conc_drop, id_meso_size_drop,&
           id_ra_meso_size_drop, id_meso_area_ice, id_meso_conc_ice, &
           id_meso_size_ice, id_meso_droplet_number, &
           id_gb_meso_conc_ice, id_gb_meso_conc_drop, id_meso_cld_amt, &
           id_meso_lwp, id_gb_meso_lwp, id_meso_iwp, id_gb_meso_iwp

! donner cell cloud microphysical properties diagnostics
integer::  id_cell_area_liq, id_cell_conc_drop, id_cell_size_drop,&
           id_ra_cell_size_drop, id_cell_area_ice, id_cell_conc_ice, &
           id_cell_size_ice, id_cell_droplet_number, &
           id_gb_cell_conc_ice, id_gb_cell_conc_drop, id_cell_cld_amt, &
           id_cell_lwp, id_gb_cell_lwp, id_cell_iwp, id_gb_cell_iwp

! uw shallow cloud microphysical properties diagnostics
integer::  id_shallow_area_liq, id_shallow_conc_drop,   &
           id_shallow_size_drop, id_ra_shallow_size_drop, &
           id_shallow_area_ice, id_shallow_conc_ice,  &
           id_shallow_size_ice, id_shallow_droplet_number, &
           id_gb_shallow_conc_ice, id_gb_shallow_conc_drop, &
           id_shallow_cld_amt, id_shallow_lwp, id_gb_shallow_lwp, &
           id_shallow_iwp, id_gb_shallow_iwp

! sum over all active cloud schemes, non-stochastic only
integer::  id_all_conc_drop, id_all_conc_ice, &
           id_predicted_cld_amt

! stochastic cloud diagnostics, avgd over all bands 
integer :: id_cldfrac_tot
integer :: id_cldfrac_ave, id_ice_conc_ave, id_drop_size_ave, &
           id_ice_size_ave, id_ra_drop_size_ave, id_drop_conc_ave, &
           id_droplet_number_ave, id_liq_col_frac_ave,  &
           id_ice_col_frac_ave, &
           id_ic_drop_conc_ave, id_ic_ice_conc_ave, id_lwp_ave,  &  
           id_ic_lwp_ave, id_iwp_ave, id_ic_iwp_ave, id_lsc_lwp_ave, &
           id_cell_lwp_ave, id_meso_lwp_ave, id_shallow_lwp_ave, &
           id_lsc_iwp_ave, id_cell_iwp_ave, id_meso_iwp_ave, &
           id_shallow_iwp_ave, id_lsc_drop_conc_ave,  &
           id_cell_drop_conc_ave, id_meso_drop_conc_ave, & 
           id_shallow_drop_conc_ave, id_lsc_ice_conc_ave,  &
           id_cell_ice_conc_ave, id_meso_ice_conc_ave, & 
           id_shallow_ice_conc_ave

!   stochastic cloud diagnostics used to show effects of extending
!   stochastic treatment to cloud types other than lsc
integer :: id_cldfrac_only_lsc, id_drop_size_only_lsc, &
           id_ice_size_only_lsc, id_ra_drop_size_only_lsc, &
           id_droplet_number_only_lsc, id_liq_col_only_lsc, &
           id_ice_col_only_lsc, id_drop_conc_only_lsc,  &
           id_ice_conc_only_lsc, id_ic_drop_conc_only_lsc, &
           id_ic_ice_conc_only_lsc, id_ic_lwp_only_lsc, &
           id_ic_iwp_only_lsc,  id_lwp_only_lsc, id_iwp_only_lsc
integer :: id_LWPr

! stochastic cloud sampling diagnostics
integer :: id_stoch_ic_cell_cf_ave, id_stoch_ic_meso_cf_ave, &
           id_stoch_ic_lsc_cf_ave, id_stoch_ic_shallow_cf_ave
integer :: id_stoch_sees_cell, id_stoch_sees_meso, &
           id_stoch_sees_lsc, id_stoch_sees_shallow
integer :: id_stoch_cell_cf_ave, id_stoch_shallow_cf_ave, &
           id_stoch_meso_cf_ave, id_stoch_lsc_cf_ave

! diagnostics for each stochastic column:
integer, dimension(:), allocatable ::    &
           id_stoch_cloud_type, &
           id_cldfrac_cols, id_ice_conc_cols, id_ice_size_cols, &
           id_drop_conc_cols, id_drop_size_cols, id_ra_drop_size_cols, &
           id_droplet_number_cols, id_lwp_cols, id_iwp_cols

! diagnostics for each stochastic column, used to show effects of 
! treating non-lsc clouds stochastically:
integer, dimension(:), allocatable ::    &
           id_cldfrac_cols_only_lsc, id_ice_conc_cols_only_lsc, &
           id_ice_size_cols_only_lsc, id_drop_conc_cols_only_lsc, &
           id_drop_size_cols_only_lsc, id_ra_drop_size_cols_only_lsc, &
           id_droplet_number_cols_only_lsc, id_lwp_cols_only_lsc, &
           id_iwp_cols_only_lsc

logical :: output_opdepth_diagnostics = .false.
logical :: module_is_initialized = .false.    ! module  initialized ?


!----------------------------------------------------------------------
!----------------------------------------------------------------------



                        contains 



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#####################################################################
! <SUBROUTINE NAME="cloudrad_diagnostics_init">
!  <OVERVIEW>
!    cloudrad_diagnostics_init is the constructor for 
!    cloudrad_diagnostics_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloudrad_diagnostics_init is the constructor for 
!    cloudrad_diagnostics_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloudrad_diagnostics_init (axes, Time)
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="real">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
! </SUBROUTINE>
!
subroutine cloudrad_diagnostics_init (min_cld_drop_rad_in,  &
                                      max_cld_drop_rad_in, &
                                      min_cld_ice_size_in, &
                                      max_cld_ice_size_in, axes, Time, &
                                      donner_meso_is_largescale_in)

!---------------------------------------------------------------------
!    cloudrad_diagnostics_init is the constructor for 
!    cloudrad_diagnostics_mod.
!------------------------------------------------------------------

real,                    intent(in)    ::   min_cld_drop_rad_in, &
                                            max_cld_drop_rad_in, &
                                            min_cld_ice_size_in, &
                                            max_cld_ice_size_in
integer, dimension(4),   intent(in)    ::   axes
type(time_type),         intent(in)    ::   Time
logical,                 intent(in)    ::   donner_meso_is_largescale_in

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       min_cld_drop_rad_in   smallest cloud droplet radius allowed by
!                             radiative properties parameterizations
!                             [ microns ]
!       max_cld_drop_rad_in   largest cloud droplet radius allowed by
!                             radiative properties parameterizations
!                             [ microns ]
!       min_cld_ice_size_in   smallest cloud ice size allowed by
!                             radiative properties parameterizations
!                             [ microns ]
!       max_cld_ice_size_in   largest cloud ice size allowed by
!                             radiative properties parameterizations
!                             [ microns ]
!       axes                  diagnostic variable axes
!       Time                  current time [time_type(days, seconds)]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer         :: unit, io, ierr, logunit

!---------------------------------------------------------------------
!   local variables:
!
!      unit     io unit for reading nml file and writing logfile
!      io       error status returned from io operation  
!      ierr     error code
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call rad_utilities_init
      call time_manager_init
      call esfsw_parameters_init
      call diag_manager_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=cloudrad_diagnostics_nml, iostat=io)
      ierr = check_nml_error(io,'cloudrad_diagnostics_nml')
#else
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=cloudrad_diagnostics_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'cloudrad_diagnostics_nml')
        enddo
10      call close_file (unit)
      endif
#endif
 
!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() )    &
                       write (logunit, nml=cloudrad_diagnostics_nml)
 
      donner_meso_is_largescale = donner_meso_is_largescale_in

!---------------------------------------------------------------------
!    define module variables to retain the smallest and largest 
!    allowable droplet and ice particle sizes which can be processed
!    by the model radiative parameterizations.
!---------------------------------------------------------------------
      min_cld_drop_rad = min_cld_drop_rad_in
      max_cld_drop_rad = max_cld_drop_rad_in
      min_cld_ice_size = min_cld_ice_size_in
      max_cld_ice_size = max_cld_ice_size_in
      mn_drp_diam    = 2.*min_cld_drop_rad
      mx_drp_diam    = 2.*max_cld_drop_rad

!---------------------------------------------------------------------
!    allocate the arrays needed to hold the diagnostics to be gener-
!    ated for each of the ncol stochastic sub-columns.
!---------------------------------------------------------------------
      ncol = Cldrad_control%nlwcldb + Solar_spect%nbands
      if (Cldrad_control%do_stochastic_clouds_iz) then
        if (Cldrad_control%do_stochastic_clouds) then
          allocate (id_stoch_cloud_type        (ncol))
          allocate (id_cldfrac_cols            (ncol))
          allocate (id_ice_conc_cols           (ncol))
          allocate (id_drop_conc_cols          (ncol))
          allocate (id_ice_size_cols           (ncol))
          allocate (id_drop_size_cols          (ncol))
          allocate (id_ra_drop_size_cols       (ncol))
          allocate (id_droplet_number_cols     (ncol))
          allocate (id_lwp_cols                (ncol))
          allocate (id_iwp_cols                (ncol))
          allocate (id_cldfrac_cols_only_lsc   (ncol))
          allocate (id_ice_conc_cols_only_lsc  (ncol))
          allocate (id_drop_conc_cols_only_lsc (ncol))
          allocate (id_ice_size_cols_only_lsc  (ncol))
          allocate (id_drop_size_cols_only_lsc      (ncol))
          allocate (id_ra_drop_size_cols_only_lsc   (ncol))
          allocate (id_droplet_number_cols_only_lsc (ncol))
          allocate (id_lwp_cols_only_lsc            (ncol))
          allocate (id_iwp_cols_only_lsc            (ncol)) 
        endif
      else
        call error_mesg ('cloudrad_diagnostics_mod',  &
         'Cldrad_control%do_stochastic_clouds not yet defined', &
                                                               FATAL)
      endif
 
!-------------------------------------------------------------------
!    initialize the netcdf diagnostics provided with this module.
!-------------------------------------------------------------------
      if (Cldrad_control%do_no_clouds_iz) then
        if (.not. Cldrad_control%do_no_clouds) then
          call diag_field_init (Time, axes)
        endif
      else
        call error_mesg ('cloudrad_diagnostics_mod',  &
         'Cldrad_control%do_no_clouds not yet defined', FATAL)
      endif

!--------------------------------------------------------------------
!    decide if isccp processing will be done here, using outdated 
!    code.
!--------------------------------------------------------------------
      if (do_isccp) then
        if (do_outdated_isccp) then
!         the outdated isccp code referenced here will be executed
        else
          call error_mesg ('cloudrad_diagnostics', &
           ' The isccp code in this module is outdated. if you REALLY &
             &want to use it, set do_outdated_isccp in this nml &
             &to .true., and resubmit; otherwise use the COSP &
             &simulator interface to obtain isccp analysis.', NOTE)
          call error_mesg ('cloudrad_diagnostics', &
             ' The yim modis output is controlled by setting &
         &do_modis_yim to .true. (default)in physics_driver_nml.', NOTE)
          do_isccp = .false.
!         if (mpp_pe() == mpp_root_pe() ) then
!           call error_mesg ('cloudrad_diagnostics', &
!            ' See above two NOTES for ways to avoid this error', FATAL)
!         endif
        endif
      else
        if (do_outdated_isccp) then
          call error_mesg ('cloudrad_diagnostics', &
           'if you REALLY want to use outdated isccp code, you must &
             &also set do_isccp in this nml to .true., and &
             &resubmit; otherwise use the COSP simulator.', NOTE)
          call error_mesg ('cloudrad_diagnostics', &
              ' To get the simple modis output, set do_modis_yim to &
               &.true. in physics_driver_nml.', NOTE)  
          if (mpp_pe() == mpp_root_pe() ) then
            call error_mesg ('cloudrad_diagnostics', &
             ' See above two NOTES for ways to avoid this error', FATAL)
          endif
        endif
      endif

!---------------------------------------------------------------------
!    initialize isccp_clouds_init 
!---------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds_iz) then
        if (Cldrad_control%do_strat_clouds) then
          if (do_isccp) call isccp_clouds_init (axes, Time)
        endif 
        if (do_isccp) then
          if (.not. Cldrad_control%do_strat_clouds) then
            call error_mesg ('cloudrad_diagnostics_mod',  &
                 'if isccp diagnostics desired, strat_clouds &
                                             &must be active', FATAL)
          endif
        endif
      else
        call error_mesg ('cloudrad_diagnostics_mod',  &
         'Cldrad_control%do_strat_clouds not yet defined', FATAL)
      endif

      nswbands = Solar_spect%nbands

!--------------------------------------------------------------------
!    define the number of shortwave bands and set integer correspond-
!    ance for diagnostics output
!
!    The understanding used in this code is that there are 2 
!    resolutions to the shortwave spectrum.  A high resolution with
!    25 bands and a low resolution with 18 bands.  The low resolution
!    is used conventional for AM2.      Here are the bands in the 
!    high and low res used for the UV, VIS, and NIR prescriptions
!    below.
!
!
!    For Low Resolution (nswbands = 18) :
!
!    Region   iband     Wavenumbers (cm-1)         Wavelength (microns)
!    ------   -----     ------------------         --------------------
!
!     UV       15          35300-36500                    0.274-0.283
!     VIS       7          16700-20000                      0.5-0.6
!     NIR       3           4200-8200                      1.22-2.38
!
!
!    For High Resolution (nswbands = 25) :
!
!    Region   iband     Wavenumbers (cm-1)         Wavelength (microns)
!    ------   -----     ------------------         --------------------
!
!     UV       22          35300-36500                    0.274-0.283
!     VIS      12          16700-20000                      0.5-0.6
!     NIR       8           6200-8200                      1.22-1.61
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Which bands to use for ISCCP cloud detection?
!
!    Note that cloud optical thickness in the visible band is sent
!    to isccp diag.  Band 6 corresponds to 14600-16700 cm-1 or 
!    0.6-0.685 microns, from 18 band structure.
!
!    If the multi-band lw emissivity is active, longwave emissivity 
!    is taken from the band closest to 10 microns (900-990 cm-1 band, 
!    10.1-11.1 microns, band 4 of 8). If the multi-band lw cloud 
!    emissivity formulation is not active, longwave emissivity is 
!    taken from band 1 (0-2200 cm-1).
!---------------------------------------------------------------------
      select case(nswbands)
        case (25) 
          isccpSwBand = 11
          iuv=22
          ivis=12
          inir=8
        case (18) 
          isccpSwBand = 6
          iuv=15
          ivis=7
          inir=3
        case default
          isccpSwBand = 6
          iuv=15
          ivis=7
          inir=3
      end select
      if (Cldrad_control%do_lw_micro_iz) then
        if (Cldrad_control%do_lw_micro) then
          isccpLwBand = 4
        else
          isccpLwBand = 1    
        end if
      else
        call error_mesg ('cloudrad_diagnostics_mod',  &
         'Cldrad_control%do_lw_micro not yet defined', FATAL)
      endif

      if (id_largescale_opdepth + id_convect_opdepth + id_strat_opdepth + &
          id_meso_opdepth + id_cell_opdepth + id_shallow_opdepth + &
          id_total_opdepth > 0) output_opdepth_diagnostics = .true.

!--------------------------------------------------------------------
!    mark the module initialized.
!--------------------------------------------------------------------
      module_is_initialized= .true.

!--------------------------------------------------------------------



end subroutine cloudrad_diagnostics_init



!--------------------------------------------------------------------

subroutine obtain_cloud_tau_and_em (is, js, Model_microphys, &
                                    Atmos_input, &
                                    Tau_stoch, Lwem_stoch)


integer,                     intent(in)     :: is, js
type(atmos_input_type),      intent(in)     :: Atmos_input
type(microphysics_type),     intent(in)     :: Model_microphys
real, dimension(:,:,:,:),    intent(inout)  :: Tau_stoch, LwEm_stoch

      integer :: n

!--------------------------------------------------------------------
!    execute the following when stochastic clouds are activated. there 
!    are separate cloud fields for each sw and lw radiative band.
!--------------------------------------------------------------------
      if (Cldrad_control%do_stochastic_clouds) then
      
!---------------------------------------------------------------------
!    after this call the Tau array is actually extinction.
!---------------------------------------------------------------------
        call isccp_microphys_sw_driver   &
                          (is, js, isccpSwBand, Model_microphys,    &
                                                 cldext=Tau_stoch) 

!---------------------------------------------------------------------
!    and to get optical thickness...
!---------------------------------------------------------------------
        do n=1,ncol
          Tau_stoch(:,:,:,n) = (Tau_stoch(:,:,:,n)*         &
                   Atmos_input%deltaz(:,:,:)/1000./isccp_scale_factor)
        end do
 
!---------------------------------------------------------------------
!    at first the LwEm array holds the absorption coefficient...
!---------------------------------------------------------------------
        call isccp_microphys_lw_driver (is, js, isccpLwBand, &
                                Model_microphys, abscoeff=LwEm_stoch)
 
!---------------------------------------------------------------------
!    and then the emissivity 
!---------------------------------------------------------------------
        do n=1,ncol
          LwEm_stoch(:,:,:,n) = 1. -   &
               exp(-1.*diffac*(LwEm_stoch(:,:,:,n)*  &
                  Atmos_input%deltaz(:,:,:)/1000.)/isccp_scale_factor)
        end do
      else
        call error_mesg ('cloudrad_diagnostics', &
              'trying to activate cosp or modis_yim without &
                                           &stochastic clouds', FATAL)
      endif ! (do_stochastic_clouds)

!-------------------------------------------------------------------

end subroutine obtain_cloud_tau_and_em 




!#####################################################################

subroutine modis_yim (is, js, Time_diag, Tau_stoch, Model_microphys, &
                      Atmos_input)

integer,                        intent(in)   :: is, js
type(time_type),                intent(in)   :: Time_diag
type(microphysics_type),        intent(in)   :: Model_microphys
type(atmos_input_type),         intent(in)   :: Atmos_input
real, dimension(:,:,:,:),       intent(in)   :: Tau_stoch           


      real, dimension(size(Atmos_input%rh2o,1),                  &
                      size(Atmos_input%rh2o,2)) :: reff_modis,   &
                                                   reff_modis2,  &
                                                   reff_modis3      
                                                  
      integer    :: i, j, n, k
      real       :: reff_n, coun_n, pres_n, Tau_m, reff_m, coun_m
      integer    :: ix, jx, kx
      real       :: min_conc
      logical    :: used

!--------------------------------------------------------------------
      ix =  size(Atmos_input%rh2o,1)
      jx =  size(Atmos_input%rh2o,2)
      kx =  size(Atmos_input%rh2o,3)

!---------------------------------------------------------------------
!     generate diagnostics related to the drop sizes which would be
!     diagnosed from MODIS satellite data. use the isccp simulator 
!     data to retrieve the drop size.
!---------------------------------------------------------------------
      if (max(id_reff_modis, id_reff_modis2, id_reff_modis3) > 0) then
        reff_modis(:,:) = 0.
        reff_modis2(:,:) = 0.
        reff_modis3(:,:) = 0.

!---------------------------------------------------------------------
!     process each model grid column. the variables ending in _n accum-
!     ulate vertical column data across the stochastic columns for a
!     given model column.
!---------------------------------------------------------------------
        do j=1,jx
          do i=1,ix
            reff_n = 0.
            coun_n = 0.
            pres_n = 0.

!---------------------------------------------------------------------
!     process each stochastic column. the variables ending in _m accu-
!     mulate data in the vertical column for a given stochastic column.
!---------------------------------------------------------------------
            do n=1,ncol
              Tau_m = 0.
              reff_m = 0.
              coun_m = 0.

!----------------------------------------------------------------------
!     scan downward in each stochastic column until the cloud optical
!     depth limit (cod_limit) is reached (the limit of the MODIS scan).
!     accumulate the effective droplet diameter for each layer the scan
!     penetrates and keep count of the number of such layers.
!----------------------------------------------------------------------
              k = 1
              do while ( k <= kx .and. Tau_m <= cod_limit)
                Tau_m = Tau_m + Tau_stoch(i,j,k,n)
                min_conc = MAX (1.0e-10, water_ice_ratio*  &
                             Model_microphys%stoch_conc_ice(i,j,k,n)) 
                if (Model_microphys%stoch_conc_drop(i,j,k,n)  &
                                                   > min_conc) then  
                  if (Model_microphys%stoch_size_drop(i,j,k,n)  &
                                                            > 1. ) then 
                    reff_m = reff_m +    &
                               Model_microphys%stoch_size_drop(i,j,k,n)
                    coun_m = coun_m + 1.
                  endif
                endif
                k = k + 1
              end do

!---------------------------------------------------------------------
!     if there were any layers in this stochastic column which are seen
!     by MODIS, add the mean droplet diameter from this column to the 
!     sum being accumulated over the stochastic columns. increment the
!     count of contributing columns, and add the pressure of maximum
!     penetration to that accumulation array. 
!---------------------------------------------------------------------
              if (coun_m >= 1.) then
                reff_n = reff_n + reff_m/coun_m
                coun_n = coun_n + 1.
                pres_n = pres_n + Atmos_input%press(i,j,k)
              endif
            end do

!---------------------------------------------------------------------
!     if there were any stochastic columns in this grid column in which
!     MODIS would have seen drops,  process the data.
!---------------------------------------------------------------------
            if (coun_n >= 1.) then

!---------------------------------------------------------------------
!     if cloud_screen is .true., then drop sizes are reported only in 
!     columns with a cloud fraction greater than cloud_cover_limit.
!     otherwise, any grid columns with cloudiness in at least one
!     stochastic column will have the drop size reported. note here that
!     droplet diameter is now converted to droplet radius, and the 
!     pressure level of maximum penetration is converted to hPa.
!---------------------------------------------------------------------
              if (cloud_screen) then
                if (coun_n/real(ncol) > cloud_cover_limit) then
                  reff_modis(i,j) = 0.5*reff_n
                  reff_modis2(i,j) = coun_n
                  reff_modis3(i,j) = pres_n*1.0e-02
                endif  
              else
                reff_modis(i,j) = 0.5*reff_n
                reff_modis2(i,j) = coun_n
                reff_modis3(i,j) = pres_n*1.0e-02
              endif
            endif
          end do
        end do

!---------------------------------------------------------------------
!     send the data to diag_manager. post-processing of these output
!     fields will be needed.
!---------------------------------------------------------------------
        used = send_data (id_reff_modis, reff_modis, Time_diag, is, js)
        used = send_data (id_reff_modis2, reff_modis2,   &
                                                     Time_diag, is, js)
        used = send_data (id_reff_modis3, reff_modis3, &
                                                     Time_diag, is, js)
      endif  ! (reff_modis)

!--------------------------------------------------------------------


end subroutine modis_yim



!#####################################################################

subroutine modis_cmip (is, js, Time_diag, Lsc_microphys, &
                      Atmos_input)

integer,                        intent(in)   :: is, js
type(time_type),                intent(in)   :: Time_diag
type(microphysics_type),        intent(in)   :: Lsc_microphys
type(atmos_input_type),         intent(in)   :: Atmos_input


      real, dimension(size(Atmos_input%rh2o,1),                  &
                      size(Atmos_input%rh2o,2)) :: cldtop_reff,  &
                                                   cldtop_dropnum,  &
                                                   cldtop_area, &      
                                                   dropnum_col
      real, dimension(size(Atmos_input%rh2o,1),                  &
                      size(Atmos_input%rh2o,2),                 &
                      size(Atmos_input%rh2o,3)) :: dpog             
                                                  
      integer    :: i, j, k
      integer    :: ix, jx, kx
      logical    :: used

!--------------------------------------------------------------------
      ix =  size(Atmos_input%rh2o,1)
      jx =  size(Atmos_input%rh2o,2)
      kx =  size(Atmos_input%rh2o,3)

      cldtop_reff = 0.0              
      cldtop_area = 0.0              
      cldtop_dropnum = 0.0              

!---------------------------------------------------------------------
!     generate diagnostics related to the drop sizes which would be
!     diagnosed from MODIS satellite data. use the isccp simulator 
!     data to retrieve the drop size.
!---------------------------------------------------------------------
      if (max(id_cldtop_reff, id_cldtop_dropnum,  &
                                    id_dropnum_col) > 0) then
!---------------------------------------------------------------------
!     process each model grid column. the variables ending in _n accum-
!     ulate vertical column data across the stochastic columns for a
!     given model column.
!---------------------------------------------------------------------
        do j=1,jx
          do i=1,ix
            do k=1, kx

!----------------------------------------------------------------------
!     scan downward in each stochastic column until the cloud optical
!     depth limit (cod_limit) is reached (the limit of the MODIS scan).
!     accumulate the effective droplet diameter for each layer the scan
!     penetrates and keep count of the number of such layers.
!----------------------------------------------------------------------
                if (Lsc_microphys%conc_drop(i,j,k) > 0.0) then
                  cldtop_reff  (i,j) =  0.5*  &
                                Lsc_microphys%cldamt(i,j,k)* &
                                         Lsc_microphys%size_drop(i,j,k)
                  cldtop_dropnum(i,j) = (Atmos_input%press(i,j,k)/  &
                                   (RDGAS*Atmos_input%temp(i,j,k)))* &
                                Lsc_microphys%cldamt(i,j,k)* &
                                    Lsc_microphys%droplet_number(i,j,k)
                  cldtop_area(i,j) =  & 
                            Lsc_microphys%cldamt(i,j,k)
                  exit
                endif
            end do
          end do
          end do

          do k=1, kx
            dpog(:,:,k) = (Atmos_input%pflux(:,:,k+1) -     &
                                     Atmos_input%pflux(:,:,k))/GRAV
          end do
          dropnum_col(:,:) = SUM  &
            (Lsc_microphys%droplet_number*dpog*Lsc_microphys%cldamt, &
                                                                 dim=3)

!---------------------------------------------------------------------
!     send the data to diag_manager. post-processing of these output
!     fields will be needed.
!---------------------------------------------------------------------
        used = send_data (id_cldtop_reff, 1.0e-06*cldtop_reff, Time_diag, is, js)
        used = send_data (id_cldtop_area, cldtop_area, Time_diag, is, js)
        used = send_data (id_cldtop_dropnum, cldtop_dropnum,   &
                                                     Time_diag, is, js)
        used = send_data (id_dropnum_col, dropnum_col, &
                                                     Time_diag, is, js)
      endif  ! (id_cldtop_reff)

!----------------------------------------------------------------------


end subroutine modis_cmip



!###################################################################
! <SUBROUTINE NAME="cloudrad_netcdf">
!  <OVERVIEW>
!    cloudrad_netcdf generates and outputs netcdf fields describing the
!    cloud radiative properties, along with isccp cloud diagnostics
!    fields. 
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloudrad_netcdf generates and outputs netcdf fields describing the
!    cloud radiative properties, along with isccp cloud diagnostics
!    fields. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloudrad_netcdf (is, js, Time, Time_diag, Atmos_input, cosz, &
!                            Lsc_microphys, Meso_microphys, &
!                            Cell_microphys, Shallow_microphys, &
!                            Lscrad_props,  Mesorad_props, &
!                            Cellrad_props, Shallowrad_props, &
!                            Cldrad_props,&
!                            Cld_spec, mask)
!  </TEMPLATE>
!  <IN NAME="is,js" TYPE="integer">
!   starting subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
!  <IN NAME="Time_diag" TYPE="time_type">
!   time on next timestep, used as stamp for 
!                        diagnostic output [ time_type (days, seconds) ]
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    atmospheric input fields on model grid,
!  </IN>
!  <IN NAME="cosz" TYPE="real">
!    cosine of solar zenith angle
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </IN>
!  <IN NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </IN>
!  <IN NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale 
!                        clouds assciated with donner convection
!  </IN>
!  <IN NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                        clouds associated with donner convection
!  </IN>
!  <IN NAME="Shallow_microphys" TYPE="microphysics_type">
!   microphysical specification for 
!                        clouds associated with uw shallow convection
!  </IN>
!  <IN NAME="Lscrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the large-scale 
!                      clouds   
!  </IN>
!  <IN NAME="Mesorad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the meso-scale
!                      clouds assciated with donner convection
!  </IN>
!  <IN NAME="Cellrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the convective cell
!                      clouds associated with donner convection 
!  </IN>
!  <IN NAME="Shallowrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the 
!                      clouds associated with uw shallow convection 
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   cloud radiative properties on model grid
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
! </SUBROUTINE>
!
subroutine cloudrad_netcdf (is, js, Time, Time_diag, Atmos_input, cosz, &
                            Lsc_microphys, Meso_microphys, &
                            Cell_microphys, Shallow_microphys, &
                            Lscrad_props,   &
                            Mesorad_props, Cellrad_props,  &
                            Shallowrad_props, Cldrad_props,&
                            Cld_spec, Model_microphys, mask)

!---------------------------------------------------------------------
!    cloudrad_netcdf generates and outputs netcdf fields describing the
!    cloud radiative properties, along with isccp cloud diagnostics
!    fields. 
!---------------------------------------------------------------------

integer,                        intent(in)      :: is, js
type(time_type),                intent(in)      :: Time, Time_diag
type(atmos_input_type),         intent(in)      :: Atmos_input
real, dimension(:,:),           intent(in)      :: cosz        
type(microphysics_type),        intent(in)      :: Lsc_microphys, &
                                                   Meso_microphys,&
                                                   Cell_microphys,&
                                                   Shallow_microphys
type(microrad_properties_type), intent(in)      :: Lscrad_props, &
                                                   Mesorad_props, &
                                                   Cellrad_props, &
                                                   Shallowrad_props
type(cldrad_properties_type),   intent(in)      :: Cldrad_props
type(cld_specification_type),   intent(in)      :: Cld_spec       
type(microphysics_type),        intent(inout)   :: Model_microphys
real, dimension(:,:,:),         intent(in),  &
                                       optional :: mask

!-------------------------------------------------------------------
!   intent(in) variables:
!
!      is,js           starting subdomain i,j indices of data 
!                      in the physics_window being integrated
!      Time            current time [time_type(days, seconds)]
!      Time_diag       time on next timestep, used as stamp for 
!                      diagnostic output [ time_type (days, seconds) ]
!      Atmos_input     atmospheric input fields on model grid,
!                      [ atmos_input_type ]
!      cosz            cosine of zenith angle --  mean value over
!                      appropriate averaging interval
!                      [ non-dimensional ]
!      Lsc_microphys   microphysical specification for large-scale 
!                      clouds
!                      [ microphysics_type ]
!      Meso_microphys  microphysical specification for meso-scale 
!                      clouds assciated with donner convection
!                      [ microphysics_type ]
!      Cell_microphys  microphysical specification for convective cell
!                      clouds associated with donner convection
!                      [ microphysics_type ]
!   Shallow_microphys  microphysical specification for 
!                      clouds associated with uw shallow convection
!                      [ microphysics_type ]
!      Lscrad_props    cloud radiative properties for the large-scale 
!                      clouds   
!                      [ microrad_properties_type ]
!      Mesorad_props   cloud radiative properties for meso-scale 
!                      clouds associated with donner convection   
!                      [ microrad_properties_type ]
!      Cellrad_props   cloud radiative properties for convective cell
!                      clouds associated with donner convection  
!                      [ microrad_properties_type ]
!    Shallowrad_props   
!                      cloud radiative properties for
!                      clouds associated with uw shallow convection  
!                      [ microrad_properties_type ]
!      Cldrad_props    total-cloud radiative properties,
!                      [ cldrad_properties_type ]
!      Cld_spec        variables on the model grid which define the 
!                      cloud location and amount     
!                      [ cld_specification_type ]
!
!   intent(in), optional variables:
!
!      mask              present when running eta vertical coordinate,
!                        mask to remove points below ground
!
!------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables:

      real, dimension(size(Atmos_input%rh2o,1),                       &
                      size(Atmos_input%rh2o,2),                       &
                      size(Atmos_input%rh2o,3))  ::    &
                    cloud, cloud2, dpog, pmass, pmass2, ptm2

      logical, dimension(size(Atmos_input%rh2o,1),                   &
                         size(Atmos_input%rh2o,2))                    &
                                                    :: tmplmask2

      logical, dimension(size(Atmos_input%rh2o,1),                    &
                         size(Atmos_input%rh2o,2),                  &
                         size(Atmos_input%rh2o,3))  :: tmplmask, &
                                                       tmplmaska, &
                                                       tmplmaskl

      logical, dimension(size(Atmos_input%rh2o,1),                    &
                         size(Atmos_input%rh2o,2),                  &
                         size(Atmos_input%rh2o,3) ,     &
         Cldrad_control%nlwcldb + Solar_spect%nbands) ::   &
                                                    tmplmask4

      real, dimension(size(Atmos_input%rh2o,1),                       &
                      size(Atmos_input%rh2o,2))   :: tca, cloud2d,  &
                                                     tca2, tca3 

      real, dimension(size(Atmos_input%rh2o,1),                       &
                      size(Atmos_input%rh2o,2),    &
                    Cldrad_control%nlwcldb + Solar_spect%nbands) ::   &
                                             cloud2n     

      real, dimension(size(Atmos_input%rh2o,1),                       &
                      size(Atmos_input%rh2o,2), 4)  :: hml_ca        

      real, dimension(size(Atmos_input%rh2o,1),                       &
                      size(Atmos_input%rh2o,2),                       &
                      size(Atmos_input%rh2o,3),     &
                    Cldrad_control%nlwcldb + Solar_spect%nbands) ::   &
                                                 Tau_stoch, LwEm_stoch

      real, dimension(size(Atmos_input%rh2o,1),                       &
                      size(Atmos_input%rh2o,2),                       &
                      size(Atmos_input%rh2o,3))   ::  Tau, LwEm,  &
                                                     tau_c, tau_s, &
                                          tau_strat, tau_meso, tau_cell, &
                                          tau_uw, tau_tot

      logical    :: used
      integer    :: ix, jx, kx
      integer    :: i, j, k, n
      integer    :: nn
      integer    :: ctr_s, ctr_c, ctr_strat, ctr_meso, ctr_cell, ctr_uw, &
                    ctr_tot
      real       :: sum_s1, sum_c1, sum_strat, sum_meso, sum_cell, sum_uw, &
                    sum_tot


      
!---------------------------------------------------------------------
!  local variables:
!
!      cloud                array used to hold the various netcdf 
!                           output fields as they are sent off to 
!                           diag_manager_mod
!      tca                  total column cloud amount [ dimensionless ]
!      hml_ca               total column cloud amount for isccp high, 
!                           middle and low clouds, individually
!      used                 flag returned from send_data indicating
!                           whether diag_manager_mod has received 
!                           data that was sent
!      kx                   number of model layers 
!      i,j,k                do-loop indices
!      Model_microphys      microphysics_type variable used to hold the
!                           cloud physical properties actuaaly seen by
!                           the model in each stochastic band (only 
!                           present when do_stochastic_clouds = .true.)
!
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('cloudrad_diagnostics_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

if (Time_diag > Time) then

!--------------------------------------------------------------------
!    define the array dimensions on the processor.
!---------------------------------------------------------------------
      ix =  size(Cld_spec%camtsw,1)
      jx =  size(Cld_spec%camtsw,2)
      kx =  size(Cld_spec%camtsw,3)

!----------------------------------------------------------------------
!    compute the depth of each model layer to be used in defining water
!    paths. include the 10**3 factor needed to produce water path units
!    of kg(h2o) / m**2 rather than  g(h2o) / m**2. 
!    pmass is consistent with the expression used for donner ans uw
!    shallow clouds, while pmass2 is consistent with that used in
!    strat_cloud.
!----------------------------------------------------------------------
      do k=1,kx
        dpog(:,:,k) = (Atmos_input%pflux(:,:,k+1) -     &
                                     Atmos_input%pflux(:,:,k))/GRAV
        ptm2(:,:,k) = (Atmos_input%pflux(:,:,k+1) -    &
                                   Atmos_input%pflux(:,:,k)) / &
                         log(Atmos_input%pflux(:,:,k+1)/  &
                                   MAX(Atmos_input%pflux(:,:,k),   &
                                        Atmos_input%press(:,:,1))) 
        pmass2(:,:,k) = dpog(:,:,k)/  &
                         (1.0e03*ptm2(:,:,k)*isccp_scale_factor/  &
                                 (RDGAS*Atmos_input%temp(:,:,k)))
        pmass(:,:,k) = dpog(:,:,k)/   &
                         (1.0e03*Atmos_input%press(:,:,k)/&
                                 (RDGAS*Atmos_input%temp(:,:,k)))
      end do

!---------------------------------------------------------------------
!    allocate and initialize the components of a microphysics_type 
!    variable Model_microphys. this variable is used to hold the 
!    combination of stochastic column cloud physical properties actually
!    used by the model when stochastic clouds are active, and the 
!    combined cloud properties of all active cloud schemes when 
!    stochastic clouds are not active.
!---------------------------------------------------------------------
      allocate (Model_microphys%conc_drop  (ix, jx, kx) )
      allocate (Model_microphys%conc_ice   (ix, jx, kx) )
      allocate (Model_microphys%conc_rain  (ix, jx, kx) )
      allocate (Model_microphys%conc_snow  (ix, jx, kx) )
      allocate (Model_microphys%size_drop  (ix, jx, kx) )
      allocate (Model_microphys%size_ice   (ix, jx, kx) )
      allocate (Model_microphys%size_rain  (ix, jx, kx) )
      allocate (Model_microphys%size_snow  (ix, jx, kx) )
      allocate (Model_microphys%cldamt     (ix, jx, kx) )
      allocate (Model_microphys%droplet_number  (ix, jx, kx) )
      allocate (Model_microphys%ice_number (ix, jx, kx) )
      Model_microphys%conc_drop = 0.
      Model_microphys%conc_ice  = 0.
      Model_microphys%conc_rain = 0.
      Model_microphys%conc_snow = 0.
      Model_microphys%size_drop = 1.0e-20
      Model_microphys%size_ice  = 1.0e-20
      Model_microphys%size_rain = 1.0e-20
      Model_microphys%size_snow = 1.0e-20
      Model_microphys%cldamt     = 0.
      Model_microphys%droplet_number = 0.              
      Model_microphys%ice_number = 0.

      if (Cldrad_control%do_stochastic_clouds) then
        allocate (Model_microphys%stoch_conc_ice (ix, jx, kx, nCol) )
        allocate (Model_microphys%stoch_conc_drop(ix, jx, kx, nCol) )
        allocate (Model_microphys%stoch_size_ice (ix, jx, kx, nCol) )
        allocate (Model_microphys%stoch_size_drop(ix, jx, kx, nCol) )
        allocate (Model_microphys%stoch_cldamt   (ix, jx, kx, nCol) )
        allocate (Model_microphys%stoch_cloud_type (ix, jx, kx, nCol) )
        allocate (Model_microphys%stoch_droplet_number   &
                                                 (ix, jx, kx, nCol) )
        allocate (Model_microphys%stoch_ice_number   &
                                                 (ix, jx, kx, nCol) )

        Model_microphys%lw_stoch_conc_drop =>    &
         Model_microphys%stoch_conc_drop(:,:,:,1:Cldrad_control%nlwcldb)
        Model_microphys%lw_stoch_conc_ice  =>    &
         Model_microphys%stoch_conc_ice (:,:,:,1:Cldrad_control%nlwcldb)
        Model_microphys%lw_stoch_size_drop =>    &
         Model_microphys%stoch_size_drop(:,:,:,1:Cldrad_control%nlwcldb)
        Model_microphys%lw_stoch_size_ice  =>    &
         Model_microphys%stoch_size_ice (:,:,:,1:Cldrad_control%nlwcldb)
        Model_microphys%sw_stoch_conc_drop =>    &
         Model_microphys%stoch_conc_drop   &
                                       (:,:,:,Cldrad_control%nlwcldb+1:)
        Model_microphys%sw_stoch_conc_ice  =>    &
         Model_microphys%stoch_conc_ice(:,:,:,Cldrad_control%nlwcldb+1:)
        Model_microphys%sw_stoch_size_drop =>    &
         Model_microphys%stoch_size_drop   &
                                       (:,:,:,Cldrad_control%nlwcldb+1:)
        Model_microphys%sw_stoch_size_ice  =>    &
         Model_microphys%stoch_size_ice(:,:,:,Cldrad_control%nlwcldb+1:)
        Model_microphys%lw_stoch_cldamt =>     &
         Model_microphys%stoch_cldamt(:,:,:,1:Cldrad_control%nlwcldb)
        Model_microphys%sw_stoch_cldamt =>     &
         Model_microphys%stoch_cldamt(:,:,:,Cldrad_control%nlwcldb+1:)
        Model_microphys%lw_stoch_droplet_number =>     &
         Model_microphys%stoch_droplet_number   &
                                        (:,:,:,1:Cldrad_control%nlwcldb)
        Model_microphys%sw_stoch_droplet_number =>     &
         Model_microphys%stoch_droplet_number  &
                                       (:,:,:,Cldrad_control%nlwcldb+1:)
        Model_microphys%stoch_conc_drop = 0.
        Model_microphys%stoch_conc_ice  = 0.
        Model_microphys%stoch_size_drop = 1.0e-20
        Model_microphys%stoch_size_ice  = 1.0e-20
        Model_microphys%stoch_cldamt = 0.0
        Model_microphys%stoch_droplet_number = 0.0
        Model_microphys%stoch_ice_number = 0.0

!---------------------------------------------------------------------
!    since the sw bands come first in Cld_spec and the lw bands come 
!    first in  Model_microphys, the band index order must be reversed 
!    in defining Model_microphys%stoch_cloud_type.  
!---------------------------------------------------------------------
        do n=1,ncol
          if ( n > Solar_spect%nbands) then
            nn = n - Solar_spect%nbands
          else
            nn = n + Cldrad_control%nlwcldb
          endif

          Model_microphys%stoch_cloud_type(:,:,:,nn)  =  &
                                     Cld_spec%stoch_cloud_type(:,:,:,n)
        end do

!---------------------------------------------------------------------
!    define the cloud properties assigned to each stochastic column,
!    based on the cloud type assignment contained in 
!    Model_microphys%stoch_cloud_type. 
!---------------------------------------------------------------------
        if (Cldrad_control%do_donner_deep_clouds .or.  &
                                   Cldrad_control%do_uw_clouds) then
          do n=1,ncol
            do k=1,kx
              do j=1,jx
                do i=1,ix
                  if (Model_microphys%stoch_cloud_type(i,j,k,n)    &
                                                             == 1) then
                    Model_microphys%stoch_conc_drop(i,j,k,n) =    &
                                 Lsc_microphys%stoch_conc_drop(i,j,k,n)
                    Model_microphys%stoch_conc_ice(i,j,k,n)  =    &
                                Lsc_microphys%stoch_conc_ice (i,j,k,n)
                    Model_microphys%stoch_size_drop(i,j,k,n) =    &
                                Lsc_microphys%stoch_size_drop(i,j,k,n)
                    Model_microphys%stoch_size_ice(i,j,k,n)  =    &
                                 Lsc_microphys%stoch_size_ice(i,j,k,n)
                    Model_microphys%stoch_cldamt(i,j,k,n) =    1.0 
                    Model_microphys%stoch_droplet_number(i,j,k,n) =   &
                            Lsc_microphys%stoch_droplet_number(i,j,k,n)
                    Model_microphys%stoch_ice_number(i,j,k,n) =   &
                            Lsc_microphys%stoch_ice_number(i,j,k,n)
                  else if (Model_microphys%stoch_cloud_type(i,j,k,n) &
                                                             == 2) then
                    Model_microphys%stoch_conc_drop(i,j,k,n) =    &
                                Meso_microphys%conc_drop(i,j,k)
                    Model_microphys%stoch_conc_ice(i,j,k,n)  =    &
                                Meso_microphys%conc_ice (i,j,k)
                    Model_microphys%stoch_size_drop(i,j,k,n) =    &
                                Meso_microphys%size_drop(i,j,k)
                    Model_microphys%stoch_size_ice(i,j,k,n)  =    &
                                Meso_microphys%size_ice(i,j,k)
                    Model_microphys%stoch_cldamt(i,j,k,n) =   1.0
                    Model_microphys%stoch_droplet_number(i,j,k,n) =   &
                           Meso_microphys%droplet_number(i,j,k) 
!cms++
!!                    Model_microphys%stoch_ice_number(i,j,k,n) =   &
!!                           Meso_microphys%ice_number(i,j,k) 
!cms--
                  else if (Model_microphys%stoch_cloud_type(i,j,k,n)   &
                                                             == 3) then
                    Model_microphys%stoch_conc_drop(i,j,k,n) =    &
                                 Cell_microphys%conc_drop(i,j,k)
                    Model_microphys%stoch_conc_ice(i,j,k,n)  =    &
                                 Cell_microphys%conc_ice (i,j,k)
                    Model_microphys%stoch_size_drop(i,j,k,n) =    &
                                 Cell_microphys%size_drop(i,j,k)
                    Model_microphys%stoch_size_ice(i,j,k,n)  =    &
                                  Cell_microphys%size_ice(i,j,k)
                    Model_microphys%stoch_cldamt(i,j,k,n) =  1.0    
                    Model_microphys%stoch_droplet_number(i,j,k,n) =   &
                           Cell_microphys%droplet_number(i,j,k) 
!cms++
!!                    Model_microphys%stoch_ice_number(i,j,k,n) =   &
!!                           Cell_microphys%ice_number(i,j,k)
!cms--
                  else if (Model_microphys%stoch_cloud_type(i,j,k,n)  &
                                                             == 4) then
                    Model_microphys%stoch_conc_drop(i,j,k,n) =    &
                              Shallow_microphys%conc_drop(i,j,k)
                    Model_microphys%stoch_conc_ice(i,j,k,n)  =    &
                              Shallow_microphys%conc_ice (i,j,k)
                    Model_microphys%stoch_size_drop(i,j,k,n) =    &
                              Shallow_microphys%size_drop(i,j,k)
                    Model_microphys%stoch_size_ice(i,j,k,n)  =    &
                              Shallow_microphys%size_ice(i,j,k)
                    Model_microphys%stoch_cldamt(i,j,k,n) =  1.0
                    Model_microphys%stoch_droplet_number(i,j,k,n) =   &
                              Shallow_microphys%droplet_number(i,j,k) 
!cms++
!!                    Model_microphys%stoch_ice_number(i,j,k,n) =   &
!!                         Shallow_microphys%ice_number(i,j,k)
!cms--
                  endif
                end do
              end do
            end do
          end do
        else  ! (do_donner_deep_clouds .or. do_uw_clouds)

!---------------------------------------------------------------------
!    if only strat cloud is active, define all column data to be that
!    coming from the strat_cloud stochasticization.
!---------------------------------------------------------------------
          Model_microphys%stoch_conc_drop =    &
                                  Lsc_microphys%stoch_conc_drop
          Model_microphys%stoch_conc_ice  =    &
                          Lsc_microphys%stoch_conc_ice
          Model_microphys%stoch_size_drop =    &
                                Lsc_microphys%stoch_size_drop
          Model_microphys%stoch_size_ice  =    & 
                                Lsc_microphys%stoch_size_ice
          Model_microphys%stoch_cldamt =       &
                                    Lsc_microphys%stoch_cldamt 
          Model_microphys%stoch_droplet_number =     &
                            Lsc_microphys%stoch_droplet_number 
          Model_microphys%stoch_ice_number =     &
                            Lsc_microphys%stoch_ice_number
        endif ! (do_donner_deep_clouds)

!---------------------------------------------------------------------
!    if stochastic clouds are not active, use Model_microphys to contain
!    the total cloud field obtained by summing the contribuutions from 
!    each active cloud scheme.
!---------------------------------------------------------------------
      else  ! (do_stochastic_clouds)

!----------------------------------------------------------------------
!    define the total cloud fraction.
!----------------------------------------------------------------------
        Model_microphys%cldamt = Lsc_microphys%cldamt
        if (Cldrad_control%do_donner_deep_clouds) then 
          Model_microphys%cldamt = Model_microphys%cldamt + &
                          Meso_microphys%cldamt + Cell_microphys%cldamt
        endif
        if (Cldrad_control%do_uw_clouds) then
          Model_microphys%cldamt = Model_microphys%cldamt + &
                                                Shallow_microphys%cldamt
        endif

!--------------------------------------------------------------------
!    define total grid box values of drop and ice cloud amounts, 
!    appropriately weighted by fractional area, and summed over all 
!    active cloud types. 
!--------------------------------------------------------------------
        Model_microphys%conc_drop = Lsc_microphys%cldamt*  &
                                              Lsc_microphys%conc_drop
        Model_microphys%conc_ice  = Lsc_microphys%cldamt*  &
                                              Lsc_microphys%conc_ice
         
        if (Cldrad_control%do_donner_deep_clouds) then 
          Model_microphys%conc_drop = Model_microphys%conc_drop +  &
                 Meso_microphys%cldamt*Meso_microphys%conc_drop +   &
                 Cell_microphys%cldamt*Cell_microphys%conc_drop
          Model_microphys%conc_ice = Model_microphys%conc_ice +  &
                 Meso_microphys%cldamt*Meso_microphys%conc_ice +   &
                 Cell_microphys%cldamt*Cell_microphys%conc_ice
        endif      
        if (Cldrad_control%do_uw_clouds) then
          Model_microphys%conc_drop = Model_microphys%conc_drop +  &
                   Shallow_microphys%cldamt*Shallow_microphys%conc_drop
          Model_microphys%conc_ice = Model_microphys%conc_ice +  &
                   Shallow_microphys%cldamt*Shallow_microphys%conc_ice
        endif

!----------------------------------------------------------------------
!    adjust the total cld fraction to be no larger than 1.0.
!----------------------------------------------------------------------
        Model_microphys%cldamt = MIN(Model_microphys%cldamt, 1.0)

      endif ! (do_stochastic_clouds)


!---------------------------------------------------------------------
!
!
!
!                   ISCCP SIMULATOR SECTION
!
!
!
!
!---------------------------------------------------------------------
 
!--------------------------------------------------------------------
!    if desired, call isccp_diag to generate isccp-relevant diagnostics
!    when running strat_clouds.
!---------------------------------------------------------------------

      if (do_isccp .or. output_opdepth_diagnostics) then  
        call obtain_cloud_tau_and_em (is, js, Model_microphys, &
                                    Atmos_input, Tau_stoch, Lwem_stoch)
        if (output_opdepth_diagnostics) then
!---------------------------------------------------------------------
!   the values of tau and lwem are available for each stochastic column.
!   here grid box mean values are obtained for the convective and
!   large-scale components. 
!---------------------------------------------------------------------
          do k=1, size(Tau_stoch,3)
            do j=1, size(Tau_stoch,2)
              do i=1, size(Tau_stoch,1)
                ctr_s = 0
                ctr_c = 0
                ctr_strat = 0
                ctr_meso = 0
                ctr_cell = 0
                ctr_uw = 0
                ctr_tot = 0
                sum_s1 = 0.
                sum_c1 = 0.
                sum_strat = 0.
                sum_meso = 0.
                sum_cell = 0.
                sum_uw = 0.
                sum_tot = 0.
                do n=1, size(Tau_stoch,4)
                  if (Model_microphys%stoch_cloud_type(i,j,k,n) == 1.) then
                    ctr_s = ctr_s + 1
                    sum_s1 = sum_s1 +  tau_stoch(i,j,k,n)
                    ctr_strat = ctr_strat + 1
                    sum_strat = sum_strat + tau_stoch(i,j,k,n)
                    ctr_tot = ctr_tot + 1
                    sum_tot = sum_tot + tau_stoch(i,j,k,n)
                  else if   &
                    (Model_microphys%stoch_cloud_type(i,j,k,n) == 2. ) then
                    ctr_meso = ctr_meso + 1
                    sum_meso = sum_meso + tau_stoch(i,j,k,n)
                    if (donner_meso_is_largescale) then
                      ctr_s = ctr_s + 1
                      sum_s1 = sum_s1 +  tau_stoch(i,j,k,n)
                    else
                      ctr_c = ctr_c + 1
                      sum_c1 = sum_c1 +  tau_stoch(i,j,k,n)
                    endif
                    ctr_tot = ctr_tot + 1
                    sum_tot = sum_tot + tau_stoch(i,j,k,n)
                  else if    &
                    (Model_microphys%stoch_cloud_type(i,j,k,n) == 3. ) then
                    ctr_cell = ctr_cell + 1
                    sum_cell = sum_cell + tau_stoch(i,j,k,n)
                    ctr_c = ctr_c + 1
                    sum_c1 = sum_c1 +  tau_stoch(i,j,k,n)
                    ctr_tot = ctr_tot + 1
                    sum_tot = sum_tot + tau_stoch(i,j,k,n)
                  else if    &
                    (Model_microphys%stoch_cloud_type(i,j,k,n) == 4. ) then
                    ctr_uw = ctr_uw + 1
                    sum_uw = sum_uw + tau_stoch(i,j,k,n)
                    ctr_c = ctr_c + 1
                    sum_c1 = sum_c1 +  tau_stoch(i,j,k,n)
                    ctr_tot = ctr_tot + 1
                    sum_tot = sum_tot + tau_stoch(i,j,k,n)
                  endif
                end do
                if (ctr_s > 0) then
                  tau_s(i,j,k) = sum_s1/ctr_s
                else
                  tau_s(i,j,k) = 0.             
                endif
                if (ctr_c > 0) then
                  tau_c(i,j,k) = sum_c1/ctr_c
                else
                  tau_c(i,j,k) = 0.             
                endif
                if (ctr_strat > 0) then
                  tau_strat(i,j,k) = sum_strat/ctr_strat
                else
                  tau_strat(i,j,k) = 0.             
                endif
                if (ctr_meso > 0) then
                  tau_meso(i,j,k) = sum_meso/ctr_meso
                else
                  tau_meso(i,j,k) = 0.             
                endif
                if (ctr_cell > 0) then
                  tau_cell(i,j,k) = sum_cell/ctr_cell
                else
                  tau_cell(i,j,k) = 0.             
                endif
                if (ctr_uw > 0) then
                  tau_uw(i,j,k) = sum_uw/ctr_uw
                else
                  tau_uw(i,j,k) = 0.             
                endif
                if (ctr_tot > 0) then
                  tau_tot(i,j,k) = sum_tot/ctr_tot
                else
                  tau_tot(i,j,k) = 0.             
                endif
              end do
            end do
          end do
          used = send_data (id_largescale_opdepth, tau_s(:,:,:), &
                            Time_diag, is, js, 1, rmask=mask)
          used = send_data (id_convect_opdepth, tau_c(:,:,:), &
                            Time_diag, is, js, 1, rmask=mask)
          used = send_data (id_strat_opdepth, tau_strat(:,:,:), &
                            Time_diag, is, js, 1, rmask=mask)
          used = send_data (id_meso_opdepth, tau_meso(:,:,:), &
                            Time_diag, is, js, 1, rmask=mask)
          used = send_data (id_cell_opdepth, tau_cell(:,:,:), &
                            Time_diag, is, js, 1, rmask=mask)
          used = send_data (id_shallow_opdepth, tau_uw(:,:,:), &
                            Time_diag, is, js, 1, rmask=mask)
          used = send_data (id_total_opdepth, tau_tot(:,:,:), &
                            Time_diag, is, js, 1, rmask=mask)
        endif
      endif  ! (do_isccp )

!--------------------------------------------------------------------
!    execute the following when stochastic clouds are activated. there 
!    are separate cloud fields for each sw and lw radiative band.
!--------------------------------------------------------------------
      if (Cldrad_control%do_stochastic_clouds) then

        if (do_isccp) then
                        
!----------------------------------------------------------------------
!    call isccp_diag_stochastic to map the stochastic clouds and cloud
!    properties to the isccp cloud categories.
!----------------------------------------------------------------------
          call isccp_diag_stochastic (is, js, Atmos_input, cosz, &
                                        Tau_stoch, LwEm_stoch,   &
                                        Model_microphys%stoch_cldamt, &
                                        Time_diag)
        endif

!--------------------------------------------------------------------
!    define the isccp properties when stochastic clouds are not active.
!    here there is only a single cloud profile for each gridbox.
!---------------------------------------------------------------------
      else  ! (do_stochastic_clouds)
        if ( do_isccp ) then 
          Tau(:,:,:) = (Lscrad_props%cldext(:,:,:,isccpSwBand)* &
                                  Atmos_input%deltaz(:,:,:)/1000.) / &
                                                     isccp_scale_factor
          LwEm(:,:,:) =  1. - exp( -1. * diffac *        &
                           (Lscrad_props%abscoeff(:,:,:,isccpLwBand)* &
                                    Atmos_input%deltaz(:,:,:)/1000.)/ &
                                                    isccp_scale_factor) 
          call isccp_diag (is, js, Cld_spec, Atmos_input, cosz, &
                                                  Tau, LwEm, Time_diag)
        endif
      endif ! (do_stochastic_clouds)

!---------------------------------------------------------------------
!
!
!
!      COMPUTE HIGH, MIDDLE, AND LOW CLOUD AMOUNTS
!
!
!---------------------------------------------------------------------
 
!---------------------------------------------------------------------
!    when stochastic clouds are active:
!---------------------------------------------------------------------

      if (Cldrad_control%do_stochastic_clouds) then

!---------------------------------------------------------------------
!    define the total cloud amount as the percentage of stochastic 
!    columns containing cloud at any model level. 
!---------------------------------------------------------------------
!       if (id_tot_cld_amt > 0  .or. id_cldfrac_tot > 0) then
        if (max(id_tot_cld_amt, id_cldfrac_tot) > 0) then
          cloud2n(:,:,:) =    &
                    SUM (Model_microphys%stoch_cldamt(:,:,:,:), dim = 3)
          cloud2n(:,:,:) = MIN (cloud2n(:,:,:), 1.0)
          tca2(:,:) = 100.0*SUM (cloud2n(:,:,:), dim = 3)/REAL (ncol)
          used = send_data (id_tot_cld_amt, tca2, Time_diag, is, js)
        endif 

!cms++
!---------------------------------------------------------------------
! fraction of stochastic columns containing lsc clouds at any model level 
!---------------------------------------------------------------------
        if ( id_lsc_cld_col > 0) then

          cloud(:,:,:) = REAL (COUNT(    &
                  Model_microphys%stoch_cloud_type(:,:,:,:) == 1,    &
                                                   dim = 4))/REAL(ncol)

          tca3(:,:) =  SUM (cloud, dim = 3)
          tca3(:,:) =  MIN (tca3(:,:), 1.0)


          used = send_data (id_lsc_cld_col, tca3, Time_diag, is, js)

        endif 

!---------------------------------------------------------------------
! stochastic columns containing anvils assuming maximum overlap
!---------------------------------------------------------------------

        if ( id_meso_cld_col > 0) then

          cloud(:,:,:) = REAL (COUNT(    &
                  Model_microphys%stoch_cloud_type(:,:,:,:) == 2,    &
                                                   dim = 4))/REAL(ncol)

          tca3(:,:) =  MAXVAL (cloud, dim = 3)
          


          used = send_data (id_meso_cld_col, tca3, Time_diag, is, js)

        endif 



        if ( id_cell_cld_col > 0) then

          cloud(:,:,:) = REAL (COUNT(    &
                  Model_microphys%stoch_cloud_type(:,:,:,:) == 3,    &
                                                   dim = 4))/REAL(ncol)

          tca3(:,:) =  MAXVAL (cloud, dim = 3)
          

          used = send_data (id_cell_cld_col, tca3, Time_diag, is, js)

        endif 


        if ( id_shallow_cld_col > 0) then

          cloud(:,:,:) = REAL (COUNT(    &
                  Model_microphys%stoch_cloud_type(:,:,:,:) == 4,    &
                                                   dim = 4))/REAL(ncol)

          tca3(:,:) =  MAXVAL (cloud, dim = 3)
    

          used = send_data (id_shallow_cld_col, tca3, Time_diag, is, js)

        endif 



        if ( id_conv_cld_col > 0) then

          cloud(:,:,:) = REAL (   &
                  COUNT( Model_microphys%stoch_cloud_type(:,:,:,:) == 2,  &
                                                  dim = 4)   +     &
                  COUNT( Model_microphys%stoch_cloud_type(:,:,:,:) == 3, &
                                                  dim = 4))/REAL(ncol)


          tca3(:,:) =  MAXVAL (cloud, dim = 3)

          cloud(:,:,:) =  REAL (   &
                  COUNT( Model_microphys%stoch_cloud_type(:,:,:,:) == 4, &
                                                dim = 4))/REAL(ncol)

          tca3(:,:) =  tca3(:,:) + MAXVAL (cloud,  dim = 3)

          tca3(:,:) = MIN(tca3(:,:), 1.)

          used = send_data (id_conv_cld_col, tca3, Time_diag, is, js)

        endif 


!cms--
!---------------------------------------------------------------------
!    define the high cloud region for each grid column. for each 
!    stochastic band, determine if any levels in the high cloud region 
!    contain cloud. if so define cloud2n to be 1.0 for that band; other-
!    wise it is 0.0. define high cloud percentage by summing over all 
!    bands.
!---------------------------------------------------------------------
        if (id_high_cld_amt > 0)  then
          nn=size(tmplmask4,3)
          do n=1,ncol
            tmplmask4(:,:,:,n) = (Atmos_input%pflux(:,:,1:nn) <= high_btm)
          end do

          cloud2n(:,:,:) =    &
               COUNT (Model_microphys%stoch_cldamt(:,:,:,:) > 0 .and. &
                                        tmplmask4(:,:,:,:), dim = 3)
          cloud2n(:,:,:) = MIN (cloud2n(:,:,:), 1.0)
          hml_ca(:,:,1) = 100.0*SUM (cloud2n(:,:,:), dim = 3)/  &
                                                            REAL (ncol)
          used = send_data (id_high_cld_amt, hml_ca(:,:,1),  &
                            Time_diag, is, js)
        endif

!---------------------------------------------------------------------
!    define the middle cloud region for each grid column. for each 
!    stochastic band, determine if any levels in the middle cloud region
!    contain cloud. if so define cloud2n to be 1.0 for that band; other-
!    wise it is 0.0. define middle cloud percentage by summing over all 
!    bands.
!---------------------------------------------------------------------
        if (id_mid_cld_amt > 0) then    
          nn=size(tmplmask4,3)
          do n=1,ncol
            tmplmask4(:,:,:,n) =     &
                       (Atmos_input%pflux(:,:,1:nn) <= mid_btm .and. &
                             Atmos_input%pflux(:,:,1:nn) > high_btm) 
          end do
                                                
          cloud2n(:,:,:) =    &
               COUNT (Model_microphys%stoch_cldamt(:,:,:,:) > 0 .and. &
                                        tmplmask4(:,:,:,:), dim = 3)
          cloud2n(:,:,:) = MIN (cloud2n(:,:,:), 1.0)
          hml_ca(:,:,2) = 100.0*SUM (cloud2n(:,:,:), dim = 3)/  &
                                                            REAL (ncol)
          used = send_data (id_mid_cld_amt, hml_ca(:,:,2),   &
                            Time_diag, is, js)
        endif

!---------------------------------------------------------------------
!    define the low cloud region for each grid column. for each 
!    stochastic band, determine if any levels in the lowe cloud region
!    contain cloud. if so define cloud2n to be 1.0 for that band; other-
!    wise it is 0.0. define low cloud percentage by summing over all 
!    bands.
!---------------------------------------------------------------------
        if (id_low_cld_amt > 0)  then            
          nn=size(tmplmask4,3)
          do n=1,ncol
            tmplmask4(:,:,:,n) = (Atmos_input%pflux(:,:,1:nn) > mid_btm)
          end do

          cloud2n(:,:,:) =    &
               COUNT (Model_microphys%stoch_cldamt(:,:,:,:) > 0 .and. &
                                       tmplmask4(:,:,:,:), dim = 3)
          cloud2n(:,:,:) = MIN (cloud2n(:,:,:), 1.0)
          hml_ca(:,:,3) = 100.0*SUM (cloud2n(:,:,:), dim = 3)/    &
                                                            REAL (ncol)
          used = send_data (id_low_cld_amt, hml_ca(:,:,3),  &
                            Time_diag, is, js)
        endif
          
!---------------------------------------------------------------------
!    define the combined low and mid cloud region for each grid column. 
!    for each stochastic band, determine if any levels in the low and
!    mid cloud region contain cloud. if so define cloud2n to be 1.0 for 
!    that band; otherwise it is 0.0. define lam cloud percentage by 
!    summing over all bands.
!---------------------------------------------------------------------
        if (id_lam_cld_amt > 0)  then            
          nn=size(tmplmask4,3)
          do n=1,ncol
            tmplmask4(:,:,:,n) = (Atmos_input%pflux(:,:,1:nn) > high_btm)
          end do

          cloud2n(:,:,:) =    &
               COUNT (Model_microphys%stoch_cldamt(:,:,:,:) > 0 .and. &
                                       tmplmask4(:,:,:,:), dim = 3)
          cloud2n(:,:,:) = MIN (cloud2n(:,:,:), 1.0)
          hml_ca(:,:,4) = 100.0*SUM (cloud2n(:,:,:), dim = 3)/    &
                                                            REAL (ncol)
          used = send_data (id_lam_cld_amt, hml_ca(:,:,4),  &
                            Time_diag, is, js)
        endif
          
!---------------------------------------------------------------------
!    when stochastic clouds are not active:
!---------------------------------------------------------------------
      else ! (do_stochastic_clouds)

!---------------------------------------------------------------------
!    define the total cloud amount. 
!---------------------------------------------------------------------
        if (id_tot_cld_amt > 0 ) then
          tca2 = 1.0  
          do k=1,kx        
            tca2(:,:) = tca2(:,:)*(1.0 - Cld_spec%camtsw(:,:,k))
          end do
          tca2 = 100.*(1. - tca2)
          used = send_data (id_tot_cld_amt, tca2, Time_diag, is, js)
        endif

!---------------------------------------------------------------------
!    if high, mid or low cloud diagnostics are desired, call 
!    compute_isccp_clds to define the amount of each. 
!---------------------------------------------------------------------
        if (max(id_high_cld_amt, id_mid_cld_amt, &
                id_low_cld_amt,  id_lam_cld_amt) > 0) then
          call compute_isccp_clds (Atmos_input%pflux, Cld_spec%camtsw, &
                                   Cld_spec%camtsw_band, hml_ca)
   
          used = send_data (id_high_cld_amt, hml_ca(:,:,1), Time_diag, is, js)
          used = send_data (id_mid_cld_amt, hml_ca(:,:,2), Time_diag, is, js)
          used = send_data (id_low_cld_amt, hml_ca(:,:,3), Time_diag, is, js)
          used = send_data (id_lam_cld_amt, hml_ca(:,:,4), Time_diag, is, js)
        endif

!---------------------------------------------------------------------
!
!
!
!                   3 DIMENSIONAL CLOUD AMOUNT
!
!
!
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    send the 3D cloud amount field to diag_manager_mod.
!----------------------------------------------------------------------
        used = send_data (id_cld_amt, 100.*Cld_spec%camtsw,   &
                                      Time_diag, is, js, 1, rmask=mask)

!----------------------------------------------------------------------
!    send the unadjusted sum of the 3D cloud amounts from all active 
!    cloud schemes to diag_manager_mod.
!----------------------------------------------------------------------
        if (id_predicted_cld_amt > 0) then
          cloud (:,:,:) = Lsc_microphys%cldamt(:,:,:)
          if (Cldrad_control%do_donner_deep_clouds) then
            cloud(:,:,:) = cloud(:,:,:) +   &
                                      Cell_microphys%cldamt(:,:,:) + &
                                           Meso_microphys%cldamt(:,:,:)
          endif 
          if (Cldrad_control%do_uw_clouds) then
            cloud(:,:,:) = cloud(:,:,:) +   &
                                      Shallow_microphys%cldamt(:,:,:) 
          endif
          used = send_data (id_predicted_cld_amt, 100.*cloud,    &
                                     Time_diag, is, js, 1, rmask=mask)

        endif 

!----------------------------------------------------------------------
!    send the grid box total cloud-fraction-weighted sums of drop and 
!    ice cloud amounts to diag_manager_mod.
!----------------------------------------------------------------------
        used = send_data (id_all_conc_ice, Model_microphys%conc_ice, &
                                                Time_diag, is, js, 1)

        used = send_data (id_all_conc_drop,   &
                            Model_microphys%conc_drop, Time_diag,   &
                                                           is, js, 1)
      endif  ! (do_stochastic_clouds)

!---------------------------------------------------------------------
!
!
!
!          SHORTWAVE RADIATIVE PROPERTIES OF STRATIFORM CLOUDS
!
!
!
!
!---------------------------------------------------------------------
 
!----------------------------------------------------------------------
!    the following diagnostics are meaningful only when strat_clouds
!    is active:
!----------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds) then

!----------------------------------------------------------------------
!    send the 3D large-scale cloud amount field to diag_manager_mod.
!----------------------------------------------------------------------
 
        used = send_data (id_lsc_cld_amt, 100.*Lsc_microphys%cldamt,   &
                      Time_diag, is, js, 1, rmask=mask)

!----------------------------------------------------------------------
!    send various large-scale cloud shortwave radiative property fields
!    to diag_manager_mod.
!----------------------------------------------------------------------
        used = send_data (id_lsc_cld_ext_uv, Lscrad_props%cldext(:,:,:,iuv), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_lsc_cld_ext_vis, Lscrad_props%cldext(:,:,:,ivis), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_lsc_cld_ext_nir, Lscrad_props%cldext(:,:,:,inir), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_lsc_cld_sct_uv, Lscrad_props%cldsct(:,:,:,iuv), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_lsc_cld_sct_vis, Lscrad_props%cldsct(:,:,:,ivis), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_lsc_cld_sct_nir, Lscrad_props%cldsct(:,:,:,inir), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_lsc_cld_asymm_uv, Lscrad_props%cldasymm(:,:,:,iuv), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_lsc_cld_asymm_vis, Lscrad_props%cldasymm(:,:,:,ivis), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_lsc_cld_asymm_nir, Lscrad_props%cldasymm(:,:,:,inir), &
                          Time_diag, is, js, 1, rmask=mask)
      endif ! (do_strat_clouds)
 
!---------------------------------------------------------------------
!
!
!
!             SHORTWAVE RADIATIVE PROPERTIES OF DONNER CLOUDS
!
!
!
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    the following diagnostics are meaningful only when 
!    donner_deep_clouds is active:
!----------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds) then

!----------------------------------------------------------------------
!    send the 3D cell-scale cloud amount field to diag_manager_mod.
!----------------------------------------------------------------------
        used = send_data (id_cell_cld_amt,100.*Cell_microphys%cldamt,&
                          Time_diag, is, js, 1, rmask=mask)

!----------------------------------------------------------------------
!    send various cell-scale cloud shortwave radiative property fields
!    to diag_manager_mod.
!----------------------------------------------------------------------
        used = send_data (id_cell_cld_ext_uv, Cellrad_props%cldext(:,:,:,iuv), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_cell_cld_ext_vis, Cellrad_props%cldext(:,:,:,ivis), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_cell_cld_ext_nir, Cellrad_props%cldext(:,:,:,inir), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_cell_cld_sct_uv, Cellrad_props%cldsct(:,:,:,iuv), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_cell_cld_sct_vis, Cellrad_props%cldsct(:,:,:,ivis), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_cell_cld_sct_nir, Cellrad_props%cldsct(:,:,:,inir), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_cell_cld_asymm_uv, Cellrad_props%cldasymm(:,:,:,iuv), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_cell_cld_asymm_vis, Cellrad_props%cldasymm(:,:,:,ivis), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_cell_cld_asymm_nir, Cellrad_props%cldasymm(:,:,:,inir), &
                          Time_diag, is, js, 1, rmask=mask )

!----------------------------------------------------------------------
!    send the 3D meso-scale cloud amount field to diag_manager_mod.
!----------------------------------------------------------------------
        used = send_data (id_meso_cld_amt,100.*Meso_microphys%cldamt,&
                          Time_diag, is, js, 1, rmask=mask)

!----------------------------------------------------------------------
!    send various meso-scale cloud shortwave radiative property fields
!    to diag_manager_mod.
!----------------------------------------------------------------------
        used = send_data (id_meso_cld_ext_uv, Mesorad_props%cldext(:,:,:,iuv), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_meso_cld_ext_vis, Mesorad_props%cldext(:,:,:,ivis), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_meso_cld_ext_nir, Mesorad_props%cldext(:,:,:,inir), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_meso_cld_sct_uv, Mesorad_props%cldsct(:,:,:,iuv), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_meso_cld_sct_vis, Mesorad_props%cldsct(:,:,:,ivis), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_meso_cld_sct_nir, Mesorad_props%cldsct(:,:,:,inir), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_meso_cld_asymm_uv, Mesorad_props%cldasymm(:,:,:,iuv), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_meso_cld_asymm_vis, Mesorad_props%cldasymm(:,:,:,ivis), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_meso_cld_asymm_nir, Mesorad_props%cldasymm(:,:,:,inir), &
                          Time_diag, is, js, 1, rmask=mask)
      endif ! (do_donner_deep_clouds)


!---------------------------------------------------------------------
!
!
!
!             SHORTWAVE RADIATIVE PROPERTIES OF UW SHALLOW CLOUDS
!
!
!
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    the following diagnostics are meaningful only when 
!    uw shallow convection is active:
!----------------------------------------------------------------------
      if (Cldrad_control%do_uw_clouds) then

!----------------------------------------------------------------------
!    send the 3D cell-scale cloud amount field to diag_manager_mod.
!----------------------------------------------------------------------
        used = send_data (id_shallow_cld_amt,   &
                               100.*Shallow_microphys%cldamt, &
                                   Time_diag, is, js, 1, rmask=mask)

!----------------------------------------------------------------------
!    send various uw shallow cloud shortwave radiative property fields
!    to diag_manager_mod.
!----------------------------------------------------------------------
        used = send_data (id_shallow_cld_ext_uv, Shallowrad_props%cldext(:,:,:,iuv), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_shallow_cld_ext_vis, Shallowrad_props%cldext(:,:,:,ivis), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_shallow_cld_ext_nir, Shallowrad_props%cldext(:,:,:,inir), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_shallow_cld_sct_uv, Shallowrad_props%cldsct(:,:,:,iuv), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_shallow_cld_sct_vis, Shallowrad_props%cldsct(:,:,:,ivis), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_shallow_cld_sct_nir, Shallowrad_props%cldsct(:,:,:,inir), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_shallow_cld_asymm_uv, Shallowrad_props%cldasymm(:,:,:,iuv), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_shallow_cld_asymm_vis, Shallowrad_props%cldasymm(:,:,:,ivis), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_shallow_cld_asymm_nir, Shallowrad_props%cldasymm(:,:,:,inir), &
                          Time_diag, is, js, 1, rmask=mask )

      endif ! (do_uw_clouds)

!---------------------------------------------------------------------
!
!
!
!             LONGWAVE RADIATIVE PROPERTIES OF CLOUDS
!
!
!
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    a total cloud field emissivity that is the weighted average
!    of the random and max overlap emissivities, over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro) then
        if (id_em_cld_10u > 0) then
          cloud(:,:,:) =    &
             (Cld_spec%crndlw(:,:,:)*Cldrad_props%emrndlw(:,:,:,5,1) + &
              Cld_spec%cmxolw(:,:,:)*Cldrad_props%emmxolw(:,:,:,5,1))/ &
             (Cld_spec%crndlw(:,:,:) + Cld_spec%cmxolw(:,:,:) + 1.0E-10)
          used = send_data (id_em_cld_10u, cloud, Time_diag,     &
                            is, js, 1, rmask=mask)
        endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define a total cloud field emissivity that is the weighted average
!    of the random and max overlap emissivities, over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
      else
        if (id_em_cld_lw > 0   ) then
          cloud(:,:,:) =      &
            (Cld_spec%crndlw(:,:,:)*Cldrad_props%emrndlw(:,:,:,1,1) +  &
             Cld_spec%cmxolw(:,:,:)*Cldrad_props%emmxolw(:,:,:,1,1))/ &
            (Cld_spec%crndlw(:,:,:) + Cld_spec%cmxolw(:,:,:) + 1.0E-10)
          used = send_data (id_em_cld_lw, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif
      endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    a large scale cloud absorption coefficient over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro) then
        used = send_data (id_abs_lsc_cld_10u,    &
                          Lscrad_props%abscoeff(:,:,:,5), Time_diag, &
                          is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define the large scale cloud absorption coefficient over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
      else
        used = send_data (id_abs_lsc_cld_lw,      &
                          Lscrad_props%abscoeff(:,:,:,1), Time_diag, &
                          is, js, 1, rmask=mask)
      endif

      if (Cldrad_control%do_donner_deep_clouds) then
!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    the cell scale cloud absorption coefficient over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
        if (Cldrad_control%do_lw_micro) then
          used = send_data (id_abs_cell_cld_10u,     &
                            Cellrad_props%abscoeff(:,:,:,5), Time_diag,&
                            is, js, 1, rmask=mask)
        else

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define the cell scale cloud absorption coefficient over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
          used = send_data (id_abs_cell_cld_lw,     &
                            Cellrad_props%abscoeff(:,:,:,1), Time_diag,&
                            is, js, 1, rmask=mask)
        endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    a meso-scale cloud absorption coefficient over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
        if (Cldrad_control%do_lw_micro) then
          used = send_data (id_abs_meso_cld_10u,    &
                            Mesorad_props%abscoeff(:,:,:,5), Time_diag,&
                            is, js, 1, rmask=mask)
        else
 
!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define the meso-scale cloud absorption coefficient over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
          used = send_data (id_abs_meso_cld_lw,    &
                            Mesorad_props%abscoeff(:,:,:,1), Time_diag,&
                            is, js, 1, rmask=mask)
        endif

      endif

      if (Cldrad_control%do_uw_clouds) then
!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    the cell scale cloud absorption coefficient over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
        if (Cldrad_control%do_lw_micro) then
          used = send_data (id_abs_shallow_cld_10u,     &
                            Shallowrad_props%abscoeff(:,:,:,5), Time_diag,&
                            is, js, 1, rmask=mask)
        else
 
!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define the cell scale cloud absorption coefficient over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
          used = send_data (id_abs_shallow_cld_lw,     &
                            Shallowrad_props%abscoeff(:,:,:,1), Time_diag,&
                            is, js, 1, rmask=mask)
        endif
      endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    the total-cloud absorption coefficient over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro) then
        used = send_data (id_abs_cld_10u,      &
                          Cldrad_props%abscoeff(:,:,:,5,1), Time_diag,&
                          is, js, 1, rmask=mask)
!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define the total-cloud absorption coefficient over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
      else
        used = send_data (id_abs_cld_lw,    &
                          Cldrad_props%abscoeff(:,:,:,1,1), Time_diag, &
                          is, js, 1, rmask=mask)
      endif

!---------------------------------------------------------------------
!
!
!
!             SHORTWAVE RADIATIVE PROPERTIES OF ALL CLOUDS COMBINED
!
!
!
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!    define total-cloud diagnostics that are associated with the micro-
!    physically-based cloud shortwave radiative properties.
!---------------------------------------------------------------------
      if (Cldrad_control%do_sw_micro) then
        used = send_data (id_ext_cld_uv, Cldrad_props%cldext(:,:,:,iuv,1), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_sct_cld_uv, Cldrad_props%cldsct(:,:,:,iuv,1), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_asymm_cld_uv, 100.0*Cldrad_props%cldasymm(:,:,:,iuv,1), &
                           Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_ext_cld_vis, Cldrad_props%cldext(:,:,:,ivis,1), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_sct_cld_vis, Cldrad_props%cldsct(:,:,:,ivis,1), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_asymm_cld_vis, 100.0*Cldrad_props%cldasymm(:,:,:,ivis,1), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_ext_cld_nir, Cldrad_props%cldext(:,:,:,inir,1), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_sct_cld_nir, Cldrad_props%cldsct(:,:,:,inir,1), &
                          Time_diag, is, js, 1, rmask=mask)
        used = send_data (id_asymm_cld_nir, 100.0*Cldrad_props%cldasymm(:,:,:,inir,1), &
                          Time_diag, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    define total-cloud diagnostics that are associated with the bulk 
!    cloud shortwave radiative properties.
!---------------------------------------------------------------------
      else

!---------------------------------------------------------------------
!    define the reflected ultra-violet.
!---------------------------------------------------------------------

        used = send_data (id_alb_uv_cld, Cldrad_props%cvisrfsw(:,:,:), &
                          Time_diag, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    define the reflected infra-red. 
!---------------------------------------------------------------------

        used = send_data (id_alb_nir_cld, Cldrad_props%cirrfsw(:,:,:), &
                          Time_diag, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    define the absorbed  ultra-violet (not implemented).
!---------------------------------------------------------------------
!       if ( id_abs_uv_cld > 0 ) then
!         cloud = 0.0
!         used = send_data (id_abs_uv_cld, cloud, Time_diag,    &
!                           is, js, 1, rmask=mask)
!       endif

!---------------------------------------------------------------------
!    define the absorbed  infra-red.
!---------------------------------------------------------------------

        used = send_data (id_abs_nir_cld, Cldrad_props%cirabsw(:,:,:), &
                          Time_diag, is, js, 1, rmask=mask)
      endif 

!---------------------------------------------------------------------
!
!
!             STRATIFORM PHYSICAL PROPERTIES
!
!
!
!
!---------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds) then

!---------------------------------------------------------------------
!    ice cloud properties: fractional area, particle size, 
!    cloud amount and path
!---------------------------------------------------------------------
!cms++


          if (id_strat_size_ice_l > 0) then

            tmplmaskl  =  Lsc_microphys%conc_ice > 1.e-3 != 1mg/m3 in-cloud

            used = send_data (id_strat_size_ice_l,   &
                              Lsc_microphys%size_ice, Time_diag,  &
                              is, js, 1, mask=tmplmaskl)
          endif





          if (id_strat_ice_number_l > 0) then

         tmplmaskl  =  Lsc_microphys%conc_ice > 1.e-3 != 1mg/m3 in-cloud

            if (Cldrad_control%do_ice_num) then
              where (Lsc_microphys%cldamt > 0.0)
                cloud = Lsc_microphys%ice_number/   &
                                                   Lsc_microphys%cldamt
!!                 cloud = Lsc_microphys%ice_number
              elsewhere
                cloud = 0.0
              end where
            else
              where (Lsc_microphys%cldamt > 0.0)
                cloud = Lsc_microphys%ice_number
              elsewhere
                cloud = 0.0
              end where
            endif


            used = send_data (id_strat_ice_number_l, cloud, &
                              Time_diag, is, js, 1, mask=tmplmaskl)
          endif


!cms--
        if (max(id_strat_area_ice, id_strat_conc_ice,   &
                id_strat_size_ice, id_lsc_iwp,   &
                                     id_strat_ice_number) > 0) then
          tmplmask = Lsc_microphys%conc_ice > 0.0
          if (id_strat_area_ice > 0) then
            cloud   = 0.
            where (tmplmask)                     
              cloud   = Lsc_microphys%cldamt
            endwhere      
            used = send_data (id_strat_area_ice, cloud, Time_diag,  &
                              is, js, 1, rmask=mask)
          endif

          used = send_data (id_strat_size_ice,   &
                            Lsc_microphys%size_ice, Time_diag,  &
                            is, js, 1, mask=tmplmask)

          used = send_data (id_strat_conc_ice,    &
                            Lsc_microphys%conc_ice, Time_diag,   &
                            is, js, 1, mask=tmplmask)

          if (id_lsc_iwp > 0) then
            cloud2d(:,:) = SUM (Lsc_microphys%conc_ice(:,:,:)*  &
                              pmass2(:,:,:), dim = 3)
            tmplmask2 = cloud2d(:,:) > 0.0
            used = send_data (id_lsc_iwp, cloud2d, Time_diag,   &
                              is, js, mask=tmplmask2)
          endif

          
!cms++

          if (id_strat_ice_number > 0) then
            if (Cldrad_control%do_ice_num) then
              where (Lsc_microphys%cldamt > 0.0)
                cloud = Lsc_microphys%ice_number/   &
                                                   Lsc_microphys%cldamt
!!                 cloud = Lsc_microphys%ice_number
              elsewhere
                cloud = 0.0
              end where
            else
              where (Lsc_microphys%cldamt > 0.0)
                cloud = Lsc_microphys%ice_number
              elsewhere
                cloud = 0.0
              end where
            endif
            used = send_data (id_strat_ice_number, cloud, &
                              Time_diag, is, js, 1, mask=tmplmask)
          endif


!cms--
        endif

        if (id_gb_strat_conc_ice > 0) then
          cloud = Lsc_microphys%conc_ice*Lsc_microphys%cldamt
          used = send_data (id_gb_strat_conc_ice, cloud, &
                            Time_diag, is, js, 1)                
        endif

        if (id_gb_lsc_iwp > 0) then
          cloud2d(:,:) = SUM (Lsc_microphys%conc_ice(:,:,:)*  &
                     Lsc_microphys%cldamt(:,:,:)*pmass2(:,:,:), dim = 3)
          used = send_data (id_gb_lsc_iwp, cloud2d, Time_diag, is, js)
        endif

!---------------------------------------------------------------------
!    water cloud properties:  fractional area, particle size, 
!    cloud amount, path and droplet number
!---------------------------------------------------------------------
!cms++
          if (id_strat_size_drop_l > 0) then

            tmplmaskl = Lsc_microphys%conc_drop > 1.e-3 != 1mg/m3 in-cloud

            used = send_data (id_strat_size_drop_l,   &
                              Lsc_microphys%size_drop, Time_diag,  &
                              is, js, 1, mask=tmplmaskl)
          endif




          if (id_strat_droplet_number_l > 0) then

            tmplmaskl = Lsc_microphys%conc_drop > 1.e-3 != 1mg/m3 in-cloud

            if (Cldrad_control%do_liq_num) then
              cloud = 0.0
              where (Lsc_microphys%cldamt > 0.0)
                cloud = Lsc_microphys%droplet_number/   &
                                                   Lsc_microphys%cldamt
              end where
            else
              cloud = 0.0
              where (Lsc_microphys%cldamt > 0.0)
                cloud = Lsc_microphys%droplet_number
              end where
            endif
            used = send_data (id_strat_droplet_number_l, cloud, &
                              Time_diag, is, js, 1, mask=tmplmaskl)
          endif



!cms--

        if (max(id_strat_area_liq,       id_strat_size_drop,  &
                id_ra_strat_size_drop,   id_strat_conc_drop,  &
                id_strat_droplet_number, id_lsc_lwp) > 0) then
          tmplmask = Lsc_microphys%conc_drop > 0.0

          if (id_strat_area_liq > 0) then
            cloud   = 0.
            where (tmplmask)                      
              cloud   = Lsc_microphys%cldamt
            endwhere      
            used = send_data (id_strat_area_liq, cloud, Time_diag,  &
                              is, js, 1, rmask=mask)
          endif

          used = send_data (id_strat_size_drop,   &
                            Lsc_microphys%size_drop, Time_diag,  &
                            is, js, 1, mask=tmplmask)

          if (id_ra_strat_size_drop > 0) then
            cloud = MAX (  &
                MIN(Lsc_microphys%size_drop, mx_drp_diam), mn_drp_diam)
            used = send_data (id_ra_strat_size_drop, cloud, Time_diag, &
                              is, js, 1, mask=tmplmask)
          endif

          used = send_data (id_strat_conc_drop,   &
                            Lsc_microphys%conc_drop, Time_diag,   &
                            is, js, 1, mask=tmplmask)

          if (id_strat_droplet_number > 0) then
            if (Cldrad_control%do_liq_num) then
              cloud = 0.0
              where (Lsc_microphys%cldamt > 0.0)
                cloud = Lsc_microphys%droplet_number/   &
                                                   Lsc_microphys%cldamt
              end where
            else
              cloud = 0.0
              where (Lsc_microphys%cldamt > 0.0)
                cloud = Lsc_microphys%droplet_number
              end where
            endif
            used = send_data (id_strat_droplet_number, cloud, &
                              Time_diag, is, js, 1, mask=tmplmask)
          endif

          if (id_lsc_lwp > 0) then
            cloud2d(:,:) = SUM (Lsc_microphys%conc_drop(:,:,:)*  &
                                pmass2(:,:,:), dim = 3)
            tmplmask2 = cloud2d(:,:) > 0.0
            used = send_data (id_lsc_lwp, cloud2d, &
                              Time_diag, is, js, mask=tmplmask2)
          endif
        endif

        if (id_gb_strat_conc_drop > 0) then
          cloud = Lsc_microphys%conc_drop*Lsc_microphys%cldamt
          used = send_data (id_gb_strat_conc_drop, cloud, &
                            Time_diag, is, js, 1)                   
        endif

        if (id_gb_lsc_lwp > 0) then
          cloud2d(:,:) = SUM (Lsc_microphys%conc_drop(:,:,:)*  &
                     Lsc_microphys%cldamt(:,:,:)*pmass2(:,:,:), dim = 3)
          used = send_data (id_gb_lsc_lwp, cloud2d, Time_diag, is, js) 
        endif

      endif ! (do_strat_clouds)

!---------------------------------------------------------------------
!
!
!             DONNER MESO PHYSICAL PROPERTIES
!
!
!---------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds) then

!--------------------------------------------------------------------
!    donner meso ice cloud properties: fractional area, size, 
!    cloud amount and path
!--------------------------------------------------------------------
        if (max(id_meso_area_ice, id_meso_size_ice, &
                id_meso_conc_ice, id_meso_iwp) > 0) then
          tmplmask = Meso_microphys%conc_ice > 0.0

          if (id_meso_area_ice > 0) then
            cloud = 0.
            where (tmplmask)                    
              cloud = Meso_microphys%cldamt
            endwhere      
            used = send_data (id_meso_area_ice, cloud, Time_diag,  &
                              is, js, 1, rmask=mask)
          endif

          used = send_data (id_meso_size_ice,   &
                            Meso_microphys%size_ice, Time_diag,  &
                            is, js, 1, mask=tmplmask)

          used = send_data (id_meso_conc_ice,   &
                            Meso_microphys%conc_ice, Time_diag,    &
                            is, js, 1, mask=tmplmask)
   
          if (id_meso_iwp > 0) then
            cloud2d(:,:) =   &
             SUM (Meso_microphys%conc_ice(:,:,:)*pmass(:,:,:), dim = 3 )
            tmplmask2 = cloud2d > 0.0
            used = send_data (id_meso_iwp, cloud2d, Time_diag,   &
                              is, js, mask=tmplmask2)
          endif
        endif

        if (id_gb_meso_conc_ice > 0) then
          cloud = Meso_microphys%conc_ice*Meso_microphys%cldamt 
          used = send_data (id_gb_meso_conc_ice, cloud, &
                            Time_diag, is, js, 1)                  
        endif

        if (id_gb_meso_iwp > 0) then
          cloud2d(:,:) =   &
             SUM (Meso_microphys%conc_ice(:,:,:)*  &
                  Meso_microphys%cldamt(:,:,:)*pmass(:,:,:), dim = 3 )
          used = send_data (id_gb_meso_iwp, cloud2d, Time_diag, is, js)
        endif

!--------------------------------------------------------------------
!    donner meso liquid cloud properties: fractional area, size, 
!    cloud amount, path and droplet number. note that the current
!    donner parameterization does not allow mesoscale liquid.
!--------------------------------------------------------------------
        if (max(id_meso_area_liq,       id_meso_size_drop,  &
                id_ra_meso_size_drop,   id_meso_conc_drop, &
                id_meso_droplet_number, id_meso_lwp) > 0) then
          tmplmask = Meso_microphys%conc_drop > 0.0

          if (id_meso_area_liq > 0) then
            cloud = 0.
            where (tmplmask)                      
              cloud = Meso_microphys%cldamt
            endwhere      
            used = send_data (id_Meso_area_liq, cloud, Time_diag,  &
                              is, js, 1, rmask=mask)
          endif

          used = send_data (id_meso_size_drop,     &
                            Meso_microphys%size_drop, Time_diag,  &
                            is, js, 1, mask=tmplmask)

          if (id_ra_meso_size_drop > 0) then
            cloud = MAX   &
               (MIN(Meso_microphys%size_drop, mx_drp_diam), mn_drp_diam)
            used = send_data (id_ra_meso_size_drop, cloud, Time_diag,  &
                              is, js, 1, mask=tmplmask)
          endif

          used = send_data (id_meso_conc_drop,   &
                            Meso_microphys%conc_drop, Time_diag,  &
                            is, js, 1, mask=tmplmask)

          used = send_data (id_meso_droplet_number,   &
                            Meso_microphys%droplet_number, &
                            Time_diag, is, js, 1, mask=tmplmask)

          if (id_meso_lwp > 0) then
            cloud2d(:,:) = SUM (Meso_microphys%conc_drop(:,:,:)*   &
                                pmass(:,:,:), dim = 3 )
            tmplmask2 = cloud2d > 0.0
            used = send_data (id_meso_lwp, cloud2d, Time_diag,   &
                              is, js, mask=tmplmask2)
          endif
        endif

        if (id_gb_meso_conc_drop > 0) then
          cloud = Meso_microphys%conc_drop*Meso_microphys%cldamt 
          used = send_data (id_gb_meso_conc_drop, cloud, &
                            Time_diag, is, js, 1)                  
        endif

        if (id_gb_meso_lwp > 0) then
          cloud2d(:,:) = SUM (Meso_microphys%conc_drop(:,:,:)*   &
                    Meso_microphys%cldamt(:,:,:)*pmass(:,:,:), dim = 3 )
          used = send_data (id_gb_meso_lwp, cloud2d, Time_diag, is, js)
        endif

!---------------------------------------------------------------------
!
!
!             DONNER CELL PHYSICAL PROPERTIES
!
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    donner cell ice cloud properties: fractional area, size, 
!    cloud amount and path. 
!--------------------------------------------------------------------
        if (max(id_cell_area_ice, id_cell_size_ice,  &
                id_cell_conc_ice, id_cell_iwp) > 0) then
          tmplmask = Cell_microphys%conc_ice > 0.0

          if (id_cell_area_ice > 0) then
            cloud = 0.
            where (tmplmask)                           
              cloud = Cell_microphys%cldamt
            endwhere      
            used = send_data (id_cell_area_ice, cloud, Time_diag,  &
                              is, js, 1, rmask=mask)
          endif

          used = send_data (id_cell_size_ice,    &
                            Cell_microphys%size_ice, Time_diag,  &
                            is, js, 1, mask=tmplmask)

          used = send_data (id_cell_conc_ice,   &
                            Cell_microphys%conc_ice, Time_diag,   &
                            is, js, 1, mask=tmplmask)

          if (id_cell_iwp > 0) then
            cloud2d(:,:) = SUM (Cell_microphys%conc_ice(:,:,:)*   &
                                pmass(:,:,:), dim = 3 )
            tmplmask2 = cloud2d  > 0.0
            used = send_data (id_cell_iwp, cloud2d, &
                              Time_diag, is, js, mask=tmplmask2)
          endif
        endif

        if (id_gb_cell_conc_ice > 0) then
          cloud = Cell_microphys%conc_ice*Cell_microphys%cldamt 
          used = send_data (id_gb_cell_conc_ice, cloud, &
                            Time_diag, is, js, 1)                
        endif

        if (id_gb_cell_iwp > 0) then
          cloud2d(:,:) = SUM (Cell_microphys%conc_ice(:,:,:)*   &
                   Cell_microphys%cldamt(:,:,:)*pmass(:,:,:), dim = 3 )
          used = send_data (id_gb_cell_iwp, cloud2d, &
                            Time_diag, is, js)                
        endif

!--------------------------------------------------------------------
!    donner cell liquid cloud properties: fractional area, size, 
!    cloud amount, path and droplet number.
!--------------------------------------------------------------------
        if (max(id_cell_area_liq,       id_cell_size_drop, &
                id_ra_cell_size_drop,   id_cell_conc_drop, &
                id_cell_droplet_number, id_cell_lwp) > 0) then
          tmplmask = Cell_microphys%conc_drop > 0.0

          if (id_cell_area_liq > 0) then
            cloud = 0.
            where (tmplmask)                       
              cloud = Cell_microphys%cldamt
            endwhere      
            used = send_data (id_cell_area_liq, cloud, Time_diag,  &
                              is, js, 1, rmask=mask)
          endif

          used = send_data (id_cell_size_drop,    &
                            Cell_microphys%size_drop, Time_diag,  &
                            is, js, 1, mask=tmplmask)

          if (id_ra_cell_size_drop > 0) then
            cloud = MAX   &
               (MIN (Cell_microphys%size_drop,mx_drp_diam), mn_drp_diam)
            used = send_data (id_ra_cell_size_drop, cloud, Time_diag,  &
                              is, js, 1, mask=tmplmask)
          endif

          used = send_data (id_cell_conc_drop,    &
                            Cell_microphys%conc_drop,Time_diag,   &
                            is, js, 1, mask=tmplmask)

          used = send_data (id_cell_droplet_number,    &
                            Cell_microphys%droplet_number, &  
                            Time_diag, is, js, 1, mask=tmplmask)

          if (id_cell_lwp > 0) then
            cloud2d = SUM (Cell_microphys%conc_drop(:,:,:)*  &
                                                pmass(:,:,:), dim = 3 )
            tmplmask2 = cloud2d  > 0.0
            used = send_data (id_cell_lwp, cloud2d, &
                              Time_diag, is, js, mask=tmplmask2)
          endif
        endif

        if (id_gb_cell_conc_drop > 0) then
          cloud = Cell_microphys%conc_drop*Cell_microphys%cldamt 
          used = send_data (id_gb_cell_conc_drop, cloud, &
                            Time_diag, is, js, 1)                
        endif

        if (id_gb_cell_lwp > 0) then
          cloud2d = SUM (Cell_microphys%conc_drop(:,:,:)*  &
                  Cell_microphys%cldamt(:,:,:)*pmass(:,:,:), dim = 3 )
          used = send_data (id_gb_cell_lwp, cloud2d, &
                            Time_diag, is, js)                
        endif
      endif ! (do_donner_deep_clouds)

!---------------------------------------------------------------------
!
!
!             UW SHALLOW PHYSICAL PROPERTIES
!
!
!---------------------------------------------------------------------

      if (Cldrad_control%do_uw_clouds) then

!--------------------------------------------------------------------
!    uw shallow ice cloud properties: fractional area, size, 
!    cloud amount and path 
!--------------------------------------------------------------------
        if (max(id_shallow_area_ice, id_shallow_size_ice, &
                id_shallow_conc_ice, id_shallow_iwp) > 0) then
          tmplmask = Shallow_microphys%conc_ice > 0.0

          if (id_shallow_area_ice > 0) then
            cloud = 0.
            where (tmplmask)                           
              cloud = Shallow_microphys%cldamt
            endwhere      
            used = send_data (id_shallow_area_ice, cloud, Time_diag,  &
                              is, js, 1, rmask=mask)
          endif

          used = send_data (id_shallow_size_ice,    &
                            Shallow_microphys%size_ice, Time_diag,  &
                            is, js, 1, mask=tmplmask)

          used = send_data (id_shallow_conc_ice,   &
                            Shallow_microphys%conc_ice, Time_diag,   &
                            is, js, 1, mask=tmplmask)

          if (id_shallow_iwp > 0) then
            cloud2d(:,:) = SUM (Shallow_microphys%conc_ice(:,:,:)*   &
                                pmass(:,:,:), dim = 3 )
            tmplmask2 = cloud2d  > 0.0
            used = send_data (id_shallow_iwp, cloud2d, &
                              Time_diag, is, js, mask=tmplmask2)
          endif
        endif

        if (id_gb_shallow_conc_ice > 0) then
          cloud = Shallow_microphys%conc_ice*Shallow_microphys%cldamt 
          used = send_data (id_gb_shallow_conc_ice, cloud, &
                            Time_diag, is, js, 1)                
        endif

        if (id_gb_shallow_iwp > 0) then
          cloud2d(:,:) = SUM (Shallow_microphys%conc_ice(:,:,:)*   &
                 Shallow_microphys%cldamt(:,:,:)*pmass(:,:,:), dim = 3 )
          used = send_data (id_gb_shallow_iwp, cloud2d, &
                            Time_diag, is, js)                
        endif

!--------------------------------------------------------------------
!    uw shallow liquid cloud properties: fractional area, size, 
!    cloud amount, path and droplet number 
!--------------------------------------------------------------------
        if (max(id_shallow_area_liq,       id_shallow_size_drop, &
                id_ra_shallow_size_drop,   id_shallow_conc_drop, &
                id_shallow_droplet_number, id_shallow_lwp) > 0) then
          tmplmask = Shallow_microphys%conc_drop > 0.0

          if (id_shallow_area_liq > 0) then
            cloud = 0.
            where (tmplmask)                       
              cloud = Shallow_microphys%cldamt
            endwhere      
            used = send_data (id_shallow_area_liq, cloud, Time_diag,  &
                              is, js, 1, rmask=mask)
          endif

          used = send_data (id_shallow_size_drop,    &
                            Shallow_microphys%size_drop, Time_diag,  &
                            is, js, 1, mask=tmplmask)

          if (id_ra_shallow_size_drop > 0) then
            cloud = MAX(              &
                      MIN(Shallow_microphys%size_drop,mx_drp_diam), &
                                                            mn_drp_diam)
            used = send_data (id_ra_shallow_size_drop, cloud,  &
                              Time_diag, is, js, 1, mask=tmplmask)
          endif

          used = send_data (id_shallow_conc_drop,    &
                            Shallow_microphys%conc_drop,Time_diag,   &
                            is, js, 1, mask=tmplmask)

          used = send_data (id_shallow_droplet_number,    &
                            Shallow_microphys%droplet_number, &  
                            Time_diag, is, js, 1, mask=tmplmask)

          if (id_shallow_lwp > 0) then
            cloud2d = SUM (Shallow_microphys%conc_drop(:,:,:)*  &
                                               pmass(:,:,:), dim = 3 )
            tmplmask2 = cloud2d  > 0.0
            used = send_data (id_shallow_lwp, cloud2d, &
                              Time_diag, is, js, mask=tmplmask2)
          endif
        endif

        if (id_gb_shallow_conc_drop > 0) then
          cloud = Shallow_microphys%conc_drop*Shallow_microphys%cldamt 
          used = send_data (id_gb_shallow_conc_drop, cloud, &
                            Time_diag, is, js, 1)                
        endif

        if (id_gb_shallow_lwp > 0) then
          cloud2d = SUM (Shallow_microphys%conc_drop(:,:,:)*  &
                 Shallow_microphys%cldamt(:,:,:)*pmass(:,:,:), dim = 3 )
          used = send_data (id_gb_shallow_lwp, cloud2d, &
                            Time_diag, is, js)                
        endif
      endif ! (do_uw_clouds)

!---------------------------------------------------------------------
!
!
!             STOCHASTIC CLOUD PROPERTIES
!
!
!---------------------------------------------------------------------

      if (Cldrad_control%do_stochastic_clouds) then

!--------------------------------------------------------------------
!
!
!                     "_ONLY_LSC" DIAGNOSTICS
!
!
!--------------------------------------------------------------------
!    these "_only_lsc" diagnostics allow one to assess the effect of
!    treating the non-lsc clouds stochastically. the difference between
!    the "_only_lsc" variables and the corresponding "_ave" variables
!    reflect the changes resulting from treating the non-lsc cloud
!    types stochastically. the cloud properties actually seen by the
!    model's radiation package are always contained in the "_ave" 
!    diagnostics.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    grid-box-mean (averaged over all stochastic bands) large-scale 
!    cloud fraction.
!--------------------------------------------------------------------
        if (id_cldfrac_only_lsc > 0) then
          cloud(:,:,:) =   &
              SUM (Lsc_microphys%stoch_cldamt(:,:,:,:), dim = 4)/ncol
          used = send_data (id_cldfrac_only_lsc, cloud, Time_diag, &
                            is, js, 1, rmask=mask)
        endif
  
!--------------------------------------------------------------------
!    grid-box-mean (averaged over all stochastic bands) large-scale 
!    cloud ice amount and icewater path.
!--------------------------------------------------------------------
        if (id_ice_conc_only_lsc > 0 .or. id_iwp_only_lsc > 0)  then   
          cloud(:,:,:) =    &
             SUM (Lsc_microphys%stoch_conc_ice(:,:,:,:), dim = 4)/ncol 
          used = send_data (id_ice_conc_only_lsc, cloud, Time_diag, &
                            is, js, 1, rmask=mask)
          if (id_iwp_only_lsc > 0)  then       
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3)
            used = send_data (id_iwp_only_lsc, cloud2d, Time_diag, &
                              is, js)
          endif
        endif

!--------------------------------------------------------------------
!    in-cloud (averaged over only cloudy stochastic bands) large-scale 
!    ice water content, ice water path, ice particle size, and fraction
!    of stochastic columns containing ice cloud.
!--------------------------------------------------------------------
        if (max(id_ic_iwp_only_lsc,   id_ic_ice_conc_only_lsc, &  
                id_ice_size_only_lsc, id_ice_col_only_lsc) > 0) then
          tmplmask4(:,:,:,:) =    &
                          Lsc_microphys%stoch_conc_ice(:,:,:,:) > 0.0 
          cloud2(:,:,:) = COUNT (tmplmask4(:,:,:,:), dim = 4)
          tmplmask = cloud2 > 0.0
          cloud(:,:,:) = 0.0
          where (tmplmask)
            cloud(:,:,:) =    &
               SUM (Lsc_microphys%stoch_conc_ice(:,:,:,:),   &
                    mask = tmplmask4, dim = 4)/(cloud2(:,:,:) + 1.0E-40)
          end where

          used = send_data (id_ic_ice_conc_only_lsc, cloud,   &
                            Time_diag, is, js, 1, mask=tmplmask)

          if (id_ic_iwp_only_lsc > 0 ) then      
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3)
            tmplmask2 = SUM (cloud2(:,:,:), dim=3) > 0.0
            used = send_data (id_ic_iwp_only_lsc, cloud2d, Time_diag, &
                              is, js, mask=tmplmask2)
          endif

          if (id_ice_size_only_lsc > 0) then
            cloud(:,:,:) = 0.0
            where (cloud2 > 0)
              cloud(:,:,:) =    &
                SUM (Lsc_microphys%stoch_size_ice(:,:,:,:),   &
                    mask = tmplmask4, dim = 4)/(cloud2(:,:,:) + 1.0E-40)
            end where

            used = send_data (id_ice_size_only_lsc, cloud, Time_diag, &
                              is, js, 1, mask=tmplmask)
          endif

          if (id_ice_col_only_lsc > 0 ) then      
            cloud2 = cloud2/float(ncol)
            used = send_data (id_ice_col_only_lsc, cloud2, Time_diag, &
                              is, js, 1, mask=tmplmask)
          endif
        endif 

!--------------------------------------------------------------------
!    grid-box-mean (averaged over all stochastic bands)large-scale 
!    cloud liquid amount and water path.
!--------------------------------------------------------------------
        if (max(id_drop_conc_only_lsc, id_lwp_only_lsc) > 0 ) then    
          cloud(:,:,:) =   &
             SUM (Lsc_microphys%stoch_conc_drop(:,:,:,:), dim = 4)/ncol

          used = send_data (id_drop_conc_only_lsc, cloud, Time_diag, &
                            is, js, 1, rmask=mask)

          if (id_lwp_only_lsc > 0 ) then    
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3)
            used = send_data (id_lwp_only_lsc, cloud2d, Time_diag, &
                              is, js)
          endif
        endif 

!--------------------------------------------------------------------
!    in-cloud (averaged only over cloudy stochastic columns) large-scale
!    cloud liq water content, liquid water path, droplet size, droplet
!    number, and fraction of stochastic columns containing cloud liquid.
!--------------------------------------------------------------------
        if (max(id_ic_drop_conc_only_lsc, id_ic_lwp_only_lsc, &
                id_liq_col_only_lsc,      id_drop_size_only_lsc, &
                id_ra_drop_size_only_lsc, id_droplet_number_only_lsc) > 0) then
          tmplmask4(:,:,:,:) =    &
                          Lsc_microphys%stoch_conc_drop(:,:,:,:) > 0.0 
          cloud2(:,:,:) = COUNT (tmplmask4(:,:,:,:), dim = 4)
          tmplmask = cloud2 > 0.0
          cloud(:,:,:) = 0.0
          where (tmplmask)   
            cloud(:,:,:) =    &
                SUM (Lsc_microphys%stoch_conc_drop(:,:,:,:),   &
                    mask = tmplmask4, dim = 4)/(cloud2(:,:,:) + 1.0E-40)
          end where

          used = send_data (id_ic_drop_conc_only_lsc, cloud,  &
                            Time_diag, is, js, 1, mask=tmplmask)

          if (id_ic_lwp_only_lsc > 0 ) then  
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3)
            tmplmask2 = SUM (cloud2(:,:,:), dim=3) > 0.0
            used = send_data (id_ic_lwp_only_lsc, cloud2d, Time_diag, &
                              is, js, mask=tmplmask2)
          endif

          if (id_drop_size_only_lsc > 0) then
            cloud(:,:,:) = 0.0
            where (tmplmask)   
              cloud(:,:,:) =    &
                 SUM (Lsc_microphys%stoch_size_drop(:,:,:,:),   &
                   mask = tmplmask4, dim = 4)/(cloud2(:,:,:) + 1.0E-40)
            end where
            used = send_data (id_drop_size_only_lsc, cloud, Time_diag, &
                              is, js, 1, mask=tmplmask)
          endif

          if (id_ra_drop_size_only_lsc > 0) then
            cloud(:,:,:) = 0.0
            where (tmplmask)   
              cloud(:,:,:) = SUM (MIN(MAX(   &
                 Lsc_microphys%stoch_size_drop(:,:,:,:), mn_drp_diam), &
                      mx_drp_diam), mask = tmplmask4, dim = 4)/ &
                                              (cloud2(:,:,:) + 1.0E-40)
            end where
            used = send_data (id_ra_drop_size_only_lsc, cloud,   &
                              Time_diag, is, js, 1, mask=tmplmask)
          endif

          if (id_droplet_number_only_lsc > 0) then
            cloud(:,:,:) = 0.0
            where (tmplmask)   
              cloud(:,:,:) =    &
                 SUM (Lsc_microphys%stoch_droplet_number(:,:,:,:),   &
                    mask = tmplmask4, dim = 4)/(cloud2(:,:,:) + 1.0E-40)
            end where
            used = send_data (id_droplet_number_only_lsc, cloud,   &
                              Time_diag, is, js, 1, mask=tmplmask)
          endif

          if (id_liq_col_only_lsc > 0 ) then  
            cloud2 = cloud2/float(ncol)
            used = send_data (id_liq_col_only_lsc, cloud2, Time_diag, &
                              is, js, 1, mask=tmplmask)
          endif
        endif 

!--------------------------------------------------------------------
!
!
!                     "_AVE" DIAGNOSTICS
!
!
!--------------------------------------------------------------------
!    these diagnostics are for the cloud properties actually seen by 
!    the model's radiation code, which includes contributions from all 
!    active cloud types, stochastically determined.
!--------------------------------------------------------------------

!------------------------------------------------------------------
!    total projected cloud fraction. note that this has been previously
!    calculated as tca2 and output via id_tot_cld_amt.
!------------------------------------------------------------------
        if (id_cldfrac_tot > 0 ) &
          used = send_data (id_cldfrac_tot, 0.01*tca2,   &
                          Time_diag, is, js)

!--------------------------------------------------------------------
!    grid-box-mean cloud fraction in each layer, averaged across all
!    stochastic columns.
!--------------------------------------------------------------------
        if (id_cldfrac_ave > 0) then
          cloud(:,:,:) =   &
              SUM (Model_microphys%stoch_cldamt(:,:,:,:), dim = 4)/ncol 
          used = send_data (id_cldfrac_ave, cloud, Time_diag, &
                            is, js, 1, rmask=mask)
        endif
  
!--------------------------------------------------------------------
!    grid-box-mean (averaged over all stochastic bands) ice cloud  
!    amount and icewater path.
!--------------------------------------------------------------------
        if (max(id_ice_conc_ave, id_iwp_ave) > 0 ) then 
          cloud(:,:,:) =    &
             sum(Model_microphys%stoch_conc_ice(:,:,:,:), dim = 4)/ncol
          used = send_data (id_ice_conc_ave, cloud, Time_diag, &
                            is, js, 1, rmask=mask)
          if (id_iwp_ave > 0 ) then 
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3)
            used = send_data (id_iwp_ave, cloud2d, Time_diag, &
                              is, js                )
          endif
        endif !(id_ice_conc_ave > 0  .or. id_iwp_ave > 0 ) 

!--------------------------------------------------------------------
!    grid-box-mean (averaged over all stochastic bands) ice water path 
!    and ice water amount contributions from large-scale, cell,
!    meso, and shallow clouds.
!--------------------------------------------------------------------
        if (max(id_lsc_iwp_ave, id_lsc_ice_conc_ave) > 0  ) then
          tmplmask4(:,:,:,:) =    &
                       (Model_microphys%stoch_cloud_type(:,:,:,:) == 1) 
          cloud(:,:,:) =   &
                  SUM (Model_microphys%stoch_conc_ice(:,:,:,:),  &
                                       mask=tmplmask4, dim = 4) / ncol
          if (id_lsc_iwp_ave  > 0  ) then
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3)
            used = send_data (id_lsc_iwp_ave, cloud2d,   &
                              Time_diag, is, js)
          endif
          used = send_data (id_lsc_ice_conc_ave, cloud, &
                            Time_diag, is, js,1)
        endif
 
        if (max(id_meso_iwp_ave, id_meso_ice_conc_ave) > 0  ) then
          tmplmask4(:,:,:,:) =    &
                      (Model_microphys%stoch_cloud_type(:,:,:,:) == 2) 
          cloud(:,:,:) =   &
                  SUM (Model_microphys%stoch_conc_ice(:,:,:,:),  &
                                     mask=tmplmask4, dim = 4) / ncol
          if (id_meso_iwp_ave > 0  ) then
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3)
            used = send_data (id_meso_iwp_ave, cloud2d,   &
                              Time_diag, is, js)
          endif
          used = send_data (id_meso_ice_conc_ave, cloud,   &
                            Time_diag, is, js,1)
        endif
        if (max(id_cell_iwp_ave, id_cell_ice_conc_ave) > 0  ) then
          tmplmask4(:,:,:,:) =   &
                   (Model_microphys%stoch_cloud_type(:,:,:,:) == 3)  
          cloud(:,:,:) =   &
              SUM (Model_microphys%stoch_conc_ice(:,:,:,:),  &
                                    mask=tmplmask4, dim = 4) / ncol
          if (id_cell_iwp_ave  > 0  ) then
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3)
            used = send_data (id_cell_iwp_ave, cloud2d,  &
                              Time_diag, is, js)
          endif
          used = send_data (id_cell_ice_conc_ave, cloud,   &
                            Time_diag, is, js,1)
        endif

        if (max(id_shallow_iwp_ave, id_shallow_ice_conc_ave) > 0  ) then
          tmplmask4(:,:,:,:) =    &
                      (Model_microphys%stoch_cloud_type(:,:,:,:) == 4)  
          cloud(:,:,:) =   &
              SUM (Model_microphys%stoch_conc_ice(:,:,:,:),  &
                                      mask=tmplmask4, dim = 4) / ncol
          if (id_shallow_iwp_ave  > 0  ) then
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3)
            used = send_data (id_shallow_iwp_ave, cloud2d,   &
                              Time_diag, is, js)
          endif
          used = send_data (id_shallow_ice_conc_ave, cloud,  &
                            Time_diag, is, js,1)
        endif


!---------------------------------------------------------------------
!    in-cloud (averaged only over cloudy stochastic columns) ice water 
!    content, ice water path, ice particle size, and fraction of 
!    stochastic columns containing cloud ice.
!---------------------------------------------------------------------
        if (max(id_ic_ice_conc_ave, id_ic_iwp_ave, &
                id_ice_size_ave,    id_ice_col_frac_ave) > 0 ) then
          tmplmask4(:,:,:,:) =     &
                   Model_microphys%stoch_conc_ice(:,:,:,:) > 0.0 
          cloud2(:,:,:) = COUNT (tmplmask4(:,:,:,:), dim = 4)
          tmplmask = cloud2 > 0.0
          cloud(:,:,:) = 0.
          where (tmplmask) 
            cloud(:,:,:) =   &
                 SUM (Model_microphys%stoch_conc_ice(:,:,:,:), &
                     mask =tmplmask4, dim = 4)/(cloud2(:,:,:) + 1.0E-40)
          end where

          used = send_data (id_ic_ice_conc_ave, cloud, Time_diag, &
                            is, js, 1, mask=tmplmask)

          if (id_ic_iwp_ave > 0 ) then 
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3)
            tmplmask2 = sum(cloud2(:,:,:), dim=3) > 0.0
            used = send_data (id_ic_iwp_ave, cloud2d, Time_diag, &
                              is, js, mask=tmplmask2)
          endif

          if (id_ice_size_ave > 0) then
            cloud(:,:,:) = 0.
            where (tmplmask) 
              cloud(:,:,:) =   &    
                 SUM (Model_microphys%stoch_size_ice(:,:,:,:), &
                    mask =tmplmask4, dim = 4)/(cloud2(:,:,:) + 1.0E-40)
            end where
            used = send_data (id_ice_size_ave, cloud, Time_diag, &
                              is, js, 1, mask=tmplmask)
          endif

          if (id_ice_col_frac_ave > 0 ) then
            cloud2 = cloud2 / float(ncol)
            used = send_data (id_ice_col_frac_ave, cloud2, Time_diag, &
                              is, js, 1, mask=tmplmask)
          endif
        endif  

!--------------------------------------------------------------------
!    grid-box-mean (averaged over all stochastic bands) liquid cloud  
!    amount and water path.
!--------------------------------------------------------------------
        if (max(id_drop_conc_ave, id_lwp_ave) > 0  ) then  
          cloud(:,:,:) =   &
             SUM (Model_microphys%stoch_conc_drop(:,:,:,:),    &
                                                        dim = 4) / ncol

          used = send_data (id_drop_conc_ave, cloud, Time_diag, &
                            is, js, 1, rmask=mask)

          if (id_lwp_ave > 0  ) then  
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3)
            used = send_data (id_lwp_ave, cloud2d, Time_diag, &
                              is, js                )
          endif
        endif 

!--------------------------------------------------------------------
!    grid-box-mean values (averaged over all stochastic bands) of the
!    contributions to total liquid cloud amount and water path 
!    from large-scale, cell, meso, and shallow clouds
!--------------------------------------------------------------------
        if (max(id_lsc_lwp_ave, id_lsc_drop_conc_ave) > 0  ) then
          tmplmask4(:,:,:,:) =      &
                       (Model_microphys%stoch_cloud_type(:,:,:,:) == 1)
          cloud(:,:,:) =   &
             SUM (Model_microphys%stoch_conc_drop(:,:,:,:), &
                                        mask=tmplmask4, dim = 4) / ncol
          if (id_lsc_lwp_ave > 0  ) then
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3)
            used = send_data (id_lsc_lwp_ave, cloud2d,    &
                              Time_diag, is, js)
          endif
          used = send_data (id_lsc_drop_conc_ave, cloud,     &
                            Time_diag, is, js,1)
        endif

        if (max(id_meso_lwp_ave, id_meso_drop_conc_ave) > 0  ) then
          tmplmask4(:,:,:,:) =    &
                       (Model_microphys%stoch_cloud_type(:,:,:,:) == 2)
          cloud(:,:,:) =   &
             SUM (Model_microphys%stoch_conc_drop(:,:,:,:), &
                                      mask=tmplmask4, dim = 4) / ncol
          if (id_meso_lwp_ave  > 0  ) then
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3)
            used = send_data (id_meso_lwp_ave, cloud2d,    &
                              Time_diag, is, js)
          endif
          used = send_data (id_meso_drop_conc_ave, cloud,   &
                            Time_diag, is, js,1)
        endif

        if (max(id_cell_lwp_ave, id_cell_drop_conc_ave) > 0  ) then
          tmplmask4(:,:,:,:) =     &
                      (Model_microphys%stoch_cloud_type(:,:,:,:) == 3)
          cloud(:,:,:) =   &
              SUM (Model_microphys%stoch_conc_drop(:,:,:,:), &
                                       mask=tmplmask4, dim = 4) / ncol
          if (id_cell_lwp_ave > 0  ) then
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3)
            used = send_data (id_cell_lwp_ave, cloud2d,   &
                              Time_diag, is, js)
          endif
          used = send_data (id_cell_drop_conc_ave, cloud,    &
                            Time_diag, is, js,1)
        endif

        if (max(id_shallow_lwp_ave, id_shallow_drop_conc_ave) > 0  ) then
          tmplmask4(:,:,:,:) =    &
                      (Model_microphys%stoch_cloud_type(:,:,:,:) == 4) 
          cloud(:,:,:) =   &
              SUM (Model_microphys%stoch_conc_drop(:,:,:,:), &
                                       mask=tmplmask4, dim = 4) / ncol
          if (id_shallow_lwp_ave  > 0  ) then
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3)
            used = send_data (id_shallow_lwp_ave, cloud2d,   &
                              Time_diag, is, js)
          endif
          used = send_data (id_shallow_drop_conc_ave, cloud,  &
                            Time_diag, is, js,1)
        endif

!--------------------------------------------------------------------
!    in-cloud (averaged only over cloudy stochastic columns) liquid 
!    water content, liquid water path, droplet size, droplet number
!    and fraction of stochastic columns containing cloud water.
!--------------------------------------------------------------------
        if (max(id_ic_drop_conc_ave, id_ic_lwp_ave,  &
                id_liq_col_frac_ave, id_drop_size_ave, &  
                id_ra_drop_size_ave, id_droplet_number_ave) > 0) then
          tmplmask4(:,:,:,:) =    &
                        Model_microphys%stoch_conc_drop(:,:,:,:) > 0.0 
          cloud2(:,:,:) = COUNT (tmplmask4(:,:,:,:), dim = 4)
          tmplmask = cloud2 > 0.0
          cloud(:,:,:) = 0.0
          where (tmplmask)   
            cloud(:,:,:) =    &
              SUM (Model_microphys%stoch_conc_drop(:,:,:,:),   &
                    mask = tmplmask4, dim = 4)/(cloud2(:,:,:) + 1.0E-40)
          end where
          
          used = send_data (id_ic_drop_conc_ave, cloud, Time_diag, &
                            is, js, 1, mask=tmplmask)

          if (id_ic_lwp_ave > 0  ) then
            cloud2d(:,:) = SUM (cloud(:,:,:)*pmass(:,:,:), dim = 3) 
            tmplmask2 = SUM (cloud2(:,:,:), dim=3) > 0.0
            used = send_data (id_ic_lwp_ave, cloud2d, Time_diag, &
                              is, js, mask=tmplmask2)
          endif

          if (id_drop_size_ave > 0) then
            cloud(:,:,:) = 0.0
            where (tmplmask)   
              cloud(:,:,:) =    &
                SUM (Model_microphys%stoch_size_drop(:,:,:,:),   &
                    mask = tmplmask4, dim = 4)/(cloud2(:,:,:) + 1.0E-40)
            end where
            used = send_data (id_drop_size_ave, cloud, Time_diag, &
                              is, js, 1, mask=tmplmask)
          endif

          if (id_ra_drop_size_ave > 0) then
            cloud(:,:,:) = 0.0
            where (tmplmask)   
              cloud(:,:,:) = SUM (  &
                  MIN(MAX(Model_microphys%stoch_size_drop(:,:,:,:),   &
                                       mn_drp_diam), mx_drp_diam),    &
                    mask = tmplmask4, dim = 4)/(cloud2(:,:,:) + 1.0E-40)
            end where
            used = send_data (id_ra_drop_size_ave, cloud, Time_diag, &
                              is, js, 1, mask=tmplmask)
          endif

          if (id_droplet_number_ave > 0) then
            cloud(:,:,:) = 0.0
            where (tmplmask)   
              cloud(:,:,:) = SUM (    &
                   Model_microphys%stoch_droplet_number(:,:,:,:),   &
                    mask = tmplmask4, dim = 4)/(cloud2(:,:,:) + 1.0E-40)
            end where
            used = send_data (id_droplet_number_ave, cloud, Time_diag,&
                              is, js, 1, mask=tmplmask)
          endif
          if (id_liq_col_frac_ave > 0 ) then   
            cloud2 = cloud2 / float(ncol)
            used = send_data (id_liq_col_frac_ave, cloud2, Time_diag, &
                              is, js, 1, mask=tmplmask)
          endif
        endif 

!cms++



    if (id_drop_size_ave_l > 0 .or. id_droplet_number_ave_l > 0) then

         tmplmask4(:,:,:,:) =    &
                        Model_microphys%stoch_conc_drop(:,:,:,:) >  1.e-3 
                
          cloud2(:,:,:) = COUNT (tmplmask4(:,:,:,:), dim = 4)
          tmplmaskl = cloud2 > 0.0

         if (id_drop_size_ave_l > 0) then
            cloud(:,:,:) = 0.0
            where (tmplmaskl)   
              cloud(:,:,:) =    &
                SUM (Model_microphys%stoch_size_drop(:,:,:,:),   &
                    mask = tmplmask4, dim = 4)/(cloud2(:,:,:) + 1.0E-40)
            end where
            used = send_data (id_drop_size_ave_l, cloud, Time_diag, &
                              is, js, 1, mask=tmplmaskl)
          endif


          if (id_droplet_number_ave_l > 0) then
            cloud(:,:,:) = 0.0
            where (tmplmaskl)   
              cloud(:,:,:) = SUM (    &
                   Model_microphys%stoch_droplet_number(:,:,:,:),   &
                    mask = tmplmask4, dim = 4)/(cloud2(:,:,:) + 1.0E-40)
            end where
            used = send_data (id_droplet_number_ave_l, cloud, Time_diag,&
                              is, js, 1, mask=tmplmaskl)
          endif


       end if



     if (id_ice_size_ave_l > 0 .or.  id_ice_number_ave_l > 0 ) then

                tmplmask4(:,:,:,:) =     &
                   Model_microphys%stoch_conc_ice(:,:,:,:)  > 1.e-3 

          cloud2(:,:,:) = COUNT (tmplmask4(:,:,:,:), dim = 4)
          tmplmaskl = cloud2 > 0.0  

        if (id_ice_size_ave_l > 0 ) then
            cloud(:,:,:) = 0.
            where (tmplmaskl) 
              cloud(:,:,:) =   &    
                 SUM (Model_microphys%stoch_size_ice(:,:,:,:), &
                    mask =tmplmask4, dim = 4)/(cloud2(:,:,:) + 1.0E-40)
            end where
            used = send_data (id_ice_size_ave_l, cloud, Time_diag, &
                              is, js, 1, mask=tmplmaskl)
          endif




          if (id_ice_number_ave_l > 0) then
            cloud(:,:,:) = 0.0
            where (tmplmaskl)   
              cloud(:,:,:) = SUM (    &
                   Model_microphys%stoch_ice_number(:,:,:,:),   &
                    mask = tmplmask4, dim = 4)/(cloud2(:,:,:) + 1.0E-40)
            end where
            used = send_data (id_ice_number_ave_l, cloud, Time_diag,&
                              is, js, 1, mask=tmplmaskl)
          endif



    end if

!cms--

!--------------------------------------------------------------------
!    special diagnostic : lwp / drop size in sw band 7 (visible band)
!--------------------------------------------------------------------
        if (id_LWPr > 0) then
          cloud(:,:,:) = Cld_spec%lwp(:,:,:)/  &
                               Lsc_microphys%stoch_size_drop(:,:,:,14)
          tca(:,:) = SUM (cloud(:,:,:), dim = 3)
          used = send_data (id_LWPr, tca, Time_diag, is, js)
        endif

!---------------------------------------------------------------------
!
!
!                   INDIVIDUAL BAND DIAGNOSTICS
!
!
!---------------------------------------------------------------------

        do n=1,ncol

!--------------------------------------------------------------------
!
!
!                     "_ONLY_LSC" DIAGNOSTICS
!
!
!--------------------------------------------------------------------
!    these "_only_lsc" diagnostics allow one to assess the effect of
!    treating the non-lsc clouds stochastically. the difference between
!    the "_only_lsc" variables and the corresponding variables without
!    that appended tag reflect the changes resulting from treating the 
!    non-lsc cloud types stochastically. the cloud properties actually 
!    seen by the model's radiation package are always contained in the 
!    non "_only_lsc" diagnostics.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    diagnostics for each stochastic band --  cloud fraction, ice water
!    content and path, liquid water content and path, ice particle size,
!    droplet size, droplet number. 
!--------------------------------------------------------------------
          used = send_data (id_cldfrac_cols_only_lsc(n),   &
                            Lsc_microphys%stoch_cldamt(:,:,:,n), &
                            Time_diag, is, js, 1, rmask=mask)

          used = send_data (id_ice_conc_cols_only_lsc(n),    &
                            Lsc_microphys%stoch_conc_ice(:,:,:,n), &
                            Time_diag, is, js, 1, rmask=mask)

          used = send_data (id_drop_conc_cols_only_lsc(n),   &
                            Lsc_microphys%stoch_conc_drop(:,:,:,n), &
                            Time_diag, is, js, 1, rmask=mask)

          if (id_iwp_cols_only_lsc(n) > 0) then
            cloud2d(:,:) =   &
               SUM (Lsc_microphys%stoch_conc_ice(:,:,:,n)*   &
                                                pmass(:,:,:), dim = 3)
               used = send_data (id_iwp_cols_only_lsc(n), cloud2d,  &
                                 Time_diag, is, js)
          endif

          if (id_lwp_cols_only_lsc(n) > 0) then
            cloud2d(:,:) =   &
              SUM (Lsc_microphys%stoch_conc_drop(:,:,:,n)*  &
                                                 pmass(:,:,:), dim = 3)
              used = send_data (id_lwp_cols_only_lsc(n), cloud2d,   &
                                Time_diag, is, js)
          endif

          if (id_ice_size_cols_only_lsc(n) > 0) then
            tmplmask = Lsc_microphys%stoch_conc_ice(:,:,:,n) > 0.0
            used = send_data (id_ice_size_cols_only_lsc(n),   &
                              Lsc_microphys%stoch_size_ice(:,:,:,n), &
                              Time_diag, is, js, 1, mask=tmplmask)
          endif

          if (max(id_drop_size_cols_only_lsc(n),    &
                  id_ra_drop_size_cols_only_lsc(n), &
                  id_droplet_number_cols_only_lsc(n)) > 0) then
            tmplmask = Lsc_microphys%stoch_conc_drop(:,:,:,n) > 0.0

            used = send_data (id_drop_size_cols_only_lsc(n),   &
                              Lsc_microphys%stoch_size_drop(:,:,:,n),&
                              Time_diag, is, js, 1, mask=tmplmask)

            if (id_ra_drop_size_cols_only_lsc(n) > 0) then
              cloud(:,:,:) = MAX (MIN   &
                       (Lsc_microphys%stoch_size_drop(:,:,:,n),   &
                                          mx_drp_diam), mn_drp_diam)
              used = send_data (id_ra_drop_size_cols_only_lsc(n),  &
                                cloud, Time_diag, is, js, 1,   &
                                mask=tmplmask)
            endif

            used = send_data (id_droplet_number_cols_only_lsc(n),   &
                        Lsc_microphys%stoch_droplet_number(:,:,:,n),&
                              Time_diag, is, js, 1, mask=tmplmask)
          endif ! (3 options)

!---------------------------------------------------------------------
!
!
!              CLOUD PROPERTIES SEEN BY RADIATION CODE 
!
!
!---------------------------------------------------------------------
          used = send_data (id_cldfrac_cols(n),   &
                            Model_microphys%stoch_cldamt(:,:,:,n), &
                            Time_diag, is, js, 1, rmask=mask)

          used = send_data (id_ice_conc_cols(n),   &
                            Model_microphys%stoch_conc_ice(:,:,:,n),&
                            Time_diag, is, js, 1, rmask=mask)

          used = send_data (id_drop_conc_cols(n),   &
                            Model_microphys%stoch_conc_drop(:,:,:,n),&
                            Time_diag, is, js, 1, rmask=mask)

          if (id_iwp_cols(n) > 0) then
            cloud2d(:,:) = SUM (  &
                 Model_microphys%stoch_conc_ice(:,:,:,n)*  &
                                                pmass(:,:,:), dim = 3)
            used = send_data (id_iwp_cols(n), cloud2d, Time_diag, &
                              is, js)
          endif

          if (id_lwp_cols(n) > 0) then
            cloud2d(:,:) = SUM (  &
                Model_microphys%stoch_conc_drop(:,:,:,n)*    &
                                                pmass(:,:,:), dim = 3)
            used = send_data (id_lwp_cols(n), cloud2d, Time_diag, &
                              is, js)
          endif

          if (id_ice_size_cols(n) > 0) then
            tmplmask = Model_microphys%stoch_conc_ice(:,:,:,n) > 0.0
            used = send_data (id_ice_size_cols(n),   &
                              Model_microphys%stoch_size_ice(:,:,:,n), &
                              Time_diag, is, js, 1, mask=tmplmask)
          endif

          if (max(id_drop_size_cols(n),    &
                  id_ra_drop_size_cols(n), &
                  id_droplet_number_cols(n)) > 0) then
            tmplmask = Model_microphys%stoch_conc_drop(:,:,:,n) > 0.0

            used = send_data (id_drop_size_cols(n),   &
                            Model_microphys%stoch_size_drop(:,:,:,n),&
                            Time_diag, is, js, 1, mask=tmplmask)

            if (id_ra_drop_size_cols(n) > 0) then
              cloud(:,:,:) = MAX(MIN(      &
                   Model_microphys%stoch_size_drop(:,:,:,n),  &
                                             mx_drp_diam), mn_drp_diam)
              used = send_data (id_ra_drop_size_cols(n), cloud,   &
                                Time_diag, is, js, 1, mask=tmplmask)
            endif

            used = send_data (id_droplet_number_cols(n),  &
                      Model_microphys%stoch_droplet_number(:,:,:,n),&
                              Time_diag, is, js, 1, mask=tmplmask)
          endif ! (3 options)
        end do

!----------------------------------------------------------------------
!    frequency of occurrence of large-scale, donner meso and cell and 
!    uw shallow clouds in cloudy stochastic columns (stoch_ic_xxx_...).
!    frequency of seeing various cloud types (largescale, donner meso 
!    and cell, uw shallow) when they are present (stoch_sees_xxx), and 
!    the grid-box-mean frequency of their being seen by the radiation 
!    package, averaged over all stochastic columns  (stoch_xxx_cf_ave). 
!----------------------------------------------------------------------
        if (max(id_stoch_ic_shallow_cf_ave, id_stoch_ic_cell_cf_ave, &
                id_stoch_ic_meso_cf_ave,    id_stoch_ic_lsc_cf_ave,  &
                id_stoch_sees_lsc,          id_stoch_sees_meso,      &
                id_stoch_sees_cell,         id_stoch_sees_shallow,   &
                id_stoch_lsc_cf_ave,        id_stoch_meso_cf_ave,    &
                id_stoch_cell_cf_ave,       id_stoch_shallow_cf_ave) > 0 ) then
          cloud(:,:,:) =   &
                  SUM (Model_microphys%stoch_cldamt(:,:,:,:), dim = 4) 
          tmplmask = cloud > 0.0
         
          if (Cldrad_control%do_strat_clouds) then
            if (max(id_stoch_ic_lsc_cf_ave, id_stoch_sees_lsc, &
                    id_stoch_lsc_cf_ave) > 0 )  then
              cloud(:,:,:) = REAL (COUNT(    &
                  Model_microphys%stoch_cloud_type(:,:,:,:) == 1,    &
                                                   dim = 4))/REAL(ncol)
              used = send_data (id_stoch_ic_lsc_cf_ave, cloud, &
                                Time_diag, is, js, 1, mask=tmplmask)
              if (id_stoch_sees_lsc > 0) then
                tmplmaska =  Lsc_microphys%cldamt > 0.0    
                used = send_data (id_stoch_sees_lsc, cloud,  &
                                  Time_diag, is, js, 1, mask=tmplmaska)
              endif
              used = send_data (id_stoch_lsc_cf_ave, cloud,  &
                                Time_diag, is, js, 1)
            endif
          endif

          if (Cldrad_control%do_donner_deep_clouds) then
            if (max(id_stoch_ic_meso_cf_ave, id_stoch_sees_meso, &
                    id_stoch_meso_cf_ave) > 0 )  then
              cloud(:,:,:) = REAL (COUNT(    &
                  Model_microphys%stoch_cloud_type(:,:,:,:) == 2,    &
                                                   dim = 4))/REAL(ncol)
              used = send_data (id_stoch_ic_meso_cf_ave, cloud, &
                                Time_diag, is, js, 1, mask=tmplmask)
              if (id_stoch_sees_meso > 0) then
                tmplmaska =  Meso_microphys%cldamt > 0.0    
                used = send_data (id_stoch_sees_meso, cloud,  &
                                  Time_diag, is, js, 1, mask=tmplmaska)
              endif
              used = send_data (id_stoch_meso_cf_ave, cloud,  &
                                Time_diag, is, js, 1)
            endif
            if (max(id_stoch_ic_cell_cf_ave, id_stoch_sees_cell, &
                    id_stoch_cell_cf_ave) > 0 )  then
              cloud(:,:,:) = REAL (COUNT(    &
                  Model_microphys%stoch_cloud_type(:,:,:,:) == 3,    &
                                                   dim = 4))/REAL(ncol)
              used = send_data (id_stoch_ic_cell_cf_ave, cloud, &
                                Time_diag, is, js, 1, mask=tmplmask)
              if (id_stoch_sees_cell > 0) then
                tmplmaska =  Cell_microphys%cldamt > 0.0    
                used = send_data (id_stoch_sees_cell, cloud,  &
                                  Time_diag, is, js, 1, mask=tmplmaska)
              endif
              used = send_data (id_stoch_cell_cf_ave, cloud,  &
                                Time_diag, is, js, 1)
            endif
          endif

          if (Cldrad_control%do_uw_clouds) then
            if (max(id_stoch_ic_shallow_cf_ave, id_stoch_sees_shallow, &
                    id_stoch_shallow_cf_ave) > 0 )  then
              cloud(:,:,:) = REAL (COUNT(    &
                  Model_microphys%stoch_cloud_type(:,:,:,:) == 4,    &
                                                   dim = 4))/REAL(ncol)
              used = send_data (id_stoch_ic_shallow_cf_ave, cloud, &
                                Time_diag, is, js, 1, mask=tmplmask)
              if (id_stoch_sees_shallow > 0) then
                tmplmaska =  Shallow_microphys%cldamt > 0.0    
                used = send_data (id_stoch_sees_shallow, cloud,  &
                                  Time_diag, is, js, 1, mask=tmplmaska)
              endif
              used = send_data (id_stoch_shallow_cf_ave, cloud,  &
                                Time_diag, is, js, 1)
            endif
          endif
        endif

!----------------------------------------------------------------------
!    cloud type assigned to each stochastic column.
!----------------------------------------------------------------------
        do n=1,ncol    
          if (id_stoch_cloud_type(n) > 0) then
            tmplmask = Model_microphys%stoch_cloud_type(:, :, :, n) /= 0
            used = send_data   &
                       (id_stoch_cloud_type(n),   &
                       REAL(Model_microphys%stoch_cloud_type(:,:,:,n)),&
                       Time_diag, is, js, 1, mask = tmplmask)
          endif
        end do 

      endif ! (do_stochastic_clouds)

endif ! (Time_Diag > Time)


!---------------------------------------------------------------------



end subroutine cloudrad_netcdf



!#####################################################################

subroutine model_micro_dealloc (Model_microphys)
 
type(microphysics_type), intent(inout) :: Model_microphys


!--------------------------------------------------------------------
!    deallocate the components of the microphysics_type derived type
!    variable.
!--------------------------------------------------------------------
     if (Cldrad_control%do_stochastic_clouds) then
        nullify (Model_microphys%lw_stoch_conc_ice)
        nullify (Model_microphys%lw_stoch_conc_drop)
        nullify (Model_microphys%lw_stoch_size_ice)
        nullify (Model_microphys%lw_stoch_size_drop)
        nullify (Model_microphys%lw_stoch_cldamt) 
        nullify (Model_microphys%lw_stoch_droplet_number)
        nullify (Model_microphys%sw_stoch_conc_ice)
        nullify (Model_microphys%sw_stoch_conc_drop)
        nullify (Model_microphys%sw_stoch_size_ice)
        nullify (Model_microphys%sw_stoch_size_drop)
        nullify (Model_microphys%sw_stoch_cldamt)
        nullify (Model_microphys%sw_stoch_droplet_number)
        deallocate (Model_microphys%stoch_conc_ice)
        deallocate (Model_microphys%stoch_conc_drop)
        deallocate (Model_microphys%stoch_size_ice)
        deallocate (Model_microphys%stoch_size_drop)
        deallocate (Model_microphys%stoch_cldamt)  
        deallocate (Model_microphys%stoch_cloud_type)
        deallocate (Model_microphys%stoch_droplet_number)
        deallocate (Model_microphys%stoch_ice_number)
      endif ! (do_stochastic_clouds)

      if ( associated(Model_microphys%conc_drop) )      deallocate (Model_microphys%conc_drop   )
      if ( associated(Model_microphys%conc_ice) )       deallocate (Model_microphys%conc_ice    )
      if ( associated(Model_microphys%conc_rain) )      deallocate (Model_microphys%conc_rain   )
      if ( associated(Model_microphys%conc_snow) )      deallocate (Model_microphys%conc_snow   )
      if ( associated(Model_microphys%size_drop) )      deallocate (Model_microphys%size_drop   )
      if ( associated(Model_microphys%size_ice) )       deallocate (Model_microphys%size_ice    )
      if ( associated(Model_microphys%size_rain) )      deallocate (Model_microphys%size_rain   )
      if ( associated(Model_microphys%size_snow) )      deallocate (Model_microphys%size_snow  )
      if ( associated(Model_microphys%cldamt) )         deallocate (Model_microphys%cldamt      )
      if ( associated(Model_microphys%droplet_number) ) deallocate (Model_microphys%droplet_number  )
      if ( associated(Model_microphys%ice_number) )     deallocate (Model_microphys%ice_number  )
       
!------------------------------------------------------------------



end subroutine model_micro_dealloc


!####################################################################
! <SUBROUTINE NAME="cloudrad_diagnostics_end">
!  <OVERVIEW>
!    cloudrad_diagnostics_end is the destructor for 
!    cloudrad_diagnostics_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloudrad_diagnostics_end is the destructor for 
!    cloudrad_diagnostics_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloudrad_diagnostics_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine cloudrad_diagnostics_end

!-------------------------------------------------------------------
!    cloudrad_diagnostics_end is the destructor for 
!    cloudrad_diagnostics_mod.
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('cloudrad_diagnostics_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!--------------------------------------------------------------------
!    close out the component modules.
!--------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds) then
        if (do_isccp) call isccp_clouds_end
      endif

!--------------------------------------------------------------------
!    mark the module as not initialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------



end subroutine cloudrad_diagnostics_end




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!####################################################################
! <SUBROUTINE NAME="diag_field_init">
!  <OVERVIEW>
!    diag_field_init registers the potential netcdf output variables
!    with diag_manager_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    diag_field_init registers the potential netcdf output variables
!    with diag_manager_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diag_field_init (axes, Time)
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="real">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
! </SUBROUTINE>
!
subroutine diag_field_init (Time, axes )

!---------------------------------------------------------------------
!    diag_field_init registers the potential netcdf output variables
!    with diag_manager_mod.
!---------------------------------------------------------------------

type(time_type), intent(in) :: Time
integer        , intent(in) :: axes(4)

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       Time      initialization time for the netcdf output fields
!       axes      diagnostic variable axes
!
!---------------------------------------------------------------------
      character(len=8) :: chvers
      integer          :: n

!---------------------------------------------------------------------
!    register the MODIS-related diagnostic fields in this module.
!---------------------------------------------------------------------
      id_reff_modis = register_diag_field    &
                         (mod_name, 'reff_modis', axes(1:2), Time, &
                          'MODIS effective radius', 'micron')

      id_reff_modis2 = register_diag_field    &
                         (mod_name, 'reff_modis2', axes(1:2), Time, &
                          'MODIS effective radius frequency', 'count')

      if ((id_reff_modis < 0 .and. id_reff_modis2 > 0) .or. &
          (id_reff_modis > 0 .and. id_reff_modis2 < 0)) then
        call error_mesg ('cloudrad_diagnostics_mod,diag_field_init:',&
          'both reff_modis and reff_modis2 must either  &
                                       &be active or inactive', FATAL)
      endif
      id_reff_modis3 = register_diag_field    &
                         (mod_name, 'reff_modis3', axes(1:2), Time, &
                          'MODIS scan pressure level', 'mbar')

      id_cldtop_reff = register_diag_field    &
                         (mod_name, 'cldtop_reff', axes(1:2), Time, &
                          'liq drop radius at cld top*cfrac', 'meters')
     
      id_cldtop_area = register_diag_field    &
                         (mod_name, 'cldtop_area', axes(1:2), Time, &
                          'liq cloud area at cld top', '1')
     
      id_cldtop_dropnum = register_diag_field    &
                         (mod_name, 'cldtop_dropnum', axes(1:2), Time, &
                          'liq droplet # at cld top*cfrac', 'm-3')
     
      id_dropnum_col = register_diag_field    &
                         (mod_name, 'dropnum_col', axes(1:2), Time, &
                    'column integrated liq droplet # * cfrac', 'm-2')
     
!---------------------------------------------------------------------
!    register various cloud fraction diagnostics.
!---------------------------------------------------------------------
      id_tot_cld_amt = register_diag_field    &
                         (mod_name, 'tot_cld_amt', axes(1:2), Time, &
                          'total cloud amount', 'percent')

      id_high_cld_amt = register_diag_field   &
                         (mod_name, 'high_cld_amt', axes(1:2), Time, &
                          'high cloud amount', 'percent')

      id_mid_cld_amt =  register_diag_field     &
                         (mod_name, 'mid_cld_amt', axes(1:2), Time, &
                          'mid cloud amount', 'percent')
  
      id_low_cld_amt = register_diag_field    &
                         (mod_name, 'low_cld_amt', axes(1:2), Time, &
                          'low cloud amount', 'percent')

      id_lam_cld_amt = register_diag_field    &
                         (mod_name, 'lam_cld_amt', axes(1:2), Time, &
                          'low and mid cloud amount', 'percent')

      if ( .not. Cldrad_control%do_stochastic_clouds) then
        id_cld_amt =  register_diag_field     &
                         (mod_name, 'cld_amt', axes(1:3), Time,      &
                          'cloud amount', 'percent',     &
                          missing_value=missing_value)

        id_predicted_cld_amt = register_diag_field    &
                         (mod_name, 'predicted_cld_amt', axes(1:3), &
                          Time, 'total raw predicted cloud amount',   &
                          'percent', missing_value=missing_value)

      endif


!cms++

      id_lsc_cld_col = register_diag_field    &
                         (mod_name, 'lsc_cld_col', axes(1:2), Time, &
                          'lsc_cld_col', 'fraction')

      id_shallow_cld_col = register_diag_field    &
                         (mod_name, 'shallow_cld_col', axes(1:2), Time, &
                          'shallow_cld_col', 'fraction')

      id_meso_cld_col = register_diag_field    &
                         (mod_name, 'meso_cld_col', axes(1:2), Time, &
                          'meso_cld_col', 'fraction')


      id_cell_cld_col = register_diag_field    &
                         (mod_name, 'cell_cld_col', axes(1:2), Time, &
                          'cell_cld_col', 'fraction')


      id_conv_cld_col = register_diag_field    &
                         (mod_name, 'conv_cld_col', axes(1:2), Time, &
                          'conv_cld_col', 'fraction')


!cms--
!---------------------------------------------------------------------
!    register lw cloud radiative properties diagnostics - all active
!    clouds, as seen by radiation package.
!---------------------------------------------------------------------
      id_em_cld_lw =  register_diag_field    &
                         (mod_name, 'em_cld_lw', axes(1:3), Time, &
                          'lw cloud emissivity', 'percent',        &
                           missing_value=missing_value)
 
      id_em_cld_10u = register_diag_field    &
                         (mod_name, 'em_cld_10u', axes(1:3), Time, &
                          'cloud emissivity 10 um band', 'percent',    &
                          missing_value=missing_value)

      id_abs_cld_lw = register_diag_field    &
                         (mod_name, 'abs_lw', axes(1:3), Time, &
                          'cloud abs coeff lw', 'percent',        &
                          missing_value=missing_value)

      id_abs_cld_10u = register_diag_field     &
                         (mod_name, 'abs_10u', axes(1:3), Time, &
                          'cloud abs coeff 10um band', 'percent',    &
                          missing_value=missing_value)

!---------------------------------------------------------------------
!    register diagnostic fields associated with the bulk shortwave
!    parameterization, all active clouds, as seen by radiation package.
!---------------------------------------------------------------------
      if (.not. Cldrad_control%do_sw_micro) then
        id_alb_uv_cld = register_diag_field     &
                         (mod_name, 'alb_uv_cld', axes(1:3), Time, &
                          'UV reflected by cloud', 'percent',       &
                          missing_value=missing_value)

        id_alb_nir_cld = register_diag_field      &
                         (mod_name, 'alb_nir_cld', axes(1:3), Time, &
                          'IR reflected by cloud', 'percent',        &
                          missing_value=missing_value)

!   --- do not output this field ---
!       id_abs_uv_cld =  register_diag_field    &
!                        (mod_name, 'abs_uv_cld', axes(1:3), Time, &
!                         'UV absorbed by cloud', 'percent',        &
!                         missing_value=missing_value)

        id_abs_nir_cld = register_diag_field     &
                         (mod_name, 'abs_nir_cld', axes(1:3), Time, &
                          'IR absorbed by cloud', 'percent',         &
                          missing_value=missing_value)

!---------------------------------------------------------------------
!    register diagnostic fields associated with the microphysically
!    based shortwave parameterization, all active clouds, as seen by
!    radiation package.
!---------------------------------------------------------------------
      else 
        id_ext_cld_uv = register_diag_field       &
                         (mod_name, 'ext_cld_uv', axes(1:3), Time, &
                          '.27um cloud extinction coeff', 'km-1',  &
                         missing_value=missing_value)

        id_sct_cld_uv = register_diag_field     &
                         (mod_name, 'sct_cld_uv', axes(1:3), Time, &
                          '.27um cloud scattering coeff', 'km-1', &
                          missing_value=missing_value)

        id_asymm_cld_uv = register_diag_field    &
                         (mod_name, 'asymm_cld_uv', axes(1:3), Time, &
                          '.27um cloud asymmetry parameter',   &
                          'percent', missing_value=missing_value) 

        id_ext_cld_vis =  register_diag_field     &
                         (mod_name, 'ext_cld_vis', axes(1:3), Time, &
                          '.55um cloud extinction coeff', 'km-1', &
                          missing_value=missing_value)

        id_sct_cld_vis = register_diag_field    &
                         (mod_name, 'sct_cld_vis', axes(1:3), Time, &
                          '.55um cloud scattering coeff', 'km-1', &
                          missing_value=missing_value)

        id_asymm_cld_vis = register_diag_field      &
                         (mod_name, 'asymm_cld_vis', axes(1:3), Time,&
                          '.55um cloud asymmetry parameter',   &
                          'percent', missing_value=missing_value)

        id_ext_cld_nir = register_diag_field    &
                         (mod_name, 'ext_cld_nir', axes(1:3), Time, &
                          '1.4um cloud extinction coeff', 'km-1', &
                          missing_value=missing_value)

        id_sct_cld_nir = register_diag_field    &
                         (mod_name, 'sct_cld_nir', axes(1:3), Time, &
                          '1.4um cloud scattering coeff', 'km-1', &
                          missing_value=missing_value)
 
        id_asymm_cld_nir = register_diag_field   &
                         (mod_name, 'asymm_cld_nir', axes(1:3), Time,&
                          '1.4um cloud asymmetry parameter',   &
                          'percent', missing_value=missing_value)

        id_largescale_opdepth = register_diag_field   &
                         (mod_name, 'largescale_opdepth', axes(1:3), &
                          Time, '.55um cloud optical depth avgd over &
                          &subcolumns with ls cloud', &
                          'none', missing_value=missing_value)
 
        id_convect_opdepth = register_diag_field   &
                         (mod_name, 'convect_opdepth', axes(1:3), &
                          Time, '.55um cloud optical depth avgd over &
                          &subcolumns with convective clouds', &
                          'none', missing_value=missing_value)
 
        id_total_opdepth = register_diag_field   &
                         (mod_name, 'total_opdepth', axes(1:3), &
                          Time, '.55um cloud optical depth avgd over &
                          &all subcolumns with cloud', &
                          'none', missing_value=missing_value)
 

!---------------------------------------------------------------------
!    register the microphysically-based cloud radiative property
!    diagnostics resulting from the large-scale clouds only.
!---------------------------------------------------------------------
        id_lsc_cld_ext_uv = register_diag_field    &
                         (mod_name, 'lsc_cld_ext_uv', axes(1:3),   &
                          Time, '.27um lsc cloud ext coeff', 'km-1',&
                          missing_value=missing_value)

        id_lsc_cld_ext_vis = register_diag_field   &
                         (mod_name, 'lsc_cld_ext_vis', axes(1:3), &
                          Time, '.55um lsc cloud ext coeff',   &
                          'km-1', missing_value=missing_value)
 
        id_strat_opdepth = register_diag_field   &
                         (mod_name, 'strat_opdepth', axes(1:3), &
                          Time, '.55um cloud optical depth avgd over &
                          &subcolumns with strat clouds',   &
                          'none', missing_value=missing_value)
 
        id_lsc_cld_ext_nir = register_diag_field    &
                         (mod_name, 'lsc_cld_ext_nir', axes(1:3),  &
                          Time, '1.4um lsc cloud ext coeff',   &
                          'km-1', missing_value=missing_value)

        id_lsc_cld_sct_uv = register_diag_field    &
                         (mod_name, 'lsc_cld_sct_uv', axes(1:3),  &
                          Time, '.27um lsc cloud sct coeff', 'km-1',&
                          missing_value=missing_value)

        id_lsc_cld_sct_vis = register_diag_field    &
                         (mod_name, 'lsc_cld_sct_vis', axes(1:3),  &
                          Time, '.55um lsc cloud sct coeff',  &
                          'km-1', missing_value=missing_value)

        id_lsc_cld_sct_nir = register_diag_field    &
                         (mod_name, 'lsc_cld_sct_nir', axes(1:3), &
                          Time, '1.4um lsc cloud sct coeff',  &
                          'km-1', missing_value=missing_value)
 
        id_lsc_cld_asymm_uv = register_diag_field   &
                         (mod_name, 'lsc_cld_asymm_uv', axes(1:3),&
                          Time, '.27um lsc cloud asymm coeff',  &
                          'percent', missing_value=missing_value)

        id_lsc_cld_asymm_vis = register_diag_field  &
                         (mod_name, 'lsc_cld_asymm_vis', axes(1:3),  &
                          Time,  '.55um lsc cloud asymm coeff',   &
                          'percent', missing_value=missing_value)

        id_lsc_cld_asymm_nir = register_diag_field   &
                         (mod_name, 'lsc_cld_asymm_nir', axes(1:3),  &
                          Time, '1.4um lsc cloud asymm coeff', &
                          'percent', missing_value=missing_value)

      endif

!---------------------------------------------------------------------
!    register the microphysically-based cloud amount and lw cloud
!    radiative properties diagnostic fields.
!---------------------------------------------------------------------
      id_lsc_cld_amt = register_diag_field    &
                         (mod_name, 'lsc_cld_amt', axes(1:3), Time, &
                          'lsc cloud amount', 'percent',             &
                          missing_value=missing_value)

      id_abs_lsc_cld_lw = register_diag_field    &
                         (mod_name, 'lsc_abs_lw', axes(1:3), Time, &
                          'lsc cloud abs coeff lw', 'percent',     &
                          missing_value=missing_value)

      id_abs_lsc_cld_10u = register_diag_field   &
                         (mod_name, 'lsc_abs_10u', axes(1:3), Time,&
                          'lsc cloud abs coeff 10um band',   &
                          'percent', missing_value=missing_value)

!---------------------------------------------------------------------
!    register the donner cell-scale cloud radiative property diagnostic 
!    fields.
!---------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds_iz) then
        if (Cldrad_control%do_donner_deep_clouds) then
          id_cell_cld_amt = register_diag_field    &
                         (mod_name, 'cell_cld_amt', axes(1:3), Time,&
                          'cell cloud amount', 'percent',           &
                          missing_value=missing_value)

          id_cell_cld_ext_uv = register_diag_field   &
                         (mod_name, 'cell_cld_ext_uv', axes(1:3),&
                          Time, '.27um cell cloud ext coeff',  &
                          'km-1', missing_value=missing_value)

          id_cell_cld_ext_vis = register_diag_field   &
                         (mod_name, 'cell_cld_ext_vis', axes(1:3),  &
                          Time, '.55um cell cloud ext coeff',  &
                          'km-1', missing_value=missing_value)

          id_cell_opdepth = register_diag_field   &
                         (mod_name, 'cell_opdepth', axes(1:3),  &
                          Time, '.55um cloud optical depth avgd over &
                          &subcolumns with cell clouds',  &
                          'none', missing_value=missing_value)

          id_cell_cld_ext_nir = register_diag_field   &
                         (mod_name, 'cell_cld_ext_nir', axes(1:3),  &
                          Time, '1.4um cell cloud ext coeff',  &
                          'km-1', missing_value=missing_value)

          id_cell_cld_sct_uv = register_diag_field    &
                         (mod_name, 'cell_cld_sct_uv', axes(1:3),&
                          Time, '.27um cell cloud sct coeff',   &
                          'km-1', missing_value=missing_value)

          id_cell_cld_sct_vis = register_diag_field    &
                         (mod_name, 'cell_cld_sct_vis', axes(1:3),  &
                          Time, '.55um cell cloud sct coeff',  &
                          'km-1', missing_value=missing_value)

          id_cell_cld_sct_nir = register_diag_field    &
                         (mod_name, 'cell_cld_sct_nir', axes(1:3), &
                          Time, '1.4um cell cloud sct coeff', &
                          'km-1', missing_value=missing_value)

          id_cell_cld_asymm_uv = register_diag_field    &
                         (mod_name, 'cell_cld_asymm_uv', axes(1:3), &
                          Time, '.27um cell cloud asymm coeff',   &
                          'percent', missing_value=missing_value)

          id_cell_cld_asymm_vis = register_diag_field     &
                         (mod_name, 'cell_cld_asymm_vis', axes(1:3), &
                          Time, '.55um cell cloud asymm coeff',   &
                          'percent', missing_value=missing_value)

          id_cell_cld_asymm_nir = register_diag_field    &
                         (mod_name, 'cell_cld_asymm_nir', axes(1:3), &
                          Time, '1.4um cell cloud asymm coeff',    &
                          'percent', missing_value=missing_value)

          id_abs_cell_cld_lw = register_diag_field    &
                         (mod_name, 'cell_abs_lw', axes(1:3), Time, &
                          'cell cloud abs coeff lw', &
                          'percent', missing_value=missing_value)

          id_abs_cell_cld_10u = register_diag_field    &
                         (mod_name, 'cell_abs_10u', axes(1:3), Time,   &
                          'cell cloud abs coeff 10um band', &
                          'percent', missing_value=missing_value)

!---------------------------------------------------------------------
!    register the donner meso-scale cloud radiative property diagnostic 
!    fields.
!---------------------------------------------------------------------
          id_meso_cld_amt = register_diag_field     &
                         (mod_name, 'meso_cld_amt', axes(1:3), Time,&
                          'meso cloud amount', 'percent',      &
                          missing_value=missing_value)

          id_meso_cld_ext_uv = register_diag_field    &
                         (mod_name, 'meso_cld_ext_uv', axes(1:3),&
                          Time, '.27um meso cloud ext coeff',   &
                          'km-1', missing_value=missing_value)

          id_meso_cld_ext_vis = register_diag_field   &
                         (mod_name, 'meso_cld_ext_vis', axes(1:3), &
                          Time, '.55um meso cloud ext coeff',  &
                          'km-1', missing_value=missing_value)

          id_meso_opdepth = register_diag_field   &
                         (mod_name, 'meso_opdepth', axes(1:3), &
                          Time, '.55um cloud optical depth avgd over &
                          &subcolumns with meso clouds',  &
                          'none', missing_value=missing_value)

          id_meso_cld_ext_nir = register_diag_field   &
                         (mod_name, 'meso_cld_ext_nir', axes(1:3), &
                          Time, '1.4um meso cloud ext coeff',  &
                          'km-1', missing_value=missing_value)

          id_meso_cld_sct_uv = register_diag_field   &
                         (mod_name, 'meso_cld_sct_uv', axes(1:3),&
                          Time, '.27um meso cloud sct coeff',  &
                          'km-1', missing_value=missing_value )

          id_meso_cld_sct_vis = register_diag_field  &
                         (mod_name, 'meso_cld_sct_vis', axes(1:3),  &
                          Time, '.55um meso cloud sct coeff',  &
                          'km-1', missing_value=missing_value)

          id_meso_cld_sct_nir = register_diag_field  &
                         (mod_name, 'meso_cld_sct_nir', axes(1:3),  &
                          Time, '1.4um meso cloud sct coeff',  &
                          'km-1', missing_value=missing_value)

          id_meso_cld_asymm_uv = register_diag_field  &
                         (mod_name, 'meso_cld_asymm_uv', axes(1:3),  &
                          Time, '.27um meso cloud asymm coeff',   &
                          'percent', missing_value=missing_value)

          id_meso_cld_asymm_vis = register_diag_field   &
                         (mod_name, 'meso_cld_asymm_vis', axes(1:3),   &
                          Time, '.55um meso cloud asymm coeff',    &
                          'percent', missing_value=missing_value)

          id_meso_cld_asymm_nir = register_diag_field    &
                         (mod_name, 'meso_cld_asymm_nir', axes(1:3), &
                          Time, '1.4um meso cloud asymm coeff',   &
                          'percent', missing_value=missing_value)

          id_abs_meso_cld_lw = register_diag_field    &
                         (mod_name, 'meso_abs_lw', axes(1:3),  &
                          Time, 'meso cloud abs coeff lw',  &
                          'percent', missing_value=missing_value)

          id_abs_meso_cld_10u = register_diag_field   &
                         (mod_name, 'meso_abs_10u', axes(1:3), Time,   &
                          'meso cloud abs coeff 10um band', &
                          'percent', missing_value=missing_value)
        endif
      else
        call error_mesg ('cloudrad_diagnostics_mod', &
            'Cldrad_control%do_donner_deep_clouds not yet defined',  &
                                                                FATAL)
      endif

!---------------------------------------------------------------------
!    register the uw shallow cloud radiative property diagnostic 
!    fields.
!---------------------------------------------------------------------
      if (Cldrad_control%do_uw_clouds_iz) then
        if (Cldrad_control%do_uw_clouds) then
          id_shallow_cld_amt = register_diag_field    &
                         (mod_name, 'shallow_cld_amt', axes(1:3),&
                          Time, 'shallow cloud amount',   &
                          'percent', missing_value=missing_value)

          id_shallow_cld_ext_uv = register_diag_field   &
                         (mod_name, 'uw_shallow_cld_ext_uv', &
                          axes(1:3), Time,   &
                          '.27um uw shallow cloud ext coeff', &
                          'km-1', missing_value=missing_value)

          id_shallow_cld_ext_vis = register_diag_field   &
                         (mod_name, 'uw_shallow_cld_ext_vis',&
                          axes(1:3), Time,   &
                          '.55um uw shallow cloud ext coeff',&
                          'km-1', missing_value=missing_value)

          id_shallow_opdepth = register_diag_field   &
                         (mod_name, 'uw_shallow_opdepth',&
                          axes(1:3), Time,   &
                          '.55um cloud optical depth avgd over &
                          &subcolumns with uw shallow clouds',&
                          'none', missing_value=missing_value)

          id_shallow_cld_ext_nir = register_diag_field   &
                         (mod_name, 'uw_shallow_cld_ext_nir',&
                          axes(1:3), Time,   &
                          '1.4um uw shallow cloud ext coeff',&
                          'km-1', missing_value=missing_value)

          id_shallow_cld_sct_uv = register_diag_field    &
                         (mod_name, 'uw_shallow_cld_sct_uv', &
                          axes(1:3), Time,  &
                          '.27um uw shallow cloud sct coeff', &
                          'km-1', missing_value=missing_value)

          id_shallow_cld_sct_vis = register_diag_field    &
                         (mod_name, 'uw_shallow_cld_sct_vis',&
                          axes(1:3), Time,    &
                          '.55um uw_shallow cloud sct coeff',&
                          'km-1', missing_value=missing_value)

          id_shallow_cld_sct_nir = register_diag_field    &
                         (mod_name, 'uw_shallow_cld_sct_nir', &
                          axes(1:3), Time,   &
                          '1.4um uw shallow cloud sct coeff',&
                          'km-1', missing_value=missing_value)

          id_shallow_cld_asymm_uv = register_diag_field    &
                         (mod_name, 'uw_shallow_cld_asymm_uv',   &
                          axes(1:3), Time,  &
                          '.27um uw shallow cloud asymm coeff',   &
                          'percent', missing_value=missing_value)

          id_shallow_cld_asymm_vis = register_diag_field     &
                         (mod_name, 'uw_shallow_cld_asymm_vis',   &
                          axes(1:3), Time,     &
                          '.55um uw shallow cloud asymm coeff',   &
                          'percent', missing_value=missing_value)

          id_shallow_cld_asymm_nir = register_diag_field    &
                         (mod_name, 'uw_shallow_cld_asymm_nir',    &
                          axes(1:3), Time,&
                          '1.4um uw shallow cloud asymm coeff',    &
                          'percent', missing_value=missing_value)

          id_abs_shallow_cld_lw = register_diag_field    &
                         (mod_name, 'uw_shallow_abs_lw', axes(1:3),  &
                          Time, 'uw shallow cloud abs coeff lw',   &
                          'percent', missing_value=missing_value)

          id_abs_shallow_cld_10u = register_diag_field    &
                         (mod_name, 'uw_shallow_abs_10u', axes(1:3), &
                          Time, 'uw shallow cloud abs coeff 10um band',&
                          'percent',  missing_value=missing_value)
        endif
      else
        call error_mesg ('cloudrad_diagnostics_mod', &
                'Cldrad_control%do_uw_clouds not yet defined', FATAL)
      endif

!--------------------------------------------------------------------
!    register total cloud condensate (non-stochastic case)   
!--------------------------------------------------------------------
      if ( .not. Cldrad_control%do_stochastic_clouds) then
        id_all_conc_drop = register_diag_field     &
                         (mod_name, 'all_conc_drop', axes(1:3), Time, &
                          'In-cloud liq water content - all clouds', &
                          'grams/m3', missing_value=missing_value)

        id_all_conc_ice = register_diag_field     &
                         (mod_name, 'all_conc_ice', axes(1:3), Time, &
                          'In-cloud ice water content - all clouds', &
                          'grams/m3', missing_value=missing_value)
      endif

!--------------------------------------------------------------------
!    register stratiform microphysical properties
!--------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds_iz) then
        if (Cldrad_control%do_strat_clouds) then
          id_strat_area_liq = register_diag_field     &
              (mod_name, 'strat_area_liq', axes(1:3), Time, &
               'Area of stratiform liquid clouds', 'fraction', &
               missing_value=missing_value)

          id_strat_conc_drop = register_diag_field     &
              (mod_name, 'strat_conc_drop', axes(1:3), Time, &
               'In-cloud liq water content of stratiform clouds', &
               'grams/m3', missing_value=missing_value,   &
               mask_variant = .true.)

          id_gb_strat_conc_drop = register_diag_field     &
              (mod_name, 'gb_strat_conc_drop', axes(1:3), Time, &
               'Grid-box-mean liq water content of stratiform clouds', &
               'grams/m3', missing_value=missing_value)      

          id_strat_size_drop = register_diag_field     &
              (mod_name, 'strat_size_drop', axes(1:3), Time, &
               'Effective diameter for stratiform liquid clouds', &
               'microns', missing_value=missing_value,   &
               mask_variant = .true.)

          id_strat_size_drop_l = register_diag_field     &
              (mod_name, 'strat_size_drop_l', axes(1:3), Time, &
               'lim Effective diameter for stratiform liquid clouds', &
               'microns', missing_value=missing_value,   &
               mask_variant = .true.)

          id_ra_strat_size_drop = register_diag_field     &
              (mod_name, 'ra_strat_size_drop', axes(1:3), Time, &
               'Adjusted effective diameter for strat liquid clouds', &
               'microns', missing_value=missing_value,   &
               mask_variant = .true.)

          id_strat_area_ice = register_diag_field     &
              (mod_name, 'strat_area_ice', axes(1:3), Time, &
               'Area of stratiform ice clouds', &
               'fraction', missing_value=missing_value)

          id_strat_conc_ice = register_diag_field     &
              (mod_name, 'strat_conc_ice', axes(1:3), Time, &
               'In-cloud ice water content of stratiform clouds', &
               'grams/m3', missing_value=missing_value,   &
               mask_variant = .true.)

          id_gb_strat_conc_ice = register_diag_field     &
              (mod_name, 'gb_strat_conc_ice', axes(1:3), Time, &
               'Grid-box-mean ice water content of stratiform clouds', &
               'grams/m3', missing_value=missing_value)   

          id_strat_size_ice = register_diag_field     &
              (mod_name, 'strat_size_ice', axes(1:3), Time, &
               'Effective diameter for stratiform ice clouds', &
               'microns', missing_value=missing_value,   &
               mask_variant = .true.)

          id_strat_size_ice_l = register_diag_field     &
              (mod_name, 'strat_size_ice_l', axes(1:3), Time, &
               'lim Effective diameter for stratiform ice clouds', &
               'microns', missing_value=missing_value,   &
               mask_variant = .true.)

          id_strat_droplet_number = register_diag_field     &
              (mod_name, 'strat_droplet_number', axes(1:3), Time, &
               'cloud droplet number for stratiform clouds', &
               '# per kg of air', missing_value=missing_value, &
               mask_variant = .true.)
     
          id_strat_droplet_number_l = register_diag_field     &
              (mod_name, 'strat_droplet_number_l', axes(1:3), Time, &
               'lim cloud droplet number for stratiform clouds', &
               '# per kg of air', missing_value=missing_value, &
               mask_variant = .true.)


          id_strat_ice_number = register_diag_field     &
              (mod_name, 'strat_ice_number', axes(1:3), Time, &
              'In-cloud Cloud ice crystal number for stratiform clouds', &
               '# per kg of air',    &
               missing_value=missing_value, mask_variant = .true.  )
 

          id_strat_ice_number_l = register_diag_field     &
             (mod_name, 'strat_ice_number_l', axes(1:3), Time, &
             'lim In-cloud Cloud ice crystal number for stratiform clouds' , &
            '# per kg of air',    &
            missing_value=missing_value, mask_variant = .true.  )

          id_lsc_lwp = register_diag_field     &
              (mod_name, 'strat_lwp', axes(1:2), Time, &
               'In-cloud liquid water path of stratiform clouds', &
               'kg/m2', missing_value=missing_value,   &
               mask_variant = .true.)

!---------------------------------------------------------------------
!RSH:
!    added to provide backward compatibility to existing diag_tables.
!    here 'LWP' and 'strat_lwp' are identical; given a choice, please 
!    use 'strat_lwp' in the diag_table.
!---------------------------------------------------------------------
          if (id_lsc_lwp <= 0) then
            id_lsc_lwp = register_diag_field     &
                (mod_name, 'LWP', axes(1:2), Time, &
                 'In-cloud liquid water path of stratiform clouds', &
                 'kg/m2', missing_value=missing_value,   &
                 mask_variant = .true.)
          endif

          id_gb_lsc_lwp = register_diag_field     &
              (mod_name, 'gb_strat_lwp', axes(1:2), Time, &
               'Grid-box-mean liquid water path of stratiform clouds', &
               'kg/m2', missing_value=missing_value)      

          id_lsc_iwp = register_diag_field     &
              (mod_name, 'strat_iwp', axes(1:2), Time, &
               'In-cloud ice water path of stratiform clouds', &
               'kg/m2', missing_value=missing_value,   &
                mask_variant = .true.)

          id_gb_lsc_iwp = register_diag_field     &
              (mod_name, 'gb_strat_iwp', axes(1:2), Time, &
               'Grid-box-mean ice water path of stratiform clouds', &
               'kg/m2', missing_value=missing_value)      
        
        endif
      else
        call error_mesg ('cloudrad_diagnostics_mod', &
               'Cldrad_control%do_strat_clouds not yet defined', FATAL)
      endif ! (do_strat_clouds_iz)

!--------------------------------------------------------------------
!    register donner meso cloud microphysical properties
!--------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds_iz) then
        if (Cldrad_control%do_donner_deep_clouds) then
          id_meso_area_liq = register_diag_field     &
              (mod_name, 'meso_area_liq', axes(1:3), Time, &
               'Area of donner meso liquid clouds', &
               'fraction', missing_value=missing_value)

          id_meso_conc_drop = register_diag_field     &
              (mod_name, 'meso_conc_drop', axes(1:3), Time, &
              'In-cloud liq water content of donner meso clouds', &
               'grams/m3', missing_value=missing_value,   &
               mask_variant = .true.)

          id_gb_meso_conc_drop = register_diag_field     &
              (mod_name, 'gb_meso_conc_drop', axes(1:3), Time, &
               'Grid-box-mean liq water content of donner meso clouds',&
               'grams/m3', missing_value=missing_value)   
       
          id_meso_size_drop = register_diag_field     &
              (mod_name, 'meso_size_drop', axes(1:3), Time, &
               'Effective diameter for donner meso liquid clouds', &
               'microns', missing_value=missing_value,   &
               mask_variant = .true.)

          id_ra_meso_size_drop = register_diag_field     &
              (mod_name, 'ra_meso_size_drop', axes(1:3), Time, &
               'Adjusted Effective diam for donner meso liq clouds',&
               'microns', missing_value=missing_value ,  &
               mask_variant = .true.)

          id_meso_area_ice = register_diag_field     &
              (mod_name, 'meso_area_ice', axes(1:3), Time, &
               'Area of donner meso ice clouds', 'fraction',    &
               missing_value=missing_value)

          id_meso_conc_ice = register_diag_field     &
              (mod_name, 'meso_conc_ice', axes(1:3), Time, &
               'In-cloud ice water content of donner meso clouds', &
               'grams/m3', missing_value=missing_value,   &
               mask_variant = .true.)

          id_gb_meso_conc_ice = register_diag_field     &
              (mod_name, 'gb_meso_conc_ice', axes(1:3), Time, &
               'Grid-box-mean ice water content of donner meso clouds',&
               'grams/m3', missing_value=missing_value)

          id_meso_size_ice = register_diag_field     &
              (mod_name, 'meso_size_ice', axes(1:3), Time, &
               'Effective diameter for donner meso ice clouds', &
               'microns', missing_value=missing_value,   &
               mask_variant = .true.)

          id_meso_droplet_number = register_diag_field     &
              (mod_name, 'meso_droplet_number', axes(1:3), Time, &
               'Cloud droplet number for donner meso clouds', &
               '# per kg of air', missing_value=missing_value,   &
               mask_variant = .true.)
     
          id_meso_lwp = register_diag_field     &
              (mod_name, 'meso_lwp', axes(1:2), Time, &
               'In-cloud liquid water path of donner meso clouds', &
               'kg/m2', missing_value=missing_value,  &
               mask_variant = .true.)

          id_gb_meso_lwp = register_diag_field     &
              (mod_name, 'gb_meso_lwp', axes(1:2), Time, &
               'Grid-box-mean liquid water path of donner meso clouds',&
               'kg/m2', missing_value=missing_value)   

          id_meso_iwp = register_diag_field     &
              (mod_name, 'meso_iwp', axes(1:2), Time, &
               'In-cloud ice water path of donner meso clouds', &
               'kg/m2',  missing_value=missing_value,   &
               mask_variant = .true.)

          id_gb_meso_iwp = register_diag_field     &
              (mod_name, 'gb_meso_iwp', axes(1:2), Time, &
               'Grid-box-mean ice water path of donner meso clouds', &
               'kg/m2', missing_value=missing_value)   

!--------------------------------------------------------------------
!    register donner cell microphysical properties
!--------------------------------------------------------------------
          id_cell_area_liq = register_diag_field     &
              (mod_name, 'cell_area_liq', axes(1:3), Time, &
               'Area of donner cell liquid clouds', 'fraction',    &
               missing_value=missing_value)

          id_cell_conc_drop = register_diag_field     &
              (mod_name, 'cell_conc_drop', axes(1:3), Time, &
               'In-cloud liquid water content of donner cell clouds', &
               'grams/m3', missing_value=missing_value,   &
               mask_variant = .true.)

          id_gb_cell_conc_drop = register_diag_field     &
              (mod_name, 'gb_cell_conc_drop', axes(1:3), Time, &
               'Grid-box-mean liq water content of donner cell clouds',&
               'grams/m3', missing_value=missing_value)   

          id_cell_size_drop = register_diag_field     &
              (mod_name, 'cell_size_drop', axes(1:3), Time, &
               'Effective diameter for donner cell liquid clouds', &
               'microns', missing_value=missing_value ,   &
               mask_variant = .true. )

          id_ra_cell_size_drop = register_diag_field     &
              (mod_name, 'ra_cell_size_drop', axes(1:3), Time, &
               'Adjusted Effective diam for donner cell liq clouds',&
               'microns', missing_value=missing_value,  &
               mask_variant = .true.)

          id_cell_area_ice = register_diag_field     &
              (mod_name, 'cell_area_ice', axes(1:3), Time, &
               'Area of donner cell ice clouds',  'fraction',    &
               missing_value=missing_value)

          id_cell_conc_ice = register_diag_field     &
              (mod_name, 'cell_conc_ice', axes(1:3), Time, &
               'In-cloud ice water content of donner cell clouds', &
               'grams/m3', missing_value=missing_value,   &
               mask_variant = .true.)

          id_gb_cell_conc_ice = register_diag_field     &
              (mod_name, 'gb_cell_conc_ice', axes(1:3), Time, &
               'Grid-box-mean ice water content of donner cell clouds',&
               'grams/m3', missing_value=missing_value)

          id_cell_size_ice = register_diag_field     &
              (mod_name, 'cell_size_ice', axes(1:3), Time, &
               'Effective diameter for donner cell ice clouds', &
               'microns', missing_value=missing_value ,  &
               mask_variant = .true. )

          id_cell_droplet_number = register_diag_field     &
              (mod_name, 'cell_droplet_number', axes(1:3), Time, &
               'Cloud droplet number for donner cell clouds', &
               '# per kg of air', missing_value=missing_value,   &
               mask_variant = .true.)
     
          id_cell_lwp = register_diag_field     &
              (mod_name, 'cell_lwp', axes(1:2), Time, &
               'In-cloud liquid water path of donner cell clouds', &
               'kg/m2', missing_value=missing_value,   &
               mask_variant = .true.)

          id_gb_cell_lwp = register_diag_field     &
              (mod_name, 'gb_cell_lwp', axes(1:2), Time, &
               'Grid-box-mean liquid water path of donner cell clouds',&
               'kg/m2', missing_value=missing_value)   

          id_cell_iwp = register_diag_field     &
              (mod_name, 'cell_iwp', axes(1:2), Time, &
               'In-cloud ice water path of donner cell clouds', &
               'kg/m2', missing_value=missing_value,   &
               mask_variant = .true.)

          id_gb_cell_iwp = register_diag_field     &
              (mod_name, 'gb_cell_iwp', axes(1:2), Time, &
               'Grid-box-mean ice water path of donner cell clouds', &
               'kg/m2', missing_value=missing_value)   
        endif
      else
        call error_mesg ('cloudrad_diagnostics_mod', &
          'Cldrad_control%do_donner_deep_clouds not yet defined', FATAL)
      endif  ! (do_donner_deep_clouds_iz)
        
!--------------------------------------------------------------------
!    register uw shallow microphysical properties
!--------------------------------------------------------------------
      if (Cldrad_control%do_uw_clouds_iz) then
        if (Cldrad_control%do_uw_clouds) then
          id_shallow_area_liq = register_diag_field     &
              (mod_name, 'shallow_area_liq', axes(1:3), Time, &
               'Area of uw shallow  liquid clouds', 'fraction',    &
                missing_value=missing_value)

          id_shallow_conc_drop = register_diag_field     &
              (mod_name, 'shallow_conc_drop', axes(1:3), Time, &
               'In-cloud liquid water content of uw shallow  clouds', &
               'grams/m3', missing_value=missing_value,   &
               mask_variant = .true.)

          id_gb_shallow_conc_drop = register_diag_field     &
              (mod_name, 'gb_shallow_conc_drop', axes(1:3), Time, &
               'Grid-box-mean liq water content of uw shallow clouds', &
               'grams/m3', missing_value=missing_value)   

          id_shallow_size_drop = register_diag_field     &
              (mod_name, 'shallow_size_drop', axes(1:3), Time, &
               'Effective diameter for uw shallow liquid clouds', &
               'microns', missing_value=missing_value ,  &
               mask_variant = .true.)

          id_ra_shallow_size_drop = register_diag_field     &
              (mod_name, 'ra_shallow_size_drop', axes(1:3), Time, &
               'Adjusted Effective diam for uw shallow  liq clouds',&
               'microns', missing_value=missing_value,   &
               mask_variant = .true.)

          id_shallow_area_ice = register_diag_field     &
              (mod_name, 'shallow_area_ice', axes(1:3), Time, &
               'Area of uw shallow  ice clouds', 'fraction',    &
               missing_value=missing_value)

          id_shallow_conc_ice = register_diag_field     &
              (mod_name, 'shallow_conc_ice', axes(1:3), Time, &
               'In-cloud ice water content of uw shallow  clouds', &
               'grams/m3', missing_value=missing_value,  &
               mask_variant = .true.)

          id_gb_shallow_conc_ice = register_diag_field     &
              (mod_name, 'gb_shallow_conc_ice', axes(1:3), Time, &
               'Grid-box-mean ice water content of uw shallow clouds', &
               'grams/m3', missing_value=missing_value)

          id_shallow_size_ice = register_diag_field     &
              (mod_name, 'shallow_size_ice', axes(1:3), Time, &
               'Effective diameter for uw shallow  ice clouds', &
               'microns', missing_value=missing_value ,  &
               mask_variant = .true. )

          id_shallow_droplet_number = register_diag_field     &
              (mod_name, 'shallow_droplet_number', axes(1:3), Time, &
               'Cloud droplet number for uw shallow  clouds', &
               '# per kg of air', missing_value=missing_value,   &
               mask_variant = .true.)
     
          id_shallow_lwp = register_diag_field     &
              (mod_name, 'shallow_lwp', axes(1:2), Time, &
               'In-cloud liquid water path of uw shallow  clouds', &
               'kg/m2', missing_value=missing_value,   &
               mask_variant = .true.)

          id_gb_shallow_lwp = register_diag_field     &
              (mod_name, 'gb_shallow_lwp', axes(1:2), Time, &
               'Grid-box-mean liquid water path of uw shallow clouds', &
               'kg/m2', missing_value=missing_value)   

          id_shallow_iwp = register_diag_field     &
              (mod_name, 'shallow_iwp', axes(1:2), Time, &
               'In-cloud ice water path of uw shallow  clouds', &
               'kg/m2', missing_value=missing_value,   &
               mask_variant = .true.)

          id_gb_shallow_iwp = register_diag_field     &
              (mod_name, 'gb_shallow_iwp', axes(1:2), Time, &
               'Grid-box-mean ice water path of uw shallow  clouds', &
               'kg/m2', missing_value=missing_value)   
        endif
      else
        call error_mesg ('cloudrad_diagnostics_mod', &
          'Cldrad_control%do_uw_clouds not yet defined', FATAL)
      endif  ! (do_uw_clouds_iz)

!---------------------------------------------------------------------
!
!    DIAGNOSTICS RELATED TO STOCHASTIC CLOUD PARAMETERIZATION
!
!---------------------------------------------------------------------
      if (Cldrad_control%do_stochastic_clouds_iz) then
        if (Cldrad_control%do_stochastic_clouds) then

!--------------------------------------------------------------------
!    the following diagnostics output the fields actually seen by the 
!    radiation code (either lsc or donner meso or donner cell or 
!    uw shallow in a given stochastic column, assuming all are 
!    activated), determined by the stochastic selection process.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    total cloud fraction summed across the stochastic subcolumns
!--------------------------------------------------------------------
          id_cldfrac_tot = register_diag_field  &
            (mod_name, 'stoch_cld_tot', axes(1:2), Time, &
             'total projected cloud fraction - stochastic clouds',&
             'fraction', missing_value=missing_value)

!--------------------------------------------------------------------
!     3d cloud fraction 
!--------------------------------------------------------------------
          id_cldfrac_ave = register_diag_field  &
            (mod_name, 'stoch_cld_ave', axes(1:3), Time, &
             'avg cloud fraction - stochastic clouds', &
             'fraction', missing_value=missing_value)

!--------------------------------------------------------------------
!     ice and water content and paths
!--------------------------------------------------------------------
          id_ice_conc_ave = register_diag_field  &
            (mod_name, 'stoch_ice_conc_ave', axes(1:3), Time, &
             'grid box avg ice water content - stochastic clouds', &
             'g/m3', missing_value=missing_value)

          id_ic_ice_conc_ave = register_diag_field  &
            (mod_name, 'stoch_incloud_ice_conc_ave', axes(1:3), Time, &
             'cloudy column avg ice water content - stochastic clouds',&
             'g/m3', missing_value=missing_value, mask_variant=.true.)

          id_drop_conc_ave = register_diag_field  &
            (mod_name, 'stoch_drop_conc_ave', axes(1:3), Time, &
             'grid box avg liquid water content - stochastic clouds', &
             'g/m3', missing_value=missing_value)

          id_ic_drop_conc_ave = register_diag_field  &
            (mod_name, 'stoch_incloud_drop_conc_ave', axes(1:3), Time, &
             'cloudy column avg liq water content - stochastic clouds',&
             'g/m3', missing_value=missing_value, mask_variant = .true.)

          id_lwp_ave = register_diag_field  &
            (mod_name, 'stoch_lwp_ave', axes(1:2), Time, &
             'grid box avg liq water path - stochastic clouds', &
             'kg/m2', missing_value=missing_value)

          id_ic_lwp_ave = register_diag_field  &
            (mod_name, 'stoch_incloud_lwp_ave', axes(1:2), Time, &
             'cloudy column avg liq water path - stochastic clouds', &
             'kg/m2', missing_value=missing_value, mask_variant=.true.)

          id_iwp_ave = register_diag_field  &
            (mod_name, 'stoch_iwp_ave', axes(1:2), Time, &
             'grid box avg ice water path - stochastic clouds', &
             'kg/m2', missing_value=missing_value)

          id_ic_iwp_ave = register_diag_field  &
            (mod_name, 'stoch_incloud_iwp_ave', axes(1:2), Time, &
             'cloudy column avg ice water path - stochastic clouds', &
             'kg/m2', missing_value=missing_value, mask_variant=.true.)

          id_LWPr = register_diag_field    &
             (mod_name, 'LWPr', axes(1:2), Time, &  
              'LWPr', 'kg m-2 micron-1')

!--------------------------------------------------------------------
!     ice and water cloud frequency distribution
!--------------------------------------------------------------------
          id_ice_col_frac_ave = register_diag_field  &
            (mod_name, 'stoch_ice_col_frac_ave', axes(1:3), Time, &
             'frctn of columns with ice clouds - stochastic clouds', &
             'fraction', missing_value=missing_value,   &
             mask_variant = .true.)

          id_liq_col_frac_ave = register_diag_field  &
            (mod_name, 'stoch_liq_col_frac_ave', axes(1:3), Time, &
             'frctn of columns with liq clouds - stochastic clouds', &
             'fraction', missing_value=missing_value, &
             mask_variant = .true.)

!--------------------------------------------------------------------
!     crystal and droplet sizes and number
!--------------------------------------------------------------------
          id_ice_size_ave = register_diag_field  &
            (mod_name, 'stoch_ice_size_ave', axes(1:3), Time, &
             'cloudy column avg ice eff diam - stochastic clouds', &
             'microns', missing_value=missing_value, &
             mask_variant = .true.)

          id_drop_size_ave = register_diag_field  &
            (mod_name, 'stoch_drop_size_ave', axes(1:3), Time, &
             'cloudy column avg droplet diam - stochastic clouds', &
             'microns', missing_value=missing_value,  &
             mask_variant = .true.)

          id_ra_drop_size_ave = register_diag_field  &
            (mod_name, 'ra_stoch_drop_size_ave', axes(1:3), Time, &
             'adjustd cloudy column avg drop diam - stochastic clouds',&
             'microns', missing_value=missing_value,  &
             mask_variant = .true. )

          id_ice_size_ave_l = register_diag_field  &
            (mod_name, 'stoch_ice_size_ave_l', axes(1:3), Time, &
             'lim cloudy column avg ice eff diam - stochastic clouds', &
             'microns', missing_value=missing_value, &
             mask_variant = .true.)

          id_drop_size_ave_l = register_diag_field  &
            (mod_name, 'stoch_drop_size_ave_l', axes(1:3), Time, &
             'lim cloudy column avg droplet diam - stochastic clouds', &
             'microns', missing_value=missing_value,  &
             mask_variant = .true.)



          id_droplet_number_ave_l = register_diag_field  &
            (mod_name, 'stoch_droplet_number_ave_l', axes(1:3), Time, &
            'lim cloudy column avg droplet number - stochastic clouds', &
            '#/kg(air)', missing_value=missing_value,  &
            mask_variant = .true. )

          id_ice_number_ave_l = register_diag_field  &
            (mod_name, 'stoch_ice_number_ave_l', axes(1:3), Time, &
              'lim cloudy column avg ice number - stochastic clouds', &
             '#/kg(air)', missing_value=missing_value,  &
             mask_variant = .true. )

          id_droplet_number_ave = register_diag_field  &
            (mod_name, 'stoch_droplet_number_ave', axes(1:3), Time, &
             'cloudy column avg droplet number - stochastic clouds', &
             '#/kg(air)', missing_value=missing_value,  &
             mask_variant = .true. )
      
!--------------------------------------------------------------------
!    the following fields are the contributions from each of the 
!    cloud types to the fields actually seen by the radiation package.
!--------------------------------------------------------------------
          id_lsc_drop_conc_ave = register_diag_field  &
            (mod_name, 'stoch_lsc_drop_conc_ave', axes(1:3), Time, &
             'grid box avg lsc liq water content - stochastic clouds', &
             'g/m3', missing_value=missing_value)

          id_lsc_ice_conc_ave = register_diag_field  &
            (mod_name, 'stoch_lsc_ice_conc_ave', axes(1:3), Time, &
             'grid box avg lsc ice water content - stochastic clouds', &
             'g/m3', missing_value=missing_value)

          id_lsc_lwp_ave = register_diag_field  &
            (mod_name, 'stoch_lsc_lwp_ave', axes(1:2), Time, &
             'grid box avg lsc liq water path  - stochastic clouds', &
             'kg/m2', missing_value=missing_value)

          id_lsc_iwp_ave = register_diag_field  &
            (mod_name, 'stoch_lsc_iwp_ave', axes(1:2), Time, &
             'grid box avg lsc ice water path - stochastic clouds', &
             'kg/m2', missing_value=missing_value)

          id_cell_drop_conc_ave = register_diag_field  &
            (mod_name, 'stoch_cell_drop_conc_ave', axes(1:3), Time, &
             'grid box avg cell liq water content - stochastic clouds',&
             'g/m3', missing_value=missing_value)

          id_cell_ice_conc_ave = register_diag_field  &
            (mod_name, 'stoch_cell_ice_conc_ave', axes(1:3), Time, &
             'grid box avg cell ice water content - stochastic clouds',&
             'g/m3', missing_value=missing_value)

          id_cell_lwp_ave = register_diag_field  &
            (mod_name, 'stoch_cell_lwp_ave', axes(1:2), Time, &
             'grid box avg cell liq water path - stochastic clouds', &
             'kg/m2', missing_value=missing_value)

          id_cell_iwp_ave = register_diag_field  &
            (mod_name, 'stoch_cell_iwp_ave', axes(1:2), Time, &
             'grid box avg cell ice water path - stochastic clouds', &
             'kg/m2', missing_value=missing_value)

          id_meso_drop_conc_ave = register_diag_field  &
            (mod_name, 'stoch_meso_drop_conc_ave', axes(1:3), Time, &
             'grid box avg meso liq water content - stochastic clouds',&
             'g/m3', missing_value=missing_value)

          id_meso_ice_conc_ave = register_diag_field  &
            (mod_name, 'stoch_meso_ice_conc_ave', axes(1:3), Time, &
             'grid box avg meso ice water content - stochastic clouds',&
             'g/m3', missing_value=missing_value)

          id_meso_lwp_ave = register_diag_field  &
            (mod_name, 'stoch_meso_lwp_ave', axes(1:2), Time, &
             'grid box avg meso liq water path - stochastic clouds', &
             'kg/m2', missing_value=missing_value)

          id_meso_iwp_ave = register_diag_field  &
            (mod_name, 'stoch_meso_iwp_ave', axes(1:2), Time, &
             'grid box avg meso ice water path - stochastic clouds', &
             'kg/m2', missing_value=missing_value)

          id_shallow_drop_conc_ave = register_diag_field  &
            (mod_name, 'stoch_shallow_drop_conc_ave', axes(1:3), Time, &
             'grid box avg shallow liq water content - stoch clouds', &
             'g/m3', missing_value=missing_value)

          id_shallow_ice_conc_ave = register_diag_field  &
            (mod_name, 'stoch_shallow_ice_conc_ave', axes(1:3), Time, &
             'grid box avg shallow ice water content - stoch clouds', &
             'g/m3', missing_value=missing_value)

          id_shallow_lwp_ave = register_diag_field  &
            (mod_name, 'stoch_shallow_lwp_ave', axes(1:2), Time, &
             'grid box avg shallow liq water path - stochastic clouds',&
             'kg/m2', missing_value=missing_value)

          id_shallow_iwp_ave = register_diag_field  &
            (mod_name, 'stoch_shallow_iwp_ave', axes(1:2), Time, &
             'grid box avg shallow ice water path - stochastic clouds',&
             'kg/m2', missing_value=missing_value)

!--------------------------------------------------------------------
!    the following fields wouild be obtained if the stochastic cloud
!    treatment was limited to just the strat_cloud component. the
!    difference between these '_only_lsc' variables and the '_ave' 
!    variables just above reflects the effect of treating the other
!    cloud types stochastically.
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!     3d cloud fraction 
!--------------------------------------------------------------------
          id_cldfrac_only_lsc = register_diag_field  &
            (mod_name, 'stoch_cld_only_lsc', axes(1:3), Time, &
             'avg cld fraction, only lsc stochastic', &
             'fraction', missing_value=missing_value)

!--------------------------------------------------------------------
!     ice and water content and paths.
!--------------------------------------------------------------------
          id_ice_conc_only_lsc = register_diag_field  &
            (mod_name, 'stoch_ice_conc_only_lsc', axes(1:3), Time, &
             'grid box avg ice water content, only lsc stochastic', &
             'g/m3', missing_value=missing_value)

          id_ic_ice_conc_only_lsc = register_diag_field  &
            (mod_name, 'stoch_incloud_ice_conc_only_lsc', axes(1:3), &
             Time, &
             'cldy column avg ice water content, only lsc stochastic',&
             'g/m3', missing_value=missing_value, mask_variant = .true.)

          id_drop_conc_only_lsc = register_diag_field  &
            (mod_name, 'stoch_drop_conc_only_lsc', axes(1:3), Time, &
             'grid box avg liq water content, only lsc stochastic', &
             'g/m3', missing_value=missing_value)

          id_ic_drop_conc_only_lsc = register_diag_field  &
            (mod_name, 'stoch_incloud_drop_conc_only_lsc', axes(1:3), &
             Time, &
             'cldy column avg liq water content, only lsc stochastic', &
             'g/m3', missing_value=missing_value, mask_variant = .true.)

          id_lwp_only_lsc = register_diag_field  &
            (mod_name, 'stoch_lwp_only_lsc', axes(1:2), Time, &
             'grid box avg liq water path, only lsc stochastic', &
             'kg/m2', missing_value=missing_value)

          id_ic_lwp_only_lsc = register_diag_field  &
            (mod_name, 'stoch_incloud_lwp_only_lsc', axes(1:2), Time, &
             'in-cloud avg liq water path, only lsc stochastic', &
             'kg/m2', missing_value=missing_value,   &
             mask_variant = .true.)

          id_iwp_only_lsc = register_diag_field  &
            (mod_name, 'stoch_iwp_only_lsc', axes(1:2), Time, &
             'grid box avg ice water path, only lsc stochastic', &
             'kg/m2', missing_value=missing_value)

          id_ic_iwp_only_lsc = register_diag_field  &
            (mod_name, 'stoch_incloud_iwp_only_lsc', axes(1:2), Time, &
             'in-cloud avg ice water path, only lsc stochastic', &
             'kg/m2', missing_value=missing_value,   &
             mask_variant = .true.)

!--------------------------------------------------------------------
!     ice and water cloud distributions
!--------------------------------------------------------------------
          id_ice_col_only_lsc = register_diag_field  &
            (mod_name, 'stoch_ice_col_frac_only_lsc', axes(1:3), Time, &
             'frctn of columns with ice clouds, only lsc stochastic',&
             'fraction', missing_value=missing_value,   &
             mask_variant=.true.)

          id_liq_col_only_lsc = register_diag_field  &
            (mod_name, 'stoch_liq_col_frac_only_lsc', axes(1:3), Time, &
             'frctn of columns with liq clouds, only lsc stochastic',&
             'fraction', missing_value=missing_value, &
             mask_variant = .true.)

!--------------------------------------------------------------------
!     crystal and droplet sizes and number
!--------------------------------------------------------------------
          id_ice_size_only_lsc = register_diag_field  &
            (mod_name, 'stoch_ice_size_only_lsc', axes(1:3), Time, &
             'cloudy column avg ice eff diam, only lsc stochastic', &
             'microns', missing_value=missing_value, &
             mask_variant = .true.)

          id_drop_size_only_lsc = register_diag_field  &
            (mod_name, 'stoch_drop_size_only_lsc', axes(1:3), Time, &
             'cloudy column avg drop diam, only lsc stochastic', &
             'microns', missing_value=missing_value,  &
             mask_variant = .true.)

          id_ra_drop_size_only_lsc = register_diag_field  &
            (mod_name, 'ra_stoch_drop_size_only_lsc', axes(1:3), Time, &
             'adjustd cldy column avg drop diam, only lsc stochastic', &
            'microns', missing_value=missing_value,  &
            mask_variant = .true.)

          id_droplet_number_only_lsc = register_diag_field  &
            (mod_name, 'stoch_droplet_number_only_lsc', axes(1:3),   &
             Time,&
             'liq cldy column avg droplet number, only lsc stochastic',&
             '#/kg(air)', missing_value=missing_value,  &
             mask_variant = .true.)

!--------------------------------------------------------------------
!    diagnostics relative to the frequency that the radiation code sees
!    lsc cloud properties.
!--------------------------------------------------------------------
          id_stoch_ic_lsc_cf_ave = register_diag_field  &
            (mod_name, 'stoch_ic_lsc_cf_ave', axes(1:3), Time, &
             'fractn of cols in cloudy grid boxes with lsc props', &
             'fraction', missing_value=missing_value,  &
             mask_variant = .true.)


          id_stoch_sees_lsc = register_diag_field  &
            (mod_name, 'stoch_sees_lsc', axes(1:3), Time, &
             'fraction of times lsc clds are used when present',&
             'fraction', missing_value=missing_value,  &
             mask_variant = .true.)
       
          id_stoch_lsc_cf_ave = register_diag_field  &
            (mod_name, 'stoch_lsc_cf_ave', axes(1:3), Time, &
             'fraction of stochastic columns assigned lsc props',&
             'fraction', missing_value=missing_value)

!--------------------------------------------------------------------
!    the following fields are only valid if the donner parameterization
!    is active.
!--------------------------------------------------------------------
          if (Cldrad_control%do_donner_deep_clouds_iz) then
            if (Cldrad_control%do_donner_deep_clouds) then

!--------------------------------------------------------------------
!    diagnostics relative to the frequency that the radiation code sees
!    donner meso and cell cloud properties.
!--------------------------------------------------------------------
              id_stoch_ic_cell_cf_ave = register_diag_field  &
                (mod_name, 'stoch_ic_cell_cf_ave', axes(1:3), Time, &
                'fractn of cols in cloudy grid boxes with cell props', &
                'fraction', missing_value=missing_value,  &
                mask_variant = .true.)

              id_stoch_ic_meso_cf_ave = register_diag_field  &
                (mod_name, 'stoch_ic_meso_cf_ave', axes(1:3), Time, &
                'fractn of cols in cloudy grid boxes with meso props', &
                'fraction', missing_value=missing_value,   &
                mask_variant = .true.)

              id_stoch_sees_cell = register_diag_field  &
                (mod_name, 'stoch_sees_cell', axes(1:3), Time, &
                 'fraction of times cell clds are seen when present', &
                 'fraction', missing_value=missing_value,  &
                 mask_variant = .true.)

              id_stoch_sees_meso = register_diag_field  &
                (mod_name, 'stoch_sees_meso', axes(1:3), Time, &
                 'fraction of times meso clds are seen when present', &
                 'fraction', missing_value=missing_value,   &
                 mask_variant = .true.)

              id_stoch_cell_cf_ave = register_diag_field  &
                (mod_name, 'stoch_cell_cf_ave', axes(1:3), Time, &
                 'fraction of stochastic columns assigned cell props',&
                 'fraction', missing_value=missing_value)

              id_stoch_meso_cf_ave = register_diag_field  &
                (mod_name, 'stoch_meso_cf_ave', axes(1:3), Time, &
                 ' fraction of stochastic columns assigned meso props',&
                 'fraction', missing_value=missing_value)
            endif
          else
            call error_mesg ('cloudrad_diagnostics_mod', &
              'Cldrad_control%do_donner_deep_clouds not yet defined', &
                                                                 FATAL)
          endif

!--------------------------------------------------------------------
!    the following fields are only valid if the uw shallow parameter-
!    ization is active.
!--------------------------------------------------------------------
          if (Cldrad_control%do_uw_clouds_iz) then
            if (Cldrad_control%do_uw_clouds) then
!--------------------------------------------------------------------
!    diagnostics indicating frequency that radiation code sees lsc, 
!    meso, cell and uw shallow clouds, and the frequency that these 
!    cloud types are seen when they exist.
!--------------------------------------------------------------------
              id_stoch_ic_shallow_cf_ave = register_diag_field  &
                (mod_name, 'stoch_ic_shallow_cf_ave', axes(1:3), Time, &
                 'frctn of cols in cldy grid boxes with shallow props',&
                 'fraction', missing_value=missing_value,  &
                 mask_variant = .true.)

              id_stoch_sees_shallow = register_diag_field  &
                (mod_name, 'stoch_sees_shallow', axes(1:3), Time, &
                 'frctn of times shallow clds are seen when present', &
                 'fraction', missing_value=missing_value,  &
                 mask_variant = .true.)

              id_stoch_shallow_cf_ave = register_diag_field  &
                (mod_name, 'stoch_shallow_cf_ave', axes(1:3), Time, &
                 'frctn of stoch columns assigned uw shallow props',&
                 'fraction', missing_value=missing_value)
            endif
          else
            call error_mesg ('cloudrad_diagnostics_mod', &
               'Cldrad_control%do_uw_clouds not yet defined', FATAL)
          endif

!---------------------------------------------------------------------
!    cloud property diagnostics for each stochastic sub-column
!---------------------------------------------------------------------
          do n=1, ncol                                          
            if (n < 10) then
              write (chvers,'(i1)') n 
            else if (n <100) then
              write (chvers,'(i2)') n 
            else
              call error_mesg ('cloudrad_diagnostics_mod', &
                'must modify code to allow writing of more than&
                 & 99 stochastic columns', FATAL)
            endif

!---------------------------------------------------------------------
!    cloud type diagnostic : 0 = no cloud, 1 = lsc, 2 = meso, 3 = cell
!                            4 = uw shallow
!----------------------------------------------------------------------
            id_stoch_cloud_type(n) = register_diag_field  &
                (mod_name, 'stoch_cloud_type_'//trim(chvers),  &
                 axes(1:3), Time, &
                 'cloud type (1-4) in stochastic col  '//trim(chvers), &
                 'none', missing_value=missing_value,   &
                 mask_variant = .true.)
 
!--------------------------------------------------------------------
!
!
!              STOCHASTIC COLUMN VALUES, "_ONLY_LSC" DIAGNOSTICS
!
!
!--------------------------------------------------------------------
!    these "_only_lsc" diagnostics allow one to assess the effect of
!    treating the non-lsc clouds stochastically. the difference between
!    the "_only_lsc" variables and the corresponding "_ave" variables
!    reflect the changes resulting from treating the non-lsc cloud
!    types stochastically. the cloud properties actually seen by the
!    model's radiation package are always contained in the "_ave" 
!    diagnostics.
!--------------------------------------------------------------------
            id_cldfrac_cols_only_lsc(n) = register_diag_field  &
              (mod_name, 'stoch_cld_col_only_lsc_'//trim(chvers),  &
               axes(1:3), Time, &
               'cloud fraction in column ' //trim(chvers) // &
                     ' when only lsc clouds treated stochastically ', &
               'fraction', missing_value=missing_value)
          
!----------------------------------------------------------------------
!    lsc ice and liquid water content and paths, and  droplet number 
!----------------------------------------------------------------------
            id_ice_conc_cols_only_lsc(n) = register_diag_field  &
              (mod_name, 'stoch_ice_conc_col_only_lsc_'//trim(chvers), &
               axes(1:3), Time, &
               'ice content in column '//trim(chvers)  // &
                    ' when only lsc clouds treated stochastically', &
               'g/m3', missing_value=missing_value)

            id_drop_conc_cols_only_lsc(n) = register_diag_field  &
              (mod_name, 'stoch_drop_conc_col_only_lsc_'//trim(chvers),&
               axes(1:3), Time, &
               'liq water content in stoch column '//trim(chvers) // &
                      ' when only lsc clouds treated stochastically', &
               'g/m3', missing_value=missing_value)

            id_droplet_number_cols_only_lsc(n) = register_diag_field  &
              (mod_name,    &
               'stoch_droplet_number_col_only_lsc_'//trim(chvers),&
               axes(1:3), Time, &
               ' droplet number in stochastic column '//trim(chvers) //&
                  ' when only lsc clouds are treated stochastically', &
               '#/kg(air)', missing_value=missing_value,  &
               mask_variant = .true.)

            id_iwp_cols_only_lsc(n) = register_diag_field  &
              (mod_name, 'stoch_iwp_col_only_lsc_'//trim(chvers),  &
               axes(1:2), Time, &
               'ice water path in stochastic column '//trim(chvers) // &
                     ' when only lsc clouds treated stochastically',&
               'kg/m2', missing_value=missing_value)

            id_lwp_cols_only_lsc(n) = register_diag_field  &
              (mod_name, 'stoch_lwp_col_only_lsc_'//trim(chvers),  &
               axes(1:2), Time, &
               'liq water path in stochastic column '//trim(chvers) // &
                    ' when only lsc clouds treated stochastically',&
               'kg/m2', missing_value=missing_value)      
  
!----------------------------------------------------------------------
!   lsc ice particle size, liquid droplet size from microphysics and as
!   adjusted for use by radiative routines 
!----------------------------------------------------------------------
            id_ice_size_cols_only_lsc(n) = register_diag_field  &
              (mod_name, 'stoch_ice_size_col_only_lsc_'//trim(chvers), &
               axes(1:3), Time, &
               'ice eff diam in stochastic column '//trim(chvers) // &
                     ' when only lsc clouds treated stochastically', &
               'microns', missing_value=missing_value, &
               mask_variant = .true.)

            id_drop_size_cols_only_lsc(n) = register_diag_field  &
              (mod_name, 'stoch_drop_size_col_only_lsc_'//trim(chvers),&
               axes(1:3), Time, &
               'droplet diam in stochastic column '//trim(chvers) // &
                     ' when only lsc clouds treated stochastically', &
               'microns', missing_value=missing_value,  &
               mask_variant = .true.)

            id_ra_drop_size_cols_only_lsc(n) = register_diag_field  &
              (mod_name,   &
               'ra_stoch_drop_size_col_only_lsc_'//trim(chvers),  &
               axes(1:3), Time, &
               'adjustd drop diam in stoch column '//trim(chvers)// &
                      ' when only lsc clouds treated stochastically', &
               'microns', missing_value=missing_value,  &
               mask_variant = .true.)
 
!--------------------------------------------------------------------
!
!
!           STOCHASTIC COLUMN VALUES, AS SEEN BY RADIATION PACKAGE
!
!
!--------------------------------------------------------------------
!    the following diagnostics show the cloud properties seen by the
!    radiation package for each of the stochastic columns. 
!--------------------------------------------------------------------

!----------------------------------------------------------------------
!    cloud fraction in column that is seen by radiation code 
!----------------------------------------------------------------------
            id_cldfrac_cols(n) = register_diag_field  &
              (mod_name, 'stoch_cld_col_'//trim(chvers),   &
               axes(1:3), Time,&
               'cloud fraction in stochastic column '//trim(chvers), &
               'fraction', missing_value=missing_value)
          
!----------------------------------------------------------------------
!    ice and liquid water content and paths, and  droplet number 
!    as seen by radiation code.
!----------------------------------------------------------------------
            id_ice_conc_cols(n) = register_diag_field  &
              (mod_name, 'stoch_ice_conc_col_'//trim(chvers),  &
               axes(1:3), Time, &
               'ice content in stochastic column '//trim(chvers), &
               'g/m3', missing_value=missing_value)

            id_drop_conc_cols(n) = register_diag_field  &
              (mod_name, 'stoch_drop_conc_col_'//trim(chvers), &
               axes(1:3), Time, &
               'water content in stochastic column '//trim(chvers), &
               'g/m3', missing_value=missing_value)
     
            id_iwp_cols(n) = register_diag_field  &
              (mod_name, 'stoch_iwp_col_'//trim(chvers), &
               axes(1:2), Time, &
               'ice water path in stochastic column '//trim(chvers), &
               'kg/m2', missing_value=missing_value)

            id_lwp_cols(n) = register_diag_field  &
              (mod_name, 'stoch_lwp_col_'//trim(chvers),  &
               axes(1:2), Time, &
               'liq water path in stochastic column '//trim(chvers), &
               'kg/m2', missing_value=missing_value)

            id_droplet_number_cols(n) = register_diag_field  &
              (mod_name, 'stoch_droplet_number_col_'//trim(chvers), &
               axes(1:3), Time, &
               'droplet number in stochastic column '//trim(chvers), &
               '#/kg(air)', missing_value=missing_value,   &
               mask_variant = .true.)

!----------------------------------------------------------------------
!    ice particle size, liquid droplet size from microphysics and as
!    adjusted for use by radiative routines 
!    (as seen by radiation code)
!----------------------------------------------------------------------
            id_ice_size_cols(n) = register_diag_field  &
              (mod_name, 'stoch_ice_size_col_'//trim(chvers), &
               axes(1:3), Time, &
               'ice particle eff diam in stoch column '//trim(chvers), &
               'microns', missing_value=missing_value,   &
               mask_variant = .true.)

            id_drop_size_cols(n) = register_diag_field  &
              (mod_name, 'stoch_drop_size_col_'//trim(chvers), &
               axes(1:3), Time, &
               'droplet diameter in stochastic column '//trim(chvers), &
               'microns', missing_value=missing_value,  &
               mask_variant = .true.)

            id_ra_drop_size_cols(n) = register_diag_field  &
              (mod_name, 'ra_stoch_drop_size_col_'//trim(chvers), &
               axes(1:3), Time, &
               'adjustd droplet diam in stoch column '//trim(chvers), &
               'microns', missing_value=missing_value,  &
               mask_variant = .true.)

          end do
        endif ! (do_stochastic_clouds)
      else
        call error_mesg ('cloudrad_diagnostics_mod', &
          'Cldrad_control%do_stochastic_clouds not yet defined', FATAL)
      endif ! (do_stochastic_clouds_iz)

!---------------------------------------------------------------------

 
end subroutine diag_field_init



!#####################################################################
! <SUBROUTINE NAME="isccp_diag">
!  <OVERVIEW>
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_diag (is, js, Cld_spec, Atmos_input, coszen, Time)
!  </TEMPLATE>
!  <IN NAME="is,js" TYPE="integer">
!   starting subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   time on next timestep, used as stamp for 
!                        diagnostic output [ time_type (days, seconds) ]
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    atmospheric input fields on model grid,
!  </IN>
!  <IN NAME="coszen" TYPE="real">
!    cosine of solar zenith angle
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </IN>
!  <IN NAME="Lsctau" TYPE="real">
!   cloud optical thickness in the visible
!  </IN>
!  <IN NAME="Lsclwem" TYPE="real">
!   10 micron cloud emissivity
!  </IN>
! </SUBROUTINE>
!
subroutine isccp_diag (is, js, Cld_spec, Atmos_input, coszen,       &
                       Lsctau, Lsclwem, Time)

!--------------------------------------------------------------------
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!---------------------------------------------------------------------
 
integer,                      intent(in)   :: is,js
type(cld_specification_type), intent(in)   :: Cld_spec
type(atmos_input_type),       intent(in)   :: Atmos_input
real, dimension(:,:),         intent(in)   :: coszen
real, dimension(:,:,:),       intent(in)   :: Lsctau, Lsclwem
type(time_type),              intent(in)   :: Time

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,js           starting/ending subdomain i,j indices of data 
!                      in the physics_window being integrated
!      Cld_spec        cloud specification properties on model grid,
!                      [ cld_specification_type ]
!      Atmos_input     atmospheric input fields on model grid,
!                      [ atmos_input_type ] 
!      coszen          cosine of zenith angle [ dimensionless ]
!      Lsctau          0.6-0.685 micron cloud optical thickness 
!                      [ dimensionless ]
!      Lsclwem         Longwave emissivity [ dimensionless ]
!                      This is from 10.1-11.1 microns if the multiband 
!                      longwave emissivity is active.
!      Time            time on next timestep, used as stamp for 
!                      diagnostic output [ time_type (days, seconds) ]
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2), &
                       size(Cld_spec%lwp,3)) :: &
                                    tau_local, em_local, cldamt_local

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2), &
                       size(Cld_spec%lwp,3) ) ::  qv, em_lw_local




      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2), &
                       size(Cld_spec%lwp,3), 4 ) ::  tau

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2), &
                      7, 7) ::       fq_isccp

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2)) :: &
                          npoints, ninhomog, inhomogeneity_parameter

      integer      :: kdim
      integer      :: max_cld
      integer      :: i, j, k
      integer, dimension(size(Cld_spec%lwp,1),size(Cld_spec%lwp,2)):: &
                        sunlit
      
!---------------------------------------------------------------------
!   local variables:
!
!      tau_local        optical depth in band 1 in the current column
!                       [ dimensionless ]
!      em_local         lw cloud emissivity in current column
!                       [ dimensionless ]
!      cldamt_local     cloud fraction in current column 
!                       [ dimensionless ]
!      qv               water vapor specific humidity
!                       [ kg vapor / kg air ]
!      em_lw_local      lw cloud emissivity [ dimensionless ]
!      rh2o             mixing ratio of water vapor at model full levels
!                       [ non-dimensional ]
!      temp             temperature at model levels (1:nlev), surface
!                       temperature is stored at value nlev+1; if eta
!                       coordinates, surface value stored in below 
!                       ground points
!                       [ deg K ]
!      tau              optical depth in 4 bands [ dimensionless ]
!      fq_isccp         matrix of fractional area covered by cloud
!                       types of a given optical depth and cloud
!                       top pressure range.  The matrix is 7x7 for
!                       7 cloud optical depths and 7 cloud top 
!                       pressure ranges
!      npoints          flag indicating whether isccp cloud is present
!                       in column (cloud + daylight needed)
!      kdim             number of model layers
!      max_cld          greatest number of clouds in any column in the
!                       current physics window
!      i,j,k            do-loop indices
!     
!      sunlit           is the given i,j point sunlit?
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!    define number of model layers.
!----------------------------------------------------------------------
      kdim = size (Cld_spec%lwp,3)
        
!---------------------------------------------------------------------
!    If optical properties are needed and no clouds exist in the 
!    window, call cloud_optical_properties_diag to define the cloud 
!    optical depth, the optical depth due to cloud ice and the 
!    longwave emissivity. If no clouds exist in the window, all the 
!    optical depths and emissivities are left are their initial
!    values of zero.to zero.
!---------------------------------------------------------------------
     
     em_lw_local = 0.
     tau = 0.

     max_cld = MAXVAL(Cld_spec%ncldsw(:,:))
     if (max_cld >= 1 .and. .not.isccp_actual_radprops)          &
       call cloud_optical_properties_diag (Cld_spec, tau, em_lw_local)

!---------------------------------------------------------------------
!    Initialize fields
!---------------------------------------------------------------------

           npoints(:,:) = 0.
           fq_isccp(:,:,:,:) = 0.
           ninhomog(:,:) = 0.
           inhomogeneity_parameter(:,:) = 0.

!---------------------------------------------------------------------
!    Compute sunlit integer flag
!---------------------------------------------------------------------
           sunlit(:,:) = 0
           where(coszen(:,:) > 1.E-06) sunlit(:,:) = 1

!--------------------------------------------------------------------
!    define the specific humidity from the mixing ratio which has been
!    input.
!--------------------------------------------------------------------
                  qv(:,:,:) = Atmos_input%cloudvapor(:,:,:)/   &
                              (1. + Atmos_input%cloudvapor(:,:,:))
                
!---------------------------------------------------------------------
!    define the column values of cloud fraction, cloud optical depth, 
!    and lw cloud emissivity. if cloud is not present, set these var-
!    iables to clear sky values.
!---------------------------------------------------------------------
           do j=1,size(Cld_spec%lwp,2)
            do i=1,size(Cld_spec%lwp,1)
             do k=1,kdim                      

                  if (Cld_spec%camtsw(i,j,k) > 0.0) then
                    cldamt_local(i,j,k) = Cld_spec%camtsw(i,j,k) 
                    if (isccp_actual_radprops) then
                      tau_local(i,j,k) = Lsctau(i,j,k)
                      em_local(i,j,k) = Lsclwem(i,j,k)
                    else 
                      tau_local(i,j,k) = tau(i,j,k,1)/ &
                                 real(Cld_spec%cld_thickness(i,j,k))
                      em_local(i,j,k) = 1.-((1.-em_lw_local(i,j,k))** &
                             (1./real(Cld_spec%cld_thickness(i,j,k))) )
                    end if          
                  else
                    cldamt_local(i,j,k) = 0.
                    tau_local(i,j,k) = 0.
                    em_local(i,j,k) = 0.
                  endif
                  
                end do
               end do
              end do
               
!---------------------------------------------------------------------
!    call isccp_cloudtypes to map each model cloud to an isccp cloud
!    type, based on its optical depth and height above the surface.
!    set a flag to indicate the presence of isccp cloud in this column.
!---------------------------------------------------------------------
                call isccp_cloudtypes (sunlit(:,:), &
                                       Atmos_input%press(:,:,1:kdim), &
                                       Atmos_input%pflux(:,:,:),&
                                       qv(:,:,:),       &
                                       Atmos_input%cloudtemp(:,:,:),  &
                                       Atmos_input%temp(:,:,kdim+1),  &
                                       cldamt_local, tau_local,   &
                                       em_local, fq_isccp(:,:,:,:),   &
                                       npoints(:,:), &
                                       inhomogeneity_parameter(:,:), &
                                       ninhomog(:,:))
 
!----------------------------------------------------------------------
!    send any desired diagnostics to the diag_manager_mod.
!----------------------------------------------------------------------
              
                
          call isccp_output (is, js, fq_isccp, npoints, &
                             inhomogeneity_parameter, ninhomog, Time)
   
!---------------------------------------------------------------------
    
    
end subroutine isccp_diag        

!#####################################################################
! <SUBROUTINE NAME="isccp_diag_stochastic">
!  <OVERVIEW>
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_diag (is, js, Cld_spec, Atmos_input, coszen, Time)
!  </TEMPLATE>
!  <IN NAME="is,js" TYPE="integer">
!   starting subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   time on next timestep, used as stamp for 
!                        diagnostic output [ time_type (days, seconds) ]
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    atmospheric input fields on model grid,
!  </IN>
!  <IN NAME="coszen" TYPE="real">
!    cosine of solar zenith angle
!  </IN>
!  <IN NAME="Lsctau" TYPE="real">
!   cloud optical thickness in the visible
!  </IN>
!  <IN NAME="Lsclwem" TYPE="real">
!   10 micron cloud emissivity
!  </IN>
! </SUBROUTINE>
!
subroutine isccp_diag_stochastic (is, js, Atmos_input, coszen,       &
                                  Lsctau, Lsclwem, LscCldAmt, Time)

!--------------------------------------------------------------------
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!---------------------------------------------------------------------
 
integer,                     intent(in)   :: is,js
type(atmos_input_type),      intent(in)   :: Atmos_input
real, dimension(:,:),        intent(in)   :: coszen
real, dimension(:,:,:,:),    intent(in)   :: Lsctau, Lsclwem, LscCldAmt
type(time_type),             intent(in)   :: Time

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,js           starting/ending subdomain i,j indices of data 
!                      in the physics_window being integrated
!      Atmos_input     atmospheric input fields on model grid,
!                      [ atmos_input_type ] 
!      coszen          cosine of zenith angle [ dimensionless ]
!      Lsctau          0.6-0.685 micron cloud optical thickness 
!                      [ dimensionless ]
!      Lsclwem         Longwave emissivity [ dimensionless ]
!                      This is from 10.1-11.1 microns if the multiband 
!                      longwave emissivity is active.
!      LsCldAmt        Cloud fraction [ dimensionless ]
!                      Values should be identically 0 or 1. 
!      Time            time on next timestep, used as stamp for 
!                      diagnostic output [ time_type (days, seconds) ]
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      real, dimension (size(Lsctau,1), size(Lsctau,2), size(Lsctau,3) ) ::  qv

      !
      ! Isccp histogram variables
      !
      real, dimension (size(Lsctau,1), size(Lsctau,2), 7, 7) &
                                                       :: fq_isccp

      real, dimension (size(Lsctau,1), size(Lsctau,2)) :: npoints, &
                       ninhomog, inhomogeneity_parameter
      

      integer      :: kdim

      integer, dimension(size(Lsctau,1),size(Lsctau,2)):: sunlit
      
!---------------------------------------------------------------------
!   local variables:
!
!      qv               water vapor specific humidity
!                       [ kg vapor / kg air ]
!      fq_isccp         matrix of fractional area covered by cloud
!                       types of a given optical depth and cloud
!                       top pressure range.  The matrix is 7x7 for
!                       7 cloud optical depths and 7 cloud top 
!                       pressure ranges
!      npoints          flag indicating whether isccp cloud is present
!                       in column (cloud + daylight needed)
!      kdim             number of model layers
!      max_cld          greatest number of clouds in any column in the
!                       current physics window
!      i,j,k            do-loop indices
!     
!      sunlit           is the given i,j point sunlit?
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!    define number of model layers.
!----------------------------------------------------------------------
      kdim = size (Lsctau,3)
        
!---------------------------------------------------------------------
!    If optical properties are needed and no clouds exist in the 
!    window, call cloud_optical_properties_diag to define the cloud 
!    optical depth, the optical depth due to cloud ice and the 
!    longwave emissivity. If no clouds exist in the window, all the 
!    optical depths and emissivities are left are their initial
!    values of zero.to zero.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    Initialize ISCCP histograms
!---------------------------------------------------------------------

     npoints(:,:) = 0.
     fq_isccp(:,:,:,:) = 0.
     ninhomog(:,:) = 0.
     inhomogeneity_parameter(:,:) = 0.

!---------------------------------------------------------------------
!    Compute sunlit integer flag
!---------------------------------------------------------------------
     sunlit(:,:) = 0
     where(coszen(:,:) > 1.E-06) sunlit(:,:) = 1

!--------------------------------------------------------------------
!    define the specific humidity from the mixing ratio which has been
!    input.
!--------------------------------------------------------------------
     qv(:,:,:) = Atmos_input%cloudvapor(:,:,:)/ (1. + Atmos_input%cloudvapor(:,:,:))
               
!---------------------------------------------------------------------
!    call isccp_cloudtypes to map each model cloud to an isccp cloud
!    type, based on its optical depth and height above the surface.
!    set a flag to indicate the presence of isccp cloud in this column.
!---------------------------------------------------------------------
     call isccp_cloudtypes_stochastic (sunlit(:,:),        &
                            Atmos_input%press(:,:,1:kdim), &
                            Atmos_input%pflux(:,:,:),      &
                            qv(:,:,:),                     &
                            Atmos_input%cloudtemp(:,:,:),  &
                            Atmos_input%temp(:,:,kdim+1),  &
                            LscCldAmt(:, :, :, :),         &
                            LscTau(:, :, :, :),            &
                            Lsclwem(:, :, :, :),           &
                            fq_isccp(:,:,:,:), npoints(:,:), &
                            inhomogeneity_parameter(:,:), &
                            ninhomog(:,:))
                                                
!----------------------------------------------------------------------
!    send any desired diagnostics to the diag_manager_mod.
!----------------------------------------------------------------------
     call isccp_output (is, js, fq_isccp, npoints, &
                             inhomogeneity_parameter, ninhomog, Time)
   
!---------------------------------------------------------------------
    
    
end subroutine isccp_diag_stochastic        


!#####################################################################
! <SUBROUTINE NAME="compute_isccp_clds">
!  <OVERVIEW>
!    subroutine compute_isccp_clds maps the model clouds into isccp
!    categories (high, middle, low) and defines the cloud fraction of
!    each.
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine compute_isccp_clds maps the model clouds into isccp
!    categories (high, middle, low) and defines the cloud fraction of
!    each.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call compute_isccp_clds (pflux, camtsw, hml_ca)
!  </TEMPLATE>
!  <IN NAME="pflux" TYPE="real">
!   average of pressure at adjacent model levels
!  </IN>
!  <IN NAME="camtsw" TYPE="real">
!   total cloud amount [ nondimensional ]
!  </IN>
!  <OUT NAME="hml_ca" TYPE="real">
!   cloud fraction for the 3 isccp cloud types
!  </OUT>
! </SUBROUTINE>
!
subroutine compute_isccp_clds (pflux, camtsw, camtsw_band, hml_ca)

!---------------------------------------------------------------------
!    subroutine compute_isccp_clds maps the model clouds into isccp
!    categories (high, middle, low) and defines the cloud fraction of
!    each.
!--------------------------------------------------------------------- 

real,  dimension(:,:,:),   intent(in)  :: pflux, camtsw
real,  dimension(:,:,:,:), intent(in)  :: camtsw_band
real,  dimension(:,:,:),   intent(out) :: hml_ca

!---------------------------------------------------------------------
!  intent(in) variables:
!
!        pflux           average of pressure at adjacent model levels
!                        [ (kg /( m s^2) ]
!        camtsw          total cloud amount [ nondimensional ]
!        camtsw_band     total cloud amount in each sw band 
!                        [ nondimensional ]
!
!  intent(out) variable:
!
!        hml_ca          cloud fraction for the 3 isccp cloud types
!                        [ nondimensional ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      integer          ::   i, j, k

!---------------------------------------------------------------------
!  local variables:
!
!         mid_btm     pressure boundary between isccp middle and 
!                     isccp low clouds [ Pa ]
!         high_btm    pressure boundary between isccp middle and
!                     isccp high clouds [ Pa ]
!         i,j,k       do-loop indices
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    initialize a column integrated hi-mid-low cloud-free area array.
!--------------------------------------------------------------------
      hml_ca = 1.0
 
!---------------------------------------------------------------------
!    define arrays giving the cloud-free area in each column within the
!    pressure regionscorresponding to the ISCCP definitions of high 
!    (10-440 hPa), middle (440-680 hPa) and low (680-1000 hPa) clouds. 
!    compute high, middle and low cloud amounts assuming that independ-
!    ent clouds overlap randomly. note that model clouds above 10 hPa 
!    and below 1000 hPa are included in the totals.
!----------------------------------------------------------------------
      do k = 1,size(pflux,3)-1
        do j=1,size(pflux,2)
          do i=1,size(pflux,1)
            if (pflux(i,j,k)  <=  high_btm) then
              hml_ca(i,j,1) = hml_ca(i,j,1) * (1. - camtsw(i,j,k))
            else if ( (pflux(i,j,k) >  high_btm) .and.  &
                      (pflux(i,j,k) <=  mid_btm) ) then
              hml_ca(i,j,2) = hml_ca(i,j,2) * (1. - camtsw(i,j,k))
            else  if ( pflux(i,j,k) > mid_btm ) then
              hml_ca(i,j,3) = hml_ca(i,j,3) * (1. - camtsw(i,j,k))
            endif
            if (pflux(i,j,k)  >  high_btm) then
              hml_ca(i,j,4) = hml_ca(i,j,4) * (1. - camtsw(i,j,k))
            endif
          end do
        end do
      end do

!--------------------------------------------------------------------
!    convert the cloud-free area to an integrated cloud fraction in 
!    the column. express the cloud area in percent.
!--------------------------------------------------------------------
      hml_ca = 1. - hml_ca
      hml_ca = 100.*hml_ca
  
!-------------------------------------------------------------------


end subroutine compute_isccp_clds



!####################################################################
! <SUBROUTINE NAME="cloud_optical_properties_diag">
!  <OVERVIEW>
!    cloud_optical_properties_diag calculates the cloud optical depth,
!    ice cloud optical depth and longwave cloud emissivity for each
!    cloudy grid box.
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloud_optical_properties_diag calculates the cloud optical depth,
!    ice cloud optical depth and longwave cloud emissivity for each
!    cloudy grid box.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_optical_properties_diag (Cld_spec, tau, em_lw)
!  </TEMPLATE>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid
!  </IN>
!  <OUT NAME="tau" TYPE="real">
!   cloud optical depth in each of the
!                     num_slingo_bands
!  </OUT>
!  <OUT NAME="em_lw" TYPE="real">
!   longwave cloud emissivity
!  </OUT>
! </SUBROUTINE>
!
subroutine cloud_optical_properties_diag (Cld_spec, tau, em_lw)

!---------------------------------------------------------------------
!    cloud_optical_properties_diag calculates the cloud optical depth,
!    ice cloud optical depth and longwave cloud emissivity for each
!    cloudy grid box.
!---------------------------------------------------------------------
                              
type(cld_specification_type), intent(in)   :: Cld_spec
real, dimension(:,:,:,:),     intent(out)  :: tau
real, dimension(:,:,:),       intent(out)  :: em_lw       

!--------------------------------------------------------------------
!   intent(in) variable:
!
!      Cld_spec       cloud specification properties on model grid,
!                     [ cld_specification_type ]
!
!   intent(out) variables:
!
!      tau            cloud optical depth in each of the
!                     num_slingo_bands [ dimensionless ]
!      em_lw          longwave cloud emissivity [ dimensionless ]  
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:


      real, dimension (size(Cld_spec%lwp,1),size(Cld_spec%lwp,2),  &
                       size(Cld_spec%lwp,3), 4) ::     &
                                                tau_liq, tau_ice

      real, dimension (size(Cld_spec%lwp,1),size(Cld_spec%lwp,2),   &
                       size(Cld_spec%lwp,3)) ::     &
                                                k_liq, k_ice, &
                                                rdrop, rice

!--------------------------------------------------------------------
!   local variables:
!
!       tau_liq    liquid cloud optical depth [ dimensionless ]
!       tau_ice    ice    cloud optical depth [ dimensionless ]
!       k_liq      liquid cloud mass absorption coefficient for longwave
!                  portion of the spectrum [ m**2 / kg condensate ]
!       k_ice      ice cloud mass absorption coefficient for longwave
!                  portion of the spectrum [ m**2 / kg condensate ]
!       rdrop      droplet radius, forced to be within the valid range 
!                  for the slingo parameterization (4.2 < rdrop < 16.6)
!                  [ microns ]
!       rice       ice particle size, forced to be within the valid 
!                  range for the ebert and curry parameterization 
!                  (10.0 < rdrop < 130.0)
!                  [ microns ]
!
!---------------------------------------------------------------------

!    NOTE: THESE SIZE LIMITS ARE INDEPENDENT OF WHAT IS USED IN THE GCM.
!          This subroutine is called only if isccp_actual_radprops =
!          .false., indicating that the actual gcm properties are not
!          being used in the isccp processing.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    compute uv cloud optical depths due to liquid. the formula for 
!    optical depth comes from: 
!    Slingo (1989), J. Atmos. Sci., vol. 46, pp. 1419-1427
!---------------------------------------------------------------------
      rdrop(:,:,:) = MAX (Cld_spec%reff_liq(:,:,:), 4.2)
      rdrop(:,:,:) = MIN (rdrop(:,:,:), 16.6)
!---------------------------------------------------------------------
!    in this program, reff_ice is limited to be between 10 microns
!    and 130 microns, which is the range of validity for the Ebert
!    and Curry (1992) radiation.
!---------------------------------------------------------------------
      rice (:,:,:) = MAX (Cld_spec%reff_ice(:,:,:), 10.0)
      rice (:,:,:) = MIN (rice (:,:,:), 130.0)
      tau_liq(:,:,:,1) = Cld_spec%lwp(:,:,:) * 1000. * &
                         (0.02817 + (1.305/rdrop(:,:,:)))
      tau_liq(:,:,:,2) = Cld_spec%lwp(:,:,:) * 1000. * &
                         (0.02682 + (1.346/rdrop(:,:,:)))
      tau_liq(:,:,:,3) = Cld_spec%lwp(:,:,:) * 1000. * &
                         (0.02264 + (1.454/rdrop(:,:,:)))
      tau_liq(:,:,:,4) = Cld_spec%lwp(:,:,:) * 1000. * &
                         (0.01281 + (1.641/rdrop(:,:,:)))
        
!---------------------------------------------------------------------
!    compute uv cloud optical depths due to ice. the formula for 
!    optical depth comes from:
!    Ebert and Curry (1992), J. Geophys. Res., vol. 97, pp. 3831-3836
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    IMPORTANT:  WE ARE CHEATING HERE BECAUSE WE ARE FORCING THE FIVE 
!                BAND MODEL OF EBERT AND CURRY INTO THE FOUR BAND MODEL 
!                OF SLINGO. THIS IS DONE BY COMBINING BANDS 3 and 4 OF 
!                EBERT AND CURRY. EVEN SO THE EXACT BAND LIMITS DO NOT 
!                MATCH.  FOR COMPLETENESS HERE ARE THE BAND LIMITS IN 
!                MICRONS:
!            BAND               SLINGO                 EBERT AND CURRY
!
!             1               0.25-0.69                0.25 - 0.7
!             2               0.69-1.19                0.7 - 1.3
!             3               1.19-2.38                1.3 - 2.5
!             4               2.38-4.00                2.5 - 3.5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!---------------------------------------------------------------------
      tau_ice(:,:,:,1) = Cld_spec%iwp(:,:,:) * 1000. * &
                         (0.003448 + (2.431/rice(:,:,:)))
      tau_ice(:,:,:,2) = tau_ice(:,:,:,1)
      tau_ice(:,:,:,3) = tau_ice(:,:,:,1)
      tau_ice(:,:,:,4) = tau_ice(:,:,:,1)

!---------------------------------------------------------------------
!     back out scaling factor

      tau_liq = tau_liq / isccp_scale_factor
      tau_ice = tau_ice / isccp_scale_factor
        
!---------------------------------------------------------------------
!    compute total cloud optical depth. the mixed phase optical prop-
!    erties are based upon equation 14 of Rockel et al. 1991, 
!    Contributions to Atmospheric Physics, volume 64, pp.1-12. 
!    thus:

!          tau = tau_liq + tau_ice
!
!    and place a minimum value on tau - taumin
!---------------------------------------------------------------------
      tau(:,:,:,:) = max(tau_liq(:,:,:,:) + tau_ice(:,:,:,:),taumin)
        
!----------------------------------------------------------------------
!    define the  mass absorption coefficient for longwave radiation 
!    for cloud ice and cloud liquid.
!
!    NOTE THAT THE NUMBERS HERE ALREADY INCLUDE THE DIFFUSIVITY
!    FACTOR!
!----------------------------------------------------------------------
      k_liq(:,:,:) = 140.
      k_ice(:,:,:) = 4.83591 + 1758.511/rice(:,:,:)
        
!----------------------------------------------------------------------
!    compute combined lw emmisivity. the mixed phase optical properties
!    are based upon equation 14 of Rockel et al. 1991, Contributions to
!    Atmospheric Physics,  volume 64, pp.1-12.  thus:
!
!    transmivvity_lw =   transmissivity_lw_ice * transmissivity_lw_liq
!
!    which can also be written as:
!
!    em_lw =  em_lw_liq + em_lw_ice -  (em_lw_liq * em_lw_ice )
!
!    which is what is solved here.
!----------------------------------------------------------------------
      em_lw(:,:,:) =  1. - exp(-1.*(k_liq(:,:,:)*Cld_spec%lwp(:,:,:) + &
                                    k_ice(:,:,:)*Cld_spec%iwp(:,:,:))/ &
                                    isccp_scale_factor)

!---------------------------------------------------------------------


end subroutine cloud_optical_properties_diag


!######################################################################



                end module cloudrad_diagnostics_mod


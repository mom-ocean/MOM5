 
               module rad_utilities_mod
!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  smf
! </REVIEWER>
! 
! <OVERVIEW>
!  Code to define the derived data types and provide table search routines.
! </OVERVIEW>
! <DESCRIPTION>
!  This code is used in the radiation code package as a helper module.
!  It defines many derived data types used in radiation calculation.
!  This code also provides table search routines and simple arithmatic
!  routines.
! </DESCRIPTION>
!
use mpp_mod,            only : input_nml_file
use fms_mod,            only : open_namelist_file, fms_init, &
                               mpp_pe, mpp_root_pe, stdlog, &
                               file_exist, write_version_number, &
                               check_nml_error, error_mesg, &
                               FATAL, close_file, lowercase
use  field_manager_mod, only : parse

use time_manager_mod,   only : time_type

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    rad_utilities_mod contains radiation table search routines,
!    some band averaging routines, and the derived-type variables 
!    used in the radiation package.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

character(len=128)  :: version =  '$Id: rad_utilities.F90,v 20.0 2013/12/13 23:20:23 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'

!---------------------------------------------------------------------
!-------  interfaces --------

public      &
          rad_utilities_init, check_derived_types, & 
          locate_in_table,  &
          looktab, table_alloc, &
          thickavg, thinavg,   &
          rad_utilities_end,   &
          get_radiative_param, assignment(=)

interface looktab
    module procedure  looktab_type1, looktab_type2, looktab_type3
end interface

interface table_alloc
   module procedure    table1_alloc, table2_alloc, table3_alloc
end interface

interface thickavg
   module procedure thickavg_3d
   module procedure thickavg_0d
   module procedure thickavg_1band
   module procedure thickavg_isccp
end interface

interface assignment(=)
  module procedure lw_output_type_eq
  module procedure sw_output_type_eq
  module procedure aerosol_props_type_eq
end interface

!------------------------------------------------------------------

public   aerosol_type              
  
!    aerosol
!    aerosol_names
!    family_members

type aerosol_type
     real, dimension(:,:,:,:),        pointer :: aerosol=>NULL()
     logical, dimension(:,:),        pointer :: family_members=>NULL()
     character(len=64), dimension(:), pointer :: aerosol_names=>NULL()
end type aerosol_type              

!------------------------------------------------------------------
 
public   aerosol_diagnostics_type
 
!    extopdep
!    absopdep
!    extopdep_vlcno
!    absopdep_vlcno
!    lw_extopdep_vlcno
!    lw_absopdep_vlcno
!    sw_heating_vlcno
 
type aerosol_diagnostics_type
!    real, dimension(:,:,:),   pointer  :: sw_heating_vlcno=>NULL()
     real, dimension(:,:,:,:),   pointer  :: sw_heating_vlcno=>NULL()
     real, dimension(:,:,:,:,:), pointer  :: extopdep=>NULL(), &
                                             absopdep=>NULL()
     real, dimension(:,:,:,:,:), pointer  :: asymdep=>NULL()

     real, dimension(:,:,:,:), pointer  :: extopdep_vlcno=>NULL(), &
                                           absopdep_vlcno=>NULL(), &
                                           lw_extopdep_vlcno=>NULL(), &
                                           lw_absopdep_vlcno=>NULL()
 
end type aerosol_diagnostics_type
 
!------------------------------------------------------------------

public   aerosol_properties_type

!    aerextband
!    aerssalbband
!    aerasymmband
!    aerextbandlw
!    aerssalbbandlw
!    sulfate_index
!    optical_index

type aerosol_properties_type
     integer, dimension(:,:,:), pointer  :: ivol=>NULL()
     real, dimension(:,:), pointer  :: aerextband=>NULL(),   &
                                       aerssalbband=>NULL(), &
                                       aerasymmband=>NULL(), &
                                       aerextbandlw=>NULL(), &
                                       aerssalbbandlw=>NULL(), &
                                       aerextbandlw_cn=>NULL(), &
                                       aerssalbbandlw_cn=>NULL()
     real, dimension(:,:,:,:), pointer :: sw_ext=>NULL(), &
                                          sw_ssa=>NULL(), &
                                          sw_asy=>NULL(), &
                                          lw_ext=>NULL(), &
                                          lw_ssa=>NULL(), &
                                          lw_asy=>NULL()
!yim
     integer, dimension(:,:), pointer :: sulfate_index=>NULL()
     integer, dimension(:), pointer :: optical_index=>NULL()
     integer, dimension(:), pointer :: omphilic_index=>NULL()
     integer, dimension(:), pointer :: bcphilic_index=>NULL()
     integer, dimension(:), pointer :: seasalt1_index=>NULL()
     integer, dimension(:), pointer :: seasalt2_index=>NULL()
     integer, dimension(:), pointer :: seasalt3_index=>NULL()
     integer, dimension(:), pointer :: seasalt4_index=>NULL()
     integer, dimension(:), pointer :: seasalt5_index=>NULL()
     integer                        :: sulfate_flag
     integer                        :: omphilic_flag
     integer                        :: bcphilic_flag
     integer                        :: seasalt1_flag
     integer                        :: seasalt2_flag
     integer                        :: seasalt3_flag
     integer                        :: seasalt4_flag
     integer                        :: seasalt5_flag
!yim
     integer                        :: bc_flag
end type aerosol_properties_type

!------------------------------------------------------------------

public astronomy_type
    
!    solar
!    cosz
!    fracday
!    rrsun

type astronomy_type
     real, dimension(:,:), pointer  :: solar=>NULL(),   &
                                       cosz=>NULL(),  &
                                       fracday=>NULL()
     real, dimension(:,:,:), pointer  :: solar_p=>NULL(),   &
                                       cosz_p=>NULL(),  &
                                       fracday_p=>NULL()
     real    :: rrsun
end type astronomy_type

!--------------------------------------------------------------------
 
public astronomy_inp_type
 
!    zenith_angle   specified zenith angles [ degrees ]
!    fracday        specified daylight fraction [ fraction ]
!    rrsun          specified earth-sun distance, normalized by mean
!                   distance 

type astronomy_inp_type
     real, dimension(:,:), pointer  :: zenith_angle=>NULL()
     real, dimension(:,:), pointer  :: fracday=>NULL()
     real                           :: rrsun
end type astronomy_inp_type

!--------------------------------------------------------------------

public atmos_input_type

!    press
!    temp
!    rh2o
!    zfull
!    pflux
!    tflux
!    deltaz
!    phalf
!    rel_hum
!    cloudtemp
!    clouddeltaz
!    cloudvapor
!    aerosoltemp
!    aerosolvapor
!    aerosolpress
!    aerosolrelhum
!    tracer_co2
!    g_rrvco2
!    tsfc
!    psfc

type atmos_input_type
     real, dimension(:,:,:), pointer :: press=>NULL(),   &
                                        temp=>NULL(), &
                                        rh2o=>NULL(),  &
                                        zfull=>NULL(),  &
                                        pflux=>NULL(), &
                                        tflux=>NULL(),  &
                                        deltaz=>NULL(),  &
                                        phalf=>NULL(),   &
                                        rel_hum=>NULL(), &
                                        cloudtemp=>NULL(),   &
                                        clouddeltaz=>NULL(), &
                                        cloudvapor=>NULL(), &
                                        aerosoltemp=>NULL(), &
                                        aerosolvapor=>NULL(), &
                                        aerosolpress=>NULL(), &
                                        aerosolrelhum=>NULL(), &
                                        tracer_co2 => NULL()
     real, dimension(:,:),   pointer :: tsfc=>NULL(),   &
                                        psfc=>NULL()              
     real                            :: g_rrvco2
end type atmos_input_type

!-------------------------------------------------------------------

public cldrad_properties_type

!    cldext
!    cldasymm
!    cldsct
!    emmxolw
!    emrndlw
!    abscoeff
!    cldemiss
!    cirabsw
!    cirrfsw
!    cvisrfsw

type cldrad_properties_type
     real, dimension(:,:,:,:,:), pointer :: cldext=>NULL(),   &
                                          cldasymm=>NULL(), &
                                          cldsct=>NULL()
     real, dimension(:,:,:,:,:), pointer ::                   &
                                          emmxolw=>NULL(), &
                                          emrndlw=>NULL(),  &
                                          abscoeff=>NULL(),  &
                                          cldemiss=>NULL()
     real, dimension(:,:,:), pointer ::   cirabsw=>NULL(), &
                                          cirrfsw=>NULL(),   &
                                          cvisrfsw=>NULL()
end type cldrad_properties_type

!------------------------------------------------------------------

public cld_space_properties_type

!    camtswkc
!    cirabswkc
!    cirrfswkc
!    cvisrfswkc
!    ktopswkc
!    kbtmswkc

type cld_space_properties_type
     real, dimension(:,:,:),    pointer :: camtswkc=>NULL()        
     real, dimension(:,:,:),    pointer :: cirabswkc=>NULL(),  &
                                           cirrfswkc=>NULL(), &
                                           cvisrfswkc=>NULL()
     integer, dimension(:,:,:), pointer :: ktopswkc=>NULL(),   &
                                           kbtmswkc=>NULL()
end type cld_space_properties_type

!------------------------------------------------------------------

public cld_specification_type

!    tau
!    lwp
!    iwp
!    reff_liq
!    reff_ice
!    liq_frac
!    cloud_water
!    cloud_ice
!    cloud_area
!    cloud_droplet
!    reff_liq_micro
!    reff_ice_micro
!    camtsw
!    cmxolw
!    crndlw
!    cld_thickness
!    ncldsw
!    nmxolw
!    nrndlw
!    hi_cloud
!    mid_cloud
!    low_cloud
!    ice_cloud

type cld_specification_type
   real, dimension(:,:,:,:),  pointer :: tau=>NULL(),  &
                                         camtsw_band=>NULL(), &
                                         crndlw_band=>NULL(), &
                                         lwp_lw_band=>NULL(), &
                                         iwp_lw_band=>NULL(), &
                                         lwp_sw_band=>NULL(), &
                                         iwp_sw_band=>NULL(), &
                                         reff_liq_lw_band=>NULL(),   &
                                         reff_ice_lw_band=>NULL(), &
                                         reff_liq_sw_band=>NULL(),   &
                                         reff_ice_sw_band=>NULL()
   real, dimension(:,:,:),    pointer :: lwp=>NULL(),   &
                                         iwp=>NULL(),  &
                                         reff_liq=>NULL(),   &
                                         reff_ice=>NULL(), &
                                         reff_liq_lim=>NULL(),   &
                                         reff_ice_lim=>NULL(), &
                                         liq_frac=>NULL(), &
                                         cloud_water=>NULL(), &
                                         cloud_ice=>NULL(),  &
                                         cloud_area=>NULL(), &
                                         cloud_droplet=>NULL(), &
                                         cloud_ice_num=>NULL(), &
  ! snow, rain
                                         rain =>NULL(), &
                                         snow =>NULL(), &
                                         rain_size =>NULL(), &
                                         snow_size =>NULL(), &

                                         reff_liq_micro=>NULL(),   &
                                         reff_ice_micro=>NULL(),&
                                         camtsw=>NULL(),   &
                                         cmxolw=>NULL(),  &
                                         crndlw=>NULL()
   integer, dimension(:,:,:), pointer :: cld_thickness=>NULL()
   integer, dimension(:,:,:,:), pointer :: stoch_cloud_type=>NULL()
   integer, dimension(:,:,:,:), pointer :: cld_thickness_lw_band=>NULL()
   integer, dimension(:,:,:,:), pointer :: cld_thickness_sw_band=>NULL()
   integer, dimension(:,:),   pointer :: ncldsw=>NULL(),   &
                                         nmxolw=>NULL(),&
                                         nrndlw=>NULL()
   integer, dimension(:,:,:), pointer :: ncldsw_band=>NULL(),   &
                                         nrndlw_band=>NULL()
   logical, dimension(:,:,:), pointer :: hi_cloud=>NULL(),   &
                                         mid_cloud=>NULL(),  &
                                         low_cloud=>NULL(),   &
                                         ice_cloud=>NULL()
end type cld_specification_type

!------------------------------------------------------------------

public cloudrad_control_type

type cloudrad_control_type
    logical :: do_pred_cld_microphys
    logical :: do_presc_cld_microphys
    logical :: do_bulk_microphys
    logical :: do_sw_micro
    logical :: do_lw_micro
    logical :: do_rh_clouds        
    logical :: do_strat_clouds        
    logical :: do_zonal_clouds        
    logical :: do_mgroup_prescribed
    logical :: do_obs_clouds        
    logical :: do_no_clouds        
    logical :: do_diag_clouds        
    logical :: do_specified_clouds        
    logical :: do_donner_deep_clouds
    logical :: do_uw_clouds
    logical :: do_zetac_clouds
    logical :: do_random_overlap
    logical :: do_max_random_overlap
    logical :: do_stochastic_clouds
    logical :: use_temp_for_seed
    logical :: do_specified_strat_clouds
    logical :: do_ica_calcs
    logical :: do_liq_num
    logical :: do_ice_num
    logical :: using_fu2007
    integer :: nlwcldb                   !   number of frequency bands 
                                         !   for which lw cloud emissiv-
                                         !   ities are defined.
    integer :: cloud_data_points
    integer :: ich
    integer :: icm
    integer :: ict
    integer :: icb
    logical :: do_pred_cld_microphys_iz
    logical :: do_presc_cld_microphys_iz
    logical :: do_bulk_microphys_iz
    logical :: do_sw_micro_iz
    logical :: do_lw_micro_iz
    logical :: do_rh_clouds_iz
    logical :: do_strat_clouds_iz
    logical :: do_zonal_clouds_iz
    logical :: do_mgroup_prescribed_iz
    logical :: do_obs_clouds_iz
    logical :: do_no_clouds_iz
    logical :: do_diag_clouds_iz
    logical :: do_specified_clouds_iz
    logical :: do_donner_deep_clouds_iz
    logical :: do_uw_clouds_iz
    logical :: do_zetac_clouds_iz
    logical :: do_random_overlap_iz
    logical :: do_max_random_overlap_iz
    logical :: do_stochastic_clouds_iz
    logical :: use_temp_for_seed_iz
    logical :: do_specified_strat_clouds_iz
    logical :: do_ica_calcs_iz
    logical :: do_liq_num_iz
    logical :: do_ice_num_iz
    logical :: using_fu2007_iz
end type cloudrad_control_type

!------------------------------------------------------------------

public fsrad_output_type

!    tdtsw
!    tdtlw
!    tdtsw_clr
!    tdtlw_clr
!    swdns
!    swups
!    lwdns
!    lwups
!    swin
!    swout
!    olr
!    swdns_clr
!    swups_clr
!    lwdns_clr
!    lwups_clr`
!    swin_clr
!    swout_clr
!    olr_clr
!    npass

type fsrad_output_type
     real, dimension(:,:,:), pointer :: tdtsw=>NULL(), &
                                        tdtlw=>NULL(),  &
                                        tdtsw_clr=>NULL(),  &
                                        tdtlw_clr=>NULL()
     real, dimension(:,:),   pointer :: swdns=>NULL(),   &
                                        swups=>NULL(),  &
                                        lwups=>NULL(), &
                                        lwdns=>NULL(), &
                                        swin=>NULL(), &
                                        swout=>NULL(), &
                                        olr=>NULL(), &
                                        swdns_clr=>NULL(),  &
                                        swups_clr=>NULL(),  &
                                        lwups_clr=>NULL(),&
                                        lwdns_clr=>NULL(),   &
                                        swin_clr=>NULL(),  &
                                        swout_clr=>NULL(), &
                                        olr_clr=>NULL()
     integer      :: npass
end type fsrad_output_type

!-------------------------------------------------------------------

public gas_tf_type

!    tdav
!    tlsqu
!    tmpdiff
!    tstdav
!    co2nbl
!    n2o9c
!    tn2o17
!    co2spnb
!    a1
!    a2

type gas_tf_type
     real, dimension(:,:,:),   pointer :: tdav=>NULL(),   &
                                          tlsqu=>NULL(),   &
                                          tmpdiff=>NULL(),   &
                                          tstdav=>NULL(),  &
                                          co2nbl=>NULL(),   &
                                          n2o9c=>NULL(),   &
                                          tn2o17=>NULL()
     real, dimension(:,:,:,:), pointer :: co2spnb=>NULL()
     real, dimension(:,:),     pointer :: a1=>NULL(),    &
                                          a2=>NULL()
end type gas_tf_type

!------------------------------------------------------------------

public longwave_control_type 

type longwave_control_type
    character(len=16) :: lw_form
    character(len=16) :: continuum_form
    character(len=16) :: linecatalog_form
    logical           :: do_cfc
    logical           :: do_lwaerosol
    logical           :: do_ch4
    logical           :: do_n2o
    logical           :: do_ch4lbltmpint
    logical           :: do_n2olbltmpint
    logical           :: do_co2
    logical           :: do_lwcldemiss
    logical           :: do_h2o
    logical           :: do_o3 
    logical           :: do_cfc_iz
    logical           :: do_lwaerosol_iz
    logical           :: do_ch4_iz
    logical           :: do_n2o_iz
    logical           :: do_ch4lbltmpint_iz
    logical           :: do_n2olbltmpint_iz
    logical           :: do_co2_iz
    logical           :: do_lwcldemiss_iz
    logical           :: do_h2o_iz
    logical           :: do_o3_iz 
end type longwave_control_type

!---------------------------------------------------------------------

public longwave_parameter_type

type longwave_parameter_type
     integer   :: offset
     integer   :: NBTRG
     integer   :: NBTRGE
     integer   :: NBLY
     integer   :: n_lwaerosol_bands
     real      :: lw_band_resolution
     logical   :: offset_iz
     logical   :: NBTRG_iz
     logical   :: NBTRGE_iz
     logical   :: NBLY_iz
     logical   :: n_lwaerosol_bands_iz
     logical   :: lw_band_resolution_iz
end type longwave_parameter_type

!--------------------------------------------------------------------

public longwave_tables1_type

!    vae
!    td
!    md
!    cd

type longwave_tables1_type
    real, dimension(:,:), pointer  ::  vae=>NULL(),   &
                                       td=>NULL(), &
                                       md=>NULL(), &
                                       cd=>NULL()
end type longwave_tables1_type

!--------------------------------------------------------------------

public longwave_tables2_type

!    vae
!    td
!    md
!    cd

type longwave_tables2_type
    real, dimension(:,:,:), pointer  ::  vae=>NULL(),  &
                                         td=>NULL(),  &
                                         md=>NULL(),   &
                                         cd=>NULL()
end type longwave_tables2_type

!---------------------------------------------------------------------

public longwave_tables3_type

!    vae
!    td

type longwave_tables3_type
     real,  dimension(:,:), pointer    ::  vae=>NULL(),   &
                                           td=>NULL()          
end type longwave_tables3_type

!---------------------------------------------------------------------

public lw_clouds_type

!    taucld_rndlw
!    taucld_mxolw
!    taunbl_mxolw

type lw_clouds_type
     real, dimension(:,:,:,:),   pointer :: taucld_rndlw=>NULL(), &
                                            taucld_mxolw=>NULL(), &
                                            taunbl_mxolw=>NULL()
end type lw_clouds_type

!------------------------------------------------------------------

public lw_diagnostics_type

!    flx1e1
!    gxcts
!    flx1e1f
!    excts
!    fctsg
!    fluxn
!    fluxncf
!    exctsn
!    cts_out
!    cts_outcf

type lw_diagnostics_type
     real, dimension(:,:),   pointer   :: flx1e1=>NULL(),  &
                                          gxcts=>NULL()
     real, dimension(:,:,:), pointer   :: flx1e1f=>NULL(),  &
                                          excts=>NULL(),&
                                          fctsg=>NULL()
     real, dimension(:,:,:,:), pointer :: fluxn=>NULL(),   &
                                          fluxncf=>NULL(),   &
                                          exctsn=>NULL(),  &
                                          cts_out=>NULL(), &
                                          cts_outcf=>NULL()
end type lw_diagnostics_type

!-------------------------------------------------------------------

public lw_output_type

!    heatra
!    flxnet
!    heatracf
!    flxnetcf
!    netlw_special
!    netlw_special_clr

type lw_output_type
     real, dimension(:,:,:), pointer :: heatra=>NULL(), &
                                        flxnet=>NULL(),  &
                                        heatracf=>NULL(), &
                                        flxnetcf=>NULL()
     real, dimension(:,:,:), pointer   :: netlw_special=>NULL(), &
                                          netlw_special_clr=>NULL(), &
                                          bdy_flx=>NULL(), &
                                          bdy_flx_clr=>NULL()
end type lw_output_type

!------------------------------------------------------------------

public lw_table_type

!    bdlocm
!    bdhicm
!    bandlo
!    bandhi
!    iband

type lw_table_type
     real, dimension(:),    pointer :: bdlocm=>NULL(),   &
                                       bdhicm=>NULL(),  &
                                       bandlo=>NULL(),  &
                                       bandhi=>NULL()
     integer, dimension(:), pointer :: iband=>NULL()
end type lw_table_type

!------------------------------------------------------------------

public microphysics_type
 
!    conc_ice
!    conc_drop
!    size_ice
!    size_drop
!    size_snow
!    conc_snow
!    size_rain
!    conc_rain
!    cldamt
!    stoch_conc_ice
!    stoch_conc_drop
!    stoch_size_ice
!    stoch_size_drop
!    stoch_cldamt
!    stoch_cloud_type
!    lw_stoch_conc_ice
!    lw_stoch_conc_drop
!    lw_stoch_size_ice
!    lw_stoch_size_drop
!    lw_stoch_cldamt
!    sw_stoch_conc_ice
!    sw_stoch_conc_drop
!    sw_stoch_size_ice
!    sw_stoch_size_drop
!    sw_stoch_cldamt

type microphysics_type
   real, dimension(:,:,:), pointer :: conc_ice=>NULL(),   &
                                      conc_drop=>NULL(),      &
                                      size_ice=>NULL(),   &
                                      size_drop=>NULL(),     &
                                      size_snow=>NULL(),   &
                                      conc_snow=>NULL(),     &
                                      size_rain=>NULL(),     &
                                      conc_rain=>NULL(),   &
                                      cldamt=>NULL(),      &
                                      droplet_number=>NULL(), &
                                      ice_number=>NULL()
real, dimension(:,:,:,:), pointer :: stoch_conc_ice=>NULL(),   &
                                     stoch_conc_drop=>NULL(),  &
                                     stoch_size_ice=>NULL(),   &
                                     stoch_size_drop=>NULL(),  &
                                     stoch_cldamt=>NULL(),     &
                                     stoch_droplet_number=>NULL(), &
                                     stoch_ice_number=>NULL()
integer, dimension(:,:,:,:), pointer ::  stoch_cloud_type=>NULL()
!
! In practice, we allocate a single set of columns for the stochastic
!   clouds, then point to sections of the larger array with the 
!   lw_ and sw_type. 
!   I.e. lw_stoch_conc_ice => stoch_conc_ice(:, :, :, 1:numLwBands)
!
real, dimension(:,:,:,:), pointer :: lw_stoch_conc_ice=>NULL(),   &
                                     lw_stoch_conc_drop=>NULL(),  &
                                     lw_stoch_size_ice=>NULL(),   &
                                     lw_stoch_size_drop=>NULL(),  &
                                     lw_stoch_cldamt=>NULL(),     &
                                     lw_stoch_droplet_number=>NULL(), &
                                     lw_stoch_ice_number=>NULL(), &
                                     sw_stoch_conc_ice=>NULL(),   &
                                     sw_stoch_conc_drop=>NULL(),  &
                                     sw_stoch_size_ice=>NULL(),   &
                                     sw_stoch_size_drop=>NULL(),  &
                                     sw_stoch_cldamt=>NULL(),     &
                                     sw_stoch_droplet_number=>NULL(), &
                                     sw_stoch_ice_number=>NULL()
end type microphysics_type

!-------------------------------------------------------------------

public microrad_properties_type
 
!    cldext
!    cldsct
!    cldasymm
!    abscoeff

type microrad_properties_type
   real, dimension(:,:,:,:), pointer :: cldext=>NULL(),  &
                                        cldsct=>NULL(), &
                                        cldasymm=>NULL(),    &
                                        abscoeff=>NULL()
end type microrad_properties_type

!--------------------------------------------------------------------

public optical_path_type

!    empl1f
!    vrpfh2o
!    xch2obd
!    tphfh2o
!    avephif
!    totaerooptdep
!    empl1
!    empl2
!    var1
!    var2
!    emx1f
!    emx2f
!    totvo2
!    avephi
!    totch2obdwd
!    xch2obdwd
!    totphi
!    cntval
!    toto3
!    tphio3
!    var3
!    var4
!    wk
!    rh2os
!    rfrgn
!    tfac
!    totaerooptdep15
!    totf11
!    totf12
!    totf113
!    totf22
!    emx1
!    emx2
!    csfah2o
!    aerooptdep_KE_15

type optical_path_type
     real, dimension (:,:,:,:), pointer :: empl1f=>NULL(),  &
                                           empl2f=>NULL(),  &
                                           vrpfh2o=>NULL(), &
                                           xch2obd=>NULL(),  &
                                           tphfh2o=>NULL(), &
                                           avephif=>NULL(), &
                                           totaerooptdep=>NULL()
     real, dimension (:,:,:),   pointer :: empl1=>NULL(), &
                                           empl2=>NULL(),  &
                                           var1=>NULL(), &
                                           var2=>NULL(), &
                                           emx1f=>NULL(),   &
                                           emx2f=>NULL(),   &
                                           totvo2=>NULL(),  &
                                           avephi=>NULL(),&
                                           totch2obdwd=>NULL(), &
                                           xch2obdwd=>NULL(), &
                                           totphi=>NULL(),   &
                                           cntval=>NULL(), &
                                           toto3=>NULL(),   &
                                           tphio3=>NULL(),  &
                                           var3=>NULL(),  &
                                           var4=>NULL(),        &
                                           wk=>NULL(),         &
                                           rh2os=>NULL(),  &
                                           rfrgn=>NULL(),  &
                                           tfac=>NULL(), &
                                           totaerooptdep_15=>NULL(), &
                                           totf11=>NULL(),   &
                                           totf12=>NULL(),  &
                                           totf113=>NULL(),   &
                                           totf22=>NULL()
      real, dimension (:,:), pointer    :: emx1=>NULL(),  &
                                           emx2=>NULL(),  &
                                           csfah2o=>NULL(), &
                                           aerooptdep_KE_15=>NULL()
end type optical_path_type

!------------------------------------------------------------------

public radiation_control_type

type radiation_control_type
    logical  :: do_totcld_forcing
    logical  :: do_aerosol
    integer  :: rad_time_step
    integer  :: lw_rad_time_step
    integer  :: sw_rad_time_step
    logical  :: do_sw_rad
    logical  :: do_lw_rad
    logical  :: hires_coszen
    integer  :: nzens
    real     :: co2_tf_calc_intrvl
    logical  :: use_current_co2_for_tf
    logical  :: calc_co2_tfs_on_first_step
    logical  :: calc_co2_tfs_monthly
    real     :: co2_tf_time_displacement
    real     :: ch4_tf_calc_intrvl
    logical  :: use_current_ch4_for_tf
    logical  :: calc_ch4_tfs_on_first_step
    logical  :: calc_ch4_tfs_monthly
    real     :: ch4_tf_time_displacement
    real     :: n2o_tf_calc_intrvl
    logical  :: use_current_n2o_for_tf
    logical  :: calc_n2o_tfs_on_first_step
    logical  :: calc_n2o_tfs_monthly
    real     :: n2o_tf_time_displacement
    integer  :: mx_spec_levs
    logical  :: time_varying_solar_constant
    logical  :: volcanic_sw_aerosols
    logical  :: volcanic_lw_aerosols
    logical  :: using_solar_timeseries_data
    logical  :: do_lwaerosol_forcing
    logical  :: do_swaerosol_forcing
    integer  :: indx_swaf
    integer  :: indx_lwaf
    logical  :: using_im_bcsul
    logical  :: do_totcld_forcing_iz
    logical  :: do_aerosol_iz
    logical  :: rad_time_step_iz
    logical  :: lw_rad_time_step_iz
    logical  :: sw_rad_time_step_iz
    logical  :: do_sw_rad_iz
    logical  :: do_lw_rad_iz
    logical  :: hires_coszen_iz
    logical  :: nzens_iz  
    logical  :: co2_tf_calc_intrvl_iz
    logical  :: use_current_co2_for_tf_iz
    logical  :: calc_co2_tfs_on_first_step_iz
    logical  :: calc_co2_tfs_monthly_iz
    logical  :: ch4_tf_calc_intrvl_iz
    logical  :: use_current_ch4_for_tf_iz
    logical  :: calc_ch4_tfs_on_first_step_iz
    logical  :: calc_ch4_tfs_monthly_iz
    logical  :: n2o_tf_calc_intrvl_iz
    logical  :: use_current_n2o_for_tf_iz
    logical  :: calc_n2o_tfs_on_first_step_iz
    logical  :: calc_n2o_tfs_monthly_iz
    logical  :: co2_tf_time_displacement_iz
    logical  :: ch4_tf_time_displacement_iz
    logical  :: n2o_tf_time_displacement_iz
    logical  :: mx_spec_levs_iz
    logical  :: time_varying_solar_constant_iz
    logical  :: volcanic_sw_aerosols_iz
    logical  :: volcanic_lw_aerosols_iz
    logical  :: using_solar_timeseries_data_iz
    logical  :: do_lwaerosol_forcing_iz
    logical  :: do_swaerosol_forcing_iz
    logical  :: indx_swaf_iz
    logical  :: indx_lwaf_iz
    logical  :: using_im_bcsul_iz
end type radiation_control_type

!------------------------------------------------------------------

public   radiative_gases_type
 
!    qo3
!    rrvch4
!    rrvn2o
!    rrvco2
!    rrvf11
!    rrvf12
!    rrvf113
!    rrvf22
!    rf11air
!    rf12air
!    rf113air
!    rf22air
!    time_varying_co2
!    time_varying_f11
!    time_varying_f12
!    time_varying_f113
!    time_varying_f22
!    time_varying_ch4
!    time_varying_n2o

type radiative_gases_type
     real, dimension(:,:,:), pointer :: qo3=>NULL()
     real                            :: rrvch4, rrvn2o, rrvco2,    &
                                        rrvf11, rrvf12, rrvf113,  &
                                        rrvf22, rf11air, rf12air,  &
                                        rf113air, rf22air, &
                                        co2_for_last_tf_calc,  &
                                        co2_tf_offset, &
                                        co2_for_next_tf_calc, &
                                        ch4_for_last_tf_calc,  &
                                        ch4_tf_offset, &
                                        ch4_for_next_tf_calc, &
                                        n2o_for_last_tf_calc,  &
                                        n2o_tf_offset, &
                                        n2o_for_next_tf_calc
     logical                         :: time_varying_co2,  &
                                        time_varying_f11, &
                                        time_varying_f12,  &
                                        time_varying_f113, &
                                        time_varying_f22,  &
                                        time_varying_ch4, &
                                        time_varying_n2o, &
                                        use_model_supplied_co2
     type(time_type)                 :: Co2_time, Ch4_time, N2o_time
end type radiative_gases_type

!------------------------------------------------------------------

public rad_output_type

!    tdt_rad
!    tdt_rad_clr
!    tdtsw
!    tdtsw_clr
!    tdtlw
!    flux_sw_surf
!    flux_lw_surf
!    coszen_angle

type rad_output_type
     real, dimension(:,:,:,:), pointer :: tdt_rad=>NULL(),  &
                                        ufsw=>NULL(),  &
                                        dfsw=>NULL(),  &
                                        tdtsw=>NULL()  
     real, dimension(:,:,:,:), pointer :: tdt_rad_clr=>NULL(), &
                                        ufsw_clr=>NULL(),  &
                                        dfsw_clr=>NULL(),  &
                                        tdtsw_clr=>NULL()
                                        
     real, dimension(:,:,:), pointer :: tdtlw=>NULL()
     real, dimension(:,:,:), pointer :: flxnet=>NULL()
     real, dimension(:,:,:), pointer :: flxnetcf=>NULL()
     real, dimension(:,:,:), pointer :: tdtlw_clr=>NULL()
     real, dimension(:,:,:),   pointer :: flux_sw_surf=>NULL(), &
                                        flux_sw_surf_refl_dir=>NULL(), &
                                        flux_sw_surf_dir=>NULL(), &
                                        flux_sw_surf_dif=>NULL(), &
                                        flux_sw_down_vis_dir=>NULL(), &
                                        flux_sw_down_vis_dif=>NULL(), &
                                       flux_sw_down_total_dir=>NULL(), &
                                       flux_sw_down_total_dif=>NULL(), &
                                        flux_sw_vis=>NULL(), &
                                        flux_sw_vis_dir=>NULL(), &
                                        flux_sw_refl_vis_dir=>NULL(), &
                                        flux_sw_vis_dif=>NULL()
     real, dimension(:,:,:),   pointer :: flux_sw_down_vis_clr=>NULL(), &
                                  flux_sw_down_total_dir_clr=>NULL(), &
                                  flux_sw_down_total_dif_clr=>NULL()
     real, dimension(:,:),   pointer :: flux_lw_surf=>NULL(), &
                                        coszen_angle=>NULL()
end type rad_output_type

!-------------------------------------------------------------------

public shortwave_control_type

type shortwave_control_type
    logical  :: do_lhsw
    logical  :: do_esfsw
    logical  :: do_swaerosol
    logical  :: do_diurnal
    logical  :: do_annual
    logical  :: do_daily_mean
    logical  :: do_cmip_diagnostics
    real     :: solar_constant
    logical  :: do_lhsw_iz
    logical  :: do_esfsw_iz
    logical  :: do_swaerosol_iz
    logical  :: do_diurnal_iz
    logical  :: do_annual_iz
    logical  :: do_daily_mean_iz
    logical  :: do_cmip_diagnostics_iz
end type shortwave_control_type

!---------------------------------------------------------------------

public solar_spectrum_type

!    solarfluxtoa      highly-resolved solar flux at toa in 
!                      Sw_control%tot_wvnums 
!                      bands [         ]
!    solflxband_lean   a time series of toa solar flux in each
!                      parameterization band
!    solflxband_lean_ann_1882   1882 average toa solar flux in each
!                               parameterization band
!    solflxband_lean_ann_2000   2000 average toa solar flux in each
!                               parameterization band
!    solflxband        toa solar flux in each parameterization band
!                      [         ]
!    endwvnbands       highest wave number in each of the solar
!                      spectral parameterization bands [ cm (-1) ]
!    tot_wvnums
!    nbands
!    nfrqpts
!    nstreams
!    nh2obands
!    visible_band_indx
!    one_micron_indx
!    eight70_band_indx
!    visible_band_indx_iz
!    one_micron_indx_iz
!    eight70_band_indx_iz

type solar_spectrum_type
    real, dimension(:),    pointer   :: solarfluxtoa=>null()
    real, dimension(:),    pointer   :: solflxband=>NULL()
    real, dimension(:),    pointer   :: solflxbandref=>NULL()
    real, dimension(:),    pointer   :: solflxband_lean_ann_1882=>NULL()
    real, dimension(:),    pointer   :: solflxband_lean_ann_2000=>NULL()
    real, dimension(:,:,:),pointer   :: solflxband_lean=>NULL()
    integer, dimension(:), pointer   :: endwvnbands=>NULL()
    integer         :: tot_wvnums
    integer         :: nbands
    integer         :: nfrqpts
    integer         :: nstreams
    integer         :: nh2obands
    integer         :: visible_band_indx, one_micron_indx
    integer         :: eight70_band_indx
    logical         :: visible_band_indx_iz, one_micron_indx_iz
    logical         :: eight70_band_indx_iz
    integer         :: w340_band_indx, w380_band_indx,  &
                       w440_band_indx, w670_band_indx
    logical         :: w340_band_iz, w380_band_iz, &
                       w440_band_iz, w670_band_iz
end type solar_spectrum_type

!---------------------------------------------------------------------

public surface_type

!    asfc
!    land

type surface_type
    real, dimension(:,:),   pointer ::  asfc=>NULL(),   &
                                        land=>NULL(),  &
                                        asfc_vis_dir=>NULL(), &
                                        asfc_nir_dir=>NULL(), &
                                        asfc_vis_dif=>NULL(), &
                                        asfc_nir_dif=>NULL()
end type surface_type
 
!-------------------------------------------------------------------

public sw_output_type

!    dfsw
!    ufsw
!    fsw
!    hsw
!    dfswcf
!    ufswcf
!    fswcf
!    hswcf
!    swdn_special
!    swup_special
!    swdn_special_clr
!    swup_special_clr

type sw_output_type
     real, dimension(:,:,:,:), pointer :: dfsw=>NULL(),   &
                                        ufsw=>NULL(),  &
                                        fsw=>NULL(),   &
                                        hsw=>NULL()   
     real, dimension(:,:,:,:), pointer :: dfswcf=>NULL(),   &
                                        ufswcf=>NULL(),&
                                        fswcf=>NULL(),  &
                                        hswcf=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_vis_sfc=>NULL(),   &
                                       ufsw_vis_sfc=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_dir_sfc=>NULL()
      real, dimension(:,:,:), pointer :: ufsw_dir_sfc=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_dir_sfc_clr=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_dif_sfc=>NULL(),   &
                                       ufsw_dif_sfc=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_dif_sfc_clr=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_vis_sfc_dir=>NULL()
      real, dimension(:,:,:), pointer :: ufsw_vis_sfc_dir=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_vis_sfc_clr=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_vis_sfc_dif=>NULL(),   &
                                       ufsw_vis_sfc_dif=>NULL()
      real, dimension(:,:,:,:), pointer   ::  bdy_flx=>NULL()
      real, dimension(:,:,:,:), pointer   ::  bdy_flx_clr=>NULL()
      real, dimension(:,:,:,:), pointer   ::                       &
                                        swup_special=>NULL(), &
                                        swup_special_clr=>NULL()
     real, dimension(:,:,:,:), pointer   :: swdn_special=>NULL(), &
                                          swdn_special_clr=>NULL()
end type sw_output_type

!-------------------------------------------------------------------

public table_axis_type

type table_axis_type
  integer :: first_col
  real    :: min_val
  real    :: max_val
  real    :: tab_inc
end type table_axis_type

!---------------------------------------------------------------------

!private      &


!--------------------------------------------------------------------
!-------- namelist  ---------

integer            ::  dummy = 0


namelist / rad_utilities_nml /   &
                                dummy


!---------------------------------------------------------------------
!------- public data ------


type (longwave_control_type),  public   ::    &
     Lw_control = longwave_control_type( '    ', '    ', '    ', &
                                         .false., .false., .false.,  &
                                         .false., .false., .false.,  &
                                         .false., .false.,  &
                                         .false., .false.,  &
                                         .false., .false., .false.,  &
                                         .false., .false.,  &
                                         .false., .false.,  &
                                         .false., .false., .false.   )

type (shortwave_control_type), public   ::  &
    Sw_control = shortwave_control_type( .false., .false., .false. , &
                                         .false., .false., .false., &
                                         .false., &
                                         0.0, &
                                         .false., .false., .false. , &
                                         .false., &
                                         .false., .false., .false.)

type (radiation_control_type), public   ::  &
   Rad_control = radiation_control_type( .false., .false., 0, 0, 0, &
                                         .false., .false., &
                                         .false., 1, &
                                         0.0,  .true.,  .false.,&
                                         .false.,  0.0, &
                                         0.0, .true., .false.,  &
                                         .false.,  0.0,  &
                                         0.0, .true., .false.,   &
                                         .false., 0.0, &
                                         0, .false., .false.,   &
                                         .false., .false.,&
                                         .false., .false.,      &
                                         0, 0, .false., &
! _iz variables:
                                         .false., .false., .false., &
                                         .false., .false., .true., &
                                         .true.,  &
                                         .false., .false., &
                                         .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., &
                                         .false., .false., &
                                         .false., .false., &
                                         .false.,          &
                                         .false., .false., .false.,  &
                                         .false., .false.,   &
                                         .false., .false., &
                                         .false., .false., .false.)

type (cloudrad_control_type), public    ::   &
 Cldrad_control = cloudrad_control_type( .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false.,                   &
                                         0,0,0,0,0,0 , &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false.   )


type (longwave_parameter_type), public  ::   &
Lw_parameters = longwave_parameter_type( 0, 0, 0, 0, 0, 10.0, &
                                         .false., .false., .false.,  &
                                         .false., .false., .true.)

type (table_axis_type),        public   ::    &
               temp_1 = table_axis_type( 1, 100.0, 370.0, 10.0), &
               mass_1 = table_axis_type( 1, -16.0,   1.9,  0.1)


!---------------------------------------------------------------------
!------- private data ------


logical :: module_is_initialized=.false.   ! module is initialized ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------



                           contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! <SUBROUTINE NAME="rad_utilities_init">
!  <OVERVIEW>
!   Subroutine to initialize radiation utility package.
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine reads the input namelist file and initializes 
!   rad_utilities_nml. It then writes out this namelist to the output
!   logfile. It also sets up the radiation calculation environment
!   and the initialization flag variable.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call rad_utilities_init
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine rad_utilities_init

!---------------------------------------------------------------------
!    rad_utilities_init is the constructor for rad_uti;lities_mod.
!---------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables:

      integer    ::  unit, ierr, io, logunit

!---------------------------------------------------------------------
!  local variables:
!
!        unit            io unit number used for namelist file
!        ierr            error code
!        io              error status returned from io operation
!                                
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return
 
!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=rad_utilities_nml, iostat=io)
      ierr = check_nml_error(io,'rad_utilities_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=rad_utilities_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'rad_utilities_nml')
        end do
10      call close_file (unit)
      endif
#endif
 
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                        write (logunit, nml=rad_utilities_nml)

!-------------------------------------------------------------------
!    mark the module as initialized.
!-------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------

end subroutine rad_utilities_init


!#####################################################################
!
! <SUBROUTINE NAME="check_derived_types">
!  <OVERVIEW>
!   Subroutine to verify that all logical elements of derived-type
!   variables were initialized during the initialization phase.
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine checks that all of the variable names ending in 
!   "_iz" in its public derived-type module variables are set to .true..
!   If additional types or variables within current public types are
!   added, it is necessary to also add the corresponding "_iz" variable,
!   initialized in this module to .false., and then set to .true. when 
!   the variable is initialized, and add a check in this routine to 
!   verify that initialization.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call check_derived_types
!  </TEMPLATE>
! </SUBROUTINE>
!

subroutine check_derived_types 

!--------------------------------------------------------------------
!    check_derived_types is called at the end of radiation package
!    initialization to verify that all logical components of public
!    derived-type module variables have been initialized.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    check the components of Lw_control.
!--------------------------------------------------------------------
      if (Lw_control%do_cfc_iz .and. &
          Lw_control%do_lwaerosol_iz .and. &
          Lw_control%do_ch4_iz .and. &
          Lw_control%do_n2o_iz .and. &
          Lw_control%do_ch4lbltmpint_iz .and. &
          Lw_control%do_n2olbltmpint_iz .and. &
          Lw_control%do_co2_iz .and. &
          Lw_control%do_h2o_iz .and. &
          Lw_control%do_o3_iz .and. &
          Lw_control%do_lwcldemiss_iz ) then  
      else
        call error_mesg ('rad_utilities_mod', &
             ' at least one component of Lw_control has not been '//&
                                     'initialized', FATAL)
      endif

!---------------------------------------------------------------------
!    check the components of Sw_control.
!--------------------------------------------------------------------
      if (Sw_control%do_lhsw_iz .and. &
          Sw_control%do_esfsw_iz .and. &
          Sw_control%do_swaerosol_iz .and. &
          Sw_control%do_diurnal_iz .and. &
          Sw_control%do_annual_iz .and. &
          Sw_control%do_cmip_diagnostics_iz .and. &
          Sw_control%do_daily_mean_iz ) then  
      else
        call error_mesg ('rad_utilities_mod', &
             ' at least one component of Sw_control has not been '//&
                                     'initialized', FATAL)
      endif

!---------------------------------------------------------------------
!    check the components of Rad_control.
!--------------------------------------------------------------------
      if (Rad_control%do_totcld_forcing_iz .and. &
          Rad_control%do_lwaerosol_forcing_iz .and.  &
          Rad_control%do_swaerosol_forcing_iz .and.  &
          Rad_control%indx_lwaf_iz .and.   &
          Rad_control%indx_swaf_iz .and.   &
          Rad_control%using_im_bcsul_iz .and. &
          Rad_control%volcanic_sw_aerosols_iz .and. &
          Rad_control%volcanic_lw_aerosols_iz .and. &
          Rad_control%time_varying_solar_constant_iz .and. &
          Rad_control%using_solar_timeseries_data_iz .and. &
          Rad_control%do_aerosol_iz .and.     &
          Rad_control%mx_spec_levs_iz .and.   &
          Rad_control%use_current_co2_for_tf_iz .and. &
          Rad_control%co2_tf_calc_intrvl_iz .and. &
          Rad_control%calc_co2_tfs_on_first_step_iz .and. &
          Rad_control%calc_co2_tfs_monthly_iz .and. &
          Rad_control%co2_tf_time_displacement_iz .and. &
          Rad_control%use_current_ch4_for_tf_iz .and. &
          Rad_control%ch4_tf_calc_intrvl_iz .and. &
          Rad_control%calc_ch4_tfs_on_first_step_iz .and. &
          Rad_control%calc_ch4_tfs_monthly_iz .and. &
          Rad_control%ch4_tf_time_displacement_iz .and. &
          Rad_control%use_current_n2o_for_tf_iz .and. &
          Rad_control%n2o_tf_calc_intrvl_iz .and. &
          Rad_control%calc_n2o_tfs_on_first_step_iz .and. &
          Rad_control%calc_n2o_tfs_monthly_iz .and. &
          Rad_control%n2o_tf_time_displacement_iz .and. &
          Rad_control%do_lw_rad_iz .and. &
          Rad_control%do_sw_rad_iz .and. &
          Rad_control%nzens_iz .and. &
          Rad_control%hires_coszen_iz .and. &
          Rad_control%lw_rad_time_step_iz .and.  &
          Rad_control%sw_rad_time_step_iz .and.  &
          Rad_control%rad_time_step_iz ) then
      else
        call error_mesg ('rad_utilities_mod', &
          ' at least one component of Rad_control has not been '//&
                                     'initialized', FATAL)
      endif

!---------------------------------------------------------------------
!    check the components of Cldrad_control.
!--------------------------------------------------------------------
      if (Cldrad_control%do_pred_cld_microphys_iz .and. &
          Cldrad_control%do_presc_cld_microphys_iz .and. &
          Cldrad_control%do_bulk_microphys_iz .and. &
          Cldrad_control%do_sw_micro_iz .and. &
          Cldrad_control%do_lw_micro_iz .and. &
          Cldrad_control%do_strat_clouds_iz .and. &
          Cldrad_control%do_rh_clouds_iz .and. &
          Cldrad_control%do_zonal_clouds_iz .and. &
          Cldrad_control%do_mgroup_prescribed_iz .and. &
          Cldrad_control%do_obs_clouds_iz .and. &
          Cldrad_control%do_no_clouds_iz .and. &
          Cldrad_control%do_diag_clouds_iz .and. &
          Cldrad_control%do_specified_clouds_iz .and. &
          Cldrad_control%do_specified_strat_clouds_iz .and. &
          Cldrad_control%do_donner_deep_clouds_iz .and. &
          Cldrad_control%do_uw_clouds_iz .and. &
          Cldrad_control%do_zetac_clouds_iz .and. &
          Cldrad_control%do_stochastic_clouds_iz .and. &
          Cldrad_control%use_temp_for_seed_iz .and. &
          Cldrad_control%do_random_overlap_iz .and. &
          Cldrad_control%do_ica_calcs_iz .and. &
          Cldrad_control%using_fu2007_iz .and.  &
          Cldrad_control%do_liq_num_iz .and.  &
          Cldrad_control%do_ice_num_iz .and.  &
          Cldrad_control%do_max_random_overlap_iz ) then     
      else
        call error_mesg ('rad_utilities_mod', &
          ' at least one component of Cldrad_control has not been '//&
                                     'initialized', FATAL)
      endif

!---------------------------------------------------------------------
!    check the components of Lw_parameters.
!--------------------------------------------------------------------
      if (Lw_parameters%offset_iz .and. &
          Lw_parameters%nbtrg_iz .and. &
          Lw_parameters%nbtrge_iz .and. &
          Lw_parameters%nbly_iz .and. &
          Lw_parameters%n_lwaerosol_bands_iz .and. &
          Lw_parameters%lw_band_resolution_iz) then
      else
        call error_mesg ('rad_utilities_mod', &
          ' at least one component of Lw_parameters has not been '//&
                                     'initialized', FATAL)
      endif

!---------------------------------------------------------------------
!    check for consistency between band structure and cfc and lwaerosol
!    effects.
!---------------------------------------------------------------------
      if (Lw_parameters%nbtrge == 0) then
        if (Lw_control%do_cfc .or. Lw_control%do_lwaerosol) then
          call error_mesg ('rad_utilities_mod', &
             'when do_cfc and / or do_lwaerosol is .true., must set &
              &sealw99_nml variable no_h2o_bands_1200_1400 > 0 ', FATAL)
        endif
      endif
!--------------------------------------------------------------------

!--------------------------------------------------------------------


end subroutine check_derived_types 



!####################################################################
! <SUBROUTINE NAME="locate_in_table">
!  <OVERVIEW>
!   Subroutine to locate index and residual value from an array provided 
!   with array and axis information
!  </OVERVIEW>
!  <DESCRIPTION>
!     given array x and an arithmetic sequence of table column headings
!     tabxmin, tabxmin+tabdeltax, ..., corresponding to column ixlow, 
!     ixlow+1, ..., ixupp, Locate returns the array ix is column 
!     indices and the array dx of residuals.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call locate_in_table(table_axis, x, dx, ix, k_min, k_max)
!  </TEMPLATE>
!  <IN NAME="table_axis" TYPE="table_axis_type">
!   table_axis contains the axis information such as, min, increment,
!   and first column values.
!  </IN>
!  <IN NAME="x" TYPE="real">
!   array from which data is to be searched
!  </IN>
!  <OUT NAME="dx" TYPE="real">
!   residual between x and x(ix+first_column)
!  </OUT>
!  <OUT NAME="ix" TYPE="integer">
!   index values of the searched domain in the array
!  </OUT>
!  <IN NAME="k_min" TYPE="integer">
!   minimum k value of the search domain 
!  </IN>
!  <IN NAME="k_max" TYPE="integer">
!   maximum k value of the search domain
!  </IN>
! </SUBROUTINE>
subroutine locate_in_table (table_axis, x, dx, ix, k_min, k_max) 

!---------------------------------------------------------------------
!    given array x and an arithmetic sequence of table column headings
!    tabxmin, tabxmin+tabdeltax, ..., corresponding to column ixlow, 
!    ixlow+1, ..., ixupp, locate_in_table returns the array ix of
!    column indices and the array dx of residuals.
!    author: c. h. goldberg
!    revised: 1/1/93
!    certified:  radiation version 1.0
!----------------------------------------------------------------------

type(table_axis_type),     intent(in)  :: table_axis
real,    dimension(:,:,:), intent(in)  :: x
integer,                   intent(in)  :: k_min, k_max
real,    dimension(:,:,:), intent(out) :: dx
integer, dimension(:,:,:), intent(out) :: ix

!--------------------------------------------------------------------
!  intent(in) variables:
!
!    table_axis
!    x
!    k_min
!    k_max
!
!  intent(out) variables:
!
!    dx
!    ix
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension (size(x,1), size(x,2), size(x,3))  ::  fx
      integer     ::  k

!---------------------------------------------------------------------
!  local variables:
!
!     fx
!     table_min
!     table_inc
!     k
!     table_col
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('rad_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do k=k_min,k_max
        fx (:,:,k) = AINT((x(:,:,k) - table_axis%min_val )/  &
                     table_axis%tab_inc)
        ix (:,:,k) = INT(fx(:,:,k)) + table_axis%first_col
        dx (:,:,k) = x(:,:,k) - fx(:,:,k)*table_axis%tab_inc - &
                     table_axis%min_val
      end do
      
!---------------------------------------------------------------------


end subroutine locate_in_table



!####################################################################
! <SUBROUTINE NAME="looktab_type1">
!  <OVERVIEW>
!   Subroutine to calculate answer from input differentials.
!  </OVERVIEW>
!  <DESCRIPTION>
!   given arrays ix(:,:,:) and iy(:,:,:) of integral subscripts and
!     arrays dx(:,:,:) and dy(:,:,:) of differences from x(:,:,:) and
!     y(:,:,:), calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:))
!     from four tables of values, f, df/dx, df/dy, and d2f/dxdy.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call looktab_type1 (tab, ix, iy, dx, dy, answer, k_min, k_max)
!  </TEMPLATE>
!  <IN NAME="tab" TYPE="longwave_tables1_type">
!   The data array that contains function values and differentials
!  </IN>
!  <IN NAME="ix" TYPE="integer">
!   x subscript of input data array
!  </IN>
!  <IN NAME="iy" TYPE="integer">
!   y subscript of input data array
!  </IN>
!  <IN NAME="dx" TYPE="real">
!   x step in the x subscript space
!  </IN>
!  <IN NAME="dy" TYPE="real">
!   y step in the y subscript space
!  </IN>
!  <OUT NAME="answer" TYPE="real">
!   the answer to be calculated
!  </OUT>
!  <IN NAME="k_min" TYPE="integer">
!   the minimum k value of the domain
!  </IN>
!  <IN NAME="k_max" TYPE="integer">
!   the maximum k value of the domain
!  </IN>
! </SUBROUTINE>
!
subroutine looktab_type1 (tab, ix, iy, dx, dy, answer, k_min, k_max)   

!----------------------------------------------------------------------
!    given arrays ix(:,:,:) and iy(:,:,:) of integral subscripts and
!    arrays dx(:,:,:) and dy(:,:,:) of differences from x(:,:,:) and
!    y(:,:,:), calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:))
!    from four tables of values, f, df/dx, df/dy, and d2f/dxdy.
!    author: c. h. goldberg
!    revised: 1/1/93
!    certified:  radiation version 1.0
!--------------------------------------------------------------------

type(longwave_tables1_type), intent(in)  :: tab
integer,dimension(:,:,:),    intent(in)  :: ix, iy
real,   dimension(:,:,:),    intent(in)  :: dx, dy   
real,   dimension(:,:,:),    intent(out) :: answer
integer,                     intent(in)  :: k_min, k_max

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    tab
!    ix
!    iy
!    dx
!    dy
!    k_min
!    k_max
!
!  intent(out) variables:
!
!    answer
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer    ::  i_min, i_max, j_min, j_max, i, j, k
  
!--------------------------------------------------------------------
!  local variables:
!
!    i_min
!    i_max
!    j_min
!    j_max
!    i,j,k
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('rad_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

      i_min = lbound(ix,1)
      i_max = ubound(ix,1)
      j_min = lbound(ix,2)
      j_max = ubound(ix,2)

      do k=k_min, k_max
        do j=j_min, j_max
          do i=i_min, i_max
            answer(i,j,k) =                                         &
                                      tab%vae (ix(i,j,k), iy(i,j,k)) + &
                            dx(i,j,k)*tab%td  (ix(i,j,k), iy(i,j,k)) + &
                            dy(i,j,k)*tab%md  (ix(i,j,k), iy(i,j,k)) + &
                  dx(i,j,k)*dy(i,j,k)*tab%cd(ix(i,j,k), iy(i,j,k))
          end do
        end do
      end do

!---------------------------------------------------------------------


end subroutine looktab_type1



!#####################################################################
! <SUBROUTINE NAME="looktab_type2">
!  <OVERVIEW>
!   Subroutine to calculate answer from input differentials.
!  </OVERVIEW>
!  <DESCRIPTION>
!   given arrays ix(:,:,:) and iy(:,:,:) of integral subscripts and
!     arrays dx(:,:,:) and dy(:,:,:) of differences from x(:,:,:) and
!     y(:,:,:), calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:))
!     from four tables of values, f, df/dx, df/dy, and d2f/dxdy.
!     The difference between this version about the version above is
!     that the differential arrays are 3 dimensional.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call looktab_type2 (tab, ix, iy, dx, dy, answer, k_min, k_max, m)
!  </TEMPLATE>
!  <IN NAME="tab" TYPE="longwave_tables2_type">
!   The data array that contains function values and differentials
!  </IN>
!  <IN NAME="ix" TYPE="integer">
!   x subscript of input data array
!  </IN>
!  <IN NAME="iy" TYPE="integer">
!   y subscript of input data array
!  </IN>
!  <IN NAME="dx" TYPE="real">
!   x step in the x subscript space
!  </IN>
!  <IN NAME="dy" TYPE="real">
!   y step in the y subscript space
!  </IN>
!  <OUT NAME="answer" TYPE="real">
!   the answer to be calculated
!  </OUT>
!  <IN NAME="k_min" TYPE="integer">
!   the minimum k value of the domain
!  </IN>
!  <IN NAME="k_max" TYPE="integer">
!   the maximum k value of the domain
!  </IN>
!  <IN NAME="m" TYPE="integer">
!   the z indice of the differential arrays
!  </IN>
! </SUBROUTINE>
!
subroutine looktab_type2 (tab, ix, iy, dx, dy, answer, k_min, k_max, m)

!-------------------------------------------------------------------
!    given arrays ix(:,:,:) and iy(:,:,:) of integral subscripts and
!    arrays dx(:,:,:) and dy(:,:,:) of differences from x(:,:,:) and
!    y(:,:,:), calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:))
!    from four tables of values, f, df/dx, df/dy, and d2f/dxdy.
!    author: c. h. goldberg
!    revised: 1/1/93
!    certified:  radiation version 1.0
!--------------------------------------------------------------------

type(longwave_tables2_type), intent(in)   :: tab
integer, dimension (:,:,:),  intent(in)   :: ix, iy
integer,                     intent(in)   :: m
real, dimension (:,:,:),     intent(in)   :: dx, dy
real, dimension (:,:,:),     intent(out)  :: answer
integer,                     intent(in)   :: k_min, k_max

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    tab
!    ix
!    iy
!    m
!    dx
!    dy
!    k_min
!    k_max
!
!  intent(out) variables:
!
!    answer
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

       integer    ::    i_min, i_max, j_min, j_max
       integer    ::    i, j, k
  
!--------------------------------------------------------------------
!  local variables:
!
!    i_min
!    i_max
!    j_min
!    j_max
!    i,j,k
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('rad_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

      i_min = lbound(ix,1)
      i_max = ubound(ix,1)
      j_min = lbound(ix,2)
      j_max = ubound(ix,2)

      do k=k_min, k_max
        do j=j_min, j_max
          do i=i_min, i_max
            answer(i,j,k) =                                           &
                                   tab%vae (ix(i,j,k), iy(i,j,k),m) + &
                         dx(i,j,k)*tab%td (ix(i,j,k), iy(i,j,k),m) + &
                         dy(i,j,k)*tab%md (ix(i,j,k), iy(i,j,k),m) + &
               dx(i,j,k)*dy(i,j,k)*tab%cd   (ix(i,j,k), iy(i,j,k),m)
           end do
        end do
      end do

!--------------------------------------------------------------------

end subroutine looktab_type2



!###################################################################
! <SUBROUTINE NAME="looktab_type3">
!  <OVERVIEW>
!   Subroutine to calculate answer from input differentials.
!  </OVERVIEW>
!  <DESCRIPTION>
!   given arrays ix(:,:,:) and iy(:,:,:) of integral subscripts and
!     arrays dx(:,:,:) and dy(:,:,:) of differences from x(:,:,:) and
!     y(:,:,:), calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:))
!     from four tables of values, f, df/dx, df/dy, and d2f/dxdy.
!   In this version, only f(x,y) and f(x,y)+dx*df/dx is used. Probably
!   the f(x,y) is homogeneous in y space.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call looktab_type3 (tab, ix, dx,  answer, k_min, k_max, n)
!  </TEMPLATE>
!  <IN NAME="tab" TYPE="longwave_tables3_type">
!   The data array that contains function values and differentials
!  </IN>
!  <IN NAME="ix" TYPE="integer">
!   x subscript of input data array
!  </IN>
!  <IN NAME="dx" TYPE="real">
!   x step in the x subscript space
!  </IN>
!  <OUT NAME="answer" TYPE="real">
!   the answer to be calculated
!  </OUT>
!  <IN NAME="k_min" TYPE="integer">
!   the minimum k value of the domain
!  </IN>
!  <IN NAME="k_max" TYPE="integer">
!   the maximum k value of the domain
!  </IN>
!  <IN NAME="n" TYPE="integer">
!   the z indice of the differential arrays
!  </IN>
! </SUBROUTINE>
!
subroutine looktab_type3 (tab, ix, dx,  answer, k_min, k_max, n)

!----------------------------------------------------------------------
!
!    given arrays ix(:,:,:) and dx(:,:,:) of integer subscripts and!
!    differences from x(:,:,:) and constant column subscript iyconst, 
!    calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:)) from four tables
!    of values f, df/dx, df/dy, and d2f/dxdy.
!    author: c. h. goldberg
!    revised: 1/1/93
!    certified:  radiation version 1.0
!-----------------------------------------------------------------------

type(longwave_tables3_type), intent(in)  :: tab
integer, dimension (:,:,:),  intent(in)  :: ix
integer,                     intent(in)  :: n
real,    dimension(:,:,:),   intent(in)  :: dx
real,    dimension(:,:,:),   intent(out) :: answer 
integer,                     intent(in)  :: k_min, k_max

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    tab
!    ix
!    n
!    dx
!    k_min
!    k_max
!
!  intent(out) variables:
!
!    answer
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer    :: i_min, i_max, j_min, j_max
      integer    :: i, j, k
   
!--------------------------------------------------------------------
!  local variables:
!
!    i_min
!    i_max
!    j_min
!    j_max
!    i,j,k
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('rad_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

!-----------------------------------------------------------------
      i_min = lbound(ix,1)
      i_max = ubound(ix,1)
      j_min = lbound(ix,2)
      j_max = ubound(ix,2)

      do k=k_min, k_max
        do j=j_min, j_max
          do i=i_min, i_max
                answer(i,j,k) =                                 &
                                      tab%vae (ix(i,j,k),n) +   &
                            dx(i,j,k)*tab%td(ix(i,j,k),n)
          end do
        end do
      end do

!------------------------------------------------------------------

end subroutine  looktab_type3


!#####################################################################
! <SUBROUTINE NAME="table1_alloc">
!  <OVERVIEW>
!   Allocate the longwave tables.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Allocate the longwave tables based on 2 dimension sizes
!  </DESCRIPTION>
!  <TEMPLATE>
!   call table1_alloc(tab, dim1, dim2)
!  </TEMPLATE>
!  <INOUT NAME="tab" TYPE="longwave_tables1_type">
!   The longwave tables
!  </INOUT>
!  <IN NAME="dim1" TYPE="integer">
!   size of the x dimension
!  </IN>
!  <IN NAME="dim2" TYPE="integer">
!   size of the y dimension
!  </IN>
! </SUBROUTINE>
!
subroutine table1_alloc (tab, dim1, dim2)

!------------------------------------------------------------------
!    table1_alloc allocates the arrays contained in a 
!    longwave_tables1_type variable.
!------------------------------------------------------------------

type(longwave_tables1_type), intent (inout) :: tab
integer,                     intent(in)     :: dim1, dim2

!-------------------------------------------------------------------
! intent(in) variables:
!
!     dim1
!     dim2
!
!  intent(inout) variables:
!
!     tab
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('rad_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    allocate the component arrays of a longwave_tables1_type variable.
!---------------------------------------------------------------------
      allocate (tab%vae(dim1, dim2))
      allocate (tab%td (dim1, dim2))
      allocate (tab%md (dim1, dim2))
      allocate (tab%cd (dim1, dim2))

!---------------------------------------------------------------------

end subroutine table1_alloc


!####################################################################
! <SUBROUTINE NAME="table2_alloc">
!  <OVERVIEW>
!   Allocate the longwave tables.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Allocate the longwave tables based on 3 dimension sizes
!  </DESCRIPTION>
!  <TEMPLATE>
!   call table2_alloc(tab, dim1, dim2)
!  </TEMPLATE>
!  <INOUT NAME="tab" TYPE="longwave_tables2_type">
!   The longwave tables
!  </INOUT>
!  <IN NAME="dim1" TYPE="integer">
!   size of the x dimension
!  </IN>
!  <IN NAME="dim2" TYPE="integer">
!   size of the y dimension
!  </IN>
!  <IN NAME="dim3" TYPE="integer">
!   size of the z dimension
!  </IN>
! </SUBROUTINE>
!
subroutine table2_alloc (tab, dim1, dim2, dim3)

!------------------------------------------------------------------
!    table2_alloc allocates the arrays contained in a 
!    longwave_tables2_type variable.
!------------------------------------------------------------------

type(longwave_tables2_type), intent (inout) :: tab
integer,                     intent(in)     :: dim1, dim2, dim3

!-------------------------------------------------------------------
! intent(in) variables:
!
!     dim1
!     dim2
!     dim3
!
!  intent(inout) variables:
!
!     tab
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('rad_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    allocate the component arrays of a longwave_tables2_type variable.
!---------------------------------------------------------------------
      allocate (tab%vae(dim1, dim2, dim3))
      allocate (tab%td (dim1, dim2, dim3))
      allocate (tab%md (dim1, dim2, dim3))
      allocate (tab%cd (dim1, dim2, dim3))

!--------------------------------------------------------------------

end subroutine table2_alloc


!#####################################################################
! <SUBROUTINE NAME="table3_alloc">
!  <OVERVIEW>
!   Allocate the longwave tables.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Allocate the longwave tables based on 2 dimension sizes
!  </DESCRIPTION>
!  <TEMPLATE>
!   call table3_alloc(tab, dim1, dim2)
!  </TEMPLATE>
!  <INOUT NAME="tab" TYPE="longwave_tables3_type">
!   The longwave tables
!  </INOUT>
!  <IN NAME="dim1" TYPE="integer">
!   size of the x dimension
!  </IN>
!  <IN NAME="dim2" TYPE="integer">
!   size of the y dimension
!  </IN>
! </SUBROUTINE>
!
subroutine table3_alloc (tab, dim1, dim2)

type(longwave_tables3_type), intent (inout) :: tab
integer,                     intent(in)     :: dim1, dim2

!-------------------------------------------------------------------
! intent(in) variables:
!
!     dim1
!     dim2
!
!  intent(inout) variables:
!
!     tab
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('rad_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    allocate the component arrays of a longwave_tables3_type variable.
!---------------------------------------------------------------------
      allocate (tab%vae(dim1, dim2))
      allocate (tab%td (dim1, dim2))

end subroutine table3_alloc


!##################################################################

! <SUBROUTINE NAME="thickavg_3d">
!  <OVERVIEW>
!   Subroutine to use thick-averaging technique to define band interval
!   single scattering properties.
!  </OVERVIEW>
!  <DESCRIPTION>
! use the thick-averaging technique to define the single-scattering    
! properties of the parameterization band spectral intervals from the  
! specified spectral intervals of the particular scatterer.            
!                                                                      
! references:                                                          
!                                                                      
! edwards,j.m. and a. slingo, studies with a flexible new radiation    
!      code I: choosing a configuration for a large-scale model.,      
!      q.j.r. meteorological society, 122, 689-719, 1996.              
!                                                                      
! note: the 1.0E-100 factor to calculate asymmband is to prevent        
!       division by zero.                                              
!  </DESCRIPTION>
!  <TEMPLATE>
!   call subroutine thickavg_3d (nivl1    , nivl2     , nivls   ,   &
!                        nbands, $
!                        extivl   , ssalbivl  , asymmivl, solflxivl, &
!                        solflxband, extband  , ssalbband , asymmband)
!  </TEMPLATE>
!  <IN NAME="nivl1" TYPE="integer">
!   interval number for the specified single-scattering                
!              properties corresponding to the first psuedo-           
!              monochromatic frequency in a given parameterization     
!              band  
!  </IN>
!  <IN NAME="nivl2" TYPE="integer">
!   interval number for the specified single-scattering     
!              properties corresponding to the last psuedo-            
!              monochromatic frequency in a given parameterization     
!              band
!  </IN>
!  <IN NAME="nivls" TYPE="integer">
!   number of specified scattering spectral intervals
!  </IN>
!  <IN NAME="nbands" TYPE="integer">
!   number of spectral bands
!  </IN>
!  <IN NAME="extivl" TYPE="real">
!   the specified spectral values of the extinction coefficient 
!  </IN>
!  <INOUT NAME="ssalbivl" TYPE="real">
!   the specified spectral values of the single-scattering albedo
!  </INOUT>
!  <IN NAME="asymmivl" TYPE="real">
!   the specified spectral values of the asymmetry factor
!  </IN>
!  <IN NAME="solflxivl" TYPE="real">
!   the solar flux in each specified scattering spectral interval
!  </IN>
!  <IN NAME="solflxband" TYPE="real">
!   the solar flux in each parameterization band
!  </IN>
!  <OUT NAME="extband" TYPE="real">
!   the parameterization band values of the extinction coefficient
!  </OUT>
!  <OUT NAME="ssalbband" TYPE="real">
!   the parameterization band values of the single-scattering albedo
!  </OUT>
!  <OUT NAME="asymmband" TYPE="real">
!   the parameterization band values of the asymmetry factor
!  </OUT>
! </SUBROUTINE>
!
subroutine thickavg_3d (nivl1, nivl2, nivls, nbands, extivl, ssalbivl,&
                        asymmivl, solflxivl, solflxband, mask, extband,&
                        ssalbband, asymmband)
 
!---------------------------------------------------------------------
!    thickavg_3d uses the thick-averaging technique to define the 
!    single-scattering properties of the parameterization band spectral
!    intervals from the  specified spectral intervals of the particular
!    scatterer, using 3d input arrays.   
!    references:                                                       
!    edwards,j.m. and a. slingo, studies with a flexible new radiation  
!      code I: choosing a configuration for a large-scale model.,   
!      q.j.r. meteorological society, 122, 689-719, 1996.            
!--------------------------------------------------------------------

integer, dimension(:),    intent(in)       :: nivl1, nivl2
integer,                  intent(in)       :: nivls
integer,                  intent(in)       :: nbands
real, dimension(:,:,:,:), intent(in)       :: extivl, asymmivl
real, dimension(:,:,:,:), intent(inout)    :: ssalbivl
real, dimension(:,:),     intent(in)       :: solflxivl             
real, dimension(:),       intent(in)       :: solflxband            
real, dimension(:,:,:,:), intent(out)      :: extband, ssalbband,   &
                                              asymmband
logical, dimension(:,:,:), intent(in)      :: mask

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    nivl1       interval number for the specified single-scattering  
!                properties corresponding to the first psuedo-         
!                monochromatic frequency in a given parameterization    
!                band                                                  
!    nivl2       interval number for the specified single-scattering 
!                properties corresponding to the last psuedo-          
!                monochromatic frequency in a given parameterization    
!                band                                                 
!    nivls       number of specified scattering spectral intervals      
!    nbands
!    extivl      specified spectral values of the extinction coefficient
!    asymmivl    the specified spectral values of the asymmetry     
!                factor                                           
!    solflxivl   the solar flux in each specified scattering spectral
!                interval                                         
!    solflxband  the solar flux in each parameterization band  
!
!  intent(inout) variables:
!
!    ssalbivl    the specified spectral values of the single-       
!                scattering albedo                                   
!
!  intent(out) variables:
!
!    extband     the parameterization band values of the extinction 
!                coefficient                                      
!    ssalbband   the parameterization band values of the single-   
!                scattering albedo                                  
!    asymmband   the parameterization band values of the asymmetry   
!                factor                                               
!    
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!  local variables:
 
      real, dimension (size(ssalbivl,1),   &
                       size(ssalbivl,2), &
                       size(ssalbivl,3), nbands)  ::   refband        

      real, dimension (size(ssalbivl,1),   &
                       size(ssalbivl,2), &
                       size(ssalbivl,3))  ::   refthick, sp, sumk,   &
                                               sumomegak, sumomegakg, &
                                               sumrefthick

      integer  :: nband
      integer  :: i, j, k, ni
 
!--------------------------------------------------------------------
!  local variables:
!
!     refband
!     refthick
!     sp
!     sumk
!     sumomegak
!     sumomegakg
!     sumrefthck
!     nband
!     i,j,k,ni
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('rad_utilities_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!------------------------------------------------------ --------------
!--------------------------------------------------------------------
      do nband = 1,nbands
        sumk(:,:,:) = 0.0
        sumomegak(:,:,:) = 0.0
        sumomegakg(:,:,:) = 0.0
        sumrefthick(:,:,:) = 0.0
        do ni = nivl1(nband),nivl2(nband)
!
          do k=1, size(ssalbivl,3)
            do j=1,size(ssalbivl,2)
              do i=1,size(ssalbivl,1)
                if (mask(i,j,k)) then
                  ssalbivl(i,j,k,ni) = MIN(ssalbivl(i,j,k,ni), 1.0)
                  sp(i,j,k) = sqrt( ( 1.0 - ssalbivl(i,j,k,ni) ) /    &
                                    ( 1.0 - ssalbivl(i,j,k,ni) *      &
                                      asymmivl(i,j,k,ni) ) )
                  refthick(i,j,k) = (1.0 - sp(i,j,k))/(1.0 + sp(i,j,k))
                  sumrefthick(i,j,k) = sumrefthick(i,j,k) +    &
                                       refthick(i,j,k)*  &
                                       solflxivl(nband,ni)
                  sumk(i,j,k) = sumk(i,j,k) + extivl(i,j,k,ni) *   &
                                solflxivl(nband,ni)
                  sumomegak(i,j,k) = sumomegak(i,j,k) +     &
                                     ssalbivl(i,j,k,ni)*   &
                                     extivl(i,j,k,ni) *   &
                                     solflxivl(nband,ni)
                  sumomegakg(i,j,k) = sumomegakg(i,j,k) +    &
                                      ssalbivl(i,j,k,ni)*&
                                      extivl(i,j,k,ni)*  &
                                      asymmivl(i,j,k,ni) * &
                                      solflxivl(nband,ni)
                endif
              end do
            end do
          end do
        end do

!---------------------------------------------------------------------
!    the 1.0E-100 factor to calculate asymmband is to prevent        
!    division by zero.                                             
!---------------------------------------------------------------------
        do k=1, size(ssalbivl,3)
          do j=1,size(ssalbivl,2)
            do i=1,size(ssalbivl,1)
              extband(i,j,k,nband) = sumk(i,j,k) / solflxband(nband)
              asymmband(i,j,k,nband) = sumomegakg(i,j,k) /         &
                                       ( sumomegak(i,j,k) + 1.0E-100)
              refband(i,j,k,nband) = sumrefthick(i,j,k)/  &
                                     solflxband(nband)
              ssalbband(i,j,k,nband) = 4.0 * refband(i,j,k,nband) / &
                                       ((1.0 +    &
                                       refband(i,j,k,nband)) ** 2 -&
                                       asymmband(i,j,k,nband) *     &
                                       (1.0 - refband(i,j,k,nband))**2 )
            end do
          end do
        end do
      end do

!---------------------------------------------------------------------

  
end subroutine thickavg_3d



!####################################################################
! <SUBROUTINE NAME="thickavg_0d">
!  <OVERVIEW>
!   Subroutine to use thick-averaging technique to define band interval
!   single scattering properties.
!  </OVERVIEW>
!  <DESCRIPTION>
! use the thick-averaging technique to define the single-scattering    
! properties of the parameterization band spectral intervals from the  
! specified spectral intervals of the particular scatterer.            
!                                                                      
! references:                                                          
!                                                                      
! edwards,j.m. and a. slingo, studies with a flexible new radiation    
!      code I: choosing a configuration for a large-scale model.,      
!      q.j.r. meteorological society, 122, 689-719, 1996.              
!                                                                      
! note: the 1.0E-100 factor to calculate asymmband is to prevent        
!       division by zero.                                              
!  </DESCRIPTION>
!  <TEMPLATE>
!   call subroutine thickavg_0d (nivl1    , nivl2     , nivls   ,   &
!                        nbands,  &
!                        extivl   , ssalbivl  , asymmivl, solflxivl, &
!                        solflxband, extband  , ssalbband , asymmband)
!  </TEMPLATE>
!  <IN NAME="nivl1" TYPE="integer">
!   interval number for the specified single-scattering                
!              properties corresponding to the first psuedo-           
!              monochromatic frequency in a given parameterization     
!              band  
!  </IN>
!  <IN NAME="nivl2" TYPE="integer">
!   interval number for the specified single-scattering     
!              properties corresponding to the last psuedo-            
!              monochromatic frequency in a given parameterization     
!              band
!  </IN>
!  <IN NAME="nivls" TYPE="integer">
!   number of specified scattering spectral intervals
!  </IN>
!  <IN NAME="nbands" TYPE="integer">
!   number of spectral bands
!  </IN>
!  <IN NAME="extivl" TYPE="real">
!   the specified spectral values of the extinction coefficient 
!  </IN>
!  <INOUT NAME="ssalbivl" TYPE="real">
!   the specified spectral values of the single-scattering albedo
!  </INOUT>
!  <IN NAME="asymmivl" TYPE="real">
!   the specified spectral values of the asymmetry factor
!  </IN>
!  <IN NAME="solflxivl" TYPE="real">
!   the solar flux in each specified scattering spectral interval
!  </IN>
!  <IN NAME="solflxband" TYPE="real">
!   the solar flux in each parameterization band
!  </IN>
!  <OUT NAME="extband" TYPE="real">
!   the parameterization band values of the extinction coefficient
!  </OUT>
!  <OUT NAME="ssalbband" TYPE="real">
!   the parameterization band values of the single-scattering albedo
!  </OUT>
!  <OUT NAME="asymmband" TYPE="real">
!   the parameterization band values of the asymmetry factor
!  </OUT>
! </SUBROUTINE>
!
subroutine thickavg_0d (nivl1, nivl2, nivls, nbands, extivl, ssalbivl,&
                        asymmivl, solflxivl, solflxband, extband,  &
                        ssalbband , asymmband)
 
!---------------------------------------------------------------------
!    thickavg_0d uses the thick-averaging technique to define the 
!    single-scattering properties of the parameterization band spectral
!    intervals from the  specified spectral intervals of the particular
!    scatterer, using 3d input arrays.   
!    references:                                                       
!    edwards,j.m. and a. slingo, studies with a flexible new radiation  
!      code I: choosing a configuration for a large-scale model.,   
!      q.j.r. meteorological society, 122, 689-719, 1996.            
!--------------------------------------------------------------------

integer, dimension(:),    intent(in)       :: nivl1, nivl2
integer,                  intent(in)       :: nivls
integer,                  intent(in)       :: nbands
real, dimension(:),       intent(in)       :: extivl, asymmivl
real, dimension(:),       intent(inout)    :: ssalbivl
real, dimension(:,:),     intent(in)       :: solflxivl             
real, dimension(:),       intent(in)       :: solflxband            
real, dimension(:),       intent(out)      :: extband, ssalbband, &
                                              asymmband

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    nivl1       interval number for the specified single-scattering  
!                properties corresponding to the first psuedo-         
!                monochromatic frequency in a given parameterization    
!                band                                                  
!    nivl2       interval number for the specified single-scattering 
!                properties corresponding to the last psuedo-          
!                monochromatic frequency in a given parameterization    
!                band                                                 
!    nivls       number of specified scattering spectral intervals      
!    nbands
!    extivl      specified spectral values of the extinction coefficient
!    asymmivl    the specified spectral values of the asymmetry     
!                factor                                           
!    solflxivl   the solar flux in each specified scattering spectral
!                interval                                         
!    solflxband  the solar flux in each parameterization band  
!
!  intent(inout) variables:
!
!    ssalbivl    the specified spectral values of the single-       
!                scattering albedo                                   
!
!  intent(out) variables:
!
!    extband     the parameterization band values of the extinction 
!                coefficient                                      
!    ssalbband   the parameterization band values of the single-   
!                scattering albedo                                  
!    asymmband   the parameterization band values of the asymmetry   
!                factor                                               
!    
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!  local variables:

      real, dimension (nbands)        ::   refband
      real                            ::   refthick, sp, sumk,   &
                                           sumomegak, sumomegakg, &
                                           sumrefthick
      integer  :: nband 
      integer  :: ni
 
 
!--------------------------------------------------------------------
!  local variables:
!
!     refband
!     refthick
!     sp
!     sumk
!     sumomegak
!     sumomegakg
!     sumrefthck
!     nband
!     i,j,k,ni
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('rad_utilities_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif
  
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      do nband = 1,nbands
        sumk        = 0.0
        sumomegak   = 0.0
        sumomegakg  = 0.0
        sumrefthick = 0.0
        do ni = nivl1(nband),nivl2(nband)
          if (extivl(ni) /= 0.0) then
            ssalbivl(ni) = MIN(ssalbivl(ni), 1.0)
            sp = sqrt( ( 1.0 - ssalbivl(ni) ) /    &
                       ( 1.0 - ssalbivl(ni) * asymmivl(ni) ) )
            refthick = (1.0 - sp)/(1.0 + sp)
            sumrefthick = sumrefthick + refthick * solflxivl(nband,ni)
            sumk = sumk + extivl(ni) * solflxivl(nband,ni)
            sumomegak = sumomegak +     &
                        ssalbivl(ni) * extivl(ni) * solflxivl(nband,ni)
            sumomegakg = sumomegakg +    &
                         ssalbivl(ni) * extivl(ni) *  &
                         asymmivl(ni) * solflxivl(nband,ni)
          endif
        end do

!--------------------------------------------------------------------- 
!
!--------------------------------------------------------------------- 
        extband(nband) = sumk / solflxband(nband)
        asymmband(nband) = sumomegakg / ( sumomegak + 1.0E-100)
        refband(nband) = sumrefthick/ solflxband(nband)
        ssalbband(nband) = 4.0 * refband(nband) / &
                           ( (1.0 + refband(nband))**2 - &
                          asymmband(nband) * (1.0 - refband(nband))**2 )
      end do

!---------------------------------------------------------------------
  

end subroutine thickavg_0d


!##################################################################

! <SUBROUTINE NAME="thickavg_isccp">
!  <OVERVIEW>
!   Subroutine to use thick-averaging technique to define band interval
!   single scattering properties for a single specified band.
!  </OVERVIEW>
!  <DESCRIPTION>
! use the thick-averaging technique to define the single-scattering    
! properties of the specified parameterization band spectral interval 
! from the specified spectral intervals of the particular scatterer.    
!                                                                      
! references:                                                          
!                                                                      
! edwards,j.m. and a. slingo, studies with a flexible new radiation    
!      code I: choosing a configuration for a large-scale model.,      
!      q.j.r. meteorological society, 122, 689-719, 1996.              
!                                                                      
! note: the 1.0E-100 factor to calculate asymmband is to prevent        
!       division by zero.                                              
!  </DESCRIPTION>
!  <TEMPLATE>
!   call subroutine thickavg (nband, nivl1, nivl2, extivl, solflxivl, &
!                             solflxband, mask, extband )
!  </TEMPLATE>
!  <IN NAME="nband" TYPE="integer">
!
!  </IN>
!  <IN NAME="nivl1" TYPE="integer">
!   interval number for the specified single-scattering                
!              properties corresponding to the first psuedo-           
!              monochromatic frequency in a given parameterization     
!              band  
!  </IN>
!  <IN NAME="nivl2" TYPE="integer">
!   interval number for the specified single-scattering     
!              properties corresponding to the last psuedo-            
!              monochromatic frequency in a given parameterization     
!              band
!  </IN>
!  <IN NAME="extivl" TYPE="real">
!   the specified spectral values of the extinction coefficient 
!  </IN>
!  <IN NAME="solflxivl" TYPE="real">
!   the solar flux in each specified scattering spectral interval
!  </IN>
!  <IN NAME="solflxband" TYPE="real">
!   the solar flux in each parameterization band
!  </IN>
!  <IN NAME="mask" TYPE="logical">
!   mask is .true. at gridpoints where extband needs to be calculated
!  </IN>
!  <OUT NAME="extband" TYPE="real">
!   the parameterization band values of the extinction coefficient
!  </OUT>
! </SUBROUTINE>
!
subroutine thickavg_isccp (nband, nivl1, nivl2, extivl,          &
                           solflxivl, solflxband, mask, extband)
 
!---------------------------------------------------------------------
!    thickavg_isccp uses the thick-averaging technique to define the 
!    solar extinction for the single specified parameterization band 
!    spectral interval (nband) from the  specified spectral intervals 
!    of the particular scatterer, using 3d input arrays.   
!    references:                                                       
!    edwards,j.m. and a. slingo, studies with a flexible new radiation  
!      code I: choosing a configuration for a large-scale model.,   
!      q.j.r. meteorological society, 122, 689-719, 1996.            
!--------------------------------------------------------------------

integer,                  intent(in)       :: nband
integer,                  intent(in)       :: nivl1, nivl2
real, dimension(:,:,:,:), intent(in)       :: extivl
real, dimension(:,:),     intent(in)       :: solflxivl             
real,                     intent(in)       :: solflxband            
logical, dimension(:,:,:),intent(in)       :: mask
real, dimension(:,:,:),   intent(out)      :: extband

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    nband       the sw parameterization band for which the optical
!                properties are being calculated
!    nivl1       interval number for the specified single-scattering  
!                properties corresponding to the first psuedo-         
!                monochromatic frequency in a given parameterization    
!                band                                                  
!    nivl2       interval number for the specified single-scattering 
!                properties corresponding to the last psuedo-          
!                monochromatic frequency in a given parameterization    
!                band                                                 
!    extivl      specified spectral values of the extinction coefficient
!    solflxband  the solar flux in each parameterization band  
!    mask        logical indicating the points at which the band values 
!                should be calculated       
!
!  intent(out) variables:
!
!    extband     the parameterization band values of the extinction 
!                coefficient                                      
!    
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!  local variables:
 
      real     ::  sumk
      integer  ::  i, j, k, ni
 
!--------------------------------------------------------------------
!  local variables:
!
!     sumk
!     i,j,k,ni
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('rad_utilities_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!--------------------------------------------------------------------
!
      do k=1, size(extivl,3)
        do j=1,size(extivl,2)
          do i=1,size(extivl,1)
            if (mask(i,j,k)) then
              sumk = 0.0
              do ni = nivl1,nivl2
                sumk = sumk + extivl(i,j,k,ni)*solflxivl(nband,ni)
              end do
              extband(i,j,k) = sumk/solflxband
            endif
          end do
        end do
      end do

!---------------------------------------------------------------------

  
end subroutine thickavg_isccp

!##################################################################

! <SUBROUTINE NAME="thickavg_1band">
!  <OVERVIEW>
!   Subroutine to use thick-averaging technique to define band interval
!   single scattering properties for a single specified band.
!  </OVERVIEW>
!  <DESCRIPTION>
! use the thick-averaging technique to define the single-scattering    
! properties of the specified parameterization band spectral interval 
! from the specified spectral intervals of the particular scatterer.    
!                                                                      
! references:                                                          
!                                                                      
! edwards,j.m. and a. slingo, studies with a flexible new radiation    
!      code I: choosing a configuration for a large-scale model.,      
!      q.j.r. meteorological society, 122, 689-719, 1996.              
!                                                                      
! note: the 1.0E-100 factor to calculate asymmband is to prevent        
!       division by zero.                                              
!  </DESCRIPTION>
!  <TEMPLATE>
!   call subroutine thickavg (nband, nivl1  , nivl2   , nivls   , &
!                        nbands, &
!                        extivl   , ssalbivl  , asymmivl, solflxivl, &
!                        solflxband,  mask, extband  , ssalbband ,  &
!                        asymmband)
!  </TEMPLATE>
!  <IN NAME="nband" TYPE="integer">
!
!  </IN>
!  <IN NAME="nivl1" TYPE="integer">
!   interval number for the specified single-scattering                
!              properties corresponding to the first psuedo-           
!              monochromatic frequency in a given parameterization     
!              band  
!  </IN>
!  <IN NAME="nivl2" TYPE="integer">
!   interval number for the specified single-scattering     
!              properties corresponding to the last psuedo-            
!              monochromatic frequency in a given parameterization     
!              band
!  </IN>
!  <IN NAME="nivls" TYPE="integer">
!   number of specified scattering spectral intervals
!  </IN>
!  <IN NAME="nbands" TYPE="integer">
!   number of spectral bands
!  </IN>
!  <IN NAME="extivl" TYPE="real">
!   the specified spectral values of the extinction coefficient 
!  </IN>
!  <INOUT NAME="ssalbivl" TYPE="real">
!   the specified spectral values of the single-scattering albedo
!  </INOUT>
!  <IN NAME="asymmivl" TYPE="real">
!   the specified spectral values of the asymmetry factor
!  </IN>
!  <IN NAME="solflxivl" TYPE="real">
!   the solar flux in each specified scattering spectral interval
!  </IN>
!  <IN NAME="solflxband" TYPE="real">
!   the solar flux in each parameterization band
!  </IN>
!  <IN NAME="mask" TYPE="logical">
!   mask is .true. at gridpoints where band calculations are needed
!  </IN>
!  <OUT NAME="extband" TYPE="real">
!   the parameterization band values of the extinction coefficient
!  </OUT>
!  <OUT NAME="ssalbband" TYPE="real">
!   the parameterization band values of the single-scattering albedo
!  </OUT>
!  <OUT NAME="asymmband" TYPE="real">
!   the parameterization band values of the asymmetry factor
!  </OUT>
! </SUBROUTINE>
!
subroutine thickavg_1band (nband, nivl1, nivl2, nivls, nbands, extivl, &
                           ssalbivl, asymmivl, solflxivl, solflxband, &
                           mask, extband, ssalbband, asymmband)
 
!---------------------------------------------------------------------
!    thickavg_1band uses the thick-averaging technique to define the 
!    single-scattering properties of the specified parameterization band
!    spectral interval from the  specified spectral intervals of the 
!    particular scatterer, using 3d input arrays.   
!    references:                                                       
!    edwards,j.m. and a. slingo, studies with a flexible new radiation  
!      code I: choosing a configuration for a large-scale model.,   
!      q.j.r. meteorological society, 122, 689-719, 1996.            
!--------------------------------------------------------------------

integer,                  intent(in)       :: nband
integer,                  intent(in)       :: nivl1, nivl2
integer,                  intent(in)       :: nivls
integer,                  intent(in)       :: nbands
real, dimension(:,:,:,:), intent(in)       :: extivl, asymmivl
real, dimension(:,:,:,:), intent(inout)    :: ssalbivl
real, dimension(:,:),     intent(in)       :: solflxivl             
real,                     intent(in)       :: solflxband            
real, dimension(:,:,:  ), intent(inout)      :: extband, ssalbband,   &
                                              asymmband
logical, dimension(:,:,:), intent(in)      :: mask

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    nband       the sw parameterization band for which the optical
!                properties are being calculated
!    nivl1       interval number for the specified single-scattering  
!                properties corresponding to the first psuedo-         
!                monochromatic frequency in a given parameterization    
!                band                                                  
!    nivl2       interval number for the specified single-scattering 
!                properties corresponding to the last psuedo-          
!                monochromatic frequency in a given parameterization    
!                band                                                 
!    nivls       number of specified scattering spectral intervals      
!    nbands
!    extivl      specified spectral values of the extinction coefficient
!    asymmivl    the specified spectral values of the asymmetry     
!                factor                                           
!    solflxivl   the solar flux in each specified scattering spectral
!                interval                                         
!    solflxband  the solar flux in each parameterization band  
!    mask        logical indicating the points at which the band values 
!                should be calculated       
!
!  intent(inout) variables:
!
!    ssalbivl    the specified spectral values of the single-       
!                scattering albedo                                   
!
!  intent(out) variables:
!
!    extband     the parameterization band values of the extinction 
!                coefficient                                      
!    ssalbband   the parameterization band values of the single-   
!                scattering albedo                                  
!    asymmband   the parameterization band values of the asymmetry   
!                factor                                               
!    
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!  local variables:
 
      real :: refband, sp, refthick
      real :: sumk, sumomegak, sumomegakg,  sumrefthick

      integer  :: i, j, k, ni
 
!--------------------------------------------------------------------
!  local variables:
!
!     refband
!     refthick
!     sp
!     sumk
!     sumomegak
!     sumomegakg
!     sumrefthck
!     nband
!     i,j,k,ni
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('rad_utilities_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      do k=1, size(ssalbivl,3)
        do j=1,size(ssalbivl,2)
          do i=1,size(ssalbivl,1)
            if (mask(i,j,k)) then
              sumk        = 0.0
              sumomegak        = 0.0
              sumomegakg        = 0.0
              sumrefthick        = 0.0
              do ni = nivl1,nivl2
                ssalbivl(i,j,k,ni) = MIN(ssalbivl(i,j,k,ni), 1.0)
                sp = sqrt((1.0 - ssalbivl(i,j,k,ni) ) /    &
                          (1.0 - ssalbivl(i,j,k,ni)*asymmivl(i,j,k,ni)))
                refthick = (1.0 - sp)/(1.0 + sp)
                sumrefthick = sumrefthick + refthick*solflxivl(nband,ni)
                sumk = sumk + extivl(i,j,k,ni)*solflxivl(nband,ni)
                sumomegak = sumomegak + ssalbivl(i,j,k,ni)*   &
                                        extivl(i,j,k,ni)*   &
                                        solflxivl(nband,ni)
                sumomegakg = sumomegakg + ssalbivl(i,j,k,ni)*&
                                          extivl(i,j,k,ni)*  &
                                          asymmivl(i,j,k,ni)* &
                                          solflxivl(nband,ni)
              end do

!---------------------------------------------------------------------
!    the 1.0E-100 factor to calculate asymmband is to prevent        
!    division by zero.                                             
!---------------------------------------------------------------------
              extband(i,j,k) = sumk/solflxband
              asymmband(i,j,k) = sumomegakg/(sumomegak + 1.0E-100)
              refband  = sumrefthick/solflxband        
              ssalbband(i,j,k) = 4.0*refband/((1.0 + refband) ** 2 - &
                                 asymmband(i,j,k)*(1.0 - refband)**2 )
            endif
          end do
        end do
      end do

!---------------------------------------------------------------------

  
end subroutine thickavg_1band

!####################################################################
! <SUBROUTINE NAME="thinavg">
!  <OVERVIEW>
!   Subroutine to use thin-averaging technique to define band interval
!   single scattering properties.
!  </OVERVIEW>
!  <DESCRIPTION>
! use the thin-averaging technique to define the single-scattering    
! properties of the parameterization band spectral intervals from the  
! specified spectral intervals of the particular scatterer.            
!                                                                      
! references:                                                          
!                                                                      
! edwards,j.m. and a. slingo, studies with a flexible new radiation    
!      code I: choosing a configuration for a large-scale model.,      
!      q.j.r. meteorological society, 122, 689-719, 1996.              
!                                                                      
! note: the 1.0E-100 factor to calculate asymmband is to prevent        
!       division by zero.                                              
!  </DESCRIPTION>
!  <TEMPLATE>
!   call subroutine thinavg (nivl1    , nivl2     , nivls   ,   &
!                        nbands, &
!                        extivl   , ssalbivl  , asymmivl, solflxivl, &
!                        solflxband, extband  , ssalbband , asymmband)
!  </TEMPLATE>
!  <IN NAME="nivl1" TYPE="integer">
!   interval number for the specified single-scattering                
!              properties corresponding to the first psuedo-           
!              monochromatic frequency in a given parameterization     
!              band  
!  </IN>
!  <IN NAME="nivl2" TYPE="integer">
!   interval number for the specified single-scattering     
!              properties corresponding to the last psuedo-            
!              monochromatic frequency in a given parameterization     
!              band
!  </IN>
!  <IN NAME="nivls" TYPE="integer">
!   number of specified scattering spectral intervals
!  </IN>
!  <IN NAME="extivl" TYPE="real">
!   the specified spectral values of the extinction coefficient 
!  </IN>
!  <IN NAME="nbands" TYPE="integer">
!   number of spectral bands
!  </IN>
!  <INOUT NAME="ssalbivl" TYPE="real">
!   the specified spectral values of the single-scattering albedo
!  </INOUT>
!  <IN NAME="asymmivl" TYPE="real">
!   the specified spectral values of the asymmetry factor
!  </IN>
!  <IN NAME="solflxivl" TYPE="real">
!   the solar flux in each specified scattering spectral interval
!  </IN>
!  <IN NAME="solflxband" TYPE="real">
!   the solar flux in each parameterization band
!  </IN>
!  <OUT NAME="extband" TYPE="real">
!   the parameterization band values of the extinction coefficient
!  </OUT>
!  <OUT NAME="ssalbband" TYPE="real">
!   the parameterization band values of the single-scattering albedo
!  </OUT>
!  <OUT NAME="asymmband" TYPE="real">
!   the parameterization band values of the asymmetry factor
!  </OUT>
! </SUBROUTINE>
!
subroutine thinavg (nivl1, nivl2, nivls, nbands, extivl, ssalbivl, &
                    asymmivl,  solflxivl, solflxband, extband,   &
                    ssalbband , asymmband)
 
!---------------------------------------------------------------------
!    thinavg uses the thin-averaging technique to define the 
!    single-scattering properties of the parameterization band spectral
!    intervals from the  specified spectral intervals of the particular
!    scatterer, using 3d input arrays.   
!    references:                                                       
!    edwards,j.m. and a. slingo, studies with a flexible new radiation  
!      code I: choosing a configuration for a large-scale model.,   
!      q.j.r. meteorological society, 122, 689-719, 1996.            
!--------------------------------------------------------------------

integer, dimension(:),    intent(in)       :: nivl1, nivl2
integer,                  intent(in)       :: nivls
integer,                  intent(in)       :: nbands
real, dimension(:,:,:,:), intent(in)       :: extivl, asymmivl
real, dimension(:,:,:,:), intent(inout)    :: ssalbivl
real, dimension(:,:),     intent(in)       :: solflxivl             
real, dimension(:),       intent(in)       :: solflxband            
real, dimension(:,:,:,:), intent(out)      :: extband, ssalbband,   &
                                              asymmband

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    nivl1       interval number for the specified single-scattering  
!                properties corresponding to the first psuedo-         
!                monochromatic frequency in a given parameterization    
!                band                                                  
!    nivl2       interval number for the specified single-scattering 
!                properties corresponding to the last psuedo-          
!                monochromatic frequency in a given parameterization    
!                band                                                 
!    nivls       number of specified scattering spectral intervals      
!    nbands
!    extivl      specified spectral values of the extinction coefficient
!    asymmivl    the specified spectral values of the asymmetry     
!                factor                                           
!    solflxivl   the solar flux in each specified scattering spectral
!                interval                                         
!    solflxband  the solar flux in each parameterization band  
!
!  intent(inout) variables:
!
!    ssalbivl    the specified spectral values of the single-       
!                scattering albedo                                   
!
!  intent(out) variables:
!
!    extband     the parameterization band values of the extinction 
!                coefficient                                      
!    ssalbband   the parameterization band values of the single-   
!                scattering albedo                                  
!    asymmband   the parameterization band values of the asymmetry   
!                factor                                               
!    
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      real, dimension (size(ssalbivl,1),   &
                       size(ssalbivl,2),  &
                       size(ssalbivl,3)) ::   sumk,  sumomegak,   &
                                              sumomegakg
 
      integer   ::   nband
      integer   ::   i, j, k, ni

!--------------------------------------------------------------------
!  local variables:
! 
!    sumk
!    sumomegak
!    sumomegakg
!    nband
!    i,j,k,ni
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('rad_utilities_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do nband = 1,nbands
        sumk(:,:,:) = 0.0
        sumomegak(:,:,:) = 0.0
        sumomegakg(:,:,:) = 0.0
        do ni = nivl1(nband),nivl2(nband)
          do k=1, size(ssalbivl,3)
            do j=1,size(ssalbivl,2)
              do i=1,size(ssalbivl,1)
                if ((ssalbivl(i,j,k,ni) +    &
                     asymmivl(i,j,k,ni)) /= 0.0) then
                  ssalbivl(i,j,k,ni) = MIN(ssalbivl(i,j,k,ni), 1.0)
                  sumk(i,j,k) = sumk(i,j,k) + extivl(i,j,k,ni) *   &
                                solflxivl(nband,ni)
                  sumomegak(i,j,k) = sumomegak(i,j,k) +    &
                                     ssalbivl(i,j,k,ni) *  &
                                     extivl(i,j,k,ni) *   &
                                     solflxivl(nband,ni)
                  sumomegakg(i,j,k) = sumomegakg(i,j,k) +    &
                                      ssalbivl(i,j,k,ni) * & 
                                      extivl(i,j,k,ni) *   &
                                      asymmivl(i,j,k,ni) *  &
                                      solflxivl(nband,ni)
                endif
              end do
            end do
          end do
        end do

!----------------------------------------------------------------------
!
!---------------------------------------------------------------------
        do k=1, size(ssalbivl,3)
          do j=1,size(ssalbivl,2)
            do i=1,size(ssalbivl,1)
              extband(i,j,k,nband) = sumk(i,j,k) / solflxband(nband)
              asymmband(i,j,k,nband) = sumomegakg(i,j,k) /    &
                                       ( sumomegak(i,j,k) + 1.0E-100 )
              ssalbband(i,j,k,nband) = sumomegak(i,j,k) /   &
                                       ( sumk(i,j,k) + 1.0E-100 )
            end do
          end do
        end do
      end do

!-------------------------------------------------------------------
  

end subroutine thinavg 


!#########################################################################
subroutine get_radiative_param(text_in_scheme,text_in_param, &
                               rad_forc_online, tr_rad_name,  &
                               tr_clim_name, tr_rad_scale_factor)

character(len=*), intent(in)    :: text_in_scheme, text_in_param
logical, intent(out)            :: rad_forc_online
character(len=*), intent(out)   :: tr_rad_name,tr_clim_name
real,             intent(out)   :: tr_rad_scale_factor
integer                         :: flag


if(lowercase(trim(text_in_scheme(1:6))) == 'online') then
       rad_forc_online = .true.
       flag=parse(text_in_param,'name_in_rad_mod', tr_rad_name)
       flag=parse(text_in_param,'name_in_clim_mod', tr_clim_name)
       tr_rad_scale_factor = 1.
       flag=parse(text_in_param,'scale_factor', tr_rad_scale_factor)
else
       rad_forc_online = .false.
       tr_rad_name  = ' '
       tr_clim_name = ' '
       tr_rad_scale_factor = 1.
endif

end subroutine get_radiative_param


!#####################################################################

! <SUBROUTINE NAME="rad_utilities_end">
!  <OVERVIEW>
!   Subroutine to close out the radiation utility package.
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine is the destructor for rad_utilies_mod. it marks
!   the module as uninitialized.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call rad_utilities_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine rad_utilities_end

!--------------------------------------------------------------------
!    rad_utilites_end is the destructor for rad_utilities_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('rad_utilites_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
       module_is_initialized = .false.

!-------------------------------------------------------------------


end subroutine rad_utilities_end


subroutine aerosol_props_type_eq(aerosol_props_out,aerosol_props_in)

   type(aerosol_properties_type), intent(inout) :: aerosol_props_out
   type(aerosol_properties_type), intent(in)    :: aerosol_props_in

!  Need to add error trap to catch unallocated aerosol_props_in
   if (ASSOCIATED(aerosol_props_in%aerextband)) then
     aerosol_props_out%aerextband    = aerosol_props_in%aerextband
     aerosol_props_out%aerssalbband    = aerosol_props_in%aerssalbband
     aerosol_props_out%aerasymmband    = aerosol_props_in%aerasymmband
    else
      call error_mesg ('=', 'extband', FATAL)
   endif
   if (ASSOCIATED(aerosol_props_in%aerextbandlw)) then
     aerosol_props_out%aerextbandlw    = aerosol_props_in%aerextbandlw
     aerosol_props_out%aerssalbbandlw  = aerosol_props_in%aerssalbbandlw
     aerosol_props_out%aerextbandlw_cn =    &
                                       aerosol_props_in%aerextbandlw_cn
     aerosol_props_out%aerssalbbandlw_cn  =    &
                                     aerosol_props_in%aerssalbbandlw_cn
    else
      call error_mesg ('=', 'extbandlw', FATAL)
   endif
  if (Rad_control%volcanic_sw_aerosols) then
   if (ASSOCIATED(aerosol_props_in%sw_ext)) then
     aerosol_props_out%sw_ext        = aerosol_props_in%sw_ext     
     aerosol_props_out%sw_ssa          = aerosol_props_in%sw_ssa       
     aerosol_props_out%sw_asy          = aerosol_props_in%sw_asy
    else
      call error_mesg ('=', 'sw volc', FATAL)
   endif
  endif
  if (Rad_control%volcanic_lw_aerosols) then
   if (ASSOCIATED(aerosol_props_in%lw_ext)) then
     aerosol_props_out%lw_ext        = aerosol_props_in%lw_ext     
     aerosol_props_out%lw_ssa          = aerosol_props_in%lw_ssa       
     aerosol_props_out%lw_asy          = aerosol_props_in%lw_asy
    else
      call error_mesg ('=', 'lw volc', FATAL)
   endif
  endif
   if (ASSOCIATED(aerosol_props_in%sulfate_index)) then
     aerosol_props_out%sulfate_index = aerosol_props_in%sulfate_index
     aerosol_props_out%optical_index = aerosol_props_in%optical_index
     aerosol_props_out%omphilic_index = aerosol_props_in%omphilic_index
     aerosol_props_out%bcphilic_index = aerosol_props_in%bcphilic_index
     aerosol_props_out%seasalt1_index = aerosol_props_in%seasalt1_index
     aerosol_props_out%seasalt2_index = aerosol_props_in%seasalt2_index
     aerosol_props_out%seasalt3_index = aerosol_props_in%seasalt3_index
     aerosol_props_out%seasalt4_index = aerosol_props_in%seasalt4_index
     aerosol_props_out%seasalt5_index = aerosol_props_in%seasalt5_index
    else
      call error_mesg ('=', 'index  ', FATAL)
   endif

   if (ASSOCIATED(aerosol_props_in%ivol)) then
     aerosol_props_out%ivol = aerosol_props_in%ivol
    else
      call error_mesg ('=', 'ivol   ', FATAL)
   endif
   
     aerosol_props_out%sulfate_flag = aerosol_props_in%sulfate_flag
     aerosol_props_out%omphilic_flag = aerosol_props_in%omphilic_flag
     aerosol_props_out%bcphilic_flag = aerosol_props_in%bcphilic_flag
     aerosol_props_out%seasalt1_flag = aerosol_props_in%seasalt1_flag
     aerosol_props_out%seasalt2_flag = aerosol_props_in%seasalt2_flag
     aerosol_props_out%seasalt3_flag = aerosol_props_in%seasalt3_flag
     aerosol_props_out%seasalt4_flag = aerosol_props_in%seasalt4_flag
     aerosol_props_out%seasalt5_flag = aerosol_props_in%seasalt5_flag
     aerosol_props_out%bc_flag = aerosol_props_in%bc_flag



end subroutine aerosol_props_type_eq



subroutine lw_output_type_eq(lw_output_out,lw_output_in)

   type(lw_output_type), intent(inout) :: lw_output_out
   type(lw_output_type), intent(in)    :: lw_output_in

!  Need to add error trap to catch unallocated lw_output_in
   lw_output_out%heatra        = lw_output_in%heatra
   lw_output_out%flxnet        = lw_output_in%flxnet
   lw_output_out%netlw_special = lw_output_in%netlw_special
   lw_output_out%bdy_flx       = lw_output_in%bdy_flx
   if (ASSOCIATED(lw_output_in%heatracf))then
       lw_output_out%heatracf          = lw_output_in%heatracf
       lw_output_out%flxnetcf          = lw_output_in%flxnetcf
       lw_output_out%netlw_special_clr = lw_output_in%netlw_special_clr
       lw_output_out%bdy_flx_clr       = lw_output_in%bdy_flx_clr
   endif
end subroutine lw_output_type_eq


subroutine sw_output_type_eq(sw_output_out,sw_output_in)

   type(sw_output_type), intent(inout) :: sw_output_out
   type(sw_output_type), intent(in)    :: sw_output_in

!  Need to add error trap to catch unallocated sw_output_in
   sw_output_out%fsw              = sw_output_in%fsw
   sw_output_out%dfsw             = sw_output_in%dfsw
   sw_output_out%ufsw             = sw_output_in%ufsw
   sw_output_out%hsw              = sw_output_in%hsw
   sw_output_out%dfsw_dir_sfc     = sw_output_in%dfsw_dir_sfc
   sw_output_out%ufsw_dir_sfc     = sw_output_in%ufsw_dir_sfc
   sw_output_out%dfsw_dif_sfc     = sw_output_in%dfsw_dif_sfc
   sw_output_out%ufsw_dif_sfc     = sw_output_in%ufsw_dif_sfc
   sw_output_out%dfsw_vis_sfc     = sw_output_in%dfsw_vis_sfc
   sw_output_out%ufsw_vis_sfc     = sw_output_in%ufsw_vis_sfc
   sw_output_out%ufsw_vis_sfc_dir = sw_output_in%ufsw_vis_sfc_dir
   sw_output_out%dfsw_vis_sfc_dir = sw_output_in%dfsw_vis_sfc_dir
   sw_output_out%dfsw_vis_sfc_dif = sw_output_in%dfsw_vis_sfc_dif
   sw_output_out%ufsw_vis_sfc_dif = sw_output_in%ufsw_vis_sfc_dif
   sw_output_out%swdn_special     = sw_output_in%swdn_special
   sw_output_out%swup_special     = sw_output_in%swup_special
   sw_output_out%bdy_flx          = sw_output_in%bdy_flx
   if (ASSOCIATED(sw_output_in%fswcf))then
       sw_output_out%fswcf            = sw_output_in%fswcf
       sw_output_out%dfswcf           = sw_output_in%dfswcf
       sw_output_out%ufswcf           = sw_output_in%ufswcf
       sw_output_out%hswcf            = sw_output_in%hswcf
       sw_output_out%dfsw_dir_sfc_clr = sw_output_in%dfsw_dir_sfc_clr
       sw_output_out%dfsw_dif_sfc_clr = sw_output_in%dfsw_dif_sfc_clr
       sw_output_out%dfsw_vis_sfc_clr = sw_output_in%dfsw_vis_sfc_clr
       sw_output_out%swdn_special_clr = sw_output_in%swdn_special_clr
       sw_output_out%swup_special_clr = sw_output_in%swup_special_clr
       sw_output_out%bdy_flx_clr      = sw_output_in%bdy_flx_clr
   endif  
end subroutine sw_output_type_eq


!####################################################################


                     end module rad_utilities_mod



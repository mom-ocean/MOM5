module flux_exchange_mod
!-----------------------------------------------------------------------
!                   GNU General Public License                        !                                                                      
! This program is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! MOM is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Bruce Wyman </CONTACT>
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> V. Balaji </CONTACT>
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Sergey Malyshev </CONTACT>


! <OVERVIEW>
!   The flux_exchange module provides interfaces to couple the following component 
!   models: atmosphere, ocean, land, and ice. All interpolation between physically 
!   distinct model grids is handled by the exchange grid (xgrid_mod) with the 
!   interpolated quantities being conserved.
! </OVERVIEW>

! <DESCRIPTION>
!  <PRE>
!  1.This version of flux_exchange_mod allows the definition of physically independent
!    grids for atmosphere, land and sea ice. Ice and ocean must share the same physical
!    grid (though the domain decomposition on parallel systems may be different). 
!    Grid information is input through the grid_spec file (URL). The masked region of the
!    land grid and ice/ocean grid must "tile" each other. The masked region of the ice grid
!    and ocean grid must be identical. 
!
!         ATMOSPHERE  |----|----|----|----|----|----|----|----|
!
!               LAND  |---|---|---|---|xxx|xxx|xxx|xxx|xxx|xxx|
!
!                ICE  |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
!
!               OCEAN |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
!
!              where  |xxx| = masked grid point
!         
!
!    The atmosphere, land, and ice grids exchange information using the exchange grid xmap_sfc.
!
!    The land and ice grids exchange runoff data using the exchange grid xmap_runoff.
!
!    Transfer of data between the ice bottom and ocean does not require an exchange 
!    grid as the grids are physically identical. The flux routines will automatically
!    detect and redistribute data if their domain decompositions are different.
!
!    To get information from the atmosphere to the ocean it must pass through the 
!    ice model, first by interpolating from the atmospheric grid to the ice grid, 
!    and then transferring from the ice grid to the ocean grid.

!  2.Each component model must have a public defined data type containing specific 
!    boundary fields. A list of these quantities is located in the NOTES of this document. 
!
!  3.The surface flux of sensible heat and surface evaporation can be implicit functions
!    of surface temperature. As a consequence, the parts of the land and sea-ice models 
!    that update the surface temperature must be called on the atmospheric time step 
!
!  4.The surface fluxes of all other tracers and of momentum are assumed to be explicit
!    functions of all surface parameters 
!
!  5.While no explicit reference is made within this module to the implicit treatment 
!    of vertical diffusion in the atmosphere and in the land or sea-ice models, the 
!    module is designed to allow for simultaneous implicit time integration on both 
!    sides of the surface interface. 
!
!  6.Due to #5, the diffusion part of the land and ice models must be called on the 
!    atmospheric time step.
  
!7. Any field passed from one component to another may be "faked" to a
!   constant value, or to data acquired from a file, using the
!   data_override feature of FMS. The fields to override are runtime
!   configurable, using the text file <tt>data_table</tt> for input.
!   See the data_override_mod documentation for more details.
!
!   We DO NOT RECOMMEND exercising the data override capabilities of
!   the FMS coupler until the user has acquired considerable
!   sophistication in running FMS.
!
!   Here is a listing of the override capabilities of the flux_exchange
!   module:
!
!   FROM the atmosphere boundary TO the exchange grid (in sfc_boundary_layer):
!  
!        t_bot, q_bot, z_bot, p_bot, u_bot, v_bot, p_surf, slp, gust
!
!   FROM the ice boundary TO the exchange grid (in sfc_boundary_layer):
!
!        t_surf, rough_mom, rough_heat, rough_moist, albedo, u_surf, v_surf
!     
!   FROM the land boundary TO the exchange grid (in sfc_boundary_layer):
!
!        t_surf, t_ca, q_ca, rough_mom, rough_heat, albedo
!
!   FROM the exchange grid TO land_ice_atmos_boundary (in
!   sfc_boundary_layer):
!
!        t, albedo, land_frac, dt_t, dt_q, u_flux, v_flux, dtaudu, dtaudv,
!        u_star, b_star, rough_mom
!   
!   FROM the atmosphere boundary TO the exchange grid (in
!    flux_down_from_atmos):
!
!        flux_sw, flux_lw, lprec, fprec, coszen, dtmass, delta_t,
!        delta_q, dflux_t, dflux_q
!        
!   FROM the exchange grid TO the land boundary (in
!    flux_down_from_atmos):
!
!    t_flux, q_flux, lw_flux, sw_flux, lprec, fprec, dhdt, dedt, dedq,
!    drdt, drag_q, p_surf
!    
!   FROM the exchange grid TO the ice boundary (in flux_down_from_atmos):
!
!        u_flux, v_flux, t_flux, q_flux, lw_flux, lw_flux_dn, sw_flux,
!        sw_flux_dn, lprec, fprec, dhdt, dedt, drdt, coszen, p 
!
!   FROM the land boundary TO the ice boundary (in flux_land_to_ice):
!
!        runoff, calving
!
!   FROM the ice boundary TO the ocean boundary (in flux_ice_to_ocean):
! 
!        u_flux, v_flux, t_flux, q_flux, salt_flux, lw_flux, sw_flux,
!        lprec, fprec, runoff, calving, p
!        
!   FROM the ocean boundary TO the ice boundary (in flux_ocean_to_ice):
!
!        u, v, t, s, frazil, sea_level
!
!   FROM the ice boundary TO the atmosphere boundary (in flux_up_to_atmos):
!
!        t_surf
!
!   FROM the land boundary TO the atmosphere boundary (in
!    flux_up_to_atmos):
!  
!        t_ca, t_surf, q_ca
!
!  See NOTES below for an explanation of the field names.
!  </PRE>
! </DESCRIPTION>

  use mpp_mod,         only: mpp_npes, mpp_pe, mpp_root_pe, &
       mpp_error, stderr, stdout, stdlog, FATAL, NOTE, mpp_set_current_pelist, &
       mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_sum, &
       CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_ROUTINE, lowercase, &
       input_nml_file
                    
  use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_compute_domains, &
                             mpp_global_sum, mpp_redistribute, operator(.EQ.)
  use mpp_domains_mod, only: mpp_get_global_domain, mpp_get_data_domain
  use mpp_domains_mod, only: mpp_set_global_domain, mpp_set_data_domain, mpp_set_compute_domain
  use mpp_domains_mod, only: mpp_deallocate_domain, mpp_copy_domain, domain2d

  use mpp_io_mod,      only: mpp_close, mpp_open, MPP_MULTI, MPP_SINGLE, MPP_OVERWR

!model_boundary_data_type contains all model fields at the boundary.
!model1_model2_boundary_type contains fields that model2 gets
!from model1, may also include fluxes. These are declared by
!flux_exchange_mod and have private components. All model fields in
!model_boundary_data_type may not be exchanged.
!will support 3 types of flux_exchange:
!REGRID: physically distinct grids, via xgrid
!REDIST: same grid, transfer in index space only
!DIRECT: same grid, same decomp, direct copy
  use atmos_model_mod, only: atmos_data_type, land_ice_atmos_boundary_type
  use ocean_model_mod, only: ocean_public_type, ice_ocean_boundary_type
  use ocean_model_mod, only: ocean_state_type
  use ice_model_mod,   only: ice_data_type, land_ice_boundary_type, &
       ocean_ice_boundary_type, atmos_ice_boundary_type, Ice_stock_pe, &
       ice_cell_area => cell_area
  use    land_model_mod, only:  land_data_type, atmos_land_boundary_type

  use  surface_flux_mod, only: surface_flux
  use monin_obukhov_mod, only: mo_profile     

  use xgrid_mod, only: xmap_type, setup_xmap, set_frac_area, &
       put_to_xgrid, get_from_xgrid, &
       xgrid_count, some, conservation_check, xgrid_init, &
       get_ocean_model_area_elements, stock_integrate_2d, &
       stock_move, stock_print


  use diag_integral_mod, only:     diag_integral_field_init, &
       sum_diag_integral_field

  use  diag_manager_mod, only: register_diag_field,  &
       register_static_field, send_data, send_tile_averaged_data

  use  time_manager_mod, only: time_type

  use sat_vapor_pres_mod, only: compute_qs

  use      constants_mod, only: rdgas, rvgas, cp_air, stefan, WTMAIR, HLV, HLF, Radius, PI, &
                                WTMCO2, WTMC

  use ocean_parameters_mod, only: cp_ocean

!Balaji
!utilities stuff into use fms_mod
  use fms_mod,                    only: clock_flag_default, check_nml_error, error_mesg
  use fms_mod,                    only: open_namelist_file, write_version_number
  use fms_mod,                    only: field_exist, field_size, read_data, get_mosaic_tile_grid

  use data_override_mod,          only: data_override
  use coupler_types_mod,          only: coupler_1d_bc_type
  use atmos_ocean_fluxes_mod,     only: atmos_ocean_fluxes_init, atmos_ocean_fluxes_calc
  use ocean_model_mod,            only: ocean_model_init_sfc, ocean_model_flux_init, ocean_model_data_get
  use coupler_types_mod,          only: coupler_type_copy
  use coupler_types_mod,          only: ind_psurf, ind_u10
  use atmos_tracer_driver_mod,    only: atmos_tracer_flux_init

  use field_manager_mod,          only: MODEL_ATMOS, MODEL_LAND, MODEL_ICE
  use tracer_manager_mod,         only: get_tracer_index
  use tracer_manager_mod,         only: get_tracer_names, get_number_tracers, NO_TRACER

  use stock_constants_mod,        only: NELEMS, ISTOCK_WATER, ISTOCK_HEAT, ISTOCK_SALT
  use stock_constants_mod,        only: ISTOCK_SIDE, ISTOCK_TOP, ISTOCK_BOTTOM , STOCK_UNITS, STOCK_NAMES
  use stock_constants_mod,        only: stocks_file, stocks_report, stocks_report_init
  use stock_constants_mod,        only: Atm_stock, Ocn_stock, Lnd_stock, Ice_stock
  use land_model_mod,             only: Lnd_stock_pe
  use ocean_model_mod,            only: Ocean_stock_pe
  use atmos_model_mod,            only: Atm_stock_pe

#ifdef SCM
! option to override various surface boundary conditions for SCM
  use scm_forc_mod,               only: do_specified_flux, scm_surface_flux,             &
                                        do_specified_tskin, TSKIN,                       &
                                        do_specified_albedo, ALBEDO_OBS,                 &
                                        do_specified_rough_leng, ROUGH_MOM, ROUGH_HEAT,  &
                                        do_specified_land
#endif

  implicit none
  include 'netcdf.inc'
private

  character(len=48), parameter :: module_name = 'flux_exchange_mod'

  public :: flux_exchange_init,   &
     sfc_boundary_layer,   &
     generate_sfc_xgrid,   &
     flux_down_from_atmos, &
     flux_up_to_atmos,     &
     flux_land_to_ice,     &
     flux_ice_to_ocean,    &
     flux_ocean_to_ice,    &
     flux_check_stocks,    &
     flux_init_stocks,     &
     flux_ice_to_ocean_stocks,&
     flux_ocean_from_ice_stocks

!-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: flux_exchange.F90,v 20.0 2013/12/13 23:27:41 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'
!-----------------------------------------------------------------------
!---- exchange grid maps -----

type(xmap_type), save :: xmap_sfc, xmap_runoff

integer         :: n_xgrid_sfc=0,  n_xgrid_runoff=0

!-----------------------------------------------------------------------
!-------- namelist (for diagnostics) ------

character(len=4), parameter :: mod_name = 'flux'

  integer :: id_drag_moist,  id_drag_heat,  id_drag_mom,     &
     id_rough_moist, id_rough_heat, id_rough_mom,    &
     id_land_mask,   id_ice_mask,     &
     id_u_star, id_b_star, id_q_star, id_u_flux, id_v_flux,   &
     id_t_surf, id_t_flux, id_r_flux, id_q_flux, id_slp,      &
     id_t_atm,  id_u_atm,  id_v_atm,  id_wind,                &
     id_t_ref,  id_rh_ref, id_u_ref,  id_v_ref, id_wind_ref,  &
     id_del_h,  id_del_m,  id_del_q,  id_rough_scale,         &
     id_t_ca,   id_q_surf, id_q_atm, id_z_atm, id_p_atm, id_gust, &
     id_t_ref_land, id_rh_ref_land, id_u_ref_land, id_v_ref_land, &
     id_q_ref,  id_q_ref_land, id_q_flux_land, id_rh_ref_cmip

integer :: id_co2_atm_dvmr, id_co2_surf_dvmr

integer, allocatable :: id_tr_atm(:), id_tr_surf(:), id_tr_flux(:), id_tr_mol_flux(:)

logical :: first_static = .true.
logical :: do_init = .true.
integer :: remap_method = 1

real, parameter :: bound_tol = 1e-7

real, parameter :: d622 = rdgas/rvgas
real, parameter :: d378 = 1.0-d622

!--- namelist interface ------------------------------------------------------
! <NAMELIST NAME="flux_exchange_nml">
!   <DATA NAME="z_ref_heat"  TYPE="real"  DEFAULT="2.0">
!    eference height (meters) for temperature and relative humidity 
!    diagnostics (t_ref,rh_ref,del_h,del_q)
!   </DATA>
!   <DATA NAME="z_ref_mom"  TYPE="real"  DEFAULT="10.0">
!    reference height (meters) for momentum diagnostics (u_ref,v_ref,del_m)
!   </DATA>
!   <DATA NAME="ex_u_star_smooth_bug"  TYPE="logical"  DEFAULT="false">
!    By default, the global exchange grid u_star will not be interpolated from 
!    atmospheric grid, this is different from Jakarta behavior and will
!    change answers. So to perserve Jakarta behavior and reproduce answers
!    explicitly set this namelist variable to .true. in input.nml.
!    Talk to mw, ens for details.
!   </DATA>
!   <DATA NAME="do_runoff"  TYPE="logical"  DEFAULT=".TRUE.">
!    Turns on/off the land runoff interpolation to the ocean.
!   </DATA>


  real ::  z_ref_heat =  2.,  &
           z_ref_mom  = 10.
  logical :: ex_u_star_smooth_bug = .false.
  logical :: sw1way_bug = .false.
  logical :: do_area_weighted_flux = .FALSE.
  logical :: debug_stocks = .FALSE.
  logical :: divert_stocks_report = .FALSE.
  logical :: do_runoff = .TRUE.
  logical :: do_forecast = .false.

namelist /flux_exchange_nml/ z_ref_heat, z_ref_mom, ex_u_star_smooth_bug, sw1way_bug, &
         do_area_weighted_flux, debug_stocks, divert_stocks_report, do_runoff, do_forecast
! </NAMELIST>

! ---- allocatable module storage --------------------------------------------
real, allocatable, dimension(:) :: &
     ! NOTE: T canopy is only differet from t_surf over vegetated land
     ex_t_surf,    &   ! surface temperature for radiation calc, degK
     ex_t_surf_miz,&   ! miz
     ex_t_ca,      &   ! near-surface (canopy) air temperature, degK
     ex_p_surf,    &   ! surface pressure
     ex_slp,       &   ! surface pressure

     ex_flux_t,    &   ! sens heat flux
     ex_flux_lw,   &   ! longwave radiation flux

     ex_dhdt_surf, &   ! d(sens.heat.flux)/d(T canopy)
     ex_dedt_surf, &   ! d(water.vap.flux)/d(T canopy)
     ex_dqsatdt_surf, &   ! d(water.vap.flux)/d(q canopy)
     ex_e_q_n,     &
     ex_drdt_surf, &   ! d(LW flux)/d(T surf)
     ex_dhdt_atm,  &   ! d(sens.heat.flux)/d(T atm)
     ex_flux_u,    &   ! u stress on atmosphere
     ex_flux_v,    &   ! v stress on atmosphere
     ex_dtaudu_atm,&   ! d(stress)/d(u)
     ex_dtaudv_atm,&   ! d(stress)/d(v)
     ex_albedo_fix,&
     ex_albedo_vis_dir_fix,&
     ex_albedo_nir_dir_fix,&
     ex_albedo_vis_dif_fix,&
     ex_albedo_nir_dif_fix,&
     ex_old_albedo,&   ! old value of albedo for downward flux calculations
     ex_drag_q,    &   ! q drag.coeff.
     ex_cd_t,      &
     ex_cd_m,      &
     ex_b_star,    &
     ex_u_star,    &
     ex_wind,      &
     ex_z_atm

#ifdef SCM
real, allocatable, dimension(:) :: &
     ex_dhdt_surf_forland, &
     ex_dedt_surf_forland, &
     ex_dedq_surf_forland
#endif

real, allocatable, dimension(:,:) :: &
     ex_tr_surf,    & ! near-surface tracer fields
     ex_flux_tr,    & ! tracer fluxes
     ex_dfdtr_surf, & ! d(tracer flux)/d(surf tracer)
     ex_dfdtr_atm,  & ! d(tracer flux)/d(atm tracer)
     ex_e_tr_n,     & ! coefficient in implicit scheme 
     ex_f_tr_delt_n   ! coefficient in implicit scheme

logical, allocatable, dimension(:) :: &
     ex_avail,     &   ! true where data on exchange grid are available
     ex_land           ! true if exchange grid cell is over land
real, allocatable, dimension(:) :: &
     ex_e_t_n,      &
     ex_f_t_delt_n

integer :: n_atm_tr  ! number of prognostic tracers in the atmos model
integer :: n_atm_tr_tot  ! number of prognostic tracers in the atmos model
integer :: n_lnd_tr  ! number of prognostic tracers in the land model 
integer :: n_lnd_tr_tot  ! number of prognostic tracers in the land model 
integer :: n_exch_tr ! number of tracers exchanged between models

type :: tracer_ind_type
   integer :: atm, ice, lnd ! indices of the tracer in the respective models
end type 
type(tracer_ind_type), allocatable :: tr_table(:) ! table of tracer indices
type :: tracer_exch_ind_type
   integer :: exch = 0  ! exchange grid index
   integer :: ice = 0   ! ice model index
   integer :: lnd = 0   ! land model index
end type tracer_exch_ind_type
type(tracer_exch_ind_type), allocatable :: tr_table_map(:) ! map atm tracers to exchange, ice and land variables
integer :: isphum = NO_TRACER       ! index of specific humidity tracer in tracer table
integer :: ico2   = NO_TRACER       ! index of co2 tracer in tracer table

type(coupler_1d_bc_type), save        :: ex_gas_fields_atm  ! gas fields in atm
                     ! Place holder for various atmospheric fields.
type(coupler_1d_bc_type), save        :: ex_gas_fields_ice  ! gas fields on ice
type(coupler_1d_bc_type), save        :: ex_gas_fluxes      ! gas flux
                     ! Place holder of intermediate calculations, such as
                     ! piston velocities etc.

integer :: ni_atm, nj_atm ! to do atmos diagnostic from flux_ocean_to_ice
real, dimension(3) :: ccc ! for conservation checks
!Balaji, sets boundary_type%xtype
!  REGRID: grids are physically different, pass via exchange grid
!  REDIST: same physical grid, different decomposition, must move data around
!  DIRECT: same physical grid, same domain decomposition, can directly copy data
integer, parameter :: REGRID=1, REDIST=2, DIRECT=3
!Balaji: clocks moved into flux_exchange
!RASF Initialise to zero irrespective  of model component
  integer :: cplClock=0, sfcClock=0, fluxAtmDnClock=0, fluxLandIceClock=0, &
             fluxIceOceanClock=0, fluxOceanIceClock=0, regenClock=0, fluxAtmUpClock=0, &
             cplOcnClock=0

  logical :: ocn_pe, ice_pe
  integer, allocatable, dimension(:) :: ocn_pelist, ice_pelist

  ! Exchange grid indices
  integer :: X1_GRID_ATM, X1_GRID_ICE, X1_GRID_LND
  integer :: X2_GRID_LND, X2_GRID_ICE
  real    :: Dt_atm, Dt_cpl
  real    :: ATM_PRECIP_NEW

integer ::  runoff_id_diag =-1 

integer :: nxc_ocn=0, nyc_ocn=0, nxc_ice=0, nyc_ice=0, nk_ice=0


contains

!#######################################################################
! <SUBROUTINE NAME="flux_exchange_init">
!  <OVERVIEW>
!   Initialization routine.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Initializes the interpolation routines,diagnostics and boundary data
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_exchange_init ( Time, Atm, Land, Ice, Ocean, &
!		atmos_ice_boundary, land_ice_atmos_boundary, &
!		land_ice_boundary, ice_ocean_boundary, ocean_ice_boundary, &
!                dt_atmos, dt_cpld )
!		
!  </TEMPLATE>
!  <IN NAME=" Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Atm" TYPE="atmos_data_type">
!   A derived data type to specify atmosphere boundary data.
!  </IN>
!  <IN NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </IN>
!  <IN NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>
!  <IN NAME="Ocean" TYPE="ocean_public_type">
!   A derived data type to specify ocean boundary data.
!  </IN>
!  <INOUT NAME="atmos_ice_boundary" TYPE="atmos_ice_boundary_type">
!   A derived data type to specify properties and fluxes passed from atmosphere to ice.
!  </INOUT>
!  <INOUT NAME="land_ice_atmos_boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere, land and ice.
!  </INOUT>
!  <INOUT NAME="land_ice_boundary" TYPE="land_ice_boundary_type">
!   A derived data type to specify properties and fluxes passed from land to ice.
!  </INOUT>
!  <INOUT NAME="ice_ocean_boundary" TYPE="ice_ocean_boundary_type">
!  A derived data type to specify properties and fluxes passed from ice to ocean.
!  </INOUT>
!  <INOUT NAME="ocean_ice_boundary" TYPE="ocean_ice_boundary_type">
!  A derived data type to specify properties and fluxes passed from ocean to ice.
!  </INOUT>
!  <IN NAME="dt_atmos" TYPE="integer">
!  Atmos time step in secs.
!  </IN>
!  <IN NAME="dt_cpld" TYPE="integer">
!  Coupled time step in secs.
!  </IN>

!
subroutine flux_exchange_init ( Time, Atm, Land, Ice, Ocean, Ocean_state,&
       atmos_ice_boundary, land_ice_atmos_boundary, &
       land_ice_boundary, ice_ocean_boundary, ocean_ice_boundary, &
       dt_atmos, dt_cpld )

  type(time_type),                   intent(in)  :: Time
  type(atmos_data_type),             intent(inout)  :: Atm
  type(land_data_type),              intent(in)  :: Land
  type(ice_data_type),               intent(inout)  :: Ice
  type(ocean_public_type),           intent(inout)  :: Ocean
  type(ocean_state_type),            pointer        :: Ocean_state
! All intent(OUT) derived types with pointer components must be 
! COMPLETELY allocated here and in subroutines called from here;
! NO pointer components should have been allocated before entry if the
! derived type has intent(OUT) otherwise they may be lost.
  type(atmos_ice_boundary_type),     intent(inout) :: atmos_ice_boundary
  type(land_ice_atmos_boundary_type),intent(inout) :: land_ice_atmos_boundary
  type(land_ice_boundary_type),      intent(inout) :: land_ice_boundary
  type(ice_ocean_boundary_type),     intent(inout) :: ice_ocean_boundary
  type(ocean_ice_boundary_type),     intent(inout) :: ocean_ice_boundary
  integer, optional,                 intent(in)    :: dt_atmos, dt_cpld

  character(len=64), parameter    :: sub_name = 'flux_exchange_init'
  character(len=256), parameter   :: error_header = '==>Error from ' // trim(module_name) //   &
                                                    '(' // trim(sub_name) // '):'
  character(len=256), parameter   :: warn_header = '==>Warning from ' // trim(module_name) //  &
                                                   '(' // trim(sub_name) // '):'
  character(len=256), parameter   :: note_header = '==>Note from ' // trim(module_name) //     &
                                                   '(' // trim(sub_name) // '):'
  character(len=64),  parameter   :: grid_file = 'INPUT/grid_spec.nc'  
  character(len=256)              :: atm_mosaic_file, tile_file 

  type(domain2d) :: domain2
  integer        :: isg, ieg, jsg, jeg
  integer        :: isc, iec, jsc, jec
  integer        :: isd, ied, jsd, jed
  integer        :: isc2, iec2, jsc2, jec2
  integer        :: nxg, nyg, ioff, joff
  integer        :: unit, ierr, io,  i, j
  integer        :: nlon, nlat, siz(4)
  integer        :: outunit, logunit
  real, dimension(:,:), allocatable :: tmpx(:,:), tmpy(:,:)
  real, dimension(:),   allocatable :: atmlonb, atmlatb
  integer :: is, ie, js, je, kd
  character(32) :: tr_name
  logical       :: found

  integer              :: n, npes_atm, npes_ocn, npes_all
  integer, allocatable :: pelist(:)

!-----------------------------------------------------------------------

!
!       initialize atmos_ocean_fluxes
! Setting up flux types, allocates the arrays.
!

!
!       ocean_tracer_flux_init is called first since it has the meaningful value to set
!       for the input/output file names for the tracer flux values used in restarts. These
!       values could be set in the field table, and this ordering allows this.
!       atmos_tracer_flux_init is called last since it will use the values set in 
!       ocean_tracer_flux_init with the exception of atm_tr_index, which can only
!       be meaningfully set from the atmospheric model (not from the field table)
!

    call ocean_model_flux_init(Ocean_state)
    call atmos_tracer_flux_init
    call atmos_ocean_fluxes_init(ex_gas_fluxes, ex_gas_fields_atm, ex_gas_fields_ice)

!-----------------------------------------------------------------------
    outunit = stdout(); logunit = stdlog()
!----- read namelist -------

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, flux_exchange_nml, iostat=io)
      ierr = check_nml_error (io, 'flux_exchange_nml')
#else
    unit = open_namelist_file()
    ierr=1; do while (ierr /= 0)
       read  (unit, nml=flux_exchange_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'flux_exchange_nml')
    enddo
10  call mpp_close(unit)
#endif

!----- write namelist to logfile -----
    call write_version_number(version, tagname)
    if( mpp_pe() == mpp_root_pe() )write( logunit, nml=flux_exchange_nml )

!----- find out number of atmospheric prognostic tracers and index of specific 
!      humidity in the tracer table
  call get_number_tracers (MODEL_ATMOS, num_tracers=n_atm_tr_tot, &
                           num_prog=n_atm_tr)
  call get_number_tracers (MODEL_LAND, num_tracers=n_lnd_tr_tot, &
                           num_prog=n_lnd_tr)

  ! assemble the table of tracer number translation by matching names of
  ! prognostic tracers in the atmosphere and surface models; skip all atmos.
  ! tracers that have no corresponding surface tracers.
  allocate(tr_table(n_atm_tr))
  allocate(tr_table_map(n_atm_tr))
  n = 1
  do i = 1,n_atm_tr
     call get_tracer_names( MODEL_ATMOS, i, tr_name )
     tr_table(n)%atm = i
     tr_table(n)%ice = get_tracer_index ( MODEL_ICE,  tr_name )
     tr_table_map(i)%ice = tr_table(n)%ice
     tr_table(n)%lnd = get_tracer_index ( MODEL_LAND, tr_name )
     tr_table_map(i)%lnd = tr_table(n)%lnd
     if(tr_table(n)%ice/=NO_TRACER.or.tr_table(n)%lnd/=NO_TRACER) then
       tr_table_map(i)%exch = n
       n = n + 1
     endif
  enddo
  n_exch_tr = n - 1
  !
  !     Set up tracer table entries for ocean-atm gas fluxes where the names of tracers in the
  !     atmosphere and ocean may not be equal
  !
  do n = 1, ex_gas_fluxes%num_bcs  !{
    if (ex_gas_fluxes%bc(n)%atm_tr_index .gt. 0) then  !{
      found = .false.
      do i = 1, n_exch_tr  !{
        if (ex_gas_fluxes%bc(n)%atm_tr_index .eq. tr_table(i)%atm) then
          found = .true.
          exit
        endif
      enddo  !} i
      if (.not. found) then
        n_exch_tr = n_exch_tr + 1
        tr_table(n_exch_tr)%atm = ex_gas_fluxes%bc(n)%atm_tr_index
        tr_table(n_exch_tr)%ice = NO_TRACER ! because ocean-atm gas fluxes are not held in the ice model as tracers
        tr_table(n_exch_tr)%lnd = NO_TRACER ! because this would have been found above
        tr_table_map(n_exch_tr)%exch = n_exch_tr
        tr_table_map(n_exch_tr)%ice = tr_table(n_exch_tr)%ice
        tr_table_map(n_exch_tr)%lnd = tr_table(n_exch_tr)%lnd
      endif
    endif  !}
  enddo  !} n
  write(outunit,*) trim(note_header), ' Number of exchanged tracers = ', n_exch_tr
  write(logunit,*) trim(note_header), ' Number of exchanged tracers = ', n_exch_tr
  do i = 1,n_exch_tr
     call get_tracer_names( MODEL_ATMOS, tr_table(i)%atm, tr_name )
     write(outunit,*)'Tracer field name :'//trim(tr_name)
     write(logunit,*)'Tracer field name :'//trim(tr_name)
  enddo

  ! find out which tracer is specific humidity

  ! +fix-me-slm+ specific humidity may not be present if we are running with
  ! dry atmosphere. Besides, model may use mixing ratio ('mix_rat') (?). However,
  ! some atmos code also assumes 'sphum' is present, so for now the following
  ! code may be good enough.

  do i = 1,n_exch_tr
     call get_tracer_names( MODEL_ATMOS, tr_table(i)%atm, tr_name )
     if(lowercase(tr_name)=='sphum') then
        isphum = i
     endif
  ! jgj: find out which exchange tracer is co2
     if(lowercase(tr_name)=='co2') then
        ico2 = i
        write(outunit,*)'Exchange tracer index for '//trim(tr_name),' : ',ico2
     endif
  enddo

  if (isphum==NO_TRACER) then
     call error_mesg('flux_exchange_mod',&
          'tracer "sphum" must be present in the atmosphere', FATAL )
  endif

  if (ico2==NO_TRACER) then
     call error_mesg('flux_exchange_mod',&
          'tracer "co2" not present in the atmosphere', NOTE )
  endif

!--------- read gridspec file ------------------
!only atmos pelists needs to do it here, ocean model will do it elsewhere

    ice_pe = Atm%pe
    ocn_pe = Ocean%is_ocean_pe
    allocate( ice_pelist(size(Atm%pelist)) ) !if ice/land become concurrent, this won't be true...
    ice_pelist(:) = Atm%pelist(:)
    allocate( ocn_pelist(size(Ocean%pelist)) )
    ocn_pelist(:) = Ocean%pelist(:)

    call get_ocean_model_area_elements(Ocean%domain, grid_file)

    if( Atm%pe )then
       call mpp_set_current_pelist(Atm%pelist)

       !
       ! check atmosphere and grid_spec.nc have same atmosphere lat/lon boundaries
       !
       call mpp_get_global_domain(Atm%domain, isg, ieg, jsg, jeg, xsize=nxg, ysize=nyg)
       call mpp_get_compute_domain(Atm%domain, isc, iec, jsc, jec)
       call mpp_get_data_domain(Atm%domain, isd, ied, jsd, jed)
       if(size(Atm%lon_bnd,1) .NE. iec-isc+2 .OR. size(Atm%lon_bnd,2) .NE. jec-jsc+2) then
          call error_mesg ('flux_exchange_mod',  &
              'size of Atm%lon_bnd does not match the Atm computational domain', FATAL)          
       endif
       ioff = lbound(Atm%lon_bnd,1) - isc
       joff = lbound(Atm%lon_bnd,2) - jsc
       if(field_exist(grid_file, "AREA_ATM" ) ) then  ! old grid
          call field_size(grid_file, "AREA_ATM", siz)
          nlon = siz(1)
          nlat = siz(2)          
          
          if (nlon /= nxg .or. nlat /= nyg) then
             if (mpp_pe()==mpp_root_pe()) then
                print *, 'grid_spec.nc has', nlon, 'longitudes,', nlat, 'latitudes; ', &
                     'atmosphere has', nxg, 'longitudes,', &
                     nyg, 'latitudes (see xba.dat and yba.dat)'
             end if
             !   <ERROR MSG="grid_spec.nc incompatible with atmosphere resolution" STATUS="FATAL">
             !      The atmosphere grid size from file grid_spec.nc is not compatible with the atmosphere 
             !      resolution from atmosphere model.
             !   </ERROR>
             call error_mesg ('flux_exchange_mod',  &
                  'grid_spec.nc incompatible with atmosphere resolution', FATAL)
          end if
          allocate( atmlonb(isg:ieg+1) )
          allocate( atmlatb(jsg:jeg+1) )
          call read_data(grid_file, 'xba', atmlonb, no_domain=.true. )
          call read_data(grid_file, 'yba', atmlatb, no_domain=.true. )

          do i=isc, iec+1
             if(abs(atmlonb(i)-Atm%lon_bnd(i+ioff,jsc+joff)*45/atan(1.0))>bound_tol) then
                print *, 'GRID_SPEC/ATMOS LONGITUDE INCONSISTENCY at i= ',i, ': ', &
                     atmlonb(i),  Atm%lon_bnd(i+ioff,jsc+joff)*45/atan(1.0)
                call error_mesg ('flux_exchange_mod', &
                     'grid_spec.nc incompatible with atmosphere longitudes (see xba.dat and yba.dat)'&
                     , FATAL)
             endif
          enddo
          !   <ERROR MSG="grid_spec.nc incompatible with atmosphere longitudes (see xba.dat and yba.dat)" STATUS="FATAL">
          !      longitude from file grid_spec.nc ( from field yba ) is different from the longitude from atmosphere model.
          !   </ERROR>
          do j=jsc, jec+1
             if(abs(atmlatb(j)-Atm%lat_bnd(isc+ioff,j+joff)*45/atan(1.0))>bound_tol) then
                print *, 'GRID_SPEC/ATMOS LATITUDE INCONSISTENCY at j= ',j, ': ', &
                     atmlatb(j),  Atm%lat_bnd(isc+ioff, j+joff)*45/atan(1.0)
                call error_mesg ('flux_exchange_mod', &
                     'grid_spec.nc incompatible with atmosphere latitudes (see xba.dat and yba.dat)'&
                     , FATAL)
             endif
          enddo
          deallocate(atmlonb, atmlatb)
        else if(field_exist(grid_file, "atm_mosaic_file" ) ) then  ! mosaic grid file.
           call read_data(grid_file, 'atm_mosaic_file', atm_mosaic_file)
           call get_mosaic_tile_grid(tile_file, 'INPUT/'//trim(atm_mosaic_file), Atm%domain)          
           call field_size(tile_file, 'area', siz)
           nlon = siz(1); nlat = siz(2)
           if( mod(nlon,2) .NE. 0) call mpp_error(FATAL,  &
                'flux_exchange_mod: atmos supergrid longitude size can not be divided by 2')
           if( mod(nlat,2) .NE. 0) call mpp_error(FATAL,  &
                'flux_exchange_mod: atmos supergrid latitude size can not be divided by 2')
           nlon = nlon/2
           nlat = nlat/2
           if (nlon /= nxg .or. nlat /= nyg) then
             if (mpp_pe()==mpp_root_pe()) then
                print *, 'atmosphere mosaic tile has', nlon, 'longitudes,', nlat, 'latitudes; ', &
                     'atmosphere has', nxg, 'longitudes,', nyg, 'latitudes'
             end if
            call error_mesg ('flux_exchange_mod',  &
                  'atmosphere mosaic tile grid file incompatible with atmosphere resolution', FATAL)
           end if

           call mpp_copy_domain(Atm%domain, domain2)
           call mpp_set_compute_domain(domain2, 2*isc-1, 2*iec+1, 2*jsc-1, 2*jec+1, 2*(iec-isc)+3, 2*(jec-jsc)+3 )
           call mpp_set_data_domain   (domain2, 2*isd-1, 2*ied+1, 2*jsd-1, 2*jed+1, 2*(ied-isd)+3, 2*(jed-jsd)+3 )   
           call mpp_set_global_domain (domain2, 2*isg-1, 2*ieg+1, 2*jsg-1, 2*jeg+1, 2*(ieg-isg)+3, 2*(jeg-jsg)+3 )   
           call mpp_get_compute_domain(domain2, isc2, iec2, jsc2, jec2)
           if(isc2 .NE. 2*isc-1 .OR. iec2 .NE. 2*iec+1 .OR. jsc2 .NE. 2*jsc-1 .OR. jec2 .NE. 2*jec+1) then
              call mpp_error(FATAL, 'flux_exchange_mod: supergrid domain is not set properly')
           endif

           allocate(tmpx(isc2:iec2,jsc2:jec2), tmpy(isc2:iec2,jsc2:jec2) )

           call read_data( tile_file, 'x', tmpx, domain2)
           call read_data( tile_file, 'y', tmpy, domain2)     
           call mpp_deallocate_domain(domain2)

           do j = jsc, jec+1
              do i = isc, iec+1
                 if (abs(tmpx(2*i-1,2*j-1)-Atm%lon_bnd(i+ioff,j+joff)*45/atan(1.0))>bound_tol) then
                    print *, 'GRID_SPEC/ATMOS LONGITUDE INCONSISTENCY at i= ',i, ', j= ', j, ': ', &
                         tmpx(2*i-1,2*j-1),  Atm%lon_bnd(i+ioff,j+joff)*45/atan(1.0)
                    !   <ERROR MSG="grid_spec.nc incompatible with atmosphere longitudes (see xba.dat and yba.dat)" STATUS="FATAL">
                    !      longitude from file grid_spec.nc ( from field xba ) is different from the longitude from atmosphere model.
                    !   </ERROR>
                    call error_mesg ('flux_exchange_mod', &
                         'grid_spec.nc incompatible with atmosphere longitudes (see '//trim(tile_file)//')'&
                         ,FATAL)
                 end if
                 if (abs(tmpy(2*i-1,2*j-1)-Atm%lat_bnd(i+ioff,j+joff)*45/atan(1.0))>bound_tol) then
                    print *, 'GRID_SPEC/ATMOS LATITUDE INCONSISTENCY at i= ',i, ', j= ', j, ': ', &
                         tmpy(2*i-1,2*j-1),  Atm%lat_bnd(i+ioff,j+joff)*45/atan(1.0)
                    !   <ERROR MSG="grid_spec.nc incompatible with atmosphere latitudes (see grid_spec.nc)" STATUS="FATAL">
                    !      latgitude from file grid_spec.nc is different from the latitude from atmosphere model.
                    !   </ERROR>
                    call error_mesg ('flux_exchange_mod', &
                         'grid_spec.nc incompatible with atmosphere latitudes (see '//trim(tile_file)//')'&
                         ,FATAL)
                 end if
              end do
           end do
           deallocate(tmpx, tmpy)
        else
           call mpp_error(FATAL, 'flux_exchange_mod: both AREA_ATMxOCN and ocn_mosaic_file does not exist in '//trim(grid_file))
        end if

 
        call xgrid_init(remap_method)

        call setup_xmap(xmap_sfc, (/ 'ATM', 'OCN', 'LND' /),   &
             (/ Atm%Domain, Ice%Domain, Land%Domain /),        &
             "INPUT/grid_spec.nc", Atm%grid)
        ! exchange grid indices
        X1_GRID_ATM = 1; X1_GRID_ICE = 2; X1_GRID_LND = 3;
        call generate_sfc_xgrid( Land, Ice )
        if (n_xgrid_sfc.eq.1) write (*,'(a,i4,6x,a)') 'PE = ', mpp_pe(), 'Surface exchange size equals one.'

        if (do_runoff) then
           call setup_xmap(xmap_runoff, (/ 'LND', 'OCN' /),       &
                (/ Land%Domain, Ice%Domain /),                    &
                "INPUT/grid_spec.nc"             )
           ! exchange grid indices
           X2_GRID_LND = 1; X2_GRID_ICE = 2;
           n_xgrid_runoff = max(xgrid_count(xmap_runoff),1)
           if (n_xgrid_runoff.eq.1) write (*,'(a,i4,6x,a)') 'PE = ', mpp_pe(), 'Runoff  exchange size equals one.'
        endif

!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!----- initialize quantities for global integral package -----

!! call diag_integral_field_init ('prec', 'f6.3')
        call diag_integral_field_init ('evap', 'f6.3')

!-----------------------------------------------------------------------
!----- initialize diagnostic fields -----
!----- all fields will be output on the atmospheric grid -----

        call diag_field_init ( Time, Atm%axes(1:2), Land%axes )
        ni_atm = size(Atm%lon_bnd,1)-1 ! to dimension "diag_atm"
        nj_atm = size(Atm%lon_bnd,2)-1 ! in flux_ocean_to_ice

!Balaji
        
!allocate atmos_ice_boundary
        call mpp_get_compute_domain( Ice%domain, is, ie, js, je )
        kd = size(Ice%ice_mask,3)
        allocate( atmos_ice_boundary%u_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%v_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%u_star(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%t_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%q_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%lw_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux_vis_dir(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux_vis_dif(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux_nir_dir(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux_nir_dif(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%lprec(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%fprec(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%dhdt(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%dedt(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%drdt(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%coszen(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%p(is:ie,js:je,kd) )
! initialize boundary values for override experiments (mjh)
        atmos_ice_boundary%u_flux=0.0
        atmos_ice_boundary%v_flux=0.0
        atmos_ice_boundary%u_star=0.0
        atmos_ice_boundary%t_flux=0.0
        atmos_ice_boundary%q_flux=0.0
        atmos_ice_boundary%lw_flux=0.0
        atmos_ice_boundary%sw_flux_vis_dir=0.0
        atmos_ice_boundary%sw_flux_vis_dif=0.0
        atmos_ice_boundary%sw_flux_nir_dir=0.0
        atmos_ice_boundary%sw_flux_nir_dif=0.0
        atmos_ice_boundary%lprec=0.0
        atmos_ice_boundary%fprec=0.0
        atmos_ice_boundary%dhdt=0.0
        atmos_ice_boundary%dedt=0.0
        atmos_ice_boundary%drdt=0.0
        atmos_ice_boundary%coszen=0.0
        atmos_ice_boundary%p=0.0

!         allocate fields for extra fluxes
! Copying initialized gas fluxes from exchange grid to atmosphere_ice boundary

        call coupler_type_copy(ex_gas_fluxes, atmos_ice_boundary%fluxes, is, ie, js, je, kd,    &
             mod_name, Ice%axes, Time, suffix = '_atm_ice')
  
!allocate land_ice_boundary
        allocate( land_ice_boundary%runoff(is:ie,js:je) )
        allocate( land_ice_boundary%calving(is:ie,js:je) )
        allocate( land_ice_boundary%runoff_hflx(is:ie,js:je) )
        allocate( land_ice_boundary%calving_hflx(is:ie,js:je) )
! initialize values for override experiments (mjh)
        land_ice_boundary%runoff=0.0
        land_ice_boundary%calving=0.0
        land_ice_boundary%runoff_hflx=0.0
        land_ice_boundary%calving_hflx=0.0
!allocate land_ice_atmos_boundary
        call mpp_get_compute_domain( Atm%domain, is, ie, js, je )
        allocate( land_ice_atmos_boundary%t(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_vis_dir(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_nir_dir(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_vis_dif(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_nir_dif(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%land_frac(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dt_t(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dt_tr(is:ie,js:je,n_atm_tr) )
        allocate( land_ice_atmos_boundary%u_flux(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%v_flux(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dtaudu(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dtaudv(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%u_star(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%b_star(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%q_star(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%rough_mom(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%frac_open_sea(is:ie,js:je) )
! initialize boundary values for override experiments (mjh)
        land_ice_atmos_boundary%t=273.0
        land_ice_atmos_boundary%albedo=0.0
        land_ice_atmos_boundary%albedo_vis_dir=0.0
        land_ice_atmos_boundary%albedo_nir_dir=0.0
        land_ice_atmos_boundary%albedo_vis_dif=0.0
        land_ice_atmos_boundary%albedo_nir_dif=0.0
        land_ice_atmos_boundary%land_frac=0.0
        land_ice_atmos_boundary%dt_t=0.0
        land_ice_atmos_boundary%dt_tr=0.0
        land_ice_atmos_boundary%u_flux=0.0
        land_ice_atmos_boundary%v_flux=0.0
        land_ice_atmos_boundary%dtaudu=0.0
        land_ice_atmos_boundary%dtaudv=0.0
        land_ice_atmos_boundary%u_star=0.0
        land_ice_atmos_boundary%b_star=0.0
        land_ice_atmos_boundary%q_star=0.0
        land_ice_atmos_boundary%rough_mom=0.01
        land_ice_atmos_boundary%frac_open_sea=0.0

! allocate fields for extra tracers
! The first call is no longer necessary, the fluxes will be passed by the land module
! The 2nd call is useful in the case of a ocean model only simulation
!
        call coupler_type_copy(ex_gas_fields_atm, Atm%fields, is, ie, js, je,                   &
             mod_name, Atm%axes(1:2), Time, suffix = '_atm')

!Balaji: clocks on atm%pe only        
    cplClock = mpp_clock_id( 'Land-ice-atm coupler', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    sfcClock = mpp_clock_id( 'SFC boundary layer', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
    fluxAtmDnClock = mpp_clock_id( 'Flux DN from atm', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    fluxLandIceClock = mpp_clock_id( 'Flux land to ice', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    regenClock = mpp_clock_id( 'XGrid generation', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    fluxAtmUpClock = mpp_clock_id( 'Flux UP to atm', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    end if

    !--- With the consideration of concurrent and series run. Also make sure pelist is monotonically increasing.
    !--- Here we can not simply call mpp_set_current_pelist() because of ensemble. The ocean_pe(n) 
    !--- should either equal to atmos_pe(n) or greater than atmos%pelist(npes_atm)
    npes_ocn = size(Ocean%pelist(:))
    npes_atm = size(Atm%pelist(:))      
    allocate(pelist(npes_ocn+npes_atm))
    pelist(1:npes_atm) = Atm%pelist(1:npes_atm)
    npes_all = npes_atm
    do n = 1, npes_ocn
       if( n <= npes_atm ) then
          if( Ocean%pelist(n) == Atm%pelist(n) ) cycle
       endif
       if( Ocean%pelist(n) < Atm%pelist(npes_atm) ) call mpp_error( FATAL, &
           'flux_exchange_init: ocean%pelist(n) should equal to atm%pelist(n) or greater than any atmos pes' )
       npes_all = npes_all + 1
       pelist(npes_all) = Ocean%pelist(n)
    enddo

    call mpp_set_current_pelist(pelist(1:npes_all) )
    deallocate(pelist)

!ocean_ice_boundary and ice_ocean_boundary must be done on all PES
!domain boundaries will assure no space is allocated on non-relevant PEs.
    call mpp_get_compute_domain( Ice%domain, is, ie, js, je )
!allocate ocean_ice_boundary
    allocate( ocean_ice_boundary%u(is:ie,js:je) )
    allocate( ocean_ice_boundary%v(is:ie,js:je) )
    allocate( ocean_ice_boundary%t(is:ie,js:je) )
    allocate( ocean_ice_boundary%s(is:ie,js:je) )
!frazil and sea_level are optional, if not present they should be nullified
    allocate( ocean_ice_boundary%frazil(is:ie,js:je) )
    allocate( ocean_ice_boundary%sea_level(is:ie,js:je) )
! initialize boundary fields for override experiments (mjh)
    ocean_ice_boundary%u=0.0
    ocean_ice_boundary%v=0.0
    ocean_ice_boundary%t=273.0
    ocean_ice_boundary%s=0.0
    ocean_ice_boundary%frazil=0.0
    ocean_ice_boundary%sea_level=0.0

!
! allocate fields for extra tracers
! Copying gas flux fields from ice to ocean_ice boundary

    call coupler_type_copy(ex_gas_fields_ice, ocean_ice_boundary%fields, is, ie, js, je,        &
         'ice_flux', Ice%axes(1:2), Time, suffix = '_ocn_ice')

!allocate ice_ocean_boundary
    call mpp_get_compute_domain( Ocean%domain, is, ie, js, je )
!ML ocean only requires t, q, lw, sw, fprec, calving
!AMIP ocean needs no input fields
!choice of fields will eventually be done at runtime
!via field_manager
    allocate( ice_ocean_boundary%u_flux   (is:ie,js:je) ) ; ice_ocean_boundary%u_flux = 0.0
    allocate( ice_ocean_boundary%v_flux   (is:ie,js:je) ) ; ice_ocean_boundary%v_flux = 0.0
    allocate( ice_ocean_boundary%t_flux   (is:ie,js:je) ) ; ice_ocean_boundary%t_flux = 0.0
    allocate( ice_ocean_boundary%q_flux   (is:ie,js:je) ) ; ice_ocean_boundary%q_flux = 0.0
    allocate( ice_ocean_boundary%salt_flux(is:ie,js:je) ) ; ice_ocean_boundary%salt_flux = 0.0
    allocate( ice_ocean_boundary%lw_flux  (is:ie,js:je) ) ; ice_ocean_boundary%lw_flux = 0.0
    allocate( ice_ocean_boundary%sw_flux_vis_dir  (is:ie,js:je) ) ; ice_ocean_boundary%sw_flux_vis_dir = 0.0
    allocate( ice_ocean_boundary%sw_flux_vis_dif  (is:ie,js:je) ) ; ice_ocean_boundary%sw_flux_vis_dif = 0.0
    allocate( ice_ocean_boundary%sw_flux_nir_dir  (is:ie,js:je) ) ; ice_ocean_boundary%sw_flux_nir_dir = 0.0
    allocate( ice_ocean_boundary%sw_flux_nir_dif  (is:ie,js:je) ) ; ice_ocean_boundary%sw_flux_nir_dif = 0.0
    allocate( ice_ocean_boundary%lprec    (is:ie,js:je) ) ; ice_ocean_boundary%lprec = 0.0
    allocate( ice_ocean_boundary%fprec    (is:ie,js:je) ) ; ice_ocean_boundary%fprec = 0.0
    allocate( ice_ocean_boundary%runoff   (is:ie,js:je) ) ; ice_ocean_boundary%runoff = 0.0
    allocate( ice_ocean_boundary%calving  (is:ie,js:je) ) ; ice_ocean_boundary%calving = 0.0
    allocate( ice_ocean_boundary%runoff_hflx   (is:ie,js:je) ) ; ice_ocean_boundary%runoff_hflx = 0.0
    allocate( ice_ocean_boundary%calving_hflx  (is:ie,js:je) ) ; ice_ocean_boundary%calving_hflx = 0.0
    allocate( ice_ocean_boundary%p        (is:ie,js:je) ) ; ice_ocean_boundary%p = 0.0
    allocate( ice_ocean_boundary%mi       (is:ie,js:je) ) ; ice_ocean_boundary%mi = 0.0
    allocate( ice_ocean_boundary%wnd      (is:ie,js:je) ) ;         ice_ocean_boundary%wnd = 0.0

!
! allocate fields for extra tracers
!

    call coupler_type_copy(ex_gas_fluxes, ice_ocean_boundary%fluxes, is, ie, js, je,    &
         'ocean_flux', Ocean%axes(1:2), Time, suffix = '_ice_ocn')

    call coupler_type_copy(ex_gas_fields_ice, Ocean%fields, is, ie, js, je,             &
         'ocean_flux', Ocean%axes(1:2), Time, suffix = '_ocn')

!pjp Why are the above not initialized to zero?
! initialize boundary values for override experiments
    ocean_ice_boundary%xtype = REDIST
    if( Ocean%domain.EQ.Ice%domain )ocean_ice_boundary%xtype = DIRECT
    ice_ocean_boundary%xtype = ocean_ice_boundary%xtype


!
! allocate fields amd fluxes for extra tracers for the Ice type
!

    call mpp_get_compute_domain( Ice%domain, is, ie, js, je )
    kd = size(Ice%ice_mask,3)
    call coupler_type_copy(ex_gas_fields_ice, Ice%ocean_fields, is, ie, js, je, kd,     &
         'ice_flux', Ice%axes, Time, suffix = '_ice')

    call coupler_type_copy(ex_gas_fluxes, Ice%ocean_fluxes, is, ie, js, je,             &
         'ice_flux', Ice%axes(1:2), Time, suffix = '_ice')

    call coupler_type_copy(ex_gas_fluxes, Ice%ocean_fluxes_top, is, ie, js, je, kd,     &
         'ice_flux', Ice%axes, Time, suffix = '_ice_top')

!       initialize the Ocean type for extra fields for surface fluxes
! Same allocation of arrays and stuff
!       (this must be done after the Ocean fields are allocated as the fields on the Ocean%fields
!       are read in in this subroutine)
!

    if ( Ocean%is_ocean_pe ) then
      call mpp_set_current_pelist(Ocean%pelist)
      call ocean_model_init_sfc(Ocean_state, Ocean)
    end if
    call mpp_set_current_pelist()

    if( Ocean%is_ocean_pe) then
       call mpp_get_compute_domain(Ocean%domain, xsize=nxc_ocn, ysize=nyc_ocn)
    endif
    if( Ice%pe) then
       call mpp_get_compute_domain(Ice%domain, xsize=nxc_ice, ysize=nyc_ice)
       nk_ice = size(Ice%part_size,3)
    endif

 

    ! required by stock_move, all fluxes used to update stocks will be zero if dt_atmos,
    ! and dt_cpld are absent
    Dt_atm = 0
    Dt_cpl = 0
    if(present(dt_atmos)) Dt_atm = dt_atmos
    if(present(dt_cpld )) Dt_cpl = dt_cpld
 
    !z1l check the flux conservation.
    if(debug_stocks) call check_flux_conservation(Ice, Ocean, Ice_Ocean_Boundary)


!Balaji
    cplOcnClock = mpp_clock_id( 'Ice-ocean coupler', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    fluxIceOceanClock = mpp_clock_id( 'Flux ice to ocean', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    fluxOceanIceClock = mpp_clock_id( 'Flux ocean to ice', flags=clock_flag_default, grain=CLOCK_ROUTINE )
!---- done ----
    do_init = .false.

  end subroutine flux_exchange_init
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="sfc_boundary_layer">
!  <OVERVIEW>
!   Computes explicit fluxes as well as derivatives that will be used to compute an implicit flux correction. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!  The following quantities in the land_ice_atmos_boundary_type are computed:
!
!     
!         t_surf_atm = surface temperature (used for radiation)    (K)
!         albedo_atm = surface albedo      (used for radiation)    (nondimensional)
!      rough_mom_atm = surface roughness for momentum (m)
!      land_frac_atm = fractional area of land beneath an atmospheric
!                      grid box 
!         dtaudu_atm, dtaudv_atm = derivatives of wind stress w.r.t. the
!                      lowest level wind speed  (Pa/(m/s))
!         flux_u_atm = zonal wind stress  (Pa)
!         flux_v_atm = meridional wind stress (Pa)
!         u_star_atm = friction velocity (m/s)
!         b_star_atm = buoyancy scale    (m2/s)
!
!         (u_star and b_star are defined so that u_star**2 = magnitude
!           of surface stress divided by density of air at the surface, 
!           and u_star*b_star = buoyancy flux at the surface)
!
!   </PRE>
!  </DESCRIPTION>

!  <TEMPLATE>
!   call sfc_boundary_layer ( dt, Time, Atm, Land, Ice, Boundary )
!		
!  </TEMPLATE>
!  <IN NAME=" dt" TYPE="real">
!   time step. 
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <INOUT NAME="Atm" TYPE="atmos_data_type">
!   A derived data type to specify atmosphere boundary data.
!  </INOUT>
!  <INOUT NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </INOUT>
!  <INOUT NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </INOUT>
!  <INOUT NAME="Boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere, land and ice.   
!  </INOUT>
!
subroutine sfc_boundary_layer ( dt, Time, Atm, Land, Ice, Land_Ice_Atmos_Boundary )

  real,                  intent(in)  :: dt
  type(time_type),       intent(in)  :: Time
  type(atmos_data_type), intent(inout)  :: Atm
  type(land_data_type),  intent(inout)  :: Land
  type(ice_data_type),   intent(inout)  :: Ice
  type(land_ice_atmos_boundary_type), intent(inout) :: Land_Ice_Atmos_Boundary

  ! ---- local vars ----------------------------------------------------------
  real, dimension(n_xgrid_sfc) :: &
       ex_albedo,     &
       ex_albedo_vis_dir,     &
       ex_albedo_nir_dir,     &
       ex_albedo_vis_dif,     &
       ex_albedo_nir_dif,     &
       ex_land_frac,  &
       ex_t_atm,      & 
       ex_p_atm,      &
       ex_u_atm, ex_v_atm,    &
       ex_gust,       &
       ex_t_surf4,    &
       ex_u_surf, ex_v_surf,  &
       ex_rough_mom, ex_rough_heat, ex_rough_moist, &
       ex_rough_scale,&
       ex_q_star,     &
       ex_cd_q,       &
       ex_ref, ex_ref_u, ex_ref_v, ex_u10, &
       ex_ref2,       &
       ex_t_ref,      &
       ex_qs_ref,     &
       ex_qs_ref_cmip,     &
       ex_del_m,      &
       ex_del_h,      &
       ex_del_q,      &
       ex_seawater,   &
       ex_frac_open_sea

  real, dimension(n_xgrid_sfc,n_exch_tr) :: ex_tr_atm
! jgj: added for co2_atm diagnostic
  real, dimension(n_xgrid_sfc)           :: ex_co2_atm_dvmr
  real, dimension(size(Land_Ice_Atmos_Boundary%t,1),size(Land_Ice_Atmos_Boundary%t,2)) :: diag_atm
  real, dimension(size(Land%t_ca, 1),size(Land%t_ca,2), size(Land%t_ca,3)) :: diag_land
  real, dimension(size(Ice%t_surf,1),size(Ice%t_surf,2),size(Ice%t_surf,3)) :: sea
  real, dimension(size(Ice%albedo,1),size(Ice%albedo,2),size(Ice%albedo,3)) ::  tmp_open_sea
  real    :: zrefm, zrefh
  logical :: used
  character(32) :: tr_name ! tracer name
  integer :: tr, n, m ! tracer indices
  integer :: i, ind_flux = 1

  ! [1] check that the module was initialized
!   <ERROR MSG="must call flux_exchange_init first " STATUS="FATAL">
!      flux_exchange_init has not been called before calling sfc_boundary_layer.
!   </ERROR>
  if (do_init) call error_mesg ('flux_exchange_mod',  &
       'must call flux_exchange_init first', FATAL)
!Balaji
  call mpp_clock_begin(cplClock)
  call mpp_clock_begin(sfcClock)
  ! [2] allocate storage for variables that are also used in flux_up_to_atmos
  allocate ( &
       ex_t_surf   (n_xgrid_sfc),  &
       ex_t_surf_miz(n_xgrid_sfc), &
       ex_p_surf   (n_xgrid_sfc),  &
       ex_slp      (n_xgrid_sfc),  &
       ex_t_ca     (n_xgrid_sfc),  &
       ex_dhdt_surf(n_xgrid_sfc),  &
       ex_dedt_surf(n_xgrid_sfc),  &
       ex_dqsatdt_surf(n_xgrid_sfc),  &
       ex_drdt_surf(n_xgrid_sfc),  &
       ex_dhdt_atm (n_xgrid_sfc),  &
       ex_flux_t   (n_xgrid_sfc),  &
       ex_flux_lw  (n_xgrid_sfc),  &
       ex_drag_q   (n_xgrid_sfc),  &
       ex_avail    (n_xgrid_sfc),  &
       ex_f_t_delt_n(n_xgrid_sfc), &

       ex_tr_surf     (n_xgrid_sfc, n_exch_tr), &
       ex_dfdtr_surf  (n_xgrid_sfc, n_exch_tr), &
       ex_dfdtr_atm   (n_xgrid_sfc, n_exch_tr), &
       ex_flux_tr     (n_xgrid_sfc, n_exch_tr), &
       ex_f_tr_delt_n (n_xgrid_sfc, n_exch_tr), &
       ex_e_tr_n      (n_xgrid_sfc, n_exch_tr), &

! MOD these were moved from local ! so they can be passed to flux down
       ex_flux_u(n_xgrid_sfc),    &
       ex_flux_v(n_xgrid_sfc),    &
       ex_dtaudu_atm(n_xgrid_sfc),&
       ex_dtaudv_atm(n_xgrid_sfc),&

! values added for LM3
       ex_cd_t     (n_xgrid_sfc),  &
       ex_cd_m     (n_xgrid_sfc),  &
       ex_b_star   (n_xgrid_sfc),  &
       ex_u_star   (n_xgrid_sfc),  &
       ex_wind     (n_xgrid_sfc),  &
       ex_z_atm    (n_xgrid_sfc),  &

       ex_e_t_n    (n_xgrid_sfc),  &
       ex_e_q_n    (n_xgrid_sfc),  &
       ex_land     (n_xgrid_sfc)   )

#ifdef SCM
  allocate ( &
       ex_dhdt_surf_forland(n_xgrid_sfc), &
       ex_dedt_surf_forland(n_xgrid_sfc), &
       ex_dedq_surf_forland(n_xgrid_sfc)  )
#endif

  ex_p_surf = 1
! Actual allocation of exchange fields for ocean_ice boundary
  do n = 1, ex_gas_fields_ice%num_bcs  !{
    do m = 1, ex_gas_fields_ice%bc(n)%num_fields  !{
      if (associated(ex_gas_fields_ice%bc(n)%field(m)%values)) then  !{
        call mpp_error( FATAL, 'sfc_boundary_layer: ex_gas_fields_ice already allocated.' )
      endif  !}
      allocate ( ex_gas_fields_ice%bc(n)%field(m)%values(n_xgrid_sfc) )
      ex_gas_fields_ice%bc(n)%field(m)%values = 0.0
    enddo  !} m
  enddo  !} n

  do n = 1, ex_gas_fields_atm%num_bcs  !{
    do m = 1, ex_gas_fields_atm%bc(n)%num_fields  !{
      if (associated(ex_gas_fields_atm%bc(n)%field(m)%values)) then  !{
        call mpp_error( FATAL, 'sfc_boundary_layer: ex_gas_fields_atm already allocated.' )
      endif  !}
      allocate ( ex_gas_fields_atm%bc(n)%field(m)%values(n_xgrid_sfc) )
      ex_gas_fields_atm%bc(n)%field(m)%values = 0.0
    enddo  !} m
  enddo  !} n

  do n = 1, ex_gas_fluxes%num_bcs  !{
    do m = 1, ex_gas_fluxes%bc(n)%num_fields  !{
      if (associated(ex_gas_fluxes%bc(n)%field(m)%values)) then  !{
        call mpp_error( FATAL, 'sfc_boundary_layer: ex_gas_fluxes already allocated.' )
      endif  !}
      allocate ( ex_gas_fluxes%bc(n)%field(m)%values(n_xgrid_sfc) )
      ex_gas_fluxes%bc(n)%field(m)%values = 0.0
    enddo  !} m
  enddo  !} n

!
!       Call the atmosphere tracer driver to gather the data needed for extra gas tracers
! For ocean only model

!  call atmos_get_fields_for_flux(Atm)

  ! [3] initialize some values on exchange grid: this is actually a safeguard
  ! against using undefined values
  ex_t_surf   = 200.
  ex_u_surf   =   0.
  ex_v_surf   =   0.
  ex_albedo = 0. ! bw 
  ex_albedo_vis_dir = 0.
  ex_albedo_nir_dir = 0.
  ex_albedo_vis_dif = 0.
  ex_albedo_nir_dif = 0.

  !---- do not use if relax time /= 0 ----
  ex_cd_t = 0.0
  ex_cd_m = 0.0
  ex_cd_q = 0.0
  ex_frac_open_sea =0.
!-----------------------------------------------------------------------
!Balaji: data_override stuff moved from coupler_main
  call data_override ('ATM', 't_bot',  Atm%t_bot , Time)
  call data_override ('ATM', 'z_bot',  Atm%z_bot , Time)
  call data_override ('ATM', 'p_bot',  Atm%p_bot , Time)
  call data_override ('ATM', 'u_bot',  Atm%u_bot , Time)
  call data_override ('ATM', 'v_bot',  Atm%v_bot , Time)
  call data_override ('ATM', 'p_surf', Atm%p_surf, Time)
  call data_override ('ATM', 'slp',    Atm%slp,    Time)
  call data_override ('ATM', 'gust',   Atm%gust,   Time)
!
! jgj: 2008/07/18 
! FV atm advects tracers in moist mass mixing ratio: kg co2 /(kg air + kg water)
! cubed sphere advects moist mass mixing ratio also (per SJ)
! data table co2 overrides for ocean (co2_flux_pcair_atm)
! and land (co2_bot) should be in dry vmr (mol/mol) units.
!  ATM: co2_flux_pcair_atm : to override atm_btm layer to send to ocean
!  ATM: co2_bot            : to override atm_btm layer to send to land

! data override for co2 to be passed to land/photosynthesis (co2_bot)
! land co2 data override is in dry_vmr units, so convert to wet_mmr for land model.
! co2mmr = (wco2/wair) * co2vmr;  wet_mmr = dry_mmr * (1-Q)
!
  do tr = 1,n_atm_tr
     call get_tracer_names( MODEL_ATMOS, tr, tr_name )
     call data_override('ATM', trim(tr_name)//'_bot', Atm%tr_bot(:,:,tr), Time, override=used)
! conversion for land co2 data override from dry vmr to moist mmr
     if (used .and. lowercase(trim(tr_name)).eq.'co2') then
       Atm%tr_bot(:,:,tr) = Atm%tr_bot(:,:,tr) * (WTMCO2/WTMAIR) *    &
                            (1.0 - Atm%tr_bot(:,:,isphum))
     end if
  enddo
! data override for co2 to be passed to ocean (co2_flux_pcair_atm) 
! atmos_co2.F90 already called: converts tr_bot passed to ocean via gas_flux   
! from moist mmr to dry vmr.
  do n = 1, atm%fields%num_bcs  !{
    do m = 1, atm%fields%bc(n)%num_fields  !{
      call data_override('ATM', atm%fields%bc(n)%field(m)%name,      &
           atm%fields%bc(n)%field(m)%values, Time, override = atm%fields%bc(n)%field(m)%override)
      ex_gas_fields_atm%bc(n)%field(m)%override = atm%fields%bc(n)%field(m)%override
    enddo  !} m
  enddo  !} n
  do n = 1, atm%fields%num_bcs  !{
     if (atm%fields%bc(n)%use_atm_pressure) then  !{
        if (.not. atm%fields%bc(n)%field(ind_psurf)%override) then  !{
           atm%fields%bc(n)%field(ind_psurf)%values = Atm%p_surf
        endif  !}
     endif  !}
  enddo  !} n
  call data_override ('ICE', 't_surf',     Ice%t_surf,      Time)
  call data_override ('ICE', 'rough_mom',  Ice%rough_mom,   Time)
  call data_override ('ICE', 'rough_heat', Ice%rough_heat,  Time)
  call data_override ('ICE', 'rough_moist',Ice%rough_moist, Time)
  call data_override ('ICE', 'albedo',     Ice%albedo,      Time)
  call data_override ('ICE', 'albedo_vis_dir', Ice%albedo_vis_dir, Time)
  call data_override ('ICE', 'albedo_nir_dir', Ice%albedo_nir_dir, Time)
  call data_override ('ICE', 'albedo_vis_dif', Ice%albedo_vis_dif, Time)
  call data_override ('ICE', 'albedo_nir_dif', Ice%albedo_nir_dif, Time)
  call data_override ('ICE', 'u_surf',     Ice%u_surf,      Time)
  call data_override ('ICE', 'v_surf',     Ice%v_surf,      Time)
  call data_override ('LND', 't_surf',     Land%t_surf,     Time)
  call data_override ('LND', 't_ca',       Land%t_ca,       Time)
  call data_override ('LND', 'rough_mom',  Land%rough_mom,  Time)
  call data_override ('LND', 'rough_heat', Land%rough_heat, Time)
  call data_override ('LND', 'albedo', Land%albedo,     Time)

! tracer data override
  do tr = 1, n_lnd_tr
     call get_tracer_names( MODEL_LAND, tr, tr_name )
     call data_override('LND', trim(tr_name)//'_surf', Land%tr(:,:,:,tr), Time)
  enddo
  do n = 1, ice%ocean_fields%num_bcs  !{
    do m = 1, ice%ocean_fields%bc(n)%num_fields  !{
      call data_override('ICE', ice%ocean_fields%bc(n)%field(m)%name, ice%ocean_fields%bc(n)%field(m)%values, Time)
      if ( Ice%ocean_fields%bc(n)%field(m)%id_diag > 0 ) then  !{
        used = send_data(Ice%ocean_fields%bc(n)%field(m)%id_diag, Ice%ocean_fields%bc(n)%field(m)%values, Time )
      endif  !}
    enddo  !} m
  enddo  !} n
  call data_override ('LND', 'albedo_vis_dir', Land%albedo_vis_dir,Time)
  call data_override ('LND', 'albedo_nir_dir', Land%albedo_nir_dir,Time)
  call data_override ('LND', 'albedo_vis_dif', Land%albedo_vis_dif,Time)
  call data_override ('LND', 'albedo_nir_dif', Land%albedo_nir_dif,Time)

!---- put atmosphere quantities onto exchange grid ----

  ! [4] put all the qantities we need onto exchange grid
  ! [4.1] put atmosphere quantities onto exchange grid
  if (do_forecast) then
    call put_to_xgrid (Atm%Surf_diff%sst_miz , 'ATM', ex_t_surf_miz, xmap_sfc, remap_method=remap_method, complete=.false.)
  endif
! put atmosphere bottom layer tracer data onto exchange grid
  do tr = 1,n_exch_tr
     call put_to_xgrid (Atm%tr_bot(:,:,tr_table(tr)%atm) , 'ATM', ex_tr_atm(:,tr), xmap_sfc, &
          remap_method=remap_method, complete=.false.)
  enddo
  do n = 1, Atm%fields%num_bcs  !{
    do m = 1, Atm%fields%bc(n)%num_fields  !{
      call put_to_xgrid (Atm%fields%bc(n)%field(m)%values, 'ATM',            &
           ex_gas_fields_atm%bc(n)%field(m)%values, xmap_sfc, remap_method=remap_method, complete=.false.)
    enddo  !} m
  enddo  !} n

  call put_to_xgrid (Atm%t_bot , 'ATM', ex_t_atm , xmap_sfc, remap_method=remap_method, complete=.false.)
  call put_to_xgrid (Atm%z_bot , 'ATM', ex_z_atm , xmap_sfc, remap_method=remap_method, complete=.false.)
  call put_to_xgrid (Atm%p_bot , 'ATM', ex_p_atm , xmap_sfc, remap_method=remap_method, complete=.false.)
  call put_to_xgrid (Atm%u_bot , 'ATM', ex_u_atm , xmap_sfc, remap_method=remap_method, complete=.false.)
  call put_to_xgrid (Atm%v_bot , 'ATM', ex_v_atm , xmap_sfc, remap_method=remap_method, complete=.false.)
  call put_to_xgrid (Atm%p_surf, 'ATM', ex_p_surf, xmap_sfc, remap_method=remap_method, complete=.false.)
  call put_to_xgrid (Atm%slp,    'ATM', ex_slp,    xmap_sfc, remap_method=remap_method, complete=.false.)
  call put_to_xgrid (Atm%gust,   'ATM', ex_gust,   xmap_sfc, remap_method=remap_method, complete=.true.)

  ! slm, Mar 20 2002: changed order in whith the data transferred from ice and land 
  ! grids, to fill t_ca first with t_surf over ocean and then with t_ca from 
  ! land, where it is different from t_surf. It is mostly to simplify 
  ! diagnostic, since surface_flux calculations distinguish between land and 
  ! not-land anyway.

  ! prefill surface values with atmospheric values before putting tracers
  ! from ice or land, so that gradient is 0 if tracers are not filled
  ex_tr_surf = ex_tr_atm

  ! [4.2] put ice quantities onto exchange grid
  ! (assume that ocean quantites are stored in no ice partition)
  ! (note: ex_avail is true at ice and ocean points)
  call put_to_xgrid (Ice%t_surf,      'OCN', ex_t_surf,      xmap_sfc)
  call put_to_xgrid (Ice%rough_mom,   'OCN', ex_rough_mom,   xmap_sfc)
  call put_to_xgrid (Ice%rough_heat,  'OCN', ex_rough_heat,  xmap_sfc)
  call put_to_xgrid (Ice%rough_moist, 'OCN', ex_rough_moist, xmap_sfc)
  call put_to_xgrid (Ice%albedo,      'OCN', ex_albedo,      xmap_sfc)
  call put_to_xgrid (Ice%albedo_vis_dir, 'OCN', ex_albedo_vis_dir, xmap_sfc)
  call put_to_xgrid (Ice%albedo_nir_dir, 'OCN', ex_albedo_nir_dir, xmap_sfc)
  call put_to_xgrid (Ice%albedo_vis_dif, 'OCN', ex_albedo_vis_dif, xmap_sfc)
  call put_to_xgrid (Ice%albedo_nir_dif, 'OCN', ex_albedo_nir_dif, xmap_sfc)
  call put_to_xgrid (Ice%u_surf,      'OCN', ex_u_surf,      xmap_sfc)
  call put_to_xgrid (Ice%v_surf,      'OCN', ex_v_surf,      xmap_sfc)

  tmp_open_sea        = 0.
  tmp_open_sea(:,:,1) = 1.
  call put_to_xgrid ( tmp_open_sea,  'OCN', ex_frac_open_sea,   xmap_sfc)

  do n = 1, ice%ocean_fields%num_bcs  !{
    do m = 1, ice%ocean_fields%bc(n)%num_fields  !{
      call put_to_xgrid (Ice%ocean_fields%bc(n)%field(m)%values, 'OCN',      &
           ex_gas_fields_ice%bc(n)%field(m)%values, xmap_sfc)
    enddo  !} m
  enddo  !} n
  sea = 0.0; sea(:,:,1) = 1.0;
  ex_seawater = 0.0
  call put_to_xgrid (sea,             'OCN', ex_seawater,    xmap_sfc)
  ex_t_ca = ex_t_surf ! slm, Mar 20 2002 to define values over the ocean

  ! [4.3] put land quantities onto exchange grid ----
  call some(xmap_sfc, ex_land, 'LND')
  if (do_forecast) then
    call put_to_xgrid (Land%t_surf,     'LND', ex_t_surf_miz,  xmap_sfc)
    ex_t_ca(:) = ex_t_surf_miz(:)
  end if

  call put_to_xgrid (Land%t_surf,     'LND', ex_t_surf,      xmap_sfc)
  call put_to_xgrid (Land%t_ca,       'LND', ex_t_ca,        xmap_sfc)
  call put_to_xgrid (Land%rough_mom,  'LND', ex_rough_mom,   xmap_sfc)
  call put_to_xgrid (Land%rough_heat, 'LND', ex_rough_heat,  xmap_sfc)
  call put_to_xgrid (Land%rough_heat, 'LND', ex_rough_moist, xmap_sfc)
  call put_to_xgrid (Land%albedo,     'LND', ex_albedo,      xmap_sfc)
  call put_to_xgrid (Land%albedo_vis_dir,     'LND', ex_albedo_vis_dir,   xmap_sfc)
  call put_to_xgrid (Land%albedo_nir_dir,     'LND', ex_albedo_nir_dir,   xmap_sfc)
  call put_to_xgrid (Land%albedo_vis_dif,     'LND', ex_albedo_vis_dif,   xmap_sfc)
  call put_to_xgrid (Land%albedo_nir_dif,     'LND', ex_albedo_nir_dif,   xmap_sfc)
  ex_rough_scale = ex_rough_mom
  call put_to_xgrid(Land%rough_scale, 'LND', ex_rough_scale, xmap_sfc)
 
  do tr = 1,n_exch_tr
     n = tr_table(tr)%lnd
     if(n /= NO_TRACER ) then
        call put_to_xgrid ( Land%tr(:,:,:,n), 'LND', ex_tr_surf(:,tr), xmap_sfc )
     else
        ! do nothing, since ex_tr_surf is prefilled with ex_tr_atm, and therefore
        ! fluxes will be 0
     endif
  enddo

  ex_land_frac = 0.0
  call put_logical_to_real (Land%mask,    'LND', ex_land_frac, xmap_sfc)

#ifdef SCM
  if (do_specified_land) then
       if (do_specified_albedo) then
            ex_albedo = ALBEDO_OBS
            ex_albedo_vis_dir = ALBEDO_OBS
            ex_albedo_nir_dir = ALBEDO_OBS
            ex_albedo_vis_dif = ALBEDO_OBS
            ex_albedo_nir_dif = ALBEDO_OBS
       endif
       if (do_specified_tskin) then
            ex_t_surf = TSKIN
            ex_t_ca   = TSKIN
            ex_tr_surf(:,isphum) = 15.e-3
       endif
       if (do_specified_rough_leng) then
            ex_rough_mom   = ROUGH_MOM
            ex_rough_heat  = ROUGH_HEAT
            ex_rough_moist = ROUGH_HEAT
       endif
  endif
#endif

  if (do_forecast) then
     ex_t_surf = ex_t_surf_miz
  end if

  ! [5] compute explicit fluxes and tendencies at all available points ---
  call some(xmap_sfc, ex_avail)
  call surface_flux (&
       ex_t_atm, ex_tr_atm(:,isphum),  ex_u_atm, ex_v_atm,  ex_p_atm,  ex_z_atm,  &
       ex_p_surf,ex_t_surf, ex_t_ca,  ex_tr_surf(:,isphum),                       &
       ex_u_surf, ex_v_surf,                                           &
       ex_rough_mom, ex_rough_heat, ex_rough_moist, ex_rough_scale,    &
       ex_gust,                                                        &
       ex_flux_t, ex_flux_tr(:,isphum), ex_flux_lw, ex_flux_u, ex_flux_v,         &
       ex_cd_m,   ex_cd_t, ex_cd_q,                                    &
       ex_wind,   ex_u_star, ex_b_star, ex_q_star,                     &
       ex_dhdt_surf, ex_dedt_surf, ex_dfdtr_surf(:,isphum),  ex_drdt_surf,        &
       ex_dhdt_atm,  ex_dfdtr_atm(:,isphum),  ex_dtaudu_atm, ex_dtaudv_atm,       &
       dt,                                                             &
       ex_land, ex_seawater .gt. 0,  ex_avail                          )

#ifdef SCM
! Option to override surface fluxes for SCM
  if (do_specified_flux) then

    call scm_surface_flux ( &
       ex_t_atm, ex_tr_atm(:,isphum),  ex_u_atm, ex_v_atm,  ex_p_atm,  ex_z_atm,  &
       ex_p_surf,ex_t_surf, ex_t_ca,  ex_tr_surf(:,isphum),                       &
       ex_u_surf, ex_v_surf,                                                      &
       ex_rough_mom, ex_rough_heat, ex_rough_moist, ex_rough_scale,               &
       ex_gust,                                                                   &
       ex_flux_t, ex_flux_tr(:,isphum), ex_flux_lw, ex_flux_u, ex_flux_v,         &
       ex_cd_m,   ex_cd_t, ex_cd_q,                                               &
       ex_wind,   ex_u_star, ex_b_star, ex_q_star,                                &
       ex_dhdt_surf, ex_dedt_surf, ex_dfdtr_surf(:,isphum),  ex_drdt_surf,        &
       ex_dhdt_atm,  ex_dfdtr_atm(:,isphum),  ex_dtaudu_atm, ex_dtaudv_atm,       &
       dt,                                                                        &
       ex_land, ex_seawater .gt. 0,  ex_avail,                                    &
       ex_dhdt_surf_forland,  ex_dedt_surf_forland,  ex_dedq_surf_forland  )

  endif
#endif

  zrefm = 10.0
  zrefh = z_ref_heat
  !      ---- optimize calculation ----
  call mo_profile ( zrefm, zrefh, ex_z_atm,   ex_rough_mom, &
       ex_rough_heat, ex_rough_moist,          &
       ex_u_star, ex_b_star, ex_q_star,        &
       ex_del_m, ex_del_h, ex_del_q, ex_avail  )
  ex_u10 = 0.
  where (ex_avail)
     ex_ref_u = ex_u_surf + (ex_u_atm-ex_u_surf) * ex_del_m 
     ex_ref_v = ex_v_surf + (ex_v_atm-ex_v_surf) * ex_del_m
     ex_u10 = sqrt(ex_ref_u**2 + ex_ref_v**2)
  endwhere
  do n = 1, ex_gas_fields_atm%num_bcs  !{
     if (atm%fields%bc(n)%use_10m_wind_speed) then  !{
        if (.not. ex_gas_fields_atm%bc(n)%field(ind_u10)%override) then  !{
           ex_gas_fields_atm%bc(n)%field(ind_u10)%values = ex_u10
        endif  !}
     endif  !}
  enddo  !} n
  ! fill derivatives for all tracers
  ! F = C0*u*rho*delta_q, C0*u*rho is the same for all tracers, copy from sphum
  do tr = 1,n_exch_tr
     if (tr==isphum) cycle
     ex_dfdtr_atm  (:,tr) = ex_dfdtr_atm  (:,isphum)
     ex_dfdtr_surf (:,tr) = ex_dfdtr_surf (:,isphum)
     ex_flux_tr    (:,tr) = ex_dfdtr_surf(:,tr)*(ex_tr_surf(:,tr)-ex_tr_atm(:,tr))
  enddo

! Combine explicit ocean flux and implicit land flux of extra flux fields.

  ! Calculate ocean explicit flux here

  call atmos_ocean_fluxes_calc(ex_gas_fields_atm, ex_gas_fields_ice, ex_gas_fluxes, ex_seawater)

  ! The following statement is a concise version of what's following and worth
  ! looking into in the future.
  ! ex_flux_tr(:,itracer) = ex_gas_fluxes%bc(itracer_ocn)%field(ind_flux)%values(:)
  ! where(ex_seawater.gt.0) ex_flux_tr(:,itracer) = F_ocn
  do n = 1, ex_gas_fluxes%num_bcs  !{
    if (ex_gas_fluxes%bc(n)%atm_tr_index .gt. 0) then  !{
      m = tr_table_map(ex_gas_fluxes%bc(n)%atm_tr_index)%exch
      do i = 1, size(ex_seawater(:))  !{
         if (ex_land(i)) cycle  ! over land, don't do anything
         ! on ocean or ice cells, flux is explicit therefore we zero derivatives. 
         ex_dfdtr_atm(i,m)  = 0.0
         ex_dfdtr_surf(i,m) = 0.0
         if (ex_seawater(i)>0) then
            ! jgj: convert to kg co2/m2/sec for atm
            ex_flux_tr(i,m)    = ex_gas_fluxes%bc(n)%field(ind_flux)%values(i) * ex_gas_fluxes%bc(n)%mol_wt * 1.0e-03
         else 
            ex_flux_tr(i,m) = 0.0 ! pure ice exchange cell
        endif  !}
      enddo  !} i
    endif  !}
  enddo  !} n

  ! [5.2] override tracer fluxes and derivatives
  do tr = 1,n_exch_tr
     if( tr_table(tr)%atm == NO_TRACER ) cycle ! it should never happen, though

     call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
     ! [5.2.1] override tracer flux. Note that "sea" and "diag_land" are repeatedly used 
     ! as temporary storage for the values we are overriding fluxes and derivative with, 
     ! over ocean and land respectively
     call data_override ( 'LND', 'ex_flux_'//trim(tr_name), diag_land, Time, override=used )
     if(used) call put_to_xgrid ( diag_land, 'LND', ex_flux_tr(:,tr), xmap_sfc )
     call data_override ( 'ICE', 'ex_flux_'//trim(tr_name), sea, Time, override=used )
     if(used) call put_to_xgrid ( sea, 'OCN', ex_flux_tr(:,tr), xmap_sfc )
     ! [5.2.2] override derivative of flux wrt surface concentration
     call data_override ( 'LND', 'ex_dfd'//trim(tr_name)//'_surf', diag_land, Time, override=used )
     if(used) call put_to_xgrid ( diag_land, 'LND', ex_dfdtr_surf(:,tr), xmap_sfc )
     call data_override ( 'ICE', 'ex_dfd'//trim(tr_name)//'_surf', sea, Time, override=used )
     if(used) call put_to_xgrid ( sea, 'OCN', ex_dfdtr_surf(:,tr), xmap_sfc )
     ! [5.2.3] override derivative of flux wrt atmospheric concentration
     call data_override ( 'LND', 'ex_dfd'//trim(tr_name)//'_atm', diag_land, Time, override=used )
     if(used) call put_to_xgrid ( diag_land, 'LND', ex_dfdtr_atm(:,tr), xmap_sfc )
     call data_override ( 'ICE', 'ex_dfd'//trim(tr_name)//'_atm', sea, Time, override=used )
     if(used) call put_to_xgrid ( sea, 'OCN', ex_dfdtr_atm(:,tr), xmap_sfc )
  enddo

  ! [5.3] override flux and derivatives for sensible heat flux
  ! [5.3.1] override flux
  call data_override ( 'LND', 'ex_flux_t', diag_land, Time, override=used )
  if (used) call put_to_xgrid ( diag_land, 'LND', ex_flux_t, xmap_sfc )
  call data_override ( 'ICE', 'ex_flux_t', sea, Time, override=used )
  if (used) call put_to_xgrid ( sea, 'OCN', ex_flux_t, xmap_sfc )
  ! [5.3.2] override derivative of flux wrt near-surface temperature
  call data_override ( 'LND', 'ex_dhdt_surf', diag_land, Time, override=used )
  if (used) call put_to_xgrid ( diag_land, 'LND', ex_dhdt_surf, xmap_sfc )
  call data_override ( 'ICE', 'ex_dhdt_surf', sea, Time, override=used )
  if (used) call put_to_xgrid ( sea, 'OCN', ex_dhdt_surf, xmap_sfc )
  ! [5.3.3] override derivative of flux wrt atmospheric temperature
  call data_override ( 'LND', 'ex_dhdt_atm', diag_land, Time,override=used )
  if (used) call put_to_xgrid ( diag_land, 'LND', ex_dhdt_atm, xmap_sfc )
  call data_override ( 'ICE', 'ex_dhdt_atm', sea, Time, override=used )
  if (used) call put_to_xgrid ( sea, 'OCN', ex_dhdt_atm, xmap_sfc )

  ! NB: names of the override fields are constructed using tracer name and certain 
  ! prefixes / suffixes. For example, for the tracer named "sphum" (specific humidity) they will be:
  ! "ex_flux_sphum", "ex_dfdsphum_surf", and "ex_dfdsphum_atm".
  ! 
  ! For sensible heat flux names are "ex_flux_t", "ex_dhdt_surf", and "ex_dhdt_atm"; 
  ! despite the name those are actually in energy units, W/m2, W/(m2 degK), and
  ! W/(m2 degK) respectively

  where (ex_avail) ex_drag_q = ex_wind*ex_cd_q
  ! [6] get mean quantities on atmosphere grid
  ! [6.1] compute t surf for radiation
  ex_t_surf4 = ex_t_surf ** 4

  ! [6.2] put relevant quantities onto atmospheric boundary
  call get_from_xgrid (Land_Ice_Atmos_Boundary%t,         'ATM', ex_t_surf4  ,  xmap_sfc, complete=.false.)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%frac_open_sea,'ATM',ex_frac_open_sea, xmap_sfc)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo,    'ATM', ex_albedo   ,  xmap_sfc, complete=.false.)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dir,    'ATM',   &
                       ex_albedo_vis_dir   ,  xmap_sfc, complete=.false.)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dir,    'ATM',   &
                       ex_albedo_nir_dir   ,  xmap_sfc, complete=.false.)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dif,    'ATM',   &
                       ex_albedo_vis_dif   ,  xmap_sfc, complete=.false.)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dif,    'ATM',   &
                       ex_albedo_nir_dif   ,  xmap_sfc, complete=.false.)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%rough_mom, 'ATM', ex_rough_mom,  xmap_sfc, complete=.false.)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%land_frac, 'ATM', ex_land_frac,  xmap_sfc, complete=.false.)

  call get_from_xgrid (Land_Ice_Atmos_Boundary%u_flux,    'ATM', ex_flux_u,     xmap_sfc, complete=.false.)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%v_flux,    'ATM', ex_flux_v,     xmap_sfc, complete=.false.)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%dtaudu,    'ATM', ex_dtaudu_atm, xmap_sfc, complete=.false.)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%dtaudv,    'ATM', ex_dtaudv_atm, xmap_sfc, complete=.false.)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%u_star,    'ATM', ex_u_star    , xmap_sfc, complete=.false.)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%b_star,    'ATM', ex_b_star    , xmap_sfc, complete=.false.)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%q_star,    'ATM', ex_q_star    , xmap_sfc, complete=.true.)

  if (do_forecast) then
     call get_from_xgrid (Ice%t_surf, 'OCN', ex_t_surf,  xmap_sfc)
  end if

  Land_Ice_Atmos_Boundary%t = Land_Ice_Atmos_Boundary%t ** 0.25
!Balaji: data_override calls moved here from coupler_main
  call data_override('ATM', 't',         Land_Ice_Atmos_Boundary%t,         Time)
  call data_override('ATM', 'albedo',    Land_Ice_Atmos_Boundary%albedo,    Time)

  call data_override('ATM', 'albedo_vis_dir',    Land_Ice_Atmos_Boundary%albedo_vis_dir,    Time)
  call data_override('ATM', 'albedo_nir_dir',    Land_Ice_Atmos_Boundary%albedo_nir_dir,    Time)
  call data_override('ATM', 'albedo_vis_dif',    Land_Ice_Atmos_Boundary%albedo_vis_dif,    Time)
  call data_override('ATM', 'albedo_nir_dif',    Land_Ice_Atmos_Boundary%albedo_nir_dif,    Time)
  call data_override('ATM', 'land_frac', Land_Ice_Atmos_Boundary%land_frac, Time)
  call data_override('ATM', 'dt_t',      Land_Ice_Atmos_Boundary%dt_t,      Time)
  do tr=1,n_atm_tr
     call get_tracer_names(MODEL_ATMOS, tr, tr_name)
     call data_override('ATM', 'dt_'//trim(tr_name), Land_Ice_Atmos_Boundary%dt_tr(:,:,tr), Time)
  enddo
  call data_override('ATM', 'u_flux',    Land_Ice_Atmos_Boundary%u_flux,    Time)
  call data_override('ATM', 'v_flux',    Land_Ice_Atmos_Boundary%v_flux,    Time)
  call data_override('ATM', 'dtaudu',    Land_Ice_Atmos_Boundary%dtaudu,    Time)
  call data_override('ATM', 'dtaudv',    Land_Ice_Atmos_Boundary%dtaudv,    Time)
  call data_override('ATM', 'u_star',    Land_Ice_Atmos_Boundary%u_star,    Time)
  call data_override('ATM', 'b_star',    Land_Ice_Atmos_Boundary%b_star,    Time)
! call data_override('ATM', 'q_star',    Land_Ice_Atmos_Boundary%q_star,    Time)
  call data_override('ATM', 'rough_mom', Land_Ice_Atmos_Boundary%rough_mom, Time)

  ! [6.3] save atmos albedo fix and old albedo (for downward SW flux calculations)
  ! on exchange grid
  ! allocate ( ex_old_albedo(n_xgrid_sfc)  )
  ! ex_old_albedo = ex_albedo
 
!!  STILL NEEDED   ????
!! IS THIS CORRECT ??
  allocate ( ex_albedo_fix(n_xgrid_sfc) )
  allocate ( ex_albedo_vis_dir_fix(n_xgrid_sfc) )
  allocate ( ex_albedo_nir_dir_fix(n_xgrid_sfc) )
  allocate ( ex_albedo_vis_dif_fix(n_xgrid_sfc) )
  allocate ( ex_albedo_nir_dif_fix(n_xgrid_sfc) )

  ex_albedo_fix = 0.
  ex_albedo_vis_dir_fix = 0.
  ex_albedo_nir_dir_fix = 0.
  ex_albedo_vis_dif_fix = 0.
  ex_albedo_nir_dif_fix = 0.


  call put_to_xgrid (Land_Ice_Atmos_Boundary%albedo, 'ATM',  ex_albedo_fix, xmap_sfc, complete=.false.)
  call put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dir, 'ATM',  &
           ex_albedo_vis_dir_fix, xmap_sfc, complete=.false.)
  call put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dir, 'ATM', &
           ex_albedo_nir_dir_fix, xmap_sfc, complete=.false.)
  call put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dif, 'ATM',   &
           ex_albedo_vis_dif_fix, xmap_sfc, complete=.false.)
  call put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dif, 'ATM',  &
           ex_albedo_nir_dif_fix, xmap_sfc, complete=.true.)
  ex_albedo_fix = (1.0-ex_albedo) / (1.0-ex_albedo_fix)
  ex_albedo_vis_dir_fix = (1.0-ex_albedo_vis_dir) / (1.0-ex_albedo_vis_dir_fix)
  ex_albedo_nir_dir_fix = (1.0-ex_albedo_nir_dir) / (1.0-ex_albedo_nir_dir_fix)
  ex_albedo_vis_dif_fix = (1.0-ex_albedo_vis_dif) / (1.0-ex_albedo_vis_dif_fix)
  ex_albedo_nir_dif_fix = (1.0-ex_albedo_nir_dif) / (1.0-ex_albedo_nir_dif_fix)

#ifdef SCM
  if (do_specified_albedo .and. do_specified_land) then
       ex_albedo_fix = 1.
       ex_albedo_vis_dir_fix = 1.
       ex_albedo_vis_dif_fix = 1.
       ex_albedo_nir_dir_fix = 1.
       ex_albedo_nir_dif_fix = 1.
  endif
#endif

! place the wind value onto the ice data type.  mac ! Check  RASF
  call get_from_xgrid (Ice%wnd, 'OCN', ex_u10, xmap_sfc)
  !=======================================================================
  ! [7] diagnostics section

  !------- save static fields first time only ------
  if (first_static) then

     !------- land fraction ------
     if ( id_land_mask > 0 ) then
        used = send_data ( id_land_mask, Land_Ice_Atmos_Boundary%land_frac, Time )
     endif

     first_static = .false.
  endif

  !------- Atm fields -----------
  do n = 1, Atm%fields%num_bcs  !{
    do m = 1, Atm%fields%bc(n)%num_fields  !{
      if ( Atm%fields%bc(n)%field(m)%id_diag > 0 ) then  !{
        if (atm%fields%bc(n)%use_10m_wind_speed .and. m .eq. ind_u10 .and. .not. Atm%fields%bc(n)%field(m)%override) then  !{
          call get_from_xgrid (Atm%fields%bc(n)%field(m)%values, 'ATM',     &
               ex_gas_fields_atm%bc(n)%field(m)%values, xmap_sfc)
        endif  !}
        if ( Atm%fields%bc(n)%field(m)%id_diag > 0 ) then  !{
           used = send_data(Atm%fields%bc(n)%field(m)%id_diag, Atm%fields%bc(n)%field(m)%values, Time )
        endif  !}
      endif  !}
    enddo  !} m
  enddo  !} n

  !------- drag coeff moisture -----------
  if ( id_wind > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_wind, xmap_sfc)
     used = send_data ( id_wind, diag_atm, Time )
  endif
  !------- drag coeff moisture -----------
  if ( id_drag_moist > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_cd_q, xmap_sfc)
     used = send_data ( id_drag_moist, diag_atm, Time )
  endif

  !------- drag coeff heat -----------
  if ( id_drag_heat > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_cd_t, xmap_sfc)
     used = send_data ( id_drag_heat, diag_atm, Time )
  endif
  
  !------- drag coeff momemtum -----------
  if ( id_drag_mom > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_cd_m, xmap_sfc)
     used = send_data ( id_drag_mom, diag_atm, Time )
  endif
  
  !------- roughness moisture -----------
  if ( id_rough_moist > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_rough_moist, xmap_sfc)
     used = send_data ( id_rough_moist, diag_atm, Time )
  endif
  
  !------- roughness heat -----------
  if ( id_rough_heat > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_rough_heat, xmap_sfc)
     used = send_data ( id_rough_heat, diag_atm, Time )
  endif
  
  !------- roughness momemtum -----------
  used = send_data ( id_rough_mom, Land_Ice_Atmos_Boundary%rough_mom, Time )
  
  !------- friction velocity -----------
  used = send_data ( id_u_star, Land_Ice_Atmos_Boundary%u_star, Time )
  
  !------- bouyancy -----------
  used = send_data ( id_b_star, Land_Ice_Atmos_Boundary%b_star, Time )

  !------- moisture scale -----------
  used = send_data ( id_q_star, Land_Ice_Atmos_Boundary%q_star, Time )

  !-----------------------------------------------------------------------
  !------ diagnostics for fields at bottom atmospheric level ------
  
  if ( id_t_atm > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_t_atm, xmap_sfc)
     used = send_data ( id_t_atm, diag_atm, Time )
  endif
  
  if ( id_u_atm > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_u_atm, xmap_sfc)
     used = send_data ( id_u_atm, diag_atm, Time )
  endif
  
  if ( id_v_atm > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_v_atm, xmap_sfc)
     used = send_data ( id_v_atm, diag_atm, Time )
  endif

  do tr = 1,n_exch_tr
     call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
     if ( id_tr_atm(tr) > 0 ) then
        call get_from_xgrid (diag_atm, 'ATM', ex_tr_atm(:,tr), xmap_sfc)
        used = send_data ( id_tr_atm(tr), diag_atm, Time )
     endif
!!jgj: add dryvmr co2_atm
! - slm Mar 25 2010: moved to resolve interdependence of diagnostic fields
     if ( id_co2_atm_dvmr > 0 .and. lowercase(trim(tr_name))=='co2') then
        ex_co2_atm_dvmr = (ex_tr_atm(:,tr) / (1.0 - ex_tr_atm(:,isphum))) * WTMAIR/WTMCO2
        call get_from_xgrid (diag_atm, 'ATM', ex_co2_atm_dvmr, xmap_sfc)
        used = send_data ( id_co2_atm_dvmr, diag_atm, Time )
     endif
  enddo

  ! - slm, Mar 25, 2002
  if ( id_p_atm > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_p_atm, xmap_sfc)
     used = send_data ( id_p_atm, diag_atm, Time )
  endif
  if ( id_z_atm > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_z_atm, xmap_sfc)
     used = send_data ( id_z_atm, diag_atm, Time )
  endif
  if ( id_gust > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_gust, xmap_sfc)
     used = send_data ( id_gust, diag_atm, Time )
  endif

  ! - bw, Sep 17, 2007
  if ( id_slp > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_slp, xmap_sfc)
     used = send_data ( id_slp, diag_atm, Time )
  endif

  !-----------------------------------------------------------------------
  !--------- diagnostics for fields at reference level ---------
  
  if ( id_t_ref > 0 .or. id_rh_ref > 0 .or. &
       id_u_ref > 0 .or. id_v_ref  > 0 .or. id_wind_ref > 0 .or. &
       id_q_ref > 0 .or. id_q_ref_land > 0 .or. &
       id_t_ref_land > 0 .or. id_rh_ref_land > 0 .or. &
       id_rh_ref_cmip >0 .or. &
       id_u_ref_land > 0 .or. id_v_ref_land  > 0 ) then
     
     zrefm = z_ref_mom
     zrefh = z_ref_heat
     !      ---- optimize calculation ----
     if ( id_t_ref <= 0 ) zrefh = zrefm
     
     call mo_profile ( zrefm, zrefh, ex_z_atm,   ex_rough_mom, &
          ex_rough_heat, ex_rough_moist,          &
          ex_u_star, ex_b_star, ex_q_star,        &
          ex_del_m, ex_del_h, ex_del_q, ex_avail  )

     !    ------- reference relative humidity -----------
     if ( id_rh_ref > 0 .or. id_rh_ref_land > 0 .or. &
          id_rh_ref_cmip > 0 .or. &
          id_q_ref > 0 .or. id_q_ref_land >0 ) then
        ex_ref = 1.0e-06
        where (ex_avail) &
           ex_ref   = ex_tr_surf(:,isphum) + (ex_tr_atm(:,isphum)-ex_tr_surf(:,isphum)) * ex_del_q
        if(id_q_ref > 0) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data(id_q_ref,diag_atm,Time)
        endif
        if(id_q_ref_land > 0) then
           call get_from_xgrid (diag_land, 'LND', ex_ref, xmap_sfc)
           used = send_tile_averaged_data(id_q_ref_land, diag_land, &
                Land%tile_size, Time, mask=Land%mask)
        endif
        ex_t_ref = 200.
        where (ex_avail) &
           ex_t_ref = ex_t_ca + (ex_t_atm-ex_t_ca) * ex_del_h
        call compute_qs (ex_t_ref, ex_p_surf, ex_qs_ref, q = ex_ref)
        call compute_qs (ex_t_ref, ex_p_surf, ex_qs_ref_cmip,  &
                         q = ex_ref, es_over_liq_and_ice = .true.)
        where (ex_avail) 
! remove cap on relative humidity -- this mod requested by cjg, ljd
!RSH       ex_ref    = MIN(100.,100.*ex_ref/ex_qs_ref)
           ex_ref2   = 100.*ex_ref/ex_qs_ref_cmip
           ex_ref    = 100.*ex_ref/ex_qs_ref
        endwhere

        if ( id_rh_ref_land > 0 ) then
           call get_from_xgrid (diag_land,'LND', ex_ref, xmap_sfc)
           used = send_tile_averaged_data ( id_rh_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if(id_rh_ref > 0) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_rh_ref, diag_atm, Time )
        endif
        if(id_rh_ref_cmip > 0) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref2, xmap_sfc)
           used = send_data ( id_rh_ref_cmip, diag_atm, Time )
        endif
     endif

     !    ------- reference temp -----------
     if ( id_t_ref > 0 .or. id_t_ref_land > 0 ) then
        where (ex_avail) &
           ex_ref = ex_t_ca + (ex_t_atm-ex_t_ca) * ex_del_h
        if (id_t_ref_land > 0) then
           call get_from_xgrid (diag_land, 'LND', ex_ref, xmap_sfc)
           used = send_tile_averaged_data ( id_t_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if ( id_t_ref > 0 ) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_t_ref, diag_atm, Time )
        endif
     endif

     !    ------- reference u comp -----------
     if ( id_u_ref > 0 .or. id_u_ref_land > 0) then
        where (ex_avail) &
           ex_ref = ex_u_surf + (ex_u_atm-ex_u_surf) * ex_del_m
        if ( id_u_ref_land > 0 ) then
           call get_from_xgrid ( diag_land, 'LND', ex_ref, xmap_sfc )
           used = send_tile_averaged_data ( id_u_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if ( id_u_ref > 0 ) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_u_ref, diag_atm, Time )
        endif
     endif

     !    ------- reference v comp -----------
     if ( id_v_ref > 0 .or. id_v_ref_land > 0 ) then
        where (ex_avail) &
           ex_ref = ex_v_surf + (ex_v_atm-ex_v_surf) * ex_del_m
        if ( id_v_ref_land > 0 ) then
           call get_from_xgrid ( diag_land, 'LND', ex_ref, xmap_sfc )
           used = send_tile_averaged_data ( id_v_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if ( id_v_ref > 0 ) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_v_ref, diag_atm, Time )
        endif
     endif

     !    ------- reference-level absolute wind -----------
     if ( id_wind_ref > 0 ) then
        where (ex_avail) &
           ex_ref = sqrt((ex_u_surf + (ex_u_atm-ex_u_surf) * ex_del_m)**2 &
                        +(ex_v_surf + (ex_v_atm-ex_v_surf) * ex_del_m)**2)
        call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
        used = send_data ( id_wind_ref, diag_atm, Time )
     endif

     !    ------- interp factor for heat ------
     if ( id_del_h > 0 ) then
        call get_from_xgrid (diag_atm, 'ATM', ex_del_h, xmap_sfc)
        used = send_data ( id_del_h, diag_atm, Time )
     endif

     !    ------- interp factor for momentum ------
     if ( id_del_m > 0 ) then
        call get_from_xgrid (diag_atm, 'ATM', ex_del_m, xmap_sfc)
        used = send_data ( id_del_m, diag_atm, Time )
     endif

     !    ------- interp factor for moisture ------
     if ( id_del_q > 0 ) then
        call get_from_xgrid (diag_atm, 'ATM', ex_del_q, xmap_sfc)
        used = send_data ( id_del_q, diag_atm, Time )
     endif

  endif
  ! topographic roughness scale
  if(id_rough_scale>0) then
     call get_from_xgrid (diag_atm, 'ATM',&
          (log(ex_z_atm/ex_rough_mom+1)/log(ex_z_atm/ex_rough_scale+1))**2, xmap_sfc)
     used = send_data(id_rough_scale, diag_atm, Time)
  endif

!Balaji
  call mpp_clock_end(sfcClock)
  call mpp_clock_end(cplClock)

!=======================================================================

end subroutine sfc_boundary_layer
! </SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="flux_down_from_atmos">
!  <OVERVIEW>
!   Returns fluxes and derivatives corrected for the implicit treatment of atmospheric 
!   diffusive fluxes, as well as the increments in the temperature and specific humidity 
!   of the lowest atmospheric layer due to all explicit processes as well as the diffusive 
!   fluxes through the top of this layer. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!    The following elements from Atmos_boundary are used as input: 
!
!        flux_u_atm = zonal wind stress (Pa)  
!        flux_v_atm = meridional wind stress (Pa)
!
!
!    The following elements of Land_boundary are output: 
!
!       flux_t_land = sensible heat flux (W/m2)
!       flux_q_land = specific humidity flux (Kg/(m2 s)
!      flux_lw_land = net longwave flux (W/m2), uncorrected for
!                     changes in surface temperature
!      flux_sw_land = net shortwave flux (W/m2)
!         dhdt_land = derivative of sensible heat flux w.r.t.
!                     surface temperature (on land model grid)  (W/(m2 K)
!         dedt_land = derivative of specific humidity flux w.r.t.
!                     surface temperature (on land model grid)  (Kg/(m2 s K)
!         drdt_land = derivative of upward longwave flux w.r.t.
!                     surface temperature (on land model grid) (W/(m2 K)
!        lprec_land = liquid precipitation, mass for one time step
!                      (Kg/m2)
!        fprec_land = frozen precipitation, mass for one time step
!                      (Kg/m2)
!
!
!    The following elements of Ice_boundary are output: 
!
!        flux_u_ice = zonal wind stress (Pa)
!        flux_v_ice = meridional wind stress (Pa)
!        coszen_ice = cosine of the zenith angle
!
!   </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_down_from_atmos (Time, Atm, Land, Ice, &
!		Atmos_boundary, Land_boundary, Ice_boundary )
!		
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <INOUT NAME="Atm" TYPE="atmos_data_type">
!   A derived data type to specify atmosphere boundary data.
!  </INOUT>
!  <IN NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </IN>
!  <IN NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>
!  <IN NAME="Atmos_boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere, land and ice.
!  </IN>
!  <INOUT NAME="Land_boundary" TYPE="atmos_land_boundary_type">
!   A derived data type to specify properties and fluxes passed from atmosphere to land.
!  </INOUT>
!  <INOUT NAME="Ice_boundary" TYPE="atmos_ice_boundary_type">
!   A derived data type to specify properties and fluxes passed from atmosphere to ice.
!  </INOUT>
!
subroutine flux_down_from_atmos (Time, Atm, Land, Ice, &
     Atmos_boundary, Land_boundary, Ice_boundary )

  type(time_type),       intent(in) :: Time
  type(atmos_data_type), intent(inout) :: Atm
  type(land_data_type),  intent(in) :: Land
  type(ice_data_type),   intent(in) :: Ice
  type(land_ice_atmos_boundary_type),intent(in) :: Atmos_boundary
  type(atmos_land_boundary_type),    intent(inout):: Land_boundary
  type(atmos_ice_boundary_type),     intent(inout):: Ice_boundary

  real, dimension(n_xgrid_sfc) :: ex_flux_sw, ex_flux_lwd, &
       ex_flux_sw_dir,  &
                                    ex_flux_sw_dif,  &
      ex_flux_sw_down_vis_dir, ex_flux_sw_down_total_dir,  &
      ex_flux_sw_down_vis_dif, ex_flux_sw_down_total_dif,  &
       ex_flux_sw_vis, &
       ex_flux_sw_vis_dir, &
       ex_flux_sw_vis_dif, &
       ex_lprec, ex_fprec,      &
       ex_tprec, & ! temperature of precipitation, currently equal to atm T
       ex_u_star_smooth,        &
       ex_coszen

  real, dimension(n_xgrid_sfc) :: ex_gamma  , ex_dtmass,  &
       ex_delta_t, ex_delta_u, ex_delta_v, ex_dflux_t

  real, dimension(n_xgrid_sfc,n_exch_tr) :: &
       ex_delta_tr, & ! tracer tendencies
       ex_dflux_tr    ! fracer flux change

  real    :: cp_inv
  logical :: used
  logical :: ov
  integer :: ier

  character(32) :: tr_name ! name of the tracer
  integer :: tr, n, m ! tracer indices

!Balaji
  call mpp_clock_begin(cplClock)
  call mpp_clock_begin(fluxAtmDnClock)
  ov = .FALSE.
!-----------------------------------------------------------------------
!Balaji: data_override calls moved here from coupler_main            
  call data_override ('ATM', 'flux_sw',  Atm%flux_sw, Time)
  call data_override ('ATM', 'flux_sw_dir',  Atm%flux_sw_dir, Time)
  call data_override ('ATM', 'flux_sw_dif',  Atm%flux_sw_dif, Time)
  call data_override ('ATM', 'flux_sw_down_vis_dir',  Atm%flux_sw_down_vis_dir, Time)
  call data_override ('ATM', 'flux_sw_down_vis_dif',  Atm%flux_sw_down_vis_dif, Time)
  call data_override ('ATM', 'flux_sw_down_total_dir',  Atm%flux_sw_down_total_dir, Time)
  call data_override ('ATM', 'flux_sw_down_total_dif',  Atm%flux_sw_down_total_dif, Time)
  call data_override ('ATM', 'flux_sw_vis',  Atm%flux_sw_vis, Time)
  call data_override ('ATM', 'flux_sw_vis_dir',  Atm%flux_sw_vis_dir, Time)
  call data_override ('ATM', 'flux_sw_vis_dif',  Atm%flux_sw_vis_dif, Time)
  call data_override ('ATM', 'flux_lw',  Atm%flux_lw, Time)
  call data_override ('ATM', 'lprec',    Atm%lprec,   Time)
  call data_override ('ATM', 'fprec',    Atm%fprec,   Time)
  call data_override ('ATM', 'coszen',   Atm%coszen,  Time)
  call data_override ('ATM', 'dtmass',   Atm%Surf_Diff%dtmass, Time)
  call data_override ('ATM', 'delta_t',  Atm%Surf_Diff%delta_t, Time)
  call data_override ('ATM', 'dflux_t',  Atm%Surf_Diff%dflux_t, Time)
  do tr = 1,n_atm_tr
     call get_tracer_names(MODEL_ATMOS,tr,tr_name)
     call data_override ('ATM', 'delta_'//trim(tr_name),  Atm%Surf_Diff%delta_tr(:,:,tr), Time)
     call data_override ('ATM', 'dflux_'//trim(tr_name),  Atm%Surf_Diff%dflux_tr(:,:,tr), Time)
  enddo

!---- put atmosphere quantities onto exchange grid ----

  if(sw1way_bug) then
     call put_to_xgrid (Atm%flux_sw, 'ATM', ex_flux_sw, xmap_sfc, complete=.false.)
     call put_to_xgrid (Atm%flux_sw_vis, 'ATM', ex_flux_sw_vis, xmap_sfc, complete=.false.)
  end if
  ex_flux_sw_dir     = 0.0
  ex_flux_sw_vis_dir = 0.0
  ex_flux_sw_dif     = 0.0
  ex_flux_sw_vis_dif = 0.0
  ex_flux_lwd        = 0.0                           
  call put_to_xgrid (Atm%flux_sw_dir, 'ATM', ex_flux_sw_dir, xmap_sfc, complete=.false.)
  call put_to_xgrid (Atm%flux_sw_vis_dir, 'ATM', ex_flux_sw_vis_dir, xmap_sfc, complete=.false.)
  call put_to_xgrid (Atm%flux_sw_dif, 'ATM', ex_flux_sw_dif, xmap_sfc, complete=.false.)
  call put_to_xgrid (Atm%flux_sw_vis_dif, 'ATM', ex_flux_sw_vis_dif, xmap_sfc, complete=.false.)
  call put_to_xgrid (Atm%flux_sw_down_vis_dir, 'ATM', ex_flux_sw_down_vis_dir, xmap_sfc, complete=.false.)
  call put_to_xgrid (Atm%flux_sw_down_total_dir, 'ATM', ex_flux_sw_down_total_dir, xmap_sfc, complete=.false.)
  call put_to_xgrid (Atm%flux_sw_down_vis_dif, 'ATM', ex_flux_sw_down_vis_dif, xmap_sfc, complete=.false.)
  call put_to_xgrid (Atm%flux_sw_down_total_dif, 'ATM', ex_flux_sw_down_total_dif, xmap_sfc, complete=.false.)

  !  ccc = conservation_check(Atm%lprec, 'ATM', xmap_sfc)
  !  if (mpp_pe()== mpp_root_pe()) print *,'LPREC', ccc

!!$  if(do_area_weighted_flux) then
!!$     call put_to_xgrid (Atm%lprec * AREA_ATM_MODEL,   'ATM', ex_lprec, xmap_sfc)
!!$     call put_to_xgrid (Atm%fprec * AREA_ATM_MODEL,   'ATM', ex_fprec, xmap_sfc)
!!$  else
     call put_to_xgrid (Atm%lprec,   'ATM', ex_lprec, xmap_sfc, complete=.false.)
     call put_to_xgrid (Atm%fprec,   'ATM', ex_fprec, xmap_sfc, complete=.false.)
     call put_to_xgrid (Atm%t_bot,   'ATM', ex_tprec, xmap_sfc, complete=.false.)
!!$  endif

  call put_to_xgrid (Atm%coszen,  'ATM', ex_coszen, xmap_sfc, complete=.true.)

  call put_to_xgrid (Atm%flux_lw, 'ATM', ex_flux_lwd, xmap_sfc, remap_method=remap_method, complete=.false.)
  if(ex_u_star_smooth_bug) then
     call put_to_xgrid (Atmos_boundary%u_star, 'ATM', ex_u_star_smooth, xmap_sfc, remap_method=remap_method, complete=.false.)
     ex_u_star = ex_u_star_smooth
  endif


! MOD changed the following two lines to put Atmos%surf_diff%delta_u and v
! on exchange grid instead of the stresses themselves so that only the 
! implicit corrections are filtered through the atmospheric grid not the
! stresses themselves
  ex_delta_u = 0.0; ex_delta_v = 0.0
  call put_to_xgrid (Atm%Surf_Diff%delta_u, 'ATM', ex_delta_u, xmap_sfc, remap_method=remap_method, complete=.false.)
  call put_to_xgrid (Atm%Surf_Diff%delta_v, 'ATM', ex_delta_v, xmap_sfc, remap_method=remap_method, complete=.true.)

  ! MOD update stresses using atmos delta's but derivatives on exchange grid
  ex_flux_u = ex_flux_u + ex_delta_u*ex_dtaudu_atm
  ex_flux_v = ex_flux_v + ex_delta_v*ex_dtaudv_atm

!-----------------------------------------------------------------------
!---- adjust sw flux for albedo variations on exch grid ----
!---- adjust 4 categories (vis/nir dir/dif) separately  ----
  if( sw1way_bug ) then ! to reproduce old results, may remove in the next major release.
!-----------------------------------------------------------------------
!---- adjust sw flux for albedo variations on exch grid ----

     ex_flux_sw = ex_flux_sw * ex_albedo_fix


     ex_flux_sw_vis = ex_flux_sw_vis * ex_albedo_vis_dir_fix
     ex_flux_sw_dir = ex_flux_sw_dir * ex_albedo_vis_dir_fix
     ex_flux_sw_dif = ex_flux_sw_dif * ex_albedo_vis_dif_fix
     ex_flux_sw_vis_dir = ex_flux_sw_vis_dir * ex_albedo_vis_dir_fix
     ex_flux_sw_vis_dif = ex_flux_sw_vis_dif * ex_albedo_vis_dif_fix
  else 
     ex_flux_sw_dir = ex_flux_sw_dir - ex_flux_sw_vis_dir     ! temporarily nir/dir
     ex_flux_sw_dir = ex_flux_sw_dir * ex_albedo_nir_dir_fix  ! fix nir/dir
     ex_flux_sw_vis_dir = ex_flux_sw_vis_dir * ex_albedo_vis_dir_fix ! fix vis/dir
     ex_flux_sw_dir = ex_flux_sw_dir + ex_flux_sw_vis_dir     ! back to total dir

     ex_flux_sw_dif = ex_flux_sw_dif - ex_flux_sw_vis_dif     ! temporarily nir/dif
     ex_flux_sw_dif = ex_flux_sw_dif * ex_albedo_nir_dif_fix  ! fix nir/dif
     ex_flux_sw_vis_dif = ex_flux_sw_vis_dif * ex_albedo_vis_dif_fix ! fix vis/dif
     ex_flux_sw_dif = ex_flux_sw_dif + ex_flux_sw_vis_dif     ! back to total dif

     ex_flux_sw_vis = ex_flux_sw_vis_dir + ex_flux_sw_vis_dif ! legacy, remove later
     ex_flux_sw     = ex_flux_sw_dir     + ex_flux_sw_dif     ! legacy, remove later
  end if

!!$  ex_flux_sw_dir = ex_flux_sw_dir - ex_flux_sw_vis_dir            ! temporarily nir/dir
!!$  ex_flux_sw_dir = ex_flux_sw_dir * ex_albedo_nir_dir_fix         ! fix nir/dir
!!$  ex_flux_sw_vis_dir = ex_flux_sw_vis_dir * ex_albedo_vis_dir_fix ! fix vis/dir
!!$  ex_flux_sw_dir = ex_flux_sw_dir + ex_flux_sw_vis_dir            ! back to total dir
!!$
!!$  ex_flux_sw_dif = ex_flux_sw_dif - ex_flux_sw_vis_dif            ! temporarily nir/dif
!!$  ex_flux_sw_dif = ex_flux_sw_dif * ex_albedo_nir_dif_fix         ! fix nir/dif
!!$  ex_flux_sw_vis_dif = ex_flux_sw_vis_dif * ex_albedo_vis_dif_fix ! fix vis/dif
!!$  ex_flux_sw_dif = ex_flux_sw_dif + ex_flux_sw_vis_dif            ! back to total dif
!!$
!!$  ex_flux_sw_vis = ex_flux_sw_vis_dir + ex_flux_sw_vis_dif        ! legacy, remove later
!!$  ex_flux_sw     = ex_flux_sw_dir     + ex_flux_sw_dif            ! legacy, remove later

  deallocate ( ex_albedo_fix )
  deallocate ( ex_albedo_vis_dir_fix )
  deallocate ( ex_albedo_nir_dir_fix )
  deallocate ( ex_albedo_vis_dif_fix )
  deallocate ( ex_albedo_nir_dif_fix )
!----- compute net longwave flux (down-up) -----
  ! (note: lw up already in ex_flux_lw)

  ex_flux_lw = ex_flux_lwd - ex_flux_lw

!-----------------------------------------------------------------------
!----- adjust fluxes for implicit dependence on atmosphere ----

  do tr = 1,n_exch_tr
     n = tr_table(tr)%atm
     call put_to_xgrid (Atm%Surf_Diff%delta_tr(:,:,n), 'ATM', ex_delta_tr(:,tr), xmap_sfc, complete=.false. )
     call put_to_xgrid (Atm%Surf_Diff%dflux_tr(:,:,n), 'ATM', ex_dflux_tr(:,tr), xmap_sfc, complete=.false. )
  enddo

  call put_to_xgrid (Atm%Surf_Diff%dtmass , 'ATM', ex_dtmass , xmap_sfc, complete=.false. )
  call put_to_xgrid (Atm%Surf_Diff%delta_t, 'ATM', ex_delta_t, xmap_sfc, complete=.false. )
  call put_to_xgrid (Atm%Surf_Diff%dflux_t, 'ATM', ex_dflux_t, xmap_sfc, complete=.true. )

  cp_inv = 1.0/cp_air

  where(ex_avail)

     ! temperature

     ex_gamma      =  1./ (1.0 - ex_dtmass*(ex_dflux_t + ex_dhdt_atm*cp_inv))
     ex_e_t_n      =  ex_dtmass*ex_dhdt_surf*cp_inv*ex_gamma
     ex_f_t_delt_n = (ex_delta_t + ex_dtmass * ex_flux_t*cp_inv) * ex_gamma    
     
     ex_flux_t     =  ex_flux_t        + ex_dhdt_atm * ex_f_t_delt_n 
     ex_dhdt_surf  =  ex_dhdt_surf     + ex_dhdt_atm * ex_e_t_n   

     ! moisture
!     ex_gamma      =  1./ (1.0 - ex_dtmass*(ex_dflux_q + ex_dedq_atm))
! here it looks like two derivatives with different units are added together,
! but in fact they are not: ex_dedt_surf and ex_dedq_surf defined in complimentary
! regions of exchange grid, so that if one of them is not zero the other is, and
! vice versa.
!     ex_e_q_n      =  ex_dtmass*(ex_dedt_surf+ex_dedq_surf) * ex_gamma
!     ex_f_q_delt_n = (ex_delta_q  + ex_dtmass * ex_flux_q) * ex_gamma    
!     ex_flux_q     =  ex_flux_q    + ex_dedq_atm * ex_f_q_delt_n 
!     ex_dedt_surf  =  ex_dedt_surf + ex_dedq_atm * ex_e_q_n
!     ex_dedq_surf  =  ex_dedq_surf + ex_dedq_atm * ex_e_q_n
     ! moisture vs. surface temperture, assuming saturation
     ex_gamma   =  1.0 / (1.0 - ex_dtmass*(ex_dflux_tr(:,isphum) + ex_dfdtr_atm(:,isphum)))
     ex_e_q_n      =  ex_dtmass * ex_dedt_surf * ex_gamma
     ex_dedt_surf  =  ex_dedt_surf + ex_dfdtr_atm(:,isphum) * ex_e_q_n
  endwhere
  do tr = 1,n_exch_tr
     where(ex_avail)
        ex_gamma   =  1.0 / (1.0 - ex_dtmass*(ex_dflux_tr(:,tr) + ex_dfdtr_atm(:,tr)))

        ex_e_tr_n(:,tr)      =  ex_dtmass*ex_dfdtr_surf(:,tr)*ex_gamma
        ex_f_tr_delt_n(:,tr) = (ex_delta_tr(:,tr)+ex_dtmass*ex_flux_tr(:,tr))*ex_gamma    
     
        ex_flux_tr(:,tr)     =  ex_flux_tr(:,tr) + ex_dfdtr_atm(:,tr)*ex_f_tr_delt_n(:,tr) 
        ex_dfdtr_surf(:,tr)  =  ex_dfdtr_surf(:,tr) + ex_dfdtr_atm(:,tr)*ex_e_tr_n(:,tr)
     endwhere
  enddo
!-----------------------------------------------------------------------
!---- output fields on the land grid -------

  call get_from_xgrid (Land_boundary%t_flux,  'LND', ex_flux_t,    xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux, 'LND', ex_flux_sw,   xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux_down_vis_dir, 'LND', ex_flux_sw_down_vis_dir,   xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux_down_total_dir, 'LND', ex_flux_sw_down_total_dir,   xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux_down_vis_dif, 'LND', ex_flux_sw_down_vis_dif,   xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux_down_total_dif, 'LND', ex_flux_sw_down_total_dif,   xmap_sfc)
  call get_from_xgrid (Land_boundary%lw_flux, 'LND', ex_flux_lw,   xmap_sfc)
#ifdef SCM
  if (do_specified_land .and. do_specified_flux) then
    call get_from_xgrid (Land_boundary%dhdt,  'LND', ex_dhdt_surf_forland, xmap_sfc)
  else
    call get_from_xgrid (Land_boundary%dhdt,  'LND', ex_dhdt_surf, xmap_sfc)
  endif
#else
  call get_from_xgrid (Land_boundary%dhdt,    'LND', ex_dhdt_surf, xmap_sfc)
#endif
  call get_from_xgrid (Land_boundary%drdt,    'LND', ex_drdt_surf, xmap_sfc)
  call get_from_xgrid (Land_boundary%p_surf,  'LND', ex_p_surf,    xmap_sfc)

  call get_from_xgrid (Land_boundary%lprec,   'LND', ex_lprec,     xmap_sfc)
  call get_from_xgrid (Land_boundary%fprec,   'LND', ex_fprec,     xmap_sfc)
  call get_from_xgrid (Land_boundary%tprec,   'LND', ex_tprec,     xmap_sfc)
!!$  if(do_area_weighted_flux) then
!!$     ! evap goes here???
!!$     do k = 1, size(Land_boundary%lprec, dim=3)
!!$        ! Note: we divide by AREA_ATM_MODEL, which should be the same as
!!$        ! AREA_LND_MODEL (but the latter may not be defined)
!!$        call divide_by_area(data=Land_boundary%lprec(:,:,k), area=AREA_ATM_MODEL)
!!$        call divide_by_area(data=Land_boundary%fprec(:,:,k), area=AREA_ATM_MODEL)
!!$     enddo
!!$  endif

  if(associated(Land_boundary%drag_q)) then
     call get_from_xgrid (Land_boundary%drag_q, 'LND', ex_drag_q,    xmap_sfc)
     call data_override('LND', 'drag_q', Land_boundary%drag_q,  Time )
  endif
  if(associated(Land_boundary%lwdn_flux)) then
     call get_from_xgrid (Land_boundary%lwdn_flux, 'LND', ex_flux_lwd, xmap_sfc)
     call data_override('LND', 'lwdn_flux', Land_boundary%lwdn_flux, Time )
  endif
  if(associated(Land_boundary%cd_m)) then
     call get_from_xgrid (Land_boundary%cd_m, 'LND', ex_cd_m, xmap_sfc)
     call data_override('LND', 'cd_m', Land_boundary%cd_m, Time )
  endif
  if(associated(Land_boundary%cd_t)) then
     call get_from_xgrid (Land_boundary%cd_t, 'LND', ex_cd_t, xmap_sfc)
     call data_override('LND', 'cd_t', Land_boundary%cd_t, Time )
  endif
  if(associated(Land_boundary%bstar)) then
     call get_from_xgrid (Land_boundary%bstar, 'LND', ex_b_star, xmap_sfc)
     call data_override('LND', 'bstar',  Land_boundary%bstar, Time )
  endif
  if(associated(Land_boundary%ustar)) then
     call get_from_xgrid (Land_boundary%ustar, 'LND', ex_u_star, xmap_sfc)
     call data_override('LND', 'ustar',  Land_boundary%ustar, Time )
  endif
  if(associated(Land_boundary%wind)) then
     call get_from_xgrid (Land_boundary%wind, 'LND', ex_wind, xmap_sfc)
     call data_override('LND', 'wind',  Land_boundary%wind, Time )
  endif
  if(associated(Land_boundary%z_bot)) then
     call get_from_xgrid (Land_boundary%z_bot, 'LND', ex_z_atm, xmap_sfc)
     call data_override('LND', 'z_bot',  Land_boundary%z_bot, Time )
  endif

  Land_boundary%tr_flux(:,:,:,:) = 0.0
  Land_boundary%dfdtr(:,:,:,:) = 0.0
  do tr = 1,n_exch_tr
     n = tr_table(tr)%lnd
     if(n /= NO_TRACER ) then
        call get_from_xgrid (Land_boundary%tr_flux(:,:,:,n), 'LND', ex_flux_tr(:,tr), xmap_sfc)
        call get_from_xgrid (Land_boundary%dfdtr(:,:,:,n),   'LND', ex_dfdtr_surf(:,tr), xmap_sfc)
#ifdef SCM
        if (do_specified_land .and. do_specified_flux .and. tr.eq.isphum) then
          call get_from_xgrid (Land_boundary%dfdtr(:,:,:,n),   'LND', ex_dedq_surf_forland(:), xmap_sfc)
        endif
#endif
     endif
  enddo

!  current time is Time: is that ok? not available in land_data_type
!Balaji: data_override calls moved here from coupler_main
  call data_override('LND', 't_flux',  Land_boundary%t_flux,  Time )
  call data_override('LND', 'lw_flux', Land_boundary%lw_flux, Time )
  call data_override('LND', 'sw_flux', Land_boundary%sw_flux, Time )
  call data_override('LND', 'sw_flux_down_vis_dir', Land_boundary%sw_flux_down_vis_dir, Time )
  call data_override('LND', 'sw_flux_down_total_dir', Land_boundary%sw_flux_down_total_dir, Time )
  call data_override('LND', 'sw_flux_down_vis_dif', Land_boundary%sw_flux_down_vis_dif, Time )
  call data_override('LND', 'sw_flux_down_total_dif', Land_boundary%sw_flux_down_total_dif, Time )
  
  call data_override('LND', 'lprec',   Land_boundary%lprec,   Time )
  call data_override('LND', 'fprec',   Land_boundary%fprec,   Time )
  call data_override('LND', 'dhdt',    Land_boundary%dhdt,    Time )
  call data_override('LND', 'drdt',    Land_boundary%drdt,    Time )
  call data_override('LND', 'p_surf',  Land_boundary%p_surf,  Time )
  do tr = 1,n_lnd_tr
     call get_tracer_names(MODEL_LAND, tr, tr_name)
     call data_override('LND', trim(tr_name)//'_flux', Land_boundary%tr_flux(:,:,:,tr), Time)
     call data_override('LND', 'dfd'//trim(tr_name),   Land_boundary%dfdtr  (:,:,:,tr), Time)
  enddo

!-----------------------------------------------------------------------
!---- output fields on the ice grid -------

  call get_from_xgrid (Ice_boundary%t_flux,   'OCN', ex_flux_t,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%q_flux,   'OCN', ex_flux_tr(:,isphum), xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux_vis_dir,  'OCN', ex_flux_sw_vis_dir,   xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux_nir_dir,  'OCN', ex_flux_sw_dir,xmap_sfc)
  Ice_boundary%sw_flux_nir_dir = Ice_boundary%sw_flux_nir_dir - Ice_boundary%sw_flux_vis_dir ! ice & ocean use these 4: dir/dif nir/vis

  call get_from_xgrid (Ice_boundary%sw_flux_vis_dif,  'OCN', ex_flux_sw_vis_dif,   xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux_nir_dif,  'OCN', ex_flux_sw_dif,xmap_sfc)
  Ice_boundary%sw_flux_nir_dif = Ice_boundary%sw_flux_nir_dif - Ice_boundary%sw_flux_vis_dif ! ice & ocean use these 4: dir/dif nir/vis

  call get_from_xgrid (Ice_boundary%lw_flux,  'OCN', ex_flux_lw,   xmap_sfc)
  call get_from_xgrid (Ice_boundary%dhdt,     'OCN', ex_dhdt_surf, xmap_sfc)
  call get_from_xgrid (Ice_boundary%dedt,     'OCN', ex_dedt_surf, xmap_sfc)
  call get_from_xgrid (Ice_boundary%drdt,     'OCN', ex_drdt_surf, xmap_sfc)
  call get_from_xgrid (Ice_boundary%u_flux,   'OCN', ex_flux_u,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%v_flux,   'OCN', ex_flux_v,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%u_star,   'OCN', ex_u_star,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%coszen,   'OCN', ex_coszen,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%p,        'OCN', ex_slp,       xmap_sfc) ! mw mod

  call get_from_xgrid (Ice_boundary%lprec,    'OCN', ex_lprec,     xmap_sfc)
  call get_from_xgrid (Ice_boundary%fprec,    'OCN', ex_fprec,     xmap_sfc)
!!$  if (do_area_weighted_flux) then
!!$     where (AREA_ATM_SPHERE /= 0)
!!$        Ice_boundary%lprec = Ice_boundary%lprec * AREA_ATM_MODEL/AREA_ATM_SPHERE
!!$        Ice_boundary%fprec = Ice_boundary%fprec * AREA_ATM_MODEL/AREA_ATM_SPHERE
!!$     end where
!!$  endif
!!$  if(do_area_weighted_flux) then
!!$     do k = 1, size(Ice_boundary%lprec, dim=3)
!!$        call divide_by_area(data=Ice_boundary%lprec(:,:,k), area=AREA_ATM_SPHERE)
!!$        call divide_by_area(data=Ice_boundary%fprec(:,:,k), area=AREA_ATM_SPHERE)
!!$     enddo
!!$  endif

! Extra fluxes
  do n = 1, Ice_boundary%fluxes%num_bcs  !{
    do m = 1, Ice_boundary%fluxes%bc(n)%num_fields  !{
      call get_from_xgrid (Ice_boundary%fluxes%bc(n)%field(m)%values, 'OCN',  &
           ex_gas_fluxes%bc(n)%field(m)%values, xmap_sfc)
    enddo  !} m
  enddo  !} n

!Balaji: data_override calls moved here from coupler_main
  call data_override('ICE', 'u_flux', Ice_boundary%u_flux,  Time)
  call data_override('ICE', 'v_flux', Ice_boundary%v_flux,  Time)
  call data_override('ICE', 't_flux', Ice_boundary%t_flux,  Time)
  call data_override('ICE', 'q_flux', Ice_boundary%q_flux,  Time)
  call data_override('ICE', 'lw_flux',Ice_boundary%lw_flux, Time)
  call data_override('ICE', 'lw_flux_dn',Ice_boundary%lw_flux, Time, override=ov)
  if (ov) then
    Ice_boundary%lw_flux = Ice_boundary%lw_flux - stefan*Ice%t_surf**4
  endif
  call data_override('ICE', 'sw_flux_nir_dir',Ice_boundary%sw_flux_nir_dir, Time)
  call data_override('ICE', 'sw_flux_vis_dir',Ice_boundary%sw_flux_vis_dir, Time)
  call data_override('ICE', 'sw_flux_nir_dif',Ice_boundary%sw_flux_nir_dif, Time, override=ov)
  call data_override('ICE', 'sw_flux_vis_dif',Ice_boundary%sw_flux_vis_dif, Time)
  call data_override('ICE', 'sw_flux_vis_dir_dn',Ice_boundary%sw_flux_vis_dir, Time, override=ov)
  if (ov) then
    Ice_boundary%sw_flux_vis_dir = Ice_boundary%sw_flux_vis_dir*(1-Ice%albedo_vis_dir)
  endif
  call data_override('ICE', 'sw_flux_vis_dif_dn',Ice_boundary%sw_flux_vis_dif, Time, override=ov)
  if (ov) then
    Ice_boundary%sw_flux_vis_dif = Ice_boundary%sw_flux_vis_dif*(1-Ice%albedo_vis_dif)
  endif
  call data_override('ICE', 'sw_flux_nir_dir_dn',Ice_boundary%sw_flux_nir_dir, Time, override=ov)
  if (ov) then
    Ice_boundary%sw_flux_nir_dir = Ice_boundary%sw_flux_nir_dir*(1-Ice%albedo_nir_dir)
  endif
  call data_override('ICE', 'sw_flux_nir_dif_dn',Ice_boundary%sw_flux_nir_dif, Time, override=ov)
  if (ov) then
    Ice_boundary%sw_flux_nir_dif = Ice_boundary%sw_flux_nir_dif*(1-Ice%albedo_nir_dif)
  endif
  call data_override('ICE', 'lprec',  Ice_boundary%lprec,   Time)
  call data_override('ICE', 'fprec',  Ice_boundary%fprec,   Time)
  call data_override('ICE', 'dhdt',   Ice_boundary%dhdt,    Time)
  call data_override('ICE', 'dedt',   Ice_boundary%dedt,    Time)
  call data_override('ICE', 'drdt',   Ice_boundary%drdt,    Time)
  call data_override('ICE', 'coszen', Ice_boundary%coszen,  Time)
  call data_override('ICE', 'p',      Ice_boundary%p,       Time)

  do n = 1, Ice_boundary%fluxes%num_bcs  !{
    do m = 1, Ice_boundary%fluxes%bc(n)%num_fields  !{
      call data_override('ICE', Ice_boundary%fluxes%bc(n)%field(m)%name,     &
           Ice_boundary%fluxes%bc(n)%field(m)%values, Time)
      if ( Ice_boundary%fluxes%bc(n)%field(m)%id_diag > 0 ) then  !{
        used = send_data(Ice_boundary%fluxes%bc(n)%field(m)%id_diag, Ice_boundary%fluxes%bc(n)%field(m)%values, Time )
      endif  !}
    enddo  !} m
  enddo  !} n

  ! compute stock changes

  ! Atm -> Lnd (precip)
  call stock_move( &
       & FROM = Atm_stock(ISTOCK_WATER),  &
       & TO   = Lnd_stock(ISTOCK_WATER), &
       & DATA = (Land_boundary%lprec + Land_boundary%fprec), &
       & grid_index=X1_GRID_LND, &
       & xmap=xmap_sfc, &
       & delta_t=Dt_atm, &
       & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
       & radius=Radius, ier=ier, verbose='stock move PRECIP (Atm->Lnd) ')

  ! Atm -> Lnd (heat)
  call stock_move( &
       & FROM = Atm_stock(ISTOCK_HEAT),  &
       & TO   = Lnd_stock(ISTOCK_HEAT), &
       & DATA = (-Land_boundary%t_flux + Land_boundary%lw_flux +  Land_boundary%sw_flux - Land_boundary%fprec*HLF), &
       & grid_index=X1_GRID_LND, &
       & xmap=xmap_sfc, &
       & delta_t=Dt_atm, &
       & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
       & radius=Radius, ier=ier, verbose='stock move HEAT (Atm->Lnd) ')

  ! Atm -> Ice (precip)
  call stock_move( &
       & FROM = Atm_stock(ISTOCK_WATER), &
       & TO   = Ice_stock(ISTOCK_WATER), &
       & DATA = (Ice_boundary%lprec + Ice_boundary%fprec), &
       & grid_index=X1_GRID_ICE, &
       & xmap=xmap_sfc, &
       & delta_t=Dt_atm, &
       & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
       & radius=Radius, ier=ier, verbose='stock move PRECIP (Atm->Ice) ')

  ! Atm -> Ice (heat)
  call stock_move( &
       & FROM = Atm_stock(ISTOCK_HEAT), &
       & TO   = Ice_stock(ISTOCK_HEAT), &
       & DATA = (-Ice_boundary%t_flux + Ice_boundary%lw_flux - Ice_boundary%fprec*HLF + Ice_boundary%sw_flux_vis_dir + &
                  Ice_boundary%sw_flux_vis_dif + Ice_boundary%sw_flux_nir_dir + Ice_boundary%sw_flux_nir_dif), &
       & grid_index=X1_GRID_ICE, &
       & xmap=xmap_sfc, &
       & delta_t=Dt_atm, &
       & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
       & radius=Radius, ier=ier, verbose='stock move HEAT (Atm->Ice) ')

  deallocate ( ex_flux_u, ex_flux_v, ex_dtaudu_atm, ex_dtaudv_atm)

  !=======================================================================
  !-------------------- diagnostics section ------------------------------

  !------- zonal wind stress -----------
  used = send_data ( id_u_flux, Atmos_boundary%u_flux, Time )

  !------- meridional wind stress -----------
  used = send_data ( id_v_flux, Atmos_boundary%v_flux, Time )

!Balaji
  call mpp_clock_end(fluxAtmDnClock)
  call mpp_clock_end(cplClock)
!=======================================================================

  end subroutine flux_down_from_atmos
! </SUBROUTINE>

!#######################################################################
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! flux_land_to_ice - translate runoff from land to ice grids                   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! <SUBROUTINE NAME="flux_land_to_ice">
!  <OVERVIEW>
!   Conservative transfer of water and snow discharge from the land model to sea ice/ocean model. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!    The following elements are transferred from the Land to the Land_ice_boundary: 
!
!        discharge --> runoff (kg/m2)
!        discharge_snow --> calving (kg/m2)
!
!  </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_land_to_ice(Time, Land, Ice, Land_Ice_Boundary )
!		
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </IN>
!  <IN NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>
!  <INOUT NAME="Land_Ice_Boundary" TYPE="land_ice_boundary_type">
!   A derived data type to specify properties and fluxes passed from land to ice.
!  </INOUT>
!
subroutine flux_land_to_ice( Time, Land, Ice, Land_Ice_Boundary )
  type(time_type),               intent(in) :: Time
  type(land_data_type),          intent(in) :: Land
  type(ice_data_type),           intent(in) :: Ice
!real, dimension(:,:), intent(out) :: runoff_ice, calving_ice
  type(land_ice_boundary_type),  intent(inout):: Land_Ice_Boundary
  
  integer :: ier
  real, dimension(n_xgrid_runoff) :: ex_runoff, ex_calving, ex_runoff_hflx, ex_calving_hflx
  real, dimension(size(Land_Ice_Boundary%runoff,1),size(Land_Ice_Boundary%runoff,2),1) :: ice_buf
!Balaji
  call mpp_clock_begin(cplClock)
  call mpp_clock_begin(fluxLandIceClock)

  ! ccc = conservation_check(Land%discharge, 'LND', xmap_runoff)
  ! if (mpp_pe()==mpp_root_pe()) print *,'RUNOFF', ccc

if (do_runoff) then
  call put_to_xgrid ( Land%discharge,      'LND', ex_runoff,  xmap_runoff)
  call put_to_xgrid ( Land%discharge_snow, 'LND', ex_calving, xmap_runoff)
  call put_to_xgrid ( Land%discharge_heat,      'LND', ex_runoff_hflx,  xmap_runoff)
  call put_to_xgrid ( Land%discharge_snow_heat, 'LND', ex_calving_hflx, xmap_runoff)
  call get_from_xgrid (ice_buf, 'OCN', ex_runoff,  xmap_runoff)
  Land_Ice_Boundary%runoff = ice_buf(:,:,1);
  call get_from_xgrid (ice_buf, 'OCN', ex_calving, xmap_runoff)
  Land_Ice_Boundary%calving = ice_buf(:,:,1);
  call get_from_xgrid (ice_buf, 'OCN', ex_runoff_hflx,  xmap_runoff)
  Land_Ice_Boundary%runoff_hflx = ice_buf(:,:,1);
  call get_from_xgrid (ice_buf, 'OCN', ex_calving_hflx, xmap_runoff)
  Land_Ice_Boundary%calving_hflx = ice_buf(:,:,1);
!Balaji
  call data_override('ICE', 'runoff' , Land_Ice_Boundary%runoff , Time)
  call data_override('ICE', 'calving', Land_Ice_Boundary%calving, Time)
  call data_override('ICE', 'runoff_hflx' , Land_Ice_Boundary%runoff_hflx , Time)
  call data_override('ICE', 'calving_hflx', Land_Ice_Boundary%calving_hflx, Time)

  ! compute stock increment
  ice_buf(:,:,1) = Land_Ice_Boundary%runoff + Land_Ice_Boundary%calving
  call stock_move(from=Lnd_stock(ISTOCK_WATER), to=Ice_stock(ISTOCK_WATER), &
              & grid_index=X2_GRID_ICE, &
              & data=ice_buf, &
              & xmap=xmap_runoff, &
              & delta_t=Dt_cpl, &
              & from_side=ISTOCK_SIDE, to_side=ISTOCK_SIDE, &
              & radius=Radius, ier=ier, verbose='stock move RUNOFF+CALVING (Lnd->Ice) ')
else   
   Land_Ice_Boundary%runoff = 0.0 
   Land_Ice_Boundary%calving = 0.0
   Land_Ice_Boundary%runoff_hflx = 0.0 
   Land_Ice_Boundary%calving_hflx = 0.0
endif

  call mpp_clock_end(fluxLandIceClock)
  call mpp_clock_end(cplClock)

end subroutine flux_land_to_ice
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="flux_ice_to_ocean">
!  <OVERVIEW>
!   Takes the ice model state (fluxes at the bottom of the ice) and interpolates it to the ocean model grid. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!   The following quantities are transferred from the Ice to the ice_ocean_boundary_type: 
!
!       flux_u = zonal wind stress (Pa)
!       flux_v = meridional wind stress (Pa)
!       flux_t = sensible heat flux (W/m2)
!       flux_q = specific humidity flux (Kg/m2/s)
!    flux_salt = salt flux (Kg/m2/s)
!      flux_sw = net (down-up) shortwave flux (W/m2)
!      flux_lw = net (down-up) longwave flux (W/m2)
!        lprec = mass of liquid precipitation since last
!                      time step (Kg/m2)
!        fprec = mass of frozen precipitation since last
!                time step (Kg/m2)
!       runoff = mass (?) of runoff since last time step
!                       (Kg/m2)
!       p_surf = surface pressure (Pa)
!  </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_ice_to_ocean ( Time, Ice, Ocean, Ice_Ocean_Boundary )
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME=" Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>
!  <IN NAME="Ocean" TYPE="ocean_public_type">
!   A derived data type to specify ocean boundary data.
!  </IN>
!  <INOUT NAME="Ice_Ocean_Boundary" TYPE="ice_ocean_boundary_type">
!   A derived data type to specify properties and fluxes passed from ice to ocean.
!  </INOUT>
!
subroutine flux_ice_to_ocean ( Time, Ice, Ocean, Ice_Ocean_Boundary )

  type(time_type),        intent(in) :: Time
  type(ice_data_type),   intent(in)  :: Ice
  type(ocean_public_type), intent(in)  :: Ocean
!  real, dimension(:,:),   intent(out) :: flux_u_ocean,  flux_v_ocean,  &
!                                         flux_t_ocean,  flux_q_ocean,  &
!                                         flux_sw_ocean, flux_lw_ocean, &
!                                         lprec_ocean,   fprec_ocean,   &
!                                         runoff_ocean,  calving_ocean, &
!                                         flux_salt_ocean, p_surf_ocean
  type(ice_ocean_boundary_type), intent(inout) :: Ice_Ocean_Boundary
  real, dimension(:,:), pointer                :: dummy_null_pointer => NULL() !  RASF hack to get around pointer association.

  integer       :: m
  integer       :: n
  logical       :: used

!Balaji
  call mpp_clock_begin(cplOcnClock)
  call mpp_clock_begin(fluxIceOceanClock)

  if(ASSOCIATED(Ice_Ocean_Boundary%u_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%flux_u, Ice_Ocean_Boundary%u_flux, Ice_Ocean_Boundary%xtype, .FALSE. )

  if(ASSOCIATED(Ice_Ocean_Boundary%v_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%flux_v, Ice_Ocean_Boundary%v_flux, Ice_Ocean_Boundary%xtype, .FALSE. )

  if(ASSOCIATED(Ice_Ocean_Boundary%p     ) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%p_surf, Ice_Ocean_Boundary%p     , Ice_Ocean_Boundary%xtype, .FALSE. )

  if(ASSOCIATED(Ice_Ocean_Boundary%mi    ) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%mi,     Ice_Ocean_Boundary%mi    , Ice_Ocean_Boundary%xtype, .FALSE. )

  ! Extra fluxes
  do n = 1, Ice_Ocean_Boundary%fluxes%num_bcs  !{
     do m = 1, Ice_Ocean_Boundary%fluxes%bc(n)%num_fields  !{
        if ( associated(Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%values) ) then  !{
           call flux_ice_to_ocean_redistribute( Ice, Ocean, Ice%ocean_fluxes%bc(n)%field(m)%values, &
                Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%values, Ice_Ocean_Boundary%xtype, .FALSE. )
        endif  !}
     enddo  !} m
  enddo  !} n

  !--- The following variables may require conserved flux exchange from ice to ocean because the 
  !--- ice area maybe different from ocean area.
  if(ASSOCIATED(Ice_Ocean_Boundary%t_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%flux_t, Ice_Ocean_Boundary%t_flux, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

  if(ASSOCIATED(Ice_Ocean_Boundary%salt_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%flux_salt, Ice_Ocean_Boundary%salt_flux, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

  if(ASSOCIATED(Ice_Ocean_Boundary%sw_flux_nir_dir) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%flux_sw_nir_dir, Ice_Ocean_Boundary%sw_flux_nir_dir, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

  if(ASSOCIATED(Ice_Ocean_Boundary%sw_flux_nir_dif) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%flux_sw_nir_dif, Ice_Ocean_Boundary%sw_flux_nir_dif, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

  if(ASSOCIATED(Ice_Ocean_Boundary%sw_flux_vis_dir) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%flux_sw_vis_dir, Ice_Ocean_Boundary%sw_flux_vis_dir, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

  if(ASSOCIATED(Ice_Ocean_Boundary%sw_flux_vis_dif) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%flux_sw_vis_dif, Ice_Ocean_Boundary%sw_flux_vis_dif, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

  if(ASSOCIATED(Ice_Ocean_Boundary%lw_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%flux_lw, Ice_Ocean_Boundary%lw_flux, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

  if(ASSOCIATED(Ice_Ocean_Boundary%lprec) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%lprec, Ice_Ocean_Boundary%lprec, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

  if(ASSOCIATED(Ice_Ocean_Boundary%fprec) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%fprec, Ice_Ocean_Boundary%fprec, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

  if(ASSOCIATED(Ice_Ocean_Boundary%runoff) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%runoff, Ice_Ocean_Boundary%runoff, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

  if(ASSOCIATED(Ice_Ocean_Boundary%calving) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%calving, Ice_Ocean_Boundary%calving, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

  if(ASSOCIATED(Ice_Ocean_Boundary%runoff_hflx) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%runoff_hflx, Ice_Ocean_Boundary%runoff_hflx, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

  if(ASSOCIATED(Ice_Ocean_Boundary%calving_hflx) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%calving_hflx, Ice_Ocean_Boundary%calving_hflx, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

  if(ASSOCIATED(Ice_Ocean_Boundary%q_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
      Ice%flux_q, Ice_Ocean_Boundary%q_flux, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

  if(ASSOCIATED(Ice_Ocean_Boundary%wnd) ) then
     if(ASSOCIATED(Ice%wnd)) then
        call flux_ice_to_ocean_redistribute( Ice, Ocean, &
           Ice%wnd(:,:,1), Ice_Ocean_Boundary%wnd, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )
     else
        call flux_ice_to_ocean_redistribute( Ice, Ocean, &
          dummy_null_pointer, Ice_Ocean_Boundary%wnd, Ice_Ocean_Boundary%xtype, do_area_weighted_flux ) !Put dummy array here.
     endif
  endif

!Balaji: moved data_override calls here from coupler_main
  if( ocn_pe )then
      call mpp_set_current_pelist(ocn_pelist)
      call data_override('OCN', 'u_flux',    Ice_Ocean_Boundary%u_flux   , Time )
      call data_override('OCN', 'v_flux',    Ice_Ocean_Boundary%v_flux   , Time )
      call data_override('OCN', 't_flux',    Ice_Ocean_Boundary%t_flux   , Time )
      call data_override('OCN', 'q_flux',    Ice_Ocean_Boundary%q_flux   , Time )
      call data_override('OCN', 'wnd',    Ice_Ocean_Boundary%wnd   , Time )
      call data_override('OCN', 'salt_flux', Ice_Ocean_Boundary%salt_flux, Time )
      call data_override('OCN', 'lw_flux',   Ice_Ocean_Boundary%lw_flux  , Time )
      call data_override('OCN', 'sw_flux_nir_dir',   Ice_Ocean_Boundary%sw_flux_nir_dir  , Time )
      call data_override('OCN', 'sw_flux_nir_dif',   Ice_Ocean_Boundary%sw_flux_nir_dif  , Time )
      call data_override('OCN', 'sw_flux_vis_dir',   Ice_Ocean_Boundary%sw_flux_vis_dir  , Time )
      call data_override('OCN', 'sw_flux_vis_dif',   Ice_Ocean_Boundary%sw_flux_vis_dif  , Time )
      call data_override('OCN', 'lprec',     Ice_Ocean_Boundary%lprec    , Time )
      call data_override('OCN', 'fprec',     Ice_Ocean_Boundary%fprec    , Time )
      call data_override('OCN', 'runoff',    Ice_Ocean_Boundary%runoff   , Time )
      call data_override('OCN', 'calving',   Ice_Ocean_Boundary%calving  , Time )
      call data_override('OCN', 'runoff_hflx',    Ice_Ocean_Boundary%runoff_hflx   , Time )
      call data_override('OCN', 'calving_hflx',   Ice_Ocean_Boundary%calving_hflx  , Time )
      call data_override('OCN', 'p',         Ice_Ocean_Boundary%p        , Time )
      call data_override('OCN', 'mi',        Ice_Ocean_Boundary%mi       , Time )


! Extra fluxes
      do n = 1, Ice_Ocean_Boundary%fluxes%num_bcs  !{
         do m = 1, Ice_Ocean_Boundary%fluxes%bc(n)%num_fields  !{
            call data_override('OCN', Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%name,   &
                  Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%values, Time)
            used = send_data(Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%id_diag,        &
                   Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%values, Time )
         enddo  !} m
      enddo  !} n

!
!       Perform diagnostic output for the fluxes
!

     do n = 1, Ice_Ocean_Boundary%fluxes%num_bcs  !{
       do m = 1, Ice_Ocean_Boundary%fluxes%bc(n)%num_fields  !{
         used = send_data(Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%id_diag,                   &
                Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%values, Time)
       enddo  !} m
     enddo  !} n
   endif
   call mpp_set_current_pelist()

!Balaji
  call mpp_clock_end(fluxIceOceanClock)
  call mpp_clock_end(cplOcnClock)
!-----------------------------------------------------------------------

  end subroutine flux_ice_to_ocean
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="flux_ocean_to_ice">
!  <OVERVIEW>
!   Takes the ocean model state and interpolates it onto the bottom of the ice. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!    The following quantities are transferred from the Ocean to the ocean_ice_boundary_type: 
!
!        t_surf = surface temperature (deg K)
!        frazil = frazil (???)
!        u_surf = zonal ocean current/ice motion (m/s)
!        v_surf = meridional ocean current/ice motion (m/s)
!  </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_ocean_to_ice ( Time, Ocean, Ice, Ocean_Ice_Boundary)
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME=" Ocean" TYPE="ocean_public_type">
!   A derived data type to specify ocean boundary data.
!  </IN>
!  <IN NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>
!  <INOUT NAME="Ocean_Ice_Boundary" TYPE="ocean_ice_boundary_type">
!   A derived data type to specify properties and fluxes passed from ocean to ice.
!  </INOUT>
!
subroutine flux_ocean_to_ice ( Time, Ocean, Ice, Ocean_Ice_Boundary )

  type(time_type),         intent(in)  :: Time
  type(ocean_public_type), intent(in)  :: Ocean
  type(ice_data_type),     intent(in)  :: Ice
!  real, dimension(:,:),   intent(out) :: t_surf_ice, u_surf_ice, v_surf_ice, &
!                                         frazil_ice, s_surf_ice, sea_lev_ice
  type(ocean_ice_boundary_type), intent(inout) :: Ocean_Ice_Boundary
  real, dimension(nxc_ice, nyc_ice, nk_ice) :: ice_frac
  real, dimension(nxc_ocn, nyc_ocn )        :: tmp
  real, dimension(n_xgrid_sfc)              :: ex_ice_frac
  real, dimension(ni_atm, nj_atm)           :: diag_atm
  logical :: used
  integer       :: m
  integer       :: n
  real          :: from_dq 


!Balaji
  call mpp_clock_begin(cplOcnClock)
  call mpp_clock_begin(fluxOceanIceClock)

  select case (Ocean_Ice_Boundary%xtype)
  case(DIRECT)
     !same grid and domain decomp for ocean and ice    
     if( ASSOCIATED(Ocean_Ice_Boundary%u) )Ocean_Ice_Boundary%u = Ocean%u_surf
     if( ASSOCIATED(Ocean_Ice_Boundary%v) )Ocean_Ice_Boundary%v = Ocean%v_surf
     if( ASSOCIATED(Ocean_Ice_Boundary%t) )Ocean_Ice_Boundary%t = Ocean%t_surf
     if( ASSOCIATED(Ocean_Ice_Boundary%s) )Ocean_Ice_Boundary%s = Ocean%s_surf
     if( ASSOCIATED(Ocean_Ice_Boundary%sea_level) )Ocean_Ice_Boundary%sea_level = Ocean%sea_lev
     if( ASSOCIATED(Ocean_Ice_Boundary%frazil) ) then
        if(do_area_weighted_flux) then
           Ocean_Ice_Boundary%frazil = Ocean%frazil * Ocean%area 
           call divide_by_area(data=Ocean_Ice_Boundary%frazil, area=Ice%area)
        else
           Ocean_Ice_Boundary%frazil = Ocean%frazil
        endif
     endif

! Extra fluxes
     do n = 1, Ocean_Ice_Boundary%fields%num_bcs  !{
       do m = 1, Ocean_Ice_Boundary%fields%bc(n)%num_fields  !{
         if ( associated(Ocean_Ice_Boundary%fields%bc(n)%field(m)%values) ) then  !{
           Ocean_Ice_Boundary%fields%bc(n)%field(m)%values = Ocean%fields%bc(n)%field(m)%values
         endif  !}
       enddo  !} m
     enddo  !} n
  case(REDIST)
     !same grid, different domain decomp for ocean and ice    
     if( ASSOCIATED(Ocean_Ice_Boundary%u) )                     &
          call mpp_redistribute(Ocean%Domain, Ocean%u_surf, Ice%Domain, Ocean_Ice_Boundary%u)
     if( ASSOCIATED(Ocean_Ice_Boundary%v) )                     &
          call mpp_redistribute(Ocean%Domain, Ocean%v_surf, Ice%Domain, Ocean_Ice_Boundary%v)
     if( ASSOCIATED(Ocean_Ice_Boundary%t) )                     &
          call mpp_redistribute(Ocean%Domain, Ocean%t_surf, Ice%Domain, Ocean_Ice_Boundary%t)
     if( ASSOCIATED(Ocean_Ice_Boundary%s) )                     &
          call mpp_redistribute(Ocean%Domain, Ocean%s_surf, Ice%Domain, Ocean_Ice_Boundary%s)

     if( ASSOCIATED(Ocean_Ice_Boundary%sea_level) )             &
          call mpp_redistribute(Ocean%Domain, Ocean%sea_lev, Ice%Domain, Ocean_Ice_Boundary%sea_level)

     if( ASSOCIATED(Ocean_Ice_Boundary%frazil) ) then
        if(do_area_weighted_flux) then
           if(Ocean%is_ocean_pe)tmp = Ocean%frazil * Ocean%area 
           call mpp_redistribute( Ocean%Domain, tmp, Ice%Domain, Ocean_Ice_Boundary%frazil)
           if(Ice%pe) call divide_by_area(data=Ocean_Ice_Boundary%frazil, area=Ice%area)
        else
           call mpp_redistribute(Ocean%Domain, Ocean%frazil, Ice%Domain, Ocean_Ice_Boundary%frazil)
        endif
     endif

! Extra fluxes
     do n = 1, Ocean_Ice_Boundary%fields%num_bcs  !{
       do m = 1, Ocean_Ice_Boundary%fields%bc(n)%num_fields  !{
         if ( associated(Ocean_Ice_Boundary%fields%bc(n)%field(m)%values) ) then  !{
           call mpp_redistribute(Ocean%Domain, Ocean%fields%bc(n)%field(m)%values,    &
                Ice%Domain, Ocean_Ice_Boundary%fields%bc(n)%field(m)%values)
         endif  !}
       enddo  !} m
     enddo  !} n
  case DEFAULT
!   <ERROR MSG="Ocean_Ice_Boundary%xtype must be DIRECT or REDIST." STATUS="FATAL">
!     The value of variable xtype of ice_ocean_boundary_type data must be DIRECT or REDIST.
!   </ERROR>
     call mpp_error( FATAL, 'FLUX_OCEAN_TO_ICE: Ocean_Ice_Boundary%xtype must be DIRECT or REDIST.' )
  end select
  if( ice_pe )then
      call mpp_set_current_pelist(ice_pelist)

!Balaji: data_override moved here from coupler_main
      call data_override('ICE', 'u',         Ocean_Ice_Boundary%u,         Time)
      call data_override('ICE', 'v',         Ocean_Ice_Boundary%v,         Time)
      call data_override('ICE', 't',         Ocean_Ice_Boundary%t,         Time)
      call data_override('ICE', 's',         Ocean_Ice_Boundary%s,         Time)
      call data_override('ICE', 'frazil',    Ocean_Ice_Boundary%frazil,    Time)
      call data_override('ICE', 'sea_level', Ocean_Ice_Boundary%sea_level, Time)

! Extra fluxes
      do n = 1, Ocean_Ice_Boundary%fields%num_bcs  !{
         do m = 1, Ocean_Ice_Boundary%fields%bc(n)%num_fields  !{
            call data_override('ICE', Ocean_Ice_Boundary%fields%bc(n)%field(m)%name,    &
                 Ocean_Ice_Boundary%fields%bc(n)%field(m)%values, Time)
         enddo  !} m
      enddo  !} n

!
!       Perform diagnostic output for the ocean_ice_boundary fields
!

     do n = 1, Ocean_Ice_Boundary%fields%num_bcs  !{
       do m = 1, Ocean_Ice_Boundary%fields%bc(n)%num_fields  !{
           used = send_data(Ocean_Ice_Boundary%fields%bc(n)%field(m)%id_diag,                   &
                Ocean_Ice_Boundary%fields%bc(n)%field(m)%values, Time)
       enddo  !} m
     enddo  !} n
   endif

  if( ocn_pe )then
      call mpp_set_current_pelist(ocn_pelist)

!
!       Perform diagnostic output for the ocean fields
!

     do n = 1, Ocean%fields%num_bcs  !{
       do m = 1, Ocean%fields%bc(n)%num_fields  !{
         used = send_data(Ocean%fields%bc(n)%field(m)%id_diag,                                &
                Ocean%fields%bc(n)%field(m)%values, Time)
       enddo  !} m
     enddo  !} n
   endif

   call mpp_set_current_pelist()
  
  if ( id_ice_mask > 0 ) then
     ice_frac        = 1.
     ice_frac(:,:,1) = 0.
     ex_ice_frac     = 0.
     call put_to_xgrid (ice_frac, 'OCN', ex_ice_frac, xmap_sfc)
     call get_from_xgrid (diag_atm, 'ATM', ex_ice_frac, xmap_sfc)
     used = send_data ( id_ice_mask, diag_atm, Time )
  endif

  if(Ice%pe) then
     ! frazil (already in J/m^2 so no need to multiply by Dt_cpl)
     from_dq = 4*PI*Radius*Radius * &
         & SUM( ice_cell_area * Ocean_Ice_Boundary%frazil )
     Ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) = Ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) - from_dq
     Ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) + from_dq
  endif

!Balaji
  call mpp_clock_end(fluxOceanIceClock)
  call mpp_clock_end(cplOcnClock)
!-----------------------------------------------------------------------

  end subroutine flux_ocean_to_ice
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="flux_check_stocks">
!  <OVERVIEW>
!   Check stock values. 
!  </OVERVIEW>
!  <DESCRIPTION>
!   Will print out any difference between the integrated flux (in time
!   and space) feeding into a component, and the stock stored in that
!   component.
!  </DESCRIPTION>

  subroutine flux_check_stocks(Time, Atm, Lnd, Ice, Ocn_state)

    type(time_type)       :: Time
    type(atmos_data_type), optional :: Atm
    type(land_data_type), optional  :: Lnd
    type(ice_data_type), optional   :: Ice
    type(ocean_state_type), optional, pointer :: Ocn_state

    real :: ref_value
    integer :: i, ier


    do i = 1, NELEMS

       if(present(Atm)) then
          ref_value = 0
          call Atm_stock_pe(Atm, index=i, value=ref_value)        
          if(i==ISTOCK_WATER .and. Atm%pe ) then
             ! decrease the Atm stock by the precip adjustment to reflect the fact that
             ! after an update_atmos_up call, the precip will be that of the future time step.
             ! Thus, the stock call will represent the (explicit ) precip at 
             ! the beginning of the preceding time step, and the (implicit) evap at the 
             ! end of the preceding time step
             call stock_integrate_2d(Atm%lprec + Atm%fprec, xmap=xmap_sfc, delta_t=Dt_atm, &
                  & radius=Radius, res=ATM_PRECIP_NEW, ier=ier)

             ref_value = ref_value + ATM_PRECIP_NEW
          endif
         
          Atm_stock(i)%q_now = ref_value
       endif

       if(present(Lnd)) then
          ref_value = 0
          call Lnd_stock_pe(Lnd, index=i, value=ref_value)
          Lnd_stock(i)%q_now = ref_value
       endif

       if(present(Ice)) then
          ref_value = 0
          call Ice_stock_pe(Ice, index=i, value=ref_value)
          Ice_stock(i)%q_now = ref_value
       endif

       if(present(Ocn_state)) then
          ref_value = 0
          call Ocean_stock_pe(Ocn_state, index=i, value=ref_value)
          Ocn_stock(i)%q_now = ref_value
       endif
    enddo

    call stocks_report(Time)


  end subroutine flux_check_stocks
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="flux_init_stocks">
!  <OVERVIEW>
!   Initialize stock values. 
!  </OVERVIEW>
!  <DESCRIPTION>
!   This will call the various component stock_pe routines to store the 
!   the initial stock values.
!  </DESCRIPTION>

subroutine flux_init_stocks(Time, Atm, Lnd, Ice, Ocn_state)
  type(time_type) , intent(in) :: Time
  type(atmos_data_type) :: Atm
  type(land_data_type)  :: Lnd
  type(ice_data_type)   :: Ice
  type(ocean_state_type), pointer :: Ocn_state

  integer i, ier

  stocks_file=stdout()
! Divert output file for stocks if requested 
  if(mpp_pe()==mpp_root_pe() .and. divert_stocks_report) then
     call mpp_open( stocks_file, 'stocks.out', action=MPP_OVERWR, threading=MPP_SINGLE, &
          fileset=MPP_SINGLE, nohdrs=.TRUE. )       
  endif
  
    ! Initialize stock values
    do i = 1, NELEMS
       call Atm_stock_pe(   Atm , index=i, value=Atm_stock(i)%q_start)

       if(i==ISTOCK_WATER .and. Atm%pe ) then
          call stock_integrate_2d(Atm%lprec + Atm%fprec, xmap=xmap_sfc, & 
               delta_t=Dt_atm, radius=Radius, res=ATM_PRECIP_NEW, ier=ier) 
          
          Atm_stock(i)%q_start = Atm_stock(i)%q_start + ATM_PRECIP_NEW
       endif

       call Lnd_stock_pe(   Lnd , index=i, value=Lnd_stock(i)%q_start)
       call Ice_stock_pe(   Ice , index=i, value=Ice_stock(i)%q_start)
       call Ocean_stock_pe( Ocn_state , index=i, value=Ocn_stock(i)%q_start)
    enddo


    call stocks_report_init(Time)


end subroutine flux_init_stocks
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="generate_sfc_xgrid">
!  <OVERVIEW>
!   Optimizes the exchange grids by eliminating land and ice partitions with no data. 
!  </OVERVIEW>
!  <DESCRIPTION>
!   Optimizes the exchange grids by eliminating land and ice partitions with no data. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call generate_sfc_xgrid( Land, Ice )
!		
!  </TEMPLATE>
!  <IN NAME=" Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </IN>
!  <IN NAME="Ice" TYPE="ice_data_type">
!  A derived data type to specify ice boundary data.
!  </IN>
!
subroutine generate_sfc_xgrid( Land, Ice )
! subroutine to regenerate exchange grid eliminating side 2 tiles with 0 frac area
    type(land_data_type), intent(in) :: Land
    type(ice_data_type),  intent(in) :: Ice

    integer :: isc, iec, jsc, jec

!Balaji
  call mpp_clock_begin(cplClock)
  call mpp_clock_begin(regenClock)

  call mpp_get_compute_domain(Ice%Domain, isc, iec, jsc, jec)

  call set_frac_area (Ice%part_size(isc:iec,jsc:jec,:) , 'OCN', xmap_sfc)
  call set_frac_area (Land%tile_size, 'LND', xmap_sfc)
  n_xgrid_sfc = max(xgrid_count(xmap_sfc),1)

!Balaji
  call mpp_clock_end(regenClock)
  call mpp_clock_end(cplClock)
  return
end subroutine generate_sfc_xgrid
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="flux_up_to_atmos">
!  <OVERVIEW>
!   Corrects the fluxes for consistency with the new surface temperatures in land 
!   and ice models.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Corrects the fluxes for consistency with the new surface temperatures in land 
!   and ice models. Final increments for temperature and specific humidity in the 
!   lowest atmospheric layer are computed and returned to the atmospheric model
!   so that it can finalize the increments in the rest of the atmosphere. 
!  <PRE>
!
!   The following elements of the land_ice_atmos_boundary_type are computed:
!        dt_t  = temperature change at the lowest
!                 atmospheric level (deg k)
!        dt_q  = specific humidity change at the lowest
!                 atmospheric level (kg/kg)
!  </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_up_to_atmos ( Time, Land, Ice, Land_Ice_Atmos_Boundary, Land_boundary, Ice_boundary )
!		
!  </TEMPLATE>
!  <IN NAME=" Time" TYPE="time_type">
!   Current time.
!  </IN>
!  <INOUT NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </INOUT>
!  <INOUT NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </INOUT>
!  <INOUT NAME="Land_Ice_Atmos_Boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere, land and ice. 
!  </INOUT>
!
subroutine flux_up_to_atmos ( Time, Land, Ice, Land_Ice_Atmos_Boundary, Land_boundary, Ice_boundary )

  type(time_type),      intent(in)  :: Time
  type(land_data_type), intent(inout)  :: Land
  type(ice_data_type),  intent(inout)  :: Ice
  type(land_ice_atmos_boundary_type), intent(inout) :: Land_Ice_Atmos_Boundary
  type(atmos_land_boundary_type) :: Land_boundary
  type(atmos_ice_boundary_type)  :: Ice_boundary

  real, dimension(n_xgrid_sfc) ::  &
       ex_t_surf_new, &
       ex_dt_t_surf,  &
       ex_delta_t_n,  &
       ex_t_ca_new,   &
       ex_dt_t_ca
  real, dimension(n_xgrid_sfc,n_exch_tr) :: &
       ex_tr_surf_new,    & ! updated tracer values at the surface
       ex_dt_tr_surf,     & ! tendency of tracers at the surface
       ex_delta_tr_n
! jgj: added for co2_surf diagnostic 
  real, dimension(n_xgrid_sfc) :: &
       ex_co2_surf_dvmr   ! updated CO2 tracer values at the surface (dry vmr)

  real, dimension(size(Land_Ice_Atmos_Boundary%dt_t,1),size(Land_Ice_Atmos_Boundary%dt_t,2)) :: diag_atm, &
       evap_atm
  real, dimension(size(Land_boundary%lprec,1), size(Land_boundary%lprec,2), size(Land_boundary%lprec,3)) :: data_lnd, diag_land
  real, dimension(size(Ice_boundary%lprec,1), size(Ice_boundary%lprec,2), size(Ice_boundary%lprec,3)) :: data_ice
  logical :: used

  integer :: tr       ! tracer index
  character(32) :: tr_name ! tracer name
  integer :: n, i, m, ier


  !Balaji
  call mpp_clock_begin(cplClock)
  call mpp_clock_begin(fluxAtmUpClock)
  !-----------------------------------------------------------------------
  !Balaji: data_override calls moved here from coupler_main
  call data_override ( 'ICE', 't_surf', Ice%t_surf,  Time)
  call data_override ( 'LND', 't_ca',   Land%t_ca,   Time)
  call data_override ( 'LND', 't_surf', Land%t_surf, Time)
  do tr = 1, n_lnd_tr
     call get_tracer_names( MODEL_LAND, tr, tr_name )
     call data_override('LND', trim(tr_name)//'_surf', Land%tr(:,:,:,tr), Time)
  enddo

  !----- compute surface temperature change -----

  ex_t_surf_new = 200.0

  call put_to_xgrid (Ice%t_surf,  'OCN', ex_t_surf_new, xmap_sfc)
  ex_t_ca_new = ex_t_surf_new  ! since it is the same thing over oceans
  call put_to_xgrid (Land%t_ca,   'LND', ex_t_ca_new,   xmap_sfc)
  call put_to_xgrid (Land%t_surf, 'LND', ex_t_surf_new, xmap_sfc)

  !  call escomp(ex_t_ca_new, ex_q_surf_new)
  !  ex_q_surf_new  = d622*ex_q_surf_new/(ex_p_surf-d378*ex_q_surf_new) 
  !  call put_to_xgrid (Land%q_ca, 'LND', ex_q_surf_new, xmap_sfc)

#ifdef SCM
  if (do_specified_flux .and. do_specified_land) then
       ex_t_surf_new = ex_t_surf
       ex_t_ca_new   = ex_t_ca
  endif
#endif

  where (ex_avail)
     ex_dt_t_ca   = ex_t_ca_new   - ex_t_ca   ! changes in near-surface T
     ex_dt_t_surf = ex_t_surf_new - ex_t_surf ! changes in radiative T
  endwhere

  if (do_forecast) then
     where (ex_avail(:) .and. (.not.ex_land(:)))
        ex_dt_t_ca  (:) = 0.
        ex_dt_t_surf(:) = 0.
     end where
  end if

  !-----------------------------------------------------------------------
  !-----  adjust fluxes and atmospheric increments for 
  !-----  implicit dependence on surface temperature -----
  do tr = 1,n_exch_tr
     ! set up updated surface tracer field so that flux to atmos for absent 
     ! tracers is zero
     do i = 1, size(ex_avail(:))
        if(.not.ex_avail(i)) cycle
        if (ex_dfdtr_surf(i,tr)/=0) then
           ex_dt_tr_surf(i,tr) = -ex_flux_tr(i,tr)/ex_dfdtr_surf(i,tr)
        else
           ex_dt_tr_surf(i,tr) = 0
        endif
        ex_tr_surf_new(i,tr) = ex_tr_surf(i,tr)+ex_dt_tr_surf(i,tr)
     enddo
     ! get all tracers available from land, and calculate changes in near-tracer field
     n = tr_table(tr)%lnd
     if(n /= NO_TRACER ) then
        call put_to_xgrid ( Land%tr(:,:,:,n), 'LND', ex_tr_surf_new(:,tr), xmap_sfc )
     endif

     ! get all tracers available from ocean here 

     ! update tracer tendencies in the atmosphere
     where (ex_avail)
        ex_dt_tr_surf(:,tr) = ex_tr_surf_new(:,tr) - ex_tr_surf(:,tr)
        ex_delta_tr_n(:,tr) = ex_f_tr_delt_n(:,tr) + ex_dt_tr_surf(:,tr) * ex_e_tr_n(:,tr)
        ex_flux_tr(:,tr)    = ex_flux_tr(:,tr)     + ex_dt_tr_surf(:,tr) * ex_dfdtr_surf(:,tr)
     endwhere
  enddo

  ! re-calculate fluxes of specific humidity over ocean
  where (ex_avail.and..not.ex_land) 
     ! note that in this region (over ocean) ex_dt_t_surf == ex_dt_t_ca
     ex_delta_tr_n(:,isphum)  = ex_f_tr_delt_n(:,isphum) + ex_dt_t_surf * ex_e_q_n
     ex_flux_tr(:,isphum)     = ex_flux_tr(:,isphum)     + ex_dt_t_surf * ex_dedt_surf
  endwhere

  do tr=1,n_exch_tr
     ! get updated tracer tendency on the atmospheic grid
     n=tr_table(tr)%atm
     call get_from_xgrid (Land_Ice_Atmos_Boundary%dt_tr(:,:,n), 'ATM', ex_delta_tr_n(:,tr), xmap_sfc)
  enddo

  ex_delta_t_n = 0.0

  where(ex_avail)
     ex_flux_t     = ex_flux_t  + ex_dt_t_ca   * ex_dhdt_surf
     ex_flux_lw    = ex_flux_lw - ex_dt_t_surf * ex_drdt_surf
     ex_delta_t_n  = ex_f_t_delt_n  + ex_dt_t_ca*ex_e_t_n
  endwhere

  !-----------------------------------------------------------------------
  !---- get mean quantites on atmospheric grid ----

  call get_from_xgrid (Land_Ice_Atmos_Boundary%dt_t, 'ATM', ex_delta_t_n, xmap_sfc)

  !=======================================================================
  !-------------------- diagnostics section ------------------------------

  !------- new surface temperature -----------
  if ( id_t_surf > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_t_surf_new, xmap_sfc)
     used = send_data ( id_t_surf, diag_atm, Time )
  endif


  ! + slm, Mar 27 2002
  ! ------ new canopy temperature --------
  !   NOTE, that in the particular case of LM2 t_ca is identical to t_surf,
  !   but this will be changed in future version of the land madel
  if ( id_t_ca > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_t_ca_new, xmap_sfc)
     used = send_data ( id_t_ca, diag_atm, Time )
  endif

  !------- updated surface tracer fields ------
  do tr=1,n_exch_tr
     call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
     if ( id_tr_surf(tr) > 0 ) then
        call get_from_xgrid (diag_atm, 'ATM', ex_tr_surf_new(:,tr), xmap_sfc)
        used = send_data ( id_tr_surf(tr), diag_atm, Time )
     endif
!!jgj:  add dryvmr co2_surf
! - slm Mar 25, 2010: moved to resolve interdependence of diagnostic fields
     if ( id_co2_surf_dvmr > 0 .and. lowercase(trim(tr_name))=='co2') then
       ex_co2_surf_dvmr = (ex_tr_surf_new(:,tr) / (1.0 - ex_tr_surf_new(:,isphum))) * WTMAIR/WTMCO2
       call get_from_xgrid (diag_atm, 'ATM', ex_co2_surf_dvmr, xmap_sfc)
       used = send_data ( id_co2_surf_dvmr, diag_atm, Time )
     endif
  enddo

  !------- sensible heat flux -----------
  if ( id_t_flux > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_flux_t, xmap_sfc)
     used = send_data ( id_t_flux, diag_atm, Time )
  endif

  !------- net longwave flux -----------
  if ( id_r_flux > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_flux_lw, xmap_sfc)
     used = send_data ( id_r_flux, diag_atm, Time )
  endif

  !------- tracer fluxes ------------
  ! tr_mol_flux diagnostic will be correct for co2 tracer only. 
  ! will need update code to use correct molar mass for tracers other than co2
  do tr=1,n_exch_tr
     if ( id_tr_flux(tr) > 0 .or. id_tr_mol_flux(tr) > 0 ) then
        call get_from_xgrid (diag_atm, 'ATM', ex_flux_tr(:,tr), xmap_sfc)
        if (id_tr_flux(tr) > 0 ) &
            used = send_data ( id_tr_flux(tr), diag_atm, Time )
        if (id_tr_mol_flux(tr) > 0 ) &
            used = send_data ( id_tr_mol_flux(tr), diag_atm*1000./WTMCO2, Time)
     endif
  enddo

  !-----------------------------------------------------------------------
  !---- accumulate global integral of evaporation (mm/day) -----
  call get_from_xgrid (evap_atm, 'ATM', ex_flux_tr(:,isphum), xmap_sfc)
  if( id_q_flux > 0 ) used = send_data ( id_q_flux, evap_atm, Time)
  if( id_q_flux_land > 0 ) then
     call get_from_xgrid (diag_land, 'LND', ex_flux_tr(:,isphum), xmap_sfc)
     used = send_tile_averaged_data(id_q_flux_land, diag_land, &
          Land%tile_size, Time, mask=Land%mask)
  endif
  call sum_diag_integral_field ('evap', evap_atm*86400.)

  ! compute stock changes

  call get_from_xgrid(data_lnd, 'LND', ex_flux_tr(:,isphum), xmap_sfc)

  ! Lnd -> Atm (evap)
  call stock_move( &
       & TO   = Atm_stock(ISTOCK_WATER), &
       & FROM = Lnd_stock(ISTOCK_WATER), &
       & DATA = data_lnd, &
       & grid_index=X1_GRID_LND, &
       & xmap=xmap_sfc, &
       & delta_t=Dt_atm, &
       & to_side=ISTOCK_SIDE, from_side=ISTOCK_TOP, &
       & radius=Radius, ier=ier, verbose='stock move EVAP (Lnd->ATm) ')

  ! Lnd -> Atm (heat lost through evap)
  call stock_move( &
       & TO   = Atm_stock(ISTOCK_HEAT), &
       & FROM = Lnd_stock(ISTOCK_HEAT), &
       & DATA = data_lnd * HLV, &
       & grid_index=X1_GRID_LND, &
       & xmap=xmap_sfc, &
       & delta_t=Dt_atm, &
       & to_side=ISTOCK_SIDE, from_side=ISTOCK_TOP, &
       & radius=Radius, ier=ier, verbose='stock move EVAP*HLV (Lnd->ATm) ')

  call get_from_xgrid(data_ice, 'OCN', ex_flux_tr(:,isphum), xmap_sfc)

  ! Ice -> Atm (evap)
  call stock_move( &
       & TO   = Atm_stock(ISTOCK_WATER), &
       & FROM = Ice_stock(ISTOCK_WATER), &
       & DATA = data_ice, &
       & grid_index=X1_GRID_ICE, &
       & xmap=xmap_sfc, &
       & delta_t=Dt_atm, &
       & to_side=ISTOCK_TOP, from_side=ISTOCK_TOP, &
       & radius=Radius, ier=ier, verbose='stock move EVAP (Ice->ATm) ')

  ! Ice -> Atm (heat lost through evap)
  call stock_move( &
       & TO   = Atm_stock(ISTOCK_HEAT), &
       & FROM = Ice_stock(ISTOCK_HEAT), &
       & DATA = data_ice * HLV, &
       & grid_index=X1_GRID_ICE, &
       & xmap=xmap_sfc, &
       & delta_t=Dt_atm, &
       & to_side=ISTOCK_TOP, from_side=ISTOCK_TOP, &
       & radius=Radius, ier=ier, verbose='stock move EVAP*HLV (Ice->ATm) ')

  !=======================================================================
  !---- deallocate module storage ----
  deallocate ( &
       ex_t_surf   ,  &
       ex_t_surf_miz, &
       ex_p_surf   ,  &
       ex_slp      ,  &
       ex_t_ca     ,  &
       ex_dhdt_surf,  &
       ex_dedt_surf,  &
       ex_dqsatdt_surf,  &
       ex_drdt_surf,  &
       ex_dhdt_atm ,  &
       ex_flux_t   ,  &
       ex_flux_lw  ,  &
       ex_drag_q   ,  &
       ex_avail    ,  &
       ex_f_t_delt_n, &
       ex_tr_surf  ,  &
       
  ex_dfdtr_surf  , &
       ex_dfdtr_atm   , &
       ex_flux_tr     , &
       ex_f_tr_delt_n , &
       ex_e_tr_n      , &
       
  ex_e_t_n    ,  &
       ex_e_q_n    ,  &
       ! values added for LM3
       ex_cd_t     ,  &
       ex_cd_m     ,  &
       ex_b_star   ,  &
       ex_u_star   ,  &
       ex_wind     ,  &
       ex_z_atm    ,  &
       
  ex_land        )

#ifdef SCM
  deallocate ( &
       ex_dhdt_surf_forland, &
       ex_dedt_surf_forland, &
       ex_dedq_surf_forland  )
#endif

! Extra fluxes
  do n = 1, ex_gas_fields_ice%num_bcs  !{
     do m = 1, ex_gas_fields_ice%bc(n)%num_fields  !{
        deallocate ( ex_gas_fields_ice%bc(n)%field(m)%values )
        nullify ( ex_gas_fields_ice%bc(n)%field(m)%values )
     enddo  !} m
  enddo  !} n

  do n = 1, ex_gas_fields_atm%num_bcs  !{
     do m = 1, ex_gas_fields_atm%bc(n)%num_fields  !{
        deallocate ( ex_gas_fields_atm%bc(n)%field(m)%values )
        nullify ( ex_gas_fields_atm%bc(n)%field(m)%values )
     enddo  !} m
  enddo  !} n

  do n = 1, ex_gas_fluxes%num_bcs  !{
     do m = 1, ex_gas_fluxes%bc(n)%num_fields  !{
        deallocate ( ex_gas_fluxes%bc(n)%field(m)%values )
        nullify ( ex_gas_fluxes%bc(n)%field(m)%values )
     enddo  !} m
  enddo  !} n

!Balaji
  call mpp_clock_end(fluxAtmUpClock)
  call mpp_clock_end(cplClock)

!-----------------------------------------------------------------------

end subroutine flux_up_to_atmos
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="flux_ice_to_ocean_stocks">
!  <OVERVIEW>
!   Updates Ice and Ocean stocks.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Integrate the fluxes over the surface and in time. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_ice_to_ocean_stocks ( Ice )
!		
!  </TEMPLATE>
!  <IN NAME=" Time" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>

  subroutine flux_ice_to_ocean_stocks(Ice)

    type(ice_data_type),   intent(in)  :: Ice

    real           :: from_dq

    ! fluxes from ice -> ocean, integrate over surface and in time 

    ! precip - evap
    from_dq = 4*PI*Radius*Radius * Dt_cpl * &
         & SUM( ice_cell_area * (Ice%lprec+Ice%fprec-Ice%flux_q) )
    Ice_stock(ISTOCK_WATER)%dq(ISTOCK_BOTTOM) = Ice_stock(ISTOCK_WATER)%dq(ISTOCK_BOTTOM) - from_dq
    Ocn_stock(ISTOCK_WATER)%dq(ISTOCK_TOP   ) = Ocn_stock(ISTOCK_WATER)%dq(ISTOCK_TOP   ) + from_dq

    ! river
    from_dq = 4*PI*Radius*Radius * Dt_cpl * &
         & SUM( ice_cell_area * (Ice%runoff + Ice%calving) )
    Ice_stock(ISTOCK_WATER)%dq(ISTOCK_BOTTOM) = Ice_stock(ISTOCK_WATER)%dq(ISTOCK_BOTTOM) - from_dq
    Ocn_stock(ISTOCK_WATER)%dq(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_WATER)%dq(ISTOCK_SIDE  ) + from_dq

    ! sensible heat + shortwave + longwave + latent heat
    from_dq = 4*PI*Radius*Radius * Dt_cpl * &
         & SUM( ice_cell_area * ( &
         &   Ice%flux_sw_vis_dir+Ice%flux_sw_vis_dif &
         & + Ice%flux_sw_nir_dir+Ice%flux_sw_nir_dif + Ice%flux_lw &
         & - (Ice%fprec + Ice%calving)*HLF - Ice%flux_t - Ice%flux_q*HLV) )
    Ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) = Ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) - from_dq
    Ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) + from_dq

    ! heat carried by river + pme (assuming reference temperature of 0 degC and river/pme temp = surface temp)
    ! Note: it does not matter what the ref temperature is but it must be consistent with that in OCN and ICE
    from_dq = 4*PI*Radius*Radius * Dt_cpl * &
         & SUM( ice_cell_area * ( &
         & (Ice%lprec+Ice%fprec-Ice%flux_q + Ice%runoff+Ice%calving)*CP_OCEAN*(Ice%t_surf(:,:,1) - 273.15)) )
    Ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) = Ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) - from_dq
    Ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) + from_dq

    !SALT flux
    from_dq = Dt_cpl* SUM( ice_cell_area * ( -Ice%flux_salt )) *4*PI*Radius*Radius
    Ice_stock(ISTOCK_SALT)%dq(ISTOCK_BOTTOM) = Ice_stock(ISTOCK_SALT)%dq(ISTOCK_BOTTOM) - from_dq
    Ocn_stock(ISTOCK_SALT)%dq(ISTOCK_TOP   ) = Ocn_stock(ISTOCK_SALT)%dq(ISTOCK_TOP   ) + from_dq


  end subroutine flux_ice_to_ocean_stocks
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="flux_ocean_from_ice_stocks">
!  <OVERVIEW>
!   Updates Ocean stocks due to input that the Ocean model gets.
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine updates the stocks of Ocean by the amount of input that the Ocean gets from Ice component.
!   Unlike subroutine flux_ice_to_ocean_stocks() that uses Ice%fluxes to update the stocks due to the amount of output from Ice
!   this subroutine uses Ice_Ocean_boundary%fluxes to calculate the amount of input to the Ocean. These fluxes are the ones
!   that Ocean model uses internally to calculate its budgets. Hence there should be no difference between this input and what 
!   Ocean model internal diagnostics uses. 
!   This bypasses the possible mismatch in cell areas between Ice and Ocean in diagnosing the stocks of Ocean
!   and should report a conserving Ocean component regardless of the glitches in fluxes.
!
!   The use of this subroutine in conjunction with  subroutine flux_ice_to_ocean_stocks() will also allow to directly
!   diagnose the amount "stocks lost in exchange" between Ice and Ocean
!
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_ocean_from_ice_stocks(ocean_state,Ocean,Ice_Ocean_boundary)
!		
!  </TEMPLATE>
  subroutine flux_ocean_from_ice_stocks(ocean_state,Ocean,Ice_Ocean_boundary)
    type(ocean_state_type),        pointer    :: ocean_state
    type(ocean_public_type),       intent(in) :: Ocean
    type(ice_ocean_boundary_type), intent(in) :: Ice_Ocean_Boundary
    real    :: from_dq, cp_ocn
    real, dimension(nxc_ocn, nyc_ocn) :: ocean_cell_area, wet, t_surf, t_pme, t_calving, t_runoff, btfHeat
    integer :: isc, iec, jsc, jec

    call mpp_get_compute_domain(Ocean%Domain, isc, iec, jsc, jec)
    call ocean_model_data_get(ocean_state,Ocean,'area'  , ocean_cell_area,isc,jsc)
    call ocean_model_data_get(ocean_state,Ocean,'mask', wet,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'t_surf', t_surf,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'t_runoff', t_runoff,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'t_pme', t_pme,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'t_calving', t_calving,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'btfHeat', btfHeat,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'c_p', cp_ocn )


    ! fluxes from ice -> ocean, integrate over surface and in time 

    ! precip - evap
    from_dq = SUM( ocean_cell_area * wet * (Ice_Ocean_Boundary%lprec+Ice_Ocean_Boundary%fprec-Ice_Ocean_Boundary%q_flux) )
    Ocn_stock(ISTOCK_WATER)%dq_IN(ISTOCK_TOP   ) = Ocn_stock(ISTOCK_WATER)%dq_IN(ISTOCK_TOP   ) + from_dq * Dt_cpl

    from_dq = SUM( ocean_cell_area * wet * (Ice_Ocean_Boundary%runoff+Ice_Ocean_Boundary%calving) )
    Ocn_stock(ISTOCK_WATER)%dq_IN(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_WATER)%dq_IN(ISTOCK_SIDE  ) + from_dq * Dt_cpl

    ! sensible heat + shortwave + longwave + latent heat

    from_dq = SUM( ocean_cell_area * wet *( Ice_Ocean_Boundary%sw_flux_vis_dir + Ice_Ocean_Boundary%sw_flux_vis_dif &
                                           +Ice_Ocean_Boundary%sw_flux_nir_dir + Ice_Ocean_Boundary%sw_flux_nir_dif &
                                           +Ice_Ocean_Boundary%lw_flux &
                                           - (Ice_Ocean_Boundary%fprec + Ice_Ocean_Boundary%calving)*HLF &
                                           - Ice_Ocean_Boundary%t_flux - Ice_Ocean_Boundary%q_flux*HLV ))

    Ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE  ) + from_dq * Dt_cpl

    ! heat carried by river + pme (assuming reference temperature of 0 degC and river/pme temp = surface temp)
    ! Note: it does not matter what the ref temperature is but it must be consistent with that in OCN and ICE

    from_dq = SUM( ocean_cell_area * wet * cp_ocn *&
                              ((Ice_Ocean_Boundary%lprec+Ice_Ocean_Boundary%fprec-Ice_Ocean_Boundary%q_flux)*t_pme &
                               +Ice_Ocean_Boundary%calving * t_calving &
                               +Ice_Ocean_Boundary%runoff  * t_runoff  ))       
    
    Ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE ) + from_dq * Dt_cpl

!   Bottom heat flux
    from_dq = - SUM( ocean_cell_area * wet * btfHeat)
    
    Ocn_stock(ISTOCK_HEAT)%dq_IN( ISTOCK_BOTTOM ) = Ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_BOTTOM ) + from_dq * Dt_cpl

!   Frazil heat

     from_dq =  SUM( ocean_cell_area *wet * Ocean%frazil )
     Ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE ) + from_dq

    !SALT flux
    from_dq = SUM( ocean_cell_area * wet * ( -Ice_Ocean_Boundary%salt_flux))  
    Ocn_stock(ISTOCK_SALT)%dq_IN(ISTOCK_TOP  ) = Ocn_stock(ISTOCK_SALT)%dq_IN(ISTOCK_TOP   ) + from_dq  * Dt_cpl


  end subroutine flux_ocean_from_ice_stocks
! </SUBROUTINE>


!#######################################################################

subroutine put_logical_to_real (mask, id, ex_mask, xmap)

  logical         , intent(in)    :: mask(:,:,:)
  character(len=3), intent(in)    :: id
  real            , intent(inout) :: ex_mask(:)
  type(xmap_type), intent(inout) :: xmap

  !-----------------------------------------------------------------------
  !    puts land or ice model masks (with partitions) onto the
  !    exchange grid as a real array (1.=true, 0.=false)
  !-----------------------------------------------------------------------

  real, dimension(size(mask,1),size(mask,2),size(mask,3)) :: rmask
  
  where (mask)
     rmask = 1.0
  elsewhere
     rmask = 0.0
  endwhere

  call put_to_xgrid(rmask, id, ex_mask, xmap)

end subroutine put_logical_to_real

!#######################################################################

subroutine diag_field_init ( Time, atmos_axes, land_axes )

  type(time_type), intent(in) :: Time
  integer,         intent(in) :: atmos_axes(2)
  integer,         intent(in) :: land_axes(2)

  integer :: iref
  character(len=6) :: label_zm, label_zh
  real, dimension(2) :: trange = (/  100., 400. /), &
       vrange = (/ -400., 400. /), &
       frange = (/ -0.01, 1.01 /)
  character(len=32)  :: name, units ! name of the tracer
  character(len=128) :: longname    ! long name of the tracer
  integer            :: tr          ! tracer index
!-----------------------------------------------------------------------
!  initializes diagnostic fields that may be output from this module
!  (the id numbers may be referenced anywhere in this module)
!-----------------------------------------------------------------------

  !------ labels for diagnostics -------
  !  (z_ref_mom, z_ref_heat are namelist variables)

  iref = int(z_ref_mom+0.5)
  if ( real(iref) == z_ref_mom ) then
     write (label_zm,105) iref
     if (iref < 10) write (label_zm,100) iref
  else
     write (label_zm,110) z_ref_mom
  endif

  iref = int(z_ref_heat+0.5)
  if ( real(iref) == z_ref_heat ) then
     write (label_zh,105) iref
     if (iref < 10) write (label_zh,100) iref
  else
     write (label_zh,110) z_ref_heat
  endif

100 format (i1,' m',3x)
105 format (i2,' m',2x)
110 format (f4.1,' m')

  !--------- initialize static diagnostic fields --------------------

  id_land_mask = &
       register_static_field ( mod_name, 'land_mask', atmos_axes,  &
       'fractional amount of land', 'none', &
       range=frange, interp_method = "conserve_order1" )
  
  !--------- initialize diagnostic fields --------------------

  id_ice_mask = &
       register_diag_field ( mod_name, 'ice_mask', atmos_axes, Time, &
       'fractional amount of sea ice', 'none',  &
       range=frange, interp_method = "conserve_order1" )
  
  id_wind = &
       register_diag_field ( mod_name, 'wind', atmos_axes, Time, &
       'wind speed for flux calculations', 'm/s', &
       range=(/0.,vrange(2)/) )
  
  id_drag_moist = &
       register_diag_field ( mod_name, 'drag_moist', atmos_axes, Time, &
       'drag coeff for moisture',    'none'     )
  
  id_drag_heat  = &
       register_diag_field ( mod_name, 'drag_heat', atmos_axes, Time, &
       'drag coeff for heat',    'none'     )
  
  id_drag_mom   = &
       register_diag_field ( mod_name, 'drag_mom',  atmos_axes, Time, &
       'drag coeff for momentum',     'none'     )
  
  id_rough_moist = &
       register_diag_field ( mod_name, 'rough_moist', atmos_axes, Time, &
       'surface roughness for moisture',  'm'  )

  id_rough_heat = &
       register_diag_field ( mod_name, 'rough_heat', atmos_axes, Time, &
       'surface roughness for heat',  'm'  )

  id_rough_mom  = &
       register_diag_field ( mod_name, 'rough_mom',  atmos_axes, Time, &
       'surface roughness for momentum',  'm'  )

  id_u_star     = &
       register_diag_field ( mod_name, 'u_star',     atmos_axes, Time, &
       'friction velocity',   'm/s'   )

  id_b_star     = &
       register_diag_field ( mod_name, 'b_star',     atmos_axes, Time, &
       'buoyancy scale',      'm/s2'   )

  id_q_star     = &
       register_diag_field ( mod_name, 'q_star',     atmos_axes, Time, &
       'moisture scale',      'kg water/kg air'   )

  id_u_flux     = &
       register_diag_field ( mod_name, 'tau_x',      atmos_axes, Time, &
       'zonal wind stress',     'pa'   )

  id_v_flux     = &
       register_diag_field ( mod_name, 'tau_y',      atmos_axes, Time, &
       'meridional wind stress',     'pa'   )

  id_t_surf     = &
       register_diag_field ( mod_name, 't_surf',     atmos_axes, Time, &
       'surface temperature',    'deg_k', &
       range=trange    )

  ! + slm, Mar 25, 2002 -- add diagnositcs for t_ca, q_ca, and q_atm
  id_t_ca       = &
       register_diag_field ( mod_name, 't_ca',     atmos_axes, Time, &
       'canopy air temperature',    'deg_k', &
       range=trange    )

  ! - slm, Mar 25, 2002
  id_z_atm      = &
       register_diag_field ( mod_name, 'z_atm',     atmos_axes, Time, &
       'height of btm level',    'm')

  id_p_atm      = &
       register_diag_field ( mod_name, 'p_atm',     atmos_axes, Time, &
       'pressure at btm level',    'pa')

  ! - bw, Mar 25, 2002 -- added diagnostic slp
  id_slp      = &
       register_diag_field ( mod_name, 'slp',      atmos_axes, Time, &
       'sea level pressure',    'pa')

  id_gust       = &
       register_diag_field ( mod_name, 'gust',     atmos_axes, Time, &
       'gust scale',    'm/s')

  id_t_flux     = &
       register_diag_field ( mod_name, 'shflx',      atmos_axes, Time, &
       'sensible heat flux',     'w/m2'    )

  id_r_flux     = &
       register_diag_field ( mod_name, 'lwflx',      atmos_axes, Time, &
       'net (down-up) longwave flux',   'w/m2'    )

  id_t_atm      = &
       register_diag_field ( mod_name, 't_atm',      atmos_axes, Time, &
       'temperature at btm level',    'deg_k', &
       range=trange     )

  id_u_atm      = &
       register_diag_field ( mod_name, 'u_atm',      atmos_axes, Time, &
       'u wind component at btm level',  'm/s', &
       range=vrange    )

  id_v_atm      = &
       register_diag_field ( mod_name, 'v_atm',      atmos_axes, Time, &
       'v wind component at btm level',  'm/s', &
       range=vrange    )

  id_t_ref      = &
       register_diag_field ( mod_name, 't_ref',      atmos_axes, Time, &
       'temperature at '//label_zh, 'deg_k' , &
       range=trange      )

  id_rh_ref     = &
       register_diag_field ( mod_name, 'rh_ref',     atmos_axes, Time,   &
       'relative humidity at '//label_zh, 'percent' )

  id_rh_ref_cmip = &
       register_diag_field ( mod_name, 'rh_ref_cmip',     atmos_axes, Time,   &
       'relative humidity at '//label_zh, 'percent' )

  id_u_ref      = &
       register_diag_field ( mod_name, 'u_ref',      atmos_axes, Time, &
       'zonal wind component at '//label_zm,  'm/s', &
       range=vrange )

  id_v_ref      = &
       register_diag_field ( mod_name, 'v_ref',      atmos_axes, Time,     &
       'meridional wind component at '//label_zm, 'm/s', &
       range=vrange )

  id_wind_ref = &
       register_diag_field ( mod_name, 'wind_ref',   atmos_axes, Time,     &
       'absolute value of wind at '//label_zm, 'm/s', &
       range=vrange )

  id_del_h      = &
       register_diag_field ( mod_name, 'del_h',      atmos_axes, Time,  &
       'ref height interp factor for heat', 'none' )
  id_del_m      = &
       register_diag_field ( mod_name, 'del_m',      atmos_axes, Time,     &
       'ref height interp factor for momentum','none' )
  id_del_q      = &
       register_diag_field ( mod_name, 'del_q',      atmos_axes, Time,     &
       'ref height interp factor for moisture','none' )

  ! + slm Jun 02, 2002 -- diagnostics of reference values over the land
  id_t_ref_land = &
       register_diag_field ( mod_name, 't_ref_land', Land_axes, Time, &
       'temperature at '//trim(label_zh)//' over land', 'deg_k' , &
       range=trange, missing_value =  -100.0)
  id_rh_ref_land= &
       register_diag_field ( mod_name, 'rh_ref_land', Land_axes, Time,   &
       'relative humidity at '//trim(label_zh)//' over land', 'percent',       &
       missing_value=-999.0)
  id_u_ref_land = &
       register_diag_field ( mod_name, 'u_ref_land',  Land_axes, Time, &
       'zonal wind component at '//trim(label_zm)//' over land',  'm/s', &
       range=vrange, missing_value=-999.0 )
  id_v_ref_land = &
       register_diag_field ( mod_name, 'v_ref_land',  Land_axes, Time,     &
       'meridional wind component at '//trim(label_zm)//' over land', 'm/s', &
       range=vrange, missing_value = -999.0 )
  ! - slm Jun 02, 2002
  id_q_ref = &
       register_diag_field ( mod_name, 'q_ref', atmos_axes, Time,     &
       'specific humidity at '//trim(label_zh), 'kg/kg', missing_value=-1.0)
  id_q_ref_land = &
       register_diag_field ( mod_name, 'q_ref_land', Land_axes, Time, &
       'specific humidity at '//trim(label_zh)//' over land', 'kg/kg',          &
       missing_value=-1.0)

  id_rough_scale = &
       register_diag_field ( mod_name, 'rough_scale', atmos_axes, Time, &
       'topographic scaling factor for momentum drag','1' )
!-----------------------------------------------------------------------

  allocate(id_tr_atm(n_exch_tr))
  allocate(id_tr_surf(n_exch_tr))
  allocate(id_tr_flux(n_exch_tr))
  allocate(id_tr_mol_flux(n_exch_tr))

  do tr = 1, n_exch_tr
     call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, name, longname, units )
     id_tr_atm(tr) = register_diag_field (mod_name, trim(name)//'_atm', atmos_axes, Time, &
          trim(longname)//' at btm level', trim(units))
     id_tr_surf(tr) = register_diag_field (mod_name, trim(name)//'_surf', atmos_axes, Time, &
          trim(longname)//' at the surface', trim(units))
     id_tr_flux(tr) = register_diag_field(mod_name, trim(name)//'_flux', atmos_axes, Time, &
          'flux of '//trim(longname), trim(units)//' kg air/(m2 s)')
!! add dryvmr co2_surf and co2_atm
     if ( lowercase(trim(name))=='co2') then
! - slm Mar 25, 2010: moved registration of mol_flux inside 'if' to disable 
! saving incorrect results (mol fluxes for other tracers computed with CO2 molar 
! mass)
       id_tr_mol_flux(tr) = register_diag_field(mod_name, trim(name)//'_mol_flux', atmos_axes, Time, &
            'flux of '//trim(longname), 'mol CO2/(m2 s)')
       id_co2_atm_dvmr = register_diag_field (mod_name, trim(name)//'_atm_dvmr', atmos_axes, Time, &
            trim(longname)//' at btm level', 'mol CO2 /mol air')
       id_co2_surf_dvmr = register_diag_field (mod_name, trim(name)//'_surf_dvmr', atmos_axes, Time, &
            trim(longname)//' at the surface', 'mol CO2 /mol air')
     else
       id_tr_mol_flux(tr) = -1
     endif
  enddo

  id_q_flux = register_diag_field( mod_name, 'evap',       atmos_axes, Time, &
         'evaporation rate',        'kg/m2/s'  )
  id_q_flux_land = register_diag_field( mod_name, 'evap_land', land_axes, Time, &
         'evaporation rate over land',        'kg/m2/s', missing_value=-1.0 )

  end subroutine diag_field_init

!#######################################################################
  subroutine flux_ice_to_ocean_redistribute(ice, ocean, ice_data, ocn_bnd_data, type, do_area_weighted )

    ! Performs a globally conservative flux redistribution across ICE/OCN.
    ! Assumes that the ice/ocn grids are the same. If ocean is present,
    ! then assume different mpp domans and redistribute

    ! should be invoked by all PEs

    type(ice_data_type),              intent(in) :: ice
    type(ocean_public_type),          intent(in) :: ocean
    real, dimension(:,:),             intent(in) :: ice_data
    real, dimension(:,:),            intent(out) :: ocn_bnd_data
    integer,                          intent(in) :: type
    logical,                          intent(in) :: do_area_weighted


    real :: tmp(nxc_ice, nyc_ice )

    select case(type)
    case(DIRECT)
       if(do_area_weighted) then
          ocn_bnd_data = ice_data * ice%area
          call divide_by_area(data=ocn_bnd_data, area=ocean%area)          
       else
          ocn_bnd_data = ice_data
       endif
    case(REDIST)
       if(do_area_weighted) then
          if( ice%pe ) tmp = ice_data  * ice%area
          call mpp_redistribute(ice%Domain, tmp, ocean%Domain, ocn_bnd_data)
          if(ocean%is_ocean_pe) call divide_by_area(ocn_bnd_data, area=ocean%area) 
       else
          call mpp_redistribute(ice%Domain, ice_data, ocean%Domain, ocn_bnd_data)
       endif
    case DEFAULT
       call mpp_error( FATAL, 'FLUX_ICE_TO_OCEAN: Ice_Ocean_Boundary%xtype must be DIRECT or REDIST.' )
    end select

  end subroutine flux_ice_to_ocean_redistribute

!######################################################################################
! Divide data by area while avoiding zero area elements
  subroutine divide_by_area(data, area)
    real, intent(inout) :: data(:,:)
    real, intent(in)    :: area(:,:)

    if(size(data, dim=1) /= size(area, dim=1) .or. size(data, dim=2) /= size(area, dim=2)) then
       ! no op
       return
    endif

    where(area /= 0) 
       data = data / area
    end where

  end subroutine divide_by_area
!#######################################################################

! This private routine will check flux conservation for routine flux_ice_to_ocean_redistribute
! when do_area_weighted_flux = false and true. 
  subroutine check_flux_conservation(Ice, Ocean, Ice_Ocean_Boundary)
  type(ice_data_type),               intent(inout)  :: Ice
  type(ocean_public_type),           intent(inout)  :: Ocean
  type(ice_ocean_boundary_type),     intent(inout) :: ice_ocean_boundary

  real, allocatable, dimension(:,:) :: ice_data, ocn_data
  real :: ice_sum, area_weighted_sum, non_area_weighted_sum
  integer :: outunit

  outunit = stdout()
  allocate(ice_data(size(Ice%flux_q,1), size(Ice%flux_q,2) ) )
  allocate(ocn_data(size(Ice_Ocean_Boundary%q_flux,1), size(Ice_Ocean_Boundary%q_flux,2) ) )
  call random_number(ice_data)
  ice_sum = sum(ice_data*ice%area)
  call mpp_sum(ice_sum)
  ocn_data = 0
  call flux_ice_to_ocean_redistribute( Ice, Ocean, ice_data, ocn_data, Ice_Ocean_Boundary%xtype, .false.)
  non_area_weighted_sum = sum(ocn_data*ocean%area)
  call mpp_sum(non_area_weighted_sum)
  ocn_data = 0
  call flux_ice_to_ocean_redistribute( Ice, Ocean, ice_data, ocn_data, Ice_Ocean_Boundary%xtype, .true.)
  area_weighted_sum = sum(ocn_data*ocean%area)
  call mpp_sum(area_weighted_sum)  
  write(outunit,*)"NOTE from flux_exchange_mod: check for flux conservation for flux_ice_to_ocean"
  write(outunit,*)"***** The global area sum of random number on ice domain (input data) is ", ice_sum
  write(outunit,*)"***** The global area sum of data after flux_ice_to_ocean_redistribute with "// &
       "do_area_weighted_flux = false is ", non_area_weighted_sum, &
       " and the difference from global input area sum = ", ice_sum - non_area_weighted_sum
  write(outunit,*)"***** The global area sum of data after flux_ice_to_ocean_redistribute with "// &
       "do_area_weighted_flux = true is ", area_weighted_sum, &
       " and the difference from global input area sum = ", ice_sum - area_weighted_sum


  end subroutine check_flux_conservation

! <DIAGFIELDS>
!   <NETCDF NAME="land_mask" UNITS="none">
!     fractional amount of land
!   </NETCDF>
!   <NETCDF NAME="wind" UNITS="m/s">
!     wind speed for flux calculations
!   </NETCDF>
!   <NETCDF NAME="drag_moist" UNITS="none">
!     drag coeff for moisture
!   </NETCDF>
!   <NETCDF NAME="drag_heat" UNITS="none">
!     drag coeff for heat
!   </NETCDF>
!   <NETCDF NAME="drag_mom" UNITS="none">
!     drag coeff for momentum
!   </NETCDF>
!   <NETCDF NAME="rough_moist" UNITS="m">
!     surface roughness for moisture
!   </NETCDF>
!   <NETCDF NAME="rough_heat" UNITS="m">
!     surface roughness for heat
!   </NETCDF>
!   <NETCDF NAME="rough_mom" UNITS="m">
!     surface roughness for momentum
!   </NETCDF>
!   <NETCDF NAME="u_star" UNITS="m/s">
!     friction velocity
!   </NETCDF>
!   <NETCDF NAME="b_star" UNITS="m/s">
!     buoyancy scale
!   </NETCDF>
!   <NETCDF NAME="q_star" UNITS="kg water/kg air">
!     moisture scale
!   </NETCDF>
!   <NETCDF NAME="t_atm" UNITS="deg_k">
!     temperature at btm level
!   </NETCDF>
!   <NETCDF NAME="u_atm" UNITS="m/s">
!     u wind component at btm level
!   </NETCDF>
!   <NETCDF NAME="v_atm" UNITS="m/s">
!     v wind component at btm level
!   </NETCDF>
!   <NETCDF NAME="q_atm" UNITS="kg/kg">
!     specific humidity at btm level
!   </NETCDF>
!   <NETCDF NAME="p_atm" UNITS="pa">
!     pressure at btm level
!   </NETCDF>
!   <NETCDF NAME="z_atm" UNITS="m">
!     height of btm level
!   </NETCDF>
!   <NETCDF NAME="gust" UNITS="m/s">
!     gust scale 
!   </NETCDF>
!   <NETCDF NAME="rh_ref" UNITS="percent">
!     relative humidity at ref height
!   </NETCDF>
!   <NETCDF NAME="t_ref" UNITS="deg_k">
!    temperature at ref height
!   </NETCDF>
!   <NETCDF NAME="u_ref" UNITS="m/s">
!    zonal wind component at ref height
!   </NETCDF>
!   <NETCDF NAME="v_ref" UNITS="m/s">
!    meridional wind component at ref height 
!   </NETCDF>
!   <NETCDF NAME="del_h" UNITS="none">
!    ref height interp factor for heat 
!   </NETCDF>
!   <NETCDF NAME="del_m" UNITS="none">
!    ref height interp factor for momentum 
!   </NETCDF>
!   <NETCDF NAME="del_q" UNITS="none">
!    ref height interp factor for moisture
!   </NETCDF>
!   <NETCDF NAME="tau_x" UNITS="pa">
!    zonal wind stress
!   </NETCDF>
!   <NETCDF NAME="tau_y" UNITS="pa">
!    meridional wind stress
!   </NETCDF>
!   <NETCDF NAME="ice_mask" UNITS="none">
!    fractional amount of sea ice 
!   </NETCDF>
!   <NETCDF NAME="t_surf" UNITS="deg_k">
!     surface temperature
!   </NETCDF>
!   <NETCDF NAME="t_ca" UNITS="deg_k">
!     canopy air temperature
!   </NETCDF>
!   <NETCDF NAME="q_surf" UNITS="kg/kg">
!     surface specific humidity 
!   </NETCDF>
!   <NETCDF NAME="shflx" UNITS="w/m2">
!     sensible heat flux
!   </NETCDF>
!   <NETCDF NAME="evap" UNITS="kg/m2/s">
!     evaporation rate 
!   </NETCDF>
!   <NETCDF NAME="lwflx" UNITS="w/m2">
!    net (down-up) longwave flux 
!   </NETCDF>

! </DIAGFIELDS>

! <INFO>


!   <NOTE>
!   <PRE>
!
!  MAIN PROGRAM EXAMPLE
!  --------------------
!
!       DO slow time steps (ocean)
!
!           call flux_ocean_to_ice
!
!           call ICE_SLOW_UP
!
!
!           DO fast time steps (atmos)
!
!                call sfc_boundary_layer
!
!                call ATMOS_DOWN
!
!                call flux_down_from_atmos
!
!                call LAND_FAST
!
!                call ICE_FAST
!
!                call flux_up_to_atmos
!
!                call ATMOS_UP
!
!           END DO
!
!           call ICE_SLOW_DN
!
!           call flux_ice_to_ocean
!
!           call OCEAN
!
!      END DO
!
!   LAND_FAST and ICE_FAST must update the surface temperature
!
! =======================================================================
!
! REQUIRED VARIABLES IN DEFINED DATA TYPES FOR COMPONENT MODELS
! --------------------------------------------------------------
!
! type (atmos_boundary_data_type) :: Atm
! type (surf_diff_type) :: Atm%Surf_Diff
!
! real, dimension(:)
!
!    Atm%lon_bnd   longitude axis grid box boundaries in radians
!                  must be monotonic
!    Atm%lat_bnd   latitude axis grid box boundaries in radians
!                  must be monotonic
!
! real, dimension(:,:)
!
!    Atm%t_bot     temperature at lowest model level
!    Atm%q_bot     specific humidity at lowest model level
!    Atm%z_bot     height above the surface for the lowest model level (m)
!    Atm%p_bot     pressure at lowest model level (pa)
!    Atm%u_bot     zonal wind component at lowest model level (m/s)
!    Atm%v_bot     meridional wind component at lowest model level (m/s)
!    Atm%p_surf    surface pressure (pa)
!    Atm%slp       sea level pressure (pa)
!    Atm%gust      gustiness factor (m/s)
!    Atm%flux_sw   net shortwave flux at the surface
!    Atm%flux_lw   downward longwave flux at the surface
!    Atm%lprec     liquid precipitation (kg/m2)
!    Atm%fprec     water equivalent frozen precipitation (kg/m2)
!    Atm%coszen    cosine of the zenith angle
!
!   (the following five fields are gathered into a data type for convenience in passing
!   this information through the different levels of the atmospheric model --
!   these fields are rlated to the simultaneous implicit time steps in the
!   atmosphere and surface models -- they are described more fully in
!   flux_exchange.tech.ps and
!   in the documntation for vert_diff_mod
!
!
!    Atm%Surf_Diff%dtmass   = dt/mass where dt = atmospheric time step ((i+1) = (i-1) for leapfrog) (s)
!                           mass = mass per unit area of lowest atmosphehic layer  (Kg/m2))
!    Atm%Surf_Diff%delta_t  increment ((i+1) = (i-1) for leapfrog) in temperature of
!                           lowest atmospheric layer  (K)
!    Atm%Surf_Diff%delta_q  increment ((i+1) = (i-1) for leapfrog) in specific humidity of
!                           lowest atmospheric layer (nondimensional -- Kg/Kg)
!    Atm%Surf_Diff%dflux_t  derivative of implicit part of downward temperature flux at top of lowest
!                           atmospheric layer with respect to temperature
!                           of lowest atmospheric layer (Kg/(m2 s))
!    Atm%Surf_Diff%dflux_q  derivative of implicit part of downward moisture flux at top of lowest
!                           atmospheric layer with respect to specific humidity of
!                           of lowest atmospheric layer (Kg/(m2 s))
!
!
! integer, dimension(4)
!
!    Atm%axes      Axis identifiers returned by diag_axis_init for the
!                  atmospheric model axes: X, Y, Z_full, Z_half.
!
! -----------------------------------------------
!
! type (land_boundary_data_type) :: Land
!
! real, dimension(:)
!
!    Land%lon_bnd     longitude axis grid box boundaries in radians
!                     must be monotonic
!    Land%lat_bnd     latitude axis grid box boundaries in radians
!                     must be monotonic
!
! logical, dimension(:,:,:)
!
!    Land%mask        land/sea mask (true for land)
!    Land%glacier     glacier mask  (true for glacier)
!
! real, dimension(:,:,:)
!
!    Land%tile_size   fractional area of each tile (partition)
!
!    Land%t_surf      surface temperature (deg k)
!    Land%albedo      surface albedo (fraction)
!    Land%rough_mom   surface roughness for momentum (m)
!    Land%rough_heat  surface roughness for heat/moisture (m)
!    Land%stomatal    stomatal resistance
!    Land%snow        snow depth (water equivalent) (kg/m2)
!    Land%water       water depth of the uppermost bucket (kg/m2)
!    Land%max_water   maximum water depth allowed in the uppermost bucket (kg/m2)
!
! -----------------------------------------------
!
!
! type (ice_boundary_data_type) :: Ice
!
! real, dimension(:)
!
!    Ice%lon_bnd       longitude axis grid box boundaries for temperature points
!                      in radians (must be monotonic)
!    Ice%lat_bnd       latitude axis grid box boundaries for temperature points
!                      in radians (must be monotonic)
!    Ice%lon_bnd_uv    longitude axis grid box boundaries for momentum points
!                      in radians (must be monotonic)
!    Ice%lat_bnd_uv    latitude axis grid box boundaries for momentum points
!                      in radians (must be monotonic)
!
! logical, dimension(:,:,:)
!
!    Ice%mask          ocean/land mask for temperature points
!                        (true for ocean, with or without ice)
!    Ice%mask_uv       ocean/land mask for momentum points
!                        (true for ocean, with or without ice)
!    Ice%ice_mask      optional ice mask (true for ice)
!
! real, dimension(:,:,:)
!
!    Ice%part_size     fractional area of each partition of a temperature grid box
!    Ice%part_size_uv  fractional area of each partition of a momentum grid box
!
!    the following fields are located on the ice top grid
!
!    Ice%t_surf        surface temperature (deg k)
!    Ice%albedo        surface albedo (fraction)
!    Ice%rough_mom     surface roughness for momentum (m)
!    Ice%rough_heat    surface roughness for heat/moisture (m)
!    Ice%u_surf        zonal (ocean/ice) current at the surface (m/s)
!    Ice%v_surf        meridional (ocean/ice) current at the surface (m/s)
!
!    the following fields are located on the ice bottom grid
!
!    Ice%flux_u        zonal wind stress (Pa)
!    Ice%flux_v        meridional wind stress (Pa)
!    Ice%flux_t        sensible heat flux (w/m2)
!    Ice%flux_q        specific humidity flux (kg/m2/s)
!    Ice%flux_sw       net (down-up) shortwave flux (w/m2)
!    Ice%flux_lw       net (down-up) longwave flux (w/m2)
!    Ice%lprec         mass of liquid precipitation since last time step (Kg/m2)
!    Ice%fprec         mass of frozen precipitation since last time step (Kg/m2)
!    Ice%runoff        mass of runoff water since last time step (Kg/m2)
!
! -----------------------------------------------
!
! type (ocean_boundary_data_type) :: Ocean
!
! real, dimension(:)
!
!    Ocean%Data%lon_bnd      longitude axis grid box boundaries for temperature
!                            points on the ocean DATA GRID (radians)
!    Ocean%Data%lat_bnd      latitude axis grid box boundaries for temperature
!                            points on the ocean DATA GRID (radians)
!    Ocean%Data%lon_bnd_uv   longitude axis grid box boundaries for momentum
!                            points on the ocean DATA GRID (radians)
!    Ocean%Data%lat_bnd_uv   latitude axis grid box boundaries for momentum
!                            points on the ocean DATA GRID (radians)
!
!    Ocean%Ocean%lon_bnd     longitude axis grid box boundaries for temperature
!                            points on the ocean MODEL GRID (radians)
!    Ocean%Ocean%lat_bnd     latitude axis grid box boundaries for temperature
!                            points on the ocean MODEL GRID (radians)
!    Ocean%Ocean%lon_bnd_uv  longitude axis grid box boundaries for momentum
!                            points on the ocean MODEL GRID (radians)
!    Ocean%Ocean%lat_bnd_uv  latitude axis grid box boundaries for momentum
!                            points on the ocean MODEL GRID (radians)
!
!      Note: The data values in all longitude and latitude grid box boundary
!            array must be monotonic.
!
! logical, dimension(:,:)
!
!    Ocean%Data%mask       ocean/land mask for temperature points on the ocean
!                          DATA GRID (true for ocean)
!    Ocean%Data%mask_uv    ocean/land mask for momentum points on the ocean
!                          DATA GRID (true for ocean)
!
!    Ocean%Ocean%mask      ocean/land mask for temperature points on the ocean
!                          MODEL GRID (true for ocean)
!    Ocean%Ocean%mask_uv   ocean/land mask for momentum points on the ocean
!                          MODEL GRID (true for ocean)
!
! real, dimension(:,:)
!
!    Ocean%t_surf_data  surface temperature on the ocean DATA GRID (deg k)
!
!    Ocean%t_surf       surface temperature on the ocean MODEL GRID (deg k)
!    Ocean%u_surf       zonal ocean current at the surface on the ocean
!                       MODEL GRID (m/s)
!    Ocean%v_surf       meridional ocean current at the surface on the
!                       ocean MODEL GRID (m/s)
!    Ocean%frazil       frazil at temperature points on the ocean MODEL GRID
!
!   </PRE>
!   </NOTE>
! </INFO>

end module flux_exchange_mod


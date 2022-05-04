module ocean_sbc_mod
!  
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
! </CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> M.J. Harrison
!</REVIEWER>
!
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> A. Rosati 
! </REVIEWER>
!
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
! V. Balaji
! </REVIEWER>
!
!<OVERVIEW>
! Set up the surface boundary conditions for MOM. 
!</OVERVIEW>
!
!<DESCRIPTION>
!
! This module sets up the surface boundary conditions for MOM. 
! Also fill Ocean_sfc derived-type used to pass information to other 
! component models.  Also write diagnostics related to surface
! boundary forcing.  
!
! The surface temperature should be the surface insitu temperature,
! which is the same as the surface potential temperature.  When the 
! model prognostic temperature variable is conservative temperature, 
! then the surface potential temperature is carried in T_diag(index_diag_temp).
! The resulting heat flux is potential enthalpy, which is the correct 
! field to be forcing the T_prog(index_temp) field when the prognostic
! temperature field is the conservative temperature.   
!
! We assume the winds passed to the ocean are on the B-grid 
! velocity point.  Likewise, we pass the currents back to the coupler 
! on the B-grid point.  Code will need to be modified if using another
! assumption.  
!
! Treatment of flux adjustments may be modified according to the FAFMIP
! experiment protocol (Gregory et al 2016). See also the design notes by
! Griffies et al http://www.met.reading.ac.uk/~jonathan/FAFMIP/GFDL_heat.pdf
! and the FAFMIP website http://www.met.reading.ac.uk/~jonathan/FAFMIP/
!
!</DESCRIPTION>
!
! <REFERENCE>
! Gregory, J. M., Bouttes, N., Griffies, S. M., Haak, H., Hurlin, W. J.,
! Jungclaus, J., Kelley, M., Lee, W. G., Marshall, J., Romanou, A., Saenko, O.
! A., Stammer, D., and Winton, M.: The Flux-Anomaly-Forced Model Intercomparison
! Project (FAFMIP) contribution to CMIP6: investigation of sea-level and ocean
! climate change in response to CO2 forcing, Geosci. Model Dev., 9, 3993-4017,
! https://doi.org/10.5194/gmd-9-3993-2016, 2016.
! </REFERENCE>
!
!<NAMELIST NAME="ocean_sbc_nml">
!
!  <DATA NAME="use_waterflux" TYPE="logical">
!  Set to true when wish to use real fresh water flux as opposed to virtual 
!  salt fluxes. This is the recommended method. The alternative virtual 
!  tracer flux method (use_waterflux=.false.) is not routinely used at 
!  GFDL, so it may suffer from poor testing.  
!  Default use_waterflux=.true. 
!  </DATA> 
!  <DATA NAME="waterflux_tavg" TYPE="logical">
!  Set to true when aiming to suppress the leap-frog computational mode
!  by setting pme and river equal to a time averaged value over the 
!  present and previous time step.  This method requires an extra
!  field in the restart file.  This method is not needed when using
!  the TWO_LEVEL time tendency.  It remains for those who wish to 
!  use the leap-frog THREE_LEVEL time stepping scheme.  
!  Note that it does not lead to simple checks of conservation across
!  model components, since there is a time averaging performed for 
!  the water flux added to the ocean model.  It is generally NOT 
!  recommended.  Default waterflux_tavg=.false. 
!  </DATA> 
!
!  <DATA NAME="use_ideal_runoff" TYPE="logical">
!  To add an idealized liquid runoff read from a file.
!  This runoff is assumed to enter the ocean with the same temperature as
!  SST, and to be liquid.  It is an additional runoff, so that any other
!  runoff remains unaltered.  The runoff coming from idealized runoff is
!  NOT subject to the global normalization realized from
!  zero_net_water_coupler=.true.   
!  Default use_ideal_runoff=.false.
!  </DATA> 
!  <DATA NAME="use_ideal_calving" TYPE="logical">
!  To add an idealized solid runoff or calving read from a file.
!  This calving runoff is assumed to require melting, so it extracts
!  latent heat of fusion from the liquid ocean
!  The runoff coming from idealized cavling is NOT subject to
!  the global normalization realized from zero_net_water_coupler=.true.  
!  Default use_ideal_calving=.false. 
!  </DATA> 
!
!  <DATA NAME="use_waterflux_override_calving" TYPE="logical">
!  Set to true will allow for model to incorporate the latent heating
!  from a calving field that comes in through coupled model instantaneous 
!  interactions, but later will over-ride the mass flux from calving with 
!  a dataset that is read in from a climatology or observations.    
!  The idea is to only modify the mass contribution from calving through 
!  the over-ride, and leave the latent heat contribution untouched.  
!  Default use_waterflux_override_calving=.false. 
!  </DATA> 
!  <DATA NAME="use_waterflux_override_fprec" TYPE="logical">
!  Set to true will allow for model to incorporate the latent heating
!  from a fprec field that comes in through coupled model instantaneous 
!  interactions, but later will over-ride the mass flux from fprec with 
!  a dataset that is read in from a climatology or observations.    
!  The idea is to only modify the mass contribution from fprec through 
!  the over-ride, and leave the latent heat contribution untouched.  
!  Default use_waterflux_override_fprec=.false. 
!  </DATA> 
!  <DATA NAME="use_waterflux_override_evap" TYPE="logical">
!  Set to true will allow for model to incorporate the latent heating
!  from an evap field that comes in through coupled model instantaneous 
!  interactions, but later will over-ride the mass flux from evap with 
!  a dataset that is read in from a climatology or observations.    
!  The idea is to only modify the mass contribution from evap through 
!  the over-ride, and leave the latent heat contribution untouched.  
!  Default use_waterflux_override_evap=.false. 
!  </DATA> 
!
!  <DATA NAME="temp_restore_tscale" UNITS="day" TYPE="real">
!  Time scale in days for restoring temperature within the top model 
!  grid cell. 
!  </DATA> 
!  <DATA NAME="salt_restore_tscale" UNITS="day" TYPE="real">
!  Time scale in days for restoring salinity within the top model 
!  grid cell. 
!  </DATA> 
!  <DATA NAME="eta_restore_tscale" UNITS="day" TYPE="real">
!  Time scale in days for restoring surface height to produce a modification to 
!  surface water flux.  This option is only available when run with 
!  use_waterflux=.true.   
!  </DATA> 
!  <DATA NAME="use_constant_sst_for_restore" TYPE="logical">
!  To over-ride the sfc_restore.nc value for temp restoring.
!  use_constant_sst_for_restore=.false. 
!  </DATA> 
!  <DATA NAME="constant_sst_for_restore" UNITS="degC"  TYPE="real">
!  The SST value used if use_constant_sst_for_restore=.true.
!  Default constant_sst_for_restore=12.0
!  </DATA> 
!  <DATA NAME="salt_restore_as_salt_flux" TYPE="logical">
!  When running a use_waterflux=.true. model, we may choose to add the 
!  salinity from a restoring condition as a salt flux or convert to 
!  a fresh water flux. The addition of salt does not alter the sea 
!  level nor does it alter the concentration of other tracers, whereas
!  converting to an implied water flux will alter sea level and other
!  concentrations.  So we generally recommend the default   
!  salt_restore_as_salt_flux=.true. 
!  </DATA> 
!  <DATA NAME="use_constant_sss_for_restore" TYPE="logical">
!  To over-ride the sfc_restore.nc value for salinity restoring.
!  use_constant_sss_for_restore=.false. 
!  </DATA> 
!  <DATA NAME="constant_sss_for_restore" UNITS="psu"  TYPE="real">
!  The SSS value used if use_constant_sss_for_restore=.true.
!  Default constant_sss_for_restore=35.0
!  </DATA> 
!  <DATA NAME="max_delta_salinity_restore" UNITS="ppt" TYPE="real">
!  When computing the restoring flux for salinity, we can define
!  a maximum absolute value for the difference between salinity(k=1)
!  and the restoring salinity from a dataset.  This approach is useful
!  especially in NAtl western boundary, where poor Gulf Stream separation
!  can lead to large salinity biases.  If restore too much the salinity
!  field, we can spuriously transport large amounts of fresh water to the 
!  subpoloar gyre, thus impacting the overturning circulation too much.  
!  If max_delta_salinity_restore < 0.0, then will NOT provide a max to the 
!  delta salinity; will instead compute an unbounded restoring flux.  
!  Default max_delta_salinity_restore=-0.50.
!  </DATA> 
!  <DATA NAME="salinity_restore_limit_lower" UNITS="ppt" TYPE="real">
!  Define a lower absolute salinity below which max_delta_salinity_restore
!  is ignored. In some coastal regions with high freshwater runoff salinity
!  can drop below zero. Set this limit to a positive value below which
!  restoring will be more aggressive.
!  Default max_delta_salinity_restore=0.0.
!  </DATA>
!  <DATA NAME="salinity_restore_limit_upper" UNITS="ppt" TYPE="real">
!  Define an upper absolute salinity above which max_delta_salinity_restore
!  is ignored. In some enclosed coastal regions with low runoff and high
!  evaporation salinity can reach physically unrealistic levels. Set this
!  limit to a positive value above which restoring will be more aggressive.
!  Default max_delta_salinity_restore=100.0.
!  </DATA>
!
!  <DATA NAME="read_restore_mask" TYPE="logical">
!  For reading in a mask that selects regions of the domain 
!  that are restored (mask=1) or not restored (mask=0).
!  Default  read_restore_mask=.false., whereby restore_mask
!  is set to tmask(k=1). 
!  </DATA> 
!  <DATA NAME="restore_mask_gfdl" TYPE="logical">
!  For modifying the restore mask based on reading in 
!  the GFDL regional mask. Default restore_mask_gfdl=.false.
!  </DATA> 
!  <DATA NAME="salinity_ref" UNITS="psu" TYPE="real">
!  Reference salinity used for converting fresh water flux
!  to salt flux. 
!  </DATA> 
!  <DATA NAME="salt_restore_under_ice" TYPE="logical">
!  Logical indicating whether to restore salinity under sea ice or not.
!  When .false. then will not restore salinity  in regions where we 
!  use a "frazil" condition as a proxy for where sea-ice is present.
!  Do not use sea ice extent from a sea ice model since we generally do 
!  not pass information regarding ice extent between the sea ice model 
!  and the ocean model.     
!  </DATA> 
!  <DATA NAME="zero_net_salt_restore" TYPE="logical">
!  Logical indicating whether to remove the area mean of the salinity 
!  restore flux so there is a net zero input of salt to the ocean
!  associated with restoring.    
!  </DATA> 
!  <DATA NAME="zero_net_salt_correction" TYPE="logical">
!  Logical indicating whether to remove the area mean of the salinity 
!  correction flux so there is a net zero input of salt to the ocean
!  associated with salt correction. 
!  </DATA> 
!  <DATA NAME="zero_net_water_restore" TYPE="logical">
!  Logical indicating whether to remove the area mean of the water 
!  restore flux so there is a net zero input of water to the ocean
!  associated with restoring.    
!  </DATA> 
!  <DATA NAME="zero_net_water_correction" TYPE="logical">
!  Logical indicating whether to remove the area mean of the water 
!  correction flux so there is a net zero input of water to the ocean
!  associated with water correction.
!  </DATA> 
!  <DATA NAME="zero_net_water_coupler" TYPE="logical">
!  Logical indicating whether to remove the area mean of the water 
!  passed through the coupler so there is a net zero input of 
!  fresh water to the ocean associated with p-e+r. Do so by removing 
!  area mean from pme--keep river values unchanged. Note that a choice
!  must be made whether to remove the area mean from rivers or pme.  
!  We choose pme since it is more evenly distributed than rivers.   
!  Also note that we DO NOT include the ice melt in this normalization.
!  The reason is that we only wish to ensure the ocean+ice system
!  has a zero net water.  When melt or form sea ice, this only transfers
!  water between liquid ocean and solid sea ice, and no normalization is
!  appropriate for this case. It is only the water exchanged with the
!  land and atmosphere that is normalized. 
!  </DATA> 
!  <DATA NAME="zero_net_water_couple_restore" TYPE="logical">
!  This logical keeps the total water forcing on the ocean+ice system
!  to a global mean of zero at each time step.  We DO NOT include
!  the ice melt in this normalization.  
!  Setting zero_net_water_couple_restore to true may be appropriate when 
!  running an ice-ocean model using a bulk formulae to compute
!  evaporation (e.g., CORE) and when only providing a weak (or zero)
!  salinity restoring.  It is not appropriate when running a coupled
!  ocean-atmosphere model, where the moisture budget should be 
!  conserved without an artificial removal of the global mean.  
!  </DATA> 
!
!  <DATA NAME="land_model_heat_fluxes" TYPE="logical">
!  For the case where land model passes through the coupler the heat flux 
!  associated with the liquid runoff and calving land ice fields.
!  This heat flux is computed relative to 0C, and takes the form 
!  heat flux = mass flux of water * temp of water * heat capacity, 
!  where the water can be either liquid or solid.  For many coupled models,
!  the water temperature is assumed to be that of the SST.  But 
!  more complete land models now carry the heat of its water relative to 0C,
!  in which case the ocean model does not need to assume anything about the 
!  heat content of the land water. 
!  Default land_model_heat_fluxes=.false.      
!  </DATA> 
!
!  <DATA NAME="debug_water_fluxes" TYPE="logical">
!  Logical for debugging water fluxes. Must be true for any of the 
!  options zero_water_fluxes, zero_calving_fluxes, zero_pme_fluxes
!  or zero_runoff_fluxes to be enabled.  
!  Default debug_water_fluxes=.false.      
!  </DATA> 
!  <DATA NAME="zero_water_fluxes" TYPE="logical">
!  Logical for debugging to zero the pme, river, and pme_taum1 into 
!  ocean, over-riding any input from Ice_ocean_boundary. 
!  Default zero_water_fluxes=.false.      
!  </DATA> 
!  <DATA NAME="zero_calving_fluxes" TYPE="logical">
!  Logical for debugging to zero the calving flux passed into the ocean.
!  Default zero_calving_fluxes=.false.      
!  </DATA> 
!  <DATA NAME="zero_pme_fluxes" TYPE="logical">
!  Logical for debugging to zero the pme flux passed into the ocean.
!  Default zero_pme_fluxes=.false.      
!  </DATA> 
!  <DATA NAME="zero_runoff_fluxes" TYPE="logical">
!  Logical for debugging to zero the runoff flux passed into the ocean.
!  Default zero_runoff_fluxes=.false.      
!  </DATA> 
!  <DATA NAME="zero_river_fluxes" TYPE="logical">
!  Logical for debugging to zero the river (calving+runoff) flux passed into the ocean.
!  Default zero_river_fluxes=.false.      
!  </DATA> 
!  <DATA NAME="convert_river_to_pme" TYPE="logical">
!  Logical for debugging.  Here we add the river water input (calving+runoff)
!  to pme, then set river=calving=runoff=0.0.
!  Default convert_river_to_pme=.false.      
!  </DATA> 
!
!  <DATA NAME="sbc_heat_fluxes_const" TYPE="logical">
!  Logical for setting the surface heat flux from the coupler 
!  to a global constant. Default is sbc_heat_fluxes_const=.false.      
!  </DATA> 
!  <DATA NAME="sbc_heat_fluxes_const_seasonal" TYPE="logical">
!  Logical for setting the surface heat flux from the coupler 
!  to a global constant, and giving it a seasonally varying amplitude. 
!  Default is sbc_heat_fluxes_const_seasonal=.false.      
!  </DATA> 
!  <DATA NAME="sbc_heat_fluxes_const_value" UNITS="W/m2" TYPE="real">
!  Value for the constant heat flux when using 
!  sbc_heat_fluxes_const=.true.  
!  Default sbc_heat_fluxes_const_value=0.0.      
!  </DATA> 

!  <DATA NAME="zero_heat_fluxes" TYPE="logical">
!  Logical for debugging to set all heat fluxes into the ocean to zero, 
!  over-riding any input from Ice_ocean_boundary.  Default is .false.      
!  </DATA> 
!  <DATA NAME="zero_surface_stress" TYPE="logical">
!  Logical for debugging to zero all surface stress applied to the ocean,
!  over-riding any input from Ice_ocean_boundary.  Default is .false.      
!  </DATA> 
!
!  <DATA NAME="rotate_winds" TYPE="logical">
!  Set to true when need to rotate the winds onto the ocean model grid.
!  This is needed for cases where the winds are on a spherical grid and 
!  the ocean model uses tripolar=.true.  If generate the wind data on 
!  the ocean model grid, then do not need to rotate, since the rotation 
!  has already been done.  
!  </DATA> 
!
!  <DATA NAME="max_ice_thickness" UNITS="m" TYPE="real">
!  When coupling MOM to an ice model, the sea ice thickness may need
!  to be restricted to prevent vanishing top-level in MOM. Set 
!  max_ice_thickness (meters) < dzt(k=1) to restrict. This truncation 
!  avoids the numerical problem but we loose mass conservation in the coupled
!  sea ice and ocean system. We also alter the pressure felt on the ocean 
!  as applied by the sea ice. Different vertical coordinates are needed 
!  to do the problem more realistically.   
!
!  Note that the problem of vanishing top layer is removed when use
!  either ZSTAR or PSTAR as vertical coordinate.  
!  </DATA> 
!
!  <DATA NAME="ice_salt_concentration" UNITS="kg salt / kg ice" TYPE="real">
!  The salt concentration of sea ice.  This is taken as a bulk value, and should 
!  be the same as that used by the ice model. Default is ice_salt_concentration=0.005,
!  as that is the value used in the GFDL coupled climate model. 
!  </DATA> 
!
!  <DATA NAME="ocean_ice_salt_limit" UNITS="kg salt / kg ice" TYPE="real">
!  The minimum salt concentration of water forming sea ice.  This shouldbe at least that used by the ice model
!  otherwise it is possible to extract more salt than physically exists.
!  Default is ocean_ice_salt_limit=0.0 in order to reproduce older results.
!  It is suggested that ocean_ice_salt_limit > ice_salt_concentration +  0.001 to prevent dropping below
!  ice_salt_concentration by mistake.
!  </DATA>
!
!  <DATA NAME="runoff_salinity" UNITS="g salt / kg runoff water (ppt)" TYPE="real">
!  The salinity of river runoff water. Default is runoff_salinity=0.0.
!  </DATA> 
!  <DATA NAME="runoff_temp_min" UNITS="DegC" TYPE="real">
!  The minimum temperature that river runoff into the ocean is assigned. 
!  Default runoff_temp_min=0.0.
!  </DATA> 
!
!  <DATA NAME="runoffspread" TYPE="logical">
!  Set to true if wish to use the spread_river_horz algorithm to spread 
!  the river runoff flux horizontally over an area into the ocean wider than 
!  set by the coupler.  This option requires the setup of a table for 
!  determining the points over which we spread. 
!  Default runoffspread=.false.
!  </DATA> 
!  <DATA NAME="calvingspread" TYPE="logical">
!  Set to true if wish to use the spread_river_horz algorithm to spread 
!  the calving flux horizontally over an area into the ocean wider than 
!  set by the coupler.  This option requires the setup of a table for 
!  determining the points over which we spread. 
!  Default calvingspread=.false.
!  </DATA> 
!
!  <DATA NAME="avg_sfc_velocity" TYPE="logical">
!  If set to true, the u and v fields passed up to the sea ice
!  are averaged over a coupling interval. TRUE by default.
!  </DATA> 
!  <DATA NAME="avg_sfc_temp_salt_eta" TYPE="logical">
!  If set to true, the t, s and sea_level fields passed up to the sea ice
!  are averaged over a coupling interval. TRUE by default.
!  </DATA> 
!
!  <DATA NAME="use_full_patm_for_sea_level" TYPE="logical">
! The option use_full_patm_for_sea_level allows for the passing 
! of the sea level including the full weight of sea ice back to
! the ice model.  This approach maintains the max weight on the liquid
! ocean according to the nml variable max_ice_thickness.  But it does 
! allow the sea ice to know when there is actually more sea ice than that
! set by max_ice_thickness.  This option then provides for a negative
! feedback on the runaway growth of sea ice, since the full pressure acting to 
! make the ice flow will be correctly felt.  This is a new option, and is not
! fully tested, So the default is use_full_patm_for_sea_level=.false
!  </DATA> 
!
!  <DATA NAME="do_flux_correction" TYPE="logical">
!  For applying surface flux correction to to a tracer or wind stress field. 
!  This code is used at GFDL for idealized perturbation experiments, such 
!  as when one wishes to artificially enhance the wind stress to test 
!  model sensitivity.  It is also appropriate for coupled models that 
!  may require a modification to the fluxes arrising from a coupled model,
!  via reading in information from a pre-defined
!  data file, 
!  Default do_flux_correction=.false.
!  </DATA> 
!  <DATA NAME="temp_correction_scale" UNITS="dimensionless" TYPE="real">
!  A scale multiplying the flux correction for temperature.  
!  Default temp_correction_scale=0.0.
!  </DATA> 
!  <DATA NAME="salt_correction_scale" UNITS="dimensionless" TYPE="real">
!  A scale multiplying the flux correction for salinity.
!  Default salt_correction_scale=0.0.
!  </DATA> 
!  <DATA NAME="tau_x_correction_scale" UNITS="dimensionless" TYPE="real">
!  A scale multiplying the flux correction for tau_x.
!  Default tau_x_correction_scale=0.0.
!  </DATA> 
!  <DATA NAME="tau_y_correction_scale" UNITS="dimensionless" TYPE="real">
!  A scale multiplying the flux correction for tau_y.
!  Default tau_y_correction_scale=0.0.
!  </DATA> 
!  <DATA NAME="do_bitwise_exact_sum" TYPE="logical">
!  Set true to do bitwise exact global sum. When it is false, the global
!  sum will be non-bitwise_exact, but will significantly increase efficiency.
!  The default value is do_bitwise_exact_sum=.true. in order to ensure answers
!  do not change when alter processors.  But if wish to enhance the efficiency
!  of coupled ocean-ice models that use one of the global normalization options
!  zero_net_salt_restore        =.true.
!  zero_net_salt_correction     =.true.
!  zero_net_water_restore       =.true.
!  zero_net_water_correction    =.true.
!  zero_net_water_coupler       =.true.
!  zero_net_water_couple_restore=.true.
!  then one may wish to consider setting do_bitwise_exact_sum=.false.
!  </DATA>
!
!  <DATA NAME="constant_hlf" TYPE="logical">
!  Treat latent heat of fusion as a constant. Otherwise, use the TEOS-10
!  approach in which hlf is function of surface salinity.
!  Note, TEOS-10 approach is only valid using Absolute Salinity and 
!  conservative temperature as the prognostic fields.  
!  Default constant_hlf = .true., which is the case for pre-TEOS-10 methods. 
!  </DATA>
!  <DATA NAME="constant_hlv" TYPE="logical">
!  Treat latent heat of vaporization as a constant. Otherwise, use the TEOS-10
!  approach in which hlf is function of surface salinity.
!  Note, TEOS-10 approach is only valid using Absolute Salinity and 
!  conservative temperature as the prognostic fields.  
!  Default constant_hlv = .true., which is the case for pre-TEOS-10 methods. 
!  </DATA>
!
!  <DATA NAME="read_stokes_drift" TYPE="logical">
!  This option is to be used when coupling to a surface wave model such as 
!  Wavewatch III that provides both the Stokes drift (m/s) velocity at the 
!  ocean surface, and a decay scale for projecting the Stokes
!  drift into the interior.  Default read_stokes_drift = .false.
!  </DATA>
!
!  <DATA NAME="do_langmuir" TYPE="logical">
!  TThis option exists in the event that boundary forcing from a future wave model is
!  provided. Not to be confused with wave mixing parameterised by supplying 10m
!  winds.
!  Default do_langmuir = .false.
!  </DATA>
!
!  <DATA NAME="do_ustar_correction" TYPE="logical">
!  Compute ustar including the adjusted/correction stress field that has now been
!  included in smf.
!  Note, however, that the FAFMIP stess experiment says we should NOT include the perturbed
!  stress when computing ustar since we do not wish to affect the mixing schemes.
!  Hence, we should set do_ustar_correction=.false. for FAFMIP in which case ustar just has
!  contributions from the unperturbed stress.
!  Default do_ustar_correction = .true. as this reproduces earlier
!  behavior.
!  </DATA>
!
!  <DATA NAME="do_frazil_redist" TYPE="logical">
!  In FAFMIP heat experiments we should be using the heat due to frazil
!  formation from the redistributed tracer. Previous code unconditionally
!  used the standard tracer. We allow the  user to override the recommended treatment to
!  recover old results by setting  do_frazil_redist=.false.. If there is no
!  frazil_redist_tracer this flag has no effect and the usual treatment of frazil proceeds.
!  This flag has no effect on the ACCESS treatment of frazil which ALWAYS uses the redistibuted
!  heat  version if it is available.
!  Note that the current approach (June 2019) the frazil heats are note quite
!  the same for the temperature and redistributed heat tracers. See the references above 
!  for details.
!  Default do_frazil_redist = .true. 
!  </DATA>
!
!</NAMELIST>
!

#include <fms_platform.h>

use constants_mod,            only: epsln, hlv, hlf, kelvin, pi
use diag_manager_mod,         only: register_diag_field, register_static_field, send_data
use fms_mod,                  only: open_namelist_file, check_nml_error, file_exist
use fms_mod,                  only: close_file, read_data, write_version_number
use fms_io_mod,               only: register_restart_field, save_restart, restore_state, restart_file_type
use mpp_domains_mod,          only: mpp_update_domains, BGRID_NE, mpp_define_domains, mpp_get_compute_domain
use mpp_domains_mod,          only: mpp_global_sum, BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
use mpp_domains_mod,          only: mpp_define_io_domain
use mpp_mod,                  only: input_nml_file, mpp_error, FATAL, stdout, stdlog
use time_interp_external_mod, only: time_interp_external, init_external_field
use time_manager_mod,         only: time_type, increment_time, get_time

use ocean_domains_mod,        only: get_local_indices
use ocean_parameters_mod,     only: missing_value, onehalf 
use ocean_parameters_mod,     only: grav, rho_cp, cp_ocean, cp_liquid_runoff, cp_solid_runoff, rho0, rho0r 
use ocean_parameters_mod,     only: CONSERVATIVE_TEMP, POTENTIAL_TEMP, PREFORMED_SALT, PRACTICAL_SALT
use ocean_parameters_mod,     only: MOM_BGRID, MOM_CGRID 
use ocean_riverspread_mod,    only: spread_river_horz
use ocean_tempsalt_mod,       only: pottemp_from_contemp
use ocean_tpm_mod,            only: ocean_tpm_sum_sfc, ocean_tpm_avg_sfc, ocean_tpm_sbc
use ocean_tpm_mod,            only: ocean_tpm_zero_sfc, ocean_tpm_sfc_end
use ocean_types_mod,          only: ocean_grid_type, ocean_domain_type, ocean_public_type
use ocean_types_mod,          only: ocean_time_type, ocean_thickness_type 
use ocean_types_mod,          only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,          only: ocean_external_mode_type, ocean_velocity_type 
use ocean_types_mod,          only: ice_ocean_boundary_type, ocean_density_type
use ocean_types_mod,          only: ocean_public_type
use ocean_workspace_mod,      only: wrk1_2d, wrk2_2d, wrk3_2d, wrk1
use ocean_util_mod,           only: diagnose_2d, diagnose_2d_u, diagnose_3d_u, diagnose_sum
use ocean_tracer_util_mod,    only: diagnose_3d_rho

#if defined(CSIRO_BGC)
use csiro_bgc_mod,            only: csiro_bgc_virtual_fluxes, do_csiro_bgc,ind_no3,ind_phy
#endif
implicit none

private

! for Bgrid or Cgrid 
integer :: horz_grid 
real, allocatable, dimension(:,:) :: wind_mask 

! for restoring input field
integer, allocatable, dimension(:) :: id_restore
integer :: id_eta_restore=-1 

! for flux corrections input fields  
integer :: id_tau_x_correction =-1
integer :: id_tau_y_correction =-1
integer, allocatable, dimension(:) :: id_correction

! for data override files associated with water contributing to latent heating
integer :: id_calving_override = -1
integer :: id_fprec_override   = -1
integer :: id_evap_override    = -1


! for non-constant latent heats 
real, allocatable, dimension(:,:) :: latent_heat_vapor
real, allocatable, dimension(:,:) :: latent_heat_fusion

integer :: index_temp         =-1
integer :: index_salt         =-1
integer :: index_diag_temp    =-1
integer :: index_frazil       =-1
integer :: prog_temp_variable =-1
integer :: prog_salt_variable =-1

! FAFMIP heat tracers 
integer :: index_added_heat  = -1
integer :: index_redist_heat = -1
integer :: index_frazil_redist =-1

integer :: memuse
integer :: num_prog_tracers
integer :: num_diag_tracers
integer :: global_sum_flag


! ids for diagnostic manager 
logical :: used
integer, allocatable, dimension(:) :: id_stf_coupler
integer, allocatable, dimension(:) :: id_stf_restore
integer, allocatable, dimension(:) :: id_stf_correct
integer, allocatable, dimension(:) :: id_stf_total
integer, allocatable, dimension(:) :: id_stf_runoff
integer, allocatable, dimension(:) :: id_stf_calving
integer, allocatable, dimension(:) :: id_stf_pme
integer, allocatable, dimension(:) :: id_stf_pme_on_nrho
integer, allocatable, dimension(:) :: id_stf_prec
integer, allocatable, dimension(:) :: id_stf_evap
integer, allocatable, dimension(:) :: id_trunoff
integer, allocatable, dimension(:) :: id_tcalving 
integer, allocatable, dimension(:) :: id_triver 

integer, allocatable, dimension(:) :: id_total_ocean_stf_coupler
integer, allocatable, dimension(:) :: id_total_ocean_stf_runoff
integer, allocatable, dimension(:) :: id_total_ocean_stf_calving 
integer, allocatable, dimension(:) :: id_total_ocean_stf_pme
integer, allocatable, dimension(:) :: id_total_ocean_stf_prec
integer, allocatable, dimension(:) :: id_total_ocean_stf_evap
integer, allocatable, dimension(:) :: id_total_ocean_stf_restore
integer, allocatable, dimension(:) :: id_total_ocean_stf_correct
integer, allocatable, dimension(:) :: id_total_ocean_stf_sum

integer :: id_tau_x_flux_correction=-1
integer :: id_tau_y_flux_correction=-1
integer :: id_tau_x_net=-1
integer :: id_tau_y_net=-1

integer :: id_latent_heat_vapor =-1
integer :: id_latent_heat_fusion=-1

integer :: id_ustokes      =-1
integer :: id_vstokes      =-1
integer :: id_stokes_depth =-1

integer :: id_ustoke          =-1
integer :: id_vstoke          =-1
integer :: id_wavlen          =-1

integer :: id_net_sfc_heating       =-1
integer :: id_total_net_sfc_heating =-1

integer :: id_net_sfc_workq       =-1
integer :: id_net_sfc_workemp     =-1

integer :: id_salt_flux_ice      =-1
integer :: id_total_salt_flux_ice=-1
integer :: id_temp_runoff_eff    =-1
integer :: id_temp_calving_eff   =-1

integer :: id_tau_x          =-1
integer :: id_tau_y          =-1
integer :: id_tau_curl       =-1
integer :: id_ekman_we       =-1
integer :: id_ekman_heat     =-1
integer :: id_swflx          =-1
integer :: id_swflx_vis      =-1

integer :: id_lw_heat            =-1
integer :: id_sens_heat          =-1
integer :: id_fprec_melt_heat    =-1
integer :: id_calving_melt_heat  =-1
integer :: id_evap_heat          =-1

integer :: id_fprec          =-1
integer :: id_lprec          =-1
integer :: id_river          =-1
integer :: id_alphasfc       =-1
integer :: id_betasfc        =-1
integer :: id_alphasfc2      =-1
integer :: id_betasfc2       =-1
integer :: id_calving        =-1
integer :: id_ideal_calving  =-1
integer :: id_runoff         =-1
integer :: id_ideal_runoff   =-1
integer :: id_melt           =-1
integer :: id_evap           =-1
integer :: id_pme_sbc        =-1
integer :: id_pme_river      =-1
integer :: id_pme_restore    =-1
integer :: id_pme_eta_restore=-1
integer :: id_pme_correct    =-1
integer :: id_pme_net        =-1
integer :: id_ice_mask       =-1
integer :: id_open_ocean_mask=-1
integer :: id_restore_mask   =-1
#if defined(ACCESS_CM) || defined(ACCESS_OM)
integer :: id_wfimelt        =-1
integer :: id_wfiform        =-1
integer :: id_aice               =-1
integer :: id_wnd                =-1
integer :: id_licefw        =-1
integer :: id_liceht        =-1
integer :: id_mh_flux            =-1
integer :: id_atm_co2            =-1
#endif
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
integer :: id_iof_nit               =-1
integer :: id_iof_alg               =-1
#endif

! ids for sea level forcing fields 
real    :: cellarea_r
logical :: diagnose_sea_level_forcing =.false.  ! internally set 
integer :: id_eta_tend_sw           = -1
integer :: id_eta_tend_lw           = -1
integer :: id_eta_tend_sens         = -1
integer :: id_eta_tend_evap_heat    = -1
integer :: id_eta_tend_fprec_melt   = -1
integer :: id_eta_tend_iceberg_melt = -1
integer :: id_eta_tend_heat_coupler = -1
integer :: id_eta_tend_heat_restore = -1

integer :: id_eta_tend_salt_coupler = -1
integer :: id_eta_tend_salt_restore = -1

integer :: id_eta_tend_evap          = -1
integer :: id_eta_tend_lprec         = -1
integer :: id_eta_tend_fprec         = -1
integer :: id_eta_tend_runoff        = -1
integer :: id_eta_tend_iceberg       = -1
integer :: id_eta_tend_water_coupler = -1
integer :: id_eta_tend_water_restore = -1


integer :: id_eta_tend_sw_glob           = -1
integer :: id_eta_tend_lw_glob           = -1
integer :: id_eta_tend_sens_glob         = -1
integer :: id_eta_tend_evap_heat_glob    = -1
integer :: id_eta_tend_fprec_melt_glob   = -1
integer :: id_eta_tend_iceberg_melt_glob = -1
integer :: id_eta_tend_heat_coupler_glob = -1
integer :: id_eta_tend_heat_restore_glob = -1

integer :: id_eta_tend_salt_coupler_glob = -1
integer :: id_eta_tend_salt_restore_glob = -1

integer :: id_eta_tend_evap_glob          = -1
integer :: id_eta_tend_lprec_glob         = -1
integer :: id_eta_tend_fprec_glob         = -1
integer :: id_eta_tend_runoff_glob        = -1
integer :: id_eta_tend_iceberg_glob       = -1
integer :: id_eta_tend_water_coupler_glob = -1
integer :: id_eta_tend_water_restore_glob = -1


! ids for scalar fields 

integer :: id_total_ocean_swflx             =-1
integer :: id_total_ocean_swflx_vis         =-1
integer :: id_total_ocean_evap_heat         =-1
integer :: id_total_ocean_lw_heat           =-1
integer :: id_total_ocean_sens_heat         =-1
integer :: id_total_ocean_river_heat        =-1
integer :: id_total_ocean_pme_heat          =-1
integer :: id_total_ocean_fprec_melt_heat   =-1
integer :: id_total_ocean_calving_melt_heat =-1

integer :: id_total_ocean_river      =-1
integer :: id_total_ocean_evap       =-1
integer :: id_total_ocean_melt       =-1
integer :: id_total_ocean_pme_sbc    =-1
integer :: id_total_ocean_pme_restore=-1
integer :: id_total_ocean_pme_correct=-1
integer :: id_total_ocean_pme_net    =-1
integer :: id_total_ocean_pme_river  =-1

integer :: id_total_ocean_fprec   =-1
integer :: id_total_ocean_lprec   =-1
integer :: id_total_ocean_calving =-1
integer :: id_total_ocean_runoff  =-1
#if defined(ACCESS_CM) || defined(ACCESS_OM)
integer :: id_total_ocean_wfimelt =-1
integer :: id_total_ocean_wfiform =-1
integer :: id_total_ocean_licefw  =-1
integer :: id_total_ocean_liceht  =-1
integer :: id_total_ocean_mh_flux =-1
#endif


! ids for rebinning mass fluxes to neutral density classes 
integer  :: id_mass_precip_on_nrho  =-1
integer  :: id_mass_evap_on_nrho    =-1
integer  :: id_mass_river_on_nrho   =-1
integer  :: id_mass_melt_on_nrho    =-1
integer  :: id_mass_pmepr_on_nrho   =-1
integer  :: id_mass_pme_adj_on_nrho =-1

! ids for rebinning temp and salt fluxes to neutral density classes 
integer  :: id_tform_rho_pbl_flux_on_nrho   =-1
integer  :: id_tform_rho_pbl_adjheat_on_nrho=-1
integer  :: id_tform_rho_pbl_adjsalt_on_nrho=-1
integer  :: id_tform_rho_pbl_heat_on_nrho   =-1
integer  :: id_tform_rho_pbl_salt_on_nrho   =-1
integer  :: id_tform_rho_pbl_sw_on_nrho     =-1
integer  :: id_tform_rho_pbl_lw_on_nrho     =-1
integer  :: id_tform_rho_pbl_sens_on_nrho   =-1
integer  :: id_tform_rho_pbl_lat_on_nrho    =-1


#include <ocean_memory.h>

#ifdef  MOM_STATIC_ARRAYS

real, dimension(isd:ied,jsd:jed) :: data
real, dimension(isd:ied,jsd:jed) :: pme_taum1     ! mass flux (kg/(m^2 sec)) of precip-evap from coupler at taum1 time step 
real, dimension(isd:ied,jsd:jed) :: river_taum1   ! mass flux of river water (liquid+solid) from coupler at taum1 time step 
real, dimension(isd:ied,jsd:jed) :: pme_river     ! mass flux of water into ocean from pme+river-melt
real, dimension(isd:ied,jsd:jed) :: restore_mask  ! mask for setting regions that are restored 
real, dimension(isd:ied,jsd:jed) :: runoff        ! mass flux of liquid river runoff 
real, dimension(isd:ied,jsd:jed) :: ideal_runoff  ! mass flux of liquid river runoff obtained from read-in file  
real, dimension(isd:ied,jsd:jed) :: calving       ! mass flux of calving land ice into ocean 
real, dimension(isd:ied,jsd:jed) :: ideal_calving ! mass flux of calving land ice obtained from read-in file 
real, dimension(isd:ied,jsd:jed) :: rhosfc_inv    ! surface ocean specific volume (m^3/kg) #
real, dimension(isd:ied,jsd:jed) :: alphasfc      ! surface thermal expansion coefficient (1/deg C) 
real, dimension(isd:ied,jsd:jed) :: betasfc       ! surface saline contraction coefficient (1/ppt) 
real, dimension(isd:ied,jsd:jed) :: alphasfc2     ! potrho surface thermal expansion coefficient (1/deg C) 
real, dimension(isd:ied,jsd:jed) :: betasfc2      ! potrho surface saline contraction coefficient (1/ppt) 

#else

real, allocatable, dimension(:,:) :: data
real, allocatable, dimension(:,:) :: pme_taum1     ! mass flux (kg/(m^2 sec)) of precip-evap from coupler at taum1 time step 
real, allocatable, dimension(:,:) :: river_taum1   ! mass flux of river water (liquid+solid) from coupler at taum1 time step
real, allocatable, dimension(:,:) :: pme_river     ! mass flux of water into ocean from pme+river-melt
real, allocatable, dimension(:,:) :: restore_mask  ! mask for setting regions that are restored 
real, allocatable, dimension(:,:) :: runoff        ! mass flux of liquid river runoff  
real, allocatable, dimension(:,:) :: ideal_runoff  ! mass flux of liquid river runoff obtained from read-in file   
real, allocatable, dimension(:,:) :: calving       ! mass flux of calving land ice into ocean 
real, allocatable, dimension(:,:) :: ideal_calving ! mass flux of calving land ice obtained from read-in file
real, allocatable, dimension(:,:) :: rhosfc_inv    ! surface ocean specific volume (m^3/kg) 
real, allocatable, dimension(:,:) :: alphasfc      ! surface thermal expansion coefficient (1/deg C) 
real, allocatable, dimension(:,:) :: betasfc       ! surface saline contraction coefficient (1/ppt) 
real, allocatable, dimension(:,:) :: alphasfc2     ! potrho surface thermal expansion coefficient (1/deg C) 
real, allocatable, dimension(:,:) :: betasfc2      ! potrho surface saline contraction coefficient (1/ppt) 

#endif
#if defined(ACCESS_CM) || defined(ACCESS_OM)
real, allocatable, dimension(:,:,:) :: sslope
real, allocatable, dimension(:,:) :: aice, iof_nit, iof_alg
#endif
#if defined(ACCESS_CM)
real, allocatable, dimension(:,:) :: co2flux
real, allocatable, dimension(:,:) :: ocn_co2
real, allocatable, dimension(:,:) :: atm_co2
#endif


! ice-ocean-boundary fields are allocated using absolute
! indices (regardless of whether ocean allocations are static)
integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd
integer :: i_shift, j_shift                      ! shift isc_bnd to isc and jsc_bnd to jsc

real :: grav_rho0_r
real :: cp_liquid_runoff_r
real :: cp_solid_runoff_r
real :: cp_ocean_r
real :: ice_salt_concentration_r
real :: dtime
real :: twopi 
real :: days_in_year_r

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()
type(restart_file_type), save    :: Sbc_restart
type(restart_file_type), save    :: Sfc_restart

public :: ocean_sbc_init
public :: sum_ocean_sfc
public :: avg_ocean_sfc
public :: zero_ocean_sfc
public :: flux_adjust
public :: get_ocean_sbc
public :: initialize_ocean_sfc
public :: ocean_sfc_end
public :: ocean_sfc_restart

private :: ocean_sbc_diag_init
private :: ocean_sbc_diag
private :: compute_latent_heat_vapor
private :: compute_latent_heat_fusion 


character(len=128) :: version=&
     '$Id: ocean_sbc.F90,v 20.0 2013/12/14 00:10:59 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

logical :: module_is_initialized          =.false.
logical :: use_waterflux                  =.true.
logical :: waterflux_tavg                 =.false.
logical :: use_waterflux_override_calving =.false.
logical :: use_waterflux_override_fprec   =.false.
logical :: use_waterflux_override_evap    =.false.
logical :: rotate_winds                   =.false.
logical :: taux_sinx                      =.false.
logical :: tauy_siny                      =.false.
logical :: runoffspread                   =.false.
logical :: calvingspread                  =.false.
logical :: salt_restore_under_ice         =.true.
logical :: salt_restore_as_salt_flux      =.true.
logical :: zero_net_salt_restore          =.false.
logical :: zero_net_salt_correction       =.false.
logical :: zero_net_water_restore         =.false.
logical :: zero_net_water_correction      =.false.
logical :: zero_net_water_coupler         =.false.
logical :: zero_net_water_couple_restore  =.false.
logical :: zero_net_pme_eta_restore       =.false.
logical :: debug_water_fluxes             =.false.
logical :: zero_water_fluxes              =.false. 
logical :: zero_pme_fluxes                =.false. 
logical :: zero_calving_fluxes            =.false. 
logical :: zero_runoff_fluxes             =.false. 
logical :: zero_river_fluxes              =.false. 
logical :: convert_river_to_pme           =.false.
logical :: zero_heat_fluxes               =.false. 
logical :: zero_surface_stress            =.false.
logical :: read_restore_mask              =.false. 
logical :: restore_mask_gfdl              =.false.
logical :: land_model_heat_fluxes         =.false. 
logical :: do_flux_correction             =.false.
logical :: sbc_heat_fluxes_const          =.false.
logical :: sbc_heat_fluxes_const_seasonal =.false.
logical :: use_constant_sss_for_restore   =.false.
logical :: use_constant_sst_for_restore   =.false.
logical :: use_ideal_runoff               =.false.
logical :: use_ideal_calving              =.false.
logical :: read_stokes_drift              =.false.
logical :: do_langmuir                    =.false.
logical :: do_ustar_correction            =.true.  ! In FAFMIP stress make this falsel
logical :: do_frazil_redist               =.true.  ! In FAFMIP heat make this false to recover old (not recommended) behaviour.

real    :: constant_sss_for_restore       = 35.0
real    :: constant_sst_for_restore       = 12.0
real    :: sbc_heat_fluxes_const_value    = 0.0    ! W/m2
real    :: ice_salt_concentration         = 0.005  ! kg/kg
real    :: ocean_ice_salt_limit           = 0.0    ! kg/kg
real    :: runoff_salinity                = 0.0    ! psu
real    :: runoff_temp_min                = 0.0    ! degC
real    :: temp_restore_tscale            = -30.
real    :: salt_restore_tscale            = -30.
real    :: eta_restore_tscale             = -30.
real    :: max_ice_thickness              = 5.0 
real    :: salinity_ref                   = 35.0
real    :: max_delta_salinity_restore     = -0.5
real    :: salinity_restore_limit_lower   = 0.0
real    :: salinity_restore_limit_upper   = 100.0
real    :: temp_damp_factor                    ! kg/(m^2*sec)
real    :: salt_damp_factor                    ! kg/(m^2*sec)
real    :: eta_damp_factor                     ! 1/sec
real    :: temp_correction_scale      = 0.0
real    :: salt_correction_scale      = 0.0
real    :: tau_x_correction_scale     = 0.0
real    :: tau_y_correction_scale     = 0.0

logical :: constant_hlf               = .true.
logical :: constant_hlv               = .true.
logical :: avg_sfc_velocity           = .true.
logical :: avg_sfc_temp_salt_eta      = .true.
logical :: use_full_patm_for_sea_level= .false. 
logical :: do_bitwise_exact_sum       = .true.

integer :: id_restore_mask_ofam = -1
logical :: restore_mask_ofam = .false.
logical :: river_temp_ofam = .false.

namelist /ocean_sbc_nml/ temp_restore_tscale, salt_restore_tscale, salt_restore_under_ice, salt_restore_as_salt_flux,        &
         eta_restore_tscale, zero_net_pme_eta_restore,                                                                       & 
         rotate_winds, taux_sinx, tauy_siny, use_waterflux, waterflux_tavg, max_ice_thickness, runoffspread, calvingspread,  &
         use_waterflux_override_calving, use_waterflux_override_evap, use_waterflux_override_fprec,                          &
         salinity_ref, zero_net_salt_restore, zero_net_water_restore, zero_net_water_coupler, zero_net_water_couple_restore, &
         zero_net_salt_correction, zero_net_water_correction,                                                                &
         debug_water_fluxes, zero_water_fluxes, zero_calving_fluxes, zero_pme_fluxes, zero_runoff_fluxes, zero_river_fluxes, &
         convert_river_to_pme, zero_heat_fluxes, zero_surface_stress, avg_sfc_velocity, avg_sfc_temp_salt_eta,               &
         ice_salt_concentration, ocean_ice_salt_limit, runoff_salinity, runoff_temp_min, read_restore_mask, restore_mask_gfdl,&
         land_model_heat_fluxes, use_full_patm_for_sea_level, max_delta_salinity_restore, do_flux_correction,                &
         salinity_restore_limit_lower, salinity_restore_limit_upper,                                                          &
         temp_correction_scale, salt_correction_scale, tau_x_correction_scale, tau_y_correction_scale, do_bitwise_exact_sum, &
         sbc_heat_fluxes_const, sbc_heat_fluxes_const_value, sbc_heat_fluxes_const_seasonal,                                 &
         use_constant_sss_for_restore, constant_sss_for_restore, use_constant_sst_for_restore, constant_sst_for_restore,     &
         use_ideal_calving, use_ideal_runoff, constant_hlf, constant_hlv, read_stokes_drift, do_langmuir,                    &
         do_ustar_correction, do_frazil_redist

namelist /ocean_sbc_ofam_nml/ restore_mask_ofam, river_temp_ofam

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_sbc_init">
!
! <DESCRIPTION>
! Initialize the ocean sbc module. 
! </DESCRIPTION>
!
subroutine ocean_sbc_init(Grid, Domain, Time, T_prog, T_diag, &
                          Ocean_sfc, Dens, time_tendency, dtime_t, hor_grid)

  type(ocean_grid_type),          intent(in),    target :: Grid
  type(ocean_domain_type),        intent(in),    target :: Domain
  type(ocean_time_type),          intent(in)            :: Time
  type(ocean_prog_tracer_type),   intent(inout), target :: T_prog(:)
  type(ocean_diag_tracer_type),   intent(inout), target :: T_diag(:)
  type(ocean_public_type),        intent(inout)         :: Ocean_sfc
  type(ocean_density_type),       intent(in)            :: Dens
  character(len=32),              intent(in)            :: time_tendency 
  real,                           intent(in)            :: dtime_t
  integer,                        intent(in)            :: hor_grid

  integer            :: ioun, ierr, io_status
  integer            :: i,j,n
  integer            :: taup1, id_field
  real               :: secday
  character(len=128) :: name, filename

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then
    call mpp_error(FATAL, '==>Error ocean_sbc_init: module has been initialized')
  endif 
  module_is_initialized = .TRUE.

  dtime     = dtime_t 
  horz_grid = hor_grid 

  call write_version_number(version, tagname)

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_sbc_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_sbc_nml')
  ierr = check_nml_error(io_status, 'ocean_sbc_nml')
  read (input_nml_file, nml=ocean_sbc_ofam_nml, iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_sbc_ofam_nml')
#else
  ioun = open_namelist_file()
  read(ioun, ocean_sbc_nml, iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_sbc_nml')
  rewind(ioun)
  read(ioun, ocean_sbc_ofam_nml, iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_sbc_ofam_nml')
  call close_file (ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_sbc_nml)  
  write (stdlogunit, ocean_sbc_nml)
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_sbc_ofam_nml)
  write (stdlogunit, ocean_sbc_ofam_nml)

  if(do_bitwise_exact_sum) then
     global_sum_flag = BITWISE_EXACT_SUM
  else
     global_sum_flag = NON_BITWISE_EXACT_SUM
  endif

#ifndef MOM_STATIC_ARRAYS  
  call get_local_indices(Domain,isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
  allocate(data(isd:ied,jsd:jed))
  allocate(pme_taum1(isd:ied,jsd:jed))
  allocate(river_taum1(isd:ied,jsd:jed))
  allocate(pme_river(isd:ied,jsd:jed))
  allocate(restore_mask(isd:ied,jsd:jed))
  allocate(runoff(isd:ied,jsd:jed))
  allocate(ideal_runoff(isd:ied,jsd:jed))
  allocate(calving(isd:ied,jsd:jed))
  allocate(ideal_calving(isd:ied,jsd:jed))
  allocate(rhosfc_inv(isd:ied,jsd:jed))
  allocate(alphasfc(isd:ied,jsd:jed))
  allocate(betasfc(isd:ied,jsd:jed))
  allocate(alphasfc2(isd:ied,jsd:jed))
  allocate(betasfc2(isd:ied,jsd:jed))
#endif
  data              = 0.0
  pme_taum1         = 0.0
  river_taum1       = 0.0
  pme_river         = 0.0
  restore_mask(:,:) = Grid%tmask(:,:,1)
  runoff            = 0.0
  calving           = 0.0
  ideal_runoff      = 0.0
  ideal_calving     = 0.0
  rhosfc_inv        = 0.0
  alphasfc          = 0.0
  betasfc           = 0.0
  alphasfc2         = 0.0
  betasfc2          = 0.0

  Dom => Domain
  Grd => Grid

  cellarea_r         = 1.0/(epsln + Grd%tcellsurf)
  grav_rho0_r        = rho0r/grav
  secday             = 1.0/(60.0*1440.0)
  taup1              = Time%taup1
  cp_liquid_runoff_r = 1.0/cp_liquid_runoff
  cp_solid_runoff_r  = 1.0/cp_solid_runoff
  cp_ocean_r         = 1.0/cp_ocean 
  days_in_year_r     = 1.0/365.25
  twopi              = 2.0*pi 

  allocate(wind_mask(isd:ied,jsd:jed))
  if(horz_grid == MOM_BGRID) then 
       wind_mask(:,:) = Grd%umask(:,:,1)
  else 
       wind_mask(:,:) = Grd%tmask(:,:,1)
  endif 

  if(ice_salt_concentration > 0.0) then
     ice_salt_concentration_r = 1.0/(epsln+ice_salt_concentration) 
  else
     ice_salt_concentration_r = 0.0
  endif 

  call mpp_define_domains((/1,Grd%ni,1,Grd%nj/),Dom%layout,Ocean_sfc%Domain,maskmap=Dom%maskmap, name='sbc', &
                          x_cyclic_offset = Domain%x_cyclic_offset, y_cyclic_offset = Domain%y_cyclic_offset )  
  call mpp_define_io_domain(Ocean_sfc%Domain, Dom%io_layout)
  call mpp_get_compute_domain(Ocean_sfc%Domain, isc_bnd, iec_bnd, &
                              jsc_bnd, jec_bnd)

  i_shift = isc - isc_bnd
  j_shift = jsc - jsc_bnd

  allocate ( Ocean_sfc%t_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%s_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%u_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%v_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%sea_lev(isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%area   (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%frazil (isc_bnd:iec_bnd,jsc_bnd:jec_bnd))
#if defined(ACCESS_CM) || defined(ACCESS_OM)
  allocate ( Ocean_sfc%gradient (isc_bnd:iec_bnd,jsc_bnd:jec_bnd,2))
  allocate ( sslope(isc:iec, jsc:jec, 2) )
  allocate ( aice(isd:ied, jsd:jed) )
#endif
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
  allocate ( iof_nit(isd:ied, jsd:jed) )
  allocate ( iof_alg(isd:ied, jsd:jed) )
  allocate ( Ocean_sfc%n_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd))
  allocate ( Ocean_sfc%alg_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd))
#endif
#if defined(ACCESS_CM)
  allocate ( Ocean_sfc%co2    (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%co2flux (isc_bnd:iec_bnd,jsc_bnd:jec_bnd)) 
  allocate ( co2flux(isd:ied,jsd:jed),ocn_co2(isd:ied,jsd:jed))
  allocate ( atm_co2(isd:ied,jsd:jed))
#endif

  Ocean_sfc%t_surf  = 0.0  ! time averaged sst (Kelvin) passed to atmosphere/ice model
  Ocean_sfc%s_surf  = 0.0  ! time averaged sss (psu) passed to atmosphere/ice models
  Ocean_sfc%u_surf  = 0.0  ! time averaged u-current (m/sec) passed to atmosphere/ice models
  Ocean_sfc%v_surf  = 0.0  ! time averaged v-current (m/sec)  passed to atmosphere/ice models 
  Ocean_sfc%sea_lev = 0.0  ! time averaged thickness of top model grid cell (m) plus patm/(grav*rho0)
                           ! minus h_geoid - h_tide 
  Ocean_sfc%frazil  = 0.0  ! time accumulated frazil (J/m^2) passed to ice model
#if defined(ACCESS_CM) || defined(ACCESS_OM)
  Ocean_sfc%gradient  = 0.0  ! gradint of ssl passed to Ice model
  sslope = 0.0
  aice = 0.0
#endif
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
  iof_nit = 0.0
  iof_alg = 0.0
  Ocean_sfc%n_surf  = 0.0 
  Ocean_sfc%alg_surf  = 0.0 
#endif
#if defined(ACCESS_CM)
  Ocean_sfc%co2       = 0.0 
  Ocean_sfc%co2flux   = 0.0 
  co2flux             = 0.0
  ocn_co2             = 0.0
  atm_co2             = 0.0
#endif

  Ocean_sfc%area    = Grid%dat(isc:iec, jsc:jec) * Grid%tmask(isc:iec, jsc:jec, 1) !grid cell area

  ! set restore_mask=0.0 in those regions where restoring is NOT applied
  ! and restore_mask=1.0 in regions where restoring is applied.  Default is 
  ! to restore everywhere (restore_mask(:,:) = Grd%tmask(:,:,1)).  
  if(read_restore_mask) then 
      if (restore_mask_ofam) then
        id_restore_mask_ofam = init_external_field('INPUT/restore_mask.nc', &
                                                   'restore_mask', &
                                                   domain=Domain%domain2d)
      else
        call read_data('INPUT/restore_mask', 'restore_mask', data, &
                       Domain%domain2d)
        do j = jsc, jec
            do i = isc, iec
                restore_mask(i,j) = data(i,j)
            end do
        end do
      end if

      ! the following is specific to the mask used at GFDL
      if(restore_mask_gfdl) then 
          do j=jsc,jec
             do i=isc,iec
                if(restore_mask(i,j)==0.0 .or. restore_mask(i,j)>=6.0) then 
                    restore_mask(i,j)=0.0
                else 
                    restore_mask(i,j)=1.0
                endif
             enddo
          enddo
      endif

      call mpp_update_domains(restore_mask(:,:), Dom%domain2d)
  endif

  ! calving > 0 means mass enters ocean.
  ! units kg/(m^2 sec)
  if(use_waterflux_override_calving) then
      filename = 'INPUT/calving_override.nc'
      if (file_exist(trim(filename))) then
          id_calving_override = init_external_field(filename, "calving", domain=Dom%domain2d)
      endif
      if (id_calving_override == -1) then 
          call mpp_error(FATAL,'==>Error in ocean_sbc_mod: did not find calving field in INPUT/calving_override.nc') 
      endif     
      write(stdoutunit,*) &
      '==>Note from ocean_sbc_mod: overriding mass contribution from calving, but not its latent heat contribution.'
  endif 

  ! frozen_precip > 0 means mass enters ocean.
  ! units kg/(m^2 sec)
  if(use_waterflux_override_fprec) then
      filename = 'INPUT/frozen_prec_override.nc'
      if (file_exist(trim(filename))) then
          id_fprec_override = init_external_field(filename, "frozen_precip", domain=Dom%domain2d)
      endif
      if (id_fprec_override == -1) then 
          call mpp_error(FATAL,'==>Error in ocean_sbc_mod: did not find fprec field in INPUT/fprec_override.nc') 
      endif     
      write(stdoutunit,*) &
      '==>Note from ocean_sbc_mod: overriding mass contribution from fprec, but not its latent heat contribution.'
  endif 

  ! evaporation > 0 means mass enters ocean.
  ! units kg/(m^2 sec)
  if(use_waterflux_override_evap) then
      filename = 'INPUT/evaporation_override.nc'
      if (file_exist(trim(filename))) then
          id_evap_override = init_external_field(filename, "evaporation", domain=Dom%domain2d)
      endif
      if (id_evap_override == -1) then 
          call mpp_error(FATAL,'==>Error in ocean_sbc_mod: did not find evap field in INPUT/evap_override.nc') 
      endif     
      write(stdoutunit,*) &
      '==>Note from ocean_sbc_mod: overriding mass contribution from evap, but not its latent heat contribution.'
  endif 

  if(use_ideal_runoff) then
      write(stdoutunit,*) &
      '==>Note from ocean_sbc_mod: adding a static prescribed liquid runoff from a file.'
     call read_data('INPUT/ideal_runoff','ideal_runoff',data,Domain%domain2d)
      do j=jsc,jec
         do i=isc,iec
            ideal_runoff(i,j) = data(i,j)*Grd%tmask(i,j,1)
         enddo
      enddo
  endif

  if(use_ideal_calving) then
      write(stdoutunit,*) &
      '==>Note from ocean_sbc_mod: adding a static prescribed calving solid land ice from a file.'
     call read_data('INPUT/ideal_calving','ideal_calving',data,Domain%domain2d)
      do j=jsc,jec
         do i=isc,iec
            ideal_calving(i,j) = data(i,j)*Grd%tmask(i,j,1)
         enddo
      enddo
  endif
  
  num_prog_tracers = size(T_prog)
  num_diag_tracers = size(T_diag)

  ! for file ids 
  allocate( id_restore    (num_prog_tracers) )
  allocate( id_correction (num_prog_tracers) )
  id_restore   (:)  = -1
  id_correction(:)  = -1

  do n=1,num_prog_tracers
     if (T_prog(n)%name == 'temp')        index_temp        = n
     if (T_prog(n)%name == 'salt')        index_salt        = n
     if (T_prog(n)%name == 'added_heat')  index_added_heat  = n
     if (T_prog(n)%name == 'redist_heat') index_redist_heat = n  
  enddo
   
  do n=1,num_diag_tracers
     if (T_diag(n)%name == 'frazil')   index_frazil    = n
     if (T_diag(n)%name == 'frazil_redist') index_frazil_redist = n
     if (T_diag(n)%name == 'con_temp') index_diag_temp = n
     if (T_diag(n)%name == 'pot_temp') index_diag_temp = n
  enddo

  if (index_temp == -1 .or. index_salt == -1) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_sbc_mod (ocean_sbc_init): temp and/or salt not identified')
  endif 
  if (index_diag_temp == -1) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_sbc_mod (ocean_sbc_init): diagnostic temp not identified')
  endif 
  if (index_frazil_redist > 0 .and. .not. do_frazil_redist) then 
    write(stdoutunit,*) &
    '==>Warning: Using standard frazil calculation despite redistributed version being available'
  endif 
  
  if(T_prog(index_temp)%longname=='Conservative temperature') prog_temp_variable = CONSERVATIVE_TEMP
  if(T_prog(index_temp)%longname=='Potential temperature')    prog_temp_variable = POTENTIAL_TEMP
  if(T_prog(index_salt)%longname=='Preformed Salinity')       prog_salt_variable = PREFORMED_SALT
  if(T_prog(index_salt)%longname=='Practical Salinity')       prog_salt_variable = PRACTICAL_SALT


  ! get file indices for wind stress flux corrections (N/m2).
  ! assume corrections are on B-grid velocity point.
  ! if using C-grid, then corrections will be averaged from B to C grid. 
  name = 'INPUT/tau_x_correction.nc'
  if (file_exist(trim(name)) .and. do_flux_correction) then
      id_tau_x_correction = init_external_field(name, "tau_x", domain=Dom%domain2d)
      write(stdoutunit,*) &
      '==>Note from ocean_sbc_mod: adjusting surface i-directed wind stress with externally read tau_x'
      if (id_tau_x_correction == -1) then 
         call mpp_error(FATAL,'==>Error in ocean_sbc_mod: failed to find tau_x field in INPUT/tau_x_correction.nc') 
      endif 
  endif
  name = 'INPUT/tau_y_correction.nc'
  if (file_exist(trim(name)) .and. do_flux_correction ) then
      id_tau_y_correction = init_external_field(name, "tau_y", domain=Dom%domain2d)
      write(stdoutunit,*) &
      '==>Note from ocean_sbc_mod:  adjusting surface j-directed wind stress with externally read tau_y'
      if (id_tau_y_correction == -1) then 
         call mpp_error(FATAL,'==>Error in ocean_sbc_mod: failed to find tau_y field in INPUT/tau_y_correction.nc') 
      endif 
  endif

  ! get file index for eta restoring (metre)
  name = 'INPUT/eta_t_restore.nc'
  if (file_exist(trim(name)) .and. eta_restore_tscale > 0.0) then
      id_eta_restore = init_external_field(name, "eta_t", domain=Dom%domain2d)
      write(stdoutunit,*) &
      '==>Note from ocean_sbc_mod: damping eta_t to dataset via a modification to the surface water fluxes'
      if (id_eta_restore == -1 .and. eta_restore_tscale > 0.0) then 
         call mpp_error(FATAL,'==>Error in ocean_sbc_mod: failed to find eta field in INPUT/eta_restore.nc') 
      endif 
  endif

  ! loop for initializing tracer fields 
  do n=1,num_prog_tracers

     ! init_external_field(file_name,field_name,domain)
#ifndef MOM_STATIC_ARRAYS         
        allocate(T_prog(n)%stf(isd:ied,jsd:jed))
        allocate(T_prog(n)%tpme(isd:ied,jsd:jed))
        allocate(T_prog(n)%triver(isd:ied,jsd:jed))
        allocate(T_prog(n)%trunoff(isd:ied,jsd:jed))
        allocate(T_prog(n)%tcalving(isd:ied,jsd:jed))
        allocate(T_prog(n)%runoff_tracer_flux(isd:ied,jsd:jed))
        allocate(T_prog(n)%calving_tracer_flux(isd:ied,jsd:jed))
        allocate(T_prog(n)%riverdiffuse(isd:ied,jsd:jed))
#endif

        ! get file indices for restoring fields on temp and salinity  
        name = 'INPUT/'//trim(T_prog(n)%name)//'_sfc_restore.nc'
        if (file_exist(trim(name))) then
            id_restore(n) = init_external_field(name, T_prog(n)%name, domain=Dom%domain2d)
            write(stdoutunit,*) &
            '==>Note from ocean_sbc_mod: applying surface restoring to '//trim(T_prog(n)%name)
            if (id_restore(n) == -1) then
               call mpp_error(FATAL,'==>ocean_sbc_mod: failure to find sfc_restore field in INPUT/_sfc_restore.nc') 
            endif
        elseif(trim(T_prog(n)%name)=='temp' .and. temp_restore_tscale > 0.0) then 
            call mpp_error(FATAL, &
            '==>ocean_sbc_mod: temp_restore_tscale > 0.0 but cannot find INPUT/temp_sfc_restore.nc') 
        elseif(trim(T_prog(n)%name)=='salt' .and. salt_restore_tscale > 0.0) then 
            call mpp_error(FATAL, &
            '==>ocean_sbc_mod: salt_restore_tscale > 0.0 but cannot find INPUT/salt_sfc_restore.nc') 
        endif

        ! get file indices for temp, added_heat, and pme flux correction 
        name = 'INPUT/'//trim(T_prog(n)%name)//'_sfc_correction.nc'
        if (file_exist(trim(name)) .and. do_flux_correction) then
            if  (n == index_temp) then
                id_correction(n) = init_external_field(name, "sfc_hflux", domain=Dom%domain2d)
                write(stdoutunit,*) '==>Note from ocean_sbc_mod: applying surface heat flux correction to '//trim(T_prog(n)%name)
                if (id_correction(n) == -1) then 
                  call mpp_error(FATAL,&
                  '==>Error in ocean_sbc_mod: failure to find temp_sfc_correction field in INPUT/temp_sfc_correction.nc') 
                endif 
            endif
            if  (n == index_added_heat) then
                id_correction(n) = init_external_field(name, "sfc_hflux", domain=Dom%domain2d)
                write(stdoutunit,*) '==>Note from ocean_sbc_mod: applying added surface heat flux from FAFMIP to '//trim(T_prog(n)%name)
                if (id_correction(n) == -1) then
                  call mpp_error(FATAL,&
                  '==>Error in ocean_sbc_mod: failure to find temp_sfc_correction field in INPUT/added_heat_sfc_correction.nc') 
                endif
            endif
            if  (n == index_salt) then
                id_correction(n) = init_external_field(name, "pme", domain=Dom%domain2d)
                write(stdoutunit,*) '==>Note from ocean_sbc_mod: applying surface pme flux correction to '//trim(T_prog(n)%name)
                if (id_correction(n) == -1) then
                  call mpp_error(FATAL, &
                   '==>Error in ocean_sbc_mod: failure to find salt_sfc_correction field in INPUT/salt_sfc_correction.nc') 
                endif 
            endif
            
        endif

        T_prog(n)%stf                 = 0.0
        T_prog(n)%tpme                = 0.0
        T_prog(n)%triver              = 0.0
        T_prog(n)%trunoff             = 0.0
        T_prog(n)%tcalving            = 0.0
        T_prog(n)%runoff_tracer_flux  = 0.0
        T_prog(n)%calving_tracer_flux = 0.0
        T_prog(n)%riverdiffuse        = 0.0
        if (n==index_salt) then
           T_prog(n)%trunoff(:,:) = runoff_salinity*Grd%tmask(:,:,1)
        endif 

  enddo ! for do n=1,num_prog_tracers

  write(stdoutunit,'(/a/)') &
  '==>Note from ocean_sbc_mod: if inputting river water, enable rivermix_mod to get river tracers into ocean.'


  !------for temperature restoring---------- 
  ! dimensions (kg/m^3)*(m/s)
  if (temp_restore_tscale > 0.0) then
     temp_damp_factor = rho0*Grd%dzt(1)*secday/temp_restore_tscale 
  else
     temp_damp_factor = -1.0
     write(stdoutunit,*) &
     '==>Note from ocean_sbc_mod: temp_restore_tscale < 0. no surface restoring for temp'
  endif


  !------for salinity or pme restoring---------- 
  ! dimensions (kg/m^3)*(m/s)
  if (salt_restore_tscale > 0.0) then
     salt_damp_factor = rho0*Grd%dzt(1)*secday/salt_restore_tscale  

     if(salt_restore_under_ice) then 
       write(stdoutunit,*) &
       '==>Note from ocean_sbc_mod: salt_restore_under_ice=.true. => sss restore even under ice.'
     else 
       write(stdoutunit,*) &
       '==>Note from ocean_sbc_mod: salt_restore_under_ice=.false. => no sss restore under ice.'
     endif 

     if(debug_water_fluxes) then 
         if(zero_water_fluxes) then 
             write(stdoutunit,*) &
                  '==>Warning: Over-riding Ice_ocean_boundary to zero the moisture fluxes: pme=river=calving=0.0.'
         endif
         if(zero_calving_fluxes) then 
             write(stdoutunit,*) &
                  '==>Warning: Over-riding Ice_ocean_boundary to zero the calving fluxes: calving=0.0.'
         endif
         if(zero_runoff_fluxes) then 
             write(stdoutunit,*) &
                  '==>Warning: Over-riding Ice_ocean_boundary to zero the river runoff fluxes: runoff=0.0.'
         endif
         if(zero_pme_fluxes) then 
             write(stdoutunit,*) &
                  '==>Warning: Over-riding Ice_ocean_boundary to zero the pme fluxes: pme=0.0.'
         endif
     endif

     if(use_waterflux) then 
       if(zero_net_water_restore) then 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_water_restore=.true.=>zero net restoring water put in ocean.'
       else 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_water_restore=.false.=>nonzero net restoring water put in ocean.'
       endif 
       if(zero_net_water_correction) then 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_water_correction=.true.=>zero net correction water put in ocean.'
       else 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_water_correction=.false.=>nonzero net correction water put in ocean.'
       endif 
       if(zero_net_water_coupler) then 
          write(stdoutunit,*) &
           '==>Note from ocean_sbc_mod: zero_net_water_coupler=.true.=>zero water into ocean via coupler (sans the sea ice).'
       else 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_water_coupler=.false.=>nonzero net water into ocean via coupler.'
       endif 
       if(zero_net_water_couple_restore) then 
          write(stdoutunit,*) &
          '==>ocean_sbc_mod: zero_net_water_couple_restore=.true.=>zero water into ocean from restore + coupler (sans the sea ice).'
       else 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_water_couple_restore=.false.'
       endif 

     else 

       if(zero_net_salt_restore) then 
          write(stdoutunit,*) &
         '==>Note from ocean_sbc_mod: zero_net_salt_restore=.true.=>zero net restoring salt put in ocean.'
       else 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_salt_restore=.false.=>nonzero net restoring salt put in ocean.'
       endif 
       if(zero_net_salt_correction) then 
          write(stdoutunit,*) &
         '==>Note from ocean_sbc_mod: zero_net_salt_correction=.true.=>zero net correction of salt put in ocean.'
       else 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_salt_correction=.false.=>nonzero net correction salt put in ocean.'
       endif 

     endif 

  else  ! no salinity or pme restoring 

     salt_damp_factor = -1.0
     write(stdoutunit,*) &
     '==>Note from ocean_sbc_mod: salt_restore_tscale < 0. no surface restoring for salt or pme.'

  endif


  ! dimensions (1/sec)
  if (eta_restore_tscale > 0.0) then
      eta_damp_factor = secday/eta_restore_tscale 
      if(.not. use_waterflux) then 
          call mpp_error(FATAL, &
               '==>Error in ocean_sbc_mod: must have use_waterflux=.true. to run with eta_restore_tscale > 0') 
      endif
      write(stdoutunit,*) &
           '==>Note from ocean_sbc_mod: eta_restore_tscale > 0. will restore eta to prescribed dataset.'
  else
      eta_damp_factor = -1.0
      write(stdoutunit,*) &
           '==>Note from ocean_sbc_mod: eta_restore_tscale < 0. no surface restoring for eta'
  endif


  if(use_full_patm_for_sea_level) then 
      write(stdoutunit,*) &
           '==>NOTE: Allowing for the ice model to feel the fully depressed sea level for purposes of its dynamics.'
  endif

  if(land_model_heat_fluxes) then 
      write(stdoutunit,*) &
           '==>NOTE: Assuming the land model carries heat of liquid runoff and solid calving.'
      write(stdoutunit,*) &
           '   Be sure to make the appropriate changes ALSO in ocean_rivermix.F90 and ocean_vert_kpp.F90'
  endif

  if(zero_heat_fluxes) then 
      write(stdoutunit,*) &
           '==>Warning: Over-riding Ice_ocean_boundary to zero the heat fluxes: stf(temp)=0.0.'
  endif

  if(sbc_heat_fluxes_const) then 
      write(stdoutunit,*) &
           '==>Warning: Over-riding Ice_ocean_boundary to set surface heat fluxes to global constant.'
  endif

  if(zero_surface_stress) then 
      write(stdoutunit,*) &
           '==>Warning: Over-riding Ice_ocean_boundary to zero the surface stress: smf=0.0.'
  endif

  if(avg_sfc_velocity) then 
      write(stdoutunit,*) &
           '==>If coupling, then avg_sfc_velocity=.true. means will pass averaged ocean velocity to ice model.'
  else 
      write(stdoutunit,*) &
           '==>If coupling, then avg_sfc_velocity=.false. means will pass most recent ocean velocity to ice model.'
  endif
  if(avg_sfc_temp_salt_eta) then 
      write(stdoutunit,*) &
           '==>If coupling, then avg_sfc_temp_salt_eta=.true. means will pass averaged sst, sss, eta to ice model.'
  else 
      write(stdoutunit,*) &
           '==>If coupling, then avg_sfc_velocity=.false. means will pass most recent sst, sss, eta to ice model.'
  endif

  if(waterflux_tavg) then 
      write(stdoutunit,*) &
           '==>Note: waterflux_tavg sets pme+river = avg of ice_ocean_boundary values over'
      write(stdoutunit,*) &
           '         tau and taum1 time steps. This may damp splitting between leap-frog modes.'
      write(stdoutunit,*) &
           '         However, it compromises conservation of mass and tracer.'
      write(stdoutunit,*) &
           '         It is NOT recommended when using the standard MOM two time level scheme.'

      if(time_tendency=='twolevel') then 
          write(stdoutunit,'(/a)') &
               '==>Warning in ocean_sbc_mod: waterflux_tavg=.true. unnecessary with time_tendency==twolevel.'
          write(stdoutunit,'(/a)') &
               '   Strongly recommend setting waterflux_tavg=.false. to conserve mass between component models.'
      endif

      filename = 'ocean_waterflux.res.nc'
      id_field = register_restart_field(Sbc_restart, filename, 'pme_taum1',   pme_taum1,  Domain%domain2d)
      id_field = register_restart_field(Sbc_restart, filename, 'river_taum1', river_taum1,Domain%domain2d)
      if (file_exist('INPUT/ocean_waterflux.res.nc')) then
          call restore_state(Sbc_restart)
      endif

  endif

  ! for latent heating  
  allocate(latent_heat_vapor(isd:ied,jsd:jed))
  allocate(latent_heat_fusion(isd:ied,jsd:jed))
  latent_heat_vapor(:,:)  = hlv*Grd%tmask(:,:,1)
  latent_heat_fusion(:,:) = hlf*Grd%tmask(:,:,1)
  if(constant_hlf) then
      write(stdoutunit,*) '==>Note from ocean_sbc_mod: Using constant latent heat of fusion at ocean surface.'
  else
      if (prog_salt_variable == PRACTICAL_SALT) then
          call mpp_error(FATAL, '==>Error in ocean_sbc_mod:' //                                                     &
               ' Must have prognostic salinity as PREFORMED SALINITY for variable latent heat of fusion')
      endif
      write(stdoutunit,*) '==>Note from ocean_sbc_mod: Using TEOS-10 non-constant latent heat of fusion at ocean surface.'
  endif

  if(constant_hlv) then
      write(stdoutunit,*) '==>Note from ocean_sbc_mod: Using constant latent heat of evaporation at ocean surface.'
  else
      if(prog_temp_variable == POTENTIAL_TEMP)  then
          call mpp_error(FATAL,                                                                                     &
               '==>Error in ocean_sbc_mod:' //                                                                          &
               ' Must have prognostic temperature as CONSERVATIVE TEMPERATURE for non-constant latent heat of evaporation.') 
      endif
      if (prog_salt_variable == PRACTICAL_SALT) then
          call mpp_error(FATAL,'==>Error in ocean_sbc_mod:' //                                                        &
               ' Must have prognostic salinity as PREFORMED SALINITY for non-constant latent heat of evaporation.')
      endif
      write(stdoutunit,*) '==>Note: Using TEOS-10 for non-constant latent heat of evaporation at ocean surface.'
  endif

  ! initialize diagnostic manager fields 
  call ocean_sbc_diag_init(Time, Dens, T_prog)


  return

end subroutine ocean_sbc_init
! </SUBROUTINE> NAME="ocean_sbc_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_sbc_diag_init">
!
! <DESCRIPTION>
! Initialize the ocean sbc diagnostics.
! Send some static diagnostics to diagnostic manager. 
! </DESCRIPTION>
!
subroutine ocean_sbc_diag_init(Time, Dens, T_prog)

  type(ocean_time_type),         intent(in)  :: Time
  type(ocean_density_type),      intent(in)  :: Dens
  type(ocean_prog_tracer_type),  intent(in)  :: T_prog(:)

  integer            :: n
  character(len=128) :: name
  integer            :: stdoutunit
  stdoutunit=stdout()
  

  allocate( id_stf_coupler           (num_prog_tracers) )
  allocate( id_stf_restore           (num_prog_tracers) )
  allocate( id_stf_correct           (num_prog_tracers) )
  allocate( id_stf_total             (num_prog_tracers) )
  allocate( id_stf_runoff            (num_prog_tracers) )
  allocate( id_stf_calving           (num_prog_tracers) )
  allocate( id_stf_pme               (num_prog_tracers) )
  allocate( id_stf_pme_on_nrho       (num_prog_tracers) )
  allocate( id_stf_prec              (num_prog_tracers) )
  allocate( id_stf_evap              (num_prog_tracers) )
  allocate( id_trunoff               (num_prog_tracers) )
  allocate( id_tcalving              (num_prog_tracers) )
  allocate( id_triver                (num_prog_tracers) )

  allocate( id_total_ocean_stf_coupler(num_prog_tracers) )
  allocate( id_total_ocean_stf_runoff (num_prog_tracers) )
  allocate( id_total_ocean_stf_calving(num_prog_tracers) )
  allocate( id_total_ocean_stf_pme    (num_prog_tracers) )
  allocate( id_total_ocean_stf_prec   (num_prog_tracers) )
  allocate( id_total_ocean_stf_evap   (num_prog_tracers) )
  allocate( id_total_ocean_stf_restore(num_prog_tracers) )
  allocate( id_total_ocean_stf_correct(num_prog_tracers) )
  allocate( id_total_ocean_stf_sum    (num_prog_tracers) )

  id_stf_coupler           (:)  = -1
  id_stf_restore           (:)  = -1
  id_stf_correct           (:)  = -1
  id_stf_total             (:)  = -1
  id_stf_runoff            (:)  = -1
  id_stf_calving           (:)  = -1
  id_stf_pme               (:)  = -1
  id_stf_pme_on_nrho       (:)  = -1
  id_stf_prec              (:)  = -1
  id_stf_evap              (:)  = -1
  id_trunoff               (:)  = -1
  id_tcalving              (:)  = -1
  id_triver                (:)  = -1

  id_total_ocean_stf_coupler(:) = -1
  id_total_ocean_stf_runoff (:) = -1
  id_total_ocean_stf_calving(:) = -1
  id_total_ocean_stf_pme    (:) = -1
  id_total_ocean_stf_prec   (:) = -1
  id_total_ocean_stf_evap   (:) = -1
  id_total_ocean_stf_restore(:) = -1
  id_total_ocean_stf_correct(:) = -1
  id_total_ocean_stf_sum    (:) = -1


  ! static fields for diagnostic manager 
  id_ideal_runoff = register_static_field('ocean_model','ideal_runoff',     & 
       Grd%tracer_axes(1:2),'prescribed static runoff flux', 'kg/(m^2 sec)',&
       missing_value=missing_value,range=(/-1e20,1e20/))
  call diagnose_2d(Time, Grd, id_ideal_runoff, ideal_runoff(:,:))

  id_ideal_calving = register_static_field('ocean_model','ideal_calving',    &
       Grd%tracer_axes(1:2),'prescribed static calving flux', 'kg/(m^2 sec)',&
       missing_value=missing_value,range=(/-1e20,1e20/))
  call diagnose_2d(Time, Grd, id_ideal_calving, ideal_calving(:,:))

  id_restore_mask = register_static_field('ocean_model','restore_mask', &
       Grd%tracer_axes(1:2),'restoring mask', 'none' ,                  &
       missing_value=missing_value,range=(/-10.,10./))
  call diagnose_2d(Time, Grd, id_restore_mask, restore_mask(:,:))


  ! register dynamic fields 

  id_ustoke = register_diag_field('ocean_model','ww3 ustoke', Grd%vel_axes_uv(1:2),&
       Time%model_time, 'i-directed stokes drift velocity', 'm/s',                   &
       missing_value=missing_value,range=(/-10.,10./),                       &
       standard_name='surface stokes drift x-velocity')

  id_vstoke = register_diag_field('ocean_model','ww3 vstoke', Grd%vel_axes_uv(1:2),&
       Time%model_time, 'j-directed stokes drift velocity', 'm/s',                   &
       missing_value=missing_value,range=(/-10.,10./),                       &
       standard_name='surface stokes drift y-velocity')

  id_wavlen = register_diag_field('ocean_model','ww3 wavlen', Grd%tracer_axes(1:2),  &
       Time%model_time, 'mean wave length', 'm')
  
  id_tau_x = register_diag_field('ocean_model','tau_x', Grd%vel_axes_u(1:2), &
       Time%model_time, 'i-directed wind stress forcing u-velocity', 'N/m^2',&
       missing_value=missing_value,range=(/-10.,10./),                       &
       standard_name='surface_downward_x_stress')

  id_tau_x_flux_correction = register_diag_field('ocean_model','tau_x_flux_correction',      &
       Grd%vel_axes_u(1:2),                                                                  &
       Time%model_time, 'i-directed wind stress flux correction forcing u-velocity', 'N/m^2',&
       missing_value=missing_value,range=(/-10.,10./),                                       &
       standard_name='surface_downward_x_stress_correction')

  id_tau_x_net = register_diag_field('ocean_model','tau_x_net',                                         &
       Grd%vel_axes_u(1:2),                                                                             &
       Time%model_time, 'net i-directed wind stress (including correction) forcing u-velocity', 'N/m^2',&
       missing_value=missing_value,range=(/-10.,10./))

  id_tau_y = register_diag_field('ocean_model','tau_y', Grd%vel_axes_u(1:2), &
       Time%model_time, 'j-directed wind stress forcing v-velocity', 'N/m^2',&
       missing_value=missing_value,range=(/-10.,10./),                       &
       standard_name='surface_downward_y_stress') 

  id_tau_y_flux_correction = register_diag_field('ocean_model','tau_y_flux_correction',     &
       Grd%vel_axes_v(1:2),                                                                 &
       Time%model_time, 'j-directed wind stress flux corretion forcing v-velocity', 'N/m^2',&
       missing_value=missing_value,range=(/-10.,10./),                                      &
       standard_name='surface_downward_y_stress_correction')

  id_tau_y_net = register_diag_field('ocean_model','tau_y_net',                                         &
       Grd%vel_axes_v(1:2),                                                                             &
       Time%model_time, 'net j-directed wind stress (including correction) forcing v-velocity', 'N/m^2',&
       missing_value=missing_value,range=(/-10.,10./))

  id_tau_curl = register_diag_field('ocean_model','tau_curl', Grd%tracer_axes_flux_y(1:2),&
       Time%model_time, 'wind stress curl averaged to U-point', 'N/m^3',                  &
       missing_value=missing_value,range=(/-10.,10./)) 

  id_ekman_we = register_diag_field('ocean_model','ekman_we', Grd%tracer_axes_flux_y(1:2),&
       Time%model_time, 'Ekman vertical velocity averaged to wt-point', 'm/s',            &
       missing_value=missing_value,range=(/-100.,100./))  

  id_ekman_heat = register_diag_field('ocean_model','ekman_heat',                        &
       Grd%tracer_axes_flux_y(1:2), Time%model_time, 'Ekman Component to heat transport',&
       'Watts', missing_value=missing_value,range=(/-1.e4,1.e4/))           

  id_stokes_depth = register_diag_field ('ocean_model', 'stokes_depth', Grd%vel_axes_uv(1:2), Time%model_time, &
     'Decady depth for Stokes drift from surface waves', 'metre', missing_value=missing_value, range=(/-10.0,1e6/))

  id_ustokes = register_diag_field ('ocean_model', 'ustokes', Grd%vel_axes_uv(1:3), Time%model_time, &
     'i-Stokes drift from surface waves', 'm/sec', missing_value=missing_value, range=(/-1e3,1e3/))

  id_vstokes = register_diag_field ('ocean_model', 'vstokes', Grd%vel_axes_uv(1:3), Time%model_time, &
     'j-Stokes drift from surface waves', 'm/sec', missing_value=missing_value, range=(/-1e3,1e3/))

  id_latent_heat_vapor = register_diag_field('ocean_model','latent_heat_vapor',&
       Grd%tracer_axes(1:2),                                                   &
       Time%model_time, 'latent heat of vaporization at ocean surface',        &
       'J/kg' , missing_value=missing_value,range=(/0.0,1.e7/))            

  id_latent_heat_fusion = register_diag_field('ocean_model','latent_heat_fusion',&
       Grd%tracer_axes(1:2),                                                     &
       Time%model_time, 'latent heat of fusion for water at ocean surface',      &
       'J/kg' , missing_value=missing_value,range=(/0.0,1.e7/))            

  id_temp_runoff_eff = register_diag_field('ocean_model','temp_runoff_eff',        &
       Grd%tracer_axes(1:2),                                                       &
       Time%model_time, 'effective temp of liquid river runoff entering the ocean',&
       'degC' , missing_value=missing_value,range=(/-1.e4,1.e4/))            

  id_temp_calving_eff = register_diag_field('ocean_model','temp_calving_eff',&
       Grd%tracer_axes(1:2),                                                 &
       Time%model_time, 'effective temp of land ice calving into ocean',     & 
       'degC' , missing_value=missing_value,range=(/-1.e4,1.e4/))            

  id_ice_mask = register_diag_field('ocean_model','ice_mask', Grd%tracer_axes(1:2),&
       Time%model_time, 'ice mask according to near-frazil condition', 'none' ,    &
       missing_value=missing_value,range=(/-10.,10./))

  id_open_ocean_mask = register_diag_field('ocean_model','open_ocean_mask', Grd%tracer_axes(1:2),&
       Time%model_time, 'open-ocean mask according to near-frazil condition', 'none' ,           &
       missing_value=missing_value,range=(/-10.,10./))

  id_river = register_diag_field('ocean_model','river', Grd%tracer_axes(1:2),   &
       Time%model_time, 'mass flux of river (runoff + calving) entering ocean', &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/))

  id_alphasfc = register_diag_field('ocean_model','alphasfc', Grd%tracer_axes(1:2),&
       Time%model_time, 'Thermal expansion coeff at surface (-1/rho drho/dT) ',    &
       '1/(deg C)', missing_value=missing_value,range=(/-1e6,1e6/))

  id_betasfc = register_diag_field('ocean_model','betasfc', Grd%tracer_axes(1:2),&
       Time%model_time, 'Saline contraction coeff at surface (1/rho drho/dS) ',  &
       '1/psu', missing_value=missing_value,range=(/-1e6,1e6/))

  id_alphasfc2 = register_diag_field('ocean_model','alphasfc2', Grd%tracer_axes(1:2),      &
       Time%model_time, 'Potrho thermal expansion coeff at surface (-1/potrho dpotrho/dT)',&
       '1/(deg C)', missing_value=missing_value,range=(/-1e6,1e6/))

  id_betasfc2 = register_diag_field('ocean_model','betasfc2', Grd%tracer_axes(1:2),        &
       Time%model_time, 'Potrho saline contraction coeff at surface (1/potrho dpotrho/dS)',&
       '1/psu', missing_value=missing_value,range=(/-1e6,1e6/))

  id_runoff = register_diag_field('ocean_model','runoff', Grd%tracer_axes(1:2),&
       Time%model_time, 'mass flux of liquid river runoff entering ocean ',    &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/),     &
       standard_name='water_flux_into_sea_water_from_rivers')

  id_calving = register_diag_field('ocean_model','ice_calving', Grd%tracer_axes(1:2),&
       Time%model_time, 'mass flux of land ice calving into ocean',                 &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/),          &
       standard_name='water_flux_into_sea_water_from_icebergs')

  id_evap = register_diag_field('ocean_model','evap', Grd%tracer_axes(1:2),         &
       Time%model_time, 'mass flux from evaporation/condensation (>0 enters ocean)',&
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/),          &
       standard_name='water_evaporation_flux' )

  id_melt = register_diag_field('ocean_model','melt', Grd%tracer_axes(1:2),               &
       Time%model_time, 'water flux transferred with sea ice form/melt (>0 enters ocean)',&
       '(kg/m^3)*(m/sec)',missing_value=missing_value,range=(/-1e6,1e6/),                 &
       standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics')

  id_pme_river= register_diag_field('ocean_model','pme_river', Grd%tracer_axes(1:2),           &
       Time%model_time, 'mass flux of precip-evap+river via sbc (liquid, frozen, evaporation)',&
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/),                     &
       standard_name='water_flux_into_sea_water')

  id_net_sfc_workemp = register_diag_field('ocean_model','net_sfc_workEmP',                        &
        Grd%tracer_axes(1:2),                                                               &
        Time%model_time, 'pme_river*g*beta2*So/rho0, beta uses pot_rho rather than neut',    &
        'm^2/s^-3' ,                                                                        &
         missing_value=missing_value,range=(/-1.e4,1.e4/))
  id_pme_sbc = register_diag_field('ocean_model','pme_sbc', Grd%tracer_axes(1:2),&
       Time%model_time, 'precip-evap via sbc (liquid, frozen, evaporation)',     &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/))

  id_pme_restore = register_diag_field('ocean_model','pme_restore', Grd%tracer_axes(1:2),&
       Time%model_time, 'precip-evap from restoring of SSS (>0 enters ocean)',           &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/))   

  id_pme_eta_restore = register_diag_field('ocean_model','pme_eta_restore', Grd%tracer_axes(1:2),&
       Time%model_time, 'precip-evap from restoring of eta (>0 enters ocean)',                   &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/))   

  id_pme_correct = register_diag_field('ocean_model','pme_correct', Grd%tracer_axes(1:2),&
       Time%model_time, 'precip-evap from flux correction (>0 enters ocean)',            &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/))   

  id_pme_net = register_diag_field('ocean_model','pme_net', Grd%tracer_axes(1:2),&
       Time%model_time, 'precip-evap into ocean (total w/ restore + normalize)', &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/))   

  id_fprec = register_diag_field('ocean_model','fprec', Grd%tracer_axes(1:2),  &
       Time%model_time, 'snow falling onto ocean (>0 enters ocean)',           &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1.e10,1.e10/), &
       standard_name='snowfall_flux')   

  id_lprec = register_diag_field('ocean_model','lprec', Grd%tracer_axes(1:2),                  &
       Time%model_time, 'liquid precip (including ice melt/form) into ocean (>0 enters ocean)',&
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1.e10,1.e10/),                 &
       standard_name='rainfall_flux')   

#if defined(ACCESS_CM) || defined(ACCESS_OM)
 id_wfimelt = register_diag_field('ocean_model','wfimelt', Grd%tracer_axes(1:2),  &
       Time%model_time, 'water into ocean due to ice melt (>0 enters ocean)',         &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1.e-1,1.e-1/),&
       standard_name='icemelt_flux')

 id_wfiform = register_diag_field('ocean_model','wfiform', Grd%tracer_axes(1:2),  &
       Time%model_time, 'water out of ocean due to ice form (>0 enters ocean)',         &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1.e-1,1.e-1/),&
       standard_name='iceform_flux')

  id_licefw = register_diag_field('ocean_model','licefw', Grd%tracer_axes(1:2), &
        Time%model_time, 'water into ocean due to land ice discharge (>0 enters ocean)', &
        '(kg/m^2/sec)', missing_value=missing_value,range=(/-1.e10,1.e10/),&
        standard_name='licefw_flux')

  id_liceht = register_diag_field('ocean_model','liceht', Grd%tracer_axes(1:2), &
        Time%model_time, 'heat into ocean due to land ice discharge-melt (>0 heats ocean)', &
        '(W/m^2)', missing_value=missing_value,range=(/-1.e10,1.e10/),&
        standard_name='liceht_flux')

  id_mh_flux = register_diag_field('ocean_model','mh_flux', Grd%tracer_axes(1:2), &
        Time%model_time, 'heat into ocean due to melting ice (>0 heats ocean)', &
        '(W/m^2)', missing_value=missing_value,range=(/-1.e10,1.e10/),&
        standard_name='mh_flux')
#endif

  id_swflx = register_diag_field('ocean_model','swflx', Grd%tracer_axes(1:2),  &
       Time%model_time, 'shortwave flux into ocean (>0 heats ocean)', 'W/m^2', &
       missing_value=missing_value,range=(/-1.e10,1.e10/),                     &
       standard_name='surface_net_downward_shortwave_flux')   

  id_swflx_vis = register_diag_field('ocean_model','swflx_vis', Grd%tracer_axes(1:2),&
       Time%model_time, 'visible shortwave into ocean (>0 heats ocean)', 'W/m^2' ,   &
       missing_value=missing_value,range=(/-1.e10,1.e10/))   

  id_evap_heat = register_diag_field('ocean_model','evap_heat', Grd%tracer_axes(1:2),&
       Time%model_time, 'latent heat flux into ocean (<0 cools ocean)', 'W/m^2',     &
       missing_value=missing_value,range=(/-1e10,1e10/),                             &
       standard_name='surface_downward_latent_heat_flux')

  id_lw_heat = register_diag_field('ocean_model','lw_heat', Grd%tracer_axes(1:2),&
       Time%model_time, 'longwave flux into ocean (<0 cools ocean)', 'W/m^2' ,  &
       missing_value=missing_value,range=(/-1.e10,1.e10/),                      &
       standard_name='surface_net_downward_longwave_flux' )   

  id_sens_heat = register_diag_field('ocean_model','sens_heat', Grd%tracer_axes(1:2),&
       Time%model_time, 'sensible heat into ocean (<0 cools ocean)', 'W/m^2' ,      &
       missing_value=missing_value,range=(/-1.e10,1.e10/),                          &
       standard_name='surface_downward_sensible_heat_flux')   

  id_fprec_melt_heat = register_diag_field('ocean_model','fprec_melt_heat', Grd%tracer_axes(1:2),&
       Time%model_time, 'heat flux to melt frozen precip (<0 cools ocean)', 'W/m^2' , &
       missing_value=missing_value,range=(/-1.e10,1.e10/),                            &
       standard_name='heat_flux_into_sea_water_due_to_snow_thermodynamics')   

  id_calving_melt_heat = register_diag_field('ocean_model','calving_melt_heat', Grd%tracer_axes(1:2),&
       Time%model_time, 'heat flux needed to melt calving ice (<0 cools ocean)', 'W/m^2',           &
       missing_value=missing_value,range=(/-1.e10,1.e10/),                                          &
       standard_name='heat_flux_into_sea_water_due_to_iceberg_thermodynamics')   

  id_total_ocean_river = register_diag_field('ocean_model','total_ocean_river',                  &
       Time%model_time, 'total liquid river water and calving ice entering ocean', 'kg/sec/1e15',&
       missing_value=missing_value,range=(/-1e6,1e6/))

  id_total_ocean_evap = register_diag_field('ocean_model','total_ocean_evap',  &
       Time%model_time, 'total evaporative ocean mass flux (>0 enters ocean)', &
       '(kg/sec)/1e15', missing_value=missing_value,range=(/-1e6,1e6/))

  id_total_ocean_melt = register_diag_field('ocean_model','total_ocean_melt',       &
       Time%model_time, 'total liquid water melted from sea ice (>0 enters ocean)', &
       '(kg/sec)/1e15', missing_value=missing_value,range=(/-1e6,1e6/))

  id_total_ocean_pme_river = register_diag_field('ocean_model','total_ocean_pme_river',        &
       Time%model_time, 'total ocean precip-evap+river via sbc (liquid, frozen, evaporation)', &
       '(kg/sec)/1e15' , missing_value=missing_value,range=(/-1e6,1e6/))

  id_total_ocean_pme_sbc = register_diag_field('ocean_model','total_ocean_pme_sbc',      &
       Time%model_time, 'total ocean precip-evap via sbc (liquid, frozen, evaporation)', &
       'kg/sec/1e15' , missing_value=missing_value,range=(/-1e6,1e6/))

  id_total_ocean_pme_restore = register_diag_field('ocean_model','total_ocean_pme_restore',&
       Time%model_time, 'total precip-evap from restore (>0 enters ocean)',                &
       '(kg/sec)/1e15', missing_value=missing_value,range=(/-1e6,1e6/))   

  id_total_ocean_pme_correct = register_diag_field('ocean_model','total_ocean_pme_correct',&
       Time%model_time, 'total precip-evap from flux correction (>0 enters ocean)',        &
       '(kg/sec)/1e15', missing_value=missing_value,range=(/-1e6,1e6/))   

  id_total_ocean_pme_net = register_diag_field('ocean_model','total_ocean_pme_net',    &
       Time%model_time, 'total precip-evap into ocean (total w/ restore + normalize)', &
       '(kg/sec)/1e15', missing_value=missing_value,range=(/-1e6,1e6/))   

  id_total_ocean_fprec = register_diag_field('ocean_model','total_ocean_fprec', &
       Time%model_time, 'total snow falling onto ocean (>0 enters ocean)',      &
       '(kg/sec)/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

  id_total_ocean_lprec = register_diag_field('ocean_model','total_ocean_lprec', &
       Time%model_time, 'total liquid precip into ocean (>0 enters ocean)',     &
       '(kg/sec)/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

#if defined(ACCESS_CM) || defined(ACCESS_OM)
  id_aice = register_diag_field('ocean_model','aice', Grd%tracer_axes(1:2),&
       Time%model_time, 'fraction of surface area covered with ice', 'm^2/m^2' ,  &
       missing_value=missing_value,range=(/-1.e1,1.e1/),                      &
       standard_name='areal_ice_concentration' )
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
  id_iof_nit = register_diag_field('ocean_model','iof_nit', Grd%tracer_axes(1:2),&
       Time%model_time, 'ice-ocean flux of nitrate', 'mmol/m^2/s^1' ,  &
       missing_value=missing_value,range=(/-1.e1,1.e1/),                      &
       standard_name='ice_ocean_nitrate_flux' )
  id_iof_alg = register_diag_field('ocean_model','iof_alg', Grd%tracer_axes(1:2),&
       Time%model_time, 'ice-ocean flux of algae', 'mmol/m^2/s^1' ,  &
       missing_value=missing_value,range=(/-1.e1,1.e1/),                      &
       standard_name='ice_ocean_algal_flux' )
#endif
  id_wnd = register_diag_field('ocean_model','wnd', Grd%tracer_axes(1:2),&
       Time%model_time, 'Wind speed', 'm/s' ,  &
       missing_value=missing_value,range=(/-1.e3,1.e3/),                      &
       standard_name='wind_speed' )
  id_total_ocean_wfimelt = register_diag_field('ocean_model','total_ocean_wfimelt',  &
       Time%model_time, 'total icemelt into ocean (>0 enters ocean)',     &
       'kg/sec/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))
  id_total_ocean_wfiform = register_diag_field('ocean_model','total_ocean_wfiform',  &
       Time%model_time, 'total iceform outof ocean (>0 enters ocean)',     &
       'kg/sec/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))
#if defined(ACCESS_CM)
  id_atm_co2 = register_diag_field('ocean_model','atm_co2', Grd%tracer_axes(1:2),&
       Time%model_time, 'Atmospheric CO2 content', 'ppm' ,  &
       missing_value=missing_value,range=(/-1.e1,1.e4/),                      &
       standard_name='atmospheric_co2' )
#endif
  id_total_ocean_licefw = register_diag_field('ocean_model','total_ocean_licefw',  &
        Time%model_time, 'total land icemelt into ocean (>0 enters ocean)',     &
        'kg/sec/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))
  id_total_ocean_liceht = register_diag_field('ocean_model','total_ocean_liceht',  &
        Time%model_time, 'total land icemelt heat flux into ocean (>0 heats ocean)', &
        'Watts/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))
  id_total_ocean_mh_flux = register_diag_field('ocean_model','total_ocean_mh_flux',  &
        Time%model_time, 'total heat flux into ocean from melting ice (>0 heats ocean)', &
        'Watts/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))
#endif

  id_total_ocean_runoff = register_diag_field('ocean_model','total_ocean_runoff',&
       Time%model_time, 'total liquid river runoff (>0 water enters ocean)',     &
       '(kg/sec)/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

  id_total_ocean_calving = register_diag_field('ocean_model','total_ocean_calving',&
       Time%model_time, 'total water entering ocean from calving land ice',        &
       '(kg/sec)/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

  id_total_ocean_swflx = register_diag_field('ocean_model','total_ocean_swflx',&
       Time%model_time, 'total shortwave flux into ocean (>0 heats ocean)',    &
       'Watts/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

  id_total_ocean_swflx_vis = register_diag_field('ocean_model','total_ocean_swflx_vis',&
       Time%model_time, 'total visible shortwave into ocean (>0 heats ocean)',         &
       'Watts/1e15' , missing_value=missing_value,range=(/-1.e10,1.e10/))   

  id_total_ocean_evap_heat = register_diag_field('ocean_model','total_ocean_evap_heat', &
       Time%model_time, 'total latent heat flux into ocean (<0 cools ocean)',           &
       'Watts/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))

  id_total_ocean_lw_heat = register_diag_field('ocean_model','total_ocean_lw_heat', &
       Time%model_time, 'total longwave flux into ocean (<0 cools ocean)',          &
       'Watts/1e15',missing_value=missing_value,range=(/-1.e10,1.e10/))   

  id_total_ocean_sens_heat = register_diag_field('ocean_model','total_ocean_sens_heat', &
       Time%model_time, 'total sensible heat into ocean (<0 cools ocean)',              &
       'Watts/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

  id_total_ocean_fprec_melt_heat = register_diag_field('ocean_model','total_ocean_fprec_melt_heat',&
       Time%model_time, 'total heat flux to melt frozen precip (<0 cools ocean)',                  &
       'Watts/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

  id_total_ocean_calving_melt_heat = register_diag_field('ocean_model','total_ocean_calving_melt_heat',&
       Time%model_time, 'total heat flux to melt frozen land ice (<0 cools ocean)',                    &
       'Watts/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

  id_total_ocean_river_heat = register_diag_field('ocean_model','total_ocean_river_heat',      &
       Time%model_time, 'total heat flux into ocean from liquid+solid runoff (<0 cools ocean)',&
       'Watts/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

  id_eta_tend_sw = register_diag_field('ocean_model','eta_tend_sw', &
       Grd%tracer_axes(1:2), Time%model_time,                       &
      'non-Bouss steric sea level tendency from shortwave heating', &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_sw > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_lw = register_diag_field('ocean_model','eta_tend_lw', &
       Grd%tracer_axes(1:2), Time%model_time,                       &
       'non-Bouss steric sea level tendency from longwave heating', &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_lw > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_sens = register_diag_field('ocean_model','eta_tend_sens', &
       Grd%tracer_axes(1:2), Time%model_time,                           &
       'non-Bouss steric sea level tendency from sensible heating',     & 
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_sens > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_evap_heat = register_diag_field('ocean_model','eta_tend_evap_heat', &
       Grd%tracer_axes(1:2), Time%model_time,                                     &
       'non-Bouss steric sea level tendency from latent heat of vaporization',    &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_evap_heat > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_fprec_melt = register_diag_field('ocean_model','eta_tend_fprec_melt',     &
       Grd%tracer_axes(1:2), Time%model_time,                                           &
       'non-Bouss steric sea level tendency from latent heat of melting frozen precip', &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_fprec_melt > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_iceberg_melt = register_diag_field('ocean_model','eta_tend_iceberg_melt', &
       Grd%tracer_axes(1:2), Time%model_time,                                           &
       'non-Bouss steric sea level tendency from latent heat of melting icebergs',      &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_iceberg_melt > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_heat_coupler = register_diag_field('ocean_model','eta_tend_heat_coupler',            &
       Grd%tracer_axes(1:2), Time%model_time,                                                      &
       'non-Bouss steric sea level tendency from total sfc heating through coupler (sans frazil)', &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_heat_coupler > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_heat_restore = register_diag_field('ocean_model','eta_tend_heat_restore', &
       Grd%tracer_axes(1:2), Time%model_time,                                           &
       'non-Bouss steric sea level tendency from total sfc heating through restore',    &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_heat_restore > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_salt_coupler = register_diag_field('ocean_model','eta_tend_salt_coupler', &
       Grd%tracer_axes(1:2), Time%model_time,                                           &
       'non-Bouss steric sea level tendency from salt flux entering through coupler',   &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_salt_coupler > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_salt_restore = register_diag_field('ocean_model','eta_tend_salt_restore', &
       Grd%tracer_axes(1:2), Time%model_time,                                           &
       'non-Bouss steric sea level tendency from salt flux entering through restoring', &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_salt_restore > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_fprec = register_diag_field('ocean_model','eta_tend_fprec', &
       Grd%tracer_axes(1:2), Time%model_time,                             &
       'non-Bouss steric sea level tendency from frozen precipitation',   &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_fprec > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_lprec = register_diag_field('ocean_model','eta_tend_lprec', &
       Grd%tracer_axes(1:2), Time%model_time,                             &
       'non-Bouss steric sea level tendency from liquid precipitation',   &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_lprec > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_evap = register_diag_field('ocean_model','eta_tend_evap', &
       Grd%tracer_axes(1:2), Time%model_time,                           &
       'non-Bouss steric sea level tendency from evaporating water',    &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_evap > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_runoff = register_diag_field('ocean_model','eta_tend_runoff', &
       Grd%tracer_axes(1:2), Time%model_time,                               &
       'non-Bouss steric sea level tendency from liquid runoff',            &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_runoff > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_iceberg = register_diag_field('ocean_model','eta_tend_iceberg', &
       Grd%tracer_axes(1:2), Time%model_time,                                 &
       'non-Bouss steric sea level tendency from icebergs', 'm/s',            &
       missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_iceberg > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_water_coupler = register_diag_field('ocean_model','eta_tend_water_coupler', &
       Grd%tracer_axes(1:2), Time%model_time,                                             &
       'non-Bouss sea level tendency from surface water fluxes passed through coupler',   &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_water_coupler > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_water_restore = register_diag_field('ocean_model','eta_tend_water_restore',               &
       Grd%tracer_axes(1:2), Time%model_time,                                                           & 
       'non-Bouss steric sea level tendency from surface water fluxes computed from surface restoring', &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_water_restore > 0) diagnose_sea_level_forcing=.true.   


  id_eta_tend_sw_glob = register_diag_field('ocean_model','eta_tend_sw_glob',   &
       Time%model_time,                                                         &
      'global mean non-Bouss steric sea level tendency from shortwave heating', &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_sw_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_lw_glob = register_diag_field('ocean_model','eta_tend_lw_glob',   &
       Time%model_time,                                                         &
       'global mean non-Bouss steric sea level tendency from longwave heating', &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_lw_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_sens_glob = register_diag_field('ocean_model','eta_tend_sens_glob', &
       Time%model_time,                                                           &
       'global mean non-Bouss steric sea level tendency from sensible heating',   & 
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_sens_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_evap_heat_glob = register_diag_field('ocean_model','eta_tend_evap_heat_glob', &
       Time%model_time,                                                                     &
       'global mean non-Bouss steric sea level tendency from latent heat of vaporization',  &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_evap_heat_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_fprec_melt_glob = register_diag_field('ocean_model','eta_tend_fprec_melt_glob',       &
       Time%model_time,                                                                             &
       'global mean non-Bouss steric sea level tendency from latent heat of melting frozen precip', &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_fprec_melt_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_iceberg_melt_glob = register_diag_field('ocean_model','eta_tend_iceberg_melt_glob', &
       Time%model_time,                                                                           &
       'global mean non-Bouss steric sea level tendency from latent heat of melting icebergs',    &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_iceberg_melt_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_heat_coupler_glob = register_diag_field('ocean_model','eta_tend_heat_coupler_glob',              &
       Time%model_time,                                                                                        &
       'global mean non-Bouss steric sea level tendency from total sfc heating through coupler (sans frazil)', &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_heat_coupler_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_heat_restore_glob = register_diag_field('ocean_model','eta_tend_heat_restore_glob', &
       Time%model_time,                                                                           &
       'global mean non-Bouss steric sea level tendency from total sfc heating through restore',  &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_heat_restore_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_salt_coupler_glob = register_diag_field('ocean_model','eta_tend_salt_coupler_glob', &
       Time%model_time,                                                                           &
       'global mean non-Bouss steric sea level tendency from salt flux entering through coupler', &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_salt_coupler_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_salt_restore_glob = register_diag_field('ocean_model','eta_tend_salt_restore_glob',   &
       Time%model_time,                                                                             &
       'global mean non-Bouss steric sea level tendency from salt flux entering through restoring', &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_salt_restore_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_fprec_glob = register_diag_field('ocean_model','eta_tend_fprec_glob', &
       Time%model_time,                                                             &
       'global mean non-Bouss steric sea level tendency from frozen precipitation', &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_fprec_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_lprec_glob = register_diag_field('ocean_model','eta_tend_lprec_glob', &
       Time%model_time,                                                             &
       'global mean non-Bouss steric sea level tendency from liquid precipitation', &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_lprec_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_evap_glob = register_diag_field('ocean_model','eta_tend_evap_glob', &
       Time%model_time,                                                           &
       'global mean non-Bouss steric sea level tendency from evaporating water',  &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_evap_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_runoff_glob = register_diag_field('ocean_model','eta_tend_runoff_glob', &
       Time%model_time,                                                               &
       'global mean non-Bouss steric sea level tendency from liquid runoff',          &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_runoff_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_iceberg_glob = register_diag_field('ocean_model','eta_tend_iceberg_glob', &
       Time%model_time,                                                                 &
       'global mean non-Bouss steric sea level tendency from icebergs', 'm/s',          &
       missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_iceberg_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_water_coupler_glob = register_diag_field('ocean_model','eta_tend_water_coupler_glob', &
       Time%model_time,                                                                             &
       'global mean non-Bouss sea level tendency from surface water fluxes passed through coupler', &
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_water_coupler_glob > 0) diagnose_sea_level_forcing=.true.   

  id_eta_tend_water_restore_glob = register_diag_field('ocean_model','eta_tend_water_restore_glob',                &
       Time%model_time,                                                                                            & 
       'global mean non-Bouss steric sea level tendency from surface water fluxes computed from surface restoring',&
       'm/s', missing_value=missing_value,range=(/-1.e10,1.e10/))   
  if(id_eta_tend_water_restore_glob > 0) diagnose_sea_level_forcing=.true.   


  id_mass_pme_adj_on_nrho = register_diag_field ('ocean_model',                                &
       'mass_pme_adj_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                     &
       'mass transport from surface restoring pme (>0 enters ocean) binned to neutral density',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_mass_precip_on_nrho = register_diag_field ('ocean_model',                                                &
       'mass_precip_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                                     &
       'mass transport from liquid+solid prec & seaice melt+form (>0 enters ocean) binned to neutral density',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_mass_evap_on_nrho = register_diag_field ('ocean_model',                                              &
       'mass_evap_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                                   &
       'mass transport from evaporation/condensation (>0 leaves ocean) binned to neutral density classes',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_mass_river_on_nrho = register_diag_field ('ocean_model',                                               & 
       'mass_river_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                                    &
       'mass transport from liquid+frozen river runoff (>0 enters ocean) binned to neutral density classes',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_mass_melt_on_nrho = register_diag_field ('ocean_model',                                             & 
       'mass_melt_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                                  &
       'mass transport from melting/forming sea ice (>0 enters ocean) binned to neutral density classes',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_mass_pmepr_on_nrho = register_diag_field ('ocean_model',                                                            & 
       'mass_pmepr_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                                                 &
       'mass transport from liquid+frozen mass and seaice melt+form (>0 enters ocean) binned to neutral density classes',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_tform_rho_pbl_flux_on_nrho = register_diag_field ('ocean_model',            & 
       'tform_rho_pbl_flux_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time, &
       'impact of surface heat and salt flux on water mass transformation',      &
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_tform_rho_pbl_adjheat_on_nrho = register_diag_field ('ocean_model',           & 
       'tform_rho_pbl_adjheat_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,&
       'impact of surface heat restoring+correction on water mass transformation', &
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_tform_rho_pbl_adjsalt_on_nrho = register_diag_field ('ocean_model',           & 
       'tform_rho_pbl_adjsalt_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,&
       'impact of surface salt restoring+correction on water mass transformation', &
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_tform_rho_pbl_heat_on_nrho = register_diag_field ('ocean_model',            & 
       'tform_rho_pbl_heat_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time, &
       'impact of surface heat flux on water mass transformation',               &
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_tform_rho_pbl_salt_on_nrho = register_diag_field ('ocean_model',            & 
       'tform_rho_pbl_salt_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time, &
       'impact of surface salt flux on water mass transformation',               &
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_tform_rho_pbl_sw_on_nrho = register_diag_field ('ocean_model',                     & 
       'tform_rho_pbl_sw_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,          &
       'impact of surface shortwave heat (>0 heats ocean) on water mass transformation',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_tform_rho_pbl_lw_on_nrho = register_diag_field ('ocean_model',                    & 
       'tform_rho_pbl_lw_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,         &
       'impact of surface longwave heat (>0 heats ocean) on water mass transformation',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_tform_rho_pbl_sens_on_nrho = register_diag_field ('ocean_model',                  & 
       'tform_rho_pbl_sens_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,       &
       'impact of surface sensible heat (>0 heats ocean) on water mass transformation',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_tform_rho_pbl_lat_on_nrho = register_diag_field ('ocean_model',                                    &  
       'tform_rho_pbl_lat_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                         &
       'impact of surface latent (vapor and solid) heat (>0 heats ocean) on water mass transformation ',&
       'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))


! register tracer-dependent dynamic fields 
  do n =1,num_prog_tracers  

     if (n == index_temp) then
         id_stf_coupler(n) = register_diag_field('ocean_model','sfc_hflux_coupler', &
              Grd%tracer_axes(1:2),                                                 &
              Time%model_time, 'surface heat flux coming through coupler',          &
              'Watts/m^2' ,                                                         &
              missing_value=missing_value,range=(/-1.e4,1.e4/))
         id_stf_restore(n) = register_diag_field('ocean_model','sfc_hflux_restore', &
              Grd%tracer_axes(1:2),                                                 &
              Time%model_time, 'surface heat flux from restoring',                  &
              'Watts/m^2' ,                                                         &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         id_stf_correct(n) = register_diag_field('ocean_model','sfc_hflux_correct', &
              Grd%tracer_axes(1:2),                                                 &
              Time%model_time, 'surface heat flux from flux correction',            &
              'Watts/m^2' ,                                                         &
              missing_value=missing_value,range=(/-1.e4,1.e4/),                     &
              standard_name='heat_flux_correction')            
         id_stf_total(n) = register_diag_field('ocean_model','sfc_hflux_total',                            &
              Grd%tracer_axes(1:2),                                                                        &
              Time%model_time, 'surface heat flux from coupler plus restore (omits mass transfer heating)',&
              'Watts/m^2' ,                                                                                &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         id_net_sfc_heating = register_diag_field('ocean_model','net_sfc_heating',                &
              Grd%tracer_axes(1:2),                                                               &
              Time%model_time, 'surface ocean heat flux coming through coupler and mass transfer',&
              'Watts/m^2' ,                                                                       &
              missing_value=missing_value,range=(/-1.e4,1.e4/))
         id_net_sfc_workq = register_diag_field('ocean_model','net_sfc_workq',                        &
              Grd%tracer_axes(1:2),                                                               &
              Time%model_time, 'net_sfc_heat*g*alpha2/Cp/rho0, alpha uses pot_rho rather than neut',    &
              'm^2/s^-3' ,                                                                        &
              missing_value=missing_value,range=(/-1.e4,1.e4/))
         id_stf_runoff(n) = register_diag_field('ocean_model','sfc_hflux_from_runoff',&
              Grd%tracer_axes(1:2),                                                   &
              Time%model_time, 'heat flux (relative to 0C) from liquid river runoff', &    
              'Watts/m^2' ,                                                           &
              missing_value=missing_value,range=(/-1.e4,1.e4/),                       &
              standard_name='temperature_flux_due_to_runoff_expressed_as_heat_flux_into_sea_water')            
         id_stf_calving(n) = register_diag_field('ocean_model','sfc_hflux_from_calving',       &
              Grd%tracer_axes(1:2),                                                            &
              Time%model_time, 'heat flux (relative to 0C) from solid land ice entering ocean',&    
              'Watts/m^2' ,                                                             &
              missing_value=missing_value,range=(/-1.e4,1.e4/),                         &
              standard_name='temperature_flux_due_to_icebergs_expressed_as_heat_flux_into_sea_water')            
         id_stf_pme(n) = register_diag_field('ocean_model','sfc_hflux_pme',                                  &
              Grd%tracer_axes(1:2),                                                                          &
              Time%model_time, 'heat flux (relative to 0C) from pme transfer of water across ocean surface', &
              'Watts/m^2' , missing_value=missing_value,range=(/-1.e4,1.e4/))            
         id_stf_pme_on_nrho(n) = register_diag_field('ocean_model','sfc_hflux_pme_on_nrho',                                  &
              Dens%neutralrho_axes(1:3),                                                                      &
              Time%model_time, 'heat flux (relative to 0C) from pme transfer of water across ocean surface binned to neutral density', &
              'Watts/m^2' , missing_value=missing_value,range=(/-1.e20,1.e20/))
         id_stf_prec(n) = register_diag_field('ocean_model','sfc_hflux_from_water_prec',      &
              Grd%tracer_axes(1:2),                                                           &
              Time%model_time, 'heat flux from precip transfer of water across ocean surface',&
              'Watts/m^2' , missing_value=missing_value,range=(/-1.e4,1.e4/),                 &
              standard_name='temperature_flux_due_to_rainfall_expressed_as_heat_flux_into_sea_water')            
         id_stf_evap(n) = register_diag_field('ocean_model','sfc_hflux_from_water_evap',    &
              Grd%tracer_axes(1:2),                                                         &
              Time%model_time, 'heat flux from evap transfer of water across ocean surface',&
              'Watts/m^2' , missing_value=missing_value,range=(/-1.e4,1.e4/),               &
              standard_name='temperature_flux_due_to_evaporation_expressed_as_heat_flux_into_sea_water')            
         id_trunoff(n) = register_diag_field('ocean_model','temp_runoff',              &
              Grd%tracer_axes(1:2),                                                    &
              Time%model_time, 'temperature of liquid river runoff entering the ocean',&
              'degC' , missing_value=missing_value,range=(/-1.e4,1.e4/))            
         id_tcalving(n) = register_diag_field('ocean_model','temp_calving', &
              Grd%tracer_axes(1:2),                                         &
              Time%model_time, 'temperature of land ice calving into ocean',&
              'degC' , missing_value=missing_value,range=(/-1.e4,1.e4/))            
         id_triver(n) = register_diag_field('ocean_model','temp_river',                      &
              Grd%tracer_axes(1:2),                                                          &
              Time%model_time, 'temperature of river water (=runoff+calving) entering ocean',&
              'degC' , missing_value=missing_value,range=(/-1.e4,1.e4/))            

         id_total_ocean_stf_coupler(n) = register_diag_field('ocean_model','total_ocean_hflux_coupler', &
              Time%model_time, 'total surface heat flux passed through coupler', 'Watts/1e15' ,         &
              missing_value=missing_value,range=(/-1.e4,1.e4/))
         id_total_ocean_stf_restore(n) = register_diag_field('ocean_model','total_ocean_hflux_restore', &
              Time%model_time, 'total surface heat flux adjustment from restoring', 'Watts/1e15',       &
              missing_value=missing_value,range=(/-1.e4,1.e4/))
         id_total_ocean_stf_correct(n) = register_diag_field('ocean_model','total_ocean_hflux_correct', &
              Time%model_time, 'total surface heat flux adjustment from correction', 'Watts/1e15',      &
              missing_value=missing_value,range=(/-1.e4,1.e4/))
         id_total_ocean_stf_sum(n) = register_diag_field('ocean_model','total_ocean_hflux_sum',   &
              Time%model_time, 'total surface heat flux from coupler and restoring', 'Watts/1e15',&
              missing_value=missing_value,range=(/-1.e4,1.e4/))
         id_total_net_sfc_heating = register_diag_field('ocean_model','total_net_sfc_heating',          &
              Time%model_time, 'total ocean surface flux from coupler and mass transfer', 'Watts/1e15' ,&
              missing_value=missing_value,range=(/-1.e4,1.e4/))
         id_total_ocean_stf_runoff(n) = register_diag_field('ocean_model','total_ocean_runoff_heat',&
              Time%model_time, 'total ocean heat flux from liquid river runoff', 'Watts/1e15',      &
              missing_value=missing_value,range=(/-1.e4,1.e4/))
         id_total_ocean_stf_calving(n) = register_diag_field('ocean_model','total_ocean_calving_heat',&
              Time%model_time, 'total ocean heat flux from calving land ice', 'Watts/1e15',           &
              missing_value=missing_value,range=(/-1.e4,1.e4/))
         id_total_ocean_stf_pme(n) = register_diag_field('ocean_model','total_ocean_hflux_pme',   &
              Time%model_time, 'total ocean heat flux from pme transferring water across surface',&
              'Watts/1e15', missing_value=missing_value,range=(/-1.e4,1.e4/))
         id_total_ocean_stf_prec(n) = register_diag_field('ocean_model','total_ocean_hflux_prec',    &
              Time%model_time, 'total ocean heat flux from precip transferring water across surface',&
              'Watts/1e15', missing_value=missing_value,range=(/-1.e4,1.e4/))
         id_total_ocean_stf_evap(n) = register_diag_field('ocean_model','total_ocean_hflux_evap',  &
              Time%model_time, 'total ocean heat flux from evap transferring water across surface',&
              'Watts/1e15', missing_value=missing_value,range=(/-1.e4,1.e4/))

     elseif(n == index_salt) then 

         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_coupler'
         id_stf_coupler(n) = register_diag_field('ocean_model',trim(name), &
              Grd%tracer_axes(1:2),                                        &
              Time%model_time, trim(name)//': flux from the coupler',      &
              'kg/(m^2*sec)' ,                                             &
              missing_value=missing_value,range=(/-1.e4,1.e4/))
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_restore'
         id_stf_restore(n) = register_diag_field('ocean_model',         &
              trim(name),                                               &
              Grd%tracer_axes(1:2),                                     &
              Time%model_time, trim(name)//': flux from restoring term',&
              'kg/(m^2*sec)' ,                                          &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_correct'
         id_stf_correct(n) = register_diag_field('ocean_model',                &
              trim(name),                                                      &
              Grd%tracer_axes(1:2),                                            &
              Time%model_time, trim(name)//': flux correction from data file', &
              'kg/(m^2*sec)' ,                                                 &
              missing_value=missing_value,range=(/-1.e4,1.e4/),                &
              standard_name='virtual_salt_flux_correction')            
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_total'
         id_stf_total(n) = register_diag_field('ocean_model',&
              trim(name),                                    &
              Grd%tracer_axes(1:2),                          &
              Time%model_time, trim(name),                   &
              'kg/(m^2*sec)' ,                               &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_runoff'
         id_stf_runoff(n) = register_diag_field('ocean_model',&
              trim(name),                                     &
              Grd%tracer_axes(1:2),                           &
              Time%model_time, trim(name),                    &
              'kg/(m^2*sec)' ,                                &
              missing_value=missing_value,range=(/-1.e4,1.e4/), &
              standard_name='salt_flux_into_sea_water_from_rivers')            
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_calving'
         id_stf_calving(n) = register_diag_field('ocean_model',&
              trim(name),                                      &
              Grd%tracer_axes(1:2),                            &
              Time%model_time, trim(name),                     &
              'kg/(m^2*sec)' ,&
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_ice'
         id_salt_flux_ice = register_diag_field('ocean_model',&
              trim(name),                                     &
              Grd%tracer_axes(1:2),                           &
              Time%model_time, trim(name),                    &
              'kg/(m^2*sec)' ,                                &
              missing_value=missing_value,range=(/-1.e4,1.e4/), &
              standard_name='downward_sea_ice_basal_salt_flux')            
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_pme'
         id_stf_pme(n) = register_diag_field('ocean_model',&
              trim(name),                                  &
              Grd%tracer_axes(1:2),                        &
              Time%model_time, trim(name),                 &
              'kg/(m^2*sec)' ,                             &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_pme_on_nrho'
         id_stf_pme_on_nrho(n) = register_diag_field('ocean_model',&
              trim(name),                                  &
              Dens%neutralrho_axes(1:3),                    &
              Time%model_time, trim(name),                 &
              'kg/(m^2*sec)' ,                             &
              missing_value=missing_value,range=(/-1.e20,1.e20/))
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_prec'
         id_stf_prec(n) = register_diag_field('ocean_model',&
              trim(name),                                   & 
              Grd%tracer_axes(1:2),                         &
              Time%model_time, trim(name),                  &
              'kg/(m^2*sec)' ,                              &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_evap'
         id_stf_evap(n) = register_diag_field('ocean_model',&
              trim(name),                                   &
              Grd%tracer_axes(1:2),                         &
              Time%model_time, trim(name),                  &
              'kg/(m^2*sec)' ,                              & 
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = trim(T_prog(n)%name)//'_runoff'
         id_trunoff(n) = register_diag_field('ocean_model',&
              trim(name),                                  &
              Grd%tracer_axes(1:2),                        &
              Time%model_time, trim(name),                 &
              trim(T_prog(n)%units),                       &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = trim(T_prog(n)%name)//'_calving'
         id_tcalving(n) = register_diag_field('ocean_model',&
              trim(name),                                   &
              Grd%tracer_axes(1:2),                         &
              Time%model_time, trim(name),                  &
              trim(T_prog(n)%units),                        &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = trim(T_prog(n)%name)//'_river'
         id_triver(n) = register_diag_field('ocean_model',&
              trim(name),                                 &
              Grd%tracer_axes(1:2),                       &
              Time%model_time, trim(name),                &
              trim(T_prog(n)%units),                      &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            

         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_ice'
         id_total_salt_flux_ice =  register_diag_field('ocean_model',    &
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',&
              missing_value=missing_value,range=(/-1.e4,1.e4/))
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_coupler'
         id_total_ocean_stf_coupler(n) =  register_diag_field('ocean_model',&
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',   &
              missing_value=missing_value,range=(/-1.e4,1.e4/))
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_restore'
         id_total_ocean_stf_restore(n) = register_diag_field('ocean_model',&
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' , &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_correct'
         id_total_ocean_stf_correct(n) = register_diag_field('ocean_model',&
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' , &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_sum'
         id_total_ocean_stf_sum(n) = register_diag_field('ocean_model',  &
              trim(name),Time%model_time, trim(name), 'kg/sec (*1e-15)', &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_runoff'
         id_total_ocean_stf_runoff(n) = register_diag_field('ocean_model',&
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' ,&
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_calving'
         id_total_ocean_stf_calving(n) = register_diag_field('ocean_model',&
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' , &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_pme'
         id_total_ocean_stf_pme(n) = register_diag_field('ocean_model',  &
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',&
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_prec'
         id_total_ocean_stf_prec(n) = register_diag_field('ocean_model', &
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',&
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_evap'
         id_total_ocean_stf_evap(n) = register_diag_field('ocean_model', &
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',&
              missing_value=missing_value,range=(/-1.e4,1.e4/))            

     else

         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_coupler'
         id_stf_coupler(n) = register_diag_field('ocean_model',trim(name), &
              Grd%tracer_axes(1:2),                                        &
              Time%model_time, trim(name)//': flux from the coupler',      &
              'kg/(m^2*sec)' ,                                             &
              missing_value=missing_value,range=(/-1.e4,1.e4/))
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_restore'
         id_stf_restore(n) = register_diag_field('ocean_model',        &
              trim(name),                                              &
              Grd%tracer_axes(1:2),                                    &
              Time%model_time, trim(name)//': flux from the restoring',&
              'kg/(m^2*sec)' ,                                         &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_correct'
         id_stf_correct(n) = register_diag_field('ocean_model',          &
              trim(name),                                                &
              Grd%tracer_axes(1:2),                                      &
              Time%model_time, trim(name)//': flux from flux correction',&
              'kg/(m^2*sec)' ,                                           &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_total'
         id_stf_total(n) = register_diag_field('ocean_model',&
              trim(name),                                    &
              Grd%tracer_axes(1:2),                          &
              Time%model_time, trim(name),                   &
              'kg/(m^2*sec)' ,&
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_runoff'
         id_stf_runoff(n) = register_diag_field('ocean_model',&
              trim(name),                                     &
              Grd%tracer_axes(1:2),                           &
              Time%model_time, trim(name),                    &
              'kg/(m^2*sec)' ,                                &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_calving'
         id_stf_calving(n) = register_diag_field('ocean_model',&
              trim(name),                                      &
              Grd%tracer_axes(1:2),                            &
              Time%model_time, trim(name),                     &
              'kg/(m^2*sec)' ,                                 &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_pme'
         id_stf_pme(n) = register_diag_field('ocean_model',&
              trim(name),                                  &
              Grd%tracer_axes(1:2),                        &
              Time%model_time, trim(name),                 &
              'kg/(m^2*sec)' ,                             &
              missing_value=missing_value,range=(/-1.e4,1.e4/))  
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_pme_on_nrho'
         id_stf_pme_on_nrho(n) = register_diag_field('ocean_model',&
              trim(name),                                  &
              Dens%neutralrho_axes(1:3),                    &
              Time%model_time, trim(name),                 &
              'kg/(m^2*sec)' ,                             &
              missing_value=missing_value,range=(/-1.e20,1.e20/))
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_prec'
         id_stf_prec(n) = register_diag_field('ocean_model',&
              trim(name),                                   &
              Grd%tracer_axes(1:2),                         &
              Time%model_time, trim(name),                  &
              'kg/(m^2*sec)' ,                              &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'sfc_'//trim(T_prog(n)%name)//'_flux_evap'
         id_stf_evap(n) = register_diag_field('ocean_model',&
              trim(name),                                   &
              Grd%tracer_axes(1:2),                         &
              Time%model_time, trim(name),                  &
              'kg/(m^2*sec)' ,                              &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = trim(T_prog(n)%name)//'_runoff'
         id_trunoff(n) = register_diag_field('ocean_model',&
              trim(name),                                  &
              Grd%tracer_axes(1:2),                        &
              Time%model_time, trim(name),                 &
              trim(T_prog(n)%units),                       &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = trim(T_prog(n)%name)//'_calving'
         id_tcalving(n) = register_diag_field('ocean_model',&
              trim(name),                                   &
              Grd%tracer_axes(1:2),                         &
              Time%model_time, trim(name),                  &
              trim(T_prog(n)%units),                        &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = trim(T_prog(n)%name)//'_river'
         id_triver(n) = register_diag_field('ocean_model',&
              trim(name),                                 &
              Grd%tracer_axes(1:2),                       &
              Time%model_time, trim(name),                &
              trim(T_prog(n)%units),                      &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            

         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_coupler'
         id_total_ocean_stf_coupler(n) =  register_diag_field('ocean_model',&
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',   &
              missing_value=missing_value,range=(/-1.e4,1.e4/))
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_restore'
         id_total_ocean_stf_restore(n) = register_diag_field('ocean_model',&
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' , &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_correct'
         id_total_ocean_stf_correct(n) = register_diag_field('ocean_model',&
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' , &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_sum'
         id_total_ocean_stf_sum(n) = register_diag_field('ocean_model', &
              trim(name),Time%model_time, trim(name), 'kg/sec (*1e-15)',&
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_runoff'
         id_total_ocean_stf_runoff(n) = register_diag_field('ocean_model',&
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' ,&
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_calving'
         id_total_ocean_stf_calving(n) = register_diag_field('ocean_model',&
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' , &
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_pme'
         id_total_ocean_stf_pme(n) = register_diag_field('ocean_model',   &
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' ,&
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_prec'
         id_total_ocean_stf_prec(n) = register_diag_field('ocean_model', &
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',&
              missing_value=missing_value,range=(/-1.e4,1.e4/))            
         name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_evap'
         id_total_ocean_stf_evap(n) = register_diag_field('ocean_model', &
              trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',&
              missing_value=missing_value,range=(/-1.e4,1.e4/))            

     endif

  enddo   ! enddo for num_prog_tracer loop 


  if(.not. diagnose_sea_level_forcing) then 
      write(stdoutunit,*) &
           '==>Note that diagnose_sea_level_forcing==.false., so no "eta_tend_" fields will be diagnosed.'
  endif


  return 

end subroutine ocean_sbc_diag_init
! </SUBROUTINE> NAME="ocean_sbc_diag_init"


!#######################################################################
! <SUBROUTINE NAME="initialize_ocean_sfc">
!
! <DESCRIPTION>
! Initialize the ocean surface type, which passes information between ocean 
! and other component models. 
!
! Note that ocean model sst passed to the atmosphere must be the surface
! potential temperature (which is equated to surface in situ temperature).
! If the ocean prognostic temperature variable is conservative temperature,
! then the sst is carried in T_diag(index_diag_temp).  If the prognostic 
! temperature is potential temperature, then the sst is carried in 
! T_prog(index_temp). 
!
! Note that we assume the winds passed to the ocean are on the B-grid 
! velocity point.  Likewise, we pass the currents back to the coupler 
! on the B-grid point.  Code will need to be modified if using another
! coupler assumption.  
!
!  Ocean_sfc%t_surf  = time averaged sst (Kelvin) passed to atmosphere/ice model
!  Ocean_sfc%s_surf  = time averaged sss (psu) passed to atmosphere/ice models
!  Ocean_sfc%u_surf  = time averaged u-current (m/sec) passed to atmosphere/ice models
!  Ocean_sfc%v_surf  = time averaged v-current (m/sec)  passed to atmosphere/ice models 
!  Ocean_sfc%sea_lev = time averaged ocean free surface height (m) plus patm/(grav*rho0) - h_geoid - h_tide 
!  Ocean_sfc%frazil  = time accumulated frazil (J/m^2) passed to ice model.  time averaging 
!                      not performed, since ice model needs the frazil accumulated over the 
!                      ocean time steps.  Note that Ocean_sfc%frazil is accumulated, whereas 
!                      T_diag%frazil (saved in diagnostic tracer restart file) is instantaneous. 
!
! </DESCRIPTION>
!
subroutine initialize_ocean_sfc(Time, Thickness, T_prog, T_diag, Velocity, Ocean_sfc)

  type(ocean_time_type),         intent(in)    :: Time
  type(ocean_thickness_type),    intent(in)    :: Thickness
  type(ocean_prog_tracer_type),  intent(in)    :: T_prog(:)
  type(ocean_diag_tracer_type),  intent(in)    :: T_diag(:)
  type(ocean_velocity_type),     intent(in)    :: Velocity
  type(ocean_public_type), target, intent(inout) :: Ocean_sfc

  integer                          :: taup1, id_field
  integer                          :: i, ii, j, jj
  real, dimension(isd:ied,jsd:jed) :: sst
  real, dimension(isd:ied,jsd:jed) :: sst_redist
  character(len=128)               :: filename

  integer :: stdoutunit
  stdoutunit=stdout()

  taup1 = Time%taup1
 
  Ocean_sfc%avg_kount = 0
  Ocean_sfc%t_surf    = kelvin
  sst                 = 0.0
  sst_redist          = 0.0

  if(prog_temp_variable==CONSERVATIVE_TEMP) then
    sst(:,:) = T_diag(index_diag_temp)%field(:,:,1)
    if(index_redist_heat > 0 ) sst_redist(:,:) = pottemp_from_contemp(T_prog(index_salt)%field(:,:,1,taup1),  &
                                                                      T_prog(index_redist_heat)%field(:,:,1,taup1))
  else
    sst(:,:) = T_prog(index_temp)%field(:,:,1,taup1)
    if(index_redist_heat > 0) sst_redist(:,:) = T_prog(index_redist_heat)%field(:,:,1,taup1)
  endif

  where (Grd%tmask(isc:iec,jsc:jec,1) == 1.0)
      Ocean_sfc%t_surf(isc_bnd:iec_bnd,jsc_bnd:jec_bnd)  = sst(isc:iec,jsc:jec) + kelvin
      Ocean_sfc%s_surf(isc_bnd:iec_bnd,jsc_bnd:jec_bnd)  = T_prog(index_salt)%field(isc:iec,jsc:jec,1,taup1)
      Ocean_sfc%u_surf(isc_bnd:iec_bnd,jsc_bnd:jec_bnd)  = Velocity%u(isc:iec,jsc:jec,1,1,taup1)
      Ocean_sfc%v_surf(isc_bnd:iec_bnd,jsc_bnd:jec_bnd)  = Velocity%u(isc:iec,jsc:jec,1,2,taup1)
      Ocean_sfc%sea_lev(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) = Thickness%sea_lev(isc:iec,jsc:jec)
      Ocean_sfc%frazil(isc_bnd:iec_bnd,jsc_bnd:jec_bnd)  = 0.0
#if defined(ACCESS_CM)
      Ocean_sfc%co2flux(isc_bnd:iec_bnd,jsc_bnd:jec_bnd)     = co2flux(isc:iec,jsc:jec)  !These should really come from a restart RASF
      Ocean_sfc%co2(isc_bnd:iec_bnd,jsc_bnd:jec_bnd)         = ocn_co2(isc:iec,jsc:jec)
#endif
  end where

#if defined(ACCESS_OM) && defined(CSIRO_BGC)
  if (ind_no3 > 0) then
   where (Grd%tmask(isc:iec,jsc:jec,1) == 1.0)
      Ocean_sfc%n_surf(isc_bnd:iec_bnd,jsc_bnd:jec_bnd)  = T_prog(ind_no3)%field(isc:iec,jsc:jec,1,taup1)
      Ocean_sfc%alg_surf(isc_bnd:iec_bnd,jsc_bnd:jec_bnd)  = T_prog(ind_phy)%field(isc:iec,jsc:jec,1,taup1)
   end where
  end if
#endif

  ! when enabled, use FAFMIP redistributed heat tracer for sst
  if(index_redist_heat > 0) then
    write(stdoutunit,*) &
       '==>Note from ocean_sbc_mod (initialize_ocean_sfc): FAFMIP - Using redistributed heat tracer for sst.'
    Ocean_sfc%t_surf    = kelvin
    where (Grd%tmask(isc:iec,jsc:jec,1) == 1.0)
      Ocean_sfc%t_surf(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) = sst_redist(isc:iec,jsc:jec) + kelvin
    end where
  endif


  ! coupler works with b-grid velocity, so we average c-grid to get b-grid
  if(horz_grid == MOM_CGRID) then 
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift
            Ocean_sfc%u_surf(i,j)  = Grd%umask(ii,jj,1)*onehalf &
                                     *(Velocity%u(ii,jj,1,1,taup1)+Velocity%u(ii,jj+1,1,1,taup1))
            Ocean_sfc%v_surf(i,j)  = Grd%umask(ii,jj,1)*onehalf &
                                     *(Velocity%u(ii,jj,1,2,taup1)+Velocity%u(ii+1,jj,1,2,taup1))
         enddo
      enddo
  endif

  filename = 'ocean_sbc.res.nc'
  id_field = register_restart_field(Sfc_restart, filename, 't_surf', Ocean_sfc%t_surf,Ocean_sfc%Domain)
  id_field = register_restart_field(Sfc_restart, filename, 's_surf', Ocean_sfc%s_surf,Ocean_sfc%Domain)
  id_field = register_restart_field(Sfc_restart, filename, 'u_surf', Ocean_sfc%u_surf,Ocean_sfc%Domain)
  id_field = register_restart_field(Sfc_restart, filename, 'v_surf', Ocean_sfc%v_surf,Ocean_sfc%Domain)
  id_field = register_restart_field(Sfc_restart, filename, 'sea_lev',Ocean_sfc%sea_lev,Ocean_sfc%Domain)
  id_field = register_restart_field(Sfc_restart, filename, 'frazil', Ocean_sfc%frazil,Ocean_sfc%Domain)
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
  id_field = register_restart_field(Sfc_restart, filename, 'n_surf', Ocean_sfc%n_surf,Ocean_sfc%Domain)
  id_field = register_restart_field(Sfc_restart, filename, 'alg_surf', Ocean_sfc%alg_surf,Ocean_sfc%Domain)
#endif
#if defined(ACCESS_CM)
!RASF Make these optional so we don't break existing runs.
  id_field = register_restart_field(Sfc_restart, filename, 'co2flux',Ocean_sfc%co2flux,Ocean_sfc%Domain, mandatory=.false.)
  id_field = register_restart_field(Sfc_restart, filename, 'ocn_co2', Ocean_sfc%co2,Ocean_sfc%Domain, mandatory=.false.)
#endif
  if (file_exist('INPUT/ocean_sbc.res.nc')) then
     call restore_state(Sfc_restart)
  endif

  return

end subroutine initialize_ocean_sfc
! </SUBROUTINE> NAME="initialize_ocean_sfc"


!#######################################################################
! <SUBROUTINE NAME="sum_ocean_sfc">
!
! <DESCRIPTION>
! Accumulate the ocean_sfc derived type over the course of the 
! ocean component sub-cycling used when coupling to other models. 
!
! Note that ocean model sst passed to the atmosphere must be the surface
! potential temperature (which is equated to surface in situ temperature).
! If the ocean prognostic temperature variable is conservative temperature,
! then the sst is carried in T_diag(index_diag_temp).  If the prognostic 
! temperature is potential temperature, then the sst is carried in 
! T_prog(index_temp). 
!
! Note that this routine is called after eta_and_pbot_diagnose,
! so Thickness%eta is eta_t(taup1).  
!
! </DESCRIPTION>
!  
subroutine sum_ocean_sfc(Time, Thickness, T_prog, T_diag, Dens, Velocity, Ocean_sfc)

#if defined(ACCESS_CM) || defined(ACCESS_OM)
  use ocean_operators_mod, only : GRAD_BAROTROPIC_P 
#endif

  type(ocean_time_type),         intent(in)    :: Time
  type(ocean_thickness_type),    intent(in)    :: Thickness
  type(ocean_prog_tracer_type),  intent(in)    :: T_prog(:)
  type(ocean_diag_tracer_type),  intent(in)    :: T_diag(:)
  type(ocean_density_type),      intent(in)    :: Dens
  type(ocean_velocity_type),     intent(in)    :: Velocity
  type(ocean_public_type), target, intent(inout) :: Ocean_sfc
  
  real, dimension(isd:ied,jsd:jed) :: sst
  real, dimension(isd:ied,jsd:jed) :: sst_redist
  integer                          :: taup1, i, j, k, ii, jj
  integer :: stdoutunit
  stdoutunit=stdout()

  if (Ocean_sfc%avg_kount == 0) call zero_ocean_sfc(Ocean_sfc)

  sst        = 0.0
  sst_redist = 0.0
  taup1 = Time%taup1
  Ocean_sfc%avg_kount = Ocean_sfc%avg_kount + 1

  if(prog_temp_variable==CONSERVATIVE_TEMP) then
    sst(:,:) = T_diag(index_diag_temp)%field(:,:,1)
    if(index_redist_heat > 0 ) sst_redist(:,:) = pottemp_from_contemp(T_prog(index_salt)%field(:,:,1,taup1),  &
                                                                      T_prog(index_redist_heat)%field(:,:,1,taup1))
  else
    sst(:,:) = T_prog(index_temp)%field(:,:,1,taup1)
    if(index_redist_heat > 0) sst_redist(:,:) = T_prog(index_redist_heat)%field(:,:,1,taup1)
  endif

#if defined(ACCESS_CM) || defined(ACCESS_OM)
  sslope(:,:,:) = GRAD_BAROTROPIC_P(Thickness%sea_lev(:,:), isc - isd, isc - isd)
  if (isc - isd /= 1) then
    call mpp_error (FATAL, '==>Error from ocean_sbc_mod (sum_ocean_sfc): grad barotropic halos')
  endif 
#endif

  if(horz_grid == MOM_BGRID) then 
  
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift
            Ocean_sfc%s_surf(i,j)  = Ocean_sfc%s_surf(i,j)  + T_prog(index_salt)%field(ii,jj,1,taup1)
            Ocean_sfc%u_surf(i,j)  = Ocean_sfc%u_surf(i,j)  + Velocity%u(ii,jj,1,1,taup1)
            Ocean_sfc%v_surf(i,j)  = Ocean_sfc%v_surf(i,j)  + Velocity%u(ii,jj,1,2,taup1)
            Ocean_sfc%sea_lev(i,j) = Ocean_sfc%sea_lev(i,j) + Thickness%sea_lev(ii,jj)
#if defined(ACCESS_CM) || defined(ACCESS_OM)
            Ocean_sfc%gradient(i,j,:) = Ocean_sfc%gradient(i,j,:) + sslope(ii,jj,:)  
#endif
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
            if (ind_no3 > 0) then
             Ocean_sfc%n_surf(i,j)  = Ocean_sfc%n_surf(i,j)  + T_prog(ind_no3)%field(ii,jj,1,taup1)
             Ocean_sfc%alg_surf(i,j)  = Ocean_sfc%alg_surf(i,j)  + T_prog(ind_phy)%field(ii,jj,1,taup1)
            end if
#endif
#if defined(ACCESS_CM)
            Ocean_sfc%co2flux(i,j) = Ocean_sfc%co2flux(i,j) + co2flux(ii,jj)
            Ocean_sfc%co2(i,j)     = Ocean_sfc%co2(i,j) + ocn_co2(ii,jj)
#endif
         enddo
      enddo

      ! for FAFMIP, sst is from redistributed heat tracer 
      if(index_redist_heat > 0) then
        do j = jsc_bnd,jec_bnd
           do i = isc_bnd,iec_bnd
              ii = i + i_shift
              jj = j + j_shift
              Ocean_sfc%t_surf(i,j)  = Ocean_sfc%t_surf(i,j) + sst_redist(ii,jj) 
           enddo
        enddo
      else! if not doing FAFMIP use the standard sst
        do j = jsc_bnd,jec_bnd
           do i = isc_bnd,iec_bnd
              ii = i + i_shift
              jj = j + j_shift
              Ocean_sfc%t_surf(i,j)  = Ocean_sfc%t_surf(i,j) + sst(ii,jj)
           enddo
        enddo
       endif

  else  ! cgrid 

      ! coupler works with b-grid velocity, so we average c-grid to get b-grid
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift
            Ocean_sfc%t_surf(i,j)  = Ocean_sfc%t_surf(i,j)  + sst(ii,jj) 
            Ocean_sfc%s_surf(i,j)  = Ocean_sfc%s_surf(i,j)  + T_prog(index_salt)%field(ii,jj,1,taup1)
            Ocean_sfc%u_surf(i,j)  = Ocean_sfc%u_surf(i,j)  + Grd%umask(ii,jj,1)*onehalf & 
                                     *(Velocity%u(ii,jj,1,1,taup1) + Velocity%u(ii,jj+1,1,1,taup1))
            Ocean_sfc%v_surf(i,j)  = Ocean_sfc%v_surf(i,j)  + Grd%umask(ii,jj,1)*onehalf &
                                     *(Velocity%u(ii,jj,1,2,taup1) + Velocity%u(ii,jj+1,1,2,taup1))
            Ocean_sfc%sea_lev(i,j) = Ocean_sfc%sea_lev(i,j) + Thickness%sea_lev(ii,jj)
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
            if (ind_no3 > 0) then
             Ocean_sfc%n_surf(i,j)  = Ocean_sfc%n_surf(i,j)  + T_prog(ind_no3)%field(ii,jj,1,taup1)
             Ocean_sfc%alg_surf(i,j)  = Ocean_sfc%alg_surf(i,j)  + T_prog(ind_phy)%field(ii,jj,1,taup1)
            end if
#endif         
         enddo
      enddo

  endif 

!  In FAFMIP heat experiments we should be using the heat from the redistributed
!  tracer. Previous code uncontionally used the standard tracer. We allow the
!  user to override the correct treatment to recover old results.

  if(index_frazil_redist > 0  .and. do_frazil_redist ) then 
     do k=1,nk
        do j = jsc_bnd,jec_bnd
           do i = isc_bnd,iec_bnd
              ii = i + i_shift
              jj = j + j_shift
              Ocean_sfc%frazil(i,j) = Ocean_sfc%frazil(i,j) + T_diag(index_frazil_redist)%field(ii,jj,k)
           enddo
        enddo
     enddo
  else if(index_frazil > 0 ) then
     do k=1,nk
        do j = jsc_bnd,jec_bnd
           do i = isc_bnd,iec_bnd
              ii = i + i_shift
              jj = j + j_shift
              Ocean_sfc%frazil(i,j) = Ocean_sfc%frazil(i,j) + T_diag(index_frazil)%field(ii,jj,k)
           enddo
        enddo
     enddo
  endif 

  call ocean_tpm_sum_sfc(Dom, T_prog(:), Dens, Ocean_sfc, Time, Grd, isc_bnd, iec_bnd, jsc_bnd, jec_bnd)
      
end subroutine sum_ocean_sfc
! </SUBROUTINE> NAME="sum_ocean_sfc"


!#######################################################################
! <SUBROUTINE NAME="zero_ocean_sfc">
!
! <DESCRIPTION>
! Zero the elements of the Ocean_sfc derived type.  
! </DESCRIPTION>
! 
subroutine zero_ocean_sfc(Ocean_sfc)

  type(ocean_public_type), intent(inout), target :: Ocean_sfc

  integer :: i, j

  Ocean_sfc%avg_kount = 0

  do j = jsc_bnd, jec_bnd
     do i = isc_bnd, iec_bnd
        Ocean_sfc%t_surf(i,j) = 0.0
        Ocean_sfc%s_surf(i,j) = 0.0
        Ocean_sfc%u_surf(i,j) = 0.0
        Ocean_sfc%v_surf(i,j) = 0.0
        Ocean_sfc%sea_lev(i,j)= 0.0
        Ocean_sfc%frazil(i,j) = 0.0
#if defined(ACCESS_CM) || defined(ACCESS_OM)
        Ocean_sfc%gradient(i,j,:)= 0.0
#endif
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
        Ocean_sfc%n_surf(i,j)= 0.0
        Ocean_sfc%alg_surf(i,j)= 0.0
#endif
#if defined(ACCESS_CM)
        Ocean_sfc%co2flux(i,j) = 0.0
        Ocean_sfc%co2(i,j)     = 0.0
#endif
     enddo
  enddo

  call ocean_tpm_zero_sfc(Ocean_sfc)

end subroutine zero_ocean_sfc
! </SUBROUTINE> NAME="zero_ocean_sfc"


!#######################################################################
! <SUBROUTINE NAME="avg_ocean_sfc">
!
! <DESCRIPTION>
! Compute average of ocean surface quantities.  This is for coupling, 
! where pass time averaged information from ocean to other component
! models. Note that Ocean_sfc%frazil is NOT time averaged.  Rather, it 
! is accumulated from T_diag(index_frazil)%field in subroutine sum_ocean_sfc.
! Doing so is necessary for heat conservation between ocean and sea 
! ice systems.  Since it is not time averaged, frazil is not part of 
! this averaging subroutine.  
!
! Note that ocean model SST passed to the atmosphere is the surface
! potential temperature (which is equal to surface in situ temperature).
! If the ocean prognostic temperature variable is conservative temperature,
! then the sst is carried in T_diag(index_diag_temp).  If the prognostic 
! temperature is potential temperature, then the sst is carried in 
! T_prog(index_temp). 
!
! Note that if one removes the averaging, then we take only the 
! latest values of the surface fields.  This approach has been 
! found useful to stabilize the "concurrent" coupling approach.  
!
! Note that this routine is called after eta_and_pbot_diagnose,
! so Thickness%eta is eta_t(taup1).  
!
! </DESCRIPTION>
!
subroutine avg_ocean_sfc(Time, Thickness, T_prog, T_diag, Velocity, Ocean_sfc)

  type(ocean_time_type),         intent(in)    :: Time
  type(ocean_thickness_type),    intent(in)    :: Thickness
  type(ocean_prog_tracer_type),  intent(in)    :: T_prog(:)
  type(ocean_diag_tracer_type),  intent(in)    :: T_diag(:)
  type(ocean_velocity_type),     intent(in)    :: Velocity
  type(ocean_public_type), target, intent(inout) :: Ocean_sfc
  
  real, dimension(isd:ied,jsd:jed) :: sst
  real, dimension(isd:ied,jsd:jed) :: sst_redist
  real                             :: divid
  integer                          :: taup1, i, j, ii, jj
  integer :: stdoutunit
  stdoutunit=stdout()

  taup1 = Time%taup1
  sst        = 0.0
  sst_redist = 0.0

  if ( Ocean_sfc%avg_kount == 0) then 
    call mpp_error (FATAL,&
    '==>Error from ocean_sbc_mod (avg_ocean_sfc): no ocean surface quantities have been time averaged')
  endif 

  divid = 1./float(Ocean_sfc%avg_kount)

  do j = jsc_bnd,jec_bnd
     do i = isc_bnd,iec_bnd
        ii = i + i_shift
        jj = j + j_shift
        if(Grd%tmask(ii,jj,1) == 1.0) then
           Ocean_sfc%t_surf(i,j)  = Ocean_sfc%t_surf(i,j)*divid + kelvin  !C --> K
           Ocean_sfc%s_surf(i,j)  = Ocean_sfc%s_surf(i,j)*divid
           Ocean_sfc%u_surf(i,j)  = Ocean_sfc%u_surf(i,j)*divid
           Ocean_sfc%v_surf(i,j)  = Ocean_sfc%v_surf(i,j)*divid
           Ocean_sfc%sea_lev(i,j) = Ocean_sfc%sea_lev(i,j)*divid 
#if defined(ACCESS_CM) || defined(ACCESS_OM)
           Ocean_sfc%gradient(i,j,:) = Ocean_sfc%gradient(i,j,:)*divid
#endif
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
           Ocean_sfc%n_surf(i,j)  = Ocean_sfc%n_surf(i,j)*divid
           Ocean_sfc%alg_surf(i,j)  = Ocean_sfc%alg_surf(i,j)*divid
#endif
#if defined(ACCESS_CM)
           Ocean_sfc%co2flux(i,j) = Ocean_sfc%co2flux(i,j)*divid
           Ocean_sfc%co2(i,j)     = Ocean_sfc%co2(i,j)*divid
#endif
        endif
     enddo
  enddo


  !replace time-averaged t,s,sealev with latest value
  if(.NOT. avg_sfc_temp_salt_eta) then 

      if(prog_temp_variable==CONSERVATIVE_TEMP) then
          sst(:,:) = T_diag(index_diag_temp)%field(:,:,1)
          if(index_redist_heat > 0 ) sst_redist(:,:) = pottemp_from_contemp(T_prog(index_salt)%field(:,:,1,taup1),  &
                                                                            T_prog(index_redist_heat)%field(:,:,1,taup1))
      else
          sst(:,:) = T_prog(index_temp)%field(:,:,1,taup1)
          if(index_redist_heat > 0 ) sst_redist(:,:) = T_prog(index_redist_heat)%field(:,:,1,taup1)
      endif
      
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift
            if(Grd%tmask(ii,jj,1) == 1.0) then
                Ocean_sfc%s_surf(i,j) = T_prog(index_salt)%field(ii,jj,1,taup1)
                Ocean_sfc%sea_lev(i,j)= Thickness%sea_lev(ii,jj)
#if defined(ACCESS_CM)
                Ocean_sfc%co2flux(i,j) = co2flux(ii,jj)
                Ocean_sfc%co2(i,j)     = ocn_co2(ii,jj)
#endif
            endif
         enddo
      enddo

      ! FAFMIP sets SST to redistributed heat tracer value 
      if(index_redist_heat > 0) then
        write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod (avg_ocean_sfc): FAFMIP - Using redistributed heat tracer for sst.'
        do j = jsc_bnd,jec_bnd
           do i = isc_bnd,iec_bnd
              ii = i + i_shift
              jj = j + j_shift
              if(Grd%tmask(ii,jj,1) == 1.0) then
                  Ocean_sfc%t_surf(i,j) = sst_redist(ii,jj) + kelvin
              endif
           enddo
        enddo
      else ! if not doing FAFMIP use the standard sst
        do j = jsc_bnd,jec_bnd
           do i = isc_bnd,iec_bnd
              ii = i + i_shift
              jj = j + j_shift
              if(Grd%tmask(ii,jj,1) == 1.0) then
                  Ocean_sfc%t_surf(i,j) = sst(ii,jj) + kelvin
              endif
           enddo
        enddo
      endif
   
      
  endif


  !replace time-averaged u,v with latest value
  if(.NOT. avg_sfc_velocity) then 

      if(horz_grid == MOM_BGRID) then 

          do j = jsc_bnd,jec_bnd
             do i = isc_bnd,iec_bnd
                ii = i + i_shift
                jj = j + j_shift
                if(Grd%tmask(ii,jj,1) == 1.0) then
                    Ocean_sfc%u_surf(i,j) = Velocity%u(ii,jj,1,1,taup1)
                    Ocean_sfc%v_surf(i,j) = Velocity%u(ii,jj,1,2,taup1)
#if defined(ACCESS_CM) || defined(ACCESS_OM)
                    Ocean_sfc%gradient(i,j,:) = sslope(i,j,:)
#endif
                endif
             enddo
          enddo

      else  ! cgrid 

          ! coupler works with b-grid velocity, so we average c-grid to get b-grid
          do j = jsc_bnd,jec_bnd
             do i = isc_bnd,iec_bnd
                ii = i + i_shift
                jj = j + j_shift
                Ocean_sfc%u_surf(i,j)  = Grd%umask(ii,jj,1)*onehalf & 
                     *(Velocity%u(ii,jj,1,1,taup1) + Velocity%u(ii,jj+1,1,1,taup1))
                Ocean_sfc%v_surf(i,j)  = Grd%umask(ii,jj,1)*onehalf &
                     *(Velocity%u(ii,jj,1,2,taup1) + Velocity%u(ii,jj+1,1,2,taup1))
             enddo
          enddo

      endif

  endif 

  call ocean_tpm_avg_sfc(Dom, Ocean_sfc, Grd, isc_bnd, iec_bnd, jsc_bnd, jec_bnd)

  !set count to zero (surface quantities will be zeroed out before next sum)
  Ocean_sfc%avg_kount = 0   

  
end subroutine avg_ocean_sfc
! </SUBROUTINE> NAME="avg_ocean_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_sfc_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_sfc_restart(time_stamp)
  character(len=*),           intent(in), optional :: time_stamp

  call save_restart(Sfc_restart, time_stamp)
  if(waterflux_tavg) call save_restart(sbc_restart, time_stamp)

end subroutine ocean_sfc_restart
! </SUBROUTINE> NAME="ocean_sbc_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_sfc_end">
!
! <DESCRIPTION>
! Save information from Ocean_sfc to restarts. Note that it is 
! important in general to distinguish the time accumulated quantity 
! Ocean_sfc%frazil, saved here, from the instantaneous quantity 
! T_diag%frazil, which is saved in the diagnostic tracer restart file.  
! </DESCRIPTION>
!
subroutine ocean_sfc_end(Ocean_sfc)
  type(ocean_public_type), intent(in), target :: Ocean_sfc

    call ocean_sfc_restart

    call ocean_tpm_sfc_end

end subroutine ocean_sfc_end
! </SUBROUTINE> NAME="ocean_sfc_end"


!#######################################################################
! <SUBROUTINE NAME="get_ocean_sbc">
!
! <DESCRIPTION>
! Subroutine to get the surface fluxes passed into the ocean from 
! other component models. 
!
! **momentum fluxes from wind stress and momentum of pme and rivers
! **stokes drift from surface wave model 
! **water fluxes and temp/salinity in water fluxes 
! **salt fluxes; real or virtual
! **heat fluxes
! **applied surface pressure
!
! </DESCRIPTION>
!

subroutine get_ocean_sbc(Time, Ice_ocean_boundary, Thickness, Dens, Ext_mode, T_prog, Velocity, &
                         pme, melt, river, runoff, calving, upme, uriver, swflx, swflx_vis, patm)


  type(ocean_time_type),          intent(in)    :: Time 
  type(ice_ocean_boundary_type),  intent(in)    :: Ice_ocean_boundary
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  type(ocean_velocity_type),      intent(inout) :: Velocity
  real, dimension(isd:,jsd:),     intent(inout) :: pme
  real, dimension(isd:,jsd:),     intent(inout) :: melt
  real, dimension(isd:,jsd:),     intent(inout) :: river
  real, dimension(isd:,jsd:),     intent(inout) :: runoff
  real, dimension(isd:,jsd:),     intent(inout) :: calving 
  real, dimension(isd:,jsd:),     intent(inout) :: swflx
  real, dimension(isd:,jsd:),     intent(inout) :: swflx_vis
  real, dimension(isd:,jsd:),     intent(inout) :: patm
  real, dimension(isd:,jsd:,:),   intent(inout) :: upme
  real, dimension(isd:,jsd:,:),   intent(inout) :: uriver

  real, dimension(isd:ied,jsd:jed) :: tmp_patm

  real, dimension(isd:ied,jsd:jed) :: liquid_precip
  real, dimension(isd:ied,jsd:jed) :: frozen_precip
  real, dimension(isd:ied,jsd:jed) :: evaporation
  real, dimension(isd:ied,jsd:jed) :: sensible 
  real, dimension(isd:ied,jsd:jed) :: longwave
  real, dimension(isd:ied,jsd:jed) :: latent 

  type(time_type)                  :: time_dtime 

  real    :: tmp_x, tmp_y
  real    :: pme_river_total, var
  real    :: tracer_input
  real    :: tmp_runoff, tmp_calving, tmp_pme
  real    :: dayreal, phase 
  real    :: stokes_factor
  integer :: tau, taup1, n, i, j, k, ii, jj
  integer :: day, sec 
  real    :: potrhosfc_inv
  real    :: active_cells, smftu, smftv

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taup1 = Time%taup1

  ! save some of the flux fields for diagnostics 
  pme_river    = 0.0   ! (kg/m^2/sec) positive when enters liquid ocean 
  melt         = 0.0   ! (kg/m^2/sec) positive when enters liquid ocean 
  liquid_precip= 0.0   ! (kg/m^2/sec) positive when enters liquid ocean 
  frozen_precip= 0.0   ! (kg/m^2/sec) positive when enters liquid ocean 
  evaporation  = 0.0   ! (kg/m^2/sec) positive when enters liquid ocean
  sensible     = 0.0   ! (W/m^2)      positive when enters liquid ocean 
  longwave     = 0.0   ! (W/m^2)      positive when enters liquid ocean 
  latent       = 0.0   ! (W/m^2)      positive when enters liquid ocean  

  ! for diagnostics related to sea level forcing
  do j=jsd,jed
     do i=isd,ied
        rhosfc_inv(i,j) = Grd%tmask(i,j,1)/(epsln+Dens%rho(i,j,1,tau))
        alphasfc(i,j)   = -rhosfc_inv(i,j)*Dens%drhodT(i,j,1) 
        betasfc(i,j)    =  rhosfc_inv(i,j)*Dens%drhodS(i,j,1) 
        
        potrhosfc_inv    = Grd%tmask(i,j,1)/(epsln+Dens%potrho(i,j,1))
        alphasfc2(i,j)   = -potrhosfc_inv*Dens%dpotrhodT(i,j,1) 
        betasfc2(i,j)    =  potrhosfc_inv*Dens%dpotrhodS(i,j,1) 
     enddo
  enddo

  !------- Get 10m winds --------------------------------
  do j = jsc_bnd,jec_bnd
     do i = isc_bnd,iec_bnd
        ii = i + i_shift
        jj = j + j_shift
        Velocity%u10(ii,jj)= Ice_ocean_boundary%wnd(i,j)*Grd%tmask(ii,jj,1) !RASF 10m winds passin ACCESS, Wombat. c.f. MOM6 approach ustar to u10. T-GRID!
     enddo
  enddo
  call mpp_update_domains(Velocity%u10(:,:)  , Dom%domain2d)

  !------- Calculate Langmuir turbulence enhancement and stokes drift if do_langmuir is true ------------------------
    if ( do_langmuir ) then 
       do j = jsc_bnd,jec_bnd
          do i = isc_bnd,iec_bnd
             ii = i + i_shift
             jj = j + j_shift
! until full wave model is implemented Ice_ocean_boundary%ustoke etc not associated 
             if(associated(Ice_ocean_boundary%ustoke)) then
                Velocity%ustoke(ii,jj) = Ice_ocean_boundary%ustoke(i,j) 
                Velocity%vstoke(ii,jj) = Ice_ocean_boundary%vstoke(i,j)
                Velocity%wavlen(ii,jj)= Ice_ocean_boundary%wavlen(i,j)
             endif
          enddo
       enddo
    endif
  !--------momentum fluxes------------------------------------- 
  !

  ! i and j momentum flux (Newton/m^2) crossing ocean surface.
  ! assume the u_flux and v_flux from coupler are on B-grid. 
  if(.not. zero_surface_stress) then 
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift
            Velocity%smf_bgrid(ii,jj,1) = Ice_ocean_boundary%u_flux(i,j)*Grd%umask(ii,jj,1)
            Velocity%smf_bgrid(ii,jj,2) = Ice_ocean_boundary%v_flux(i,j)*Grd%umask(ii,jj,1)
         enddo
      enddo
  endif

  ! for case when winds are not aligned properly with the ocean grid 
  if (rotate_winds) then
      do j=jsc,jec
         do i=isc,iec
            tmp_x =  Grd%cos_rot(i,j)*Velocity%smf_bgrid(i,j,1) + Grd%sin_rot(i,j)*Velocity%smf_bgrid(i,j,2)
            tmp_y = -Grd%sin_rot(i,j)*Velocity%smf_bgrid(i,j,1) + Grd%cos_rot(i,j)*Velocity%smf_bgrid(i,j,2)
            Velocity%smf_bgrid(i,j,1) = tmp_x
            Velocity%smf_bgrid(i,j,2) = tmp_y
         enddo
      enddo
  endif
 
  ! for idealized tests 
  if (taux_sinx) then
      do j=jsc,jec
         do i=isc,iec
            Velocity%smf_bgrid(i,j,1) = 0.1*sin(Grd%xu(i,j)/(2.0*pi))*Grd%umask(i,j,1)
            Velocity%smf_bgrid(i,j,2) = 0.0
         enddo
      enddo
  endif
 
  ! for idealized tests 
  if (tauy_siny) then
      do j=jsc,jec
         do i=isc,iec
            Velocity%smf_bgrid(i,j,1) = 0.0
            Velocity%smf_bgrid(i,j,2) = 0.1*sin(Grd%yu(i,j)/(2.0*pi))*Grd%umask(i,j,1)
         enddo
      enddo
  endif
 
  ! mpp update domain is needed for vertical mixing schemes such
  ! as KPP and GOTM, as well as to get c-grid version of smf. 
  call mpp_update_domains(Velocity%smf_bgrid(:,:,1),Velocity%smf_bgrid(:,:,2),Dom%domain2d,gridtype=BGRID_NE)
  
  ! average the b-grid smf to get c-grid version.
  ! if running MOM_BGRID, then smf_cgrid is purely diagnostic.   
  ! if running MOM_CGRID, then smf_cgrid is used for momentum equation.
  ! note there is no mpp update domains needed for smf_cgrid,   
  ! since we use smf_bgrid for vertical mixing schemes such as KPP and GOTM
  ! and it is for these schemes we need to reach into halos with smf_bgrid.
  Velocity%smf_cgrid(:,:,:) = 0.0
  do j=jsc,jec
     do i=isc,iec
        if(Grd%umask(i,j,1) + Grd%umask(i,j-1,1) > 0.0) then 
            Velocity%smf_cgrid(i,j,1) =                                                                   &
            (Grd%umask(i,j,1)*Velocity%smf_bgrid(i,j,1) + Grd%umask(i,j-1,1)*Velocity%smf_bgrid(i,j-1,1)) &
            / (Grd%umask(i,j,1) + Grd%umask(i,j-1,1)) 
        endif
        if(Grd%umask(i,j,1) + Grd%umask(i-1,j,1) > 0.0) then 
            Velocity%smf_cgrid(i,j,2) =                                                                   &
            (Grd%umask(i,j,1)*Velocity%smf_bgrid(i,j,2) + Grd%umask(i-1,j,1)*Velocity%smf_bgrid(i-1,j,2)) &
            / (Grd%umask(i,j,1) + Grd%umask(i-1,j,1)) 
        endif
     enddo
  enddo

  ! for use in forcing momentum  
  if(horz_grid == MOM_BGRID) then 
     do n=1,2
        do j=jsc,jec
           do i=isc,iec
              Velocity%smf(i,j,n) = Velocity%smf_bgrid(i,j,n)
           enddo
        enddo
     enddo
  else
     do n=1,2
        do j=jsc,jec
           do i=isc,iec
              Velocity%smf(i,j,n) = Velocity%smf_cgrid(i,j,n)
           enddo
        enddo
     enddo
  endif 

  ! Calculate surface friction velocity
    do j=jsc,jec
        do i=isc,iec

!         ustar is needed on the "T-grid".  It is assumed that masking of 
!         smf over land was performed inside of the ocean_sbc module. 
!         smf has units of N/m^2 and so we need rho0r to get ustar in m/s.   
!         swflx has units W/m^2 so needs 1/rho_cp to get to C*m/s units.
!         these are proper units for buoyancy fluxes. 
          active_cells = Grd%umask(i,j,1)   + Grd%umask(i-1,j,1)   &
                        +Grd%umask(i,j-1,1) + Grd%umask(i-1,j-1,1) + epsln
          smftu = rho0r*(Velocity%smf_bgrid(i,j,1)   + Velocity%smf_bgrid(i-1,j,1)     &
                        +Velocity%smf_bgrid(i,j-1,1) + Velocity%smf_bgrid(i-1,j-1,1))  &
                  /active_cells
          smftv = rho0r*(Velocity%smf_bgrid(i,j,2)   + Velocity%smf_bgrid(i-1,j,2)    &
                        +Velocity%smf_bgrid(i,j-1,2) + Velocity%smf_bgrid(i-1,j-1,2)) &
                   /active_cells
          Velocity%ustar(i,j) = sqrt( sqrt(smftu**2 + smftv**2) )
     enddo
  enddo
  

  !--------stokes drift from surface wave model------------------------------- 
  ! smg: place holder until get code updates to Ice_ocean_boundary.
  ! Also to be determined: will stokes drift velocity be on A,B, or C grid? 
  ! Perhaps it should be same grid as smf. Or perhaps A grid is best...?
  ! i and j stokes drift (m/s) and stokes decay depth 
  if(read_stokes_drift) then 
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift
!            Velocity%stokes_drift(ii,jj,1,1) = Ice_ocean_boundary%ustokes(i,j)*Grd%umask(ii,jj,1)
!            Velocity%stokes_drift(ii,jj,2,1) = Ice_ocean_boundary%vstokes(i,j)*Grd%umask(ii,jj,1)
!            Velocity%stokes_depth(ii,jj)     = Ice_ocean_boundary%stokes_depth(i,j)*Grd%umask(ii,jj,1)
         enddo
      enddo
      do k=1,nk
          do j=jsc,jec
              do i=isc,iec
                  stokes_factor = -Grd%umask(i,j,k)*Thickness%depth_zu(i,j,k)/(epsln+Velocity%stokes_depth(i,j))
                  Velocity%stokes_drift(i,j,k,1) = Velocity%stokes_drift(i,j,1,1)*exp(stokes_factor)    
                  Velocity%stokes_drift(i,j,k,2) = Velocity%stokes_drift(i,j,1,2)*exp(stokes_factor)    
              enddo
          enddo
      enddo
  endif 


  !--------water and salt fluxes--------------------------------- 
  !
  ! Set temperature for water in evaporation and precipitation equal to the 
  ! ocean surface value. Default for other tracers is to have zero concentration 
  ! in evaporation and precipitation.
  do j=jsc,jec
     do i=isc,iec
       T_prog(index_temp)%tpme(i,j) = T_prog(index_temp)%field(i,j,1,tau)
     enddo
  enddo
  if(index_redist_heat > 0) then
    do j=jsc,jec
       do i=isc,iec
         T_prog(index_redist_heat)%tpme(i,j) = T_prog(index_redist_heat)%field(i,j,1,tau)
       enddo
    enddo
  endif
  if(index_added_heat > 0) then
    do j=jsc,jec
       do i=isc,iec
         T_prog(index_added_heat)%tpme(i,j) = T_prog(index_added_heat)%field(i,j,1,tau)
       enddo
    enddo
  endif

  ! set velocity of pme and river water to that of upper ocean cell.   
  ! generalizations may be suitable with refined component models.   
  do n = 1, size(upme,3)
     do j = jsc, jec
        do i = isc, iec
           upme(i,j,n)   = Velocity%u(i,j,1,n,tau)   ! velocity of precip-evap water
           uriver(i,j,n) = Velocity%u(i,j,1,n,tau)   ! velocity of river water 
        enddo
     enddo
  enddo

  ! start of long if-block for use_waterflux true or false. 
  if (use_waterflux) then

      ! water flux in (kg/m^3)*(m/s) from liquid, frozen, and q_flux (evap/condense).  
      ! note sign switch on evaporation(i,j), so that evaporation(i,j) > 0 when 
      ! water enters ocean (as when condensation occurs from fog).  for most regions,
      ! evaporation(i,j) < 0.  
      do j = jsc_bnd,jec_bnd
          do i = isc_bnd,iec_bnd
              ii = i + i_shift
              jj = j + j_shift
              pme(ii,jj) = (Ice_ocean_boundary%lprec(i,j) + Ice_ocean_boundary%fprec(i,j) &
#if defined(ACCESS_CM) || defined(ACCESS_OM)
                          ! PME is meant to include "melt", in a MOM+SIS configuration it
                          ! is added by the coupler. We add it here. 
                          + Ice_ocean_boundary%wfimelt(i,j) & 
                          + Ice_ocean_boundary%wfiform(i,j) &
                          + Ice_ocean_boundary%licefw(i,j) &
#endif
                          - Ice_ocean_boundary%q_flux(i,j))*Grd%tmask(ii,jj,1) 
              liquid_precip(ii,jj) =  Ice_ocean_boundary%lprec(i,j)*Grd%tmask(ii,jj,1) 
              frozen_precip(ii,jj) =  Ice_ocean_boundary%fprec(i,j)*Grd%tmask(ii,jj,1) 
              evaporation(ii,jj)   = -Ice_ocean_boundary%q_flux(i,j)*Grd%tmask(ii,jj,1) 
          enddo
      enddo
      if(waterflux_tavg) then 
        do j = jsc_bnd,jec_bnd
          do i = isc_bnd,iec_bnd
              ii = i + i_shift
              jj = j + j_shift
              ! Take half the waterflux from the previous timestep and add to half 
              ! the flux from the current timestep calculated above. Set pme_taum1 
              ! to calculated pme for this timestep
              tmp_pme = pme(ii,jj)
              pme(ii,jj) = 0.5*pme(ii,jj) + 0.5*pme_taum1(ii,jj)
              pme_taum1(ii,jj) = tmp_pme
          end do
        end do
      endif 

      ! water flux in (kg/m^3)*(m/s) from liquid runoff and calving land ice  
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift
            runoff(ii,jj)  = Ice_ocean_boundary%runoff(i,j)*Grd%tmask(ii,jj,1)
            calving(ii,jj) = Ice_ocean_boundary%calving(i,j)*Grd%tmask(ii,jj,1)
         enddo
      enddo
      if(use_ideal_runoff) then
          do j=jsc,jec
             do i=isc,iec
                runoff(i,j)  = runoff(i,j)  + ideal_runoff(i,j)
             enddo
          enddo
      endif
      if(use_ideal_calving) then
          do j=jsc,jec
             do i=isc,iec
                calving(i,j) = calving(i,j) + ideal_calving(i,j)
             enddo
          enddo
      endif
      if(runoffspread)  call spread_river_horz(runoff)
      if(calvingspread) call spread_river_horz(calving)

      ! river=runoff+calving is the water mass flux 
      ! entering ocean, other than through pme. 
      ! calving is immediately melted, so it is appropriate
      ! to add it to runoff for purposes of getting a mass
      ! flux of liquid, from land, into the ocean.   
      if(waterflux_tavg) then 
          do j=jsc,jec
             do i=isc,iec
                river(i,j)       = 0.5*river_taum1(i,j)
                river_taum1(i,j) = runoff(i,j) + calving(i,j)
                river(i,j)       = river(i,j) + 0.5*river_taum1(i,j)
             enddo
          enddo
      else 
          do j=jsc,jec
             do i=isc,iec
                river(i,j) = runoff(i,j) + calving(i,j)
             enddo
          enddo
      endif


      ! Set the temperature flux associated with the water 
      ! entering ocean from land. This flux equals to the mass 
      ! flux of water times the temperature of water, whether 
      ! this be liquid or solid water.  
      if(land_model_heat_fluxes) then  

          ! for cases where the land model computes the heat flux (W/m^2)
          ! (relative to 0degC) associated with the liquid runoff and 
          ! solid calving ice.  In this case, we have heat_land/cp_ocean
          ! as the temperature flux going into evolving the ocean temp. 
          ! Notably, we should NOT be dividing by cp_liquid_runoff nor 
          ! cp_solid_runoff to get this temperature flux.  
          do j = jsc_bnd,jec_bnd
             do i = isc_bnd,iec_bnd
                ii = i + i_shift
                jj = j + j_shift
                T_prog(index_temp)%runoff_tracer_flux(ii,jj) = cp_ocean_r &
                     *Ice_ocean_boundary%runoff_hflx(i,j)*Grd%tmask(ii,jj,1)
                T_prog(index_temp)%calving_tracer_flux(ii,jj) = cp_ocean_r &
                     *Ice_ocean_boundary%calving_hflx(i,j)*Grd%tmask(ii,jj,1)
             enddo
          enddo
          if(index_redist_heat > 0) then
            do j=jsc,jec
              do i=isc,iec
                T_prog(index_redist_heat)%runoff_tracer_flux(i,j)  = T_prog(index_temp)%runoff_tracer_flux(i,j)
                T_prog(index_redist_heat)%calving_tracer_flux(i,j) = T_prog(index_temp)%calving_tracer_flux(i,j)
             enddo
           enddo
          endif
          
          ! For diagnostic purposes, compute the effective temperatures 
          ! trunoff and tcalving.  These are NOT the temperatures that the
          ! land model may have, as those temperatures are obtained by using
          ! cp_liquid_runoff and cp_solid_runoff. Instead, this is an effective
          ! temperature resulting from dividing the land heat flux by the ocean
          ! heat capacity. We do so for purposes of diagnostics transparency 
          ! in ocean_tracer_diag...
          !
          ! also note the code to avoid non-representable numbers in 
          ! trunoff/tcalving diagnostics when runoff/calving are small.
          do j=jsc,jec
             do i=isc,iec
                T_prog(index_temp)%trunoff(i,j)  = 0.0
                T_prog(index_temp)%tcalving(i,j) = 0.0

                if(runoff(i,j) > epsln) then
                    T_prog(index_temp)%trunoff(i,j) = &
                    min(100.0, max(-273.15,T_prog(index_temp)%runoff_tracer_flux(i,j)/runoff(i,j)))
                endif
                if(calving(i,j) > epsln) then
                    T_prog(index_temp)%tcalving(i,j)= &
                    min(100.0, max(-273.15,T_prog(index_temp)%calving_tracer_flux(i,j)/calving(i,j)))
                endif
             enddo
          enddo

          ! for FAFMIP redistributed heat tracer 
          if(index_redist_heat > 0) then
            do j=jsc,jec
               do i=isc,iec
                  T_prog(index_redist_heat)%trunoff(i,j)  = 0.0
                  T_prog(index_redist_heat)%tcalving(i,j) = 0.0

                  if(runoff(i,j) > epsln) then
                      T_prog(index_redist_heat)%trunoff(i,j) = &
                      min(100.0, max(-273.15,T_prog(index_redist_heat)%runoff_tracer_flux(i,j)/runoff(i,j)))
                  endif
                  if(calving(i,j) > epsln) then
                      T_prog(index_redist_heat)%tcalving(i,j)= &
                      min(100.0, max(-273.15,T_prog(index_redist_heat)%calving_tracer_flux(i,j)/calving(i,j)))
                  endif
               enddo
            enddo
          endif           

          ! for FAFMIP added heat tracer 
          if(index_added_heat > 0) then
            do j=jsc,jec
               do i=isc,iec
                  T_prog(index_added_heat)%trunoff(i,j)  = 0.0
                  T_prog(index_added_heat)%tcalving(i,j) = 0.0

                  if(runoff(i,j) > epsln) then
                      T_prog(index_added_heat)%trunoff(i,j) = &
                      min(100.0, max(-273.15,T_prog(index_added_heat)%runoff_tracer_flux(i,j)/runoff(i,j)))
                  endif
                  if(calving(i,j) > epsln) then
                      T_prog(index_added_heat)%tcalving(i,j)= &
                      min(100.0, max(-273.15,T_prog(index_added_heat)%calving_tracer_flux(i,j)/calving(i,j)))
                  endif
               enddo
            enddo
          endif           
          
      else if (river_temp_ofam) then
          ! OFAM sometimes specifies the temperatures of the rivers according
          ! to climatological SST.
          call time_interp_external(id_restore(index_temp),Time%model_time, &
                                    data)
          do j = jsc, jec
              do i = isc, iec
                  T_prog(index_temp)%trunoff(i, j) &
                    = max(runoff_temp_min, data(i, j))
                  T_prog(index_temp)%tcalving(i, j) &
                    = T_prog(index_temp)%field(i, j, 1, tau)
                  T_prog(index_temp)%runoff_tracer_flux(i, j) &
                    = Grd%tmask(i, j, 1) * T_prog(index_temp)%trunoff(i,j) &
                        * runoff(i, j)
                  T_prog(index_temp)%calving_tracer_flux(i,j) &
                    = Grd%tmask(i, j, 1) * T_prog(index_temp)%tcalving(i,j) &
                        * calving(i, j)
              end do
          end do

      else 

          ! for cases where the land model does not carry heat flux associated with 
          ! the liquid runoff and solid calving ice. We assign a temperature to the 
          ! liquid and solid runoff.  
          ! Set temperature for liquid runoff water to ocean surface value, but no 
          ! less than runoff_temp_min. For other tracers, by default keep them   
          ! set to their initial concentration.
          ! For calving temperature, set this equal to the SST.  
          do j=jsc,jec
             do i=isc,iec
                T_prog(index_temp)%trunoff(i,j)  = &
                   max(runoff_temp_min,T_prog(index_temp)%field(i,j,1,tau))
                T_prog(index_temp)%tcalving(i,j) = &
                   T_prog(index_temp)%field(i,j,1,tau)
                T_prog(index_temp)%runoff_tracer_flux(i,j) = &
                   Grd%tmask(i,j,1)*T_prog(index_temp)%trunoff(i,j)*runoff(i,j)
                T_prog(index_temp)%calving_tracer_flux(i,j)= &
                   Grd%tmask(i,j,1)*T_prog(index_temp)%tcalving(i,j)*calving(i,j)
             enddo
          enddo

          ! for FAFMIP redistributed heat tracer 
          if(index_redist_heat > 0) then
            do j=jsc,jec
               do i=isc,iec
                 T_prog(index_redist_heat)%trunoff(i,j)  = &
                   max(runoff_temp_min,T_prog(index_redist_heat)%field(i,j,1,tau))
                 T_prog(index_redist_heat)%tcalving(i,j) = &
                   T_prog(index_redist_heat)%field(i,j,1,tau)
                 T_prog(index_redist_heat)%runoff_tracer_flux(i,j) = &
                   Grd%tmask(i,j,1)*T_prog(index_redist_heat)%trunoff(i,j)*runoff(i,j)
                 T_prog(index_redist_heat)%calving_tracer_flux(i,j)= &
                   Grd%tmask(i,j,1)*T_prog(index_redist_heat)%tcalving(i,j)*calving(i,j)
               enddo
            enddo
          endif           

          ! for FAFMIP added heat tracer 
          if(index_added_heat > 0) then
            do j=jsc,jec
               do i=isc,iec
                 T_prog(index_added_heat)%trunoff(i,j)  = &
                   max(runoff_temp_min,T_prog(index_added_heat)%field(i,j,1,tau))
                 T_prog(index_added_heat)%tcalving(i,j) = &
                   T_prog(index_added_heat)%field(i,j,1,tau)
                 T_prog(index_added_heat)%runoff_tracer_flux(i,j) = &
                   Grd%tmask(i,j,1)*T_prog(index_added_heat)%trunoff(i,j)*runoff(i,j)
                 T_prog(index_added_heat)%calving_tracer_flux(i,j)= &
                   Grd%tmask(i,j,1)*T_prog(index_added_heat)%tcalving(i,j)*calving(i,j)
               enddo
            enddo
          endif           

      endif     ! land_model_heat_fluxes

      ! compute temperature and salinity of "river" according to 
      ! temperature and salinity of calving and runoff. 
      ! allow for possibility of runoff and/or calving to be negative,
      ! which may exist for certain land models that "suck" water from 
      ! the oceans to address limitations in their formulation. 
      ! We do not include such negative points when computing triver. 
      ! Also compute the salinity flux associated with liquid and solid runoff. 
      do j=jsc,jec
         do i=isc,iec
            tmp_runoff = max(0.0,runoff(i,j))
            tmp_calving= max(0.0,calving(i,j))
            T_prog(index_temp)%triver(i,j) = Grd%tmask(i,j,1)    &
                 *(tmp_runoff *T_prog(index_temp)%trunoff(i,j)   &
                  +tmp_calving*T_prog(index_temp)%tcalving(i,j)) &
                 /(epsln + tmp_runoff + tmp_calving)

            T_prog(index_salt)%triver(i,j) = Grd%tmask(i,j,1)    &
                 *(tmp_runoff *T_prog(index_salt)%trunoff(i,j)   &
                  +tmp_calving*T_prog(index_salt)%tcalving(i,j)) &
                 /(epsln + tmp_runoff + tmp_calving)
         enddo
      enddo

      if(index_redist_heat > 0) then 
        do j=jsc,jec
           do i=isc,iec
              tmp_runoff = max(0.0,runoff(i,j))
              tmp_calving= max(0.0,calving(i,j))
              T_prog(index_redist_heat)%triver(i,j) = Grd%tmask(i,j,1)    &
                   *(tmp_runoff *T_prog(index_redist_heat)%trunoff(i,j)   &
                    +tmp_calving*T_prog(index_redist_heat)%tcalving(i,j)) &
                  /(epsln + tmp_runoff + tmp_calving)
           enddo
        enddo
      endif       

      if(index_added_heat > 0) then 
        do j=jsc,jec
           do i=isc,iec
              tmp_runoff = max(0.0,runoff(i,j))
              tmp_calving= max(0.0,calving(i,j))
              T_prog(index_added_heat)%triver(i,j) = Grd%tmask(i,j,1)    &
                   *(tmp_runoff *T_prog(index_added_heat)%trunoff(i,j)   &
                    +tmp_calving*T_prog(index_added_heat)%tcalving(i,j)) &
                  /(epsln + tmp_runoff + tmp_calving)
           enddo
        enddo
      endif       
      
      ! for debugging the water fluxes 
      if(debug_water_fluxes) then 

          write(stdoutunit,'(a)')

          wrk1_2d(:,:) = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*cp_ocean &
                               *T_prog(index_temp)%runoff_tracer_flux(i,j)*dtime
             enddo
          enddo
          tracer_input = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), global_sum_flag)

          ! truncate out the lower order bits if using non_bitwise_reproducible sums
          if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) tracer_input = real(tracer_input, kind=FLOAT_KIND)
          write (stdoutunit,'(a,es24.17,a)') ' Heat input via liquid runoff in ocean_sbc          = ',&
                                             tracer_input,' Joule'
          wrk1_2d(:,:) = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*cp_ocean &
                               *T_prog(index_temp)%calving_tracer_flux(i,j)*dtime
             enddo
          enddo
          tracer_input = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), global_sum_flag)

          ! truncate out the lower order bits if using non_bitwise_reproducible sums
          if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) tracer_input = real(tracer_input, kind=FLOAT_KIND)
          write (stdoutunit,'(a,es24.17,a)') ' Heat input via calving land ice in ocean_sbc       = ',&
                                           tracer_input,' Joule'

          write(stdoutunit,'(a)')

          if(zero_water_fluxes) then 
              river     = 0.0
              pme       = 0.0
              pme_taum1 = 0.0
              runoff    = 0.0
              calving   = 0.0
          endif
          if(zero_river_fluxes) then 
              river     = 0.0
          endif
          if(zero_pme_fluxes) then 
              pme       = 0.0
              pme_taum1 = 0.0
          endif
          if(zero_runoff_fluxes) then 
              runoff    = 0.0
          endif
          if(zero_calving_fluxes) then 
              calving   = 0.0
          endif
          if(convert_river_to_pme) then 
              do j=jsc,jec
                 do i=isc,iec 
                    pme(i,j) = pme(i,j) + river(i,j)
                 enddo
              enddo
              river   = 0.0
              calving = 0.0
              runoff  = 0.0
          endif

      endif  ! for debug_water_fluxes


      ! set riverdiffuse to determine where to enhance diff_cbt 
      ! inside ocean_rivermix_mod.  this option is sometimes used
      ! to help mix tracers vertically near river mouths.  
      do n=1,num_prog_tracers
         do j=jsc,jec
            do i=isc,iec
               T_prog(n)%riverdiffuse(i,j) = river(i,j) 
            enddo
         enddo
      enddo

      ! Allow for nonzero salt flux from, say, ice melt.  
      ! flux from ice-melt is in units of (kg/m^3)*(m/sec).
      ! to convert to rho*dz*salinity units for the model, 
      ! need to multiply by 1000.0.
      ! 
      ! The melt field is deduced from knowing that 
      ! salt transferred between the ice and ocean
      ! occurs in the presence of liquid water transport
      ! between the ice and ocean. When salt is added to 
      ! the ocean from the ice melt, fresh liquid water is 
      ! also added.
      ! 
      ! The minus sign arises since the ice model produces a 
      ! salt_flux > 0 when there is salt added to the ice model,
      ! and hence taken away from the ocean model.  
      melt = 0.0
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift      
            T_prog(index_salt)%stf(ii,jj) = -Ice_ocean_boundary%salt_flux(i,j)*1000.0
            melt(ii,jj)                   = -Ice_ocean_boundary%salt_flux(i,j)*ice_salt_concentration_r
         enddo
      enddo

      ! produce a zero area average of pme + river. 
      ! note that we remove the ice melt water  
      ! since the coupler has already included the 
      ! ice melt within the pme field.
      ! Also remove ideal_runoff and ideal_calving, since
      ! wish to have these fields contribute to a net change
      ! in the water content of the ocean.   
      pme_river    = 0.0
      wrk1_2d(:,:) = 0.0
      if(use_ideal_runoff .or. use_ideal_calving) then 
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = ideal_runoff(i,j) + ideal_calving(i,j)
             enddo
          enddo
      endif
      if(zero_net_water_coupler) then 
         do j=jsc,jec
            do i=isc,iec
               pme_river(i,j) = pme(i,j) + river(i,j) - melt(i,j) - wrk1_2d(i,j)
            enddo
         enddo
         pme_river_total = mpp_global_sum(Dom%domain2d,pme_river(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1),&
                                          global_sum_flag)/Grd%tcellsurf

         ! truncate out the lower order bits if using non_bitwise_reproducible sums
         if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) pme_river_total = real(pme_river_total, kind=FLOAT_KIND)

         ! for diagnostic purposes, assume that the global adjustment occurs solely
         ! within liquid_precip(i,j), so that we maintain the diagnostic equality 
         ! pme(i,j) = liquid_precip(i,j)+frozen_precip(i,j)+evaporation(i,j). 
         ! (Remember that each of the terms liquid_precip, frozen_precip, and evaporation
         ! are > 0 when water enters ocean and < 0 when water leaves).
         do j=jsc,jec
            do i=isc,iec
               pme(i,j)           = pme(i,j)           - pme_river_total*Grd%tmask(i,j,1)
               liquid_precip(i,j) = liquid_precip(i,j) - pme_river_total*Grd%tmask(i,j,1)
            enddo
         enddo
      endif 
      
      ! when have a nonzero salinity of runoff, then there is a 
      ! nonzero salt flux (kg/(m^2*sec)) into the ocean with the runoff. 
      do j=jsc,jec
         do i=isc,iec
            T_prog(index_salt)%runoff_tracer_flux(i,j) = &
            T_prog(index_salt)%trunoff(i,j)*runoff(i,j)*T_prog(index_salt)%conversion 
         enddo
      enddo

  else     
  ! now enter code block for use_waterflux=.false. 
  ! this option is becoming obsolete and is maintained
  ! at a rudimenatary level solely for legacy purposes. 
  ! It is generally not recommended.  

     
      ! water flux in (kg/m^3)*(m/s) from rivers and calving land glaciers.
      ! no ideal_runoff nor ideal_calving when use_waterflux=.false. 
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift
            runoff(ii,jj)  = Ice_ocean_boundary%runoff(i,j)*Grd%tmask(ii,jj,1)
            calving(ii,jj) = Ice_ocean_boundary%calving(i,j)*Grd%tmask(ii,jj,1)
         enddo
      enddo
      if(runoffspread)  call spread_river_horz(runoff)
      if(calvingspread) call spread_river_horz(calving)     
      do j=jsc,jec
         do i=isc,iec
            river(i,j) = runoff(i,j) + calving(i,j)
         enddo
      enddo

      ! convert salt flux from ice to psu with multiplication by 1000.0. 
      ! convert freshwater mass fluxes (kg/m^3)*(m/sec) into virtual salt
      ! fluxes with multiplication by salinity_ref. 
         do j = jsc_bnd, jec_bnd
            do i = isc_bnd, iec_bnd
               ii = i + i_shift
               jj = j + j_shift
               T_prog(index_salt)%stf(ii,jj) = -Ice_ocean_boundary%salt_flux(i,j)*1000.0     -&
                                               (Ice_ocean_boundary%lprec(i,j)                +&
#if defined(ACCESS_CM) || defined(ACCESS_OM)
                                               (Ice_ocean_boundary%wfimelt(i,j)              +&
                                               Ice_ocean_boundary%wfiform(i,j))              +&
#endif
                                                Ice_ocean_boundary%fprec(i,j) + river(ii,jj) -&
                                                Ice_ocean_boundary%q_flux(i,j))*salinity_ref*Grd%tmask(ii,jj,1) 
         enddo
      enddo
      ! set the riverdiffuse "mask" to determine where to enhance diff_cbt
      do n=1,num_prog_tracers
         do j = jsc,jec
            do i = isc, iec
               T_prog(n)%riverdiffuse(i,j) = river(i,j) 
            enddo
         enddo
      enddo

      ! set pme and river to zero since use_waterflux=.false. 
      pme     = 0.0
      river   = 0.0
      runoff  = 0.0
      calving = 0.0

  endif    ! end of long if-block for if(use_waterflux) 

  call mpp_update_domains(pme(:,:)  , Dom%domain2d)
  call mpp_update_domains(river(:,:), Dom%domain2d)


  !---------------surface pressure from ice and atmos-----------
  !
  ! apply atmosphere and ice pressure to ocean only when there is
  ! mass transfer allowed between ocean and atmos/ice via fresh water. 
  ! if do not use fresh water, then cannot allow ocean to feel 
  ! weight of atmosphere and ice.  
  !
  ! NOTE: the option use_full_patm_for_sea_level allows for passing 
  ! the sea level, including the full weight of sea ice, back to
  ! the ice model.  This approach maintains the max weight on the liquid
  ! ocean according to the nml variable max_ice_thickness.  But it does 
  ! allow the sea ice to know when there is actually more sea ice than that
  ! set by max_ice_thickness.  This option then provides for a negative
  ! feedback on the runaway growth of sea ice, since the full pressure acting to 
  ! make the ice flow will be correctly felt.  This is a new option, and is not
  ! fully tested, so the default is use_full_patm_for_sea_level=.false.
  ! 
  tmp_patm(:,:) = 0.0
  if(use_waterflux) then 
      var = grav*rho0*max_ice_thickness
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift
            patm(ii,jj)     = min(Ice_ocean_boundary%p(i,j),var)
            tmp_patm(ii,jj) = Ice_ocean_boundary%p(i,j)
         enddo
      enddo
      call mpp_update_domains(patm(:,:), Dom%domain2d)
      call mpp_update_domains(tmp_patm(:,:), Dom%domain2d)
      if(use_full_patm_for_sea_level) then 
          do j=jsd,jed
             do i=isd,ied
                Ext_mode%patm_t(i,j,taup1)     = patm(i,j) 
                Ext_mode%patm_for_sea_lev(i,j) = tmp_patm(i,j)
             enddo
          enddo
      else
          do j=jsd,jed
             do i=isd,ied
                Ext_mode%patm_t(i,j,taup1)     = patm(i,j) 
                Ext_mode%patm_for_sea_lev(i,j) = patm(i,j) 
             enddo
          enddo
      endif
  endif


  !--------compute surface heat fluxes-------------------------- 
  !
  ! set latent heat of fusion and vaporization 
  call compute_latent_heat_fusion(Dens%rho_salinity(:,:,1,tau)) 
  call compute_latent_heat_vapor(Dens%rho_salinity(:,:,1,tau), T_prog(index_temp)%field(:,:,1,tau))
  
  if(.not. zero_heat_fluxes) then
     
    do j = jsc_bnd, jec_bnd
       do i = isc_bnd, iec_bnd
          ii = i + i_shift
          jj = j + j_shift  
          T_prog(index_temp)%stf(ii,jj) = (Ice_ocean_boundary%sw_flux_vis_dir(i,j)               + &
                                           Ice_ocean_boundary%sw_flux_vis_dif(i,j)               + &
                                           Ice_ocean_boundary%sw_flux_nir_dir(i,j)               + &
                                           Ice_ocean_boundary%sw_flux_nir_dif(i,j)               + &
                                           Ice_ocean_boundary%lw_flux(i,j)                       - &
                                           (Ice_ocean_boundary%fprec(i,j)                        + &
                                            calving(ii,jj))*latent_heat_fusion(ii,jj)            - &
                                           Ice_ocean_boundary%t_flux(i,j)                        - &
                                           Ice_ocean_boundary%q_flux(i,j)*latent_heat_vapor(ii,jj) &
#if defined(ACCESS_CM) || defined(ACCESS_OM)
                                           + Ice_ocean_boundary%mh_flux(i,j)  &
                                           + Ice_ocean_boundary%liceht(i,j) &
#endif
                                           )/cp_ocean*Grd%tmask(ii,jj,1)
          
       enddo
    enddo

    ! to ensure bitwise reproducibility with result prior to saving these 
    ! diagnostic fields, keep the following loop outside of the main loop above. 
    do j = jsc_bnd, jec_bnd
       do i = isc_bnd, iec_bnd
          ii = i + i_shift
          jj = j + j_shift  
          sensible(ii,jj) = -Ice_ocean_boundary%t_flux(i,j)*Grd%tmask(ii,jj,1)
          longwave(ii,jj) =  Ice_ocean_boundary%lw_flux(i,j)*Grd%tmask(ii,jj,1)
          latent(ii,jj)   =  latent_heat_vapor(ii,jj)*evaporation(ii,jj) &
#if defined(ACCESS_CM) || defined(ACCESS_OM)
          +Ice_ocean_boundary%liceht(ii,jj) &   
#endif
                            -latent_heat_fusion(ii,jj)*(frozen_precip(ii,jj)+calving(ii,jj))
       enddo
    enddo

  endif  ! .not. zero_heat_fluxes


  ! over-ride the boundary heat fluxes with a constant value. 
  ! this code is only of use for idealized tests.  
  if(sbc_heat_fluxes_const) then 
      if(sbc_heat_fluxes_const_seasonal) then
          time_dtime = increment_time(Time%model_time, int(dtime),0)  
          call get_time(time_dtime, sec, day)                                 
          dayreal = day + sec/86400.0 
          phase   = sin(twopi*dayreal*days_in_year_r) 
          do j=jsc,jec
             do i=isc,iec
                T_prog(index_temp)%stf(i,j) = phase*sbc_heat_fluxes_const_value*Grd%tmask(i,j,1)/cp_ocean
             enddo
          enddo
      else 
          do j=jsc,jec
             do i=isc,iec
                T_prog(index_temp)%stf(i,j) = sbc_heat_fluxes_const_value*Grd%tmask(i,j,1)/cp_ocean
             enddo
          enddo
      endif
  endif 


  ! shortwave flux (W/m^2) into ocean
  if(.not. zero_heat_fluxes) then 
      swflx(isc:iec,jsc:jec) = Grd%tmask(isc:iec,jsc:jec,1)                  &
       *(Ice_ocean_boundary%sw_flux_nir_dir(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) &
       + Ice_ocean_boundary%sw_flux_nir_dif(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) &
       + Ice_ocean_boundary%sw_flux_vis_dir(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) &
       + Ice_ocean_boundary%sw_flux_vis_dif(isc_bnd:iec_bnd,jsc_bnd:jec_bnd))

      swflx_vis(isc:iec,jsc:jec) = Grd%tmask(isc:iec,jsc:jec,1)               &
        *(Ice_ocean_boundary%sw_flux_vis_dir(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) &
       +  Ice_ocean_boundary%sw_flux_vis_dif(isc_bnd:iec_bnd,jsc_bnd:jec_bnd))
  endif 


  ! Boundary condition for FAFMIP redistributed tracer is same as for regular temperature field
  ! plus the Qfaf added heat which is added inside the "flux_adjust" routine (called after get_ocean_sbc). 
  if(index_redist_heat > 0) then 
    do j=jsc,jec
      do i=isc,iec
        T_prog(index_redist_heat)%stf(i,j) = T_prog(index_temp)%stf(i,j)
      enddo
    enddo
 endif
 
  ! FAFMIP added heat has zero surface flux at this point.
  ! It will have "correction" flux added later in "flux_adjust". 
  if(index_added_heat > 0) then 
    do j=jsc,jec
      do i=isc,iec
        T_prog(index_added_heat)%stf(i,j) = 0.0
      enddo
    enddo
  endif
 
  !------------------------------------------------------------------
  ! waterflux override code to override mass flux from 
  ! terms contributing to latent heat fluxes.  
  ! mass flux is over-ridden whereas latent (computed earlier) is not.  
  ! assume that use_waterflux=.true. for the over-ride to make sense. 

  ! water mass flux from calving.
  if(use_waterflux_override_calving .and. id_calving_override > 0) then 

     ! remove from river the calving input through Ice_ocean_boundary%calving
     do j=jsc,jec
         do i=isc,iec
            river(i,j) = Grd%tmask(i,j,1)*(river(i,j) - calving(i,j))           
         enddo
     enddo

     ! fill calving with value from input file 
     ! (calving > 0 means water enters ocean)
     calving(:,:) = 0.0
     wrk1_2d(:,:) = 0.0
     call time_interp_external(id_calving_override, Time%model_time, wrk1_2d)
     do j=jsc,jec
         do i=isc,iec
            calving(i,j) = wrk1_2d(i,j)*Grd%tmask(i,j,1)
            river(i,j)   = river(i,j) + calving(i,j)
            latent(i,j)  = latent_heat_vapor(i,j)*evaporation(i,j) &
                          -latent_heat_fusion(i,j)*(frozen_precip(i,j)+calving(i,j))
         enddo
     enddo
     call mpp_update_domains(river(:,:), Dom%domain2d)

  endif 

  ! water mass flux from fprec = frozen_precip.
  if(use_waterflux_override_fprec .and. id_fprec_override > 0) then 

     ! remove from pme the mass input from Ice_ocean_boundary%fprec 
     do j=jsc,jec
         do i=isc,iec
            pme(i,j) = Grd%tmask(i,j,1)*(pme(i,j) - frozen_precip(i,j))           
         enddo
     enddo

     ! fill frozen_precip with value from input file 
     ! (frozen_precip > 0 means water enters ocean)
     frozen_precip(:,:) = 0.0
     wrk1_2d(:,:)       = 0.0
     call time_interp_external(id_fprec_override, Time%model_time, wrk1_2d)
     do j=jsc,jec
         do i=isc,iec
            frozen_precip(i,j) = wrk1_2d(i,j)*Grd%tmask(i,j,1)
            pme(i,j)           = pme(i,j) + frozen_precip(i,j)
            latent(i,j)        = latent_heat_vapor(i,j)*evaporation(i,j) &
                                -latent_heat_fusion(i,j)*(frozen_precip(i,j)+calving(i,j))
         enddo
     enddo

  endif 


  ! water mass flux from evaporation = -q_flux.
  if(use_waterflux_override_evap .and. id_evap_override > 0) then 

     ! remove from pme the mass input from Ice_ocean_boundary
     do j=jsc,jec
         do i=isc,iec
            pme(i,j) = Grd%tmask(i,j,1)*(pme(i,j) - evaporation(i,j))           
         enddo
     enddo

     ! fill evaporation with value from input file
     ! (evaporation > 0 means liquid water enters ocean)
     evaporation(:,:) = 0.0
     wrk1_2d(:,:) = 0.0
     call time_interp_external(id_evap_override, Time%model_time, wrk1_2d)
     do j=jsc,jec
         do i=isc,iec
            evaporation(i,j) = wrk1_2d(i,j)*Grd%tmask(i,j,1)
            pme(i,j)         = pme(i,j) + evaporation(i,j)
            latent(i,j)      = latent_heat_vapor(i,j)*evaporation(i,j) &
                              -latent_heat_fusion(i,j)*(frozen_precip(i,j)+calving(i,j))
         enddo
     enddo

  endif 

  if(use_waterflux_override_fprec .or. use_waterflux_override_evap) then 
     call mpp_update_domains(pme(:,:), Dom%domain2d)
  endif 

#if defined(ACCESS_CM)
  do j = jsc_bnd,jec_bnd
     do i = isc_bnd,iec_bnd
        ii = i + i_shift
        jj = j + j_shift
        atm_co2(ii,jj) = Ice_ocean_boundary%co2(i,j)*Grd%tmask(ii,jj,1)
     enddo
  enddo
#endif

#if defined(ACCESS_CM) || defined(ACCESS_OM)
  do j = jsc_bnd,jec_bnd
     do i = isc_bnd,iec_bnd
        ii = i + i_shift
        jj = j + j_shift
        aice(ii,jj) = Ice_ocean_boundary%aice(i,j)*Grd%tmask(ii,jj,1)
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
        iof_nit(ii,jj) = Ice_ocean_boundary%iof_nit(i,j)*Grd%tmask(ii,jj,1)
        iof_alg(ii,jj) = Ice_ocean_boundary%iof_alg(i,j)*Grd%tmask(ii,jj,1)
#endif
     enddo
  enddo
#endif


  !--------compute surface tracer fluxes from tracer packages------------------- 
  !
#if defined(ACCESS_CM) && defined(CSIRO_BGC)
  call ocean_tpm_sbc(Dom, Grd, T_prog(:), Time, Ice_ocean_boundary%fluxes, runoff, &
                     isc_bnd, iec_bnd, jsc_bnd, jec_bnd,aice, Velocity%u10, &
     use_waterflux, salt_restore_as_salt_flux, atm_co2, co2flux, ocn_co2)
#elif defined(ACCESS_OM) && defined(CSIRO_BGC)
! Do not pass co2flux, ocn_co2 or atm_co2
  call ocean_tpm_sbc(Dom, Grd, T_prog(:), Time, Ice_ocean_boundary%fluxes, runoff, &
                     isc_bnd, iec_bnd, jsc_bnd, jec_bnd,aice=aice, iof_nit=iof_nit, iof_alg=iof_alg, wnd=Velocity%u10, &
                     use_waterflux=use_waterflux, salt_restore_as_salt_flux=salt_restore_as_salt_flux)
#else 
  call ocean_tpm_sbc(Dom, Grd, T_prog(:), Time, Ice_ocean_boundary%fluxes, runoff, &
                     isc_bnd, iec_bnd, jsc_bnd, jec_bnd)
! Leave this case blank at the moment. We want to use the bgc with SIS etc.
#endif

  !--------send diagnostics------------------- 
  !
  call ocean_sbc_diag (Time, Velocity, Thickness, Dens, T_prog, Ice_ocean_boundary,        &
                      pme, runoff, calving, river, alphasfc, betasfc, alphasfc2, betasfc2, &
                      melt, liquid_precip, frozen_precip, evaporation, sensible, longwave, &
                      latent, swflx, swflx_vis)

end subroutine get_ocean_sbc
! </SUBROUTINE> NAME="get_ocean_sbc"



!#######################################################################
! <SUBROUTINE NAME="flux_adjust">
!
! <DESCRIPTION>
! Subroutine to compute the surface fluxes derived from a 
! restoring condition and/or correction from an input file. 
!
! We use a convention whereby a positive 
! flux enters the ocean:  (+) down convention. 
!
! When restoring salinity, one may choose to convert this
! flux to an implied water flux, or keep it a salt flux.
! Converting to a water flux will alter the sea level, and
! so alter the concentration of other tracers.  
! The default is to keep it as a salt flux. 
!
! </DESCRIPTION>
!

subroutine flux_adjust(Time, T_diag, Dens, Ext_mode, T_prog, Velocity, river, melt, pme)
#if defined(ACCESS_CM) || defined(ACCESS_OM)

  use auscom_ice_parameters_mod, only : use_ioaice, aice_cutoff

#endif

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_diag_tracer_type),   intent(in)    :: T_diag(:)
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  type(ocean_velocity_type),      intent(inout) :: Velocity
  real, dimension(isd:,jsd:),     intent(in)    :: river
  real, dimension(isd:,jsd:),     intent(in)    :: melt 
  real, dimension(isd:,jsd:),     intent(inout) :: pme

  real, dimension(isd:ied,jsd:jed) :: open_ocean_mask
  real, dimension(isd:ied,jsd:jed) :: pme_restore, flx_restore
  real, dimension(isd:ied,jsd:jed) :: pme_correct
  real, dimension(isd:ied,jsd:jed) :: flx_correct
  real, dimension(isd:ied,jsd:jed) :: flx_added_heat
  real, dimension(isd:ied,jsd:jed) :: tau_x_correction 
  real, dimension(isd:ied,jsd:jed) :: tau_y_correction 
  real, dimension(isd:ied,jsd:jed) :: pme_eta_restore
  real                             :: pme_eta_restore_total
  real                             :: pme_river_total
  real                             :: pme_restore_total, flx_restore_total 
  real                             :: pme_correct_total, flx_correct_total 
  real                             :: tmp_delta_salinity 
  integer                          :: i, j, k, n, tau, taum1
  logical                          :: used
  logical                          :: ice_present
  real                             :: active_cells, smftu, smftv

#if defined(ACCESS_CM)
  ! Changed in CM2. Make parameter to isolate change
  real, parameter                  :: ice_detection_parameter = 0.054
#else
  real, parameter                  :: ice_detection_parameter = 0.0539
#endif
  
  tau      = Time%tau
  taum1    = Time%taum1  

  pme_eta_restore = 0.0
  pme_correct     = 0.0
  pme_restore     = 0.0
  flx_correct     = 0.0 
  flx_restore     = 0.0
  flx_added_heat  = 0.0 

  open_ocean_mask(isd:ied,jsd:jed) = Grd%tmask(isd:ied,jsd:jed,1)

  ! Update restore_mask if present
  if (id_restore_mask_ofam > 0) then
      call time_interp_external(id_restore_mask_ofam, Time%model_time, &
                                restore_mask)
  end if

  ! add restoring to fluxes from coupled model, or
  ! add flux correction to fluxes from coupled model.
  ! NOTE: (+) down convention


  !-----salinity or pme restoring----------------------------------
  !
  if (id_restore(index_salt) > 0 .and. salt_damp_factor > 0.0) then

      call time_interp_external(id_restore(index_salt), Time%model_time, data)
      if(use_constant_sss_for_restore) then 
        data(:,:) = Grd%tmask(:,:,1)*constant_sss_for_restore
      endif

      ! initialization of restoring fields
      pme_restore = 0.0
      flx_restore = 0.0 

      ! use near-frazil condition as a proxy for where sea-ice is present 
      if(.not. salt_restore_under_ice) then 
          do j=jsc,jec
             do i=isc,iec
                if(Grd%tmask(i,j,1) == 1.0) then 
                    ice_present = T_prog(index_temp)%field(i,j,1,tau) <= &
                        -ice_detection_parameter*T_prog(index_salt)%field(i,j,1,tau)
#if defined(ACCESS_CM) || defined(ACCESS_OM)
                    if (use_ioaice) then
                      ice_present = aice(i,j) >= aice_cutoff
                    endif
#endif
                    if (ice_present) then
                        open_ocean_mask(i,j) = 0.0
                    endif
                endif
             enddo
          enddo
      endif

      if (use_waterflux .and. .not. salt_restore_as_salt_flux) then

          ! put all water restore into pme_restore (no river_restore has been coded) 
          if(max_delta_salinity_restore < 0.0) then 
              do j=jsc,jec
                 do i=isc,iec
                    pme_restore(i,j) = salt_damp_factor*open_ocean_mask(i,j)*restore_mask(i,j) &
                         *(T_prog(index_salt)%field(i,j,1,taum1)-data(i,j))                    &
                         /(T_prog(index_salt)%field(i,j,1,taum1)+epsln)
                 enddo
              enddo
          else 
              do j=jsc,jec
                 do i=isc,iec
                    tmp_delta_salinity = T_prog(index_salt)%field(i,j,1,taum1)-data(i,j)
                    tmp_delta_salinity = sign(1.0,tmp_delta_salinity) &
                                        *min(abs(tmp_delta_salinity),max_delta_salinity_restore)
                    pme_restore(i,j) = salt_damp_factor*open_ocean_mask(i,j)*restore_mask(i,j) &
                         *tmp_delta_salinity                                                   &
                         /(T_prog(index_salt)%field(i,j,1,taum1)+epsln)
                 enddo
              enddo
          endif

          ! produce a zero area average so there is no net input 
          ! of pme mass to the ocean associated with the restoring 
          if(zero_net_water_restore) then 
             pme_restore_total = mpp_global_sum(Dom%domain2d,pme_restore(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1), &
                             global_sum_flag)/Grd%tcellsurf

             ! truncate out the lower order bits if using non_bitwise_reproducible sums
             if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) pme_restore_total = real(pme_restore_total, kind=FLOAT_KIND)
             pme_restore(isc:iec,jsc:jec) = pme_restore(isc:iec,jsc:jec) - pme_restore_total*Grd%tmask(isc:iec,jsc:jec,1)
          endif 

          ! add pme_restore to pme
          pme(isc:iec,jsc:jec) = pme(isc:iec,jsc:jec) + pme_restore(isc:iec,jsc:jec)

          ! produce a zero area average of pme_restore + pme + river.
          ! remove the melt term from sea ice since the coupler 
          ! has already included the ice melt within the pme field. 
          pme_river(:,:) = 0.0
          if(zero_net_water_couple_restore) then 
              do j=jsc,jec
                 do i=isc,iec
                    pme_river(i,j) = pme(i,j) + river(i,j) - melt(i,j)
                 enddo
              enddo
              pme_river_total = mpp_global_sum(Dom%domain2d,pme_river(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1), &
                                global_sum_flag)/Grd%tcellsurf

              ! truncate out the lower order bits if using non_bitwise_reproducible sums
              if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) pme_river_total = real(pme_river_total, kind=FLOAT_KIND)
              do j=jsc,jec
                 do i=isc,iec
                    pme(i,j)  = pme(i,j) - pme_river_total*Grd%tmask(i,j,1)
                 enddo
              enddo
          endif
          call mpp_update_domains(pme(:,:), Dom%domain2d)


      else  ! for opposite of if (use_waterflux .and. .not. salt_restore_as_salt_flux) then

          if(max_delta_salinity_restore < 0.0) then 
              do j=jsc,jec
                 do i=isc,iec
                    flx_restore(i,j) = salt_damp_factor*open_ocean_mask(i,j)*restore_mask(i,j) &
                                       *(data(i,j) - T_prog(index_salt)%field(i,j,1,taum1))
                 enddo
              enddo
          else 
              do j=jsc,jec
                 do i=isc,iec
                    tmp_delta_salinity = data(i,j) - T_prog(index_salt)%field(i,j,1,taum1)
                    ! Only limit salinity restoring within specified salinity tolerances
                    if (T_prog(index_salt)%field(i,j,1,taum1) > salinity_restore_limit_lower &
                         .and.  T_prog(index_salt)%field(i,j,1,taum1) < salinity_restore_limit_upper ) then
                       tmp_delta_salinity = sign(1.0,tmp_delta_salinity) &
                            *min(abs(tmp_delta_salinity),max_delta_salinity_restore)
                    end if
                    flx_restore(i,j) = salt_damp_factor*open_ocean_mask(i,j)*restore_mask(i,j) &
                                       *tmp_delta_salinity
                 enddo
              enddo
          endif

          ! Restrict outgoing salt flux (ice being formed) in cases where salinity falls below
          ! ocean_ice_salt_limit. If zero_net_salt_restore=.true. then the mismatch will be
          ! spread globally by the following block of code.

          if( ocean_ice_salt_limit > 0.0 ) then
             do j=jsc,jec
                do i=isc,iec
                   if (T_prog(index_salt)%stf(i,j) < 0.0 .and. T_prog(index_salt)%field(i,j,1,taum1) < ocean_ice_salt_limit * 1000.0) then
                        flx_restore(i,j) = flx_restore(i,j) - T_prog(index_salt)%stf(i,j)
                   endif
                enddo
             enddo
          endif

          ! produce a zero area average so there is no net input
          ! of salt to the ocean associated with the restoring 
          if(zero_net_salt_restore) then 
            flx_restore_total =  &
            mpp_global_sum(Dom%domain2d,flx_restore(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1), global_sum_flag)/Grd%tcellsurf

            ! truncate out the lower order bits if using non_bitwise_reproducible sums
            if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) flx_restore_total = real(flx_restore_total, kind=FLOAT_KIND)
            flx_restore(isc:iec,jsc:jec) = flx_restore(isc:iec,jsc:jec) - flx_restore_total*Grd%tmask(isc:iec,jsc:jec,1)
          endif 

          ! add flx_restore to stf
          do j = jsc,jec
             do i = isc,iec
                T_prog(index_salt)%stf(i,j) = T_prog(index_salt)%stf(i,j) + flx_restore(i,j)
             enddo
          enddo

#if defined(CSIRO_BGC)
          ! if salt fluxes are used to restore salinity, then virtual fluxes are
          ! needed for csiro BGC tracers. mac, dec12.
          if (do_csiro_bgc) then
            call csiro_bgc_virtual_fluxes(isc, iec, jsc, jec, isd, ied, jsd, jed, flx_restore, T_prog)
          endif
#endif

      endif  ! endif for if (use_waterflux .and. .not. salt_restore_as_salt_flux) then

  endif  ! endif for if(id_restore(index_salt) > 0 .and. salt_damp_factor > 0.0)



  !-----salinity or pme flux correction--------------------------
  ! 
  ! add salt fluxes from a data file to perform a flux correction
  if ( index_salt > 0 ) then
    if (id_correction(index_salt) > 0 ) then

      call time_interp_external(id_correction(index_salt), Time%model_time, data)

      ! initialization of correction fields 
      pme_correct = 0.0
      flx_correct = 0.0 

      ! salt flux correction is assumed to be a pme_mass field (units kg/(m^2*sec))
      if (use_waterflux) then

          do j = jsc,jec
             do i = isc,iec
                pme_correct(i,j) = Grd%tmask(i,j,1)*salt_correction_scale*data(i,j)
             enddo
          enddo     

          ! produce a zero area average so there is no net input 
          ! of water to the ocean associated with the flux correction 
          if(zero_net_water_correction) then 
             pme_correct_total = mpp_global_sum(Dom%domain2d,pme_correct(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1), &
                                 global_sum_flag)/Grd%tcellsurf

             ! truncate out the lower order bits if using non_bitwise_reproducible sums
             if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) pme_correct_total = real(pme_correct_total, kind=FLOAT_KIND)
             pme_correct(isc:iec,jsc:jec) = pme_correct(isc:iec,jsc:jec) - pme_correct_total*Grd%tmask(isc:iec,jsc:jec,1)
          endif 

          ! add pme_correct to pme
          pme(isc:iec,jsc:jec) = pme(isc:iec,jsc:jec) + pme_correct(isc:iec,jsc:jec)
          call mpp_update_domains(pme(:,:), Dom%domain2d)

      else  ! salt fluxes 

          ! need to convert the pme correction to an implied salt flux  
          do j = jsc,jec
             do i = isc,iec
                flx_correct(i,j) = Grd%tmask(i,j,1)*salt_correction_scale*data(i,j)/(rho0*0.001) 
             enddo
          enddo
         
          ! produce a zero area average so there is no net input
          ! of salt to the ocean associated with the flux correction 
          if(zero_net_salt_correction) then 
             flx_correct_total =  &
             mpp_global_sum(Dom%domain2d,flx_correct(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1), global_sum_flag)/Grd%tcellsurf

             ! truncate out the lower order bits if using non_bitwise_reproducible sums
             if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) flx_correct_total = real(flx_correct_total, kind=FLOAT_KIND)
             flx_correct(isc:iec,jsc:jec) = flx_correct(isc:iec,jsc:jec) - flx_correct_total*Grd%tmask(isc:iec,jsc:jec,1)
          endif 

          ! add flx_correct to stf
          do j = jsc,jec
             do i = isc,iec
                T_prog(index_salt)%stf(i,j) = T_prog(index_salt)%stf(i,j) + flx_correct(i,j)
             enddo
          enddo

      endif 

   endif

  endif   ! endif for if (id_correction(index_salt) > 0 )


  ! diagnostics for salinity or pme restoring and correction 

  ! salt from restoring
  if ( index_salt > 0 ) then
    if (id_stf_restore(index_salt) > 0) then
     call diagnose_2d(Time, Grd, id_stf_restore(index_salt), flx_restore(:,:)*T_prog(index_salt)%conversion)
    endif
  endif
  ! total salt from restoring
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_restore(index_salt), flx_restore, 1e-15)

  ! salt from correction 
  if ( index_salt > 0 ) then
    if (id_stf_correct(index_salt) > 0) then
     call diagnose_2d(Time, Grd, id_stf_correct(index_salt), flx_correct(:,:)*T_prog(index_salt)%conversion)
    endif
  endif
  ! total salt from correction 
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_correct(index_salt), flx_correct, 1e-15)

  if(id_tform_rho_pbl_adjsalt_on_nrho > 0) then
      wrk1(:,:,:)      = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = (flx_restore(i,j)+flx_correct(i,j))*Dens%drhodS(i,j,k) &
                 *Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_tform_rho_pbl_adjsalt_on_nrho, wrk1)
  endif

  if(id_mass_pme_adj_on_nrho > 0) then
      wrk1(:,:,:)      = 0.0
      k=1 
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%dat(i,j)*(pme_restore(i,j)+pme_correct(i,j))
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_mass_pme_adj_on_nrho, wrk1)
  endif

  ! salt from all surface fluxes 
  if (id_stf_total(index_salt) > 0) then
     call diagnose_2d(Time, Grd, id_stf_total(index_salt), T_prog(index_salt)%stf(:,:)*T_prog(index_salt)%conversion)
  endif
  ! total salt from all fluxes 
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_sum(index_salt), T_prog(index_salt)%stf, 1e-15*T_prog(index_salt)%conversion)

  ! pme from salt restoring
  call diagnose_2d(Time, Grd, id_pme_restore, pme_restore(:,:))
  ! total pme from salt restoring
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_pme_restore, pme_restore, 1e-15)

  ! pme from salt correction 
  call diagnose_2d(Time, Grd, id_pme_correct, pme_correct(:,:))
  ! total pme from salt correction 
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_pme_correct, pme_correct, 1e-15)

  ! pme from all surface terms 
  call diagnose_2d(Time, Grd, id_pme_net, pme(:,:))
  ! total pme from all surface terms 
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_pme_net, pme, 1e-15)

  ! heat input from net pme relative to 0 degrees C (W/m2)
  if (id_stf_pme(index_temp) > 0) then
     call diagnose_2d(Time, Grd, id_stf_pme(index_temp),       &
             pme(:,:)*T_prog(index_temp)%tpme(:,:)*T_prog(index_temp)%conversion)
  endif
  ! heat input from net pme relative to 0 degrees C (W/m2) binned to
  ! neutral density
  if (id_stf_pme_on_nrho(index_temp) > 0) then
      wrk1(:,:,:)      = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = pme(i,j)*T_prog(index_temp)%tpme(i,j)*T_prog(index_temp)%conversion
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_stf_pme_on_nrho(index_temp), wrk1)
  endif
  ! total heat flux from net pme (Watts)
  if(id_total_ocean_stf_pme(index_temp) > 0) then 
     call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_pme(index_temp), pme(:,:)*T_prog(index_temp)%tpme(:,:), 1e-15*T_prog(index_temp)%conversion)
  endif

  if (id_ice_mask > 0) then
     call diagnose_2d(Time, Grd, id_ice_mask, (1.0-open_ocean_mask(:,:))*Grd%tmask(:,:,1))
  endif

  call diagnose_2d(Time, Grd, id_open_ocean_mask, open_ocean_mask(:,:))
  
  ! contribution to sea level from salt restoring flux 
  if(id_eta_tend_salt_restore > 0 .or. id_eta_tend_salt_restore_glob > 0) then 
      wrk1_2d = 0.0
      do j=jsc,jec
         do i=isc,iec
            wrk1_2d(i,j) = -betasfc(i,j)*flx_restore(i,j)*rhosfc_inv(i,j) 
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_salt_restore, wrk1_2d(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_salt_restore_glob, wrk1_2d, cellarea_r)
  endif

  ! contribution to sea level from water restoring flux 
  if(id_eta_tend_water_restore > 0 .or. id_eta_tend_water_restore_glob > 0) then 
      wrk1_2d = 0.0
      do j=jsc,jec
         do i=isc,iec
            wrk1_2d(i,j) = pme_restore(i,j)*rhosfc_inv(i,j) 
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_water_restore, wrk1_2d(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_water_restore_glob, wrk1_2d, cellarea_r)
  endif



  !-------temperature restoring and flux correction----------------
  !
  ! initialization of flux fields for temperature restoring and correction
  flx_correct = 0.0 
  flx_restore = 0.0 

  ! temperature restoring   
  if (id_restore(index_temp) > 0 .and. temp_damp_factor > 0.0) then

     call time_interp_external(id_restore(index_temp), Time%model_time, data)
     if(use_constant_sst_for_restore) then 
       data(:,:) = Grd%tmask(:,:,1)*constant_sst_for_restore
     endif


     if(prog_temp_variable==CONSERVATIVE_TEMP) then 
         do j=jsc,jec
            do i=isc,iec
               flx_restore(i,j) = temp_damp_factor*restore_mask(i,j) &
                                  *(data(i,j)-T_diag(index_diag_temp)%field(i,j,1))
               T_prog(index_temp)%stf(i,j) = T_prog(index_temp)%stf(i,j) + flx_restore(i,j)
            enddo
         enddo
     else 
         do j=jsc,jec
            do i=isc,iec
               flx_restore(i,j) = temp_damp_factor*restore_mask(i,j) &
                                  *(data(i,j)-T_prog(index_temp)%field(i,j,1,taum1))
               T_prog(index_temp)%stf(i,j) = T_prog(index_temp)%stf(i,j) + flx_restore(i,j)
            enddo
         enddo
     endif

  endif ! (id_restore(index_temp) > 0 .and. temp_damp_factor > 0.0)


  ! temperature flux correction 
  if ( index_temp > 0 ) then

    if (id_correction(index_temp) > 0) then

       call time_interp_external(id_correction(index_temp), Time%model_time, data)

       ! assume flux correction in units of W/m2, so divide by cp_ocean to 
       ! convert to a temperature flux in units of (rho*dz/time)*(deg C) 
       do j=jsc,jec
          do i=isc,iec
             flx_correct(i,j) = Grd%tmask(i,j,1)*data(i,j)*temp_correction_scale*cp_ocean_r
             T_prog(index_temp)%stf(i,j) = T_prog(index_temp)%stf(i,j) + flx_correct(i,j)
          enddo
       enddo

    endif

  endif

  ! FAFMIP added heat flux correction 
  if ( index_added_heat > 0 ) then
    if (id_correction(index_added_heat) > 0) then

       call time_interp_external(id_correction(index_added_heat), Time%model_time, data)

       ! assume FAFMIP added heat in units of W/m2, so divide by cp_ocean to 
       ! convert to a temperature flux in units of (rho*dz/time)*(deg C) 
       do j=jsc,jec
          do i=isc,iec
             flx_added_heat(i,j) = Grd%tmask(i,j,1)*data(i,j)*temp_correction_scale*cp_ocean_r
             T_prog(index_added_heat)%stf(i,j) = T_prog(index_added_heat)%stf(i,j) + flx_added_heat(i,j)
          enddo
       enddo
    endif

  endif


  ! diagnostics for temperature restoring and flux correction

  ! restoring heat flux 
  if (id_stf_restore(index_temp) > 0) then
     call diagnose_2d(Time, Grd, id_stf_restore(index_temp), flx_restore(:,:)*T_prog(index_temp)%conversion)
  endif

  ! total of restoring heat flux 
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_restore(index_temp), flx_restore(:,:), 1e-15*T_prog(index_temp)%conversion)

  ! flux correction heat flux 
  if ( index_temp > 0 ) then
    if (id_stf_correct(index_temp) > 0) then
      call diagnose_2d(Time, Grd, id_stf_correct(index_temp), flx_correct(:,:)*T_prog(index_temp)%conversion)
    endif
  endif

  ! FAFMIP added heat  
  if ( index_added_heat > 0 ) then
    if (id_stf_correct(index_added_heat) > 0) then
      call diagnose_2d(Time, Grd, id_stf_correct(index_added_heat), flx_added_heat(:,:)*T_prog(index_added_heat)%conversion)
    endif
  endif

  ! total of flux correction heat flux 
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_correct(index_temp), flx_correct(:,:), 1e-15*T_prog(index_temp)%conversion)

  ! total of added heat for FAFMIP
  if ( index_added_heat > 0 ) then
    call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_correct(index_added_heat), flx_added_heat(:,:), 1e-15*T_prog(index_added_heat)%conversion)
  endif

  if(id_tform_rho_pbl_adjheat_on_nrho > 0) then
      wrk1(:,:,:) = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = (flx_restore(i,j)+flx_correct(i,j))*Dens%drhodT(i,j,k) &
                          *Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_tform_rho_pbl_adjheat_on_nrho, wrk1)
  endif

  ! net heat from stf   
  if (id_stf_total(index_temp) > 0) then
     call diagnose_2d(Time, Grd, id_stf_total(index_temp),              &
             T_prog(index_temp)%stf(:,:)*T_prog(index_temp)%conversion)
  endif

  ! total of net heat flux 
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_sum(index_temp), T_prog(index_temp)%stf(:,:), 1e-15*T_prog(index_temp)%conversion)


  ! contribution to sea level from temperature restoring flux 
  if(id_eta_tend_heat_restore > 0 .or. id_eta_tend_heat_restore_glob > 0) then 
      wrk1_2d = 0.0
      do j=jsc,jec
         do i=isc,iec
            wrk1_2d(i,j) = alphasfc(i,j)*flx_restore(i,j)*rhosfc_inv(i,j) 
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_heat_restore, wrk1_2d(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_heat_restore_glob, wrk1_2d, cellarea_r)
  endif



  !---------wind flux correction-----------------------------------
  tau_x_correction = 0.0
  tau_y_correction = 0.0
  wrk1_2d          = 0.0  
  wrk2_2d          = 0.0

  ! tau_x flux correction (N/m2), assumed to be on B-grid
  if (id_tau_x_correction > 0 ) then
     call time_interp_external(id_tau_x_correction, Time%model_time, data)
     do j=jsc,jec
        do i=isc,iec
           tau_x_correction(i,j) = tau_x_correction_scale*Grd%umask(i,j,1)*data(i,j)
        enddo
     enddo
  endif

  ! tau_y flux correction (N/m2), assumed to be on B-grid
  if (id_tau_y_correction > 0 ) then
     call time_interp_external(id_tau_y_correction, Time%model_time, data)
     do j=jsc,jec
        do i=isc,iec
           tau_y_correction(i,j) = tau_y_correction_scale*Grd%umask(i,j,1)*data(i,j)
        enddo
     enddo
  endif

  ! mpp update domain is needed for vertical mixing schemes such
  ! as KPP and GOTM, as well as to get c-grid version of smf.
  if(id_tau_x_correction > 0 .or. id_tau_y_correction > 0) then 
     call mpp_update_domains(tau_x_correction(:,:),tau_y_correction(:,:),Dom%domain2d,gridtype=BGRID_NE)

     do j=jsd,jed
        do i=isd,ied
           Velocity%smf_bgrid(i,j,1) = Velocity%smf_bgrid(i,j,1) + tau_x_correction(i,j)
           Velocity%smf_bgrid(i,j,2) = Velocity%smf_bgrid(i,j,2) + tau_y_correction(i,j)
        enddo
     enddo

     ! c-grid version of correction 
     do j=jsc,jec
        do i=isc,iec
           if(Grd%umask(i,j,1) + Grd%umask(i,j-1,1) > 0.0) then
               wrk1_2d(i,j) =                                                                        &
               (Grd%umask(i,j,1)*tau_x_correction(i,j) + Grd%umask(i,j-1,1)*tau_x_correction(i,j-1)) &
               / (Grd%umask(i,j,1) + Grd%umask(i,j-1,1))
           endif
           if(Grd%umask(i,j,1) + Grd%umask(i-1,j,1) > 0.0) then
               wrk2_2d(i,j) =                                                                        &
               (Grd%umask(i,j,1)*tau_y_correction(i,j) + Grd%umask(i-1,j,1)*tau_y_correction(i-1,j)) &
               / (Grd%umask(i,j,1) + Grd%umask(i-1,j,1))
           endif
        enddo
     enddo

     do j=jsc,jec
        do i=isc,iec
           Velocity%smf_cgrid(i,j,1) = Velocity%smf_cgrid(i,j,1) + wrk1_2d(i,j)
           Velocity%smf_cgrid(i,j,2) = Velocity%smf_cgrid(i,j,2) + wrk2_2d(i,j)
        enddo
     enddo

     ! for use in forcing momentum
     if(horz_grid == MOM_BGRID) then
        do n=1,2
           do j=jsd,jed
              do i=isd,ied
                 Velocity%smf(i,j,n) = Velocity%smf_bgrid(i,j,n)
              enddo
           enddo
        enddo
     else
        do n=1,2
           do j=jsc,jec
              do i=isc,iec
                 Velocity%smf(i,j,n) = Velocity%smf_cgrid(i,j,n)
              enddo
           enddo
        enddo
     endif

  ! Calculate surface friction velocity
  ! FAFMIP stess experiment says do not do this. Default is .TRUE. See notes for
  ! ocean_sbc_nml namelist above

     if ( do_ustar_correction ) then
          do j=jsc,jec
              do i=isc,iec

!         ustar is needed on the "T-grid".  It is assumed that masking of 
!         smf over land was performed inside of the ocean_sbc module. 
!         smf has units of N/m^2 and so we need rho0r to get ustar in m/s.   
!         swflx has units W/m^2 so needs 1/rho_cp to get to C*m/s units.
!         these are proper units for buoyancy fluxes. 
             active_cells = Grd%umask(i,j,1)   + Grd%umask(i-1,j,1)   &
                        +Grd%umask(i,j-1,1) + Grd%umask(i-1,j-1,1) + epsln
             smftu = rho0r*(Velocity%smf_bgrid(i,j,1)   + Velocity%smf_bgrid(i-1,j,1)     &
                        +Velocity%smf_bgrid(i,j-1,1) + Velocity%smf_bgrid(i-1,j-1,1))  &
                  /active_cells
             smftv = rho0r*(Velocity%smf_bgrid(i,j,2)   + Velocity%smf_bgrid(i-1,j,2)    &
                        +Velocity%smf_bgrid(i,j-1,2) + Velocity%smf_bgrid(i-1,j-1,2)) &
                   /active_cells
             Velocity%ustar(i,j) = sqrt( sqrt(smftu**2 + smftv**2) )
           enddo
        enddo
     endif
  

  endif

  if(horz_grid == MOM_BGRID) then 

     ! net i-directed wind stress getting into the ocean (N/m2)
     if (id_tau_x_net > 0) then 
          used =  send_data(id_tau_x_net, Velocity%smf_bgrid(:,:,1),&
                  Time%model_time, rmask=wind_mask(:,:),            &
                  is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif 

     ! net j-directed wind stress getting into the ocean (N/m2)
     if (id_tau_y_net > 0) then 
          used =  send_data(id_tau_y_net, Velocity%smf_bgrid(:,:,2),&
                  Time%model_time, rmask=wind_mask(:,:),            &
                  is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif 

  else ! c-grid 

     ! net i-directed wind stress getting into the ocean (N/m2)
     if (id_tau_x_net > 0) then 
          used =  send_data(id_tau_x_net, Velocity%smf_cgrid(:,:,1),&
                  Time%model_time, rmask=wind_mask(:,:),            &
                  is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif 

     ! net j-directed wind stress getting into the ocean (N/m2)
     if (id_tau_y_net > 0) then 
          used =  send_data(id_tau_y_net, Velocity%smf_cgrid(:,:,2),&
                  Time%model_time, rmask=wind_mask(:,:),            &
                  is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif 

     tau_x_correction(:,:) = wrk1_2d(:,:)
     tau_y_correction(:,:) = wrk2_2d(:,:)

  endif 

  ! diagnostics for i-directed wind stress correction (N/m^2)
  if (id_tau_x_flux_correction > 0) then 
       used =  send_data(id_tau_x_flux_correction, tau_x_correction(:,:),&
               Time%model_time, rmask=wind_mask(:,:),                    &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 

  ! diagnostics for j-directed wind stress correction (N/m^2)
  if (id_tau_y_flux_correction > 0) then 
       used =  send_data(id_tau_y_flux_correction, tau_y_correction(:,:),&
               Time%model_time, rmask=wind_mask(:,:),                    &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 



  !-------------sea level restoring----------------------------
  !
  ! flux fields for eta restoring to modify water flux
  pme_eta_restore = 0.0
  if (id_eta_restore > 0 .and. eta_damp_factor > 0.0) then

     ! rho0 factor converts to mass flux 
     call time_interp_external(id_eta_restore, Time%model_time, data)
     do j=jsc,jec
        do i=isc,iec
           pme_eta_restore(i,j) = rho0*eta_damp_factor*restore_mask(i,j) &
                                  *(data(i,j)-Ext_mode%eta_t(i,j,taum1))
        enddo
     enddo

     if(zero_net_pme_eta_restore) then 
         pme_eta_restore_total = mpp_global_sum(Dom%domain2d,pme_eta_restore(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1), &
                                 global_sum_flag)/Grd%tcellsurf
         ! truncate out the lower order bits if using non_bitwise_reproducible sums
         if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) pme_eta_restore_total = real(pme_eta_restore_total, kind=FLOAT_KIND)
         pme_eta_restore(isc:iec,jsc:jec) = pme_eta_restore(isc:iec,jsc:jec) - pme_eta_restore_total*Grd%tmask(isc:iec,jsc:jec,1)
    endif 

    ! add pme_restore to pme
    pme(isc:iec,jsc:jec) = pme(isc:iec,jsc:jec) + pme_eta_restore(isc:iec,jsc:jec)

    call diagnose_2d(Time, Grd, id_pme_eta_restore, pme_eta_restore(:,:))

  endif ! (id_eta_restore > 0 .and. eta_damp_factor > 0.0)



  return

end subroutine flux_adjust
! </SUBROUTINE> NAME="flux_adjust"



!#######################################################################
! <SUBROUTINE NAME="ocean_sbc_diag">
!
! <DESCRIPTION>
! Compute and send diagnostics from get_ocean_sbc. 
! </DESCRIPTION>
!
subroutine ocean_sbc_diag(Time, Velocity, Thickness, Dens, T_prog, Ice_ocean_boundary,     &
                      pme, runoff, calving, river, alphasfc, betasfc, alphasfc2, betasfc2, &
                      melt, liquid_precip,  frozen_precip, evaporation, sensible, longwave,&
                      latent, swflx, swflx_vis)

  type(ocean_time_type),          intent(in) :: Time 
  type(ocean_velocity_type),      intent(in) :: Velocity
  type(ocean_thickness_type),     intent(in) :: Thickness
  type(ocean_density_type),       intent(in) :: Dens
  type(ocean_prog_tracer_type),   intent(in) :: T_prog(:)
  type(ice_ocean_boundary_type),  intent(in) :: Ice_ocean_boundary
  real, dimension(isd:,jsd:),     intent(in) :: pme
  real, dimension(isd:,jsd:),     intent(in) :: runoff
  real, dimension(isd:,jsd:),     intent(in) :: calving 
  real, dimension(isd:,jsd:),     intent(in) :: river
  real, dimension(isd:,jsd:),     intent(in) :: alphasfc
  real, dimension(isd:,jsd:),     intent(in) :: betasfc
  real, dimension(isd:,jsd:),     intent(in) :: alphasfc2
  real, dimension(isd:,jsd:),     intent(in) :: betasfc2
  real, dimension(isd:,jsd:),     intent(in) :: melt
  real, dimension(isd:,jsd:),     intent(in) :: liquid_precip
  real, dimension(isd:,jsd:),     intent(in) :: frozen_precip
  real, dimension(isd:,jsd:),     intent(in) :: evaporation 
  real, dimension(isd:,jsd:),     intent(in) :: sensible 
  real, dimension(isd:,jsd:),     intent(in) :: longwave
  real, dimension(isd:,jsd:),     intent(in) :: latent
  real, dimension(isd:,jsd:),     intent(in) :: swflx
  real, dimension(isd:,jsd:),     intent(in) :: swflx_vis

  real, dimension(isd:ied,jsd:jed) :: tmp_flux

  integer :: i,j,k 
  integer :: ii, jj
  integer :: tau
  real    :: total_stuff 
  real    :: umask_norm
  real    :: totz, tbrz, maskt, factor 

  tau = Time%tau 


  !----momentum related diagnostics--------------------------------

  if(horz_grid == MOM_BGRID) then 

     ! i-directed wind stress (N/m^2)
     if (id_tau_x > 0) used =  send_data(id_tau_x, Velocity%smf_bgrid(:,:,1), &
                               Time%model_time, rmask=wind_mask(:,:),         &
                               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

     ! j-directed wind stress (N/m^2)
     if (id_tau_y > 0) used =  send_data(id_tau_y, Velocity%smf_bgrid(:,:,2), &
                               Time%model_time, rmask=wind_mask(:,:),         &
                               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  

  else 

     ! i-directed wind stress (N/m^2)
     if (id_tau_x > 0) used =  send_data(id_tau_x, Velocity%smf_cgrid(:,:,1), &
                               Time%model_time, rmask=wind_mask(:,:),         &
                               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

     ! j-directed wind stress (N/m^2)
     if (id_tau_y > 0) used =  send_data(id_tau_y, Velocity%smf_cgrid(:,:,2), &
                               Time%model_time, rmask=wind_mask(:,:),         &
                               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  

  endif 

  !----Langmuir turbulence related diagnostics-------------------------------
  ! i-directed stokes drift velocity (m/s)
  if (do_langmuir) then
     if (id_ustoke > 0) used =  send_data(id_ustoke, Velocity%ustoke(:,:), &
                            Time%model_time, rmask=Grd%umask(:,:,1), &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

     ! j-directed stokes drift velocity (m/s)
     if (id_vstoke > 0) used =  send_data(id_vstoke, Velocity%vstoke(:,:), &
                            Time%model_time, rmask=Grd%umask(:,:,1), &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

     ! mean wave length (m)
     if (id_wavlen > 0) used =  send_data(id_wavlen, Velocity%wavlen(:,:), &
                            Time%model_time, rmask=Grd%tmask(:,:,1), &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! wind stress curl (N/m^3) averaged to U-point
  ! Ekman pumping velocity averaged to U-point 
  if (id_tau_curl > 0 .or. id_ekman_we > 0) then 
      wrk2_2d(:,:) = 0.0
      wrk3_2d(:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec
            umask_norm   = Grd%umask(i,j,1)/(epsln+Grd%umask(i+1,j,1)+Grd%umask(i-1,j,1))
            wrk2_2d(i,j) = umask_norm*(                                                                         &
                  Grd%umask(i+1,j,1)*( Velocity%smf_bgrid(i+1,j,2)-Velocity%smf_bgrid(i,j,2))  /Grd%dxtn(i+1,j) & 
                 +Grd%umask(i-1,j,1)*( Velocity%smf_bgrid(i,j,2)  -Velocity%smf_bgrid(i-1,j,2))/Grd%dxtn(i,j)   &
                 )
            umask_norm   = Grd%umask(i,j,1)/(epsln+Grd%umask(i,j+1,1)+Grd%umask(i,j-1,1))
            wrk2_2d(i,j) = wrk2_2d(i,j)                                                                         &
                 -umask_norm*(                                                                                  &
                  Grd%umask(i,j+1,1)*( Velocity%smf_bgrid(i,j+1,1)-Velocity%smf_bgrid(i,j,1))  /Grd%dyte(i,j+1) & 
                 +Grd%umask(i,j-1,1)*( Velocity%smf_bgrid(i,j,1)  -Velocity%smf_bgrid(i,j-1,1))/Grd%dyte(i,j)   &
                 )
            if(abs(Grd%yu(i,j)) > 1.0) then 
              wrk3_2d(i,j) = rho0r*wrk2_2d(i,j)/Grd%f(i,j) 
            endif 
         enddo
      enddo
      call diagnose_2d_u(Time, Grd, id_tau_curl, wrk2_2d(:,:))
      call diagnose_2d_u(Time, Grd, id_ekman_we, wrk3_2d(:,:))
  endif 


  ! Ekman Heat Transport component
  ! Defined as: Hek = -rho_cp \Int[ Taux/(f rho0) (T1 - <T>) ]dx
  ! where T1=first layer Temp; < >=vertical mean.
  ! Zonal sum done off-line.
  ! TAUx = smf*rho0, we want TAUx/rho0, therefore we use smf.  
  ! Algorithm from riccardo.farneti
  if (id_ekman_heat > 0) then

      wrk1_2d(:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec

            factor = 0.0
            if(Grd%f(i,j)==0.0 .and. j>1) then
                factor = 4.*Grd%f(i,j-1)
            else
                factor = 4.*Grd%f(i,j)
            endif

            totz = 1.0
            tbrz = 0.0
            maskt= 0.0
            do k=1,nk 
               maskt = Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
               tbrz  = tbrz + min(Thickness%dzt(i,j,k),Thickness%dzt(i,j+1,k))*Thickness%dzt(i,j,k)*maskt &
                       *(T_prog(index_temp)%field(i,j,k,tau)+T_prog(index_temp)%field(i,j+1,k,tau))
               totz  = totz + min(Thickness%dzt(i,j,k),Thickness%dzt(i,j+1,k))*Thickness%dzt(i,j,k)*maskt 
            enddo
            if (totz /= 0.0) then
                tbrz = tbrz/totz
                wrk1_2d(i,j) = wrk1_2d(i,j)                                                                 &
                   -( Velocity%smf_bgrid(i,j,1)*Grd%dxu(i,j) + Velocity%smf_bgrid(i-1,j,1)*Grd%dxu(i-1,j) ) &
                  *( T_prog(index_temp)%field(i,j,1,tau) + T_prog(index_temp)%field(i,j+1,1,tau) -tbrz )    &
                  *cos(Grd%phiu(i,j)) / factor    
            endif

         enddo
      enddo
      call diagnose_2d(Time, Grd, id_ekman_heat, wrk1_2d(:,:)*rho_cp)

  endif


  !--------stokes drift velocity and decay depth----------------------
  call diagnose_3d_u(Time, Grd, id_ustokes, Velocity%stokes_drift(:,:,:,1))
  call diagnose_3d_u(Time, Grd, id_vstokes, Velocity%stokes_drift(:,:,:,2))
  call diagnose_2d_u(Time, Grd, id_stokes_depth, Velocity%stokes_depth(:,:))


  !--------runoff/calving/river related diagnostics----------------------
  !
  ! temperature of runoff
  if (id_trunoff(index_temp) > 0) then

      wrk1_2d=0.0
      if(land_model_heat_fluxes) then 
          ! temp_runoff = heat_runoff/(cp_liquid_runoff*runoff).
          ! this is the temperature that the land model would 
          ! call its runoff temperature. 
          do j=jsc,jec
             do i=isc,iec
                if(runoff(i,j) /= 0.0) then 
                    wrk1_2d(i,j) = T_prog(index_temp)%trunoff(i,j)*cp_liquid_runoff_r*cp_ocean
                endif
             enddo
          enddo

      else 

          do j=jsc,jec
             do i=isc,iec
                if(runoff(i,j) /= 0.0) then 
                    wrk1_2d(i,j) = T_prog(index_temp)%trunoff(i,j)
                endif
             enddo
          enddo

      endif

      call diagnose_2d(Time, Grd, id_trunoff(index_temp), wrk1_2d(:,:))
  endif


  ! salinity of runoff
  if (id_trunoff(index_salt) > 0) then
      wrk1_2d=0.0
      do j=jsc,jec
         do i=isc,iec
            if(runoff(i,j) /= 0.0) wrk1_2d(i,j) = 1.0
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_trunoff(index_salt),  &
             wrk1_2d(:,:)*T_prog(index_salt)%trunoff(:,:))
  endif

  ! temp of calving solid runoff 
  if (id_tcalving(index_temp) > 0) then

      wrk1_2d=0.0
      if(land_model_heat_fluxes) then 
          ! temp_calving = heat_calving/(cp_solid_runoff*calving).
          ! this is the temperature that the land model would 
          ! call its calving land ice temperature. 
          do j=jsc,jec
             do i=isc,iec
                if(calving(i,j) /= 0.0) then 
                    wrk1_2d(i,j) = T_prog(index_temp)%tcalving(i,j)*cp_solid_runoff_r*cp_ocean
                endif
             enddo
          enddo

      else 

          do j=jsc,jec
             do i=isc,iec
                if(calving(i,j) /= 0.0) then 
                    wrk1_2d(i,j) = T_prog(index_temp)%tcalving(i,j)
                endif
             enddo
          enddo

      endif

      call diagnose_2d(Time, Grd, id_tcalving(index_temp), wrk1_2d(:,:))
  endif

  ! effective temp of calving, computed as 
  ! temp_calve = heat_calve/(cp_ocean*calve)
  if (id_temp_calving_eff > 0) then
      wrk1_2d=0.0
      do j=jsc,jec
         do i=isc,iec
            if(calving(i,j) /= 0.0) then 
                wrk1_2d(i,j) = 1.0
            endif 
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_temp_calving_eff,      &
             wrk1_2d(:,:)*T_prog(index_temp)%tcalving(:,:))
  endif

  ! salinity of calving
  if (id_tcalving(index_salt) > 0) then
      wrk1_2d=0.0
      do j=jsc,jec
         do i=isc,iec
            if(calving(i,j) /= 0.0) wrk1_2d(i,j) = 1.0
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_tcalving(index_salt), &
             wrk1_2d(:,:)*T_prog(index_salt)%tcalving(:,:))
  endif

  ! temperature of river 
  if (id_triver(index_temp) > 0) then
      wrk1_2d=0.0
      do j=jsc,jec
         do i=isc,iec
            if(river(i,j) /= 0.0) wrk1_2d(i,j) = 1.0
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_triver(index_temp), &
             wrk1_2d(:,:)*T_prog(index_temp)%triver(:,:))
  endif

  ! salinity of river 
  if (id_triver(index_salt) > 0) then
      wrk1_2d=0.0
      do j=jsc,jec
         do i=isc,iec
            if(river(i,j) /= 0.0) wrk1_2d(i,j) = 1.0
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_triver(index_salt), &
             wrk1_2d(:,:)*T_prog(index_salt)%triver(:,:))
  endif


  !--------heat related diagnostics ------------------------------------
  !
  ! latent heat of vaporization
  call diagnose_2d(Time, Grd, id_latent_heat_vapor, latent_heat_vapor(:,:))
  
  ! latent heat of fusion 
  call diagnose_2d(Time, Grd, id_latent_heat_fusion, latent_heat_fusion(:,:))
  
  ! surface heat flux (W/m2) passed through the coupler 
  if (id_stf_coupler(index_temp) > 0) then
     call diagnose_2d(Time, Grd, id_stf_coupler(index_temp),           &
             T_prog(index_temp)%stf(:,:)*T_prog(index_temp)%conversion)
  endif

  ! FAFMIP sft passed through the coupler and/or through flux correction 
  if ( index_added_heat > 0 ) then
    if (id_stf_coupler(index_added_heat) > 0) then
       call diagnose_2d(Time, Grd, id_stf_coupler(index_added_heat),           &
             T_prog(index_added_heat)%stf(:,:)*T_prog(index_added_heat)%conversion)
    endif
  endif
  
  ! total surface heat flux (Watts) passed through coupler 
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_coupler(index_temp), T_prog(index_temp)%stf(:,:), 1e-15*T_prog(index_temp)%conversion)
  
  ! total FAFMIP added heat flux (Watts) 
  if ( index_added_heat > 0 ) then
    call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_coupler(index_added_heat), T_prog(index_added_heat)%stf(:,:), 1e-15*T_prog(index_added_heat)%conversion)  
  endif

  ! heat input from liquid river runoff relative to 0 degrees C (W/m2)
  if (id_stf_runoff(index_temp) > 0) then
     call diagnose_2d(Time, Grd, id_stf_runoff(index_temp),                           &
             T_prog(index_temp)%runoff_tracer_flux(:,:)*T_prog(index_temp)%conversion)
  endif
  ! total heat flux from liquid river runoff (Watts), relative to 0C. 
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_runoff(index_temp), T_prog(index_temp)%runoff_tracer_flux(:,:), 1e-15*T_prog(index_temp)%conversion)

  ! heat input from solid calving land ice relative to 0 degrees C (W/m2)
  if (id_stf_calving(index_temp) > 0) then
     call diagnose_2d(Time, Grd, id_stf_calving(index_temp),                           &
             T_prog(index_temp)%calving_tracer_flux(:,:)*T_prog(index_temp)%conversion)
  endif
  ! total heat flux from solid calving land ice (Watts), relative to 0C. 
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_calving(index_temp), T_prog(index_temp)%calving_tracer_flux(:,:), 1e-15*T_prog(index_temp)%conversion)

  ! total heat flux from liquid runoff + solid calving land ice (Watts), relative to 0C. 

  if(id_total_ocean_river_heat > 0) then 
     call diagnose_sum(Time, Grd, Dom, id_total_ocean_river_heat, T_prog(index_temp)%calving_tracer_flux(:,:) + T_prog(index_temp)%runoff_tracer_flux, 1e-15*T_prog(index_temp)%conversion)
  endif

  ! heat input from liquid precip relative to 0 degrees C (W/m2).
  ! note that frozen precip arrives at 0C, so contributes no heat 
  ! relative to 0C.  Assume temp of liquid precip same as tpme
  if (id_stf_prec(index_temp) > 0) then
     call diagnose_2d(Time, Grd, id_stf_prec(index_temp),                                   &
             liquid_precip(:,:)*T_prog(index_temp)%tpme(:,:)*T_prog(index_temp)%conversion)
  endif
  ! total heat flux from liquid precip (Watts)
  if(id_total_ocean_stf_prec(index_temp) > 0) then 
     call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_prec(index_temp), liquid_precip(:,:)*T_prog(index_temp)%tpme(:,:), 1e-15*T_prog(index_temp)%conversion)
  endif

  ! heat sent away from ocean due to water mass leaving ocean
  ! via evaporation, measured relative to 0 degrees C (W/m2).
  ! Assume temp of evaporating water is same as tpme
  if (id_stf_evap(index_temp) > 0) then
     call diagnose_2d(Time, Grd, id_stf_evap(index_temp),                                &
             evaporation(:,:)*T_prog(index_temp)%tpme(:,:)*T_prog(index_temp)%conversion)
  endif
  ! total heat flux from evaporating water carrying heat away from ocean (Watts)
  if(id_total_ocean_stf_evap(index_temp) > 0) then 
     call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_evap(index_temp), evaporation(:,:)*T_prog(index_temp)%tpme(:,:), 1e-15*T_prog(index_temp)%conversion)
  endif

  ! net heat flux from radiation+latent+sensible (as passed through coupler) + mass transport 
  ! note that the addition of frozen_precip is operationally how the model computes the 
  ! heat contribution from mass transport.  however, it is arguably not correct, since 
  ! the frozen precip is typically best approximated to be at 0C, rather than SST. But
  ! we diagnose the contribution in this manner in order to agree with the prognostic model
  ! methods.  
  if(id_net_sfc_heating > 0) then 
      wrk1_2d(:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec
            wrk1_2d(i,j) =   T_prog(index_temp)%conversion*(                                              &
                   T_prog(index_temp)%stf(i,j)                                                            &
                 + T_prog(index_temp)%runoff_tracer_flux(i,j)                                             &
                 + T_prog(index_temp)%calving_tracer_flux(i,j)                                            &
                 + (frozen_precip(i,j)+liquid_precip(i,j)+evaporation(i,j)                                &
#if defined(ACCESS_CM) || defined(ACCESS_OM)
                 + Ice_ocean_boundary%wfimelt(i,j) &
                 + Ice_ocean_boundary%wfiform(i,j) &
                 + Ice_ocean_boundary%licefw(i,j) &
#endif
                   )*T_prog(index_temp)%tpme(i,j))
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_net_sfc_heating, wrk1_2d(:,:))
  endif
  
  ! net_surface_heating*g*alphasfc2/rho0/Cp  
  if(id_net_sfc_workq > 0) then 
      wrk1_2d(:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec
            wrk1_2d(i,j) =   grav/rho0/cp_ocean*alphasfc2(i,j)*                                           &
                   T_prog(index_temp)%conversion*(                                                        &
                   T_prog(index_temp)%stf(i,j)                                                            &
                 + T_prog(index_temp)%runoff_tracer_flux(i,j)                                             &
                 + T_prog(index_temp)%calving_tracer_flux(i,j)                                            &
                 + (frozen_precip(i,j)+liquid_precip(i,j)+evaporation(i,j)                               &
#if defined(ACCESS_CM) || defined(ACCESS_OM)
                 + Ice_ocean_boundary%wfimelt(i,j)                                                       & 
                 + Ice_ocean_boundary%wfiform(i,j)                                                       &
                 + Ice_ocean_boundary%licefw(i,j)                                                        &
#endif
                   )*T_prog(index_temp)%tpme(i,j) )
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_net_sfc_workq, wrk1_2d(:,:))
  endif
  
  ! area integrated total net heat flux 
  if(id_total_net_sfc_heating > 0) then 
      wrk1_2d(:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec
            wrk1_2d(i,j) = T_prog(index_temp)%conversion*(                                               &
                   T_prog(index_temp)%stf(i,j)                                                           &
                 + T_prog(index_temp)%runoff_tracer_flux(i,j)                                            &
                 + T_prog(index_temp)%calving_tracer_flux(i,j)                                           &
                 + (frozen_precip(i,j)+liquid_precip(i,j)+evaporation(i,j)                               &
#if defined(ACCESS_CM) || defined(ACCESS_OM)
                 + Ice_ocean_boundary%wfimelt(i,j)                                                       & 
                 + Ice_ocean_boundary%wfiform(i,j)                                                       &
                 + Ice_ocean_boundary%licefw(i,j)                                                        &
#endif
                   )*T_prog(index_temp)%tpme(i,j) )
         enddo
      enddo
      call diagnose_sum(Time, Grd, Dom, id_total_net_sfc_heating, wrk1_2d, 1e-15)
  endif

  ! shortwave flux (W/m2)
  call diagnose_2d(Time, Grd, id_swflx, swflx(:,:))
  ! total shortwave heat transport (Watts)
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_swflx, swflx, 1e-15)
  ! swflx impacts on water mass transformation in neutral density classes 
  if(id_tform_rho_pbl_sw_on_nrho > 0) then
      wrk1(:,:,:) = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = swflx(i,j)*cp_ocean_r*Dens%drhodT(i,j,k) &
                         *Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_tform_rho_pbl_sw_on_nrho, wrk1)
  endif

  ! visible shortwave flux (W/m2)
  call diagnose_2d(Time, Grd, id_swflx_vis, swflx_vis(:,:))
  ! total visible shortwave (Watts)
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_swflx_vis, swflx_vis, 1e-15)


  ! evaporative heat flux (W/m2) (<0 cools ocean)
  if (id_evap_heat > 0) then
      tmp_flux=0.0
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = -latent_heat_vapor(ii,jj)*Ice_ocean_boundary%q_flux(i,j)
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_evap_heat, tmp_flux(:,:))
  endif
  ! total evaporative heating (Watts) 
  if (id_total_ocean_evap_heat > 0) then
      tmp_flux=0.0
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = -latent_heat_vapor(ii,jj)*Ice_ocean_boundary%q_flux(i,j)
         enddo
      enddo
      call diagnose_sum(Time, Grd, Dom, id_total_ocean_evap_heat, tmp_flux, 1e-15)
  endif

  ! latent heat (liquid-vapor and solid-liquid) 
  ! impacts on water mass transformation in neutral density classes 
  if(id_tform_rho_pbl_lat_on_nrho > 0) then
      wrk1(:,:,:)      = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = latent(i,j)*cp_ocean_r*Dens%drhodT(i,j,k) &
                         *Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_tform_rho_pbl_lat_on_nrho, wrk1)
  endif


  ! longwave heat flux (W/m2)
  call diagnose_2d(Time, Grd, id_lw_heat, longwave(:,:))
  ! total longwave heating (Watts) 
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_lw_heat, longwave, 1e-15)
  ! longwave impacts on water mass transformation in neutral density classes 
  if(id_tform_rho_pbl_lw_on_nrho > 0) then
      wrk1(:,:,:)      = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = longwave(i,j)*cp_ocean_r*Dens%drhodT(i,j,k) &
                         *Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_tform_rho_pbl_lw_on_nrho, wrk1)
  endif

  ! heat flux from melting the frozen precip (W/m2)
  if (id_fprec_melt_heat > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = -Ice_ocean_boundary%fprec(i,j)*latent_heat_fusion(ii,jj)
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_fprec_melt_heat, tmp_flux(:,:))
  endif
  ! total heating from melting the frozen precip (Watts) 
  if (id_total_ocean_fprec_melt_heat > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = -Ice_ocean_boundary%fprec(i,j)*latent_heat_fusion(ii,jj)
         enddo
      enddo
      call diagnose_sum(Time, Grd, Dom, id_total_ocean_fprec_melt_heat, tmp_flux, 1e-15)
  endif

  ! heat flux from the melting of calved land ice (W/m2)
  if (id_calving_melt_heat > 0) then
      wrk1_2d(:,:) = -calving(:,:)*latent_heat_fusion(:,:)
      call diagnose_2d(Time, Grd, id_calving_melt_heat, wrk1_2d(:,:))
  endif
  ! total heating from the melting of calved land ice (Watts)  
  if (id_total_ocean_calving_melt_heat > 0) then
     call diagnose_sum(Time, Grd, Dom, id_total_ocean_calving_melt_heat, -calving(:,:)*latent_heat_fusion(:,:), 1e-15)
  endif

  ! sensible heat flux (W/m2)
  call diagnose_2d(Time, Grd, id_sens_heat, sensible(:,:))
  ! total sensible heat transport (Watts) 
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_sens_heat, sensible, 1e-15)
  ! sensible heat impacts on water mass transformation in neutral density classes 
  if(id_tform_rho_pbl_sens_on_nrho > 0) then
      wrk1(:,:,:)      = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = sensible(i,j)*cp_ocean_r*Dens%drhodT(i,j,k) &
                         *Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_tform_rho_pbl_sens_on_nrho, wrk1)
  endif

  ! contribution from total pbl heat fluxes on 
  ! water mass transformation in neutral density classes 
  if(id_tform_rho_pbl_heat_on_nrho > 0) then
      wrk1(:,:,:)      = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = (longwave(i,j)+sensible(i,j)+latent(i,j)+swflx(i,j)) &
                          *cp_ocean_r*Dens%drhodT(i,j,k)                       &
                          *Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_tform_rho_pbl_heat_on_nrho, wrk1)
  endif

  ! contribution from total pbl heat and salt fluxes on 
  ! water mass transformation in neutral density classes 
  if(id_tform_rho_pbl_flux_on_nrho > 0) then
      wrk1(:,:,:)      = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) =  (T_prog(index_salt)%stf(i,j)*Dens%drhodS(i,j,k)      &
                          +(longwave(i,j)+sensible(i,j)+latent(i,j)+swflx(i,j)) &
                           *cp_ocean_r*Dens%drhodT(i,j,k))                      & 
                          *Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_tform_rho_pbl_flux_on_nrho, wrk1)
  endif

#if defined(ACCESS_CM) || defined(ACCESS_OM)
   ! Heat into ocean due to land ice discharge-melt (>0 heats ocean)
   if (id_liceht > 0) then
       do j=jsc_bnd,jec_bnd
          do i=isc_bnd,iec_bnd
             ii=i+i_shift
             jj=j+j_shift
             tmp_flux(ii,jj) = Ice_ocean_boundary%liceht(i,j)
          enddo
       enddo
       used = send_data(id_liceht, tmp_flux(:,:),        &
              Time%model_time, rmask=Grd%tmask(:,:,1),  &
              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      call diagnose_2d(Time, Grd, id_liceht, tmp_flux(:,:))
   endif
   if (id_total_ocean_liceht > 0) then
       do j=jsc_bnd,jec_bnd
          do i=isc_bnd,iec_bnd
             ii=i+i_shift
             jj=j+j_shift
             tmp_flux(ii,jj) = Ice_ocean_boundary%liceht(i,j)
          enddo
       enddo
       call diagnose_sum(Time, Grd, Dom, id_total_ocean_liceht, tmp_flux, 1e-15)
    endif

   ! Heat into ocean due to melting ice (>0 heats ocean)
   if (id_mh_flux > 0) then
       do j=jsc_bnd,jec_bnd
          do i=isc_bnd,iec_bnd
             ii=i+i_shift
             jj=j+j_shift
             tmp_flux(ii,jj) = Ice_ocean_boundary%mh_flux(i,j)
          enddo
       enddo
      call diagnose_2d(Time, Grd, id_mh_flux, tmp_flux(:,:))
   endif
   if (id_total_ocean_mh_flux > 0) then
       do j=jsc_bnd,jec_bnd
          do i=isc_bnd,iec_bnd
             ii=i+i_shift
             jj=j+j_shift
             tmp_flux(ii,jj) = Ice_ocean_boundary%mh_flux(i,j)
          enddo
       enddo
       call diagnose_sum(Time, Grd, Dom, id_total_ocean_mh_flux, tmp_flux, 1e-15)
    endif
#endif

  !--------salt related diagnostics ------------------------------------
  !
  ! salt flux (kg/(m2*sec)) passed through the coupler 
  if (id_stf_coupler(index_salt) > 0) then
     call diagnose_2d(Time, Grd, id_stf_coupler(index_salt),            &
             T_prog(index_salt)%stf(:,:)*T_prog(index_salt)%conversion)
  endif

  ! total salt flux (kg/sec) passed through coupler 
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_coupler(index_salt), T_prog(index_salt)%stf(:,:), 1e-15*T_prog(index_salt)%conversion)

  ! salt input from liquid river runoff (kg/(m2*sec))
  if (id_stf_runoff(index_salt) > 0) then
     call diagnose_2d(Time, Grd, id_stf_runoff(index_salt),                           &
             T_prog(index_salt)%runoff_tracer_flux(:,:)*T_prog(index_salt)%conversion)
  endif

  ! salt input from calving land ice (kg/(m2*sec))
  if (id_stf_runoff(index_salt) > 0) then
     call diagnose_2d(Time, Grd, id_stf_runoff(index_salt),                            &
             T_prog(index_salt)%calving_tracer_flux(:,:)*T_prog(index_salt)%conversion)
  endif

  ! total salt flux from liquid river runoff (kg/sec)
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_runoff(index_salt), T_prog(index_salt)%runoff_tracer_flux(:,:), 1e-15*T_prog(index_salt)%conversion)

  ! total salt flux from calving land ice (kg/sec)
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_stf_calving(index_salt), T_prog(index_salt)%calving_Tracer_flux(:,:), 1e-15*T_prog(index_salt)%conversion)

  ! salt input from ice (kg/(m2*sec))
  if (id_salt_flux_ice > 0) then
     call diagnose_2d(Time, Grd, id_salt_flux_ice, melt(:,:)*ice_salt_concentration)
  endif
  ! total salt flux from ice (kg/sec)
  call diagnose_sum(Time, Grd, Dom, id_total_salt_flux_ice, melt, 1e-15*ice_salt_concentration)

  ! salt flux impacts on water mass transformation in neutral density classes 
  if(id_tform_rho_pbl_salt_on_nrho > 0) then
      wrk1(:,:,:)      = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = T_prog(index_salt)%stf(i,j)*Dens%drhodS(i,j,k) &
                          *Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_tform_rho_pbl_salt_on_nrho, wrk1)
  endif


  !--------mass flux related diagnostics-----------------------------
  !
  ! total mass flux per area from pme and river (kg/(m2*sec))  
  if (id_pme_river > 0) then
     call diagnose_2d(Time, Grd, id_pme_river, pme(:,:) + river(:,:))
  endif
  
  !ape work from E-P+R: pme_river*g*betasfc2*So/rho0  
  if(id_net_sfc_workemp > 0) then
     call diagnose_2d(Time, Grd, id_net_sfc_workemp,        &
        grav*salinity_ref/rho0*betasfc2(:,:)*(pme(:,:) + river(:,:)) )
  endif
  
  ! total mass flux from pme+river (kg/sec)
  if(id_total_ocean_pme_river > 0) then 
     call diagnose_sum(Time, Grd, Dom,id_total_ocean_pme_river, pme(:,:) + river(:,:), 1e-15)
  endif
  ! bin pme+river into neutral density classes 
  if(id_mass_pmepr_on_nrho > 0) then
      wrk1(:,:,:)      = 0.0
      k=1 
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%dat(i,j)*(pme(i,j)+river(i,j))
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_mass_pmepr_on_nrho, wrk1)
  endif

  ! mass flux per area from pme_sbc (kg/(m2*sec))  
  call diagnose_2d(Time, Grd, id_pme_sbc, pme(:,:))

  ! total mass flux from pme_sbc (kg/sec)
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_pme_sbc, pme, 1e-15)

  ! mass flux per area from ice melt (kg/(m2*sec))  
  call diagnose_2d(Time, Grd, id_melt, melt(:,:))

  ! total mass flux from ice melt (kg/sec)
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_melt, melt, 1e-15)

  ! bin ice melt/form into neutral density classes 
  if(id_mass_melt_on_nrho > 0) then
      wrk1(:,:,:)      = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%dat(i,j)*melt(i,j)
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_mass_melt_on_nrho, wrk1)
  endif


  ! evaporative mass flux (kg/(m2*sec))
  ! evaporation > 0 means liquid water enters ocean. 
  call diagnose_2d(Time, Grd, id_evap, evaporation(:,:))
  ! total mass transport from evap (kg/sec)
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_evap, evaporation, 1e-15)
  ! bin evap/condense mass transport into neutral density classes 
  if(id_mass_evap_on_nrho > 0) then
      wrk1(:,:,:)      = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%dat(i,j)*evaporation(i,j)
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_mass_evap_on_nrho, wrk1)
  endif

  ! frozen precip (kg/(m2*sec))
  call diagnose_2d(Time, Grd, id_fprec, frozen_precip(:,:))
  ! total frozen precip (kg/sec)
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_fprec, frozen_precip, 1e-15)


  ! liquid precip (kg/(m2*sec))
  call diagnose_2d(Time, Grd, id_lprec, liquid_precip(:,:))
  ! total liquid precip (kg/sec)
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_lprec, liquid_precip, 1e-15)
  ! bin precip (liquid and frozen) mass transport into neutral density classes 
  if(id_mass_precip_on_nrho > 0) then
      wrk1(:,:,:)      = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%dat(i,j)*(liquid_precip(i,j) + frozen_precip(i,j))
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_mass_precip_on_nrho, wrk1)
  endif

#if defined(ACCESS_CM) || defined(ACCESS_OM)
  ! waterflux associated with ice melt (kg/(m2*sec))
  if (id_wfimelt > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift
            tmp_flux(ii,jj) = Ice_ocean_boundary%wfimelt(i,j)
         enddo
      enddo
      used = send_data(id_wfimelt, tmp_flux(:,:),        &
             Time%model_time, rmask=Grd%tmask(:,:,1),  &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total water into ocean associated ice melt (kg/sec)
  if (id_total_ocean_wfimelt > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift
            tmp_flux(ii,jj) = Ice_ocean_boundary%wfimelt(i,j)
         enddo
      enddo
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*tmp_flux(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_wfimelt, total_stuff*1e-15, Time%model_time)
  endif
  ! waterflux associated with ice form (kg/(m2*sec))
  if (id_wfiform > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift
            tmp_flux(ii,jj) = Ice_ocean_boundary%wfiform(i,j)
         enddo
      enddo
      used = send_data(id_wfiform, tmp_flux(:,:),        &
             Time%model_time, rmask=Grd%tmask(:,:,1),  &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total water outof ocean associated ice melt (kg/sec)
  if (id_total_ocean_wfiform > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift
            tmp_flux(ii,jj) = Ice_ocean_boundary%wfiform(i,j)
         enddo
      enddo
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*tmp_flux(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_wfiform, total_stuff*1e-15, Time%model_time)
  endif
  if (id_aice > 0) used = send_data(id_aice, aice(:,:),    &
                 Time%model_time, rmask=Grd%tmask(:,:,1),  &
                 is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
  if (id_iof_nit > 0) used = send_data(id_iof_nit, iof_nit(:,:),    &
                 Time%model_time, rmask=Grd%tmask(:,:,1),  &
                 is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_iof_alg > 0) used = send_data(id_iof_alg, iof_alg(:,:),    &
                 Time%model_time, rmask=Grd%tmask(:,:,1),  &
                 is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
#endif
  if (id_wnd > 0) used = send_data(id_wnd, Velocity%u10(:,:),    &
                 Time%model_time, rmask=Grd%tmask(:,:,1),  &
                 is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
   ! waterflux associated with land ice melt into ocean (kg/(m2*sec))
   if (id_licefw > 0) then
       do j=jsc_bnd,jec_bnd
          do i=isc_bnd,iec_bnd
             ii=i+i_shift
             jj=j+j_shift
             tmp_flux(ii,jj) = Ice_ocean_boundary%licefw(i,j)
          enddo
       enddo
       used = send_data(id_licefw, tmp_flux(:,:),        &
              Time%model_time, rmask=Grd%tmask(:,:,1),  &
              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
   endif
#if defined(ACCESS_CM)
  if (id_atm_co2 > 0) used = send_data(id_atm_co2, atm_co2(:,:),    &
                 Time%model_time, rmask=Grd%tmask(:,:,1),  &
                 is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
#endif
#endif


  ! output saline contraction coeff (1/rho drho/dS) at the ocean surface (1/psu)
  call diagnose_2d(Time, Grd, id_betasfc, betasfc(:,:))
  ! output thermal expansion coeff (-1/rho drho/dT) at the ocean surface (1/(deg C))
  call diagnose_2d(Time, Grd, id_alphasfc, alphasfc(:,:))
  ! output saline contraction coeff (1/potrho dpotrho/dS) at the ocean surface (1/psu)
  call diagnose_2d(Time, Grd, id_betasfc2, betasfc2(:,:))
  ! output thermal expansion coeff (-1/potrho dpotrho/dT) at the ocean surface (1/(deg C))
  call diagnose_2d(Time, Grd, id_alphasfc2, alphasfc2(:,:))
  ! river (mass flux of land water (liquid+solid) ) entering ocean (kg/m^3)*(m/s)
  call diagnose_2d(Time, Grd, id_river, river(:,:))
  ! river (mass flux of land water (liquid+solid) ) entering ocean (kg/m^3)*(m/s)
  call diagnose_2d(Time, Grd, id_river, river(:,:))
  ! global sum of river input (kg/sec)
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_river, river, 1e-15)
  ! bin river (liquid and frozen) runoff into neutral density classes 
  if(id_mass_river_on_nrho > 0) then
      wrk1(:,:,:)      = 0.0
      k=1
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Grd%dat(i,j)*river(i,j)
         enddo
      enddo
      call diagnose_3d_rho(Time, Dens, id_mass_river_on_nrho, wrk1)
  endif

  ! calving land ice (kg/(m2*sec)) entering the ocean 
  call diagnose_2d(Time, Grd, id_calving, calving(:,:))
  ! total mass of calving (kg/sec)
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_calving, calving, 1e-15)

  ! liquid river runoff entering the ocean (kg/m^3)*(m/s)
  call diagnose_2d(Time, Grd, id_runoff, runoff(:,:))
  ! total liquid river runoff (kg/sec)
  call diagnose_sum(Time, Grd, Dom, id_total_ocean_runoff, runoff, 1e-15)


  !----------------------------------------------------------------------
  ! contributions to sea level forcing 
  if(diagnose_sea_level_forcing) then 

      ! for holding sum of sfc heat fluxes  
      wrk2_2d(:,:) = 0.0

      ! shortwave contribution to sea level 
      do j=jsc,jec
         do i=isc,iec
            tmp_flux(i,j) = cp_ocean_r*swflx(i,j)*Grd%tmask(i,j,1) 
            wrk2_2d(i,j)  = tmp_flux(i,j)
         enddo
      enddo
      if(id_eta_tend_sw > 0 .or. id_eta_tend_sw_glob > 0) then 
          wrk1_2d = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = alphasfc(i,j)*tmp_flux(i,j)*rhosfc_inv(i,j) 
             enddo
          enddo
          call diagnose_2d(Time, Grd, id_eta_tend_sw, wrk1_2d(:,:))
          call diagnose_sum(Time, Grd, Dom, id_eta_tend_sw_glob, wrk1_2d, cellarea_r)
      endif


      ! longwave contribution to sea level 
      do j=jsc,jec
         do i=isc,iec
            tmp_flux(i,j) = cp_ocean_r*longwave(i,j)*Grd%tmask(i,j,1)
            wrk2_2d(i,j)  = wrk2_2d(i,j) + tmp_flux(i,j)
         enddo
      enddo
      if(id_eta_tend_lw > 0 .or. id_eta_tend_lw_glob > 0) then 
          wrk1_2d = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = alphasfc(i,j)*tmp_flux(i,j)*rhosfc_inv(i,j) 
             enddo
          enddo
          call diagnose_2d(Time, Grd, id_eta_tend_lw, wrk1_2d(:,:))
          call diagnose_sum(Time, Grd, Dom, id_eta_tend_lw_glob, wrk1_2d, cellarea_r)
      endif

      ! sensible heat contribution to sea level 
      do j=jsc,jec
         do i=isc,iec
            tmp_flux(i,j) = cp_ocean_r*sensible(i,j)*Grd%tmask(i,j,1)
            wrk2_2d(i,j)  = wrk2_2d(i,j) + tmp_flux(i,j)
         enddo
      enddo
      if(id_eta_tend_sens > 0 .or. id_eta_tend_sens_glob > 0) then 
          wrk1_2d = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = alphasfc(i,j)*tmp_flux(i,j)*rhosfc_inv(i,j) 
             enddo
          enddo
          call diagnose_2d(Time, Grd, id_eta_tend_sens, wrk1_2d(:,:))
          call diagnose_sum(Time, Grd, Dom, id_eta_tend_sens_glob, wrk1_2d, cellarea_r)
      endif

      ! latent heat from vaporization contribution to sea level 
      ! evaporation > 0 means liquid water enters ocean. 
      do j=jsc,jec
         do i=isc,iec
            tmp_flux(i,j) = cp_ocean_r*latent_heat_vapor(i,j)*evaporation(i,j)
            wrk2_2d(i,j)  = wrk2_2d(i,j) + tmp_flux(i,j)*Grd%tmask(i,j,1)
         enddo
      enddo
      if(id_eta_tend_evap_heat > 0 .or. id_eta_tend_evap_heat_glob > 0) then 
          wrk1_2d = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = alphasfc(i,j)*tmp_flux(i,j)*rhosfc_inv(i,j) 
             enddo
          enddo
          call diagnose_2d(Time, Grd, id_eta_tend_evap_heat, wrk1_2d(:,:))
          call diagnose_sum(Time, Grd, Dom, id_eta_tend_evap_heat, wrk1_2d, cellarea_r)
      endif

      ! latent heat from melting frozen precip contribution to sea level.
      ! frozen_precip > 0 means frozen precipitation enters liquid ocean 
      ! and so must be melted by liquid ocean. 
      do j=jsc,jec
         do i=isc,iec
            tmp_flux(i,j) = -cp_ocean_r*frozen_precip(i,j)*latent_heat_fusion(i,j)
            wrk2_2d(i,j)  = wrk2_2d(i,j) + tmp_flux(i,j)*Grd%tmask(i,j,1)
         enddo
      enddo
      if(id_eta_tend_fprec_melt > 0 .or. id_eta_tend_fprec_melt_glob > 0) then 
          wrk1_2d = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = alphasfc(i,j)*tmp_flux(i,j)*rhosfc_inv(i,j) 
             enddo
          enddo
          call diagnose_2d(Time, Grd, id_eta_tend_fprec_melt, wrk1_2d(:,:))
          call diagnose_sum(Time, Grd, Dom, id_eta_tend_fprec_melt_glob, wrk1_2d, cellarea_r)
      endif

      ! latent heat from melting icebergs contribution to sea level  
      ! calving > 0 means solid land icen enters liquid ocean 
      ! and so must be melted by liquid ocean. 
      do j=jsc,jec
         do i=isc,iec
            tmp_flux(i,j) = -cp_ocean_r*calving(i,j)*latent_heat_fusion(i,j)
         enddo
      enddo
      do j=jsc,jec
         do i=isc,iec
            wrk2_2d(i,j) = wrk2_2d(i,j) + tmp_flux(i,j)*Grd%tmask(i,j,1)
         enddo
      enddo
      if(id_eta_tend_iceberg_melt > 0 .or. id_eta_tend_iceberg_melt_glob > 0) then 
          wrk1_2d = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = alphasfc(i,j)*tmp_flux(i,j)*rhosfc_inv(i,j) 
             enddo
          enddo
          call diagnose_2d(Time, Grd, id_eta_tend_iceberg_melt, wrk1_2d(:,:))
          call diagnose_sum(Time, Grd, Dom, id_eta_tend_iceberg_melt_glob, wrk1_2d, cellarea_r)
      endif

      ! sum of sfc heat contributions to sea level from heat passed through coupler
      if(id_eta_tend_heat_coupler > 0 .or. id_eta_tend_heat_coupler_glob > 0) then 
          wrk1_2d = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = alphasfc(i,j)*wrk2_2d(i,j)*rhosfc_inv(i,j) 
             enddo
          enddo
          call diagnose_2d(Time, Grd, id_eta_tend_heat_coupler, wrk1_2d(:,:))
          call diagnose_sum(Time, Grd, Dom, id_eta_tend_heat_coupler_glob, wrk1_2d, cellarea_r)
      endif


      ! salt flux contribution to sea level  
      do j=jsc,jec
         do i=isc,iec
            tmp_flux(i,j) = T_prog(index_salt)%stf(i,j)*Grd%tmask(i,j,1)
         enddo
      enddo
      if(id_eta_tend_salt_coupler > 0 .or. id_eta_tend_salt_coupler_glob > 0) then 
          wrk1_2d = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = -betasfc(i,j)*tmp_flux(i,j)*rhosfc_inv(i,j) 
             enddo
          enddo
          call diagnose_2d(Time, Grd, id_eta_tend_salt_coupler, wrk1_2d(:,:))
          call diagnose_sum(Time, Grd, Dom, id_eta_tend_salt_coupler_glob, wrk1_2d, cellarea_r)
      endif



      ! for holding sum of water fluxes   
      wrk3_2d(:,:) = 0.0 

      ! evaporative water flux contribution to sea level 
      do j=jsc,jec
         do i=isc,iec
            tmp_flux(i,j) = evaporation(i,j)*Grd%tmask(i,j,1)
            wrk3_2d(i,j)  = tmp_flux(i,j)
         enddo
      enddo
      if(id_eta_tend_evap > 0 .or. id_eta_tend_evap_glob > 0) then 
          wrk1_2d = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = tmp_flux(i,j)*rhosfc_inv(i,j) 
             enddo
          enddo
          call diagnose_2d(Time, Grd, id_eta_tend_evap, wrk1_2d(:,:))
          call diagnose_sum(Time, Grd, Dom, id_eta_tend_evap_glob, wrk1_2d, cellarea_r)
      endif


      ! liquid precip contribution to sea level 
      do j=jsc,jec
         do i=isc,iec
            tmp_flux(i,j) = liquid_precip(i,j)*Grd%tmask(i,j,1)
            wrk3_2d(i,j)  = wrk3_2d(i,j) + tmp_flux(i,j)
         enddo
      enddo
      if(id_eta_tend_lprec > 0 .or. id_eta_tend_lprec_glob > 0) then 
          wrk1_2d = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = tmp_flux(i,j)*rhosfc_inv(i,j) 
             enddo
          enddo
          call diagnose_2d(Time, Grd, id_eta_tend_lprec, wrk1_2d(:,:))
          call diagnose_sum(Time, Grd, Dom, id_eta_tend_lprec_glob, wrk1_2d, cellarea_r)
      endif


      ! frozen precip contribution to sea level 
      do j=jsc,jec
         do i=isc,iec
            tmp_flux(i,j) = frozen_precip(i,j)*Grd%tmask(i,j,1)
            wrk3_2d(i,j)  = wrk3_2d(i,j) + tmp_flux(i,j)
         enddo
      enddo
      if(id_eta_tend_fprec > 0 .or. id_eta_tend_fprec_glob > 0) then 
          wrk1_2d = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = tmp_flux(i,j)*rhosfc_inv(i,j) 
             enddo
          enddo
          call diagnose_2d(Time, Grd, id_eta_tend_fprec, wrk1_2d(:,:))
          call diagnose_sum(Time, Grd, Dom, id_eta_tend_fprec_glob, wrk1_2d, cellarea_r)
      endif


      ! liquid runoff contribution to sea level 
      do j=jsc,jec
         do i=isc,iec
            tmp_flux(i,j) = runoff(i,j)*Grd%tmask(i,j,1)
            wrk3_2d(i,j)  = wrk3_2d(i,j) + tmp_flux(i,j)
         enddo
      enddo
      if(id_eta_tend_runoff > 0 .or. id_eta_tend_runoff_glob > 0) then 
          wrk1_2d = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = tmp_flux(i,j)*rhosfc_inv(i,j) 
             enddo
          enddo
          call diagnose_2d(Time, Grd, id_eta_tend_runoff, wrk1_2d(:,:))
          call diagnose_sum(Time, Grd, Dom, id_eta_tend_runoff_glob, wrk1_2d, cellarea_r)
      endif

      ! iceberg contribution to sea level 
      do j=jsc,jec
         do i=isc,iec
            tmp_flux(i,j) = calving(i,j)*Grd%tmask(i,j,1)
            wrk3_2d(i,j)  = wrk3_2d(i,j) + tmp_flux(i,j)
         enddo
      enddo
      if(id_eta_tend_iceberg > 0 .or. id_eta_tend_iceberg_glob > 0) then 
          wrk1_2d = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = tmp_flux(i,j)*rhosfc_inv(i,j) 
             enddo
          enddo
          call diagnose_2d(Time, Grd, id_eta_tend_iceberg, wrk1_2d(:,:))
          call diagnose_sum(Time, Grd, Dom, id_eta_tend_iceberg_glob, wrk1_2d, cellarea_r)
      endif

      ! sum of the water fluxes passed through coupler contribution to sea level 
      if(id_eta_tend_water_coupler > 0 .or. id_eta_tend_water_coupler_glob > 0) then
          do j=jsc,jec
             do i=isc,iec
                wrk3_2d(i,j) = wrk3_2d(i,j)*rhosfc_inv(i,j) 
             enddo
          enddo
          call diagnose_2d(Time, Grd, id_eta_tend_water_coupler, wrk3_2d(:,:))
          call diagnose_sum(Time, Grd, Dom, id_eta_tend_water_coupler_glob, wrk1_2d, cellarea_r)
      endif


  endif   ! endif for diagnose_sea_level_forcing
  !----------------------------------------------------------------------


end subroutine ocean_sbc_diag
! </SUBROUTINE> NAME="ocean_sbc_diag"


!#######################################################################
! <SUBROUTINE> NAME="compute_latent_heat_fusion"
!
! <DESCRIPTION>
! TEOS-10 expression for latent heat of fusion at the sea surface (p=0dbar)
!
! The following is from the from the matlab routine due to 
! McDougall and Barker with respect to the full pressure dependent 
! formulation.
!
! Note that the computed latent heat of fusion from this function has 
! errors that range between -0.4 and 0.3 J kg^-1, when compared with the 
! latent heats of melting derived from the Gibbs functions of ice and of 
! seawater (using the SIA code of TEOS-10).  However, the underlying data to
! the Gibbs function contains uncertainities of 200 J kg^-1 (IOC et al., 2010).  
!
! The reference for this routine is 
!
!  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
!  seawater - 2010: Calculation and use of thermodynamic properties.  
!  Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
!  UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
!  See section 3.34 of the TEOS-10 Manual.
! </DESCRIPTION>

subroutine compute_latent_heat_fusion(salinity)

real, dimension(isd:,jsd:), intent(in) :: salinity

real, parameter  :: c0  =  3.334265169240710d5
real, parameter  :: c1  = -2.789444646733159d0
real, parameter  :: c3  = -4.984585692734338d3
real, parameter  :: c6  =  1.195857305019339d3
real, parameter  :: c10 = -5.792068522727968d2
real, parameter  :: c15 =  6.836527214265952d2
real, parameter  :: c21 = -2.371103254714944d2

real    :: x 
integer :: i,j

if (constant_hlf) return 

do j=jsd,jed
   do i=isd,ied
      x = salinity(i,j)/(40.0*(35.16504/35.0))
      latent_heat_fusion(i,j) =  c0 + x*(c1 +  x*(c3 + x*(c6 + &
                                 x*(c10  +  x*(c15 + c21*x))))) 
   enddo
enddo

return

end subroutine compute_latent_heat_fusion
! </SUBROUTINE> NAME="compute_latent_heat_fusion"



!#######################################################################
! <SUBROUTINE> NAME="compute_latent_heat_vapor"
!
! <DESCRIPTION>
! TEOS-10 expression for latent heat of vaporization at the sea surface 
! (p=0dbar).

! The following is from the From the original matlab routine due to
! Barker, McDougall and Feistel
!
!  Calculates latent heat, or enthalpy, of evaporation at p = 0 (the 
!  surface).  It is defined as a function of Absolute Salinity, SA, and
!  Conservative Temperature, CT, and is valid in the ranges 
!  0 < SA < 42 g/kg and 0 < CT < 40 deg C.  The errors range between 
!  -0.4 and 0.6 J/kg
!
!  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
!   seawater - 2010: Calculation and use of thermodynamic properties.  
!   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
!   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
!   See section 3.39 of the TEOS-10 Manual.
!
! </DESCRIPTION>
subroutine compute_latent_heat_vapor(salinity, theta)

real, dimension(isd:,jsd:), intent(in) :: salinity
real, dimension(isd:,jsd:), intent(in) :: theta
 
real, parameter :: c0 =   2.499065844825125e6;
real, parameter :: c1 =  -1.544590633515099e-1;
real, parameter :: c2 =  -9.096800915831875e4;
real, parameter :: c3 =   1.665513670736000e2;
real, parameter :: c4 =   4.589984751248335e1;
real, parameter :: c5 =   1.894281502222415e1;
real, parameter :: c6 =   1.192559661490269e3;
real, parameter :: c7 =  -6.631757848479068e3;
real, parameter :: c8 =  -1.104989199195898e2;
real, parameter :: c9 =  -1.207006482532330e3;
real, parameter :: c10 = -3.148710097513822e3;
real, parameter :: c11 =  7.437431482069087e2;
real, parameter :: c12 =  2.519335841663499e3;
real, parameter :: c13 =  1.186568375570869e1;
real, parameter :: c14 =  5.731307337366114e2;
real, parameter :: c15 =  1.213387273240204e3;
real, parameter :: c16 =  1.062383995581363e3;
real, parameter :: c17 = -6.399956483223386e2;
real, parameter :: c18 = -1.541083032068263e3;
real, parameter :: c19 =  8.460780175632090e1;
real, parameter :: c20 = -3.233571307223379e2;
real, parameter :: c21 = -2.031538422351553e2;
real, parameter :: c22 =  4.351585544019463e1;
real, parameter :: c23 = -8.062279018001309e2;
real, parameter :: c24 =  7.510134932437941e2;
real, parameter :: c25 =  1.797443329095446e2;
real, parameter :: c26 = -2.389853928747630e1;
real, parameter :: c27 =  1.021046205356775e2;

real    :: x,y
integer :: i,j

if (constant_hlv) return 

do j=jsd,jed
   do i=isd,ied

      x = salinity(i,j)/(40*(35.16504/35))
      y=theta(i,j)/40.0

      latent_heat_vapor(i,j) = c0 + x*(c1 + c4*y + x*(c3               &
           + y*(c7 + c12*y) + x*(c6 + y*(c11 + y*(c17 + c24*y))        &
           + x*(c10 + y*(c16 + c23*y) + x*(c15 + c22*y + c21*x)))))    &
           + y*(c2 + y*(c5 + c8*x + y*(c9 + x*(c13 + c18*x)            &
           + y*(c14 + x*(c19 + c25*x) + y*(c20 + c26*x + c27*y)))))

   enddo
enddo

return

end subroutine compute_latent_heat_vapor
! </SUBROUTINE> NAME="compute_latent_heat_vapor"


end module ocean_sbc_mod

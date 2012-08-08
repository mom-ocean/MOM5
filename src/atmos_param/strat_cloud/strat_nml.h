! $Id: strat_nml.h,v 19.0 2012/01/06 20:27:25 fms Exp $
! $Name: siena_201207 $

!------------------------------------------------------------------------
!---namelist------

!------------------------------------------------------------------------
!   the namelist variables are stored in a strat_nml_type derived type 
!   variable so that they may be conveniently passed to other modules.
!------------------------------------------------------------------------

! <NAMELIST NAME="strat_cloud_nml">
!  <DATA NAME="do_netcdf_restart" TYPE="logical" DEFAULT=".true.">
!   a netcdf restart file should be written ? 
!  </DATA>
!  <DATA NAME="U00" UNITS="fraction" TYPE="real" DEFAULT="0.80">
!   threshold rel hum for cloud formation by large-scale condensation.  
!  </DATA>
!  <DATA NAME="u00_profile" TYPE="logical"  DEFAULT=".false.">
!   the low-level u00 ECMWF profile should be applied? 
!  </DATA>
!  <DATA NAME="rthresh" UNITS="microns" TYPE="real" DEFAULT="10.">
!   liquid cloud drop radius threshold for autoconversion. 
!  </DATA>
!  <DATA NAME="use_kk_auto" TYPE="logical"  DEFAULT=".false.">
!   the Khairoutdinov and Kogan (2000) autoconversion should be used? 
!  </DATA>
!  <DATA NAME="U_evap" UNITS="fraction" TYPE="real" DEFAULT="1.0">
!   critical relative humidity above which rain does not evaporate. 
!  </DATA>
!  <DATA NAME="eros_scale" UNITS="1/sec" TYPE="real" DEFAULT="1.E-06">
!   normal erosion rate constant for cloud destruction (default = 1.E-06) 
!  </DATA>
!  <DATA NAME="eros_choice" TYPE="logical"  DEFAULT=".false.">
!   enhanced erosion should occur in turbulent or convective conditions ?
!  </DATA>
!  <DATA NAME="eros_scale_c" UNITS="1/sec" TYPE="real"  DEFAULT="8.E-06">
!   erosion rate constant for cloud destruction for convective conditions. 
!  </DATA>
!  <DATA NAME="eros_scale_t" UNITS="1/sec" TYPE="real" DEFAULT="5.E-05">
!   erosion rate constant for cloud destruction for turbulent conditions. 
!  </DATA>
!  <DATA NAME="mc_thresh" UNITS="kg/m2/sec" TYPE="real" DEFAULT="0.001">
!   convective mass-flux threshold for enhanced erosion to turn on.  
!  </DATA>
!  <DATA NAME="diff_thresh" UNITS="m2/s" TYPE="real" DEFAULT="1.0">
!   diffusion coefficient threshold for enhanced erosion to turn on. 
!  </DATA>
!  <DATA NAME="super_choice" TYPE="logical" DEFAULT=".false.">
!   excess vapor in supersaturated conditions should be put into 
!   cloud water (true) or precipitation fluxes (false) ? 
!  </DATA>
!  <DATA NAME="tracer_advec" TYPE="logical" DEFAULT=".false.">
!   cloud liquid, ice and fraction are advected by the grid resolved motion?
!  </DATA>
!  <DATA NAME="qmin" UNITS="kg h2o/kg air" TYPE="real"  DEFAULT="1.E-10">
!   minimum permissible value of cloud liquid, cloud ice, saturated volume 
!   fraction, or rain and snow areas.
!   NOTE: qmin should be chosen such that the range of {qmin, max(qa,ql,qi)}
!   is resolved by the precision of the numbers used. 
!  </DATA>
!  <DATA NAME="Dmin" UNITS="dimensionless" TYPE="real"  DEFAULT="1.0e-08">
!   minimum permissible dissipation in analytic integration of qa, ql, qi 
!   equations. This constant only affects the method by which the 
!   prognostic equations are integrated.
!   NOTE: Dmin will be MACHINE DEPENDENT and occur when
!   a) 1. -exp(-Dmin) = 0. instead of Dmin in the limit of very small Dmin.
!   AND
!   b) 1. - exp(-D) < D for all D > Dmin
!  </DATA>
!  <DATA NAME="efact" UNITS="dimensionless" TYPE="real" DEFAULT="0.0">
!   factor used in formula to enhance cloud erosion as a function of height
!  </DATA>
!  <DATA NAME="vfact" UNITS="" TYPE="real" DEFAULT="1.0">
!   factor used to enhance ice fall velocity.
!  </DATA>
!  <DATA NAME="cfact" UNITS="" TYPE="real" DEFAULT="1.0">
!   factor in bergeron process calculation when the effect of dust
!   is not considered; values > 1 enhance the conversion to ice. 
!  </DATA>
!  <DATA NAME="do_old_snowmelt" TYPE="logical" DEFAULT=".false.">
!   the old version of snow melting, which has a bug, should be run? 
!  </DATA>
!  <DATA NAME="iwc_crit" UNITS="kg(ice)/m**3" TYPE="real"  DEFAULT="0.0">
!   critical ice-water content below which to apply alternate fall 
!   speed formula  
!  </DATA>
!  <DATA NAME="vfall_const2" UNITS="" TYPE="real"  DEFAULT="3.29">
!   factor for alternate fall speed formula 
!  </DATA>
!  <DATA NAME="vfall_exp2" UNITS="" TYPE="real" DEFAULT="0.16">
!   exponent for alternate fall speed formula 
!  </DATA>
!  <DATA NAME="num_mass_ratio1" UNITS="" TYPE="real" DEFAULT="1.0">
!                 
!  </DATA>
!  <DATA NAME="num_mass_ratio2" UNITS="" TYPE="real" DEFAULT="1.0">
!                  
!  </DATA>
!  <DATA NAME="microphys_scheme" TYPE="character" DEFAULT="rotstayn_klein">
!   the microphysics scheme being used (currently either 
!   "morrison_gettelman" or "rotstayn_klein")
!  </DATA>
!  <DATA NAME="super_ice_opt" TYPE="integer" DEFAULT="0">
!   flag to indicate how to treat supersaturation; 0 => don't allow.  
!  </DATA>
!  <DATA NAME="do_ice_nucl_wpdf" TYPE="logical" DEFAULT=".false.">
!   should we use activate ice nuclei by assuming a pdf ?
!  </DATA>
!  <DATA NAME="use_online_aerosol" TYPE="logical"  DEFAULT=".false.">
!   the online aerosol fields should be used for nucleation source 
!   (rather than climo values) ?  
!  </DATA>
!  <DATA NAME="use_sub_seasalt" TYPE="logical"  DEFAULT=".true.">
!   only sub-micron seasalt particles should be used for nucleation ?
!  </DATA>
!  <DATA NAME="sea_salt_scale" UNITS="" TYPE="real" DEFAULT="0.1">
!   scaling factor used when offline seasalt aerosol is used for 
!   nucleation
!  </DATA>
!  <DATA NAME="om_to_oc" UNITS="" TYPE="real"  DEFAULT="1.67">
!   scaling factor used to convert offline organic matter to organic 
!   carbon for nucleation    
!  </DATA>
!  <DATA NAME="N_land" UNITS="1/(m*m*m)" TYPE="real" DEFAULT="250.E+06">
!   assumed number of cloud drops per unit volume in liquid clouds 
!   over land when droplet number is not prdicted. 
!  </DATA>
!  <DATA NAME="N_ocean" UNITS="1/(m*m*m)" TYPE="real" DEFAULT="100.E+06">
!   assumed number of cloud drops per unit volume in liquid clouds 
!   over ocean when droplet number is not predicted.
!  </DATA>
!  <DATA NAME="var_limit" UNITS=" (m**2)/(s**2)" TYPE="real" DEFAULT="0.0">
!   minimum value of the variance in the vertical velocity pdf 
!  </DATA>
!  <DATA NAME="do_liq_num" TYPE="logical"  DEFAULT=".false.">
!   the prognostic droplet number option is activated ?
!  </DATA>
!  <DATA NAME="do_dust_berg" TYPE="logical" DEFAULT=".false.">
!   sub-micron dust particles are used as ice nuclei for bergeron 
!   process evaluation?
!  </DATA>
!  <DATA NAME="N_min" TYPE="real" DEFAULT="1.0E6">
!   minimum number of droplets allowed in a grid box when predicted
!   droplet number code is activated
!  </DATA>
!  <DATA NAME="do_pdf_clouds" TYPE="logical" DEFAULT=".false.">
!   the statistical cloud scheme should be run?
!  </DATA>
!  <DATA NAME="betaP" UNITS="" TYPE="integer" DEFAULT="5">
!   p-parameter to the beta distribution - used when do_pdf_clouds is true 
!  </DATA>
!  <DATA NAME="qthalfwidth" UNITS="" TYPE="real" DEFAULT="0.1">
!   half-width to the qt PDF - used when do_pdf_clouds is true and 
!   diagnostic variance
!   The fraction of qtbar (mean total water in the grid box) that the 
!   maximum and minimum of the distribution differ from qtbar. That is, 
!   total water at the sub-grid scale may take on values anywhere between 
!   (1.-qthalfwidth)*qtbar and (1.+qthalfwidth)*qtbar
!  </DATA>
!  <DATA NAME="nsublevels" UNITS="" TYPE="integer" DEFAULT="1">
!   number of sublevels to vertical sub-grid cloud structure - used when   
!   do_pdf_cloud is true 
!  </DATA>
!  <DATA NAME="kmap" UNITS="" TYPE="integer" DIM="" DEFAULT="1">
!   PPM partial remap integer - used when do_pdf_cloud is true and if 
!   vertical subgrid structure is used
!  </DATA>
!  <DATA NAME="kord" UNITS="" TYPE="integer"  DEFAULT="7">
!   PPM method number - used when do_pdf_cloud is true and if vertical 
!   subgrid structure is used 
!  </DATA>
!  <DATA NAME="pdf_org" TYPE="logical" DEFAULT=".true.">
!   define pdf based on input water fields (rather than updated ones) ?
!  </DATA>
!  <DATA NAME="num_strat_pts" TYPE="integer"DEFAULT="0">
!   number of grid points where instantaneous output will be saved to 
!   file strat.data.
!   num_strat_pts must be <= max_strat_pts.
!   NOTE: this option is currently not supported.
!  </DATA>
!  <DATA NAME="strat_pts" TYPE="integer" DEFAULT="0,0">
!   "num_strat_pts" pairs of grid indices, i.e., the global indices for i,j.
!   NOTE: this option is currently not supported.
!  </DATA>
!  <DATA NAME="debugo" TYPE="logical" DEFAULT=".false.">
!   should we activated debug diagnostics ?
!  </DATA>
!  <DATA NAME="isamp" TYPE="integer" DEFAULT="1">
!   i coordinate of grid point for which debug data is to be written 
!  </DATA>
!  <DATA NAME="jsamp" TYPE="integer" DEFAULT="1">
!   j coordinate of grid point for which debug data is to be written 
!  </DATA>
!  <DATA NAME="ksamp" TYPE="integer" DEFAULT="1">
!   k coordinate of grid point for which debug data is to be written 
!  </DATA>
! </NAMELIST>

INTEGER, PARAMETER  :: max_strat_pts = 5

  logical           :: do_netcdf_restart = .true.
  real              :: U00            =  0.80
  logical           :: u00_profile    =  .false.
  real              :: rthresh        =  10.
  logical           :: use_kk_auto    =  .false.
  real              :: U_evap         =  1.0
  real              :: eros_scale     =  1.E-06
  logical           :: eros_choice    =  .false.
  real              :: eros_scale_c   =  8.E-06
  real              :: eros_scale_t   =  5.E-05
  real              :: mc_thresh      =  0.001
  real              :: diff_thresh    =  1.0
  logical           :: super_choice   =  .false.
  logical           :: tracer_advec   =  .false.
  real              :: qmin           =  1.E-10
  real              :: Dmin           =  1.E-08
  real              :: efact          = 0.0
  real              :: vfact          = 1.0
  real              :: cfact          = 1.0
  logical           :: do_old_snowmelt= .false.
  real              :: iwc_crit       = 0.
  real              :: vfall_const2   = 3.29
  real              :: vfall_exp2     = 0.16
  real              :: num_mass_ratio1= 1.
  real              :: num_mass_ratio2= 1.
  character(len=64) :: microphys_scheme = 'rotstayn_klein'
  integer           :: super_ice_opt = 0
  logical           :: do_ice_nucl_wpdf = .false.

  logical           :: use_online_aerosol = .false.
  logical           :: use_sub_seasalt = .true.
  real              :: sea_salt_scale =  0.1
  real              :: om_to_oc       =  1.67
  real              :: N_land         =  250.E+06
  real              :: N_ocean        =  100.E+06
  real              :: var_limit      = 0.0
  logical           :: do_liq_num   = .false.
  logical           :: do_dust_berg   = .false.
  real              :: N_min          = 1.E6

  logical           :: do_pdf_clouds  = .false.
  integer           :: betaP          = 5
  real              :: qthalfwidth    = 0.1
  integer           :: nsublevels     = 1
  integer           :: kmap           = 1
  integer           :: kord           = 7
  logical           :: pdf_org  = .true.

  integer                             :: num_strat_pts  =  0
  integer,dimension(2,max_strat_pts)  :: strat_pts = 0
  logical                             :: debugo = .false.
  integer                             :: isamp = 1 
  integer                             :: jsamp = 1 
  integer                             :: ksamp = 1 


namelist / strat_cloud_nml /   &
       do_netcdf_restart, U00, u00_profile, rthresh, use_kk_auto, &
       U_evap, eros_scale, eros_choice, eros_scale_c, eros_scale_t, &
       mc_thresh, diff_thresh, super_choice, tracer_advec, qmin, Dmin, &
       efact, vfact, cfact, do_old_snowmelt, iwc_crit, vfall_const2,  &
       vfall_exp2, num_mass_ratio1, num_mass_ratio2,   &
       microphys_scheme, super_ice_opt, do_ice_nucl_wpdf,  &

       use_online_aerosol, use_sub_seasalt, sea_salt_scale, om_to_oc, &
       N_land, N_ocean, var_limit, do_liq_num, do_dust_berg, N_min,  & 

       do_pdf_clouds, betaP, qthalfwidth, nsublevels, kmap, kord, pdf_org, &

       num_strat_pts, strat_pts, debugo, isamp, jsamp, ksamp



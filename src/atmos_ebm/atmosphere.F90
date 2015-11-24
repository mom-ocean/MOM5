module atmosphere_mod
! 
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> B.L. Samuels
!</CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  Zhi Liang
! </REVIEWER>

!<OVERVIEW>
!         interface for spherical grid dynamical core and physics
!                   for the energy balance model (EBM)
!</OVERVIEW>
!
!<DESCRIPTION>
!    The Energy Balance Model is a simple spectral atmospheric model
!     using diffusion and radiative balance. The EBM solves two
!     prognostic equations for atmospheric temperature and specific
!     humidity 
!</DESCRIPTION>
!
!<INFO>
!     The Energy Balance Model was originally developed at GFDL
!     by I. Held and extensively modified by J. Russell to be
!     MOM4 compatible. Zhi Liang made the appropriate modifications
!     for the code to be fms_io compatible and meet RTS standards.
!      No current reference of this code exists.
!</INFO>
!<REFERENCE>
!    Zhang,Rong and Geoffrey Vallis, " The Role of the Deep Western Boundary Current in the
!     Separation of the Gulf Stream", Journal of Physical Oceanography,submitted 2004
!     see Appendix A  
!</REFERENCE>
!<NOTE>
!     Note by J. Russell: the number of processors you use must be less than the number of
!     atmospheric latitude boxes and must divide into that number equally
!</NOTE>
!
!-----------------------------------------------------------------------
!
!         interface for spherical grid dynamical core and physics
!                   for the energy balance model (EBM)
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------

use mpp_mod,             only: mpp_pe, stdout, stdlog, mpp_error
use fms_mod,             only: file_exist, close_file, read_data, write_data
use fms_mod,             only: error_mesg, FATAL, check_nml_error, open_namelist_file
use fms_mod,             only: write_version_number
use constants_mod,       only: radius, rdgas, rvgas, cp_air, hlv, hlf, tfreeze, stefan, PI
use transforms_mod,      only: transforms_init, grid_domain, spectral_domain
use transforms_mod,      only: trans_grid_to_spherical, trans_spherical_to_grid
use transforms_mod,      only: get_deg_lat, get_eigen_laplacian, compute_gradient_cos
use transforms_mod,      only: compute_div, get_grid_boundaries, divide_by_cos2
use transforms_mod,      only: triangular_truncation, get_grid_domain, get_spec_domain
use transforms_mod,      only: get_wts_lat
use time_manager_mod,    only: time_type, get_time, get_date, operator(+)
use ebm_diagnostics_mod, only: ebm_diagnostics_init, ebm_diagnostics_up, ebm_diagnostics_down
use astronomy_mod,       only: astronomy_init, daily_mean_solar, annual_mean_solar
use sat_vapor_pres_mod,  only: escomp, descomp
use mpp_domains_mod,     only: domain2D
use tracer_manager_mod,  only: get_number_tracers
use xgrid_mod,           only: grid_box_type

!==================================================================================
implicit none
private
!==================================================================================

! version information 

character(len=128), parameter :: version = &
'$Id: atmosphere.F90,v 19.0 2012/01/06 20:00:07 fms Exp $'

character(len=128), parameter :: tagname = &
'$Name: tikal $'

!==================================================================================

! public interfaces

public :: atmosphere_init, atmosphere_down, atmosphere_up
public :: atmosphere_end, atmosphere_resolution, atmosphere_boundary
public :: get_bottom_mass, get_bottom_wind, get_atmosphere_axes
public :: surf_diff_type, radiation, precipitation
public :: diffusion, atmosphere_domain
public :: get_stock_pe, atmosphere_cell_area
public :: atmosphere_restart

!==================================================================================

! module variables

!<PUBLICTYPE>
type surf_diff_type
  real, pointer, dimension(:,:) :: dtmass  => NULL() ! dt/mass, where dt = atmospheric time step (sec)
                                                     ! mass = mass of lowest atmospheric layer (Kg/m2)
  real, pointer, dimension(:,:) :: dflux_t => NULL() ! derivative of the temperature flux at the top of the lowest
                                                     ! atmospheric layer with respect to the temperature
                                                     ! of that layer  (J/(m2 K))
  real, pointer, dimension(:,:,:) :: dflux_tr => NULL() ! derivative of the flux of specific humidity
                                                     ! at the top of the lowest atmospheric layer with respect to
                                                     ! the specific humidity of that layer  (--/(m2 K))
  real, pointer, dimension(:,:) :: delta_t => NULL() ! the increment in temperature in the lowest atmospheric
                                                     ! layer (((i+1)-(i-1) if atmos model is leapfrog) (K)
                                                     ! (defined in  gcm_vert_diff_down as the increment computed up
                                                     ! to this point in model, including effect of vertical
                                                     ! diffusive flux at top of lowest model level, presumed
                                                     ! to be modified to include effects of surface fluxes
                                                     ! outside of this module, then used to start upward
                                                     ! tridiagonal sweep,
  real, pointer, dimension(:,:,:):: delta_tr => NULL() ! similarly for the increment in specific humidity
                                                     ! (non-dimensional  = Kg/Kg)
  real, pointer, dimension(:,:) :: delta_u => NULL()
  real, pointer, dimension(:,:) :: delta_v => NULL()
  real, pointer, dimension(:,:) :: sst_miz => NULL()
end type surf_diff_type
!</PUBLICTYPE>

complex, allocatable, dimension(:,:) :: ts    ! atmospheric temperature  deg_k
complex, allocatable, dimension(:,:) :: qs    ! specific humidity  no units
complex, allocatable, dimension(:,:) :: dt_ts ! rate of change in atmospheric temperature (degK/s)
complex, allocatable, dimension(:,:) :: dt_qs ! rate of change in specific humidity (1/s)
real   , allocatable, dimension(:,:) :: tg    ! atmospheric temperature  ( deg_k )
real   , allocatable, dimension(:,:) :: qg    ! specific humidity ( no units )
real   , allocatable, dimension(:,:) :: dt_tg ! rate of change in atmospheric temperature (degK/s)
real   , allocatable, dimension(:,:) :: dt_qg ! rate of change in specific humidity (1/s)
real   , allocatable, dimension(:,:) :: u_bot ! zonal wind component at lowest model level (m/s)
real   , allocatable, dimension(:,:) :: v_bot ! meridional wind component at lowest model level (m/s)
real   , allocatable, dimension(:,:) :: t_bot ! temperature at lowest model level (deg k)
real   , allocatable, dimension(:,:) :: q_bot ! specific humidity at lowest model level (kg/kg)

real   , allocatable, dimension(:) :: deg_lat, deg_lon, rad_lat, annual_solar, annual_cosz

integer, dimension(4) :: axis_id         ! axes identifiers provided by ebm_diagnostics_mod
real                  :: dt              ! time step (secs)
type(time_type)       :: Time, Time_next, Time_step

integer :: is,ie,js,je,ms,me,ns,ne

real, parameter :: d622 = rdgas/rvgas
real, parameter :: d378 = 1.-d622
logical         :: module_is_initialized =.false.
integer         :: q_ind = 1
!-----------------------------------------------------------------------
!---- namelist (saved in file input.nml) ----
!
! namelist parameters with default values
! Default values are based on ebm run with 1p5 ocean model
!   120 x 65 as run by Joellen Russell

integer  :: lon_max        = 180
integer  :: lat_max        = 120
integer  :: num_fourier    = 42

real    :: diff            = 1.e04
real    :: nu              = 4.0e16
real    :: solar_constant  = 1360.0
real    :: atm_abs         = 0.15
real    :: atm_ref         = 0.05
real    :: mass            = 8.e03
real    :: rh              = 0.8
real    :: t_bot_atm       = 30.0
logical :: seasonal_solar  = .true.
logical :: diff_is_uniform = .false.

real    :: tg_init         = 273.0
real    :: gustiness       = 1.0

logical :: debug_ebm_atm   = .false.

!<NAMELIST NAME="atmosphere_nml">
!  namelist/atmosphere_nml/ lon_max, lat_max, num_fourier, &
!                        mass, rh, t_bot_atm,          &
!                        solar_constant, atm_abs, atm_ref, seasonal_solar, &
!                        diff_is_uniform, nu, diff, tg_init, debug_ebm_atm
!  <DATA NAME="lon_max" TYPE="integer"  DEFAULT="180">
!  This is the longitude of the spectral atmospheric grid
!  </DATA>
!  <DATA NAME="lat_max" TYPE="integer" DEFAULT="120">
!  This is the latitude of the spectral atmospheric grid
!  Note that the number of processors you use must be less than the number of
!  atmospheric latitude boxes and must divide into that number equally
!  </DATA>
!  <DATA NAME="num_fourier" TYPE="integer" DEFAULT="42">
!   num_fourier is used along with lon_max, lat_max in the spectral transforms of grids
!  </DATA>
!  <DATA NAME="diff" UNITS="m^2/sec" TYPE="real" DEFAULT="1.e04">
!   background atmospheric lateral eddy diffusion coefficient for the EBM
!  </DATA>
!  <DATA NAME="nu"  TYPE="real" DEFAULT="4.0e16">
!    scaling factor for spherical grid analysis used in diffusion subroutine
!  </DATA>
!  <DATA NAME="solar_constant" UNITS="Watts/m^2" TYPE="real" DEFAULT="1360.0">
!   standard value for the annual mean solar flux at the top of atmosphere.
!    Default value if seasonal_solar=.false.
!  </DATA>
!  <DATA NAME="atm_abs" TYPE="real" DEFAULT="0.15">
!   atmospheric absorption (fraction)
!  </DATA>
!  <DATA NAME="atm_ref" TYPE="real" DEFAULT="0.05">
!   atmosphere reflectance  (fraction)
!  </DATA>
!  <DATA NAME="mass" UNITS="Kg/m2" TYPE="real" DEFAULT="8.e03">
!   mass = mass of lowest atmospheric layer (Kg/m2)
!  </DATA>
!  <DATA NAME="rh" TYPE="real" DEFAULT="0.8">
!     relative humidity (fraction)
!  </DATA>
!  <DATA NAME="t_bot_atm" UNITS="deg K"  TYPE="real" DEFAULT="30.0">
!    temperature at the bottom of the atmosphere
!  </DATA>
!  <DATA NAME="seasonal_solar" UNITS="Watts/m^2" TYPE="logical" DEFAULT=".true.">
!     If .true., seasonal solar value is computed
!  </DATA>
!  <DATA NAME="diff_is_uniform" TYPE="logical" DEFAULT=".false.">
!  </DATA>
!  <DATA NAME="tg_init" UNITS="deg K" TYPE="real" DEFAULT="273.0">
!   constant value for inital atmospheric temperature  (degK)
!  </DATA>
!  <DATA NAME="debug_ebm_atm" TYPE="logical" DEFAULT=".false.">
!  For debugging purposes.
!  </DATA>

namelist/atmosphere_nml/ lon_max, lat_max, num_fourier, mass, rh, t_bot_atm, &
                         solar_constant, atm_abs, atm_ref, seasonal_solar,   &
                         diff_is_uniform, nu, diff, tg_init, debug_ebm_atm, gustiness
!</NAMELIST>

!==================================================================================

contains



! ==================================================================================
! ==================================================================================

! <SUBROUTINE NAME="atmosphere_init">
!
! <OVERVIEW>
!  public routine required for atmospheric components of coupled models
!  Read in restart files and initialize arrays
! </OVERVIEW>
!
! <DESCRIPTION>
!  public routine required for atmospheric components of coupled models
!  Read in restart files and initialize arrays
! </DESCRIPTION>
!
! <TEMPLATE>
!   call atmosphere_init(Time_init, Time_in, Time_step_in, Surf_diff, Grid_box)
! </TEMPLATE>
!
! <IN NAME="Time_init" TYPE="type(time_type)" >
!   The base (or initial) time of the experiment.
! </IN>

! <IN NAME="Time_in" TYPE="type(time_type)" >
!   The current time.
! </IN>

! <IN NAME="Time_step_in" TYPE="type(time_type)" >
!   The atmospheric model/physics time step.
! </IN>

! <INOUT NAME="Surf_diff" TYPE="type(surf_diff_type)">
!   The surface terms for vertical diffusion that are exchanged
!   with other component models. On input fields have not been
!   allocated, on output all arrays are allocated.
! </INOUT>

! <INOUT NAME="Grid_box" TYPE="type(grid_box_type)">
!   Contains information about the grid cells 
! </INOUT>

subroutine atmosphere_init(Time_init, Time_in, Time_step_in, Surf_diff, Grid_box)

type(time_type),      intent(in)    :: Time_init, Time_in, Time_step_in
type(surf_diff_type), intent(inout) :: Surf_diff
type(grid_box_type),  intent(inout) :: Grid_box

integer                           :: ierr, io, unit, seconds, days, m, n
character(len=64)                 :: filename
real, allocatable, dimension(:,:) :: real_part, imag_part
!-----------------------------------------------------------------------------------------
! managing time

Time      = Time_in
Time_step = Time_step_in

call get_time(Time_step, seconds, days)
dt      = float(86400*days + seconds)

!-----------------------------------------------------------------------------------------
! read namelist and copy to logfile

unit = open_namelist_file()
read  (unit, nml=atmosphere_nml,iostat=io)
write (stdout(),'(/)')
write (stdout(), atmosphere_nml)  
write (stdlog(), atmosphere_nml)
ierr = check_nml_error(io, 'atmosphere_nml')
call close_file(unit)

!--- write out version information
call write_version_number(version, tagname)


!-----------------------------------------------------------------------------------------
! Initialize spherical fourier
call transforms_init(radius, lat_max, lon_max, num_fourier, 1, num_fourier+1)

!-----------------------------------------------------------------------------------------

call get_grid_domain(is,ie,js,je)
call get_spec_domain(ms,me,ns,ne)

! allocate module variables
 
allocate ( tg    (is:ie, js:je) )
allocate ( qg    (is:ie, js:je) )
allocate ( dt_tg (is:ie, js:je) )
allocate ( dt_qg (is:ie, js:je) )

allocate ( u_bot (is:ie, js:je) )
allocate ( v_bot (is:ie, js:je) )
allocate ( t_bot (is:ie, js:je) )
allocate ( q_bot (is:ie, js:je) )

allocate ( ts    (ms:me, ns:ne) )
allocate ( qs    (ms:me, ns:ne) )
allocate ( dt_ts (ms:me, ns:ne) )
allocate ( dt_qs (ms:me, ns:ne) )

allocate ( deg_lon      (is:ie) )
allocate ( deg_lat      (js:je) )
allocate ( rad_lat      (js:je) )

allocate ( annual_cosz  (js:je) )
allocate ( annual_solar (js:je) )

allocate (real_part( ms:me, ns:ne) )
allocate (imag_part( ms:me, ns:ne) )

! hack: in case EBM use more than 1 tracer other than water vapor, revamp this code.
! fil
allocate (Surf_diff%dtmass  (is:ie, js:je) )
allocate (Surf_diff%dflux_t (is:ie, js:je) )
allocate (Surf_diff%dflux_tr (is:ie, js:je, q_ind) )
allocate (Surf_diff%delta_t (is:ie, js:je) )
allocate (Surf_diff%delta_tr (is:ie, js:je, q_ind) )
allocate (Surf_diff%delta_u (is:ie, js:je) )
allocate (Surf_diff%delta_v (is:ie, js:je) )

Surf_diff%dtmass  = 0.0
Surf_diff%dflux_t = 0.0
Surf_diff%dflux_tr = 0.0
Surf_diff%delta_t = 0.0
Surf_diff%delta_tr = 0.0
Surf_diff%delta_u = 0.0
Surf_diff%delta_v = 0.0

!-----------------------------------------------------------------------------------------

! these elements of Surf_diff are static

Surf_diff%dtmass  = dt/mass
Surf_diff%dflux_t = 0.0
Surf_diff%dflux_tr = 0.0

!-----------------------------------------------------------------------------------------

! read restart

filename = 'INPUT/atmosphere.res.nc'

if(file_exist(filename)) then

  call read_data(filename, 'ts_real', real_part, spectral_domain)
  call read_data(filename, 'ts_imag', imag_part, spectral_domain)
  do n = ns, ne
     do m = ms, me
        ts(m,n) = cmplx(real_part(m,n),imag_part(m,n))
     enddo
  enddo
  call read_data(filename, 'qs_real', real_part, spectral_domain)
  call read_data(filename, 'qs_imag', imag_part, spectral_domain)

  do n = ns, ne
     do m = ms, me
        qs(m,n) = cmplx(real_part(m,n),imag_part(m,n))
     enddo
  enddo

  call read_data(filename, 'tg', tg, grid_domain)
  call read_data(filename, 'qg', qg, grid_domain)
else
! tg-initial atmospheric temperature default value 273 deg Kelvin
! qg-initial atmospheric specific humidity
  tg = tg_init
  qg = rh*sat_mixing_ratio(tg)

  call trans_grid_to_spherical(tg, ts)
  call trans_grid_to_spherical(qg, qs)

endif

v_bot = 0.0
t_bot = tg + t_bot_atm
q_bot = rh*sat_mixing_ratio(t_bot)  

! initalize other modules

call ebm_diagnostics_init(Time, lon_max, lat_max, axis_id)
call astronomy_init()

! get latitudes of grid points in degrees and radians
call get_deg_lat(deg_lat)
rad_lat = deg_lat*PI/180.

! if annual mean option is used, compute solar radiation here
!  Use version for  spectral 2-layer model 
if(.not.seasonal_solar) call annual_mean_solar ( rad_lat, annual_cosz, annual_solar ) 

module_is_initialized = .true.

return
end subroutine atmosphere_init
! </SUBROUTINE>

! ==================================================================================
! <SUBROUTINE NAME="atmosphere_down">
! <OVERVIEW>
!  public routine required for atmospheric components of coupled models
! </OVERVIEW>

! <DESCRIPTION>
!   atmosphere_down
!  This routine calls the dynamical core and the
!  "downward pass" of the atmospheric physics.
!  It should only be called once per time step and before
!  calling atmosphere_up.
! </DESCRIPTION>

! <TEMPLATE>
!     call atmosphere_down(Time, frac_land, t_surf, albedo, rough_mom, &
!                          u_star, b_star, q_star, dtau_du, dtau_dv, tau_x, tau_y,  &
!                          gust, coszen, net_surf_sw_down, surf_lw_down, &
!                          Surf_diff)
! </TEMPLATE>

! <IN NAME="Time" TYPE="type(time_type)">
!  time at the current time level (tau)
! </IN>

! <IN NAME="frac_land" TYPE="real, dimension(is:ie,js:je)">
!  fraction (0. to 1.) of underlying surface which covered by land
! </IN>

! <IN NAME="t_surf" TYPE="real, dimension(is:ie,js:je)">
!  surface (skin) temperature (in deg k)
! </IN>

! <IN NAME="albedo" TYPE="real, dimension(is:ie,js:je)">
!  surface albedo
! </IN>

! <IN NAME="rough_mom" TYPE="real, dimension(is:ie,js:je)">
!  momentum roughness length (units=m)
! </IN>

! <IN NAME="u_star" TYPE=" real, dimension(is:ie,js:je)">
!  friction velocity
! </IN>

! <IN NAME="b_star" TYPE="real, dimension(is:ie,js:je) ">
!  buoyancy scale
! </IN>

! <IN NAME="q_star" TYPE="real, dimension(is:ie,js:je) ">
!  moisture scale
! </IN>

! <IN NAME="dtau_du" TYPE=" real, dimension(is:ie,js:je)">
!  derivative of zonal wind stress w.r.t. the lowest level wind speed
! </IN>

! <IN NAME="dtau_dv" TYPE=" real, dimension(is:ie,js:je)">
!  derivative of meridional wind stress w.r.t. the lowest level wind speed
! </IN>

! <OUT NAME="gust" TYPE="real, dimension(is:ie,js:je)" >
!  wind gustiness
! </OUT>

! <OUT NAME="coszen" TYPE="real, dimension(is:ie,js:je)" >
!  cosine of the zenith angle
! </OUT>

! <OUT NAME="net_surf_sw_down" TYPE="real, dimension(is:ie,js:je)" >
!  net shortwave surface flux (down minus up) (in watts/m**2)
! </OUT>

! <OUT NAME="surf_lw_down" TYPE="real, dimension(is:ie,js:je)" >
!  downward longwave surface flux (in watts/m**2)
! </OUT>

! <INOUT NAME=" Surf_diff" TYPE="type(surf_diff_type)" >
!  Surface diffusion terms computed by the vertical diffusion scheme
! </INOUT>

subroutine atmosphere_down(Time, frac_land, t_surf, albedo, albedo_vis_dir, &
                           albedo_nir_dir, albedo_vis_dif,  albedo_nir_dif, &
                           rough_mom, u_star, b_star, q_star, dtau_du, dtau_dv, &
                           tau_x, tau_y,frac_open_sea, gust, coszen, net_surf_sw_down,    &
                           flux_sw_dir, flux_sw_dif, flux_sw_down_vis_dir,  &
                           flux_sw_down_vis_dif, flux_sw_down_total_dir,    &
                           flux_sw_down_total_dif, flux_sw_vis,             &
                           flux_sw_vis_dir, flux_sw_vis_dif,surf_lw_down,   &
                           Surf_diff)

!  public routine required for atmospheric components of coupled models
type(time_type),  intent(in)                 :: Time
real, intent(in),     dimension(is:ie,js:je) :: frac_land, t_surf, albedo, rough_mom
real, intent(in),     dimension(is:ie,js:je) :: albedo_vis_dir, albedo_nir_dir
real, intent(in),     dimension(is:ie,js:je) :: albedo_vis_dif, albedo_nir_dif
real, intent(in),     dimension(is:ie,js:je) :: u_star, b_star, q_star, dtau_du, dtau_dv, frac_open_sea
real,  intent(inout), dimension(is:ie,js:je) :: tau_x, tau_y
real,  intent(out),   dimension(is:ie,js:je) :: gust, coszen, net_surf_sw_down, surf_lw_down
real, intent(out),    dimension(is:ie,js:je) :: flux_sw_dir, flux_sw_dif, flux_sw_down_vis_dir
real, intent(out),    dimension(is:ie,js:je) :: flux_sw_down_vis_dif, flux_sw_down_total_dir
real, intent(out),    dimension(is:ie,js:je) :: flux_sw_down_total_dif, flux_sw_vis
real, intent(out),    dimension(is:ie,js:je) :: flux_sw_vis_dir, flux_sw_vis_dif
type(surf_diff_type),   intent(inout)        :: Surf_diff

! initialize tendencies of temperature and mixing ratio to zero 
dt_tg = 0.0
dt_qg = 0.0

! add temperature tendency due to radiation, and compute radiative fluxes and
!    zenith angle needed by flux module

call radiation (Time, t_surf, albedo, coszen, net_surf_sw_down, surf_lw_down)

!
!           VIS                NIR          VIS+NIR
!-----------------------------------------------------
!       |                 |              |
! DIR   |flux_sw_vis_dir  |              | flux_sw_dir
!       |                 |              |
!-----------------------------------------------------
!       |                 |              |
! DIF   |flux_sw_vis_dif  |              | flux_sw_dif
!-----------------------------------------------------
!       |                 |              |
!DIR+DIF|flux_sw_vis      |              | flux_sw        
!       |                 |              |(diagnostic only)
!-----------------------------------------------------

! bls-total Short wave is equally divided into VIS and NIR
! flux_sw_dif - total Diffusive Short wave : VIS + NIR
flux_sw_dif = .5*net_surf_sw_down
! flux_sw_dir - total Direct Short wave : VIS + NIR
flux_sw_dir = net_surf_sw_down - flux_sw_dif

! bls- Visible Direct and Diffusive short wave will be assumed to be
! half of the total VIS+DIR components
flux_sw_vis_dir = .5*flux_sw_dir
flux_sw_vis_dif = .5*flux_sw_dif
! Total Visible short wave (DIR + DIFFUSE)
 flux_sw_vis = flux_sw_vis_dir + flux_sw_vis_dif
! NIR components are determined in the ice model

! for surface flux computation 
gust    = gustiness
v_bot   = 0.0
t_bot   = tg + t_bot_atm  
q_bot   = rh*sat_mixing_ratio(t_bot)  

Surf_diff%delta_t = dt*dt_tg   ! needed in this form for surface models
Surf_diff%delta_tr(:,:,q_ind) = dt*dt_qg

return
end subroutine atmosphere_down
! </SUBROUTINE>


! ==================================================================================
! <SUBROUTINE NAME="atmosphere_up">
! <OVERVIEW>
!  public routine required for atmospheric components of coupled models
!</OVERVIEW>
!
!<DESCRIPTION>
!   atmosphere_up
!  This routine calls the "upward pass" of the atmospheric physics,
!  spherical diagnostics, and time differencing.  The prognostic
!  variables are advanced to the next time step.  It should only be
!  called once per time step and after calling atmosphere_down.
!</DESCRIPTION>

!   <TEMPLATE>
!     call  atmosphere_up(Time, frac_land, Surf_diff, lprec, fprec, gust)
!   </TEMPLATE>

! <IN NAME="Time" TYPE="type(time_type)">
!  time at the current time level (tau)
! </IN>

! <IN NAME="frac_land" TYPE="real, dimension(is:ie,js:je)">
!  fraction (0. to 1.) of underlying surface which covered by land
! </IN>

! <OUT NAME="lprec" TYPE="real, dimension(is:ie,js:je)" >
!  liquid precipitation rate (rain) in kg/m2/s
! </OUT>

! <OUT NAME="fprec" TYPE="real, dimension(is:ie,js:je)" >
!  frozen precipitation rate (snow) in kg/m2/s
! </OUT>

! <OUT NAME="gust" TYPE="real, dimension(is:ie,js:je)" >
!  wind gustiness
! </OUT>

! <INOUT NAME=" Surf_diff" TYPE="type(surf_diff_type)" >
!  urface diffusion terms computed by the vertical diffusion scheme
! </INOUT>

! <IN NAME="u_star" TYPE="real,dimension(:,:)">
!  NOT USED in this routine. Dummy argument added for compilation success.
! </IN>
! <IN NAME="b_star" TYPE="real,dimension(:,:)">
!  NOT USED in this routine. Dummy argument added for compilation success.
! </IN>
! <IN NAME="q_star" TYPE="real,dimension(:,:)">
!  NOT USED in this routine. Dummy argument added for compilation success.
! </IN>

subroutine atmosphere_up(Time,  frac_land, Surf_diff, lprec, fprec, gust, u_star, b_star, q_star )

    
type(time_type),      intent(in)                          :: Time
real,                 intent(in),  dimension(is:ie,js:je) :: frac_land
type(surf_diff_type), intent(in)                          :: Surf_diff
real,                 intent(out), dimension(is:ie,js:je) :: lprec, fprec, gust
real,                 intent(in),  dimension(:,:)         :: u_star, b_star, q_star

complex, dimension(ms:me,ns:ne) :: lap, spec
real,    dimension(ms:me,ns:ne) :: eigen

integer :: year, month, day, hour, minute, second

dt_tg = Surf_diff%delta_t/dt
dt_qg = Surf_diff%delta_tr(:,:,q_ind)/dt

!initialize arrays to zero
 gust = gustiness

call precipitation(lprec, fprec)
call diffusion

Time_next = Time + Time_step
     if(debug_ebm_atm) then
 call get_date(Time_next, year, month, day, hour, minute, second)
 if (mpp_pe() == 0 .and. hour+minute == 0 ) then
    print *,"atmosphere_up:Current time: year =",year," month = ",month," day = ",day
    write (stdout(),1000) tg(ie,je), qg(ie,je), lprec(ie,je), fprec(ie,je)
 endif
     endif
call ebm_diagnostics_up(Time_next, tg, qg, lprec, fprec, dt_qg)
     if(debug_ebm_atm) then
! print out a few things for degbugging
 if (mpp_pe() == 0 .and. hour+minute == 0 ) then 
    print *,"atmosphere_up aft diag:Current time: year =",year," month = ",month," day = ",day
    write (stdout(),1000) tg(ie,je), qg(ie,je), lprec(ie,je), fprec(ie,je)
 endif
 1000 format(4(1x,e11.4))
     endif


return
end subroutine atmosphere_up
! </SUBROUTINE>

! ==================================================================================
! <SUBROUTINE NAME="radiation">
! <DESCRIPTION>
!    radiation
!    Radiation module for prognostic equation of atmospheric temperature
!    The atmospheric tendency term (dt_tg) includes the balance of shortwave
!    and longwave radiation terms at the surface
! </DESCRIPTION>
!
! <TEMPLATE>
!  call radiation (Time, t_surf, albedo, coszen, net_surf_sw_down, surf_lw_down)
! </TEMPLATE>

! <IN NAME="Time" TYPE="type(time_type)">
!  time at the current time level (tau)
! </IN>

! <IN NAME="t_surf" TYPE="real, dimension(is:ie,js:je)">
!  surface (skin) temperature (in deg k)
! </IN>

! <IN NAME="albedo" TYPE="real, dimension(is:ie,js:je)">
!  surface albedo
! </IN>

! <OUT NAME="coszen" TYPE="real, dimension(is:ie,js:je)" >
!  cosine of the zenith angle
! </OUT>

! <OUT NAME="net_surf_sw_down" TYPE="real, dimension(is:ie,js:je)" >
!  net shortwave surface flux (down minus up) (in watts/m**2)
! </OUT>

! <OUT NAME="surf_lw_down" TYPE="real, dimension(is:ie,js:je)" >
!  downward longwave surface flux (in watts/m**2)
! </OUT>

subroutine radiation (Time, t_surf, albedo, coszen, net_surf_sw_down, surf_lw_down)
!
type(time_type),      intent(in)                    :: Time
real,                 intent(in),    dimension(is:ie,js:je) :: t_surf, albedo
real,                 intent(out),   dimension(is:ie,js:je) :: coszen, net_surf_sw_down, surf_lw_down

integer :: j

real, dimension(is:ie,js:je) :: net_top_sw_down, top_lw_up, surf_lw_up
real, dimension(js:je) :: cosz, solar, incident_solar

integer :: year, month, day, hour, minute, second

if(seasonal_solar) then 
! Use spectral 2-layer model version
  call daily_mean_solar  (rad_lat, Time, cosz, solar) 
else
  cosz  = annual_cosz
  solar = annual_solar
endif
incident_solar =  solar_constant*solar

     if(debug_ebm_atm) then
 call get_date(Time, year, month, day, hour, minute, second)
     endif

do j = js, je

  net_top_sw_down  (:,j) = incident_solar(j)*(1.0 - atm_ref - (1.0 - atm_ref - atm_abs)*albedo(:,j))  
  net_surf_sw_down (:,j) = incident_solar(j)*(1.0 - atm_ref- atm_abs)*(1. - albedo(:,j))              
  top_lw_up        (:,j) = stefan*( tg     (:,j)      **4)                                           
  surf_lw_down     (:,j) = stefan*((tg     (:,j)+20.0)**4)                                         
  surf_lw_up       (:,j) = stefan*( t_surf (:,j)      **4)                                         

  dt_tg            (:,j) = dt_tg(:,j) +                                     &
                       ( + surf_lw_up     (:,j)                             &
                         - top_lw_up      (:,j) - surf_lw_down    (:,j)     &
                         + net_top_sw_down(:,j) - net_surf_sw_down(:,j) )/(cp_air*mass)

  coszen           (:,j) = cosz(j) 

     if(debug_ebm_atm) then
 if (mpp_pe() == 0 .and. hour+minute == 0 ) then
    print *,"radiation:Current time: year =",year," month = ",month," day = ",day, " j=",j
    write (stdout(),1000) incident_solar(j),albedo(ie,j),atm_ref,atm_abs
    write (stdout(),1000) tg(ie,j), qg(ie,j), net_top_sw_down(ie,j), net_surf_sw_down(ie,j)
    write (stdout(),1000) top_lw_up(ie,j), surf_lw_down(ie,j), surf_lw_up(ie,j), dt_tg(ie,j)
 endif
 1000 format(4(1x,e11.4))
      endif
enddo

call ebm_diagnostics_down(Time, top_lw_up, surf_lw_down, surf_lw_up, net_top_sw_down, net_surf_sw_down)

     if(debug_ebm_atm) then
 if (mpp_pe() == 0 .and. hour+minute == 0 ) then
    print *,"radiation aft diag:Current time: year =",year," month = ",month," day = ",day, " j=",j-1
    write (stdout(),1000) incident_solar(j-1),albedo(ie,j-1),atm_ref,atm_abs
    write (stdout(),1000) tg(ie,j-1), qg(ie,j-1), net_top_sw_down(ie,j-1), net_surf_sw_down(ie,j-1)
    write (stdout(),1000) top_lw_up(ie,j-1), surf_lw_down(ie,j-1), surf_lw_up(ie,j-1), dt_tg(ie,j-1)
 endif
      endif

return
end subroutine radiation
! </SUBROUTINE>


! ==================================================================================
! <SUBROUTINE NAME="precipitation">
!
! <OVERVIEW>
!    Precipiatation module for prognostic equation of atmospheric temperature
!     and specific humidity
! </OVERVIEW>

! <DESCRIPTION>
!    Precipiatation module for prognostic equation of atmospheric temperature
!     and specific humidity
!    The atmospheric tendency term (dt_tg) includes the latent heat flux released
!     during precipitation (tdel)
!    The specific humidity is the balance of evaporation and precipitation (dt_qg)
! </DESCRIPTION>

!   <TEMPLATE>
!     call   precipitation(lprec, fprec)
!   </TEMPLATE>

! <OUT NAME="lprec" TYPE="real, dimension(is:ie,js:je)" >
!  liquid precipitation rate (rain) in kg/m2/s
! </OUT>

! <OUT NAME="fprec" TYPE="real, dimension(is:ie,js:je)" >
!  frozen precipitation rate (snow) in kg/m2/s
! </OUT>

subroutine precipitation(lprec, fprec)

real, intent(out), dimension(is:ie,js:je) :: lprec, fprec
real,              dimension(is:ie,js:je) :: qsat, dqsat, tdel, qdel, &
                                             tg_temp, qg_temp
real :: hlcp, hfcp
integer :: i,j

integer :: year, month, day, hour, minute, second

tg_temp = tg + dt_tg*dt + t_bot_atm  
qg_temp = qg + dt_qg*dt

qsat   =       sat_mixing_ratio(tg_temp)   
dqsat  = deriv_sat_mixing_ratio(tg_temp)
! tfreeze - temperature where fresh water freezes 273.16 degK 
! hlv - Latent heat of evaporation J/kg
! hlf - latent heat of fusion      J/kg condensate
! cp_air  - specific heat of air at         J/kg air/K
hlcp   =  hlv       /cp_air
hfcp   = (hlv + hlf)/cp_air
! initialize precipiation to zero
lprec = 0.0
fprec = 0.0

where (qg_temp > qsat .and. tg_temp > tfreeze)
  qdel = (qsat - qg_temp)/(1.0 + hlcp*dqsat)
  tdel = -hlcp*qdel
  lprec = -qdel*mass/dt
elsewhere (qg_temp > qsat .and. tg_temp <= tfreeze)
  qdel = (qsat - qg_temp)/(1.0 + hfcp*dqsat)
  tdel = -hfcp*qdel
  fprec = -qdel*mass/dt
elsewhere
  qdel = 0.0
  tdel = 0.0
endwhere
dt_tg = dt_tg + tdel/dt
dt_qg = dt_qg + qdel/dt

     if(debug_ebm_atm) then
 if (mpp_pe() == 0 ) then
    print *,"in precip:"
    write (stdout(),1000)  qg_temp(ie,je), qsat(ie,je), tg_temp(ie,je), dqsat(ie,je)
    write (stdout(),1000) tg(ie,je), qg(ie,je), lprec(ie,je), fprec(ie,je)
    write (stdout(),1000) dt_tg(ie,je), dt_qg(ie,je), tdel(ie,je), qdel(ie,je)
 endif
 1000 format(4(1x,e11.4))
     endif

return
end subroutine precipitation
! </SUBROUTINE>

! ==================================================================================
! <SUBROUTINE NAME="diffusion">
!
! <OVERVIEW>
!    Lateral eddy diffusion module for prognostic equation of atmospheric temperature
!     and specific humidity
! </OVERVIEW>

! <DESCRIPTION>
!    Lateral eddy diffusion module for prognostic equation of atmospheric temperature
!     and specific humidity
!    There are two options: uniform or variable diffusion
! </DESCRIPTION>

!   <TEMPLATE>
!     call diffusion()
!   </TEMPLATE>

subroutine diffusion  

real   , dimension (ms:me,ns:ne) :: eigen, denom 
complex, dimension (ms:me,ns:ne) :: dxs, dys, div
real   , dimension (is:ie,js:je) :: dxg, dyg, d

integer :: m, n

! real :: nu = 4.0e16

call trans_grid_to_spherical(dt_tg, dt_ts)
call trans_grid_to_spherical(dt_qg, dt_qs)

call get_eigen_laplacian(eigen)

if(diff_is_uniform) then

  dt_ts = dt_ts - diff*eigen*ts
  dt_qs = dt_qs - diff*eigen*qs

else

  d = diffusivity()  

! temperature

  call compute_gradient_cos(ts, dxs, dys)

  call trans_spherical_to_grid(dxs,dxg)
  call trans_spherical_to_grid(dys,dyg)

  dxg = d*dxg
  dyg = d*dyg

  call divide_by_cos2(dxg)
  call divide_by_cos2(dyg)

  call trans_grid_to_spherical(dxg,dxs)
  call trans_grid_to_spherical(dyg,dys)

  div = compute_div(dxs, dys)
  call triangular_truncation(div)

  dt_ts = dt_ts + div

! water vapor

  call compute_gradient_cos(qs, dxs, dys)

  call trans_spherical_to_grid(dxs,dxg)
  call trans_spherical_to_grid(dys,dyg)

  dxg = d*dxg
  dyg = d*dyg

  call divide_by_cos2(dxg)
  call divide_by_cos2(dyg)

  call trans_grid_to_spherical(dxg,dxs)
  call trans_grid_to_spherical(dyg,dys)

  div = compute_div(dxs, dys)
  call triangular_truncation(div)

  dt_qs = dt_qs + div

endif


dt_ts = dt_ts - nu*eigen*eigen*ts
dt_qs = dt_qs - nu*eigen*eigen*qs


denom = 1./(1.0 + (diff*eigen + nu*eigen*eigen)*dt)

ts = ts + denom*dt*dt_ts
qs = qs + denom*dt*dt_qs

call trans_spherical_to_grid(ts,tg)
call trans_spherical_to_grid(qs,qg)

     if(debug_ebm_atm) then
 if (mpp_pe() == 0 ) then
    print *,"diffusion:"
    write (stdout(),1000) tg(ie,je), qg(ie,je), denom(me,ne), eigen(me,ne)
    write (stdout(),1000) ts(me,ne), qs(me,ne), dt_qs(me,ne), dt_ts(me,ne)
 endif
 1000 format(4(1x,e11.4))
      endif

return
end subroutine diffusion
! </SUBROUTINE>

! ==================================================================================
! <FUNCTION NAME="diffusivity">
!
! <OVERVIEW>
!  For the a given scalar diffusion coefficient the routine returns an array
! </OVERVIEW>

! <DESCRIPTION>
!  For the a given scalar diffusion coefficient the routine returns an array
! </DESCRIPTION>

!   <TEMPLATE>
!     d = diffusivity()
!   </TEMPLATE>

function diffusivity() result(d)

real, dimension(is:ie,js:je) :: d
d = diff

return
end function diffusivity
! </FUNCTION>


!#######################################################################
! <FUNCTION NAME="sat_mixing_ratio">
!
! <OVERVIEW>
!  For the given temperatures routine returns the saturation mixing ratio
!   esat - saturation vapor pressure
! </OVERVIEW>
!   <TEMPLATE>
!     q = sat_mixing_ratio(t)
!   </TEMPLATE>
! <DESCRIPTION>
!  For the given temperatures routine returns the saturation mixing ratio
!   esat - saturation vapor pressure
! </DESCRIPTION>

function sat_mixing_ratio(t) result(q)

real, dimension(is:ie,js:je) :: t, q, esat

call escomp (t, esat)
q = d622*esat/1.e05   

return
end function sat_mixing_ratio
! </FUNCTION>

! ==================================================================================
! <FUNCTION NAME="deriv_sat_mixing_ratio">
!
! <OVERVIEW>
!  For the given temperatures, routine returns the derivative of the  saturation mixing ratio
! </OVERVIEW>

! <DESCRIPTION>
!  For the given temperatures, routine returns the derivative of the  saturation mixing ratio
! </DESCRIPTION>

!   <TEMPLATE>
!     q = deriv_sat_mixing_ratio(t)
!   </TEMPLATE>
function deriv_sat_mixing_ratio(t) result(q)

real, dimension(is:ie,js:je) :: t, q, desat

call descomp (t, desat)
q = d622*desat/1.e05   

return
end function deriv_sat_mixing_ratio
! </FUNCTION>

! ==================================================================================
! <SUBROUTINE NAME="get_bottom_mass">
!
! <OVERVIEW>
! returns temp, sphum, pres, height at the lowest model level
!         and surface pressure
! </OVERVIEW>

! <DESCRIPTION>
!  public routine required for atmospheric components of coupled models
!  returns temp, sphum, pres, height at the lowest model level
!         and surface pressure
! </DESCRIPTION>

!   <TEMPLATE>
!     call get_bottom_mass (t_bot_out, q_bot_out, p_bot, z_bot, p_surf)
!   </TEMPLATE>

!   <OUT NAME="t_bot_out" TYPE="real">
!                near surface temperature in degrees Kelvin
!   </OUT>

!   <OUT NAME="q_bot_out" TYPE="real">
!                near surface mixing  ratio
!   </OUT>

!   <OUT NAME="p_bot" TYPE="real">
!                 pressure at which atmos near usrface values are assumed to be defined
!   </OUT>

!   <OUT NAME="z_bot" TYPE="real">
!                 height at which atmos near usrface values are assumed to be defined
!   </OUT>

!   <OUT NAME="p_surf" TYPE="real">
!                 surface pressure
!   </OUT>

subroutine get_bottom_mass (t_bot_out, tr_bot_out, p_bot, z_bot, p_surf, slp)
real, intent(out), dimension(is:ie,js:je) :: t_bot_out, p_bot, z_bot, p_surf, slp
real, intent(out), dimension(:,:,:) :: tr_bot_out

p_surf    = 1000.0e02   ! surface pressure
p_bot     = 990.0e02    ! pressure at which atmos near usrface values are assumed to be defined
z_bot     = 100.0       !   height at which atmos near surface values are assumed to be defined
t_bot_out = t_bot(is:ie,js:je)       ! near surface temperature
tr_bot_out(:,:,1) = q_bot(is:ie,js:je)       ! near surface mixing ratio
slp = p_surf

return
end subroutine get_bottom_mass
! </SUBROUTINE>

! ==================================================================================
! <SUBROUTINE NAME="get_bottom_wind">
!
! <OVERVIEW>
! returns u and v on the mass grid at the lowest model level
! </OVERVIEW>

! <DESCRIPTION>
!  public routine required for atmospheric components of coupled models
!  returns u and v on the mass grid at the lowest model level
! </DESCRIPTION>

!   <TEMPLATE>
!     call get_bottom_wind (u_bot_out, v_bot_out)
!   </TEMPLATE>

!   <OUT NAME="u_bot"  TYPE="real" >
!    near surface zonal      wind
!   </OUT>

!   <OUT NAME="v_bot"  TYPE="real" >
!    near surface meridional wind
!   </OUT>

subroutine get_bottom_wind (u_bot_out, v_bot_out)
real, intent(out), dimension(is:ie,js:je) :: u_bot_out, v_bot_out

u_bot_out = u_bot(is:ie,js:je)   ! near surface zonal      wind
v_bot_out = v_bot(is:ie,js:je)   ! near surface meridional wind

return
end subroutine get_bottom_wind
! </SUBROUTINE>

! ==================================================================================
! <SUBROUTINE NAME="atmosphere_resolution">
! <OVERVIEW>
!    returns the number of longitude and latitude grid points
!    for either the local PEs grid (default) or the global grid
! </OVERVIEW>

! <DESCRIPTION>
!  public routine required for atmospheric components of coupled models
!    returns the number of longitude and latitude grid points
!    for either the local PEs grid (default) or the global grid
! </DESCRIPTION>

! <TEMPLATE>
!   call  atmosphere_resolution(num_lon_out, num_lat_out, global)
! </TEMPLATE>

! <OUT NAME="num_lon_out" TYPE="integer">
!  The number of longitude points in the compute domain.
! </OUT>

! <OUT NAME="num_lat_out" TYPE="integer">
!  The number of latitude points in the compute domain.
! </OUT>

! <IN NAME="global" TYPE="logical,optional">
!  Flag that specifies whether the returned compute domain size is
!  for the global grid (TRUE) or for the current processor (FALSE).
! </IN>

subroutine atmosphere_resolution(num_lon_out, num_lat_out, global)
 
integer, intent(out)          :: num_lon_out, num_lat_out
logical, intent(in), optional :: global
logical :: global_tmp

if (present(global)) then
  global_tmp = global
else
  global_tmp = .false.
endif

if(global_tmp) then
  num_lon_out = lon_max
  num_lat_out = lat_max
else
  num_lon_out = ie+1-is
  num_lat_out = je+1-js
endif

return
end subroutine atmosphere_resolution
! </SUBROUTINE>

! ==================================================================================
! <SUBROUTINE NAME="get_atmosphere_axes">
!
! <OVERVIEW>
!    returns the axis indices associated with the coupling grid
! </OVERVIEW>

! <DESCRIPTION>
!  public routine required for atmospheric components of coupled models
!    returns the axis indices associated with the coupling grid
! </DESCRIPTION>

! <TEMPLATE>
!     call  get_atmosphere_axes(axes_out)
! </TEMPLATE>

! <OUT NAME="axes_out" TYPE="integer,dimension(:)" >
!  The axis identifiers for the atmospheric grids.
!  The size of axes must be least 1 but not greater than 4.
!  The axes are returned in the order (/ x, y, p_full, p_half /)
!</OUT>

subroutine get_atmosphere_axes(axes_out)
integer, intent(out), dimension(:) :: axes_out

axes_out = axis_id

return
end subroutine get_atmosphere_axes
! </SUBROUTINE>

! ==================================================================================
! <SUBROUTINE NAME="atmosphere_boundary">
! <OVERVIEW>
!    returns the longitude and latitude grid box edges
!    for either the local PEs grid (default) or the global grid
! </OVERVIEW>
! <DESCRIPTION>
!  public routine required for atmospheric components of coupled models
!    returns the longitude and latitude grid box edges
!    for either the local PEs grid (default) or the global grid
! </DESCRIPTION>

! <TEMPLATE>
!   call atmosphere_boundary(lon_boundaries, lat_boundaries, global)
! </TEMPLATE>

! <OUT NAME = "lon_boundaries" TYPE = "real,dimension(:,:)"> 
!  The west-to-east longitude edges of grid boxes (in radians).
! </OUT>

! <OUT NAME = "lat_boundaries" TYPE = "real,dimension(:,:)"> 
!  The south-to-north latitude edges of grid boxes (in radians).
! </OUT>

! <IN NAME="global" TYPE="logical, optional">
!  Flag that specifies whether the returned grid box edges are
!  for the global grid (TRUE) or for the current processor (FALSE).
! </IN>
subroutine atmosphere_boundary(lon_boundaries, lat_boundaries, global)

real,    intent(out), dimension(:,:) :: lon_boundaries, lat_boundaries
logical, intent(in),  optional     :: global
real, dimension(size(lon_boundaries,1)) :: tmpx
real, dimension(size(lon_boundaries,2)) :: tmpy
integer :: i
logical :: global_tmp


if(present(global)) then
  global_tmp = global
else
  global_tmp = .false.
endif
call get_grid_boundaries(tmpx, tmpy, global_tmp)

do i = 1, size(lon_boundaries,2)
  lon_boundaries(:,i) = tmpx(:)
end do

do i = 1, size(lon_boundaries,1)
  lat_boundaries(i,:) = tmpy(:)
end do

return
end subroutine atmosphere_boundary
! </SUBROUTINE>

! ==================================================================================
! <SUBROUTINE NAME="atmosphere_cell_area">
!
! <OVERVIEW>
!   returns grid cell areas for the domain
! </OVERVIEW>
!
! <DESCRIPTION>
!   returns grid cell areas for the domain
! </DESCRIPTION>

 subroutine atmosphere_cell_area(area_out)
 real, dimension(:,:),  intent(out) :: area_out
 integer :: xsize, ysize, j
 real, dimension(size(area_out,2)) :: wts_lat

 xsize = ie-is+1
 ysize = je-js+1
 if(any(shape(area_out) /= (/xsize,ysize/))) then
   call mpp_error(FATAL,'atmosphere_cell_area: argument has wrong shape. Its shape is ', &
     shape(area_out),'  It should be ',(/xsize,ysize/))
 endif

 call get_wts_lat(wts_lat)
 do j=1,ysize
   area_out(:,j) = 2*PI*radius*radius*wts_lat(j)/lon_max
 enddo

 end subroutine atmosphere_cell_area

! </SUBROUTINE>

! ==================================================================================
! <SUBROUTINE NAME="atmosphere_domain">
!
! <OVERVIEW>
!    returns the domain2d variable associated with the coupling grid
! </OVERVIEW>
!
! <DESCRIPTION>
!  public routine required for atmospheric components of coupled models
!    returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos
!
!OUTPUT
!  Domain   The domain2d variable describing the grid used for coupling.
!           For the B-grid, this corresponds to the temperature grid
!           without halos.
! </DESCRIPTION>

! <TEMPLATE>
!     call atmosphere_domain(domain)
! </TEMPLATE>

! <INOUT NAME="domain" TYPE="type (domain2d)">
!  The domain2d variable describing the grid used for coupling.
!  For the B-grid, this corresponds to the temperature grid
!  without halos.
! </INOUT>

subroutine atmosphere_domain(domain)

type (domain2d), intent(inout) :: domain

   
domain = grid_domain;
   
end subroutine atmosphere_domain
! </SUBROUTINE>

! =====================================================================
! <SUBROUTINE NAME="atmosphere_end">
!
! <OVERVIEW>
!  write out restart file
!  Termination routine for atmosphere_mod.
! </OVERVIEW>
! <DESCRIPTION>
!  public routine required for atmospheric components of coupled models
!  write out restart file
!  Termination routine for atmosphere_mod.
! </DESCRIPTION>

! <TEMPLATE>
!     call atmosphere_end(Time, Grid_box)
! </TEMPLATE>

! <IN NAME="Time" TYPE="type(time_type)">
!  time at the current time level (tau)
! </IN>

! <INOUT NAME="Grid_box" TYPE="type(grid_box_type)">
!   Contains information about the grid cells 
! </INOUT>

subroutine atmosphere_end(Time, Grid_box)

type(time_type), intent(in) :: Time
type(grid_box_type), intent(inout) :: Grid_box
integer :: time_level,seconds,days,unit
character(len=64)            :: filename

if(.not.module_is_initialized) then
  call error_mesg('atmosphere_end',' atmosphere_init has not been called.', FATAL)
end if

!--- write out restart file
filename = 'RESTART/atmosphere.res.nc'

call write_data(filename, 'ts_real', real(ts), spectral_domain)
call write_data(filename, 'ts_imag', aimag(ts), spectral_domain)
call write_data(filename, 'qs_real', real(qs), spectral_domain)
call write_data(filename, 'qs_imag', aimag(qs), spectral_domain)
call write_data(filename, 'tg', tg, grid_domain)
call write_data(filename, 'qg', qg, grid_domain)

return
end subroutine atmosphere_end
! </SUBROUTINE>

  !#######################################################################
  ! <SUBROUTINE NAME="atmosphere_restart">
  ! <DESCRIPTION>
  !  dummy routine.
  ! </DESCRIPTION>
  subroutine atmosphere_restart(timestamp)
    character(len=*),  intent(in) :: timestamp

    call error_mesg ('atmosphere_restart in atmosphere_mod', &
                     'intermediate restart capability is not implemented for this model', FATAL)

  end subroutine atmosphere_restart
  ! </SUBROUTINE>

! =====================================================================
! <SUBROUTINE NAME="get_stock_pe">
!
! <OVERVIEW>
!  Not implemented properly. Puts value=0
! </OVERVIEW>
! <DESCRIPTION>
!  This is a stub routine needed for compilation of mom4_test6.
!  It puts value=0 for any index.
! </DESCRIPTION>

! <TEMPLATE>
!     call get_stock_pe(index, value)
! </TEMPLATE>

! <IN NAME="index" TYPE="integer">
! </IN>
! <OUT NAME="value" TYPE="real">
! </OUT>

subroutine get_stock_pe(index, value)
  use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT
  integer, intent(in) :: index
  real, intent(out)   :: value

  select case (index)

  case (ISTOCK_WATER)
     ! don't know how to compute water stock
     value = 0

  case (ISTOCK_HEAT)
     ! don't know how to compute water stock
     value = 0

  case default
     value = 0.0
  end select

end subroutine get_stock_pe

! ==================================================================================
! ==================================================================================

end module atmosphere_mod

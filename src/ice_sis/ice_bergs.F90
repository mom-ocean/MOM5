module ice_bergs

use constants_mod, only: radius, pi, omega, HLF
use fms_mod, only: open_namelist_file, check_nml_error, close_file
use fms_mod, only: field_exist, get_global_att_value
use fms_mod, only: stdlog, stderr, error_mesg, FATAL, WARNING
use fms_mod, only: write_version_number, read_data, write_data, file_exist
use mosaic_mod, only: get_mosaic_ntiles, get_mosaic_ncontacts
use mpp_mod, only: mpp_npes, mpp_pe, mpp_root_pe, mpp_sum, mpp_min, mpp_max, NULL_PE
use mpp_mod, only: mpp_send, mpp_recv, mpp_sync_self, mpp_chksum, input_nml_file
use mpp_mod, only: mpp_clock_begin, mpp_clock_end, mpp_clock_id
use mpp_mod, only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_LOOP
use mpp_mod, only: COMM_TAG_1, COMM_TAG_2, COMM_TAG_3, COMM_TAG_4
use mpp_mod, only: COMM_TAG_5, COMM_TAG_6, COMM_TAG_7, COMM_TAG_8
use mpp_mod, only: COMM_TAG_9, COMM_TAG_10
use mpp_mod, only: mpp_gather
use fms_mod, only: clock_flag_default
use fms_io_mod, only: get_instance_filename
use mpp_domains_mod, only: domain2D, mpp_update_domains, mpp_define_domains
use mpp_parameter_mod, only: SCALAR_PAIR, CGRID_NE, BGRID_NE, CORNER, AGRID
use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain
use mpp_domains_mod, only: CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE
use mpp_domains_mod, only: mpp_get_neighbor_pe, NORTH, SOUTH, EAST, WEST
use mpp_domains_mod, only: mpp_define_io_domain
use time_manager_mod, only: time_type, get_date, get_time, set_date, operator(-)
use diag_manager_mod, only: register_diag_field, register_static_field, send_data
use diag_manager_mod, only: diag_axis_init


implicit none ; private

include 'netcdf.inc'

public icebergs_init, icebergs_end, icebergs_run, icebergs_stock_pe
public icebergs_incr_mass, icebergs_save_restart

type :: icebergs_gridded
  type(domain2D), pointer :: domain ! MPP domain
  integer :: halo ! Nominal halo width
  integer :: isc, iec, jsc, jec ! Indices of computational domain
  integer :: isd, ied, jsd, jed ! Indices of data domain
  integer :: my_pe, pe_N, pe_S, pe_E, pe_W ! MPI PE identifiers
  real, dimension(:,:), pointer :: lon=>null() ! Longitude of cell corners
  real, dimension(:,:), pointer :: lat=>null() ! Latitude of cell corners
  real, dimension(:,:), pointer :: lonc=>null() ! Longitude of cell centers
  real, dimension(:,:), pointer :: latc=>null() ! Latitude of cell centers
  real, dimension(:,:), pointer :: dx=>null() ! Length of cell edge (m)
  real, dimension(:,:), pointer :: dy=>null() ! Length of cell edge (m)
  real, dimension(:,:), pointer :: area=>null() ! Area of cell (m^2)
  real, dimension(:,:), pointer :: msk=>null() ! Ocean-land mask (1=ocean)
  real, dimension(:,:), pointer :: cos=>null() ! Cosine from rotation matrix to lat-lon coords
  real, dimension(:,:), pointer :: sin=>null() ! Sine from rotation matrix to lat-lon coords
  real, dimension(:,:), pointer :: uo=>null() ! Ocean zonal flow (m/s)
  real, dimension(:,:), pointer :: vo=>null() ! Ocean meridional flow (m/s)
  real, dimension(:,:), pointer :: ui=>null() ! Ice zonal flow (m/s)
  real, dimension(:,:), pointer :: vi=>null() ! Ice meridional flow (m/s)
  real, dimension(:,:), pointer :: ua=>null() ! Atmosphere zonal flow (m/s)
  real, dimension(:,:), pointer :: va=>null() ! Atmosphere meridional flow (m/s)
  real, dimension(:,:), pointer :: ssh=>null() ! Sea surface height (m)
  real, dimension(:,:), pointer :: sst=>null() ! Sea surface temperature (oC)
  real, dimension(:,:), pointer :: cn=>null() ! Sea-ice concentration (0 to 1)
  real, dimension(:,:), pointer :: hi=>null() ! Sea-ice thickness (m)
  real, dimension(:,:), pointer :: calving=>null() ! Calving mass rate [frozen runoff] (kg/s) (into stored ice)
  real, dimension(:,:), pointer :: calving_hflx=>null() ! Calving heat flux [heat content of calving] (W/m2) (into stored ice)
  real, dimension(:,:), pointer :: floating_melt=>null() ! Net melting rate to icebergs + bits (kg/s/m^2)
  real, dimension(:,:), pointer :: berg_melt=>null() ! Melting+erosion rate of icebergs (kg/s/m^2)
  real, dimension(:,:), pointer :: melt_buoy=>null() ! Buoyancy componenet of melting rate (kg/s/m^2)
  real, dimension(:,:), pointer :: melt_eros=>null() ! Erosion component of melting rate (kg/s/m^2)
  real, dimension(:,:), pointer :: melt_conv=>null() ! Convective component of melting rate (kg/s/m^2)
  real, dimension(:,:), pointer :: bergy_src=>null() ! Mass flux from berg erosion into bergy bits (kg/s/m^2)
  real, dimension(:,:), pointer :: bergy_melt=>null() ! Melting rate of bergy bits (kg/s/m^2)
  real, dimension(:,:), pointer :: bergy_mass=>null() ! Mass distribution of bergy bits (kg/s/m^2)
  real, dimension(:,:), pointer :: virtual_area=>null() ! Virtual surface coverage by icebergs (m^2)
  real, dimension(:,:), pointer :: mass=>null() ! Mass distribution (kg/m^2)
  real, dimension(:,:,:), pointer :: mass_on_ocean=>null() ! Mass distribution partitioned by neighbor (kg/m^2)
  real, dimension(:,:), pointer :: tmp=>null() ! Temporary work space
  real, dimension(:,:), pointer :: tmpc=>null() ! Temporary work space
  real, dimension(:,:,:), pointer :: stored_ice=>null() ! Accumulated ice mass flux at calving locations (kg)
  real, dimension(:,:), pointer :: stored_heat=>null() ! Heat content of stored ice (J)
  real, dimension(:,:,:), pointer :: real_calving=>null() ! Calving rate into iceberg class at calving locations (kg/s)
  real, dimension(:,:), pointer :: iceberg_heat_content=>null() ! Distributed heat content of bergs (J/m^2)
  real, dimension(:,:), pointer :: parity_x=>null() ! X component of vector point from i,j to i+1,j+1 (for detecting tri-polar fold)
  real, dimension(:,:), pointer :: parity_y=>null() ! Y component of vector point from i,j to i+1,j+1 (for detecting tri-polar fold)
  ! Diagnostics handles
  integer :: id_uo=-1, id_vo=-1, id_calving=-1, id_stored_ice=-1, id_accum=-1, id_unused=-1, id_floating_melt=-1
  integer :: id_melt_buoy=-1, id_melt_eros=-1, id_melt_conv=-1, id_virtual_area=-1, id_real_calving=-1
  integer :: id_calving_hflx_in=-1, id_stored_heat=-1, id_melt_hflx=-1, id_heat_content=-1
  integer :: id_mass=-1, id_ui=-1, id_vi=-1, id_ua=-1, id_va=-1, id_sst=-1, id_cn=-1, id_hi=-1
  integer :: id_bergy_src=-1, id_bergy_melt=-1, id_bergy_mass=-1, id_berg_melt=-1
end type icebergs_gridded

type :: xyt
  real :: lon, lat, day
  real :: mass, thickness, width, length, uvel, vvel
  real :: uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, cn, hi
  real :: mass_of_bits, heat_density
  integer :: year
  type(xyt), pointer :: next=>null()
end type xyt

type :: iceberg
  type(iceberg), pointer :: prev=>null(), next=>null()
  ! State variables (specific to the iceberg, needed for restarts)
  real :: lon, lat, uvel, vvel, mass, thickness, width, length
  real :: start_lon, start_lat, start_day, start_mass, mass_scaling
  real :: mass_of_bits, heat_density
  integer :: start_year
  integer :: ine, jne ! nearest index in NE direction (for convenience)
  real :: xi, yj ! Non-dimensional coords within current cell (0..1)
  ! Environment variables (as seen by the iceberg)
  real :: uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, cn, hi
  type(xyt), pointer :: trajectory=>null()
end type iceberg

type :: buffer
  integer :: size=0
  real, dimension(:,:), pointer :: data
end type buffer

type, public :: icebergs ; private
  type(icebergs_gridded), pointer :: grd
  type(iceberg), pointer :: first=>null()
  type(xyt), pointer :: trajectories=>null()
  real :: dt           ! Time-step between iceberg calls (should make adaptive?)
  integer :: current_year
  real :: current_yearday ! 1.00-365.99
  integer :: traj_sample_hrs
  integer :: verbose_hrs
  integer :: clock, clock_mom, clock_the, clock_int, clock_cal, clock_com, clock_ini, clock_ior, clock_iow, clock_dia ! ids for fms timers
  real :: rho_bergs ! Density of icebergs
  real :: LoW_ratio ! Initial ratio L/W for newly calved icebergs
  real :: bergy_bit_erosion_fraction ! Fraction of erosion melt flux to divert to bergy bits
  real :: sicn_shift ! Shift of sea-ice concentration in erosion flux modulation (0<sicn_shift<1)
  real, dimension(:), pointer :: initial_mass, distribution, mass_scaling
  real, dimension(:), pointer :: initial_thickness, initial_width, initial_length
  logical :: restarted=.false. ! Indicate whether we read state from a restart or not
  logical :: use_operator_splitting=.true. ! Use first order operator splitting for thermodynamics
  logical :: add_weight_to_ocean=.true. ! Add weight of bergs to ocean
  logical :: passive_mode=.false. ! Add weight of icebergs + bits to ocean
  logical :: time_average_weight=.false. ! Time average the weight on the ocean
  real :: speed_limit=0. ! CFL speed limit for a berg
  type(buffer), pointer :: obuffer_n=>null(), ibuffer_n=>null()
  type(buffer), pointer :: obuffer_s=>null(), ibuffer_s=>null()
  type(buffer), pointer :: obuffer_e=>null(), ibuffer_e=>null()
  type(buffer), pointer :: obuffer_w=>null(), ibuffer_w=>null()
  ! Budgets
  real :: net_calving_received=0., net_calving_returned=0.
  real :: net_incoming_calving=0., net_outgoing_calving=0.
  real :: net_incoming_calving_heat=0., net_outgoing_calving_heat=0.
  real :: net_incoming_calving_heat_used=0., net_heat_to_bergs=0.
  real :: stored_start=0., stored_end=0.
  real :: stored_heat_start=0., stored_heat_end=0., net_heat_to_ocean=0.
  real :: net_calving_used=0., net_calving_to_bergs=0.
  real :: floating_mass_start=0., floating_mass_end=0.
  real :: floating_heat_start=0., floating_heat_end=0.
  real :: icebergs_mass_start=0., icebergs_mass_end=0.
  real :: bergy_mass_start=0., bergy_mass_end=0.
  real :: returned_mass_on_ocean=0.
  real :: net_melt=0., berg_melt=0., bergy_src=0., bergy_melt=0.
  integer :: nbergs_calved=0, nbergs_melted=0, nbergs_start=0, nbergs_end=0
  integer :: nspeeding_tickets=0
  integer, dimension(:), pointer :: nbergs_calved_by_class=>null()
end type icebergs

! Global constants
character(len=*), parameter :: version = '$Id: ice_bergs.F90,v 20.0 2013/12/13 23:28:21 fms Exp $'
character(len=*), parameter :: tagname = '$Name: tikal $'

integer, parameter :: nclasses=10 ! Number of ice bergs classes
integer, parameter :: file_format_major_version=0
integer, parameter :: file_format_minor_version=1
integer, parameter :: delta_buf=25 ! Size by which to increment buffers
real, parameter :: pi_180=pi/180. ! Converts degrees to radians
real, parameter :: rho_ice=916.7 ! Density of fresh ice @ 0oC (kg/m^3)
real, parameter :: rho_water=999.8 ! Density of fresh water @ 0oC (kg/m^3)
real, parameter :: rho_air=1.1 ! Density of air @ 0oC (kg/m^3) ???
real, parameter :: rho_seawater=1025. ! Approx. density of surface sea water @ 0oC (kg/m^3)
real, parameter :: gravity=9.8 ! Gravitational acceleratio (m/s^2)
real, parameter :: Cd_av=1.3 ! (Vertical) Drag coefficient between bergs and atmos (?)
real, parameter :: Cd_ah=0.0055 ! (Horizontal) Drag coefficient between bergs and atmos (?)
real, parameter :: Cd_wv=0.9 ! (Vertical) Drag coefficient between bergs and ocean (?)
real, parameter :: Cd_wh=0.0012 ! (Horizontal) Drag coefficient between bergs and ocean (?)
real, parameter :: Cd_iv=0.9 ! (Vertical) Drag coefficient between bergs and sea-ice (?)
!TOM> no horizontal drag for sea ice! real, parameter :: Cd_ih=0.0012 ! (Horizontal) Drag coefficient between bergs and sea-ice (?)

! Global data (minimal for debugging)
logical :: verbose=.false. ! Be verbose to stderr
logical :: budget=.true. ! Calculate budgets
logical :: debug=.false. ! Turn on debugging
logical :: really_debug=.false. ! Turn on debugging
logical :: parallel_reprod=.true. ! Reproduce across different PE decompositions
logical :: use_slow_find=.true. ! Use really slow (but robust) find_cell for reading restarts
logical :: ignore_ij_restart=.false. ! Read i,j location from restart if available (needed to use restarts on different grids)
logical :: generate_test_icebergs=.false. ! Create icebergs in absence of a restart file
logical :: use_roundoff_fix=.true. ! Use a "fix" for the round-off discrepancy between is_point_in_cell() and pos_within_cell()
logical :: old_bug_rotated_weights=.false. ! Skip the rotation of off-center weights for rotated halo updates
logical :: make_calving_reproduce=.false. ! Make the calving.res.nc file reproduce across pe count changes.
character(len=10) :: restart_input_dir = 'INPUT/'

logical :: folded_north_on_pe = .false.

contains

! ##############################################################################

subroutine accel(bergs, berg, i, j, xi, yj, lat, uvel, vvel, uvel0, vvel0, dt, ax, ay, debug_flag)
! Arguments
type(icebergs), pointer :: bergs
type(iceberg), pointer :: berg
integer, intent(in) :: i, j
real, intent(in) :: xi, yj, lat, uvel, vvel, uvel0, vvel0, dt
real, intent(inout) :: ax, ay
logical, optional :: debug_flag
! Local variables
type(icebergs_gridded), pointer :: grd
real :: uo, vo, ui, vi, ua, va, uwave, vwave, ssh_x, ssh_y, sst, cn, hi
real :: f_cori, T, D, W, L, M, F
real :: drag_ocn, drag_atm, drag_ice, wave_rad
real :: c_ocn, c_atm, c_ice
real :: ampl, wmod, Cr, Lwavelength, Lcutoff, Ltop
real, parameter :: alpha=0.0, beta=1.0, accel_lim=1.e-2, Cr0=0.06, vel_lim=15.
real :: lambda, detA, A11, A12, axe, aye, D_hi
real :: uveln, vveln, us, vs, speed, loc_dx, new_speed
logical :: dumpit
integer :: itloop
integer :: stderrunit

  ! Get the stderr unit number.
  stderrunit = stderr()

  ! For convenience
  grd=>bergs%grd

  ! Interpolate gridded fields to berg
  call interp_flds(grd, i, j, xi, yj, uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, cn, hi)

  f_cori=(2.*omega)*sin(pi_180*lat)

  M=berg%mass
  T=berg%thickness ! total thickness
  D=(bergs%rho_bergs/rho_seawater)*T ! draught (keel depth)
  F=T-D ! freeboard
  W=berg%width
  L=berg%length

  hi=min(hi,D)
  D_hi=max(0.,D-hi)

  ! Wave radiation
  uwave=ua-uo; vwave=va-vo  ! Use wind speed rel. to ocean for wave model (aja)?
  wmod=uwave*uwave+vwave*vwave ! The wave amplitude and length depend on the wind speed relative to the ocean current;
                               !  actually wmod is wmod**2 here.
  ampl=0.5*0.02025*wmod ! This is "a", the wave amplitude
  Lwavelength=0.32*wmod ! Surface wave length fitted to data in table at
  !      http://www4.ncsu.edu/eos/users/c/ceknowle/public/chapter10/part2.html
  Lcutoff=0.125*Lwavelength
  Ltop=0.25*Lwavelength
  Cr=Cr0*min(max(0.,(L-Lcutoff)/((Ltop-Lcutoff)+1.e-30)),1.) ! Wave radiation coefficient
  !     fitted to graph from Carrieres et al.,  POAC Drift Model.
  wave_rad=0.5*rho_seawater/M*Cr*gravity*ampl*min(ampl,F)*(2.*W*L)/(W+L)
  wmod = sqrt(ua*ua+va*va) ! Wind speed
  if (wmod.ne.0.) then
    uwave=ua/wmod ! Wave radiation force acts in wind direction ...
    vwave=va/wmod
  else
    uwave=0.; vwave=0.; wave_rad=0. ! ... and only when wind is present.
  endif

  ! Weighted drag coefficients
  c_ocn=rho_seawater/M*(0.5*Cd_wv*W*(D_hi)+Cd_wh*W*L)
  c_atm=rho_air     /M*(0.5*Cd_av*W*F     +Cd_ah*W*L)
  c_ice=rho_ice     /M*(0.5*Cd_iv*W*hi              )
  if (abs(ui)+abs(vi).eq.0.) c_ice=0.

  uveln=uvel; vveln=vvel ! Copy starting uvel, vvel
  do itloop=1,2 ! Iterate on drag coefficients

    us=0.5*(uveln+uvel); vs=0.5*(vveln+vvel)
    drag_ocn=c_ocn*sqrt( (us-uo)**2+(vs-vo)**2 )
    drag_atm=c_atm*sqrt( (us-ua)**2+(vs-va)**2 )
    drag_ice=c_ice*sqrt( (us-ui)**2+(vs-vi)**2 )

    ! Explicit accelerations
   !axe= f_cori*vvel -gravity*ssh_x +wave_rad*uwave &
   !    -drag_ocn*(uvel-uo) -drag_atm*(uvel-ua) -drag_ice*(uvel-ui)
   !aye=-f_cori*uvel -gravity*ssh_y +wave_rad*vwave &
   !    -drag_ocn*(vvel-vo) -drag_atm*(vvel-va) -drag_ice*(vvel-vi)
    axe=-gravity*ssh_x +wave_rad*uwave
    aye=-gravity*ssh_y +wave_rad*vwave
    if (alpha>0.) then ! If implicit, use time-level (n) rather than RK4 latest
      axe=axe+f_cori*vvel0
      aye=aye-f_cori*uvel0
    else
      axe=axe+f_cori*vvel
      aye=aye-f_cori*uvel
    endif
    if (beta>0.) then ! If implicit, use time-level (n) rather than RK4 latest
      axe=axe-drag_ocn*(uvel0-uo) -drag_atm*(uvel0-ua) -drag_ice*(uvel0-ui)
      aye=aye-drag_ocn*(vvel0-vo) -drag_atm*(vvel0-va) -drag_ice*(vvel0-vi)
    else
      axe=axe-drag_ocn*(uvel-uo) -drag_atm*(uvel-ua) -drag_ice*(uvel-ui)
      aye=aye-drag_ocn*(vvel-vo) -drag_atm*(vvel-va) -drag_ice*(vvel-vi)
    endif

    ! Solve for implicit accelerations
    if (alpha+beta.gt.0.) then
      lambda=drag_ocn+drag_atm+drag_ice
      A11=1.+beta*dt*lambda
      A12=alpha*dt*f_cori
      detA=1./(A11**2+A12**2)
      ax=detA*(A11*axe+A12*aye)
      ay=detA*(A11*aye-A12*axe)
    else
      ax=axe; ay=aye
    endif

    uveln=uvel0+dt*ax
    vveln=vvel0+dt*ay

  enddo ! itloop

  ! Limit speed of bergs based on a CFL criteria
  if (bergs%speed_limit>0.) then
    speed=sqrt(uveln*uveln+vveln*vveln) ! Speed of berg
    if (speed>0.) then
      loc_dx=min(0.5*(grd%dx(i,j)+grd%dx(i,j-1)),0.5*(grd%dy(i,j)+grd%dy(i-1,j))) ! min(dx,dy)
     !new_speed=min(loc_dx/dt*bergs%speed_limit,speed) ! Restrict speed to dx/dt x factor
      new_speed=loc_dx/dt*bergs%speed_limit ! Speed limit as a factor of dx / dt
      if (new_speed<speed) then
        uveln=uveln*(new_speed/speed) ! Scale velocity to reduce speed
        vveln=vveln*(new_speed/speed) ! without changing the direction
        bergs%nspeeding_tickets=bergs%nspeeding_tickets+1
      endif
    endif
  endif

  dumpit=.false.
  if (abs(uveln)>vel_lim.or.abs(vveln)>vel_lim) then
    dumpit=.true.
    write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Dump triggered by excessive velocity'
  endif
  if (abs(ax)>accel_lim.or.abs(ay)>accel_lim) then
    dumpit=.true.
    write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Dump triggered by excessive acceleration'
  endif
  if (present(debug_flag)) then
    if (debug_flag) dumpit=.true.
    write(stderrunit,'("pe=",i3,x,a)') mpp_pe(),'Debug dump flagged by arguments'
  endif
  if (dumpit) then
 100 format('pe=',i3,a15,9(x,a8,es12.3))
 200 format('pe=',i3,a15,(x,a8,i12),9(x,a8,es12.3))
    write(stderrunit,200) mpp_pe(),'Starting pars:', &
      'yr0=',berg%start_year, 'day0=',berg%start_day, &
      'lon0=',berg%start_lon, 'lat0=',berg%start_lat, 'mass0=',berg%start_mass, &
      'sclng=',berg%mass_scaling
    write(stderrunit,100) mpp_pe(),'Geometry:', &
      'M=',M, 'T=',T, 'D=',D, 'F=',F, 'W=',W, 'L=',L
    write(stderrunit,100) mpp_pe(),'delta U:', &
      'u(n)=',uvel0, 'u(*)=', uvel, 'u(n+1)=',uvel+dt*ax, 'del u=',dt*ax
    write(stderrunit,100) mpp_pe(),'U terms', &
      'f*v=',f_cori*vvel, &
      'g*H_x=',-gravity*ssh_x, &
      'wave*ua=',wave_rad*uwave, &
      'd*(u-uo)=',-drag_ocn*(uvel-uo), &
      'd*(u-ua)=',-drag_atm*(uvel-ua), &
      'd*(u-ui)=',-drag_ice*(uvel-ui)
    write(stderrunit,100) mpp_pe(),'U accel.', &
      'axe=',axe, &
      'ax=',ax, &
      'ax(cori)=',detA*(A11*(f_cori*vvel)+A12*(-f_cori*uvel)), &
      'ax(grav)=',detA*(A11*(-gravity*ssh_x)+A12*(-gravity*ssh_y)), &
      'ax(wave)=',detA*(A11*(wave_rad*uwave)+A12*(wave_rad*vwave)), &
      'ax(ocn)=',detA*(A11*(-drag_ocn*(uvel-uo))+A12*(-drag_ocn*(vvel-vo))), &
      'ax(atm)=',detA*(A11*(-drag_atm*(uvel-ua))+A12*(-drag_atm*(vvel-va))), &
      'ax(ice)=',detA*(A11*(-drag_ice*(uvel-ui))+A12*(-drag_ice*(vvel-vi)))
    write(stderrunit,100) mpp_pe(),'delta V:', &
      'v(n)=',vvel0, 'v(*)=', vvel, 'v(n+1)=',vvel+dt*ay, 'del v=',dt*ay
    write(stderrunit,100) mpp_pe(),'V terms', &
      'f*u=',-f_cori*uvel, &
      'g*H_y=',-gravity*ssh_y, &
      'wave*va=',wave_rad*vwave, &
      'd*(v-vo)=',-drag_ocn*(vvel-vo), &
      'd*(v-va)=',-drag_atm*(vvel-va), &
      'd*(v-vi)=',-drag_ice*(vvel-vi)
    write(stderrunit,100) mpp_pe(),'V accel. pe=', &
      'aye=',aye, &
      'ay=',ay, &
      'ay(cori)=',detA*(-A12*(f_cori*vvel)+A11*(-f_cori*uvel)), &
      'ay(grav)=',detA*(-A12*(-gravity*ssh_x)+A11*(-gravity*ssh_y)), &
      'ay(wave)=',detA*(-A12*(wave_rad*uwave)+A11*(wave_rad*vwave)), &
      'ay(ocn)=',detA*(-A12*(-drag_ocn*(uvel-uo))+A11*(-drag_ocn*(vvel-vo))), &
      'ay(atm)=',detA*(-A12*(-drag_atm*(uvel-ua))+A11*(-drag_atm*(vvel-va))), &
      'ay(ice)=',detA*(-A12*(-drag_ice*(uvel-ui))+A11*(-drag_ice*(vvel-vi)))
    write(stderrunit,100) mpp_pe(),'Vel scales', &
      '|va-vo|=',sqrt((ua-uo)**2+(va-vo)**2), &
      '|vo-vb|=',sqrt((uvel-uo)**2+(vvel-vo)**2), &
      '|va-vb|=',sqrt((uvel-ua)**2+(vvel-va)**2), &
      '|vi-vb|=',sqrt((uvel-ui)**2+(vvel-vi)**2), &
      '|vb|=',sqrt((uvel)**2+(vvel)**2), &
      '|va|=',sqrt((ua)**2+(va)**2), &
      '|vo|=',sqrt((uo)**2+(vo)**2), &
      '|vi|=',sqrt((ui)**2+(vi)**2)
    write(stderrunit,100) mpp_pe(),'Time scales', &
      'f=',f_cori, 'wave_rad=',wave_rad, 'do=',drag_ocn, 'da=',drag_atm, 'di=',drag_ice
    write(stderrunit,100) mpp_pe(),'u*', &
      'd*=',lambda, &
      'u*=',(drag_ocn*uo+drag_atm*ua+drag_ice*ui)/lambda, &
      'uo*=',(drag_ocn*uo)/lambda, &
      'ua*=',(drag_atm*ua)/lambda, &
      'ui*=',(drag_ice*ui)/lambda
    write(stderrunit,100) mpp_pe(),'v*', &
      'd*=',lambda, &
      'v*=',(drag_ocn*vo+drag_atm*va+drag_ice*vi)/lambda, &
      'vo*=',(drag_ocn*vo)/lambda, &
      'va*=',(drag_atm*va)/lambda, &
      'vi*=',(drag_ice*vi)/lambda
    write(stderrunit,100) mpp_pe(),'params', &
      'a=',ampl, 'Lwl=',Lwavelength, 'Lcut=',Lcutoff, 'Ltop=',Ltop, 'hi=',hi, 'Cr=',Cr
    write(stderrunit,100) mpp_pe(),'Position', &
      'xi=',xi, 'yj=',yj, 'lat=',lat
    call dump_locfld(grd,i,j,grd%msk,'MSK')
    call dump_locfld(grd,i,j,grd%ssh,'SSH')
    call dump_locfld(grd,i,j,grd%sst,'SST')
    call dump_locvel(grd,i,j,grd%uo,'Uo')
    call dump_locvel(grd,i,j,grd%vo,'Vo')
    call dump_locvel(grd,i,j,grd%ua,'Ua')
    call dump_locvel(grd,i,j,grd%va,'Va')
    call dump_locvel(grd,i,j,grd%ui,'Ui')
    call dump_locvel(grd,i,j,grd%vi,'Vi')
    call dump_locfld(grd,i,j,grd%hi,'HI')
    call dump_locfld(grd,i,j,grd%cn,'CN')
    call dump_locvel(grd,i,j,grd%lon,'Lon')
    call dump_locvel(grd,i,j,grd%lat,'Lat')
    call print_berg(stderrunit,berg,'diamonds, accel, large accel')
  endif

  contains

  subroutine dump_locfld(grd,i0,j0,A,lbl)
  ! Arguments
  type(icebergs_gridded), pointer :: grd
  integer, intent(in) :: i0, j0
  real, dimension(grd%isd:grd%ied,grd%jsd:grd%jed), intent(in) :: A
  character(len=*) :: lbl
! Local variables
  integer :: i, j, ii, jj
  real :: B(-1:1,-1:1), fac

  do jj=-1,1
    j=max(grd%jsd,min(grd%jed,jj+j0))
    do ii=-1,1
      i=max(grd%isd,min(grd%ied,ii+i0))
      B(ii,jj)=A(i,j)
      if ((i.ne.ii+i0).or.(j.ne.jj+j0)) B(ii,jj)=-9.999999e-99
    enddo
  enddo
  write(stderrunit,'("pe=",i3,x,a8,3i12)') mpp_pe(),lbl,(i0+ii,ii=-1,1)
  do jj=1,-1,-1
    write(stderrunit,'("pe=",i3,x,i8,3es12.4)') mpp_pe(),j0+jj,(B(ii,jj),ii=-1,1)
  enddo
  end subroutine dump_locfld

  subroutine dump_locvel(grd,i0,j0,A,lbl)
  ! Arguments
  type(icebergs_gridded), pointer :: grd
  integer, intent(in) :: i0, j0
  real, dimension(grd%isd:grd%ied,grd%jsd:grd%jed), intent(in) :: A
  character(len=*) :: lbl
! Local variables
  integer :: i, j, ii, jj
  real :: B(-1:0,-1:0), fac

  do jj=-1,0
    j=max(grd%jsd,min(grd%jed,jj+j0))
    do ii=-1,0
      i=max(grd%isd,min(grd%ied,ii+i0))
      B(ii,jj)=A(i,j)
      if ((i.ne.ii+i0).or.(j.ne.jj+j0)) B(ii,jj)=-9.999999e-99
    enddo
  enddo
  write(stderrunit,'("pe=",i3,x,a8,3i12)') mpp_pe(),lbl,(i0+ii,ii=-1,0)
  do jj=0,-1,-1
    write(stderrunit,'("pe=",i3,x,i8,3es12.4)') mpp_pe(),j0+jj,(B(ii,jj),ii=-1,0)
  enddo
  end subroutine dump_locvel

end subroutine accel

! ##############################################################################

subroutine thermodynamics(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(icebergs_gridded), pointer :: grd
real :: M, T, W, L, SST, Vol, Ln, Wn, Tn, nVol, IC, Dn
real :: Mv, Me, Mb, melt, dvo, dva, dM, Ss, dMe, dMb, dMv
real :: Mnew, Mnew1, Mnew2
real :: Mbits, nMbits, dMbitsE, dMbitsM, Lbits, Abits, Mbb
integer :: i,j, stderrunit
type(iceberg), pointer :: this, next
real, parameter :: perday=1./86400.

  ! For convenience
  grd=>bergs%grd

  this=>bergs%first
  do while(associated(this))
    if (debug) call check_position(grd, this, 'thermodynamics (top)')

    call interp_flds(grd, this%ine, this%jne, this%xi, this%yj, this%uo, this%vo, &
            this%ui, this%vi, this%ua, this%va, this%ssh_x, this%ssh_y, this%sst, &
            this%cn, this%hi)
    SST=this%sst
    IC=min(1.,this%cn+bergs%sicn_shift) ! Shift sea-ice concentration
    M=this%mass
    T=this%thickness ! total thickness
  ! D=(bergs%rho_bergs/rho_seawater)*T ! draught (keel depth)
  ! F=T-D ! freeboard
    W=this%width
    L=this%length
    i=this%ine
    j=this%jne
    Vol=T*W*L

    ! Environment
    dvo=sqrt((this%uvel-this%uo)**2+(this%vvel-this%vo)**2)
    dva=sqrt((this%ua-this%uo)**2+(this%va-this%vo)**2)
    Ss=1.5*(dva**0.5)+0.1*dva ! Sea state

    ! Melt rates in m/s
    Mv=max( 7.62e-3*SST+1.29e-3*(SST**2), 0.) &! Buoyant convection at sides
        *perday ! convert to m/s
    Mb=max( 0.58*(dvo**0.8)*(SST+4.0)/(L**0.2), 0.) &! Basal turbulent melting
        *perday ! convert to m/s
    Me=max( 1./12.*(SST+2.)*Ss*(1+cos(pi*(IC**3))) ,0.) &! Wave erosion
        *perday ! convert to m/s

    if (bergs%use_operator_splitting) then
      ! Operator split update of volume/mass
      Tn=max(T-Mb*bergs%dt,0.) ! new total thickness (m)
      nVol=Tn*W*L ! new volume (m^3)
      Mnew1=(nVol/Vol)*M ! new mass (kg)
      dMb=M-Mnew1 ! mass lost to basal melting (>0) (kg)

      Ln=max(L-Mv*bergs%dt,0.) ! new length (m)
      Wn=max(W-Mv*bergs%dt,0.) ! new width (m)
      nVol=Tn*Wn*Ln ! new volume (m^3)
      Mnew2=(nVol/Vol)*M ! new mass (kg)
      dMv=Mnew1-Mnew2 ! mass lost to buoyant convection (>0) (kg)

      Ln=max(Ln-Me*bergs%dt,0.) ! new length (m)
      Wn=max(Wn-Me*bergs%dt,0.) ! new width (m)
      nVol=Tn*Wn*Ln ! new volume (m^3)
      Mnew=(nVol/Vol)*M ! new mass (kg)
      dMe=Mnew2-Mnew ! mass lost to erosion (>0) (kg)
      dM=M-Mnew ! mass lost to all erosion and melting (>0) (kg)
    else
      ! Update dimensions of berg
      Ln=max(L-(Mv+Me)*(bergs%dt),0.) ! (m)
      Wn=max(W-(Mv+Me)*(bergs%dt),0.) ! (m)
      Tn=max(T-Mb*(bergs%dt),0.) ! (m)
      ! Update volume and mass of berg
      nVol=Tn*Wn*Ln ! (m^3)
      Mnew=(nVol/Vol)*M ! (kg)
      dM=M-Mnew ! (kg)
      dMb=(M/Vol)*(W*L)*Mb*bergs%dt ! approx. mass loss to basal melting (kg)
      dMe=(M/Vol)*(T*(W+L))*Me*bergs%dt ! approx. mass lost to erosion (kg)
      dMv=(M/Vol)*(T*(W+L))*Mv*bergs%dt ! approx. mass loss to buoyant convection (kg)
    endif

    ! Bergy bits
    if (bergs%bergy_bit_erosion_fraction>0.) then
      Mbits=this%mass_of_bits ! mass of bergy bits (kg)
      dMbitsE=bergs%bergy_bit_erosion_fraction*dMe ! change in mass of bits (kg)
      nMbits=Mbits+dMbitsE ! add new bergy bits to mass (kg)
      Lbits=min(L,W,T,40.) ! assume bergy bits are smallest dimension or 40 meters
      Abits=(Mbits/bergs%rho_bergs)/Lbits ! Effective bottom area (assuming T=Lbits)
      Mbb=max( 0.58*(dvo**0.8)*(SST+2.0)/(Lbits**0.2), 0.) &! Basal turbulent melting (for bits)
           *perday ! convert to m/s
      Mbb=bergs%rho_bergs*Abits*Mbb ! in kg/s
      dMbitsM=min(Mbb*bergs%dt,nMbits) ! bergy bits mass lost to melting (kg)
      nMbits=nMbits-dMbitsM ! remove mass lost to bergy bits melt
      if (Mnew==0.) then ! if parent berg has completely melted then
        dMbitsM=dMbitsM+nMbits ! instantly melt all the bergy bits
        nMbits=0.
      endif
    else
      Abits=0.
      dMbitsE=0.
      dMbitsM=0.
      nMbits=this%mass_of_bits ! retain previous value incase non-zero
    endif

    ! Add melting to the grid and field diagnostics
    if (grd%area(i,j).ne.0.) then
      melt=(dM-(dMbitsE-dMbitsM))/bergs%dt ! kg/s
      grd%floating_melt(i,j)=grd%floating_melt(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
      melt=melt*this%heat_density ! kg/s x J/kg = J/s
      grd%calving_hflx(i,j)=grd%calving_hflx(i,j)+melt/grd%area(i,j)*this%mass_scaling ! W/m2
      bergs%net_heat_to_ocean=bergs%net_heat_to_ocean+melt*this%mass_scaling*bergs%dt ! J
      melt=dM/bergs%dt ! kg/s
      grd%berg_melt(i,j)=grd%berg_melt(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
      melt=dMbitsE/bergs%dt ! mass flux into bergy bits in kg/s
      grd%bergy_src(i,j)=grd%bergy_src(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
      melt=dMbitsM/bergs%dt ! melt rate of bergy bits in kg/s
      grd%bergy_melt(i,j)=grd%bergy_melt(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
      if(grd%id_melt_buoy>0) then
        melt=dMb/bergs%dt ! melt rate due to buoyancy term in kg/s
        grd%melt_buoy(i,j)=grd%melt_buoy(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
      endif
      if(grd%id_melt_eros>0) then
        melt=dMe/bergs%dt ! erosion rate in kg/s
        grd%melt_eros(i,j)=grd%melt_eros(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
      endif
      if(grd%id_melt_conv>0) then
        melt=dMv/bergs%dt ! melt rate due to convection term in kg/s
        grd%melt_conv(i,j)=grd%melt_conv(i,j)+melt/grd%area(i,j)*this%mass_scaling ! kg/m2/s
      endif
    else
      stderrunit = stderr()
      write(stderrunit,*) 'diamonds, thermodynamics: berg appears to have grounded!!!! PE=',mpp_pe(),i,j
      call print_berg(stderrunit,this,'thermodynamics, grounded')
      if (associated(this%trajectory)) &
        write(stderrunit,*) 'traj=',this%trajectory%lon,this%trajectory%lat
      write(stderrunit,*) 'msk=',grd%msk(i,j),grd%area(i,j)
      call error_mesg('diamonds, thermodynamics', 'berg appears to have grounded!', FATAL)
    endif

    ! Rolling
    Dn=(bergs%rho_bergs/rho_seawater)*Tn ! draught (keel depth)
    if ( Dn>0. ) then
       if ( max(Wn,Ln)<sqrt(0.92*(Dn**2)+58.32*Dn) ) then
          T=Tn
          Tn=Wn
          Wn=T
      end if
    endif

    ! Store the new state of iceberg (with L>W)
    this%mass=Mnew
    this%mass_of_bits=nMbits
    this%thickness=Tn
    this%width=min(Wn,Ln)
    this%length=max(Wn,Ln)

    next=>this%next

    ! Did berg completely melt?
    if (Mnew<=0.) then ! Delete the berg
      call move_trajectory(bergs, this)
      call delete_iceberg_from_list(bergs%first, this)
      bergs%nbergs_melted=bergs%nbergs_melted+1
    else ! Diagnose mass distribution on grid
      if (grd%id_virtual_area>0)&
           & grd%virtual_area(i,j)=grd%virtual_area(i,j)+(Wn*Ln+Abits)*this%mass_scaling ! m^2
      if (grd%id_mass>0 .or. bergs%add_weight_to_ocean)&
           & grd%mass(i,j)=grd%mass(i,j)+Mnew/grd%area(i,j)*this%mass_scaling ! kg/m2
      if (grd%id_bergy_mass>0 .or. bergs%add_weight_to_ocean)&
           & grd%bergy_mass(i,j)=grd%bergy_mass(i,j)+nMbits/grd%area(i,j)*this%mass_scaling ! kg/m2
      if (bergs%add_weight_to_ocean .and. .not. bergs%time_average_weight) &
         call spread_mass_across_ocean_cells(grd, i, j, this%xi, this%yj, Mnew, nMbits, this%mass_scaling)
    endif

    this=>next
  enddo

end subroutine thermodynamics

! ##############################################################################

subroutine spread_mass_across_ocean_cells(grd, i, j, x, y, Mberg, Mbits, scaling)
  ! Arguments
  type(icebergs_gridded), pointer :: grd
  integer, intent(in) :: i, j
  real, intent(in) :: x, y, Mberg, Mbits, scaling
  ! Local variables
  real :: xL, xC, xR, yD, yC, yU, Mass
  real :: yDxL, yDxC, yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR

  Mass=(Mberg+Mbits)*scaling
  xL=min(0.5, max(0., 0.5-x))
  xR=min(0.5, max(0., x-0.5))
  xC=max(0., 1.-(xL+xR))
  yD=min(0.5, max(0., 0.5-y))
  yU=min(0.5, max(0., y-0.5))
  yC=max(0., 1.-(yD+yU))

  yDxL=yD*xL*grd%msk(i-1,j-1)
  yDxC=yD*xC*grd%msk(i  ,j-1)
  yDxR=yD*xR*grd%msk(i+1,j-1)
  yCxL=yC*xL*grd%msk(i-1,j  )
  yCxR=yC*xR*grd%msk(i+1,j  )
  yUxL=yU*xL*grd%msk(i-1,j+1)
  yUxC=yU*xC*grd%msk(i  ,j+1)
  yUxR=yU*xR*grd%msk(i+1,j+1)
  yCxC=1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )

  grd%mass_on_ocean(i,j,1)=grd%mass_on_ocean(i,j,1)+yDxL*Mass
  grd%mass_on_ocean(i,j,2)=grd%mass_on_ocean(i,j,2)+yDxC*Mass
  grd%mass_on_ocean(i,j,3)=grd%mass_on_ocean(i,j,3)+yDxR*Mass
  grd%mass_on_ocean(i,j,4)=grd%mass_on_ocean(i,j,4)+yCxL*Mass
  grd%mass_on_ocean(i,j,5)=grd%mass_on_ocean(i,j,5)+yCxC*Mass
  grd%mass_on_ocean(i,j,6)=grd%mass_on_ocean(i,j,6)+yCxR*Mass
  grd%mass_on_ocean(i,j,7)=grd%mass_on_ocean(i,j,7)+yUxL*Mass
  grd%mass_on_ocean(i,j,8)=grd%mass_on_ocean(i,j,8)+yUxC*Mass
  grd%mass_on_ocean(i,j,9)=grd%mass_on_ocean(i,j,9)+yUxR*Mass

end subroutine spread_mass_across_ocean_cells

! ##############################################################################

subroutine interp_flds(grd, i, j, xi, yj, uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, cn, hi)
! Arguments
type(icebergs_gridded), pointer :: grd
integer, intent(in) :: i, j
real, intent(in) :: xi, yj
real, intent(out) :: uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, cn, hi
! Local variables
real :: cos_rot, sin_rot
#ifdef USE_OLD_SSH_GRADIENT
real :: dxm, dx0, dxp
#endif
real :: hxm, hxp
real, parameter :: ssh_coast=0.00

  cos_rot=bilin(grd, grd%cos, i, j, xi, yj)
  sin_rot=bilin(grd, grd%sin, i, j, xi, yj)

  uo=bilin(grd, grd%uo, i, j, xi, yj)
  vo=bilin(grd, grd%vo, i, j, xi, yj)
  ui=bilin(grd, grd%ui, i, j, xi, yj)
  vi=bilin(grd, grd%vi, i, j, xi, yj)
  ua=bilin(grd, grd%ua, i, j, xi, yj)
  va=bilin(grd, grd%va, i, j, xi, yj)
  ! These fields are cell centered (A-grid) and would
  ! best be interpolated using PLM. For now we use PCM!
  sst=grd%sst(i,j) ! A-grid
  cn=grd%cn(i,j) ! A-grid
  hi=grd%hi(i,j) ! A-grid

  ! Estimate SSH gradient in X direction
#ifdef USE_OLD_SSH_GRADIENT
  dxp=0.5*(grd%dx(i+1,j)+grd%dx(i+1,j-1))
  dx0=0.5*(grd%dx(i,j)+grd%dx(i,j-1))
  dxm=0.5*(grd%dx(i-1,j)+grd%dx(i-1,j-1))
  hxm=2.*(grd%ssh(i,j)-grd%ssh(i-1,j))/(dx0+dxm)*grd%msk(i-1,j) &
        +(-ssh_coast)/(dx0+dxm)*(1.-grd%msk(i-1,j)) ! force to drive bergs away from coasts
  hxp=2.*(grd%ssh(i+1,j)-grd%ssh(i,j))/(dx0+dxp)*grd%msk(i+1,j) &
        +(+ssh_coast)/(dx0+dxp)*(1.-grd%msk(i+1,j)) ! force to drive bergs away from coasts
#else
  if (yj>=0.5) then
    hxp=(yj-0.5)*ddx_ssh(grd,i  ,j+1)+(1.5-yj)*ddx_ssh(grd,i  ,j  )
    hxm=(yj-0.5)*ddx_ssh(grd,i-1,j+1)+(1.5-yj)*ddx_ssh(grd,i-1,j  )
  else
    hxp=(yj+0.5)*ddx_ssh(grd,i  ,j  )+(0.5-yj)*ddx_ssh(grd,i  ,j-1)
    hxm=(yj+0.5)*ddx_ssh(grd,i-1,j  )+(0.5-yj)*ddx_ssh(grd,i-1,j-1)
  endif
#endif
  ! ssh_x is at the u-point on a C-grid
  ssh_x=xi*hxp+(1.-xi)*hxm

  ! Estimate SSH gradient in Y direction
#ifdef USE_OLD_SSH_GRADIENT
  dxp=0.5*(grd%dy(i,j+1)+grd%dy(i-1,j+1))
  dx0=0.5*(grd%dy(i,j)+grd%dy(i-1,j))
  dxm=0.5*(grd%dy(i,j-1)+grd%dy(i-1,j-1))
  hxm=2.*(grd%ssh(i,j)-grd%ssh(i,j-1))/(dx0+dxm)*grd%msk(i,j-1) &
        +(-ssh_coast)/(dx0+dxm)*(1.-grd%msk(i,j-1)) ! force to drive bergs away from coasts
  hxp=2.*(grd%ssh(i,j+1)-grd%ssh(i,j))/(dx0+dxp)*grd%msk(i,j+1) &
        +(+ssh_coast)/(dx0+dxp)*(1.-grd%msk(i,j+1)) ! force to drive bergs away from coasts
#else
  if (xi>=0.5) then
    hxp=(xi-0.5)*ddy_ssh(grd,i+1,j  )+(1.5-xi)*ddy_ssh(grd,i  ,j  )
    hxm=(xi-0.5)*ddy_ssh(grd,i+1,j-1)+(1.5-xi)*ddy_ssh(grd,i  ,j-1)
  else
    hxp=(xi+0.5)*ddy_ssh(grd,i  ,j  )+(0.5-xi)*ddy_ssh(grd,i-1,j  )
    hxm=(xi+0.5)*ddy_ssh(grd,i  ,j-1)+(0.5-xi)*ddy_ssh(grd,i-1,j-1)
  endif
#endif
  ! ssh_y is at the v-point on a C-grid
  ssh_y=yj*hxp+(1.-yj)*hxm

  ! Rotate vectors from local grid to lat/lon coordinates
  call rotate(uo, vo, cos_rot, sin_rot)
  call rotate(ui, vi, cos_rot, sin_rot)
  call rotate(ua, va, cos_rot, sin_rot)
  call rotate(ssh_x, ssh_y, cos_rot, sin_rot)

  contains

  real function ddx_ssh(grd,i,j)
  ! Arguments
  type(icebergs_gridded), pointer :: grd
  integer, intent(in) :: i, j
  ! Local variables
  real :: dxp,dx0
    dxp=0.5*(grd%dx(i+1,j)+grd%dx(i+1,j-1))
    dx0=0.5*(grd%dx(i,j)+grd%dx(i,j-1))
    ddx_ssh=2.*(grd%ssh(i+1,j)-grd%ssh(i,j))/(dx0+dxp)*grd%msk(i+1,j)*grd%msk(i,j)
  end function ddx_ssh

  real function ddy_ssh(grd,i,j)
  ! Arguments
  type(icebergs_gridded), pointer :: grd
  integer, intent(in) :: i, j
  ! Local variables
  real :: dyp,dy0
    dyp=0.5*(grd%dy(i,j+1)+grd%dy(i-1,j+1))
    dy0=0.5*(grd%dy(i,j)+grd%dy(i-1,j))
    ddy_ssh=2.*(grd%ssh(i,j+1)-grd%ssh(i,j))/(dy0+dyp)*grd%msk(i,j+1)*grd%msk(i,j)
  end function ddy_ssh

  subroutine rotate(u, v, cos_rot, sin_rot)
  ! Arguments
  real, intent(inout) :: u, v
  real, intent(in) :: cos_rot, sin_rot
  ! Local variables
  real :: u_old, v_old

    u_old=u
    v_old=v
    u=cos_rot*u_old+sin_rot*v_old
    v=cos_rot*v_old-sin_rot*u_old

  end subroutine rotate

end subroutine interp_flds

! ##############################################################################

real function bilin(grd, fld, i, j, xi, yj)
! Arguments
type(icebergs_gridded), pointer :: grd
real, intent(in) :: fld(grd%isd:grd%ied,grd%jsd:grd%jed), xi, yj
integer, intent(in) :: i, j
! Local variables

  bilin=(fld(i,j  )*(1.-xi)+fld(i-1,j  )*xi)*(1.-yj) &
       +(fld(i,j-1)*(1.-xi)+fld(i-1,j-1)*xi)*yj

end function bilin

! ##############################################################################

subroutine icebergs_run(bergs, time, calving, uo, vo, ui, vi, tauxa, tauya, ssh, sst, calving_hflx, cn, hi)
! Arguments
type(icebergs), pointer :: bergs
type(time_type), intent(in) :: time
real, dimension(:,:), intent(inout) :: calving, calving_hflx
real, dimension(:,:), intent(in) :: uo, vo, ui, vi, tauxa, tauya, ssh, sst, cn, hi
! Local variables
integer :: iyr, imon, iday, ihr, imin, isec, k
type(icebergs_gridded), pointer :: grd
logical :: lerr, sample_traj, lbudget, lverbose
real :: unused_calving, tmpsum, grdd_berg_mass, grdd_bergy_mass
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()

  call mpp_clock_begin(bergs%clock)
  call mpp_clock_begin(bergs%clock_int)

  ! For convenience
  grd=>bergs%grd

  grd%floating_melt(:,:)=0.
  grd%berg_melt(:,:)=0.
  grd%melt_buoy(:,:)=0.
  grd%melt_eros(:,:)=0.
  grd%melt_conv(:,:)=0.
  grd%bergy_src(:,:)=0.
  grd%bergy_melt(:,:)=0.
  grd%bergy_mass(:,:)=0.
  grd%mass(:,:)=0.
  if (bergs%add_weight_to_ocean) grd%mass_on_ocean(:,:,:)=0.
  grd%virtual_area(:,:)=0.

  ! Manage time
  call get_date(time, iyr, imon, iday, ihr, imin, isec)
  bergs%current_year=iyr
  bergs%current_yearday=yearday(imon, iday, ihr, imin, isec)
  ! Turn on sampling of trajectories, verbosity, budgets
  sample_traj=.false.
  if (bergs%traj_sample_hrs>0) then
     if (mod(24*iday+ihr,bergs%traj_sample_hrs).eq.0) sample_traj=.true.
  end if
  lverbose=.false.
  if (bergs%verbose_hrs>0) then
     if (mod(24*iday+ihr,bergs%verbose_hrs).eq.0) lverbose=verbose
  end if
  lbudget=.false.
  if (bergs%verbose_hrs>0) then
     if (mod(24*iday+ihr,bergs%verbose_hrs).eq.0) lbudget=budget
  end if
  if (mpp_pe()==mpp_root_pe().and.lverbose) write(*,'(a,3i5,a,3i5,a,i5,f8.3)') &
       'diamonds: y,m,d=',iyr, imon, iday,' h,m,s=', ihr, imin, isec, &
       ' yr,yrdy=', bergs%current_year, bergs%current_yearday

  ! Adapt calving flux from coupler for use here
 !call sanitize_field(grd%calving,1.e20)
  tmpsum=sum( calving(:,:)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%net_calving_received=bergs%net_calving_received+tmpsum*bergs%dt
  grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec)=calving(:,:) ! Units of kg/m2/s
  grd%calving(:,:)=grd%calving(:,:)*grd%msk(:,:)*grd%area(:,:) ! Convert to kg/s
  tmpsum=sum( grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%net_incoming_calving=bergs%net_incoming_calving+tmpsum*bergs%dt
  if (grd%id_calving>0) &
    lerr=send_data(grd%id_calving, grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec), Time)

  ! Adapt calving heat flux from coupler
  grd%calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec)=calving_hflx(:,:) ! Units of W/m2
  grd%calving_hflx(:,:)=grd%calving_hflx(:,:)*grd%msk(:,:) ! Mask (just in case)
  if (grd%id_calving_hflx_in>0) &
    lerr=send_data(grd%id_calving_hflx_in, grd%calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  tmpsum=sum( grd%calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%net_incoming_calving_heat=bergs%net_incoming_calving_heat+tmpsum*bergs%dt ! Units of J

  ! Copy ocean flow (resides on B grid)
  grd%uo(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=uo(:,:)
  grd%vo(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=vo(:,:)
  call mpp_update_domains(grd%uo, grd%vo, grd%domain, gridtype=BGRID_NE)
  ! Copy ice flow (resides on B grid)
  grd%ui(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=ui(:,:)
  grd%vi(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=vi(:,:)
  call mpp_update_domains(grd%ui, grd%vi, grd%domain, gridtype=BGRID_NE)
  ! Copy atmospheric stress (resides on A grid)
  grd%ua(grd%isc:grd%iec,grd%jsc:grd%jec)=tauxa(:,:) ! Note rough conversion from stress to speed
  grd%va(grd%isc:grd%iec,grd%jsc:grd%jec)=tauya(:,:) ! Note rough conversion from stress to speed
  call invert_tau_for_du(grd%ua, grd%va) ! Note rough conversion from stress to speed
 !grd%ua(grd%isc:grd%iec,grd%jsc:grd%jec)=sign(sqrt(abs(tauxa(:,:))/0.01),tauxa(:,:))  ! Note rough conversion from stress to speed
 !grd%va(grd%isc:grd%iec,grd%jsc:grd%jec)=sign(sqrt(abs(tauya(:,:))/0.01),tauya(:,:))  ! Note rough conversion from stress to speed
  call mpp_update_domains(grd%ua, grd%va, grd%domain, gridtype=BGRID_NE)
  ! Copy sea surface height and temperature(resides on A grid)
  grd%ssh(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=ssh(:,:)
  call mpp_update_domains(grd%ssh, grd%domain)
  grd%sst(grd%isc:grd%iec,grd%jsc:grd%jec)=sst(:,:)-273.15 ! Note convert from Kelvin to Celsius
  call mpp_update_domains(grd%sst, grd%domain)
  ! Copy sea-ice concentration and thickness (resides on A grid)
  grd%cn(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=cn(:,:)
  call mpp_update_domains(grd%cn, grd%domain)
  grd%hi(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1)=hi(:,:)
  call mpp_update_domains(grd%hi, grd%domain)

  if (debug) call bergs_chksum(bergs, 'run bergs (top)')
  if (debug) call checksum_gridded(bergs%grd, 'top of s/r run')

  ! Accumulate ice from calving
  call accumulate_calving(bergs)
  if (grd%id_accum>0) then
    grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec)=calving(:,:)
    grd%tmp(:,:)=grd%tmp(:,:)*grd%msk(:,:)*grd%area(:,:)-grd%calving(:,:)
    lerr=send_data(grd%id_accum, grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  endif
  if (grd%id_unused>0) &
    lerr=send_data(grd%id_unused, grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  unused_calving=sum( grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec) )
  call mpp_clock_end(bergs%clock_int)

  call mpp_clock_begin(bergs%clock_cal)
  ! Calve excess stored ice into icebergs
  call calve_icebergs(bergs)
  if (debug) call bergs_chksum(bergs, 'run bergs (calved)')
  if (debug) call checksum_gridded(bergs%grd, 's/r run after calving')
  call mpp_clock_end(bergs%clock_cal)

  ! For each berg, evolve
  call mpp_clock_begin(bergs%clock_mom)
  if (associated(bergs%first)) call evolve_icebergs(bergs)
  if (debug) call bergs_chksum(bergs, 'run bergs (evolved)',ignore_halo_violation=.true.)
  if (debug) call checksum_gridded(bergs%grd, 's/r run after evolve')
  call mpp_clock_end(bergs%clock_mom)

  ! Send bergs to other PEs
  call mpp_clock_begin(bergs%clock_com)
  call send_bergs_to_other_pes(bergs)
  if (debug) call bergs_chksum(bergs, 'run bergs (exchanged)')
  if (debug) call checksum_gridded(bergs%grd, 's/r run after exchange')
  call mpp_clock_end(bergs%clock_com)

  ! Ice berg thermodynamics (melting) + rolling
  call mpp_clock_begin(bergs%clock_the)
  if (associated(bergs%first)) call thermodynamics(bergs)
  if (debug) call bergs_chksum(bergs, 'run bergs (thermo)')
  if (debug) call checksum_gridded(bergs%grd, 's/r run after thermodynamics')
  call mpp_clock_end(bergs%clock_the)

  ! For each berg, record
  call mpp_clock_begin(bergs%clock_dia)
  if (sample_traj.and.associated(bergs%first)) call record_posn(bergs)

  ! Gridded diagnostics
  if (grd%id_uo>0) &
    lerr=send_data(grd%id_uo, grd%uo(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_vo>0) &
    lerr=send_data(grd%id_vo, grd%vo(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_ui>0) &
    lerr=send_data(grd%id_ui, grd%ui(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_vi>0) &
    lerr=send_data(grd%id_vi, grd%vi(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_ua>0) &
    lerr=send_data(grd%id_ua, grd%ua(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_va>0) &
    lerr=send_data(grd%id_va, grd%va(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_sst>0) &
    lerr=send_data(grd%id_sst, grd%sst(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_cn>0) &
    lerr=send_data(grd%id_cn, grd%cn(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_hi>0) &
    lerr=send_data(grd%id_hi, grd%hi(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_floating_melt>0) &
    lerr=send_data(grd%id_floating_melt, grd%floating_melt(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_berg_melt>0) &
    lerr=send_data(grd%id_berg_melt, grd%berg_melt(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_melt_buoy>0) &
    lerr=send_data(grd%id_melt_buoy, grd%melt_buoy(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_melt_eros>0) &
    lerr=send_data(grd%id_melt_eros, grd%melt_eros(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_melt_conv>0) &
    lerr=send_data(grd%id_melt_conv, grd%melt_conv(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_virtual_area>0) &
    lerr=send_data(grd%id_virtual_area, grd%virtual_area(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_bergy_src>0) &
    lerr=send_data(grd%id_bergy_src, grd%bergy_src(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_bergy_melt>0) &
    lerr=send_data(grd%id_bergy_melt, grd%bergy_melt(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_bergy_mass>0) &
    lerr=send_data(grd%id_bergy_mass, grd%bergy_mass(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_mass>0) &
    lerr=send_data(grd%id_mass, grd%mass(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_stored_ice>0) &
    lerr=send_data(grd%id_stored_ice, grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:), Time)
  if (grd%id_real_calving>0) &
    lerr=send_data(grd%id_real_calving, grd%real_calving(grd%isc:grd%iec,grd%jsc:grd%jec,:), Time)

  ! Dump icebergs to screen
  if (really_debug) call print_bergs(stderrunit,bergs,'icebergs_run, status')
  call mpp_clock_end(bergs%clock_dia)

  ! Return what ever calving we did not use and additional icebergs melt
  call mpp_clock_begin(bergs%clock_int)
  if (.not. bergs%passive_mode) then
    where (grd%area(grd%isc:grd%iec,grd%jsc:grd%jec)>0.)
      calving(:,:)=grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec)/grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) &
                  +grd%floating_melt(grd%isc:grd%iec,grd%jsc:grd%jec)
    elsewhere
      calving(:,:)=0.
    end where
    calving_hflx(:,:)=grd%calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec)
  endif
  call mpp_clock_end(bergs%clock_int)

  ! Diagnose budgets
  call mpp_clock_begin(bergs%clock_dia)
  tmpsum=sum( grd%floating_melt(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%net_melt=bergs%net_melt+tmpsum*bergs%dt
  bergs%net_outgoing_calving=bergs%net_outgoing_calving+(unused_calving+tmpsum)*bergs%dt
  tmpsum=sum( grd%berg_melt(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%berg_melt=bergs%berg_melt+tmpsum*bergs%dt
  tmpsum=sum( grd%bergy_src(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%bergy_src=bergs%bergy_src+tmpsum*bergs%dt
  tmpsum=sum( grd%bergy_melt(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%bergy_melt=bergs%bergy_melt+tmpsum*bergs%dt
  tmpsum=sum( calving(:,:)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%net_calving_returned=bergs%net_calving_returned+tmpsum*bergs%dt
  tmpsum=sum( grd%calving_hflx(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
  bergs%net_outgoing_calving_heat=bergs%net_outgoing_calving_heat+tmpsum*bergs%dt ! Units of J
  if (lbudget) then
    bergs%stored_end=sum( grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:) )
    bergs%stored_heat_end=sum( grd%stored_heat(grd%isc:grd%iec,grd%jsc:grd%jec) )
    bergs%floating_mass_end=sum_mass(bergs%first)
    bergs%icebergs_mass_end=sum_mass(bergs%first,justbergs=.true.)
    bergs%bergy_mass_end=sum_mass(bergs%first,justbits=.true.)
    bergs%floating_heat_end=sum_heat(bergs%first)
    grd%tmpc(:,:)=0.;
    call mpp_clock_end(bergs%clock); call mpp_clock_end(bergs%clock_dia) ! To enable calling of public s/r
    call icebergs_incr_mass(bergs, grd%tmpc)
    call mpp_clock_begin(bergs%clock_dia); call mpp_clock_begin(bergs%clock) ! To enable calling of public s/r
    bergs%returned_mass_on_ocean=sum( grd%tmpc(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    bergs%nbergs_end=count_bergs(bergs)
    call mpp_sum(bergs%stored_end)
    call mpp_sum(bergs%stored_heat_end)
    call mpp_sum(bergs%floating_mass_end)
    call mpp_sum(bergs%icebergs_mass_end)
    call mpp_sum(bergs%bergy_mass_end)
    call mpp_sum(bergs%floating_heat_end)
    call mpp_sum(bergs%returned_mass_on_ocean)
    call mpp_sum(bergs%nbergs_end)
    call mpp_sum(bergs%nbergs_calved)
    do k=1,nclasses; call mpp_sum(bergs%nbergs_calved_by_class(k)); enddo
    call mpp_sum(bergs%nbergs_melted)
    call mpp_sum(bergs%nspeeding_tickets)
    call mpp_sum(bergs%net_calving_returned)
    call mpp_sum(bergs%net_outgoing_calving)
    call mpp_sum(bergs%net_calving_received)
    call mpp_sum(bergs%net_incoming_calving)
    call mpp_sum(bergs%net_incoming_calving_heat)
    call mpp_sum(bergs%net_incoming_calving_heat_used)
    call mpp_sum(bergs%net_outgoing_calving_heat)
    call mpp_sum(bergs%net_calving_used)
    call mpp_sum(bergs%net_calving_to_bergs)
    call mpp_sum(bergs%net_heat_to_bergs)
    call mpp_sum(bergs%net_heat_to_ocean)
    call mpp_sum(bergs%net_melt)
    call mpp_sum(bergs%berg_melt)
    call mpp_sum(bergs%bergy_src)
    call mpp_sum(bergs%bergy_melt)
    grdd_berg_mass=sum( grd%mass(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    call mpp_sum(grdd_berg_mass)
    grdd_bergy_mass=sum( grd%bergy_mass(grd%isc:grd%iec,grd%jsc:grd%jec)*grd%area(grd%isc:grd%iec,grd%jsc:grd%jec) )
    call mpp_sum(grdd_bergy_mass)
    if (mpp_pe().eq.mpp_root_pe()) then
 100 format("diamonds: ",a19,3(a18,"=",es14.7,x,a2,:,","),a12,i8)
 200 format("diamonds: ",a19,10(a18,"=",es14.7,x,a2,:,","))
      call report_state('stored ice','kg','',bergs%stored_start,'',bergs%stored_end,'')
      call report_state('floating','kg','',bergs%floating_mass_start,'',bergs%floating_mass_end,'',bergs%nbergs_end)
      call report_state('icebergs','kg','',bergs%icebergs_mass_start,'',bergs%icebergs_mass_end,'')
      call report_state('bits','kg','',bergs%bergy_mass_start,'',bergs%bergy_mass_end,'')
      call report_istate('berg #','',bergs%nbergs_start,'',bergs%nbergs_end,'')
      call report_ibudget('berg #','calved',bergs%nbergs_calved, &
                                   'melted',bergs%nbergs_melted, &
                                   '#',bergs%nbergs_start,bergs%nbergs_end)
      call report_budget('stored mass','kg','calving used',bergs%net_calving_used, &
                                            'bergs',bergs%net_calving_to_bergs, &
                                            'stored mass',bergs%stored_start,bergs%stored_end)
      call report_budget('floating mass','kg','calving used',bergs%net_calving_to_bergs, &
                                              'bergs',bergs%net_melt, &
                                              'stored mass',bergs%floating_mass_start,bergs%floating_mass_end)
      call report_budget('berg mass','kg','calving',bergs%net_calving_to_bergs, &
                                          'melt+eros',bergs%berg_melt, &
                                          'berg mass',bergs%icebergs_mass_start,bergs%icebergs_mass_end)
      call report_budget('bits mass','kg','eros used',bergs%bergy_src, &
                                          'bergs',bergs%bergy_melt, &
                                          'stored mass',bergs%bergy_mass_start,bergs%bergy_mass_end)
      call report_budget('net mass','kg','recvd',bergs%net_calving_received, &
                                         'rtrnd',bergs%net_calving_returned, &
                                         'net mass',bergs%stored_start+bergs%floating_mass_start, &
                                                    bergs%stored_end+bergs%floating_mass_end)
      call report_consistant('iceberg mass','kg','gridded',grdd_berg_mass,'bergs',bergs%icebergs_mass_end)
      call report_consistant('bits mass','kg','gridded',grdd_bergy_mass,'bits',bergs%bergy_mass_end)
      call report_consistant('wieght','kg','returned',bergs%returned_mass_on_ocean,'floating',bergs%floating_mass_end)
      call report_state('net heat','J','',bergs%stored_heat_start+bergs%floating_heat_start,'',&
           & bergs%stored_heat_end+bergs%floating_heat_end,'')
      call report_state('stored heat','J','',bergs%stored_heat_start,'',bergs%stored_heat_end,'')
      call report_state('floating heat','J','',bergs%floating_heat_start,'',bergs%floating_heat_end,'')
      call report_budget('net heat','J','net heat',bergs%net_incoming_calving_heat, &
                                        'net heat',bergs%net_outgoing_calving_heat, &
                                        'net heat',bergs%stored_heat_start+bergs%floating_heat_start, &
                                                   bergs%stored_heat_end+bergs%floating_heat_end)
      call report_budget('stored heat','J','calving used',bergs%net_incoming_calving_heat_used, &
                                           'bergs',bergs%net_heat_to_bergs, &
                                           'net heat',bergs%stored_heat_start,bergs%stored_heat_end)
      call report_budget('flting heat','J','calved',bergs%net_heat_to_bergs, &
                                           'melt',bergs%net_heat_to_ocean, &
                                           'net heat',bergs%floating_heat_start,bergs%floating_heat_end)
      if (debug) then
        call report_consistant('top interface','kg','from SIS',bergs%net_incoming_calving,'seen by diamonds',&
             & bergs%net_calving_received)
        call report_consistant('bot interface','kg','sent',bergs%net_outgoing_calving,'seen by SIS',bergs%net_calving_returned)
      endif
      write(*,'("diamonds: calved by class = ",i4,20(",",i4))') (bergs%nbergs_calved_by_class(k),k=1,nclasses)
      if (bergs%nspeeding_tickets>0) write(*,'("diamonds: speeding tickets issued = ",i4)') bergs%nspeeding_tickets
    endif
    bergs%nbergs_start=bergs%nbergs_end
    bergs%stored_start=bergs%stored_end
    bergs%nbergs_melted=0
    bergs%nbergs_calved=0
    bergs%nbergs_calved_by_class(:)=0
    bergs%nspeeding_tickets=0
    bergs%stored_heat_start=bergs%stored_heat_end
    bergs%floating_heat_start=bergs%floating_heat_end
    bergs%floating_mass_start=bergs%floating_mass_end
    bergs%icebergs_mass_start=bergs%icebergs_mass_end
    bergs%bergy_mass_start=bergs%bergy_mass_end
    bergs%net_calving_used=0.
    bergs%net_calving_to_bergs=0.
    bergs%net_heat_to_bergs=0.
    bergs%net_heat_to_ocean=0.
    bergs%net_calving_received=0.
    bergs%net_calving_returned=0.
    bergs%net_incoming_calving=0.
    bergs%net_outgoing_calving=0.
    bergs%net_incoming_calving_heat=0.
    bergs%net_incoming_calving_heat_used=0.
    bergs%net_outgoing_calving_heat=0.
    bergs%net_melt=0.
    bergs%berg_melt=0.
    bergs%bergy_melt=0.
    bergs%bergy_src=0.
  endif

  if (debug) call bergs_chksum(bergs, 'run bergs (bot)')
  if (debug) call checksum_gridded(bergs%grd, 'end of s/r run')
  call mpp_clock_end(bergs%clock_dia)

  call mpp_clock_end(bergs%clock)

  contains

  subroutine report_state(budgetstr,budgetunits,startstr,startval,endstr,endval,delstr,nbergs)
  ! Arguments
  character*(*), intent(in) :: budgetstr, budgetunits, startstr, endstr, delstr
  real, intent(in) :: startval, endval
  integer, intent(in), optional :: nbergs
  ! Local variables
  if (present(nbergs)) then
    write(*,100) budgetstr//' state:', &
                        startstr//' start',startval,budgetunits, &
                        endstr//' end',endval,budgetunits, &
                        'Delta '//delstr,endval-startval,budgetunits, &
                        '# of bergs',nbergs
  else
    write(*,100) budgetstr//' state:', &
                        startstr//' start',startval,budgetunits, &
                        endstr//' end',endval,budgetunits, &
                        delstr//'Delta',endval-startval,budgetunits
  endif
  100 format("diamonds: ",a19,3(a18,"=",es14.7,x,a2,:,","),a12,i8)
  end subroutine report_state

  subroutine report_consistant(budgetstr,budgetunits,startstr,startval,endstr,endval)
  ! Arguments
  character*(*), intent(in) :: budgetstr, budgetunits, startstr, endstr
  real, intent(in) :: startval, endval
  ! Local variables
  write(*,200) budgetstr//' check:', &
                      startstr,startval,budgetunits, &
                      endstr,endval,budgetunits, &
                      'error',(endval-startval)/((endval+startval)+1e-30),'nd'
  200 format("diamonds: ",a19,10(a18,"=",es14.7,x,a2,:,","))
  end subroutine report_consistant

  subroutine report_budget(budgetstr,budgetunits,instr,inval,outstr,outval,delstr,startval,endval)
  ! Arguments
  character*(*), intent(in) :: budgetstr, budgetunits, instr, outstr, delstr
  real, intent(in) :: inval, outval, startval, endval
  ! Local variables
  write(*,200) budgetstr//' budget:', &
                      instr//' in',inval,budgetunits, &
                      outstr//' out',outval,budgetunits, &
                      'Delta '//delstr,inval-outval,budgetunits, &
                      'error',((endval-startval)-(inval-outval))/max(1.e-30,max(abs(endval-startval),abs(inval-outval))),'nd'
  200 format("diamonds: ",a19,3(a18,"=",es14.7,x,a2,:,","),a8,"=",es10.3,x,a2)
  end subroutine report_budget

  subroutine report_istate(budgetstr,startstr,startval,endstr,endval,delstr)
  ! Arguments
  character*(*), intent(in) :: budgetstr, startstr, endstr, delstr
  integer, intent(in) :: startval, endval
  ! Local variables
  write(*,100) budgetstr//' state:', &
                        startstr//' start',startval, &
                        endstr//' end',endval, &
                        delstr//'Delta',endval-startval
  100 format("diamonds: ",a19,3(a18,"=",i14,x,:,","))
  end subroutine report_istate

  subroutine report_ibudget(budgetstr,instr,inval,outstr,outval,delstr,startval,endval)
  ! Arguments
  character*(*), intent(in) :: budgetstr, instr, outstr, delstr
  integer, intent(in) :: inval, outval, startval, endval
  ! Local variables
  write(*,200) budgetstr//' budget:', &
                      instr//' in',inval, &
                      outstr//' out',outval, &
                      'Delta '//delstr,inval-outval, &
                      'error',((endval-startval)-(inval-outval))
  200 format("diamonds: ",a19,10(a18,"=",i14,x,:,","))
  end subroutine report_ibudget

end subroutine icebergs_run

! ##############################################################################

subroutine icebergs_incr_mass(bergs, mass)
! Arguments
type(icebergs), pointer :: bergs
real, dimension(bergs%grd%isc:bergs%grd%iec,bergs%grd%jsc:bergs%grd%jec), intent(inout) :: mass
! Local variables
integer :: i, j
type(icebergs_gridded), pointer :: grd
real :: dmda

  if (.not. associated(bergs)) return

  if (.not. bergs%add_weight_to_ocean) return

  call mpp_clock_begin(bergs%clock)
  call mpp_clock_begin(bergs%clock_int)

  ! For convenience
  grd=>bergs%grd

  ! Add iceberg+bits mass field to non-haloed SIS field (kg/m^2)
 !mass(:,:)=mass(:,:)+( grd%mass(grd%isc:grd%iec,grd%jsc:grd%jec) &
 !                    + grd%bergy_mass(grd%isc:grd%iec,grd%jsc:grd%jec) )


  if (debug) then
    bergs%grd%tmp(:,:)=0.; bergs%grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec)=mass
    call grd_chksum2(bergs%grd, bergs%grd%tmp, 'mass in (incr)')
  endif

  call mpp_update_domains(grd%mass_on_ocean, grd%domain)
  if (.not. old_bug_rotated_weights) then
    do j=grd%jsd, grd%jed; do i=grd%isd, grd%ied
      if (grd%parity_x(i,j)<0.) then
        ! This block assumes both parity_x and parity_y are negative
        ! (i.e. a 180 degree rotation). In general, we should handle
        ! +/- 90 degree rotations as well but in CM2*-class models
        ! this is not necessary. -aja
        dmda=grd%mass_on_ocean(i,j,9); grd%mass_on_ocean(i,j,9)=grd%mass_on_ocean(i,j,1); grd%mass_on_ocean(i,j,1)=dmda
        dmda=grd%mass_on_ocean(i,j,8); grd%mass_on_ocean(i,j,8)=grd%mass_on_ocean(i,j,2); grd%mass_on_ocean(i,j,2)=dmda
        dmda=grd%mass_on_ocean(i,j,7); grd%mass_on_ocean(i,j,7)=grd%mass_on_ocean(i,j,3); grd%mass_on_ocean(i,j,3)=dmda
        dmda=grd%mass_on_ocean(i,j,6); grd%mass_on_ocean(i,j,6)=grd%mass_on_ocean(i,j,4); grd%mass_on_ocean(i,j,4)=dmda
      endif
    enddo; enddo
  endif
  do j=grd%jsc, grd%jec; do i=grd%isc, grd%iec
    dmda=grd%mass_on_ocean(i,j,5) &
         + ( ( (grd%mass_on_ocean(i-1,j-1,9)+grd%mass_on_ocean(i+1,j+1,1))   &
         +     (grd%mass_on_ocean(i+1,j-1,7)+grd%mass_on_ocean(i-1,j+1,3)) ) &
         +   ( (grd%mass_on_ocean(i-1,j  ,6)+grd%mass_on_ocean(i+1,j  ,4))   &
         +     (grd%mass_on_ocean(i  ,j-1,8)+grd%mass_on_ocean(i  ,j+1,2)) ) )
    if (grd%area(i,j)>0) dmda=dmda/grd%area(i,j)*grd%msk(i,j)
    if (.not. bergs%passive_mode) mass(i,j)=mass(i,j)+dmda
  enddo; enddo

  if (debug) then
    call grd_chksum3(bergs%grd, bergs%grd%mass_on_ocean, 'mass bergs (incr)')
    bergs%grd%tmp(:,:)=0.; bergs%grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec)=mass
    call grd_chksum2(bergs%grd, bergs%grd%tmp, 'mass out (incr)')
  endif

  call mpp_clock_end(bergs%clock_int)
  call mpp_clock_end(bergs%clock)

end subroutine icebergs_incr_mass

! ##############################################################################

subroutine accumulate_calving(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(icebergs_gridded), pointer :: grd
real :: remaining_dist, net_calving_used
integer :: k, i, j
logical, save :: first_call=.true.
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>bergs%grd

  ! This is a hack to simplify initialization
  if (first_call.and..not.bergs%restarted) then
    first_call=.false.
   !do k=1, nclasses
   !  where (grd%calving==0.) grd%stored_ice(:,:,k)=0.
   !enddo
    bergs%stored_start=sum( grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:) )
    call mpp_sum( bergs%stored_start )
    if (mpp_pe().eq.mpp_root_pe()) write(*,'(a,es13.6,a)') &
        'diamonds, accumulate_calving: initial stored mass=',bergs%stored_start,' kg'
    do j=grd%jsc,grd%jec; do i=grd%isc,grd%iec
      if (grd%calving(i,j).ne.0.) grd%stored_heat(i,j)= & ! Need units of J
            sum(grd%stored_ice(i,j,:)) & ! initial stored ice in kg
           *grd%calving_hflx(i,j)*grd%area(i,j) & ! J/s/m2 x m^2 = J/s
           /grd%calving(i,j) ! /calving in kg/s
    enddo; enddo
    bergs%stored_heat_start=sum( grd%stored_heat(grd%isc:grd%iec,grd%jsc:grd%jec) )
    call mpp_sum( bergs%stored_heat_start )
    if (mpp_pe().eq.mpp_root_pe()) write(*,'(a,es13.6,a)') &
        'diamonds, accumulate_calving: initial stored heat=',bergs%stored_heat_start,' J'
   endif

  remaining_dist=1.
  do k=1, nclasses
    grd%stored_ice(:,:,k)=grd%stored_ice(:,:,k)+bergs%dt*grd%calving(:,:)*bergs%distribution(k)
    remaining_dist=remaining_dist-bergs%distribution(k)
  enddo
  if (remaining_dist.lt.0.) then
    write(stderrunit,*) 'diamonds, accumulate_calving: sum(distribution)>1!!!',remaining_dist
    call error_mesg('diamonds, accumulate_calving', 'calving is OVER distributed!', WARNING)
  endif
  net_calving_used=sum( grd%calving(grd%isc:grd%iec,grd%jsc:grd%jec) )*(1.-remaining_dist)
  bergs%net_calving_used=bergs%net_calving_used+net_calving_used*bergs%dt
  ! Remove the calving accounted for by accumulation
  grd%calving(:,:)=grd%calving(:,:)*remaining_dist

  ! Do the same for heat (no separate classes needed)
  grd%tmp(:,:)=bergs%dt*grd%calving_hflx(:,:)*grd%area(:,:)*(1.-remaining_dist)
  bergs%net_incoming_calving_heat_used=bergs%net_incoming_calving_heat_used+sum( grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec) )
  grd%stored_heat(:,:)=grd%stored_heat(:,:)+grd%tmp(:,:) ! +=bergs%dt*grd%calving_hflx(:,:)*grd%area(:,:)*(1.-remaining_dist)
  grd%calving_hflx(:,:)=grd%calving_hflx(:,:)*remaining_dist

end subroutine accumulate_calving

! ##############################################################################

subroutine calve_icebergs(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(icebergs_gridded), pointer :: grd
integer :: i,j,k,icnt,icntmax
type(iceberg) :: newberg
logical :: lret
real :: xi, yj, ddt, calving_to_bergs, calved_to_berg, heat_to_bergs, heat_to_berg
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>bergs%grd

  grd%real_calving(:,:,:)=0.
  calving_to_bergs=0.
  heat_to_bergs=0.
  icntmax=0

  do k=1, nclasses
    do j=grd%jsc, grd%jec
      do i=grd%isc, grd%iec
        ddt=0.; icnt=0
        do while (grd%stored_ice(i,j,k).ge.bergs%initial_mass(k)*bergs%mass_scaling(k))
          newberg%lon=0.25*((grd%lon(i,j)+grd%lon(i-1,j-1))+(grd%lon(i-1,j)+grd%lon(i,j-1)))
          newberg%lat=0.25*((grd%lat(i,j)+grd%lat(i-1,j-1))+(grd%lat(i-1,j)+grd%lat(i,j-1)))
         !write(stderr(),*) 'diamonds, calve_icebergs: creating new iceberg at ',newberg%lon,newberg%lat
          lret=pos_within_cell(grd, newberg%lon, newberg%lat, i, j, xi, yj)
          if (.not.lret) then
            write(stderrunit,*) 'diamonds, calve_icebergs: something went very wrong!',i,j,xi,yj
            call error_mesg('diamonds, calve_icebergs', 'berg is not in the correct cell!', FATAL)
          endif
          if (debug.and.(xi<0..or.xi>1..or.yj<0..or.yj>1.)) then
            write(stderrunit,*) 'diamonds, calve_icebergs: something went very wrong!',i,j,xi,yj
            call error_mesg('diamonds, calve_icebergs', 'berg xi,yj is not correct!', FATAL)
          endif
          newberg%ine=i
          newberg%jne=j
          newberg%xi=xi
          newberg%yj=yj
          newberg%uvel=0.
          newberg%vvel=0.
          newberg%mass=bergs%initial_mass(k)
          newberg%thickness=bergs%initial_thickness(k)
          newberg%width=bergs%initial_width(k)
          newberg%length=bergs%initial_length(k)
          newberg%start_lon=newberg%lon
          newberg%start_lat=newberg%lat
          newberg%start_year=bergs%current_year
          newberg%start_day=bergs%current_yearday+ddt/86400.
          newberg%start_mass=bergs%initial_mass(k)
          newberg%mass_scaling=bergs%mass_scaling(k)
          newberg%mass_of_bits=0.
          newberg%heat_density=grd%stored_heat(i,j)/grd%stored_ice(i,j,k) ! This is in J/kg
          call add_new_berg_to_list(bergs%first, newberg)
          calved_to_berg=bergs%initial_mass(k)*bergs%mass_scaling(k) ! Units of kg
          ! Heat content
          heat_to_berg=calved_to_berg*newberg%heat_density ! Units of J
          grd%stored_heat(i,j)=grd%stored_heat(i,j)-heat_to_berg
          heat_to_bergs=heat_to_bergs+heat_to_berg
          ! Stored mass
          grd%stored_ice(i,j,k)=grd%stored_ice(i,j,k)-calved_to_berg
          calving_to_bergs=calving_to_bergs+calved_to_berg
          grd%real_calving(i,j,k)=grd%real_calving(i,j,k)+calved_to_berg/bergs%dt
          ddt=ddt+bergs%dt*2./17. ! Minor offset to start day
          icnt=icnt+1
          bergs%nbergs_calved=bergs%nbergs_calved+1
          bergs%nbergs_calved_by_class(k)=bergs%nbergs_calved_by_class(k)+1
        enddo
        icntmax=max(icntmax,icnt)
      enddo
    enddo
  enddo

  if (debug.and.icntmax>1) write(stderrunit,*) 'calve_icebergs: icnt=',icnt,' on',mpp_pe()

  bergs%net_calving_to_bergs=bergs%net_calving_to_bergs+calving_to_bergs
  bergs%net_heat_to_bergs=bergs%net_heat_to_bergs+heat_to_bergs

end subroutine calve_icebergs

! ##############################################################################

subroutine evolve_icebergs(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(icebergs_gridded), pointer :: grd
real :: uvel1, vvel1, lon1, lat1, u1, v1, dxdl1, ax1, ay1
real :: uvel2, vvel2, lon2, lat2, u2, v2, dxdl2, ax2, ay2
real :: uvel3, vvel3, lon3, lat3, u3, v3, dxdl3, ax3, ay3
real :: uvel4, vvel4, lon4, lat4, u4, v4, dxdl4, ax4, ay4
real :: uveln, vveln, lonn, latn
real :: x1, xdot1, xddot1, y1, ydot1, yddot1
real :: x2, xdot2, xddot2, y2, ydot2, yddot2
real :: x3, xdot3, xddot3, y3, ydot3, yddot3
real :: x4, xdot4, xddot4, y4, ydot4, yddot4
real :: xn, xdotn, yn, ydotn
real :: r180_pi, dt, dt_2, dt_6, dydl, Rearth
integer :: i, j
integer :: i1,j1,i2,j2,i3,j3,i4,j4
real :: xi, yj
logical :: bounced, on_tangential_plane, error_flag
type(iceberg), pointer :: berg
integer :: stderrunit

  ! 4th order Runge-Kutta to solve:
  !    d/dt X = V,  d/dt V = A
  ! with I.C.'s:
  !    X=X1 and V=V1
  !
  !  A1 = A(X1)
  !  X2 = X1+dt/2*V1 ; V2 = V1+dt/2*A1; A2=A(X2)
  !  X3 = X1+dt/2*V2 ; V3 = V1+dt/2*A2; A3=A(X3)
  !  X4 = X1+  dt*V3 ; V4 = V1+  dt*A3; A4=A(X4)
  !
  !  Xn = X1+dt*(V1+2*V2+2*V3+V4)/6
  !  Vn = V1+dt*(A1+2*A2+2*A3+A4)/6

  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>bergs%grd

  ! Common constants
  r180_pi=1./pi_180
  dt=bergs%dt
  dt_2=0.5*dt
  dt_6=dt/6.
  Rearth=6360.e3

  berg=>bergs%first
  do while (associated(berg)) ! loop over all bergs

  if (.not. is_point_in_cell(bergs%grd, berg%lon, berg%lat, berg%ine, berg%jne) ) then
    write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lon',(i,i=grd%isd,grd%ied)
    do j=grd%jed,grd%jsd,-1
      write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lon(i,j),i=grd%isd,grd%ied)
    enddo
    write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lat',(i,i=grd%isd,grd%ied)
    do j=grd%jed,grd%jsd,-1
      write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lat(i,j),i=grd%isd,grd%ied)
    enddo
    call print_berg(stderrunit, berg, 'evolve_iceberg, berg is not in proper starting cell')
    write(stderrunit,'(a,i3,2(i4,3f8.2))') 'evolve_iceberg: pe,lon/lat(i,j)=', mpp_pe(), &
             berg%ine,berg%lon,grd%lon(berg%ine-1,berg%jne-1),grd%lon(berg%ine,berg%jne), &
             berg%jne,berg%lat,grd%lat(berg%ine-1,berg%jne-1),grd%lat(berg%ine,berg%jne)
    if (debug) call error_mesg('diamonds, evolve_iceberg','berg is in wrong starting cell!',FATAL)
  endif

  if (debug) call check_position(grd, berg, 'evolve_iceberg (top)')

  i=berg%ine
  j=berg%jne
  xi=berg%xi
  yj=berg%yj
  bounced=.false.
  on_tangential_plane=.false.
  if (berg%lat>89.) on_tangential_plane=.true.
  i1=i;j1=j
  if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
    call spread_mass_across_ocean_cells(grd, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling)

  ! A1 = A(X1)
  lon1=berg%lon; lat1=berg%lat
  if (on_tangential_plane) call rotpos_to_tang(lon1,lat1,x1,y1)
  dxdl1=r180_pi/(Rearth*cos(lat1*pi_180))
  dydl=r180_pi/Rearth
  uvel1=berg%uvel; vvel1=berg%vvel
  if (on_tangential_plane) call rotvec_to_tang(lon1,uvel1,vvel1,xdot1,ydot1)
  u1=uvel1*dxdl1; v1=vvel1*dydl
  call accel(bergs, berg, i, j, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, ax1, ay1)
  if (on_tangential_plane) call rotvec_to_tang(lon1,ax1,ay1,xddot1,yddot1)

  !  X2 = X1+dt/2*V1 ; V2 = V1+dt/2*A1; A2=A(X2)
 !if (debug) write(stderr(),*) 'diamonds, evolve: x2=...'
  if (on_tangential_plane) then
    x2=x1+dt_2*xdot1; y2=y1+dt_2*ydot1
    xdot2=xdot1+dt_2*xddot1; ydot2=ydot1+dt_2*yddot1
    call rotpos_from_tang(x2,y2,lon2,lat2)
    call rotvec_from_tang(lon2,xdot2,ydot2,uvel2,vvel2)
  else
    lon2=lon1+dt_2*u1; lat2=lat1+dt_2*v1
    uvel2=uvel1+dt_2*ax1; vvel2=vvel1+dt_2*ay1
  endif
  i=i1;j=j1;xi=berg%xi;yj=berg%yj
  call adjust_index_and_ground(grd, lon2, lat2, uvel2, vvel2, i, j, xi, yj, bounced, error_flag)
  i2=i; j2=j
  if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
    call spread_mass_across_ocean_cells(grd, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling)
  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon2,lat2,x2,y2)
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(bergs%grd, lon2, lat2, i, j)) error_flag=.true.
  endif
  if (error_flag) then
   call print_fld(grd, grd%msk, 'msk')
   call print_fld(grd, grd%ssh, 'ssh')
   call print_fld(grd, grd%sst, 'sst')
   call print_fld(grd, grd%hi, 'hi')
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: i1,i2,i=',i1,i2,i
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: j1,j2,j=',j1,j2,j
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lon1,lon2=',lon1,lon2,berg%lon
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lat1,lat2=',lat1,lat2,berg%lat
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: u1,u2,u0=',uvel1,uvel2,berg%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: v1,v2,v0=',vvel1,vvel2,berg%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ax1,ax2=',dt*ax1,dt*ax2
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ay1,ay2=',dt*ay1,dt*ay2
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2,u0=',dt*uvel1,dt*uvel2,dt*berg%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2,v0=',dt*vvel1,dt*vvel2,dt*berg%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2 (deg)=',dt*u1,dt*u2
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2 (deg)=',dt*v1,dt*v2
   write(stderrunit,*) 'Acceleration terms for position 1'
   error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
   call accel(bergs, berg, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, ax1, ay1, debug_flag=.true.)
    call print_berg(stderrunit, berg, 'evolve_iceberg, out of position at 2')
    write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'pos2 i,j,lon,lat,xi,yj=',i,j,lon2,lat2,xi,yj
    write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'pos2 box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
    bounced=is_point_in_cell(bergs%grd, lon2, lat2, i, j, explain=.true.)
    call error_mesg('diamonds, evolve_iceberg','berg is out of posn at 2!',FATAL)
  endif
  dxdl2=r180_pi/(Rearth*cos(lat2*pi_180))
  u2=uvel2*dxdl2; v2=vvel2*dydl
  call accel(bergs, berg, i, j, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, ax2, ay2)
  if (on_tangential_plane) call rotvec_to_tang(lon2,ax2,ay2,xddot2,yddot2)

  !  X3 = X1+dt/2*V2 ; V3 = V1+dt/2*A2; A3=A(X3)
 !if (debug) write(stderr(),*) 'diamonds, evolve: x3=...'
  if (on_tangential_plane) then
    x3=x1+dt_2*xdot2; y3=y1+dt_2*ydot2
    xdot3=xdot1+dt_2*xddot2; ydot3=ydot1+dt_2*yddot2
    call rotpos_from_tang(x3,y3,lon3,lat3)
    call rotvec_from_tang(lon3,xdot3,ydot3,uvel3,vvel3)
  else
    lon3=lon1+dt_2*u2; lat3=lat1+dt_2*v2
    uvel3=uvel1+dt_2*ax2; vvel3=vvel1+dt_2*ay2
  endif
  i=i1;j=j1;xi=berg%xi;yj=berg%yj
  call adjust_index_and_ground(grd, lon3, lat3, uvel3, vvel3, i, j, xi, yj, bounced, error_flag)
  i3=i; j3=j
  if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
    call spread_mass_across_ocean_cells(grd, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling)
  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon3,lat3,x3,y3)
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(bergs%grd, lon3, lat3, i, j)) error_flag=.true.
  endif
  if (error_flag) then
   call print_fld(grd, grd%msk, 'msk')
   call print_fld(grd, grd%ssh, 'ssh')
   call print_fld(grd, grd%sst, 'sst')
   call print_fld(grd, grd%hi, 'hi')
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: i1,i2,i3,i=',i1,i2,i3,i
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: j1,j2,j3,j=',j1,j2,j3,j
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lon1,lon2,lon3=',lon1,lon2,lon3,berg%lon
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lat1,lat2,lat3=',lat1,lat2,lat3,berg%lat
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: u1,u2,u3,u0=',uvel1,uvel2,uvel3,berg%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: v1,v2,v3,v0=',vvel1,vvel2,vvel3,berg%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ax1,ax2,ax3=',dt*ax1,dt*ax2,dt*ax3
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ay1,ay2,ay3=',dt*ay1,dt*ay2,dt*ay3
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2,u3,u0=',dt*uvel1,dt*uvel2,dt*uvel3,dt*berg%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2,v3,v0=',dt*vvel1,dt*vvel2,dt*vvel3,dt*berg%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2,u3 (deg)=',dt*u1,dt*u2,dt*u3
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2,v3 (deg)=',dt*v1,dt*v2,dt*v3
   write(stderrunit,*) 'Acceleration terms for position 1'
   error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
   call accel(bergs, berg, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, ax1, ay1, debug_flag=.true.)
   write(stderrunit,*) 'Acceleration terms for position 2'
   error_flag=pos_within_cell(grd, lon2, lat2, i2, j2, xi, yj)
   call accel(bergs, berg, i2, j2, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, ax2, ay2, debug_flag=.true.)
    call print_berg(stderrunit, berg, 'evolve_iceberg, out of position at 3')
    write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'pos3 i,j,lon,lat,xi,yj=',i,j,lon3,lat3,xi,yj
    write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'pos3 box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
    bounced=is_point_in_cell(bergs%grd, lon2, lat2, i, j, explain=.true.)
    call error_mesg('diamonds, evolve_iceberg','berg is out of posn at 3!',FATAL)
  endif
  dxdl3=r180_pi/(Rearth*cos(lat3*pi_180))
  u3=uvel3*dxdl3; v3=vvel3*dydl
  call accel(bergs, berg, i, j, xi, yj, lat3, uvel3, vvel3, uvel1, vvel1, dt, ax3, ay3)
  if (on_tangential_plane) call rotvec_to_tang(lon3,ax3,ay3,xddot3,yddot3)

  !  X4 = X1+dt*V3 ; V4 = V1+dt*A3; A4=A(X4)
 !if (debug) write(stderr(),*) 'diamonds, evolve: x4=...'
  if (on_tangential_plane) then
    x4=x1+dt*xdot3; y4=y1+dt*ydot3
    xdot4=xdot1+dt*xddot3; ydot4=ydot1+dt*yddot3
    call rotpos_from_tang(x4,y4,lon4,lat4)
    call rotvec_from_tang(lon4,xdot4,ydot4,uvel4,vvel4)
  else
    lon4=lon1+dt*u3; lat4=lat1+dt*v3
    uvel4=uvel1+dt*ax3; vvel4=vvel1+dt*ay3
  endif
  i=i1;j=j1;xi=berg%xi;yj=berg%yj
  call adjust_index_and_ground(grd, lon4, lat4, uvel4, vvel4, i, j, xi, yj, bounced, error_flag)
  i4=i; j4=j
  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon4,lat4,x4,y4)
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(bergs%grd, lon4, lat4, i, j)) error_flag=.true.
  endif
  if (error_flag) then
   call print_fld(grd, grd%msk, 'msk')
   call print_fld(grd, grd%ssh, 'ssh')
   call print_fld(grd, grd%sst, 'sst')
   call print_fld(grd, grd%hi, 'hi')
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: i1,i2,i3,i4,i=',i1,i2,i3,i4,i
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: j1,j2,j3,j4,j=',j1,j2,j3,j4,j
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lon1,lon2,lon3,lon4=',lon1,lon2,lon3,lon4,berg%lon
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lat1,lat2,lat3,lat4=',lat1,lat2,lat3,lat4,berg%lat
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: u1,u2,u3,u4,u0=',uvel1,uvel2,uvel3,uvel4,berg%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: v1,v2,v3,v4,v0=',vvel1,vvel2,vvel3,vvel4,berg%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ax1,ax2,ax3,ax4=',dt*ax1,dt*ax2,dt*ax3,dt*ax4
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ay1,ay2,ay3,ay4=',dt*ay1,dt*ay2,dt*ay3,dt*ay4
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2,u3,u4,u0=',dt*uvel1,dt*uvel2,dt*uvel3,dt*uvel4,dt*berg%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2,v3,v4,v0=',dt*vvel1,dt*vvel2,dt*vvel3,dt*vvel4,dt*berg%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2,u3,u4 (deg)=',dt*u1,dt*u2,dt*u3,dt*u4
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2,v3,v4 (deg)=',dt*v1,dt*v2,dt*v3,dt*v4
   write(stderrunit,*) 'Acceleration terms for position 1'
   error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
   call accel(bergs, berg, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, ax1, ay1, debug_flag=.true.)
   write(stderrunit,*) 'Acceleration terms for position 2'
   error_flag=pos_within_cell(grd, lon2, lat2, i2, j2, xi, yj)
   call accel(bergs, berg, i2, j2, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, ax2, ay2, debug_flag=.true.)
   write(stderrunit,*) 'Acceleration terms for position 3'
   error_flag=pos_within_cell(grd, lon3, lat3, i3, j3, xi, yj)
   call accel(bergs, berg, i3, j3, xi, yj, lat3, uvel3, vvel3, uvel1, vvel1, dt, ax3, ay3, debug_flag=.true.)
    call print_berg(stderrunit, berg, 'evolve_iceberg, out of position at 4')
    write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'pos4 i,j,lon,lat,xi,yj=',i,j,lon4,lat4,xi,yj
    write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'pos4 box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
    bounced=is_point_in_cell(bergs%grd, lon2, lat2, i, j, explain=.true.)
    call error_mesg('diamonds, evolve_iceberg','berg is out of posn at 4!',FATAL)
  endif
  dxdl4=r180_pi/(Rearth*cos(lat4*pi_180))
  u4=uvel4*dxdl4; v4=vvel4*dydl
  call accel(bergs, berg, i, j, xi, yj, lat4, uvel4, vvel4, uvel1, vvel1, dt, ax4, ay4)
  if (on_tangential_plane) call rotvec_to_tang(lon4,ax4,ay4,xddot4,yddot4)

  !  Xn = X1+dt*(V1+2*V2+2*V3+V4)/6
  !  Vn = V1+dt*(A1+2*A2+2*A3+A4)/6
  if (on_tangential_plane) then
    xn=x1+dt_6*( (xdot1+xdot4)+2.*(xdot2+xdot3) )
    yn=y1+dt_6*( (ydot1+ydot4)+2.*(ydot2+ydot3) )
    xdotn=xdot1+dt_6*( (xddot1+xddot4)+2.*(xddot2+xddot3) )
    ydotn=ydot1+dt_6*( (yddot1+yddot4)+2.*(yddot2+yddot3) )
    call rotpos_from_tang(xn,yn,lonn,latn)
    call rotvec_from_tang(lonn,xdotn,ydotn,uveln,vveln)
  else
    lonn=berg%lon+dt_6*( (u1+u4)+2.*(u2+u3) )
    latn=berg%lat+dt_6*( (v1+v4)+2.*(v2+v3) )
    uveln=berg%uvel+dt_6*( (ax1+ax4)+2.*(ax2+ax3) )
    vveln=berg%vvel+dt_6*( (ay1+ay4)+2.*(ay2+ay3) )
  endif
  i=i1;j=j1;xi=berg%xi;yj=berg%yj
  call adjust_index_and_ground(grd, lonn, latn, uveln, vveln, i, j, xi, yj, bounced, error_flag)
  if (bergs%add_weight_to_ocean .and. bergs%time_average_weight) &
    call spread_mass_across_ocean_cells(grd, i, j, xi, yj, berg%mass, berg%mass_of_bits, 0.25*berg%mass_scaling)

  if (.not.error_flag) then
    if (.not. is_point_in_cell(bergs%grd, lonn, latn, i, j)) error_flag=.true.
  endif
  if (error_flag) then
   call print_fld(grd, grd%msk, 'msk')
   call print_fld(grd, grd%ssh, 'ssh')
   call print_fld(grd, grd%sst, 'sst')
   call print_fld(grd, grd%hi, 'hi')
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: i1,i2,i3,i4,i=',i1,i2,i3,i4,i
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_iceberg: j1,j2,j3,j4,j=',j1,j2,j3,j4,j
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lon1,lon2,lon3,lon4,lonn=',lon1,lon2,lon3,lon4,lonn,berg%lon
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: lat1,lat2,lat3,lat4,latn=',lat1,lat2,lat3,lat4,latn,berg%lat
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: u1,u2,u3,u4,un,u0=',uvel1,uvel2,uvel3,uvel4,uveln,berg%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: v1,v2,v3,v4,vn,v0=',vvel1,vvel2,vvel3,vvel4,vveln,berg%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ax1,ax2,ax3,ax4,axn=',&
        & dt*ax1,dt*ax2,dt*ax3,dt*ax4,dt_6*( (ax1+ax4)+2.*(ax2+ax3) )
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* ay1,ay2,ay3,ay4,ayn=',&
        & dt*ay1,dt*ay2,dt*ay3,dt*ay4,dt_6*( (ay1+ay4)+2.*(ay2+ay3) )
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2,u3,u4,un,u0=',&
        & dt*uvel1,dt*uvel2,dt*uvel3,dt*uvel4,dt*uveln,dt*berg%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2,v3,v4,vn,v0=',&
        & dt*vvel1,dt*vvel2,dt*vvel3,dt*vvel4,dt*vveln,dt*berg%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* u1,u2,u3,u4,u_rk (deg)=',&
        & dt*u1,dt*u2,dt*u3,dt*u4,dt_6*( (u1+u4)+2.*(u2+u3) )
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_iceberg: dt* v1,v2,v3,v4,v_rk (deg)=',&
        & dt*v1,dt*v2,dt*v3,dt*v4,dt_6*( (v1+v4)+2.*(v2+v3) )
   write(stderrunit,*) 'diamonds, evolve_iceberg: on_tangential_plane=',on_tangential_plane
   write(stderrunit,*) 'Acceleration terms for position 1'
   error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
   call accel(bergs, berg, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, ax1, ay1, debug_flag=.true.)
   write(stderrunit,*) 'Acceleration terms for position 2'
   error_flag=pos_within_cell(grd, lon2, lat2, i2, j2, xi, yj)
   call accel(bergs, berg, i2, j2, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, ax2, ay2, debug_flag=.true.)
   write(stderrunit,*) 'Acceleration terms for position 3'
   error_flag=pos_within_cell(grd, lon3, lat3, i3, j3, xi, yj)
   call accel(bergs, berg, i3, j3, xi, yj, lat3, uvel3, vvel3, uvel1, vvel1, dt, ax3, ay3, debug_flag=.true.)
   write(stderrunit,*) 'Acceleration terms for position 4'
   error_flag=pos_within_cell(grd, lon4, lat4, i4, j4, xi, yj)
   call accel(bergs, berg, i4, j4, xi, yj, lat4, uvel4, vvel4, uvel1, vvel1, dt, ax4, ay4, debug_flag=.true.)
    write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'posn i,j,lon,lat,xi,yj=',i,j,lonn,latn,xi,yj
    write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'posn box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
    call print_berg(stderrunit, berg, 'evolve_iceberg, out of cell at end!')
    bounced=is_point_in_cell(bergs%grd, lonn, latn, i, j, explain=.true.)
    if (debug) call error_mesg('diamonds, evolve_iceberg','berg is out of posn at end!',FATAL)
    write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lon',(i,i=grd%isd,grd%ied)
    do j=grd%jed,grd%jsd,-1
      write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lon(i,j),i=grd%isd,grd%ied)
    enddo
    write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lat',(i,i=grd%isd,grd%ied)
    do j=grd%jed,grd%jsd,-1
      write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lat(i,j),i=grd%isd,grd%ied)
    enddo
  endif

  berg%lon=lonn
  berg%lat=latn
  berg%uvel=uveln
  berg%vvel=vveln
  berg%ine=i
  berg%jne=j
  berg%xi=xi
  berg%yj=yj
 !call interp_flds(grd, i, j, xi, yj, berg%uo, berg%vo, berg%ui, berg%vi, berg%ua, berg%va, berg%ssh_x, berg%ssh_y, berg%sst)

 !if (debug) call print_berg(stderr(), berg, 'evolve_iceberg, final posn.')
  if (debug) call check_position(grd, berg, 'evolve_iceberg (bot)')

  berg=>berg%next
  enddo ! loop over all bergs

  contains

  subroutine rotpos_to_tang(lon, lat, x, y)
  ! Arguments
  real, intent(in) :: lon, lat
  real, intent(out) :: x, y
  ! Local variables
  real :: r,colat,clon,slon

    if (lat>90.) then
      write(stderrunit,*) 'diamonds, rotpos_to_tang: lat>90 already!',lat
      call error_mesg('diamonds, rotpos_to_tang','Something went very wrong!',FATAL)
    endif
    if (lat==90.) then
      write(stderrunit,*) 'diamonds, rotpos_to_tang: lat==90 already!',lat
      call error_mesg('diamonds, rotpos_to_tang','Something went wrong!',FATAL)
    endif

    colat=90.-lat
    r=Rearth*(colat*pi_180)
    clon=cos(lon*pi_180)
    slon=sin(lon*pi_180)
    x=r*clon
    y=r*slon

  end subroutine rotpos_to_tang

  subroutine rotpos_from_tang(x, y, lon, lat)
  ! Arguments
  real, intent(in) :: x, y
  real, intent(out) :: lon, lat
  ! Local variables
  real :: r

    r=sqrt(x**2+y**2)
    lat=90.-(r180_pi*r/Rearth)
    lon=r180_pi*acos(x/r)*sign(1.,y)

  end subroutine rotpos_from_tang

  subroutine rotvec_to_tang(lon, uvel, vvel, xdot, ydot)
  ! Arguments
  real, intent(in) :: lon, uvel, vvel
  real, intent(out) :: xdot, ydot
  ! Local variables
  real :: clon,slon

    clon=cos(lon*pi_180)
    slon=sin(lon*pi_180)
    xdot=-slon*uvel-clon*vvel
    ydot=clon*uvel-slon*vvel

  end subroutine rotvec_to_tang

  subroutine rotvec_from_tang(lon, xdot, ydot, uvel, vvel)
  ! Arguments
  real, intent(in) :: lon, xdot, ydot
  real, intent(out) :: uvel, vvel
  ! Local variables
  real :: clon,slon

    clon=cos(lon*pi_180)
    slon=sin(lon*pi_180)
    uvel=-slon*xdot+clon*ydot
    vvel=-clon*xdot-slon*ydot

  end subroutine rotvec_from_tang

! ##############################################################################

subroutine adjust_index_and_ground(grd, lon, lat, uvel, vvel, i, j, xi, yj, bounced, error)
! Arguments
type(icebergs_gridded), pointer :: grd
real, intent(inout) :: lon, lat, uvel, vvel, xi, yj
integer, intent(inout) :: i,j
logical, intent(out) :: bounced, error
! Local variables
logical lret, lpos
real, parameter :: posn_eps=0.05
integer :: icount, i0, j0, inm, jnm
real :: xi0, yj0, lon0, lat0

  bounced=.false.
  error=.false.
  lon0=lon; lat0=lat ! original position
  i0=i; j0=j ! original i,j
  lret=pos_within_cell(grd, lon, lat, i, j, xi, yj)
  xi0=xi; yj0=yj ! original xi,yj
  if (debug) then
    !Sanity check lret, xi and yj
    lret=is_point_in_cell(grd, lon, lat, i, j)
    if (xi<0. .or. xi>1. .or. yj<0. .or. yj>1.) then
      if (lret) then
        write(stderrunit,*) 'diamonds, adjust: WARNING!!! lret=T but |xi,yj|>1',mpp_pe()
        write(stderrunit,*) 'diamonds, adjust: xi=',xi,' lon=',lon
        write(stderrunit,*) 'diamonds, adjust: x3 x2=',grd%lon(i-1,j),grd%lon(i,j)
        write(stderrunit,*) 'diamonds, adjust: x0 x1=',grd%lon(i-1,j-1),grd%lon(i,j-1)
        write(stderrunit,*) 'diamonds, adjust: yi=',yj,' lat=',lat
        write(stderrunit,*) 'diamonds, adjust: y3 y2=',grd%lat(i-1,j),grd%lat(i,j)
        write(stderrunit,*) 'diamonds, adjust: y0 y1=',grd%lat(i-1,j-1),grd%lat(i,j-1)
        lret=is_point_in_cell(grd, lon, lat, i, j, explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn is_point_in_cell=',lret
        lret=pos_within_cell(grd, lon, lat, i, j, xi, yj, explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn pos_within_cell=',lret
        write(0,*) 'This should never happen!'
        error=.true.; return
     endif
    else
      if (.not.lret) then
        write(stderrunit,*) 'diamonds, adjust: WARNING!!! lret=F but |xi,yj|<1',mpp_pe()
        write(stderrunit,*) 'diamonds, adjust: xi=',xi,' lon=',lon
        write(stderrunit,*) 'diamonds, adjust: x3 x2=',grd%lon(i-1,j),grd%lon(i,j)
        write(stderrunit,*) 'diamonds, adjust: x0 x1=',grd%lon(i-1,j-1),grd%lon(i,j-1)
        write(stderrunit,*) 'diamonds, adjust: yi=',yj,' lat=',lat
        write(stderrunit,*) 'diamonds, adjust: y3 y2=',grd%lat(i-1,j),grd%lat(i,j)
        write(stderrunit,*) 'diamonds, adjust: y0 y1=',grd%lat(i-1,j-1),grd%lat(i,j-1)
        lret=is_point_in_cell(grd, lon, lat, i, j, explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn is_point_in_cell=',lret
        lret=pos_within_cell(grd, lon, lat, i, j, xi, yj, explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn pos_within_cell=',lret
        write(0,*) 'This should never happen!'
        error=.true.; return
      endif
    endif
    lret=pos_within_cell(grd, lon, lat, i, j, xi, yj)
  endif ! debug
  if (lret) return ! Berg was already in cell

  ! Find inm, jnm (as if adjusting i,j) based on xi,yj
  ! ignoring the mand mask.
  ! NOTE:  This search appears to have *NO* active role
  ! in the algorithm other than to flag a warning.
  icount=0
  inm=i0; jnm=j0 ! original i,j
  do while (debug .and. .not.lret .and. icount<4)
    icount=icount+1
    if (xi.lt.0.) then
      if (inm>grd%isd) then
        inm=inm-1
      endif
    elseif (xi.gt.1.) then
      if (inm<grd%ied) then
        inm=inm+1
      endif
    endif
    if (yj.lt.0.) then
      if (jnm>grd%jsd) then
        jnm=jnm-1
      endif
    elseif (yj.gt.1.) then
      if (jnm<grd%jed) then
        jnm=jnm+1
      endif
    endif
    lret=pos_within_cell(grd, lon, lat, inm, jnm, xi, yj) ! Update xi and yj
  enddo
  if (abs(inm-i0)>1) then
    write(stderrunit,*) 'pe=',mpp_pe(),'diamonds, adjust: inm,i0,inm-i0=',inm,i0,inm-i0
   !stop 'Moved too far in i without mask!'
  endif
  if (abs(jnm-j0)>1) then
    write(stderrunit,*) 'pe=',mpp_pe(),'diamonds, adjust: jnm,i0,jnm-j0=',jnm,j0,inm-j0
   !stop 'Moved too far in j without mask!'
  endif

  ! Adjust i,j based on xi,yj while bouncing off of masked land cells
  icount=0
  lret=pos_within_cell(grd, lon, lat, i0, j0, xi, yj)
  do while ( .not.lret.and. icount<4 )
    icount=icount+1
    if (xi.lt.0.) then
      if (i>grd%isd) then
        if (grd%msk(i-1,j)>0.) then
          if (i>grd%isd+1) i=i-1
        else
         !write(stderr(),'(a,6f8.3,i)') 'diamonds, adjust: bouncing berg from west',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
      endif
    elseif (xi.gt.1.) then
      if (i<grd%ied) then
        if (grd%msk(i+1,j)>0.) then
          if (i<grd%ied) i=i+1
        else
         !write(stderr(),'(a,6f8.3,i)') 'diamonds, adjust: bouncing berg from east',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
      endif
    endif
    if (yj.lt.0.) then
      if (j>grd%jsd) then
        if (grd%msk(i,j-1)>0.) then
          if (j>grd%jsd+1) j=j-1
        else
         !write(stderr(),'(a,6f8.3,i)') 'diamonds, adjust: bouncing berg from south',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
      endif
    elseif (yj.gt.1.) then
      if (j<grd%jed) then
        if (grd%msk(i,j+1)>0.) then
          if (j<grd%jed) j=j+1
        else
         !write(stderr(),'(a,6f8.3,i)') 'diamonds, adjust: bouncing berg from north',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
      endif
    endif
    if (bounced) then
      if (xi>1.) xi=1.-posn_eps
      if (xi<0.) xi=posn_eps
      if (yj>1.) yj=1.-posn_eps
      if (yj<0.) yj=posn_eps
      lon=bilin(grd, grd%lon, i, j, xi, yj)
      lat=bilin(grd, grd%lat, i, j, xi, yj)
    endif
    if (debug) then
      if (grd%msk(i,j)==0.) stop 'diamonds, adjust: Berg is in land! This should not happen...'
    endif
    lret=pos_within_cell(grd, lon, lat, i, j, xi, yj) ! Update xi and yj
  enddo

 !if (debug) then
 !  if (abs(i-i0)>2) then
 !    stop 'diamonds, adjust: Moved too far in i!'
 !  endif
 !  if (abs(j-j0)>2) then
 !    stop 'diamonds, adjust: Moved too far in j!'
 !  endif
 !endif

  if (.not.bounced.and.lret.and.grd%msk(i,j)>0.) return ! Landed in ocean without bouncing so all is well
  if (.not.bounced.and..not.lret) then ! This implies the berg traveled many cells without getting far enough
                                       ! OR that it did not move at all (round-off problem)
    if (debug) then
      write(stderrunit,*) 'diamonds, adjust: lon0, lat0=',lon0,lat0
      write(stderrunit,*) 'diamonds, adjust: xi0, yj0=',xi0,yj0
      write(stderrunit,*) 'diamonds, adjust: i0,j0=',i0,j0
      write(stderrunit,*) 'diamonds, adjust: lon, lat=',lon,lat
      write(stderrunit,*) 'diamonds, adjust: xi,yj=',xi,yj
      write(stderrunit,*) 'diamonds, adjust: i,j=',i,j
      write(stderrunit,*) 'diamonds, adjust: inm,jnm=',inm,jnm
      write(stderrunit,*) 'diamonds, adjust: icount=',icount
      lret=pos_within_cell(grd, lon, lat, i, j, xi, yj, explain=.true.)
      write(stderrunit,*) 'diamonds, adjust: lret=',lret
    endif
    if (abs(i-i0)+abs(j-j0)==0) then
      if (use_roundoff_fix) then
        ! This is a special case due to round off where is_point_in_cell()
        ! returns false but xi and yj are between 0 and 1.
        ! It occurs very rarely but often enough to have brought down
        ! ESM2G four times since the spin-up began. (as of 8/10/2010)
        ! This temporary fix arbitrarily moves the berg toward the
        ! center of the current cell.
        xi=(xi-0.5)*(1.-posn_eps)+0.5
        yj=(yj-0.5)*(1.-posn_eps)+0.5
      endif
      call error_mesg('diamonds, adjust', 'Berg did not move or bounce during iterations AND was not in cell. Adjusting!', WARNING)
    else
      call error_mesg('diamonds, adjust', 'Berg iterated many times without bouncing!', WARNING)
    endif
  endif
  if (xi>1.) xi=1.-posn_eps
  if (xi<0.) xi=posn_eps
  if (yj>1.) yj=1.-posn_eps
  if (yj<0.) yj=posn_eps
  lon=bilin(grd, grd%lon, i, j, xi, yj)
  lat=bilin(grd, grd%lat, i, j, xi, yj)
  lret=pos_within_cell(grd, lon, lat, i, j, xi, yj) ! Update xi and yj

  if (.not. lret) then
    write(stderrunit,*) 'diamonds, adjust: Should not get here! Berg is not in cell after adjustment'
    if (debug) error=.true.
  endif

end subroutine adjust_index_and_ground

end subroutine evolve_icebergs

! ##############################################################################

subroutine send_bergs_to_other_pes(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(iceberg), pointer :: kick_the_bucket, this
integer :: nbergs_to_send_e, nbergs_to_send_w
integer :: nbergs_to_send_n, nbergs_to_send_s
integer :: nbergs_rcvd_from_e, nbergs_rcvd_from_w
integer :: nbergs_rcvd_from_n, nbergs_rcvd_from_s
integer, parameter :: buffer_width=18
type(icebergs_gridded), pointer :: grd
integer :: i, nbergs_start, nbergs_end
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>bergs%grd

  if (debug) then
    nbergs_start=count_bergs(bergs)
  endif

  ! Find number of bergs that headed east/west
  nbergs_to_send_e=0
  nbergs_to_send_w=0
  if (associated(bergs%first)) then
    this=>bergs%first
    do while (associated(this))
      if (this%ine.gt.bergs%grd%iec) then
        kick_the_bucket=>this
        this=>this%next
        nbergs_to_send_e=nbergs_to_send_e+1
        call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_e, nbergs_to_send_e)
        call move_trajectory(bergs, kick_the_bucket)
        call delete_iceberg_from_list(bergs%first,kick_the_bucket)
      elseif (this%ine.lt.bergs%grd%isc) then
        kick_the_bucket=>this
        this=>this%next
        nbergs_to_send_w=nbergs_to_send_w+1
        call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_w, nbergs_to_send_w)
        call move_trajectory(bergs, kick_the_bucket)
        call delete_iceberg_from_list(bergs%first,kick_the_bucket)
      else
        this=>this%next
      endif
    enddo
  endif

  ! Send bergs east
  if (grd%pe_E.ne.NULL_PE) then
    call mpp_send(nbergs_to_send_e, plen=1, to_pe=grd%pe_E, tag=COMM_TAG_1)
    if (nbergs_to_send_e.gt.0) then
      call mpp_send(bergs%obuffer_e%data, nbergs_to_send_e*buffer_width, grd%pe_E, tag=COMM_TAG_2)
    endif
  endif

  ! Send bergs west
  if (grd%pe_W.ne.NULL_PE) then
    call mpp_send(nbergs_to_send_w, plen=1, to_pe=grd%pe_W, tag=COMM_TAG_3)
    if (nbergs_to_send_w.gt.0) then
      call mpp_send(bergs%obuffer_w%data, nbergs_to_send_w*buffer_width, grd%pe_W, tag=COMM_TAG_4)
    endif
  endif

  ! Receive bergs from west
  if (grd%pe_W.ne.NULL_PE) then
    nbergs_rcvd_from_w=-999
    call mpp_recv(nbergs_rcvd_from_w, glen=1, from_pe=grd%pe_W, tag=COMM_TAG_1)
    if (nbergs_rcvd_from_w.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_w,' from',grd%pe_W,' (W) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nbergs_rcvd_from_w.gt.0) then
      call increase_ibuffer(bergs%ibuffer_w, nbergs_rcvd_from_w)
      call mpp_recv(bergs%ibuffer_w%data, nbergs_rcvd_from_w*buffer_width, grd%pe_W, tag=COMM_TAG_2)
      do i=1, nbergs_rcvd_from_w
        call unpack_berg_from_buffer2(bergs%first, bergs%ibuffer_w, i)
      enddo
    endif
  else
    nbergs_rcvd_from_w=0
  endif

  ! Receive bergs from east
  if (grd%pe_E.ne.NULL_PE) then
    nbergs_rcvd_from_e=-999
    call mpp_recv(nbergs_rcvd_from_e, glen=1, from_pe=grd%pe_E, tag=COMM_TAG_3)
    if (nbergs_rcvd_from_e.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_e,' from',grd%pe_E,' (E) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nbergs_rcvd_from_e.gt.0) then
      call increase_ibuffer(bergs%ibuffer_e, nbergs_rcvd_from_e)
      call mpp_recv(bergs%ibuffer_e%data, nbergs_rcvd_from_e*buffer_width, grd%pe_E, tag=COMM_TAG_4)
      do i=1, nbergs_rcvd_from_e
        call unpack_berg_from_buffer2(bergs%first, bergs%ibuffer_e, i)
      enddo
    endif
  else
    nbergs_rcvd_from_e=0
  endif

  ! Find number of bergs that headed north/south
  ! (note: this block should technically go ahead of the E/W recv block above
  !  to handle arbitrary orientation of PEs. But for simplicity, it is
  !  here to accomodate diagonal transfer of bergs between PEs -AJA)
  nbergs_to_send_n=0
  nbergs_to_send_s=0
  if (associated(bergs%first)) then
    this=>bergs%first
    do while (associated(this))
      if (this%jne.gt.bergs%grd%jec) then
        kick_the_bucket=>this
        this=>this%next
        nbergs_to_send_n=nbergs_to_send_n+1
        call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_n, nbergs_to_send_n)
        call move_trajectory(bergs, kick_the_bucket)
        call delete_iceberg_from_list(bergs%first,kick_the_bucket)
      elseif (this%jne.lt.bergs%grd%jsc) then
        kick_the_bucket=>this
        this=>this%next
        nbergs_to_send_s=nbergs_to_send_s+1
        call pack_berg_into_buffer2(kick_the_bucket, bergs%obuffer_s, nbergs_to_send_s)
        call move_trajectory(bergs, kick_the_bucket)
        call delete_iceberg_from_list(bergs%first,kick_the_bucket)
      else
        this=>this%next
      endif
    enddo
  endif

  ! Send bergs north
  if (grd%pe_N.ne.NULL_PE) then
    if(folded_north_on_pe) then
       call mpp_send(nbergs_to_send_n, plen=1, to_pe=grd%pe_N, tag=COMM_TAG_9)
    else
       call mpp_send(nbergs_to_send_n, plen=1, to_pe=grd%pe_N, tag=COMM_TAG_5)
    endif
    if (nbergs_to_send_n.gt.0) then
       if(folded_north_on_pe) then
          call mpp_send(bergs%obuffer_n%data, nbergs_to_send_n*buffer_width, grd%pe_N, tag=COMM_TAG_10)
       else
          call mpp_send(bergs%obuffer_n%data, nbergs_to_send_n*buffer_width, grd%pe_N, tag=COMM_TAG_6)
       endif
    endif
  endif

  ! Send bergs south
  if (grd%pe_S.ne.NULL_PE) then
    call mpp_send(nbergs_to_send_s, plen=1, to_pe=grd%pe_S, tag=COMM_TAG_7)
    if (nbergs_to_send_s.gt.0) then
      call mpp_send(bergs%obuffer_s%data, nbergs_to_send_s*buffer_width, grd%pe_S, tag=COMM_TAG_8)
    endif
  endif

  ! Receive bergs from south
  if (grd%pe_S.ne.NULL_PE) then
    nbergs_rcvd_from_s=-999
    call mpp_recv(nbergs_rcvd_from_s, glen=1, from_pe=grd%pe_S, tag=COMM_TAG_5)
    if (nbergs_rcvd_from_s.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_s,' from',grd%pe_S,' (S) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nbergs_rcvd_from_s.gt.0) then
      call increase_ibuffer(bergs%ibuffer_s, nbergs_rcvd_from_s)
      call mpp_recv(bergs%ibuffer_s%data, nbergs_rcvd_from_s*buffer_width, grd%pe_S, tag=COMM_TAG_6)
      do i=1, nbergs_rcvd_from_s
        call unpack_berg_from_buffer2(bergs%first, bergs%ibuffer_s, i)
      enddo
    endif
  else
    nbergs_rcvd_from_s=0
  endif

  ! Receive bergs from north
  if (grd%pe_N.ne.NULL_PE) then
    nbergs_rcvd_from_n=-999
    if(folded_north_on_pe) then
       call mpp_recv(nbergs_rcvd_from_n, glen=1, from_pe=grd%pe_N, tag=COMM_TAG_9)
    else
       call mpp_recv(nbergs_rcvd_from_n, glen=1, from_pe=grd%pe_N, tag=COMM_TAG_7)
    endif
    if (nbergs_rcvd_from_n.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nbergs_rcvd_from_n,' from',grd%pe_N,' (N) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nbergs_rcvd_from_n.gt.0) then
      call increase_ibuffer(bergs%ibuffer_n, nbergs_rcvd_from_n)
      if(folded_north_on_pe) then
         call mpp_recv(bergs%ibuffer_n%data, nbergs_rcvd_from_n*buffer_width, grd%pe_N, tag=COMM_TAG_10)
      else
         call mpp_recv(bergs%ibuffer_n%data, nbergs_rcvd_from_n*buffer_width, grd%pe_N, tag=COMM_TAG_8)
      endif
      do i=1, nbergs_rcvd_from_n
        call unpack_berg_from_buffer2(bergs%first, bergs%ibuffer_n, i)
      enddo
    endif
  else
    nbergs_rcvd_from_n=0
  endif

  if (debug) then
    nbergs_end=count_bergs(bergs)
    i=nbergs_rcvd_from_n+nbergs_rcvd_from_s+nbergs_rcvd_from_e+nbergs_rcvd_from_w &
     -nbergs_to_send_n-nbergs_to_send_s-nbergs_to_send_e-nbergs_to_send_w
    if (nbergs_end-(nbergs_start+i).ne.0) then
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_end=',nbergs_end,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_start=',nbergs_start,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: delta=',i,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: error=',nbergs_end-(nbergs_start+i),' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_to_send_n=',nbergs_to_send_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_to_send_s=',nbergs_to_send_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_to_send_e=',nbergs_to_send_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_to_send_w=',nbergs_to_send_w,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_rcvd_from_n=',nbergs_rcvd_from_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_rcvd_from_s=',nbergs_rcvd_from_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_rcvd_from_e=',nbergs_rcvd_from_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_bergs_to_other_pes: nbergs_rcvd_from_w=',nbergs_rcvd_from_w,' on PE',mpp_pe()
      call error_mesg('diamonds, send_bergs_to_other_pes:', 'We lost some bergs!', FATAL)
    endif
  endif

  if (debug) then
    i=0
    this=>bergs%first
    do while (associated(this))
      call check_position(grd, this, 'exchange (bot)')
      if (this%ine.lt.bergs%grd%isc .or. &
          this%ine.gt.bergs%grd%iec .or. &
          this%jne.lt.bergs%grd%jsc .or. &
          this%jne.gt.bergs%grd%jec) i=i+1
      this=>this%next
    enddo ! while
    call mpp_sum(i)
    if (i>0 .and. mpp_pe()==mpp_root_pe()) then
      write(stderrunit,'(a,i4)') 'diamonds, send_bergs_to_other_pes: # of bergs outside computational domain = ',i
      call error_mesg('diamonds, send_bergs_to_other_pes:', 'there are bergs still in halos!', FATAL)
    endif ! root_pe
  endif ! debug

  call mpp_sync_self()

contains

  subroutine pack_berg_into_buffer2(berg, buff, n)
  ! Arguments
  type(iceberg), pointer :: berg
  type(buffer), pointer :: buff
  integer, intent(in) :: n
  ! Local variables

    if (.not.associated(buff)) call increase_buffer(buff,delta_buf)
    if (n>buff%size) call increase_buffer(buff,delta_buf)

    buff%data(1,n)=berg%lon
    buff%data(2,n)=berg%lat
    buff%data(3,n)=berg%uvel
    buff%data(4,n)=berg%vvel
    buff%data(5,n)=berg%xi
    buff%data(6,n)=berg%yj
    buff%data(7,n)=berg%start_lon
    buff%data(8,n)=berg%start_lat
    buff%data(9,n)=float(berg%start_year)
    buff%data(10,n)=berg%start_day
    buff%data(11,n)=berg%start_mass
    buff%data(12,n)=berg%mass
    buff%data(13,n)=berg%thickness
    buff%data(14,n)=berg%width
    buff%data(15,n)=berg%length
    buff%data(16,n)=berg%mass_scaling
    buff%data(17,n)=berg%mass_of_bits
    buff%data(18,n)=berg%heat_density

  end subroutine pack_berg_into_buffer2

  subroutine increase_buffer(old,delta)
  ! Arguments
  type(buffer), pointer :: old
  integer, intent(in) :: delta
  ! Local variables
  type(buffer), pointer :: new
  integer :: new_size

    if (.not.associated(old)) then
      new_size=delta
    else
      new_size=old%size+delta
    endif
    allocate(new)
    allocate(new%data(buffer_width,new_size))
    new%size=new_size
    if (associated(old)) then
      new%data(:,1:old%size)=old%data(:,1:old%size)
      deallocate(old%data)
      deallocate(old)
    endif
    old=>new
   !write(stderr(),*) 'diamonds, increase_buffer',mpp_pe(),' increased to',new_size

  end subroutine increase_buffer

  subroutine unpack_berg_from_buffer2(first, buff, n)
  ! Arguments
  type(iceberg), pointer :: first
  type(buffer), pointer :: buff
  integer, intent(in) :: n
  ! Local variables
 !real :: lon, lat, uvel, vvel, xi, yj
 !real :: start_lon, start_lat, start_day, start_mass
 !integer :: ine, jne, start_year
  logical :: lres
  type(iceberg) :: localberg

    localberg%lon=buff%data(1,n)
    localberg%lat=buff%data(2,n)
    localberg%uvel=buff%data(3,n)
    localberg%vvel=buff%data(4,n)
    localberg%xi=buff%data(5,n)
    localberg%yj=buff%data(6,n)
    localberg%start_lon=buff%data(7,n)
    localberg%start_lat=buff%data(8,n)
    localberg%start_year=nint(buff%data(9,n))
    localberg%start_day=buff%data(10,n)
    localberg%start_mass=buff%data(11,n)
    localberg%mass=buff%data(12,n)
    localberg%thickness=buff%data(13,n)
    localberg%width=buff%data(14,n)
    localberg%length=buff%data(15,n)
    localberg%mass_scaling=buff%data(16,n)
    localberg%mass_of_bits=buff%data(17,n)
    localberg%heat_density=buff%data(18,n)
    lres=find_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
    if (lres) then
      lres=pos_within_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne, localberg%xi, localberg%yj)
      call add_new_berg_to_list(first, localberg)
    else
      lres=find_cell_wide(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
      if (lres) then
        lres=pos_within_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne, localberg%xi, localberg%yj)
        call add_new_berg_to_list(first, localberg)
      else
        write(stderrunit,'("diamonds, unpack_berg_from_buffer pe=(",i3,a,2i4,a,2f8.2)')&
             & mpp_pe(),') Failed to find i,j=',localberg%ine,localberg%jne,' for lon,lat=',localberg%lon,localberg%lat
        write(stderrunit,*) localberg%lon,localberg%lat
        write(stderrunit,*) localberg%uvel,localberg%vvel
        write(stderrunit,*) grd%isc,grd%iec,grd%jsc,grd%jec
        write(stderrunit,*) grd%isd,grd%ied,grd%jsd,grd%jed
        write(stderrunit,*) grd%lon(grd%isc-1,grd%jsc-1),grd%lon(grd%iec,grd%jsc)
        write(stderrunit,*) grd%lat(grd%isc-1,grd%jsc-1),grd%lat(grd%iec,grd%jec)
        write(stderrunit,*) grd%lon(grd%isd,grd%jsd),grd%lon(grd%ied,grd%jsd)
        write(stderrunit,*) grd%lat(grd%isd,grd%jsd),grd%lat(grd%ied,grd%jed)
        write(stderrunit,*) lres
        call error_mesg('diamonds, unpack_berg_from_buffer', 'can not find a cell to place berg in!', FATAL)
      endif
    endif

  end subroutine unpack_berg_from_buffer2

  subroutine increase_ibuffer(old,delta)
  ! Arguments
  type(buffer), pointer :: old
  integer, intent(in) :: delta
  ! Local variables
  type(buffer), pointer :: new
  integer :: new_size, old_size

    if (.not.associated(old)) then
      new_size=delta+delta_buf
      old_size=0
    else
      old_size=old%size
      if (delta<old%size) then
        new_size=old%size+delta
      else
        new_size=delta+delta_buf
      endif
    endif

    if (old_size.ne.new_size) then
      allocate(new)
      allocate(new%data(buffer_width,new_size))
      new%size=new_size
      if (associated(old)) then
        new%data(:,1:old%size)=old%data(:,1:old%size)
        deallocate(old%data)
        deallocate(old)
      endif
      old=>new
     !write(stderr(),*) 'diamonds, increase_ibuffer',mpp_pe(),' increased to',new_size
    endif

  end subroutine increase_ibuffer

end subroutine send_bergs_to_other_pes

! ##############################################################################

subroutine icebergs_init(bergs, &
             gni, gnj, layout, io_layout, axes, maskmap, x_cyclic, tripolar_grid, &
             dt, Time, ice_lon, ice_lat, ice_wet, ice_dx, ice_dy, ice_area, &
             cos_rot, sin_rot)
! Arguments
type(icebergs), pointer :: bergs
integer, intent(in) :: gni, gnj, layout(2), io_layout(2), axes(2)
logical, intent(in), optional :: maskmap(:,:)
logical, intent(in) :: x_cyclic, tripolar_grid
real, intent(in) :: dt
type (time_type), intent(in) :: Time ! current time
real, dimension(:,:), intent(in) :: ice_lon, ice_lat, ice_wet
real, dimension(:,:), intent(in) :: ice_dx, ice_dy, ice_area
real, dimension(:,:), intent(in) :: cos_rot, sin_rot
! Namelist parameters (and defaults)
integer :: halo=4 ! Width of halo region
integer :: traj_sample_hrs=24 ! Period between sampling of position for trajectory storage
integer :: verbose_hrs=24 ! Period between verbose messages
real :: rho_bergs=850. ! Density of icebergs
real :: LoW_ratio=1.5 ! Initial ratio L/W for newly calved icebergs
real :: bergy_bit_erosion_fraction=0. ! Fraction of erosion melt flux to divert to bergy bits
real :: sicn_shift=0. ! Shift of sea-ice concentration in erosion flux modulation (0<sicn_shift<1)
logical :: use_operator_splitting=.true. ! Use first order operator splitting for thermodynamics
logical :: add_weight_to_ocean=.true. ! Add weight of icebergs + bits to ocean
logical :: passive_mode=.false. ! Add weight of icebergs + bits to ocean
logical :: time_average_weight=.false. ! Time average the weight on the ocean
logical :: fix_restart_dates=.true. ! After a restart, check that bergs were created before the current model date
logical :: reproduce_siena=.false. !To reproduce siena answers which change across PE layout change set to .true.
real :: speed_limit=0. ! CFL speed limit for a berg
real, dimension(nclasses) :: initial_mass=(/8.8e7, 4.1e8, 3.3e9, 1.8e10, 3.8e10, 7.5e10, 1.2e11, 2.2e11, 3.9e11, 7.4e11/) ! Mass thresholds between iceberg classes (kg)
real, dimension(nclasses) :: distribution=(/0.24, 0.12, 0.15, 0.18, 0.12, 0.07, 0.03, 0.03, 0.03, 0.02/) ! Fraction of calving to apply to this class (non-dim)
real, dimension(nclasses) :: mass_scaling=(/2000, 200, 50, 20, 10, 5, 2, 1, 1, 1/) ! Ratio between effective and real iceberg mass (non-dim)
real, dimension(nclasses) :: initial_thickness=(/40., 67., 133., 175., 250., 250., 250., 250., 250., 250./) ! Total thickness of newly calved bergs (m)
namelist /icebergs_nml/ verbose, budget, halo, traj_sample_hrs, initial_mass, &
         distribution, mass_scaling, initial_thickness, verbose_hrs, &
         rho_bergs, LoW_ratio, debug, really_debug, use_operator_splitting, bergy_bit_erosion_fraction, &
         parallel_reprod, use_slow_find, sicn_shift, add_weight_to_ocean, passive_mode, ignore_ij_restart, &
         time_average_weight, generate_test_icebergs, speed_limit, fix_restart_dates, use_roundoff_fix, &
         old_bug_rotated_weights, make_calving_reproduce, restart_input_dir, reproduce_siena
! Local variables
integer :: ierr, iunit, i, j, id_class, axes3d(3), is,ie,js,je
type(icebergs_gridded), pointer :: grd
real :: minl
logical :: lerr
integer :: stdlogunit, stderrunit

  ! Get the stderr and stdlog unit numbers
  stderrunit=stderr()
  stdlogunit=stdlog()

! Read namelist parameters
 !write(stderrunit,*) 'diamonds: reading namelist'
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=icebergs_nml, iostat=ierr)
#else
  iunit = open_namelist_file()
  read  (iunit, icebergs_nml,iostat=ierr)
  call close_file (iunit)
#endif
  ierr = check_nml_error(ierr,'icebergs_nml')

  if (really_debug) debug=.true. ! One implies the other...

! Log version and parameters
  call write_version_number(version, tagname)
  write (stdlogunit, icebergs_nml)

  if( reproduce_siena ) then
     if( mpp_pe() == mpp_root_pe() ) then
        call error_mesg("ice_bergs: You have overridden the default value of reproduce_siena " // &
                        "and set it to .true. in icebergs_nml. This is a temporary workaround to " // &
                        "allow for consistency in continuing experiments.",  "Please use the default " //&
                        "value (.false.) as this option will be removed in a future release. ", WARNING)
     endif
  endif

! Allocate overall structure
 !write(stderrunit,*) 'diamonds: allocating bergs'
  allocate(bergs)
  allocate(bergs%grd)
  grd=>bergs%grd ! For convenience to avoid bergs%grd%X
 !write(stderrunit,*) 'diamonds: allocating domain'
  allocate(grd%domain)

! Clocks
  bergs%clock=mpp_clock_id( 'Icebergs', flags=clock_flag_default, grain=CLOCK_COMPONENT )
  bergs%clock_mom=mpp_clock_id( 'Icebergs-momentum', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_the=mpp_clock_id( 'Icebergs-thermodyn', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_int=mpp_clock_id( 'Icebergs-interface', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_cal=mpp_clock_id( 'Icebergs-calving', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_com=mpp_clock_id( 'Icebergs-communication', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_ini=mpp_clock_id( 'Icebergs-initialization', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_ior=mpp_clock_id( 'Icebergs-I/O read', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_iow=mpp_clock_id( 'Icebergs-I/O write', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  bergs%clock_dia=mpp_clock_id( 'Icebergs-diagnostics', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  call mpp_clock_begin(bergs%clock)
  call mpp_clock_begin(bergs%clock_ini)

! Set up iceberg domain
 !write(stderrunit,*) 'diamonds: defining domain'
  if(tripolar_grid) then
    call mpp_define_domains( (/1,gni,1,gnj/), layout, grd%domain, &
!                            maskmap=maskmap, &
                             xflags=CYCLIC_GLOBAL_DOMAIN, xhalo=halo,  &
                             yflags=FOLD_NORTH_EDGE, yhalo=halo, name='diamond')
  else if(x_cyclic) then
    call mpp_define_domains( (/1,gni,1,gnj/), layout, grd%domain, &
!                            maskmap=maskmap, &
                             xflags=CYCLIC_GLOBAL_DOMAIN, &
                             xhalo=halo, yhalo=halo, name='diamond')
  else
    call mpp_define_domains( (/1,gni,1,gnj/), layout, grd%domain, &
!                            maskmap=maskmap, &
                             xhalo=halo, yhalo=halo, name='diamond')
  endif

  call mpp_define_io_domain(grd%domain, io_layout)

 !write(stderrunit,*) 'diamond: get compute domain'
  call mpp_get_compute_domain( grd%domain, grd%isc, grd%iec, grd%jsc, grd%jec )
  call mpp_get_data_domain( grd%domain, grd%isd, grd%ied, grd%jsd, grd%jed )

  call mpp_get_neighbor_pe(grd%domain, NORTH, grd%pe_N)
  call mpp_get_neighbor_pe(grd%domain, SOUTH, grd%pe_S)
  call mpp_get_neighbor_pe(grd%domain, EAST, grd%pe_E)
  call mpp_get_neighbor_pe(grd%domain, WEST, grd%pe_W)

  folded_north_on_pe = .false.
  if(tripolar_grid .and. grd%jec == gnj) folded_north_on_pe = .true.
 !write(stderrunit,'(a,6i4)') 'diamonds, icebergs_init: pe,n,s,e,w =',mpp_pe(),grd%pe_N,grd%pe_S,grd%pe_E,grd%pe_W, NULL_PE

 !if (verbose) &
 !write(stderrunit,'(a,i3,a,4i4,a,4f8.2)') 'diamonds, icebergs_init: (',mpp_pe(),') [ij][se]c=', &
 !     grd%isc,grd%iec,grd%jsc,grd%jec, &
 !     ' [lon|lat][min|max]=', minval(ice_lon),maxval(ice_lon),minval(ice_lat),maxval(ice_lat)
 !write(stderrunit,*) 'diamonds, int args = ', mpp_pe(),gni, gnj, layout, axes


 !write(stderrunit,*) 'diamonds: allocating grid'
  allocate( grd%lon(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%lon(:,:)=999.
  allocate( grd%lat(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%lat(:,:)=999.
  allocate( grd%lonc(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%lon(:,:)=999.
  allocate( grd%latc(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%lat(:,:)=999.
  allocate( grd%dx(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%dx(:,:)=0.
  allocate( grd%dy(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%dy(:,:)=0.
  allocate( grd%area(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%area(:,:)=0.
  allocate( grd%msk(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%msk(:,:)=0.
  allocate( grd%cos(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%cos(:,:)=1.
  allocate( grd%sin(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%sin(:,:)=0.
  allocate( grd%calving(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%calving(:,:)=0.
  allocate( grd%calving_hflx(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%calving_hflx(:,:)=0.
  allocate( grd%stored_heat(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%stored_heat(:,:)=0.
  allocate( grd%floating_melt(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%floating_melt(:,:)=0.
  allocate( grd%berg_melt(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%berg_melt(:,:)=0.
  allocate( grd%melt_buoy(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%melt_buoy(:,:)=0.
  allocate( grd%melt_eros(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%melt_eros(:,:)=0.
  allocate( grd%melt_conv(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%melt_conv(:,:)=0.
  allocate( grd%bergy_src(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%bergy_src(:,:)=0.
  allocate( grd%bergy_melt(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%bergy_melt(:,:)=0.
  allocate( grd%bergy_mass(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%bergy_mass(:,:)=0.
  allocate( grd%virtual_area(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%virtual_area(:,:)=0.
  allocate( grd%mass(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%mass(:,:)=0.
  allocate( grd%mass_on_ocean(grd%isd:grd%ied, grd%jsd:grd%jed, 9) ); grd%mass_on_ocean(:,:,:)=0.
  allocate( grd%stored_ice(grd%isd:grd%ied, grd%jsd:grd%jed, nclasses) ); grd%stored_ice(:,:,:)=0.
  allocate( grd%real_calving(grd%isd:grd%ied, grd%jsd:grd%jed, nclasses) ); grd%real_calving(:,:,:)=0.
  allocate( grd%uo(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%uo(:,:)=0.
  allocate( grd%vo(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%vo(:,:)=0.
  allocate( grd%ui(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%ui(:,:)=0.
  allocate( grd%vi(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%vi(:,:)=0.
  allocate( grd%ua(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%ua(:,:)=0.
  allocate( grd%va(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%va(:,:)=0.
  allocate( grd%ssh(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%ssh(:,:)=0.
  allocate( grd%sst(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%sst(:,:)=0.
  allocate( grd%cn(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%cn(:,:)=0.
  allocate( grd%hi(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%hi(:,:)=0.
  allocate( grd%tmp(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%tmp(:,:)=0.
  allocate( grd%tmpc(grd%isc:grd%iec, grd%jsc:grd%jec) ); grd%tmpc(:,:)=0.
  allocate( bergs%nbergs_calved_by_class(nclasses) ); bergs%nbergs_calved_by_class(:)=0
  allocate( grd%parity_x(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%parity_x(:,:)=1.
  allocate( grd%parity_y(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%parity_y(:,:)=1.

 !write(stderrunit,*) 'diamonds: copying grid'
  ! Copy data declared on ice model computational domain
  is=grd%isc; ie=grd%iec; js=grd%jsc; je=grd%jec
  grd%lon(is:ie,js:je)=ice_lon(:,:)
  grd%lat(is:ie,js:je)=ice_lat(:,:)
  grd%area(is:ie,js:je)=ice_area(:,:)*(4.*pi*radius*radius)
  ! Copy data declared on ice model data domain
  is=grd%isc-1; ie=grd%iec+1; js=grd%jsc-1; je=grd%jec+1
  grd%dx(is:ie,js:je)=ice_dx(:,:)
  grd%dy(is:ie,js:je)=ice_dy(:,:)
  grd%msk(is:ie,js:je)=ice_wet(:,:)
  grd%cos(is:ie,js:je)=cos_rot(:,:)
  grd%sin(is:ie,js:je)=sin_rot(:,:)

  call mpp_update_domains(grd%lon, grd%domain, position=CORNER)
  call mpp_update_domains(grd%lat, grd%domain, position=CORNER)
  call mpp_update_domains(grd%dy, grd%dx, grd%domain, gridtype=CGRID_NE, flags=SCALAR_PAIR)
  call mpp_update_domains(grd%area, grd%domain)
  call mpp_update_domains(grd%msk, grd%domain)
  call mpp_update_domains(grd%cos, grd%domain, position=CORNER)
  call mpp_update_domains(grd%sin, grd%domain, position=CORNER)
  call mpp_update_domains(grd%parity_x, grd%parity_y, grd%domain, gridtype=AGRID) ! If either parity_x/y is -ve, we need rotation of vectors

  ! Sanitize lon and lat at the SW edges
  do j=grd%jsc-1,grd%jsd,-1; do i=grd%isd,grd%ied
      if (grd%lon(i,j).gt.900.) grd%lon(i,j)=grd%lon(i,j+1)
      if (grd%lat(i,j).gt.900.) grd%lat(i,j)=2.*grd%lat(i,j+1)-grd%lat(i,j+2)
  enddo; enddo

  do j=grd%jsd,grd%jed; do i=grd%isd,grd%ied
      if (grd%lon(i,j).gt.900.) write(stderrunit,*) 'bad lon: ',mpp_pe(),i-grd%isc+1,j-grd%jsc+1
      if (grd%lat(i,j).gt.900.) write(stderrunit,*) 'bad lat: ',mpp_pe(),i-grd%isc+1,j-grd%jsc+1
  enddo; enddo

  ! Sanitize lon for the tile (need continuous longitudes within one tile)
  if(reproduce_siena) then
  j=grd%jsc; do i=grd%isc+1,grd%ied
    minl=grd%lon(i-1,j)-180.
    grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
  enddo
  j=grd%jsc; do i=grd%isc-1,grd%isd,-1
    minl=grd%lon(i+1,j)-180.
    grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
  enddo
  do j=grd%jsc+1,grd%jed; do i=grd%isd,grd%ied
      minl=grd%lon(i,j-1)-180.
      grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
  enddo; enddo
  do j=grd%jsc-1,grd%jsd,-1; do i=grd%isd,grd%ied
      minl=grd%lon(i,j+1)-180.
      grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
  enddo; enddo
  else  !The fix to reproduce across PE layout change, from AJA
  j=grd%jsc; do i=grd%isc+1,grd%ied
    minl=grd%lon(i-1,j)-180.
    if (abs(grd%lon(i,j)-(modulo(grd%lon(i,j)-minl,360.)+minl))>180.) &
       grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
  enddo
  j=grd%jsc; do i=grd%isc-1,grd%isd,-1
    minl=grd%lon(i+1,j)-180.
    if (abs(grd%lon(i,j)-(modulo(grd%lon(i,j)-minl,360.)+minl))>180.) &
       grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
  enddo
  do j=grd%jsc+1,grd%jed; do i=grd%isd,grd%ied
      minl=grd%lon(i,j-1)-180.
      if (abs(grd%lon(i,j)-(modulo(grd%lon(i,j)-minl,360.)+minl))>180.) &
         grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
  enddo; enddo
  do j=grd%jsc-1,grd%jsd,-1; do i=grd%isd,grd%ied
      minl=grd%lon(i,j+1)-180.
      if (abs(grd%lon(i,j)-(modulo(grd%lon(i,j)-minl,360.)+minl))>180.) &
         grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
  enddo; enddo
  endif
  ! lonc, latc used for searches
  do j=grd%jsd+1,grd%jed; do i=grd%isd+1,grd%ied
    grd%lonc(i,j)=0.25*( (grd%lon(i,j)+grd%lon(i-1,j-1)) &
                        +(grd%lon(i-1,j)+grd%lon(i,j-1)) )
    grd%latc(i,j)=0.25*( (grd%lat(i,j)+grd%lat(i-1,j-1)) &
                        +(grd%lat(i-1,j)+grd%lat(i,j-1)) )
  enddo; enddo

  if (debug) then
    write(stderrunit,'(a,i3,a,4i4,a,4f8.2)') 'diamonds, icebergs_init: (',mpp_pe(),') [ij][se]c=', &
         grd%isc,grd%iec,grd%jsc,grd%jec, &
         ' [lon|lat][min|max]=', minval(grd%lon),maxval(grd%lon),minval(grd%lat),maxval(grd%lat)
  endif

 !if (mpp_pe().eq.3) then
 !  write(stderrunit,'(a3,32i7)') 'Lon',(i,i=grd%isd,grd%ied)
 !  do j=grd%jed,grd%jsd,-1
 !    write(stderrunit,'(i3,32f7.1)') j,(grd%lon(i,j),i=grd%isd,grd%ied)
 !  enddo
 !  write(stderrunit,'(a3,32i7)') 'Lat',(i,i=grd%isd,grd%ied)
 !  do j=grd%jed,grd%jsd,-1
 !    write(stderrunit,'(i3,32f7.1)') j,(grd%lat(i,j),i=grd%isd,grd%ied)
 !  enddo
 !  write(stderrunit,'(a3,32i7)') 'Msk',(i,i=grd%isd,grd%ied)
 !  do j=grd%jed,grd%jsd,-1
 !    write(stderrunit,'(i3,32f7.1)') j,(grd%msk(i,j),i=grd%isd,grd%ied)
 !  enddo
 !endif

 ! Parameters
  bergs%dt=dt
  bergs%traj_sample_hrs=traj_sample_hrs
  bergs%verbose_hrs=verbose_hrs
  bergs%grd%halo=halo
  bergs%rho_bergs=rho_bergs
  bergs%LoW_ratio=LoW_ratio
  bergs%use_operator_splitting=use_operator_splitting
  bergs%bergy_bit_erosion_fraction=bergy_bit_erosion_fraction
  bergs%sicn_shift=sicn_shift
  bergs%passive_mode=passive_mode
  bergs%time_average_weight=time_average_weight
  bergs%speed_limit=speed_limit
  bergs%add_weight_to_ocean=add_weight_to_ocean
  allocate( bergs%initial_mass(nclasses) ); bergs%initial_mass(:)=initial_mass(:)
  allocate( bergs%distribution(nclasses) ); bergs%distribution(:)=distribution(:)
  allocate( bergs%mass_scaling(nclasses) ); bergs%mass_scaling(:)=mass_scaling(:)
  allocate( bergs%initial_thickness(nclasses) ); bergs%initial_thickness(:)=initial_thickness(:)
  allocate( bergs%initial_width(nclasses) )
  allocate( bergs%initial_length(nclasses) )
  bergs%initial_width(:)=sqrt(initial_mass(:)/(LoW_ratio*rho_bergs*initial_thickness(:)))
  bergs%initial_length(:)=LoW_ratio*bergs%initial_width(:)

  call mpp_clock_end(bergs%clock_ini)
  call mpp_clock_begin(bergs%clock_ior)
  call read_restart_bergs(bergs,Time)
  call bergs_chksum(bergs, 'read_restart bergs')
  if (fix_restart_dates) call offset_berg_dates(bergs,Time)
  call read_restart_calving(bergs)
  call mpp_clock_end(bergs%clock_ior)
  call mpp_clock_begin(bergs%clock_ini)

  if (really_debug) call print_bergs(stderrunit,bergs,'icebergs_init, initial status')

  ! Diagnostics
  id_class = diag_axis_init('mass_class', initial_mass, 'kg','Z', 'iceberg mass')
  axes3d(1:2)=axes
  axes3d(3)=id_class
  grd%id_calving=register_diag_field('icebergs', 'calving', axes, Time, &
     'Incoming Calving mass rate', 'kg/s')
  grd%id_calving_hflx_in=register_diag_field('icebergs', 'calving_hflx_in', axes, Time, &
     'Incoming Calving heat flux', 'J/s')
  grd%id_accum=register_diag_field('icebergs', 'accum_calving', axes, Time, &
     'Accumulated calving mass rate', 'kg/s')
  grd%id_unused=register_diag_field('icebergs', 'unused_calving', axes, Time, &
     'Unused calving mass rate', 'kg/s')
  grd%id_floating_melt=register_diag_field('icebergs', 'melt', axes, Time, &
     'Melt rate of icebergs + bits', 'kg/(m^2*s)')
  grd%id_berg_melt=register_diag_field('icebergs', 'berg_melt', axes, Time, &
     'Melt rate of icebergs', 'kg/(m^2*s)')
  grd%id_melt_buoy=register_diag_field('icebergs', 'melt_buoy', axes, Time, &
     'Buoyancy component of iceberg melt rate', 'kg/(m^2*s)')
  grd%id_melt_eros=register_diag_field('icebergs', 'melt_eros', axes, Time, &
     'Erosion component of iceberg melt rate', 'kg/(m^2*s)')
  grd%id_melt_conv=register_diag_field('icebergs', 'melt_conv', axes, Time, &
     'Convective component of iceberg melt rate', 'kg/(m^2*s)')
  grd%id_bergy_src=register_diag_field('icebergs', 'bergy_src', axes, Time, &
     'Mass source of bergy bits', 'kg/(m^2*s)')
  grd%id_bergy_melt=register_diag_field('icebergs', 'bergy_melt', axes, Time, &
     'Melt rate of bergy bits', 'kg/(m^2*s)')
  grd%id_bergy_mass=register_diag_field('icebergs', 'bergy_mass', axes, Time, &
     'Bergy bit density field', 'kg/(m^2)')
  grd%id_virtual_area=register_diag_field('icebergs', 'virtual_area', axes, Time, &
     'Virtual coverage by icebergs', 'm^2')
  grd%id_mass=register_diag_field('icebergs', 'mass', axes, Time, &
     'Iceberg density field', 'kg/(m^2)')
  grd%id_stored_ice=register_diag_field('icebergs', 'stored_ice', axes3d, Time, &
     'Accumulated ice mass by class', 'kg')
  grd%id_real_calving=register_diag_field('icebergs', 'real_calving', axes3d, Time, &
     'Calving into iceberg class', 'kg/s')
  grd%id_uo=register_diag_field('icebergs', 'uo', axes, Time, &
     'Ocean zonal component of velocity', 'm s^-1')
  grd%id_vo=register_diag_field('icebergs', 'vo', axes, Time, &
     'Ocean meridional component of velocity', 'm s^-1')
  grd%id_ui=register_diag_field('icebergs', 'ui', axes, Time, &
     'Ice zonal component of velocity', 'm s^-1')
  grd%id_vi=register_diag_field('icebergs', 'vi', axes, Time, &
     'Ice meridional component of velocity', 'm s^-1')
  grd%id_ua=register_diag_field('icebergs', 'ua', axes, Time, &
     'Atmos zonal component of velocity', 'm s^-1')
  grd%id_va=register_diag_field('icebergs', 'va', axes, Time, &
     'Atmos meridional component of velocity', 'm s^-1')
  grd%id_sst=register_diag_field('icebergs', 'sst', axes, Time, &
     'Sea surface temperature', 'degrees_C')
  grd%id_cn=register_diag_field('icebergs', 'cn', axes, Time, &
     'Sea ice concentration', '(fraction)')
  grd%id_hi=register_diag_field('icebergs', 'hi', axes, Time, &
     'Sea ice thickness', 'm')

  ! Static fields
  id_class=register_static_field('icebergs', 'lon', axes, &
               'longitude (corners)', 'degrees_E')
  if (id_class>0) lerr=send_data(id_class, grd%lon(grd%isc:grd%iec,grd%jsc:grd%jec))
  id_class=register_static_field('icebergs', 'lat', axes, &
               'latitude (corners)', 'degrees_N')
  if (id_class>0) lerr=send_data(id_class, grd%lat(grd%isc:grd%iec,grd%jsc:grd%jec))
  id_class=register_static_field('icebergs', 'area', axes, &
               'cell area', 'm^2')
  if (id_class>0) lerr=send_data(id_class, grd%area(grd%isc:grd%iec,grd%jsc:grd%jec))
  id_class=register_static_field('icebergs', 'mask', axes, &
               'wet point mask', 'none')
  if (id_class>0) lerr=send_data(id_class, grd%msk(grd%isc:grd%iec,grd%jsc:grd%jec))

  if (debug) then
    call grd_chksum2(grd, grd%lon, 'init lon')
    call grd_chksum2(grd, grd%lat, 'init lat')
    call grd_chksum2(grd, grd%lonc, 'init lonc')
    call grd_chksum2(grd, grd%latc, 'init latc')
    call grd_chksum2(grd, grd%area, 'init area')
    call grd_chksum2(grd, grd%msk, 'init msk')
    call grd_chksum2(grd, grd%cos, 'init cos')
    call grd_chksum2(grd, grd%sin, 'init sin')
  endif

 !write(stderrunit,*) 'diamonds: done'
  call mpp_clock_end(bergs%clock_ini)
  call mpp_clock_end(bergs%clock)
end subroutine icebergs_init

! ##############################################################################

subroutine read_restart_bergs(bergs,Time)
! Arguments
type(icebergs), pointer :: bergs
type(time_type), intent(in) :: Time
! Local variables
integer, dimension(:), allocatable :: found_restart_int
integer :: k, ierr, ncid, dimid, nbergs_in_file
integer :: lonid, latid, uvelid, vvelid, ineid, jneid
integer :: massid, thicknessid, widthid, lengthid
integer :: start_lonid, start_latid, start_yearid, start_dayid, start_massid
integer :: scaling_id, mass_of_bits_id, heat_density_id
logical :: lres, found_restart, multiPErestart
real :: lon0, lon1, lat0, lat1
character(len=33) :: filename, filename_base
type(icebergs_gridded), pointer :: grd
type(iceberg) :: localberg ! NOT a pointer but an actual local variable
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  ! For convenience
  grd=>bergs%grd

  ! Find a restart file
  multiPErestart=.false.

  ! Zero out nbergs_in_file
  nbergs_in_file = 0

  filename_base=trim(restart_input_dir)//'icebergs.res.nc'

  found_restart = find_restart_file(filename_base, filename, multiPErestart)

  ! Check if no restart found on any pe
  allocate(found_restart_int(mpp_npes()))
  if (found_restart .eqv. .true.) then
     k=1
  else
     k=0
  endif
  call mpp_gather((/k/),found_restart_int)
  if (sum(found_restart_int)==0.and.mpp_pe()==mpp_root_pe())&
       & write(*,'(a)') 'diamonds, read_restart_bergs: no restart file found'
  deallocate(found_restart_int)

  if (.not.found_restart) then

  multiPErestart=.true. ! This is to force sanity checking in a mulit-PE mode if no file was found on this PE

  elseif (found_restart) then ! if (.not.found_restart)
  ! only do the following if a file was found

  if (verbose.and.mpp_pe()==mpp_root_pe()) write(*,'(2a)') 'diamonds, read_restart_bergs: found restart file = ',filename

  ierr=nf_open(filename, NF_NOWRITE, ncid)
  if (ierr .ne. NF_NOERR) write(stderrunit,*) 'diamonds, read_restart_bergs: nf_open failed'

  ierr=nf_inq_unlimdim(ncid, dimid)
  if (ierr .ne. NF_NOERR) write(stderrunit,*) 'diamonds, read_restart_bergs: nf_inq_unlimdim failed'

  ierr=nf_inq_dimlen(ncid, dimid, nbergs_in_file)
  if (ierr .ne. NF_NOERR) write(stderrunit,*) 'diamonds, read_restart_bergs: nf_inq_dimlen failed'
 !write(stderrunit,*) 'diamonds, read_restart_bergs: nbergs in file =', nbergs_in_file

  lonid=inq_var(ncid, 'lon')
  latid=inq_var(ncid, 'lat')
  uvelid=inq_var(ncid, 'uvel')
  vvelid=inq_var(ncid, 'vvel')
  massid=inq_var(ncid, 'mass')
  thicknessid=inq_var(ncid, 'thickness')
  widthid=inq_var(ncid, 'width')
  lengthid=inq_var(ncid, 'length')
  start_lonid=inq_var(ncid, 'start_lon')
  start_latid=inq_var(ncid, 'start_lat')
  start_yearid=inq_var(ncid, 'start_year')
  start_dayid=inq_var(ncid, 'start_day')
  start_massid=inq_var(ncid, 'start_mass')
  scaling_id=inq_var(ncid, 'mass_scaling')
  mass_of_bits_id=inq_var(ncid, 'mass_of_bits',unsafe=.true.)
  heat_density_id=inq_var(ncid, 'heat_density',unsafe=.true.)
  ineid=inq_var(ncid, 'ine',unsafe=.true.)
  jneid=inq_var(ncid, 'jne',unsafe=.true.)

  ! Find approx outer bounds for tile
  lon0=minval( grd%lon(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lon1=maxval( grd%lon(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lat0=minval( grd%lat(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lat1=maxval( grd%lat(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
 !write(stderrunit,'(a,i3,a,4f10.3)') 'diamonds, read_restart_bergs: (',mpp_pe(),') ',lon0,lon1,lat0,lat1
  do k=1, nbergs_in_file
   !write(stderrunit,*) 'diamonds, read_restart_bergs: reading berg ',k
    localberg%lon=get_double(ncid, lonid, k)
    localberg%lat=get_double(ncid, latid, k)
    if (ineid>0 .and. jneid>0 .and. .not. ignore_ij_restart) then ! read i,j position and avoid the "find" step
      localberg%ine=get_int(ncid, ineid, k)
      localberg%jne=get_int(ncid, jneid, k)
      if ( localberg%ine>=grd%isc .and. localberg%ine<=grd%iec .and. &
           localberg%jne>=grd%jsc .and.localberg%jne<=grd%jec ) then
        lres=.true.
      else
        lres=.false.
      endif
    else ! i,j are not available from the file so we search the grid to find out if we reside on this PE
      if (use_slow_find) then
        lres=find_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
      else
        lres=find_cell_by_search(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne)
      endif
    endif
    if (really_debug) then
      write(stderrunit,'(a,i8,a,2f9.4,a,i8)') 'diamonds, read_restart_bergs: berg ',k,' is at ',localberg%lon,localberg%lat,&
           & ' on PE ',mpp_pe()
      write(stderrunit,*) 'diamonds, read_restart_bergs: lres = ',lres
    endif
    if (lres) then ! true if we reside on this PE grid
      localberg%uvel=get_double(ncid, uvelid, k)
      localberg%vvel=get_double(ncid, vvelid, k)
      localberg%mass=get_double(ncid, massid, k)
      localberg%thickness=get_double(ncid, thicknessid, k)
      localberg%width=get_double(ncid, widthid, k)
      localberg%length=get_double(ncid, lengthid, k)
      localberg%start_lon=get_double(ncid, start_lonid, k)
      localberg%start_lat=get_double(ncid, start_latid, k)
      localberg%start_year=get_int(ncid, start_yearid, k)
      localberg%start_day=get_double(ncid, start_dayid, k)
      localberg%start_mass=get_double(ncid, start_massid, k)
      localberg%mass_scaling=get_double(ncid, scaling_id, k)
      if (mass_of_bits_id>0) then ! Allow reading of older restart with no bergy bits
        localberg%mass_of_bits=get_double(ncid, mass_of_bits_id, k)
      else
        localberg%mass_of_bits=0.
      endif
      if (heat_density_id>0) then ! Allow reading of older restart with no heat content
        localberg%heat_density=get_double(ncid, heat_density_id, k)
      else
        localberg%heat_density=0.
      endif
      if (really_debug) lres=is_point_in_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne, explain=.true.)
      lres=pos_within_cell(grd, localberg%lon, localberg%lat, localberg%ine, localberg%jne, localberg%xi, localberg%yj)
     !call add_new_berg_to_list(bergs%first, localberg, quick=.true.)
      call add_new_berg_to_list(bergs%first, localberg)
      if (really_debug) call print_berg(stderrunit, bergs%first, 'read_restart_bergs, add_new_berg_to_list')
    elseif (multiPErestart) then
      call error_mesg('diamonds, read_restart_bergs', 'berg in PE file was not on PE!', FATAL)
    endif
  enddo

  else ! if no restart file was read on this PE
    nbergs_in_file=0
  endif ! if (.not.found_restart)

  ! Sanity check
  k=count_bergs(bergs)
  if (verbose) write(*,'(2(a,i8))') 'diamonds, read_restart_bergs: # bergs =',k,' on PE',mpp_pe()
  if (multiPErestart) call mpp_sum(nbergs_in_file) ! In case PE 0 didn't open a file
  call mpp_sum(k)
  bergs%nbergs_start=k
  if (mpp_pe().eq.mpp_root_pe()) then
    write(*,'(a,i8,a,i8,a)') 'diamonds, read_restart_bergs: there were',nbergs_in_file,' bergs in the restart file and', &
     k,' bergs have been read'
  endif
  if (k.ne.nbergs_in_file) call error_mesg('diamonds, read_restart_bergs', 'wrong number of bergs read!', FATAL)

  if (.not. found_restart .and. bergs%nbergs_start==0 .and. generate_test_icebergs) call generate_bergs(bergs,Time)

  bergs%floating_mass_start=sum_mass(bergs%first)
  call mpp_sum( bergs%floating_mass_start )
  bergs%icebergs_mass_start=sum_mass(bergs%first,justbergs=.true.)
  call mpp_sum( bergs%icebergs_mass_start )
  bergs%bergy_mass_start=sum_mass(bergs%first,justbits=.true.)
  call mpp_sum( bergs%bergy_mass_start )
  if (mpp_pe().eq.mpp_root_pe().and.verbose) write(*,'(a)') 'diamonds, read_restart_bergs: completed'

contains

  subroutine generate_bergs(bergs,Time)
  ! Arguments
  type(icebergs), pointer :: bergs
  type(time_type), intent(in) :: Time
  ! Local variables
  integer :: i,j
  type(iceberg) :: localberg ! NOT a pointer but an actual local variable
  integer :: iyr, imon, iday, ihr, imin, isec

    ! For convenience
    grd=>bergs%grd

    call get_date(Time, iyr, imon, iday, ihr, imin, isec)

    do j=grd%jsc,grd%jec; do i=grd%isc,grd%iec
      if (grd%msk(i,j)>0. .and. abs(grd%latc(i,j))>60.) then
        localberg%xi=0.5
        localberg%yj=0.5
        localberg%ine=i
        localberg%jne=j
        localberg%lon=bilin(grd, grd%lon, i, j, localberg%xi, localberg%yj)
        localberg%lat=bilin(grd, grd%lat, i, j, localberg%xi, localberg%yj)
        localberg%mass=bergs%initial_mass(1)
        localberg%thickness=bergs%initial_thickness(1)
        localberg%width=bergs%initial_width(1)
        localberg%length=bergs%initial_length(1)
        localberg%start_lon=localberg%lon
        localberg%start_lat=localberg%lat
        localberg%start_year=iyr
        localberg%start_day=float(iday)+(float(ihr)+float(imin)/60.)/24.
        localberg%start_mass=localberg%mass
        localberg%mass_scaling=bergs%mass_scaling(1)
        localberg%mass_of_bits=0.
        localberg%heat_density=0.
        localberg%uvel=1.
        localberg%vvel=0.
        call add_new_berg_to_list(bergs%first, localberg)
        localberg%uvel=-1.
        localberg%vvel=0.
        call add_new_berg_to_list(bergs%first, localberg)
        localberg%uvel=0.
        localberg%vvel=1.
        call add_new_berg_to_list(bergs%first, localberg)
        localberg%uvel=0.
        localberg%vvel=-1.
        call add_new_berg_to_list(bergs%first, localberg)
      endif
    enddo; enddo

    bergs%nbergs_start=count_bergs(bergs)
    call mpp_sum(bergs%nbergs_start)
    if (mpp_pe().eq.mpp_root_pe()) &
      write(*,'(a,i8,a)') 'diamonds, generate_bergs: ',bergs%nbergs_start,' were generated'

  end subroutine generate_bergs

end subroutine read_restart_bergs

! ##############################################################################

subroutine read_restart_calving(bergs)
use random_numbers_mod, only: initializeRandomNumberStream, getRandomNumbers, randomNumberStream
! Arguments
type(icebergs), pointer :: bergs
! Local variables
integer :: k,i,j
character(len=37) :: filename, actual_filename
type(icebergs_gridded), pointer :: grd
real, allocatable, dimension(:,:) :: randnum
type(randomNumberStream) :: rns

  ! For convenience
  grd=>bergs%grd

  ! Read stored ice
  filename=trim(restart_input_dir)//'calving.res.nc'
  if (file_exist(filename)) then
    if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(2a)') &
     'diamonds, read_restart_calving: reading ',filename
    call read_data(filename, 'stored_ice', grd%stored_ice, grd%domain)
    if (field_exist(filename, 'stored_heat')) then
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
       'diamonds, read_restart_calving: reading stored_heat from restart file.'
      call read_data(filename, 'stored_heat', grd%stored_heat, grd%domain)
    else
      if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'diamonds, read_restart_calving: stored_heat WAS NOT FOUND in the file. Setting to 0.'
      grd%stored_heat(:,:)=0.
    endif
    bergs%restarted=.true.
  else
    if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(*,'(a)') &
     'diamonds, read_restart_calving: initializing stored ice to random numbers'
    if ( make_calving_reproduce ) then
       allocate(randnum(1,nclasses))
       do j=grd%jsc, grd%jec
          do i=grd%isc, grd%iec
             rns=initializeRandomNumberStream(i+10000*j)
             call getRandomNumbers(rns,randnum(1,:))
             do k=1, nclasses
                grd%stored_ice(i,j,k)=randnum(1,k) * grd%msk(i,j) * bergs%initial_mass(k) * bergs%mass_scaling(k)
             end do
          end do
       end do
    else
       allocate(randnum(grd%jsc:grd%jec,nclasses))
       do i=grd%isc, grd%iec
          rns = initializeRandomNumberStream(i)
          call getRandomNumbers(rns,randnum)
          do k=1, nclasses
             grd%stored_ice(i,grd%jsc:grd%jec,k) = randnum(:,k) * grd%msk(i,grd%jsc:grd%jec) * &
                  & bergs%initial_mass(k) * bergs%mass_scaling(k)
          end do
       end do
    end if
    deallocate(randnum)
  endif

  call grd_chksum3(bergs%grd, bergs%grd%stored_ice, 'read_restart_calving, stored_ice')
  call grd_chksum2(bergs%grd, bergs%grd%stored_heat, 'read_restart_calving, stored_heat')

  bergs%stored_start=sum( grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:) )
  call mpp_sum( bergs%stored_start )
  bergs%stored_heat_start=sum( grd%stored_heat(grd%isc:grd%iec,grd%jsc:grd%jec) )
  call mpp_sum( bergs%stored_heat_start )
  bergs%floating_heat_start=sum_heat(bergs%first)
  call mpp_sum( bergs%floating_heat_start )

end subroutine read_restart_calving

! ##############################################################################

subroutine offset_berg_dates(bergs,Time)
! Arguments
type(icebergs), pointer :: bergs
type(time_type), intent(in) :: Time
! Local variables
type(iceberg), pointer :: this
integer :: iyr, imon, iday, ihr, imin, isec, yr_offset
real :: latest_start_year, berg_start_year

  call get_date(Time, iyr, imon, iday, ihr, imin, isec)
  latest_start_year=iyr-99999

  this=>bergs%first
  if (associated(this)) then
   latest_start_year=float(this%start_year)+this%start_day/367.
  endif
  do while (associated(this))
    berg_start_year=float(this%start_year)+this%start_day/367.
    if (berg_start_year>latest_start_year) latest_start_year=berg_start_year
    this=>this%next
  enddo
  call mpp_max(latest_start_year)

  if (latest_start_year<=float(iyr)+yearday(imon, iday, ihr, imin, isec)/367.) return ! No conflicts!

  yr_offset=int(latest_start_year+1.)-iyr
  if (mpp_pe().eq.mpp_root_pe()) write(*,'(a,i8,a)') &
    'diamonds: Bergs found with creation dates after model date! Adjusting berg dates by ',yr_offset,' years'
  call bergs_chksum(bergs, 'before adjusting start dates')
  this=>bergs%first
  do while (associated(this))
    this%start_year=this%start_year-yr_offset
    this=>this%next
  enddo
  call bergs_chksum(bergs, 'after adjusting start dates')

end subroutine offset_berg_dates

! ##############################################################################

subroutine add_new_berg_to_list(first, bergvals, quick)
! Arguments
type(iceberg), pointer :: first
type(iceberg), intent(in) :: bergvals
logical, intent(in), optional :: quick
! Local variables
type(iceberg), pointer :: new=>null()

  new=>null()
  call create_iceberg(new, bergvals)

  if (present(quick)) then
    if(quick) call insert_berg_into_list(first, new, quick=.true.)
  else
    call insert_berg_into_list(first, new)
  endif

  !Clear new
  new=>null()

end subroutine add_new_berg_to_list

! ##############################################################################

subroutine count_out_of_order(bergs,label)
! Arguments
type(icebergs), pointer :: bergs
character(len=*) :: label
! Local variables
type(iceberg), pointer :: this, next
integer :: i, icnt1, icnt2, icnt3

  icnt1=0; icnt3=0
  this=>bergs%first
  next=>null()
  if (associated(this)) then
    if (associated(this%next)) next=>this%next
  endif
  do while (associated(next))
    if (.not. inorder(this,next)) icnt1=icnt1+1
    if (inorder(this,next).and.inorder(next,this)) icnt3=icnt3+1
    this=>next
    next=>next%next
  enddo
  call mpp_sum(icnt1)

  i=0;
  icnt2=0
  this=>bergs%first
  do while (associated(this))
    i=1
    if (this%ine<bergs%grd%isc .or. &
        this%ine>bergs%grd%iec .or. &
        this%jne<bergs%grd%jsc .or. &
        this%jne>bergs%grd%jec) icnt2=icnt2+1
    this=>this%next
    if (i>1.and..not.associated(this%prev)) then
      call error_mesg('diamonds, count_out_of_order', 'Pointer %prev is unassociated. This should not happen!', FATAL)
    endif
  enddo
  call mpp_sum(icnt2)

  if ((debug.or.icnt1.ne.0).and.mpp_pe().eq.mpp_root_pe()) then
    write(*,'(a,3(x,a,i6),x,a)') 'diamonds, count_out_of_order:', &
      '# out of order=', icnt1,'# in halo=',icnt2,'# identicals=',icnt3,label
  endif

  call check_for_duplicates(bergs,label)

end subroutine count_out_of_order

! ##############################################################################

subroutine check_for_duplicates(bergs,label)
! Arguments
type(icebergs), pointer :: bergs
character(len=*) :: label
! Local variables
type(iceberg), pointer :: this1, next1, this2, next2
integer :: icnt_id, icnt_same

  icnt_id=0
  icnt_same=0
  this1=>bergs%first
  do while (associated(this1))
    this2=>this1%next
    do while (associated(this2))
      if (sameid(this1,this2)) icnt_id=icnt_id+1
      if (sameberg(this1,this2)) icnt_same=icnt_same+1
      this2=>this2%next
    enddo
    this1=>this1%next
  enddo
  call mpp_sum(icnt_id)
  call mpp_sum(icnt_same)

  if ((debug.or.icnt_id>0.or.icnt_same>0).and.mpp_pe().eq.mpp_root_pe()) then
    write(*,'(a,2(x,a,i9),x,a)') 'diamonds, check_for_duplicates:', &
      '# with same id=', icnt_id,'# identical bergs=',icnt_same,label
  endif

end subroutine check_for_duplicates

! ##############################################################################

subroutine insert_berg_into_list(first, newberg, quick)
! Arguments
type(iceberg), pointer :: first, newberg
logical, intent(in), optional :: quick
! Local variables
type(iceberg), pointer :: this, prev
logical :: quickly = .false.

if(present(quick)) quickly = quick

  if (associated(first)) then
    if (.not. parallel_reprod .or. quickly) then
      newberg%next=>first
      first%prev=>newberg
      first=>newberg
    else
      if (inorder(newberg,first)) then
        ! Insert at front of list
        newberg%next=>first
        first%prev=>newberg
        first=>newberg
      else
        this=>first
        prev=>null()
        do while( associated(this) )
          if (inorder(newberg,this) ) then
            exit
          endif
          prev=>this
          this=>this%next
        enddo
        prev%next=>newberg
        newberg%prev=>prev
        if (associated(this)) this%prev=>newberg
        newberg%next=>this
      endif
    endif
  else
    ! list is empty so create it
    first=>newberg
  endif

end subroutine insert_berg_into_list

! ##############################################################################

logical function inorder(berg1, berg2)
! Arguments
type(iceberg), pointer :: berg1, berg2
! Local variables
  if (berg1%start_year<berg2%start_year) then ! want newer first
    inorder=.true.
    return
  else if (berg1%start_year>berg2%start_year) then
    inorder=.false.
    return
  endif
  if (berg1%start_day<berg2%start_day) then ! want newer first
    inorder=.true.
    return
  else if (berg1%start_day>berg2%start_day) then
    inorder=.false.
    return
  endif
  if (berg1%start_mass<berg2%start_mass) then ! want lightest first
    inorder=.true.
    return
  else if (berg1%start_mass>berg2%start_mass) then
    inorder=.false.
    return
  endif
  if (berg1%start_lon<berg2%start_lon) then ! want eastward first
    inorder=.true.
    return
  else if (berg1%start_lon>berg2%start_lon) then
    inorder=.false.
    return
  endif
  if (berg1%start_lat<berg2%start_lat) then ! want southern first
    inorder=.true.
    return
  else if (berg1%start_lat>berg2%start_lat) then
    inorder=.false.
    return
  endif
  inorder=.true. ! passing the above tests mean the bergs 1 and 2 are identical?
end function inorder

! ##############################################################################

  real function time_hash(berg)
  ! Arguments
  type(iceberg), pointer :: berg
    time_hash=berg%start_day+366.*float(berg%start_year)
  end function time_hash

! ##############################################################################

  real function pos_hash(berg)
  ! Arguments
  type(iceberg), pointer :: berg
    pos_hash=berg%start_lon+360.*(berg%start_lat+90.)
  end function pos_hash

! ##############################################################################

logical function sameid(berg1, berg2)
! Arguments
type(iceberg), pointer :: berg1, berg2
! Local variables
  sameid=.false.
  if (berg1%start_year.ne.berg2%start_year) return
  if (berg1%start_day.ne.berg2%start_day) return
  if (berg1%start_mass.ne.berg2%start_mass) return
  if (berg1%start_lon.ne.berg2%start_lon) return
  if (berg1%start_lat.ne.berg2%start_lat) return
  sameid=.true. ! passing the above tests means that bergs 1 and 2 have the same id
end function sameid

! ##############################################################################

logical function sameberg(berg1, berg2)
! Arguments
type(iceberg), pointer :: berg1, berg2
! Local variables
  sameberg=.false.
  if (.not. sameid(berg1, berg2)) return
  if (berg1%lon.ne.berg2%lon) return
  if (berg1%lat.ne.berg2%lat) return
  if (berg1%mass.ne.berg2%mass) return
  if (berg1%uvel.ne.berg2%uvel) return
  if (berg1%vvel.ne.berg2%vvel) return
  if (berg1%thickness.ne.berg2%thickness) return
  if (berg1%width.ne.berg2%width) return
  if (berg1%length.ne.berg2%length) return
  sameberg=.true. ! passing the above tests mean that bergs 1 and 2 are identical
end function sameberg

! ##############################################################################

real function yearday(imon, iday, ihr, imin, isec)
! Arguments
integer, intent(in) :: imon, iday, ihr, imin, isec

  yearday=float(imon-1)*31.+float(iday-1)+(float(ihr)+(float(imin)+float(isec)/60.)/60.)/24.

end function yearday

! ##############################################################################

subroutine create_iceberg(berg, bergvals)
! Arguments
type(iceberg), pointer :: berg
type(iceberg), intent(in) :: bergvals
! Local variables
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  if (associated(berg)) then
    write(stderrunit,*) 'diamonds, create_iceberg: berg already associated!!!!',mpp_pe()
    call error_mesg('diamonds, create_iceberg', 'berg already associated. This should not happen!', FATAL)
  endif
  allocate(berg)
  berg=bergvals
  berg%prev=>null()
  berg%next=>null()

end subroutine create_iceberg

! ##############################################################################

subroutine delete_iceberg_from_list(first, berg)
! Arguments
type(iceberg), pointer :: first, berg
! Local variables

  ! Connect neighbors to each other
  if (associated(berg%prev)) then
    berg%prev%next=>berg%next
  else
    first=>berg%next
  endif
  if (associated(berg%next)) berg%next%prev=>berg%prev

  ! Bye-bye berg
  call destroy_iceberg(berg)

end subroutine delete_iceberg_from_list

! ##############################################################################

subroutine destroy_iceberg(berg)
! Arguments
type(iceberg), pointer :: berg
! Local variables

  ! Bye-bye berg
  deallocate(berg)

end subroutine destroy_iceberg

! ##############################################################################

subroutine print_berg(iochan, berg, label)
! Arguments
integer, intent(in) :: iochan
type(iceberg), pointer :: berg
character(len=*) :: label
! Local variables

  write(iochan,'("diamonds, print_berg: ",a," pe=(",i3,") start lon,lat,yr,day,mass=",2f10.4,i5,f7.2,es12.4)') &
    label, mpp_pe(), berg%start_lon, berg%start_lat, &
    berg%start_year, berg%start_day, berg%start_mass
  write(iochan,'("diamonds, print_berg: ",a," pe=(",i3,a,2i5,3(a,2f10.4),a,2l2)') &
    label, mpp_pe(), ') i,j=',berg%ine, berg%jne, &
    ' xi,yj=', berg%xi, berg%yj, &
    ' lon,lat=', berg%lon, berg%lat, &
    ' u,v=', berg%uvel, berg%vvel, &
    ' p,n=', associated(berg%prev), associated(berg%next)
  write(iochan,'("diamonds, print_berg: ",a," pe=(",i3,") ",6(a,2f10.4))') &
    label, mpp_pe(), 'uo,vo=', berg%uo, berg%vo, 'ua,va=', berg%ua, berg%va, 'ui,vi=', berg%ui, berg%vi

end subroutine print_berg

! ##############################################################################

subroutine print_bergs(iochan, bergs, label)
! Arguments
integer, intent(in) :: iochan
type(icebergs), pointer :: bergs
character(len=*) :: label
! Local variables
integer :: nbergs, nnbergs
type(iceberg), pointer :: this

  this=>bergs%first
  do while(associated(this))
    call print_berg(iochan, this, label)
    this=>this%next
  enddo
  nbergs=count_bergs(bergs)
  nnbergs=nbergs
  call mpp_sum(nnbergs)
  if (nbergs.gt.0) write(iochan,'("diamonds, ",a," there are",i5," bergs out of",i6," on PE ",i4)') label, nbergs, nnbergs, mpp_pe()

end subroutine print_bergs

! ##############################################################################

integer function count_bergs(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(iceberg), pointer :: this

  count_bergs=0
  this=>bergs%first
  do while(associated(this))
    count_bergs=count_bergs+1
    this=>this%next
  enddo

end function count_bergs

! ##############################################################################

subroutine record_posn(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(xyt) :: posn
type(iceberg), pointer :: this

  this=>bergs%first
  do while (associated(this))
    posn%lon=this%lon
    posn%lat=this%lat
    posn%year=bergs%current_year
    posn%day=bergs%current_yearday
    posn%uvel=this%uvel
    posn%vvel=this%vvel
    posn%mass=this%mass
    posn%mass_of_bits=this%mass_of_bits
    posn%heat_density=this%heat_density
    posn%thickness=this%thickness
    posn%width=this%width
    posn%length=this%length
    posn%uo=this%uo
    posn%vo=this%vo
    posn%ui=this%ui
    posn%vi=this%vi
    posn%ua=this%ua
    posn%va=this%va
    posn%ssh_x=this%ssh_x
    posn%ssh_y=this%ssh_y
    posn%sst=this%sst
    posn%cn=this%cn
    posn%hi=this%hi

    call push_posn(this%trajectory, posn)

    this=>this%next
  enddo

end subroutine record_posn

! ##############################################################################

subroutine push_posn(trajectory, posn_vals)
! Arguments
type(xyt), pointer :: trajectory
type(xyt) :: posn_vals
! Local variables
type(xyt), pointer :: new_posn

  allocate(new_posn)
  new_posn=posn_vals
  new_posn%next=>trajectory
  trajectory=>new_posn

end subroutine push_posn

! ##############################################################################

subroutine move_trajectory(bergs, berg)
! Arguments
type(icebergs), pointer :: bergs
type(iceberg), pointer :: berg
! Local variables
type(xyt), pointer :: next, last
type(xyt) :: vals

  ! If the trajectory is empty, ignore it
  if (.not.associated(berg%trajectory)) return

  ! Push identifying info into first posn (note reverse order due to stack)
  vals%lon=berg%start_lon
  vals%lat=berg%start_lat
  vals%year=berg%start_year
  vals%day=berg%start_day
  vals%mass=berg%start_mass
  call push_posn(berg%trajectory, vals)
  vals%lon=0.
  vals%lat=99.
  vals%year=0
  vals%day=0.
  vals%mass=berg%mass_scaling
  call push_posn(berg%trajectory, vals)

  ! Find end of berg trajectory and point it to start of existing trajectories
  next=>berg%trajectory
  do while (associated(next))
    last=>next
    next=>next%next
  enddo
  last%next=>bergs%trajectories

  bergs%trajectories=>berg%trajectory
  berg%trajectory=>null()

end subroutine move_trajectory

! ##############################################################################

logical function find_cell_by_search(grd, x, y, i, j)
! Arguments
type(icebergs_gridded), pointer :: grd
real, intent(in) :: x, y
integer, intent(inout) :: i, j
! Local variables
integer :: is,ie,js,je,di,dj,io,jo,icnt
real :: d0,d1,d2,d3,d4,d5,d6,d7,d8,dmin
logical :: explain=.false.

911 continue

  find_cell_by_search=.false.
  is=grd%isc; ie=grd%iec; js=grd%jsc; je=grd%jec

  ! Start at nearest corner
  d1=dcost(x,y,grd%lonc(is+1,js+1),grd%latc(is+1,js+1))
  d2=dcost(x,y,grd%lonc(ie-1,js+1),grd%latc(ie-1,js+1))
  d3=dcost(x,y,grd%lonc(ie-1,je-1),grd%latc(ie-1,je-1))
  d4=dcost(x,y,grd%lonc(is+1,je-1),grd%latc(is+1,je-1))
  dmin=min(d1,d2,d3,d4)
  if (d1==dmin) then; i=is+1; j=js+1
  elseif (d2==dmin) then; i=ie-1; j=js+1
  elseif (d3==dmin) then; i=ie-1; j=je-1
  elseif (d4==dmin) then; i=is+1; j=je-1
  else
    call error_mesg('diamonds, find_cell_by_search:', 'This should never EVER happen! (1)', FATAL)
  endif

  if (explain) then
    write(0,'(i3,a,2i4,f9.4)') mpp_pe(),'Initial corner i-is,j-js,cost=',i-is,j-js,dmin
    write(0,'(i3,a,3f9.4)') mpp_pe(),'cost ',d4,d3
    write(0,'(i3,a,3f9.4)') mpp_pe(),'cost ',d1,d2
  endif

  if (is_point_in_cell(grd, x, y, i, j)) then
    find_cell_by_search=.true.
    return
  endif

  do icnt=1, 1*(ie-is+je-js)
    io=i; jo=j

    d0=dcost(x,y,grd%lonc(io,jo),grd%latc(io,jo))
    d1=dcost(x,y,grd%lonc(io,jo+1),grd%latc(io,jo+1))
    d2=dcost(x,y,grd%lonc(io-1,jo+1),grd%latc(io-1,jo+1))
    d3=dcost(x,y,grd%lonc(io-1,jo),grd%latc(io-1,jo))
    d4=dcost(x,y,grd%lonc(io-1,jo-1),grd%latc(io-1,jo-1))
    d5=dcost(x,y,grd%lonc(io,jo-1),grd%latc(io,jo-1))
    d6=dcost(x,y,grd%lonc(io+1,jo-1),grd%latc(io+1,jo-1))
    d7=dcost(x,y,grd%lonc(io+1,jo),grd%latc(io+1,jo))
    d8=dcost(x,y,grd%lonc(io+1,jo+1),grd%latc(io+1,jo+1))

  ! dmin=min(d0,d1,d3,d5,d7)
    dmin=min(d0,d1,d2,d3,d4,d5,d6,d7,d8)
    if (d0==dmin) then; di=0; dj=0
    elseif (d2==dmin) then; di=-1; dj=1
    elseif (d4==dmin) then; di=-1; dj=-1
    elseif (d6==dmin) then; di=1; dj=-1
    elseif (d8==dmin) then; di=1; dj=1
    elseif (d1==dmin) then; di=0; dj=1
    elseif (d3==dmin) then; di=-1; dj=0
    elseif (d5==dmin) then; di=0; dj=-1
    elseif (d7==dmin) then; di=1; dj=0
    else
      call error_mesg('diamonds, find_cell_by_search:', 'This should never EVER happen!', FATAL)
    endif

    i=min(ie, max(is, io+di))
    j=min(je, max(js, jo+dj))

    if (explain) then
      write(0,'(i3,a,2i4,f9.5,a,2i4,a,2i4)') mpp_pe(),'Current position i,j,cost=',i,j,dmin,' di,dj=',di,dj,' old io,jo=',io,jo
      write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d2,d1,d8
      write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d3,d0,d7
      write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d4,d5,d6
    endif

    if (is_point_in_cell(grd, x, y, i, j)) then
      find_cell_by_search=.true.
      return
    endif

    if ((i==io.and.j==jo) &
        .and. .not.find_better_min(grd, x, y, 3, i, j) &
       ) then
      ! Stagnated
      find_cell_by_search=find_cell_loc(grd, x, y, is, ie, js, je, 1, i, j)
      if (.not. find_cell_by_search) find_cell_by_search=find_cell_loc(grd, x, y, is, ie, js, je, 3, i, j)
      i=min(ie, max(is, i))
      j=min(je, max(js, j))
      if (is_point_in_cell(grd, x, y, i, j)) then
        find_cell_by_search=.true.
      else
  !     find_cell_by_search=find_cell(grd, x, y, io, jo)
  !     if (find_cell_by_search) then
  !       if (explain) then
  !         call print_fld(grd, grd%lat, 'Lat')
  !         call print_fld(grd, grd%lon, 'Lat')
  !         do j=grd%jsd, grd%jed; do i=grd%isd, grd%ied
  !           grd%tmp(i,j) = dcost(x,y,grd%lonc(i,j),grd%latc(i,j))
  !         enddo; enddo
  !         call print_fld(grd, grd%tmp, 'Cost')
  !         stop 'Avoid recursing'
  !       endif
  !       write(0,'(i3,a,2i5,a,2i3,a,2f8.3)') mpp_pe(),'diamonds, find_cell_by_search: false negative io,jo=',io,jo,' di,dj=',di,dj,' targ=',x,y
  !       explain=.true.; goto 911
  !     endif
      endif
      return
    endif

  enddo

  find_cell_by_search=find_cell(grd, x, y, i, j)
  if (find_cell_by_search) then
    write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d2,d1,d8
    write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d3,d0,d7
    write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d4,d5,d6
    write(0,'(i3,a,2f9.5)') mpp_pe(),'x,y ',x,y
    write(0,'(i3,a,4i5)') mpp_pe(),'io,jo ',io,jo,di,dj
    write(0,'(i3,a,2i5,a,2i3)') mpp_pe(),'diamonds, find_cell_by_search: false negative 2 i,j=',i-is,j-js,' di,dj=',di,dj
    write(0,'(i3,a,2i5,a,2f8.3)') mpp_pe(),'diamonds, find_cell_by_search: false negative 2 io,jo=',io,jo
    write(0,'(i3,a,2i5,a,2f8.3)') mpp_pe(),'diamonds, find_cell_by_search: false negative 2 i,j=',i,j,' targ=',x,y
    return
  endif
  find_cell_by_search=.false.

  contains

! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  real function dcost(x1, y1, x2, y2)
  ! Arguments
  real, intent(in) :: x1, x2, y1, y2
  ! Local variables
  real :: x1m

    x1m=modulo(x1-(x2-180.),360.)+(x2-180.)
  ! dcost=(x2-x1)**2+(y2-y1)**2
    dcost=(x2-x1m)**2+(y2-y1)**2
  end function dcost

! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  logical function find_better_min(grd, x, y, w, oi, oj)
  ! Arguments
  type(icebergs_gridded), intent(in) :: grd
  real, intent(in) :: x, y
  integer, intent(in) :: w
  integer, intent(inout) :: oi, oj
  ! Local variables
  integer :: i,j,xs,xe,ys,ye
  real :: dmin, dcst

  xs=max(grd%isc, oi-w)
  xe=min(grd%iec, oi+w)
  ys=max(grd%jsc, oj-w)
  ye=min(grd%jec, oj+w)

  find_better_min=.false.
  dmin=dcost(x,y,grd%lonc(oi,oj),grd%latc(oi,oj))
  do j=ys,ye; do i=xs,xe
      dcst=dcost(x,y,grd%lonc(i,j),grd%latc(i,j))
      if (dcst<dmin) then
        find_better_min=.true.
        dmin=dcst
        oi=i;oj=j
      endif
  enddo; enddo

  end function find_better_min

! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  logical function find_cell_loc(grd, x, y, is, ie, js, je, w, oi, oj)
  ! Arguments
  type(icebergs_gridded), intent(in) :: grd
  real, intent(in) :: x, y
  integer, intent(in) :: is, ie, js, je, w
  integer, intent(inout) :: oi, oj
  ! Local variables
  integer :: i,j,xs,xe,ys,ye

    xs=max(is, oi-w)
    xe=min(ie, oi+w)
    ys=max(js, oj-w)
    ye=min(je, oj+w)

    find_cell_loc=.false.
    do j=ys,ye; do i=xs,xe
        if (is_point_in_cell(grd, x, y, i, j)) then
          oi=i; oj=j; find_cell_loc=.true.
          return
        endif
    enddo; enddo

  end function find_cell_loc

end function find_cell_by_search

! ##############################################################################

logical function find_cell(grd, x, y, oi, oj)
! Arguments
type(icebergs_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(out) :: oi, oj
! Local variables
integer :: i,j

  find_cell=.false.; oi=-999; oj=-999

  do j=grd%jsc,grd%jec; do i=grd%isc,grd%iec
      if (is_point_in_cell(grd, x, y, i, j)) then
        oi=i; oj=j; find_cell=.true.
        return
      endif
  enddo; enddo

end function find_cell

! ##############################################################################

logical function find_cell_wide(grd, x, y, oi, oj)
! Arguments
type(icebergs_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(out) :: oi, oj
! Local variables
integer :: i,j

  find_cell_wide=.false.; oi=-999; oj=-999

  do j=grd%jsd+1,grd%jed; do i=grd%isd+1,grd%ied
      if (is_point_in_cell(grd, x, y, i, j)) then
        oi=i; oj=j; find_cell_wide=.true.
        return
      endif
  enddo; enddo

end function find_cell_wide

! ##############################################################################

logical function is_point_in_cell(grd, x, y, i, j, explain)
! Arguments
type(icebergs_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(in) :: i, j
logical, intent(in), optional :: explain
! Local variables
real :: xlo, xhi, ylo, yhi
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  ! Safety check index bounds
  if (i-1.lt.grd%isd.or.i.gt.grd%ied.or.j-1.lt.grd%jsd.or.j.gt.grd%jed) then
    write(stderrunit,'(a,i3,(a,3i4))') &
                     'diamonds, is_point_in_cell: pe=(',mpp_pe(),') i,s,e=', &
                     i,grd%isd,grd%ied,' j,s,e=', j,grd%jsd,grd%jed
    call error_mesg('diamonds, is_point_in_cell', 'test is off the PE!', FATAL)
  endif

  is_point_in_cell=.false.

  ! Test crude bounds
  xlo=min( modulo(grd%lon(i-1,j-1)-(x-180.),360.)+(x-180.), &
           modulo(grd%lon(i  ,j-1)-(x-180.),360.)+(x-180.), &
           modulo(grd%lon(i-1,j  )-(x-180.),360.)+(x-180.), &
           modulo(grd%lon(i  ,j  )-(x-180.),360.)+(x-180.) )
  xhi=max( modulo(grd%lon(i-1,j-1)-(x-180.),360.)+(x-180.), &
           modulo(grd%lon(i  ,j-1)-(x-180.),360.)+(x-180.), &
           modulo(grd%lon(i-1,j  )-(x-180.),360.)+(x-180.), &
           modulo(grd%lon(i  ,j  )-(x-180.),360.)+(x-180.) )
  if (x.lt.xlo .or. x.gt.xhi) return
  ylo=min( grd%lat(i-1,j-1), grd%lat(i,j-1), grd%lat(i-1,j), grd%lat(i,j) )
  yhi=max( grd%lat(i-1,j-1), grd%lat(i,j-1), grd%lat(i-1,j), grd%lat(i,j) )
  if (y.lt.ylo .or. y.gt.yhi) return

  if (grd%lat(i,j).gt.89.999) then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                        x, y, explain=explain)
  elseif (grd%lat(i-1,j).gt.89.999) then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i  ,j  ),grd%lat(i-1,j  ), &
                                        grd%lon(i-1,j-1),grd%lat(i-1,j  ), &
                                        x, y, explain=explain)
  elseif (grd%lat(i-1,j-1).gt.89.999) then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j  ),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                        x, y, explain=explain)
  elseif (grd%lat(i,j-1).gt.89.999) then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i-1,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                        x, y, explain=explain)
  else
  is_point_in_cell=sum_sign_dot_prod4(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                      grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                      grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                      grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                      x, y, explain=explain)
  endif

end function is_point_in_cell

! ##############################################################################

logical function sum_sign_dot_prod4(x0, y0, x1, y1, x2, y2, x3, y3, x, y, explain)
! Arguments
real, intent(in) :: x0, y0, x1, y1, x2, y2, x3, y3, x, y
logical, intent(in), optional :: explain
! Local variables
real :: p0,p1,p2,p3,xx
real :: l0,l1,l2,l3
real :: xx0,xx1,xx2,xx3
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  sum_sign_dot_prod4=.false.
  xx=modulo(x-(x0-180.),360.)+(x0-180.) ! Reference x to within 180 of x0
  xx0=modulo(x0-(x0-180.),360.)+(x0-180.) ! Reference x0 to within 180 of xx
  xx1=modulo(x1-(x0-180.),360.)+(x0-180.) ! Reference x1 to within 180 of xx
  xx2=modulo(x2-(x0-180.),360.)+(x0-180.) ! Reference x2 to within 180 of xx
  xx3=modulo(x3-(x0-180.),360.)+(x0-180.) ! Reference x3 to within 180 of xx

  l0=(xx-xx0)*(y1-y0)-(y-y0)*(xx1-xx0)
  l1=(xx-xx1)*(y2-y1)-(y-y1)*(xx2-xx1)
  l2=(xx-xx2)*(y3-y2)-(y-y2)*(xx3-xx2)
  l3=(xx-xx3)*(y0-y3)-(y-y3)*(xx0-xx3)

  p0=sign(1., l0); if (l0.eq.0.) p0=0.
  p1=sign(1., l1); if (l1.eq.0.) p1=0.
  p2=sign(1., l2); if (l2.eq.0.) p2=0.
  p3=sign(1., l3); if (l3.eq.0.) p3=0.

  if ( (abs(p0)+abs(p2))+(abs(p1)+abs(p3)) .eq. abs((p0+p2)+(p1+p3)) ) then
    sum_sign_dot_prod4=.true.
  endif


  if (present(explain)) then
   if(explain) then
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: x=',mpp_pe(),':', &
                           x0,x1,x2,x3, x
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: X=',mpp_pe(),':', &
                           xx0,xx1,xx2,xx3, xx
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: y=',mpp_pe(),':', &
                           y0,y1,y2,y3, y
   write(stderrunit,'(a,i3,a,1p10e12.4)') 'sum_sign_dot_prod4: l=',mpp_pe(),':', &
                           l0,l1,l2,l3
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: p=',mpp_pe(),':', &
                           p0,p1,p2,p3, abs( (abs(p0)+abs(p2))+(abs(p1)+abs(p3)) - abs((p0+p2)+(p1+p3)) )
   endif
  endif

end function sum_sign_dot_prod4

! ##############################################################################

logical function sum_sign_dot_prod5(x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x, y, explain)
! Arguments
real, intent(in) :: x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x, y
logical, intent(in), optional :: explain
! Local variables
real :: p0,p1,p2,p3,p4,xx
real :: l0,l1,l2,l3,l4
real :: xx0,xx1,xx2,xx3,xx4
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  sum_sign_dot_prod5=.false.
  xx=modulo(x-(x0-180.),360.)+(x0-180.) ! Reference x to within 180 of x0
  xx0=modulo(x0-(x0-180.),360.)+(x0-180.) ! Reference x0 to within 180 of xx
  xx1=modulo(x1-(x0-180.),360.)+(x0-180.) ! Reference x1 to within 180 of xx
  xx2=modulo(x2-(x0-180.),360.)+(x0-180.) ! Reference x2 to within 180 of xx
  xx3=modulo(x3-(x0-180.),360.)+(x0-180.) ! Reference x3 to within 180 of xx
  xx4=modulo(x4-(x0-180.),360.)+(x0-180.) ! Reference x4 to within 180 of xx

  l0=(xx-xx0)*(y1-y0)-(y-y0)*(xx1-xx0)
  l1=(xx-xx1)*(y2-y1)-(y-y1)*(xx2-xx1)
  l2=(xx-xx2)*(y3-y2)-(y-y2)*(xx3-xx2)
  l3=(xx-xx3)*(y4-y3)-(y-y3)*(xx4-xx3)
  l4=(xx-xx4)*(y0-y4)-(y-y4)*(xx0-xx4)

  p0=sign(1., l0); if (l0.eq.0.) p0=0.
  p1=sign(1., l1); if (l1.eq.0.) p1=0.
  p2=sign(1., l2); if (l2.eq.0.) p2=0.
  p3=sign(1., l3); if (l3.eq.0.) p3=0.
  p4=sign(1., l4); if (l4.eq.0.) p4=0.

  if ( ((abs(p0)+abs(p2))+(abs(p1)+abs(p3)))+abs(p4) - abs(((p0+p2)+(p1+p3))+p4) .lt. 0.5 ) then
    sum_sign_dot_prod5=.true.
  endif

  if (present(explain)) then
   if(explain) then
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: x=',mpp_pe(),':', &
                           x0,x1,x2,x3,x4, x
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: X=',mpp_pe(),':', &
                           xx0,xx1,xx2,xx3,xx4, xx
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: y=',mpp_pe(),':', &
                           y0,y1,y2,y3,y4, y
   write(stderrunit,'(a,i3,a,1p10e12.4)') 'sum_sign_dot_prod5: l=',mpp_pe(),':', &
                           l0,l1,l2,l3,l4
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: p=',mpp_pe(),':', &
                           p0,p1,p2,p3,p4
   endif
  endif

end function sum_sign_dot_prod5

! ##############################################################################

logical function pos_within_cell(grd, x, y, i, j, xi, yj, explain)
! Arguments
type(icebergs_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(in) :: i, j
real, intent(out) :: xi, yj
logical, intent(in), optional :: explain
! Local variables
real :: x1,y1,x2,y2,x3,y3,x4,y4,xx,yy,fac
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  pos_within_cell=.false.; xi=-999.; yj=-999.
  if (i-1<grd%isd) return
  if (j-1<grd%jsd) return
  if (i>grd%ied) return
  if (j>grd%jed) return

  x1=grd%lon(i-1,j-1)
  y1=grd%lat(i-1,j-1)
  x2=grd%lon(i  ,j-1)
  y2=grd%lat(i  ,j-1)
  x3=grd%lon(i  ,j  )
  y3=grd%lat(i  ,j  )
  x4=grd%lon(i-1,j  )
  y4=grd%lat(i-1,j  )

  if (present(explain)) then
    if(explain) then
    write(stderrunit,'(a,4f12.6)') 'pos_within_cell: x1..x4 ',x1,x2,x3,x4
    write(stderrunit,'(a,2f12.6)') 'pos_within_cell: x',x
    write(stderrunit,'(a,4f12.6)') 'pos_within_cell: y1..y4 ',y1,y2,y3,y4
    write(stderrunit,'(a,2f12.6)') 'pos_within_cell: y',y
    endif
  endif

  if (max(y1,y2,y3,y4)<89.999) then
    call calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, x, y, xi, yj, explain=explain)
  else
    if (debug) write(stderrunit,*) 'diamonds, pos_within_cell: working in tangential plane!'
    xx=(90.-y)*cos(x*pi_180)
    yy=(90.-y)*sin(x*pi_180)
    x1=(90.-y1)*cos(grd%lon(i-1,j-1)*pi_180)
    y1=(90.-y1)*sin(grd%lon(i-1,j-1)*pi_180)
    x2=(90.-y2)*cos(grd%lon(i  ,j-1)*pi_180)
    y2=(90.-y2)*sin(grd%lon(i  ,j-1)*pi_180)
    x3=(90.-y3)*cos(grd%lon(i  ,j  )*pi_180)
    y3=(90.-y3)*sin(grd%lon(i  ,j  )*pi_180)
    x4=(90.-y4)*cos(grd%lon(i-1,j  )*pi_180)
    y4=(90.-y4)*sin(grd%lon(i-1,j  )*pi_180)
    if (present(explain)) then
      if(explain) then
      write(stderrunit,'(a,4f12.6)') 'pos_within_cell: x1..x4 ',x1,x2,x3,x4
      write(stderrunit,'(a,2f12.6)') 'pos_within_cell: x',xx
      write(stderrunit,'(a,4f12.6)') 'pos_within_cell: y1..y4 ',y1,y2,y3,y4
      write(stderrunit,'(a,2f12.6)') 'pos_within_cell: y',yy
      endif
    endif
    call calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, xx, yy, xi, yj, explain=explain)
    if (is_point_in_cell(grd, x, y, i, j)) then
      if (abs(xi-0.5)>0.5.or.abs(yj-0.5)>0.5) then
        ! Scale internal coordinates to be consistent with is_point_in_cell()
        ! Note: this is intended to fix the inconsistency between the tangent plane
        ! and lat-lon calculations
        fac=2.*max( abs(xi-0.5), abs(yj-0.5) ); fac=max(1., fac)
        xi=0.5+(xi-0.5)/fac
        yj=0.5+(yj-0.5)/fac
        if (debug) call error_mesg('diamonds, pos_within_cell', 'in cell so scaling internal coordinates!', WARNING)
      endif
    else
      if (abs(xi-0.5)<=0.5.and.abs(yj-0.5)<=0.5) then
        if (debug) call error_mesg('diamonds, pos_within_cell', 'out of cell but coordinates <=0.5!', WARNING)
      endif
    endif
  endif

  if (present(explain)) then
	if(explain) write(stderrunit,'(a,2f12.6)') 'pos_within_cell: xi,yj=',xi,yj
  endif

 !if (.not. is_point_in_cell(grd, x, y, i, j) ) then
 !   write(stderrunit,'(a,i3,a,8f8.2,a)') 'diamonds, pos_within_cell: (',mpp_pe(),') ', &
 !                   x1, y1, x2, y2, x3, y3, x4, y4, ' NOT IN CELL!'
 !endif

  if (xi.ge.0. .and. xi.le.1. .and. yj.ge.0. .and. yj.le.1.) then
    pos_within_cell=is_point_in_cell(grd, x, y, i, j, explain=explain)
    if (.not. pos_within_cell .and. verbose) then
      if (debug) call error_mesg('diamonds, pos_within_cell', 'pos_within_cell is in cell BUT is_point_in_cell disagrees!', WARNING)
    endif
   !pos_within_cell=.true. ! commenting this out makes pos_within_cell agree with is_point_in_cell
  endif

  contains

  subroutine calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, x, y, xi, yj, explain)
  ! Arguments
  real,  intent(in) :: x1, x2, x3, x4, y1, y2, y3, y4, x, y
  real, intent(out) :: xi, yj
  logical, intent(in), optional :: explain
  ! Local variables
  real :: alpha, beta, gamma, delta, epsilon, kappa, a, b, c, d, dx, dy, yy1, yy2
  logical :: expl=.false.
  integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  expl=.false.
  if (present(explain)) then
     if(explain) expl=.true.
  endif
  alpha=x2-x1
  delta=y2-y1
  beta=x4-x1
  epsilon=y4-y1
  gamma=(x3-x1)-(alpha+beta)
  kappa=(y3-y1)-(delta+epsilon)
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs alpha,beta,gamma',alpha,beta,gamma,delta,epsilon,kappa
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs delta,epsilon,kappa',alpha,beta,gamma,delta,epsilon,kappa

  a=(kappa*beta-gamma*epsilon)
  dx=modulo(x-(x1-180.),360.)+(x1-180.)-x1
  dy=y-y1
  b=(delta*beta-alpha*epsilon)-(kappa*dx-gamma*dy)
  c=(alpha*dy-delta*dx)
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs dx,dy=',dx,dy
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs A,B,C=',a,b,c
  if (abs(a)>1.e-12) then
    d=0.25*(b**2)-a*c
    if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs D=',d
    if (d.ge.0.) then
      if (expl) write(stderrunit,'(a,1p3e12.4)') 'Roots for b/2a, sqrt(d) = ',-0.5*b/a,sqrt(d)/a
      yy1=-(0.5*b+sqrt(d))/a
      yy2=-(0.5*b-sqrt(d))/a
      if (abs(yy1-0.5).lt.abs(yy2-0.5)) then; yj=yy1; else; yj=yy2; endif
      if (expl) write(stderrunit,'(a,1p3e12.4)') 'Roots for y = ',yy1,yy2,yj
    else
      write(stderrunit,'(a,i3,4f8.2)') 'calc_xiyj: x1..x4 ',mpp_pe(),x1,x2,x3,x4
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: x2..x4 - x1',mpp_pe(),x2-x1,x3-x1,x4-x1
      write(stderrunit,'(a,i3,4f8.2)') 'calc_xiyj: y1..y4 ',mpp_pe(),y1,y2,y3,y4
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: y2..y4 - x1',mpp_pe(),y2-y1,y3-y1,y4-y1
      write(stderrunit,'(a,i3,1p6e12.4)') 'calc_xiyj: coeffs alpha..kappa',mpp_pe(),alpha,beta,gamma,delta,epsilon,kappa
      write(stderrunit,'(a,i3)') 'calc_xiyj: b<0 in quadratic root solver!!!!',mpp_pe()
      write(stderrunit,'(a,i3,1p6e12.4)') 'calc_xiyj: coeffs a,b,c,d,dx,dy',mpp_pe(),a,b,c,d,dx,dy
      call error_mesg('diamonds, calc_xiyj', 'We have complex roots. The grid must be very distorted!', FATAL)
    endif
  else
    if (b.ne.0.) then
      yj=-c/b
    else
      yj=0.
    endif
  endif

  a=(alpha+gamma*yj)
  b=(delta+kappa*yj)
  if (a.ne.0.) then
    xi=(dx-beta*yj)/a
  elseif (b.ne.0.) then
    xi=(dy-epsilon*yj)/b
  else
    c=(epsilon*alpha-beta*delta)+(epsilon*gamma-beta*kappa)*yj
    if (c.ne.0.) then
      xi=(epsilon*dx-beta*dy)/c
    else
      write(stderrunit,'(a,i3,4f8.2)') 'calc_xiyj: x1..x4 ',mpp_pe(),x1,x2,x3,x4
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: x2..x4 - x1',mpp_pe(),x2-x1,x3-x1,x4-x1
      write(stderrunit,'(a,i3,4f8.2)') 'calc_xiyj: y1..y4 ',mpp_pe(),y1,y2,y3,y4
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: y2..y4 - x1',mpp_pe(),y2-y1,y3-y1,y4-y1
      write(stderrunit,'(a,i3,1p6e12.4)') 'calc_xiyj: coeffs alpha..kappa',mpp_pe(),alpha,beta,gamma,delta,epsilon,kappa
      write(stderrunit,'(a,i3,1p2e12.4)') 'calc_xiyj: coeffs a,b',mpp_pe(),a,b
      call error_mesg('diamonds, calc_xiyj', 'Can not invert either linear equaton for xi! This should not happen!', FATAL)
    endif
  endif
  if (expl) write(stderrunit,'(a,2e12.4)') 'calc_xiyj: xi,yj=',xi,yj

  end subroutine calc_xiyj

end function pos_within_cell

! ##############################################################################

subroutine check_position(grd, berg, label)
! Arguments
type(icebergs_gridded), pointer :: grd
type(iceberg), pointer :: berg
character(len=*) :: label
! Local variables
real :: xi, yj
logical :: lret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  lret=pos_within_cell(grd, berg%lon, berg%lat, berg%ine, berg%jne, xi, yj)
  if (xi.ne.berg%xi.or.yj.ne.berg%yj) then
    write(stderrunit,'("diamonds: check_position (",i4,") b%x,x,-=",3(es12.4,x),a)') mpp_pe(),berg%xi,xi,berg%xi-xi,label
    write(stderrunit,'("diamonds: check_position (",i4,") b%y,y,-=",3(es12.4,x),a)') mpp_pe(),berg%yj,yj,berg%yj-yj,label
    call print_berg(stderrunit, berg, 'check_position')
    call error_mesg('diamonds, check_position','berg has inconsistent xi,yj!',FATAL)
  endif

end subroutine check_position

! ##############################################################################

real function sum_mass(first,justbits,justbergs)
! Arguments
type(iceberg), pointer :: first
logical, intent(in), optional :: justbits, justbergs
! Local variables
type(iceberg), pointer :: this

  sum_mass=0.
  this=>first
  do while(associated(this))
    if (present(justbergs)) then
      sum_mass=sum_mass+this%mass*this%mass_scaling
    elseif (present(justbits)) then
      sum_mass=sum_mass+this%mass_of_bits*this%mass_scaling
    else
      sum_mass=sum_mass+(this%mass+this%mass_of_bits)*this%mass_scaling
    endif
    this=>this%next
  enddo

end function sum_mass

! ##############################################################################

real function sum_heat(first,justbits,justbergs)
! Arguments
type(iceberg), pointer :: first
logical, intent(in), optional :: justbits, justbergs
! Local variables
type(iceberg), pointer :: this
real :: dm

  sum_heat=0.
  this=>first
  do while(associated(this))
    dm=0.
    if (present(justbergs)) then
      dm=this%mass*this%mass_scaling
    elseif (present(justbits)) then
      dm=this%mass_of_bits*this%mass_scaling
    else
      dm=(this%mass+this%mass_of_bits)*this%mass_scaling
    endif
    sum_heat=sum_heat+dm*this%heat_density
    this=>this%next
  enddo

end function sum_heat

! ##############################################################################

subroutine icebergs_stock_pe(bergs, index, value)
! Modules
use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT
! Arguments
type(icebergs), pointer :: bergs
integer, intent(in) :: index
real, intent(out) :: value
! Local variables
type(icebergs_gridded), pointer :: grd
real :: berg_mass, stored_mass

! For convenience
grd=>bergs%grd

select case (index)

  case (ISTOCK_WATER)
    berg_mass=sum_mass(bergs%first)
    stored_mass=sum( grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:) )
    value=stored_mass+berg_mass

  case (ISTOCK_HEAT)
    berg_mass=sum_mass(bergs%first)
    stored_mass=sum( grd%stored_ice(grd%isc:grd%iec,grd%jsc:grd%jec,:) )
    value=-(stored_mass+berg_mass)*HLF ! HLF is in (J/kg) from constants_mod

  case default
    value = 0.0

end select

end subroutine icebergs_stock_pe

! ##############################################################################

subroutine icebergs_save_restart(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables

  if (.not.associated(bergs)) return

  call mpp_clock_begin(bergs%clock_iow)
  call bergs_chksum(bergs, 'write_restart bergs')
  call write_restart(bergs)
  call mpp_clock_end(bergs%clock_iow)

end subroutine icebergs_save_restart

! ##############################################################################

subroutine icebergs_end(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
type(iceberg), pointer :: this, next

  if (.not.associated(bergs)) return

  call icebergs_save_restart(bergs)

  call mpp_clock_begin(bergs%clock_ini)
  ! Delete bergs and structures
  this=>bergs%first
  do while (associated(this))
    next=>this%next
    call move_trajectory(bergs, this)
    call destroy_iceberg(this)
    this=>next
  enddo

  if (associated(bergs%trajectories)) then
    call write_trajectory(bergs%trajectories)
  endif

  deallocate(bergs%grd%lon)
  deallocate(bergs%grd%lat)
  deallocate(bergs%grd%lonc)
  deallocate(bergs%grd%latc)
  deallocate(bergs%grd%dx)
  deallocate(bergs%grd%dy)
  deallocate(bergs%grd%area)
  deallocate(bergs%grd%msk)
  deallocate(bergs%grd%cos)
  deallocate(bergs%grd%sin)
  deallocate(bergs%grd%calving)
  deallocate(bergs%grd%calving_hflx)
  deallocate(bergs%grd%stored_heat)
  deallocate(bergs%grd%floating_melt)
  deallocate(bergs%grd%berg_melt)
  deallocate(bergs%grd%melt_buoy)
  deallocate(bergs%grd%melt_eros)
  deallocate(bergs%grd%melt_conv)
  deallocate(bergs%grd%bergy_src)
  deallocate(bergs%grd%bergy_melt)
  deallocate(bergs%grd%bergy_mass)
  deallocate(bergs%grd%virtual_area)
  deallocate(bergs%grd%mass)
  deallocate(bergs%grd%mass_on_ocean)
  deallocate(bergs%grd%tmp)
  deallocate(bergs%grd%tmpc)
  deallocate(bergs%grd%stored_ice)
  deallocate(bergs%grd%real_calving)
  deallocate(bergs%grd%uo)
  deallocate(bergs%grd%vo)
  deallocate(bergs%grd%ui)
  deallocate(bergs%grd%vi)
  deallocate(bergs%grd%ua)
  deallocate(bergs%grd%va)
  deallocate(bergs%grd%ssh)
  deallocate(bergs%grd%sst)
  deallocate(bergs%grd%cn)
  deallocate(bergs%grd%hi)
  deallocate(bergs%grd%domain)
  deallocate(bergs%grd)
  deallocate(bergs%initial_mass)
  deallocate(bergs%distribution)
  deallocate(bergs%initial_thickness)
  deallocate(bergs%initial_width)
  deallocate(bergs%initial_length)
  call dealloc_buffer(bergs%obuffer_n)
  call dealloc_buffer(bergs%obuffer_s)
  call dealloc_buffer(bergs%obuffer_e)
  call dealloc_buffer(bergs%obuffer_w)
  call dealloc_buffer(bergs%ibuffer_n)
  call dealloc_buffer(bergs%ibuffer_s)
  call dealloc_buffer(bergs%ibuffer_e)
  call dealloc_buffer(bergs%ibuffer_w)
  call mpp_clock_end(bergs%clock_ini)
  deallocate(bergs)

  if (mpp_pe()==mpp_root_pe()) write(*,'(a,i8)') 'diamonds: icebergs_end complete',mpp_pe()

  contains

  subroutine dealloc_buffer(buff)
  ! Arguments
  type(buffer), pointer :: buff
  ! Local variables
    if (associated(buff)) then
      if (associated(buff%data)) deallocate(buff%data)
      deallocate(buff)
    endif
  end subroutine dealloc_buffer

end subroutine icebergs_end

! ##############################################################################

subroutine invert_tau_for_du(u,v)
! Arguments
real, dimension(:,:),intent(inout) :: u, v
! Local variables
integer :: i, j
real :: cd, cddvmod, tau2

  cd=0.0015

  do j=lbound(u,2), ubound(u,2)
    do i=lbound(u,1), ubound(u,1)
      tau2=u(i,j)*u(i,j)+v(i,j)*v(i,j)
      cddvmod=sqrt(cd*sqrt(tau2))
      if (cddvmod.ne.0.) then
        u(i,j)=u(i,j)/cddvmod
        v(i,j)=v(i,j)/cddvmod
      else
        u(i,j)=0.
        v(i,j)=0.
      endif
    enddo
  enddo

end subroutine invert_tau_for_du

! ##############################################################################

subroutine sanitize_field(arr,val)
! Arguments
real, dimension(:,:),intent(inout) :: arr
real, intent(in) :: val
! Local variables
integer :: i, j

  do j=lbound(arr,2), ubound(arr,2)
    do i=lbound(arr,1), ubound(arr,1)
      if (abs(arr(i,j)).ge.val) arr(i,j)=0.
    enddo
  enddo

end subroutine sanitize_field

! ##############################################################################

subroutine write_restart(bergs)
! Arguments
type(icebergs), pointer :: bergs
! Local variables
integer :: iret, ncid, i_dim, i
integer :: lonid, latid, uvelid, vvelid, ineid, jneid
integer :: massid, thicknessid, lengthid, widthid
integer :: start_lonid, start_latid, start_yearid, start_dayid, start_massid
integer :: scaling_id, mass_of_bits_id, heat_density_id
character(len=35) :: filename
type(iceberg), pointer :: this
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  ! Only create a restart file for this PE if we have anything to say
  if (associated(bergs%first)) then

    call get_instance_filename("RESTART/icebergs.res.nc", filename)
    write(filename,'(A,".",I4.4)') trim(filename), mpp_pe()
    if (verbose) write(*,'(2a)') 'diamonds, write_restart: creating ',filename

    iret = nf_create(filename, NF_CLOBBER, ncid)
    if (iret .ne. NF_NOERR) write(stderrunit,*) 'diamonds, write_restart: nf_create failed'

    ! Dimensions
    iret = nf_def_dim(ncid, 'i', NF_UNLIMITED, i_dim)
    if (iret .ne. NF_NOERR) write(stderrunit,*) 'diamonds, write_restart: nf_def_dim i failed'

    ! Variables
    lonid = def_var(ncid, 'lon', NF_DOUBLE, i_dim)
    latid = def_var(ncid, 'lat', NF_DOUBLE, i_dim)
    uvelid = def_var(ncid, 'uvel', NF_DOUBLE, i_dim)
    vvelid = def_var(ncid, 'vvel', NF_DOUBLE, i_dim)
    massid = def_var(ncid, 'mass', NF_DOUBLE, i_dim)
    ineid = def_var(ncid, 'ine', NF_INT, i_dim)
    jneid = def_var(ncid, 'jne', NF_INT, i_dim)
    thicknessid = def_var(ncid, 'thickness', NF_DOUBLE, i_dim)
    widthid = def_var(ncid, 'width', NF_DOUBLE, i_dim)
    lengthid = def_var(ncid, 'length', NF_DOUBLE, i_dim)
    start_lonid = def_var(ncid, 'start_lon', NF_DOUBLE, i_dim)
    start_latid = def_var(ncid, 'start_lat', NF_DOUBLE, i_dim)
    start_yearid = def_var(ncid, 'start_year', NF_INT, i_dim)
    start_dayid = def_var(ncid, 'start_day', NF_DOUBLE, i_dim)
    start_massid = def_var(ncid, 'start_mass', NF_DOUBLE, i_dim)
    scaling_id = def_var(ncid, 'mass_scaling', NF_DOUBLE, i_dim)
    mass_of_bits_id = def_var(ncid, 'mass_of_bits', NF_DOUBLE, i_dim)
    heat_density_id = def_var(ncid, 'heat_density', NF_DOUBLE, i_dim)

    ! Attributes
    call put_att(ncid, lonid, 'long_name', 'longitude')
    call put_att(ncid, lonid, 'units', 'degrees_E')
    call put_att(ncid, latid, 'long_name', 'latitude')
    call put_att(ncid, latid, 'units', 'degrees_N')
    call put_att(ncid, uvelid, 'long_name', 'zonal velocity')
    call put_att(ncid, uvelid, 'units', 'm/s')
    call put_att(ncid, vvelid, 'long_name', 'meridional velocity')
    call put_att(ncid, vvelid, 'units', 'm/s')
    call put_att(ncid, ineid, 'long_name', 'i index')
    call put_att(ncid, ineid, 'units', 'none')
    call put_att(ncid, jneid, 'long_name', 'j index')
    call put_att(ncid, jneid, 'units', 'none')
    call put_att(ncid, massid, 'long_name', 'mass')
    call put_att(ncid, massid, 'units', 'kg')
    call put_att(ncid, thicknessid, 'long_name', 'thickness')
    call put_att(ncid, thicknessid, 'units', 'm')
    call put_att(ncid, widthid, 'long_name', 'width')
    call put_att(ncid, widthid, 'units', 'm')
    call put_att(ncid, lengthid, 'long_name', 'length')
    call put_att(ncid, lengthid, 'units', 'm')
    call put_att(ncid, start_lonid, 'long_name', 'longitude of calving location')
    call put_att(ncid, start_lonid, 'units', 'degrees_E')
    call put_att(ncid, start_latid, 'long_name', 'latitude of calving location')
    call put_att(ncid, start_latid, 'units', 'degrees_N')
    call put_att(ncid, start_yearid, 'long_name', 'calendar year of calving event')
    call put_att(ncid, start_yearid, 'units', 'years')
    call put_att(ncid, start_dayid, 'long_name', 'year day of calving event')
    call put_att(ncid, start_dayid, 'units', 'days')
    call put_att(ncid, start_massid, 'long_name', 'initial mass of calving berg')
    call put_att(ncid, start_massid, 'units', 'kg')
    call put_att(ncid, scaling_id, 'long_name', 'scaling factor for mass of calving berg')
    call put_att(ncid, scaling_id, 'units', 'none')
    call put_att(ncid, mass_of_bits_id, 'long_name', 'mass of bergy bits')
    call put_att(ncid, mass_of_bits_id, 'units', 'kg')
    call put_att(ncid, heat_density_id, 'long_name', 'heat density')
    call put_att(ncid, heat_density_id, 'units', 'J/kg')
    iret = nf_put_att_int(ncid, NF_GLOBAL, 'file_format_major_version', NF_INT, 1, file_format_major_version)
    iret = nf_put_att_int(ncid, NF_GLOBAL, 'file_format_minor_version', NF_INT, 3, file_format_minor_version)

    ! End define mode
    iret = nf_enddef(ncid)

    ! Write variables
   !this=>bergs%first; i=0
    this=>last_berg(bergs%first); i=0
    do while (associated(this))
      i=i+1
      call put_double(ncid, lonid, i, this%lon)
      call put_double(ncid, latid, i, this%lat)
      call put_double(ncid, uvelid, i, this%uvel)
      call put_double(ncid, vvelid, i, this%vvel)
      call put_int(ncid, ineid, i, this%ine)
      call put_int(ncid, jneid, i, this%jne)
      call put_double(ncid, massid, i, this%mass)
      call put_double(ncid, thicknessid, i, this%thickness)
      call put_double(ncid, widthid, i, this%width)
      call put_double(ncid, lengthid, i, this%length)
      call put_double(ncid, start_lonid, i, this%start_lon)
      call put_double(ncid, start_latid, i, this%start_lat)
      call put_int(ncid, start_yearid, i, this%start_year)
      call put_double(ncid, start_dayid, i, this%start_day)
      call put_double(ncid, start_massid, i, this%start_mass)
      call put_double(ncid, scaling_id, i, this%mass_scaling)
      call put_double(ncid, mass_of_bits_id, i, this%mass_of_bits)
      call put_double(ncid, heat_density_id, i, this%heat_density)
     !this=>this%next
      this=>this%prev
    enddo

    ! Finish up
    iret = nf_close(ncid)
    if (iret .ne. NF_NOERR) write(stderrunit,*) 'diamonds, write_restart: nf_close failed'

  endif ! associated(bergs%first)

  ! Write stored ice
  filename='RESTART/calving.res.nc'
  if (verbose.and.mpp_pe().eq.mpp_root_pe()) write(stderrunit,'(2a)') 'diamonds, write_restart: writing ',filename
  call grd_chksum3(bergs%grd, bergs%grd%stored_ice, 'write stored_ice')
  call write_data(filename, 'stored_ice', bergs%grd%stored_ice, bergs%grd%domain)
  call grd_chksum2(bergs%grd, bergs%grd%stored_heat, 'write stored_heat')
  call write_data(filename, 'stored_heat', bergs%grd%stored_heat, bergs%grd%domain)

  contains

  function last_berg(berg)
  ! Arguments
  type(iceberg), pointer :: last_berg, berg
  ! Local variables

    last_berg=>berg
    do while (associated(last_berg%next))
      last_berg=>last_berg%next
    enddo

  end function last_berg

end subroutine write_restart

! ##############################################################################

subroutine write_trajectory(trajectory)
! Arguments
type(xyt), pointer :: trajectory
! Local variables
integer :: iret, ncid, i_dim, i
integer :: lonid, latid, yearid, dayid, uvelid, vvelid
integer :: uoid, void, uiid, viid, uaid, vaid, sshxid, sshyid, sstid
integer :: cnid, hiid
integer :: mid, did, wid, lid, mbid, hdid
character(len=37) :: filename
character(len=7) :: pe_name
type(xyt), pointer :: this, next
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  call get_instance_filename("iceberg_trajectories.nc", filename)
  if (mpp_npes()>10000) then
     write(pe_name,'(a,i6.6)' )'.', mpp_pe()
  else
     write(pe_name,'(a,i4.4)' )'.', mpp_pe()
  endif
  filename=trim(filename)//trim(pe_name)
  if (debug) write(stderrunit,*) 'diamonds, write_trajectory: creating ',filename

  iret = nf_create(filename, NF_CLOBBER, ncid)
  if (iret .ne. NF_NOERR) write(stderrunit,*) 'diamonds, write_trajectory: nf_create failed'

  ! Dimensions
  iret = nf_def_dim(ncid, 'i', NF_UNLIMITED, i_dim)
  if (iret .ne. NF_NOERR) write(stderrunit,*) 'diamonds, write_trajectory: nf_def_dim i failed'

  ! Variables
  lonid = def_var(ncid, 'lon', NF_DOUBLE, i_dim)
  latid = def_var(ncid, 'lat', NF_DOUBLE, i_dim)
  yearid = def_var(ncid, 'year', NF_INT, i_dim)
  dayid = def_var(ncid, 'day', NF_DOUBLE, i_dim)
  uvelid = def_var(ncid, 'uvel', NF_DOUBLE, i_dim)
  vvelid = def_var(ncid, 'vvel', NF_DOUBLE, i_dim)
  uoid = def_var(ncid, 'uo', NF_DOUBLE, i_dim)
  void = def_var(ncid, 'vo', NF_DOUBLE, i_dim)
  uiid = def_var(ncid, 'ui', NF_DOUBLE, i_dim)
  viid = def_var(ncid, 'vi', NF_DOUBLE, i_dim)
  uaid = def_var(ncid, 'ua', NF_DOUBLE, i_dim)
  vaid = def_var(ncid, 'va', NF_DOUBLE, i_dim)
  mid = def_var(ncid, 'mass', NF_DOUBLE, i_dim)
  mbid = def_var(ncid, 'mass_of_bits', NF_DOUBLE, i_dim)
  hdid = def_var(ncid, 'heat_density', NF_DOUBLE, i_dim)
  did = def_var(ncid, 'thickness', NF_DOUBLE, i_dim)
  wid = def_var(ncid, 'width', NF_DOUBLE, i_dim)
  lid = def_var(ncid, 'length', NF_DOUBLE, i_dim)
  sshxid = def_var(ncid, 'ssh_x', NF_DOUBLE, i_dim)
  sshyid = def_var(ncid, 'ssh_y', NF_DOUBLE, i_dim)
  sstid = def_var(ncid, 'sst', NF_DOUBLE, i_dim)
  cnid = def_var(ncid, 'cn', NF_DOUBLE, i_dim)
  hiid = def_var(ncid, 'hi', NF_DOUBLE, i_dim)

  ! Attributes
  iret = nf_put_att_int(ncid, NF_GLOBAL, 'file_format_major_version', NF_INT, 1, 0)
  iret = nf_put_att_int(ncid, NF_GLOBAL, 'file_format_minor_version', NF_INT, 1, 1)
  call put_att(ncid, lonid, 'long_name', 'longitude')
  call put_att(ncid, lonid, 'units', 'degrees_E')
  call put_att(ncid, latid, 'long_name', 'latitude')
  call put_att(ncid, latid, 'units', 'degrees_N')
  call put_att(ncid, yearid, 'long_name', 'year')
  call put_att(ncid, yearid, 'units', 'years')
  call put_att(ncid, dayid, 'long_name', 'year day')
  call put_att(ncid, dayid, 'units', 'days')
  call put_att(ncid, uvelid, 'long_name', 'zonal spped')
  call put_att(ncid, uvelid, 'units', 'm/s')
  call put_att(ncid, vvelid, 'long_name', 'meridional spped')
  call put_att(ncid, vvelid, 'units', 'm/s')
  call put_att(ncid, uoid, 'long_name', 'ocean zonal spped')
  call put_att(ncid, uoid, 'units', 'm/s')
  call put_att(ncid, void, 'long_name', 'ocean meridional spped')
  call put_att(ncid, void, 'units', 'm/s')
  call put_att(ncid, uiid, 'long_name', 'ice zonal spped')
  call put_att(ncid, uiid, 'units', 'm/s')
  call put_att(ncid, viid, 'long_name', 'ice meridional spped')
  call put_att(ncid, viid, 'units', 'm/s')
  call put_att(ncid, uaid, 'long_name', 'atmos zonal spped')
  call put_att(ncid, uaid, 'units', 'm/s')
  call put_att(ncid, vaid, 'long_name', 'atmos meridional spped')
  call put_att(ncid, vaid, 'units', 'm/s')
  call put_att(ncid, mid, 'long_name', 'mass')
  call put_att(ncid, mid, 'units', 'kg')
  call put_att(ncid, mbid, 'long_name', 'mass_of_bits')
  call put_att(ncid, mbid, 'units', 'kg')
  call put_att(ncid, hdid, 'long_name', 'heat_density')
  call put_att(ncid, hdid, 'units', 'J/kg')
  call put_att(ncid, did, 'long_name', 'thickness')
  call put_att(ncid, did, 'units', 'm')
  call put_att(ncid, wid, 'long_name', 'width')
  call put_att(ncid, wid, 'units', 'm')
  call put_att(ncid, lid, 'long_name', 'length')
  call put_att(ncid, lid, 'units', 'm')
  call put_att(ncid, sshxid, 'long_name', 'sea surface height gradient_x')
  call put_att(ncid, sshxid, 'units', 'non-dim')
  call put_att(ncid, sshyid, 'long_name', 'sea surface height gradient_y')
  call put_att(ncid, sshyid, 'units', 'non-dim')
  call put_att(ncid, sstid, 'long_name', 'sea surface temperature')
  call put_att(ncid, sstid, 'units', 'degrees_C')
  call put_att(ncid, cnid, 'long_name', 'sea ice concentration')
  call put_att(ncid, cnid, 'units', 'none')
  call put_att(ncid, hiid, 'long_name', 'sea ice thickness')
  call put_att(ncid, hiid, 'units', 'm')

  ! End define mode
  iret = nf_enddef(ncid)

  ! Write variables
  this=>trajectory; i=0
  do while (associated(this))
    i=i+1
    call put_double(ncid, lonid, i, this%lon)
    call put_double(ncid, latid, i, this%lat)
    call put_int(ncid, yearid, i, this%year)
    call put_double(ncid, dayid, i, this%day)
    call put_double(ncid, uvelid, i, this%uvel)
    call put_double(ncid, vvelid, i, this%vvel)
    call put_double(ncid, uoid, i, this%uo)
    call put_double(ncid, void, i, this%vo)
    call put_double(ncid, uiid, i, this%ui)
    call put_double(ncid, viid, i, this%vi)
    call put_double(ncid, uaid, i, this%ua)
    call put_double(ncid, vaid, i, this%va)
    call put_double(ncid, mid, i, this%mass)
    call put_double(ncid, hdid, i, this%heat_density)
    call put_double(ncid, did, i, this%thickness)
    call put_double(ncid, wid, i, this%width)
    call put_double(ncid, lid, i, this%length)
    call put_double(ncid, sshxid, i, this%ssh_x)
    call put_double(ncid, sshyid, i, this%ssh_y)
    call put_double(ncid, sstid, i, this%sst)
    call put_double(ncid, cnid, i, this%cn)
    call put_double(ncid, hiid, i, this%hi)
    next=>this%next
    deallocate(this)
    this=>next
  enddo
  trajectory=>null()

  ! Finish up
  iret = nf_close(ncid)
  if (iret .ne. NF_NOERR) write(stderrunit,*) 'diamonds, write_trajectory: nf_close failed',mpp_pe(),filename

end subroutine write_trajectory


! ##############################################################################

integer function inq_var(ncid, var, unsafe)
! Arguments
integer, intent(in) :: ncid
character(len=*), intent(in) :: var
logical, optional, intent(in) :: unsafe
! Local variables
integer :: iret
integer :: stderrunit
logical :: unsafely=.false.

if(present(unsafe)) unsafely=unsafe
  ! Get the stderr unit number
  stderrunit=stderr()

  iret=nf_inq_varid(ncid, var, inq_var)
  if (iret .ne. NF_NOERR) then
    if (.not. unsafely) then
      write(stderrunit,*) 'diamonds, inq_var: nf_inq_varid ',var,' failed'
      call error_mesg('diamonds, inq_var', 'netcdf function returned a failure!', FATAL)
    else
      inq_var=-1
    endif
  endif

end function inq_var

! ##############################################################################

integer function def_var(ncid, var, ntype, idim)
! Arguments
integer, intent(in) :: ncid, ntype, idim
character(len=*), intent(in) :: var
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_def_var(ncid, var, ntype, 1, idim, def_var)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, def_var: nf_def_var failed for ',trim(var)
    call error_mesg('diamonds, def_var', 'netcdf function returned a failure!', FATAL)
  endif

end function def_var

! ##############################################################################

subroutine put_att(ncid, id, att, attval)
! Arguments
integer, intent(in) :: ncid, id
character(len=*), intent(in) :: att, attval
! Local variables
integer :: vallen, iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  vallen=len_trim(attval)
  iret = nf_put_att_text(ncid, id, att, vallen, attval)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, put_att: nf_put_att_text failed adding', &
      trim(att),' = ',trim(attval)
    call error_mesg('diamonds, put_att', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_att

! ##############################################################################

real function get_double(ncid, id, i)
! Arguments
integer, intent(in) :: ncid, id, i
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret=nf_get_var1_double(ncid, id, i, get_double)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, get_double: nf_get_var1_double failed reading'
    call error_mesg('diamonds, get_double', 'netcdf function returned a failure!', FATAL)
  endif

end function get_double

! ##############################################################################

integer function get_int(ncid, id, i)
! Arguments
integer, intent(in) :: ncid, id, i
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret=nf_get_var1_int(ncid, id, i, get_int)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, get_int: nf_get_var1_int failed reading'
    call error_mesg('diamonds, get_int', 'netcdf function returned a failure!', FATAL)
  endif

end function get_int

! ##############################################################################

subroutine put_double(ncid, id, i, val)
! Arguments
integer, intent(in) :: ncid, id, i
real, intent(in) :: val
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_put_vara_double(ncid, id, i, 1, val)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, put_double: nf_put_vara_double failed writing'
    call error_mesg('diamonds, put_double', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_double

! ##############################################################################

subroutine put_int(ncid, id, i, val)
! Arguments
integer, intent(in) :: ncid, id, i, val
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_put_vara_int(ncid, id, i, 1, val)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, put_int: nf_put_vara_int failed writing'
    call error_mesg('diamonds, put_int', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_int

! ##############################################################################

subroutine print_fld(grd, fld, label)
! Arguments
type(icebergs_gridded), pointer :: grd
real, intent(in) :: fld(grd%isd:grd%ied,grd%jsd:grd%jed)
character(len=*) :: label
! Local variables
integer :: i, j
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  write(stderrunit,'("pe=",i3,x,a8,32i10)') mpp_pe(),label,(i,i=grd%isd,grd%ied)
  do j=grd%jed,grd%jsd,-1
    write(stderrunit,'("pe=",i3,x,i8,32es10.2)') mpp_pe(),j,(fld(i,j),i=grd%isd,grd%ied)
  enddo

end subroutine print_fld

! ##############################################################################

subroutine checksum_gridded(grd, label)
! Arguments
type(icebergs_gridded), pointer :: grd
character(len=*) :: label
! Local variables

  if (mpp_pe().eq.mpp_root_pe()) write(*,'(2a)') 'diamonds: checksumming gridded data @ ',trim(label)

  ! external forcing
  call grd_chksum2(grd, grd%uo, 'uo')
  call grd_chksum2(grd, grd%vo, 'vo')
  call grd_chksum2(grd, grd%ua, 'ua')
  call grd_chksum2(grd, grd%va, 'va')
  call grd_chksum2(grd, grd%ui, 'ui')
  call grd_chksum2(grd, grd%vi, 'vi')
  call grd_chksum2(grd, grd%ssh, 'ssh')
  call grd_chksum2(grd, grd%sst, 'sst')
  call grd_chksum2(grd, grd%hi, 'hi')
  call grd_chksum2(grd, grd%cn, 'cn')
  call grd_chksum2(grd, grd%calving, 'calving')
  call grd_chksum2(grd, grd%calving_hflx, 'calving_hflx')

  ! state
  call grd_chksum2(grd, grd%mass, 'mass')
  call grd_chksum3(grd, grd%mass_on_ocean, 'mass_on_ocean')
  call grd_chksum3(grd, grd%stored_ice, 'stored_ice')
  call grd_chksum2(grd, grd%stored_heat, 'stored_heat')
  call grd_chksum2(grd, grd%melt_buoy, 'melt_b')
  call grd_chksum2(grd, grd%melt_eros, 'melt_e')
  call grd_chksum2(grd, grd%melt_conv, 'melt_v')
  call grd_chksum2(grd, grd%bergy_src, 'bergy_src')
  call grd_chksum2(grd, grd%bergy_melt, 'bergy_melt')
  call grd_chksum2(grd, grd%bergy_mass, 'bergy_mass')
  call grd_chksum2(grd, grd%virtual_area, 'varea')
  call grd_chksum2(grd, grd%floating_melt, 'floating_melt')
  call grd_chksum2(grd, grd%berg_melt, 'berg_melt')

  ! static
  call grd_chksum2(grd, grd%lon, 'lon')
  call grd_chksum2(grd, grd%lat, 'lat')
  call grd_chksum2(grd, grd%lonc, 'lonc')
  call grd_chksum2(grd, grd%latc, 'latc')
  call grd_chksum2(grd, grd%dx, 'dx')
  call grd_chksum2(grd, grd%dy, 'dy')
  call grd_chksum2(grd, grd%area, 'area')
  call grd_chksum2(grd, grd%msk, 'msk')
  call grd_chksum2(grd, grd%cos, 'cos')
  call grd_chksum2(grd, grd%sin, 'sin')

end subroutine checksum_gridded

! ##############################################################################

subroutine grd_chksum3(grd, fld, txt)
! Arguments
type(icebergs_gridded), pointer :: grd
real, dimension(:,:,:), intent(in) :: fld
character(len=*), intent(in) :: txt
! Local variables
integer :: i, j, k, halo, icount, io, jo
real :: mean, rms, SD, minv, maxv
real, dimension(lbound(fld,1):ubound(fld,1), lbound(fld,2):ubound(fld,2), lbound(fld,3):ubound(fld,3)) :: tmp

  halo=grd%halo
  mean=0.
  rms=0.
  sd=0.
  icount=0
  i=lbound(fld,1)+halo
  j=lbound(fld,2)+halo
  k=lbound(fld,3)
  minv=fld(i,j,k)
  maxv=fld(i,j,k)
  tmp(:,:,:)=0.
  io=grd%isd-lbound(fld,1)
  jo=grd%jsd-lbound(fld,2)
  do k=lbound(fld,3), ubound(fld,3)
    do j=lbound(fld,2)+halo, ubound(fld,2)-halo
      do i=lbound(fld,1)+halo, ubound(fld,1)-halo
        icount=icount+1
        mean=mean+fld(i,j,k)
        rms=rms+fld(i,j,k)**2
        minv=min(minv,fld(i,j,k))
        maxv=max(maxv,fld(i,j,k))
        tmp(i,j,k)=fld(i,j,k)*float(i+io+2*(j+jo)+3*(k-1))
      enddo
    enddo
  enddo
  call mpp_sum(icount)
  call mpp_sum(mean)
  call mpp_sum(rms)
  call mpp_min(minv)
  call mpp_max(maxv)
  mean=mean/float(icount)
  rms=sqrt(rms/float(icount))
  do k=lbound(fld,3), ubound(fld,3)
    do j=lbound(fld,2)+halo, ubound(fld,2)-halo
      do i=lbound(fld,1)+halo, ubound(fld,1)-halo
        sd=sd+(fld(i,j,k)-mean)**2
      enddo
    enddo
  enddo
  call mpp_sum(sd)
  sd=sqrt(sd/float(icount))
  i=mpp_chksum( fld(lbound(fld,1)+halo:ubound(fld,1)-halo, &
                    lbound(fld,2)+halo:ubound(fld,2)-halo,:) )
  j=mpp_chksum( tmp(lbound(fld,1)+halo:ubound(fld,1)-halo, &
                    lbound(fld,2)+halo:ubound(fld,2)-halo,:) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("diamonds, grd_chksum3: ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9))') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd
#ifdef CHECKSUM_HALOS
  i=mpp_chksum( fld(lbound(fld,1):ubound(fld,1), &
                    lbound(fld,2):ubound(fld,2),:) )
  j=mpp_chksum( tmp(lbound(fld,1):ubound(fld,1), &
                    lbound(fld,2):ubound(fld,2),:) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("diamonds, grd_chksum3* ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9))') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd
#endif

end subroutine grd_chksum3

! ##############################################################################

subroutine grd_chksum2(grd, fld, txt)
! Arguments
type(icebergs_gridded), pointer :: grd
real, dimension(grd%isd:grd%ied,grd%jsd:grd%jed), intent(in) :: fld
character(len=*), intent(in) :: txt
! Local variables
integer :: i, j, icount
real :: mean, rms, SD, minv, maxv

  grd%tmp(:,:)=0.

  mean=0.
  rms=0.
  sd=0.
  icount=0
  minv=fld(grd%isc,grd%jsc)
  maxv=fld(grd%isc,grd%jsc)
  do j=grd%jsc, grd%jec
    do i=grd%isc, grd%iec
      icount=icount+1
      mean=mean+fld(i,j)
      rms=rms+fld(i,j)**2
      minv=min(minv,fld(i,j))
      maxv=max(maxv,fld(i,j))
      grd%tmp(i,j)=fld(i,j)*float(i+2*j)
    enddo
  enddo
  call mpp_sum(icount)
  call mpp_sum(mean)
  call mpp_sum(rms)
  call mpp_min(minv)
  call mpp_max(maxv)
  mean=mean/float(icount)
  rms=sqrt(rms/float(icount))
  do j=grd%jsc, grd%jec
    do i=grd%isc, grd%iec
      sd=sd+(fld(i,j)-mean)**2
    enddo
  enddo
  call mpp_sum(sd)
  sd=sqrt(sd/float(icount))
  i=mpp_chksum( fld(grd%isc:grd%iec,grd%jsc:grd%jec) )
  j=mpp_chksum( grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("diamonds, grd_chksum2: ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9),x,a,"=",i8)') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd!, '#', icount
#ifdef CHECKSUM_HALOS
  i=mpp_chksum( fld(grd%isd:grd%ied,grd%jsd:grd%jed) )
  j=mpp_chksum( grd%tmp(grd%isd:grd%ied,grd%jsd:grd%jed) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("diamonds, grd_chksum2* ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9),x,a,"=",i)') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd!, '#', icount
#endif

end subroutine grd_chksum2

! ##############################################################################

subroutine bergs_chksum(bergs, txt, ignore_halo_violation)
! Arguments
type(icebergs), pointer :: bergs
character(len=*), intent(in) :: txt
logical, optional :: ignore_halo_violation
! Local variables
integer :: i, nbergs, ichk1, ichk2, ichk3, ichk4, ichk5, iberg
real, allocatable :: fld(:,:), fld2(:,:)
integer, allocatable :: icnt(:,:)
type(iceberg), pointer :: this
type(icebergs_gridded), pointer :: grd
logical :: check_halo

! For convenience
  grd=>bergs%grd

  nbergs=count_bergs(bergs)
  call mpp_max(nbergs)
  allocate( fld( nbergs, 11 ) )
  allocate( fld2( nbergs, 11 ) )
  allocate( icnt( grd%isd:grd%ied, grd%jsd:grd%jed ) )
  fld(:,:)=0.
  fld2(:,:)=0.
  icnt(:,:)=0
  grd%tmp(:,:)=0.

  this=>bergs%first
  i=0; ichk5=0
  do while(associated(this))
    i=i+1
    iberg=berg_chksum(this)
    fld(i,1) = this%lon
    fld(i,2) = this%lat
    fld(i,3) = this%uvel
    fld(i,4) = this%vvel
    fld(i,5) = this%mass
    fld(i,6) = this%thickness
    fld(i,7) = this%width
    fld(i,8) = this%length
    fld(i,9) = time_hash(this)
    fld(i,10) = pos_hash(this)
    fld(i,11) = float(iberg)
    icnt(this%ine,this%jne)=icnt(this%ine,this%jne)+1
    fld2(i,:) = fld(i,:)*float( icnt(this%ine,this%jne) ) !*float( i )
    grd%tmp(this%ine,this%jne)=grd%tmp(this%ine,this%jne)+time_hash(this)*pos_hash(this)+log(this%mass)
    ichk5=ichk5+iberg
    this=>this%next
  enddo

  ichk1=mpp_chksum( fld )
  ichk2=mpp_chksum( fld2 )
  ichk3=mpp_chksum( grd%tmp )
  ichk4=mpp_chksum( grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec) )
  call mpp_sum( ichk5 )
  nbergs=count_bergs(bergs)

  if (nbergs.ne.sum(icnt(:,:))) then
    write(*,'("diamonds, bergs_chksum: ",2(a,i8))') &
      '# bergs =', nbergs, ' sum(icnt) =',sum(icnt(:,:))
    call error_mesg('diamonds, bergs_chksum:', 'mismatch in berg count!', FATAL)
  endif

  check_halo=.true.
  if (present(ignore_halo_violation)) then
    if (ignore_halo_violation) check_halo=.false.
  endif
  if (check_halo.and.nbergs.ne.sum(icnt(grd%isc:grd%iec, grd%jsc:grd%jec))) then
    write(*,'("diamonds, bergs_chksum: ",2(a,i8))') &
      '# bergs =', nbergs, ' sum(icnt(comp_dom)) =',sum(icnt(:,:))
    call error_mesg('diamonds, bergs_chksum:', 'mismatch in berg count on computational domain!', FATAL)
  endif

  call mpp_sum(nbergs)
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("diamonds, bergs_chksum: ",a18,6(x,a,"=",i22))') &
      txt, 'chksum', ichk1, 'chksum2', ichk2, 'chksum3', ichk3, 'chksum4', ichk4, 'chksum5', ichk5, '#', nbergs

  grd%tmp(:,:)=real(icnt(:,:))
  call grd_chksum2(grd,grd%tmp,'# of bergs/cell')

  deallocate( fld )
  deallocate( fld2 )
  deallocate( icnt )

  if (debug) call count_out_of_order(bergs,txt)

end subroutine bergs_chksum

! ##############################################################################

integer function berg_chksum(berg )
! Arguments
type(iceberg), pointer :: berg
! Local variables
real :: rtmp(28)
integer :: itmp(28+3), i8=0, ichk1, ichk2, ichk3
integer :: i

  rtmp(:)=0.
  rtmp(1)=berg%lon
  rtmp(2)=berg%lat
  rtmp(3)=berg%uvel
  rtmp(4)=berg%vvel
  rtmp(5)=berg%mass
  rtmp(6)=berg%thickness
  rtmp(7)=berg%width
  rtmp(8)=berg%length
  rtmp(9)=berg%start_lon
  rtmp(10)=berg%start_lat
  rtmp(11)=berg%start_day
  rtmp(12)=berg%start_mass
  rtmp(13)=berg%mass_scaling
  rtmp(14)=berg%mass_of_bits
  rtmp(15)=berg%heat_density
  rtmp(16)=berg%xi
  rtmp(17)=berg%yj
  rtmp(19)=berg%uo
  rtmp(20)=berg%vo
  rtmp(21)=berg%ui
  rtmp(22)=berg%vi
  rtmp(23)=berg%ua
  rtmp(24)=berg%va
  rtmp(25)=berg%ssh_x
  rtmp(26)=berg%ssh_y
  rtmp(27)=berg%cn
  rtmp(28)=berg%hi

  itmp(1:28)=transfer(rtmp,i8)
  itmp(29)=berg%start_year
  itmp(30)=berg%ine
  itmp(31)=berg%jne

  ichk1=0; ichk2=0; ichk3=0
  do i=1,28+3
   ichk1=ichk1+itmp(i)
   ichk2=ichk2+itmp(i)*i
   ichk3=ichk3+itmp(i)*i*i
  enddo
  berg_chksum=ichk1+ichk2+ichk3

end function berg_chksum

! ##############################################################################

logical function find_restart_file(filename, actual_file, multiPErestart)
  character(len=*), intent(in) :: filename
  character(len=*), intent(out) :: actual_file
  logical, intent(out) :: multiPErestart

  character(len=6) :: pe_name

  find_restart_file = .false.

  ! If running as ensemble, add the ensemble id string to the filename
  call get_instance_filename(filename, actual_file)

  ! Prefer combined restart files.
  inquire(file=actual_file,exist=find_restart_file)
  if (find_restart_file) return

  ! Uncombined restart
  if (mpp_npes()>10000) then
     write(pe_name,'(a,i6.6)' )'.', mpp_pe()
  else
     write(pe_name,'(a,i4.4)' )'.', mpp_pe()
  endif
  actual_file=trim(actual_file)//trim(pe_name)
  inquire(file=actual_file,exist=find_restart_file)
  if (find_restart_file) then
     multiPErestart=.true.
     return
  endif

  ! No file found, Reset all return parameters
  find_restart_file=.false.
  actual_file = ''
  multiPErestart=.false.

end function find_restart_file

end module

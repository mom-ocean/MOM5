
module bgrid_diagnostics_mod

!-----------------------------------------------------------------------

use       bgrid_horiz_mod, only: horiz_grid_type
use        bgrid_vert_mod, only: vert_grid_type, compute_pres_full,  &
                                 compute_pres_half, compute_pres_depth
use       bgrid_masks_mod, only: grid_mask_type
use    bgrid_prog_var_mod, only: prog_var_type
use bgrid_change_grid_mod, only: change_grid, TEMP_GRID, WIND_GRID

use      diag_manager_mod, only: diag_axis_init, register_diag_field, &
                                 register_static_field, send_data
use      time_manager_mod, only: time_type

use            fms_mod, only: file_exist, open_namelist_file,    &
                              error_mesg, NOTE, check_nml_error, &
                              mpp_pe, mpp_root_pe, stdlog,       &
                              close_file, write_version_number
use      constants_mod, only: GRAV, KAPPA, RDGAS
use  field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_names, get_number_tracers



implicit none
private

public :: bgrid_diagnostics_init, &
          bgrid_diagnostics,      &
          bgrid_diagnostics_tend

!-----------------------------------------------------------------------
!------------------------- axis names ----------------------------------

character(len=8) :: axiset = 'dynamics'
character(len=8) :: mod_name = 'dynamics'

!-----------------------------------------------------------------------

   integer, parameter :: MXTR = 10
   real,    parameter :: GINV = 1./GRAV

!-----------------------------------------------------------------------
! axis and field identifiers for the diag manager

integer :: id_hlonb, id_hlon , id_hlatb, id_hlat , &
           id_vlonb, id_vlon , id_vlatb, id_vlat , &
           id_phalf, id_pfull, id_hlat_wgt, id_vlat_wgt

integer :: id_bk   , id_pk   , id_zsurf, id_res  , id_wspd,          &
           id_ps   , id_ucomp, id_vcomp, id_temp , id_pres_full,     &
           id_omega, id_div  , id_vor  , id_pgfx , id_pgfy,          &
           id_udt  , id_vdt  , id_tdt  , id_pres_half,               &
           id_alm  , id_aph  , id_theta, id_mfew,  id_mfns,  id_slp
integer, allocatable :: id_tracer(:), id_tracer_tend(:)

integer :: id_ucomp_sq, id_vcomp_sq, id_temp_sq, id_omega_sq, &
           id_ucomp_vcomp, id_omega_temp

!-----------------------------------------------------------------------
! need to save surface geopotential height argument to initialization call
! the surface height will be needed for computing sea level pressure
real, pointer, dimension(:,:) :: zsurfg
!-----------------------------------------------------------------------

 character(len=128) :: version = '$Id: bgrid_diagnostics.F90,v 19.0 2012/01/06 19:54:38 fms Exp $'
 character(len=128) :: tag = '$Name: tikal $'

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine bgrid_diagnostics_init ( Time, Hgrid, Vgrid, Var, &
                                    fis, res,                &
                                    mass_axes, vel_axes      )

!-----------------------------------------------------------------------
!             setup/write netcdf metadata and static fields
!-----------------------------------------------------------------------
! Time      = current/initial time
! Hgrid     = horizontal grid constants
! Vgrid     = vertical grid constants
! fis       = geopotential height of the surface
! res       = reciprocal of eta at the surface
! mass_axes = axis identifiers for the temperature (mass) grid
! vel_axes  = axis identifiers for the velocity grid
!   NOTE:  The axes identifiers are for the lon, lat, pfull, and
!       phalf axes, respectively. They are returned by the diag_manager.
!-----------------------------------------------------------------------

   type(time_type),       intent(in)  :: Time
   type(horiz_grid_type), intent(in)  :: Hgrid
   type (vert_grid_type), intent(in)  :: Vgrid
   type  (prog_var_type), intent(in)  :: Var
   real, intent(in), dimension(Hgrid%ilb:,Hgrid%jlb:), target  :: fis, res
   integer, dimension(4), intent(out) :: mass_axes,  vel_axes

!-----------------------------------------------------------------------
      real, dimension(Hgrid%Tmp%isg:Hgrid%Tmp%ieg+1) :: hlonb
      real, dimension(Hgrid%Vel%isg:Hgrid%Vel%ieg+1) :: vlonb
      real, dimension(Hgrid%Tmp%jsg:Hgrid%Tmp%jeg+1) :: hlatb
      real, dimension(Hgrid%Vel%jsg:Hgrid%Vel%jeg+1) :: vlatb

      real, dimension(Hgrid%Tmp%isg:Hgrid%Tmp%ieg) :: hlon
      real, dimension(Hgrid%Vel%isg:Hgrid%Vel%ieg) :: vlon
      real, dimension(Hgrid%Tmp%jsg:Hgrid%Tmp%jeg) :: hlat
      real, dimension(Hgrid%Vel%jsg:Hgrid%Vel%jeg) :: vlat

      real, dimension(Hgrid%Tmp%js:Hgrid%Tmp%je) :: hlat_wgt
      real, dimension(Hgrid%Vel%js:Hgrid%Vel%je) :: vlat_wgt

      real, dimension(1,1)              :: psurf
      real, dimension(1,1,Vgrid%nlev)   :: pfull
      real, dimension(1,1,Vgrid%nlev+1) :: phalf


      real    :: vrange(2), trange(2), prange(2)
      real    :: rad2deg
      integer :: i, j, n, unit, io, ierr, ntprog
      integer :: isg, ieg, hsg, heg, vsg, veg
      integer :: is, ie, hs, he, vs, ve
      integer :: uflx_axes(4), vflx_axes(4)
      integer :: logunit
      logical :: used
      character(len=128) :: tname
      character(len=256) :: longname, units

!--------------------------- set up axes -------------------------------

      ! compute grid indices
      is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie
      hs = Hgrid % Tmp % js;  he = Hgrid % Tmp % je
      vs = Hgrid % Vel % js;  ve = Hgrid % Vel % je

      ! global grid indices
      isg = Hgrid % Tmp % isg;  ieg = Hgrid % Tmp % ieg
      hsg = Hgrid % Tmp % jsg;  heg = Hgrid % Tmp % jeg
      vsg = Hgrid % Vel % jsg;  veg = Hgrid % Vel % jeg

      ! grid box boundaries in degrees
      rad2deg = 90./acos(0.0)
      hlonb(isg:ieg+1) = Hgrid % Tmp % blong(isg:ieg+1) * rad2deg
      hlatb(hsg:heg+1) = Hgrid % Tmp % blatg(hsg:heg+1) * rad2deg
      vlonb(isg:ieg+1) = Hgrid % Vel % blong(isg:ieg+1) * rad2deg
      vlatb(vsg:veg+1) = Hgrid % Vel % blatg(vsg:veg+1) * rad2deg

      logunit = stdlog()
      ! grid box centers in degrees
   do i = isg, ieg
      hlon(i) = 0.5*(hlonb(i)+hlonb(i+1))
      vlon(i) = 0.5*(vlonb(i)+vlonb(i+1))
   enddo

   do j = hsg, heg
      hlat(j) = 0.5*(hlatb(j)+hlatb(j+1))
   enddo

   do j = vsg, veg
      vlat(j) = 0.5*(vlatb(j)+vlatb(j+1))
   enddo

  ! compute a reference profile of pressure based on psurf = 1000 hPa
    psurf = reshape ( (/ 100000. /), (/ 1, 1 /) )
    call compute_pres_full (Vgrid, psurf, pfull)
    call compute_pres_half (Vgrid, psurf, phalf)
  ! in units of hPa
    pfull = pfull*0.01
    phalf = phalf*0.01


!----- initialize mass axes ------

 id_hlonb = diag_axis_init ( 'lonb', hlonb, 'degrees_E', 'x',     &
                             'longitude edges', set_name='atmos', &
                             Domain2=Hgrid%Tmp%Domain_nohalo      )

 id_hlon  = diag_axis_init ( 'lon', hlon, 'degrees_E', 'x',       &
                             'longitude', set_name='atmos',       &
                             edges=id_hlonb,                      &
                             Domain2=Hgrid%Tmp%Domain_nohalo      )

 id_hlatb = diag_axis_init ( 'latb', hlatb, 'degrees_N', 'y',     &
                             'latitude edges', set_name='atmos',  &
                             Domain2=Hgrid%Tmp%Domain_nohalo      )

 id_hlat  = diag_axis_init ( 'lat', hlat, 'degrees_N', 'y',       &
                             'latitude', set_name='atmos',        &
                             edges=id_hlatb,                      &
                             Domain2=Hgrid%Tmp%Domain_nohalo      )

!----- initialize velocity axes ------

 id_vlonb = diag_axis_init ( 'vlonb', vlonb, 'degrees_E', 'x',    &
                             'longitude edges', set_name='atmos', &
                             Domain2=Hgrid%Vel%Domain_nohalo      )

 id_vlon  = diag_axis_init ( 'vlon', vlon, 'degrees_E', 'x',      &
                             'longitude', set_name='atmos',       &
                             edges=id_vlonb,                      &
                             Domain2=Hgrid%Vel%Domain_nohalo      )

 id_vlatb = diag_axis_init ( 'vlatb', vlatb, 'degrees_N', 'y',    &
                             'latitude edges', set_name='atmos',  &
                             Domain2=Hgrid%Vel%Domain_nohalo      )

 id_vlat  = diag_axis_init ( 'vlat', vlat, 'degrees_N', 'y',      &
                             'latitude', set_name='atmos',        &
                             edges=id_vlatb,                      &
                             Domain2=Hgrid%Vel%Domain_nohalo      )

!----- initialize vertical axes -----

 id_phalf = diag_axis_init ( 'phalf', phalf(1,1,:), 'hPa', 'z', &
                             'approx half pressure level',      &
                             direction=-1, set_name='atmos'     )

 id_pfull = diag_axis_init ( 'pfull', pfull(1,1,:), 'hPa', 'z', &
                             'approx full pressure level',      &
                             direction=-1, edges=id_phalf,      &
                             set_name='atmos'                   )

!-----------------------------------------------------------------------
!-------- initialize and output variables with no time axis ------------

    mass_axes = (/ id_hlon, id_hlat, id_pfull, id_phalf /)
     vel_axes = (/ id_vlon, id_vlat, id_pfull, id_phalf /)
    uflx_axes = (/ id_vlon, id_hlat, id_pfull, id_phalf /)
    vflx_axes = (/ id_hlon, id_vlat, id_pfull, id_phalf /)

   ! valid range for some fields
     vrange = (/ -400., +400. /)  ! momentum
     trange = (/  100.,  400. /)  ! temperature
     prange = (/ -1., 107500. /)  ! pressure

!-----------------------------------------------------------------------
!---- register static fields -------

  id_bk    = register_static_field ( mod_name, 'bk', (/id_phalf/), &
                        'vertical coordinate sigma value', 'none' )

  id_pk    = register_static_field ( mod_name, 'pk', (/id_phalf/), &
   'vertical coordinate reference pressure value (ak*pref)', 'pascals' )

  id_zsurf = register_static_field ( mod_name, 'zsurf', mass_axes(1:2),&
                                       'surface height', 'm' )

  id_res   = register_static_field ( mod_name, 'res', mass_axes(1:2), &
                     'reciprocal of sigma/eta at the surface', 'none' )

  id_alm   = register_static_field ( mod_name, 'alm', mass_axes(1:2), &
                'actual longitudes for temperature grid', 'degrees_E' )

  id_aph   = register_static_field ( mod_name, 'aph', mass_axes(1:2), &
                 'actual latitudes for temperature grid', 'degrees_N' )

! these changes cannot be implemented until changes to diag_manager
! initialize fields useful for computing offline global averages
!
! id_hlat_wgt = register_static_field ( mod_name, 'lat_wgt',  &
!                 (/id_hlat/), 'latitude weight for mass grid', 'none' )
!
! id_vlat_wgt = register_static_field ( mod_name, 'vlat_wgt',  &
!             (/id_vlat/), 'latitude weight for momentum grid', 'none' )

      if ( id_bk > 0 ) &
      used = send_data ( id_bk, Vgrid%eta, Time )

      if ( id_pk > 0 ) &
      used = send_data ( id_pk, Vgrid%peta, Time )

      if ( id_zsurf > 0 ) &
      used = send_data ( id_zsurf, fis(is:ie,hs:he)*GINV, Time )

      if ( id_res > 0 ) &
      used = send_data ( id_res, res(is:ie,hs:he), Time )

      if ( id_alm > 0 ) &
      used = send_data ( id_alm, Hgrid%Tmp%alm(is:ie,hs:he)*rad2deg, Time )

      if ( id_aph > 0 ) &
      used = send_data ( id_aph, Hgrid%Tmp%aph(is:ie,hs:he)*rad2deg, Time )

!     if ( id_hlat_wgt > 0 ) then
!        hlat_wgt = sin(Hgrid%Tmp%blatg(hs+1:he+1))-sin(Hgrid%Tmp%blatg(hs:he))
!        used = send_data ( id_hlat_wgt, hlat_wgt, Time )
!     endif
!
!     if ( id_vlat_wgt > 0 ) then
!        vlat_wgt = sin(Hgrid%Vel%blatg(vs+1:ve+1))-sin(Hgrid%Vel%blatg(vs:ve))
!        used = send_data ( id_vlat_wgt, vlat_wgt, Time )
!     endif

!---- register non-static fields -------

   id_ps   = register_diag_field ( mod_name, 'ps', mass_axes(1:2), &
                               Time, 'surface pressure', 'pascals' )

   id_slp  = register_diag_field ( mod_name, 'slp', mass_axes(1:2), &
                               Time, 'sea level pressure', 'pascals' )

   id_ucomp = register_diag_field ( mod_name, 'ucomp', vel_axes(1:3), &
                           Time, 'zonal wind component', 'm/sec',     &
                           missing_value=vrange(1), range=vrange      )

   id_vcomp = register_diag_field ( mod_name, 'vcomp', vel_axes(1:3), &
                        Time, 'meridional wind component', 'm/sec',   &
                        missing_value=vrange(1), range=vrange         )

   id_temp = register_diag_field ( mod_name, 'temp', mass_axes(1:3), &
                             Time, 'temperature', 'deg_k',           &
                             missing_value=trange(1), range=trange   )

   id_pres_full = register_diag_field ( mod_name, 'pres_full', mass_axes(1:3), &
                  Time, 'pressure at full model levels', 'pascals',  &
                              missing_value=prange(1), range=prange  )

   ! pressure at half levels
   id_pres_half = register_diag_field ( mod_name, 'pres_half',       &
                             (/ id_hlon, id_hlat, id_phalf /), Time, &
                        'pressure at half model levels', 'pascals',  &
                              missing_value=prange(1), range=prange  )

   id_omega = register_diag_field ( mod_name, 'omega', mass_axes(1:3),&
                                 Time, 'omega vertical velocity',     &
                                 'pascals/sec',                       &
                                 missing_value=-999.                  )

   id_theta = register_diag_field ( mod_name, 'theta', mass_axes(1:3), &
                             Time, 'potential temperature', 'deg_k',   &
                             missing_value=-999.                       )

   id_mfew  = register_diag_field ( mod_name, 'mfew', uflx_axes(1:3),  &
                             Time, 'Zonal mass flux', 'Pa-m2/s',       &
                             missing_value=-1.e30                      )

   id_mfns  = register_diag_field ( mod_name, 'mfns', vflx_axes(1:3),  &
                             Time, 'Meridional mass flux', 'Pa-m2/s',  &
                             missing_value=-1.e30                      )

 !  write version (to log file) 
    call write_version_number (version,tag)

 ! register diagnostics for all tracers
   allocate (id_tracer(Var%ntrace))
   if (mpp_pe() == mpp_root_pe()) write(logunit,100) trim(mod_name)
   do n = 1, Var%ntrace
     call get_tracer_names ( MODEL_ATMOS, n, tname, longname, units )
     if (mpp_pe() == mpp_root_pe()) write(logunit,110) trim(tname),trim(longname),trim(units)
     id_tracer(n) = register_diag_field ( mod_name, trim(tname),  &
                            mass_axes(1:3), Time, trim(longname), &
                            trim(units), missing_value=-999.      )
   enddo
100 format ('Diagnostics for the following tracer fields are available for module name = ',a)
110 format (3x,a,' (',a,'; ',a,')')

!-------- register second-moment quantities -------
! (for now we are only saving fields on the same grids)

   id_ucomp_sq = register_diag_field ( mod_name, 'ucomp_sq', vel_axes(1:3), &
                           Time, 'zonal wind component squared', 'm2/s2',   &
                           missing_value=-1., range=(/0.,vrange(2)**2/)   )

   id_vcomp_sq = register_diag_field ( mod_name, 'vcomp_sq', vel_axes(1:3), &
                      Time, 'meridional wind component squared', 'm2/s2',   &
                       missing_value=-1., range=(/0.,vrange(2)**2/)   )

   id_temp_sq = register_diag_field ( mod_name, 'temp_sq', mass_axes(1:3), &
                             Time, 'temperature squared', 'deg_K**2',      &
                             missing_value=-1., range=(/0.,trange(2)**2/)  )

   id_omega_sq = register_diag_field ( mod_name, 'omega_sq', mass_axes(1:3),&
                                 Time, 'omega vertical velocity squared',   &
                                 'Pa**2/s**2', missing_value=-999.          )

   id_ucomp_vcomp = register_diag_field ( mod_name, 'ucomp_vcomp', vel_axes(1:3),&
                       Time, 'zonal times meridional wind components', 'm2/s2',  &
                       missing_value=-1. )

   id_omega_temp = register_diag_field ( mod_name, 'omega_temp', mass_axes(1:3),&
                               Time, 'omega vertical velocity time temperature',&
                                 'Pascals*deg_K/sec', missing_value=-999.       )

!-------- wind speed, divergence, vorticity ----------------------------

   id_wspd = register_diag_field ( mod_name, 'wspd', vel_axes(1:3),   &
                       Time, 'wind speed', 'm/s', missing_value=-999.,&
                       range=(/0.,vrange(2)/) )

   id_div  = register_diag_field ( mod_name, 'div', mass_axes(1:3),   &
                       Time, 'divergence', '1/s', missing_value=-999. )

   id_vor  = register_diag_field ( mod_name, 'vor', mass_axes(1:3),   &
               Time, 'relative vorticity', '1/s', missing_value=-999. )

!--------- pressure gradient components (NOT USED) ---------------------

!  id_pgfx = register_diag_field ( mod_name, 'pgfx', vel_axes(1:3), &
!                            Time, 'zonal pressure gradient force', &
!                            'm/s2', missing_value=-999.            )

!  id_pgfy = register_diag_field ( mod_name, 'pgfy', vel_axes(1:3), &
!                       Time, 'meridional pressure gradient force', &
!                       'm/s2', missing_value=-999.                 )

!-----------------------------------------------------------------------
!         -------- tendencies ---------

   id_udt = register_diag_field ( mod_name, 'udt_dyn', vel_axes(1:3),  &
                            Time, 'zonal wind tendency for dynamics', &
                            'm/s2', missing_value=-999. )

   id_vdt = register_diag_field ( mod_name, 'vdt_dyn', vel_axes(1:3),  &
                        Time, 'meridional wind tendency for dynamics', &
                        'm/s2', missing_value=-999. )

   id_tdt = register_diag_field ( mod_name, 'tdt_dyn', mass_axes(1:3), &
                            Time, 'temperature tendency for dynamics', &
                            'deg_k/sec', missing_value=-999. )

 ! tendencies for prognostic tracers only
   call get_number_tracers ( MODEL_ATMOS, num_prog=ntprog )
   allocate (id_tracer_tend(ntprog))
   do n = 1, ntprog
     call get_tracer_names ( MODEL_ATMOS, n, tname, longname, units )
     tname    = trim(tname)   //'_dt_dyn'
     longname = trim(longname)//' tendency for dynamics'
     units    = trim(units)   //'/s'
     if (units == 'none') units = '1/sec'
     if (mpp_pe() == mpp_root_pe()) write(logunit,110) trim(tname),trim(longname),trim(units)
     id_tracer_tend(n) = register_diag_field ( mod_name, trim(tname),  &
                                 mass_axes(1:3), Time, trim(longname), &
                                 trim(units), missing_value=-999.      )
   enddo

 ! save surface geopotential height for computing sea level pressure
   zsurfg => fis

!-----------------------------------------------------------------------

   end subroutine bgrid_diagnostics_init

!#######################################################################

 subroutine bgrid_diagnostics ( Hgrid, Vgrid, Var, Masks, Time, &
                                omega, div, mfew, mfns )

!-----------------------------------------------------------------------
!           write netcdf fields
!-----------------------------------------------------------------------
!  Hgrid  = horizontal grid constants
!  Vgrid  = vertical grid constants
!  Var    = prognostic variables at diagnostics Time
!  Masks  = grid box masks for step-mountain topography
!  Time   = diagnostics time
!  omega  = omega (vertical velocity) diagnostic
!-----------------------------------------------------------------------
type(horiz_grid_type), intent(in) :: Hgrid
type (vert_grid_type), intent(in) :: Vgrid
type (prog_var_type),  intent(in) :: Var
type(grid_mask_type),  intent(in) :: Masks
type(time_type),       intent(in) :: Time
   real, intent(in), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: &
                            omega, div, mfew, mfns

!  real, intent(in), dimension(Hgrid%ilb:,Hgrid%jlb:,:), optional ::&
!                                               div, pgfx, pgfy
                                                       
   real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: slp
   real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub, &
                   Vgrid%nlev) :: wspd, vor, dp, udp, vdp
   real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub, &
                   Vgrid%nlev+1) :: ph
logical, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub, &
                   Vgrid%nlev+1) :: lmask
!-----------------------------------------------------------------------
   integer :: is, ie, hs, he, vs, ve, n, j, k
   logical :: used
!-----------------------------------------------------------------------

      is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie
      hs = Hgrid % Tmp % js;  he = Hgrid % Tmp % je
      vs = Hgrid % Vel % js;  ve = Hgrid % Vel % je

!-----------------------------------------------------------------------
!---------------- surface fields ---------------------------------------

      if ( id_ps > 0 ) &
      used = send_data ( id_ps , Var%ps(is:ie,hs:he), Time )

!---------------- 3d momentum fields (u & v) ---------------------------

      if ( id_ucomp > 0 ) &
      used = send_data ( id_ucomp, Var%u(is:ie,vs:ve,:), Time,    &
                         mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

      if ( id_vcomp > 0 ) &
      used = send_data ( id_vcomp, Var%v(is:ie,vs:ve,:), Time,    &
                         mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

      if ( id_temp > 0 ) &
      used = send_data ( id_temp, Var%t(is:ie,hs:he,:), Time,     &
                         mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )

      do n = 1, Var%ntrace
        if ( id_tracer(n) > 0 ) &
        used = send_data ( id_tracer(n), Var%r(is:ie,hs:he,:,n), Time, &
                           mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5    )
      enddo

      if ( id_omega > 0 ) &
      used = send_data ( id_omega, omega(is:ie,hs:he,:), Time,    &
                         mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )

      ! pressure at full levels
      ! Note: not computational efficient to recompute pfull
      if ( id_pres_full > 0 ) then
         call compute_pres_full (Vgrid, Var%pssl, dp)
         used = send_data ( id_pres_full, dp(is:ie,hs:he,:), Time,   &
                            mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )
      endif

      ! pressure at half (interface) levels
      ! Note: not computational efficient to recompute phalf
      if ( id_pres_half > 0 ) then
         call compute_pres_half (Vgrid, Var%pssl, ph)
         lmask(is:ie,hs:he,1) = .true.
         lmask(is:ie,hs:he,2:Vgrid%nlev+1) = Masks%Tmp%mask(is:ie,hs:he,:) > 0.5
         used = send_data ( id_pres_half, ph(is:ie,hs:he,:), Time,   &
                            mask=lmask(is:ie,hs:he,:) )
      endif

     !--- sea level pressure ---
      if ( id_slp > 0 ) then
         if ( id_pres_full <= 0 ) call compute_pres_full (Vgrid, Var%pssl, dp)
         call sea_level_pressure ( Var%ps, zsurfg, dp, Var%t, slp )
         used = send_data ( id_slp, slp(is:ie,hs:he), Time )
      endif

     ! potential temperature (compute pfull if necessary)
      if ( id_theta > 0 ) then
          if ( id_pres_full <= 0 .and. id_slp <= 0 ) call compute_pres_full (Vgrid, Var%pssl, dp)
          dp = Var%t * (1000.e2/dp)**KAPPA
          used = send_data ( id_theta, dp(is:ie,hs:he,:), Time,       &
                             mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )
      endif

!--------- second moment quantities ----------

      if ( id_ucomp_sq > 0 ) &
      used = send_data ( id_ucomp_sq, Var%u(is:ie,vs:ve,:)**2, Time, &
                         mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

      if ( id_vcomp_sq > 0 ) &
      used = send_data ( id_vcomp_sq, Var%v(is:ie,vs:ve,:)**2, Time, &
                         mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

      if ( id_temp_sq > 0 ) &
      used = send_data ( id_temp_sq, Var%t(is:ie,hs:he,:)**2, Time, &
                         mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )

      if ( id_omega_sq > 0 ) &
      used = send_data ( id_omega_sq, omega(is:ie,hs:he,:)**2, Time, &
                         mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )

      if ( id_ucomp_vcomp > 0 ) used = send_data ( id_ucomp_vcomp, &
                  Var%u(is:ie,vs:ve,:)*Var%v(is:ie,vs:ve,:), Time, &
                         mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

      if ( id_omega_temp > 0 ) used = send_data ( id_omega_temp, &
                omega(is:ie,hs:he,:)*Var%t(is:ie,hs:he,:), Time, &
                        mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )

!------------ wind speed, divergence, vorticity ------------------------

      if ( id_wspd > 0 ) then
          wspd(is:ie,vs:ve,:) = sqrt &
                      ( Var%u(is:ie,vs:ve,:)*Var%u(is:ie,vs:ve,:) + &
                        Var%v(is:ie,vs:ve,:)*Var%v(is:ie,vs:ve,:) )
          used = send_data ( id_wspd, wspd(is:ie,vs:ve,:), Time,      &
                             mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )
      endif

      if ( id_div > 0 ) then
          used = send_data ( id_div, div(is:ie,hs:he,:), Time,      &
                           mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )
      endif

      if ( id_vor > 0 ) then
     !if ( id_vor > 0 .or. id_div > 0 ) then
     !--- precompute quantities common to both vor and div ---
         call compute_pres_depth (Vgrid, Var%pssl, dp)
         call change_grid (Hgrid, TEMP_GRID, WIND_GRID, dp, udp)
         vdp = Var%v * udp  ! note: using udp to store dp at vel pts
         udp = Var%u * udp
        !if ( id_vor > 0 ) then
             call compute_vorticity (Hgrid, dp, udp, vdp, vor )
             used = send_data ( id_vor, vor(is:ie,hs:he,:), Time,      &
                              mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )
        !endif
        !if ( id_div > 0 ) then
        !    call compute_divergence (Hgrid, dp, udp, vdp, div )
        !    used = send_data ( id_div, div(is:ie,hs:he,:), Time,      &
        !                     mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )
        !endif
      endif

!------- mass fluxes (without topography masks) -----------
      if ( id_mfew > 0 ) then
          used = send_data ( id_mfew, mfew(is:ie,hs:he,:), Time )
      endif
      if ( id_mfns > 0 ) then
          used = send_data ( id_mfns, mfns(is:ie,vs:ve,:), Time )
      endif

!--------------- pressure gradient components --------------------------
!------------------------ NOT USED -------------------------------------

!     if ( id_pgfx > 0 .and. present(pgfx) )  &
!     used = send_data ( id_pgfx, pgfx(is:ie,vs:ve,:), Time,      &
!                        mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

!     if ( id_pgfy > 0 .and. present(pgfy) )  &
!     used = send_data ( id_pgfy, pgfy(is:ie,vs:ve,:), Time,      &
!                        mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

!-----------------------------------------------------------------------

 end subroutine bgrid_diagnostics

!#######################################################################

 subroutine bgrid_diagnostics_tend ( Hgrid, Var_dt, Masks, Time )

!-----------------------------------------------------------------------
!  Hgrid  = horizontal grid constants
!  Var_dt = prognostic variables tendencies FROM ONLY THE DYNAMICS
!  Masks  = grid box masks for step-mountain topography
!  Time   = diagnostics time
!-----------------------------------------------------------------------
type(horiz_grid_type), intent(in) :: Hgrid
type (prog_var_type),  intent(in) :: Var_dt
type(grid_mask_type),  intent(in) :: Masks
type(time_type),       intent(in) :: Time
!-----------------------------------------------------------------------
   integer :: is, ie, hs, he, vs, ve, n
   logical :: used
!-----------------------------------------------------------------------

      ! compute domain indices
      is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie
      hs = Hgrid % Tmp % js;  he = Hgrid % Tmp % je
      vs = Hgrid % Vel % js;  ve = Hgrid % Vel % je

!-----------------------------------------------------------------------
!---------------- 3d prognostic fields ---------------------------

      if ( id_udt > 0 ) &
      used = send_data ( id_udt, Var_dt%u(is:ie,vs:ve,:), Time,    &
                          mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

      if ( id_vdt > 0 ) &
      used = send_data ( id_vdt, Var_dt%v(is:ie,vs:ve,:), Time,   &
                         mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

      if ( id_tdt > 0 ) &
      used = send_data ( id_tdt, Var_dt%t(is:ie,hs:he,:), Time,   &
                         mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )

      do n = 1, Var_dt%ntrace
       if ( id_tracer_tend(n) > 0 ) &
       used = send_data ( id_tracer_tend(n), Var_dt%r(is:ie,hs:he,:,n),  &
                          Time, mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )
      enddo

!-----------------------------------------------------------------------

 end subroutine bgrid_diagnostics_tend

!#######################################################################

 subroutine compute_vorticity ( Hgrid, dp, udp, vdp, vor )

!-----------------------------------------------------------------------
!  Computes relative vorticity on B-grid
!     Hgrid = horizontal grid constants
!     dp    = pressure thickness of model layers at temperature points
!     udp   = zonal wind, u * dp, at velocity points
!     vdp   = meridional wind, v * dp, at velocity points
!     vor   = relative vorticity (1/s) at temperature points
!-----------------------------------------------------------------------
type(horiz_grid_type), intent(in) :: Hgrid
real, intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: dp, udp, vdp
real, intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: vor

real,dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: &
                          vdy, udx, few, fns
integer :: i, j, k, is, ie, js, je

   is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie
   js = Hgrid % Tmp % js;  je = Hgrid % Tmp % je

   do k = 1, size(dp,3)
      do j = js-1, je
         vdy(:,j) = vdp(:,j,k)*Hgrid%Vel%dy
         udx(:,j) = udp(:,j,k)*Hgrid%Vel%dx(j)
      enddo

      do j = js,   je
      do i = is-1, ie
         fns(i,j) = (vdy(i,j-1)+vdy(i,j))*0.5
      enddo
      enddo
      do j = js-1, je
      do i = is,   ie
         few(i,j) = (udx(i-1,j)+udx(i,j))*0.5
      enddo
      enddo

!  ------ vorticity ------
      do j = js, je
      do i = is, ie
         vor(i,j,k)=((fns(i,j)-fns(i-1,j))-(few(i,j)-few(i,j-1))) &
                    /(dp(i,j,k)*Hgrid%Tmp%area(j))
      enddo
      enddo
   enddo

 end subroutine compute_vorticity

!#######################################################################

 subroutine compute_divergence ( Hgrid, dp, udp, vdp, div )

!-----------------------------------------------------------------------
!  Computes divergence on B-grid
!     Hgrid = horizontal grid constants
!     dp    = pressure thickness of model layers at temperature points
!     udp   = zonal wind, u * dp, at velocity points
!     vdp   = meridional wind, v * dp, at velocity points
!     div   = divergence (1/s) at temperature points
!-----------------------------------------------------------------------
type(horiz_grid_type), intent(in) :: Hgrid
real, intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: dp, udp, vdp
real, intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: div

real,dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) ::  &
                          udy, vdx, few, fns
integer :: i, j, k, is, ie, js, je

   is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie
   js = Hgrid % Tmp % js;  je = Hgrid % Tmp % je

   do k = 1, size(dp,3)
      do j = js-1, je
         udy(:,j) = udp(:,j,k)*Hgrid%Vel%dy
         vdx(:,j) = vdp(:,j,k)*Hgrid%Vel%dx(j)
      enddo

      do j = js,   je
      do i = is-1, ie
         few(i,j) = (udy(i,j-1)+udy(i,j))*0.5
      enddo
      enddo
      do j = js-1, je
      do i = is,   ie
         fns(i,j) = (vdx(i-1,j)+vdx(i,j))*0.5
      enddo
      enddo

!  ------ divergence ------
      do j = js, je
      do i = is, ie
         div(i,j,k)=((few(i,j)+fns(i,j))-(few(i-1,j)+fns(i,j-1))) &
                    /(dp(i,j,k)*Hgrid%Tmp%area(j))
      enddo
      enddo
   enddo

 end subroutine compute_divergence

!#######################################################################

 subroutine sea_level_pressure ( psurf, zsurf, pfull, tfull, slp )

 real, intent(in),  dimension(:,:)   :: psurf, zsurf
 real, intent(in),  dimension(:,:,:) :: pfull, tfull
 real, intent(out), dimension(:,:)   :: slp

! psurf = surface pressure
! zsurf = surface geopotential height in meters^2/sec^2
! pfull = pressure at full model levels
! tfull = temperature at full model levels
! slp   = sea level pressure in pascals

 integer :: i, j, k, kr
 real    :: sig, tbot

 real, parameter :: TLAPSE = 6.5e-3
 real, parameter :: GORG = GRAV/(RDGAS*TLAPSE)
 real, parameter :: MRGOG = -1./GORG

     do j = 1, size(psurf,2)
     do i = 1, size(psurf,1)

        if ( abs(zsurf(i,j)) > 0.0001 ) then

            !---- get ref level for temp ----
             do k = 1, size(tfull,3)
                sig = pfull(i,j,k)/psurf(i,j)
                if ( sig > 0.8 ) then
                     kr = k
                     exit
                endif
             enddo

             tbot = tfull(i,j,kr) * sig ** MRGOG
             slp(i,j) =  psurf(i,j) * ( 1.0 + TLAPSE * zsurf(i,j) / (tbot*GRAV) ) ** GORG
        else
             slp(i,j) = psurf(i,j)
        endif
     enddo
     enddo

 end subroutine sea_level_pressure

!#######################################################################

end module bgrid_diagnostics_mod


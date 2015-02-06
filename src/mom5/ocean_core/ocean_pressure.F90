module ocean_pressure_mod
#define COMP isc:iec,jsc:jec
!  
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
! S.M. Griffies 
! </CONTACT>
!
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
! A. Rosati 
! </REVIEWER>
!
!<OVERVIEW>
! Compute the hydrostatic pressure and forces from pressure. 
! Includes methods for either Bgrid or Cgrid. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes hydrostatic pressure and the pressure force
! acting at a velocity point (traditional finite difference approach).
! This force is used for the linear momentum equation.  
!
! This module takes account of the vertical coordinate,
! which determines details of the calculation.  
!
! This module allows for either Bgrid or Cgrid calculation. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
!  S.M. Griffies, 2012: Elements of MOM 
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_pressure_nml">
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging. 
!  </DATA> 
!
!  <DATA NAME="zero_correction_term_grad" TYPE="logical">
!  For debugging it is often useful to zero the contribution to the 
!  pressure gradient that arises from the "correction" term. 
!  Implemented only for depth based vertical coordinate models.
!  </DATA> 
!  <DATA NAME="zero_diagonal_press_grad" TYPE="logical">
!  For debugging it is often useful to zero the contribution to the 
!  pressure gradient that arises from the along k-level gradient.
!  Implemented only for depth based vertical coordinate models.
!  </DATA> 
!  <DATA NAME="zero_pressure_force" TYPE="logical">
!  For debugging it is often useful to zero the pressure force
!  to zero.  
!  </DATA> 
!  <DATA NAME="zero_eta_over_h_zstar_pressure" TYPE="logical">
!  For debugging zstar, we drop any eta/H contribution to 
!  the hydrostatic pressure.  This is wrong physically, but 
!  useful for certain tests.   
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: c2dbars, epsln 
use diag_manager_mod, only: register_diag_field, send_data
use fms_mod,          only: open_namelist_file, check_nml_error, close_file, write_version_number
use mpp_io_mod,       only: mpp_open, mpp_close, MPP_RDONLY, MPP_ASCII
use mpp_mod,          only: input_nml_file, mpp_error, FATAL, stdout, stdlog
use mpp_domains_mod,  only: mpp_update_domains

use ocean_domains_mod,    only: get_local_indices
use ocean_operators_mod,  only: FAY, FAX, FDX_NT, FDY_ET 
use ocean_parameters_mod, only: MOM_BGRID, MOM_CGRID 
use ocean_parameters_mod, only: GEOPOTENTIAL, ZSTAR, ZSIGMA
use ocean_parameters_mod, only: PRESSURE, PSTAR, PSIGMA
use ocean_parameters_mod, only: DEPTH_BASED, PRESSURE_BASED
use ocean_parameters_mod, only: ENERGETIC, FINITEVOLUME
use ocean_parameters_mod, only: missing_value, onehalf, rho0, rho0r, grav
use ocean_types_mod,      only: ocean_domain_type, ocean_grid_type
use ocean_types_mod,      only: ocean_time_type, ocean_velocity_type
use ocean_types_mod,      only: ocean_thickness_type, ocean_density_type
use ocean_types_mod,      only: ocean_lagrangian_type
use ocean_util_mod,       only: write_timestamp, diagnose_3d, diagnose_3d_u, write_chksum_3d
use ocean_workspace_mod,  only: wrk1, wrk2, wrk1_zw, wrk2_zw, wrk1_v, wrk2_v, wrk1_2d
use ocean_workspace_mod,  only: wrk3, wrk4
use ocean_obc_mod,        only: store_ocean_obc_pressure_grad

implicit none

private

public ocean_pressure_init
public pressure_in_dbars
public pressure_in_dbars_blob
public pressure_force
public hydrostatic_pressure
public hydrostatic_pressure_blob
public geopotential_anomaly 
public geopotential_anomaly_blob

private press_grad_force_depth_bgrid
private press_grad_force_depth_cgrid
private press_grad_force_depth_blob

private press_grad_force_press_bgrid
private press_grad_force_press_cgrid
private press_grad_force_press_blob

#include <ocean_memory.h>

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()

character(len=128) :: version = &
     '$Id: ocean_pressure.F90,v 20.0 2013/12/14 00:10:57 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

! for vertical coordinate
integer :: vert_coordinate 
integer :: vert_coordinate_class

! for Bgrid or Cgrid
integer :: horz_grid

! useful constants
real :: grav_rho0r
real :: grav_rho0
real :: p5grav
real :: p5gravr
real :: p5grav_rho0r
real :: dbars2pa

real, allocatable, dimension(:,:,:) :: gridmask

logical :: use_blobs

! for computing geopotential at T-cell bottom from depth_zwt
real :: geopotential_switch=0.0

! for diagnostics
integer :: id_geopotential   =-1
integer :: id_anomgeopot     =-1
integer :: id_anomgeopot_zwt =-1
integer :: id_anompress      =-1
integer :: id_anompress_zwt  =-1

integer :: id_press_force(2)           =-1
integer :: id_pgrad_klev(2)            =-1
integer :: id_rhoprime_geograd_klev(2) =-1
integer :: id_rhoprime_pgrad_klev(2)   =-1
integer :: id_rho_geogradprime_klev(2) =-1

logical :: used

logical :: debug_this_module         =.false. 
logical :: zero_pressure_force       =.false.
logical :: zero_correction_term_grad =.false.
logical :: zero_diagonal_press_grad  =.false.
logical :: module_is_initialized     =.false.
logical :: have_obc                  =.false. 
logical :: zero_eta_over_h_zstar_pressure = .false. 

namelist /ocean_pressure_nml/ debug_this_module, zero_pressure_force, &
         zero_correction_term_grad, zero_diagonal_press_grad,         &
         zero_eta_over_h_zstar_pressure

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_pressure_init">
!
! <DESCRIPTION>
! Initialize the pressure module
! </DESCRIPTION>
!
subroutine ocean_pressure_init(Grid, Domain, Time, ver_coordinate, &
                               ver_coordinate_class, hor_grid, obc, blobs, debug)

  type(ocean_grid_type),   target, intent(in) :: Grid
  type(ocean_domain_type), target, intent(in) :: Domain
  type(ocean_time_type),           intent(in) :: Time
  integer,                         intent(in) :: ver_coordinate
  integer,                         intent(in) :: ver_coordinate_class
  integer,                         intent(in) :: hor_grid
  logical,                         intent(in) :: obc
  logical,                         intent(in) :: blobs
  logical,               optional, intent(in) :: debug

  integer :: ioun, io_status, ierr

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error( FATAL, &
    '==>Error: ocean_pressure_mod (ocean_pressure_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  vert_coordinate       = ver_coordinate
  vert_coordinate_class = ver_coordinate_class
  horz_grid             = hor_grid 
  have_obc              = obc
  use_blobs             = blobs

  grav_rho0    = grav*rho0
  grav_rho0r   = grav*rho0r
  p5grav       = 0.5*grav
  p5gravr      = 0.5/grav
  p5grav_rho0r = 0.5*grav*rho0r 
  dbars2pa     = 1.0/c2dbars

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_pressure_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_pressure_nml')
#else
  ioun = open_namelist_file()
  read  (ioun, ocean_pressure_nml,iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_pressure_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_pressure_nml)
  write (stdlogunit, ocean_pressure_nml)

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain,isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif

  Dom => Domain
  Grd => Grid

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_pressure_mod with debug_this_module=.true.'  
  endif 

  if(zero_pressure_force) then 
    write(stdoutunit,*) &
    '==>Warning: Running MOM with zero horizontal pressure force; meant for debugging.'
  endif 

  if(use_blobs .and. horz_grid == MOM_CGRID) then 
    write(stdoutunit,*) &
    '==>Warning from pressure mod: MOM C-grid with blobs is incomplete; blobs have only been tested with Bgrid.'
  endif 

  if(vert_coordinate==GEOPOTENTIAL) then 
     geopotential_switch=0.0
  else
     geopotential_switch=1.0
  endif


  allocate (gridmask(isd:ied,jsd:jed,nk))
  if(horz_grid == MOM_BGRID) then 
     gridmask(:,:,:) = Grd%umask(:,:,:)
  else
     gridmask(:,:,:) = Grd%tmask(:,:,:)
  endif 


  write(stdoutunit,*) &
       '==>NOTE: Running MOM with finite difference formulation of pressure force.'

  if(zero_correction_term_grad) then 
      write(stdoutunit,*) &
           '==>Warning: Running MOM with zero_correction_term_grad=.true. Unrealistic simulation.'
  endif

  if(zero_diagonal_press_grad) then 
      write(stdoutunit,*) &
           '==>Warning: Running MOM with zero_diagonal_press_grad=.true. Unrealistic simulation.'
  endif


  ! register fields

  id_geopotential = register_diag_field ('ocean_model', 'geopotential', Grd%tracer_axes(1:3),&
   Time%model_time, 'geopotential at T-point', 'm^2/s^2',                                    &
   missing_value=missing_value, range=(/-1e8,1e8/))

  id_anompress = register_diag_field ('ocean_model', 'anompress', Grd%tracer_axes(1:3),&
   Time%model_time, 'anomalous hydrostatic pressure at T-point', 'Pa',                 &
   missing_value=missing_value, range=(/-1e8,1e8/))

  id_anompress_zwt = register_diag_field ('ocean_model', 'anompress_zwt', Grd%tracer_axes(1:3),&
   Time%model_time, 'anomalous hydrostatic pressure at T-cell bottom', 'Pa',                   &
   missing_value=missing_value, range=(/-1e8,1e8/))

  id_anomgeopot = register_diag_field ('ocean_model', 'anomgeopot', Grd%tracer_axes(1:3),&
   Time%model_time, 'anomalous geopotential at T-point', 'm^2/s^2',                      &
   missing_value=missing_value, range=(/-1e8,1e8/))

  id_anomgeopot_zwt = register_diag_field ('ocean_model', 'anomgeopot_zwt', Grd%tracer_axes(1:3),&
   Time%model_time, 'anomalous geopotential at T-cell bottom', 'm^2/s^2',                        &
   missing_value=missing_value, range=(/-1e8,1e8/))

  id_pgrad_klev(1) = register_diag_field ('ocean_model', 'dzupgrad_x_klev', Grd%vel_axes_u(1:3), &
     Time%model_time, 'dz*dp_dx on k-level', 'Pa', missing_value=missing_value, range=(/-1e9,1e9/))
  id_pgrad_klev(2) = register_diag_field ('ocean_model', 'dzupgrad_y_klev', Grd%vel_axes_v(1:3), &
     Time%model_time, 'dz*dp_dy on k-level', 'Pa', missing_value=missing_value, range=(/-1e9,1e9/))

  id_rhoprime_geograd_klev(1) = register_diag_field ('ocean_model', 'dzurhoprime_geograd_x_klev',&
     Grd%vel_axes_u(1:3), Time%model_time, 'dz*rhoprime*d(geopot)/dx on k-level', 'Pa',          &
     missing_value=missing_value, range=(/-1e9,1e9/))
  id_rhoprime_geograd_klev(2) = register_diag_field ('ocean_model', 'dzurhoprime_geograd_y_klev',&
    Grd%vel_axes_v(1:3), Time%model_time, 'dz*rhoprime*d(geopot)/dx on k-level', 'Pa',           &
    missing_value=missing_value, range=(/-1e9,1e9/))

  id_rhoprime_pgrad_klev(1) = register_diag_field ('ocean_model', 'dzurhoprime_pgrad_x_klev',&
     Grd%vel_axes_u(1:3), Time%model_time, 'dz*(rhoprime/rho0)*dp/dx on k-level', 'Pa',      &
     missing_value=missing_value, range=(/-1e9,1e9/))
  id_rhoprime_pgrad_klev(2) = register_diag_field ('ocean_model', 'dzurhoprime_pgrad_y_klev',&
     Grd%vel_axes_v(1:3), Time%model_time, 'dz*(rhoprime/rho0)*dp/dy on k-level', 'Pa',      &
     missing_value=missing_value, range=(/-1e9,1e9/))

  id_rho_geogradprime_klev(1) = register_diag_field ('ocean_model', 'dzurho_geogradprime_x_klev',&
    Grd%vel_axes_u(1:3), Time%model_time, 'dzu*rho*d(geopot_prime)/dx on k-level', 'Pa',         &
    missing_value=missing_value, range=(/-1e9,1e9/))
  id_rho_geogradprime_klev(2) = register_diag_field ('ocean_model', 'dzurho_geogradprime_y_klev',&
    Grd%vel_axes_v(1:3), Time%model_time, 'dzu*rho*d(geopot_prime)/dy on k-level', 'Pa/m',       &
    missing_value=missing_value, range=(/-1e9,1e9/))

  id_press_force(1) = register_diag_field ('ocean_model', 'press_force_u',&
     Grd%vel_axes_u(1:3), Time%model_time, 'i-directed pressure force',   &
     'Pa', missing_value=missing_value, range=(/-1e9,1e9/))
  id_press_force(2) = register_diag_field ('ocean_model', 'press_force_v',&
     Grd%vel_axes_v(1:3), Time%model_time, 'j-directed pressure force',   &
     'Pa', missing_value=missing_value, range=(/-1e9,1e9/))
 

end subroutine ocean_pressure_init
! </SUBROUTINE> NAME="ocean_pressure_init"


!#######################################################################
! <SUBROUTINE NAME="pressure_force">
!
! <DESCRIPTION>
! Compute the horizontal force [Pa=N/m^2] from pressure. 
!
! Use the traditional approach whereby the pressure force
! is computed as a finite difference gradient centred 
! at the U-cell point. 
!
! </DESCRIPTION>
!
subroutine pressure_force(Time, Thickness, Dens, Velocity, L_system, rho)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_velocity_type),      intent(inout) :: Velocity
  type(ocean_lagrangian_type),    intent(inout) :: L_system
  real, dimension(isd:,jsd:,:),   intent(in)    :: rho

  integer :: tau, taup1
 
  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau 
  taup1 = Time%taup1

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_pressure_mod (pressure_force): module must be initialized')
  endif 

  if(horz_grid == MOM_BGRID) then 

     if(use_blobs) then
        if(vert_coordinate_class==DEPTH_BASED) then 
           call press_grad_force_depth_blob(Time, Thickness, Velocity, L_system, rho)
        else 
           call press_grad_force_press_blob(Time, Thickness, Dens, Velocity, L_system, rho)
        endif
     else
        if(vert_coordinate_class==DEPTH_BASED) then 
           call press_grad_force_depth_bgrid(Time, Thickness, Velocity, rho)
        else 
           call press_grad_force_press_bgrid(Time, Thickness, Dens, Velocity, rho)
        endif
     endif

  else ! MOM_CGRID (smg: blobs not yet coded for Cgrid)

     if(vert_coordinate_class==DEPTH_BASED) then 
        call press_grad_force_depth_cgrid(Time, Thickness, Velocity, rho)
     else
        call press_grad_force_press_cgrid(Time, Thickness, Dens, Velocity, rho)
     endif

  endif 

  ! for debugging 
  if(zero_pressure_force) then 
      Velocity%press_force = 0.0
  endif

  ! for open boundaries 
  if (have_obc) call store_ocean_obc_pressure_grad(Thickness, Velocity%press_force) 

  ! send some diagnostics 

  if(id_geopotential > 0) call diagnose_3d(Time, Grd, id_geopotential, -grav*Thickness%geodepth_zt(:,:,:))

  if (id_press_force(1) > 0) then 
       used = send_data( id_press_force(1), Velocity%press_force(:,:,:,1), &
       Time%model_time, rmask=gridmask(:,:,:),                             &
       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_press_force(2) > 0) then 
       used = send_data( id_press_force(2), Velocity%press_force(:,:,:,2), &
       Time%model_time, rmask=gridmask(:,:,:),                             &
       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 


  ! more debugging 
  if(debug_this_module) then 
      write(stdoutunit,*)' '
      write(stdoutunit,*) 'From ocean_pressure_mod: pressure chksums'
      call write_timestamp(Time%model_time)
      call write_chksum_3d('rho_dzu(tau)', Thickness%rho_dzu(COMP,:,tau)*Grd%umask(COMP,:))
      call write_chksum_3d('press_force(1)', Velocity%press_force(COMP,:,1)*Grd%umask(COMP,:))
      call write_chksum_3d('press_force(2)', Velocity%press_force(COMP,:,2)*Grd%umask(COMP,:))
  endif 


end subroutine pressure_force
! </SUBROUTINE> NAME="pressure_force"



!#######################################################################
! <SUBROUTINE NAME="press_grad_force_depth_bgrid">
!
! <DESCRIPTION>
! Compute the force from pressure using a finite difference method
! to compute the thickness weighted pressure gradient at the 
! velocity cell point. 
!
! Assume B-grid arrangement here.
!
! For depth-like vertical coordinates, we exclude surface and applied 
! pressures (i.e., we are computing here the gradient of the baroclinic 
! pressure).  The surface and applied pressures are accounted for in 
! the barotropic module. 
!
! Account is taken of variable partial cell thickness.
!
! 1 = dp/dx; 2 = dp/dy
!
! Thickness weight since this is what we wish to use in update of 
! the velocity. Resulting thickness weighted pressure gradient has 
! dimensions of Pa = N/m^2 = kg/(m*s^2).
!
! Thickness%dzu should be at tau.
!
! </DESCRIPTION>
!
subroutine press_grad_force_depth_bgrid(Time, Thickness, Velocity, rho)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_velocity_type),      intent(inout) :: Velocity
  real, dimension(isd:,jsd:,:),   intent(in)    :: rho

  real, dimension(isd:ied,jsd:jed) :: diff_geo_x
  real, dimension(isd:ied,jsd:jed) :: diff_geo_y
  real, dimension(isd:ied,jsd:jed) :: tmp1
  real, dimension(isd:ied,jsd:jed) :: tmp2
  integer :: i, j, k
  integer :: tau

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_pressure_mod (press_grad_force_depth_bgrid): module must be initialized')
  endif 

  tau    = Time%tau
  tmp1   = 0.0
  tmp2   = 0.0
  wrk1_v = 0.0
  wrk2_v = 0.0


  ! use density anomaly to improve accuracy 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk2(i,j,k) = Grd%tmask(i,j,k)*(rho(i,j,k)-rho0)
        enddo
     enddo
  enddo
  
  ! compute hydrostatic pressure anomaly at T-point 
  ! ( Pa=N/m^2=kg/(m*s^2) )
  wrk1(:,:,:) = hydrostatic_pressure(Thickness, wrk2(:,:,:))  

  diff_geo_x(:,:) = 0.0
  diff_geo_y(:,:) = 0.0

  do k=1,nk

     ! geopotential = -grav*geodepth, hence the minus sign below 
     do j=jsd,jed
        do i=isd,iec
           diff_geo_x(i,j) = -grav*(Thickness%geodepth_zt(i+1,j,k)-Thickness%geodepth_zt(i,j,k))
        enddo
     enddo
     do j=jsd,jec
        do i=isd,ied
           diff_geo_y(i,j) = -grav*(Thickness%geodepth_zt(i,j+1,k)-Thickness%geodepth_zt(i,j,k))
        enddo
     enddo

     ! density anomaly times geopotential gradient on k-levels 
     ! (dzurhoprime_geograd_x_klev,dzurhoprime_geograd_y_klev)
     tmp1(:,:) = FAY(FAX(wrk2(:,:,k))*diff_geo_x(:,:))
     tmp2(:,:) = FAX(FAY(wrk2(:,:,k))*diff_geo_y(:,:))
     do j=jsc,jec
        do i=isc,iec
           wrk1_v(i,j,k,1) = tmp1(i,j)*Grd%dxur(i,j)*Thickness%dzu(i,j,k)
           wrk1_v(i,j,k,2) = tmp2(i,j)*Grd%dyur(i,j)*Thickness%dzu(i,j,k)
        enddo
     enddo
     if(zero_correction_term_grad) then 
         wrk1_v(:,:,:,:) = 0.0
     endif

     ! gradient of anomalous baroclinic pressure on k-levels 
     ! (dzupgrad_x_klev,dzupgrad_y_klev)
     tmp1(:,:) = FDX_NT(FAY(wrk1(:,:,k)))
     tmp2(:,:) = FDY_ET(FAX(wrk1(:,:,k)))
     do j=jsc,jec
        do i=isc,iec
           wrk2_v(i,j,k,1) = tmp1(i,j)*Thickness%dzu(i,j,k)
           wrk2_v(i,j,k,2) = tmp2(i,j)*Thickness%dzu(i,j,k)
        enddo
     enddo

     if(zero_diagonal_press_grad) then 
         wrk2_v(:,:,:,:) = 0.0
     endif

     ! thickness weighted baroclinic pressure gradient on z-levels
     ! (dzu_pgrad_u,dzu_pgrad_v)
     do j=jsc,jec
        do i=isc,iec    
           Velocity%press_force(i,j,k,1) = -Grd%umask(i,j,k)*( wrk1_v(i,j,k,1) + wrk2_v(i,j,k,1) )
           Velocity%press_force(i,j,k,2) = -Grd%umask(i,j,k)*( wrk1_v(i,j,k,2) + wrk2_v(i,j,k,2) )
        enddo
     enddo

  enddo   ! k do-loop finish 


  ! diagnostics 

  call diagnose_3d(Time, Grd, id_anompress, wrk1(:,:,:))

  ! (dzurhoprime_geograd_x_klev,dzurhoprime_geograd_y_klev)
  call diagnose_3d_u(Time, Grd, id_rhoprime_geograd_klev(1),wrk1_v(:,:,:,1))
  call diagnose_3d_u(Time, Grd, id_rhoprime_geograd_klev(2),wrk1_v(:,:,:,2))

  ! (dzupgrad_x_klev,dzupgrad_y_klev)
  call diagnose_3d_u(Time, Grd, id_pgrad_klev(1), wrk2_v(:,:,:,1))
  call diagnose_3d_u(TIme, Grd, id_pgrad_klev(2), wrk2_v(:,:,:,2))

end subroutine press_grad_force_depth_bgrid
! </SUBROUTINE> NAME="press_grad_force_depth_bgrid"



!#######################################################################
! <SUBROUTINE NAME="press_grad_force_press_bgrid">
!
! <DESCRIPTION>
! Compute the force from pressure using a finite difference method
! to compute the thickness weighted pressure gradient at the 
! velocity cell corner point. 
!
! Assume B-grid arrangement here.
!
! For pressure-like vertical coordinates, we omit the bottom pressure
! and bottom geopotential. These pressures are accounted for in the 
! barotropic module.   
!
! Account is taken of variable partial cell thickness. 
!
! 1 = dp/dx; 2 = dp/dy
!
! Thickness weight since this is what we wish to use in update of 
! the velocity. Resulting thickness weighted pressure gradient has 
! dimensions of Pa = N/m^2 = kg/(m*s^2).
!
! Thickness%dzu should be at tau.
!
! </DESCRIPTION>
!
subroutine press_grad_force_press_bgrid(Time, Thickness, Dens, Velocity, rho)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_velocity_type),      intent(inout) :: Velocity
  real, dimension(isd:,jsd:,:),   intent(in)    :: rho

  real, dimension(isd:ied,jsd:jed) :: diff_press_x
  real, dimension(isd:ied,jsd:jed) :: diff_press_y
  real, dimension(isd:ied,jsd:jed) :: tmp1
  real, dimension(isd:ied,jsd:jed) :: tmp2
  integer :: i, j, k
  integer :: tau

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_pressure_mod (press_grad_force_press_bgrid): module must be initialized')
  endif 

  tau    = Time%tau
  tmp1   = 0.0
  tmp2   = 0.0
  wrk1_v = 0.0
  wrk2_v = 0.0

  ! use density anomaly to improve accuracy 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk2(i,j,k) = Grd%tmask(i,j,k)*(rho(i,j,k)-rho0)
        enddo
     enddo
  enddo

  ! compute geopotential anomaly (m^2/s^2) at T-cell point 
  wrk1(:,:,:) = geopotential_anomaly(Thickness, wrk2(:,:,:)) 

  diff_press_x(:,:) = 0.0
  diff_press_y(:,:) = 0.0

  do k=1,nk

     do j=jsd,jed
        do i=isd,iec
           diff_press_x(i,j) = -dbars2pa*(Dens%pressure_at_depth(i+1,j,k)-Dens%pressure_at_depth(i,j,k))
        enddo
     enddo
     do j=jsd,jec
        do i=isd,ied
           diff_press_y(i,j) = -dbars2pa*(Dens%pressure_at_depth(i,j+1,k)-Dens%pressure_at_depth(i,j,k))
        enddo
     enddo

     ! anomalous density times gradient of pressure on k-level               
     ! (dzurhoprime_pgrad_x_klev,dzurhoprime_pgrad_y_klev)
     tmp1(:,:) = rho0r*FAY(FAX(wrk2(:,:,k))*diff_press_x(:,:))
     tmp2(:,:) = rho0r*FAX(FAY(wrk2(:,:,k))*diff_press_y(:,:))
     do j=jsc,jec
        do i=isc,iec
           wrk1_v(i,j,k,1) = tmp1(i,j)*Grd%dxur(i,j)*Thickness%dzu(i,j,k)
           wrk1_v(i,j,k,2) = tmp2(i,j)*Grd%dyur(i,j)*Thickness%dzu(i,j,k)
        enddo
     enddo

     ! density times gradient of anomalous geopotential on k-level
     ! (dzurho_geogradprime_x_klev,dzurho_geogradprime_y_klev) 
     tmp1(:,:) = FDX_NT(FAY(wrk1(:,:,k)))
     tmp2(:,:) = FDY_ET(FAX(wrk1(:,:,k)))
     do j=jsc,jec
        do i=isc,iec
           wrk2_v(i,j,k,1) = tmp1(i,j)*Thickness%rho_dzu(i,j,k,tau)
           wrk2_v(i,j,k,2) = tmp2(i,j)*Thickness%rho_dzu(i,j,k,tau)
        enddo
     enddo

     ! thickness weighted baroclinic pressure gradient on z-levels
     ! (dzu_pgrad_u,dzu_pgrad_v)
     do j=jsc,jec
        do i=isc,iec    
           Velocity%press_force(i,j,k,1) = -Grd%umask(i,j,k)*( wrk1_v(i,j,k,1) + wrk2_v(i,j,k,1) )
           Velocity%press_force(i,j,k,2) = -Grd%umask(i,j,k)*( wrk1_v(i,j,k,2) + wrk2_v(i,j,k,2) )
        enddo
     enddo

  enddo  ! k do-loop finish 


  ! diagnostics 

  call diagnose_3d(Time, Grd, id_anomgeopot, wrk1(:,:,:))

  ! (dzurhoprime_pgrad_x_klev,dzurhoprime_pgrad_y_klev)
  call diagnose_3d_u(Time, Grd, id_rhoprime_pgrad_klev(1),wrk1_v(:,:,:,1))
  call diagnose_3d_u(Time, Grd, id_rhoprime_pgrad_klev(2),wrk1_v(:,:,:,2))

  ! (dzurho_geogradprime_x_klev,dzurho_geogradprime_y_klev)
  call diagnose_3d_u(Time, Grd, id_rho_geogradprime_klev(1),wrk2_v(:,:,:,1))
  call diagnose_3d_u(Time, Grd, id_rho_geogradprime_klev(2),wrk2_v(:,:,:,2))

end subroutine press_grad_force_press_bgrid
! </SUBROUTINE> NAME="press_grad_force_press_bgrid"


!#######################################################################
! <SUBROUTINE NAME="press_grad_force_depth_cgrid">
!
! <DESCRIPTION>
! Compute the force from pressure using a finite difference method
! to compute the thickness weighted pressure gradient at the 
! T-cell face, thus acting on the C-grid velocity components.  
!
! Assume C-grid arrangement here.
!
! For depth-like vertical coordinates, we exclude surface and applied 
! pressures (i.e., we are computing here the gradient of the baroclinic 
! pressure).  The surface and applied pressures are included in 
! ocean_barotropic.
!
! Account is taken of variable partial cell thickness.
!
! 1 = dp/dx; 2 = dp/dy
!
! Thickness weight since this is what we wish to use in update of 
! the velocity. Resulting thickness weighted pressure gradient has 
! dimensions of Pa = N/m^2 = kg/(m*s^2).
!
! Thickness%dzten should be at time tau.
!
! </DESCRIPTION>
!
subroutine press_grad_force_depth_cgrid(Time, Thickness, Velocity, rho)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_velocity_type),    intent(inout) :: Velocity
  real, dimension(isd:,jsd:,:), intent(in)    :: rho

  real    :: rhoave 
  integer :: i, j, k
  integer :: tau

  tau    = Time%tau
  wrk1   = 0.0
  wrk2   = 0.0
  wrk1_v = 0.0
  wrk2_v = 0.0

  ! density anomaly to improve accuracy 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk2(i,j,k) = Grd%tmask(i,j,k)*(rho(i,j,k)-rho0)
        enddo
     enddo
  enddo
  
  ! hydrostatic pressure anomaly at T-point ( Pa = N/m^2 = kg/(m*s^2) )
  wrk1(:,:,:) = hydrostatic_pressure(Thickness, wrk2(:,:,:))  

  ! wrk1_v = horizontal gradient of geopotential on T-cell faces (geopotential = -grav*geodepth);
  ! wrk2_v = horizontal gradient of anomalous baroclinic pressure on T-cell faces;
  ! press_force = thickness weighted baroclinic pressure gradient on z-levels at T-cell faces. 
  do k=1,nk
     do j=jsd,jec
        do i=isd,iec
           rhoave          = onehalf*(wrk2(i+1,j,k)+wrk2(i,j,k))
           wrk1_v(i,j,k,1) = -grav*rhoave*(Thickness%geodepth_zt(i+1,j,k)-Thickness%geodepth_zt(i,j,k))
           rhoave          = onehalf*(wrk2(i,j+1,k)+wrk2(i,j,k))
           wrk1_v(i,j,k,2) = -grav*rhoave*(Thickness%geodepth_zt(i,j+1,k)-Thickness%geodepth_zt(i,j,k))
           wrk2_v(i,j,k,1) = wrk1(i+1,j,k)-wrk1(i,j,k)
           wrk2_v(i,j,k,2) = wrk1(i,j+1,k)-wrk1(i,j,k)
        enddo
     enddo
  enddo
  if(zero_correction_term_grad) then 
     wrk1_v(:,:,:,:) = 0.0
  endif
  if(zero_diagonal_press_grad) then 
      wrk2_v(:,:,:,:) = 0.0
  endif
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           Velocity%press_force(i,j,k,1) = -Grd%tmasken(i,j,k,1)*Thickness%dzten(i,j,k,1)*Grd%dxter(i,j) * ( wrk1_v(i,j,k,1) + wrk2_v(i,j,k,1) )
           Velocity%press_force(i,j,k,2) = -Grd%tmasken(i,j,k,2)*Thickness%dzten(i,j,k,2)*Grd%dytnr(i,j) * ( wrk1_v(i,j,k,2) + wrk2_v(i,j,k,2) )
        enddo
     enddo
  enddo

  ! diagnostics 
  call diagnose_3d(Time, Grd, id_anompress, wrk1(:,:,:))

  if (id_rhoprime_geograd_klev(1) > 0) then 
     do k=1,nk
        do j=jsd,jec
           do i=isd,iec
              wrk1_v(i,j,k,1) = wrk1_v(i,j,k,1)*Thickness%dzten(i,j,k,1)*Grd%dxter(i,j)
           enddo
        enddo
     enddo
     call diagnose_3d(Time, Grd, id_rhoprime_geograd_klev(1), wrk1_v(:,:,:,1))
  endif
  if (id_rhoprime_geograd_klev(2) > 0) then 
     do k=1,nk
        do j=jsd,jec
           do i=isd,iec
              wrk1_v(i,j,k,2) = wrk1_v(i,j,k,2)*Thickness%dzten(i,j,k,2)*Grd%dytnr(i,j)
           enddo
        enddo
     enddo
     call diagnose_3d(Time, Grd, id_rhoprime_geograd_klev(2), wrk1_v(:,:,:,2))
  endif
  if (id_pgrad_klev(1) > 0) then 
     do k=1,nk
        do j=jsd,jec
           do i=isd,iec
              wrk2_v(i,j,k,1) = wrk2_v(i,j,k,1)*Thickness%dzten(i,j,k,1)*Grd%dxter(i,j)
           enddo
        enddo
     enddo
     call diagnose_3d(Time, Grd, id_pgrad_klev(1), wrk2_v(:,:,:,1))
  endif
  if (id_pgrad_klev(2) > 0) then 
     do k=1,nk
        do j=jsd,jec
           do i=isd,iec
              wrk2_v(i,j,k,2) = wrk2_v(i,j,k,2)*Thickness%dzten(i,j,k,2)*Grd%dytnr(i,j)
           enddo
        enddo
     enddo
     call diagnose_3d(Time, Grd, id_pgrad_klev(2), wrk2_v(:,:,:,2))
  endif


end subroutine press_grad_force_depth_cgrid
! </SUBROUTINE> NAME="press_grad_force_depth_cgrid"


!#######################################################################
! <SUBROUTINE NAME="press_grad_force_press_cgrid">
!
! <DESCRIPTION>
! Compute the force from pressure using a finite difference method
! to compute the thickness weighted pressure gradient at the
! T-cell face, thus acting on the C-grid velocity components.
!
! Assume C-grid arrangement here.
!
! For pressure-like vertical coordinates, we omit the bottom pressure
! and bottom geopotential; these are handled in ocean_barotropic.  
!
! Account is taken of variable partial cell thickness.  
!
! 1 = dp/dx; 2 = dp/dy
!
! Thickness weight since this is what we wish to use in update of 
! the velocity. Resulting thickness weighted pressure gradient has 
! dimensions of Pa = N/m^2 = kg/(m*s^2).
!
! Thickness%dzu should be at tau.
!
! </DESCRIPTION>
!
subroutine press_grad_force_press_cgrid(Time, Thickness, Dens, Velocity, rho)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_velocity_type),    intent(inout) :: Velocity
  real, dimension(isd:,jsd:,:), intent(in)    :: rho

  integer :: i, j, k
  integer :: tau
  real    :: rhoave 

  tau    = Time%tau
  wrk1   = 0.0
  wrk2   = 0.0
  wrk1_v = 0.0
  wrk2_v = 0.0

  ! use density anomaly to improve accuracy 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk2(i,j,k) = Grd%tmask(i,j,k)*(rho(i,j,k)-rho0)
        enddo
     enddo
  enddo

  ! geopotential anomaly (m^2/s^2) at T-cell point 
  wrk1(:,:,:) = geopotential_anomaly(Thickness, wrk2(:,:,:)) 

  ! wrk1_v = rho0r*rhoave * dz * horizontal gradient of pressure on T-cell faces;
  ! wrk2_v = rho_dz * horizontal gradient of anomalous geopotential on T-cell faces;
  ! press_force = thickness weighted baroclinic pressure gradient on z-levels at T-cell faces. 
  do k=1,nk
     do j=jsd,jec
        do i=isd,iec
           rhoave          = onehalf*(wrk2(i+1,j,k)+wrk2(i,j,k))
           wrk1_v(i,j,k,1) = -rhoave*rho0r*Thickness%dzten(i,j,k,1)*dbars2pa &
                             *(Dens%pressure_at_depth(i+1,j,k)-Dens%pressure_at_depth(i,j,k))
           rhoave          = onehalf*(wrk2(i,j+1,k)+wrk2(i,j,k))
           wrk1_v(i,j,k,2) = -rhoave*rho0r*Thickness%dzten(i,j,k,2)*dbars2pa &
                             *(Dens%pressure_at_depth(i,j+1,k)-Dens%pressure_at_depth(i,j,k))
           wrk2_v(i,j,k,1) = Thickness%rho_dzten(i,j,k,1)*(wrk1(i+1,j,k)-wrk1(i,j,k))
           wrk2_v(i,j,k,2) = Thickness%rho_dzten(i,j,k,2)*(wrk1(i,j+1,k)-wrk1(i,j,k))
           Velocity%press_force(i,j,k,1) = -Grd%dxter(i,j) * ( wrk1_v(i,j,k,1) + wrk2_v(i,j,k,1) )
           Velocity%press_force(i,j,k,2) = -Grd%dytnr(i,j) * ( wrk1_v(i,j,k,2) + wrk2_v(i,j,k,2) )
        enddo
     enddo
  enddo

  ! diagnostics 
  call diagnose_3d(Time, Grd, id_anomgeopot, wrk1(:,:,:))

  call diagnose_3d(Time, Grd, id_rhoprime_pgrad_klev(1), wrk1_v(:,:,:,1))
  call diagnose_3d(Time, Grd, id_rhoprime_pgrad_klev(2), wrk1_v(:,:,:,2))

  call diagnose_3d(Time, Grd, id_rho_geogradprime_klev(1), wrk2_v(:,:,:,1))
  call diagnose_3d(Time, Grd, id_rho_geogradprime_klev(2), wrk2_v(:,:,:,2))


end subroutine press_grad_force_press_cgrid
! </SUBROUTINE> NAME="press_grad_force_press_cgrid"


!#######################################################################
! <FUNCTION NAME="pressure_in_dbars">
!
! <DESCRIPTION>
! Compute pressure (dbars) exerted at T cell grid point by weight of
! water column above the grid point. 
!
! rho = density in kg/m^3
!
! psurf = surface pressure in Pa= kg/(m*s^2) = hydrostatic pressure 
! at z=0 associated with fluid between z=0 and z=eta_t.
! Also include pressure from atmosphere and ice, both of which 
! are part of the patm array. 
!
! This routine is used by ocean_density to compute the pressure 
! used in the equation of state.  It is only called when the 
! vertical coordinate is DEPTH_BASED.
!
! </DESCRIPTION>
function pressure_in_dbars (Thickness, rho, psurf)

  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: rho   
  real, dimension(isd:,jsd:),   intent(in) :: psurf 
  real, dimension(isd:ied,jsd:jed,nk)      :: pressure_in_dbars
  integer :: k

  if ( .not. module_is_initialized ) then   
    call mpp_error(FATAL, &
    '==>Error in ocean_pressure_mod (pressure_in_dbars): module must be initialized')
  endif 

  pressure_in_dbars(:,:,:) = hydrostatic_pressure (Thickness, rho(:,:,:))*c2dbars

  ! add pressure from free surface height and loading from atmosphere and/or sea ice
  do k=1,nk
     pressure_in_dbars(:,:,k) = pressure_in_dbars(:,:,k) + psurf(:,:)*c2dbars 
  enddo


end function pressure_in_dbars
! </FUNCTION> NAME="pressure_in_dbars"


!#######################################################################
! <FUNCTION NAME="hydrostatic_pressure">
!
! <DESCRIPTION>
! Hydrostatic pressure [Pa=N/m^2=kg/(m*s^2)] at T cell grid points.
!
! For GEOPOTENTIAL vertical coordinate, integration is 
! from z=0 to depth of grid point.  This integration results in  
! the so-called "baroclinic" pressure. 
!
! For ZSTAR or ZSIGMA, vertical coordinate, integration is from z=eta to 
! depth of grid point.  This is allowed because ZSTAR and ZSIGMA 
! absorb the undulations of the surface height into their definition.
!
! If the input density "rho" is an anomoly, the resulting presure 
! will be a hydrostatic pressure anomoly. If "rho" is full density, 
! the presure will be a full hydrostatic pressure.
!
! </DESCRIPTION>
!
function hydrostatic_pressure (Thickness, rho)

  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: rho

  integer :: i, j, k
  real, dimension(isd:ied,jsd:jed,nk) :: hydrostatic_pressure

  if(vert_coordinate==GEOPOTENTIAL) then 
      do j=jsd,jed
         do i=isd,ied
            hydrostatic_pressure(i,j,1) = rho(i,j,1)*grav*Grd%dzw(0)
         enddo
      enddo
  elseif(vert_coordinate==ZSTAR .or. vert_coordinate==ZSIGMA) then 
      do j=jsd,jed
         do i=isd,ied
            hydrostatic_pressure(i,j,1) = rho(i,j,1)*grav*Thickness%dzwt(i,j,0)
         enddo
      enddo
  endif

  if(Thickness%method==ENERGETIC) then 

      do k=2,nk
         do j=jsd,jed
            do i=isd,ied
               hydrostatic_pressure(i,j,k) = hydrostatic_pressure(i,j,k-1) &
               + Grd%tmask(i,j,k)*p5grav                                   &
                *(rho(i,j,k-1)+rho(i,j,k))*Thickness%dzwt(i,j,k-1) 
            enddo
         enddo
      enddo

  elseif(Thickness%method==FINITEVOLUME) then 

      do k=2,nk
         do j=jsd,jed
            do i=isd,ied
               hydrostatic_pressure(i,j,k) = hydrostatic_pressure(i,j,k-1) &
               + Grd%tmask(i,j,k)*grav                                     &
                *(Thickness%dztlo(i,j,k-1)*rho(i,j,k-1)+Thickness%dztup(i,j,k)*rho(i,j,k)) 
            enddo
         enddo
      enddo

  endif

  ! For some debugging...we here drop the contribution from a 
  ! nonzero eta in the calculation of hydrostatic_pressure. 
  ! This is wrong physically, but it is useful for certain tests
  ! of the zstar vertical coordinate.  
  if(vert_coordinate==ZSTAR .and. zero_eta_over_h_zstar_pressure) then 
      do j=jsd,jed
         do i=isd,ied
            hydrostatic_pressure(i,j,1) = rho(i,j,1)*grav*Thickness%dswt(i,j,0)
         enddo
      enddo
      do k=2,nk
         do j=jsd,jed
            do i=isd,ied
               hydrostatic_pressure(i,j,k) = hydrostatic_pressure(i,j,k-1) &
               + Grd%tmask(i,j,k)*p5grav                                   &
                *(rho(i,j,k-1)+rho(i,j,k))*Thickness%dswt(i,j,k-1) 
            enddo
         enddo
      enddo
  endif 


end function hydrostatic_pressure
! </FUNCTION> NAME="hydrostatic_pressure"



!#######################################################################
! <FUNCTION NAME="geopotential_anomaly">
!
! <DESCRIPTION>
! Geopotential anomaly [m^2/s^2] at T cell grid points.
! Integration here is from z=-H to depth of grid point.  
!
! Input should be density anomaly rhoprime = rho-rho0.
!
! This function is needed when computing pressure gradient
! for PRESSURE_BASED vertical coordinates. 
!
! WARNING: Thickness%method==FINITEVOLUME has been found to be 
! problematic.  It remains under development.  It is NOT 
! supported for general use. 
! 
! </DESCRIPTION>
!
function geopotential_anomaly (Thickness, rhoprime)

  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: rhoprime

  integer :: i, j, k, kbot
  real, dimension(isd:ied,jsd:jed,nk) :: geopotential_anomaly 
  
  geopotential_anomaly(:,:,:) = 0.0

  if(Thickness%method==ENERGETIC) then 
      do j=jsd,jed
         do i=isd,ied
            kbot=Grd%kmt(i,j)
            if(kbot > 1) then  
                k=kbot
                geopotential_anomaly(i,j,k) = -grav_rho0r*rhoprime(i,j,k)*Thickness%dzwt(i,j,k) 
                do k=kbot-1,1,-1
                   geopotential_anomaly(i,j,k) = geopotential_anomaly(i,j,k+1) &
                   -p5grav_rho0r*(rhoprime(i,j,k+1)+rhoprime(i,j,k))*Thickness%dzwt(i,j,k)
                enddo
            endif
         enddo
      enddo
  elseif(Thickness%method==FINITEVOLUME) then 
      do j=jsd,jed
         do i=isd,ied
            kbot=Grd%kmt(i,j)
            if(kbot > 1) then
                k=kbot
                geopotential_anomaly(i,j,k) = -grav_rho0r*rhoprime(i,j,k)*Thickness%dztlo(i,j,k)   
                do k=kbot-1,1,-1
                   geopotential_anomaly(i,j,k) = geopotential_anomaly(i,j,k+1) &
                   -grav_rho0r*(Thickness%dztlo(i,j,k)*rhoprime(i,j,k)+Thickness%dztup(i,j,k+1)*rhoprime(i,j,k+1)) 
                enddo
            endif
         enddo
      enddo
  endif

end function geopotential_anomaly
! </FUNCTION> NAME="geopotential_anomaly"


!#######################################################################
! <SUBROUTINE NAME="press_grad_force_depth_blob">
!
! <DESCRIPTION>
! This routine respects the partition between Eulerian system mass
! and the Lagrangian system mass associated with the Lagrangian blobs
! model.  The pressure gradient from the total (combined Eulerian and 
! Lagrangian systems) is calculated.  Aside from that, the routine is 
! the same as for press_grad_force_depth.
!
! Compute the force from pressure using a finite difference method
! to compute the thickness weighted pressure gradient at the 
! velocity cell point. 
!
! For depth-like vertical coordinates, we exclude surface and applied 
! pressures (i.e., we are computing here the gradient of the baroclinic 
! pressure).  Account is taken of variable partial cell thickness.
! 1 = dp/dx; 2 = dp/dy
!
! Thickness weight since this is what we wish to use in update of 
! the velocity. Resulting thickness weighted pressure gradient has 
! dimensions of Pa = N/m^2 = kg/(m*s^2).
!
! Thickness%dzu should be at tau.
! </DESCRIPTION>
!
subroutine press_grad_force_depth_blob(Time, Thickness, Velocity, L_system, rho)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_velocity_type),      intent(inout) :: Velocity
  type(ocean_lagrangian_type),    intent(inout) :: L_system
  real, dimension(isd:,jsd:,:),   intent(in)    :: rho

  real, dimension(isd:ied,jsd:jed,nk) :: rho_anomup
  real, dimension(isd:ied,jsd:jed,nk) :: rho_anomlo
  real, dimension(isd:ied,jsd:jed) :: diff_geo_x
  real, dimension(isd:ied,jsd:jed) :: diff_geo_y
  real, dimension(isd:ied,jsd:jed) :: tmp1
  real, dimension(isd:ied,jsd:jed) :: tmp2
  real    :: rhoT
  integer :: i, j, k
  integer :: tau, taup1

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_pressure_mod (press_grad_force_depth): module must be initialized')
  endif 

  tau    = Time%tau
  taup1  = Time%taup1
  tmp1   = 0.0
  tmp2   = 0.0
  wrk1_v = 0.0
  wrk2_v = 0.0
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  rho_anomup(:,:,:) = 0.0
  rho_anomlo(:,:,:) = 0.0

  ! use density anomaly to improve accuracy 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           ! Thickness weighted density of upper and lower part of grid box
           wrk2(i,j,k) = Thickness%dztup(i,j,k)*rho(i,j,k) + L_system%rho_dztup(i,j,k)
           wrk3(i,j,k) = Thickness%dztlo(i,j,k)*rho(i,j,k) + L_system%rho_dztlo(i,j,k)
           
           ! Density anomaly of full grid box
           rhoT        = ( wrk2(i,j,k)+wrk3(i,j,k) )/(Thickness%dztT(i,j,k,taup1)+epsln)
           wrk4(i,j,k) = Grd%tmask(i,j,k)*(rhoT - rho0)

           ! Density anomaly of upper and lower part of grid box
           rho_anomup(i,j,k) = Grd%tmask(i,j,k)*( -rho0 + (wrk2(i,j,k)/(Thickness%dztupT(i,j,k)+epsln)) )
           rho_anomlo(i,j,k) = Grd%tmask(i,j,k)*( -rho0 + (wrk3(i,j,k)/(Thickness%dztloT(i,j,k)+epsln)) )
        enddo
     enddo
  enddo
  
  ! compute hydrostatic pressure anomaly at T-point 
  ! ( Pa=N/m^2=kg/(m*s^2) )
  wrk1(:,:,:) = hydrostatic_pressure_blob(Thickness, rho_anomup(:,:,:), rho_anomlo(:,:,:))  

  call diagnose_3d(Time, Grd, id_anompress, wrk1(:,:,:))

  diff_geo_x(:,:) = 0.0
  diff_geo_y(:,:) = 0.0

  do k=1,nk

     do j=jsd,jed
        do i=isd,iec
           diff_geo_x(i,j) = -grav*(Thickness%geodepth_zt(i+1,j,k)-Thickness%geodepth_zt(i,j,k))
        enddo
     enddo
     do j=jsd,jec
        do i=isd,ied
           diff_geo_y(i,j) = -grav*(Thickness%geodepth_zt(i,j+1,k)-Thickness%geodepth_zt(i,j,k))
        enddo
     enddo

     ! density anomaly times geopotential gradient on k-levels 
     ! (dzurhoprime_geograd_x_klev,dzurhoprime_geograd_y_klev)
     tmp1(:,:) = FAY(FAX(wrk4(:,:,k))*diff_geo_x(:,:))
     tmp2(:,:) = FAX(FAY(wrk4(:,:,k))*diff_geo_y(:,:))
     do j=jsc,jec
        do i=isc,iec
           wrk1_v(i,j,k,1) = tmp1(i,j)*Grd%dxur(i,j)*Thickness%dzu(i,j,k)
           wrk1_v(i,j,k,2) = tmp2(i,j)*Grd%dyur(i,j)*Thickness%dzu(i,j,k)
        enddo
     enddo
     if(zero_correction_term_grad) then 
         wrk1_v(:,:,:,:) = 0.0
     endif

     ! gradient of anomalous baroclinic pressure on k-levels 
     ! (dzupgrad_x_klev,dzupgrad_y_klev)
     tmp1(:,:) = FDX_NT(FAY(wrk1(:,:,k)))
     tmp2(:,:) = FDY_ET(FAX(wrk1(:,:,k)))
     do j=jsc,jec
        do i=isc,iec
           wrk2_v(i,j,k,1) = tmp1(i,j)*Thickness%dzu(i,j,k)
           wrk2_v(i,j,k,2) = tmp2(i,j)*Thickness%dzu(i,j,k)
        enddo
     enddo

     if(zero_diagonal_press_grad) then 
         wrk2_v(:,:,:,:) = 0.0
     endif

     ! thickness weighted baroclinic pressure gradient on z-levels
     ! (dzu_pgrad_u,dzu_pgrad_v)
     do j=jsc,jec
        do i=isc,iec    
           Velocity%press_force(i,j,k,1) = -Grd%umask(i,j,k)*( wrk1_v(i,j,k,1) + wrk2_v(i,j,k,1) )
           Velocity%press_force(i,j,k,2) = -Grd%umask(i,j,k)*( wrk1_v(i,j,k,2) + wrk2_v(i,j,k,2) )
        enddo
     enddo

  enddo   ! k do-loop finish 

  ! (dzurhoprime_geograd_x_klev,dzurhoprime_geograd_y_klev)
  call diagnose_3d_u(Time, Grd, id_rhoprime_geograd_klev(1),wrk1_v(:,:,:,1))
  call diagnose_3d_u(Time, Grd, id_rhoprime_geograd_klev(2),wrk1_v(:,:,:,2))

  ! (dzupgrad_x_klev,dzupgrad_y_klev)
  call diagnose_3d_u(Time, Grd, id_pgrad_klev(1), wrk2_v(:,:,:,1))
  call diagnose_3d_u(Time, Grd, id_pgrad_klev(2), wrk2_v(:,:,:,2))

end subroutine press_grad_force_depth_blob
! </SUBROUTINE> NAME="press_grad_force_depth_blob"


!#######################################################################
! <SUBROUTINE NAME="press_grad_force_press_blob">
!
! <DESCRIPTION>
! This routine respects the partition between Eulerian system mass
! and the Lagrangian system mass associated with the Lagrangian blobs
! model.  The pressure gradient from the total (combined Eulerian and 
! Lagrangian systems) system is calculated. Aside from that, the routine 
! is the same as for press_grad_force_press.
!
! Compute the force from pressure using a finite difference method
! to compute the thickness weighted pressure gradient at the 
! velocity cell point. 
!
! For pressure-like vertical coordinates, we omit the bottom pressure
! and bottom geopotential.  Account is taken of variable partial cell
! thickness.  1 = dp/dx; 2 = dp/dy
!
! Thickness weight since this is what we wish to use in update of 
! the velocity. Resulting thickness weighted pressure gradient has 
! dimensions of Pa = N/m^2 = kg/(m*s^2).
!
! Thickness%dzu should be at tau.
!
! </DESCRIPTION>
!
subroutine press_grad_force_press_blob(Time, Thickness, Dens, Velocity, L_system, rho)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_velocity_type),      intent(inout) :: Velocity
  type(ocean_lagrangian_type),    intent(inout) :: L_system
  real, dimension(isd:,jsd:,:),   intent(in)    :: rho

  real, dimension(isd:ied,jsd:jed,nk) :: rho_anomup
  real, dimension(isd:ied,jsd:jed,nk) :: rho_anomlo
  real, dimension(isd:ied,jsd:jed) :: diff_press_x
  real, dimension(isd:ied,jsd:jed) :: diff_press_y
  real, dimension(isd:ied,jsd:jed) :: tmp1
  real, dimension(isd:ied,jsd:jed) :: tmp2
  real    :: rhoT
  integer :: i, j, k
  integer :: tau, taup1

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_pressure_mod (press_grad_force_press): module must be initialized')
  endif 

  tau    = Time%tau
  taup1  = Time%taup1
  tmp1   = 0.0
  tmp2   = 0.0
  wrk1_v = 0.0
  wrk2_v = 0.0
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  rho_anomup(:,:,:) = 0.0
  rho_anomlo(:,:,:) = 0.0

  ! use density anomaly to improve accuracy 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           ! Thickness weighted density of upper and lower part of grid box
           wrk2(i,j,k) = Thickness%dztup(i,j,k)*rho(i,j,k) + L_system%rho_dztup(i,j,k)
           wrk3(i,j,k) = Thickness%dztlo(i,j,k)*rho(i,j,k) + L_system%rho_dztlo(i,j,k)
           
           ! Density anomaly of full grid box
           rhoT        = ( wrk2(i,j,k)+wrk3(i,j,k) )/(Thickness%dztT(i,j,k,taup1)+epsln)
           wrk4(i,j,k) = Grd%tmask(i,j,k)*(rhoT - rho0)

           ! Density anomaly of upper and lower part of grid box
           rho_anomup(i,j,k) = Grd%tmask(i,j,k)*( -rho0 + wrk2(i,j,k)/(Thickness%dztupT(i,j,k)+epsln) )
           rho_anomlo(i,j,k) = Grd%tmask(i,j,k)*( -rho0 + wrk3(i,j,k)/(Thickness%dztloT(i,j,k)+epsln) )
        enddo
     enddo
  enddo

  ! compute geopotential anomaly (m^2/s^2) at T-cell point 
  wrk1(:,:,:) = geopotential_anomaly_blob(Thickness, rho_anomup(:,:,:), rho_anomlo(:,:,:)) 

  call diagnose_3d(Time, Grd, id_anomgeopot, wrk1(:,:,:))

  diff_press_x(:,:) = 0.0
  diff_press_y(:,:) = 0.0

  do k=1,nk

     do j=jsd,jed
        do i=isd,iec
           diff_press_x(i,j) = -dbars2pa*(Dens%pressure_at_depth(i+1,j,k)-Dens%pressure_at_depth(i,j,k))
        enddo
     enddo
     do j=jsd,jec
        do i=isd,ied
           diff_press_y(i,j) = -dbars2pa*(Dens%pressure_at_depth(i,j+1,k)-Dens%pressure_at_depth(i,j,k))
        enddo
     enddo

     ! anomalous density times gradient of pressure on k-level               
     ! (dzurhoprime_pgrad_x_klev,dzurhoprime_pgrad_y_klev)
     tmp1(:,:) = rho0r*FAY(FAX(wrk4(:,:,k))*diff_press_x(:,:))
     tmp2(:,:) = rho0r*FAX(FAY(wrk4(:,:,k))*diff_press_y(:,:))
     do j=jsc,jec
        do i=isc,iec
           wrk1_v(i,j,k,1) = tmp1(i,j)*Grd%dxur(i,j)*Thickness%dzu(i,j,k)
           wrk1_v(i,j,k,2) = tmp2(i,j)*Grd%dyur(i,j)*Thickness%dzu(i,j,k)
        enddo
     enddo

     ! density times gradient of anomalous geopotential on k-level
     ! (dzurho_geogradprime_x_klev,dzurho_geogradprime_y_klev) 
     tmp1(:,:) = FDX_NT(FAY(wrk1(:,:,k)))
     tmp2(:,:) = FDY_ET(FAX(wrk1(:,:,k)))
     do j=jsc,jec
        do i=isc,iec
           wrk2_v(i,j,k,1) = tmp1(i,j)*Thickness%rho_dzu(i,j,k,tau)
           wrk2_v(i,j,k,2) = tmp2(i,j)*Thickness%rho_dzu(i,j,k,tau)
        enddo
     enddo

     ! thickness weighted baroclinic pressure gradient on z-levels
     ! (dzu_pgrad_u,dzu_pgrad_v)
     do j=jsc,jec
        do i=isc,iec    
           Velocity%press_force(i,j,k,1) = -Grd%umask(i,j,k)*( wrk1_v(i,j,k,1) + wrk2_v(i,j,k,1) )
           Velocity%press_force(i,j,k,2) = -Grd%umask(i,j,k)*( wrk1_v(i,j,k,2) + wrk2_v(i,j,k,2) )
        enddo
     enddo

  enddo  ! k do-loop finish 


  ! (dzurhoprime_pgrad_x_klev,dzurhoprime_pgrad_y_klev)
  call diagnose_3d_u(Time, Grd, id_rhoprime_pgrad_klev(1),wrk1_v(:,:,:,1))
  call diagnose_3d_u(Time, Grd, id_rhoprime_pgrad_klev(2),wrk1_v(:,:,:,2))

  ! (dzurho_geogradprime_x_klev,dzurho_geogradprime_y_klev)
  call diagnose_3d_u(Time, Grd, id_rho_geogradprime_klev(1),wrk2_v(:,:,:,1))
  call diagnose_3d_u(Time, Grd, id_rho_geogradprime_klev(2),wrk2_v(:,:,:,2))

end subroutine press_grad_force_press_blob
! </SUBROUTINE> NAME="press_grad_force_press_blob"


!#######################################################################
! <FUNCTION NAME="pressure_in_dbars_blob">
!
! <DESCRIPTION>
! This routine respects the partition between Eulerian system mass
! and the Lagrangian system mass associated with the Lagrangian blobs
! model.  The pressure from the combined system is calculated.
! Aside from that, the routine is the same as for pressure_in_dbars.
!
! Compute pressure (dbars) exerted at T cell grid point by weight of
! water column above the grid point. 
!
! rho = density in kg/m^3
!
! psurf = surface pressure in Pa= kg/(m*s^2) = hydrostatic pressure 
! at z=0 associated with fluid between z=0 and z=eta_t.
! Also include pressure from atmosphere and ice, both of which 
! are part of the patm array. 
!
! This routine is used by ocean_density to compute the pressure 
! used in the equation of state.  It is only called when the 
! vertical coordinate is DEPTH_BASED.
!
! </DESCRIPTION>
function pressure_in_dbars_blob (Thickness, rho, Ldzrhoup, Ldzrholo, psurf)

  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: rho   
  real, dimension(isd:,jsd:,:), intent(in) :: Ldzrhoup
  real, dimension(isd:,jsd:,:), intent(in) :: Ldzrholo
  real, dimension(isd:,jsd:),   intent(in) :: psurf 
  real, dimension(isd:ied,jsd:jed,nk)      :: pressure_in_dbars_blob
  integer :: k

  if ( .not. module_is_initialized ) then   
    call mpp_error(FATAL, &
    '==>Error in ocean_pressure_mod (pressure_in_dbars): module must be initialized')
  endif 

  wrk1(:,:,:) = ( Thickness%dztup(:,:,:)*rho(:,:,:) + Ldzrhoup(:,:,:) )/(Thickness%dztupT(:,:,:) + epsln)
  wrk2(:,:,:) = ( Thickness%dztlo(:,:,:)*rho(:,:,:) + Ldzrholo(:,:,:) )/(Thickness%dztloT(:,:,:) + epsln)

  pressure_in_dbars_blob(:,:,:) = hydrostatic_pressure_blob (Thickness, wrk1(:,:,:), wrk2(:,:,:))*c2dbars

  ! add pressure from free surface height and loading from atmosphere and/or sea ice
  do k=1,nk
     pressure_in_dbars_blob(:,:,k) = pressure_in_dbars_blob(:,:,k) + psurf(:,:)*c2dbars 
  enddo


end function pressure_in_dbars_blob
! </FUNCTION> NAME="pressure_in_dbars_blob"


!#######################################################################
! <FUNCTION NAME="hydrostatic_pressure_blob">
!
! <DESCRIPTION>
! This routine respects the partition between Eulerian system mass
! and the Lagrangian system mass associated with the Lagrangian blobs
! model.  The hydrostatic pressure from the combined system is 
! calculated. Aside from that, the routine is the same as for 
! Hydrostatic_pressure.
!
! Hydrostatic pressure [Pa=N/m^2=kg/(m*s^2)] at T cell grid points.
!
! For GEOPOTENTIAL vertical coordinate, integration is 
! from z=0 to depth of grid point.  This integration results in  
! the so-called "baroclinic" pressure. 
!
! For ZSTAR, vertical coordinate, integration is from z=eta to 
! depth of grid point.  This is allowed because ZSTAR
! absorbs the undulations of the surface height into their definition.
!
! If the input density "rho" is an anomoly, the resulting presure 
! will be a hydrostatic pressure anomoly. If "rho" is full density, 
! the presure will be a full hydrostatic pressure.
!
! </DESCRIPTION>
!
function hydrostatic_pressure_blob (Thickness, rhoup, rholo)

  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: rhoup
  real, dimension(isd:,jsd:,:), intent(in) :: rholo

  integer :: i, j, k
  real, dimension(isd:ied,jsd:jed,nk) :: hydrostatic_pressure_blob

  if(vert_coordinate==GEOPOTENTIAL) then 
      do j=jsd,jed
         do i=isd,ied
            hydrostatic_pressure_blob(i,j,1) = rhoup(i,j,1)*grav*Grd%dzw(0)
         enddo
      enddo
  else !ZSTAR
      do j=jsd,jed
         do i=isd,ied
            hydrostatic_pressure_blob(i,j,1) = rhoup(i,j,1)*grav*Thickness%dzwtT(i,j,0)
         enddo
      enddo
  endif

  if(Thickness%method==ENERGETIC) then 

      do k=2,nk
         do j=jsd,jed
            do i=isd,ied
               hydrostatic_pressure_blob(i,j,k) = hydrostatic_pressure_blob(i,j,k-1) &
               + Grd%tmask(i,j,k)*p5grav                                   &
                *(rholo(i,j,k-1)+rhoup(i,j,k))*Thickness%dzwtT(i,j,k-1) 
            enddo
         enddo
      enddo

  elseif(Thickness%method==FINITEVOLUME) then 

      do k=2,nk
         do j=jsd,jed
            do i=isd,ied
               hydrostatic_pressure_blob(i,j,k) = hydrostatic_pressure_blob(i,j,k-1) &
               + Grd%tmask(i,j,k)*grav                                     &
                *(Thickness%dztloT(i,j,k-1)*rholo(i,j,k-1)+Thickness%dztupT(i,j,k)*rhoup(i,j,k)) 
            enddo
         enddo
      enddo

  endif

  ! For some debugging...we here drop the contribution from a 
  ! nonzero eta in the calculation of hydrostatic_pressure. 
  ! This is wrong physically, but it is useful for certain tests
  ! of the zstar vertical coordinate.  
  if(vert_coordinate==ZSTAR .and. zero_eta_over_h_zstar_pressure) then 
      do j=jsd,jed
         do i=isd,ied
            hydrostatic_pressure_blob(i,j,1) = rhoup(i,j,1)*grav*Thickness%dswt(i,j,0)
         enddo
      enddo
      do k=2,nk
         do j=jsd,jed
            do i=isd,ied
               hydrostatic_pressure_blob(i,j,k) = hydrostatic_pressure_blob(i,j,k-1) &
               + Grd%tmask(i,j,k)*p5grav                                   &
                *(rholo(i,j,k-1)+rhoup(i,j,k))*Thickness%dswt(i,j,k-1) 
            enddo
         enddo
      enddo
  endif 


end function hydrostatic_pressure_blob
! </FUNCTION> NAME="hydrostatic_pressure_blob"



!#######################################################################
! <FUNCTION NAME="geopotential_anomaly_blob">
!
! <DESCRIPTION>
! This routine respects the partition between Eulerian system mass
! and the Lagrangian system mass associated with the Lagrangian blobs
! model.  The geopotential from the combined system is calculated.
! Aside from that, the routine is the same as for geopotential_anomaly.
!
! Geopotential anomaly [m^2/s^2] at T cell grid points.
! Integration here is from z=-H to depth of grid point.  
!
! Input should be density anomaly rhoprime = rho-rho0.
!
! This function is needed when computing pressure gradient
! for PRESSURE_BASED vertical coordinates. 
!
! WARNING: Thickness%method==FINITEVOLUME has been found to be 
! problematic.  It remains under development.  It is NOT 
! supported for general use. 
! 
! </DESCRIPTION>
!
function geopotential_anomaly_blob (Thickness, rhoprimeup, rhoprimelo)

  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: rhoprimeup
  real, dimension(isd:,jsd:,:), intent(in) :: rhoprimelo

  integer :: i, j, k, kbot, kp1
  real, dimension(isd:ied,jsd:jed,nk) :: geopotential_anomaly_blob
  
  geopotential_anomaly_blob(:,:,:) = 0.0

  if(Thickness%method==ENERGETIC) then 
      do j=jsd,jed
         do i=isd,ied
            kbot=Grd%kmt(i,j)
            if(kbot > 1) then  
                k=kbot
                geopotential_anomaly_blob(i,j,k) = -grav_rho0r*rhoprimelo(i,j,k)*Thickness%dztloT(i,j,k)
                do k=kbot-1,1,-1
                   kp1 = k+1
                   geopotential_anomaly_blob(i,j,k) = geopotential_anomaly_blob(i,j,kp1) &
                   -p5grav_rho0r*(rhoprimeup(i,j,kp1)+rhoprimelo(i,j,k))*Thickness%dzwtT(i,j,k)
                enddo
            endif
         enddo
      enddo
  elseif(Thickness%method==FINITEVOLUME) then 
      do j=jsd,jed
         do i=isd,ied
            kbot=Grd%kmt(i,j)
            if(kbot > 1) then
                k=kbot
                geopotential_anomaly_blob(i,j,k) = -grav_rho0r*rhoprimelo(i,j,k)*Thickness%dztloT(i,j,k)
                do k=kbot-1,1,-1
                   kp1 = k+1
                   geopotential_anomaly_blob(i,j,k) = geopotential_anomaly_blob(i,j,kp1) &
                   -grav_rho0r*(Thickness%dztloT(i,j,k)*rhoprimelo(i,j,k)      &
                              + Thickness%dztupT(i,j,kp1)*rhoprimeup(i,j,kp1)) 
                enddo
            endif
         enddo
      enddo
  endif

end function geopotential_anomaly_blob
! </FUNCTION> NAME="geopotential_anomaly_blob"

end module ocean_pressure_mod

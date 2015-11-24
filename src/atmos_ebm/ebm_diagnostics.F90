!
!    Diagnostics for use in the Energy Balance Model
!     diagnostic tag is 'ebm'

! Add geopotential height to diagnostic fields

! ==========================================================================================

module ebm_diagnostics_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> B.L. Samuels
!</CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  Zhi Liang
! </REVIEWER>
!<OVERVIEW>
!  Diagnostics for use in the Energy Balance Model
!</OVERVIEW>
!<DESCRIPTION>
!  Diagnostics for use in the Energy Balance Model
!</DESCRIPTION>
! <INFO>
!   diagnostic tag is 'ebm'
! </INFO>

use          fms_mod, only: write_version_number
use    constants_mod, only: RADIAN
use   transforms_mod, only: get_grid_boundaries, get_deg_lon 
use   transforms_mod, only: get_deg_lat, grid_domain
use diag_manager_mod, only: diag_axis_init, register_diag_field
use diag_manager_mod, only: register_static_field, send_data
use time_manager_mod, only: time_type

implicit none
private

public :: ebm_diagnostics_init, ebm_diagnostics_up, ebm_diagnostics_down

character(len=84), parameter :: version = '$Id: ebm_diagnostics.F90,v 11.0 2004/09/28 19:10:27 fms Exp $'
character(len=84), parameter :: tagname='$Name: tikal $'
!-----------------------------------------------------------------------
!------------------------- axis names ----------------------------------
character(len=8) :: axiset   = 'ebm'
character(len=8) :: mod_name = 'ebm'
!-----------------------------------------------------------------------
! axis and field identifiers for the diag manager
integer :: id_t, id_q, id_prec, id_lprec, id_fprec, id_dt_qg
integer :: id_olr_toa, id_olr_boa, id_ilr_surf, id_nsw_toa, id_nsw_surf

contains

!-----------------------------------------------------------------------------------------------------------------
! <SUBROUTINE NAME="ebm_diagnostics_init">
!
! <OVERVIEW>
!             setup/write netcdf metadata and static fields
! </OVERVIEW>

! <DESCRIPTION>
!             setup/write netcdf metadata and static fields
! </DESCRIPTION>

!   <TEMPLATE>
!     call  ebm_diagnostics_init(Time, lon_max, lat_max,  axis_id)
!   </TEMPLATE>

! <IN NAME="Time" TYPE="type(time_type)">
!  Current/Init time.
! </IN>

! <IN NAME="lon_max" TYPE="integer">
!  longitude of spectral atmopsheric grid
! </IN>

! <IN NAME="lat_max" TYPE="integer">
!  latitude of spectral atmopsheric grid
! </IN>

! <OUT NAME="axis_id" TYPE="integer, dimension(:)">
!  axes identifiers
! </OUT>

subroutine ebm_diagnostics_init(Time, lon_max, lat_max,  axis_id)

type(time_type),         intent(in) :: Time
integer,                 intent(in) :: lon_max, lat_max
integer, intent(out),dimension(:)   :: axis_id

real, dimension(lon_max  ) :: lon
real, dimension(lon_max+1) :: lonb
real, dimension(lat_max  ) :: lat
real, dimension(lat_max+1) :: latb

integer :: axes_2d(2)
integer :: log_unit, id_lonb, id_lon, id_latb, id_lat
integer :: namelist_unit, ierr, io
logical :: used

!--- write out version information -------------------------------------
call write_version_number(version, tagname)

call get_grid_boundaries(lonb,latb,global=.true.)
call get_deg_lon(lon)
call get_deg_lat(lat)

id_lonb=diag_axis_init('lonb', RADIAN*lonb, 'degrees_E', 'x', 'longitude edges', set_name=axiset, Domain2=grid_domain)
id_latb=diag_axis_init('latb', RADIAN*latb, 'degrees_N', 'y', 'latitude edges',  set_name=axiset, Domain2=grid_domain)
id_lon =diag_axis_init('lon', lon, 'degrees_E', 'x', 'longitude', set_name=axiset, Domain2=grid_domain, edges=id_lonb)
id_lat =diag_axis_init('lat', lat, 'degrees_N', 'y', 'latitude',  set_name=axiset, Domain2=grid_domain, edges=id_latb)

axis_id(1) = id_lon
axis_id(2) = id_lat

axes_2d = axis_id(1:2)

id_t     = register_diag_field(mod_name, 'temp'  , axes_2d, Time, 'temperature'      , 'deg_k'    ) 
id_q     = register_diag_field(mod_name, 'sphum' , axes_2d, Time, 'specific humidity', 'none'     )
id_prec  = register_diag_field(mod_name, 'prec'  , axes_2d, Time, 'precipitation'    , 'Kg/(m2-s)')
id_lprec = register_diag_field(mod_name, 'lprec',axes_2d, Time, 'rain'             , 'Kg/(m2-s)')
id_fprec = register_diag_field(mod_name, 'fprec' , axes_2d, Time, 'snow'             ,' Kg/(m2-s)')

id_dt_qg = register_diag_field(mod_name, 'dt_qg' , axes_2d, Time, 'advection'             ,' 1/s')

id_olr_toa = register_diag_field(mod_name, 'olr_toa', axes_2d, Time, 'outgoing longwave (to space)', 'W/m2')
id_olr_boa = register_diag_field(mod_name, 'olr_boa', axes_2d, Time, 'outgoing longwave (to ocean)', 'W/m2')
id_ilr_surf = register_diag_field(mod_name, 'ilr_surf', axes_2d, Time, 'incoming longwave (from ocean)', 'W/m2')
id_nsw_toa = register_diag_field(mod_name, 'nsw_toa', axes_2d, Time, 'net shortwave (toa)', 'W/m2')
id_nsw_surf = register_diag_field(mod_name, 'nsw_surf', axes_2d, Time, 'net shortwave (to ocean)', 'W/m2')

return
end subroutine ebm_diagnostics_init

! </SUBROUTINE>

! <DIAGFIELDS>
!   <NETCDF NAME="temp" UNITS="deg_K">
!     atmospheric temperature
!   </NETCDF>
!   <NETCDF NAME="shpum" UNITS="none">
!     specific humidity
!   </NETCDF>
!   <NETCDF NAME="prec" UNITS="kg/(m2 s)">
!     total precipitation per unit ocean area
!   </NETCDF>
!   <NETCDF NAME="lprec" UNITS="kg/(m2 s)">
!     rain  per unit ocean area
!   </NETCDF>
!   <NETCDF NAME="fprec" UNITS="kg/(m2 s)">
!     snow  per unit ocean area
!   </NETCDF>
!   <NETCDF NAME="dt_qg" UNITS="1/s">
!     advection 
!   </NETCDF>
!   <NETCDF NAME="olr_toa" UNITS="W/m2">
!     outgoing longwave (to space)
!   </NETCDF>
!   <NETCDF NAME="olr_boa" UNITS="W/m2">
!     outgoing longwave (to ocean)
!   </NETCDF>
!   <NETCDF NAME="ilr_surf" UNITS="W/m2">
!     incoming longwave (from ocean)
!   </NETCDF>
!   <NETCDF NAME="nsw_toa" UNITS="W/m2">
!     net shortwave (toa)
!   </NETCDF>
!   <NETCDF NAME="nsw_surf" UNITS="W/m2">
!     net shortwave (to ocean)
!   </NETCDF>
! </DIAGFIELDS>

!--------------------------------------------------------------------------------------------
! <SUBROUTINE NAME="ebm_diagnostics_up">
!
! <OVERVIEW>
!  Diagnostics passed from the "ocean" surface up thru the atmosphere
! </OVERVIEW>

! <DESCRIPTION>
!  Diagnostics passed from the "ocean" surface up thru the atmosphere
! </DESCRIPTION>

!   <TEMPLATE>
!     call  ebm_diagnostics_up(Time,tg, qg, lprec, fprec, dt_qg )
!   </TEMPLATE>
!
! <IN NAME="Time" TYPE="type(time_type)">
!  time at the current time level (tau)
! </IN>

! <IN NAME="tg" TYPE="real, dimension(:,:)" >
!  atmospheric temperature  ( deg_k )
! </IN>

! <IN NAME="qg" TYPE="real, dimension(:,:)" >
!  specific humidity ( no units )
! </IN>

! <IN NAME="lprec" TYPE="real, dimension(:,:)" >
!  liquid precipitation rate (rain) in kg/m2/s
! </IN>

! <IN NAME="fprec" TYPE="real, dimension(:,:)" >
!  frozen precipitation rate (snow) in kg/m2/s
! </IN>

! <IN NAME="dt_qg" TYPE="real, dimension(:,:)" >
!  rate of change in specific humidity (1/s)
! </IN>

subroutine ebm_diagnostics_up(Time,tg, qg, lprec, fprec, dt_qg )

type(time_type),        intent(in) :: Time
real, intent(in), dimension(:,:)   :: tg, qg, lprec, fprec, dt_qg

logical :: used


if(id_t      > 0) used = send_data(id_t      , tg          , time)
if(id_q      > 0) used = send_data(id_q      , qg          , time)
if(id_prec   > 0) used = send_data(id_prec   , lprec+fprec , time)
if(id_lprec  > 0) used = send_data(id_lprec  , lprec       , time)
if(id_fprec  > 0) used = send_data(id_fprec  , fprec       , time)
if(id_dt_qg  > 0) used = send_data(id_dt_qg  , dt_qg       , time)

return
end subroutine ebm_diagnostics_up
! </SUBROUTINE>

!--------------------------------------------------------------------------------------------
! <SUBROUTINE NAME="ebm_diagnostics_down">
!
! <OVERVIEW>
!   Diagnostics passed down from the atmosphere to the "ocean" surface
! </OVERVIEW>

! <DESCRIPTION>
!   Diagnostics passed down from the atmosphere to the "ocean" surface
! </DESCRIPTION>

!   <TEMPLATE>
!     call  ebm_diagnostics_down(Time,olr_toa,olr_boa,ilr_surf,nsw_toa,nsw_surf)
!   </TEMPLATE>
!
! <IN NAME="Time" TYPE="type(time_type)">
!  time at the current time level (tau)
! </IN>


subroutine ebm_diagnostics_down(Time,olr_toa,olr_boa,ilr_surf,nsw_toa,nsw_surf)

type(time_type), intent(in) :: Time
real, intent(in), dimension(:,:)   :: olr_toa,olr_boa,ilr_surf,nsw_toa,nsw_surf

logical :: used

if(id_olr_toa  > 0) used = send_data(id_olr_toa ,olr_toa  ,time)
if(id_olr_boa  > 0) used = send_data(id_olr_boa ,olr_boa  ,time)
if(id_ilr_surf > 0) used = send_data(id_ilr_surf,ilr_surf ,time)
if(id_nsw_toa  > 0) used = send_data(id_nsw_toa ,nsw_toa  ,time)
if(id_nsw_surf > 0) used = send_data(id_nsw_surf,nsw_surf ,time)

return
end subroutine ebm_diagnostics_down
! </SUBROUTINE>

!--------------------------------------------------------------------------------------------

end module ebm_diagnostics_mod








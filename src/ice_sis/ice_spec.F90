!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_spec_mod - sea ice and SST specified from data as per GFDL climate group !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_spec_mod

use fms_mod, only: open_namelist_file, check_nml_error, close_file, &
                   stdlog, stdout, mpp_pe, mpp_root_pe, write_version_number
use mpp_mod, only: input_nml_file

use time_manager_mod, only: time_type, get_date, set_date
use data_override_mod,only: data_override
use constants_mod,    only: Tfreeze

implicit none
include 'netcdf.inc'
private
public :: get_sea_surface

character(len=128), parameter :: version = '$Id: ice_spec.F90,v 20.0 2013/12/13 23:28:29 fms Exp $'
character(len=128), parameter :: tagname = '$Name: tikal $'

logical :: module_is_initialized = .false.

logical :: mcm_ice = .false. ! When mcm_ice=.true., ice is handled as in supersource
real    :: sst_pert = 0.     ! global temperature perturbation used for sensitivity experiments

real    :: minimum_ice_concentration = 0.2
real    :: minimum_ice_thickness     = 1.0
logical :: do_leads = .true.   ! when do_leads=false there is no fractional ice concentration
                               ! also you should set the minimum_ice_concentration = 0.5
logical :: sst_degk = .false.  ! when sst_degk=true the input sst data is in degrees Kelvin
                               ! otherwise it is assumed to be in degrees Celsius

!amip date for repeating single day (rsd) option
 integer :: repeat_date(3)=(/-1,-1,-1/)

namelist /ice_spec_nml/ mcm_ice, do_leads, minimum_ice_concentration, &
                          minimum_ice_thickness, sst_degk, sst_pert, repeat_date

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! get_sea_surface - get SST, ice concentration and thickness from data         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine get_sea_surface(Time, ts, cn, iceh)
type (time_type),                         intent(in)  :: Time
real, dimension(:, :),                    intent(out) :: ts
real, dimension(size(ts,1),size(ts,2),2), intent(out) :: cn
real, dimension(size(ts,1),size(ts,2)),   intent(out) :: iceh

real, dimension(size(ts,1),size(ts,2))                :: sst, icec

real ::  t_sw_freeze0 = -1.8
real ::  t_sw_freeze
integer :: ierr, io, unit
type(time_type) :: Spec_Time
integer :: tod(3),dum
integer :: stdoutunit,stdlogunit

stdoutunit=stdout()
stdlogunit=stdlog()

if(.not.module_is_initialized) then
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ice_spec_nml, iostat=io)
#else
    unit = open_namelist_file()
    read  (unit, ice_spec_nml,iostat=io)
    call close_file (unit)
#endif
    ierr = check_nml_error(io,'ice_spec_nml')
    write (stdoutunit,'(/)')
    write (stdoutunit, ice_spec_nml)
    write (stdlogunit, ice_spec_nml)

   call write_version_number(version, tagname)
   module_is_initialized = .true.
endif

! modify time repeating single day option
  if (all(repeat_date>0)) then
     call get_date(Time,dum,dum,dum,tod(1),tod(2),tod(3))
     Spec_Time = set_date(repeat_date(1),repeat_date(2),repeat_date(3),tod(1),tod(2),tod(3))
  else
     Spec_Time = Time
  endif

  t_sw_freeze = t_sw_freeze0
  if (sst_degk) then
     t_sw_freeze = t_sw_freeze0 + Tfreeze ! convert sea water freeze point to degK
  endif
  icec = 0.0; iceh = 0.0; sst = t_sw_freeze;
  call data_override('ICE', 'sic_obs', icec, Spec_Time)
  call data_override('ICE', 'sit_obs', iceh, Spec_Time)
  call data_override('ICE', 'sst_obs', sst,  Spec_Time)

  if(mcm_ice) then
    icec = 0.0
!   TK Mod: Limit minimum non-zero sea ice thickness to 0.01m.
!           This is to eliminate some very thin but non-zero
!           sea ice thickness values, where they really should be zero
!           but have become nonzero due to spatial interpolation
!           where the input grid and model grid are not
!           EXACTLY the same.  0.01 was obtained by trial and
!           error to roughly match supersource behavior.
!           5/22/01; 8/23/01

    where (iceh < 0.01) iceh=0.0
    where (iceh>0.0)
      icec = 1.0
      sst  = t_sw_freeze
    end where
  else
    where (icec >= minimum_ice_concentration)
      iceh = max(iceh, minimum_ice_thickness)
      sst = t_sw_freeze
    elsewhere
      icec = 0.0
      iceh = 0.0
    end where
    if (.not.do_leads) then
       where (icec >= minimum_ice_concentration) icec = 1.0
    endif
  endif

  where (icec==0.0 .and. sst<=t_sw_freeze) sst = t_sw_freeze+1e-10

! add on non-zero sea surface temperature perturbation (namelist option)
! this perturbation may be useful in accessing model sensitivities

   if ( abs(sst_pert) > 0.0001 ) then
      sst = sst + sst_pert
   endif

  cn(:,:,2) = icec
  cn(:,:,1) = 1-cn(:,:,2)

  if (sst_degk) then
     ts = sst
  else
     ts = sst+Tfreeze
  endif

  return
end subroutine get_sea_surface

end module ice_spec_mod

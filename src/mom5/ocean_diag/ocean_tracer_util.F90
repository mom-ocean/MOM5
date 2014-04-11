module ocean_tracer_util_mod
#define COMP isc:iec,jsc:jec
!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
! S. M. Griffies 
! </CONTACT>
!
!<OVERVIEW>
! This module contains many routines of use for tracer diagnostics in MOM. 
!</OVERVIEW>
!
!<DESCRIPTION>
! Tracer utility module for MOM. Of use for tracer diagnostics. 
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_tracer_util_nml">
!  <DATA NAME="rebin_onto_rho_all_values" TYPE="logical">
!  Set true to if wish to bin all values into density classes, even 
!  those cells whose density is outside the max and min range of the 
!  density bins.  The default is rebin_onto_rho_all_values=.true.,
!  which means those cells with extreme density values will be included. 
!  This default is consistent with the default computation of 
!  transport_on_nrho.  
!  </DATA> 
!
!  <DATA NAME="debug_diagnose_mass_of_layer" TYPE="logical">
!  To help debug the algorithm to diagnose mass of fluid within
!  a neutral density layer. 
!  Default:  debug_diagnose_mass_of_layer=.false.
!  </DATA> 
!  <DATA NAME="epsln_diagnose_mass_of_layer" TYPE="real">
!  Relative mass difference allowable between layer and level total mass. 
!  Default: epsln_diagnose_mass_of_layer=1e-4.
!  </DATA> 
!
!</NAMELIST>
!
use constants_mod,        only: epsln
use fms_mod,              only: open_namelist_file, check_nml_error, close_file
use mpp_mod,              only: input_nml_file, stdout, stdlog, FATAL
use mpp_mod,              only: mpp_error, mpp_pe, mpp_min, mpp_max
use platform_mod,         only: i8_kind
use diag_manager_mod,     only: send_data
  
use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: ADVECT_PSOM
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,      only: ocean_time_type, ocean_thickness_type 
use ocean_types_mod,      only: ocean_density_type
use ocean_util_mod,       only: write_timestamp, write_chksum_3d
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk4
use ocean_workspace_mod,  only: wrk1_v, wrk2_v, wrk3_v
use ocean_workspace_mod,  only: wrk1_2d, wrk2_2d , wrk3_2d 

implicit none

private

character(len=256) :: version='CVS $Id'
character(len=256) :: tagname='Tag $Name'

! for output
integer :: unit=6

logical :: use_blobs

#include <ocean_memory.h>

logical :: rebin_onto_rho_all_values    = .true.
logical :: module_is_initialized        = .false.
logical :: debug_diagnose_mass_of_layer = .false.
real    :: epsln_diagnose_mass_of_layer = 1e-5

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()

public ocean_tracer_util_init
public tracer_min_max
public dzt_min_max
public tracer_prog_chksum 
public tracer_diag_chksum 
public tracer_psom_chksum
public sort_pick_array
public sort_shell_array
public rebin_onto_rho
public diagnose_mass_of_layer
public diagnose_3d_rho

namelist /ocean_tracer_util_nml/ rebin_onto_rho_all_values, &
          debug_diagnose_mass_of_layer, epsln_diagnose_mass_of_layer


contains


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_util_init">
!
! <DESCRIPTION>
! Initialize MOM tracer utilities.
! </DESCRIPTION>
!
subroutine ocean_tracer_util_init (Grid, Domain, blobs)

  type(ocean_grid_type),   intent(in), target :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  logical,                 intent(in)         :: blobs

  integer :: ioun, io_status
#ifdef INTERNAL_FILE_NML
  integer :: ierr
#endif
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 


  if (module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error in ocean_tracer_util_mod (ocean_tracer_util_init): module already initialized')
  endif 

  module_is_initialized = .true.
  stdlogunit=stdlog()
  write( stdlogunit,'(/a/)') trim(version)

  use_blobs = blobs

  Dom => Domain
  Grd => Grid

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk=Grd%nk
#endif

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_tracer_util_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_tracer_util_nml')
#else
  ioun = open_namelist_file()
  read(ioun, ocean_tracer_util_nml, iostat=io_status)
  call close_file(ioun)
#endif
  write (stdlogunit, ocean_tracer_util_nml)
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_tracer_util_nml)

  if(rebin_onto_rho_all_values) then
      write (stdoutunit,'(/a)') &
      'Note: ocean_tracer_util: rebin_onto_rho will include density values outside bounds range.'
  else
      write (stdoutunit,'(/a)') &
      'Note: ocean_tracer_util: rebin_onto_rho will NOT include density values outside bounds range.'
  endif


end subroutine ocean_tracer_util_init
! </SUBROUTINE> NAME="ocean_tracer_util_init">



!#######################################################################
! <SUBROUTINE NAME="tracer_min_max">
!
! <DESCRIPTION>
! Compute the global min and max for tracers.  
!
! Vectorized using maxloc() and minloc() intrinsic functions by 
! Russell.Fiedler@csiro.au (May 2005).
!
! Modified by Zhi.Liang (July 2005)
!          
! </DESCRIPTION>
!
subroutine tracer_min_max(Time, Thickness, Tracer)
  
  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: Tracer

  real    :: tmax, tmin, tmax0, tmin0
  integer :: itmax, jtmax, ktmax, itmin, jtmin, ktmin
  integer :: tau 
  real    :: fudge

  ! arrays to enable vectorization
  integer :: iminarr(3),imaxarr(3)

  if (.not. module_is_initialized) then
     call mpp_error(FATAL,&
     '==>Error from ocean_tracer_util_mod (tracer_min_max): module not initialized')
  endif 

  tmax=-1.e10;tmin=1.e10
  itmax=0;jtmax=0;ktmax=0
  itmin=0;jtmin=0;ktmin=0

  tau = Time%tau

  call write_timestamp(Time%model_time)
  
  wrk1(isc:iec,jsc:jec,:) = Tracer%field(isc:iec,jsc:jec,:,tau)
  
  if(ANY(Grd%tmask(isc:iec,jsc:jec,:) > 0.)) then
     iminarr=minloc(wrk1(isc:iec,jsc:jec,:),Grd%tmask(isc:iec,jsc:jec,:) > 0.)
     imaxarr=maxloc(wrk1(isc:iec,jsc:jec,:),Grd%tmask(isc:iec,jsc:jec,:) > 0.)
     itmin=iminarr(1)+isc-1
     jtmin=iminarr(2)+jsc-1
     ktmin=iminarr(3)
     itmax=imaxarr(1)+isc-1
     jtmax=imaxarr(2)+jsc-1 
     ktmax=imaxarr(3)
     tmin=wrk1(itmin,jtmin,ktmin)
     tmax=wrk1(itmax,jtmax,ktmax)
  end if

  ! use "fudge" to distinguish processors when tracer extreme is independent of processor
  fudge = 1.0 + 1.e-12*mpp_pe() 
  tmax = tmax*fudge
  tmin = tmin*fudge
  if(tmax == 0.0) then 
    tmax = tmax + 1.e-12*mpp_pe() 
  endif 
  if(tmin == 0.0) then 
    tmin = tmin + 1.e-12*mpp_pe() 
  endif 
  

  tmax0=tmax;tmin0=tmin

  call mpp_max(tmax)
  call mpp_min(tmin)

  if (tmax0 == tmax) then
      if (trim(Tracer%name) == 'temp') then
          write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
               ' The maximum T is ', tmax,&
               ' deg C at (i,j,k) = (',itmax+Dom%ioff,',',jtmax+Dom%joff,',',ktmax,&
               '),  (lon,lat,dpt) = (',Grd%xt(itmax,jtmax),',', &
               Grd%yt(itmax,jtmax),',', Grd%zt(ktmax),' m)'
          if (use_blobs) then
             write(unit,'(a,es22.12,a,es22.12,a,es22.12,a,a,es22.12,a)')      &
                  ' The grid dimensions are (dxt, dyt, dztT), dzt = (', &
                  Grd%dxt(itmax,jtmax),',', Grd%dyt(itmax,jtmax),',', &
                  Thickness%dztT(itmax,jtmax,ktmax,tau),' m)',',', Thickness%dzt(itmax,jtmax,ktmax),' m'
             write(unit,'(a,es22.12,a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dztT,rho_dzt) = ', &
                  Thickness%dst(itmax,jtmax,ktmax),', ',Thickness%rho_dztT(itmax,jtmax,ktmax,tau),', ',&
                  Thickness%rho_dzt(itmax,jtmax,ktmax,tau)
          else
             write(unit,'(a,es22.12,a,es22.12,a,es22.12,a)')      &
                  ' The grid dimensions are (dxt, dyt, dzt) = (', &
                  Grd%dxt(itmax,jtmax),',', Grd%dyt(itmax,jtmax),',', Thickness%dzt(itmax,jtmax,ktmax),' m)'
             write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
                  Thickness%dst(itmax,jtmax,ktmax),', ',Thickness%rho_dzt(itmax,jtmax,ktmax,tau)
          endif
          write(unit,'(a,i6)') ' And the number of cells in the column are kmt = ', Grd%kmt(itmax,jtmax)
          if (tmax > Tracer%max_tracer) then
              call mpp_error(FATAL, &
               '==>Error from ocean_tracer_util_mod: The maximum temperature is outside allowable range')
          endif
      else if (trim(Tracer%name) == 'salt') then
          write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
               ' The maximum S is ', tmax,&
               ' psu at (i,j,k) = (',itmax+Dom%ioff,',',jtmax+Dom%joff,',',ktmax,&
               '),  (lon,lat,dpt) = (',Grd%xt(itmax,jtmax),',', &
               Grd%yt(itmax,jtmax),',', Grd%zt(ktmax),' m)'
          if (use_blobs) then
             write(unit,'(a,es22.12,a,es22.12,a,es22.12,a,a,es22.12,a)')      &
                  ' The grid dimensions are (dxt, dyt, dztT), dzt = (', &
                  Grd%dxt(itmax,jtmax),',', Grd%dyt(itmax,jtmax),',', Thickness%dztT(itmax,jtmax,ktmax,tau),' m)',&
                  ',', Thickness%dzt(itmax,jtmax,ktmax),' m'
             write(unit,'(a,es22.12,a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dztT,rho_dzt) = ', &
                  Thickness%dst(itmax,jtmax,ktmax),', ',Thickness%rho_dztT(itmax,jtmax,ktmax,tau),', '&
                  ,Thickness%rho_dzt(itmax,jtmax,ktmax,tau)
          else
             write(unit,'(a,es22.12,a,es22.12,a,es22.12,a)')      &
                  ' The grid dimensions are (dxt, dyt, dzt) = (', &
                  Grd%dxt(itmax,jtmax),',', Grd%dyt(itmax,jtmax),',', Thickness%dzt(itmax,jtmax,ktmax),' m)'
             write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
                  Thickness%dst(itmax,jtmax,ktmax),', ',Thickness%rho_dzt(itmax,jtmax,ktmax,tau)
          endif
          write(unit,'(a,i6)') ' And the number of cells in the column are kmt = ', Grd%kmt(itmax,jtmax)
          if (tmax > Tracer%max_tracer) then
              call mpp_error(FATAL,&
              '==>Error from ocean_tracer_util_mod: The maximum salinity is outside allowable range')
          endif
      else
          write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
               ' The maximum '//trim(Tracer%name)// ' is ', tmax,&
                ' '// trim(Tracer%units)// ' at (i,j,k) = (',itmax+Dom%ioff,',',jtmax+Dom%joff,','&
               ,ktmax,'),  (lon,lat,dpt) = (',Grd%xt(itmax,jtmax),',', &
               Grd%yt(itmax,jtmax),',', Grd%zt(ktmax),' m)'
          if (use_blobs) then
             write(unit,'(a,es22.12,a,es22.12,a,es22.12,a,a,es22.12,a)')      &
                  ' The grid dimensions are (dxt, dyt, dztT), dzt = (', &
                  Grd%dxt(itmax,jtmax),',', Grd%dyt(itmax,jtmax),',', Thickness%dztT(itmax,jtmax,ktmax,tau),' m)',&
                  ',', Thickness%dzt(itmax,jtmax,ktmax),' m'
             write(unit,'(a,es22.12,a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dztT,rho_dzt) = ', &
                  Thickness%dst(itmin,jtmin,ktmin),', ',Thickness%rho_dztT(itmin,jtmin,ktmin,tau),', ',&
                  Thickness%rho_dzt(itmin,jtmin,ktmin,tau)
          else
             write(unit,'(a,es22.12,a,es22.12,a,es22.12,a)')      &
                  ' The grid dimensions are (dxt, dyt, dzt) = (', &
                  Grd%dxt(itmax,jtmax),',', Grd%dyt(itmax,jtmax),',', Thickness%dzt(itmax,jtmax,ktmax),' m)'
             write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
                  Thickness%dst(itmin,jtmin,ktmin),', ',Thickness%rho_dzt(itmin,jtmin,ktmin,tau)
          endif
          write(unit,'(a,i6)') ' And the number of cells in the column are kmt = ', Grd%kmt(itmax,jtmax)
          if (tmax > Tracer%max_tracer) then
              call mpp_error(FATAL, &
               '==>Error from ocean_tracer_util_mod: The maximum tracer is outside allowable range')
          endif 
      endif
  endif
  
  if (tmin0 == tmin) then
      if (trim(Tracer%name) == 'temp') then
          write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
               ' The minimum T is ', tmin,&
               ' deg C at (i,j,k) = (',itmin+Dom%ioff,',',jtmin+Dom%joff,',',ktmin,&
               '),  (lon,lat,dpt) = (',Grd%xt(itmin,jtmin),',', &
               Grd%yt(itmin,jtmin),',', Grd%zt(ktmin),' m)'
          if (use_blobs) then
             write(unit,'(a,es22.12,a,es22.12,a,es22.12,a,a,es22.12,a)')      &
                  ' The grid dimensions are (dxt, dyt, dztT), dzt = (', &
                  Grd%dxt(itmin,jtmin),',', Grd%dyt(itmin,jtmin),',', Thickness%dztT(itmin,jtmin,ktmin,tau),' m)',','&
                  , Thickness%dzt(itmin,jtmin,ktmin),' m'
             write(unit,'(a,es22.12,a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dztT,rho_dzt) = ', &
                  Thickness%dst(itmin,jtmin,ktmin),', ',Thickness%rho_dztT(itmin,jtmin,ktmin,tau),', ',&
                  Thickness%rho_dzt(itmin,jtmin,ktmin,tau)
          else
             write(unit,'(a,es22.12,a,es22.12,a,es22.12,a)')      &
                  ' The grid dimensions are (dxt, dyt, dzt) = (', &
                  Grd%dxt(itmin,jtmin),',', Grd%dyt(itmin,jtmin),',', Thickness%dzt(itmin,jtmin,ktmin),' m)'
             write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
                  Thickness%dst(itmin,jtmin,ktmin),', ',Thickness%rho_dzt(itmin,jtmin,ktmin,tau)
          endif
          write(unit,'(a,i6)') ' And the number of cells in the column are kmt = ', Grd%kmt(itmin,jtmin)
          if (tmin < Tracer%min_tracer) then
              call mpp_error(FATAL, &
               '==>Error from ocean_tracer_util_mod (tracer_min_max): minimum temp outside allowable range')
          endif
      else if (trim(Tracer%name) == 'salt') then
          write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
               ' The minimum S is ', tmin,&
               ' psu at (i,j,k) = (',itmin+Dom%ioff,',',jtmin+Dom%joff,',',ktmin,&
               '),  (lon,lat,dpt) = (',Grd%xt(itmin,jtmin),',', &
               Grd%yt(itmin,jtmin),',', Grd%zt(ktmin),' m)'
          if (use_blobs) then
             write(unit,'(a,es22.12,a,es22.12,a,es22.12,a,a,es22.12,a)')      &
                  ' The grid dimensions are (dxt, dyt, dztT) = (', &
                  Grd%dxt(itmin,jtmin),',', Grd%dyt(itmin,jtmin),',', Thickness%dztT(itmin,jtmin,ktmin,tau),' m)',','&
                  , Thickness%dzt(itmin,jtmin,ktmin),' m'
             write(unit,'(a,es22.12,a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dztT,rho_dzt) = ', &
                  Thickness%dst(itmin,jtmin,ktmin),', ',Thickness%rho_dztT(itmin,jtmin,ktmin,tau),', ',&
                  Thickness%rho_dzt(itmin,jtmin,ktmin,tau)
          else
             write(unit,'(a,es22.12,a,es22.12,a,es22.12,a)')      &
                  ' The grid dimensions are (dxt, dyt, dzt) = (', &
                  Grd%dxt(itmin,jtmin),',', Grd%dyt(itmin,jtmin),',', Thickness%dzt(itmin,jtmin,ktmin),' m)'
             write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
                  Thickness%dst(itmin,jtmin,ktmin),', ',Thickness%rho_dzt(itmin,jtmin,ktmin,tau)
          endif
          write(unit,'(a,i6)') ' And the number of cells in the column are kmt = ', Grd%kmt(itmin,jtmin)
          if (tmin < Tracer%min_tracer) then
              call mpp_error(FATAL, &
               '==>Error from ocean_tracer_util_mod: The minimum salinity is outside allowable range')
          endif
      else
          write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
               ' The minimum '//trim(Tracer%name)// ' is ', tmin,&
                ' '// trim(Tracer%units)// ' at (i,j,k) = (',itmin+Dom%ioff,',',jtmin+Dom%joff,','&
               ,ktmin,'),  (lon,lat,dpt) = (',Grd%xt(itmin,jtmin),',', &
               Grd%yt(itmin,jtmin),',', Grd%zt(ktmin),' m)'
          if (use_blobs) then
             write(unit,'(a,es22.12,a,es22.12,a,es22.12,a,a,es22.12,a)')      &
                  ' The grid dimensions are (dxt, dyt, dztT), dzt = (', &
                  Grd%dxt(itmin,jtmin),',', Grd%dyt(itmin,jtmin),',', Thickness%dztT(itmin,jtmin,ktmin,tau),' m)',',',&
                  Thickness%dzt(itmin,jtmin,ktmin),' m'
             write(unit,'(a,es22.12,a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dztT,rho_dzt) = ', &
                  Thickness%dst(itmin,jtmin,ktmin),', ',Thickness%rho_dztT(itmin,jtmin,ktmin,tau),', ',&
                  Thickness%rho_dzt(itmin,jtmin,ktmin,tau)
          else
             write(unit,'(a,es22.12,a,es22.12,a,es22.12,a)')      &
                  ' The grid dimensions are (dxt, dyt, dzt) = (', &
                  Grd%dxt(itmin,jtmin),',', Grd%dyt(itmin,jtmin),',', Thickness%dzt(itmin,jtmin,ktmin),' m)'
             write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
                  Thickness%dst(itmin,jtmin,ktmin),', ',Thickness%rho_dzt(itmin,jtmin,ktmin,tau)
          endif
          write(unit,'(a,i6)') ' And the number of cells in the column are kmt = ', Grd%kmt(itmin,jtmin)
          if (tmin < Tracer%min_tracer) then
              call mpp_error(FATAL, &
               '==>Error from ocean_tracer_util_mod (tracer_min_max): minimum tracer outside allowable range')
          endif
      endif
  endif

  return


end subroutine tracer_min_max
! </SUBROUTINE>  NAME="tracer_min_max"


!#######################################################################
! <SUBROUTINE NAME="dzt_min_max">
!
! <DESCRIPTION>
! Compute the global min and max for dzt.  
!
! Modified by Stephen.Griffies from subroutine tracer_min_max
!          
! </DESCRIPTION>
!
subroutine dzt_min_max(Time, Thickness, filecaller)
  
  type(ocean_time_type),      intent(in) :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  character(len=*),           intent(in) :: filecaller

  real    :: dztmax, dztmin, dztmax0, dztmin0
  integer :: itmax, jtmax, ktmax, itmin, jtmin, ktmin
  integer :: tau
  real    :: fudge

  ! arrays to enable vectorization
  integer :: iminarr(3),imaxarr(3)

  tau = Time%tau

  dztmax=-1.e10;dztmin=1.e10
  itmax=0;jtmax=0;ktmax=0
  itmin=0;jtmin=0;ktmin=0

  call write_timestamp(Time%model_time)

  if (use_blobs) then
     wrk1(isc:iec,jsc:jec,:) = Thickness%dztT(isc:iec,jsc:jec,:,tau)
  else
     wrk1(isc:iec,jsc:jec,:) = Thickness%dzt(isc:iec,jsc:jec,:)
  endif
  
  if(ANY(Grd%tmask(isc:iec,jsc:jec,:) > 0.)) then
     iminarr=minloc(wrk1(isc:iec,jsc:jec,:),Grd%tmask(isc:iec,jsc:jec,:) > 0.)
     imaxarr=maxloc(wrk1(isc:iec,jsc:jec,:),Grd%tmask(isc:iec,jsc:jec,:) > 0.)
     itmin=iminarr(1)+isc-1
     jtmin=iminarr(2)+jsc-1
     ktmin=iminarr(3)
     itmax=imaxarr(1)+isc-1
     jtmax=imaxarr(2)+jsc-1 
     ktmax=imaxarr(3)
     dztmin=wrk1(itmin,jtmin,ktmin)
     dztmax=wrk1(itmax,jtmax,ktmax)
  end if

  ! use "fudge" to distinguish processors when extreme is independent of processor
  fudge = 1.0 + 1.e-12*mpp_pe() 
  dztmax = dztmax*fudge
  dztmin = dztmin*fudge
  if(dztmax == 0.0) then 
    dztmax = dztmax + 1.e-12*mpp_pe() 
  endif 
  if(dztmin == 0.0) then 
    dztmin = dztmin + 1.e-12*mpp_pe() 
  endif 

  dztmax0=dztmax;dztmin0=dztmin

  call mpp_max(dztmax)
  call mpp_min(dztmin)

  if (dztmax0 == dztmax) then

      write(unit,'(/a)') trim(filecaller)
      call write_timestamp(Time%model_time)

      write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
           ' The maximum dzt is ', dztmax,&
           ' metre at (i,j,k) = (',itmax+Dom%ioff,',',jtmax+Dom%joff,',',ktmax,&
           '),  (lon,lat,dpt) = (',Grd%xt(itmax,jtmax),',', &
           Grd%yt(itmax,jtmax),',', Thickness%depth_zt(itmax,jtmax,ktmax),' m)'
      if (use_blobs) then
         write(unit,'(a,es22.12,a,es22.12,a,es22.12)')      &
              ' The grid dimensions are (dxt, dyt, dztT) = (', &
              Grd%dxt(itmax,jtmax),',', Grd%dyt(itmax,jtmax),',', Thickness%dztT(itmax,jtmax,ktmax,tau),' m)'
         write(unit,'(a,i6)') ' The number of cells in the column are kmt = ', Grd%kmt(itmax,jtmax)
         write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dztT) = ', &
              Thickness%dst(itmax,jtmax,ktmax),', ',Thickness%rho_dztT(itmax,jtmax,ktmax,tau)
      else
         write(unit,'(a,es22.12,a,es22.12,a,es22.12)')      &
              ' The grid dimensions are (dxt, dyt, dzt) = (', &
              Grd%dxt(itmax,jtmax),',', Grd%dyt(itmax,jtmax),',', Thickness%dzt(itmax,jtmax,ktmax),' m)'
         write(unit,'(a,i6)') ' The number of cells in the column are kmt = ', Grd%kmt(itmax,jtmax)
         write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
              Thickness%dst(itmax,jtmax,ktmax),', ',Thickness%rho_dzt(itmax,jtmax,ktmax,tau)
      endif
  endif

  if (dztmin0 == dztmin) then
      write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
           ' The minimum dzt is ', dztmin,&
           ' metre at (i,j,k) = (',itmin+Dom%ioff,',',jtmin+Dom%joff,',',ktmin,&
           '),  (lon,lat,dpt) = (',Grd%xt(itmin,jtmin),',', &
           Grd%yt(itmin,jtmin),',', Thickness%depth_zt(itmin,jtmin,ktmin),' m)'
      if (use_blobs) then
         write(unit,'(a,es22.12,a,es22.12,a,es22.12,a)')      &
              ' The grid dimensions are (dxt, dyt, dztT) = (', &
              Grd%dxt(itmin,jtmin),',', Grd%dyt(itmin,jtmin),',', Thickness%dztT(itmin,jtmin,ktmin,tau),' m)'
         write(unit,'(a,i6)') ' And the number of cells in the column are kmt = ', Grd%kmt(itmin,jtmin)
         write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dztT) = ', &
              Thickness%dst(itmin,jtmin,ktmin),', ',Thickness%rho_dztT(itmin,jtmin,ktmin,tau)
      else
         write(unit,'(a,es22.12,a,es22.12,a,es22.12,a)')      &
              ' The grid dimensions are (dxt, dyt, dzt) = (', &
              Grd%dxt(itmin,jtmin),',', Grd%dyt(itmin,jtmin),',', Thickness%dzt(itmin,jtmin,ktmin),' m)'
         write(unit,'(a,i6)') ' And the number of cells in the column are kmt = ', Grd%kmt(itmin,jtmin)
         write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
              Thickness%dst(itmin,jtmin,ktmin),', ',Thickness%rho_dzt(itmin,jtmin,ktmin,tau)         
      endif

      if(dztmin < 0.0) then 
          write(unit,'(a)') '==>Error in ocean_tracer_util_mod (dzt_min_max): dztT < 0 detected.'
          call mpp_error(FATAL,&
               '==>Error in ocean_tracer_util_mod (dzt_min_max): dzt < 0 detected.')
      endif

  endif


  return


end subroutine dzt_min_max
! </SUBROUTINE>  NAME="dzt_min_max"


!#######################################################################
! <SUBROUTINE NAME="tracer_prog_chksum">
!
! <DESCRIPTION>
! Compute checksums for prognostic tracers 
! </DESCRIPTION>
subroutine tracer_prog_chksum(Time, Tracer, index, chksum)

  type(ocean_time_type),        intent(in)  :: Time
  type(ocean_prog_tracer_type), intent(in)  :: Tracer
  integer,                      intent(in)  :: index
  integer(i8_kind), optional, intent(inout) :: chksum
  
  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not. module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_tracer_util_mod (tracer_prog_chksum): module not initialized ')
  endif 

  write(stdoutunit,*) ' '
  write(stdoutunit,*) '=== Prognostic tracer checksum follows ==='
  call write_timestamp(Time%model_time)
  if (PRESENT(chksum)) then
     call write_chksum_3d(Tracer%name, Tracer%field(COMP,:,index)*Grd%tmask(COMP,:), chksum)
  else
     call write_chksum_3d(Tracer%name, Tracer%field(COMP,:,index)*Grd%tmask(COMP,:))
  endif

end subroutine tracer_prog_chksum
! </SUBROUTINE>  NAME="tracer_prog_chksum"


!#######################################################################
! <SUBROUTINE NAME="tracer_diag_chksum">
!
! <DESCRIPTION>
! Compute checksums for diagnostic tracers 
! </DESCRIPTION>
subroutine tracer_diag_chksum(Time, Tracer, chksum)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_diag_tracer_type), intent(in)    :: Tracer
  integer(i8_kind), optional,   intent(inout) :: chksum
  
  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not. module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_tracer_util_mod (tracer_diag_chksum): module not initialized ')
  endif 

  write(stdoutunit,*) '=== Diagnostic tracer checksum follows ==='
  call write_timestamp(Time%model_time)
  if (PRESENT(chksum)) then
     call write_chksum_3d(Tracer%name, Tracer%field(COMP,:)*Grd%tmask(COMP,:), chksum)
  else
     call write_chksum_3d(Tracer%name, Tracer%field(COMP,:)*Grd%tmask(COMP,:))
  endif

end subroutine tracer_diag_chksum
! </SUBROUTINE>  NAME="tracer_diag_chksum"


!#######################################################################
! <SUBROUTINE NAME="tracer_psom_chksum">
!
! <DESCRIPTION>
! Compute checksums for PSOM advection second order moments. 
! </DESCRIPTION>
subroutine tracer_psom_chksum(Time, Tracer)

  type(ocean_time_type),        intent(in)  :: Time
  type(ocean_prog_tracer_type), intent(in)  :: Tracer
  
  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not. module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_tracer_util_mod (tracer_psom_chksum): module not initialized ')
  endif 

  write(stdoutunit,*) ' '
  write (stdoutunit,*) 'Writing psom moments for tracer ',trim(Tracer%name)
  call write_timestamp(Time%model_time)

  call write_chksum_3d('Tracer%s0', Tracer%s0(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_3d('Tracer%sx', Tracer%sx(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_3d('Tracer%sxx', Tracer%sxx(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_3d('Tracer%sy', Tracer%sy(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_3d('Tracer%syy', Tracer%syy(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_3d('Tracer%sz', Tracer%sz(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_3d('Tracer%szz', Tracer%szz(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_3d('Tracer%sxy', Tracer%sxy(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_3d('Tracer%sxz', Tracer%sxy(COMP,:)*Grd%tmask(COMP,:))
  call write_chksum_3d('Tracer%syz', Tracer%syz(COMP,:)*Grd%tmask(COMP,:))

end subroutine tracer_psom_chksum
! </SUBROUTINE>  NAME="tracer_psom_chksum"



!#######################################################################
! <SUBROUTINE NAME="sort_pick_array">
!
! <DESCRIPTION>
! Simplest, and slowest, sorting algorithm from Numerical Recipes.
! Called "sort_pick" in Numerical Recipes.  
!
! Input are two arrays, first array defines the ascending sort 
! and second is a slave to the sort. 
!
! Typical example is sorting a vector of water parcels lightest
! to densest, with slave being volume of the parcels.  
!
! More sophisticated sorting algorithms exist, and may need to 
! be coded should this method prove too slow. 
!
! This scheme has order N^2 operations, which is a lot. 
!
! output has array(1) smallest and a(nsortpts) largest 
! with corresponding slave array.  
!
! coded Stephen.Griffies June 2005
!
! </DESCRIPTION>
!
subroutine sort_pick_array(array, slave)

  real, dimension(:), intent(inout) :: array
  real, dimension(:), intent(inout) :: slave

  real    :: tmp_a, tmp_s
  integer :: m, n, nsortpts

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_tracer_util_mod (sort_pick_array): module not initialized')
  endif 

  nsortpts = size(array)

  do n=2,nsortpts
      tmp_a = array(n)
      tmp_s = slave(n)
     do m=n-1,1,-1
        if(array(m) <= tmp_a) exit
        array(m+1)=array(m)
        slave(m+1)=slave(m)
     enddo
     array(m+1) = tmp_a
     slave(m+1) = tmp_s    
  enddo


end subroutine sort_pick_array
! </SUBROUTINE>  NAME="sort_pick_array"



!#######################################################################
! <SUBROUTINE NAME="sort_shell_array">
!
! <DESCRIPTION>
! Shell (or diminishing increment) sort from Numerical Recipes.
! Called "sort_shell" in Numerical Recipes.  
!
! Input are two arrays, first array defines the ascending sort 
! and second is a slave to the sort array. 
!
! Typical example is sorting a vector of water parcels lightest  
! to densest, with slave being volume of the parcels.  
!
! More sophisticated sorting algorithms exist, and may need to 
! be coded should this method prove too slow. 
!
! This scheme has order N^(5/4) operations. 
!
! output has array(1) smallest and a(nsortpts) largest,
! with corresponding ordering for slave array.  
!
! coded Stephen.Griffies June 2005
!
! </DESCRIPTION>
!
subroutine sort_shell_array(array, slave)

  real, dimension(:), intent(inout) :: array
  real, dimension(:), intent(inout) :: slave

  real    :: tmp_a, tmp_s
  integer :: m, n, inc, nsortpts

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_tracer_util_mod (sort_shell_array): module not initialized')
  endif 

  nsortpts = size(array)

  inc=1
  do 
     inc=3*inc+1
     if(inc > nsortpts) exit
  enddo

  do 
     inc=inc/3

     do m=inc+1,nsortpts
        tmp_a = array(m)
        tmp_s = slave(m)

        n=m
        do
           if(array(n-inc) <= tmp_a) exit
           array(n) = array(n-inc)
           slave(n) = slave(n-inc)
           n=n-inc
           if(n <= inc) exit
        enddo

        array(n)=tmp_a
        slave(n)=tmp_s

     enddo

     if(inc <=1) exit
  enddo

end subroutine sort_shell_array
! </SUBROUTINE>  NAME="sort_shell_array"


!#######################################################################
! <SUBROUTINE NAME="rebin_onto_rho">
!
! <DESCRIPTION>
! Bin a level input tendency field according to density classes. 
! The binning is meant for tendencies and transports, as used in particular
! for the neut_rho, wdian_rho, and tform_rho diagnostics.  
!
! Note that if use rebin_onto_rho_all_values=.false. then will 
! not be consistent with transport_on_nrho calculation, which includes
! bins all grid cells, including those outside of range for the bounds.  
!
! Stephen.Griffies 
! April 2012: algorithm made identical to transport_on_nrho as 
! computed in ocean_adv_vel_diag.
!
! </DESCRIPTION>
!
subroutine rebin_onto_rho (rho_bounds, rho_level, infield_level, outfield_rho)

  real, dimension(:),           intent(in)    :: rho_bounds
  real, dimension(isd:,jsd:,:), intent(in)    :: rho_level
  real, dimension(isd:,jsd:,:), intent(in)    :: infield_level
  real, dimension(isd:,jsd:,:), intent(inout) :: outfield_rho

  integer :: i, j, k, k_rho, rho_nk

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_neutral_util_mod (rebin_onto_rho): module needs initialization ')
  endif 

  rho_nk = size(rho_bounds(:)) - 1 
  outfield_rho(:,:,:) = 0.0

  ! rebin infield_level to get outfield_rho

  ! include even those cells with density outside of the 
  ! maximum and minimum of rho_bounds.
  if(rebin_onto_rho_all_values) then 

      do k_rho=1,rho_nk
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec

                  ! light waters 
                  if (k_rho == 1) then
                      if(rho_level(i,j,k) < rho_bounds(k_rho)) then 
                          outfield_rho(i,j,k_rho) = outfield_rho(i,j,k_rho) + infield_level(i,j,k)
                      endif

                  ! within the predefined range 
                  elseif(k_rho < rho_nk) then 
                      if( (rho_bounds(k_rho) <= rho_level(i,j,k)) .and.  &
                           (rho_level(i,j,k)  <  rho_bounds(k_rho+1)) ) then 
                          outfield_rho(i,j,k_rho) = outfield_rho(i,j,k_rho) + infield_level(i,j,k)
                      endif

                  ! denser than the last interface 
                  else   ! if (k_rho == rho_nk) then
                      if(rho_bounds(k_rho) <= rho_level(i,j,k)) then 
                          outfield_rho(i,j,k_rho) = outfield_rho(i,j,k_rho) + infield_level(i,j,k)             
                      endif
                  endif

               enddo
            enddo
         enddo
      enddo

  else 

      ! ignore those cells with density outside of the 
      ! maximum and minimum of rho_bounds.
      do k_rho=1,rho_nk
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  if( (rho_bounds(k_rho) <= rho_level(i,j,k)) .and.  &
                      (rho_level(i,j,k)  <  rho_bounds(k_rho+1)) ) then 
                      outfield_rho(i,j,k_rho) = outfield_rho(i,j,k_rho) + infield_level(i,j,k)
                  endif
               enddo
            enddo
         enddo
      enddo

  endif
 
  ! apply land-sea masks for level-mask at k=1 
  do k_rho=1,rho_nk
     do j=jsc,jec
        do i=isc,iec
           outfield_rho(i,j,k_rho) = outfield_rho(i,j,k_rho)*Grd%tmask(i,j,1)
        enddo
     enddo
  enddo

end subroutine rebin_onto_rho
! </SUBROUTINE> NAME="rebin_onto_rho"



!#######################################################################
! <SUBROUTINE NAME="diagnose_3d_rho">
!
! <DESCRIPTION>
! Helper function for diagnosting 3D data mapped onto density levels.
! </DESCRIPTION>
!
subroutine diagnose_3d_rho(Time, Dens, id_name, data)
  type(ocean_time_type), intent(in) :: Time
  type(ocean_density_type), intent(in) :: Dens
  integer, intent(in) :: id_name
  real, dimension(isd:,jsd:,:), intent(in) :: data

  real :: nrho_work(isd:ied,jsd:jed,size(Dens%neutralrho_ref(:)))
  logical :: used

  integer :: neutralrho_nk

  if (id_name > 0) then
     neutralrho_nk = size(Dens%neutralrho_ref(:))
     nrho_work(:,:,:) = 0.0
     call rebin_onto_rho(Dens%neutralrho_bounds, Dens%neutralrho, data, nrho_work)
     used = send_data (id_name, nrho_work(:,:,:), &
          Time%model_time,                        &
          is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=neutralrho_nk)
  endif

end subroutine diagnose_3d_rho
! </SUBROUTINE> NAME="diagnose_3d_rho"


!#######################################################################
! <SUBROUTINE NAME="diagnose_mass_of_layer_orig">
! <DESCRIPTION>
!
! Diagnose the mass of a layer as a function of layer and (i,j). 
!
! Method uses linear interpolation to find the mass per area of
! layer boundaries, which are then used to get mass per area of
! layer, and then mass through multipication by area.  
!
! Scheme currently does not forward (backwards) interpolate if 
! layer boundary lies within lowest (uppermost) grid cell.
!
! Diagnostic only makes sense when layer_level is monotonically
! increasing as go deeper in water column. 
!
! layer_bounds = neutral density of layer interfaces 
! mass_level   = mass per area at bottom of tracer cell
! layer_level  = neutral density of model grid levels   
! mass_layer   = diagnosed mass of a neutral density layer
!
! </DESCRIPTION>
!
subroutine diagnose_mass_of_layer_orig(layer_bounds, mass_level, layer_level, mass_layer)

  real, dimension(:),           intent(in)    :: layer_bounds
  real, dimension(isd:,jsd:,:), intent(in)    :: mass_level
  real, dimension(isd:,jsd:,:), intent(in)    :: layer_level
  real, dimension(isd:,jsd:,:), intent(inout) :: mass_layer

  real      :: w1,w2
  integer   :: layer_nk
  integer   :: i, j, k, n

  if (.not.module_is_initialized) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_tracer_util (diagnose_mass_of_layer): module needs initialization ')
  endif

  ! find mass per area at the layer bounds 
  layer_nk = size(layer_bounds(:)) - 1 
  mass_layer(:,:,:) = 0.0
  do n=1,layer_nk
     do j=jsc,jec
        do i=isc,iec
kloop:     do k=nk-1,1,-1
              if(    layer_level(i,j,k) <= layer_bounds(n)) then
                  if(layer_bounds(n)    < layer_level(i,j,k+1)) then
                      if(Grd%tmask(i,j,k+1) > 0) then 
                          W1= layer_bounds(n)      - layer_level(i,j,k)
                          W2= layer_level(i,j,k+1) - layer_bounds(n)
                          mass_layer(i,j,n) = ( mass_level(i,j,k+1)*W1  &
                                               +mass_level(i,j,k)  *W2) &
                                               /(W1 + W2 + epsln)
                          exit kloop 
                      endif
                  endif
              endif
           enddo kloop
        enddo
     enddo
  enddo

  ! take differences of mass per area at layer boundaries to then 
  ! get the mass within a layer. 
  ! if mass < 0, then this means we have reached edge of the density classes, 
  ! so set mass to zero. 
  do n=1,layer_nk-1
     do j=jsc,jec
        do i=isc,iec
            if(mass_layer(i,j,n+1) >= mass_layer(i,j,n)) then 
               mass_layer(i,j,n) = Grd%dat(i,j)*(mass_layer(i,j,n+1)-mass_layer(i,j,n))
            else
               mass_layer(i,j,n) = 0.0
            endif 
        enddo
     enddo
  enddo

  ! bottom layer mass equated to that above (cludge)
  n=layer_nk
  do j=jsc,jec
     do i=isc,iec
         mass_layer(i,j,n) = Grd%dat(i,j)*mass_layer(i,j,n-1)
     enddo
  enddo


end subroutine diagnose_mass_of_layer_orig
! </SUBROUTINE>  NAME="diagnose_mass_of_layer_orig"




!#######################################################################
! <SUBROUTINE NAME="diagnose_mass_of_layer">
! <DESCRIPTION>
!
! Diagnose the mass of a layer as a function of layer and (i,j). 
!
! Method uses linear interpolation to find the mass per area of
! layer boundaries, which are then used to get mass per area of
! layer, and then mass through multipication by area.  
!
! Diagnostic only makes sense when layer_level is monotonically
! increasing as go deeper in water column. 
!
! </DESCRIPTION>
!
subroutine diagnose_mass_of_layer(area_t, dzt, dztlo, dztup, rho_dzt, &
                                  nrho_level, nrho_nk, nrho_bdy, nrho_layer_mass)

  real, dimension(isd:,jsd:),   intent(in)    :: area_t 
  real, dimension(isd:,jsd:,:), intent(in)    :: dzt 
  real, dimension(isd:,jsd:,:), intent(in)    :: dztlo 
  real, dimension(isd:,jsd:,:), intent(in)    :: dztup
  real, dimension(isd:,jsd:,:), intent(in)    :: rho_dzt 
  real, dimension(isd:,jsd:,:), intent(in)    :: nrho_level
  integer,                      intent(in)    :: nrho_nk
  real, dimension(:),           intent(in)    :: nrho_bdy
  real, dimension(isd:,jsd:,:), intent(inout) :: nrho_layer_mass

  real    :: W1, W2
  real    :: tmp_mass, relative_mass_diff
  real, allocatable, dimension(:) :: nrho_bdy_mass
  integer :: i, j, k, n, kbot
  integer :: nrho_bdy_nk
  integer :: num_nonzero_bdy
  
  if (.not.module_is_initialized) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_tracer_util (diagnose_mass_of_layer): module needs initialization ')
  endif
 
  nrho_bdy_nk = nrho_nk + 1 
  allocate(nrho_bdy_mass(nrho_nk+1))
  nrho_bdy_mass(:) = 0.0

  ! compute mass per area for levels 
  wrk1(:,:,:)  = 0.0  ! mass per area sitting above a T-cell center 
  wrk1_2d(:,:) = 0.0  ! mass per area in full T-cell column 
  do j=jsc,jec
     do i=isc,iec

        kbot = Grd%kmt(i,j)
        if(kbot > 1) then 

            k=1 
            wrk1(i,j,k) = rho_dzt(i,j,k)*dztup(i,j,k)/dzt(i,j,k)
            do k=2,kbot
               wrk1(i,j,k) = wrk1(i,j,k-1) + rho_dzt(i,j,k-1)*dztlo(i,j,k-1)/dzt(i,j,k-1) &
                                           + rho_dzt(i,j,k)  *dztup(i,j,k)  /dzt(i,j,k)
            enddo
            k=kbot
            wrk1_2d(i,j)  = wrk1(i,j,k) + rho_dzt(i,j,k)*dztlo(i,j,k)/dzt(i,j,k)

        endif

     enddo
  enddo


  ! determine mass within neutral rho layers 
  wrk2_2d(:,:) = 0.0  ! number of nonzero mass layers 
  do j=jsc,jec
     do i=isc,iec

        kbot = Grd%kmt(i,j)
        if(kbot > 1) then 

            wrk2_2d(i,j)     = 0.0
            num_nonzero_bdy  = 0

            ! loop to determine mass per area sitting above depth of layer interface 
            do n=1,nrho_bdy_nk
               nrho_bdy_mass(n) = 0.0

               ! linear interpolate for mass sitting above interior layer interfaces 
               do k=2,kbot
                  if(nrho_level(i,j,k) > nrho_bdy(n) .and. nrho_level(i,j,k-1) <= nrho_bdy(n)) then 
                      W1 = nrho_level(i,j,k) - nrho_bdy(n)
                      W2 = nrho_bdy(n)       - nrho_level(i,j,k-1) 
                      nrho_bdy_mass(n) = wrk1(i,j,k-1) + (                                                          &
                        W1*rho_dzt(i,j,k-1)*dztlo(i,j,k-1)/dzt(i,j,k-1) + W2*rho_dzt(i,j,k)*dztup(i,j,k)/dzt(i,j,k) &
                        )/(W1 + W2 + epsln)
                      wrk2_2d(i,j)     = wrk2_2d(i,j) + 1.0
                      num_nonzero_bdy  = num_nonzero_bdy + 1
                  endif
               enddo

               ! interfaces heavier than nrho_level(k=kbot) assumed to 
               ! have full column mass sitting above the interfaces.  
               k=kbot
               if(nrho_level(i,j,k) < nrho_bdy(n)) then 
                   nrho_bdy_mass(n) = wrk1_2d(i,j)
                   wrk2_2d(i,j)     = wrk2_2d(i,j) + 1.0 
                   num_nonzero_bdy  = num_nonzero_bdy + 1
               endif

            enddo

            ! determine mass within each neutral density layer. 
            ! linear interpolation can create negatives, so need to 
            ! apply an abs. The abs will cause some nonconsevation of mass,
            ! but hopefully will be negligible for purposes of the diagnostic. 
            do n=1,nrho_nk
               nrho_layer_mass(i,j,n) = area_t(i,j)*abs(nrho_bdy_mass(n+1)-nrho_bdy_mass(n))*Grd%tmask(i,j,1)
            enddo                  
                   
            ! if all nrho_bdy_mass containt full mass, we place 
            ! all mass in lightest of the nrho_layer_mass in order
            ! to conserve the full column mass. 
            if(num_nonzero_bdy == nrho_bdy_nk) then 
               nrho_layer_mass(i,j,1) = area_t(i,j)*wrk1_2d(i,j)*Grd%tmask(i,j,1)
            endif 


        endif ! endif for kbot > 1 

     enddo   ! enddo for i
  enddo   ! enddo for j


  if(debug_diagnose_mass_of_layer) then 
      wrk3_2d(:,:) = 0.0

      do j=jsc,jec
         do i=isc,iec
            relative_mass_diff = 0.0
            kbot = Grd%kmt(i,j)
            if(kbot > 0) then 
                do n=1,nrho_nk
                   wrk3_2d(i,j) = wrk3_2d(i,j) + nrho_layer_mass(i,j,n)  
                enddo
                tmp_mass           = area_t(i,j)*wrk1_2d(i,j)
                relative_mass_diff = (tmp_mass-wrk3_2d(i,j))/tmp_mass
            endif
            if(abs(relative_mass_diff) > epsln_diagnose_mass_of_layer) then 
                write(unit,'(/a,i5,a,i5,a,e20.10,a,e20.10,a,e20.10)')                               &
                'relative mass difference and masses at (i,j) = (',i+Dom%ioff,',',j+Dom%joff,') = ',&
                relative_mass_diff,',',tmp_mass,',',wrk3_2d(i,j)
                write(unit,'(a,f12.1)') &
                'number of nonzero interface masses in this column are ',wrk2_2d(i,j)  
            endif
            do n=1,nrho_nk
               if(nrho_layer_mass(i,j,n) < 0.0) then 
                   write(unit,'(/a,i5,a,i5,a,i5,e20.10)')                                      &
                   'negative layer mass at (i,j,n) = (',i+Dom%ioff,',',j+Dom%joff,',',n,') = ',&
                   nrho_layer_mass(i,j,n) 
               endif
            enddo

         enddo
      enddo
  endif


end subroutine diagnose_mass_of_layer
! </SUBROUTINE>  NAME="diagnose_mass_of_layer"


end module ocean_tracer_util_mod


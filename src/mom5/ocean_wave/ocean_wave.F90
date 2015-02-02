module ocean_wave_mod
#define COMP isc:iec,jsc:jec
!
! <CONTACT EMAIL="Martin.Schmidt@io-warnemuende.de"> M. Schmidt
! </CONTACT>
!
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
! </REVIEWER>
!
!<OVERVIEW>
! Idealized wave model for delivering wave number and wave height
! for coupled current-wave action on sediment. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This idealized wave model delivers wave number and wave height 
! for calculation of the coupled current-wave action on sediment.
! Swell is not included in this model.
!
! All fields are defined at tracer grid, for later use in sediment dynamics
! in such modules as ocean_shared/generic_tracers/generic_ERGOM.F90
! </DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! P.C. Liu, D.J. Schwab and J.R. Bennett, 
!      Journ. of Phys. Oceanography 14, 1514 (1984).
!      D.J. Schwab, J.R. Bennett, P.C. Liu and M.A. Donelan,
!      Journ. of Geophysical Research  89 (C3), 3586 (1984).  
!</REFERENCE>
!
!<REFERENCE>
! Hughes, S. A. 1984. "TMA Shallow-Water Spectrum:
!      Description and Application," Technical Report CERC-84-7, 
!      US Army Engineer Waterways Experiment Station, Vicksburg, Miss.
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST>
!</NAMELIST>

use constants_mod,          only: pi, grav
use diag_manager_mod,       only: register_diag_field
use fms_mod,                only: open_namelist_file, check_nml_error, close_file, file_exist
use fms_mod,                only: write_version_number, FATAL, WARNING, NOTE
use fms_mod,                only: clock_flag_default
use mpp_mod,                only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
use fms_io_mod,             only: reset_field_pointer, restart_file_type
use fms_io_mod,             only: register_restart_field, save_restart, restore_state
use mpp_mod,                only: input_nml_file, mpp_error, stdout, stdlog
use mpp_mod,                only: mpp_max, mpp_min, mpp_pe
use mpp_mod,                only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE, CLOCK_ROUTINE
use ocean_domains_mod,      only: get_local_indices, get_global_indices
use mpp_domains_mod,        only: mpp_update_domains
use ocean_types_mod,        only: ocean_domain_type, ocean_grid_type, ocean_options_type
use ocean_types_mod,        only: ocean_time_type, ocean_time_steps_type
use ocean_parameters_mod,   only: missing_value
use ocean_workspace_mod,    only: wrk1_2d, wrk2_2d, wrk3_2d, wrk4_2d
use ocean_util_mod,         only: write_timestamp, diagnose_2d, write_chksum_2d
use ocean_types_mod,        only: ice_ocean_boundary_type
use wave_types_mod,         only: ocean_wave_type
use data_override_mod,      only: data_override

implicit none
private

public ocean_wave_init
public ocean_wave_model
public ocean_wave_end
public wave_model_is_initialised

type(ocean_domain_type), pointer   :: Dom =>NULL()
type(ocean_grid_type), pointer     :: Grd =>NULL()

character(len=128) :: &
     version='$Id: ocean_wave.F90,v 20.0 2013/12/14 00:17:26 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

real,parameter:: gtpi=grav/2.0/pi, epsln=1e-20 
real,parameter:: sqrt_2 = sqrt(2.)
real,parameter:: rhoair=1.2233, rhoh2o=1000., gamma=0.028  
real,parameter:: fac2=0.5*gamma*rhoair/rhoh2o/grav
real,parameter:: twopi=2.*pi, third=1./3., seventh=1./7.
real,parameter:: rho_ice=905.

#include <ocean_memory.h>

#ifndef MOM_STATIC_ARRAYS
real, dimension(:,:),  allocatable:: windx, windy      !wind fields for surface stress
real, dimension(:,:),  allocatable:: sn,cs,c,s         !direction fields for wave propagation
real, dimension(:,:),  allocatable:: wrk1,  wrk2       !work space in the compute domain
#else
real, dimension(isd:ied,jsd:jed)  :: windx, windy      !wind fields for surface stress
real, dimension(isd:ied,jsd:jed)  :: sn,cs,c,s         !direction fields for wave propagation
real, dimension(isc:iec,jsc:jec)  :: wrk1, wrk2        !work space in the compute domain
#endif
real     :: wavedamp = -10.
logical  :: module_is_initialized = .FALSE.
logical  :: first_call            = .TRUE.
logical  :: damp_where_ice        = .TRUE.
logical  :: filter_wave_mom       = .TRUE.
logical  :: write_a_restart       = .TRUE.
logical  :: use_this_module       = .FALSE.
logical  :: use_TMA               = .TRUE.
logical  :: debug_this_module     = .FALSE. ! for debugging--prints out a lot of checksums 
integer  :: tau_w, taup1_w
real     :: cgmax, gridmin, dttw, dtts, wave_damp

integer  :: id_windx  =-1 
integer  :: id_windy  =-1
integer  :: id_xmom   =-1
integer  :: id_ymom   =-1
integer  :: id_wave_k =-1
integer  :: id_height =-1
integer  :: id_wave_p =-1
integer  :: id_init
integer  :: id_diag

! for write statements 
integer :: stdoutunit,stdlogunit 


! for restart
integer                       :: id_restart(2) = 0
type(restart_file_type), save :: wave_restart



namelist /ocean_wave_nml/ wavedamp, damp_where_ice, write_a_restart, debug_this_module &
                          ,use_TMA, filter_wave_mom, use_this_module

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_wave_init">
!
! <DESCRIPTION>
! Initialize the wave module
! </DESCRIPTION>
!

subroutine ocean_wave_init(Grid, Domain, Waves, Time, Time_steps, Ocean_options, debug)
  type(ocean_grid_type),          intent(in), target   :: Grid
  type(ocean_domain_type),        intent(in), target   :: Domain
  type(ocean_wave_type),          intent(inout)        :: Waves
  type(ocean_time_type),          intent(in)           :: Time
  type(ocean_time_steps_type),    intent(inout)        :: Time_steps 
  type(ocean_options_type),       intent(inout)        :: Ocean_options 
  logical,                        intent(in), optional :: debug
  
  real    :: gridsp
  integer :: ioun, io_status, ierr
  integer :: i, j

  stdoutunit=stdout();stdlogunit=stdlog() 
  
  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_wave_mod (ocean_wave_init): module already initialized.')
  endif 
  
  id_init = mpp_clock_id('(Wave model initialization) ',grain=CLOCK_ROUTINE)
  id_diag = mpp_clock_id('(Wave diagnostics) '         ,grain=CLOCK_ROUTINE)
  
  call mpp_clock_begin(id_init)

  call write_version_number( version, tagname )
  
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, ocean_wave_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_wave_nml')
#else
  ioun = open_namelist_file()
  read  (ioun, ocean_wave_nml,iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_wave_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_wave_nml)  
  write (stdlogunit, ocean_wave_nml)


  if(use_this_module) then
    write(stdoutunit,'(a)') '==>Note: Using the MOM simple ocean surface wave model'  
    module_is_initialized = .TRUE.
    Ocean_options%ocean_ideal_surf_wave = 'Used idealized ocean surface wave module.'
  else
    write(stdoutunit,'(a)') '==>Note: Not using the idealized ocean surface wave module.'  
    module_is_initialized = .FALSE.
    Ocean_options%ocean_ideal_surf_wave = 'Not using the idealized ocean surface wave module.'
    return
  endif

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_wave with debug_this_module=.true.'  
  endif 
  if(use_TMA) then 
    write(stdoutunit,'(a)') '==>Note: using TMA-spectrum for shallow areas'  
  else
    write(stdoutunit,'(a)') '==>Note: not using TMA-spectrum for shallow areas'  
  endif 

  if(.not. write_a_restart) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_wave with write_a_restart=.false.'
    write(stdoutunit,'(a)') '   Will NOT write restart file, and so cannot restart the run.'
  endif 

  Dom => Domain
  Grd => Grid

#ifndef MOM_STATIC_ARRAYS

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  call get_global_indices(Domain, isg, ieg, jsg, jeg)

  allocate(Waves%xmom(isd:ied,jsd:jed,0:1))
  allocate(Waves%ymom(isd:ied,jsd:jed,0:1))
  allocate(Waves%wave_k(isd:ied,jsd:jed))
  allocate(Waves%height(isd:ied,jsd:jed))
  allocate(Waves%wave_p(isd:ied,jsd:jed))

  allocate(windx(isd:ied,jsd:jed))
  allocate(windy(isd:ied,jsd:jed))
  allocate(sn(isd:ied,jsd:jed))
  allocate(cs(isd:ied,jsd:jed))
  allocate(s(isd:ied,jsd:jed))
  allocate(c(isd:ied,jsd:jed))
  allocate(wrk1(isc:iec,jsc:jec))
  allocate(wrk2(isc:iec,jsc:jec))
#endif

  Waves%xmom   = 0.0
  Waves%ymom   = 0.0
  Waves%wave_k = 0.0
  Waves%height = 0.0
  Waves%wave_p = 0.0
  windx        = 0.0
  windy        = 0.0
  sn           = 0.0
  cs           = 0.0
  s            = 0.0
  c            = 0.0
  wrk1         = 0.0 
  wrk2         = 0.0 


  id_windx  = register_diag_field ('ocean_model', 'windx_wave', Grd%tracer_axes(1:2), &
              Time%model_time, 'zonal wind speed for waves', 'm/s',                   &
              missing_value=missing_value, range=(/-1e3,1e3/))
  id_windy  = register_diag_field ('ocean_model', 'windy_wave', Grd%tracer_axes(1:2), &
              Time%model_time, 'meridional wind speed for waves', 'm/s',              &
              missing_value=missing_value, range=(/-1e3,1e3/))
  id_xmom   = register_diag_field ('ocean_model', 'xmom', Grd%tracer_axes(1:2), &
              Time%model_time, 'zonal wave momentum', 'm^2/s',                  &
              missing_value=missing_value, range=(/-1e3,1e3/))
  id_ymom   = register_diag_field ('ocean_model', 'ymom', Grd%tracer_axes(1:2), &
              Time%model_time, 'meridional wave momentum', 'm^2/s',             &
              missing_value=missing_value, range=(/-1e3,1e3/))
  id_height = register_diag_field ('ocean_model', 'height', Grd%tracer_axes(1:2), &
              Time%model_time, 'significant wave height', 'm',                    &
              missing_value=missing_value, range=(/-1e3,1e3/))
  id_wave_p = register_diag_field ('ocean_model', 'wave_p', Grd%tracer_axes(1:2), &
              Time%model_time, 'peak frequency', '1/s',                           &
              missing_value=missing_value, range=(/-1e3,1e3/))
  id_wave_k = register_diag_field ('ocean_model', 'wave_k', Grd%tracer_axes(1:2), &
              Time%model_time, 'wave number', '1/m',                              &
              missing_value=missing_value, range=(/-1e3,1e3/))

  ! initialise the time stepping
  tau_w   = 0
  taup1_w = 1
  dtts     = Time_steps%dtts

! find the minimum and maximum grid spacing to organise substepping
  gridmin          = 1.0e20
  do j=jsc,jec
     do i=isc,iec
        if (Grd%kmt(i,j) > 0) then
           gridsp = 1.0/(Grd%dxu(i,j)*Grd%dxu(i,j)) + 1.0/(Grd%dyu(i,j)*Grd%dyu(i,j))
           gridsp = sqrt(1.0/gridsp) 
           gridmin = min(gridsp,gridmin)
        endif
     enddo
  enddo
  call mpp_min (gridmin)
! readthe initial field or the restart file
  call read_wave(Waves)

  wave_damp = wavedamp/rho_ice  

  call mpp_clock_end(id_init) 

  return
end subroutine ocean_wave_init
! </SUBROUTINE> NAME="ocean_wave_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_wave_model">
!
! <DESCRIPTION>
! time step the wave model
! </DESCRIPTION>
!
subroutine ocean_wave_model(Time, Waves, Ice_ocean_boundary)
  type(ocean_time_type),          intent(in)    :: Time 
  type(ocean_wave_type),          intent(inout) :: Waves
  type(ice_ocean_boundary_type),  intent(in)    :: Ice_ocean_boundary

  integer :: i,j
  integer :: ndtt, nww
  real    :: wmax, cspeed, dtwmax

  if ( .not.module_is_initialized ) return  

  if(debug_this_module) then 
     write(stdoutunit,*) 'starting ocean_wave_model'
  endif 

  ! If we do not initialize them they will have junk in haloes (huge values) after data_overrides below.
  ! The subsequent mpp_update_domains calls are supposed to update these haloes with sensible values.
  ! But, that does not happen in the Baltic testcase and the haloes will contain junk even after mpp_update_domains calls!
  ! This would cause crash later when those huge values in the haloes are accessed for doing arithmatic 
  ! (e.g., in use_TMA block below).
  ! Zeroing out windx and windy on data domain at the onset avoids such crash. 
  ! But the root cause is not clear and there may be a bug somewhere. 
  ! Niki.Zadeh 
  windx = 0.0
  windy = 0.0

  call data_override('OCN', 'u_bot', wrk1, Time%model_time )
  call data_override('OCN', 'v_bot', wrk2, Time%model_time )

  do j=jsc,jec
     do i=isc,iec
       windx(i,j) = wrk1(i,j)
       windy(i,j) = wrk2(i,j)
     enddo
  enddo

  call mpp_update_domains (windx, Dom%domain2d)
  call mpp_update_domains (windy, Dom%domain2d)

  ! wavediag is needed to initialize wave_p      
  if (first_call) then
    call ocean_wave_diag(Waves) 
    first_call=.false.
    if(debug_this_module) write(stdoutunit,*) 'first start of ocean_wave_model'
  endif
! find the maximum speed for stability of time stepping
  wmax=epsln
  do j=jsc,jec
     do i=isc,iec
        if (Grd%kmt(i,j) > 0) then
           wmax=max(wmax, abs(windx(i,j)), abs(windy(i,j)))
!phase speed
           cspeed=gtpi/(Waves%wave_p(i,j)+epsln)
           wmax=max(wmax,abs(cspeed))        
        endif
     enddo
  enddo
  call mpp_max(wmax )
  dtwmax=0.5*gridmin/sqrt_2/wmax     
  ndtt=int(dtts/dtwmax)+1   
  
  if(debug_this_module) write(stdoutunit,*) 'ocean_wave_model: wmax    = ', wmax
  if(debug_this_module) write(stdoutunit,*) 'ocean_wave_model: gridmin = ', gridmin
  if(debug_this_module) write(stdoutunit,*) 'ocean_wave_model: ndtt    = ', ndtt
  
  dttw=dtts/float(ndtt)
!loop for waves       
  do nww=1,ndtt      
    tau_w=abs(tau_w-1)
    taup1_w=abs(taup1_w-1)  
    call ocean_wave_diag(Waves)  
    call ocean_wave_prop(Waves, Ice_ocean_boundary)  
!    call ocean_wave_prop(Time)  
  enddo      
! filtering xmom and ymom       
!  if (mod(itt,10).eq.0) then
  if(filter_wave_mom) call ocean_wave_filter(Waves, dtts)
!  endif

  call diagnose_2d(Time, Grd, id_windx, windx(:,:))
  call diagnose_2d(Time, Grd, id_windy, windy(:,:))
  if (id_xmom  > 0) call diagnose_2d(Time, Grd, id_xmom, Waves%xmom(:,:,taup1_w)*grav)
  if (id_ymom  > 0) call diagnose_2d(Time, Grd, id_ymom, Waves%ymom(:,:,taup1_w)*grav)
  call diagnose_2d(Time, Grd, id_height, Waves%height(:,:))
  call diagnose_2d(Time, Grd, id_wave_p, Waves%wave_p(:,:))
  call diagnose_2d(Time, Grd, id_wave_k, Waves%wave_k(:,:))

  if(debug_this_module) write(stdoutunit,*) 'ending ocean_wave_model'
    
  return
end subroutine ocean_wave_model
! </SUBROUTINE> NAME="ocean_wave_model"


!#######################################################################
! <SUBROUTINE NAME="ocean_wave_prop">
!
! <DESCRIPTION>
! wave propagation
! </DESCRIPTION>
!
subroutine ocean_wave_prop(Waves, Ice_ocean_boundary)
  type(ocean_wave_type),          intent(inout) :: Waves
  type(ice_ocean_boundary_type),  intent(in)    :: Ice_ocean_boundary

  integer :: i,j
  real    :: cm, adv_Txy, adv_Txx, adv_Tyy, adv_Tyx, SIP, SIM, SJP, SJM
  real    :: w1, w2, wspd, wspdr, tmp, cthe, sg, damp, diff

  if(debug_this_module) then 
     write(stdoutunit,*) 'ocean_wave_model: start wave_prop'
  endif 

  do j=jsd,jed
    do i=isd,ied      
      cm        = sqrt(Waves%xmom(i,j,tau_w)**2+Waves%ymom(i,j,tau_w)**2+epsln)   ! m^2/s
      cs(i,j)   = Waves%xmom(i,j,tau_w)/cm
      sn(i,j)   = Waves%ymom(i,j,tau_w)/cm  
      tmp       = gtpi/(Waves%wave_p(i,j)+epsln)
      c(i,j)    = tmp*Grd%tmask(i,j,1)                                ! m/s
      tmp       = tmp*cm    
      wrk1_2d(i,j) = tmp*Grd%tmask(i,j,1)
      s(i,j)    = sqrt(tmp)*Grd%tmask(i,j,1)
    enddo
  enddo

  do j=jsc,jec
    do i=isc,iec      
      adv_Txy = 0.5*( wrk1_2d(i,j)  *cs(i,j)  *(sn(i,j)  +abs(sn(i,j  )))  &
           + wrk1_2d(i,j+1)*cs(i,j+1)*(sn(i,j+1)-abs(sn(i,j+1)))) &
           -0.5*( wrk1_2d(i,j-1)*cs(i,j-1)*(sn(i,j-1)+abs(sn(i,j-1)))  &
           + wrk1_2d(i,j)  *cs(i,j)  *(sn(i,j)  -abs(sn(i,j  ))))

      adv_Txx = 0.5*( wrk1_2d(i,j)  *cs(i,j)  *(cs(i,j)  +abs(cs(i,j  )))  &
           + wrk1_2d(i+1,j)*cs(i+1,j)*(cs(i+1,j)-abs(cs(i+1,j)))) &
           -0.5*( wrk1_2d(i-1,j)*cs(i-1,j)*(cs(i-1,j)+abs(cs(i-1,j)))  &
           + wrk1_2d(i,j)  *cs(i,j)  *(cs(i,j)  -abs(cs(i,j  ))))
            
      SIP = Grd%tmask(i+1,j,1)*s(i+1,j)+(1.0-Grd%tmask(i+1,j,1))*&
            (2.0*s(i,j)-s(i-1,j))
      SIM = Grd%tmask(i-1,j,1)*s(i-1,j)+(1.0-Grd%tmask(i-1,j,1))*&
            (2.0*s(i,j)-s(i+1,j))            
      Waves%xmom(i,j,taup1_w) = Waves%xmom(i,j,tau_w) &
                          - dttw*( 0.25*adv_Txx*Grd%dxtr(i,j)+0.25*adv_Txy*Grd%dytr(i,j) &
                          + (SIP-SIM)*(SIP+SIM)/16.*Grd%dxtr(i,j))*Grd%tmask(i,j,1)
    enddo
  enddo

  do j=jsc,jec
    do i=isc,iec      
      adv_Tyy = 0.5*( wrk1_2d(i,j)  *sn(i,j)  *(sn(i,j)  +abs(sn(i,j  )))  &
           +wrk1_2d(i,j+1)*sn(i,j+1)*(sn(i,j+1)-abs(sn(i,j+1)))) &
           -0.5*( wrk1_2d(i,j-1)*sn(i,j-1)*(sn(i,j-1)+abs(sn(i,j-1)))  &
           +wrk1_2d(i,j)  *sn(i,j)  *(sn(i,j)  -abs(sn(i,j  ))))
      adv_Tyx = 0.5*( wrk1_2d(i,j)  *sn(i,j)  *(cs(i,j)  +abs(cs(i,j  )))  &
           + wrk1_2d(i+1,j)*sn(i+1,j)*(cs(i+1,j)-abs(cs(i+1,j)))) &
           -0.5*( wrk1_2d(i-1,j)*sn(i-1,j)*(cs(i-1,j)+abs(cs(i-1,j)))  &
           + wrk1_2d(i,j)  *sn(i,j)  *(cs(i,j)  -abs(cs(i,j  )))) 
            
      SJP = Grd%tmask(i,j+1,1)*s(i,j+1)+(1.0-Grd%tmask(i,j+1,1))*&
            (2.0*s(i,j)-s(i,j-1))
      SJM = Grd%tmask(i,j-1,1)*s(i,j-1)+(1.0-Grd%tmask(i,j-1,1))*&
            (2.0*s(i,j)-s(i,j+1))                
      Waves%ymom(i,j,taup1_w) = Waves%ymom(i,j,tau_w)       &
           - dttw*( 0.25*adv_Tyy*Grd%dytr(i,j)+0.25*adv_Tyx*Grd%dxtr(i,j)         &
           + (SJP-SJM)*(SJP+SJM)/16.*Grd%dytr(i,j))*Grd%tmask(i,j,1)
    enddo
  enddo    

  !  calculate the wind input to momentum flux for two directions:
  ! the wind direction (w1) and the wave field direction (w2)
  ! drag coefficient depends on direction, maximal for w2, reduced for w1
  do j=jsc,jec
    do i=isc,iec      
      wspd   = sqrt(windx(i,j)**2+windy(i,j)**2+epsln)
      wspdr  = 1/wspd
      cthe   = (windx(i,j)*cs(i,j) + windy(i,j)*sn(i,j))*wspdr

      sg     = max(5.0E-03,s(i,j) * abs(cthe))
      tmp    = wspd-0.83*c(i,j)*cthe
      w1     = tmp*abs(tmp)*wspdr*(0.4/log(50./sg))**2

      sg     = max(5.0E-06,s(i,j))
      tmp    = cthe*wspd-0.83*c(i,j)
      w2     = tmp*abs(tmp)*(0.4/log(50./sg))**2
      wrk3_2d(i,j) = dttw*(w1*windx(i,j) + w2*cs(i,j))*fac2 *Grd%tmask(i,j,1) 
      wrk4_2d(i,j) = dttw*(w1*windy(i,j) + w2*sn(i,j))*fac2 *Grd%tmask(i,j,1)
    enddo
  enddo    

  if (damp_where_ice) then
! Damp the time tendency
     do j=jsc,jec
        do i=isc,iec      
           if (Ice_ocean_boundary%mi(i,j) > 0) then
              damp = exp(wave_damp*Ice_ocean_boundary%mi(i,j)) 
              wrk3_2d(i,j) = wrk3_2d(i,j) * damp
              wrk4_2d(i,j) = wrk4_2d(i,j) * damp
           endif 
        enddo
     enddo 
  endif
   
  do j=jsc,jec
    do i=isc,iec      
      Waves%xmom(i,j,taup1_w) = Waves%xmom(i,j,taup1_w) + wrk3_2d(i,j)
      Waves%ymom(i,j,taup1_w) = Waves%ymom(i,j,taup1_w) + wrk4_2d(i,j)
    enddo
  enddo 

  if (damp_where_ice) then
! This is engeneering needs improvement with resolved ice classes ...
    do j=jsc,jec
      do i=isc,iec
        diff = Ice_ocean_boundary%mi(i,j)-Waves%height(i,j)*rho_ice
        if(diff >= 0.) then
           damp = exp(wave_damp*diff)
           Waves%xmom(i,j,taup1_w)=Waves%xmom(i,j,taup1_w)*damp   
           Waves%ymom(i,j,taup1_w)=Waves%ymom(i,j,taup1_w)*damp   
        endif
      enddo
    enddo
  endif

  call mpp_update_domains (Waves%xmom(:,:,taup1_w), Dom%domain2d)
  call mpp_update_domains (Waves%ymom(:,:,taup1_w), Dom%domain2d)

  if(debug_this_module) then 
     write(stdoutunit,*) 'ocean_wave_model: end wave_prop'
  endif 

  return

end subroutine ocean_wave_prop
! </SUBROUTINE> NAME="ocean_wave_prop"


!#######################################################################
! <SUBROUTINE NAME="ocean_wave_diag">
!
! <DESCRIPTION>
! wave diagnostics
! </DESCRIPTION>
!
subroutine ocean_wave_diag(Waves)
  type(ocean_wave_type), intent(inout) :: Waves

  real, parameter   :: const1=0.01788735, const2=14343.09    
  real, parameter   :: flimit=0.760545
  real              :: fphilf, uhilf
  real              :: cosm, sinm, depth, omega, omh, cm
  real, parameter          :: a0=1., a1=0.666, a2=0.445, a3=-0.105, a4=0.272
  integer :: i, j
    
    
  call mpp_clock_begin(id_diag)

  if(debug_this_module) then 
      write(stdoutunit,*) 'ocean_wave_model: start wave_diag'
  endif 


  if (use_TMA) then
    do j=jsd,jed 
      do i=isd,ied 
         cm=sqrt(Waves%xmom(i,j,tau_w)**2+Waves%ymom(i,j,tau_w)**2+epsln)
         cosm=Waves%xmom(i,j,tau_w)/cm
         sinm=Waves%ymom(i,j,tau_w)/cm
! The relevant 10m wind is the projection at the wave propagation direction
! This is not stated in the paper.
         uhilf=windx(i,j)*cosm+windy(i,j)*sinm
         depth=Grd%ht(i,j)
! calculate  omega=2*pi*f_p
         fphilf=const1*(max(epsln,uhilf)**2/(cm+epsln)**3)**seventh
         if (fphilf*uhilf.le.flimit) fphilf=(const2*cm)**(-third)      
         Waves%wave_p(i,j)=fphilf
         omega = twopi*Waves%wave_p(i,j)  
         omh   = omega**2*depth/grav
      
! using Pade approximation for wave number
        Waves%wave_k(i,j)=sqrt(omh**2+omh/((((a4*omh+a3)*omh+a2)*omh+a1)*omh+a0)) &
                   /(depth+epsln)
! calculate significant wave height
! height=4*sigma
        Waves%height(i,j)=4.0*sqrt(grav*cm/(omega+epsln))   ! cm ist scaled with 1/grav
! using TMA approximation, Hughes 1984
        if (omh.le.4.) then
          if (omh.le.1) then
            Waves%height(i,j)=Waves%height(i,j)*sqrt(0.5*omh)
          else
            Waves%height(i,j)=Waves%height(i,j)*sqrt(1.0-0.5*(2.0-sqrt(omh))**2)
          endif
        endif
 
      enddo
    enddo    

  else    ! use_TMA=.false. in following block 

    do j=jsd,jed 
      do i=isd,ied 
         cm=sqrt(Waves%xmom(i,j,tau_w)**2+Waves%ymom(i,j,tau_w)**2+epsln)
         cosm=Waves%xmom(i,j,tau_w)/cm
         sinm=Waves%ymom(i,j,tau_w)/cm
! The relevant 10m wind is the projection at the wave propagation direction
! This is not stated in the paper.
         uhilf=windx(i,j)*cosm+windy(i,j)*sinm
         depth=Grd%ht(i,j)
! calculate  omega=2*pi*f_p
         fphilf=const1*(max(epsln,uhilf)**2/(cm+epsln)**3)**seventh
         if (fphilf*uhilf.le.flimit) fphilf=(const2*cm)**(-third)      
         Waves%wave_p(i,j)=fphilf
         omega = twopi*Waves%wave_p(i,j)  
         omh   = omega**2*depth/grav
      
! using Pade approximation for wave number
        Waves%wave_k(i,j)=sqrt(omh**2+omh/((((a4*omh+a3)*omh+a2)*omh+a1)*omh+a0)) &
                   /(depth+epsln)
! calculate significant wave height
! height=4*sigma
        Waves%height(i,j)=4.0*sqrt(grav*cm/(omega+epsln))
 
      enddo
    enddo    
  
  endif

  if(debug_this_module) then 
      write(stdoutunit,*) 'ocean_wave_model: end wave_diag'
  endif 

  
  call mpp_clock_end(id_diag)

  return
end subroutine ocean_wave_diag
! </SUBROUTINE> NAME="ocean_wave_diag"


!#######################################################################
! <SUBROUTINE NAME="ocean_wave_filter">
!
! <DESCRIPTION>
! wave filter
! </DESCRIPTION>
!
subroutine ocean_wave_filter(Waves, dttw)
  type(ocean_wave_type), intent(inout) :: Waves
  real,                  intent(in)    :: dttw

  real, parameter :: wi=100.0
  integer :: i, j, n, ium, jum, im
  real    :: wich, rm


  if(debug_this_module) then 
     write(stdoutunit,*) 'ocean_wave_model: start wave_filter'
  endif 

  wich=min(0.9,wi/dttw)

  do j=jsc,jec
    do i=isc,iec
      n=0
      wrk1_2d(i,j)=0.0
      wrk2_2d(i,j)=0.0
      if (Grd%kmt(i,j).gt.0) then
         do ium=-1,1
            do jum=-1,1
               rm=Grd%tmask(i+ium,j+jum,1)
               im=nint(rm)
               n=n+im
               wrk1_2d(i,j)=wrk1_2d(i,j)+Waves%xmom(i+ium,j+jum,taup1_w)*rm
               wrk2_2d(i,j)=wrk2_2d(i,j)+Waves%ymom(i+ium,j+jum,taup1_w)*rm
            enddo
         enddo
         wrk1_2d(i,j)=wrk1_2d(i,j)/float(n)
         wrk2_2d(i,j)=wrk2_2d(i,j)/float(n)
      endif
   enddo
  enddo

  do j=jsc,jec
    do i=isc,iec
      if (Grd%kmt(i,j).gt.0) then
       Waves%xmom(i,j,taup1_w)=wich*Waves%xmom(i,j,taup1_w)+(1.-wich)*wrk1_2d(i,j)
       Waves%ymom(i,j,taup1_w)=wich*Waves%ymom(i,j,taup1_w)+(1.-wich)*wrk2_2d(i,j)
      endif
    enddo
  enddo

  call mpp_update_domains (Waves%xmom(:,:,taup1_w), Dom%domain2d)
  call mpp_update_domains (Waves%ymom(:,:,taup1_w), Dom%domain2d)

  if(debug_this_module) then 
      write(stdoutunit,*) 'ocean_wave_model: ending wave_filter'
  endif 
 
end subroutine ocean_wave_filter
! </SUBROUTINE> NAME="ocean_wave_filter"


!#######################################################################
! <SUBROUTINE NAME="read_wave">
!
! <DESCRIPTION>
!  Read wave restart information. 
! </DESCRIPTION>

subroutine read_wave(Waves)
  type(ocean_wave_type), intent(inout) :: Waves

  character*128 file_name

  file_name = 'ocean_wave.res.nc'
  id_restart(1) = register_restart_field(wave_restart, file_name, 'xmom', &
       Waves%xmom(:,:,taup1_w), domain=Dom%domain2d)
  id_restart(2) = register_restart_field(wave_restart, file_name, 'ymom', &
       Waves%ymom(:,:,taup1_w), domain=Dom%domain2d)       
  
  file_name = 'INPUT/ocean_wave.res.nc'
  if(.NOT. file_exist(trim(file_name)) ) return
  
  call restore_state(wave_restart)

  write (stdoutunit,'(/a)') &
  '  Reading TWO_LEVEL restart from INPUT/ocean_wave.res.nc'
  write (stdoutunit,'(a)')  &
  '  Expecting only one time record for each restart field.'
  
  call mpp_update_domains(Waves%xmom(:,:,taup1_w), Dom%domain2d)
  call mpp_update_domains(Waves%ymom(:,:,taup1_w), Dom%domain2d)
  Waves%xmom(:,:,tau_w) = Waves%xmom(:,:,taup1_w)
  Waves%ymom(:,:,tau_w) = Waves%ymom(:,:,taup1_w)
  
  return
end subroutine read_wave
! </SUBROUTINE> NAME="read_wave"


!#######################################################################
! <SUBROUTINE NAME="ocean_wave_restart">
!
! <DESCRIPTION>
!  Save wave restart information. 
! </DESCRIPTION>

subroutine ocean_wave_restart(Waves, time_stamp)
  type(ocean_wave_type), intent(inout)        :: Waves
  character(len=*),      intent(in), optional :: time_stamp
  
  call reset_field_pointer(wave_restart, id_restart(1),  Waves%xmom(:,:,taup1_w) )
  call reset_field_pointer(wave_restart, id_restart(2),  Waves%ymom(:,:,taup1_w) )

  call save_restart(wave_restart, time_stamp)
  return
end subroutine ocean_wave_restart
! </SUBROUTINE> NAME="ocean_wave_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_wave_end">
!
! <DESCRIPTION>
!  Write out external mode fields to restart file. 
! </DESCRIPTION>
subroutine ocean_wave_end(Time, Waves)
  type(ocean_time_type),  intent(in)    :: Time
  type(ocean_wave_type),  intent(inout) :: Waves

  if ( .not.use_this_module ) return

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_wave_mod (ocean_wave_end): module must be initialized')
  endif 

  if(.not. write_a_restart) then
    write(stdoutunit,'(/a)') &
    '==>Warning from ocean_wave_mod (ocean_wave_end): NO restart written.'
    call mpp_error(WARNING, &
    '==>Warning from ocean_wave_mod (ocean_wave_end): NO restart written.')
    return
  endif 

  call ocean_wave_restart(Waves)

  write(stdoutunit,*) ' '
  write(stdoutunit,*) &
  ' From ocean_wave_mod: ending external mode chksums at taup1'
  call write_timestamp(Time%model_time)
  call wave_chksum(Waves, taup1_w)

  module_is_initialized = .FALSE.

  nullify(Grd)
  nullify(Dom)

end subroutine ocean_wave_end
! </SUBROUTINE> NAME="ocean_wave_end"


!#######################################################################
! <SUBROUTINE NAME="wave_chksum">
!
! <DESCRIPTION>
!  Compute checksum for external mode fields.
! </DESCRIPTION>
!
subroutine wave_chksum(Waves, index)
  type(ocean_wave_type),  intent(inout) :: Waves
  integer,                intent(in)    :: index

  call write_chksum_2d('xmom', Waves%xmom(COMP,index))
  call write_chksum_2d('ymom', Waves%ymom(COMP,index))
  call write_chksum_2d('wave_k', Waves%wave_k(COMP))
  call write_chksum_2d('wave_p', Waves%wave_p(COMP))

end subroutine wave_chksum 
! </SUBROUTINE> NAME="wave_chksum"


!#######################################################################
! <FUNCTION NAME="wave_model_is_initialised">
!   <OVERVIEW>
!     Returns .true. if the wave model is initialised   
!   </OVERVIEW>
!   <DESCRIPTION>
!     This function returns .true. if the wave model is initialised
!     It is needed, because the wave model may be initialised after some module that requires a wave model
!   </DESCRIPTION>
!   <TEMPLATE>
!     use ocean_wave_mod,       only: wave_model_is_initialised     
!     if (wave_model_is_initialised() ) then
!   </TEMPLATE>
!   <IN NAME=""  TYPE="" >
!     No inputs needed.
!   </IN>
!   <OUT NAME=""  TYPE="logical" >
!     This function returns a logical. 
!   </OUT>

function wave_model_is_initialised()
  logical :: wave_model_is_initialised

  wave_model_is_initialised = .false.
  if(module_is_initialized) wave_model_is_initialised = .true.

end function wave_model_is_initialised
! </FUNCTION> NAME="wave_model_is_initialised"

end module ocean_wave_mod 

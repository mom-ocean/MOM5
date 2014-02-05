module ocean_coriolis_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> A. Rosati
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Compute the Coriolis acceleration for either Bgrid or Cgrid. 
!</OVERVIEW>

!<DESCRIPTION>
! This module computes Coriolis acceleration on either a Bgrid or Cgrid. 
! Coriolis and beta parameters are located at B-grid
! velocity point, which equals the C-grid vorticity point. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! R.C. Pacanowski and S.M. Griffies
! The MOM3 Manual (1999)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R.C. Pacanowski, and A. Rosati 
! A Technical Guide to MOM4 (2003)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies: Elements of MOM (2012)
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_coriolis_nml">
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging.  
!  </DATA> 
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to add contributions from Coriolis force.   
!  </DATA> 
!
!  <DATA NAME="acor" TYPE="real">
!  acor=0.0 means explicit Coriolis force.  0.5 < = acor < 1.0 means semi-implicit,
!  and acor = 1.0 is implicit.  This option is only relevant for the Bgrid, since
!  the C-grid compute Coriolis using 3rd order Adams-Bashforth scheme. For the Bgrid, the 
!  semi-implicit method removes dtuv time step constraint associated with inertial oscillations,
!  but it leads to Coriolis force affecting energy balances.  
!  If use two-level tendency discretization, then acor=0 is NOT allowed since the 
!  model will be linearly unstable with growth rate going as f*(delta time). 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: radius, radian, pi, epsln, omega
use diag_manager_mod, only: register_static_field, register_diag_field
use fms_mod,          only: write_version_number, mpp_error, FATAL, WARNING
use fms_mod,          only: check_nml_error, close_file, open_namelist_file
use mpp_mod,          only: input_nml_file, mpp_max, stdout, stdlog

use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: TWO_LEVEL
use ocean_parameters_mod, only: missing_value, onefourth
use ocean_parameters_mod, only: MOM_BGRID, MOM_CGRID 
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,      only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,      only: ocean_velocity_type, ocean_adv_vel_type
use ocean_types_mod,      only: ocean_options_type, ocean_thickness_type
use ocean_util_mod,       only: write_timestamp, diagnose_2d_u, diagnose_3d_u, diagnose_3d_en, write_chksum_3d
use ocean_workspace_mod,  only: wrk1, wrk2, wrk1_v  

implicit none

private

! for diagnostics 
logical :: used
integer :: id_coriolis  =-1
integer :: id_beta      =-1
integer :: id_beta_eff  =-1
integer :: id_cor_u     =-1
integer :: id_cor_v     =-1
integer :: id_hrho_cor_u=-1
integer :: id_hrho_cor_v=-1
integer :: id_ucori_impl=-1
integer :: id_vcori_impl=-1

#include <ocean_memory.h>

public ocean_coriolis_init
public coriolis_force_bgrid
public coriolis_force_bgrid_implicit
public coriolis_force_cgrid

type(ocean_grid_type), pointer :: Grd =>NULL()

character(len=128) :: version = &
     '$Id: ocean_coriolis.F90,v 20.0 2013/12/14 00:10:38 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

! for Bgrid or Cgrid
integer :: horz_grid

real    :: dtime 
logical :: module_is_initialized = .FALSE.

logical :: debug_this_module = .FALSE.
logical :: use_this_module   = .true.
real    :: acor              = 0.5

namelist /ocean_coriolis_nml/  debug_this_module, use_this_module, acor 

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_coriolis_init">
!
! <DESCRIPTION>
! Initialize the Coriolis module.
! </DESCRIPTION>
!
subroutine ocean_coriolis_init(Grid, Domain, Time, Time_steps, Ocean_options, hor_grid, debug)

  type(ocean_grid_type),       intent(inout),target :: Grid
  type(ocean_domain_type),     intent(in),   target :: Domain  
  type(ocean_time_type),       intent(in)           :: Time
  type(ocean_time_steps_type), intent(inout)        :: Time_steps
  type(ocean_options_type),    intent(inout)        :: Ocean_options
  integer,                     intent(in)           :: hor_grid
  logical,                     intent(in), optional :: debug

  real    :: deg2m, sin1, cos1, y, fmax
  real    :: max_dt_for_inertial_oscillation
  integer :: i, j
  integer :: ioun, io_status, ierr

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_coriolis_mod (ocean_coriolis_init): module already initialized.')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_coriolis_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_coriolis_nml')
#else
  ioun = open_namelist_file()
  read  (ioun, ocean_coriolis_nml,iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_coriolis_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_coriolis_nml)  
  write (stdlogunit, ocean_coriolis_nml)

  if (PRESENT(debug) .and. .not. debug_this_module) then 
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') &
    '==>Note: running ocean_coriolis_mod with debug_this_module=.true.'  
  endif 
  if(.not. use_this_module) then 
    write(stdoutunit,'(a)') &
    '==>Warning: NOT adding Coriolis to acceleration. Simulation has no rotational effects.'  
  endif 

  Grd=> Grid
  dtime     = Time_steps%dtime_u
  horz_grid = hor_grid 

  ! consistency checks for acor parameter 
  if(horz_grid == MOM_CGRID) then 

     acor = 0.0
     write(stdoutunit,'(a)') &
     '==>Note: Cgrid MOM computes Coriolis force using 3rd order Adams-Bashforth, just as for momentum advection.'
     write(stdoutunit,'(a)') &
     '==>Note: Cgrid MOM sets acor=0 since cannot compute C-grid Coriolis force implicitly in time.'
     Ocean_options%coriolis = 'Computed Coriolis force on Cgrid using 3rd order Adams-Bashforth time stepping.'

  else ! Bgrid 

     if(acor == 0.0) then 
        write(stdoutunit,'(a)')'==>Note: Coriolis force on Bgrid computed explicitly in time.'
        Ocean_options%coriolis = 'Computed Coriolis force on Bgrid using explicit time stepping.'
     elseif(acor == 1.0) then 
        write(stdoutunit,'(a)')'==>Note: Coriolis on Bgrid computed implicitly so to remove inertial time step constraint.'
        Ocean_options%coriolis = 'Computed Coriolis force on Bgrid using implicit time stepping.'
     elseif(acor >= 0.5 .and. acor < 1.0) then 
        write(stdoutunit,'(a)')'==>Note: Coriolis on Bgrid computed semi-implicitly to remove inertial time step constraint.'
        Ocean_options%coriolis = 'Computed Coriolis force on Bgrid using semi-implicit time stepping.'
     elseif(acor > 0.0 .and. acor < 0.5) then 
        acor = 0.5 
        write(stdoutunit,'(a)')'==>Warning: Coriolis "acor" parameter spuriously set outside relevant range. Resetting acor to MOM default 0.5.'
        Ocean_options%coriolis = 'Computed Coriolis force on Bgrid using semi-implicit time stepping to remove inertial time step constraint.'
     elseif(acor < 0.0 .or. acor > 1.0) then 
        acor = 0.5 
        write(stdoutunit,'(a)')'==>Warning: Coriolis "acor" parameter spuriously set outside relevant range. Resetting acor to MOM default 0.5.'
        Ocean_options%coriolis = 'Computed Coriolis force on Bgrid using semi-implicit time stepping.'
     endif 

     if(Time_steps%tendency==TWO_LEVEL .and. acor==0.0) then  
        call mpp_error(WARNING,&
        '==>ocean_coriolis_mod: acor=0.0 w/ 2-level tendency is unstable. MOM reseting acor to default acor=0.5')
        write(stdoutunit,'(a)')'==>Warning: acor=0.0 w/ 2-level tendency is unstable. MOM reseting acor to default acor=0.5.'
        acor = 0.5 
     endif 

  endif 
  Time_steps%acor = acor


#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
  allocate (Grid%f(isd:ied,jsd:jed))
  allocate (Grid%fstar(isd:ied,jsd:jed))
  allocate (Grid%beta(isd:ied,jsd:jed))
  allocate (Grid%beta_eff(isd:ied,jsd:jed))
#endif

  ! coriolis and beta parameters are located at B-grid
  ! velocity point, which equals the C-grid vorticity point. 
  if (Grid%beta_plane .or. Grid%f_plane) then

    ! beta plane with f     = f0     + beta*y where f0    is at f_plane_latitude
    ! beta plane with fstar = fstar0 + beta*y where fstar is at f_plane_latitude
    ! if f_plane then beta  = 0

    if (Grid%f_plane) then
      Grid%beta(:,:) = 0.0
    else
      Grid%beta(:,:) = 2.0*omega*cos(Grid%f_plane_latitude*pi/180.0)/radius
    endif

    deg2m = radius/radian
    sin1  = sin(Grid%f_plane_latitude*pi/180.0)
    cos1  = cos(Grid%f_plane_latitude*pi/180.0)
    do j=jsd,jed
      do i=isd,ied
        y = (Grid%yu(i,j)-Grid%f_plane_latitude)*deg2m
        Grid%f(i,j)     = 2.0*omega*sin1  + y*Grid%beta(i,j)
        Grid%fstar(i,j) = 2.0*omega*cos1  + y*Grid%beta(i,j)
      enddo
    enddo
    if (Grid%f_plane) then
      write (stdoutunit,'(//,a,f6.2,a//)') &
      ' Note: "f plane" set using f0 at latitude =', Grid%f_plane_latitude,'deg'
    else
      write (stdoutunit,'(//,a,f6.2,a,g14.7//)') &
      ' Note: "beta plane" set using f0 at latitude=',&
       Grid%f_plane_latitude,'deg and beta =',Grid%beta(isc,jsc)
    endif

  else

    Grid%f(:,:)    = 2.0*omega*sin(Grid%phiu(:,:))
    Grid%fstar(:,:)= 2.0*omega*cos(Grid%phiu(:,:))
    Grid%beta(:,:) = 2.0*omega*cos(Grid%phiu(:,:))/radius

  endif

  ! beta_eff centered on u-cell
  ! beta_eff = H*|grad(f/H)| --> |(f/hu)*dht_dx| + |beta-(f/hu)*dht_dy|
  Grid%beta_eff(:,:) = 0.0
  do j=jsc,jec
     do i=isc,iec
        Grid%beta_eff(i,j) = sqrt( (Grid%f(i,j)/(epsln+Grid%hu(i,j))*Grid%dht_dx(i,j))**2 &
                                  +(Grid%beta(i,j)-Grid%f(i,j)/(epsln+Grid%hu(i,j))*Grid%dht_dy(i,j))**2 )
     enddo
  enddo 

  ! check for marginally resolved inertial oscillation.
  ! assume 2*pi timesteps to resolve oscillation
  fmax = 2.0*omega*0.0002 ! ~ 0.01 deg latitude
  do j=jsc,jec
     do i=isc,iec
        if (Grid%kmu(i,j) /= 0) then
           fmax = max(fmax,abs(Grid%f(i,j)))
        endif
     enddo
  enddo
  call mpp_max (fmax)
  max_dt_for_inertial_oscillation = int(1.0/fmax) 
  write (stdoutunit,'(/1x, a,f8.0,a)')&
  ' ==> Note: 2*pi timesteps/(min inertial period) implies a maximum dtuv for time-explicit Coriolis =',&
  max_dt_for_inertial_oscillation,' sec.'

  if(horz_grid == MOM_BGRID) then 
     if (Time_steps%dtuv > max_dt_for_inertial_oscillation .and. acor==0.0) then
       call mpp_error(FATAL,&
       '==> Error in ocean_coriolis_mod: inertial oscillation not resolved. Reduce time step or set 0.5<=acor<=1.0.')
     endif  
  endif

  ! diagnostic manager registers and sends for static fields 

  id_coriolis  = register_static_field ('ocean_model', 'f_coriolis', Grd%vel_axes_uv(1:2), &
                                        'Coriolis frequency on U-cell', '1/s',             &
                                         missing_value=missing_value, range=(/-10.0,10.0/))
  call diagnose_2d_u(Time, Grd, id_coriolis, Grid%f(:,:))

  id_beta  = register_static_field ('ocean_model', 'beta', Grd%vel_axes_uv(1:2), &
                                    'planetary beta', '1/(m*s)',&
                                    missing_value=missing_value, range=(/-10.0,10.0/))
  call diagnose_2d_u(Time, Grd, id_beta, Grd%beta(:,:))

  id_beta_eff  = register_static_field ('ocean_model', 'beta_eff', Grd%vel_axes_uv(1:2), &
                                        'effective beta', '1/(m*s)',&
                                        missing_value=missing_value, range=(/-10.0,10.0/))
  call diagnose_2d_u(Time, Grd, id_beta_eff, Grd%beta_eff(:,:))

  ! diagnostic manager registers for dynamic fields 
  id_cor_u =  register_diag_field ('ocean_model', 'cor_u', Grd%vel_axes_u(1:3), Time%model_time, &
     'explicit coriolis accel in i-direct', 'm/s^2', missing_value=missing_value, range=(/-1e9,1e9/))
  id_cor_v =  register_diag_field ('ocean_model', 'cor_v', Grd%vel_axes_v(1:3), Time%model_time, &
     'explicit coriolis accel in j-direct', 'm/s^2', missing_value=missing_value, range=(/-1e9,1e9/))

  id_hrho_cor_u =  register_diag_field ('ocean_model', 'hrho_cor_u', Grd%vel_axes_u(1:3), Time%model_time,&
                  'rho*dz weighted explicit coriolis accel in i-direct', 'N/m^2',                         &
                   missing_value=missing_value, range=(/-1e9,1e9/))
  id_hrho_cor_v =  register_diag_field ('ocean_model', 'hrho_cor_v', Grd%vel_axes_v(1:3), Time%model_time,&
                  'rho*dz weighted explicit coriolis accel in j-direct', 'N/m^2',                         &
                   missing_value=missing_value, range=(/-1e9,1e9/))

  id_ucori_impl = register_diag_field ('ocean_model', 'ucori_impl', Grd%vel_axes_u(1:3), Time%model_time, &
     'implicit Coriolis force in i-direction', 'N/m^2', missing_value=missing_value, range=(/-1e9,1e9/))
  id_vcori_impl = register_diag_field ('ocean_model', 'vcori_impl', Grd%vel_axes_v(1:3), Time%model_time, &
     'implicit Coriolis force in j-direction', 'N/m^2', missing_value=missing_value, range=(/-1e9,1e9/))


end subroutine ocean_coriolis_init
! </SUBROUTINE> NAME="ocean_coriolis_init"


!#######################################################################
! <SUBROUTINE NAME="coriolis_force_bgrid">
!
! <DESCRIPTION>
! Compute thickness and density weighted acceleration due to Coriolis
! force on a B-grid. 
! </DESCRIPTION>
!
subroutine coriolis_force_bgrid(Time, Thickness, Velocity, energy_analysis_step)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_velocity_type),  intent(inout) :: Velocity
  logical,                    intent(in)    :: energy_analysis_step

  integer :: itime, tau, taum1, taup1
  integer :: i, j, k, n
  
  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_coriolis_mod (coriolis_force): module must be initialized')
  endif 

  if(.not. use_this_module) then 
      if(energy_analysis_step) then 
          Velocity%wrkv=0.0
      endif
      return
  endif

  wrk1_v = 0.0

  tau    = Time%tau
  taum1  = Time%taum1
  taup1  = Time%taup1

  if(acor == 0.0) then 
    itime = tau
  else 
    itime = taum1
  endif 


  if(energy_analysis_step) then 

      if(acor > 0) then 

          ! Coriolis force alters kinetic energy when acor>0.  To determine 
          ! effects via the energy analysis, compute Coriolis force as here.  
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   Velocity%wrkv(i,j,k,1) = Grd%f(i,j)  &
                        *(Velocity%u(i,j,k,2,taum1)*(1.0-acor) + acor*Velocity%u(i,j,k,2,taup1)) &
                        *Thickness%rho_dzu(i,j,k,tau)
                   Velocity%wrkv(i,j,k,2) =-Grd%f(i,j) &
                        *(Velocity%u(i,j,k,1,taum1)*(1.0-acor) + acor*Velocity%u(i,j,k,1,taup1)) &
                        *Thickness%rho_dzu(i,j,k,tau)
                enddo
             enddo
          enddo

      else 

          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   Velocity%wrkv(i,j,k,1) =  Grd%f(i,j)*Velocity%u(i,j,k,2,itime)*Thickness%rho_dzu(i,j,k,tau)
                   Velocity%wrkv(i,j,k,2) = -Grd%f(i,j)*Velocity%u(i,j,k,1,itime)*Thickness%rho_dzu(i,j,k,tau)
                enddo
             enddo
          enddo
      endif

  else 

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,1) =  Grd%f(i,j)*Velocity%u(i,j,k,2,itime)
               wrk1_v(i,j,k,2) = -Grd%f(i,j)*Velocity%u(i,j,k,1,itime)
            enddo
         enddo
      enddo

      call diagnose_3d_u(Time, Grd, id_cor_u, wrk1_v(:,:,:,1))
      call diagnose_3d_u(Time, Grd, id_cor_v, wrk1_v(:,:,:,2))

      ! weight acceleration with thickness and density of velocity cell
      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1_v(i,j,k,n)         = wrk1_v(i,j,k,n)*Thickness%rho_dzu(i,j,k,tau)
                  Velocity%accel(i,j,k,n) = Velocity%accel(i,j,k,n) + wrk1_v(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

      call diagnose_3d_u(Time, Grd, id_hrho_cor_u, wrk1_v(:,:,:,1))
      call diagnose_3d_u(TIme, Grd, id_hrho_cor_v, wrk1_v(:,:,:,2))
  endif


end subroutine coriolis_force_bgrid
! </SUBROUTINE> NAME="coriolis_force_bgrid"


!#######################################################################
! <SUBROUTINE NAME="coriolis_force_bgrid_implicit">
!
! <DESCRIPTION>
! Contributions to thickness weighted and density weighted 
! acceleration from time-implicit Coriolis force.
! </DESCRIPTION>
!
subroutine coriolis_force_bgrid_implicit(Time, Velocity)

  type(ocean_time_type),     intent(in)    :: Time
  type(ocean_velocity_type), intent(inout) :: Velocity

  integer :: i, j, k
  real    :: dtimeacor
  real    :: lambda,factor

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) then 
      return
  endif


  wrk1_v = 0.0
  
  if(acor > 0.0) then 

      dtimeacor = dtime*acor 

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,:) = Velocity%accel(i,j,k,:)  
               lambda          = dtimeacor*Grd%f(i,j)
               factor          = 1.0/(1.0 + lambda*lambda) 
               wrk1(i,j,k)     = (Velocity%accel(i,j,k,1) + lambda*Velocity%accel(i,j,k,2))*factor
               wrk2(i,j,k)     = (Velocity%accel(i,j,k,2) - lambda*Velocity%accel(i,j,k,1))*factor
            enddo
         enddo
      enddo

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec       
               Velocity%accel(i,j,k,1) = wrk1(i,j,k) 
               Velocity%accel(i,j,k,2) = wrk2(i,j,k) 
               wrk1_v(i,j,k,:)         = Velocity%accel(i,j,k,:) - wrk1_v(i,j,k,:) 
            enddo
         enddo
      enddo

      if (debug_this_module) then
          write(stdoutunit,*) ' ' 
          write(stdoutunit,*) 'From ocean_coriolis_mod: chksums after acor>0 Coriolis' 
          call write_timestamp(Time%model_time)
          call write_chksum_3d('accel(1)', Velocity%accel(COMP,:,1))
          call write_chksum_3d('accel(2)', Velocity%accel(COMP,:,2))
          call write_chksum_3d('vel(1)', wrk1_v(COMP,:,1))
          call write_chksum_3d('vel(2)', wrk1_v(COMP,:,2))
      endif

  endif

  ! send to diagnostics manager  
  call diagnose_3d_u(Time, Grd, id_ucori_impl, wrk1_v(:,:,:,1))
  call diagnose_3d_u(Time, Grd, id_vcori_impl, wrk1_v(:,:,:,2))

end subroutine coriolis_force_bgrid_implicit
! </SUBROUTINE> NAME="coriolis_force_bgrid_implicit"


!#######################################################################
! <SUBROUTINE NAME="coriolis_force_cgrid">
!
! <DESCRIPTION>
!
! Compute thickness and density weighted acceleration due to Coriolis
! force on a C-grid. 
! 
! </DESCRIPTION>
!
subroutine coriolis_force_cgrid(Time, Adv_vel, Velocity, abtau_m0, abtau_m1, abtau_m2, energy_analysis_step)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_adv_vel_type),   intent(in)    :: Adv_vel
  type(ocean_velocity_type),  intent(inout) :: Velocity
  real,                       intent(in)    :: abtau_m0
  real,                       intent(in)    :: abtau_m1
  real,                       intent(in)    :: abtau_m2
  logical,                    intent(in)    :: energy_analysis_step

  integer :: tau
  integer :: tau_m0, tau_m1, tau_m2
  integer :: i, j, k, n
  
  if(.not. use_this_module) then 
      if(energy_analysis_step) then 
          Velocity%wrkv=0.0
      endif
      return
  endif

  wrk1_v = 0.0
  tau    = Time%tau
  tau_m0 = Time%tau_m0
  tau_m1 = Time%tau_m1
  tau_m2 = Time%tau_m2

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1_v(i,j,k,1) =  onefourth*Grd%tmasken(i,j,k,1)                                      &
                            *( Grd%f(i,j)  *(Adv_vel%vhrho_nt(i,j,k)  +Adv_vel%vhrho_nt(i+1,j,k)) &
                              +Grd%f(i,j-1)*(Adv_vel%vhrho_nt(i,j-1,k)+Adv_vel%vhrho_nt(i+1,j-1,k)) )
           wrk1_v(i,j,k,2) = -onefourth*Grd%tmasken(i,j,k,2)                                      &
                            *( Grd%f(i,j)  *(Adv_vel%uhrho_et(i,j,k)  +Adv_vel%uhrho_et(i,j+1,k)) &
                              +Grd%f(i-1,j)*(Adv_vel%uhrho_et(i-1,j,k)+Adv_vel%uhrho_et(i-1,j+1,k)) )
        enddo
     enddo
  enddo

  if(energy_analysis_step) then 
     do n=1,2 
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 Velocity%wrkv(i,j,k,n) = wrk1_v(i,j,k,n)
              enddo
           enddo
        enddo
     enddo
  else 
     do n=1,2
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 Velocity%coriolis(i,j,k,n,tau_m0) = Velocity%coriolis(i,j,k,n,tau_m0) + wrk1_v(i,j,k,n)
                 Velocity%accel(i,j,k,n) = Velocity%accel(i,j,k,n)                     + &
                                           abtau_m0*Velocity%coriolis(i,j,k,n,tau_m0)  + &
                                           abtau_m1*Velocity%coriolis(i,j,k,n,tau_m1)  + &
                                           abtau_m2*Velocity%coriolis(i,j,k,n,tau_m2)
              enddo
           enddo
        enddo
     enddo


     call diagnose_3d_en(Time, Grd, id_hrho_cor_u, id_hrho_cor_v, wrk1_v(:,:,:,:))
  endif

end subroutine coriolis_force_cgrid
! </SUBROUTINE> NAME="coriolis_force_cgrid"


end module ocean_coriolis_mod

module atmos_model_mod
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Bruce Wyman  
!</CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  Zhi Liang
! </REVIEWER>
!-----------------------------------------------------------------------
!<OVERVIEW>
!  Driver for the atmospheric model, contains routines to advance the
!  atmospheric model state by one time step.
!</OVERVIEW>

!<DESCRIPTION>
!     This version of atmos_model_mod has been designed around the implicit
!     version diffusion scheme of the GCM. It requires two routines to advance
!     the atmospheric model one time step into the future. These two routines
!     correspond to the down and up sweeps of the standard tridiagonal solver.
!     Most atmospheric processes (dynamics,radiation,etc.) are performed
!     in the down routine. The up routine finishes the vertical diffusion
!     and computes moisture related terms (convection,large-scale condensation,
!     and precipitation).

!     The boundary variables needed by other component models for coupling
!     are contained in a derived data type. A variable of this derived type
!     is returned when initializing the atmospheric model. It is used by other
!     routines in this module and by coupling routines. The contents of
!     this derived type should only be modified by the atmospheric model.

!</DESCRIPTION>

use mpp_mod,            only: mpp_pe, mpp_root_pe, mpp_clock_id, mpp_clock_begin
use mpp_mod,            only: mpp_clock_end, CLOCK_COMPONENT, mpp_error, mpp_chksum
use mpp_domains_mod,    only: domain2d
#ifdef INTERNAL_FILE_NML
use mpp_mod,            only: input_nml_file
#else
use fms_mod,            only: open_namelist_file
#endif
use fms_mod,            only: file_exist, error_mesg, field_size, FATAL, NOTE, WARNING
use fms_mod,            only: close_file,  write_version_number, stdlog, stdout
use fms_mod,            only: read_data, write_data, clock_flag_default
use fms_mod,            only: open_restart_file, check_nml_error
use fms_io_mod,         only: get_restart_io_mode
use fms_io_mod,         only: restart_file_type, register_restart_field
use fms_io_mod,         only: save_restart, restore_state, get_mosaic_tile_file
use time_manager_mod,   only: time_type, operator(+), get_time
use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_number_tracers, get_tracer_index, NO_TRACER
use diag_integral_mod,  only: diag_integral_init, diag_integral_end
use diag_integral_mod,  only: diag_integral_output
use atmosphere_mod,     only: atmosphere_cell_area
use xgrid_mod,          only: grid_box_type
use atmosphere_mod,     only: atmosphere_up, atmosphere_down, atmosphere_init
use atmosphere_mod,     only: atmosphere_end, get_bottom_mass, get_bottom_wind
use atmosphere_mod,     only: atmosphere_resolution, atmosphere_domain
use atmosphere_mod,     only: atmosphere_boundary, get_atmosphere_axes
use atmosphere_mod,     only: get_stock_pe
use atmosphere_mod,     only: surf_diff_type
use atmosphere_mod,     only: atmosphere_restart
use coupler_types_mod,  only: coupler_2d_bc_type


!-----------------------------------------------------------------------

implicit none
private

public update_atmos_model_down, update_atmos_model_up
public atmos_model_init, atmos_model_end, atmos_data_type
public land_ice_atmos_boundary_type, land_atmos_boundary_type
public atm_stock_pe
public ice_atmos_boundary_type
public atmos_model_restart
public atmos_data_type_chksum
public lnd_ice_atm_bnd_type_chksum, lnd_atm_bnd_type_chksum
public ice_atm_bnd_type_chksum
!-----------------------------------------------------------------------

!<PUBLICTYPE >
 type atmos_data_type
     type (domain2d)               :: domain             ! domain decomposition
     integer                       :: axes(4)            ! axis indices (returned by diag_manager) for the atmospheric grid 
                                                         ! (they correspond to the x, y, pfull, phalf axes)
     real, pointer, dimension(:,:) :: lon_bnd  => null() ! local longitude axis grid box corners in radians.
     real, pointer, dimension(:,:) :: lat_bnd  => null() ! local latitude axis grid box corners in radians.
     real, pointer, dimension(:,:) :: t_bot    => null() ! temperature at lowest model level
     real, pointer, dimension(:,:,:) :: tr_bot   => null() ! tracers at lowest model level
     real, pointer, dimension(:,:) :: z_bot    => null() ! height above the surface for the lowest model level
     real, pointer, dimension(:,:) :: p_bot    => null() ! pressure at lowest model level
     real, pointer, dimension(:,:) :: u_bot    => null() ! zonal wind component at lowest model level
     real, pointer, dimension(:,:) :: v_bot    => null() ! meridional wind component at lowest model level
     real, pointer, dimension(:,:) :: p_surf   => null() ! surface pressure 
     real, pointer, dimension(:,:) :: slp      => null() ! sea level pressure 
     real, pointer, dimension(:,:) :: gust     => null() ! gustiness factor
     real, pointer, dimension(:,:) :: coszen   => null() ! cosine of the zenith angle
     real, pointer, dimension(:,:) :: flux_sw  => null() ! net shortwave flux (W/m2) at the surface
     real, pointer, dimension(:,:) :: flux_sw_dir            =>null()
     real, pointer, dimension(:,:) :: flux_sw_dif            =>null()
     real, pointer, dimension(:,:) :: flux_sw_down_vis_dir   =>null()
     real, pointer, dimension(:,:) :: flux_sw_down_vis_dif   =>null()
     real, pointer, dimension(:,:) :: flux_sw_down_total_dir =>null()
     real, pointer, dimension(:,:) :: flux_sw_down_total_dif =>null()
     real, pointer, dimension(:,:) :: flux_sw_vis            =>null()
     real, pointer, dimension(:,:) :: flux_sw_vis_dir        =>null()
     real, pointer, dimension(:,:) :: flux_sw_vis_dif        =>null()
     real, pointer, dimension(:,:) :: flux_lw  => null() ! net longwave flux (W/m2) at the surface
     real, pointer, dimension(:,:) :: lprec    => null() ! mass of liquid precipitation since last time step (Kg/m2)
     real, pointer, dimension(:,:) :: fprec    => null() ! ass of frozen precipitation since last time step (Kg/m2)
     logical, pointer, dimension(:,:) :: maskmap =>null()! A pointer to an array indicating which
                                                         ! logical processors are actually used for
                                                         ! the ocean code. The other logical
                                                         ! processors would be all land points and
                                                         ! are not assigned to actual processors.
                                                         ! This need not be assigned if all logical
                                                         ! processors are used. This variable is dummy and need 
                                                         ! not to be set, but it is needed to pass compilation.
     type (surf_diff_type)         :: Surf_diff          ! store data needed by the multi-step version of the diffusion algorithm
     type (time_type)              :: Time               ! current time
     type (time_type)              :: Time_step          ! atmospheric time step.
     type (time_type)              :: Time_init          ! reference time.
     integer, pointer              :: pelist(:) =>null() ! pelist where atmosphere is running.
     logical                       :: pe                 ! current pe.
     type(coupler_2d_bc_type)      :: fields             ! array of fields used for additional tracers
     type(grid_box_type)           :: grid               ! hold grid information needed for 2nd order conservative flux exchange 
                                                         ! to calculate gradient on cubic sphere grid.
 end type atmos_data_type
!</PUBLICTYPE >

!<PUBLICTYPE >
type land_ice_atmos_boundary_type
   ! variables of this type are declared by coupler_main, allocated by flux_exchange_init.
!quantities going from land+ice to atmos
   real, dimension(:,:),   pointer :: t              =>null() ! surface temperature for radiation calculations
   real, dimension(:,:),   pointer :: albedo         =>null() ! surface albedo for radiation calculations
   real, dimension(:,:),   pointer :: albedo_vis_dir =>null()
   real, dimension(:,:),   pointer :: albedo_nir_dir =>null()
   real, dimension(:,:),   pointer :: albedo_vis_dif =>null()
   real, dimension(:,:),   pointer :: albedo_nir_dif =>null()
   real, dimension(:,:),   pointer :: land_frac      =>null() ! fraction amount of land in a grid box 
   real, dimension(:,:),   pointer :: dt_t           =>null() ! temperature tendency at the lowest level
   real, dimension(:,:,:), pointer :: dt_tr          =>null() ! tracer tendency at the lowest level
   real, dimension(:,:),   pointer :: u_flux         =>null() ! zonal wind stress
   real, dimension(:,:),   pointer :: v_flux         =>null() ! meridional wind stress
   real, dimension(:,:),   pointer :: dtaudu         =>null() ! derivative of zonal wind stress w.r.t. the lowest zonal level wind speed
   real, dimension(:,:),   pointer :: dtaudv         =>null() ! derivative of meridional wind stress w.r.t. the lowest meridional level wind speed
   real, dimension(:,:),   pointer :: u_star         =>null() ! friction velocity
   real, dimension(:,:),   pointer :: b_star         =>null() ! bouyancy scale
   real, dimension(:,:),   pointer :: q_star         =>null() ! moisture scale
   real, dimension(:,:),   pointer :: rough_mom      =>null() ! surface roughness (used for momentum)
   real, dimension(:,:),   pointer :: frac_open_sea  =>null() ! non-seaice fraction (%)
   real, dimension(:,:,:), pointer :: data           =>null() !collective field for "named" fields above
   integer                         :: xtype                   !REGRID, REDIST or DIRECT
end type land_ice_atmos_boundary_type
!</PUBLICTYPE >

!<PUBLICTYPE >
type :: land_atmos_boundary_type
   real, dimension(:,:), pointer :: data =>null() ! quantities going from land alone to atmos (none at present)
end type land_atmos_boundary_type
!</PUBLICTYPE >

!<PUBLICTYPE >
!quantities going from ice alone to atmos (none at present)
type :: ice_atmos_boundary_type
   real, dimension(:,:), pointer :: data =>null() ! quantities going from ice alone to atmos (none at present)
end type ice_atmos_boundary_type
!</PUBLICTYPE >

!Balaji
integer :: atmClock

!for restart
integer                 :: ipts, jpts, dto
type(restart_file_type), pointer, save :: Atm_restart => null()
type(restart_file_type), pointer, save :: Til_restart => null()
logical                                :: in_different_file = .false.

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: atmos_model.F90,v 20.0 2013/12/13 23:08:05 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

integer :: ivapor = NO_TRACER ! index of water vapor tracer

!-----------------------------------------------------------------------
character(len=80) :: restart_format = 'atmos_coupled_mod restart format 01'
!-----------------------------------------------------------------------
logical           :: do_netcdf_restart = .true.
logical           :: restart_tbot_qbot = .false.
namelist /atmos_model_nml/ do_netcdf_restart, restart_tbot_qbot  

contains

!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_down">
!
! <OVERVIEW>
!   compute the atmospheric tendencies for dynamics, radiation, 
!   vertical diffusion of momentum, tracers, and heat/moisture.
! </OVERVIEW>
!
!<DESCRIPTION>
!   Called every time step as the atmospheric driver to compute the
!   atmospheric tendencies for dynamics, radiation, vertical diffusion of
!   momentum, tracers, and heat/moisture.  For heat/moisture only the
!   downward sweep of the tridiagonal elimination is performed, hence
!   the name "_down". 
!</DESCRIPTION>

!   <TEMPLATE>
!     call  update_atmos_model_down( Surface_boundary, Atmos )
!   </TEMPLATE>

! <IN NAME = "Surface_boundary" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains quantities going from land+ice to atmos.  
! </IN>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </INOUT>

subroutine update_atmos_model_down( Surface_boundary, Atmos )
!
!-----------------------------------------------------------------------
  type(land_ice_atmos_boundary_type), intent(inout) :: Surface_boundary
  type (atmos_data_type), intent(inout) :: Atmos
                                      
!-----------------------------------------------------------------------
  call mpp_clock_begin(atmClock)

    call atmosphere_down (Atmos%Time, Surface_boundary%land_frac,        &
                          Surface_boundary%t,  Surface_boundary%albedo,  &
                          Surface_boundary%albedo_vis_dir,   &
                          Surface_boundary%albedo_nir_dir,   &
                          Surface_boundary%albedo_vis_dif,   &
                          Surface_boundary%albedo_nir_dif,   &
                          Surface_boundary%rough_mom,   &
                          Surface_boundary%u_star,      &
                          Surface_boundary%b_star,      &
                          Surface_boundary%q_star, &
                          Surface_boundary%dtaudu,      &
                          Surface_boundary%dtaudv,      &
                          Surface_boundary%u_flux,      &
                          Surface_boundary%v_flux,      &
                          Surface_boundary%frac_open_sea, &
                          Atmos%gust,                   &
                          Atmos%coszen,                 &
                          Atmos%flux_sw,                &
                          Atmos%flux_sw_dir,            &
                          Atmos%flux_sw_dif,            &
                          Atmos%flux_sw_down_vis_dir,   &
                          Atmos%flux_sw_down_vis_dif,   &
                          Atmos%flux_sw_down_total_dir, &
                          Atmos%flux_sw_down_total_dif, &
                          Atmos%flux_sw_vis,            &
                          Atmos%flux_sw_vis_dir,        &
                          Atmos%flux_sw_vis_dif,        &
                          Atmos%flux_lw,                &
                          Atmos%Surf_diff               )

!-----------------------------------------------------------------------

  call mpp_clock_end(atmClock)
 end subroutine update_atmos_model_down
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_up">
!
!-----------------------------------------------------------------------
! <OVERVIEW>
!   upward vertical diffusion of heat/moisture and moisture processes
! </OVERVIEW>

!<DESCRIPTION>
!   Called every time step as the atmospheric driver to finish the upward
!   sweep of the tridiagonal elimination for heat/moisture and compute the
!   convective and large-scale tendencies.  The atmospheric variables are
!   advanced one time step and tendencies set back to zero. 
!</DESCRIPTION>

! <TEMPLATE>
!     call  update_atmos_model_up( Surface_boundary, Atmos )
! </TEMPLATE>

! <IN NAME = "Surface_boundary" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains quantities going from land+ice to atmos.  
! </IN>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </INOUT>

 subroutine update_atmos_model_up( Surface_boundary, Atmos )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

type(land_ice_atmos_boundary_type), intent(in) :: Surface_boundary
type (atmos_data_type), intent(inout) :: Atmos
                                      
!-----------------------------------------------------------------------
  call mpp_clock_begin(atmClock)


    Atmos%Surf_diff%delta_t = Surface_boundary%dt_t
    Atmos%Surf_diff%delta_tr = Surface_boundary%dt_tr

    call atmosphere_up (Atmos%Time,  Surface_boundary%land_frac, Atmos%Surf_diff, &
                        Atmos%lprec, Atmos%fprec, Atmos%gust, &
                        Surface_boundary%u_star, Surface_boundary%b_star, Surface_boundary%q_star)

!   --- advance time ---

    Atmos % Time = Atmos % Time + Atmos % Time_step


    call get_bottom_mass (Atmos % t_bot,  Atmos % tr_bot, &
                          Atmos % p_bot,  Atmos % z_bot,  &
                          Atmos % p_surf, Atmos % slp     )

    call get_bottom_wind (Atmos % u_bot, Atmos % v_bot)


!------ global integrals ------

    call diag_integral_output (Atmos % Time)

!-----------------------------------------------------------------------
  call mpp_clock_end(atmClock)

end subroutine update_atmos_model_up
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_model_init">
!
! <OVERVIEW>
! Routine to initialize the atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!     This routine allocates storage and returns a variable of type
!     atmos_boundary_data_type, and also reads a namelist input and restart file. 
! </DESCRIPTION>

! <TEMPLATE>
!     call atmos_model_init (Atmos, Time_init, Time, Time_step)
! </TEMPLATE>

! <IN NAME="Time_init" TYPE="type(time_type)" >
!   The base (or initial) time of the experiment.
! </IN>

! <IN NAME="Time" TYPE="type(time_type)" >
!   The current time.
! </IN>

! <IN NAME="Time_step" TYPE="type(time_type)" >
!   The atmospheric model/physics time step.
! </IN>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>

subroutine atmos_model_init (Atmos, Time_init, Time, Time_step)


type (atmos_data_type), intent(inout) :: Atmos
type (time_type), intent(in) :: Time_init, Time, Time_step

  integer :: unit, ntrace, ntprog, ntdiag, ntfamily, i, j
  integer :: mlon, mlat, nlon, nlat, sec, day, dt
  character(len=80) :: control
  real, dimension(:,:), allocatable :: area
  integer :: ierr, io, logunit
  character(len=64) :: filename, filename2
  integer           :: id_restart
!-----------------------------------------------------------------------

!---- set the atmospheric model time ------

   Atmos % Time_init = Time_init
   Atmos % Time      = Time
   Atmos % Time_step = Time_step
   logunit = stdlog()

   IF ( file_exist('input.nml')) THEN
#ifdef INTERNAL_FILE_NML
      read(input_nml_file, nml=atmos_model_nml, iostat=io)
      ierr = check_nml_error(io, 'atmos_model_nml')
#else
      unit = open_namelist_file ( )
      ierr=1
      do while (ierr /= 0)
         read  (unit, nml=atmos_model_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'atmos_model_nml')
      enddo
 10     call close_file (unit)
#endif
   endif
   call get_restart_io_mode(do_netcdf_restart)

!-----------------------------------------------------------------------
! how many tracers have been registered?
!  (will print number below)
   call get_number_tracers ( MODEL_ATMOS, ntrace, ntprog, ntdiag, ntfamily )
   if ( ntfamily > 0 ) call error_mesg ('atmos_model', 'ntfamily > 0', FATAL)
   ivapor = get_tracer_index( MODEL_ATMOS, 'sphum' )
   if (ivapor==NO_TRACER) &
        ivapor = get_tracer_index( MODEL_ATMOS, 'mix_rat' )
   if (ivapor==NO_TRACER) &
        call error_mesg('atmos_model_init', 'Cannot find water vapor in ATM tracer table', FATAL)

!-----------------------------------------------------------------------
!  ----- initialize atmospheric model -----

    call atmosphere_init (Atmos%Time_init, Atmos%Time, Atmos%Time_step,&
                          Atmos%Surf_diff, Atmos%grid )
                           
!-----------------------------------------------------------------------
!---- allocate space ----

    call atmosphere_resolution (mlon, mlat, global=.true.)
    call atmosphere_resolution (nlon, nlat, global=.false.)
    call atmosphere_domain     (Atmos%domain)

    allocate ( Atmos %  lon_bnd (nlon+1,nlat+1), &
               Atmos %  lat_bnd (nlon+1,nlat+1), &
               Atmos % t_bot    (nlon,nlat), &
               Atmos % tr_bot    (nlon,nlat, ntprog), &
               Atmos % z_bot    (nlon,nlat), &
               Atmos % p_bot    (nlon,nlat), &
               Atmos % u_bot    (nlon,nlat), &
               Atmos % v_bot    (nlon,nlat), &
               Atmos % p_surf   (nlon,nlat), &
               Atmos % slp      (nlon,nlat), &
               Atmos % gust     (nlon,nlat), &
               Atmos % flux_sw  (nlon,nlat), &
               Atmos % flux_sw_dir (nlon,nlat), &
               Atmos % flux_sw_dif (nlon,nlat), &
               Atmos % flux_sw_down_vis_dir (nlon,nlat), &
               Atmos % flux_sw_down_vis_dif (nlon,nlat), &
               Atmos % flux_sw_down_total_dir (nlon,nlat), &
               Atmos % flux_sw_down_total_dif (nlon,nlat), &
               Atmos % flux_sw_vis (nlon,nlat), &
               Atmos % flux_sw_vis_dir (nlon,nlat), &
               Atmos % flux_sw_vis_dif(nlon,nlat), &
               Atmos % flux_lw  (nlon,nlat), &
               Atmos % coszen   (nlon,nlat), &
               Atmos % lprec    (nlon,nlat), &
               Atmos % fprec    (nlon,nlat)  )

    do j = 1, nlat
       do i = 1, nlon    
          Atmos % flux_sw(i,j)                 = 0.0
          Atmos % flux_lw(i,j)                 = 0.0    
          Atmos % flux_sw_dir (i,j)            = 0.0
          Atmos % flux_sw_dif (i,j)            = 0.0 
          Atmos % flux_sw_down_vis_dir (i,j)   = 0.0 
          Atmos % flux_sw_down_vis_dif (i,j)   = 0.0 
          Atmos % flux_sw_down_total_dir (i,j) = 0.0 
          Atmos % flux_sw_down_total_dif (i,j) = 0.0 
          Atmos % flux_sw_vis (i,j)            = 0.0 
          Atmos % flux_sw_vis_dir (i,j)        = 0.0 
          Atmos % flux_sw_vis_dif(i,j)         = 0.0 
          Atmos % coszen(i,j)                  = 0.0 
       enddo
    enddo
!-----------------------------------------------------------------------
!------ get initial state for dynamics -------

    call get_atmosphere_axes ( Atmos % axes )

    call atmosphere_boundary ( Atmos %  lon_bnd, Atmos %  lat_bnd, &
                               global=.false. )

    call get_bottom_mass (Atmos % t_bot,  Atmos % tr_bot, &
                          Atmos % p_bot,  Atmos % z_bot,  &
                          Atmos % p_surf, Atmos % slp     )

    call get_bottom_wind (Atmos % u_bot, Atmos % v_bot)

!-----------------------------------------------------------------------
!---- print version number to logfile ----

   call write_version_number(version, tagname)
!  write the namelist to a log file
   if (mpp_pe() == mpp_root_pe()) then
      unit = stdlog( )
      write (unit, nml=atmos_model_nml)
      call close_file (unit)
!  number of tracers
      write (unit, '(a,i3)') 'Number of tracers =', ntrace
      write (unit, '(a,i3)') 'Number of prognostic tracers =', ntprog
      write (unit, '(a,i3)') 'Number of diagnostic tracers =', ntdiag
   endif

!------ read initial state for several atmospheric fields ------
   filename = 'atmos_coupled.res.nc'
   call get_mosaic_tile_file(filename, filename2, .false., Atmos%domain ) 
   allocate(Atm_restart)
   if(trim(filename2) == trim(filename)) then
      Til_restart => Atm_restart
      in_different_file = .false.
      id_restart = register_restart_field(Atm_restart, filename, 'glon_bnd', ipts, domain=Atmos%domain)
      id_restart = register_restart_field(Atm_restart, filename, 'glat_bnd', jpts, domain=Atmos%domain)
      id_restart = register_restart_field(Atm_restart, filename, 'dt', dto, domain=Atmos%domain)
   else
      in_different_file = .true.
      allocate(Til_restart)
      id_restart = register_restart_field(Atm_restart, filename, 'glon_bnd', ipts, no_domain=.true.)
      id_restart = register_restart_field(Atm_restart, filename, 'glat_bnd', jpts, no_domain=.true.)
      id_restart = register_restart_field(Atm_restart, filename, 'dt', dto, no_domain=.true.)
   endif

   id_restart = register_restart_field(Til_restart, filename, 'lprec', Atmos % lprec, domain=Atmos%domain)
   id_restart = register_restart_field(Til_restart, filename, 'fprec', Atmos % fprec, domain=Atmos%domain)
   id_restart = register_restart_field(Til_restart, filename, 'gust',  Atmos % gust,  domain=Atmos%domain)
   if (restart_tbot_qbot) then
      id_restart = register_restart_field(Til_restart, filename, 't_bot', Atmos%t_bot, domain=Atmos%domain)
      id_restart = register_restart_field(Til_restart, filename, 'q_bot', Atmos%tr_bot(:,:,ivapor), domain=Atmos%domain)
   end if

   call get_time (Atmos % Time_step, sec, day)
   dt = sec + 86400*day  ! integer seconds

   filename = 'INPUT/atmos_coupled.res.nc'
   if ( file_exist(filename, domain=Atmos%domain) ) then
       if(mpp_pe() == mpp_root_pe() ) call mpp_error ('atmos_model_mod', &
                   'Reading netCDF formatted restart file: INPUT/atmos_coupled.res.nc', NOTE)
       call restore_state(Atm_restart)
       if(in_different_file)call restore_state(Til_restart)
       if (ipts /= mlon .or. jpts /= mlat) call error_mesg &
               ('coupled_atmos_init', 'incorrect resolution on restart file', WARNING)

!---- if the time step has changed then convert ----
!        tendency to conserve mass of water
       if (dto /= dt) then
          Atmos % lprec = Atmos % lprec * real(dto)/real(dt)
          Atmos % fprec = Atmos % fprec * real(dto)/real(dt)
          if (mpp_pe() == mpp_root_pe()) write (logunit,50)
       endif
   else if (file_exist('INPUT/atmos_coupled.res')) then
          if(mpp_pe() == mpp_root_pe() ) call mpp_error ('atmos_model_mod', &
                   'Reading native formatted restart file: INPUT/atmos_coupled.res', NOTE)
          unit = open_restart_file ('INPUT/atmos_coupled.res', 'read')
          !--- check version number (format) of restart file ---
          read  (unit) control
          if (trim(control) /= trim(restart_format)) call error_mesg &
               ('coupled_atmos_init', 'invalid restart format', FATAL)
          !--- check resolution and time step ---
          read  (unit) ipts,jpts,dto
          if (ipts /= mlon .or. jpts /= mlat) call error_mesg &
               ('coupled_atmos_init', 'incorrect resolution on restart file', FATAL)

          !--- read data ---
          call read_data ( unit, Atmos % lprec )
          call read_data ( unit, Atmos % fprec )
          call read_data ( unit, Atmos % gust  )
          if (restart_tbot_qbot) then
             call read_data ( unit, Atmos % t_bot  )
             call read_data ( unit, Atmos % tr_bot(:,:,ivapor) )
          endif
          call close_file (unit)

!---- if the time step has changed then convert ----
!        tendency to conserve mass of water
          if (dto /= dt) then
             Atmos % lprec = Atmos % lprec * real(dto)/real(dt)
             Atmos % fprec = Atmos % fprec * real(dto)/real(dt)
             if (mpp_pe() == mpp_root_pe()) write (logunit,50)
 50         format (/,'The model time step changed .... &
                      &modifying precipitation tendencies')
          endif
   else
        Atmos % lprec = 0.0
        Atmos % fprec = 0.0
        Atmos % gust  = 1.0
   endif

   ! to be written to restart file
   ipts = mlon
   jpts = mlat  
   dto  = dt 

!------ initialize global integral package ------
!**** TEMPORARY FIX FOR GRID CELL CORNER PROBLEM ****

    allocate (area (nlon, nlat))
! call atmosphere_cell_area to obtain array of grid cell areas needed
! by diag_integral_init
    call atmosphere_cell_area (area)
    call diag_integral_init (Atmos % Time_init, Atmos % Time,  &
                             Atmos % lon_bnd(:,:),  &
                             Atmos % lat_bnd(:,:), area)
    deallocate (area)

!-----------------------------------------------------------------------
atmClock = mpp_clock_id( 'Atmosphere', flags=clock_flag_default, grain=CLOCK_COMPONENT )
end subroutine atmos_model_init
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_model_end">
!
! <OVERVIEW>
!  termination routine for atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!  Call once to terminate this module and any other modules used.
!  This routine writes a restart file and deallocates storage
!  used by the derived-type variable atmos_boundary_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_model_end (Atmos)
! </TEMPLATE>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>

subroutine atmos_model_end (Atmos)

type (atmos_data_type), intent(inout) :: Atmos

!-----------------------------------------------------------------------
!---- termination routine for atmospheric model ----
                                              
  call atmosphere_end (Atmos % Time, Atmos%grid)

!------ global integrals ------

  call diag_integral_end (Atmos % Time)

!------ write several atmospheric fields ------
!        also resolution and time step
  call atmos_model_local_restart(Atmos)

!-------- deallocate space --------

  deallocate ( Atmos %  lon_bnd , &
               Atmos %  lat_bnd , &
               Atmos % t_bot    , &
               Atmos % tr_bot   , &
               Atmos % z_bot    , &
               Atmos % p_bot    , &
               Atmos % u_bot    , &
               Atmos % v_bot    , &
               Atmos % p_surf   , &
               Atmos % slp      , &
               Atmos % gust     , &
               Atmos % flux_sw  , &
               Atmos % flux_sw_dir  , &
               Atmos % flux_sw_dif  , &
               Atmos % flux_sw_down_vis_dir  , &
               Atmos % flux_sw_down_vis_dif  , &
               Atmos % flux_sw_down_total_dir  , &
               Atmos % flux_sw_down_total_dif  , &
               Atmos % flux_sw_vis  , &
               Atmos % flux_sw_vis_dir  , &
               Atmos % flux_sw_vis_dif  , &
               Atmos % flux_lw  , &
               Atmos % coszen   , &
               Atmos % lprec    , &
               Atmos % fprec      )

!-----------------------------------------------------------------------

end subroutine atmos_model_end
! </SUBROUTINE>
  !#######################################################################
  ! <SUBROUTINE NAME="atmos_model_restart">
  ! <DESCRIPTION>
  !  Write out restart files registered through register_restart_file
  ! </DESCRIPTION>
  subroutine atmos_model_restart(Atmos, timestamp)
    type (atmos_data_type),   intent(inout) :: Atmos
    character(len=*),  intent(in)           :: timestamp

    call atmosphere_restart(timestamp)
    call atmos_model_local_restart(Atmos, timestamp)

  end subroutine atmos_model_restart
  ! </SUBROUTINE>

  !#######################################################################
  ! <SUBROUTINE NAME="atmos_model_local_restart">
  ! <DESCRIPTION>
  !  Write out restart files registered through register_restart_file
  ! </DESCRIPTION>
  subroutine atmos_model_local_restart(Atmos, timestamp)
    type (atmos_data_type),   intent(inout) :: Atmos
    character(len=*),  intent(in), optional :: timestamp
    integer :: unit
    if( do_netcdf_restart) then
       if(mpp_pe() == mpp_root_pe()) then
          call mpp_error ('atmos_model_mod', 'Writing netCDF formatted restart file.', NOTE)
       endif
       call save_restart(Atm_restart, timestamp)
       if(in_different_file) call save_restart(Til_restart, timestamp)
    else
       if(present(timestamp)) call mpp_error ('atmos_model_mod',  &
            'intermediate restart capability is not implemented for non-netcdf file', FATAL)
       unit = open_restart_file ('RESTART/atmos_coupled.res', 'write')
       if (mpp_pe() == mpp_root_pe()) then
          write (unit) restart_format
          write (unit) ipts, jpts, dto
       endif
       call write_data ( unit, Atmos % lprec )
       call write_data ( unit, Atmos % fprec )
       call write_data ( unit, Atmos % gust  )
       if(restart_tbot_qbot) then
          call write_data ( unit, Atmos % t_bot  )
          call write_data ( unit, Atmos % tr_bot(:,:,ivapor)  )
       endif
       call close_file (unit)
    endif

  end subroutine atmos_model_local_restart
  ! </SUBROUTINE>
!#######################################################################
! <SUBROUTINE NAME="atm_stock_pe">
!
! <OVERVIEW>
!  returns the total stock in atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!  Called to compute and return the total stock (e.g., water, heat, etc.)
! in the atmospheric on the current PE.
! </DESCRIPTION>

! <TEMPLATE>
!   call atm_stock_pe (Atmos, index, value)
! </TEMPLATE>

! <INOUT NAME="Atm" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>
!
! <IN NAME="index" TYPE="integer">
!   Index of stock to be computed.
! </IN>
!
! <OUT NAME="value" TYPE="real">
!   Value of stock on the current processor.
! </OUT>

subroutine atm_stock_pe (Atm, index, value)

type (atmos_data_type), intent(inout) :: Atm
integer,                intent(in)    :: index
real,                   intent(out)   :: value

   value = 0.0
   if(Atm%pe) call get_stock_pe (index, value)

end subroutine atm_stock_pe

! </SUBROUTINE>

!#######################################################################
!#######################################################################
! <SUBROUTINE NAME="atmos_data_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the atmos_data_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the atmos_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_data_type_chksum(id, timestep, atm)
! </TEMPLATE>

! <IN NAME="Atm" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields in the atmos_data_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!
subroutine atmos_data_type_chksum(id, timestep, atm)
type(atmos_data_type), intent(in) :: atm 
    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    integer :: n, outunit

100 format("CHECKSUM::",A32," = ",Z20)
101 format("CHECKSUM::",A16,a,'%',a," = ",Z20)

  outunit = stdout()
  write(outunit,*) 'BEGIN CHECKSUM(Atmos_data_type):: ', id, timestep
  write(outunit,100) ' atm%lon_bnd                ', mpp_chksum(atm%lon_bnd               )
  write(outunit,100) ' atm%lat_bnd                ', mpp_chksum(atm%lat_bnd               )
  write(outunit,100) ' atm%t_bot                  ', mpp_chksum(atm%t_bot                 )
  do n = 1, size(atm%tr_bot,3)
  write(outunit,100) ' atm%tr_bot(:,:,n)          ', mpp_chksum(atm%tr_bot(:,:,n)         )
  enddo
  write(outunit,100) ' atm%z_bot                  ', mpp_chksum(atm%z_bot                 )
  write(outunit,100) ' atm%p_bot                  ', mpp_chksum(atm%p_bot                 )
  write(outunit,100) ' atm%u_bot                  ', mpp_chksum(atm%u_bot                 )
  write(outunit,100) ' atm%v_bot                  ', mpp_chksum(atm%v_bot                 )
  write(outunit,100) ' atm%p_surf                 ', mpp_chksum(atm%p_surf                )
  write(outunit,100) ' atm%slp                    ', mpp_chksum(atm%slp                   )
  write(outunit,100) ' atm%gust                   ', mpp_chksum(atm%gust                  )
  write(outunit,100) ' atm%coszen                 ', mpp_chksum(atm%coszen                )
  write(outunit,100) ' atm%flux_sw                ', mpp_chksum(atm%flux_sw               )
  write(outunit,100) ' atm%flux_sw_dir            ', mpp_chksum(atm%flux_sw_dir           )
  write(outunit,100) ' atm%flux_sw_dif            ', mpp_chksum(atm%flux_sw_dif           )
  write(outunit,100) ' atm%flux_sw_down_vis_dir   ', mpp_chksum(atm%flux_sw_down_vis_dir  )
  write(outunit,100) ' atm%flux_sw_down_vis_dif   ', mpp_chksum(atm%flux_sw_down_vis_dif  )
  write(outunit,100) ' atm%flux_sw_down_total_dir ', mpp_chksum(atm%flux_sw_down_total_dir)
  write(outunit,100) ' atm%flux_sw_down_total_dif ', mpp_chksum(atm%flux_sw_down_total_dif)
  write(outunit,100) ' atm%flux_sw_vis            ', mpp_chksum(atm%flux_sw_vis           )
  write(outunit,100) ' atm%flux_sw_vis_dir        ', mpp_chksum(atm%flux_sw_vis_dir       )
  write(outunit,100) ' atm%flux_sw_vis_dif        ', mpp_chksum(atm%flux_sw_vis_dif       )
  write(outunit,100) ' atm%flux_lw                ', mpp_chksum(atm%flux_lw               )
  write(outunit,100) ' atm%lprec                  ', mpp_chksum(atm%lprec                 )
  write(outunit,100) ' atm%fprec                  ', mpp_chksum(atm%fprec                 )
!  call surf_diff_type_chksum(id, timestep, atm%surf_diff)

end subroutine atmos_data_type_chksum

! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="lnd_ice_atm_bnd_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the land_ice_atmos_boundary_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the land_ice_atmos_boundary_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_data_type_chksum(id, timestep, bnd_type)
! </TEMPLATE>

! <IN NAME="bnd_type" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains fields in the land_ice_atmos_boundary_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!


subroutine lnd_ice_atm_bnd_type_chksum(id, timestep, bnd_type)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(land_ice_atmos_boundary_type), intent(in) :: bnd_type
 integer ::   n, outunit

    outunit = stdout()
    write(outunit,*) 'BEGIN CHECKSUM(lnd_ice_Atm_bnd_type):: ', id, timestep
100 format("CHECKSUM::",A32," = ",Z20)
    write(outunit,100) 'lnd_ice_atm_bnd_type%t             ',mpp_chksum(bnd_type%t              )
    write(outunit,100) 'lnd_ice_atm_bnd_type%albedo        ',mpp_chksum(bnd_type%albedo         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%albedo_vis_dir',mpp_chksum(bnd_type%albedo_vis_dir )
    write(outunit,100) 'lnd_ice_atm_bnd_type%albedo_nir_dir',mpp_chksum(bnd_type%albedo_nir_dir )
    write(outunit,100) 'lnd_ice_atm_bnd_type%albedo_vis_dif',mpp_chksum(bnd_type%albedo_vis_dif )
    write(outunit,100) 'lnd_ice_atm_bnd_type%albedo_nir_dif',mpp_chksum(bnd_type%albedo_nir_dif )
    write(outunit,100) 'lnd_ice_atm_bnd_type%land_frac     ',mpp_chksum(bnd_type%land_frac      )
    write(outunit,100) 'lnd_ice_atm_bnd_type%dt_t          ',mpp_chksum(bnd_type%dt_t           )
    do n = 1, size(bnd_type%dt_tr,3)
    write(outunit,100) 'lnd_ice_atm_bnd_type%dt_tr(:,:,n)  ',mpp_chksum(bnd_type%dt_tr(:,:,n)   )
    enddo
    write(outunit,100) 'lnd_ice_atm_bnd_type%u_flux        ',mpp_chksum(bnd_type%u_flux         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%v_flux        ',mpp_chksum(bnd_type%v_flux         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%dtaudu        ',mpp_chksum(bnd_type%dtaudu         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%dtaudv        ',mpp_chksum(bnd_type%dtaudv         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%u_star        ',mpp_chksum(bnd_type%u_star         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%b_star        ',mpp_chksum(bnd_type%b_star         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%q_star        ',mpp_chksum(bnd_type%q_star         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%rough_mom     ',mpp_chksum(bnd_type%rough_mom      )
!    write(outunit,100) 'lnd_ice_atm_bnd_type%data          ',mpp_chksum(bnd_type%data           )

end subroutine lnd_ice_atm_bnd_type_chksum
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="lnd_atm_bnd_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the land_atmos_boundary_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the land_atmos_boundary_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call lnd_atm_bnd_type_chksum(id, timestep, bnd_type)
! </TEMPLATE>

! <IN NAME="bnd_type" TYPE="type(land_atmos_boundary_type)">
!   Derived-type variable that contains fields in the land_atmos_boundary_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!


subroutine lnd_atm_bnd_type_chksum(id, timestep, bnd_type)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(land_atmos_boundary_type), intent(in) :: bnd_type
 integer ::   n, outunit

    outunit = stdout()
    write(outunit,*) 'BEGIN CHECKSUM(lnd_atmos_boundary_type):: ', id, timestep
!    write(outunit,100) 'lnd_atm_bnd_type%data',mpp_chksum(bnd_type%data)

100 format("CHECKSUM::",A32," = ",Z20)

end subroutine lnd_atm_bnd_type_chksum
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="ice_atm_bnd_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the ice_atmos_boundary_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the ice_atmos_boundary_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call ice_atm_bnd_type_chksum(id, timestep, bnd_type)
! </TEMPLATE>

! <IN NAME="bnd_type" TYPE="type(ice_atmos_boundary_type)">
!   Derived-type variable that contains fields in the ice_atmos_boundary_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!


subroutine ice_atm_bnd_type_chksum(id, timestep, bnd_type)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ice_atmos_boundary_type), intent(in) :: bnd_type
 integer ::   n, outunit

    outunit = stdout()
    write(outunit,*) 'BEGIN CHECKSUM(ice_atmos_boundary_type):: ', id, timestep
!    write(outunit,100) 'ice_atm_bnd_type%data',mpp_chksum(data_type%data)

100 format("CHECKSUM::",A32," = ",Z20)


end subroutine ice_atm_bnd_type_chksum
! </SUBROUTINE>


end module atmos_model_mod


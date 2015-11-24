module ocean_nphysics_new_mod

!<CONTACT EMAIL="tim.leslie@gmail.com"> Tim Leslie
!</CONTACT>
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Compute the effects of neutral physics processes 
! (neutral diffusion and neutral skew-diffusion).
!</OVERVIEW>
!
!<DESCRIPTION>
! Compute the effects of neutral physics processes 
! (neutral diffusion and neutral skew-diffusion).
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies 
! Fundamentals of Ocean Climate Models (FOCM) (2004)
! Princeton University Press 
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, Elements of MOM (2012)
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_nphysics_new_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be set .true. to use this module.
!  Default use_this_module = .false.
!  </DATA> 
!
!  <DATA NAME="drhodz_smooth_vert" TYPE="logical">
!  For smoothing vertical density gradient before computing 
!  neutral slope.  Exercise caution if using this option.   
!  Default drhodz_smooth_vert=.false.
!  </DATA> 
!
!  <DATA NAME="drhodz_smooth_horz" TYPE="logical">
!  For smoothing vertical density gradient before computing 
!  neutral slope.  Exercise caution if using this option.   
!  Default drhodz_smooth_horz=.false.
!  </DATA> 
!
!  <DATA NAME="smax" TYPE="real">
!  Slope maximum parameter for setting behaviour of neutral 
!  physics. Default smax=0.01.
!  </DATA> 
!
!  <DATA NAME="vel_micom_smooth" TYPE="real" UNITS="m/s">
!  For horizontal smoothing of drhodz before computing neutral
!  slopes. Default vel_micom_smooth=0.2.
!  </DATA> 
!
!</NAMELIST>

#define COMP    isc:iec,jsc:jec
#define COMPXLL isc-1:iec-1,jsc:jec
#define COMPYLL isc:iec,jsc-1:jec-1

  use constants_mod,   only: epsln, pi
  use fms_mod,         only: FATAL, WARNING, open_namelist_file, check_nml_error 
  use fms_mod,         only: close_file, file_exist, write_version_number
  use fms_io_mod,      only: restart_file_type
  use fms_io_mod,      only: register_restart_field, save_restart, restore_state
  use mpp_mod,         only: mpp_error, stdout, stdlog, mpp_clock_id, input_nml_file
  use mpp_mod,         only: CLOCK_ROUTINE, mpp_clock_begin, mpp_clock_end
  use mpp_domains_mod, only: mpp_update_domains

  use ocean_domains_mod,           only: get_local_indices, set_ocean_domain
  use ocean_nphysics_util_new_mod, only: ocean_nphysics_util_new_init, vert_smooth, compute_rossby_radius
  use ocean_nphysics_flux_mod,     only: ocean_nphysics_flux_init, flux_calculations
  use ocean_nphysics_tensor_mod,   only: ocean_nphysics_tensor_init, compute_tensors
  use ocean_nphysics_skew_mod,     only: ocean_nphysics_skew_init
  use ocean_nphysics_diff_mod,     only: ocean_nphysics_diff_init, compute_diffusivity
  use ocean_operators_mod,         only: FDX_T, FDY_T, FDX_ZT, FDY_ZT, FMX, FMY, LAP_T
  use ocean_parameters_mod,        only: rho0, grav, TERRAIN_FOLLOWING
  use ocean_types_mod,             only: ocean_grid_type, ocean_domain_type, ocean_time_type
  use ocean_types_mod,             only: ocean_time_steps_type, ocean_options_type
  use ocean_types_mod,             only: ocean_thickness_type, ocean_prog_tracer_type, ocean_density_type
  use ocean_util_mod,              only: write_note
  use ocean_util_mod,              only: diagnose_2d, diagnose_2d_int, diagnose_3d, diagnose_3d_int
  use ocean_util_mod,              only: register_2d_t_field, register_3d_t_field
  use ocean_util_mod,              only: register_3d_xte_field, register_3d_ytn_field, register_3d_ztb_field

  implicit none

  public ocean_nphysics_new_init
  public neutral_physics_new
  public ocean_nphysics_new_restart
  public ocean_nphysics_new_end

  private tracer_gradients

  private gradrho
  private adjust_drhodz

  private neutral_slopes
  private neutral_blayer

  private

#include <ocean_memory.h>

character(len=128) :: version=&
     '$Id: ocean_nphysics_new.F90,v 20.0 2013/12/14 00:14:46 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

character(len=*), parameter :: FILENAME=&
     __FILE__


  logical :: module_is_initialized = .FALSE.

  integer :: id_clock_tracer_gradients
  integer :: id_clock_density
  integer :: id_clock_diffusivity
  integer :: id_clock_tensors
  integer :: id_clock_flux

  integer :: id_drhodx
  integer :: id_drhody
  integer :: id_drhodz
  integer :: id_drhodz_up
  integer :: id_drhodz_lo

  integer :: id_slopex
  integer :: id_slopey
  integer :: id_absslope
  integer :: id_N2
  integer :: id_gravity_wave_speed
  integer :: id_rossby_radius
  integer :: id_ksurf_blayer

  integer, dimension(:), allocatable :: id_dTdx
  integer, dimension(:), allocatable :: id_dTdy
  integer, dimension(:), allocatable :: id_dTdz

  integer :: id_neut_rho_nphysics
  integer :: id_wdian_rho_nphysics

  type(restart_file_type), save :: nphysics_new_restart

  real, dimension(:,:), allocatable :: grid_length ! grid length array (metre)
  real, dimension(:,:), allocatable :: smooth_diff ! (m^2/sec) for smoothing

  ! Diffusivity arrays (m^2/s)
  real, dimension(:,:), allocatable :: ah_array
  real, dimension(:,:,:), allocatable :: agm_array
  real, dimension(:,:,:), allocatable :: aredi_array

  real :: dtime
  integer :: index_temp, index_salt

  logical :: use_this_module

  ! for Tim%init=.true. and ocean_neutral.res.nc exists, but still wish
  ! to start from initial fields.
  logical :: nphysics_new_zero_init=.false.

  real :: smax = 0.01

  ! adjust_drhodz()
  logical :: drhodz_smooth_vert = .false.
  logical :: drhodz_smooth_horz = .false.
  real    :: vel_micom_smooth   = 0.2 ! m/sec

  namelist /ocean_nphysics_new_nml/ use_this_module, &
       drhodz_smooth_vert, drhodz_smooth_horz, smax, vel_micom_smooth

contains

  !#######################################################################
  ! <SUBROUTINE NAME="ocean_nphysics_new_init">
  !
  ! <DESCRIPTION>
  ! Initialises diagnostics, namelists and constants.
  ! </DESCRIPTION>
  !
  subroutine ocean_nphysics_new_init(Grid, Domain, Time, Time_steps, Thickness, &
                                     T_prog, Ocean_options, vert_coordinate_type)

    type(ocean_grid_type),                      intent(in)    :: Grid
    type(ocean_domain_type),                    intent(inout) :: Domain
    type(ocean_time_type),                      intent(in)    :: Time
    type(ocean_time_steps_type),                intent(in)    :: Time_steps
    type(ocean_thickness_type),                 intent(in)    :: Thickness
    type(ocean_prog_tracer_type), dimension(:), intent(in)    :: T_prog
    type(ocean_options_type),                   intent(inout) :: Ocean_options
    integer,                                    intent(in)    :: vert_coordinate_type

    logical :: agm_z_dep
    integer :: ioun, ierr, io_status
    integer :: num_prog_tracers, n
    integer :: id_restart

    character(len=*), parameter :: restart_filename = 'ocean_neutral_new.res.nc'

    integer :: stdoutunit, stdlogunit
    stdoutunit = stdout()
    stdlogunit = stdlog()

    num_prog_tracers = size(T_prog)

    if ( module_is_initialized ) then 
       call mpp_error(FATAL, &
            '==>Error from ocean_nphysics_mod (ocean_nphysics_init):already initialized')
    endif

    module_is_initialized = .TRUE.

    call write_version_number(version, tagname)

    ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_nphysics_new_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_nphysics_new_nml')
#else
    ioun =  open_namelist_file()
    read (ioun,ocean_nphysics_new_nml,IOSTAT=io_status)
    ierr = check_nml_error(io_status,'ocean_nphysics_new_nml')
    call close_file (ioun)
#endif
    write (stdoutunit,'(/)')
    write (stdoutunit,ocean_nphysics_new_nml)  
    write (stdlogunit,ocean_nphysics_new_nml)

  call ocean_nphysics_util_new_init(Domain, Grid)

  if(use_this_module) then 
     call write_note(FILENAME,&
     'USING ocean_nphysics_new.')
     if(vert_coordinate_type==TERRAIN_FOLLOWING) then 
        call mpp_error(WARNING, &
             '==>Warning: ocean_nphysics_new is NOT supported with TERRRAIN_FOLLOWING vertical coordinates.')
     endif
     Ocean_options%neutral_physics_new = 'Used neutral physics new.'
  else 
     call write_note(FILENAME,&
     'NOT using ocean_nphysics_new.')
     Ocean_options%neutral_physics_new = 'Did NOT use neutral physics new option.'
     return
  endif 


#ifndef MOM_STATIC_ARRAYS
    call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
    nk = Grid%nk
#endif

    index_temp=-1;index_salt=-1
    do n=1,num_prog_tracers
       if (T_prog(n)%name == 'temp') index_temp = n
       if (T_prog(n)%name == 'salt') index_salt = n
    enddo
    if (index_temp == -1 .or. index_salt == -1) then
       call mpp_error(FATAL, &
            '==>Error: temp and/or salt not identified in call to ocean_nphysics_util_init')
    endif

    id_clock_tracer_gradients = mpp_clock_id('(Ocean neutral: tracer gradients)' ,grain=CLOCK_ROUTINE)
    id_clock_density          = mpp_clock_id('(Ocean neutral: density calcs)'    ,grain=CLOCK_ROUTINE)
    id_clock_diffusivity      = mpp_clock_id('(Ocean neutral: diffusivity)'    ,grain=CLOCK_ROUTINE)
    id_clock_tensors          = mpp_clock_id('(Ocean neutral: tensors)'    ,grain=CLOCK_ROUTINE)
    id_clock_flux             = mpp_clock_id('(Ocean neutral: flux)'    ,grain=CLOCK_ROUTINE)

    ! some useful constants
    dtime = Time_steps%dtime_t

    allocate (grid_length(isd:ied,jsd:jed))
    allocate (smooth_diff(isd:ied,jsd:jed))
    allocate (ah_array(isd:ied,jsd:jed))
    allocate (agm_array(isd:ied,jsd:jed,nk))
    allocate (aredi_array(isd:ied,jsd:jed,nk))

    grid_length(:,:) = 2.0*Grid%dxt(:,:)*Grid%dyt(:,:)/(Grid%dxt(:,:)+Grid%dyt(:,:)) 
    smooth_diff(:,:) = vel_micom_smooth*grid_length(:,:)

    call ocean_nphysics_diff_init(Domain, Grid, Time, Thickness, grid_length, dtime, smax, &
         ah_array, agm_array, aredi_array, agm_z_dep)
    call ocean_nphysics_flux_init(Domain, Grid, Time, T_prog)
    call ocean_nphysics_tensor_init(Domain, Grid, Time)
    call ocean_nphysics_skew_init(Domain, Grid, Time, agm_z_dep)

    id_restart = register_restart_field(nphysics_new_restart, restart_filename, 'agm_array', agm_array, &
         domain=Domain%domain2d, mandatory=.false.)
    id_restart = register_restart_field(nphysics_new_restart, restart_filename, 'aredi_array', aredi_array, &
         domain=Domain%domain2d, mandatory=.false.)

    if (Time%init .and. nphysics_new_zero_init) then
       !write(stdoutunit,'(1x,a)') ' Starting ocean_nphysics_util fields from raw initialization.'
    elseif (.not. file_exist('INPUT/' // restart_filename)) then
       if (.not. Time%init) call mpp_error(FATAL,'Expecting file INPUT/ocean_neutral_new.res.nc to exist.&
            &This file was not found and Time%init=.false.')
    else
       call restore_state(nphysics_new_restart)

       call mpp_update_domains(agm_array, Domain%domain2d)
       call mpp_update_domains(aredi_array, Domain%domain2d)
    endif

    call register_fields(Grid, Time, T_prog)

  end subroutine ocean_nphysics_new_init
  ! </SUBROUTINE> NAME="ocean_nphysics_new_init"

  !#######################################################################
  ! <SUBROUTINE NAME="register fields">
  !
  ! <DESCRIPTION>
  ! Register diagnostic fields.
  ! </DESCRIPTION>
  !
  subroutine register_fields(Grid, Time, T_prog)
    type(ocean_grid_type), intent(in) :: Grid
    type(ocean_time_type), intent(in) :: Time
    type(ocean_prog_tracer_type), dimension(:), intent(in) :: T_prog

    integer n, num_prog_tracers

    num_prog_tracers = size(T_prog(:))

    id_drhodx = register_3d_t_field(Grid, Time, 'drhodx',                           &
         'drho/dx at tracer point calculated as average of left/right components ', &
         'rho/m', (/-1e-4,1.e-4/))
    id_drhody = register_3d_t_field(Grid, Time, 'drhody',                           &
         'drho/dy at tracer point calculated as average of left/right components ', &
         'rho/m', (/-1e-4,1.e-4/))
    id_drhodz = register_3d_t_field(Grid, Time, 'drhodz',                        &
         'drho/dz at tracer point calculated as average of up/down components ', &
         'rho/m', (/-1.e-1, 0.0/))
    id_drhodz_up = register_3d_t_field(Grid, Time, 'drhodz_up', &
         'drho/dz at tracer point up component ', 'rho/m', (/-1.e-1, 0.0/))
    id_drhodz_lo = register_3d_t_field(Grid, Time, 'drhodz_lo', &
         'drho/dz at tracer point lo component ', 'rho/m', (/-1.e-1, 0.0/))

    id_slopex = register_3d_t_field(Grid, Time, 'slopex', &
         'x component of neutral slope averaged over triads', '', (/-1e-1,1.e-1/))
    id_slopey = register_3d_t_field(Grid, Time, 'slopey', &
         'y component of neutral slope averaged over triads', '', (/-1e-1,1.e-1/))
    id_absslope = register_3d_t_field(Grid, Time, 'absslope', &
         'magnitude of the neutral slope vector', '', (/0.0,1.-1/))

    id_N2 = register_3d_t_field(Grid, Time, 'N2', &
         'Squared neutral buoyancy frequency', '1/s^2', (/0.0,1.e8/))
    id_gravity_wave_speed = register_2d_t_field(Grid, Time, 'gravity_wave_speed', &
         'Gravity wave speed.', 'm/s', (/-1e-8,1.e8/))

    id_rossby_radius = register_2d_t_field(Grid, Time, 'rossby_radius', &
         'Rossby radius', 'm', (/-1e-8,1.e8/))

    id_ksurf_blayer = register_2d_t_field(Grid, Time, 'ksurf_blayer', &
         'K-index of surface boundary layer.', 'k', (/0.0,100.0/))

    allocate(id_dTdx(num_prog_tracers))
    allocate(id_dTdy(num_prog_tracers))
    allocate(id_dTdz(num_prog_tracers))

    do n = 1, num_prog_tracers
       id_dTdx(n) = register_3d_xte_field(Grid, Time, 'dTdx_'//trim(T_prog(n)%name), &
            'Gradient of '//trim(T_prog(n)%name)//' in x direction', &
            '['//trim(T_prog(n)%name)//']/m', (/-1.e18,1.e18/))
       id_dTdy(n) = register_3d_ytn_field(Grid, Time, 'dTdy_'//trim(T_prog(n)%name), &
            'Gradient of '//trim(T_prog(n)%name)//' in y direction', &
            '['//trim(T_prog(n)%name)//']/m', (/-1.e18,1.e18/))
       id_dTdz(n) = register_3d_ztb_field(Grid, Time, 'dTdz_'//trim(T_prog(n)%name), &
            'Gradient of '//trim(T_prog(n)%name)//' in z direction', &
            '['//trim(T_prog(n)%name)//']/m', (/-1.e18,1.e18/))
    enddo

    id_neut_rho_nphysics = register_3d_t_field(Grid, Time, 'neut_rho_nphysics',&
         'update of neutral density from explicit in time nphysics',           &
         'rho*rho_dz/sec', (/-1.e10,1.e10/))
    id_wdian_rho_nphysics = register_3d_t_field(Grid, Time, 'wdian_rho_nphysics',&
         'dianeutral velocity component due to explicit in time nphysics',       &
         'm/sec', (/-1.e10,1.e10/))

  end subroutine register_fields
  ! </SUBROUTINE> NAME="register_fields"


  !#######################################################################
  ! <SUBROUTINE NAME="neutral_physics_new">
  !
  ! <DESCRIPTION>
  ! This public interface computes the effects of neutral physics
  ! processes (neutral diffusion and neutral skew-diffusion).
  !
  ! --Returns rossby_radius for use in lapgen_friction module.
  ! --Returns gm_diffusivity for use in visc_form_drag module
  ! --Updates the value of T_prog(n)%th_tendency with the explicit thickness
  !   weighted tendency due to neutral physics and also T_prog(n)%K33_implicit
  !   with the implicit component.
  ! </DESCRIPTION>
  !
  subroutine neutral_physics_new(Time, Domain, Thickness, Dens, Grid, surf_blthick, &
       T_prog, rossby_radius, gm_diffusivity)

    type(ocean_time_type),      intent(in)    :: Time
    type(ocean_domain_type),    intent(inout) :: Domain
    type(ocean_thickness_type), intent(in)    :: Thickness
    type(ocean_density_type),   intent(in)    :: Dens
    type(ocean_grid_type),      intent(in)    :: Grid
    real, dimension(isd:,jsd:), intent(in)    :: surf_blthick

    type(ocean_prog_tracer_type), dimension(:), intent(inout) :: T_prog
    real, dimension(isd:,jsd:),                 intent(out) :: rossby_radius
    real, dimension(isd:,jsd:,:),               intent(out) :: gm_diffusivity

    real, dimension(isd:ied,jsd:jed,nk,size(T_prog)) :: dTdx, dTdy
    real, dimension(isd:ied,jsd:jed,0:nk,size(T_prog)) :: dTdz

    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: slopex_xz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: slopey_yz
    real, dimension(isd:ied,jsd:jed,nk,0:1) :: Ndz_z
    real, dimension(isd:ied,jsd:jed,nk,0:1) :: N2_z
    real, dimension(isd:ied,jsd:jed) :: gravity_wave_speed

    real, dimension(isd:ied,jsd:jed,nk,0:1) :: absslope_z
    real, dimension(isd:ied,jsd:jed,nk) :: absslope
    integer, dimension(isd:ied,jsd:jed) :: ksurf_blayer   ! k-value at base of surface nblayer

    ! Neutral diffusion tensor components [m^2/s]
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: symm_tensor13_xz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: symm_tensor23_yz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: symm_tensor33_xz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: symm_tensor33_yz
    real, dimension(isd:ied,jsd:jed,nk,0:1)     :: symm_tensor1122_z

    ! Skew diffusion tensor components (transport vector) [m^2/s]
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: upsilonx_xz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: upsilony_yz

    real, dimension(isd:ied,jsd:jed,nk,size(T_prog)) :: total_th_tendency

    integer :: n, num_prog_tracers

    if (.not. use_this_module) return

    if ( .not. module_is_initialized ) then 
       call mpp_error(FATAL, &
            '==>Error from ocean_nphysics_new (neutral_physics_new): needs initialization')
    endif

    num_prog_tracers = size(T_prog)

    call mpp_clock_begin(id_clock_tracer_gradients)
    call tracer_gradients(Time, Grid, Time%taum1, T_prog, Thickness, &
         dTdx, dTdy, dTdz)
    call mpp_clock_end(id_clock_tracer_gradients)   

    call mpp_clock_begin(id_clock_density)
    call density_calculations(Grid, Time, Thickness, Dens, dTdx, dTdy, dTdz, smax,    &
         slopex_xz, slopey_yz, absslope_z, absslope, N2_z, Ndz_z, gravity_wave_speed, &
         rossby_radius, ksurf_blayer)
    call mpp_clock_end(id_clock_density)

    call mpp_clock_begin(id_clock_diffusivity)
    call compute_diffusivity(Domain, Grid, Time, Thickness, Dens, T_prog(index_temp), &
         T_prog(index_salt), absslope_z, N2_z, gravity_wave_speed, rossby_radius, &
         ksurf_blayer, grid_length, agm_array, aredi_array)
    call mpp_clock_end(id_clock_diffusivity)

    call mpp_clock_begin(id_clock_tensors)
    call compute_tensors(Time, Domain, Grid, Thickness, slopex_xz, slopey_yz, Ndz_z, &
         N2_z, gravity_wave_speed, rossby_radius, absslope, absslope_z,              &
         aredi_array, ah_array, agm_array, ksurf_blayer, smax, surf_blthick,         &
         symm_tensor13_xz, symm_tensor23_yz, symm_tensor33_xz, symm_tensor33_yz,     &
         symm_tensor1122_z, upsilonx_xz, upsilony_yz)
    call mpp_clock_end(id_clock_tensors)

    call mpp_clock_begin(id_clock_flux)
    call flux_calculations(Domain, Grid, Time, Thickness, dtime,        &
         symm_tensor13_xz, symm_tensor23_yz, symm_tensor33_xz,          &
         symm_tensor33_yz, symm_tensor1122_z, upsilonx_xz, upsilony_yz, &
         dTdx, dTdy, dTdz, aredi_array,                                 &
         total_th_tendency, T_prog)
    call mpp_clock_end(id_clock_flux)

    do n = 1, num_prog_tracers
       T_prog(n)%th_tendency(COMP,:) = T_prog(n)%th_tendency(COMP,:) + total_th_tendency(COMP,:,n)
    enddo

    call diagnose_3d(Time, Grid, id_neut_rho_nphysics,              &
         (Dens%drhodT(COMP,:)*total_th_tendency(COMP,:,index_temp)  &
         +Dens%drhodS(COMP,:)*total_th_tendency(COMP,:,index_salt)))
    call diagnose_3d(Time, Grid, id_wdian_rho_nphysics,             &
         (Dens%drhodT(COMP,:)*total_th_tendency(COMP,:,index_temp)  &
         +Dens%drhodS(COMP,:)*total_th_tendency(COMP,:,index_salt)) &
         /(epsln+Dens%drhodz_zt(COMP,:)*Thickness%rho_dzt(COMP,:,Time%tau)))


    ! gm_diffusivity passed to neutral_physics for computing form drag viscosity
    gm_diffusivity(:,:,:) = agm_array(:,:,:)

  end subroutine neutral_physics_new
  ! </SUBROUTINE> NAME="neutral_physics_new"

  !#######################################################################
  ! <SUBROUTINE NAME="tracer_gradients">
  !
  ! <DESCRIPTION>
  ! Compute the tracer derivatives.
  !
  ! G(2004) 16.58
  ! dTdx(i) = (T(i+1) - T(i))/dxte
  ! dTdT(j) = (T(j+1) - T(j))/dytn
  ! dTdz(k) = (T(k) - T(k+1))/dzwt
  !
  ! Horizontal derivatives are taken along surfaces of 
  ! constant vertical coordinate (constant k-level)
  !
  ! This approach ensures that when neutral physics defaults to "horizontal" physics
  ! next to boundaries, it will do so as horizontal, defined along surfaces of constant 
  ! s-surfaces, and so will not generate spurious extrema.  
  !
  ! Additionally, when using generalized vertical coordinates, the neutral diffusion
  ! slope should be computed relative to the s-surfaces.  The skew diffusion slope 
  ! should ideally be computed with respect to z-surfaces, as z-surfaces define
  ! available potential energy. However, when s and z surfaces are reasonably close, 
  ! as they are in the interior for zstar and pstar vertical coordinates, then we 
  ! choose to to dissipate thickness as defined relative to the zstar or pstar surfaces. 
  ! This should not be such a big deal, and it is certainly easier computationally than
  ! worrying about computing two separate sets of slopes.  More on this detail is 
  ! discussed in "Elements of MOM".
  ! 
  ! NOTE: This approach is not appropriate for sigma-models. Indeed, many assumptions
  ! in the neutral physics modules need to be rethought for terrain following vertical
  ! coordinates.  
  !
  ! </DESCRIPTION>
  !
  subroutine tracer_gradients(Time, Grid, taum1, T_prog, Thickness, &
       dTdx, dTdy, dTdz)

    type(ocean_time_type),                      intent(in) :: Time
    type(ocean_grid_type),                      intent(in) :: Grid
    integer,                                    intent(in) :: taum1
    type(ocean_thickness_type),                 intent(in) :: Thickness
    type(ocean_prog_tracer_type), dimension(:), intent(in) :: T_prog

    real, dimension(isd:,jsd:,:,:),  intent(out) :: dTdx
    real, dimension(isd:,jsd:,:,:),  intent(out) :: dTdy
    real, dimension(isd:,jsd:,0:,:), intent(out) :: dTdz

    integer :: k, n, num_prog_tracers
    real, dimension(isd:ied,jsd:jed,2) :: field

    num_prog_tracers = size(T_prog)

    ! lateral derivatives taken along surfaces of
    ! constant vertical coordinate (constant k-level)
    do n=1,num_prog_tracers
       do k = 1, nk
          field(:,:,1) = T_prog(n)%field(:,:,k,taum1)
          dTdx(:,:,k,n) = FDX_T(field(:,:,1))*FMX(Grid%tmask(:,:,k))
          dTdy(:,:,k,n) = FDY_T(field(:,:,1))*FMY(Grid%tmask(:,:,k))
       enddo
    enddo
 
    ! vertical derivative
    do n=1,num_prog_tracers
       dTdz(:,:,0,n) = 0.0
       do k = 1, nk-1
          dTdz(:,:,k,n) = Grid%tmask(:,:,k+1)*(T_prog(n)%field(:,:,k,taum1) - T_prog(n)%field(:,:,k+1,taum1)) &
                          /(Thickness%dzwt(:,:,k) + epsln)
       enddo
       dTdz(:,:,nk,n) = 0.0
    enddo

    do n=1,num_prog_tracers
       call diagnose_3d(Time, Grid, id_dTdx(n), dTdx(:,:,:,n))
       call diagnose_3d(Time, Grid, id_dTdy(n), dTdy(:,:,:,n))
       call diagnose_3d(Time, Grid, id_dTdz(n), dTdz(:,:,1:nk,n))
    enddo

  end subroutine tracer_gradients
  ! </SUBROUTINE> NAME="tracer_gradients"

  !#######################################################################
  ! <SUBROUTINE NAME="density_calculations">
  !
  ! <DESCRIPTION>
  ! This subroutine computes a number of values based on the density gradient.
  !
  ! These are the neutral slope vector, the neutral buoyancy frequency, 
  ! gravity wave speed, Rossby radius and boundary layer depth.
  ! </DESCRIPTION>
  !
  subroutine density_calculations(Grid, Time, Thickness, Dens, dTdx, dTdy, dTdz, smax, &
       slopex_xz, slopey_yz, absslope_z, absslope, N2_z, Ndz_z, gravity_wave_speed,    &
       rossby_radius, ksurf_blayer)

    type(ocean_density_type),        intent(in) :: Dens
    type(ocean_time_type),           intent(in) :: Time
    type(ocean_grid_type),           intent(in) :: Grid
    type(ocean_thickness_type),      intent(in) :: Thickness
    real, dimension(isd:,jsd:,:,:),  intent(in) :: dTdx
    real, dimension(isd:,jsd:,:,:),  intent(in) :: dTdy
    real, dimension(isd:,jsd:,0:,:), intent(in) :: dTdz
    real,                            intent(in) :: smax

    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: slopex_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: slopey_yz
    real, dimension(isd:,jsd:,:,0:),    intent(out) :: absslope_z
    real, dimension(isd:,jsd:,:),       intent(out) :: absslope
    real, dimension(isd:,jsd:,:,0:),    intent(out) :: N2_z
    real, dimension(isd:,jsd:,:,0:),    intent(out) :: Ndz_z
    real, dimension(isd:,jsd:),         intent(out) :: gravity_wave_speed
    real, dimension(isd:,jsd:),         intent(out) :: rossby_radius

    integer, dimension(isd:,jsd:), intent(out) :: ksurf_blayer

    ! neutral density derivatives 
    real, dimension(isd:ied,jsd:jed,nk,0:1) :: drhodx_x
    real, dimension(isd:ied,jsd:jed,nk,0:1) :: drhody_y
    real, dimension(isd:ied,jsd:jed,nk,0:1) :: drhodz_z 

    call gradrho(Dens, dTdx, dTdy, dTdz, drhodx_x, drhody_y, drhodz_z)
    call adjust_drhodz(Grid, drhodz_z)
    call diagnose_3d(Time, Grid, id_drhodx, 0.5*sum(drhodx_x, dim=4))
    call diagnose_3d(Time, Grid, id_drhody, 0.5*sum(drhody_y, dim=4))
    call diagnose_3d(Time, Grid, id_drhodz, 0.5*sum(drhodz_z, dim=4), abs_min=1.e-30)
    call diagnose_3d(Time, Grid, id_drhodz_up, drhodz_z(:,:,:,0),     abs_min=1.e-30)
    call diagnose_3d(Time, Grid, id_drhodz_lo, drhodz_z(:,:,:,1),     abs_min=1.e-30)

    call neutral_slopes(Time, Grid, drhodx_x, drhody_y, drhodz_z, &
         slopex_xz, slopey_yz, absslope_z, absslope)

    N2_z(:,:,:,:) = -(grav/rho0)*drhodz_z(:,:,:,:) ! EoM: 16.65

    Ndz_z(COMP,:,0) = sqrt(N2_z(COMP,:,0))*Thickness%dztup(COMP,:)
    Ndz_z(COMP,:,1) = sqrt(N2_z(COMP,:,1))*Thickness%dztlo(COMP,:)

    ! EoM: 16.64, FCOM(2004) 14.83
    gravity_wave_speed(COMP) = sum(sum(Ndz_z(COMP,:,:), dim=4), dim=3)/pi 

    call compute_rossby_radius(gravity_wave_speed, rossby_radius)

    call neutral_blayer(smax, absslope, ksurf_blayer)

    call diagnose_2d_int(Time, Grid, id_ksurf_blayer, ksurf_blayer)
    call diagnose_3d(Time, Grid, id_absslope, absslope, abs_max=1.e4)
    call diagnose_3d(Time, Grid, id_N2, 0.5*sum(N2_z, dim=4), abs_min=1.e-30)
    call diagnose_2d(Time, Grid, id_gravity_wave_speed, gravity_wave_speed)
    call diagnose_2d(Time, Grid, id_rossby_radius, rossby_radius)

  end subroutine density_calculations
  ! </SUBROUTINE> NAME="density_calculations"

  !#######################################################################
  ! <SUBROUTINE NAME="gradrho">
  !
  ! <DESCRIPTION>
  ! Calculate the raw density gradients. No smoothing or limiting  
  ! applied here.
  ! </DESCRIPTION>
  !
  subroutine gradrho(Dens, dTdx, dTdy, dTdz, drhodx_x, drhody_y, drhodz_z)
    type(ocean_density_type),        intent(in) :: Dens
    real, dimension(isd:,jsd:,:,:),  intent(in) :: dTdx
    real, dimension(isd:,jsd:,:,:),  intent(in) :: dTdy
    real, dimension(isd:,jsd:,0:,:), intent(in) :: dTdz

    real, dimension(isd:,jsd:,:,0:), intent(out) :: drhodx_x
    real, dimension(isd:,jsd:,:,0:), intent(out) :: drhody_y
    real, dimension(isd:,jsd:,:,0:), intent(out) :: drhodz_z

    real, dimension(isd:ied,jsd:jed,nk) :: drhodT, drhodS
    real, dimension(isd:ied,jsd:jed,nk) :: dSdx, dSdy, dTdx_, dTdy_
    real, dimension(isd:ied,jsd:jed,0:nk) :: dSdz, dTdz_

    drhodT = Dens%drhodT
    drhodS = Dens%drhodS
    dSdx = dTdx(:,:,:,index_salt)
    dSdy = dTdy(:,:,:,index_salt)
    dSdz = dTdz(:,:,:,index_salt)
    dTdx_ = dTdx(:,:,:,index_temp)
    dTdy_ = dTdy(:,:,:,index_temp)
    dTdz_ = dTdz(:,:,:,index_temp)

    drhodz_z(:,:,:,0) = drhodT(:,:,:)*dTdz_(:,:,0:nk-1) + drhodS(:,:,:)*dSdz(:,:,0:nk-1)
    drhodz_z(:,:,:,1) = drhodT(:,:,:)*dTdz_(:,:,1:nk)   + drhodS(:,:,:)*dSdz(:,:,1:nk)

    drhodx_x(COMP,:,0) = drhodT(COMP,:)*dTdx_(COMPXLL,:) + drhodS(COMP,:)*dSdx(COMPXLL,:) ! G(2004) 16.64
    drhodx_x(COMP,:,1) = drhodT(COMP,:)*dTdx_(COMP,:)    + drhodS(COMP,:)*dSdx(COMP,:)

    drhody_y(COMP,:,0) = drhodT(COMP,:)*dTdy_(COMPYLL,:) + drhodS(COMP,:)*dSdy(COMPYLL,:)
    drhody_y(COMP,:,1) = drhodT(COMP,:)*dTdy_(COMP,:)    + drhodS(COMP,:)*dSdy(COMP,:)

  end subroutine gradrho
  ! </SUBROUTINE> NAME="gradrho"


  !#######################################################################
  ! <SUBROUTINE NAME="adjust_drhodz">
  !
  ! <DESCRIPTION>
  ! Comments about smoothing drhodz:
  !
  ! 1/ Tests in coupled 1-degree model showed extreme sensitivity
  ! of MOC to smoothing.  GFDL users generally do NOT smooth, hence
  ! the default drhodz_smooth_vert=drhodz_smooth_horz=.false.
  !
  ! 2/ Smoothing the vertical derivative of drhodz helps
  ! produce a regularized (i.e., well behaved) neutral slope vector.
  !
  ! 3/ An attempt was made to smooth dTdz and dSdz rather
  ! than drhodz.  The resulting slope was smooth, but not as
  ! smooth as when acting on drhodz itself.
  ! </DESCRIPTION>
  !
  subroutine adjust_drhodz(Grid, drhodz_z)
    type(ocean_grid_type),           intent(in) :: Grid
    real, dimension(isd:,jsd:,:,0:), intent(inout) :: drhodz_z

    integer :: k

    ! vertical neutral density derivative for use in fz_terms
    ! and fz_flux, and for use in fx_flux and fy_flux.
    ! note that the derivative at k=nk vanishes by definition
    ! since these derivatives are at the bottom of tracer cell.
    ! also note the use of -epsln_drhodz ensures the vertical
    ! derivative is always < 0.
    drhodz_z(:,:,:,:) = min(drhodz_z(:,:,:,:), -epsln)

    ! vertically smooth the vertical derivative of density to
    ! produce a smooth neutral slope vector for all flux components.
    if (drhodz_smooth_vert) then
       call vert_smooth(Grid%kmt(:,:), drhodz_z(:,:,:,0))
       call vert_smooth(Grid%kmt(:,:), drhodz_z(:,:,:,1))
    endif

    ! horizontally smooth (via a diffusivity) the vertical derivative of 
    ! density to produce a smooth neutral slope vector for all flux components.
    if (drhodz_smooth_horz) then
       do k = 1, nk
          drhodz_z(:,:,k,0) = drhodz_z(:,:,k,0) + dtime*LAP_T(drhodz_z(:,:,k,0), smooth_diff(:,:))
          drhodz_z(:,:,k,1) = drhodz_z(:,:,k,1) + dtime*LAP_T(drhodz_z(:,:,k,1), smooth_diff(:,:))
       enddo
    endif

    drhodz_z(:,:,:,0) = Grid%tmask(:,:,:)*drhodz_z(:,:,:,0)
    do k = 1, nk-1
       drhodz_z(:,:,k,1) = Grid%tmask(:,:,k+1)*drhodz_z(:,:,k,1)
    enddo

  end subroutine adjust_drhodz
  ! </SUBROUTINE> NAME="adjust_drhodz"


  !#######################################################################
  ! <SUBROUTINE NAME="neutral_slopes">
  !
  ! <DESCRIPTION>
  ! Compute the neutral slope vector along with its magnitude.
  ! The neutral slope vector is defined as -grad_h(rho)/(drho/dz).
  ! </DESCRIPTION>
  !
  subroutine neutral_slopes(Time, Grid, drhodx_x, drhody_y, drhodz_z, &
       slopex_xz, slopey_yz, absslope_z, absslope)

    type(ocean_time_type),           intent(in) :: Time
    type(ocean_grid_type),           intent(in) :: Grid
    real, dimension(isd:,jsd:,:,0:), intent(in) :: drhodx_x
    real, dimension(isd:,jsd:,:,0:), intent(in) :: drhody_y
    real, dimension(isd:,jsd:,:,0:), intent(in) :: drhodz_z

    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: slopex_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: slopey_yz
    real, dimension(isd:,jsd:,:,0:),    intent(out) :: absslope_z
    real, dimension(isd:,jsd:,:),       intent(out) :: absslope

    real, dimension(isd:ied,jsd:jed,nk)     :: slopex
    real, dimension(isd:ied,jsd:jed,nk)     :: slopey
    real, dimension(isd:ied,jsd:jed,nk,0:1) :: slopex_z
    real, dimension(isd:ied,jsd:jed,nk,0:1) :: slopey_z
    integer :: k

    slopex_xz(COMP,:,:,:) = 0.0
    slopey_yz(COMP,:,:,:) = 0.0

    ! Upper half cell slopes
    do k = 2, nk       
       where (Grid%tmask(COMP,k) == 1.0)
          slopex_xz(COMP,k,0,0) = -drhodx_x(COMP,k,0)/drhodz_z(COMP,k,0)
          slopex_xz(COMP,k,1,0) = -drhodx_x(COMP,k,1)/drhodz_z(COMP,k,0)
       
          slopey_yz(COMP,k,0,0) = -drhody_y(COMP,k,0)/drhodz_z(COMP,k,0)
          slopey_yz(COMP,k,1,0) = -drhody_y(COMP,k,1)/drhodz_z(COMP,k,0)
       endwhere
    enddo

    ! Lower half cell slopes
    do k = 1, nk-1
       where (Grid%tmask(COMP,k+1) == 1.0)
          slopex_xz(COMP,k,0,1) = -drhodx_x(COMP,k,0)/drhodz_z(COMP,k,1)
          slopex_xz(COMP,k,1,1) = -drhodx_x(COMP,k,1)/drhodz_z(COMP,k,1)

          slopey_yz(COMP,k,0,1) = -drhody_y(COMP,k,0)/drhodz_z(COMP,k,1)
          slopey_yz(COMP,k,1,1) = -drhody_y(COMP,k,1)/drhodz_z(COMP,k,1)
       endwhere
    enddo

    ! Set the very top level slope
    slopex_xz(COMP,1,:,0) = slopex_xz(COMP,1,:,1)
    slopey_yz(COMP,1,:,0) = slopey_yz(COMP,1,:,1)

    ! Set the very bottom level slope
    do k = 1, nk
       where (k == Grid%kmt(COMP))
          slopex_xz(COMP,k,0,1) = slopex_xz(COMP,k,0,0)
          slopex_xz(COMP,k,1,1) = slopex_xz(COMP,k,1,0)
          slopey_yz(COMP,k,0,1) = slopey_yz(COMP,k,0,0)
          slopey_yz(COMP,k,1,1) = slopey_yz(COMP,k,1,0)
       endwhere
    enddo

    slopex_z(COMP,:,:)   = sum(slopex_xz(COMP,:,:,:), dim=4)/2.0
    slopey_z(COMP,:,:)   = sum(slopey_yz(COMP,:,:,:), dim=4)/2.0
    absslope_z(COMP,:,:) = sqrt(slopex_z(COMP,:,:)**2 + slopey_z(COMP,:,:)**2)

    slopex(COMP,:)   = sum(slopex_z(COMP,:,:), dim=4)/2.0
    slopey(COMP,:)   = sum(slopey_z(COMP,:,:), dim=4)/2.0
    absslope(COMP,:) = sqrt(slopex(COMP,:)**2 + slopey(COMP,:)**2)

    call diagnose_3d(Time, Grid, id_slopex, slopex, abs_max=1.e10)
    call diagnose_3d(Time, Grid, id_slopey, slopey, abs_max=1.e10)

  end subroutine neutral_slopes
  ! </SUBROUTINE> NAME="neutral_slopes"

  !#######################################################################
  ! <SUBROUTINE NAME="neutral_blayer">
  !
  ! <DESCRIPTION>
  ! Locate the vertical index of the neutral boundary layer. This layer is
  ! defined as the point where the magnitude of the neutral slope vector
  ! first drops below smax, when searching down from the surface.
  ! </DESCRIPTION>
  !
  subroutine neutral_blayer(smax, absslope, ksurf_blayer)
    real,                          intent(in) :: smax
    real, dimension(isd:,jsd:,:),  intent(in) :: absslope
    integer, dimension(isd:,jsd:), intent(out) :: ksurf_blayer

    logical, dimension(isc:iec,jsc:jec) :: slopeok
    integer :: k

    slopeok(COMP) = .false.
    ksurf_blayer(COMP) = nk
    do k = nk, 1, -1
       where (absslope(COMP,k) < smax .and. .not. slopeok(COMP))
          ksurf_blayer(COMP) = k - 1
       elsewhere
          slopeok(COMP) = .true.
       endwhere
       if (all(slopeok(COMP))) exit
    enddo

  end subroutine neutral_blayer
  ! </SUBROUTINE> NAME="neutral_blayer"

  !#######################################################################
  ! <SUBROUTINE NAME="ocean_nphysics_new_restart">
  !
  ! <DESCRIPTION>
  ! Write out the restart data for this module
  ! </DESCRIPTION>
  !
  subroutine ocean_nphysics_new_restart(time_stamp)
    character(len=*), intent(in), optional :: time_stamp

    if(.not. use_this_module) return

    call save_restart(nphysics_new_restart, time_stamp)

  end subroutine ocean_nphysics_new_restart
  ! </SUBROUTINE> NAME="ocean_nphysics_new_restart"

  !#######################################################################
  ! <SUBROUTINE NAME="ocean_nphysics_new_end">
  !
  ! <DESCRIPTION>
  ! Writes out the restart data.
  ! </DESCRIPTION>
  !
  subroutine ocean_nphysics_new_end()

    if(.not. use_this_module) return

    if (.not. module_is_initialized ) then 
       call mpp_error(FATAL, &
            '==>Error from ocean_nphysics_new (ocean_nphysics_new_end): needs initialization')
    endif

    call ocean_nphysics_new_restart()

  end subroutine ocean_nphysics_new_end
  ! </SUBROUTINE> NAME="ocean_nphysics_new_end"

end module ocean_nphysics_new_mod

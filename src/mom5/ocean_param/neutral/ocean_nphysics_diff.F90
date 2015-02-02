module ocean_nphysics_diff_mod
! 
!<CONTACT EMAIL="tim.leslie@gmail.com"> Tim Leslie
!</CONTACT>
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Compute the neutral diffusivity in this module. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the neutral diffusivity. There are many 
! methods available. 
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
!<NAMELIST NAME="ocean_nphysics_diff_nml">
!
!  <DATA NAME="" TYPE="logical">
!  </DATA> 
!
!</NAMELIST>

#define COMP isc:iec,jsc:jec
#define HARM_MEAN(a,b) (2.0*(a)*(b)/((a)+(b)+epsln))

  use constants_mod,               only: epsln
  use fms_mod,                     only: open_namelist_file, check_nml_error, close_file, FATAL
  use mpp_domains_mod,             only: mpp_update_domains
  use mpp_mod,                     only: mpp_pe, mpp_min, stdout, stdlog, mpp_error, input_nml_file
  
  use ocean_domains_mod,           only: get_local_indices, set_ocean_domain
  use ocean_nphysics_util_new_mod, only: vert_smooth, horz_smooth, compute_rossby_radius
  use ocean_operators_mod,         only: S2D
  use ocean_parameters_mod,        only: omega_earth
  use ocean_tracer_diag_mod,       only: calc_mixed_layer_depth
  use ocean_types_mod,             only: ocean_grid_type, ocean_domain_type, ocean_time_type, ocean_thickness_type
  use ocean_types_mod,             only: ocean_prog_tracer_type, ocean_density_type
  use ocean_util_mod,              only: write_note
  use ocean_util_mod,              only: diagnose_2d, diagnose_3d, register_3d_t_field, register_2d_t_field

  implicit none

  public ocean_nphysics_diff_init
  public compute_diffusivity

  private check_stability
  private check_nml_options
  private compute_agm
  private compute_raw_growth_rate
  private compute_growth_rate
  private vertical_average
  private compute_length
  private compute_bczone_radius
  private compute_aredi
  private apply_grid_scaling

  private

#include <ocean_memory.h>

character(len=128) :: version=&
     '$Id: ocean_nphysics_diff.F90,v 20.0 2013/12/14 00:14:42 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

character(len=*), parameter :: FILENAME=&
     __FILE__


  integer :: id_agm_array
  integer :: id_aredi_array

  integer :: id_N2_ref
  integer :: id_agm_length
  integer :: id_raw_growth_rate
  integer :: id_growth_rate
  integer :: id_growth_rate_zave
  integer :: id_agm_fast

  integer :: id_ml_depth
  integer :: id_ave_ml_rate
  integer :: id_eg_rate

  integer :: id_bczone_radius
  integer :: id_rhines_length


  ! 2D array of micom gm diffusivities (m^2/sec)
  real, dimension(:,:), allocatable :: agm_micom

  ! 2D array of gm growth rate max (sec^-1)
  real, dimension(:,:), allocatable :: growth_rate_max

  ! beta = d(Coriolis)/dy (m^-1 sec^-1)
  real, dimension(:,:), allocatable :: beta_param

  ! for determining baroclinic zone radius (metre)
  real,    dimension(:,:), allocatable :: bczone_dxt
  real,    dimension(:,:), allocatable :: bczone_dyt
  integer, dimension(:,:), allocatable :: bczone_over

  type(ocean_domain_type), save :: BCzone_domain
  real :: agm_gamma_damp


  ! -- GM Diffusivity parameters set via namelist --
  
  ! globally constant diffusivities
  logical :: agm_const = .true.
  real    :: agm = 1.0e3 ! constant skew-diffusion diffusivity (m^2/sec)
  
  ! specify agm according to local horizontal area
  ! of grid and a specified velocity scale (m/s)
  logical :: agm_closure_micom = .false.
  real    :: agm_micom_vel     = 0.0 ! m/s

  ! specify agm according to latitude zones
  logical :: agm_lat_bands          = .false.
  real    :: agm_lat_bands_boundary = -999.   ! boundary between agm in the south and north zones in degrees
  real    :: agm_lat_bands_ratio    = 1.0     ! ratio agm(south)/agm(north)

  ! for computing time dependent agm according to flow properties
  logical :: agm_closure = .false.
  real    :: agm_max     = 2.e3   ! maximum diffusivity allowed when agm_closure=.true.
  real    :: agm_min     = 2.e2   ! minimum diffusivity allowed when agm_closure=.true.
  real    :: agm_scaling = 2.0    ! dimensionless tuning parameter for flow dependent agm

  ! agm computed as growth_rate*length_scale**2
  logical :: agm_closure_rate_len2 = .false.

  ! scaling flow dependent agm according to Rossby radius and grid scale 
  logical :: agm_grid_scaling       = .false.
  real    :: agm_grid_scaling_power = 2.0 ! no units

  ! for smoothing agm_array in space and time
  logical :: agm_smooth_space = .false.
  logical :: agm_smooth_time  = .false.
  real    :: agm_damping_time = 10.0 ! days

  ! for scaling agm by buoyancy ratio (N/Nref)^2
  logical :: agm_closure_n2_scale   = .false.
  real    :: agm_n2_scale_coeff     = 1e3   ! m^2/s
  logical :: agm_n2_scale_nref_cst  = .false.
  real    :: agm_n2_scale_buoy_freq = 0.004 ! 1/sec

  ! depths over which compute the Eady growth and/or horizontal density gradient
  ! defaults of 100-2000 taken from Treguier, Held, and Larichev (1997)
  real    :: agm_rate_upper_depth = 100.0  ! metre
  real    :: agm_rate_lower_depth = 2000.0 ! metre

  logical :: agm_rate_eady       = .false.
  logical :: agm_rate_baroclinic = .false.

  ! sec^-1 buoyancy frequency for use with agm_rate_baroclinic
  real    :: agm_rate_baro_buoy_freq  = 0.004 ! 1/sec

  ! for smoothing and capping the growth rate
  logical :: agm_rate_smooth_vert  = .false.
  logical :: agm_rate_smooth_horiz = .false.
  logical :: agm_rate_cap          = .false.

  ! dimensionless number to set maximum value of agm_growth.
  ! agm_rate_cap_scale = 0.5 yields a maximum agm_growth rate of coriolis parameter.
  real    :: agm_rate_cap_scale = 0.5 ! no units

  logical :: agm_rate_ave_mixed = .false.
  logical :: agm_rate_zave      = .false.

  ! compute agm according to Eden and Greatbatch (2008) and Eden(2007)
  logical :: agm_rate_eden_greatbatch   = .false.
  real    :: agm_rate_eg_alpha          = 0.07 ! no units
  logical :: agm_length_eden_greatbatch = .false.

  ! for capping the agm_length scale
  logical :: agm_length_cap = .false.
  real    :: agm_length_max = 50e3 ! metres

  ! for fixed length scale (metres) set by agm_length
  logical :: agm_length_fixed = .false.
  real    :: agm_length       = 50.e3  ! metres

  ! for length scale set according to estimate of first baroclinic Rossby radius
  logical :: agm_length_rossby = .false.

  ! for length scale set according to radius of baroclinic zone
  logical :: agm_length_bczone    = .false.
  integer :: agm_bczone_max_pts   = 10 ! max # points searched for determining baroclinic zone width
  real    :: agm_bczone_crit_rate = 1.4e-6  ! critical growth rate for determining baroclinic zone (sec^-1)


  ! -- Diffusivity for Redi diffusion tensor -- 
  ! will set aredi_array=agm_array
  logical :: aredi_equal_agm =.true.

  ! constant neutral diffusion tracer diffusivity (m^2/sec)
  real    :: aredi = 1.0e3  ! m^2/s
  logical :: aredi_fixed = .false.

  ! for scaling aredi_array according to size of grid and Rossby radius
  logical :: aredi_grid_scaling       = .false.
  real    :: aredi_grid_scaling_power = 2.0 ! no units



  ! -- Horizontal boundary layer diffusivity --
  ! for applying lateral diffusivity within neutral boundary layer
  logical :: neutral_horiz_mix_bdy =.false.
  logical :: ah_micom     = .false.
  real    :: ah_micom_vel = 0.0  ! velocity scale (m/s) for horizontal diffusivity w/i boundary
  logical :: ah_const     = .false.
  real    :: ah           = 0.0  ! constant horiz diffusivity in neutral bdy



  namelist /ocean_nphysics_diff_nml/                                                            &
       neutral_horiz_mix_bdy, ah_micom, ah_micom_vel, ah_const, ah,                             &       
       aredi_equal_agm, aredi_fixed, aredi, aredi_grid_scaling, aredi_grid_scaling_power,       &
       agm_const, agm, agm_lat_bands, agm_lat_bands_boundary, agm_lat_bands_ratio,              &
       agm_closure, agm_max, agm_min, agm_grid_scaling, agm_grid_scaling_power,                 &
       agm_smooth_space, agm_smooth_time, agm_damping_time,                                     &
       agm_closure_micom, agm_micom_vel,                                                        &
       agm_closure_n2_scale, agm_n2_scale_coeff, agm_n2_scale_nref_cst, agm_n2_scale_buoy_freq, &
       agm_closure_rate_len2, agm_scaling, agm_rate_upper_depth, agm_rate_lower_depth,          &
       agm_rate_eady, agm_rate_baroclinic, agm_rate_baro_buoy_freq,                             &
       agm_rate_smooth_vert, agm_rate_smooth_horiz, agm_rate_cap, agm_rate_cap_scale,           &
       agm_rate_ave_mixed, agm_rate_zave, agm_rate_eden_greatbatch, agm_rate_eg_alpha,          &
       agm_length_cap, agm_length_max, agm_length_fixed, agm_length,                            &
       agm_length_rossby, agm_length_bczone, agm_bczone_max_pts, agm_bczone_crit_rate,          &
       agm_length_eden_greatbatch

contains

  subroutine ocean_nphysics_diff_init(Domain, Grid, Time, Thickness, grid_length, dtime, smax, &
       ah_array, agm_array, aredi_array, agm_z_dep)

    type(ocean_domain_type),    intent(in) :: Domain
    type(ocean_grid_type),      intent(in) :: Grid
    type(ocean_time_type),      intent(in) :: Time
    type(ocean_thickness_type), intent(in) :: Thickness
    real, dimension(isd:,jsd:), intent(in) :: grid_length
    real,                       intent(in) :: dtime
    real,                       intent(in) :: smax

    real, dimension(isd:,jsd:),   intent(out) :: ah_array
    real, dimension(isd:,jsd:,:), intent(out) :: agm_array
    real, dimension(isd:,jsd:,:), intent(out) :: aredi_array
    logical,                      intent(out) :: agm_z_dep

    call check_nml_options(agm_z_dep)

#ifndef MOM_STATIC_ARRAYS
    call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
    nk = Grid%nk
#endif

    call init_globals(Domain, Grid, grid_length, dtime)

    call register_fields(Grid, Time)

    call diffusivity_init(Grid, grid_length, ah_array, agm_array, aredi_array)
    
    call check_stability(Grid, Thickness, aredi_array, smax, dtime)

  end subroutine ocean_nphysics_diff_init

  !#######################################################################
  ! <SUBROUTINE NAME="check_nml_options">
  !
  ! <DESCRIPTION>
  ! Read in the namelist parameters and ensure that valid values have been
  ! choosen.
  ! Also determine whether agm is z-dependent.
  ! </DESCRIPTION>
  !
  subroutine check_nml_options(agm_z_dep)

    logical, intent(out) :: agm_z_dep

    integer :: ierr, ioun, io_status
    integer :: stdoutunit, stdlogunit
    stdoutunit = stdout()
    stdlogunit = stdlog()

    ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_nphysics_diff_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_nphysics_diff_nml')
#else
    ioun =  open_namelist_file()
    read (ioun,ocean_nphysics_diff_nml,IOSTAT=io_status)
    ierr = check_nml_error(io_status,'ocean_nphysics_diff_nml')
    call close_file (ioun)
#endif
    write (stdoutunit,'(/)')
    write (stdoutunit,ocean_nphysics_diff_nml)  
    write (stdlogunit,ocean_nphysics_diff_nml)

    ! Check for horizontal boundary layer mixing
    if (neutral_horiz_mix_bdy) then
       if (count((/ah_micom, ah_const /)) /= 1) then
          call mpp_error(FATAL, &
               '==>Error: with neutral_horiz_mix_bdy=.true., must choose exactly 1 horizontal diffusivity closure')
       endif
       call write_note(FILENAME,&
       'neutral_horiz_mix_bdy=.true. => Adding horizontal diffusion in the boundary layer')
    endif

    ! Check neutral diffusion method
    if (count((/aredi_fixed, aredi_equal_agm /)) /= 1) then
          call mpp_error(FATAL, &
               '==>Error: must choose exactly 1 neutral diffusivity closure')
    endif
    if (aredi_fixed) then
       call write_note(FILENAME,&
       'Using fixed constant neutral diffusion diffusivity')
    else
       call write_note(FILENAME,&
       'Using same diffusion for neutral diffusion and GM')
    endif

    if (count((/agm_closure, agm_const /)) /= 1) then
          call mpp_error(FATAL, &
               '==>Error: must choose to use either constant GM diffusivity or a closure')
    endif

    ! Check GM method
    agm_z_dep = .false.
    if (agm_closure) then
       if (count((/agm_closure_n2_scale, agm_closure_micom, agm_closure_rate_len2 /)) /= 1) then
          call mpp_error(FATAL, &
               '==>Error: with agm_closure=.true., must choose exactly 1 GM diffusivity closure')
       endif
       
       if (agm_closure_rate_len2) then
          if (count((/agm_rate_eady, agm_rate_baroclinic /)) /= 1) then
             call mpp_error(FATAL, &
                  '==>Error: with agm_closure_rate_len2=.true., must choose exactly 1 growth rate')
          endif

          if (count((/agm_length_fixed, agm_length_rossby, agm_length_bczone, agm_length_eden_greatbatch /)) /= 1) then
             call mpp_error(FATAL, &
                  '==>Error: with agm_closure_rate_len2=.true., must choose exactly 1 length scale')
          endif
          call write_note(FILENAME,&
          'GM diffusivity = rate*length^2')
       else if (agm_closure_n2_scale) then
          call write_note(FILENAME,&
          'GM diffusivity = const*N^2')
       else if (agm_closure_micom) then
          call write_note(FILENAME,&
          'GM diffusivity = grid length * const velocity')
       endif

       ! Check for z-dependence.
       if (agm_closure_micom) then
       else if (agm_closure_n2_scale) then
          agm_z_dep = .true.
       else if (agm_closure_rate_len2) then
          if (.not. agm_rate_zave) agm_z_dep = .true.
          if (agm_length_eden_greatbatch) agm_z_dep = .true.
       endif
    else
       call write_note(FILENAME,&
       'Using constant GM diffusivity')
    endif

  end subroutine check_nml_options
  ! </SUBROUTINE> NAME="check_nml_options"

  !#######################################################################
  ! <SUBROUTINE NAME="init_globals">
  !
  ! <DESCRIPTION>
  ! Allocate and initialise all (non-namelist) global variables in the module.
  ! </DESCRIPTION>
  !
  subroutine init_globals(Domain, Grid, grid_length, dtime)

    type(ocean_domain_type),    intent(in) :: Domain
    type(ocean_grid_type),      intent(in) :: Grid
    real, dimension(isd:,jsd:), intent(in) :: grid_length
    real,                       intent(in) :: dtime

    real, dimension(isd:ied,jsd:jed) :: coriolis_param

    integer :: max_pts

    ! Coriolis parameter and beta parameter
    allocate (beta_param(isd:ied,jsd:jed))
    beta_param(COMP) = max(2.28e-11*abs(cos(Grid%phit(COMP))), epsln)

    if (agm_closure) then
       if (agm_smooth_time) then
          agm_gamma_damp = dtime/(24*60*60*agm_damping_time) !for damping time dependent agm_array
       endif
       
       if (agm_closure_micom) then
          allocate (agm_micom(isd:ied,jsd:jed))
          agm_micom(:,:) = agm_micom_vel*grid_length(:,:)
       endif

       if (agm_closure_rate_len2) then
          if (agm_rate_cap) then
             allocate (growth_rate_max(isd:ied,jsd:jed))
             coriolis_param(COMP) = 2.0*omega_earth*abs(sin(Grid%phit(COMP)))
             growth_rate_max(COMP) = agm_rate_cap_scale*coriolis_param(COMP)
          endif
         
          ! for computing the baroclinic zone radius using Hadley Centre search algorithm
          if (agm_length_bczone) then
             max_pts = agm_bczone_max_pts
             call set_ocean_domain(BCzone_domain,Grid,xhalo=max_pts,yhalo=max_pts,name='bczone',maskmap=Domain%maskmap)
             allocate (bczone_over(isc-max_pts:iec+max_pts,jsc-max_pts:jec+max_pts))
             allocate (bczone_dxt (isc-max_pts:iec+max_pts,jsc-max_pts:jec+max_pts))
             allocate (bczone_dyt (isc-max_pts:iec+max_pts,jsc-max_pts:jec+max_pts))
             bczone_dxt(COMP) = Grid%dxt(COMP)
             bczone_dyt(COMP) = Grid%dyt(COMP)
             call mpp_update_domains(bczone_dxt(:,:), BCzone_domain%domain2d)
             call mpp_update_domains(bczone_dyt(:,:), BCzone_domain%domain2d)
          endif
       endif
    endif

  end subroutine init_globals
  ! </SUBROUTINE> NAME="init_globals"

  !#######################################################################
  ! <SUBROUTINE NAME="register fields">
  !
  ! <DESCRIPTION>
  ! Register diagnostic fields.
  ! </DESCRIPTION>
  !
  subroutine register_fields(Grid, Time)
    type(ocean_grid_type), intent(in) :: Grid
    type(ocean_time_type), intent(in) :: Time

    id_agm_array = register_3d_t_field(Grid, Time, 'agm_array', &
         'Skew diffusivity', 'm^2/s', (/-1e-8,1.e8/))
    id_aredi_array = register_3d_t_field(Grid, Time, 'aredi_array', &
         'Neutral diffusivity', 'm^2/s', (/-1e-8,1.e8/))

    id_N2_ref = register_2d_t_field(Grid, Time, 'N2_ref', &
         'Reference squared buoyancy frequency', '1/s^2', (/-1e-8,1.e8/))
    id_agm_length = register_3d_t_field(Grid, Time, 'agm_length', &
         'Characteristic length in diffusivity closure calculation', 'm', (/-1e-8,1.e8/))

    id_raw_growth_rate = register_3d_t_field(Grid, Time, 'raw_growth_rate', &
         'Raw growth rate from either baroclinicity or eady rate', '1/s', (/-1e-8,1.e8/))
    id_growth_rate = register_3d_t_field(Grid, Time, 'growth_rate', &
         'Characteristic growth rate in diffusivity closure calculation', '1/s', (/-1e-8,1.e8/))
    id_growth_rate = register_2d_t_field(Grid, Time, 'growth_rate_zave', &
         'Verical average of growth rate used in diffusivity closure', '1/s', (/-1e-8,1.e8/))

    id_agm_fast = register_3d_t_field(Grid, Time, 'agm_fast', &
         'Skew diffusivity closure value before temporal smoothing is applied', 'm^2/s', (/-1e-8,1.e8/))

    id_ml_depth = register_2d_t_field(Grid, Time, 'ml_depth', &
         'Mixed layer depth used to compute constant growth rate region', 'm', (/0.0, 8000.0/))
    id_ave_ml_rate = register_2d_t_field(Grid, Time, 'ave_ml_rate', &
         'Average growth rate within the mixed layer', '1/s', (/-1e-8,1.e8/))
    id_eg_rate = register_2d_t_field(Grid, Time, 'eg_rate', &
         'Growth rate as calculated by the method of Eden & Greatbatch', '1/s', (/-1e-8,1.e8/))

    id_bczone_radius = register_2d_t_field(Grid, Time, 'bczone_radius', &
         'Radius of the baroclinic zone', 'm', (/-1e-8,1.e8/))
    id_rhines_length = register_2d_t_field(Grid, Time, 'rhines_length', &
         'The Rhines length as used in the Eden & Greatbatch closure', 'm', (/-1e-8,1.e8/))

  end subroutine register_fields
  ! </SUBROUTINE> NAME="register_fields"

  !#######################################################################
  ! <SUBROUTINE NAME="diffusivity_init">
  !
  ! <DESCRIPTION>
  ! Initialise the three diffusivity arrays.
  ! </DESCRIPTION>
  !
  subroutine diffusivity_init(Grid, grid_length, &
       ah_array, agm_array, aredi_array)

    type(ocean_grid_type),      intent(in) :: Grid
    real, dimension(isd:,jsd:), intent(in) :: grid_length

    real, dimension(isd:,jsd:),   intent(out) :: ah_array
    real, dimension(isd:,jsd:,:), intent(out) :: agm_array
    real, dimension(isd:,jsd:,:), intent(out) :: aredi_array

    real, dimension(isd:ied,jsd:jed) :: gravity_wave_speed
    real, dimension(isd:ied,jsd:jed) :: rossby_radius

    gravity_wave_speed(COMP) = 0.0
    call compute_rossby_radius(gravity_wave_speed, rossby_radius)

    ! horizontal diffusivity in neutral boundary layer region
    if (neutral_horiz_mix_bdy) then
       if (ah_micom) then
          ah_array(COMP) = ah_micom_vel*grid_length(COMP)*Grid%tmask(COMP,1)
       elseif (ah_const) then
          ah_array(:,:) = ah*Grid%tmask(:,:,1)
       endif
    else
       ah_array(:,:) = 0.0
    endif

    if (agm_const) then
       ! initialization based on constant coefficients
       agm_array(:,:,:) = agm

       ! Set diffusivity according to agm_lat_bands
       if (agm_lat_bands) then
          where (spread(Grid%yt(:,:), 3, nk) <= agm_lat_bands_boundary)
             agm_array(:,:,:) = agm*agm_lat_bands_ratio
          endwhere
       endif
    elseif (agm_closure) then
       ! grid-scale dependent diffusivity suggested by that commonly used in MICOM
       ! vel_micom (m/s) sets velocity scale.
       ! space scale is set by grid size.
       if (agm_closure_micom) then
          agm_array(:,:,:) = spread(agm_micom(:,:), 3, nk)
       else
          ! initialization based on constant coefficients
          agm_array(:,:,:) = agm
       endif
    endif

    ! now apply tracer mask
    agm_array(:,:,:) = agm_array(:,:,:)*Grid%tmask(:,:,:)

    ! Initialise Redi diffusivity
    if (aredi_fixed) then
       aredi_array(:,:,:) = aredi*Grid%tmask(:,:,:)
    endif

    call compute_aredi(Grid%tmask, grid_length, rossby_radius, agm_array, &
         aredi_array)

  end subroutine diffusivity_init
  ! </SUBROUTINE> NAME="diffusivity_init"

  !#######################################################################
  ! <SUBROUTINE NAME="check_stability">
  !
  ! <DESCRIPTION>
  ! Check the stability assumptions and print details of the limits of stability.
  ! </DESCRIPTION>
  !
  subroutine check_stability(Grid, Thickness, aredi_array, smax, dtime)

    type(ocean_grid_type),        intent(in) :: Grid
    type(ocean_thickness_type),   intent(in) :: Thickness
    real, dimension(isd:,jsd:,:), intent(in) :: aredi_array
    real,                         intent(in) :: smax
    real,                         intent(in) :: dtime

    real :: delta_iso, delta_iso0, delta_min, A_max, A_max0
    real, dimension(isd:ied,jsd:jed,nk) :: min_delta, min_xy_dz
    real, dimension(isd:ied,jsd:jed) :: min_xy
    integer :: i_delta, j_delta, k_delta, k
    integer, dimension(3) :: index_delta
    integer :: unit = 6 ! FIXME: This is a hack to have a non-root PE write to stdout
                        ! but will fail on SGICRAY where stdout != 6.

    ! FOCM - 15.1 Eqn 15.4

    ! compute maximum stable neutral slope available for neutral diffusion
    min_xy(COMP) = min(Grid%dxt(COMP), Grid%dyt(COMP))
    do k = 1, nk
       min_xy_dz(COMP,k) = min_xy(COMP)*Thickness%dzt(COMP,k)
    enddo
    min_delta(COMP,:) = min_xy_dz(COMP,:)/(2*aredi_array(COMP,:)*dtime + epsln)

    delta_iso = minval(min_delta(COMP,:), mask=Grid%tmask(COMP,:)==1.0)
    delta_iso = delta_iso + 1.e-6*mpp_pe() ! to separate redundancies
    delta_iso0 = delta_iso
    call mpp_min (delta_iso)

    ! show most unstable location
    if (delta_iso == delta_iso0) then
       index_delta = minloc(min_delta(COMP,:), mask=Grid%tmask(COMP,:)==1.0)
       i_delta = index_delta(1) - 1 + isc
       j_delta = index_delta(2) - 1 + jsc
       k_delta = index_delta(3)
       write(unit,'(a)')
       write(unit,'(a)')'---Neutral direction slope check I for linear stability of neutral diffusion---'
       write(unit,'(a,e14.7)')'With a neutral physics time step (secs) of ', dtime
       write(unit,'(a)')'the most stringent linear stability constraint was found at the following ocean cell:'
       write(unit,'(a,i4,a,i4,a,e14.7)')'long(',i_delta,',',j_delta,')                   = ',Grid%xt(i_delta,j_delta)
       write(unit,'(a,i4,a,i4,a,e14.7)')'lat (',i_delta,',',j_delta,')                   = ',Grid%yt(i_delta,j_delta)
       write(unit,'(a,i4,a,i4,a,i4,a,e14.7)')'thick(',i_delta,',',j_delta,',',k_delta,') = ', &
            Thickness%dzt(i_delta,j_delta,k_delta)
       write(unit,'(a,i4,a,i4,a,i4,a,e14.7)')'aredi(',i_delta,',',j_delta,',',k_delta,') = ', &
            aredi_array(i_delta,j_delta,k_delta)
       write(unit,'(a,e14.7,a)')'delta_iso           = ',delta_iso,' is the maximum neutral direction slope'
       write(unit,'(a)')'available for linear stability of the neutral diffusion scheme.'
       write(unit,'(a)')'The namelist parameter smax should conservatively be <= delta_iso.'
       if (smax >= delta_iso) then
          write(unit,'(a,f10.5,a)')'==> Warning: The namelist parameter smax= ',smax, ' is >= to delta_iso.'
          write(unit,'(a)')'Linear stability of the neutral diffusion scheme may be compromised.'
       endif
       write(unit,'(a)')
    endif

    ! Compute maximum diffusivity available given a maximum slope of smax
    delta_min = minval(min_xy_dz(COMP,:), mask=Grid%tmask(COMP,:)==1.0)
    A_max = delta_min/(2*smax*dtime + epsln) + 1.e-6*mpp_pe() ! to separate redundancies
    A_max0 = A_max
    call mpp_min (A_max)    

    ! show most unstable location
    if (A_max == A_max0) then
       index_delta = minloc(min_xy_dz(COMP,:), mask=Grid%tmask(COMP,:)==1.0)
       i_delta = index_delta(1) - 1 + isc
       j_delta = index_delta(2) - 1 + jsc
       k_delta = index_delta(3)

       write(unit,'(a)')
       write(unit,'(a)')'---Neutral direction slope check II for linear stability of neutral diffusion---'
       write(unit,'(a,e14.7)')'Assuming maximum Redi neutral diffusion slope of ', smax
       write(unit,'(a,e14.7)')'and neutral physics time step (secs) of ', dtime
       write(unit,'(a)')'the most stringent linear stability constraint was found at the following ocean cell:'
       write(unit,'(a,i4,a,i4,a,e14.7)')'long(',i_delta,',',j_delta,')                  = ',&
            Grid%xt(i_delta,j_delta)
       write(unit,'(a,i4,a,i4,a,e14.7)')'lat (',i_delta,',',j_delta,')                  = ',&
            Grid%yt(i_delta,j_delta)
       write(unit,'(a,i4,a,i4,a,i4,a,e14.7)')'thick(',i_delta,',',j_delta,',',k_delta,')= ',&
            Thickness%dzt(i_delta,j_delta,k_delta)
       write(unit,'(a,e14.7,a)')'A_max      = ',A_max,' (m^2/sec) is the maximum neutral diffusivity'
       write(unit,'(a)')'available for linear stability of the neutral diffusion scheme.'
       write(unit,'(a)')'Conservatively, neutral diffusivities used in the model should be less than A_max.'
       write(unit,'(a)')'--------------------------------------------------------------------------------'
       write(unit,'(a)')
    endif

  end subroutine check_stability
  ! </SUBROUTINE> NAME="check_stability"


  !#######################################################################
  ! <SUBROUTINE NAME="compute_diffusivity">
  !
  ! <DESCRIPTION>
  ! </DESCRIPTION>
  !
  subroutine compute_diffusivity(Domain, Grid, Time, Thickness, Dens, T_temp, T_salt, &
       absslope_z, N2_z, gravity_wave_speed, rossby_radius, &
       ksurf_blayer, grid_length, &
       agm_array, aredi_array)

    type(ocean_domain_type),      intent(inout) :: Domain
    type(ocean_grid_type),        intent(in)    :: Grid
    type(ocean_time_type),        intent(in)    :: Time
    type(ocean_thickness_type),   intent(in)    :: Thickness
    type(ocean_density_type),     intent(in)    :: Dens
    type(ocean_prog_tracer_type), intent(in)    :: T_temp
    type(ocean_prog_tracer_type), intent(in)    :: T_salt

    real, dimension(isd:,jsd:,:,0:), intent(in) :: absslope_z
    real, dimension(isd:,jsd:,:,0:), intent(in) :: N2_z
    real, dimension(isd:,jsd:),      intent(in) :: gravity_wave_speed
    real, dimension(isd:,jsd:),      intent(in) :: rossby_radius

    integer, dimension(isd:,jsd:), intent(in) :: ksurf_blayer
    real, dimension(isd:,jsd:),    intent(in) :: grid_length

    real, dimension(isd:,jsd:,:), intent(inout) :: agm_array
    real, dimension(isd:,jsd:,:), intent(inout) :: aredi_array

    ! update closure-based diffusivity for next time step
    if (agm_closure) then
       call compute_agm(Domain, Grid, Time, Thickness, Dens, T_temp, T_salt, &
            absslope_z, N2_z, gravity_wave_speed, rossby_radius, &
            ksurf_blayer, grid_length, &
            agm_array)
    endif

    call compute_aredi(Grid%tmask, grid_length, rossby_radius, agm_array, &
         aredi_array)

    call diagnose_3d(Time, Grid, id_agm_array, agm_array)
    call diagnose_3d(Time, Grid, id_aredi_array, aredi_array)

  end subroutine compute_diffusivity
  ! </SUBROUTINE> NAME="compute_diffusivity"


  !#######################################################################
  ! <SUBROUTINE NAME="compute_agm">
  !
  ! <DESCRIPTION>
  ! Compute the flow-dependent GM diffusivity.
  ! </DESCRIPTION>
  !
  subroutine compute_agm(Domain, Grid, Time, Thickness, Dens, T_temp, T_salt, &
       absslope_z, N2_z, gravity_wave_speed, rossby_radius, &
       ksurf_blayer, grid_length, &
       agm_array)
 
    type(ocean_domain_type),      intent(inout) :: Domain
    type(ocean_grid_type),        intent(in)    :: Grid
    type(ocean_time_type),        intent(in)    :: Time
    type(ocean_thickness_type),   intent(in)    :: Thickness
    type(ocean_density_type),     intent(in)    :: Dens
    type(ocean_prog_tracer_type), intent(in)    :: T_temp
    type(ocean_prog_tracer_type), intent(in)    :: T_salt

    real, dimension(isd:,jsd:,:,0:), intent(in) :: absslope_z
    real, dimension(isd:,jsd:,:,0:), intent(in) :: N2_z
    real, dimension(isd:,jsd:),      intent(in) :: gravity_wave_speed
    real, dimension(isd:,jsd:),      intent(in) :: rossby_radius

    integer, dimension(isd:,jsd:), intent(in) :: ksurf_blayer
    real, dimension(isd:,jsd:),    intent(in) :: grid_length

    real, dimension(isd:,jsd:,:), intent(inout) :: agm_array

    integer :: i, j, k
    real, dimension(isd:ied,jsd:jed,nk) :: agm_fast

    real, dimension(isd:ied,jsd:jed) :: N2_ref
    real, dimension(isd:ied,jsd:jed,nk) :: raw_growth_rate, growth_rate
    real, dimension(isd:ied,jsd:jed) :: growth_rate_zave
    real, dimension(isd:ied,jsd:jed,nk) :: agm_length

    agm_fast = 0.0
    if (agm_closure_micom) then
       do k = 1, nk
          agm_fast(COMP,k) = agm_micom(COMP)
       enddo
    elseif (agm_closure_n2_scale) then
       ! diffusivity computed as coeff*(N/Nref)**2
       if (agm_n2_scale_nref_cst) then
          N2_ref(COMP) = agm_n2_scale_buoy_freq**2
       else
          ! reference buoyancy taken one level beneath blayer base,
          ! and no deeper than one cell from bottom.
          N2_ref(COMP) = 0.0
          do j = jsc, jec
             do i = isc, iec
                if (ksurf_blayer(i,j) > 0 .and. Grid%kmt(i,j) > 1) then
                   k = min(Grid%kmt(i,j)-1, ksurf_blayer(i,j)+1)
                   N2_ref(i,j) = 0.5*(sum(N2_z(i,j,k,:)))
                endif
             enddo
          enddo
       endif
       call diagnose_2d(Time, Grid, id_N2_ref, N2_ref)

       do k = 1, nk
          agm_fast(COMP,k) = agm_n2_scale_coeff*(0.5*sum(N2_z(COMP,k,:))/(epsln+N2_ref(COMP)))
       enddo
    elseif (agm_closure_rate_len2) then
       ! diffusivity computed as growth_rate*length_scale**2
       call compute_raw_growth_rate(absslope_z, N2_z, raw_growth_rate)
       call compute_growth_rate(Domain, Time, Grid, Thickness, T_temp, T_salt, Dens, &
            gravity_wave_speed, rossby_radius, raw_growth_rate, &
            growth_rate, growth_rate_zave)
       call compute_length(Time, Grid, rossby_radius, growth_rate, growth_rate_zave, grid_length, agm_length)

       call diagnose_3d(Time, Grid, id_agm_length, agm_length)
       call diagnose_3d(Time, Grid, id_raw_growth_rate, raw_growth_rate)
       call diagnose_3d(Time, Grid, id_growth_rate, growth_rate)
       call diagnose_2d(Time, Grid, id_growth_rate_zave, growth_rate_zave)

       agm_fast(COMP,:) = agm_scaling*growth_rate(COMP,:)*agm_length(COMP,:)**2
    endif

    if (agm_grid_scaling) then
       call apply_grid_scaling(grid_length, rossby_radius, agm_grid_scaling_power, Grid%tmask, &
            agm_fast)
    endif

    ! Place upper and lower bounds on agm.
    agm_fast(COMP,:) = max(agm_min, min(agm_max, agm_fast(COMP,:)))

    call diagnose_3d(Time, Grid, id_agm_fast, agm_fast)

    ! time damping to get slowly evolving diffusivity
    if (agm_smooth_time) then
       agm_array(COMP,:) = agm_array(COMP,:) - agm_gamma_damp*(agm_array(COMP,:) - agm_fast(COMP,:))
    else
       agm_array(COMP,:) = agm_fast(COMP,:)
    endif

    ! spatial smoothing
    if (agm_smooth_space) then
       call horz_smooth(Domain, Grid%tmask, agm_array)
    endif

    ! need agm_array on full data domain
    agm_array(COMP,:) = agm_array(COMP,:)*Grid%tmask(COMP,:)
    call mpp_update_domains(agm_array(:,:,:), Domain%domain2d)

  end subroutine compute_agm
  ! </SUBROUTINE> NAME="compute_agm"

  !#######################################################################
  ! <SUBROUTINE NAME="compute_raw_growth_rate">
  !
  ! <DESCRIPTION>
  ! Compute the raw growth rate at each grid point, using either the eady
  ! growth rate or the baroclinicity to obtain raw_growth_rate = NS
  ! </DESCRIPTION>
  !
  subroutine compute_raw_growth_rate(absslope_z, N2_z, raw_growth_rate)

    real, dimension(isd:,jsd:,:,0:), intent(in) :: absslope_z
    real, dimension(isd:,jsd:,:,0:), intent(in) :: N2_z
    real, dimension(isd:,jsd:,:),    intent(out) :: raw_growth_rate

    ! average over vertical quarter cells to get raw Eady growth rate,
    ! NS, with regularization not yet applied.
    if (agm_rate_baroclinic) then
       raw_growth_rate(COMP,:) = sum(N2_z(COMP,:,:)*absslope_z(COMP,:,:), dim=4)/(epsln + 2.0*agm_rate_baro_buoy_freq)
    elseif (agm_rate_eady) then
       raw_growth_rate(COMP,:) = sum(sqrt(N2_z(COMP,:,:))*absslope_z(COMP,:,:), dim=4)/2.0
    endif
    
  end subroutine compute_raw_growth_rate
  ! </SUBROUTINE> NAME="compute_raw_growth_rate"

  !#######################################################################
  ! <SUBROUTINE NAME="compute_growth_rate">
  !
  ! <DESCRIPTION>
  ! Take the raw growth rate and convert it to a final growth rate to be
  ! used in the diffusivity calculations.
  ! </DESCRIPTION>
  !
  subroutine compute_growth_rate(Domain, Time, Grid, Thickness, T_temp, T_salt, Dens, &
       gravity_wave_speed, rossby_radius, raw_growth_rate, &
       growth_rate, growth_rate_zave)

    type(ocean_domain_type),    intent(inout) :: Domain
    type(ocean_time_type),      intent(in)    :: Time
    type(ocean_grid_type),      intent(in)    :: Grid
    type(ocean_thickness_type), intent(in)    :: Thickness
    type(ocean_density_type),   intent(in)    :: Dens

    real, dimension(isd:,jsd:),   intent(in) :: gravity_wave_speed
    real, dimension(isd:,jsd:),   intent(in) :: rossby_radius
    type(ocean_prog_tracer_type), intent(in) :: T_temp
    type(ocean_prog_tracer_type), intent(in) :: T_salt
    real, dimension(isd:,jsd:,:), intent(in) :: raw_growth_rate
    real, dimension(isd:,jsd:,:), intent(out) :: growth_rate
    real, dimension(isd:,jsd:),   intent(out) :: growth_rate_zave

    integer :: k, tau
    real, dimension(isd:ied,jsd:jed) :: ml_depth
    real, dimension(isd:ied,jsd:jed) :: int_dz
    real, dimension(isd:ied,jsd:jed) :: int_rate_dz
    real, dimension(isd:ied,jsd:jed) :: ave_ml_rate

    real, dimension(isd:ied,jsd:jed) :: eg_rate

    growth_rate(COMP,:) = raw_growth_rate(COMP,:)

    ! apply cap to the growth growth rate
    if (agm_rate_cap) then
       do k = 1, nk
          growth_rate(COMP,k) = HARM_MEAN(growth_rate(COMP,k), growth_rate_max(COMP))
       enddo
    endif

    ! vertically average growth rate within surface mixed layer
    growth_rate(COMP,:) = Grid%tmask(COMP,:)*growth_rate(COMP,:)
    if (agm_rate_ave_mixed) then
       tau = Time%tau
       ml_depth(COMP) = 0.0
       call calc_mixed_layer_depth(Thickness,    &
            T_salt%field(:,:,:,tau), T_temp%field(:,:,:,tau), &
            Dens%rho(:,:,:,tau), Dens%pressure_at_depth(:,:,:), &
            ml_depth)
       ml_depth(COMP) = Grid%tmask(COMP,1)*min(ml_depth(COMP), Grid%ht(COMP))

       ! do not believe growth_rate at k=1, so always skip its contribution
       int_dz(COMP) = 0.0
       int_rate_dz(COMP) = 0.0
       do k = 2, nk
          where (Thickness%depth_zt(COMP,k) <= ml_depth(COMP))
             int_dz(COMP) = int_dz(COMP) + Thickness%dzt(COMP,k)
             int_rate_dz(COMP) = int_rate_dz(COMP) + growth_rate(COMP,k)*Thickness%dzt(COMP,k)
          endwhere
       enddo
       ave_ml_rate(COMP) = int_rate_dz(COMP)/(int_dz(COMP)+epsln)

       call diagnose_2d(Time, Grid, id_ml_depth, ml_depth)
       call diagnose_2d(Time, Grid, id_ave_ml_rate, ave_ml_rate)

       do k = 1, nk
          where (Thickness%depth_zt(COMP,k) <= ml_depth(COMP))
             growth_rate(COMP,k) = Grid%tmask(COMP,k)*ave_ml_rate(COMP)
          endwhere
       enddo
    endif

    ! apply vertical 1-2-1 smoothing
    if (agm_rate_smooth_vert) then
       call vert_smooth(Grid%kmt(:,:), growth_rate)
    endif

    ! apply horizontal 1-2-1 smoothing
    if (agm_rate_smooth_horiz) then
       call mpp_update_domains(growth_rate(:,:,:), Domain%domain2d)
       do k = 1, nk
          growth_rate(:,:,k) = S2D(growth_rate(:,:,k))
       enddo
    endif

    ! FIXME: need to find some documentation of this method.
    if (agm_rate_eden_greatbatch) then
       eg_rate(COMP) = agm_rate_eg_alpha*gravity_wave_speed(COMP)/(rossby_radius(COMP) + epsln)
       call diagnose_2d(Time, Grid, id_eg_rate, eg_rate)
       do k = 1, nk
          growth_rate(COMP,k) = HARM_MEAN(growth_rate(COMP,k), eg_rate(COMP))
       enddo
    endif

    ! compute vertical average over specified depth
    call vertical_average(Thickness, Grid%tmask, growth_rate, growth_rate_zave)

    if (agm_rate_zave) then
       do k = 1, nk
          growth_rate(COMP,k) = growth_rate_zave(COMP)
       enddo
    endif

  end subroutine compute_growth_rate
  ! </SUBROUTINE> NAME="compute_growth_rate"  

  !#######################################################################
  ! <SUBROUTINE NAME="vertical_average">
  !
  ! <DESCRIPTION>
  ! Compute the vertical average of the given array between D_t and D_b, 
  ! as specified by agm_rate_upper_depth and agm_rate_lower_depth.
  ! </DESCRIPTION>
  !
  subroutine vertical_average(Thickness, tmask, array, array_ave)
    type(ocean_thickness_type),   intent(in) :: Thickness
    real, dimension(isd:,jsd:,:), intent(in) :: tmask
    real, dimension(isd:,jsd:,:), intent(in) :: array
    real, dimension(isd:,jsd:),   intent(out) :: array_ave

    integer :: k 
    real, dimension(isd:ied,jsd:jed) :: int_dz
    real, dimension(isd:ied,jsd:jed) :: int_array_dz

    int_array_dz(COMP) = 0.0
    int_dz(COMP) = 0.0
    do k = 1, nk-1
       where (Thickness%depth_zt(COMP,k) >= agm_rate_upper_depth .and. &
            Thickness%depth_zt(COMP,k) <= agm_rate_lower_depth)
          int_array_dz(COMP) = int_array_dz(COMP) + array(COMP,k)*Thickness%dzt(COMP,k)
          int_dz(COMP) = int_dz(COMP) + Thickness%dzt(COMP,k)
       endwhere
    enddo
    array_ave(COMP) = tmask(COMP,1)*int_array_dz(COMP)/(int_dz(COMP)+epsln)

  end subroutine vertical_average
  ! </SUBROUTINE> NAME="vertical_average"  


  !#######################################################################
  ! <SUBROUTINE NAME="compute_length">
  !
  ! <DESCRIPTION>
  ! Compute the flow-dependent length scale involved in the GM diffusivity
  ! calculations.
  ! </DESCRIPTION>
  !
  subroutine compute_length(Time, Grid, rossby_radius, growth_rate, growth_rate_zave, grid_length, agm_length)

    type(ocean_time_type),        intent(in) :: Time
    type(ocean_grid_type),        intent(in) :: Grid
    real, dimension(isd:,jsd:),   intent(in) :: rossby_radius
    real, dimension(isd:,jsd:,:), intent(in) :: growth_rate
    real, dimension(isd:,jsd:),   intent(in) :: growth_rate_zave
    real, dimension(isd:,jsd:),   intent(in) :: grid_length
    real, dimension(isd:,jsd:,:), intent(out) :: agm_length

    integer :: k
    real, dimension(isd:ied,jsd:jed,nk) :: rhines_length
    real, dimension(isd:ied,jsd:jed) :: bczone_radius

    ! set the length scale
    agm_length = 0.0
    if (agm_length_fixed) then
       agm_length(COMP,:) = agm_length
    elseif (agm_length_rossby) then
       agm_length(COMP,:) = spread(HARM_MEAN(grid_length(COMP), rossby_radius(COMP)), 3, nk)
    elseif (agm_length_bczone) then
       call compute_bczone_radius(growth_rate_zave, grid_length, bczone_radius)
       call diagnose_2d(Time, Grid, id_bczone_radius, bczone_radius)
       agm_length(COMP,:) = spread(HARM_MEAN(grid_length(COMP), bczone_radius(COMP)), 3, nk)
    elseif (agm_length_eden_greatbatch) then
       do k = 1, nk
          rhines_length(COMP,k) = growth_rate(COMP,k)/beta_param(COMP)
          agm_length(COMP,k) = min(rhines_length(COMP,k), rossby_radius(COMP))
          agm_length(COMP,k) = HARM_MEAN(grid_length(COMP), agm_length(COMP,k))
       enddo
       call diagnose_3d(Time, Grid, id_rhines_length, rhines_length)
    endif

    if (agm_length_cap) then
       agm_length(COMP,:) = min(agm_length_max, agm_length(COMP,:))
    endif

  end subroutine compute_length
  ! </SUBROUTINE> NAME="compute_length"

  !#######################################################################
  ! <SUBROUTINE NAME="compute_bczone_radius">
  !
  ! <DESCRIPTION>
  ! Subroutine computes the radius of the baroclinic zone in a manner
  ! suggested by the Hadley Centre approach (Malcolm Roberts, personal
  ! communication).
  !
  ! Algorithm is used in MOM3 and documented in the MOM3 Manual.
  ! </DESCRIPTION>
  !
  subroutine compute_bczone_radius(growth_rate_zave, grid_length, bczone_radius)

    real, dimension(isd:,jsd:), intent(in) :: growth_rate_zave
    real, dimension(isd:,jsd:), intent(in) :: grid_length
    real, dimension(isd:,jsd:), intent(out) :: bczone_radius

    integer :: i, j
    real :: n_zone, e_zone, s_zone, w_zone, fract, nstot, ewtot
    integer :: ip, jq, max_pts

    max_pts = agm_bczone_max_pts
    bczone_over = 0
    where (growth_rate_zave(COMP) > agm_bczone_crit_rate) bczone_over(COMP) = 1
    call mpp_update_domains(bczone_over(:,:), BCzone_domain%domain2d)

    ! Set default length for points outside baroclinic region.
    where (bczone_over(COMP) == 0) bczone_radius(COMP) = grid_length(COMP)

    do j = jsc, jec
       do i = isc, iec
          if (bczone_over(i,j)==1) then
             ! search northward
             n_zone = bczone_dyt(i,j)
             do jq = j+1, j+max_pts
                if (bczone_over(i,jq)==1) then
                   n_zone = n_zone + bczone_dyt(i,jq)
                else
                   exit
                endif
             enddo

             ! search southward
             s_zone = bczone_dyt(i,j)
             do jq = j-1, j-max_pts, -1
                if (bczone_over(i,jq)==1) then
                   s_zone = s_zone + bczone_dyt(i,jq)
                else
                   exit
                endif
             enddo

             ! search eastward
             e_zone = bczone_dxt(i,j)
             do ip = i+1, i+max_pts
                if (bczone_over(ip,j)==1) then
                   e_zone = e_zone + bczone_dxt(ip,j)
                else
                   exit
                endif
             enddo

             ! search westward
             w_zone = bczone_dxt(i,j)
             do ip = i-1, i-max_pts, -1
                if (bczone_over(ip,j)==1) then
                   w_zone = w_zone + bczone_dxt(ip,j)
                else
                   exit
                endif
             enddo

             ! total radius (subtraction accounts for double-counting central point)
             nstot = n_zone + s_zone - bczone_dyt(i,j)
             ewtot = e_zone + w_zone - bczone_dxt(i,j)

             if (nstot < ewtot) then
                fract = min(n_zone,s_zone)/max(n_zone,s_zone)
                bczone_radius(i,j) = fract*nstot
             else
                fract = min(e_zone,w_zone)/max(e_zone,w_zone)
                bczone_radius(i,j) = fract*ewtot
             endif
          endif
       enddo
    enddo

  end subroutine compute_bczone_radius
  ! </SUBROUTINE> NAME="compute_bczone_radius"

  !#######################################################################
  ! <SUBROUTINE NAME="compute_aredi">
  !
  ! <DESCRIPTION>
  ! Compute the flow-dependent neutral (Redi) diffusivity.
  ! </DESCRIPTION>
  !
  subroutine compute_aredi(tmask, grid_length, rossby_radius, agm_array, &
       aredi_array)

    real, dimension(isd:,jsd:,:), intent(in) :: tmask
    real, dimension(isd:,jsd:),   intent(in) :: grid_length
    real, dimension(isd:,jsd:),   intent(in) :: rossby_radius
    real, dimension(isd:,jsd:,:), intent(in) :: agm_array

    real, dimension(isd:,jsd:,:), intent(inout) :: aredi_array

    if (aredi_equal_agm) then
       aredi_array(COMP,:) = agm_array(COMP,:)
    elseif (aredi_grid_scaling) then
       aredi_array(COMP,:) = aredi*tmask(COMP,:)
       call apply_grid_scaling(grid_length, rossby_radius, aredi_grid_scaling_power, tmask, &
            aredi_array)
    endif

  end subroutine compute_aredi
  ! </SUBROUTINE> NAME="compute_aredi"

  !#######################################################################
  ! <SUBROUTINE NAME="apply_grid_scaling">
  !
  ! <DESCRIPTION>
  ! Scale the supplied array as a function of the grid length and the
  ! Rossby radius. The scaling factor takes a value between zero and one.
  !
  ! Interesting values of the scaling factor are
  ! 1.0 if Rossby radius = 0.0
  ! 0.5 if Rossby radius = grid length
  ! -> 0 as rossby radius -> inf
  ! </DESCRIPTION>
  !
  subroutine apply_grid_scaling(grid_length, rossby_radius, power, tmask, &
       array)

    real, dimension(isd:,jsd:),   intent(in) :: grid_length
    real, dimension(isd:,jsd:),   intent(in) :: rossby_radius
    real,                         intent(in) :: power
    real, dimension(isd:,jsd:,:), intent(in) :: tmask
    real, dimension(isd:,jsd:,:), intent(inout) :: array

    real, dimension(isd:ied,jsd:jed) :: grid_scaling
    integer :: k

    grid_scaling(COMP) = 1/(1 + (rossby_radius(COMP)/grid_length(COMP))**power)
    grid_scaling(COMP) = tmask(COMP,1)*grid_scaling(COMP)
    do k = 1, nk
       array(COMP,k) = array(COMP,k)*grid_scaling(COMP)
    enddo

  end subroutine apply_grid_scaling
  ! </SUBROUTINE> NAME="apply_grid_scaling"

end module ocean_nphysics_diff_mod

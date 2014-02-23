module ocean_nphysics_skew_mod
  !
  !<CONTACT EMAIL="t.leslie@unsw.edu.au"> Tim Leslie
  !</CONTACT>
  !
  !<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
  !</CONTACT>
  ! 
  !<OVERVIEW>
  ! Thickness weighted and density weighted time tendency for tracer
  ! from Laplacian neutral diffusion + Laplacian skew-diffusion.
  !</OVERVIEW>
  !
  !<DESCRIPTION>
  ! This module computes the cell thickness weighted and density
  ! weighted tracer tendency from small angle Laplacian neutral diffusion
  ! plus Laplacian skew-diffusion.  The algorithms for neutral diffusion
  ! are based on mom4p0d methods.  The algorithm for neutral skewsion
  ! are based on a projection onto a few of the lowest baroclinic
  ! modes. This module is experimental, and should be used with caution.
  !</DESCRIPTION>
  !
  ! <INFO>
  !
  ! <REFERENCE>
  ! S.M. Griffies
  ! The Gent-McWilliams Skew-flux
  ! Journal of Physical Oceanography (1998) vol 28 pages 831-841
  ! </REFERENCE>
  !
  ! <REFERENCE>
  ! R. Ferrari, S.M. Griffies, A.J.G. Nurser, and G.K. Vallis
  ! A boundary value problem for the parameterized mesoscale eddy transport
  ! Ocean Modelling, 2009.
  ! </REFERENCE>
  !
  ! <REFERENCE>
  ! S.M. Griffies
  ! Fundamentals of Ocean Climate Models (2004)
  ! Princeton University Press
  ! </REFERENCE>
  !
  ! <REFERENCE>
  ! S.M. Griffies: Elements of MOM (2012)
  ! </REFERENCE>
  !
  ! </INFO>
  !
  !<NAMELIST NAME="ocean_nphysics_skew_nml">
  !  <DATA NAME="gm_transport" TYPE="logical">
  !  To compute tendency from GM skewsion. Default gm_transport=.false.
  !  </DATA>
  !  <DATA NAME="bc_modes_transport" TYPE="logical">
  !  To compute tendency from GM skewsion using streamfunction established
  !  by baroclinic modes. Default bc_modes_transport=.false.
  !  </DATA>
  !  <DATA NAME="bvp_transport" TYPE="logical">
  !  To compute tendency from GM skewsion using streamfunction established
  !  by a boundary value problem. Default bvp_transport=.false.
  !  </DATA>
  !
  !  <DATA NAME="number_bc_modes" TYPE="integer">
  !  The number of baroclinic modes used to construct the eddy induced
  !  streamfunction when bc_modes_transport. Default number_bc_modes=1.
  !  </DATA>
  !  <DATA NAME="min_bc_speed"  UNITS="m/s"  TYPE="real">
  !  The minimum speed used for computing the baroclinic modes.
  !  Default min_bc_speed=1e-6
  !  </DATA>
  !  <DATA NAME="regularize_transport" TYPE="logical">
  !  To reduce the magnitude of psi in regions of weak stratification,
  !  using the slope = smax_psi to set the overall scale of the max allowed
  !  for psi. Default regularize_transport=.true.
  !  </DATA>
  !  <DATA NAME="smax_psi" TYPE="real">
  !  Maximum slope used for setting the overall scale of a modal
  !  contribution to the parameterized transport.
  !  Default smax_psi=0.1.
  !  </DATA>
  !
  !  <DATA NAME="bvp_constant_speed" TYPE="logical">
  !  For taking a constant speed to be used for the calculation
  !  of the BVP streamfunction. Default bvp_constant_speed=.false.
  !  </DATA>
  !  <DATA NAME="bvp_speed" UNITS="m/s" TYPE="real">
  !  For setting the speed weighting the second order derivative operator
  !  in the BVP streamfunction method:
  !  c^2 = max[bvp_min_speed, (bvp_speed-c_mode)^2].
  !  If bvp_constant_speed, then  c^2 = bvp_speed^2.
  !  Default bvp_speed=0.0, in which case c^2 = c_mode^2.
  !  </DATA>
  !
  !  <DATA NAME="bvp_bc_mode" TYPE="integer">
  !  The particular baroclinic mode used to construct the BVP streamfunction.
  !  If bvp_bc_mode=0, then will set bc_speed=0 when computing the BVP streamfunction.
  !  Default bvp_bc_mode=1.
  !  </DATA>
  !  <DATA NAME="bvp_min_speed" UNITS="m/s" TYPE="real">
  !  For setting a minimum speed for use with the calculation
  !  of the BVP streamfunction. We need  bvp_min_speed>0 to ensure
  !  that the second order derivative operator contributes to the
  !  calculation of the streamfunction.
  !  Default bvp_min_speed=0.1.
  !  </DATA>
  !
  !  <DATA NAME="bv_freq_smooth_vert" TYPE="logical">
  !  To smooth the buoyancy frequency for use in
  !  computing the baroclinic modes. Generally this field has already
  !  been smooted in ocean_density_mod, but we maintain the possibility of
  !  further smoothing here.  Default bv_freq_smooth_vert=.false.
  !  </DATA>
  !
  !  <DATA NAME="smooth_transport" TYPE="logical">
  !  For doing a horizontal 1-2-1 smoothing on the psix and psiy fields.
  !  This is useful to reduce noise. Default smooth_psi=.true.
  !  </DATA>
  !
  !</NAMELIST>

#define COMP isc:iec,jsc:jec
#define COMPXRR isc+1:iec+1,jsc:jec
#define COMPYRR isc:iec,jsc+1:jec+1

  use constants_mod, only: epsln, pi
  use fms_mod,       only: open_namelist_file, check_nml_error, close_file, FATAL
  use mpp_mod,       only: mpp_error, stdout, stdlog
  use mpp_mod,       only: input_nml_file

  use ocean_domains_mod,           only: get_local_indices
  use ocean_nphysics_util_new_mod, only: horz_smooth
  use ocean_parameters_mod,        only: grav
  use ocean_types_mod,             only: ocean_grid_type, ocean_domain_type
  use ocean_types_mod,             only: ocean_thickness_type, ocean_time_type
  use ocean_util_mod,              only: diagnose_2d, diagnose_3d
  use ocean_util_mod,              only: register_2d_t_field, register_3d_t_field

  implicit none

  public ocean_nphysics_skew_init
  public skew_transport

  private

#include <ocean_memory.h>

character(len=128) :: version=&
     '$Id: ocean_nphysics_skew.F90,v 20.0 2013/12/14 00:14:48 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

character(len=*), parameter :: FILENAME=&
     __FILE__


  integer :: id_bc_speed
  integer :: id_sine_arg_up
  integer :: id_sine_arg_lo
  integer :: id_sine_arg_pre_up
  integer :: id_sine_arg_pre_lo
  integer :: id_bc_mode_1
  integer :: id_norm_1
  integer :: id_depth_blayer_base
  integer :: id_agm_slopex_max
  integer :: id_agm_slopey_max
  integer :: id_N2_Sx_agm_dz
  integer :: id_N2_Sy_agm_dz

  !------Namelist parameters 

  logical :: gm_transport = .false.

  ! for computing skewsion via baroclinic mode decomposition  
  logical :: bc_modes_transport = .false.
  ! number of modes used to construct GM transport when bc_modes_transport=.true.
  integer :: number_bc_modes = 1

  ! for regularizing the transport when computed from modes
  logical :: regularize_transport = .true.
  real    :: smax_psi = 0.01

  ! for computing skewsion via boundary value problem (BVP)
  logical :: bvp_transport = .false.

  ! for setting the speed to constant with the BVP calculation.
  logical :: bvp_constant_speed = .false.
  real    :: bvp_speed          = 0.0

  ! the particular baroclinic mode speed used to construct the BVP transport.
  integer :: bvp_bc_mode = 1
  real    :: bvp_min_speed = 1.0

  ! for smoothing upsilonx and upsilony
  logical :: smooth_transport = .true.

  namelist /ocean_nphysics_skew_nml/                                      & 
       gm_transport, bc_modes_transport, number_bc_modes,                 &
       regularize_transport, smax_psi, bvp_transport, bvp_constant_speed, &
       bvp_speed, bvp_bc_mode, bvp_min_speed, smooth_transport

  contains

  !#######################################################################
  ! <SUBROUTINE NAME="ocean_nphysics_skew_init">
  !
  ! <DESCRIPTION>
  ! Initialization routine. 
  ! </DESCRIPTION>
  !
  subroutine ocean_nphysics_skew_init(Domain, Grid, Time, agm_z_dep)

    type(ocean_domain_type), intent(in) :: Domain
    type(ocean_grid_type),   intent(in) :: Grid
    type(ocean_time_type),   intent(in) :: Time
    logical,                 intent(in) :: agm_z_dep

    integer :: stdoutunit, stdlogunit, ierr, ioun, io_status
    stdoutunit=stdout();stdlogunit=stdlog() 

#ifndef MOM_STATIC_ARRAYS
    call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
    nk = Grid%nk
#endif

    ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_nphysics_skew_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_nphysics_skew_nml')
#else
    ioun =  open_namelist_file()
    read (ioun,ocean_nphysics_skew_nml,IOSTAT=io_status)
    ierr = check_nml_error(io_status,'ocean_nphysics_skew_nml')
    call close_file (ioun)
#endif
    write (stdoutunit,'(/)')
    write (stdoutunit,ocean_nphysics_skew_nml)  
    write (stdlogunit,ocean_nphysics_skew_nml)

    if (count((/gm_transport, bc_modes_transport, bvp_transport /)) /= 1) then
       call mpp_error(FATAL, &
            '==>Error from ocean_nphysics_skew: must choose exactly one neutral skew transport method.')
    endif

    if (agm_z_dep .and. (bc_modes_transport .or. bvp_transport)) then
       call mpp_error(FATAL, &
            '==>Error from ocean_nphysics_skew (ocean_nphysics_skew_init): Cannot use a z-dependent ' // &
            'skew diffusivity with bc_modes_transport or bvp_transport.')
    endif

    call register_fields(Grid, Time)

  end subroutine ocean_nphysics_skew_init
  ! </SUBROUTINE> NAME="ocean_nphysics_skew_init"

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

    id_norm_1 = register_2d_t_field(Grid, Time, 'norm_1', &
         'Normalisation factor for 1st baroclinic mode.', 'm/s', (/-1e20,1.e20/))
    id_sine_arg_up = register_3d_t_field(Grid, Time, 'sine_arg_up', &
         'Argument to sine used in computing baroclinic modes', 'rad', (/-2*pi, 2*pi/))
    id_sine_arg_lo = register_3d_t_field(Grid, Time, 'sine_arg_lo', &
         'Argument to sine used in computing baroclinic modes', 'rad', (/-2*pi, 2*pi/))
    id_sine_arg_pre_up = register_3d_t_field(Grid, Time, 'sine_arg_pre_up', &
         'Argument to sine used in computing baroclinic modes', 'rad', (/-2*pi, 2*pi/))
    id_sine_arg_pre_lo = register_3d_t_field(Grid, Time, 'sine_arg_pre_lo', &
         'Argument to sine used in computing baroclinic modes', 'rad', (/-2*pi, 2*pi/))
    id_bc_mode_1 = register_3d_t_field(Grid, Time, 'bc_mode_1', &
         '1st baroclinic mode', '', (/-1e8,1.e8/)) ! Dimensionless

    id_depth_blayer_base = register_2d_t_field(Grid, Time, 'depth_blayer_base', &
         'Depth of the boundary layer base used to taper GM transport', 'm', (/0.0, 1.e5/))
    id_agm_slopex_max = register_2d_t_field(Grid, Time, 'agm_slopex_max', &
         'x component of (kS)_max at boundary layer', 'm^2/s', (/-1e8,1.e8/))
    id_agm_slopey_max = register_2d_t_field(Grid, Time, 'agm_slopey_max', &
         'y component of (kS)_max at boundary layer', 'm^2/s', (/-1e8,1.e8/))

    id_N2_Sx_agm_dz = register_3d_t_field(Grid, Time, 'N2_Sx_agm_dz', &
         'x component of \kappa \grad\rho dz', 'rad', (/-1e20, 1e20/))
    id_N2_Sy_agm_dz = register_3d_t_field(Grid, Time, 'N2_Sy_agm_dz', &
         'x component of \kappa \grad\rho dz', 'rad', (/-1e20, 1e20/))

  end subroutine register_fields
  ! </SUBROUTINE> NAME="register_fields"

  subroutine skew_transport(Domain, Time, Grid, Thickness,               &
       agm_array, slopex_xz, slopey_yz, Ndz_z, N2_z, gravity_wave_speed, &
       depth_taper, ksurf_blayer,                                        &
       upsilonx_xz, upsilony_yz)

    type(ocean_domain_type),    intent(inout) :: Domain
    type(ocean_time_type),      intent(in)    :: Time
    type(ocean_grid_type),      intent(in)    :: Grid
    type(ocean_thickness_type), intent(in)    :: Thickness

    real, dimension(isd:,jsd:,:),       intent(in) :: agm_array
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: slopex_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: slopey_yz
    real, dimension(isd:,jsd:,:,0:),    intent(in) :: Ndz_z
    real, dimension(isd:,jsd:,:,0:),    intent(in) :: N2_z
    real, dimension(isd:,jsd:),         intent(in) :: gravity_wave_speed

    real, dimension(isd:,jsd:,:),  intent(in) :: depth_taper
    integer, dimension(isd:,jsd:), intent(in) :: ksurf_blayer

    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: upsilonx_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: upsilony_yz

    integer :: ip, kr

    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: agm_slopex_xz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: agm_slopey_yz
    real, dimension(isd:ied,jsd:jed,nk,0:1)     :: N2dz_z
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: N2_Sx_agm_dz_xz, N2_Sy_agm_dz_yz

    upsilonx_xz(:,:,:,:,:) = 0.0
    upsilony_yz(:,:,:,:,:) = 0.0

    do ip=0,1
       do kr=0,1
          agm_slopex_xz(COMP,:,ip,kr) = agm_array(COMP,:)*slopex_xz(COMP,:,ip,kr)
          agm_slopey_yz(COMP,:,ip,kr) = agm_array(COMP,:)*slopey_yz(COMP,:,ip,kr)
       enddo
    enddo
    if(gm_transport) then
       call gm_tensor(Time, Grid, Thickness, agm_slopex_xz, agm_slopey_yz, &
            depth_taper, ksurf_blayer, &
            upsilonx_xz, upsilony_yz)
    else
       ! We need the following quantities for the vertical integrals which follow
       N2dz_z(COMP,:,0) = N2_z(COMP,:,0)*Thickness%dztup(COMP,:)
       N2dz_z(COMP,:,1) = N2_z(COMP,:,1)*Thickness%dztlo(COMP,:)
       do ip=0,1
          N2_Sx_agm_dz_xz(COMP,:,ip,:) = N2dz_z(COMP,:,:)*agm_slopex_xz(COMP,:,ip,:)
          N2_Sy_agm_dz_yz(COMP,:,ip,:) = N2dz_z(COMP,:,:)*agm_slopey_yz(COMP,:,ip,:)
       enddo
       call diagnose_3d(Time, Grid, id_N2_Sx_agm_dz, 0.25*sum(sum(N2_Sx_agm_dz_xz, dim=5), dim=4))
       call diagnose_3d(Time, Grid, id_N2_Sy_agm_dz, 0.25*sum(sum(N2_Sy_agm_dz_yz, dim=5), dim=4))
       
       if (bc_modes_transport) then
          call compute_transport_modes(Time, Grid, agm_array, Ndz_z, N2dz_z, &
               gravity_wave_speed, N2_Sx_agm_dz_xz, N2_Sy_agm_dz_yz,         &
               upsilonx_xz, upsilony_yz)
       elseif (bvp_transport) then
          call compute_transport_bvp(Grid%tmask, Thickness, gravity_wave_speed, &
               N2dz_z, N2_Sx_agm_dz_xz, N2_Sy_agm_dz_yz, upsilonx_xz, upsilony_yz)
       endif
    endif

    if (smooth_transport) then
       call do_smooth_transport(Domain, Grid%tmask, upsilonx_xz, upsilony_yz)
    endif

  end subroutine skew_transport

  !#######################################################################
  ! <SUBROUTINE NAME="gm_tensor">
  !
  ! <DESCRIPTION>
  ! Compute the skew tensor terms using the GM method.
  ! </DESCRIPTION>
  !
  subroutine gm_tensor(Time, Grid, Thickness, agm_slopex_xz, agm_slopey_yz, &
       depth_taper, ksurf_blayer, upsilonx_xz, upsilony_yz)

    type(ocean_time_type),              intent(in) :: Time
    type(ocean_grid_type),              intent(in) :: Grid
    type(ocean_thickness_type),         intent(in) :: Thickness
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: agm_slopex_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: agm_slopey_yz

    real, dimension(isd:,jsd:,:),  intent(in) :: depth_taper
    integer, dimension(isd:,jsd:), intent(in) :: ksurf_blayer

    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: upsilonx_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: upsilony_yz

    real, dimension(isd:ied,jsd:jed,nk) :: agm_slope_taper_x
    real, dimension(isd:ied,jsd:jed,nk) :: agm_slope_taper_y
    real, dimension(isd:ied,jsd:jed)    :: agm_slopex_max
    real, dimension(isd:ied,jsd:jed)    :: agm_slopey_max
    real, dimension(isd:ied,jsd:jed)    :: depth_blayer_base
    real, dimension(isd:ied,jsd:jed)    :: agm_slopex_blayer
    real, dimension(isd:ied,jsd:jed)    :: agm_slopey_blayer
    real, dimension(isd:ied,jsd:jed)    :: depth_taper_blayer

    integer :: i, j, k, k_blayer

    ! Calculate boundary layer values
    do j=jsc,jec
       do i=isc,iec
          k_blayer = ksurf_blayer(i,j)
          if (k_blayer > 0) then
             depth_blayer_base(i,j)  = Thickness%depth_zt(i,j,k_blayer)
             agm_slopex_blayer(i,j)  = sum(agm_slopex_xz(i,j,k_blayer,:,:))/4.0
             agm_slopey_blayer(i,j)  = sum(agm_slopey_yz(i,j,k_blayer,:,:))/4.0
             depth_taper_blayer(i,j) = depth_taper(i,j,k_blayer)
          else
             depth_blayer_base(i,j)  = 0.0
             agm_slopex_blayer(i,j)  = 0.0
             agm_slopey_blayer(i,j)  = 0.0
             depth_taper_blayer(i,j) = 0.0
          endif
       enddo
    enddo

    agm_slopex_max(COMP) = 0.0
    agm_slopey_max(COMP) = 0.0
    where (ksurf_blayer(COMP) > 0)
       agm_slopex_max(COMP) = depth_taper_blayer(COMP)*agm_slopex_blayer(COMP)
       agm_slopey_max(COMP) = depth_taper_blayer(COMP)*agm_slopey_blayer(COMP)
    endwhere
    
    call diagnose_2d(Time, Grid, id_depth_blayer_base, depth_blayer_base)
    call diagnose_2d(Time, Grid, id_agm_slopex_max, agm_slopex_max, abs_max=1.e5)
    call diagnose_2d(Time, Grid, id_agm_slopey_max, agm_slopey_max, abs_max=1.e5)

    do k = 1, nk
       ! taper times slope for use with GM skewsion
       where (k <= ksurf_blayer(COMP))
          ! Linear taper of (kS)_max in the steep slope region
          agm_slope_taper_x(COMP,k) = agm_slopex_max(COMP)*Thickness%depth_zt(COMP,k)/depth_blayer_base(COMP)
          agm_slope_taper_y(COMP,k) = agm_slopey_max(COMP)*Thickness%depth_zt(COMP,k)/depth_blayer_base(COMP)
          upsilonx_xz(COMP,k,0,0)   = -agm_slope_taper_x(COMP,k)
          upsilonx_xz(COMP,k,0,1)   = -agm_slope_taper_x(COMP,k)
          upsilonx_xz(COMP,k,1,0)   = -agm_slope_taper_x(COMP,k)
          upsilonx_xz(COMP,k,1,1)   = -agm_slope_taper_x(COMP,k)

          upsilony_yz(COMP,k,0,0) = -agm_slope_taper_y(COMP,k)
          upsilony_yz(COMP,k,0,1) = -agm_slope_taper_y(COMP,k)
          upsilony_yz(COMP,k,1,0) = -agm_slope_taper_y(COMP,k)
          upsilony_yz(COMP,k,1,1) = -agm_slope_taper_y(COMP,k)
       elsewhere
          ! Sine taper of (kS) elsewhere
          upsilonx_xz(COMP,k,0,0) = -depth_taper(COMP,k)*agm_slopex_xz(COMP,k,0,0)
          upsilonx_xz(COMP,k,0,1) = -depth_taper(COMP,k)*agm_slopex_xz(COMP,k,0,1)
          upsilonx_xz(COMP,k,1,0) = -depth_taper(COMP,k)*agm_slopex_xz(COMP,k,1,0)
          upsilonx_xz(COMP,k,1,1) = -depth_taper(COMP,k)*agm_slopex_xz(COMP,k,1,1)

          upsilony_yz(COMP,k,0,0) = -depth_taper(COMP,k)*agm_slopey_yz(COMP,k,0,0)
          upsilony_yz(COMP,k,0,1) = -depth_taper(COMP,k)*agm_slopey_yz(COMP,k,0,1)
          upsilony_yz(COMP,k,1,0) = -depth_taper(COMP,k)*agm_slopey_yz(COMP,k,1,0)
          upsilony_yz(COMP,k,1,1) = -depth_taper(COMP,k)*agm_slopey_yz(COMP,k,1,1)
       endwhere
    enddo

  end subroutine gm_tensor
  ! </SUBROUTINE> NAME="gm_tensor"

  !#######################################################################
  ! <SUBROUTINE NAME="compute_transport_modes">
  !
  ! <DESCRIPTION>
  ! Compute transport as projection onto baroclinic modes.
  !
  ! Units of upsilon are m^2/sec
  !
  ! Subroutine computes the baroclinic wave speeds and the dimensionless
  ! baroclinic mode eigenfunction for the vertical velocity baroclinic
  ! modes.  These modes vanish at the surface and the bottom.  We use
  ! the Chelton etal WKB analytic formulae for the speeds and modes.
  !
  ! The baroclinic modes are dimensionless, and normalized over the
  ! depth of the ocean, from free surface to bottom.
  !
  ! The speeds are m/sec.
  !
  ! </DESCRIPTION>
  !
  subroutine compute_transport_modes(Time, Grid, agm_array, Ndz_z, N2dz_z, &
       gravity_wave_speed, N2_Sx_agm_dz_xz, N2_Sy_agm_dz_yz, upsilonx_xz, upsilony_yz)
    type(ocean_time_type),              intent(in) :: Time
    type(ocean_grid_type),              intent(in) :: Grid
    real, dimension(isd:,jsd:,:),       intent(in) :: agm_array
    real, dimension(isd:,jsd:,:,0:),    intent(in) :: Ndz_z
    real, dimension(isd:,jsd:,:,0:),    intent(in) :: N2dz_z
    real, dimension(isd:,jsd:),         intent(in) :: gravity_wave_speed
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: N2_Sx_agm_dz_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: N2_Sy_agm_dz_yz

    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: upsilonx_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: upsilony_yz

    integer :: k, n, ip, jq, kr
    real, dimension(isd:ied,jsd:jed,nk,0:1)                 :: sine_arg_z
    real, dimension(isd:ied,jsd:jed,nk,0:1,number_bc_modes) :: bc_modes_zn
    real, dimension(isd:ied,jsd:jed,number_bc_modes)        :: norm_n
    real, dimension(isd:ied,jsd:jed,0:1,number_bc_modes)    :: upsilon_xn, upsilon_yn

    ! compute argument of sine used in defining the baroclinic modes
    sine_arg_z(COMP,:,:)  = 0.0
    sine_arg_z(COMP,nk,1) = 0.5*Ndz_z(COMP,nk,1)
    sine_arg_z(COMP,nk,0) = sine_arg_z(COMP,nk,1) + 0.5*(Ndz_z(COMP,nk,1) + Ndz_z(COMP,nk,0))
    do k = nk-1, 1, -1
       sine_arg_z(COMP,k,1) = sine_arg_z(COMP,k+1,0) + 0.5*(Ndz_z(COMP,k+1,0) + Ndz_z(COMP,k,1))
       sine_arg_z(COMP,k,0) = sine_arg_z(COMP,k,1) + 0.5*(Ndz_z(COMP,k,1) + Ndz_z(COMP,k,0))
    enddo

    call diagnose_3d(Time, Grid, id_sine_arg_pre_up, sine_arg_z(:,:,:,0))
    call diagnose_3d(Time, Grid, id_sine_arg_pre_lo, sine_arg_z(:,:,:,1))

    do k = 1, nk
       do kr = 0, 1
          sine_arg_z(COMP,k,kr) = sine_arg_z(COMP,k,kr)/(epsln + gravity_wave_speed(COMP))
       enddo
    enddo

    ! compute baroclinic modes
    bc_modes_zn(COMP,:,:,:) = 0.0
    do n = 1, number_bc_modes
       ! dimensionless baroclinic mod
       bc_modes_zn(COMP,:,:,n) = sin(sine_arg_z(COMP,:,:)*n)
       ! normalize the baroclinic mode
       norm_n(COMP,n) = sum(sum(N2dz_z(COMP,:,:)*bc_modes_zn(COMP,:,:,n)**2, dim=4), dim=3)
    enddo
    norm_n(COMP,:) = sqrt(grav/(epsln+norm_n(COMP,:)))

    call diagnose_2d(Time, Grid, id_norm_1, norm_n(:,:,1))
    call diagnose_3d(Time, Grid, id_sine_arg_up, sine_arg_z(:,:,:,0))
    call diagnose_3d(Time, Grid, id_sine_arg_lo, sine_arg_z(:,:,:,1))
    call diagnose_3d(Time, Grid, id_bc_mode_1, 0.5*sum(bc_modes_zn(:,:,:,:,1), dim=4))

    do k = 1, nk
       do kr = 0,1
          bc_modes_zn(COMP,k,kr,:) = bc_modes_zn(COMP,k,kr,:)*norm_n(COMP,:)
       enddo
    enddo

    do n = 1, number_bc_modes
       ! projection of agm*grad_rho onto baroclinic modes,
       ! integrated over the depth of an ocean column.
       ! upsilon_xn and upsilon_yn have dimensions kg/(m*sec)
       do ip = 0, 1
          jq = ip
          upsilon_xn(COMP,ip,n) = -sum(sum(bc_modes_zn(COMP,:,:,n)*N2_Sx_agm_dz_xz(COMP,:,ip,:), dim=4), dim=3)/grav
          upsilon_yn(COMP,jq,n) = -sum(sum(bc_modes_zn(COMP,:,:,n)*N2_Sy_agm_dz_yz(COMP,:,jq,:), dim=4), dim=3)/grav
       enddo
    enddo

    do k = 1, nk
       do ip = 0, 1
          jq = ip
          do kr = 0, 1
             upsilonx_xz(COMP,k,ip,kr) = sum(bc_modes_zn(COMP,k,kr,:)*upsilon_xn(COMP,ip,:), dim=3)
             upsilony_yz(COMP,k,jq,kr) = sum(bc_modes_zn(COMP,k,kr,:)*upsilon_yn(COMP,jq,:), dim=3)
          enddo
       enddo
    enddo

    if (regularize_transport) then
       call do_regularize_transport(agm_array, upsilonx_xz, upsilony_yz)
    endif

  end subroutine compute_transport_modes
  ! </SUBROUTINE> NAME="compute_transport_modes"

  !#######################################################################
  ! <SUBROUTINE NAME="do_regularize_transport">
  !
  ! <DESCRIPTION>
  ! regularize the transport to keep magnitude
  ! under control in regions of weak vertical stratification.
  ! </DESCRIPTION>
  !
  subroutine do_regularize_transport(agm_array, upsilonx_xz, upsilony_yz)
    real, dimension(isd:,jsd:,:),       intent(in) :: agm_array
    real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: upsilonx_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: upsilony_yz

    integer :: k, ip, jq, kr
    real, dimension(isd:ied,jsd:jed) :: max_agm
    real, dimension(isd:ied,jsd:jed,0:1) :: col_max_x, col_max_y
    real, dimension(isd:ied,jsd:jed,0:1) :: reg_psix_y, reg_psiy_x

    ! agm_array is ot z-dependent, so the maximum value over the column
    ! is equal to its value at the surface
    max_agm(COMP) = agm_array(COMP,1)*smax_psi
    col_max_x(COMP,:) = maxval(maxval(abs(upsilonx_xz(COMP,:,:,:)), dim=5), dim=3)
    col_max_y(COMP,:) = maxval(maxval(abs(upsilony_yz(COMP,:,:,:)), dim=5), dim=3)

    reg_psix_y(:,:,:) = 1.0
    reg_psiy_x(:,:,:) = 1.0
    do ip=0,1
       jq=ip
       where (col_max_y(COMP,jq) > max_agm(COMP))
          reg_psix_y(COMP,jq) = max_agm(COMP)/col_max_y(COMP,jq)
       endwhere
       where (col_max_x(COMP,ip) > max_agm(COMP))
          reg_psiy_x(COMP,ip) = max_agm(COMP)/col_max_x(COMP,ip)
       endwhere
    enddo
    do k = 1, nk
       do kr = 0, 1
          upsilonx_xz(COMP,k,:,kr) = upsilonx_xz(COMP,k,:,kr)*reg_psiy_x(COMP,:)
          upsilony_yz(COMP,k,:,kr) = upsilony_yz(COMP,k,:,kr)*reg_psix_y(COMP,:)
       enddo
    enddo

  end subroutine do_regularize_transport
  ! </SUBROUTINE> NAME="do_regularize_transport"


  !#######################################################################
  ! <SUBROUTINE NAME="compute_transport_bvp">
  !
  ! <DESCRIPTION>
  ! Compute transport by solving a boundary value problem.
  !
  ! psi is centered on bottom of tracer cell; for example,
  ! psi(k=1)=psi at bottom of tracer cell k=1.
  ! psi vanishes at the ocean surface: psi(k=0)=0
  ! and ocean bottom: psi(k=kmt)=0.
  !
  ! We solve for psi(k=1,kmt-1) using a tridiagonal solver from
  ! Section 2.4 of Press etal 1986.
  !
  ! Units of psi are m^2/sec
  !
  ! </DESCRIPTION>
  !
  subroutine compute_transport_bvp(tmask, Thickness, gravity_wave_speed, &
       N2dz_z, N2_Sx_agm_dz_xz, N2_Sy_agm_dz_yz, upsilonx_xz, upsilony_yz)

    real, dimension(isd:,jsd:,:),       intent(in) :: tmask
    type(ocean_thickness_type),         intent(in) :: Thickness
    real, dimension(isd:,jsd:),         intent(in) :: gravity_wave_speed
    real, dimension(isd:,jsd:,:,0:),    intent(in) :: N2dz_z
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: N2_Sx_agm_dz_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: N2_Sy_agm_dz_yz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: upsilonx_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: upsilony_yz

    integer :: k
    real, dimension(isd:ied,jsd:jed,nk,0:1) :: rhs_dz_ztb_x, rhs_dz_ztb_y
    real, dimension(isd:ied,jsd:jed,nk) :: c2_dzr_a, n2_dz_b, c2_dzr_c, N2dz_ztb
    real, dimension(isd:ied,jsd:jed,0:nk,0:1) :: upsilonx_ztb_x, upsilony_ztb_y

    ! squared wave speed (m/s) for chosen baroclinic mode
    real, dimension(isd:ied,jsd:jed) :: bc_speed2 

    ! compute RHS source terms
    do k = 1, nk-1
       rhs_dz_ztb_x(COMP,k,:) = N2_Sx_agm_dz_xz(COMP,k,:,1) + N2_Sx_agm_dz_xz(COMP,k+1,:,0)
       rhs_dz_ztb_y(COMP,k,:) = N2_Sy_agm_dz_yz(COMP,k,:,1) + N2_Sy_agm_dz_yz(COMP,k+1,:,0)
    enddo

    ! compute c^2
    if (bvp_constant_speed) then
       bc_speed2(COMP) = bvp_speed**2
    else
       bc_speed2(COMP) = max(bvp_min_speed, gravity_wave_speed(COMP)/bvp_bc_mode)**2
    endif

    ! a,b,c fields (using Press etal nomenclature) to invert tridiagonal matrix
    c2_dzr_a(COMP,:) = 0.0 ! a
    n2_dz_b(COMP,:)  = 0.0 ! b
    c2_dzr_c(COMP,:) = 0.0 ! c
    do k = 1, nk-1
       c2_dzr_a(COMP,k) = bc_speed2(COMP)/(epsln+Thickness%dzt(COMP,k))
       c2_dzr_c(COMP,k) = bc_speed2(COMP)/(epsln+Thickness%dzt(COMP,k+1))

       N2dz_ztb(COMP,k) = (N2dz_z(COMP,k,1) + N2dz_z(COMP,k+1,0))
       n2_dz_b(COMP,k)  = -(N2dz_ztb(COMP,k) + c2_dzr_a(COMP,k) + c2_dzr_c(COMP,k))
    enddo

    ! Mask things out
    do k = 1, nk-1
       c2_dzr_a(COMP,k)       = tmask(COMP,k+1)*c2_dzr_a(COMP,k)
       n2_dz_b(COMP,k)        = tmask(COMP,k+1)*n2_dz_b(COMP,k)
       c2_dzr_c(COMP,k)       = tmask(COMP,k+1)*c2_dzr_c(COMP,k)
       rhs_dz_ztb_x(COMP,k,0) = tmask(COMP,k+1)*rhs_dz_ztb_x(COMP,k,0)
       rhs_dz_ztb_x(COMP,k,1) = tmask(COMP,k+1)*rhs_dz_ztb_x(COMP,k,1)
       rhs_dz_ztb_y(COMP,k,0) = tmask(COMP,k+1)*rhs_dz_ztb_y(COMP,k,0)
       rhs_dz_ztb_y(COMP,k,1) = tmask(COMP,k+1)*rhs_dz_ztb_y(COMP,k,1)
    enddo

    ! invert to solve for upsilonx, upsilony...
    call invtri_bvp(c2_dzr_a, n2_dz_b, c2_dzr_c, rhs_dz_ztb_x, &
                    rhs_dz_ztb_y, upsilonx_ztb_x, upsilony_ztb_y)

    do k = 1, nk
       upsilonx_xz(COMP,k,:,0) = upsilonx_ztb_x(COMP,k-1,:)
       upsilonx_xz(COMP,k,:,1) = upsilonx_ztb_x(COMP,k,:)
       upsilony_yz(COMP,k,:,0) = upsilony_ztb_y(COMP,k-1,:)
       upsilony_yz(COMP,k,:,1) = upsilony_ztb_y(COMP,k,:)
    enddo

  end subroutine compute_transport_bvp
  ! </SUBROUTINE> NAME="compute_transport_bvp"


  !#######################################################################
  ! <SUBROUTINE NAME="invtri_bvp">
  !
  ! <DESCRIPTION>
  ! Solve the vertical diffusion equation implicitly using the
  ! method of inverting a tridiagonal matrix as described in
  ! Numerical Recipes in Fortran, The art of Scientific Computing,
  ! Second Edition, Press, Teukolsky, Vetterling, Flannery, 1992
  ! pages 42,43.
  !
  ! enforce upsilon(k=kmt) = 0 via use of mask(k+1).
  !
  ! </DESCRIPTION>
  !
  subroutine invtri_bvp (a, b, c, r_x, r_y, upsilon_x, upsilon_y)

    real, dimension(isd:,jsd:,:),     intent(in) :: a
    real, dimension(isd:,jsd:,:),     intent(in) :: b
    real, dimension(isd:,jsd:,:),     intent(in) :: c
    real, dimension(isd:,jsd:,:,0:),  intent(in) :: r_x
    real, dimension(isd:,jsd:,:,0:),  intent(in) :: r_y
    real, dimension(isd:,jsd:,0:,0:), intent(out) :: upsilon_x
    real, dimension(isd:,jsd:,0:,0:), intent(out) :: upsilon_y

    real, dimension(isd:ied,jsd:jed,nk) :: bet
    real, dimension(isd:ied,jsd:jed,nk) :: gam

    integer :: k, ip, jq

    bet = 0.0
    gam = 0.0

    ! decomposition and forward substitution
    bet(COMP,1) = 1.0/(b(COMP,1) + epsln)
    do k = 2, nk-1
       gam(COMP,k) = c(COMP,k-1)*bet(COMP,k-1)
       bet(COMP,k) = 1.0/(b(COMP,k) - a(COMP,k)*gam(COMP,k) + epsln)
    enddo

    upsilon_x(COMP,0,:) = 0.0
    upsilon_y(COMP,0,:) = 0.0
    do ip=0,1
       jq=ip
       upsilon_x(COMP,1,ip) = r_x(COMP,1,ip)*bet(COMP,1)
       upsilon_y(COMP,1,jq) = r_y(COMP,1,jq)*bet(COMP,1)
       do k = 2, nk-1
          upsilon_x(COMP,k,ip) = (r_x(COMP,k,ip) - a(COMP,k)*upsilon_x(COMP,k-1,ip))*bet(COMP,k)
          upsilon_y(COMP,k,jq) = (r_y(COMP,k,jq) - a(COMP,k)*upsilon_y(COMP,k-1,jq))*bet(COMP,k)
       enddo
       ! back substitution
       do k = nk-2, 1, -1
          upsilon_x(COMP,k,ip) = upsilon_x(COMP,k,ip) - gam(COMP,k+1)*upsilon_x(COMP,k+1,ip)
          upsilon_y(COMP,k,jq) = upsilon_y(COMP,k,jq) - gam(COMP,k+1)*upsilon_y(COMP,k+1,jq)
       enddo
    enddo
    upsilon_x(COMP,nk,:) = 0.0
    upsilon_y(COMP,nk,:) = 0.0

  end subroutine invtri_bvp
  ! </SUBROUTINE> NAME="invtri_bvp">

  !#######################################################################
  ! <SUBROUTINE NAME="do_smooth_transport">
  !
  ! <DESCRIPTION>
  ! </DESCRIPTION>
  !
  subroutine do_smooth_transport(Domain, tmask, upsilonx_xz, upsilony_yz)

    type(ocean_domain_type),            intent(inout) :: Domain
    real, dimension(isd:,jsd:,:),       intent(in)    :: tmask
    real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: upsilonx_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: upsilony_yz

    integer :: kr, ip, jq

    ! smooth transport to reduce potentials for checkerboard noise
    do ip=0,1
       jq = ip
       do kr = 0, 1
          call horz_smooth(Domain, tmask, upsilonx_xz(:,:,:,ip,kr))
          call horz_smooth(Domain, tmask, upsilony_yz(:,:,:,jq,kr))
          upsilonx_xz(COMP,:,ip,kr) = upsilonx_xz(COMP,:,ip,kr)*tmask(COMP,:)*tmask(COMPXRR,:)
          upsilony_yz(COMP,:,jq,kr) = upsilony_yz(COMP,:,jq,kr)*tmask(COMP,:)*tmask(COMPYRR,:)
       enddo
    enddo

  end subroutine do_smooth_transport
  ! </SUBROUTINE> NAME="do_smooth_transport"


end module ocean_nphysics_skew_mod

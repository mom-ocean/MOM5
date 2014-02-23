module ocean_nphysics_tensor_mod
  !
  !<CONTACT EMAIL="t.leslie@unsw.edu.au"> Tim Leslie
  !</CONTACT>
  !
  !<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
  !</CONTACT>
  ! 
  !<OVERVIEW>
  !</OVERVIEW>
  !
  !<DESCRIPTION>
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
  !<NAMELIST NAME="ocean_nphysics_tensor_nml">
  !  <DATA NAME="" TYPE="logical">
  !  </DATA>
  !
  !</NAMELIST>

#define COMP isc:iec,jsc:jec

  use constants_mod, only: epsln, pi
  use fms_mod,       only: open_namelist_file, check_nml_error, close_file
  use mpp_mod,       only: stdout, stdlog, input_nml_file 

  use ocean_domains_mod,           only: get_local_indices
  use ocean_nphysics_skew_mod,     only: skew_transport
  use ocean_sigma_transport_mod,   only: tmask_sigma_on, tmask_sigma
  use ocean_types_mod,             only: ocean_grid_type, ocean_domain_type
  use ocean_types_mod,             only: ocean_time_type, ocean_thickness_type
  use ocean_util_mod,              only: register_3d_t_field, diagnose_3d

  implicit none

  public ocean_nphysics_tensor_init
  public compute_tensors

  private sine_taper
  private compute_diffusion_tapers
  private compute_masked_tensor

  private

#include <ocean_memory.h>

character(len=128) :: version=&
     '$Id: ocean_nphysics_tensor.F90,v 20.0 2013/12/14 00:14:50 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

character(len=*), parameter :: FILENAME=&
     __FILE__


  ! Diagnostics

  ! Tensors
  integer :: id_symm_tensor1122
  integer :: id_symm_tensor33x
  integer :: id_symm_tensor33y
  integer :: id_symm_tensor13
  integer :: id_symm_tensor23
  integer :: id_upsilonx
  integer :: id_upsilony

  ! Tapers
  integer :: id_slope_taper
  integer :: id_depth_taper
  integer :: id_neutral_taper
  integer :: id_displacement_depth

  ! Namelist Values

  ! for tapering neutral physics over penetration depth of
  ! eddies as determined by slope and Rossby radius
  logical :: neutral_sine_taper = .true.
  real :: turb_blayer_min   = 0.0 ! metres (was always set to zero in mom4p0)
  real :: rossby_radius_max = 100e3  ! metres
  real :: rossby_radius_min = 15e3   ! metres

  ! for setting the slope tapering methods
  logical :: dm_taper = .true. ! true -> tanh tapering scheme of Danabasoglu and McWilliams
                               ! false -> quadratic tapering of Gerdes, Koberle, and Willebrand
  ! slope width over which fluxes tapered using tanh function using dm_taper scheme
  real :: swidth = 0.05*.01    ! [dimensionless]

  ! to remove diagonal elements to neutral diffusion within the neutral boundar layer
  logical :: neutral_taper_diagonal = .false.

  ! to reduce neutral fluxes to horz/vert diffusion next to model boundaries
  logical :: tmask_neutral_on = .false.

  namelist /ocean_nphysics_tensor_nml/ neutral_sine_taper, turb_blayer_min, &
       rossby_radius_max, rossby_radius_min, dm_taper, swidth,              &
       neutral_taper_diagonal, tmask_neutral_on

contains

  !#######################################################################
  ! <SUBROUTINE NAME="ocean_nphysics_tensor_init">
  !
  ! <DESCRIPTION>
  ! Initialise this module.
  ! </DESCRIPTION>
  !
  subroutine ocean_nphysics_tensor_init(Domain, Grid, Time)

    type(ocean_domain_type), intent(in) :: Domain
    type(ocean_grid_type),   intent(in) :: Grid
    type(ocean_time_type),   intent(in) :: Time

    integer :: ierr, ioun, io_status
    integer :: stdoutunit, stdlogunit
    stdoutunit = stdout()
    stdlogunit = stdlog()

    ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_nphysics_tensor_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_nphysics_tensor_nml')
#else
    ioun =  open_namelist_file()
    read (ioun,ocean_nphysics_tensor_nml,IOSTAT=io_status)
    ierr = check_nml_error(io_status,'ocean_nphysics_tensor_nml')
    call close_file (ioun)
#endif
    write (stdoutunit,'(/)')
    write (stdoutunit,ocean_nphysics_tensor_nml)  
    write (stdlogunit,ocean_nphysics_tensor_nml)

#ifndef MOM_STATIC_ARRAYS
    call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
    nk = Grid%nk
#endif

    call register_fields(Grid, Time)
   
  end subroutine ocean_nphysics_tensor_init
  ! </SUBROUTINE> NAME="ocean_nphysics_tensor_init"

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

    id_symm_tensor1122 = register_3d_t_field(Grid, Time, 'symm_tensor1122', &
         'Horizontal diagonal diffusion tensor component', 'm^2/s', (/-1e8,1.e8/))
    id_symm_tensor33x = register_3d_t_field(Grid, Time, 'symm_tensor33x',      &
         'x-component of 33 term in neutral diffusion tensor', 'm^2/s', range=(/-1e8,1.e8/))
    id_symm_tensor33y = register_3d_t_field(Grid, Time, 'symm_tensor33y',      &
         'y-component of 33 term in neutral diffusion tensor', 'm^2/s', range=(/-1e8,1.e8/))
    id_symm_tensor13 = register_3d_t_field(Grid, Time, 'symm_tensor13',      &
         '13 term in neutral diffusion tensor', 'm^2/s', range=(/-1e8,1.e8/))
    id_symm_tensor23 = register_3d_t_field(Grid, Time, 'symm_tensor23',      &
         '23 term in neutral diffusion tensor', 'm^2/s', range=(/-1e8,1.e8/))
    id_upsilonx = register_3d_t_field(Grid, Time, 'upsilonx',      &
         'x component of GM transport vector', 'm^2/s', range=(/-1e8,1.e8/))
    id_upsilony = register_3d_t_field(Grid, Time, 'upsilony',      &
         'y component of GM transport vector', 'm^2/s', range=(/-1e8,1.e8/))

    id_slope_taper = register_3d_t_field(Grid, Time, 'slope_taper',      &
         'Taper function for handling steep neutral slope', '-', range=(/0.0, 1.0/))
    id_depth_taper = register_3d_t_field(Grid, Time, 'depth_taper',      &
         'Taper function for handling surface region', '-', range=(/0.0, 1.0/))
    id_neutral_taper = register_3d_t_field(Grid, Time, 'neutral_taper',      &
         'Taper function product of slope_taper and depth_taper', '-', range=(/0.0, 1.0/))
    id_displacement_depth = register_3d_t_field(Grid, Time, 'displacement_depth',      &
         'Effective displacement depth due to undulating density surfaces', '-', range=(/0.0, 1.e5/))

  end subroutine register_fields
  ! </SUBROUTINE> NAME="register_fields"

  !#######################################################################
  ! <SUBROUTINE NAME="compute_tensors">
  !
  ! <DESCRIPTION>
  ! Compute the tensor components for the diffusion and skew-diffusion tensor.
  ! </DESCRIPTION>
  !
  subroutine compute_tensors(Time, Domain, Grid, Thickness,                    &
       slopex_xz, slopey_yz, Ndz_z, N2_z, gravity_wave_speed,                  &
       rossby_radius, absslope, absslope_z,                                    &
       aredi_array, ah_array, agm_array, ksurf_blayer, smax, surf_blthick,     &
       symm_tensor13_xz, symm_tensor23_yz, symm_tensor33_xz, symm_tensor33_yz, &
       symm_tensor1122_z, upsilonx_xz, upsilony_yz)

    type(ocean_time_type),      intent(in)    :: Time
    type(ocean_domain_type),    intent(inout) :: Domain
    type(ocean_grid_type),      intent(in)    :: Grid
    type(ocean_thickness_type), intent(in)    :: Thickness

    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: slopex_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: slopey_yz
    real, dimension(isd:,jsd:,:,0:),    intent(in) :: Ndz_z
    real, dimension(isd:,jsd:,:,0:),    intent(in) :: N2_z
    real, dimension(isd:,jsd:),         intent(in) :: gravity_wave_speed
    real, dimension(isd:,jsd:),         intent(in) :: rossby_radius
    real, dimension(isd:,jsd:,:),       intent(in) :: absslope
    real, dimension(isd:,jsd:,:,0:),    intent(in) :: absslope_z

    real, dimension(isd:,jsd:,:), intent(in) :: aredi_array
    real, dimension(isd:,jsd:),   intent(in) :: ah_array
    real, dimension(isd:,jsd:,:), intent(in) :: agm_array

    integer, dimension(isd:,jsd:), intent(in) :: ksurf_blayer
    real,                          intent(in) :: smax
    real, dimension(isd:,jsd:),    intent(in) :: surf_blthick

    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: symm_tensor13_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: symm_tensor23_yz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: symm_tensor33_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: symm_tensor33_yz
    real, dimension(isd:,jsd:,:,0:),    intent(out) :: symm_tensor1122_z
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: upsilonx_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: upsilony_yz

    ! Taper function based on depth, values in [0..1]
    real, dimension(isd:ied,jsd:jed,nk) :: depth_taper

    ! Taper over the eddy depth
    call sine_taper(Time, Grid, Thickness, absslope, surf_blthick, rossby_radius, depth_taper)

    call diffusion_tensor(Time, Grid, slopex_xz, slopey_yz, absslope_z, &
         aredi_array, ah_array, smax, depth_taper,                      &
         symm_tensor1122_z, symm_tensor13_xz, symm_tensor23_yz,         &
         symm_tensor33_xz, symm_tensor33_yz)

    call skew_transport(Domain, Time, Grid, Thickness,                     &
         agm_array, slopex_xz, slopey_yz, Ndz_z, N2_z, gravity_wave_speed, &
         depth_taper, ksurf_blayer, upsilonx_xz, upsilony_yz)

    call compute_masked_tensor(Grid%kmt, aredi_array, symm_tensor1122_z, &
         symm_tensor33_xz, symm_tensor33_yz, symm_tensor13_xz,           &
         symm_tensor23_yz, upsilonx_xz, upsilony_yz)

    call diagnose_3d(Time, Grid, id_symm_tensor1122, 0.5*sum(symm_tensor1122_z, dim=4))
    call diagnose_3d(Time, Grid, id_symm_tensor33x, 0.25*sum(sum(symm_tensor33_xz, dim=5), dim=4))
    call diagnose_3d(Time, Grid, id_symm_tensor33y, 0.25*sum(sum(symm_tensor33_yz, dim=5), dim=4))
    call diagnose_3d(Time, Grid, id_symm_tensor13, 0.25*sum(sum(symm_tensor13_xz, dim=5), dim=4))
    call diagnose_3d(Time, Grid, id_symm_tensor23, 0.25*sum(sum(symm_tensor23_yz, dim=5), dim=4))
    call diagnose_3d(Time, Grid, id_upsilonx, 0.25*sum(sum(upsilonx_xz, dim=5), dim=4))
    call diagnose_3d(Time, Grid, id_upsilony, 0.25*sum(sum(upsilony_yz, dim=5), dim=4))

  end subroutine compute_tensors
  ! </SUBROUTINE> NAME="compute_tensors"

  !#######################################################################
  ! <SUBROUTINE NAME="sine_taper">
  !
  ! <DESCRIPTION>
  ! Calculate a sine taper for those points shallower than the given threshold.
  !
  ! The depth_function gives the depth (in metres) at each (i,j,k) point.
  ! The threshold_depth is the depth (in metres) at each (i,j) point.
  ! The returned taper gives a value in the range [0, 1] at each (i,j,k) point,
  ! where those points at a depth below the threshold are set to 1.0 and those
  ! above approach 0.0 as they approach the surface according to a sine
  ! function.
  !
  ! </DESCRIPTION>
  !
  subroutine sine_taper(Time, Grid, Thickness, absslope, surf_blthick, rossby_radius, taper)
    type(ocean_time_type),        intent(in) :: Time
    type(ocean_grid_type),        intent(in) :: Grid
    type(ocean_thickness_type),   intent(in) :: Thickness
    real, dimension(isd:,jsd:,:), intent(in) :: absslope
    real, dimension(isd:,jsd:),   intent(in) :: surf_blthick
    real, dimension(isd:,jsd:),   intent(in) :: rossby_radius
    real, dimension(isd:,jsd:,:), intent(out) :: taper

    real, dimension(isd:ied,jsd:jed) :: capped_rossby_radius  ! [m]
    real, dimension(isd:ied,jsd:jed,nk) :: displacement_depth ! [m]
    real, dimension(isd:ied,jsd:jed,nk) :: depth_ratio        ! [1]

    integer :: k

    displacement_depth(:,:,:) = 0.0
    if (neutral_sine_taper) then
       ! EoM: 16.68
       capped_rossby_radius(COMP) = min(rossby_radius_max, max(rossby_radius_min, rossby_radius(COMP))) 
       do k = 1, nk
          ! FOCM 15.22
          displacement_depth(COMP,k) = max(turb_blayer_min, surf_blthick(COMP), capped_rossby_radius(COMP)*absslope(COMP,k))
       enddo

       where (Thickness%depth_zt(COMP,:) <= displacement_depth(COMP,:))
          depth_ratio(COMP,:) = Thickness%depth_zt(COMP,:)/(epsln+displacement_depth(COMP,:))
          taper(COMP,:) = 0.5*(1.0 + sin(pi*(depth_ratio(COMP,:) - 0.5))) ! FOCM 15.24
       elsewhere
          taper(COMP,:) = 1.0
       endwhere
       call diagnose_3d(Time, Grid, id_displacement_depth, displacement_depth, abs_max=1.5e4)
    else
       taper(COMP,:) = 1.0
    endif

    call diagnose_3d(Time, Grid, id_depth_taper, taper)

  end subroutine sine_taper
  ! </SUBROUTINE> NAME="sine_taper"


  !#######################################################################
  ! <SUBROUTINE NAME="diffusion_tensor">
  !
  ! <DESCRIPTION>
  ! Compute the components of the neutral diffusion tensor.
  ! </DESCRIPTION>
  !
  subroutine diffusion_tensor(Time, Grid, slopex_xz, slopey_yz, absslope_z,           &
       aredi_array, ah_array, smax, depth_taper, symm_tensor1122_z, symm_tensor13_xz, &
       symm_tensor23_yz, symm_tensor33_xz, symm_tensor33_yz)

    type(ocean_time_type),              intent(in) :: Time
    type(ocean_grid_type),              intent(in) :: Grid
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: slopex_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: slopey_yz
    real, dimension(isd:,jsd:,:,0:),    intent(in) :: absslope_z
    real, dimension(isd:,jsd:,:),       intent(in) :: aredi_array
    real, dimension(isd:,jsd:),         intent(in) :: ah_array
    real,                               intent(in) :: smax
    real, dimension(isd:,jsd:,:),       intent(in) :: depth_taper

    real, dimension(isd:,jsd:,:,0:),    intent(out) :: symm_tensor1122_z
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: symm_tensor13_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: symm_tensor23_yz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: symm_tensor33_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: symm_tensor33_yz

    ! Taper functions, values in [0..1]
    real, dimension(isd:ied,jsd:jed,nk,0:1) :: neutral_taper_z
    real, dimension(isd:ied,jsd:jed,nk,0:1) :: diagonal_taper_z
    real, dimension(isd:ied,jsd:jed,nk)     :: horiz_diff_taper

    ! Tapered Redi diffusivity [m^2/s]
    real, dimension(isd:ied,jsd:jed,nk,0:1)     :: tapered_aredi_z
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: full_tapered_aredi

    ! Tapered ah diffusivity [m^2/s]
    real, dimension(isd:ied,jsd:jed,nk) :: ah_term

    call compute_diffusion_tapers(Time, Grid, absslope_z, smax, depth_taper, &
         neutral_taper_z, diagonal_taper_z, horiz_diff_taper)
   
    tapered_aredi_z(COMP,:,0)      = aredi_array(COMP,:)*neutral_taper_z(COMP,:,0)
    tapered_aredi_z(COMP,:,1)      = aredi_array(COMP,:)*neutral_taper_z(COMP,:,1)
    full_tapered_aredi(COMP,:,:,:) = spread(tapered_aredi_z(COMP,:,:), 5, 2)

    ! Off diagonal terms
    symm_tensor13_xz(COMP,:,:,:) = full_tapered_aredi(COMP,:,:,:)*slopex_xz(COMP,:,:,:)
    symm_tensor23_yz(COMP,:,:,:) = full_tapered_aredi(COMP,:,:,:)*slopey_yz(COMP,:,:,:)

    ! Diagonal terms
    ah_term(COMP,:)             = spread(ah_array(COMP), 3, nk)*horiz_diff_taper(COMP,:)
    symm_tensor1122_z(COMP,:,0) = aredi_array(COMP,:)*diagonal_taper_z(COMP,:,0) + ah_term(COMP,:)
    symm_tensor1122_z(COMP,:,1) = aredi_array(COMP,:)*diagonal_taper_z(COMP,:,1) + ah_term(COMP,:)

    symm_tensor33_xz(COMP,:,:,:) = full_tapered_aredi(COMP,:,:,:)*slopex_xz(COMP,:,:,:)**2
    symm_tensor33_yz(COMP,:,:,:) = full_tapered_aredi(COMP,:,:,:)*slopey_yz(COMP,:,:,:)**2

  end subroutine diffusion_tensor
  ! </SUBROUTINE> NAME="diffusion_tensor"

  !#######################################################################
  ! <SUBROUTINE NAME="compute_diffusion_tapers">
  !
  ! <DESCRIPTION>
  ! Compute the taper functions to be applied to the components of the diffusion
  ! tensor.
  ! </DESCRIPTION>
  !
  subroutine compute_diffusion_tapers(Time, Grid, absslope_z, smax, depth_taper, &
       neutral_taper_z, diagonal_taper_z, horiz_diff_taper)
    
    type(ocean_time_type),           intent(in) :: Time
    type(ocean_grid_type),           intent(in) :: Grid
    real, dimension(isd:,jsd:,:,0:), intent(in) :: absslope_z
    real,                            intent(in) :: smax
    real, dimension(isd:,jsd:,:),    intent(in) :: depth_taper

    real, dimension(isd:,jsd:,:,0:), intent(out) :: neutral_taper_z
    real, dimension(isd:,jsd:,:,0:), intent(out) :: diagonal_taper_z
    real, dimension(isd:,jsd:,:),    intent(out) :: horiz_diff_taper

    ! Taper function, values in [0..1]
    real, dimension(isd:ied,jsd:jed,nk,0:1) :: slope_taper_z

    ! Taper for steep slope regions
    if (dm_taper) then
       slope_taper_z(COMP,:,:) = 0.5*(1.0 + tanh((smax-absslope_z(COMP,:,:))/(swidth+epsln)))
    else
       where (absslope_z(COMP,:,:) > smax)
          slope_taper_z(COMP,:,:) = (smax/absslope_z(COMP,:,:))**2
       elsewhere
          slope_taper_z(COMP,:,:) = 1.0
       endwhere
    endif

    neutral_taper_z(COMP,:,0) = slope_taper_z(COMP,:,0)*depth_taper(COMP,:)
    neutral_taper_z(COMP,:,1) = slope_taper_z(COMP,:,1)*depth_taper(COMP,:)

    if (neutral_taper_diagonal) then
       diagonal_taper_z(COMP,:,:) = neutral_taper_z(COMP,:,:)
    else
       diagonal_taper_z(COMP,:,:) = 1.0
    endif

    horiz_diff_taper(COMP,:) = 1.0 - depth_taper(COMP,:)

    call diagnose_3d(Time, Grid, id_slope_taper, 0.5*sum(slope_taper_z(:,:,:,:), dim=4))
    call diagnose_3d(Time, Grid, id_neutral_taper, 0.5*sum(neutral_taper_z(:,:,:,:), dim=4))

  end subroutine compute_diffusion_tapers
  ! </SUBROUTINE> NAME="compute_diffusion_tapers"


  !#######################################################################
  ! <SUBROUTINE NAME="compute_masked_tensor">
  !
  ! <DESCRIPTION>
  ! This function applies the boundary conditions relevent to the flags
  ! tmask_neutral_on and tmask_sigma_on.
  ! </DESCRIPTION>
  !
  subroutine compute_masked_tensor(kmt, aredi_array,          &
       symm_tensor1122_z, symm_tensor33_xz, symm_tensor33_yz, &
       symm_tensor13_xz, symm_tensor23_yz, upsilonx_xz, upsilony_yz)
    
    integer, dimension(isd:,jsd:),      intent(in)    :: kmt
    real, dimension(isd:,jsd:,:),       intent(in)    :: aredi_array
    real, dimension(isd:,jsd:,:,0:),    intent(inout) :: symm_tensor1122_z
    real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: symm_tensor33_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: symm_tensor33_yz

    real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: symm_tensor13_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: symm_tensor23_yz
    real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: upsilonx_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: upsilony_yz

    integer :: k, ip, kr

    ! Apply masks to use horizonal diffusion at boundaries
    if (tmask_neutral_on) then
       ! Surface
       symm_tensor1122_z(COMP,1,0) = aredi_array(COMP,1)
       symm_tensor1122_z(COMP,1,1) = aredi_array(COMP,1)
       symm_tensor33_xz(COMP,1,:,:) = 0.0
       symm_tensor33_yz(COMP,1,:,:) = 0.0
       symm_tensor13_xz(COMP,1,:,:) = 0.0
       symm_tensor23_yz(COMP,1,:,:) = 0.0
       upsilonx_xz(COMP,1,:,:) = 0.0
       upsilony_yz(COMP,1,:,:) = 0.0

       ! Bottom
       do k = 2, nk
          where (k == kmt(COMP))
             symm_tensor1122_z(COMP,k,0) = aredi_array(COMP,k)
             symm_tensor1122_z(COMP,k,1) = aredi_array(COMP,k)
          endwhere
          do ip = 0, 1
             do kr = 0, 1
                where (k == kmt(COMP))
                   symm_tensor33_xz(COMP,k,ip,kr) = 0.0
                   symm_tensor33_yz(COMP,k,ip,kr) = 0.0
                   symm_tensor13_xz(COMP,k,ip,kr) = 0.0
                   symm_tensor23_yz(COMP,k,ip,kr) = 0.0
                   upsilonx_xz(COMP,k,ip,kr) = 0.0
                   upsilony_yz(COMP,k,ip,kr) = 0.0
                endwhere
             enddo
          enddo
       enddo
    endif

    if (tmask_sigma_on) then
       do k = 2, nk
          do ip = 0, 1
             do kr = 0, 1
                where (k == kmt(COMP) .and. tmask_sigma(COMP) == 1.0) 
                   symm_tensor1122_z(COMP,k,kr) = 0.0
                   symm_tensor33_xz(COMP,k,ip,kr) = 0.0
                   symm_tensor33_yz(COMP,k,ip,kr) = 0.0
                   symm_tensor13_xz(COMP,k,ip,kr) = 0.0
                   symm_tensor23_yz(COMP,k,ip,kr) = 0.0
                   upsilonx_xz(COMP,k,ip,kr) = 0.0
                   upsilony_yz(COMP,k,ip,kr) = 0.0
                endwhere
             enddo
          enddo
       enddo
    endif

  end subroutine compute_masked_tensor
  ! </SUBROUTINE> NAME="compute_masked_tensor"

end module ocean_nphysics_tensor_mod

module ocean_nphysics_flux_mod

!<CONTACT EMAIL="tim.leslie@gmail.com"> Tim Leslie
!</CONTACT>
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Compute the neutral physics fluxes. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the neutral physics fluxes. 
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
! S.M. Griffies: Elements of MOM (2012)
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_nphysics_flux_nml">
!
!  <DATA NAME="diffusion_all_explicit" TYPE="logical">
! To compute K33 explicitly in time.  This setting is meant
! only for debugging tests, since in general the simulation 
! will go unstable. 
! Default diffusion_all_explicit=.false.
!  </DATA> 
!
!  <DATA NAME="neutral_physics_limit" TYPE="logical">
! Revert to horizontal diffusion when tracer falls outside specified range.
! Default neutral_physics_limit=.true., so to keep tracers from going
! too far outside of physical range. 
!  </DATA> 
!
!</NAMELIST>

#define COMP isc:iec,jsc:jec
#define COMPXL isc-1:iec,jsc:jec
#define COMPXR isc:iec+1,jsc:jec
#define COMPXLL isc-1:iec-1,jsc:jec
#define COMPYL isc:iec,jsc-1:jec
#define COMPYR isc:iec,jsc:jec+1
#define COMPYLL isc:iec,jsc-1:jec-1
#define COMPXLYL isc-1:iec,jsc-1:jec

  use constants_mod,   only: epsln
  use fms_mod,         only: open_namelist_file, check_nml_error, close_file
  use mpp_domains_mod, only: mpp_update_domains, CGRID_NE
  use mpp_mod,         only: stdout, stdlog, mpp_clock_id, CLOCK_ROUTINE
  use mpp_mod,         only: mpp_clock_begin, mpp_clock_end
  use mpp_mod,         only: input_nml_file

  use ocean_domains_mod,           only: get_local_indices, set_ocean_domain
  use ocean_nphysics_util_new_mod, only: stencil_centre_to_vert, stencil_centre_to_horiz
  use ocean_nphysics_util_new_mod, only: vert_smooth
  use ocean_operators_mod,         only: FMX, FMY
  use ocean_types_mod,             only: ocean_grid_type, ocean_domain_type, ocean_time_type
  use ocean_types_mod,             only: ocean_thickness_type, ocean_prog_tracer_type
  use ocean_util_mod,              only: write_note, diagnose_3d 
  use ocean_util_mod,              only: register_3d_t_field, register_3d_xte_field
  use ocean_util_mod,              only: register_3d_ytn_field, register_3d_ztb_field

  implicit none

  public flux_calculations
  public ocean_nphysics_flux_init

  private compute_mass_diff
  private geometric_terms
  private compute_33_term
  private compute_fluxes
  private apply_tracer_limits
  private update_tendencies

  private

#include <ocean_memory.h>

character(len=128) :: version=&
     '$Id: ocean_nphysics_flux.F90,v 20.0 2013/12/14 00:14:44 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

character(len=*), parameter :: FILENAME=&
     __FILE__


  type(ocean_domain_type), save :: Dom_flux

  ! Diagnostics

  integer :: id_clock_flux_mass_diff
  integer :: id_clock_flux_33_term
  integer :: id_clock_flux_fluxes
  integer :: id_clock_flux_tracer_limits
  integer :: id_clock_flux_misc
  integer :: id_clock_flux_update

  ! Mass weighted tensor components
  integer :: id_symm_mass_diff_11_xte
  integer :: id_symm_mass_diff_22_ytn

  integer :: id_symm_mass_diff_13_h_xz
  integer :: id_symm_mass_diff_23_h_yz
  integer :: id_symm_mass_diff_31_v_xz
  integer :: id_symm_mass_diff_32_v_yz
  
  integer :: id_skew_mass_diff_13_h_xz
  integer :: id_skew_mass_diff_23_h_yz
  integer :: id_skew_mass_diff_31_v_xz
  integer :: id_skew_mass_diff_32_v_yz
  
  integer :: id_m33_explicit_ztb
  integer :: id_m33_implicit

  ! Fluxes
  integer, dimension(:), allocatable :: id_diff_flux_x_xte
  integer, dimension(:), allocatable :: id_diff_flux_y_ytn
  integer, dimension(:), allocatable :: id_diff_flux_z_ztb
  integer, dimension(:), allocatable :: id_skew_flux_x_xte
  integer, dimension(:), allocatable :: id_skew_flux_y_ytn
  integer, dimension(:), allocatable :: id_skew_flux_z_ztb

  ! Tendencies
  integer, dimension(:), allocatable :: id_diff_th_tendency
  integer, dimension(:), allocatable :: id_skew_th_tendency

  ! Namelist

  ! to compute K33 explicitly in time.
  ! diffusion_all_explicit=.false. for realistic simulations.
  logical :: diffusion_all_explicit = .false.

  ! revert to horizontal diffusion when tracer falls outside specified range
  logical :: neutral_physics_limit = .true.

  namelist /ocean_nphysics_flux_nml/ diffusion_all_explicit, neutral_physics_limit

contains

  !#######################################################################
  ! <SUBROUTINE NAME="ocean_nphysics_flux_init">
  !
  ! <DESCRIPTION>
  ! Initialise namelist variables and prepare diagnostics.
  ! </DESCRIPTION>
  !
  subroutine ocean_nphysics_flux_init(Domain, Grid, Time, T_prog)

    type(ocean_domain_type), intent(in) :: Domain
    type(ocean_grid_type),   intent(in) :: Grid
    type(ocean_time_type),   intent(in) :: Time
    type(ocean_prog_tracer_type), dimension(:), intent(in) :: T_prog

    integer :: ierr, ioun, io_status
    integer :: stdoutunit, stdlogunit
    stdoutunit = stdout()
    stdlogunit = stdlog()

    ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_nphysics_flux_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_nphysics_flux_nml')
#else
    ioun =  open_namelist_file()
    read (ioun,ocean_nphysics_flux_nml,IOSTAT=io_status)
    ierr = check_nml_error(io_status,'ocean_nphysics_flux_nml')
    call close_file (ioun)
#endif
    write (stdoutunit,'(/)')
    write (stdoutunit,ocean_nphysics_flux_nml)  
    write (stdlogunit,ocean_nphysics_flux_nml)

#ifndef MOM_STATIC_ARRAYS
    call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
    nk = Grid%nk
#endif

    call set_ocean_domain(Dom_flux,Grid,xhalo=Domain%xhalo,yhalo=Domain%yhalo,name='flux dom neutral',maskmap=Domain%maskmap)
    
    id_clock_flux_mass_diff     = mpp_clock_id('(Ocean neutral: flux: mass diff)'     ,grain=CLOCK_ROUTINE)
    id_clock_flux_33_term       = mpp_clock_id('(Ocean neutral: flux: 33 term)'       ,grain=CLOCK_ROUTINE)
    id_clock_flux_fluxes        = mpp_clock_id('(Ocean neutral: flux: fluxes)'        ,grain=CLOCK_ROUTINE)
    id_clock_flux_tracer_limits = mpp_clock_id('(Ocean neutral: flux: tracer limits)' ,grain=CLOCK_ROUTINE)
    id_clock_flux_misc          = mpp_clock_id('(Ocean neutral: flux: misc)'          ,grain=CLOCK_ROUTINE)
    id_clock_flux_update        = mpp_clock_id('(Ocean neutral: flux: update)'        ,grain=CLOCK_ROUTINE)

    call register_fields(Grid, Time, T_prog)

  end subroutine ocean_nphysics_flux_init
  ! </SUBROUTINE> NAME="ocean_nphysics_flux_init"


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

    id_symm_mass_diff_11_xte = register_3d_xte_field(Grid, Time, 'symm_mass_diff_11_xte', &
         '11 component of mass weighted neutral diffusivity at Eastern tracer-cell face', &
         'kg m^2/s', range=(/-1e8,1.e8/))
    id_symm_mass_diff_22_ytn = register_3d_ytn_field(Grid, Time, 'symm_mass_diff_22_ytn',  &
         '22 component of mass weighted neutral diffusivity at Northern tracer-cell face', &
         'kg m^2/s', range=(/-1e8,1.e8/))
    id_symm_mass_diff_13_h_xz = register_3d_xte_field(Grid, Time,  'symm_mass_diff_13_h_xz',             &
         'Triad-averaged 13 component of mass weighted neutral diffusivity at Eastern tracer-cell face', &
         'kg m^2/s', range=(/-1e8,1.e8/))
    id_symm_mass_diff_23_h_yz = register_3d_ytn_field(Grid, Time, 'symm_mass_diff_23_h_yz',               &
         'Triad-averaged 13 component of mass weighted neutral diffusivity at Northern tracer-cell face', &
         'kg m^2/s', range=(/-1e8,1.e8/))
    id_symm_mass_diff_31_v_xz = register_3d_ztb_field(Grid, Time, 'symm_mass_diff_31_v_xz',             &
         'Triad-averaged 31 component of mass weighted neutral diffusivity at bottom tracer-cell face', &
         'kg m^2/s', range=(/-1e8,1.e8/))
    id_symm_mass_diff_32_v_yz = register_3d_ztb_field(Grid, Time, 'symm_mass_diff_32_v_yz',             &
         'Triad-averaged 32 component of mass weighted neutral diffusivity at bottom tracer-cell face', &
         'kg m^2/s', range=(/-1e8,1.e8/))

    id_skew_mass_diff_13_h_xz = register_3d_xte_field(Grid, Time, 'skew_mass_diff_13_h_xz',          &
         'Triad-averaged 13 component of mass weighted skew diffusivity at Eastern tracer-cell face',&
         'kg m^2/s', range=(/-1e8,1.e8/))
    id_skew_mass_diff_23_h_yz = register_3d_ytn_field(Grid, Time, 'skew_mass_diff_23_h_yz',            &
         'Triad-averaged 13 component of mass weighted skew diffusivity at Northern tracer-cell face', &
         'kg m^2/s', range=(/-1e8,1.e8/))
    id_skew_mass_diff_31_v_xz = register_3d_ztb_field(Grid, Time, 'skew_mass_diff_31_v_xz',          &
         'Triad-averaged 31 component of mass weighted skew diffusivity at bottom tracer-cell face', &
         'kg m^2/s', range=(/-1e8,1.e8/))
    id_skew_mass_diff_32_v_yz = register_3d_ztb_field(Grid, Time, 'skew_mass_diff_32_v_yz',          &
         'Triad-averaged 32 component of mass weighted skew diffusivity at bottom tracer-cell face', &
         'kg m^2/s', range=(/-1e8,1.e8/))

    id_m33_explicit_ztb = register_3d_ztb_field(Grid, Time, 'm33_explicit_ztb',                   &
         'Explicit 33 component of mass weighted neutral diffusivity at bottom tracer-cell face', &
         'kg m^2/s', range=(/-1e8,1.e8/))
    id_m33_implicit = register_3d_t_field(Grid, Time, 'm33_implicit',                        &
         'Implicit 33 component of mass weighted neutral diffusivity at tracer-cell center', &
         'kg m^2/s', range=(/-1e8,1.e8/))

    allocate(id_diff_flux_x_xte(num_prog_tracers))
    allocate(id_diff_flux_y_ytn(num_prog_tracers))
    allocate(id_diff_flux_z_ztb(num_prog_tracers))
    allocate(id_skew_flux_x_xte(num_prog_tracers))
    allocate(id_skew_flux_y_ytn(num_prog_tracers))
    allocate(id_skew_flux_z_ztb(num_prog_tracers))

    do n = 1, num_prog_tracers
       id_diff_flux_x_xte(n) = register_3d_xte_field(Grid, Time, 'diff_flux_x_xte_'//trim(T_prog(n)%name),          &
            'Density-area weighted flux of '//trim(T_prog(n)%name)//' through Eastern face from neutral diffusion', &
            'rho_dzt*dyt*flux of '//trim(T_prog(n)%name), (/-1.e18,1.e18/))
       id_diff_flux_y_ytn(n) = register_3d_ytn_field(Grid, Time, 'diff_flux_y_ytn_'//trim(T_prog(n)%name),           &
            'Density-area weighted flux of '//trim(T_prog(n)%name)//' through Northern face from neutral diffusion', &
            'rho_dzt*dxt*flux of '//trim(T_prog(n)%name), (/-1.e18,1.e18/))
       id_diff_flux_z_ztb(n) = register_3d_ztb_field(Grid, Time, 'diff_flux_z_ztb_'//trim(T_prog(n)%name),           &
            'Density-area weighted flux of '//trim(T_prog(n)%name)//' through bottom face from neutral diffusivion', &
            'rho*dxt*dyt*flux of '//trim(T_prog(n)%name), (/-1.e18,1.e18/))

       id_skew_flux_x_xte(n) = register_3d_xte_field(Grid, Time, 'skew_flux_x_xte_'//trim(T_prog(n)%name),       &
            'Density-area weighted flux of '//trim(T_prog(n)%name)//' through Eastern face from skew diffusion', &
            'rho_dzt*dyt*flux of '//trim(T_prog(n)%name), (/-1.e18,1.e18/))
       id_skew_flux_y_ytn(n) = register_3d_ytn_field(Grid, Time, 'skew_flux_y_ytn_'//trim(T_prog(n)%name),        &
            'Density-area weighted flux of '//trim(T_prog(n)%name)//' through Northern face from skew diffusion', &
            'rho_dzt*dxt*flux of '//trim(T_prog(n)%name), (/-1.e18,1.e18/))
       id_skew_flux_z_ztb(n) = register_3d_ztb_field(Grid, Time, 'skew_flux_z_ztb_'//trim(T_prog(n)%name),        &
            'Density-area weighted flux of '//trim(T_prog(n)%name)//' through bottom face from skew diffusivion', &
            'rho*dxt*dyt*flux of '//trim(T_prog(n)%name), (/-1.e18,1.e18/))
    enddo

    allocate(id_diff_th_tendency(num_prog_tracers))
    allocate(id_skew_th_tendency(num_prog_tracers))

    do n = 1, num_prog_tracers
       id_diff_th_tendency(n) = register_3d_t_field(Grid, Time, 'diff_th_tendency_'//trim(T_prog(n)%name), &
            'Thickness density weighted tendency of '//trim(T_prog(n)%name)//' due to neutral diffusion',  &
            'rho_dzt*'//trim(T_prog(n)%units)//'/s', (/-1.e18,1.e18/))
       id_skew_th_tendency(n) = register_3d_t_field(Grid, Time, 'skew_th_tendency_'//trim(T_prog(n)%name), &
            'Thickness density weighted tendency of '//trim(T_prog(n)%name)//' due to skew diffusion',     &
            'rho_dzt*'//trim(T_prog(n)%units)//'/s', (/-1.e18,1.e18/))
    enddo

  end subroutine register_fields
  ! </SUBROUTINE> NAME="register_fields"

  !#######################################################################
  ! <SUBROUTINE NAME="flux_calculations">
  !
  ! <DESCRIPTION>
  ! This function computes the thickness weighted tendency of tracers
  ! due to neutral physics as well as the implicit vertical diffusivity
  ! term.
  ! </DESCRIPTION>
  !
  subroutine flux_calculations(Domain, Grid, Time, Thickness, dtime,               &
       symm_tensor13_xz, symm_tensor23_yz, symm_tensor33_xz, symm_tensor33_yz,     &
       symm_tensor1122_z, upsilonx_xz, upsilony_yz, dTdx, dTdy, dTdz, aredi_array, &
       total_th_tendency, T_prog)

    type(ocean_domain_type),     intent(inout) :: Domain
    type(ocean_grid_type),       intent(in)    :: Grid
    type(ocean_time_type),       intent(in)    :: Time
    type(ocean_thickness_type),  intent(in)    :: Thickness
    real,                        intent(in)    :: dtime

    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: symm_tensor13_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: symm_tensor23_yz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: symm_tensor33_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: symm_tensor33_yz
    real, dimension(isd:,jsd:,:,0:),    intent(in) :: symm_tensor1122_z
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: upsilonx_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: upsilony_yz

    real, dimension(isd:,jsd:,:,:),  intent(in) :: dTdx
    real, dimension(isd:,jsd:,:,:),  intent(in) :: dTdy
    real, dimension(isd:,jsd:,0:,:), intent(in) :: dTdz

    real, dimension(isd:,jsd:,:), intent(in) :: aredi_array

    real, dimension(isd:,jsd:,:,:),             intent(out)   :: total_th_tendency
    type(ocean_prog_tracer_type), dimension(:), intent(inout) :: T_prog

    ! Mass weighted neutral diffusion tensor components [kgm^2/s]
    real, dimension(isd:ied,jsd:jed,nk)         :: symm_mass_diff_11_xte
    real, dimension(isd:ied,jsd:jed,nk)         :: symm_mass_diff_22_ytn
    real, dimension(isd:ied,jsd:jed,nk,0:1)     :: symm_mass_diff_33_z
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: symm_mass_diff_13_h_xz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: symm_mass_diff_23_h_yz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: symm_mass_diff_31_v_xz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: symm_mass_diff_32_v_yz

    ! Mass weighted skew diffusion tensor components [kgm^2/s]
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: skew_mass_diff_13_h_xz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: skew_mass_diff_23_h_yz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: skew_mass_diff_31_v_xz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: skeW_mass_diff_32_v_yz

    ! non-tensor based value of (kS) for used horizontal on diffusion [m^2/s]
    real, dimension(isd:ied,jsd:jed,nk) :: simple_mass_diff_11_xte 
    real, dimension(isd:ied,jsd:jed,nk) :: simple_mass_diff_22_ytn

    ! [kg*T/s]
    real, dimension(isd:ied,jsd:jed,nk,size(T_prog))   :: diff_flux_x_xte, diff_flux_y_ytn
    real, dimension(isd:ied,jsd:jed,0:nk,size(T_prog)) :: diff_flux_z_ztb
    real, dimension(isd:ied,jsd:jed,nk,size(T_prog))   :: skew_flux_x_xte, skew_flux_y_ytn
    real, dimension(isd:ied,jsd:jed,0:nk,size(T_prog)) :: skew_flux_z_ztb

    real, dimension(isd:ied,jsd:jed,nk) :: m33_explicit_ztb
    integer :: n, num_prog_tracers

    num_prog_tracers = size(T_prog)

    call mpp_clock_begin(id_clock_flux_mass_diff)
    call compute_mass_diff(Time, symm_tensor13_xz, symm_tensor23_yz, symm_tensor33_xz, &
         symm_tensor33_yz, symm_tensor1122_z, upsilonx_xz, upsilony_yz,                &
         aredi_array, Time%tau, Thickness, Grid, Domain,                               &
         simple_mass_diff_11_xte, simple_mass_diff_22_ytn,                             &
         symm_mass_diff_11_xte, symm_mass_diff_22_ytn,                                 &
         symm_mass_diff_13_h_xz, symm_mass_diff_23_h_yz,                               &
         symm_mass_diff_31_v_xz, symm_mass_diff_32_v_yz,                               &
         symm_mass_diff_33_z, skew_mass_diff_13_h_xz,                                  &
         skew_mass_diff_23_h_yz, skew_mass_diff_31_v_xz, skew_mass_diff_32_v_yz)
    call mpp_clock_end(id_clock_flux_mass_diff)

    call mpp_clock_begin(id_clock_flux_33_term)
    call compute_33_term(Time, Grid, symm_mass_diff_33_z, dtime, Time%tau, Thickness, &
         m33_explicit_ztb, T_prog)
    call mpp_clock_end(id_clock_flux_33_term)

    call mpp_clock_begin(id_clock_flux_fluxes)
    call compute_fluxes(symm_mass_diff_11_xte, symm_mass_diff_13_h_xz,     &
         symm_mass_diff_22_ytn, symm_mass_diff_23_h_yz,                    &
         symm_mass_diff_31_v_xz, symm_mass_diff_32_v_yz, m33_explicit_ztb, &
         skew_mass_diff_13_h_xz, skew_mass_diff_23_h_yz,                   & 
         skew_mass_diff_31_v_xz, skew_mass_diff_32_v_yz,                   &
         dTdx, dTdy, dTdz, Grid, Thickness,                                &
         diff_flux_x_xte, diff_flux_y_ytn, diff_flux_z_ztb,                &
         skew_flux_x_xte, skew_flux_y_ytn, skew_flux_z_ztb, T_prog)
    call mpp_clock_end(id_clock_flux_fluxes)

    call mpp_clock_begin(id_clock_flux_tracer_limits)
    call apply_tracer_limits(simple_mass_diff_11_xte,       &
         simple_mass_diff_22_ytn, Grid, dTdx, dTdy,         &
         diff_flux_x_xte, diff_flux_y_ytn, diff_flux_z_ztb, &
         skew_flux_x_xte, skew_flux_y_ytn, skew_flux_z_ztb, T_prog)
    call mpp_clock_end(id_clock_flux_tracer_limits)

    call mpp_clock_begin(id_clock_flux_misc)
    if (Grid%tripolar) then
       do n = 1, num_prog_tracers
          call mpp_update_domains(diff_flux_x_xte(:,:,:,n), diff_flux_y_ytn(:,:,:,n), &
               Dom_flux%domain2d, gridtype=CGRID_NE, complete=T_prog(n)%complete)
          call mpp_update_domains(skew_flux_x_xte(:,:,:,n), skew_flux_y_ytn(:,:,:,n), &
               Dom_flux%domain2d, gridtype=CGRID_NE, complete=T_prog(n)%complete)
       enddo
    endif

    ! minus sign is a MOM convention.
    do n = 1, num_prog_tracers
       call diagnose_3d(Time, Grid, id_diff_flux_x_xte(n), -diff_flux_x_xte(:,:,:,n))
       call diagnose_3d(Time, Grid, id_diff_flux_y_ytn(n), -diff_flux_y_ytn(:,:,:,n))
       call diagnose_3d(Time, Grid, id_diff_flux_z_ztb(n), -diff_flux_z_ztb(:,:,1:,n))
       call diagnose_3d(Time, Grid, id_skew_flux_x_xte(n), -skew_flux_x_xte(:,:,:,n))
       call diagnose_3d(Time, Grid, id_skew_flux_y_ytn(n), -skew_flux_y_ytn(:,:,:,n))
       call diagnose_3d(Time, Grid, id_skew_flux_z_ztb(n), -skew_flux_z_ztb(:,:,1:,n))
    enddo
    call mpp_clock_end(id_clock_flux_misc)

    call mpp_clock_begin(id_clock_flux_update)
    call update_tendencies(Time, Grid, diff_flux_x_xte, &
         diff_flux_y_ytn, diff_flux_z_ztb,              &
         skew_flux_x_xte, skew_flux_y_ytn,              &
         skew_flux_z_ztb, total_th_tendency)
    call mpp_clock_end(id_clock_flux_update)

  end subroutine flux_calculations
  ! </SUBROUTINE> NAME="flux_calculations"

  !#######################################################################
  ! <SUBROUTINE NAME="compute_mass_diff">
  !
  ! <DESCRIPTION>
  ! Subroutine computes the vertical neutral diffusion tracer flux component.
  ! Compute this component for all tracers at level k.
  ! Surface and bottom boundary condition fz(k=0)=fz(k=kmt(i,j))=0
  !
  ! fz has physical dimensions (density*diffusivity*tracer gradient)
  !
  ! This is nearly the same as the subroutine in ocean_nphysicsA.
  !
  ! </DESCRIPTION>
  !
  subroutine compute_mass_diff(Time, symm_tensor13_xz, symm_tensor23_yz, &
         symm_tensor33_xz, symm_tensor33_yz, symm_tensor1122_z,          &
         upsilonx_xz, upsilony_yz,                                       &
         aredi_array, tau, Thickness, Grid, Domain,                      &
         simple_mass_diff_11_xte, simple_mass_diff_22_ytn,               &
         symm_mass_diff_11_xte,  symm_mass_diff_22_ytn,                  &
         symm_mass_diff_13_h_xz, symm_mass_diff_23_h_yz,                 &
         symm_mass_diff_31_v_xz, symm_mass_diff_32_v_yz,                 &
         symm_mass_diff_33_z, skew_mass_diff_13_h_xz,                    &
         skew_mass_diff_23_h_yz, skew_mass_diff_31_v_xz, skew_mass_diff_32_v_yz)
   
    type(ocean_time_type),              intent(in) :: Time
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: symm_tensor13_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: symm_tensor23_yz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: symm_tensor33_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: symm_tensor33_yz
    real, dimension(isd:,jsd:,:,0:),    intent(in) :: symm_tensor1122_z
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: upsilonx_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: upsilony_yz

    real, dimension(isd:,jsd:,:), intent(in)    :: aredi_array
    integer,                      intent(in)    :: tau
    type(ocean_thickness_type),   intent(in)    :: Thickness
    type(ocean_grid_type),        intent(in)    :: Grid
    type(ocean_domain_type),      intent(inout) :: Domain

    real, dimension(isd:,jsd:,:),       intent(out) :: simple_mass_diff_11_xte
    real, dimension(isd:,jsd:,:),       intent(out) :: simple_mass_diff_22_ytn
    real, dimension(isd:,jsd:,:),       intent(out) :: symm_mass_diff_11_xte
    real, dimension(isd:,jsd:,:),       intent(out) :: symm_mass_diff_22_ytn
    real, dimension(isd:,jsd:,:,0:),    intent(out) :: symm_mass_diff_33_z

    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: symm_mass_diff_13_h_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: symm_mass_diff_23_h_yz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: symm_mass_diff_31_v_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: symm_mass_diff_32_v_yz

    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: skew_mass_diff_13_h_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: skew_mass_diff_23_h_yz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: skew_mass_diff_31_v_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: skew_mass_diff_32_v_yz

    ! Density weighted quarter cell volumes [kg]
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: rho_qcv_xz, rho_qcv_yz

    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: simple_mass_diff_11_h_xz, simple_mass_diff_22_h_yz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: symm_mass_diff_11_h_xz, symm_mass_diff_22_h_yz

    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: simple_mass_diff_11_xz, simple_mass_diff_22_yz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: symm_mass_diff_11_xz, symm_mass_diff_22_yz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: symm_mass_diff_13_xz, symm_mass_diff_23_yz
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: skew_mass_diff_13_xz, skew_mass_diff_23_yz

    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: full_tensor1122 ! [m^2/s]

    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: full_aredi ! [m^2/s]

    call geometric_terms(Grid, Thickness%rho_dzt(:,:,:,tau), Thickness%dzt(:,:,:), rho_qcv_xz, rho_qcv_yz)

    full_aredi(COMP,:,:,:)      = spread(spread(aredi_array(COMP,:), 4, 2), 5, 2)
    full_tensor1122(COMP,:,:,:) = spread(symm_tensor1122_z(COMP,:,:), 5, 2)

    simple_mass_diff_11_xz(COMP,:,:,:) = full_aredi(COMP,:,:,:)*rho_qcv_xz(COMP,:,:,:)
    simple_mass_diff_22_yz(COMP,:,:,:) = full_aredi(COMP,:,:,:)*rho_qcv_yz(COMP,:,:,:)
    symm_mass_diff_11_xz(COMP,:,:,:)   = full_tensor1122(COMP,:,:,:)*rho_qcv_xz(COMP,:,:,:)
    symm_mass_diff_22_yz(COMP,:,:,:)   = full_tensor1122(COMP,:,:,:)*rho_qcv_yz(COMP,:,:,:)
    symm_mass_diff_13_xz(COMP,:,:,:)   = symm_tensor13_xz(COMP,:,:,:)*rho_qcv_xz(COMP,:,:,:)
    symm_mass_diff_23_yz(COMP,:,:,:)   = symm_tensor23_yz(COMP,:,:,:)*rho_qcv_yz(COMP,:,:,:)
    skew_mass_diff_13_xz(COMP,:,:,:)   = upsilonx_xz(COMP,:,:,:)*rho_qcv_xz(COMP,:,:,:)
    skew_mass_diff_23_yz(COMP,:,:,:)   = upsilony_yz(COMP,:,:,:)*rho_qcv_yz(COMP,:,:,:)
    call mpp_update_domains(simple_mass_diff_11_xz, Domain%domain2d)
    call mpp_update_domains(simple_mass_diff_22_yz, Domain%domain2d)
    call mpp_update_domains(symm_mass_diff_11_xz, Domain%domain2d)
    call mpp_update_domains(symm_mass_diff_22_yz, Domain%domain2d)
    call mpp_update_domains(symm_mass_diff_13_xz, Domain%domain2d)
    call mpp_update_domains(symm_mass_diff_23_yz, Domain%domain2d)
    call mpp_update_domains(skew_mass_diff_13_xz, Domain%domain2d)
    call mpp_update_domains(skew_mass_diff_23_yz, Domain%domain2d)

    ! off diagonal terms
    call stencil_centre_to_horiz(symm_mass_diff_13_xz,   symm_mass_diff_23_yz, &
                                 symm_mass_diff_13_h_xz, symm_mass_diff_23_h_yz)
    call stencil_centre_to_vert( symm_mass_diff_13_xz,   symm_mass_diff_23_yz, &
                                 symm_mass_diff_31_v_xz, symm_mass_diff_32_v_yz)
    call stencil_centre_to_horiz(skew_mass_diff_13_xz,   skew_mass_diff_23_yz,  &
                                 skew_mass_diff_13_h_xz, skew_mass_diff_23_h_yz)
    call stencil_centre_to_vert(-skew_mass_diff_13_xz,  -skew_mass_diff_23_yz,  &
                                 skew_mass_diff_31_v_xz, skew_mass_diff_32_v_yz)

    ! horizontal diagonal terms    
    call stencil_centre_to_horiz(symm_mass_diff_11_xz, symm_mass_diff_22_yz, &
                                 symm_mass_diff_11_h_xz, symm_mass_diff_22_h_yz)
    symm_mass_diff_11_xte(COMPXL,:) = sum(sum(symm_mass_diff_11_h_xz(COMPXL,:,:,:), dim=5), dim=4)
    symm_mass_diff_22_ytn(COMPYL,:) = sum(sum(symm_mass_diff_22_h_yz(COMPYL,:,:,:), dim=5), dim=4)

    ! Simple horizontal
    call stencil_centre_to_horiz(simple_mass_diff_11_xz, simple_mass_diff_22_yz, &
                                 simple_mass_diff_11_h_xz, simple_mass_diff_22_h_yz)
    simple_mass_diff_11_xte(COMPXLYL,:) = sum(sum(simple_mass_diff_11_h_xz(COMPXLYL,:,:,:), dim=5), dim=4)
    simple_mass_diff_22_ytn(COMPXLYL,:) = sum(sum(simple_mass_diff_22_h_yz(COMPXLYL,:,:,:), dim=5), dim=4)

    ! Sum over each half cell
    symm_mass_diff_33_z(COMP,:,:) = sum(symm_tensor33_xz(COMP,:,:,:)*rho_qcv_xz(COMP,:,:,:) &
         +                              symm_tensor33_yz(COMP,:,:,:)*rho_qcv_yz(COMP,:,:,:), dim=4)

    call diagnose_3d(Time, Grid, id_symm_mass_diff_11_xte, symm_mass_diff_11_xte)
    call diagnose_3d(Time, Grid, id_symm_mass_diff_22_ytn, symm_mass_diff_22_ytn)

    call diagnose_3d(Time, Grid, id_symm_mass_diff_13_h_xz, sum(sum(symm_mass_diff_13_h_xz, dim=5), dim=4))
    call diagnose_3d(Time, Grid, id_symm_mass_diff_23_h_yz, sum(sum(symm_mass_diff_23_h_yz, dim=5), dim=4))
    call diagnose_3d(Time, Grid, id_symm_mass_diff_31_v_xz, sum(sum(symm_mass_diff_31_v_xz, dim=5), dim=4))
    call diagnose_3d(Time, Grid, id_symm_mass_diff_32_v_yz, sum(sum(symm_mass_diff_32_v_yz, dim=5), dim=4))
    
    call diagnose_3d(Time, Grid, id_skew_mass_diff_13_h_xz, sum(sum(skew_mass_diff_13_h_xz, dim=5), dim=4))
    call diagnose_3d(Time, Grid, id_skew_mass_diff_23_h_yz, sum(sum(skew_mass_diff_23_h_yz, dim=5), dim=4))
    call diagnose_3d(Time, Grid, id_skew_mass_diff_31_v_xz, sum(sum(skew_mass_diff_31_v_xz, dim=5), dim=4))
    call diagnose_3d(Time, Grid, id_skew_mass_diff_32_v_yz, sum(sum(skew_mass_diff_32_v_yz, dim=5), dim=4))

  end subroutine compute_mass_diff
  ! </SUBROUTINE> NAME="compute_mass_diff"


  !#######################################################################
  ! <SUBROUTINE NAME="geometric_terms">
  !
  ! <DESCRIPTION>
  ! Calculate the density weighted quarter cell volumes of the triads.
  ! </DESCRIPTION>
  !
  subroutine geometric_terms(Grid, rho_dzt, dzt, &
       rho_qcv_xz, rho_qcv_yz)

    type(ocean_grid_type),        intent(in) :: Grid
    real, dimension(isd:,jsd:,:), intent(in) :: rho_dzt
    real, dimension(isd:,jsd:,:), intent(in) :: dzt

    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: rho_qcv_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: rho_qcv_yz

    integer :: k, ip, kr

    ! density calculated by removing dzt from rho_dzt [kg/m^3]
    real, dimension(isd:ied,jsd:jed,nk) :: rho

    ! minimum half cell thickness between adjacent tracer cells [m]
    real, dimension(isd:ied,jsd:jed,nk,0:1) :: min_dzt_xte_z
    real, dimension(isd:ied,jsd:jed,nk,0:1) :: min_dzt_ytn_z

    ! quarter cell area [m^2]
    real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: qca_xz, qca_yz

    do k = 1, nk
       do kr = 0, 1
          min_dzt_xte_z(:,:,k,kr) = FMX(Grid%fracdz(k,kr)*dzt(:,:,k)) ! G(2004): 16.53
          min_dzt_ytn_z(:,:,k,kr) = FMY(Grid%fracdz(k,kr)*dzt(:,:,k)) ! G(2004): 16.54
       enddo
    enddo

    do k = 1, nk
       qca_xz(COMP,k,0,0) = Grid%dtw(COMP)*min_dzt_xte_z(COMPXLL,k,0) ! ~ V(1)
       qca_xz(COMP,k,0,1) = Grid%dtw(COMP)*min_dzt_xte_z(COMPXLL,k,1) ! ~ V(3)
       qca_xz(COMP,k,1,0) = Grid%dte(COMP)*min_dzt_xte_z(COMP,k,0)    ! ~ V(2)
       qca_xz(COMP,k,1,1) = Grid%dte(COMP)*min_dzt_xte_z(COMP,k,1)    ! ~ V(4)

       qca_yz(COMP,k,0,0) = Grid%dts(COMP)*min_dzt_ytn_z(COMPYLL,k,0)
       qca_yz(COMP,k,0,1) = Grid%dts(COMP)*min_dzt_ytn_z(COMPYLL,k,1)
       qca_yz(COMP,k,1,0) = Grid%dtn(COMP)*min_dzt_ytn_z(COMP,k,0)
       qca_yz(COMP,k,1,1) = Grid%dtn(COMP)*min_dzt_ytn_z(COMP,k,1)
    enddo

    rho(COMP,:) = (rho_dzt(COMP,:)/(dzt(COMP,:) + epsln))
    do k = 1, nk
       do ip = 0, 1
          do kr = 0, 1
             rho_qcv_xz(COMP,k,ip,kr) = rho(COMP,k)*qca_xz(COMP,k,ip,kr)*Grid%dyt(COMP)
             rho_qcv_yz(COMP,k,ip,kr) = rho(COMP,k)*qca_yz(COMP,k,ip,kr)*Grid%dxt(COMP)
          enddo
       enddo
    enddo

  end subroutine geometric_terms
  ! </SUBROUTINE> NAME="geometric_terms"

  !#######################################################################
  ! <SUBROUTINE NAME="compute_33_term">
  !
  ! <DESCRIPTION>
  ! K33 is the (3,3) term in small angle Redi diffusion tensor.
  ! It is broken into an explicit in time piece and implicit
  ! in time piece.  It is weighted by density for non-Boussinesq
  ! and rho0 for Boussinesq.
  !
  ! K33 has units (kg/m^3)*m^2/sec.
  ! </DESCRIPTION>
  subroutine compute_33_term(Time, Grid, symm_mass_diff_33_z, dtime, tau, Thickness, &
       m33_explicit_ztb, T_prog)

    type(ocean_time_type),              intent(in) :: Time
    type(ocean_grid_type),              intent(in) :: Grid
    real, dimension(isd:,jsd:,:,0:),    intent(in) :: symm_mass_diff_33_z
    real,                               intent(in) :: dtime
    integer,                            intent(in) :: tau
    type(ocean_thickness_type),         intent(in) :: Thickness

    real, dimension(isd:,jsd:,:), intent(out) :: m33_explicit_ztb
    type(ocean_prog_tracer_type), dimension(:), intent(inout) :: T_prog

    integer :: k, kr
    integer :: n, num_prog_tracers

    real, dimension(isd:ied,jsd:jed,nk,0:1) :: m33_explicit_z, m33_crit_z
    real, dimension(isd:ied,jsd:jed,nk) :: m33_implicit

    ! [(kg/m^3)*(m^2/s)]
    real, dimension(isd:ied,jsd:jed,nk) :: rho_k33_implicit

    num_prog_tracers = size(T_prog)

    ! Compute the explicit component of m33
    if (.not. diffusion_all_explicit) then
       ! Critical value at the tracer point
       do k = 1, nk
          do kr = 0, 1
             m33_crit_z(COMP,k,kr) = &
             (0.5/dtime)*(Grid%fracdz(k,kr)*Thickness%dzt(COMP,k))**2 &
             *(Grid%fracdz(k,kr)*Thickness%rho_dzt(COMP,k,tau)*Grid%dat(COMP))
          enddo
       enddo
       m33_explicit_z(COMP,:,:) = min(symm_mass_diff_33_z(COMP,:,:), m33_crit_z(COMP,:,:))
    else
       m33_explicit_z(COMP,:,:) = symm_mass_diff_33_z(COMP,:,:)
    endif

    ! Tracer-centred implicit component for use in other modules
    m33_implicit(COMP,:) = sum((symm_mass_diff_33_z(COMP,:,:) - m33_explicit_z(COMP,:,:)), dim=4)
    do k = 1, nk
       rho_k33_implicit(COMP,k) = &
       m33_implicit(COMP,k)*Grid%datr(COMP)/(Thickness%dzt(COMP,k) + epsln)
    enddo
    do n = 1, num_prog_tracers
       T_prog(n)%K33_implicit(:,:,:) = rho_k33_implicit(:,:,:)
    enddo

    ! Cell-bottom centred explicit component for use in this module.    
    do k = 1, nk-1
       m33_explicit_ztb(COMP,k) = m33_explicit_z(COMP,k,1) + m33_explicit_z(COMP,k+1,0)
    enddo

    call diagnose_3d(Time, Grid, id_m33_explicit_ztb, m33_explicit_ztb)
    call diagnose_3d(Time, Grid, id_m33_implicit, m33_implicit)
  
  end subroutine compute_33_term
  ! </SUBROUTINE> NAME="compute_33_term"

  
  !#######################################################################
  ! <SUBROUTINE NAME="compute_fluxes">
  !
  ! <DESCRIPTION>
  ! Computes the tracer fluxes due to neutral diffusion and skew diffusion.
  ! Fluxes are computed at the tracer cell faces and have units of [kg*T/s].
  ! </DESCRIPTION>
  !
  subroutine compute_fluxes(symm_mass_diff_11_xte, symm_mass_diff_13_h_xz,     &
       symm_mass_diff_22_ytn,  symm_mass_diff_23_h_yz,                         &
       symm_mass_diff_31_v_xz, symm_mass_diff_32_v_yz, m33_explicit_ztb,       &
       skew_mass_diff_13_h_xz, skew_mass_diff_23_h_yz, skew_mass_diff_31_v_xz, &
       skew_mass_diff_32_v_yz, dTdx, dTdy, dTdz, Grid, Thickness,              &
       diff_flux_x_xte, diff_flux_y_ytn, diff_flux_z_ztb, skew_flux_x_xte,     &
       skew_flux_y_ytn, skew_flux_z_ztb, T_prog)

    real, dimension(isd:,jsd:,:),       intent(in) :: symm_mass_diff_11_xte
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: symm_mass_diff_13_h_xz
    real, dimension(isd:,jsd:,:),       intent(in) :: symm_mass_diff_22_ytn
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: symm_mass_diff_23_h_yz

    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: symm_mass_diff_31_v_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: symm_mass_diff_32_v_yz
    real, dimension(isd:,jsd:,:),       intent(in) :: m33_explicit_ztb

    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: skew_mass_diff_13_h_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: skew_mass_diff_23_h_yz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: skew_mass_diff_31_v_xz
    real, dimension(isd:,jsd:,:,0:,0:), intent(in) :: skew_mass_diff_32_v_yz

    real, dimension(isd:,jsd:,:,:),  intent(in) :: dTdx
    real, dimension(isd:,jsd:,:,:),  intent(in) :: dTdy
    real, dimension(isd:,jsd:,0:,:), intent(in) :: dTdz
    type(ocean_grid_type),           intent(in) :: Grid
    type(ocean_thickness_type),      intent(in) :: Thickness

    real, dimension(isd:,jsd:,:,:),  intent(out) :: diff_flux_x_xte
    real, dimension(isd:,jsd:,:,:),  intent(out) :: diff_flux_y_ytn
    real, dimension(isd:,jsd:,0:,:), intent(out) :: diff_flux_z_ztb
    real, dimension(isd:,jsd:,:,:),  intent(out) :: skew_flux_x_xte
    real, dimension(isd:,jsd:,:,:),  intent(out) :: skew_flux_y_ytn
    real, dimension(isd:,jsd:,0:,:), intent(out) :: skew_flux_z_ztb

    type(ocean_prog_tracer_type), dimension(:), intent(inout) :: T_prog

    integer :: k, n, num_prog_tracers
    ! Components of z-flux from dTdz and dTdz respectivly [(rho*dz*T/s)/m^2]
    real, dimension(isd:ied,jsd:jed,nk) :: flux_x_ztb, flux_y_ztb

    num_prog_tracers = size(T_prog)

    ! Compute the flux due to neutral diffusion
    diff_flux_x_xte(:,:,:,:) = 0.0
    diff_flux_y_ytn(:,:,:,:) = 0.0
    diff_flux_z_ztb(:,:,:,:) = 0.0
    do n = 1, num_prog_tracers
       do k = 1, nk
          diff_flux_x_xte(COMPXL,k,n) = -(symm_mass_diff_11_xte(COMPXL,k)*dTdx(COMPXL,k,n) &
               + symm_mass_diff_13_h_xz(COMPXL,k,0,0)*dTdz(COMPXL,k-1,n)                   &
               + symm_mass_diff_13_h_xz(COMPXL,k,0,1)*dTdz(COMPXL,k  ,n)                   &
               + symm_mass_diff_13_h_xz(COMPXL,k,1,0)*dTdz(COMPXR,k-1,n)                   &
               + symm_mass_diff_13_h_xz(COMPXL,k,1,1)*dTdz(COMPXR,k  ,n))*Grid%dxter(COMPXL)
          diff_flux_y_ytn(COMPYL,k,n) = -(symm_mass_diff_22_ytn(COMPYL,k)*dTdy(COMPYL,k,n) &
               + symm_mass_diff_23_h_yz(COMPYL,k,0,0)*dTdz(COMPYL,k-1,n)                   &
               + symm_mass_diff_23_h_yz(COMPYL,k,0,1)*dTdz(COMPYL,k,n)                     &
               + symm_mass_diff_23_h_yz(COMPYL,k,1,0)*dTdz(COMPYR,k-1,n)                   &
               + symm_mass_diff_23_h_yz(COMPYL,k,1,1)*dTdz(COMPYR,k,n))*Grid%dytnr(COMPYL)
       enddo

       flux_x_ztb(:,:,:) = 0.0
       flux_y_ztb(:,:,:) = 0.0
       do k = 1, nk-1
          flux_x_ztb(COMP,k) = symm_mass_diff_31_v_xz(COMP,k,0,0)*dTdx(COMPXLL,k,n)   &
               +               symm_mass_diff_31_v_xz(COMP,k,0,1)*dTdx(COMPXLL,k+1,n) &
               +               symm_mass_diff_31_v_xz(COMP,k,1,0)*dTdx(COMP,k,n)      &
               +               symm_mass_diff_31_v_xz(COMP,k,1,1)*dTdx(COMP,k+1,n)
          flux_y_ztb(COMP,k) = symm_mass_diff_32_v_yz(COMP,k,0,0)*dTdy(COMPYLL,k,n)   &
               +               symm_mass_diff_32_v_yz(COMP,k,0,1)*dTdy(COMPYLL,k+1,n) &
               +               symm_mass_diff_32_v_yz(COMP,k,1,0)*dTdy(COMP,k,n)      &
               +               symm_mass_diff_32_v_yz(COMP,k,1,1)*dTdy(COMP,k+1,n)
       enddo
       diff_flux_z_ztb(COMP,1:nk,n) = -((flux_x_ztb(COMP,:) + flux_y_ztb(COMP,:)  &
           + m33_explicit_ztb(COMP,:)*dTdz(COMP,1:nk,n))/(Thickness%dzwt(COMP,1:nk)+epsln))
    enddo

    ! Compute the flux due to skew diffusion (same algorithm as above, but without diagonal terms)
    skew_flux_x_xte(:,:,:,:) = 0.0
    skew_flux_y_ytn(:,:,:,:) = 0.0
    skew_flux_z_ztb(:,:,:,:) = 0.0
    do n = 1, num_prog_tracers
       do k = 1, nk
          skew_flux_x_xte(COMPXL,k,n) = -(                               &
               + skew_mass_diff_13_h_xz(COMPXL,k,0,0)*dTdz(COMPXL,k-1,n) &
               + skew_mass_diff_13_h_xz(COMPXL,k,0,1)*dTdz(COMPXL,k  ,n) &
               + skew_mass_diff_13_h_xz(COMPXL,k,1,0)*dTdz(COMPXR,k-1,n) &
               + skew_mass_diff_13_h_xz(COMPXL,k,1,1)*dTdz(COMPXR,k  ,n))*Grid%dxter(COMPXL)
          skew_flux_y_ytn(COMPYL,k,n) = -(                               &
               + skew_mass_diff_23_h_yz(COMPYL,k,0,0)*dTdz(COMPYL,k-1,n) &
               + skew_mass_diff_23_h_yz(COMPYL,k,0,1)*dTdz(COMPYL,k,n)   &
               + skew_mass_diff_23_h_yz(COMPYL,k,1,0)*dTdz(COMPYR,k-1,n) &
               + skew_mass_diff_23_h_yz(COMPYL,k,1,1)*dTdz(COMPYR,k,n))*Grid%dytnr(COMPYL)
       enddo

       flux_x_ztb(:,:,:) = 0.0
       flux_y_ztb(:,:,:) = 0.0
       do k = 1, nk-1
          flux_x_ztb(COMP,k) = skew_mass_diff_31_v_xz(COMP,k,0,0)*dTdx(COMPXLL,k,n)   &
               +               skew_mass_diff_31_v_xz(COMP,k,0,1)*dTdx(COMPXLL,k+1,n) &
               +               skew_mass_diff_31_v_xz(COMP,k,1,0)*dTdx(COMP,k,n)      &
               +               skew_mass_diff_31_v_xz(COMP,k,1,1)*dTdx(COMP,k+1,n)
          flux_y_ztb(COMP,k) = skew_mass_diff_32_v_yz(COMP,k,0,0)*dTdy(COMPYLL,k,n)   &
               +               skew_mass_diff_32_v_yz(COMP,k,0,1)*dTdy(COMPYLL,k+1,n) &
               +               skew_mass_diff_32_v_yz(COMP,k,1,0)*dTdy(COMP,k,n)      &
               +               skew_mass_diff_32_v_yz(COMP,k,1,1)*dTdy(COMP,k+1,n)
       enddo
       skew_flux_z_ztb(COMP,1:nk,n) = -((flux_x_ztb(COMP,:) + flux_y_ztb(COMP,:))/(Thickness%dzwt(COMP,1:nk)+epsln))
    enddo

  end subroutine compute_fluxes
  ! </SUBROUTINE> NAME="compute_fluxes"

  !#######################################################################
  ! <SUBROUTINE NAME="apply_tracer_limits">
  !
  ! <DESCRIPTION>
  ! If the neutral_physics_limit flag is set, then the flux used in regions of large
  ! tracer gradients (as defined by T_prog(n)tmask_limit) are set to have purely
  ! horizontal diffusion, with no vertical or skew terms.
  ! </DESCRIPTION>
  !
  subroutine apply_tracer_limits(simple_mass_diff_11_xte, simple_mass_diff_22_ytn, &
       Grid, dTdx, dTdy, diff_flux_x_xte, diff_flux_y_ytn, diff_flux_z_ztb,        &
       skew_flux_x_xte, skew_flux_y_ytn, skew_flux_z_ztb, T_prog)
    
    real, dimension(isd:,jsd:,:),    intent(in) :: simple_mass_diff_11_xte
    real, dimension(isd:,jsd:,:),    intent(in) :: simple_mass_diff_22_ytn
    type(ocean_grid_type),           intent(in) :: Grid
    real, dimension(isd:,jsd:,:,:),  intent(in) :: dTdx
    real, dimension(isd:,jsd:,:,:),  intent(in) :: dTdy

    real, dimension(isd:,jsd:,:,:),  intent(inout) :: diff_flux_x_xte
    real, dimension(isd:,jsd:,:,:),  intent(inout) :: diff_flux_y_ytn
    real, dimension(isd:,jsd:,0:,:), intent(inout) :: diff_flux_z_ztb
    real, dimension(isd:,jsd:,:,:),  intent(inout) :: skew_flux_x_xte
    real, dimension(isd:,jsd:,:,:),  intent(inout) :: skew_flux_y_ytn
    real, dimension(isd:,jsd:,0:,:), intent(inout) :: skew_flux_z_ztb
    type(ocean_prog_tracer_type), dimension(:), intent(inout) :: T_prog

    integer :: k, n, num_prog_tracers

    num_prog_tracers = size(T_prog)

    if (neutral_physics_limit) then
       do n = 1, num_prog_tracers
          do k = 1, nk
             where (T_prog(n)%tmask_limit(COMPXLYL,k) == 1.0)
                diff_flux_x_xte(COMPXLYL,k,n) =  &
                  -simple_mass_diff_11_xte(COMPXLYL,k)*dTdx(COMPXLYL,k,n)*Grid%dxter(COMPXLYL)
                diff_flux_y_ytn(COMPXLYL,k,n) =  &
                  -simple_mass_diff_22_ytn(COMPXLYL,k)*dTdy(COMPXLYL,k,n)*Grid%dytnr(COMPXLYL)
                skew_flux_x_xte(COMPXLYL,k,n) = 0.0
                skew_flux_y_ytn(COMPXLYL,k,n) = 0.0
             endwhere
          enddo

          where (T_prog(n)%tmask_limit(COMP,:) == 1.0) 
             diff_flux_z_ztb(COMP,1:nk,n)   = 0.0
             skew_flux_z_ztb(COMP,1:nk,n)   = 0.0
             T_prog(n)%K33_implicit(COMP,:) = 0.0
          endwhere
       enddo
    endif

  end subroutine apply_tracer_limits
  ! </SUBROUTINE> NAME="apply_tracer_limits"

  !#######################################################################
  ! <SUBROUTINE NAME="update_tendencies">
  !
  ! <DESCRIPTION>
  ! Update the tendency for each tracer in each cell based on the total flux
  ! flowing through each of the six cell faces. The tendency is calculated 
  ! separately for the flux due to neutral diffusion and due skew diffusion.
  ! </DESCRIPTION>
  !
  subroutine update_tendencies(Time, Grid, diff_flux_x_xte, diff_flux_y_ytn, diff_flux_z_ztb, &
       skew_flux_x_xte, skew_flux_y_ytn, skew_flux_z_ztb, total_th_tendency)

    type(ocean_time_type),             intent(in) :: Time
    type(ocean_grid_type),            intent(in)  :: Grid
    real, dimension(isd:,jsd:,:,:),   intent(in)  :: difF_flux_x_xte
    real, dimension(isd:,jsd:,:,:),   intent(in)  :: diff_flux_y_ytn
    real, dimension(isd:,jsd:,0:,:),  intent(in)  :: diff_flux_z_ztb
    real, dimension(isd:,jsd:,:,:),   intent(in)  :: skew_flux_x_xte
    real, dimension(isd:,jsd:,:,:),   intent(in)  :: skew_flux_y_ytn
    real, dimension(isd:,jsd:,0:,:),  intent(in)  :: skew_flux_z_ztb
    real, dimension(isd:,jsd:,:,:),   intent(out) :: total_th_tendency

    integer :: k, n, num_prog_tracers
    
    ! density thickness weighted tracer tendency due to neutral and skew diffusion [rho*dz*T/s]
    real, dimension(isd:ied,jsd:jed,nk) :: diff_th_tendency
    real, dimension(isd:ied,jsd:jed,nk) :: skew_th_tendency

    num_prog_tracers = size(total_th_tendency, dim=4)

    do n = 1, num_prog_tracers
       diff_th_tendency(:,:,:) = 0.0
       skew_th_tendency(:,:,:) = 0.0
       do k = 1, nk
          diff_th_tendency(COMP,k) =                                        &
               -(diff_flux_x_xte(COMP,k,n)   - diff_flux_x_xte(COMPXLL,k,n) &
               + diff_flux_y_ytn(COMP,k,n)   - diff_flux_y_ytn(COMPYLL,k,n) &
               + diff_flux_z_ztb(COMP,k-1,n) - diff_flux_z_ztb(COMP,k,n) )*Grid%datr(COMP)*Grid%tmask(COMP,k)
          skew_th_tendency(COMP,k) =                                        &
               -(skew_flux_x_xte(COMP,k,n)   - skew_flux_x_xte(COMPXLL,k,n) &
               + skew_flux_y_ytn(COMP,k,n)   - skew_flux_y_ytn(COMPYLL,k,n) &
               + skew_flux_z_ztb(COMP,k-1,n) - skew_flux_z_ztb(COMP,k,n) )*Grid%datr(COMP)*Grid%tmask(COMP,k)
       enddo
       total_th_tendency(COMP,:,n) = diff_th_tendency(COMP,:) + skew_th_tendency(COMP,:)
       call diagnose_3d(Time, Grid, id_diff_th_tendency(n), diff_th_tendency)
       call diagnose_3d(Time, Grid, id_skew_th_tendency(n), skew_th_tendency)
    enddo

  end subroutine update_tendencies
  ! </SUBROUTINE> NAME="update_tendencies"

end module ocean_nphysics_flux_mod

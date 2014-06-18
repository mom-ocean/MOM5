module ocean_tracer_advect_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Matt Harrison 
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> John Dunne 
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Alistair Adcroft 
!</CONTACT>
!
!<OVERVIEW>
! This module computes thickness weighted tracer advection tendencies. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes thickness weighted tracer advection tendencies
! using a variety of advection schemes.  
!</DESCRIPTION>
!
!<INFO>
!
!<REFERENCE>
! S.-J. Lin
! A "Vertically Lagrangian" Finite-Volume Dynamical Core for Global Models
! Month. Weather Rev. (2004) 132, 2293-2307
! (Appendix B)
!</REFERENCE>
!
!<REFERENCE>
! H.T. Huynh
! Schemes and Constraints for advection
! 15th Intern. Conf. on Numeric. Meth. in Fluid Mech., Springer (1997)
!</REFERENCE>
!
!<REFERENCE>
! Colella P. and P.R. Woodward
! The piecewise parabloic method (PPM) for gasdynamical simulations
! J. Comput. Phys. (1984) 54, 174-201
!</REFERENCE>
!
!<REFERENCE>
! A. Suresh and H.T. Huynh
! Accurate Monotonicity-Preserving Schemes with Runge-Kutta Time Splitting
! J. Comput. Phys. (1997) 136, 83-99
!</REFERENCE>
!
!<REFERENCE>
! V. Daru and  C. Tenaud
! High order one-step monotonicity-preserving schemes for unsteady
! compressible flow calculations.
! J. Comp. Phys. (2004) 193, 563-594
!</REFERENCE>
!
!<REFERENCE>
! R.C. Easter
! Two modified versions of Botts positive-definite numerical advection scheme.
! Month. Weath. Rev. (1993) 121, 297-304
!</REFERENCE>
!
!<REFERENCE>
! Prather, M. J.,"Numerical Advection by Conservation of Second-Order Moments"
! JGR, Vol 91, NO. D6, p 6671-6681, May 20, 1986
!</REFERENCE>
!
!<REFERENCE>
! Merryfield and Holloway (2003), "Application of an accurate advection
! algorithm to sea-ice modelling". Ocean Modelling, Vol 5, p 1-15.
!</REFERENCE>
!
!<REFERENCE>
! Hundsdorder and Trompert (1994), "Method of lines and 
! direct discretization: a comparison for linear 
! advection", Applied Numerical Mathematics,
! pages 469--490.
!</REFERENCE>
!
!<REFERENCE>
! Sweby (1984): "High-resolution schemes using flux 
! limiters for hyperbolic conservation laws", 
! SIAM Journal of Numerical Analysis, vol. 21
! pages 995-1011. 
!</REFERENCE>
!
!</INFO>
!
!<NOTE>
! Contains a version of quicker with MOM3 masking.
!</NOTE>
!
!<NAMELIST NAME="ocean_tracer_advect_nml">
!  <DATA NAME="limit_with_upwind" TYPE="logical">
!  If true, will compute tracer fluxes entering a cell using upwind 
!  if the tracer value is outside a specified range. Implemented
!  only for quick at this time. This is an ad hoc and incomplete attempt
!  to maintain monotonicity with the quicker scheme.  
!  </DATA> 
!  <DATA NAME="advect_sweby_all" TYPE="logical">
!  For running all tracers with sweby, thereby utilizing a bitwise same
!  routine that reorganizes loops and can be faster for certain configurations.  
!  Default advect_sweby_all=.false.
!  </DATA> 
!  <DATA NAME="zero_tracer_advect_horz" TYPE="logical">
!  For debugging.  Set to .true. to turn off horizontal advection.  
!  </DATA> 
!  <DATA NAME="zero_tracer_advect_vert" TYPE="logical">
!  For debugging.  Set to .true. to turn off vertical advection.  
!  </DATA> 
!  <DATA NAME="psom_limit_prather" TYPE="logical">
!  For running with the original Prather limiter for the PSOM scheme.
!  The limiter is positive definite, but not monotonic.  This limiter
!  is NOT recommended for most applications.  The default is 
!  psom_limit_prather=.false., since we prefer to use the limiter 
!  from Merryfield and Holloway (2003).
!  </DATA> 
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging 
!  </DATA> 
!
!  <DATA NAME="write_a_restart" TYPE="logical">
!  Set true to write a restart.  False setting only for rare 
!  cases where wish to benchmark model without measuring the cost
!  of writing restarts and associated chksums.  
!  Default is write_a_restart=.true. 
!  </DATA> 
!
!  <DATA NAME="read_basin_mask" TYPE="logical">
!  For reading in a mask that selects regions of the domain 
!  for performing gyre and overturning diagnostics. 
!  The basin-mask convention used at GFDL has 
!  Southern=1.0,Atlantic=2.0,Pacific=3.0,Arctic=4.0,Indian=5.0
!  Default read_basin_mask=.false., whereby basin_mask
!  is set to tmask(k=1). 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,       only: epsln
use diag_manager_mod,    only: register_diag_field, send_data
use fms_mod,             only: write_version_number, check_nml_error, close_file, open_namelist_file
use fms_mod,             only: read_data, field_exist, file_exist
use fms_io_mod,          only: register_restart_field, save_restart, restore_state
use fms_io_mod,          only: restart_file_type 
use mpp_domains_mod,     only: mpp_update_domains, mpp_define_domains, mpp_global_field
use mpp_domains_mod,     only: mpp_start_update_domains, mpp_complete_update_domains
use mpp_domains_mod,     only: BGRID_NE, CGRID_NE, domain2d, XUPDATE, YUPDATE
use mpp_domains_mod,     only: mpp_get_domain_components, mpp_get_layout, domain1d, mpp_get_pelist
use mpp_mod,             only: input_nml_file, mpp_error, FATAL, WARNING, NOTE, stdout, stdlog
use mpp_mod,             only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use mpp_mod,             only: mpp_sum

use ocean_domains_mod,     only: get_local_indices, set_ocean_domain
use ocean_obc_mod,         only: store_ocean_obc_tracer_flux, ocean_obc_zero_boundary
use ocean_parameters_mod,  only: ADVECT_UPWIND, ADVECT_2ND_ORDER, ADVECT_4TH_ORDER, ADVECT_6TH_ORDER
use ocean_parameters_mod,  only: ADVECT_QUICKER, ADVECT_QUICKMOM3, ADVECT_MDFL_SUP_B
use ocean_parameters_mod,  only: ADVECT_MDFL_SWEBY_TEST, ADVECT_MDFL_SWEBY
use ocean_parameters_mod,  only: ADVECT_MDPPM_TEST, ADVECT_MDPPM 
use ocean_parameters_mod,  only: ADVECT_DST_LINEAR_TEST, ADVECT_DST_LINEAR
use ocean_parameters_mod,  only: ADVECT_PSOM, ADVECT_MDMDT_TEST
use ocean_parameters_mod,  only: TWO_LEVEL
use ocean_parameters_mod,  only: missing_value, onesixth, rho0, rho0r
use ocean_topog_mod,       only: ocean_topog_init
use ocean_tracer_util_mod, only: tracer_psom_chksum
use ocean_types_mod,       only: ocean_domain_type, ocean_grid_type, ocean_density_type, ocean_time_type
use ocean_types_mod,       only: ocean_prog_tracer_type, ocean_thickness_type, ocean_adv_vel_type
use ocean_workspace_mod,   only: wrk1, wrk2, wrk3, wrk4, wrk5, wrk1_2d 
use ocean_util_mod,        only: diagnose_2d, diagnose_3d
use ocean_tracer_util_mod, only: diagnose_3d_rho

implicit none

private
integer :: id0, id1, now
integer :: id_clock_up_horz
integer :: id_clock_up_vert
integer :: id_clock_2nd_horz
integer :: id_clock_2nd_vert
integer :: id_clock_4th_horz
integer :: id_clock_4th_vert
integer :: id_clock_6th_horz
integer :: id_clock_6th_vert
integer :: id_clock_quick_horz
integer :: id_clock_quick_vert
integer :: id_clock_quickmom3_horz
integer :: id_clock_quickmom3_vert
integer :: id_clock_mdfl_sup_b
integer :: id_clock_dst_linear
integer :: id_clock_dst_linear_test
integer :: id_clock_mdfl_sweby
integer :: id_clock_mdfl_sweby_all
integer :: id_clock_mdfl_sweby_test
integer :: id_clock_psom
integer :: id_clock_psom_x
integer :: id_clock_psom_y
integer :: id_clock_psom_z
integer :: id_clock_mdppm 
integer :: id_clock_mdppm_test
integer :: id_clock_mdmdt_test 
integer :: id_clock_gyre_over
integer :: id_clock_adv_diss
integer :: id_clock_mdfl_sweby_mpi
integer :: id_clock_mdfl_sweby_cu1
integer :: id_clock_mdfl_sweby_cu2
integer :: id_clock_mdfl_sweby_cu3
integer :: id_clock_mdfl_sweby_dia

!for restart file
type(restart_file_type), save :: Adv_restart

#include <ocean_memory.h>

#ifdef MOM_STATIC_ARRAYS

real, dimension(isd:ied,jsd:jed,nk) :: flux_x
real, dimension(isd:ied,jsd:jed,nk) :: flux_y
real, dimension(isd:ied,jsd:jed,nk) :: flux_z

real, dimension(isd:ied,jsd:jed,nk) :: advect_tendency
real, dimension(isd:ied,jsd:jed,nk) :: neutral_temp_advect
real, dimension(isd:ied,jsd:jed,nk) :: neutral_salt_advect

! some fields for quicker 
real, dimension(isd:ied,jsd:jed,nk)         :: quick_x
real, dimension(isd:ied,jsd:jed,nk)         :: quick_y                    
real, dimension(nk,2)                       :: quick_z                             
real, dimension(isd:ied,jsd:jed,3)          :: curv_xp
real, dimension(isd:ied,jsd:jed,3)          :: curv_xn
real, dimension(isd:ied,jsd:jed,3)          :: curv_yp
real, dimension(isd:ied,jsd:jed,3)          :: curv_yn  
real, dimension(nk,3)                       :: curv_zp                   
real, dimension(nk,3)                       :: curv_zn                    
real, dimension(isc-2:iec+2,jsc-2:jec+2)    :: dxt_quick
real, dimension(isc-2:iec+2,jsc-2:jec+2)    :: dyt_quick
real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: tmask_quick
real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: tracer_quick

! fields for 4th adn 6th order 
real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: tmask_fourth
real, dimension(isc-3:iec+3,jsc-3:jec+3,nk) :: tmask_sixth
real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: tracer_fourth
real, dimension(isc-3:iec+3,jsc-3:jec+3,nk) :: tracer_sixth

! fields for mdfl 
real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: tracer_mdfl
real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: tmask_mdfl
real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: mass_mdfl
real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: tracermass_mdfl

! fields for mdppm
real, dimension(isc-4:iec+4,jsc-4:jec+4,nk) :: tmask_mdppm
real, dimension(isc-4:iec+4,jsc-4:jec+4,nk) :: tracer_mdppm
real, dimension(isc-4:iec+4,jsc-4:jec+4,nk) :: mass_mdppm
real, dimension(isc-4:iec+4,jsc-4:jec+4,nk) :: tracermass_mdppm

! fields for mdmdt 
real, dimension(isc-4:iec+4,jsc-4:jec+4,nk) :: tmask_mdmdt
real, dimension(isc-4:iec+4,jsc-4:jec+4,nk) :: tracer_mdmdt
real, dimension(isc-4:iec+4,jsc-4:jec+4,nk) :: mass_mdmdt
real, dimension(isc-4:iec+4,jsc-4:jec+4,nk) :: tracermass_mdmdt

! for doing all tracers with sweby 
type  :: tracer_mdfl_type
  real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: field
end type tracer_mdfl_type

#else

real, dimension(:,:,:), allocatable :: flux_x
real, dimension(:,:,:), allocatable :: flux_y
real, dimension(:,:,:), allocatable :: flux_z

real, dimension(:,:,:), allocatable :: advect_tendency
real, dimension(:,:,:), allocatable :: neutral_temp_advect
real, dimension(:,:,:), allocatable :: neutral_salt_advect

! some flow independent quantities for quicker 
real, dimension(:,:,:), allocatable :: quick_x
real, dimension(:,:,:), allocatable :: quick_y 
real, dimension(:,:),   allocatable :: quick_z 
real, dimension(:,:,:), allocatable :: curv_xp
real, dimension(:,:,:), allocatable :: curv_xn
real, dimension(:,:,:), allocatable :: curv_yp
real, dimension(:,:,:), allocatable :: curv_yn
real, dimension(:,:),   allocatable :: curv_zp
real, dimension(:,:),   allocatable :: curv_zn
real, dimension(:,:),   allocatable :: dxt_quick
real, dimension(:,:),   allocatable :: dyt_quick
real, dimension(:,:,:), allocatable :: tmask_quick
real, dimension(:,:,:), allocatable :: tracer_quick

! masks for other advection schemes 
real, dimension(:,:,:), allocatable :: tmask_fourth
real, dimension(:,:,:), allocatable :: tmask_sixth
real, dimension(:,:,:), allocatable :: tracer_fourth
real, dimension(:,:,:), allocatable :: tracer_sixth

real, dimension(:,:,:), allocatable :: tracer_mdfl
real, dimension(:,:,:), allocatable :: tmask_mdfl
real, dimension(:,:,:), allocatable :: mass_mdfl
real, dimension(:,:,:), allocatable :: tracermass_mdfl

real, dimension(:,:,:), allocatable :: tracer_mdppm
real, dimension(:,:,:), allocatable :: tmask_mdppm
real, dimension(:,:,:), allocatable :: mass_mdppm
real, dimension(:,:,:), allocatable :: tracermass_mdppm

real, dimension(:,:,:), allocatable :: tmask_mdmdt
real, dimension(:,:,:), allocatable :: tracer_mdmdt
real, dimension(:,:,:), allocatable :: mass_mdmdt
real, dimension(:,:,:), allocatable :: tracermass_mdmdt

! for doing all tracers with sweby 
type  :: tracer_mdfl_type
  real, dimension(:,:,:), pointer   :: field => NULL()
end type tracer_mdfl_type

#endif

! for psom, the mass of seawater in a tracer cell 
! and the tendency for use with diagnostics 
real, dimension(:,:,:),  allocatable :: sm

integer :: num_prog_tracers = 0
integer :: index_temp=-1
integer :: index_salt=-1
integer :: index_temp_sq=-1
integer :: index_salt_sq=-1

character(len=256) :: version='CVS $Id: ocean_tracer_advect.F90,v 20.0 2013/12/14 00:17:22 fms Exp $'
character(len=256) :: tagname='Tag $Name: tikal $'


type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type)  , pointer :: Grd =>NULL()
type(ocean_domain_type), save    :: Dom_quicker
type(ocean_domain_type), save    :: Dom_fourth
type(ocean_domain_type), save    :: Dom_sixth
type(ocean_domain_type), save    :: Dom_mdfl
type(ocean_domain_type), save    :: Dom_flux
type(ocean_domain_type), save    :: Dom_mdppm
type(ocean_domain_type), save    :: Dom_mdmdt

type(tracer_mdfl_type), dimension(:), allocatable :: tracer_mdfl_all  

real, parameter :: a4=7.0/12.0, b4=-1.0/12.0
real, parameter :: a6=37.0/60.0, b6=-2.0/15.0, c6=1.0/60.0

logical :: used

integer  :: id_neut_rho_advect
integer  :: id_neut_rho_advect_on_nrho
integer  :: id_wdian_rho_advect
integer  :: id_wdian_rho_advect_on_nrho
integer  :: id_tform_rho_advect
integer  :: id_tform_rho_advect_on_nrho

integer  :: id_neut_temp_advect
integer  :: id_neut_temp_advect_on_nrho
integer  :: id_wdian_temp_advect
integer  :: id_wdian_temp_advect_on_nrho
integer  :: id_tform_temp_advect
integer  :: id_tform_temp_advect_on_nrho

integer  :: id_neut_salt_advect
integer  :: id_neut_salt_advect_on_nrho
integer  :: id_wdian_salt_advect
integer  :: id_wdian_salt_advect_on_nrho
integer  :: id_tform_salt_advect
integer  :: id_tform_salt_advect_on_nrho

integer, dimension(:), allocatable :: id_tracer_advection
integer, dimension(:), allocatable :: id_tracer2_advection
integer, dimension(:), allocatable :: id_tracer_adv_diss

integer, dimension(:), allocatable :: id_sweby_advect
integer, dimension(:), allocatable :: id_psom_advect
integer, dimension(:), allocatable :: id_vert_advect
integer, dimension(:), allocatable :: id_horz_advect
integer, dimension(:), allocatable :: id_advection_x
integer, dimension(:), allocatable :: id_advection_y
integer, dimension(:), allocatable :: id_advection_z
integer, dimension(:), allocatable :: id_xflux_adv
integer, dimension(:), allocatable :: id_yflux_adv
integer, dimension(:), allocatable :: id_zflux_adv
integer, dimension(:), allocatable :: id_xflux_adv_int_z
integer, dimension(:), allocatable :: id_yflux_adv_int_z

! for gyre/overturning diagnostics 
integer                              :: id_basin_mask=-1
integer, dimension(:,:), allocatable :: id_merid_flux_advect          
integer, dimension(:,:), allocatable :: id_merid_flux_over
integer, dimension(:,:), allocatable :: id_merid_flux_gyre
real, dimension(:,:), allocatable    :: basin_mask     
real, dimension(:,:), allocatable    :: basin_mask_jp1 
real, dimension(:,:), allocatable    :: data           

real, dimension(:,:,:), allocatable :: global_basin_mask
real, dimension(:,:),   allocatable :: global_dxtn
real, dimension(:,:),   allocatable :: global_dxt
real, dimension(:,:,:), allocatable :: global_tmask
real, dimension(:,:,:), allocatable :: global_rho_dzt
real, dimension(:,:,:), allocatable :: global_tracer
real, dimension(:,:,:), allocatable :: global_vhrho_nt
real, dimension(:,:,:), allocatable :: global_flux_y

real, dimension(:,:,:), allocatable :: global_basin_mask_jp1
real, dimension(:,:),   allocatable :: global_dxt_jp1
real, dimension(:,:,:), allocatable :: global_tracer_jp1
real, dimension(:,:,:), allocatable :: global_tmask_jp1

public  horz_advect_tracer
private horz_advect_tracer_upwind
private horz_advect_tracer_2nd_order
private horz_advect_tracer_4th_order
private horz_advect_tracer_6th_order
private horz_advect_tracer_quicker
private horz_advect_tracer_quickmom3

public  vert_advect_tracer
private vert_advect_tracer_upwind
private vert_advect_tracer_2nd_order
private vert_advect_tracer_4th_order
private vert_advect_tracer_6th_order
private vert_advect_tracer_quicker
private vert_advect_tracer_quickmom3

private advect_tracer_mdfl_sup_b
private advect_tracer_mdfl_sweby
private advect_tracer_mdfl_sweby_test
private advect_tracer_sweby_all
private advect_tracer_mdppm
private advect_tracer_mdppm_test
private advect_tracer_mdmdt_test

private ppm_limit_ifc
private ppm_limit_cw84
private ppm_limit_sh

private advect_tracer_psom
private psom_x
private psom_y
private psom_z

private quicker_init
private mdfl_init
private fourth_sixth_init
private mdppm_init
private psom_init
private mdmdt_init

private advection_diag_init
private watermass_diag_init

private gyre_overturn_diagnose_init
private gyre_overturn_diagnose
private compute_adv_diss
private watermass_diag 

public ocean_tracer_advect_init
public ocean_tracer_advect_end
public ocean_tracer_advect_restart

logical :: module_is_initialized   = .false.
logical :: have_obc                = .false.
logical :: compute_watermass_diag  = .false. 

! for gyre/overturning diagnostics 
type(domain1d),save :: Domx,Domy
logical :: read_basin_mask                = .false.
logical :: compute_gyre_overturn_diagnose = .false.
integer, allocatable                   :: pelist_x(:)
integer                                :: layout(2)
real,    dimension(:,:,:), allocatable :: local_basin_mask
real,    dimension(:,:,:), allocatable :: local_basin_mask_jp1
logical, dimension(0:5)                :: do_this_basin

logical :: limit_with_upwind       = .false.
logical :: debug_this_module       = .false.
logical :: advect_sweby_all        = .false.
logical :: zero_tracer_advect_horz = .false.
logical :: zero_tracer_advect_vert = .false. 
logical :: write_a_restart         = .true. 
logical :: psom_limit_prather      = .false. 
logical :: async_domain_update     = .false. 

namelist /ocean_tracer_advect_nml/ debug_this_module, limit_with_upwind, advect_sweby_all,   &
                                   zero_tracer_advect_horz, zero_tracer_advect_vert,         &
                                   write_a_restart, psom_limit_prather, read_basin_mask,     &
                                   async_domain_update

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_advect_init">
!
! <DESCRIPTION>
! Initialize the tracer advection module.
! </DESCRIPTION>
!
subroutine ocean_tracer_advect_init (Grid, Domain, Time, Dens, T_prog, obc, debug)

  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_density_type),     intent(in)           :: Dens
  type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
  logical,                      intent(in)           :: obc
  logical,                      intent(in), optional :: debug
  
  integer :: n
  integer :: ioun, io_status, ierr
  
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  write( stdlogunit,'(/a/)') trim(version)

  module_is_initialized = .true.

  have_obc = obc

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_tracer_advect_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_tracer_advect_nml')
#else
  ioun = open_namelist_file()
  read  (ioun, ocean_tracer_advect_nml,iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_tracer_advect_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_tracer_advect_nml)  
  write (stdlogunit, ocean_tracer_advect_nml)

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk  = Grid%nk
  allocate(flux_x(isd:ied,jsd:jed,nk))
  allocate(flux_y(isd:ied,jsd:jed,nk))
  allocate(flux_z(isd:ied,jsd:jed,nk))  
  allocate(advect_tendency(isd:ied,jsd:jed,nk))  
  allocate(neutral_temp_advect(isd:ied,jsd:jed,nk))  
  allocate(neutral_salt_advect(isd:ied,jsd:jed,nk))  
#endif

  Dom => Domain
  Grd => Grid

  flux_x = 0.0
  flux_y = 0.0
  flux_z = 0.0

  advect_tendency(:,:,:)     = 0.0
  neutral_temp_advect(:,:,:) = 0.0
  neutral_salt_advect(:,:,:) = 0.0

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') &
    '==>Note: running ocean_tracer_advect_mod with debug_this_module=.true.'  
  endif 

  num_prog_tracers = size(T_prog(:))
  do n=1,num_prog_tracers
     if (trim(T_prog(n)%name) == 'temp') index_temp = n
     if (trim(T_prog(n)%name) == 'salt') index_salt = n
     if (trim(T_prog(n)%name) == 'passive_temp_sq') index_temp_sq = n
     if (trim(T_prog(n)%name) == 'passive_salt_sq') index_salt_sq = n
  enddo

  if(.not. write_a_restart) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_tracer_advect with write_a_restart=.false.'
    write(stdoutunit,'(a)') '   Will NOT write restart file, so cannot restart if using PSOM advection.'
  endif

  if(zero_tracer_advect_horz) then 
     call mpp_error(WARNING,'==>ocean_tracer_advect_mod: have turned OFF horizontal tracer advection')
     write(stdoutunit,*)'==>WARNING: have turned off horizontal tracer advection. Unrealistic simulation.'
  endif 
  if(zero_tracer_advect_vert) then 
     call mpp_error(WARNING,'==>ocean_tracer_advect_mod: have turned OFF vertical tracer advection')
     write(stdoutunit,*)'==>WARNING: have turned off vertical tracer advection. Unrealistic simulation.'
  endif 

  if(advect_sweby_all) then
      if(async_domain_update) then 
        write(stdoutunit,'(a)') &
        '==>Note: using asynchrnous domain update for MDFL SWEBY_all'  
      endif 

      write(stdoutunit,'(/a)')'==>ocean_tracer_advect_mod: advect_sweby_all=.true. so all tracers advected'
      write(stdoutunit,'(a)') '   with mdfl_sweby, regardless the settings in field_table.'
      write(stdoutunit,'(a/)')'   This method exploits mpp_update_domain capabilities and is faster in some cases.'

  else

      write(stdoutunit,'(a)') ' ' 
      write(stdoutunit,'(a)') ' From ocean_tracer_advect_init: SUMMARY OF TRACER ADVECTION SCHEMES' 

      do n=1,num_prog_tracers

         if(T_prog(n)%horz_advect_scheme==ADVECT_UPWIND) then 
             write(stdoutunit,'(1x,a)') &
             trim(T_prog(n)%name) //'is using first order upwind for horz advection.'
         endif
         if(T_prog(n)%vert_advect_scheme==ADVECT_UPWIND) then 
             write(stdoutunit,'(1x,a)') &
                  trim(T_prog(n)%name) //'is using first order upwind for vert advection.'
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_2ND_ORDER) then 
             write(stdoutunit,'(1x,a)') &
                  trim(T_prog(n)%name) //'is using 2nd order centred for horz advection.'
         endif
         if(T_prog(n)%vert_advect_scheme==ADVECT_2ND_ORDER) then 
             write(stdoutunit,'(1x,a)') &
                  trim(T_prog(n)%name) //'is using 2nd order centred for vert advection.'
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_4TH_ORDER) then 
             write(stdoutunit,'(1x,a)') &
                  trim(T_prog(n)%name) //' is using 4th order centred for horz advection.'
         endif
         if(T_prog(n)%vert_advect_scheme==ADVECT_4TH_ORDER) then 
             write(stdoutunit,'(1x,a)') &
                  trim(T_prog(n)%name) //'is using 4th order centred for vert advection.'
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_6TH_ORDER) then 
             write(stdoutunit,'(1x,a)') &
                  trim(T_prog(n)%name) //' is using 6th order centred for horz advection.'
         endif
         if(T_prog(n)%vert_advect_scheme==ADVECT_6TH_ORDER) then 
             write(stdoutunit,'(1x,a)') &
               trim(T_prog(n)%name) //' is using 6th order centred for vert advection.'
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_QUICKER) then 
             write(stdoutunit,'(1x,a)') &
                  trim(T_prog(n)%name) //' is using Quicker for horz advection.'
             if(limit_with_upwind) then 
                 write(stdoutunit,'(1x,a)')'limit_with_upwind reverts quicker to upwind if tracer outside limits'
             endif
         endif
         if(T_prog(n)%vert_advect_scheme==ADVECT_QUICKER) then 
             write(stdoutunit,'(1x,a)') &
                 trim(T_prog(n)%name) //' is using Quicker for vert advection.'
             if(limit_with_upwind) then 
                 write(stdoutunit,'(1x,a)')'limit_with_upwind reverts quicker to upwind if tracer outside limits'
             endif
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_QUICKMOM3) then 
             write(stdoutunit,'(1x,a)') &
                  trim(T_prog(n)%name) //' is using MOM3-Quicker for horz advection.'
         endif
         if(T_prog(n)%vert_advect_scheme==ADVECT_QUICKMOM3) then 
             write(stdoutunit,'(1x,a)') 'From ocean_tracer_advect_init: tracer ',trim(T_prog(n)%name), &
                  ' is using MOM3-Quicker for vertical advection.'
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_MDFL_SUP_B) then 
             write(stdoutunit,'(1x,a)') &
                  trim(T_prog(n)%name) //' is using MDFL Super-B for horz/vert advection.'
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_MDFL_SWEBY) then 
             write(stdoutunit,'(1x,a)') &
                  trim(T_prog(n)%name) //' is using MDFL Sweby for horz/vert advection.'
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_DST_LINEAR) then 
             write(stdoutunit,'(1x,a)') &
                  trim(T_prog(n)%name) //' is using linear direct space-time for horiz/vert advection.'
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_PSOM) then 
             write(stdoutunit,'(1x,a)') &
                    trim(T_prog(n)%name) //' is using Prather SOM for horz/vert advection.'
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_MDPPM) then 
             write(stdoutunit,'(1x,a)') &
                    trim(T_prog(n)%name) //' is using multi-dim piecewise parabolic for horz/vert advection.'
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_MDFL_SWEBY_TEST) then 
             write(stdoutunit,'(1x,a)') &
                  trim(T_prog(n)%name) //' is using unsupported test MDFL Sweby for horz/vert advection.'
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_DST_LINEAR_TEST) then 
             write(stdoutunit,'(1x,a)') &
                  trim(T_prog(n)%name) //' is using unsupported test linear direct space-time for horiz/vert advection.'
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_MDPPM_TEST) then 
             write(stdoutunit,'(1x,a)') &
                    trim(T_prog(n)%name) //' is using unsupported test multi-dim piecewise parabolic for horz/vert advection.'
         endif

         if(T_prog(n)%horz_advect_scheme==ADVECT_MDMDT_TEST) then 
             write(stdoutunit,'(1x,a)') &
                    trim(T_prog(n)%name) //' is using unsupported test multi-dim modified Daru & Tenaud for horz/vert advection.'
         endif

      enddo
      write(stdoutunit,'(a)') ' ' 

  endif

  call mpp_define_domains( (/1,Grid%ni,1,Grid%nj/), Domain%layout, Dom_flux%domain2d, maskmap=Domain%maskmap&
             , xflags = Domain%xflags, yflags = Domain%yflags, xhalo=1, yhalo=1,name='flux'&
             , x_cyclic_offset = Domain%x_cyclic_offset, y_cyclic_offset = Domain%y_cyclic_offset)


  ! for the gyre/overturning diagnostic
  call gyre_overturn_diagnose_init(Time, T_prog(:))

  ! initialize schemes other than 2nd order centered and first order upwind 
  call fourth_sixth_init
  call quicker_init
  call mdfl_init
  call mdppm_init
  call psom_init(Time, T_prog(:) )
  call mdmdt_init

  ! array needed for psom advection
  allocate( sm(isd:ied,jsd:jed,nk) )
  sm(:,:,:)            = 0.0

  call advection_diag_init(Time, T_prog(:))
  call watermass_diag_init(Time, Dens)

  ! initialize clock ids 
  id_clock_up_horz        = mpp_clock_id('(Ocean advect: horz up)          ',grain=CLOCK_ROUTINE)
  id_clock_up_vert        = mpp_clock_id('(Ocean advect: vert up)          ',grain=CLOCK_ROUTINE)
  id_clock_2nd_horz       = mpp_clock_id('(Ocean advect: horz 2nd)         ',grain=CLOCK_ROUTINE)
  id_clock_2nd_vert       = mpp_clock_id('(Ocean advect: vert 2nd)         ',grain=CLOCK_ROUTINE)
  id_clock_gyre_over      = mpp_clock_id('(Ocean advect: gyre_overturn)    ',grain=CLOCK_ROUTINE)
  id_clock_adv_diss       = mpp_clock_id('(Ocean advect: advect diss)      ',grain=CLOCK_ROUTINE)


end subroutine ocean_tracer_advect_init
! </SUBROUTINE> NAME="ocean_tracer_advect_init"


!#######################################################################
! <SUBROUTINE NAME="advection_diag_init">
!
! <DESCRIPTION>
! Initialize the main tracer advection diagnostics. 
! </DESCRIPTION>
!
subroutine advection_diag_init (Time, T_prog)

  type(ocean_time_type),        intent(in)  :: Time
  type(ocean_prog_tracer_type), intent(in)  :: T_prog(:)

  integer :: n  
  integer :: stdoutunit
  stdoutunit=stdout()

  allocate (id_tracer_advection(num_prog_tracers))
  allocate (id_tracer2_advection(num_prog_tracers))
  allocate (id_tracer_adv_diss(num_prog_tracers))
  allocate (id_sweby_advect(num_prog_tracers))
  allocate (id_psom_advect(num_prog_tracers))
  allocate (id_horz_advect(num_prog_tracers))
  allocate (id_vert_advect(num_prog_tracers))
  allocate (id_advection_x(num_prog_tracers))
  allocate (id_advection_y(num_prog_tracers))
  allocate (id_advection_z(num_prog_tracers))
  allocate (id_xflux_adv(num_prog_tracers))
  allocate (id_yflux_adv(num_prog_tracers))
  allocate (id_zflux_adv(num_prog_tracers))
  allocate (id_xflux_adv_int_z(num_prog_tracers))
  allocate (id_yflux_adv_int_z(num_prog_tracers))

  id_tracer_advection =-1
  id_tracer2_advection=-1
  id_tracer_adv_diss  =-1
  id_sweby_advect     =-1
  id_psom_advect      =-1
  id_horz_advect      =-1
  id_vert_advect      =-1
  id_advection_x      =-1
  id_advection_y      =-1
  id_advection_z      =-1
  id_xflux_adv        =-1
  id_yflux_adv        =-1
  id_zflux_adv        =-1
  id_xflux_adv_int_z  =-1
  id_yflux_adv_int_z  =-1

  do n=1,num_prog_tracers 

    if(n==index_temp) then 
      id_tracer_advection(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_advection', &
                   Grd%tracer_axes(1:3), Time%model_time, 'cp*rho*dzt*advection tendency',             &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_tracer2_advection(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sq_advection', &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*advection tendency for (cp*temp)^2',    &
                   'kg/(m^2*sec)*(J/kg)^2', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_sweby_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sweby_advect', &
                   Grd%tracer_axes(1:3), Time%model_time, 'cp*rho*dzt*sweby advect tendency',         &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_psom_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_psom_advect', &
                   Grd%tracer_axes(1:3), Time%model_time, 'cp*rho*dzt*psom advect tendency',        &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_horz_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_horz_advect', &
                   Grd%tracer_axes(1:3), Time%model_time, 'cp*rho*dzt*horz advect tendency',        &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_vert_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_vert_advect', &
                   Grd%tracer_axes(1:3), Time%model_time, 'cp*rho*dzt*vert advect tendency',        &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_advection_x(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_x_adv', &
                   Grd%tracer_axes(1:3), Time%model_time, 'cp*rho*dzt*x-advective heating',   &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))      
      id_advection_y(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_y_adv', &
                   Grd%tracer_axes(1:3), Time%model_time, 'cp*rho*dzt*y-advective heating',   &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))      
      id_advection_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_z_adv', &
                   Grd%tracer_axes(1:3), Time%model_time, 'cp*rho*dzt*z-advective heating',   &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))      
      id_xflux_adv(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_adv', &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time, 'cp*rho*dzt*dyt*u*temp',        &
                   'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_yflux_adv(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_adv', &
                   Grd%tracer_axes_flux_y(1:3), Time%model_time, 'cp*rho*dzt*dxt*v*temp',        &
                   'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_zflux_adv(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_zflux_adv', &
                   Grd%tracer_axes_wt(1:3), Time%model_time, 'cp*rho*dxt*dyt*wt*temp',           &
                   'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))      
      id_xflux_adv_int_z(n) = register_diag_field ('ocean_model', &
                   trim(T_prog(n)%name)//'_xflux_adv_int_z',      &
                   Grd%tracer_axes_flux_x(1:2), Time%model_time,   &
                   'z-integral of cp*rho*dyt*u*temp',              &
                   'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_yflux_adv_int_z(n) = register_diag_field ('ocean_model', &
                   trim(T_prog(n)%name)//'_yflux_adv_int_z',      &
                   Grd%tracer_axes_flux_y(1:2), Time%model_time,   &
                   'z-integral of cp*rho*dxt*v*temp',              &
                   'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_tracer_adv_diss(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_adv_diss', &
      Grd%tracer_axes(1:3), Time%model_time,                                                         &
      'dissipation of squared '//trim(T_prog(n)%name)//' from advection errors',                     &
      '(Watt/m^2)^2', missing_value=missing_value,                                                   &
       range=(/-1.e18,1.e18/))

    else

      id_tracer_advection(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_advection', &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*advection tendency',         &
                   'kg/(sec*m^2)', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_tracer2_advection(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sq_advection',    &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*advection tendency for tracer sqrd',&
                   'kg/(m^2*s)*(tracer squared)', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_sweby_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sweby_advect', &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*sweby advect tendency',     &
                   'kg/(sec*m^2)', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_psom_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_psom_advect',   &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*psom advect tendency',      &
                   'kg/(sec*m^2)', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_horz_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_horz_advect', &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*horz advect tendency',    &
                   'kg/(sec*m^2)', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_vert_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_vert_advect', &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*vert advect tendency',    &
                   'kg/(sec*m^2)', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_advection_x(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_x_adv', &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*x-advective tendency',     &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))            
      id_advection_y(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_y_adv', &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*y-advective tendency',     &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))            
      id_advection_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_z_adv', &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*z-advective tendency',     &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))            
      id_xflux_adv(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_adv', &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time, 'rho*dzt*dyt*u*tracer',         &
                   'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_yflux_adv(n)   = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_adv', &
                   Grd%tracer_axes_flux_y(1:3), Time%model_time, 'rho*dzt*dxt*v*tracer',           &
                   'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_zflux_adv(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_zflux_adv', &
                   Grd%tracer_axes_wt(1:3), Time%model_time, 'rho*dxt*dyt*wt*tracer',            &
                   'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))            
      id_xflux_adv_int_z(n) = register_diag_field ('ocean_model', &
                   trim(T_prog(n)%name)//'_xflux_adv_int_z',      &
                   Grd%tracer_axes_flux_x(1:2), Time%model_time,   &
                   'z-integral of rho*dzt*dyt*u*tracer',           &
                   'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_yflux_adv_int_z(n) = register_diag_field ('ocean_model', &
                   trim(T_prog(n)%name)//'_yflux_adv_int_z',      &
                   Grd%tracer_axes_flux_y(1:2), Time%model_time,   &
                   'z-integral of rho*dzt*dxt*v*tracer',           &
                   'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_tracer_adv_diss(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_adv_diss', &
      Grd%tracer_axes(1:3), Time%model_time,                                                         &
      'dissipation of squared '//trim(T_prog(n)%name)//' from advection errors',                     &
      '[kg/(m^2*sec)]^2', missing_value=missing_value,                                               &
       range=(/-1.e18,1.e18/))


    endif 

  enddo 


end subroutine advection_diag_init
! </SUBROUTINE> NAME="advection_diag_init"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag_init">
!
! <DESCRIPTION>
! Initialize the watermass diagnostics. 
! </DESCRIPTION>
!
subroutine watermass_diag_init (Time, Dens)
  type(ocean_time_type),    intent(in)  :: Time
  type(ocean_density_type), intent(in)  :: Dens

  integer :: stdoutunit
  stdoutunit=stdout()

  id_neut_rho_advect = register_diag_field ('ocean_model', 'neut_rho_advect',    &
     Grd%tracer_axes(1:3), Time%model_time, 'advection tendency for neutral rho',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_advect > 0) compute_watermass_diag = .true. 

  id_wdian_rho_advect = register_diag_field ('ocean_model', 'wdian_rho_advect',&
    Grd%tracer_axes(1:3), Time%model_time,                                     &
    'dianeutral mass transport from advection tendency for neutral rho',       &
    'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_advect > 0) compute_watermass_diag = .true. 

  id_tform_rho_advect = register_diag_field ('ocean_model',          &
   'tform_rho_advect', Grd%tracer_axes(1:3), Time%model_time,        &
   'water mass transformation from advection on level (pre-binning)',&
   'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_advect > 0) compute_watermass_diag = .true. 

  id_neut_rho_advect_on_nrho = register_diag_field ('ocean_model',            &
   'neut_rho_advect_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,     &
   'update of locally ref potrho from advection as binned to neutral density',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_rho_advect_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_advect_on_nrho = register_diag_field ('ocean_model',            &
    'wdian_rho_advect_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,    &
    'dianeutral mass transport due to advection as binned to neutral density', &
    'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_rho_advect_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_advect_on_nrho = register_diag_field ('ocean_model',          &
   'tform_rho_advect_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time, &
   'water mass transformation from advection as binned to neutral density',&
   'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_rho_advect_on_nrho > 0) compute_watermass_diag = .true. 
  

  ! temperature contributions 
  id_neut_temp_advect = register_diag_field ('ocean_model', 'neut_temp_advect',               &
     Grd%tracer_axes(1:3), Time%model_time, 'temp related advection tendency for neutral rho',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_temp_advect > 0) compute_watermass_diag = .true. 

  id_wdian_temp_advect = register_diag_field ('ocean_model', 'wdian_temp_advect',    &
    Grd%tracer_axes(1:3), Time%model_time,                                           &
    'temp related dianeutral mass transport from advection tendency for neutral rho',&
    'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_temp_advect > 0) compute_watermass_diag = .true. 

  id_tform_temp_advect = register_diag_field ('ocean_model',        &
   'tform_temp_advect', Grd%tracer_axes(1:3), Time%model_time,      &
   'temp related water mass transformation from advection on level',&
   'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_temp_advect > 0) compute_watermass_diag = .true. 

  id_neut_temp_advect_on_nrho = register_diag_field ('ocean_model',                     &
   'neut_temp_advect_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,              &
   'temp related update of locally ref potrho from advection binned to neutral density',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_temp_advect_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_temp_advect_on_nrho = register_diag_field ('ocean_model',                    &
    'wdian_temp_advect_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,            &
    'temp related dianeutral mass transport due to advection binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_temp_advect_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_temp_advect_on_nrho = register_diag_field ('ocean_model',                 &
   'tform_temp_advect_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,          &
   'temp related water mass transformation from advection binned to neutral density',&
   'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_temp_advect_on_nrho > 0) compute_watermass_diag = .true. 
  

  ! salinity contributions 
  id_neut_salt_advect = register_diag_field ('ocean_model', 'neut_salt_advect',               &
     Grd%tracer_axes(1:3), Time%model_time, 'salt related advection tendency for neutral rho',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_salt_advect > 0) compute_watermass_diag = .true. 

  id_wdian_salt_advect = register_diag_field ('ocean_model', 'wdian_salt_advect',    &
    Grd%tracer_axes(1:3), Time%model_time,                                           &
    'salt related dianeutral mass transport from advection tendency for neutral rho',&
    'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_salt_advect > 0) compute_watermass_diag = .true. 

  id_tform_salt_advect = register_diag_field ('ocean_model',        &
   'tform_salt_advect', Grd%tracer_axes(1:3), Time%model_time,      &
   'salt related water mass transformation from advection on level',&
   'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_salt_advect > 0) compute_watermass_diag = .true. 

  id_neut_salt_advect_on_nrho = register_diag_field ('ocean_model',                     &
   'neut_salt_advect_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,              &
   'salt related update of locally ref potrho from advection binned to neutral density',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_neut_salt_advect_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_salt_advect_on_nrho = register_diag_field ('ocean_model',                    &
    'wdian_salt_advect_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,            &
    'salt related dianeutral mass transport due to advection binned to neutral density',&
    'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_wdian_salt_advect_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_salt_advect_on_nrho = register_diag_field ('ocean_model',                 &
   'tform_salt_advect_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,          &
   'salt related water mass transformation from advection binned to neutral density',&
   'kg/sec', missing_value=missing_value, range=(/-1e20,1e20/))
  if(id_tform_salt_advect_on_nrho > 0) compute_watermass_diag = .true. 
  


  if(compute_watermass_diag) then 
    write(stdoutunit,'(/a/)') &
    '==>Note: running ocean_tracer_advect w/ compute_watermass_diag=.true.'  
  endif 

end subroutine watermass_diag_init
! </SUBROUTINE> NAME="watermass_diag_init"


!#######################################################################
! <SUBROUTINE NAME="gyre_overturn_diagnose_init">
!
! <DESCRIPTION>
! Initialize diagnostics and fields for gyre/overturning.
!
! March 2012 for faster approach:
! russell.fiedler@csiro.au
!
! Some reorganization
! Stephen.Griffies
!
! </DESCRIPTION>
!
subroutine gyre_overturn_diagnose_init(Time, T_prog)
  type(ocean_time_type),        intent(in)  :: Time
  type(ocean_prog_tracer_type), intent(in)  :: T_prog(:)

  integer :: i, j, n, nbasin

  allocate (id_merid_flux_advect(0:5,num_prog_tracers))        
  allocate (id_merid_flux_over(0:5,num_prog_tracers))        
  allocate (id_merid_flux_gyre(0:5,num_prog_tracers))        

  id_merid_flux_advect=-1
  id_merid_flux_over  =-1
  id_merid_flux_gyre  =-1

  ! register the diagnostic fields 
  do n=1,num_prog_tracers 

    if(n==index_temp) then 

      id_merid_flux_advect(0,n)=register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_advect_global',            &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'Global integral(cp*rho*dzt*dxt*v*temp)',                     &
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_over(0,n) = register_diag_field ('ocean_model',     &
      trim(T_prog(n)%name)//'_merid_flux_over_global',                  &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                     &
      'overturn contribution to global integral(cp*rho*dzt*dxt*v*temp)',&
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_gyre(0,n) = register_diag_field ('ocean_model',  &
      trim(T_prog(n)%name)//'_merid_flux_gyre_global',               &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                  &
      'gyre contribution to global integral(cp*rho*dzt*dxt*v*temp)', &
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_advect(1,n)=register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_advect_southern',          &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'ACC integral(cp*rho*dzt*dxt*v*temp)',                        &
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_over(1,n) = register_diag_field ('ocean_model',  &
      trim(T_prog(n)%name)//'_merid_flux_over_southern',             &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                  &
      'overturn contribution to ACC integral(cp*rho*dzt*dxt*v*temp)',&
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_gyre(1,n) = register_diag_field ('ocean_model',  &
      trim(T_prog(n)%name)//'_merid_flux_gyre_southern',             &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                  &
      'gyre contribution to ACC integral(cp*rho*dzt*dxt*v*temp)',    &
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_advect(2,n)=register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_advect_atlantic',          &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'Atlantic integral(cp*rho*dzt*dxt*v*temp)',                   &
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_over(2,n) = register_diag_field ('ocean_model',       &
      trim(T_prog(n)%name)//'_merid_flux_over_atlantic',                  &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                       &
      'overturn contribution to Atlantic integral(cp*rho*dzt*dxt*v*temp)',&
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_gyre(2,n) = register_diag_field ('ocean_model',   &
      trim(T_prog(n)%name)//'_merid_flux_gyre_atlantic',              &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                   &
      'gyre contribution to Atlantic integral(cp*rho*dzt*dxt*v*temp)',&
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_advect(3,n)=register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_advect_pacific',           &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'Pacific integral(cp*rho*dzt*dxt*v*temp)',                    &
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_over(3,n) = register_diag_field ('ocean_model',      &
      trim(T_prog(n)%name)//'_merid_flux_over_pacific',                  &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                      &
      'overturn contribution to Pacific integral(cp*rho*dzt*dxt*v*temp)',&
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_gyre(3,n) = register_diag_field ('ocean_model',  &
      trim(T_prog(n)%name)//'_merid_flux_gyre_pacific',              &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                  &
      'gyre contribution to Pacific integral(cp*rho*dzt*dxt*v*temp)',&
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_advect(4,n)=register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_advect_arctic',            &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'Arctic integral(cp*rho*dzt*dxt*v*temp)',                     &
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_over(4,n) = register_diag_field ('ocean_model',     &
      trim(T_prog(n)%name)//'_merid_flux_over_arctic',                  &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                     &
      'overturn contribution to Arctic integral(cp*rho*dzt*dxt*v*temp)',&
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_gyre(4,n) = register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_gyre_arctic',              &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'gyre contribution to Arctic integral(cp*rho*dzt*dxt*v*temp)',&
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_advect(5,n)=register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_advect_indian',            &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'Indian integral(cp*rho*dzt*dxt*v*temp)',                     &
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_over(5,n) = register_diag_field ('ocean_model',     &
      trim(T_prog(n)%name)//'_merid_flux_over_indian',                  &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                     &
      'overturn contribution to Indian integral(cp*rho*dzt*dxt*v*temp)',&
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_gyre(5,n) = register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_gyre_indian',              &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'gyre contribution to Indian integral(cp*rho*dzt*dxt*v*temp)',&
      'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

    else

      id_merid_flux_advect(0,n)=register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_advect_global',            &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'Global integral(rho*dzt*dxt*v*temp)',                        &
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_over(0,n) = register_diag_field ('ocean_model',  &
      trim(T_prog(n)%name)//'_merid_flux_over_global',               &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                  &
      'overturn contribution to global integral(rho*dzt*dxt*v*temp)',&
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_gyre(0,n) = register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_gyre_global',              &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'gyre contribution to global integral(rho*dzt*dxt*v*temp)',   &
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_advect(1,n)=register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_advect_southern',          &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'ACC integral(rho*dzt*dxt*v*temp)',                           &
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_over(1,n) = register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_over_southern',            &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'overturn contribution to ACC integral(rho*dzt*dxt*v*temp)',  &
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_gyre(1,n) = register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_gyre_southern',            &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'gyre contribution to ACC integral(rho*dzt*dxt*v*temp)',      &
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_advect(2,n)=register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_advect_atlantic',          &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'Atlantic integral(rho*dzt*dxt*v*temp)',                      &
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_over(2,n) = register_diag_field ('ocean_model',    &
      trim(T_prog(n)%name)//'_merid_flux_over_atlantic',               &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                    &
      'overturn contribution to Atlantic integral(rho*dzt*dxt*v*temp)',&
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_gyre(2,n) = register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_gyre_atlantic',            &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'gyre contribution to Atlantic integral(rho*dzt*dxt*v*temp)', &
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_advect(3,n)=register_diag_field ('ocean_model', &
      trim(T_prog(n)%name)//'_merid_flux_advect_pacific',           &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                 &
      'Pacific integral(rho*dzt*dxt*v*temp)',                       &
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_over(3,n) = register_diag_field ('ocean_model',   &
      trim(T_prog(n)%name)//'_merid_flux_over_pacific',               &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                   &
      'overturn contribution to Pacific integral(rho*dzt*dxt*v*temp)',&
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_gyre(3,n) = register_diag_field ('ocean_model',&
      trim(T_prog(n)%name)//'_merid_flux_gyre_pacific',            &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                &
      'gyre contribution to Pacific integral(rho*dzt*dxt*v*temp)', &
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_advect(4,n)=register_diag_field ('ocean_model',&
      trim(T_prog(n)%name)//'_merid_flux_advect_arctic',           &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                &
      'Arctic integral(rho*dzt*dxt*v*temp)',                       &
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_over(4,n) = register_diag_field ('ocean_model',  &
      trim(T_prog(n)%name)//'_merid_flux_over_arctic',               &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                  &
      'overturn contribution to Arctic integral(rho*dzt*dxt*v*temp)',&
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_gyre(4,n) = register_diag_field ('ocean_model',&
      trim(T_prog(n)%name)//'_merid_flux_gyre_arctic',             &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                &
      'gyre contribution to Arctic integral(rho*dzt*dxt*v*temp)',  &
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_advect(5,n)=register_diag_field ('ocean_model',&
      trim(T_prog(n)%name)//'_merid_flux_advect_indian',           &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                &
      'Indian integral(rho*dzt*dxt*v*temp)',                       &
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_over(5,n) = register_diag_field ('ocean_model',  &
      trim(T_prog(n)%name)//'_merid_flux_over_indian',               &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                  &
      'overturn contribution to Indian integral(rho*dzt*dxt*v*temp)',&
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_merid_flux_gyre(5,n) = register_diag_field ('ocean_model',&
      trim(T_prog(n)%name)//'_merid_flux_gyre_indian',             &
      Grd%tracer_axes_flux_y(2:2), Time%model_time,                &
      'gyre contribution to Indian integral(rho*dzt*dxt*v*temp)',  &
      'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

    endif ! endif for index_temp 

    ! check to see if any of the fields have been requested in diagnostic table 
    do nbasin=0,5
      if(id_merid_flux_advect(nbasin,n) > 0) compute_gyre_overturn_diagnose = .true.
      if(id_merid_flux_over(nbasin,n)   > 0) compute_gyre_overturn_diagnose = .true.
      if(id_merid_flux_gyre(nbasin,n)   > 0) compute_gyre_overturn_diagnose = .true.
    enddo

  enddo   ! enddo for num_prog_tracers 

  ! return if not needing to compute the diagnostic 
  if(.not. compute_gyre_overturn_diagnose) then
      return
  endif



  allocate(basin_mask(isd:ied,jsd:jed))
  allocate(data(isd:ied,jsd:jed))
  basin_mask(:,:) = Grd%tmask(:,:,1)
  data(:,:)       = 0.0

  allocate(basin_mask_jp1(isd:ied,jsd:jed))
  do j=jsc,jec
     do i=isc,iec
        basin_mask_jp1(i,j) = Grd%tmask(i,j+1,1)
     enddo
  enddo
  call mpp_update_domains(basin_mask_jp1(:,:), Dom%domain2d)


  ! define basin mask, with default basin_mask(:,:)=Grd%tmask(:,:,1). 
  if(read_basin_mask) then 
      call read_data('INPUT/basin_mask','basin_mask', data, Dom%domain2d)

      ! use max function to eliminate missing values < 0.0
      do j=jsc,jec
         do i=isc,iec
            basin_mask(i,j)     = max(0.0,data(i,j)*Grd%tmask(i,j,1))
            basin_mask_jp1(i,j) = max(0.0,data(i,j+1)*Grd%tmask(i,j+1,1))
         enddo
      enddo
      call mpp_update_domains(basin_mask(:,:)    , Dom%domain2d)
      call mpp_update_domains(basin_mask_jp1(:,:), Dom%domain2d)
  endif

  ! This approach only requires those PEs in the horizontal to perform sums. 
  ! Partial sums are computed locally.

      ! get pelist
      call mpp_get_layout (Dom%domain2d, layout)
      allocate ( pelist_x(layout(1)))
      call mpp_get_domain_components (Dom%domain2d, Domx,Domy)
      call mpp_get_pelist( Domx, pelist_x )
      allocate(local_basin_mask(isd:ied,jsd:jed,0:5))    ; local_basin_mask     = 0.0
      allocate(local_basin_mask_jp1(isd:ied,jsd:jed,0:5)); local_basin_mask_jp1 = 0.0
      if(read_basin_mask) then 
          do j=jsd,jed
             do i=isc,iec
                if(basin_mask(i,j)==1.0) local_basin_mask(i,j,1) = 1.0
                if(basin_mask(i,j)==2.0) local_basin_mask(i,j,2) = 1.0
                if(basin_mask(i,j)==3.0) local_basin_mask(i,j,3) = 1.0
                if(basin_mask(i,j)==4.0) local_basin_mask(i,j,4) = 1.0
                if(basin_mask(i,j)==5.0) local_basin_mask(i,j,5) = 1.0

                if(basin_mask_jp1(i,j)==1.0) local_basin_mask_jp1(i,j,1) = 1.0
                if(basin_mask_jp1(i,j)==2.0) local_basin_mask_jp1(i,j,2) = 1.0
                if(basin_mask_jp1(i,j)==3.0) local_basin_mask_jp1(i,j,3) = 1.0
                if(basin_mask_jp1(i,j)==4.0) local_basin_mask_jp1(i,j,4) = 1.0
                if(basin_mask_jp1(i,j)==5.0) local_basin_mask_jp1(i,j,5) = 1.0
             enddo
          enddo
          ! fill basin=0 mask with tmask
          do nbasin=0,0
             local_basin_mask(:,:,nbasin)           = Grd%tmask(:,:,1)
             local_basin_mask_jp1(:,jsc:jec,nbasin) = Grd%tmask(:,jsc+1:jec+1,1)
          enddo

      else 
      ! fill basin mask with tmask as default  
          do nbasin=0,5
             local_basin_mask(isc:iec,jsc:jec,nbasin) = Grd%tmask(isc:iec,jsc:jec,1)
             local_basin_mask_jp1(:,jsc:jec,nbasin)   = Grd%tmask(:,jsc+1:jec+1,1)
          enddo

      endif
      call mpp_update_domains(local_basin_mask(:,:,:), Dom%domain2d)
      call mpp_update_domains(local_basin_mask_jp1(:,:,:), Dom%domain2d)


      ! basin masks do not change with time.  
      ! we precalculate a mask to decide if want to do
      ! anything on this PE.
      do nbasin=0,5
         if(any(local_basin_mask(:,:,nbasin)==1.0) )  then
             do_this_basin(nbasin) = .true.
         else
             do_this_basin(nbasin) = .false.
         endif
      enddo

end subroutine gyre_overturn_diagnose_init
! </SUBROUTINE> NAME="gyre_overturn_diagnose_init



!#######################################################################
! <SUBROUTINE NAME="quicker_init">
!
! <DESCRIPTION>
! Initialize quicker specific fields. 
! </DESCRIPTION>
!
subroutine quicker_init

  integer :: i, j, k, kp1, km1, kp2

#ifndef MOM_STATIC_ARRAYS
  allocate(quick_x(isd:ied,jsd:jed,2))
  allocate(quick_y(isd:ied,jsd:jed, 2))
  allocate(quick_z(nk,2))
  allocate(curv_xp(isd:ied,jsd:jed,3))
  allocate(curv_xn(isd:ied,jsd:jed, 3))
  allocate(curv_yp(isd:ied,jsd:jed,3))
  allocate(curv_yn(isd:ied,jsd:jed, 3))
  allocate(curv_zp(nk,3))
  allocate(curv_zn(nk,3))
  allocate(dxt_quick(isc-2:iec+2,jsc-2:jec+2))
  allocate(dyt_quick(isc-2:iec+2,jsc-2:jec+2))
  allocate(tmask_quick(isc-2:iec+2,jsc-2:jec+2,nk))    
  allocate(tracer_quick(isc-2:iec+2,jsc-2:jec+2,nk))
#endif
 
  quick_x       = 0.0
  quick_y       = 0.0
  quick_z       = 0.0
  curv_xp       = 0.0
  curv_xn       = 0.0
  curv_yp       = 0.0
  curv_yn       = 0.0
  curv_zp       = 0.0
  curv_zn       = 0.0
  dxt_quick     = 0.0
  dyt_quick     = 0.0
  tmask_quick   = 0.0
  tracer_quick  = 0.0

  call set_ocean_domain(Dom_quicker,Grd,xhalo=2,yhalo=2,name='quicker',maskmap=Dom%maskmap)

  do i=isc,iec
     do j=jsc,jec
        dxt_quick(i,j) = Grd%dxt(i,j)
        dyt_quick(i,j) = Grd%dyt(i,j)
        do k=1,nk
           tmask_quick(i,j,k)  = Grd%tmask(i,j,k)
        enddo
     enddo
  enddo
  call mpp_update_domains(tmask_quick,Dom_quicker%domain2d)  


  ! fill in boundaries (needed for non-cyclic case to prevent blow-ups)
  do i=isc-2,isc-1
     dxt_quick(i,jsc:jec) = dxt_quick(isc,jsc:jec)
     dyt_quick(i,:) = dyt_quick(isc,:)
  enddo

  do i=iec+1,iec+2
     dxt_quick(i,jsc:jec) = dxt_quick(iec,jsc:jec)
     dyt_quick(i,:) = dyt_quick(iec,:)
  enddo

  do j=jsc-2,jsc-1
     dxt_quick(:,j) = dxt_quick(:,jsc)
     dyt_quick(:,j) = dyt_quick(:,jsc)
  enddo

  do j=jec+1,jec+2
     dxt_quick(:,j) = dxt_quick(:,jec)
     dyt_quick(:,j) = dyt_quick(:,jec)
  enddo
  call mpp_update_domains(dxt_quick,Dom_quicker%domain2d)
  call mpp_update_domains(dyt_quick,Dom_quicker%domain2d)

  
! calculate quicker weights on computational domain (not filled in halos)
  
  do i=isc-1,iec
     do j=jsc-1,jec
        quick_x(i,j,1) = dxt_quick(i+1,j)/(dxt_quick(i+1,j)+dxt_quick(i,j))
        quick_x(i,j,2) = dxt_quick(i,j)/(dxt_quick(i+1,j)+dxt_quick(i,j))
        quick_y(i,j,1) = dyt_quick(i,j+1)/(dyt_quick(i,j+1)+dyt_quick(i,j))
        quick_y(i,j,2) = dyt_quick(i,j)/(dyt_quick(i,j+1)+dyt_quick(i,j))
        
        curv_xp(i,j,1) = dxt_quick(i,j)*dxt_quick(i+1,j)/&
             ((dxt_quick(i-1,j)+2.0*dxt_quick(i,j)+dxt_quick(i+1,j))&
             *(dxt_quick(i,j)+dxt_quick(i+1,j)))
        curv_xp(i,j,2) = -dxt_quick(i,j)*dxt_quick(i+1,j)/&
             ((dxt_quick(i,j)+dxt_quick(i+1,j))&
             *(dxt_quick(i-1,j)+dxt_quick(i,j)))        
        curv_xp(i,j,3) = dxt_quick(i,j)*dxt_quick(i+1,j)/&
             ((dxt_quick(i-1,j)+2.0*dxt_quick(i,j)+dxt_quick(i+1,j))&
             *(dxt_quick(i-1,j)+dxt_quick(i,j)))

        curv_xn(i,j,1) = dxt_quick(i,j)*dxt_quick(i+1,j)/&
             ((dxt_quick(i,j)+2.0*dxt_quick(i+1,j)+dxt_quick(i+2,j))&
             *(dxt_quick(i+1,j)+dxt_quick(i+2,j)))
        curv_xn(i,j,2) = -dxt_quick(i,j)*dxt_quick(i+1,j)/&
             ((dxt_quick(i+1,j)+dxt_quick(i+2,j))&
             *(dxt_quick(i,j)+dxt_quick(i+1,j)))        
        curv_xn(i,j,3) = dxt_quick(i,j)*dxt_quick(i+1,j)/&
             ((dxt_quick(i,j)+2.0*dxt_quick(i+1,j)+dxt_quick(i+2,j))&
             *(dxt_quick(i,j)+dxt_quick(i+1,j)))

        curv_yp(i,j,1) = dyt_quick(i,j)*dyt_quick(i,j+1)/&
             ((dyt_quick(i,j-1)+2.0*dyt_quick(i,j)+dyt_quick(i,j+1))&
             *(dyt_quick(i,j)+dyt_quick(i,j+1)))
        curv_yp(i,j,2) = -dyt_quick(i,j)*dyt_quick(i,j+1)/&
             ((dyt_quick(i,j)+dyt_quick(i,j+1))&
             *(dyt_quick(i,j-1)+dyt_quick(i,j)))        
        curv_yp(i,j,3) = dyt_quick(i,j)*dyt_quick(i,j+1)/&
             ((dyt_quick(i,j-1)+2.0*dyt_quick(i,j)+dyt_quick(i,j+1))&
             *(dyt_quick(i,j-1)+dyt_quick(i,j)))

        curv_yn(i,j,1) = dyt_quick(i,j)*dyt_quick(i,j+1)/&
             ((dyt_quick(i,j)+2.0*dyt_quick(i,j+1)+dyt_quick(i,j+2))&
             *(dyt_quick(i,j+1)+dyt_quick(i,j+2)))
        curv_yn(i,j,2) = -dyt_quick(i,j)*dyt_quick(i,j+1)/&
             ((dyt_quick(i,j+1)+dyt_quick(i,j+2))&
             *(dyt_quick(i,j)+dyt_quick(i,j+1)))        
        curv_yn(i,j,3) = dyt_quick(i,j)*dyt_quick(i,j+1)/&
             ((dyt_quick(i,j)+2.0*dyt_quick(i,j+1)+dyt_quick(i,j+2))&
             *(dyt_quick(i,j)+dyt_quick(i,j+1)))
        
     enddo
  enddo

  do k= 1,nk
     kp2 = min(k+2,nk)
     kp1 = min(k+1,nk)
     km1 = max(k-1,1)
     quick_z(k,1) = Grd%dzt(kp1)/(Grd%dzt(kp1)+Grd%dzt(k))
     quick_z(k,2) = Grd%dzt(k  )/(Grd%dzt(kp1)+Grd%dzt(k))
     curv_zp(k,1) = Grd%dzt(k)*Grd%dzt(kp1)/((Grd%dzt(km1)+2.0*Grd%dzt(k)+Grd%dzt(kp1))*(Grd%dzt(k)+Grd%dzt(kp1)))
     curv_zp(k,2) =-Grd%dzt(k)*Grd%dzt(kp1)/((Grd%dzt(k)+Grd%dzt(kp1))*(Grd%dzt(km1)+Grd%dzt(k)))
     curv_zp(k,3) = Grd%dzt(k)*Grd%dzt(kp1)/((Grd%dzt(km1)+2.0*Grd%dzt(k)+Grd%dzt(kp1))*(Grd%dzt(km1)+Grd%dzt(k)))
     curv_zn(k,1) = Grd%dzt(k)*Grd%dzt(kp1)/((Grd%dzt(k)+2.0*Grd%dzt(kp1)+Grd%dzt(kp2))*(Grd%dzt(kp1)+Grd%dzt(kp2)))
     curv_zn(k,2) =-Grd%dzt(k)*Grd%dzt(kp1)/((Grd%dzt(kp1)+Grd%dzt(kp2))*(Grd%dzt(k)+Grd%dzt(kp1)))
     curv_zn(k,3) = Grd%dzt(k)*Grd%dzt(kp1)/((Grd%dzt(k)+2.0*Grd%dzt(kp1)+Grd%dzt(kp2))*(Grd%dzt(k)+Grd%dzt(kp1)))
  enddo

  ! initialize clock ids 
  id_clock_quick_horz     = mpp_clock_id('(Ocean advect: horz quk)     ',grain=CLOCK_ROUTINE)
  id_clock_quick_vert     = mpp_clock_id('(Ocean advect: vert quk)     ',grain=CLOCK_ROUTINE)
  id_clock_quickmom3_horz = mpp_clock_id('(Ocean advect: horz qukmom3) ',grain=CLOCK_ROUTINE)
  id_clock_quickmom3_vert = mpp_clock_id('(Ocean advect: vert qukmom3) ',grain=CLOCK_ROUTINE)

end subroutine quicker_init
! </SUBROUTINE> NAME="quicker_init"


!#######################################################################
! <SUBROUTINE NAME="fourth_sixth_init">
!
! <DESCRIPTION>
! Initialize the fourth order and sixth order advection fields.
! </DESCRIPTION>
!
subroutine fourth_sixth_init

  integer :: i, j, k

#ifndef MOM_STATIC_ARRAYS
  allocate(tmask_fourth (isc-2:iec+2,jsc-2:jec+2,nk))    
  allocate(tracer_fourth(isc-2:iec+2,jsc-2:jec+2,nk))    
  allocate(tmask_sixth  (isc-3:iec+3,jsc-3:jec+3,nk))    
  allocate(tracer_sixth (isc-3:iec+3,jsc-3:jec+3,nk))    
#endif
  tmask_fourth  = 0.0
  tracer_fourth = 0.0
  tmask_sixth   = 0.0
  tracer_sixth  = 0.0
  
  call set_ocean_domain(Dom_fourth,Grd,xhalo=2,yhalo=2,name='fourth',maskmap=Dom%maskmap)
  call set_ocean_domain(Dom_sixth ,Grd,xhalo=3,yhalo=3,name='sixth' ,maskmap=Dom%maskmap)

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tmask_fourth(i,j,k) = Grd%tmask(i,j,k)
           tmask_sixth(i,j,k)  = Grd%tmask(i,j,k)
        enddo
     enddo
  enddo
  call mpp_update_domains(tmask_fourth,Dom_fourth%domain2d)  
  call mpp_update_domains(tmask_sixth ,Dom_sixth%domain2d)  

  ! initialize clock ids 
  id_clock_4th_horz       = mpp_clock_id('(Ocean advect: horz 4th)     ',grain=CLOCK_ROUTINE)
  id_clock_4th_vert       = mpp_clock_id('(Ocean advect: vert 4th)     ',grain=CLOCK_ROUTINE)
  id_clock_6th_horz       = mpp_clock_id('(Ocean advect: horz 6th)     ',grain=CLOCK_ROUTINE)
  id_clock_6th_vert       = mpp_clock_id('(Ocean advect: vert 6th)     ',grain=CLOCK_ROUTINE)


end subroutine fourth_sixth_init
! </SUBROUTINE> NAME="fourth_sixth_init"


!#######################################################################
! <SUBROUTINE NAME="mdfl_init">
!
! <DESCRIPTION>
! Initialize mdfl and dst_linear specific fields.
! </DESCRIPTION>
!
subroutine mdfl_init 

  integer :: i, j, k, n 

  call set_ocean_domain(Dom_mdfl,Grd,xhalo=2,yhalo=2,name='mdfl',maskmap=Dom%maskmap)

  allocate(tracer_mdfl_all(num_prog_tracers))
#ifndef MOM_STATIC_ARRAYS
  allocate(tmask_mdfl (isc-2:iec+2,jsc-2:jec+2,nk))    
  allocate(tracer_mdfl(isc-2:iec+2,jsc-2:jec+2,nk))    
  allocate(mass_mdfl(isc-2:iec+2,jsc-2:jec+2,nk))    
  allocate(tracermass_mdfl(isc-2:iec+2,jsc-2:jec+2,nk))    
  do n=1,num_prog_tracers
     allocate (tracer_mdfl_all(n)%field(isc-2:iec+2,jsc-2:jec+2,nk))
  enddo
#endif
  tmask_mdfl  = 0.0
  tracer_mdfl = 0.0 
  mass_mdfl = 0.0 
  tracermass_mdfl = 0.0 
  do n=1,num_prog_tracers
     tracer_mdfl_all(n)%field(:,:,:) = 0.0
  enddo

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tmask_mdfl(i,j,k) = Grd%tmask(i,j,k)
        enddo
     enddo
  enddo
  call mpp_update_domains(tmask_mdfl,Dom_mdfl%domain2d)  

  ! initialize clock ids 
  id_clock_mdfl_sup_b       = mpp_clock_id('(Ocean advect: MDFL-sup-b)       ',grain=CLOCK_ROUTINE)
  id_clock_mdfl_sweby       = mpp_clock_id('(Ocean advect: MDFL-sweby)       ',grain=CLOCK_ROUTINE)
  id_clock_mdfl_sweby_all   = mpp_clock_id('(Ocean advect: MDFL-sweby-all)   ',grain=CLOCK_ROUTINE)
  id_clock_mdfl_sweby_test  = mpp_clock_id('(Ocean advect: MDFL-sweby-test)  ',grain=CLOCK_ROUTINE)
  id_clock_dst_linear       = mpp_clock_id('(Ocean advect: DST-linear)       ',grain=CLOCK_ROUTINE)
  id_clock_dst_linear_test  = mpp_clock_id('(Ocean advect: DST-linear-test)  ',grain=CLOCK_ROUTINE)
  id_clock_mdfl_sweby_mpi = mpp_clock_id('(Ocean advect: MDFL-sweby-mpi)   ',grain=CLOCK_ROUTINE)
  id_clock_mdfl_sweby_cu1 = mpp_clock_id('(Ocean advect: MDFL-sweby-cuk)   ',grain=CLOCK_ROUTINE)
  id_clock_mdfl_sweby_cu2 = mpp_clock_id('(Ocean advect: MDFL-sweby-cui)   ',grain=CLOCK_ROUTINE)
  id_clock_mdfl_sweby_cu3 = mpp_clock_id('(Ocean advect: MDFL-sweby-cuj)   ',grain=CLOCK_ROUTINE)
  id_clock_mdfl_sweby_dia = mpp_clock_id('(Ocean advect: MDFL-sweby-diag)  ',grain=CLOCK_ROUTINE)


end subroutine mdfl_init
! </SUBROUTINE> NAME="mdfl_init"


!#######################################################################
! <SUBROUTINE NAME="mdppm_init">
!
! <DESCRIPTION>
! Initialize mdppm specific fields.
! </DESCRIPTION>
!
subroutine mdppm_init

  integer :: i, j, k

  call set_ocean_domain(Dom_mdppm,Grd,xhalo=4,yhalo=4,name='mdppm',maskmap=Dom%maskmap)         
                                                              
#ifndef MOM_STATIC_ARRAYS                                         
  allocate(tmask_mdppm (isc-4:iec+4,jsc-4:jec+4,nk))          
  allocate(tracer_mdppm(isc-4:iec+4,jsc-4:jec+4,nk))          
  allocate(tracermass_mdppm(isc-4:iec+4,jsc-4:jec+4,nk))          
  allocate(mass_mdppm(isc-4:iec+4,jsc-4:jec+4,nk))          
#endif                                                        
  tmask_mdppm      = 0.0                                           
  tracer_mdppm     = 0.0                                          
  tracermass_mdppm = 0.0                                          
  mass_mdppm       = 0.0                                          
                                                              
  do k=1,nk                                                   
     do j=jsc,jec                                             
        do i=isc,iec                                          
           tmask_mdppm(i,j,k) = Grd%tmask(i,j,k)               
        enddo                                                 
     enddo                                                    
  enddo                                                       
  call mpp_update_domains(tmask_mdppm,Dom_mdppm%domain2d)       
                                                              
  ! initialize clock ids                                      
  id_clock_mdppm      = mpp_clock_id('(Ocean advect: MDPPM)      ',grain=CLOCK_ROUTINE)
  id_clock_mdppm_test = mpp_clock_id('(Ocean advect: MDPPM-TEST) ',grain=CLOCK_ROUTINE)
                                             
end subroutine mdppm_init                                     
! </SUBROUTINE> NAME="mdppm_init"


!#######################################################################
! <SUBROUTINE NAME="mdmdt_init">
!
! <DESCRIPTION>
! Initialize mdmdt specific fields.
! </DESCRIPTION>
!
subroutine mdmdt_init

  integer :: i, j, k

  call set_ocean_domain(Dom_mdmdt,Grd,xhalo=4,yhalo=4,name='mdmdt',maskmap=Dom%maskmap)         
                                                              
#ifndef MOM_STATIC_ARRAYS                                         
  allocate(tmask_mdmdt (isc-4:iec+4,jsc-4:jec+4,nk))          
  allocate(tracer_mdmdt(isc-4:iec+4,jsc-4:jec+4,nk))          
  allocate(tracermass_mdmdt(isc-4:iec+4,jsc-4:jec+4,nk))          
  allocate(mass_mdmdt(isc-4:iec+4,jsc-4:jec+4,nk))          
#endif                                                        
  tmask_mdmdt      = 0.0                                           
  tracer_mdmdt     = 0.0                                          
  tracermass_mdmdt = 0.0                                          
  mass_mdmdt       = 0.0                                          
                                                              
  do k=1,nk                                                   
     do j=jsc,jec                                             
        do i=isc,iec                                          
           tmask_mdmdt(i,j,k) = Grd%tmask(i,j,k)               
        enddo                                                 
     enddo                                                    
  enddo                                                       
  call mpp_update_domains(tmask_mdmdt,Dom_mdmdt%domain2d)       
                                                              
  ! initialize clock ids                                      
  id_clock_mdmdt_test = mpp_clock_id('(Ocean advect: MDMDT-TEST) ',grain=CLOCK_ROUTINE)
                                             
end subroutine mdmdt_init                                     
! </SUBROUTINE> NAME="mdmdt_init"


!#######################################################################
! <SUBROUTINE NAME="psom_init">
!
! <DESCRIPTION>
! Read restart or initialize moments for prather advection
! </DESCRIPTION>
!
subroutine psom_init(Time, T_prog)

  type(ocean_time_type),           intent(in)    :: Time
  type(ocean_prog_tracer_type),    intent(inout) :: T_prog(:)

  logical :: mandatory  
  integer :: n, id_restart
  character(len=128) :: filename

  integer :: stdoutunit 
  stdoutunit=stdout() 

  filename = 'ocean_psom_moments.res.nc'

  do n=1,num_prog_tracers

     if (T_prog(n)%horz_advect_scheme == ADVECT_PSOM) then

        allocate( T_prog(n)%s0(isd:ied,jsd:jed,nk) ) 
        allocate( T_prog(n)%sx(isd:ied,jsd:jed,nk) )
        allocate( T_prog(n)%sxx(isd:ied,jsd:jed,nk) )
        allocate( T_prog(n)%sy(isd:ied,jsd:jed,nk) )
        allocate( T_prog(n)%syy(isd:ied,jsd:jed,nk) )
        allocate( T_prog(n)%sz(isd:ied,jsd:jed,nk) )
        allocate( T_prog(n)%szz(isd:ied,jsd:jed,nk) )
        allocate( T_prog(n)%sxy(isd:ied,jsd:jed,nk) )
        allocate( T_prog(n)%sxz(isd:ied,jsd:jed,nk) )
        allocate( T_prog(n)%syz(isd:ied,jsd:jed,nk) )

        T_prog(n)%s0(:,:,:)  = 0.0
        T_prog(n)%sx(:,:,:)  = 0.0
        T_prog(n)%sxx(:,:,:) = 0.0
        T_prog(n)%sy(:,:,:)  = 0.0
        T_prog(n)%syy(:,:,:) = 0.0
        T_prog(n)%sz(:,:,:)  = 0.0
        T_prog(n)%szz(:,:,:) = 0.0
        T_prog(n)%sxy(:,:,:) = 0.0
        T_prog(n)%sxz(:,:,:) = 0.0
        T_prog(n)%syz(:,:,:) = 0.0
        mandatory = .false.
        if(field_exist('INPUT/'//trim(filename), trim(T_prog(n)%name)//'_prather_s0')) mandatory = .true.
        id_restart = register_restart_field(Adv_restart, filename, trim(T_prog(n)%name)//'_prather_s0', &
             T_prog(n)%s0(:,:,:), Dom%domain2d, mandatory = mandatory)
        id_restart = register_restart_field(Adv_restart, filename, trim(T_prog(n)%name)//'_prather_sx', &
             T_prog(n)%sx(:,:,:), Dom%domain2d, mandatory = mandatory)
        id_restart = register_restart_field(Adv_restart, filename, trim(T_prog(n)%name)//'_prather_sxx', &
             T_prog(n)%sxx(:,:,:), Dom%domain2d, mandatory = mandatory)
        id_restart = register_restart_field(Adv_restart, filename, trim(T_prog(n)%name)//'_prather_sy', &
             T_prog(n)%sy(:,:,:), Dom%domain2d, mandatory = mandatory)
        id_restart = register_restart_field(Adv_restart, filename, trim(T_prog(n)%name)//'_prather_syy', &
             T_prog(n)%syy(:,:,:), Dom%domain2d, mandatory = mandatory)
        id_restart = register_restart_field(Adv_restart, filename, trim(T_prog(n)%name)//'_prather_sz', &
             T_prog(n)%sz(:,:,:), Dom%domain2d, mandatory = mandatory)
        id_restart = register_restart_field(Adv_restart, filename, trim(T_prog(n)%name)//'_prather_szz', &
             T_prog(n)%szz(:,:,:), Dom%domain2d, mandatory = mandatory)
        id_restart = register_restart_field(Adv_restart, filename, trim(T_prog(n)%name)//'_prather_sxy', &
             T_prog(n)%sxy(:,:,:), Dom%domain2d, mandatory = mandatory)
        id_restart = register_restart_field(Adv_restart, filename, trim(T_prog(n)%name)//'_prather_sxz', &
             T_prog(n)%sxz(:,:,:), Dom%domain2d, mandatory = mandatory)
        id_restart = register_restart_field(Adv_restart, filename, trim(T_prog(n)%name)//'_prather_syz', &
             T_prog(n)%syz(:,:,:), Dom%domain2d, mandatory = mandatory)
     endif
  enddo

  filename = 'INPUT/ocean_psom_moments.res.nc'
  do n=1,num_prog_tracers
     if(ANY(T_prog(1:num_prog_tracers)%horz_advect_scheme == ADVECT_PSOM)) then 
        if(file_exist(filename)) call restore_state(Adv_restart)
     endif
  enddo

  do n=1,num_prog_tracers

     if (T_prog(n)%horz_advect_scheme == ADVECT_PSOM) then   
         if (.not. field_exist(trim(filename), trim(T_prog(n)%name)//'_prather_s0')) then
             write (stdoutunit,*) 'Initializing psom moments for tracer ',trim(T_prog(n)%name)
         else
             write (stdoutunit,*) 'Read psom moments for tracer ',trim(T_prog(n)%name), ' from file INPUT/',filename
             write (stdoutunit,*) 'Completed read of psom restart for tracer ', trim(T_prog(n)%name)
             call mpp_update_domains(T_prog(n)%s0(:,:,:),  Dom%domain2d, complete = .false.)
             call mpp_update_domains(T_prog(n)%sx(:,:,:),  Dom%domain2d, complete = .false.)
             call mpp_update_domains(T_prog(n)%sxx(:,:,:), Dom%domain2d, complete = .false.)
             call mpp_update_domains(T_prog(n)%sy(:,:,:),  Dom%domain2d, complete = .false.)
             call mpp_update_domains(T_prog(n)%syy(:,:,:), Dom%domain2d, complete = .false.)
             call mpp_update_domains(T_prog(n)%sz(:,:,:),  Dom%domain2d, complete = .false.)
             call mpp_update_domains(T_prog(n)%szz(:,:,:), Dom%domain2d, complete = .false.)
             call mpp_update_domains(T_prog(n)%sxy(:,:,:), Dom%domain2d, complete = .false.)
             call mpp_update_domains(T_prog(n)%sxz(:,:,:), Dom%domain2d, complete = .false.)
             call mpp_update_domains(T_prog(n)%syz(:,:,:), Dom%domain2d, complete = .true.)

         endif

         call tracer_psom_chksum(Time, T_prog(n))

     endif ! endif for ADVECT_PSOM

  enddo    ! enddo for num_prog_tracers 

  id_clock_psom    = mpp_clock_id('(Ocean advect: psom total)   ',grain=CLOCK_ROUTINE)
  id_clock_psom_x  = mpp_clock_id('(Ocean psom advect: psom_x)  ',grain=CLOCK_ROUTINE)
  id_clock_psom_y  = mpp_clock_id('(Ocean psom advect: psom_y)  ',grain=CLOCK_ROUTINE)
  id_clock_psom_z  = mpp_clock_id('(Ocean psom advect: psom_z)  ',grain=CLOCK_ROUTINE)


end subroutine psom_init
! </SUBROUTINE> NAME="psom_init"


!#######################################################################
! <SUBROUTINE NAME="horz_advect_tracer">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers 
! </DESCRIPTION>
!
subroutine horz_advect_tracer(Time, Adv_vel, Thickness, Dens, T_prog, Tracer, ntracer, dtime, store_flux)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_prog_tracer_type), intent(inout) :: Tracer
  integer,                      intent(in)    :: ntracer
  real,                         intent(in)    :: dtime
  logical,            optional, intent(in)    :: store_flux

  real,dimension(isd:ied,jsd:jed)             :: tmp_flux
  integer                                     :: i, j, k
  integer                                     :: taum1, tau
  logical                                     :: store

  if(zero_tracer_advect_horz) return 

  store = .TRUE.
  if(present(store_flux)) store = store_flux

  taum1 = Time%taum1
  tau   = Time%tau

  if(.not. advect_sweby_all) then 

      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Tracer%wrk1(i,j,k) = 0.0
            enddo
         enddo
      enddo

      select case (Tracer%horz_advect_scheme)

      case (ADVECT_UPWIND)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  &
          -horz_advect_tracer_upwind(Adv_vel, Tracer%field(:,:,:,taum1))
      case (ADVECT_2ND_ORDER)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  &
          -horz_advect_tracer_2nd_order(Adv_vel, Tracer%field(:,:,:,tau))
      case (ADVECT_4TH_ORDER)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  &
          -horz_advect_tracer_4th_order(Adv_vel, Tracer%field(:,:,:,tau))
      case (ADVECT_6TH_ORDER)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  &
          -horz_advect_tracer_6th_order(Adv_vel, Tracer%field(:,:,:,tau))
      case (ADVECT_QUICKER)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  &
          -horz_advect_tracer_quicker(Adv_vel, Tracer, Tracer%field(:,:,:,taum1), Tracer%field(:,:,:,tau))
      case (ADVECT_QUICKMOM3)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  &
          -horz_advect_tracer_quickmom3(Adv_vel, Tracer,Tracer%field(:,:,:,taum1), Tracer%field(:,:,:,tau))

      case (ADVECT_MDFL_SUP_B)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  &
          -advect_tracer_mdfl_sup_b(Time, Adv_vel, Thickness, Tracer%field(:,:,:,taum1), dtime)
      case (ADVECT_MDFL_SWEBY)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  &
          -advect_tracer_mdfl_sweby(Time, Adv_vel, Thickness, Tracer%field(:,:,:,taum1), dtime, sweby_limiter=1.0)
      case (ADVECT_DST_LINEAR)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  &
          -advect_tracer_mdfl_sweby(Time, Adv_vel, Thickness, Tracer%field(:,:,:,taum1), dtime, sweby_limiter=0.0)
      case (ADVECT_PSOM)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  &
          -advect_tracer_psom(Time, Adv_vel, Tracer, Thickness, dtime, ntracer, do_passive_sq=.false.)
      case (ADVECT_MDPPM)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  &
          -advect_tracer_mdppm(Time, Adv_vel, Tracer, Thickness, Tracer%field(:,:,:,taum1), dtime)

      case (ADVECT_MDFL_SWEBY_TEST)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  &
          -advect_tracer_mdfl_sweby_test(Time, Adv_vel, Thickness, Tracer%field(:,:,:,taum1), dtime, sweby_limiter=1.0)
      case (ADVECT_DST_LINEAR_TEST)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  &
          -advect_tracer_mdfl_sweby_test(Time, Adv_vel, Thickness, Tracer%field(:,:,:,taum1), dtime, sweby_limiter=0.0)
      case (ADVECT_MDPPM_TEST)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  &
          -advect_tracer_mdppm_test(Time, Adv_vel, Tracer, Thickness, Tracer%field(:,:,:,taum1), dtime)
      case (ADVECT_MDMDT_TEST)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  &
          -advect_tracer_mdmdt_test(Time, Adv_vel, Tracer, Thickness, Tracer%field(:,:,:,taum1), dtime)

      case default
          call mpp_error(FATAL,&
          '==>Error from ocean_tracer_advect_mod (horz_advect_tracer): chose invalid horz advection scheme')
      end select
      
      if (have_obc) call ocean_obc_zero_boundary(Tracer%wrk1, "T")

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               Tracer%th_tendency(i,j,k) = Tracer%th_tendency(i,j,k) + Tracer%wrk1(i,j,k)
            enddo
         enddo
      enddo


      ! fill some diagnostic fields       
      advect_tendency(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               advect_tendency(i,j,k) = Tracer%wrk1(i,j,k)
            enddo
         enddo
      enddo

      if(ntracer==index_temp) then 
          neutral_temp_advect(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   neutral_temp_advect(i,j,k) = Tracer%wrk1(i,j,k)
                enddo
             enddo
          enddo
      endif
      if(ntracer==index_salt) then 
          neutral_salt_advect(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   neutral_salt_advect(i,j,k) = Tracer%wrk1(i,j,k)
                enddo
             enddo
          enddo
      endif


      ! send some diagnostics 

      if (id_horz_advect(ntracer) > 0 .and. Tracer%horz_advect_scheme /= ADVECT_PSOM) then 
         call diagnose_3d(Time, Grd, id_horz_advect(ntracer), Tracer%conversion*Tracer%wrk1(:,:,:))
      endif 
      if (id_psom_advect(ntracer) > 0 .and. Tracer%horz_advect_scheme == ADVECT_PSOM) then 
         call diagnose_3d(Time, Grd, id_psom_advect(ntracer), Tracer%conversion*Tracer%wrk1(:,:,:))
      endif 
      if (id_xflux_adv(ntracer) > 0) then 
         call diagnose_3d(Time, Grd, id_xflux_adv(ntracer), Tracer%conversion*flux_x(:,:,:))
      endif 
      if (id_yflux_adv(ntracer) > 0) then 
         call diagnose_3d(Time, Grd, id_yflux_adv(ntracer), Tracer%conversion*flux_y(:,:,:))
      endif 

      if (id_xflux_adv_int_z(ntracer) > 0) then 
          tmp_flux(:,:) = 0.0
          do k=1,nk
             tmp_flux(isc:iec,jsc:jec) = tmp_flux(isc:iec,jsc:jec) +  flux_x(isc:iec,jsc:jec,k) 
          enddo
          call diagnose_2d(Time, Grd, id_xflux_adv_int_z(ntracer), Tracer%conversion*tmp_flux(:,:))
      endif

      if (id_yflux_adv_int_z(ntracer) > 0) then 
          tmp_flux(:,:) = 0.0
          do k=1,nk
             tmp_flux(isc:iec,jsc:jec) = tmp_flux(isc:iec,jsc:jec) +  flux_y(isc:iec,jsc:jec,k) 
          enddo
          call diagnose_2d(Time, Grd, id_yflux_adv_int_z(ntracer), Tracer%conversion*tmp_flux(:,:))
      endif

      if (have_obc) then
        call store_ocean_obc_tracer_flux(Time,Tracer,flux_x,ntracer,'z','adv')
        call store_ocean_obc_tracer_flux(Time,Tracer,flux_y,ntracer,'m','adv')
      endif

      ! This is only needed for diagnostics of the flux through open boundaries.
      ! It must be checked, if still correct with the new scheme
      ! for gyre/overturning diagnostics
      ! Note we need to call this from within advect_tracer_sweby_all
      if(compute_gyre_overturn_diagnose) then
         call gyre_overturn_diagnose(Time, Adv_vel, Tracer, ntracer)
      endif

  endif   ! endif for (.not. advect_sweby_all)


  if(advect_sweby_all .and. ntracer==1) then 
    call advect_tracer_sweby_all(Time, Adv_vel, Dens, T_prog, Thickness, dtime)
  endif 

 
end subroutine horz_advect_tracer
! </SUBROUTINE> NAME="horz_advect_tracer"


!#######################################################################
! <SUBROUTINE NAME="vert_advect_tracer">
!
! <DESCRIPTION>
! Compute vertical advection of tracers for case 
! with advect_sweby_all=.false.
! </DESCRIPTION>
!
subroutine vert_advect_tracer(Time, Adv_vel, Dens, Thickness, T_prog, Tracer, ntracer, dtime)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_prog_tracer_type), intent(inout) :: Tracer
  integer,                      intent(in)    :: ntracer
  real,                         intent(in)    :: dtime

  integer :: tau, taum1
  integer :: i,j,k

  if(zero_tracer_advect_vert) return 

  tau   = Time%tau
  taum1 = Time%taum1

  if(.not. advect_sweby_all) then 

      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Tracer%wrk1(i,j,k) = 0.0
            enddo
         enddo
      enddo

      select case (Tracer%vert_advect_scheme)

      case (ADVECT_UPWIND)
          Tracer%wrk1(isc:iec,jsc:jec,:) = &
          -vert_advect_tracer_upwind(Adv_vel, Tracer%field(:,:,:,taum1))
      case (ADVECT_2ND_ORDER)
          Tracer%wrk1(isc:iec,jsc:jec,:) = &
          -vert_advect_tracer_2nd_order(Adv_vel, Tracer%field(:,:,:,tau))
      case (ADVECT_4TH_ORDER)
          Tracer%wrk1(isc:iec,jsc:jec,:) = &
          -vert_advect_tracer_4th_order(Adv_vel, Tracer%field(:,:,:,tau))
      case (ADVECT_6TH_ORDER)
          Tracer%wrk1(isc:iec,jsc:jec,:) = &
          -vert_advect_tracer_6th_order(Adv_vel, Tracer%field(:,:,:,tau))
      case (ADVECT_QUICKER)
          Tracer%wrk1(isc:iec,jsc:jec,:) = &
          -vert_advect_tracer_quicker(Adv_vel, Tracer, Tracer%field(:,:,:,taum1), Tracer%field(:,:,:,tau))
      case (ADVECT_QUICKMOM3)
          Tracer%wrk1(isc:iec,jsc:jec,:) = &
          -vert_advect_tracer_quickmom3(Adv_vel, Tracer, Tracer%field(:,:,:,taum1), Tracer%field(:,:,:,tau))

      ! the mdfl, mdppm, dst_linear, mdmdt, and PSOM schemes are three-dimensional, 
      ! with 3-dimensional advection tendency computed in horz_advect_tracer
      case (ADVECT_MDFL_SUP_B)
      case (ADVECT_MDFL_SWEBY)
      case (ADVECT_DST_LINEAR)
      case (ADVECT_MDPPM)
      case (ADVECT_PSOM)
      case (ADVECT_MDFL_SWEBY_TEST)
      case (ADVECT_DST_LINEAR_TEST)
      case (ADVECT_MDPPM_TEST)
      case (ADVECT_MDMDT_TEST)

      case default
        call mpp_error(FATAL,&
        '==>Error from ocean_tracer_advect_mod (vert_advect_tracer): invalid advection scheme chosen')
      end select

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               Tracer%th_tendency(i,j,k) = Tracer%th_tendency(i,j,k) + Tracer%wrk1(i,j,k)
            enddo
         enddo
      enddo

      ! for diagnostics, here and in compute_adv_diss 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               advect_tendency(i,j,k) = advect_tendency(i,j,k) + Tracer%wrk1(i,j,k)
            enddo
         enddo
      enddo

      if(id_tracer_advection(ntracer) > 0) then 
         call diagnose_3d(Time, Grd, id_tracer_advection(ntracer), Tracer%conversion*advect_tendency(:,:,:))
      endif

      ! for watermass_diag 
          if(ntracer==index_temp) then 
              do k=1,nk
                 do j=jsc,jec
                    do i=isc,iec
                       neutral_temp_advect(i,j,k) = neutral_temp_advect(i,j,k) + Tracer%wrk1(i,j,k)
                        enddo
                     enddo
                  enddo
              endif
          if(ntracer==index_salt) then 
              do k=1,nk
                 do j=jsc,jec
                    do i=isc,iec
                       neutral_salt_advect(i,j,k) = neutral_salt_advect(i,j,k) + Tracer%wrk1(i,j,k)
                           enddo
                        enddo
                     enddo
                 endif

      ! send some more diagnostics 
      if (id_vert_advect(ntracer) > 0) then 
         call diagnose_3d(Time, Grd, id_vert_advect(ntracer), Tracer%conversion*Tracer%wrk1(:,:,:))
      endif 
      if (id_zflux_adv(ntracer) > 0)  then 
         call diagnose_3d(Time, Grd, id_zflux_adv(ntracer), Tracer%conversion*flux_z(:,:,:))
      endif 

      if(ntracer == num_prog_tracers) then 
         call watermass_diag(Time, Dens)
      endif 

  endif  ! endif for (.not. advect_sweby_all)


  if(id_tracer_adv_diss(ntracer) > 0) then
     call compute_adv_diss(Time, Adv_vel, Thickness, T_prog(:), Tracer, ntracer, dtime)
  endif 

     
end subroutine vert_advect_tracer
! </SUBROUTINE> NAME="vert_advect_tracer"


!#######################################################################
! <FUNCTION NAME="horz_advect_tracer_upwind">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers from first order upwind.  
! This scheme is positive definite but very diffusive.  
! </DESCRIPTION>
!
function horz_advect_tracer_upwind(Adv_vel, Tracer_field)

  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field 

  real,dimension(isc:iec,jsc:jec,nk) :: horz_advect_tracer_upwind
  
  real, dimension(isd:ied,jsd:jed) :: fe, fn 
  real                             :: velocity, upos, uneg
  integer                          :: i, j, k

  call mpp_clock_begin(id_clock_up_horz)

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
     '==>Error from ocean_tracer_advect (horz_advect_tracer_upwind): ocean_tracer_advect_mod not yet initialized')
  endif 

  do k=1,nk 

     ! i-flux
     do j=jsc,jec
        do i=isc-1,iec
           velocity = 0.5*Adv_vel%uhrho_et(i,j,k)
           upos     = velocity + abs(velocity)
           uneg     = velocity - abs(velocity)
           fe(i,j)  = Grd%dyte(i,j)*(upos*Tracer_field(i,j,k) + uneg*Tracer_field(i+1,j,k)) &
                      *Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)  
           flux_x(i,j,k) = fe(i,j)
        enddo
     enddo

     ! j-flux
     do j=jsc-1,jec
        do i=isc,iec
           velocity = 0.5*Adv_vel%vhrho_nt(i,j,k)
           upos     = velocity + abs(velocity)
           uneg     = velocity - abs(velocity)
           fn(i,j)  = Grd%dxtn(i,j)*(upos*Tracer_field(i,j,k) + uneg*Tracer_field(i,j+1,k)) &
                      *Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
           flux_y(i,j,k) = fn(i,j)
        enddo
     enddo

     do j=jsc,jec
        do i=isc,iec
           horz_advect_tracer_upwind(i,j,k) = &
            Grd%tmask(i,j,k)*(fe(i,j)-fe(i-1,j)+fn(i,j)-fn(i,j-1))*Grd%datr(i,j)
        enddo
     enddo

  enddo

  call mpp_clock_end(id_clock_up_horz)
  
  
end function horz_advect_tracer_upwind
! </FUNCTION> NAME="horz_advect_tracer_upwind"


!#######################################################################
! <FUNCTION NAME="horz_advect_tracer_2nd_order">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers from second order 
! centered differences.  
! </DESCRIPTION>
!
function horz_advect_tracer_2nd_order(Adv_vel, Tracer_field)

  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field 

  real,dimension(isc:iec,jsc:jec,nk) :: horz_advect_tracer_2nd_order
  
  real, dimension(isd:ied) :: fs, fn
  integer                  :: i, j, k
  real                     :: fe, fw

  call mpp_clock_begin(id_clock_2nd_horz)

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
  '==>Error from ocean_tracer_advect (horz_advect_tracer_2nd_order): ocean_tracer_advect_mod not initialized')
  endif 
  
  fs(:) = 0.0; fn(:) = 0.0
  do k=1,nk 
    j     = jsc-1
    fs(:) = Grd%dxtn(:,j)*Adv_vel%vhrho_nt(:,j,k)*0.5*(Tracer_field(:,j+1,k)+Tracer_field(:,j,k))
    do j=jsc,jec
      i = isc-1
      fw = Grd%dyte(i,j)*Adv_vel%uhrho_et(i,j,k)*0.5*(Tracer_field(i+1,j,k)+Tracer_field(i,j,k))
      do i=isc,iec
        fe    = Grd%dyte(i,j)*Adv_vel%uhrho_et(i,j,k)*0.5*(Tracer_field(i+1,j,k)+Tracer_field(i,j,k))
        fn(i) = Grd%dxtn(i,j)*Adv_vel%vhrho_nt(i,j,k)*0.5*(Tracer_field(i,j+1,k)+Tracer_field(i,j,k))
        flux_x(i,j,k) = fe
        flux_y(i,j,k) = fn(i)
        horz_advect_tracer_2nd_order(i,j,k) = Grd%tmask(i,j,k)*(fe - fw + fn(i) - fs(i))*Grd%datr(i,j)
        fw = fe
      enddo
      fs(:) = fn(:)
    enddo
  enddo

  call mpp_clock_end(id_clock_2nd_horz)


end function horz_advect_tracer_2nd_order
! </FUNCTION> NAME="horz_advect_tracer_2nd_order"


!#######################################################################
! <FUNCTION NAME="horz_advect_tracer_4th_order">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers from fourth order centered
! differences. 
!
! WARNING: This code does NOT account for non-uniform grids.   
!
! </DESCRIPTION>
!
function horz_advect_tracer_4th_order(Adv_vel, Tracer_field)

  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field 

  real,dimension(isc:iec,jsc:jec,nk) :: horz_advect_tracer_4th_order

  integer :: i, j, k
  integer :: im1, ip2, jm1, jp2

  call mpp_clock_begin(id_clock_4th_horz)

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
     '==>Error from ocean_tracer_advect (horz_advect_tracer_4th_order): ocean_tracer_advect_mod not yet initialized')
  endif 

  tracer_fourth = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tracer_fourth(i,j,k) = Tracer_field(i,j,k)
        enddo
     enddo
  enddo
  call mpp_update_domains (tracer_fourth, Dom_fourth%domain2d)

  flux_x = 0.0
  flux_y = 0.0

  do k=1,nk

     do j=jsc,jec
        do i=isc-1,iec
           im1   = int(tmask_fourth(i-1,j,k)*(i-1) + (1.0-tmask_fourth(i-1,j,k))*i)
           ip2   = int(tmask_fourth(i+2,j,k)*(i+2) + (1.0-tmask_fourth(i+2,j,k))*(i+1))
           flux_x(i,j,k) = Grd%dyte(i,j)*Adv_vel%uhrho_et(i,j,k)*  &
                (a4*(tracer_fourth(i+1,j,k)+tracer_fourth(i,j,k)) + &
                 b4*(tracer_fourth(ip2,j,k)+tracer_fourth(im1,j,k)))
        enddo
     enddo

     do j=jsc-1,jec
        do i=isc,iec
           jm1   = int(tmask_fourth(i,j-1,k)*(j-1) + (1.0-tmask_fourth(i,j-1,k))*j)
           jp2   = int(tmask_fourth(i,j+2,k)*(j+2) + (1.0-tmask_fourth(i,j+2,k))*(j+1))
           flux_y(i,j,k) = Grd%dxtn(i,j)*Adv_vel%vhrho_nt(i,j,k)*&
                (a4*(tracer_fourth(i,j+1,k)+tracer_fourth(i,j,k)) + &
                 b4*(tracer_fourth(i,jp2,k)+tracer_fourth(i,jm1,k)))
        enddo
     enddo

  enddo

  ! to handle redundancies across the bipolar fold 
  if (Grd%tripolar) call mpp_update_domains(flux_x, flux_y, Dom_flux%domain2d, gridtype=CGRID_NE) 

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
          horz_advect_tracer_4th_order(i,j,k)  &
          = Grd%tmask(i,j,k)*(flux_x(i,j,k) - flux_x(i-1,j,k) + flux_y(i,j,k) - flux_y(i,j-1,k))*Grd%datr(i,j)
        enddo
     enddo
  enddo

  call mpp_clock_end(id_clock_4th_horz)

end function horz_advect_tracer_4th_order
! </FUNCTION> NAME="horz_advect_tracer_4th_order"


!#######################################################################
! <FUNCTION NAME="horz_advect_tracer_6th_order">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers from sixth order centered
! differences. 
!
! WARNING: this code does NOT account for non-uniform grids.   
!
! </DESCRIPTION>
!
function horz_advect_tracer_6th_order(Adv_vel, Tracer_field)

  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field 

  real,dimension(isc:iec,jsc:jec,nk) :: horz_advect_tracer_6th_order

  integer :: i, j, k
  integer :: im1, ip2, jm1, jp2

  call mpp_clock_begin(id_clock_6th_horz)

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
     '==>Error from ocean_tracer_advect (horz_advect_tracer_6th_order): ocean_tracer_advect_mod not yet initialized')
  endif 

  tracer_sixth=0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tracer_sixth(i,j,k) = Tracer_field(i,j,k)
        enddo
     enddo
  enddo
  call mpp_update_domains (tracer_sixth, Dom_sixth%domain2d)

  flux_x = 0.0
  flux_y = 0.0
  do k=1,nk

     do j=jsc,jec
        do i=isc-1,iec
           if (tmask_sixth(i-2,j,k)*tmask_sixth(i-1,j,k) == 0 .or. &
               tmask_sixth(i+3,j,k)*tmask_sixth(i+2,j,k) == 0) then
               im1   = int(tmask_sixth(i-1,j,k)*(i-1) + (1.0-tmask_sixth(i-1,j,k))*i)
               ip2   = int(tmask_sixth(i+2,j,k)*(i+2) + (1.0-tmask_sixth(i+2,j,k))*(i+1))
               flux_x(i,j,k) = Grd%dyte(i,j)*Adv_vel%uhrho_et(i,j,k)   &
                    *(  a4*(tracer_sixth(i+1,j,k)+tracer_sixth(i,j,k)) &
                      + b4*(tracer_sixth(ip2,j,k)+tracer_sixth(im1,j,k)))
           else
               flux_x(i,j,k) = Grd%dyte(i,j)*Adv_vel%uhrho_et(i,j,k)   &
                    *( a6*(tracer_sixth(i+1,j,k)+tracer_sixth(i,j,k))  &
                     + b6*(tracer_sixth(i+2,j,k)+tracer_sixth(i-1,j,k))&
                     + c6*(tracer_sixth(i+3,j,k)+tracer_sixth(i-2,j,k)))
           endif
        enddo
     enddo

     do j=jsc-1,jec
        do i=isc,iec
           if (tmask_sixth(i,j-2,k)*tmask_sixth(i,j-1,k) == 0 .or. &
               tmask_sixth(i,j+3,k)*tmask_sixth(i,j+2,k) == 0) then
               jm1   = int(tmask_sixth(i,j-1,k)*(j-1) + (1.0-tmask_sixth(i,j-1,k))*j)
               jp2   = int(tmask_sixth(i,j+2,k)*(j+2) + (1.0-tmask_sixth(i,j+2,k))*(j+1))
               flux_y(i,j,k) = Grd%dxtn(i,j)*Adv_vel%vhrho_nt(i,j,k) &
               *(  a4*(tracer_sixth(i,j+1,k)+tracer_sixth(i,j,k))    &
                 + b4*(tracer_sixth(i,jp2,k)+tracer_sixth(i,jm1,k)))
           else
               flux_y(i,j,k) = Grd%dxtn(i,j)*Adv_vel%vhrho_nt(i,j,k) &
               *(  a6*(tracer_sixth(i,j+1,k)+tracer_sixth(i,j,k))    &
                 + b6*(tracer_sixth(i,j+2,k)+tracer_sixth(i,j-1,k))  &
                 + c6*(tracer_sixth(i,j+3,k)+tracer_sixth(i,j-2,k)))
           endif
        enddo
     enddo

  enddo

  ! to handle redundancies across the bipolar fold 
  if (Grd%tripolar) call mpp_update_domains(flux_x, flux_y, Dom_flux%domain2d, gridtype=CGRID_NE) 

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
          horz_advect_tracer_6th_order(i,j,k) &
          = Grd%tmask(i,j,k)*(flux_x(i,j,k) - flux_x(i-1,j,k) + flux_y(i,j,k) - flux_y(i,j-1,k))*Grd%datr(i,j)
        enddo
     enddo
  enddo

  call mpp_clock_end(id_clock_6th_horz)

end function horz_advect_tracer_6th_order
! </FUNCTION> NAME="horz_advect_tracer_6th_order"


!#######################################################################
! <FUNCTION NAME="horz_advect_tracer_quicker">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers from quicker. 
! </DESCRIPTION>
!
function horz_advect_tracer_quicker(Adv_vel, Tracer, Tracer_field_taum1, Tracer_field_tau)

  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field_taum1 
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field_tau

  real,dimension(isc:iec,jsc:jec,nk) :: horz_advect_tracer_quicker

  integer :: i, j, k
  real    :: rnormsk,soutmsk, eastmsk, westmsk
  real    :: upos, uneg, vel

  call mpp_clock_begin(id_clock_quick_horz)

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
     '==>Error from ocean_tracer_advect (horz_advect_tracer_quicker): ocean_tracer_advect_mod not yet initialized')
  endif 

  tracer_quick = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tracer_quick(i,j,k) = Tracer_field_taum1(i,j,k)
        enddo
     enddo
  enddo
  call mpp_update_domains (tracer_quick, Dom_quicker%domain2d)

  flux_x = 0.0
  flux_y = 0.0
  do k=1,nk

     do j=jsc-1,jec
        do i=isc,iec
           vel  = Grd%dxtn(i,j)*Adv_vel%vhrho_nt(i,j,k)
           upos = 0.5*(vel + abs(vel))
           uneg = 0.5*(vel - abs(vel))
           rnormsk = tmask_quick(i,j,k)*(1-tmask_quick(i,j-1,k))
           soutmsk = tmask_quick(i,j+1,k)*(1-tmask_quick(i,j+2,k))
           flux_y(i,j,k) = &
               vel*(quick_y(i,j,1)*Tracer_field_tau(i,j,k)+quick_y(i,j,2)*Tracer_field_tau(i,j+1,k))&
                - upos*(curv_yp(i,j,1)*tracer_quick(i,j+1,k) + &
                curv_yp(i,j,2)*tracer_quick(i,j,k)   + &
                curv_yp(i,j,3)*(tracer_quick(i,j-1,k)*(1-rnormsk) + tracer_quick(i,j,k)*rnormsk)) &
                - uneg*(curv_yn(i,j,1)*(tracer_quick(i,j+2,k)*(1-soutmsk)+tracer_quick(i,j+1,k)*soutmsk) + &
                curv_yn(i,j,2)*tracer_quick(i,j+1,k)   + &
                curv_yn(i,j,3)*tracer_quick(i,j,k))        
        enddo
     enddo

     do j=jsc,jec
        do i=isc-1,iec
           vel  = Grd%dyte(i,j)*Adv_vel%uhrho_et(i,j,k)
           upos = 0.5*(vel + abs(vel))
           uneg = 0.5*(vel - abs(vel))
           eastmsk = tmask_quick(i,j,k)*(1-tmask_quick(i-1,j,k))
           westmsk = tmask_quick(i+1,j,k)*(1-tmask_quick(i+2,j,k))
           flux_x(i,j,k) =  &
            vel*(quick_x(i,j,1)*Tracer_field_tau(i,j,k)+quick_x(i,j,2)*Tracer_field_tau(i+1,j,k))&
                - upos*(curv_xp(i,j,1)*tracer_quick(i+1,j,k) &
                + curv_xp(i,j,2)*tracer_quick(i,j,k) + &
                curv_xp(i,j,3)*(tracer_quick(i-1,j,k)*(1-eastmsk)+tracer_quick(i,j,k)*eastmsk)) &
                - uneg*(curv_xn(i,j,1)*(tracer_quick(i+2,j,k)*(1-westmsk)+tracer_quick(i+1,j,k)*westmsk) &
                +       curv_xn(i,j,2)*tracer_quick(i+1,j,k) &
                +       curv_xn(i,j,3)*tracer_quick(i,j,k))
        enddo
     enddo


     ! recompute fluxes as upwind if tmask_limit flag satisfied 
     if(limit_with_upwind) then 
         do j=jsc,jec
            do i=isc-1,iec
               if(Tracer%tmask_limit(i,j,k)==1.0) then   
                   vel  = Adv_vel%uhrho_et(i,j,k)
                   upos = 0.5*(vel + abs(vel))
                   uneg = 0.5*(vel - abs(vel))
                   flux_x(i,j,k)  =  &
                     Grd%dyte(i,j)*(upos*Tracer_field_taum1(i,j,k) + uneg*Tracer_field_taum1(i+1,j,k)) &
                     *Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)  
               endif
            enddo
         enddo
         do j=jsc-1,jec
            do i=isc,iec
               if(Tracer%tmask_limit(i,j,k)==1.0) then   
                   vel  = Adv_vel%vhrho_nt(i,j,k)
                   upos = 0.5*(vel + abs(vel))
                   uneg = 0.5*(vel - abs(vel))
                   flux_y(i,j,k) =  &
                     Grd%dxtn(i,j)*(upos*Tracer_field_taum1(i,j,k) + uneg*Tracer_field_taum1(i,j+1,k)) &
                     *Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k) 
               endif
            enddo
         enddo
     endif

  enddo  ! end of k-loop 

  ! to handle redundancies across the bipolar fold 
  if (Grd%tripolar) call mpp_update_domains(flux_x, flux_y, Dom_flux%domain2d, gridtype=CGRID_NE) 

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           horz_advect_tracer_quicker(i,j,k) = &
              Grd%tmask(i,j,k)*(flux_x(i,j,k) - flux_x(i-1,j,k) + flux_y(i,j,k) - flux_y(i,j-1,k))*Grd%datr(i,j)
        enddo
     enddo
  enddo

  call mpp_clock_end(id_clock_quick_horz)

end function horz_advect_tracer_quicker
! </FUNCTION> NAME="horz_advect_tracer_quicker"


!#######################################################################
! <FUNCTION NAME="horz_advect_tracer_quickmom3">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers from quicker using MOM3 masking. 
! This method has proven useful for reproducing MOM3 results with later MOM.  
! </DESCRIPTION>
!
function horz_advect_tracer_quickmom3(Adv_vel, Tracer, Tracer_field_taum1, Tracer_field_tau)

  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field_taum1 
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field_tau

  real,dimension(isc:iec,jsc:jec,nk) :: horz_advect_tracer_quickmom3

  integer  :: i, j, k, jp1, jp2, jp2_mod
  real     :: upos, uneg, vel

  call mpp_clock_begin(id_clock_quickmom3_horz)

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
     '==>Error from ocean_tracer_advect (horz_advect_tracer_quickmom3): ocean_tracer_advect_mod not yet initialized')
  endif 

  tracer_quick = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tracer_quick(i,j,k) = Tracer_field_taum1(i,j,k)
        enddo
     enddo
  enddo
  call mpp_update_domains (tracer_quick, Dom_quicker%domain2d)

  flux_x = 0.0
  flux_y = 0.0

  do k=1,nk

     do j=jsc-1,jec
        do i=isc,iec
           jp1 = min(j+1, jec)
           jp2 = min(j+2, jec)
           jp2_mod = nint(tmask_quick(i,jp2,k)*jp2 + (1.0-tmask_quick(i,jp2,k))*jp1)

           vel  = Grd%dxtn(i,j)*Adv_vel%vhrho_nt(i,j,k)
           upos = 0.5*(vel + abs(vel))*tmask_quick(i,j-1,k)*tmask_quick(i,j,k)*tmask_quick(i,j+1,k)
           uneg = 0.5*(vel - abs(vel))*tmask_quick(i,jp2,k)*tmask_quick(i,j+1,k)*tmask_quick(i,j,k)

           flux_y(i,j,k) = vel*(quick_y(i,j,1)*Tracer_field_tau(i,j,k)+  &
                quick_y(i,j,2)*Tracer_field_tau(i,j+1,k))                &
                - upos*(curv_yp(i,j,1)*tracer_quick(i,j+1,k) +           &
                curv_yp(i,j,2)*tracer_quick(i,j,k)   +                   &
                curv_yp(i,j,3)*tracer_quick(i,j-1,k))                    &
                - uneg*(curv_yn(i,j,1)*tracer_quick(i,jp2,k) +           &
                curv_yn(i,j,2)*tracer_quick(i,j+1,k)   +                 &
                curv_yn(i,j,3)*tracer_quick(i,j,k))        
        enddo
     enddo

     do j=jsc,jec
        do i=isc-1,iec
           vel  = Grd%dyte(i,j)*Adv_vel%uhrho_et(i,j,k)
           upos = 0.5*(vel + abs(vel))*tmask_quick(i-1,j,k)*tmask_quick(i,j,k)*tmask_quick(i+1,j,k)
           uneg = 0.5*(vel - abs(vel))*tmask_quick(i+2,j,k)*tmask_quick(i+1,j,k)*tmask_quick(i,j,k)

           flux_x(i,j,k) = vel*(quick_x(i,j,1)*Tracer_field_tau(i,j,k)+ &
                quick_x(i,j,2)*Tracer_field_tau(i+1,j,k))               &
                - upos*(curv_xp(i,j,1)*tracer_quick(i+1,j,k)            &
                + curv_xp(i,j,2)*tracer_quick(i,j,k) +                  &
                curv_xp(i,j,3)*tracer_quick(i-1,j,k))                   &
                - uneg*(curv_xn(i,j,1)*tracer_quick(i+2,j,k)            &
                +       curv_xn(i,j,2)*tracer_quick(i+1,j,k)            &
                +       curv_xn(i,j,3)*tracer_quick(i,j,k))
        enddo
     enddo


     ! recompute fluxes as upwind if tmask_limit flag satisfied 
     if(limit_with_upwind) then 
         do j=jsc,jec
            do i=isc-1,iec
               if(Tracer%tmask_limit(i,j,k)==1.0) then   
                   vel  = Adv_vel%uhrho_et(i,j,k)
                   upos = 0.5*(vel + abs(vel))
                   uneg = 0.5*(vel - abs(vel))
                   flux_x(i,j,k) = Grd%dyte(i,j)*(upos*Tracer_field_taum1(i,j,k) + uneg*Tracer_field_taum1(i+1,j,k)) &
                        *Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)  
               endif
            enddo
         enddo
         do j=jsc-1,jec
            do i=isc,iec
               if(Tracer%tmask_limit(i,j,k)==1.0) then   
                   vel  = Adv_vel%vhrho_nt(i,j,k)
                   upos = 0.5*(vel + abs(vel))
                   uneg = 0.5*(vel - abs(vel))
                   flux_y(i,j,k) = Grd%dxtn(i,j)*(upos*Tracer_field_taum1(i,j,k) + uneg*Tracer_field_taum1(i,j+1,k)) &
                        *Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k) 
               endif
            enddo
         enddo
     endif

  enddo ! end of k-loop

  ! to handle redundancies across the bipolar fold 
  if (Grd%tripolar) call mpp_update_domains(flux_x, flux_y, Dom_flux%domain2d, gridtype=CGRID_NE) 

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           horz_advect_tracer_quickmom3(i,j,k) = &
           Grd%tmask(i,j,k)*(flux_x(i,j,k) - flux_x(i-1,j,k) + flux_y(i,j,k) - flux_y(i,j-1,k))*Grd%datr(i,j)
        enddo
     enddo
  enddo

  call mpp_clock_end(id_clock_quickmom3_horz)

end function horz_advect_tracer_quickmom3
! </FUNCTION> NAME="horz_advect_tracer_quickmom3"


!#######################################################################
! <FUNCTION NAME="vert_advect_tracer_upwind">
!
! <DESCRIPTION>
! Compute vertical advection of tracers from first order upwind. 
! This scheme is positive definite, but very diffusive. 
! </DESCRIPTION>
!
function vert_advect_tracer_upwind(Adv_vel, Tracer_field)

  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field

  real,dimension(isc:iec,jsc:jec,nk) :: vert_advect_tracer_upwind

  real,dimension(isc:iec,jsc:jec)    :: ft1, ft2 
  real                               :: velocity, wpos, wneg
  integer                            :: k, i, j, kp1

  call mpp_clock_begin(id_clock_up_vert)
 
  ft1  = 0.0
  do k=1,nk
    kp1 = min(k+1,nk)
    do j=jsc,jec
      do i=isc,iec
        velocity = 0.5*Adv_vel%wrho_bt(i,j,k)
        wpos     = velocity + abs(velocity) 
        wneg     = velocity - abs(velocity) 
        ft2(i,j) = (wneg*Tracer_field(i,j,k) + wpos*Tracer_field(i,j,kp1)) &
                   *Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1) 
        flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
        vert_advect_tracer_upwind(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
        ft1(i,j) = ft2(i,j)
      enddo
    enddo
  enddo

  call mpp_clock_end(id_clock_up_vert)

end function vert_advect_tracer_upwind
! </FUNCTION> NAME="vert_advect_tracer_upwind"


!#######################################################################
! <FUNCTION NAME="vert_advect_tracer_2nd_order">
!
! <DESCRIPTION>
! Compute vertical advection of tracers from second order centered
! differences. 
! </DESCRIPTION>
!
function vert_advect_tracer_2nd_order(Adv_vel, Tracer_field)

  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field

  real,dimension(isc:iec, jsc:jec, nk)     :: vert_advect_tracer_2nd_order 
  integer                                  :: k, kp
  real,dimension(isc:iec,jsc:jec)          :: ft1, ft2

  call mpp_clock_begin(id_clock_2nd_vert)

  ft1 = 0.0
  
  do k=1,nk
    kp = min(k+1,nk)
    ft2(isc:iec,jsc:jec) = Adv_vel%wrho_bt(isc:iec,jsc:jec,k)*0.5*(Tracer_field(isc:iec,jsc:jec,k)+&
                           Tracer_field(isc:iec,jsc:jec,kp))
    flux_z(isc:iec,jsc:jec,k) = Grd%dat(isc:iec,jsc:jec)*ft2(isc:iec,jsc:jec)
    vert_advect_tracer_2nd_order(isc:iec,jsc:jec,k) = Grd%tmask(isc:iec,jsc:jec,k)* &
                                                      (ft1(isc:iec,jsc:jec)-ft2(isc:iec,jsc:jec))
    ft1(:,:) = ft2(:,:)
  enddo

  call mpp_clock_end(id_clock_2nd_vert)

end function vert_advect_tracer_2nd_order
! </FUNCTION> NAME="vert_advect_tracer_2nd_order"


!#######################################################################
! <FUNCTION NAME="vert_advect_tracer_4th_order">
!
! <DESCRIPTION>
! Compute vertical advection of tracers from fourth order centered
! differences. NOTE: this code does not account for non-uniform grids.   
! </DESCRIPTION>
!
function vert_advect_tracer_4th_order(Adv_vel, Tracer_field)

  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field

  real,dimension(isc:iec, jsc:jec, nk)     :: vert_advect_tracer_4th_order 
  integer                                  :: k, i, j, kp2, km1, p1, p2
  real,dimension(isc:iec,jsc:jec)          :: ft1, ft2

  call mpp_clock_begin(id_clock_4th_vert)

  ft1(isc:iec,jsc:jec) = 0.0
  do k=1,nk
    km1 = max(k-1,1)
    p1  = min(k+1,nk)
    p2  = min(k+2,nk)
    do j=jsc,jec
      do i=isc,iec
        kp2   = nint(Grd%tmask(i,j,p2)*p2 + (1.0-Grd%tmask(i,j,p2))*p1)    
        ft2(i,j) = Adv_vel%wrho_bt(i,j,k)*(a4*(Tracer_field(i,j,k)+&
                   Tracer_field(i,j,p1)) + b4*(Tracer_field(i,j,km1)+&
                   Tracer_field(i,j,kp2)))
        flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
        vert_advect_tracer_4th_order(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
      enddo
    enddo
    ft1(isc:iec,jsc:jec) = ft2(isc:iec,jsc:jec)
  enddo

  call mpp_clock_end(id_clock_4th_vert)

end function vert_advect_tracer_4th_order
! </FUNCTION> NAME="vert_advect_tracer_4th_order"


!#######################################################################
! <FUNCTION NAME="vert_advect_tracer_6th_order">
!
! <DESCRIPTION>
! Compute vertical advection of tracers from sixth order centered
! differences. NOTE: this code does not account for non-uniform grids.   
! </DESCRIPTION>
!
function vert_advect_tracer_6th_order(Adv_vel, Tracer_field)

  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field

  real,dimension(isc:iec,jsc:jec)    :: ft1, ft2
  real,dimension(isc:iec,jsc:jec,nk) :: vert_advect_tracer_6th_order 
  integer                            :: i, j, k, km1, kp2, p1, p2

  call mpp_clock_begin(id_clock_6th_vert)

  ft1(isc:iec,jsc:jec) = 0.0
  do k=1,nk
    km1 = max(k-1,1)
    p1 = min(k+1,nk)
    p2 = min(k+2,nk)
    if (k <= 2 .or. k >= nk-1) then
      do j=jsc,jec
        do i=isc,iec
          kp2   = nint(Grd%tmask(i,j,p2)*p2 + (1.0-Grd%tmask(i,j,p2))*p1)    
          ft2(i,j) = Adv_vel%wrho_bt(i,j,k)*(a4*(Tracer_field(i,j,k)+&
                     Tracer_field(i,j,p1)) + &
                     b4*(Tracer_field(i,j,km1)+Tracer_field(i,j,kp2)))
          flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)          
          vert_advect_tracer_6th_order(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
        enddo
      enddo
    else
      do j=jsc,jec
        do i=isc,iec
          kp2   = nint(Grd%tmask(i,j,p2)*p2 + (1.0-Grd%tmask(i,j,p2))*p1)    
          if (k < Grd%kmt(i,j)-2) then
            ft2(i,j) = Adv_vel%wrho_bt(i,j,k)*(a6*(Tracer_field(i,j,k)+&
                       Tracer_field(i,j,p1)) + &
                       b6*(Tracer_field(i,j,km1)+ &
                       Tracer_field(i,j,kp2))&
                       + c6*(Tracer_field(i,j,k-2)+ &
                       Tracer_field(i,j,k+3)))
          else
            ft2(i,j) = Adv_vel%wrho_bt(i,j,k)*(a4*(Tracer_field(i,j,k)+&
                       Tracer_field(i,j,p1)) +&
                       b4*(Tracer_field(i,j,km1)+&
                       Tracer_field(i,j,kp2)))
        endif
          flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)        
          vert_advect_tracer_6th_order(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
        enddo
      enddo
    endif
    ft1(isc:iec,jsc:jec) = ft2(isc:iec,jsc:jec)
  enddo

  call mpp_clock_end(id_clock_6th_vert)

end function vert_advect_tracer_6th_order
! </FUNCTION> NAME="vert_advect_tracer_6th_order"


!#######################################################################
! <FUNCTION NAME="vert_advect_tracer_quicker">
!
! <DESCRIPTION>
! Compute vertical advection of tracers from quicker.
! </DESCRIPTION>
!
function vert_advect_tracer_quicker(Adv_vel, Tracer, Tracer_field_taum1, Tracer_field_tau)

  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field_taum1 
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field_tau

  real,dimension(isc:iec,jsc:jec,nk)       :: vert_advect_tracer_quicker 
  integer                                  :: k, i, j
  real,dimension(isc:iec,jsc:jec)          :: ft1, ft2
  real                                     :: upos, uneg, vel, upmsk
  integer                                  :: km1, kp2, kp1, p2
  
  call mpp_clock_begin(id_clock_quick_vert)

  ft1(isc:iec,jsc:jec) = 0.0
  do k=1,nk
    km1 = max(k-1,1)
    kp1 = min(k+1,nk)
    p2  = min(k+2,nk)
    do j=jsc,jec
      do i=isc,iec
        vel  = Adv_vel%wrho_bt(i,j,k)
        upos = 0.5*(vel + abs(vel))
        uneg = 0.5*(vel - abs(vel))
        if(Tracer%tmask_limit(i,j,k)==1.0) then 
          ft2(i,j) = (uneg*Tracer_field_taum1(i,j,k) + upos*Tracer_field_taum1(i,j,kp1)) &
                          *Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1)
        else 
          kp2   = nint(Grd%tmask(i,j,p2)*p2 + (1.0-Grd%tmask(i,j,p2))*kp1)
          upmsk = Grd%tmask(i,j,k)*(1-Grd%tmask(i,j,km1))
          ft2(i,j) = vel*(quick_z(k,1)*Tracer_field_tau(i,j,k) +      &
                        quick_z(k,2)*Tracer_field_tau(i,j,kp1))       &
                    -uneg*(curv_zp(k,1)*Tracer_field_taum1(i,j,kp1)   &
                          +curv_zp(k,2)*Tracer_field_taum1(i,j,k)     &
                          +curv_zp(k,3)*Tracer_field_taum1(i,j,km1))  &
                    -upos*(curv_zn(k,1)*Tracer_field_taum1(i,j,kp2)   &
                          +curv_zn(k,2)*Tracer_field_taum1(i,j,kp1)   &
                          +curv_zn(k,3)*(Tracer_field_taum1(i,j,k)*(1-upmsk)+Tracer_field_taum1(i,j,kp1)*upmsk))

        endif 
        flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
        vert_advect_tracer_quicker(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
      enddo
    enddo
    ft1(isc:iec,jsc:jec) = ft2(isc:iec,jsc:jec)
  enddo

  call mpp_clock_end(id_clock_quick_vert)
  
end function vert_advect_tracer_quicker
! </FUNCTION> NAME="vert_advect_tracer_quicker"


!#######################################################################
! <FUNCTION NAME="vert_advect_tracer_quickmom3">
!
! <DESCRIPTION>
! Compute vertical advection of tracers from quicker using MOM3 masking.
! This method has proven useful for reproducing MOM3 results with later MOM.  
! </DESCRIPTION>
!
function vert_advect_tracer_quickmom3(Adv_vel, Tracer, Tracer_field_taum1, Tracer_field_tau)

  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field_taum1 
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field_tau

  real,dimension(isc:iec,jsc:jec,nk)       :: vert_advect_tracer_quickmom3 
  integer                                  :: k, i, j
  real,dimension(isc:iec,jsc:jec)          :: ft1, ft2
  real                                     :: upos, uneg, vel
  integer                                  :: kp1
  
  call mpp_clock_begin(id_clock_quickmom3_vert)

  ft1(isc:iec,jsc:jec) = 0.0
  do k=1,nk

     if (k == 1) then
      do j=jsc,jec
        do i=isc,iec
          vel  = Adv_vel%wrho_bt(i,j,k)
          upos = 0.5*(vel + abs(vel))*Grd%tmask(i,j,k+2)*Grd%tmask(i,j,k+1)*Grd%tmask(i,j,k)
          if(Tracer%tmask_limit(i,j,k)==1.0) then 
            kp1 = min(k+1,nk)
          ft2(i,j) = (uneg*Tracer_field_taum1(i,j,k) +               &
                      upos*Tracer_field_taum1(i,j,kp1))              &
                     *Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1)
          else 
            ft2(i,j) = vel*(quick_z(k,1)*Tracer_field_tau(i,j,k) +   &
                        quick_z(k,2)*Tracer_field_tau(i,j,k+1))      &
                    -upos*(curv_zn(k,1)*Tracer_field_taum1(i,j,k+2)  &
                          +curv_zn(k,2)*Tracer_field_taum1(i,j,k+1)  &
                          +curv_zn(k,3)*Tracer_field_taum1(i,j,k))
          endif 
         flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
         vert_advect_tracer_quickmom3(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
       enddo
     enddo

    elseif (k == nk-1) then
      do j=jsc,jec
        do i=isc,iec
         vel  = Adv_vel%wrho_bt(i,j,k)
         uneg = 0.5*(vel - abs(vel))*Grd%tmask(i,j,k-1)*Grd%tmask(i,j,k)*Grd%tmask(i,j,k+1)
          if(Tracer%tmask_limit(i,j,k)==1.0) then 
            kp1 = min(k+1,nk)
            ft2(i,j) = (uneg*Tracer_field_taum1(i,j,k) +     &
                        upos*Tracer_field_taum1(i,j,kp1))    &
                        *Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1)
          else
            ft2(i,j) = vel*(quick_z(k,1)*Tracer_field_tau(i,j,k) +   &
                        quick_z(k,2)*Tracer_field_tau(i,j,k+1))      &
                    -uneg*(curv_zp(k,1)*Tracer_field_taum1(i,j,k+1)  &
                          +curv_zp(k,2)*Tracer_field_taum1(i,j,k)    &
                          +curv_zp(k,3)*Tracer_field_taum1(i,j,k-1))
          endif 
         flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
         vert_advect_tracer_quickmom3(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
        enddo
      enddo

    else
      do j=jsc,jec
        do i=isc,iec
         vel  = Adv_vel%wrho_bt(i,j,k)
         upos = 0.5*(vel + abs(vel))*Grd%tmask(i,j,k+2)*Grd%tmask(i,j,k+1)*Grd%tmask(i,j,k)
         uneg = 0.5*(vel - abs(vel))*Grd%tmask(i,j,k-1)*Grd%tmask(i,j,k)*Grd%tmask(i,j,k+1)
          if(Tracer%tmask_limit(i,j,k)==1.0) then 
            kp1 = min(k+1,nk)
            ft2(i,j) = (uneg*Tracer_field_taum1(i,j,k) +     &
                        upos*Tracer_field_taum1(i,j,kp1))    &
                        *Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1)
          else


            ft2(i,j) = vel*(quick_z(k,1)*Tracer_field_tau(i,j,k) +   &
                        quick_z(k,2)*Tracer_field_tau(i,j,k+1))      &
                    -upos*(curv_zn(k,1)*Tracer_field_taum1(i,j,k+2)  &
                          +curv_zn(k,2)*Tracer_field_taum1(i,j,k+1)  &
                          +curv_zn(k,3)*Tracer_field_taum1(i,j,k))   &
                    -uneg*(curv_zp(k,1)*Tracer_field_taum1(i,j,k+1)  &
                          +curv_zp(k,2)*Tracer_field_taum1(i,j,k)    &
                          +curv_zp(k,3)*Tracer_field_taum1(i,j,k-1))
          endif 
         flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
         vert_advect_tracer_quickmom3(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
        enddo
      enddo
     endif


    ft1(isc:iec,jsc:jec) = ft2(isc:iec,jsc:jec)
  enddo

  call mpp_clock_end(id_clock_quickmom3_vert)
  
end function vert_advect_tracer_quickmom3
! </FUNCTION> NAME="vert_advect_tracer_quickmom3"


!#######################################################################
! <FUNCTION NAME="advect_tracer_mdfl_sup_b">
!
! <DESCRIPTION>
!-------
! NOTE: This is legacy code. It may exhibit problems with limiters 
! that are improperly implemented in the case of generalized vertical level 
! coordinates of MOM. Under certain circumstances, the resulting tracer 
! concentration can exhibit non-monotonic behaviour (i.e., produce extrema).  
!
! Sept 2011 
! Stephen.Griffies
! Alistair.Adcroft
!-------
!
! Compute tendency due to 3D advection of tracers using a 
! multi-dimensional flux-limited method.  This method differs from 
! other MOM advection methods in the following ways:
!
! 1) Horizontal and vertical advection are combined
! 2) Calculations of the three coordinates (Z, X, Y) are performed sequentially as updates
!        to the tracer, so that the advection components for X and Y depend on Z and Z,Y
!    respectively. This helps limit the flux.
! 3) During the update for each direction, the 2nd order Super-B flux limiter is applied.
! 4) Flux divergence is included within the calculation to also help limit the flux:
!          - During the update for each direction, the divergence in each direction is added.
!          - During the overall tendency calculated, the divergence in all three directions 
!                 is removed.
! 5) All fluxes are functions of the tracer field at the taum1 time step.  This 
!    means that this method is ideally suited for the twolevel time stepping scheme,
!    in which "taum1=tau", thus enabling twice the tracer time step available for the 
!    threelevel scheme. 
!
!The calculation proceeds as follows:
!
! IMPORTANT NOTE: If this scheme is used at all, it must be used as the option for BOTH
! horizontal and vertical advection.  In the the tracer tendency, it is applied as the 
! horizontal term, but applies to vertical as well, for which case the vertical term in
! the tracer tendency equation is set to zero.
!
! This scheme was ported to mom4 from the MIT-GCM by John Dunne and Alistair Adcroft  
! during Summer 2003 
!
! </DESCRIPTION>
!
function advect_tracer_mdfl_sup_b(Time, Adv_vel, Thickness, Tracer_field, dtime)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field
  real,                         intent(in) :: dtime

  real,dimension(isc  :iec,jsc-1:jec)      :: fn
  real,dimension(isc-1:iec,jsc  :jec)      :: fe
  real,dimension(isc:iec,jsc:jec,nk)       :: advect_tracer_mdfl_sup_b 
  integer                                  :: k, i, j, kp1, kp2, km1, tau, taum1
  real,dimension(isc:iec,jsc:jec)          :: ftp, fbt, wkm1
  real,dimension(isc-2:iec+2,jsc-2:jec+2)  :: tmp_updated_tracer
  real                                     :: Rjm, Rj, Rjp, Cr, cfl, massflux

  call mpp_clock_begin(id_clock_mdfl_sup_b)

  ! initialize to zero local arrays 
  tmp_updated_tracer = 0.0
  ftp                = 0.0
  fbt                = 0.0
  wkm1               = 0.0
  fn                 = 0.0
  fe                 = 0.0

  ! set time indices 
  tau   = Time%tau
  taum1 = Time%taum1

!
!Boundary condition at surface
!
    do j=jsc,jec
      do i=isc,iec
        ftp(i,j) = 0.0
        wkm1(i,j) = 0.0
      enddo !i
    enddo !j
!
! Main loop over all depths
!
  do k=1,nk

!
! Get tracer field on computational grid
!
    do j=jsc,jec
      do i=isc,iec
        tmp_updated_tracer(i,j) = Tracer_field(i,j,k)
      enddo !i
    enddo !j
!
! Account for boundaries
!
    kp1 = min(k+1,nk)
    kp2 = min(k+2,nk)
    km1 = max(k-1,1)
!
! Calculate flux at bottom of box and update the tracer
!
    do j=jsc,jec
      do i=isc,iec
        Rjp = ( Tracer_field(i,j,km1) - Tracer_field(i,j,k) )  &
                  * tmask_mdfl(i,j,km1) *       tmask_mdfl(i,j,k)
        Rj  = ( Tracer_field(i,j,k) - Tracer_field(i,j,kp1) )    &
                  * tmask_mdfl(i,j,k) *       tmask_mdfl(i,j,kp1)
        Rjm  = ( Tracer_field(i,j,kp1) - Tracer_field(i,j,kp2) )   &
                  * tmask_mdfl(i,j,kp1) *        tmask_mdfl(i,j,kp2)
        massflux = Grd%dat(i,j) * Adv_vel%wrho_bt(i,j,k)
        cfl = abs(Adv_vel%wrho_bt(i,j,k) * dtime  / Thickness%rho_dzt(i,j,k,tau))
        if ( Rj .NE. 0.0) then
         if ( massflux .GT. 0.0) then
           Cr = Rjm / Rj
         else
           Cr = Rjp / Rj
         endif
        else
         if ( massflux .GT. 0.0) then
           Cr = Rjm * 1.0e20
         else
           Cr = Rjp * 1.0e20
         endif
        endif
        Cr = max(0.0, max(min(1.0, 2.0 * Cr), min(2.0, Cr)))
        fbt(i,j) = 0.5 * ( massflux                                     &
                 * ( Tracer_field(i,j,k) + Tracer_field(i,j,kp1) )      &
                 - abs(massflux) * ( 1.0 + (cfl - 1 ) * Cr ) * Rj       &
                  ) * tmask_mdfl(i,j,kp1) * tmask_mdfl(i,j,k)
        tmp_updated_tracer(i,j) = tmp_updated_tracer(i,j)               &
               + dtime / Thickness%rho_dzt(i,j,k,tau) * (               &
                 Grd%datr(i,j) * ( fbt(i,j) - ftp(i,j))                 &
                 + Tracer_field(i,j,k) *                                &
                 (wkm1(i,j) - Adv_vel%wrho_bt(i,j,k)) )
        flux_z(i,j,k) = fbt(i,j)
        ftp(i,j) = fbt(i,j)
      enddo !i
    enddo !j
!
! Update the two-point tracer halo
!
    call mpp_update_domains (tmp_updated_tracer, Dom_mdfl%domain2d)
!
! Calculate flux at the eastern wall of the boxes
!
    do j=jsc,jec
      do i=isc-1,iec
        Rjp = ( tmp_updated_tracer(i+2,j) - tmp_updated_tracer(i+1,j) )       &
                      * tmask_mdfl(i+2,j,k) *       tmask_mdfl(i+1,j,k)
        Rj =  ( tmp_updated_tracer(i+1,j) - tmp_updated_tracer( i ,j) )       &
                      * tmask_mdfl(i+1,j,k) *       tmask_mdfl( i ,j,k)
        Rjm = ( tmp_updated_tracer( i ,j) - tmp_updated_tracer(i-1,j) )       &
                      * tmask_mdfl( i ,j,k) *       tmask_mdfl(i-1,j,k)
        massflux = Grd%dyte(i,j) * Adv_vel%uhrho_et(i,j,k)
        cfl = abs( Adv_vel%uhrho_et(i,j,k) * dtime * 2.0                      &
        / ((Thickness%rho_dzt(i,j,k,tau) + Thickness%rho_dzt(i+1,j,k,tau)) * Grd%dxte(i,j)))
        if ( Rj .NE. 0.0) then
         if ( massflux .GT. 0.0) then
           Cr = Rjm / Rj
         else
           Cr = Rjp / Rj
         endif
        else
         if ( massflux .GT. 0.0) then
           Cr = Rjm * 1.0e20
         else
           Cr = Rjp * 1.0e20
         endif
        endif
        Cr = max(0.0, max(min(1.0, 2.0 * Cr), min(2.0, Cr)))
        fe(i,j) = 0.5 * ( massflux                                           &
                   * (tmp_updated_tracer(i+1,j) + tmp_updated_tracer(i,j) )   &
                 - abs(massflux) * ( 1.0 + (cfl - 1 ) * Cr ) * Rj            &
                  ) * tmask_mdfl(i,j,k) * tmask_mdfl(i+1,j,k)
      enddo !i
    enddo !j
!
! Update the tracer
!
    do j=jsc,jec
      do i=isc,iec
        tmp_updated_tracer(i,j) = tmp_updated_tracer(i,j)                           &
            + dtime * tmask_mdfl(i,j,k) * Grd%datr(i,j) / Thickness%rho_dzt(i,j,k,tau)  &
            * (fe(i-1,j) - fe(i,j)                                                  &
            + Tracer_field(i,j,k)* (                                          &
                      Grd%dyte(i,j) *   Adv_vel%uhrho_et( i ,j,k)                       &
                    - Grd%dyte(i-1,j) * Adv_vel%uhrho_et(i-1,j,k) ) )
        flux_x(i,j,k) = fe(i,j)
      enddo !i
    enddo !j
!
! Update the two-point tracer halo
!
    call mpp_update_domains (tmp_updated_tracer, Dom_mdfl%domain2d)
!
! Calculate flux at the northern wall of the boxes
!
    do j=jsc-1,jec
      do i=isc,iec
        Rjp = ( tmp_updated_tracer(i,j+2) - tmp_updated_tracer(i,j+1) )       &
                      * tmask_mdfl(i,j+2,k) *       tmask_mdfl(i,j+1,k)
        Rj =  ( tmp_updated_tracer(i,j+1) - tmp_updated_tracer(i, j ) )       &
                      * tmask_mdfl(i,j+1,k) *       tmask_mdfl(i, j ,k)
        Rjm = ( tmp_updated_tracer(i, j ) - tmp_updated_tracer(i,j-1) )       &
                      * tmask_mdfl(i, j ,k) *       tmask_mdfl(i,j-1,k)
        massflux = Grd%dxtn(i,j) * Adv_vel%vhrho_nt(i,j,k)
        cfl = abs(Adv_vel%vhrho_nt(i,j,k) * dtime * 2.0                       &
        / ((Thickness%rho_dzt(i,j,k,tau) + Thickness%rho_dzt(i,j+1,k,tau)) * Grd%dytn(i,j)))
        if ( Rj .NE. 0.0) then
         if ( massflux .GT. 0.0) then
           Cr = Rjm / Rj
         else
           Cr = Rjp / Rj
         endif
        else
         if ( massflux .GT. 0.0) then
           Cr = Rjm * 1.0e20
         else
           Cr = Rjp * 1.0e20
         endif
        endif
        Cr = max(0.0, max(min(1.0, 2.0 * Cr), min(2.0, Cr)))
        fn(i,j) = 0.5 * ( massflux                                           &
                   * (tmp_updated_tracer(i,j+1) + tmp_updated_tracer(i,j) )   &
                 - abs(massflux) * ( 1.0 + (cfl - 1.0 ) * Cr ) * Rj          &
                  ) * tmask_mdfl(i,j,k) * tmask_mdfl(i,j+1,k)
      enddo !i
    enddo !j
!
! Update the tracer
!
    do j=jsc,jec
      do i=isc,iec
        tmp_updated_tracer(i,j) = tmp_updated_tracer(i,j)                           &
         + dtime * Grd%tmask(i,j,k) * Grd%datr(i,j) / Thickness%rho_dzt(i,j,k,tau)  &
         * (fn(i,j-1) - fn(i,j))
        flux_y(i,j,k) = fn(i,j)
      enddo !i
    enddo !j
!
! Calculate the overall tendency
!
    do j=jsc,jec
      do i=isc,iec
        tmp_updated_tracer(i,j) = tmp_updated_tracer(i,j)                    &
           + dtime * Tracer_field(i,j,k) / Thickness%rho_dzt(i,j,k,tau) * (  &
                      Adv_vel%wrho_bt(i,j,k) - wkm1(i,j)                     &
             + Grd%datr(i,j)*(                                               &
                         ( Grd%dyte(i-1,j) * Adv_vel%uhrho_et(i-1,j,k)       &
                         - Grd%dyte( i ,j) * Adv_vel%uhrho_et( i ,j,k))))
!
! By convention, advection is applied as a negative tendency (i.e., convergence)
!
        advect_tracer_mdfl_sup_b(i,j,k) = - Thickness%rho_dzt(i,j,k,tau) *   &
          ( tmp_updated_tracer(i,j) - Tracer_field(i,j,k) )                  &
            / dtime * tmask_mdfl(i,j,k)
      enddo !i
    enddo !j
!
! Update vertical velocity for next level
!
    do j=jsc,jec
      do i=isc,iec
        wkm1(i,j) = Adv_vel%wrho_bt(i,j,k)
      enddo !i
    enddo !j

  enddo !k
  
  call mpp_clock_end(id_clock_mdfl_sup_b)

end function advect_tracer_mdfl_sup_b
! </FUNCTION> NAME="advect_tracer_mdfl_sup_b"



!#######################################################################
! <FUNCTION NAME="advect_tracer_mdfl_sweby_test">
!
! <DESCRIPTION>
! Compute tendency due to 3D advection of tracers using a multi-dimensional flux-limited
! method.  This method differs from other methods in the following ways:
!
! 1) Horizontal and vertical advection are combined
! 2) Calculations of the three coordinates (Z, X, Y) are performed sequentially as updates
!        to the tracer, so that the advection components for X and Y depend on Z and Z,Y
!        respectively... This helps limit the flux.
! 3) During the update for each direction, the 3rd order Sweby flux limiter is applied.
! 4) Flux divergence is included within the calculation to also help limit the flux:
!          - During the update for each direction, the divergence in each direction is added.
!          - During the overall tendency calculated, the divergence in all three directions
!                is removed.
! 5) All fluxes are functions of the tracer field at the taum1 time step.  This 
!    means that this method is ideally suited for the twolevel time stepping scheme,
!    in which "taum1=tau", thus enabling twice the tracer time step available for the 
!    threelevel scheme. 
!
! The calculation proceeds as follows:
!
! IMPORTANT NOTE: If this scheme is used at all, it must be used as the option for BOTH
! horizontal and vertical advection.  In the tracer tendency, it is applied as the 
! horizontal term, but applies to vertical as well, for which case the vertical term in
! the tracer tendency equation is set to zero.
!
! This scheme was ported to mom4 from the MIT-GCM by John Dunne and Alistair Adcroft  
! during Summer 2003 
!
! Griffies: 5/27/04
! Optimized by filling 3d arrays prior to sending for mpp_update_domains
!
! 07/11/2007 by Stephen.Griffies
! When the nonlinear limiters are removed (sweby_limiter=0.0), 
! this scheme reduces to a linear direct space-time method 
! (ADVECT_DST_LINEAR).  
! 
! This is a test version of the mdfl_sweby scheme. It is not supported for 
! general use. 
!
! </DESCRIPTION>
!
function advect_tracer_mdfl_sweby_test(Time, Adv_vel, Thickness, Tracer_field, dtime, sweby_limiter)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field
  real,                         intent(in) :: dtime
  real,                         intent(in) :: sweby_limiter

  real,dimension(isc:iec,jsc:jec)   :: ftp
  real,dimension(isc:iec,jsc:jec)   :: fbt
  real,dimension(isc:iec,jsc:jec)   :: wkm1

  real,dimension(isc:iec,jsc:jec,nk) :: advect_tracer_mdfl_sweby_test 

  integer :: i, j, k
  integer :: kp1, kp2, km1
  integer :: tau, taum1
  real    :: Rjm, Rj, Rjp, cfl, massflux
  real    :: d0, d1, thetaP, psiP 
  real    :: thetaM, psiM

!aja real :: tmin0,tmax0

  if(sweby_limiter==1.0) then
     call mpp_clock_begin(id_clock_mdfl_sweby_test)
  else 
     call mpp_clock_begin(id_clock_dst_linear_test)
  endif 

  ftp             = 0.0
  fbt             = 0.0
  wkm1            = 0.0
  tracer_mdfl     = 0.0
  mass_mdfl       = 0.0
  tracermass_mdfl = 0.0
  flux_x          = 0.0
  flux_y          = 0.0

  tau   = Time%tau
  taum1 = Time%taum1

!aja call get_tracer_stats(Tracer_field(isc:iec,jsc:jec,:),tmin0,tmax0)

  do k=1,nk

     do j=jsc,jec
        do i=isc,iec
           tracer_mdfl(i,j,k) = Tracer_field(i,j,k)
           mass_mdfl(i,j,k) = Thickness%rho_dzt(i,j,k,tau) * Grd%dat(i,j)
           tracermass_mdfl(i,j,k) = mass_mdfl(i,j,k) * Tracer_field(i,j,k)
        enddo
     enddo

     kp1 = min(k+1,nk)
     kp2 = min(k+2,nk)
     km1 = max(k-1,1)

     ! calculate flux at bottom of box and update the tracer
     do j=jsc,jec
        do i=isc,iec

           Rjp = (Tracer_field(i,j,km1) - Tracer_field(i,j,k))    &
                 *tmask_mdfl(i,j,km1)*tmask_mdfl(i,j,k)
           Rj  = (Tracer_field(i,j,k) - Tracer_field(i,j,kp1))    &
                 *tmask_mdfl(i,j,k)*tmask_mdfl(i,j,kp1)
           Rjm = (Tracer_field(i,j,kp1) - Tracer_field(i,j,kp2))  &
                 *tmask_mdfl(i,j,kp1)*tmask_mdfl(i,j,kp2)

           massflux = Grd%dat(i,j) * Adv_vel%wrho_bt(i,j,k)

           if (massflux*Thickness%rho_dzt(i,j,kp1,tau)>0.0) then
              cfl = abs(Adv_vel%wrho_bt(i,j,k)) * dtime / Thickness%rho_dzt(i,j,kp1,tau)
           elseif (massflux*Thickness%rho_dzt(i,j,k,tau)<0.0) then
              cfl = abs(Adv_vel%wrho_bt(i,j,k)) * dtime / Thickness%rho_dzt(i,j,k,tau)
           else
              cfl = 0.0
           endif

           d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
           d1 = (1.0 - cfl * cfl) * onesixth

           thetaP = Rjm / (sign(1.0e-30,Rj) + Rj)
           thetaM = Rjp / (sign(1.0e-30,Rj) + Rj)

           psiP = d0 + d1 * thetaP
           psiP = psiP*(1.0-sweby_limiter) +                 &
                  max(0.0, min(min(1.0, d0 + d1 * thetaP),   &
                  (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) ) &
                  *sweby_limiter

           psiM = d0 + d1 * thetaM
           psiM = psiM*(1.0-sweby_limiter) +                 &
                  max(0.0, min(min(1.0, d0 + d1 * thetaM),   &
                  (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) ) &
                  *sweby_limiter  

           fbt(i,j) =  0.5 * ( ( massflux + abs(massflux) )      &
                * ( Tracer_field(i,j,kp1) + psiP * Rj )          &  
                + ( massflux - abs(massflux) )                   &
                * ( Tracer_field(i,j, k ) - psiM * Rj ) )        &
                * tmask_mdfl(i,j,kp1) * tmask_mdfl(i,j,k)

           mass_mdfl(i,j,k) = mass_mdfl(i,j,k) + dtime * Grd%dat(i,j) * ( &
                    Adv_vel%wrho_bt(i,j,k) - wkm1(i,j)  )
 
           tracermass_mdfl(i,j,k) = tracermass_mdfl(i,j,k) + dtime * ( fbt(i,j) - ftp(i,j) )

           if (mass_mdfl(i,j,k)>0.) tracer_mdfl(i,j,k) = tracermass_mdfl(i,j,k) / mass_mdfl(i,j,k)

           flux_z(i,j,k) = fbt(i,j)
           ftp(i,j)      = fbt(i,j)
           wkm1(i,j) = Adv_vel%wrho_bt(i,j,k)

        enddo
     enddo

  enddo ! end of k-loop
  call mpp_update_domains (tracer_mdfl, Dom_mdfl%domain2d, flags=XUPDATE)
  call mpp_update_domains (tracermass_mdfl, Dom_mdfl%domain2d, flags=XUPDATE)
  call mpp_update_domains (mass_mdfl, Dom_mdfl%domain2d, flags=XUPDATE)

  ! calculate flux at the eastern wall of the boxes
  do k=1,nk

     do j=jsc,jec
        do i=isc-1,iec

           Rjp = (tracer_mdfl(i+2,j,k) - tracer_mdfl(i+1,j,k))    &
                 *tmask_mdfl(i+2,j,k)*tmask_mdfl(i+1,j,k)
           Rj  = (tracer_mdfl(i+1,j,k) - tracer_mdfl( i ,j,k) )   &
                 *tmask_mdfl(i+1,j,k)*tmask_mdfl(i,j,k)
           Rjm = (tracer_mdfl(i,j,k) - tracer_mdfl(i-1,j,k))      &
                 *tmask_mdfl(i,j,k)*tmask_mdfl(i-1,j,k)

           massflux = Grd%dyte(i,j) * Adv_vel%uhrho_et(i,j,k)

           if (massflux*mass_mdfl(i,j,k)>0.0) then                              
              cfl = abs(massflux) * dtime / mass_mdfl(i,j,k) 
           elseif (massflux*mass_mdfl(i+1,j,k)<0.0) then
              cfl = abs(massflux) * dtime / mass_mdfl(i+1,j,k) 
           else
              cfl = 0.0
           endif

           d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
           d1 = (1.0 - cfl * cfl) * onesixth

           thetaP = Rjm / (sign(1.0e-30,Rj) + Rj)
           thetaM = Rjp / (sign(1.0e-30,Rj) + Rj)

           psiP = d0 + d1 * thetaP
           psiP = psiP*(1.0-sweby_limiter) +                 & 
                  max(0.0, min(min(1.0, d0 + d1 * thetaP),   &
                  (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) ) &
                  *sweby_limiter

           psiM = d0 + d1 * thetaM
           psiM = psiM*(1.0-sweby_limiter) +                 & 
                  max(0.0, min(min(1.0, d0 + d1 * thetaM),   &
                  (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) ) &
                  *sweby_limiter

           flux_x(i,j,k) =  0.5 * ( ( massflux + abs(massflux) ) &
                * ( tracer_mdfl( i ,j,k) + psiP * Rj )           &
                + ( massflux - abs(massflux) )                   &
                * ( tracer_mdfl(i+1,j,k) - psiM * Rj ) )         &
                * tmask_mdfl(i,j,k) * tmask_mdfl(i+1,j,k)

        enddo
     enddo

     ! update the tracer
     do j=jsc,jec
        do i=isc,iec

           mass_mdfl(i,j,k) = mass_mdfl(i,j,k) + dtime * (         &
                  Grd%dyte(i-1,j) * Adv_vel%uhrho_et(i-1,j,k)      &
                - Grd%dyte(i,j) * Adv_vel%uhrho_et(i,j,k) )

           tracermass_mdfl(i,j,k) = tracermass_mdfl(i,j,k) + dtime * ( flux_x(i-1,j,k) - flux_x(i,j,k) )

           if (mass_mdfl(i,j,k)>0.) tracer_mdfl(i,j,k) = tracermass_mdfl(i,j,k) / mass_mdfl(i,j,k)

        enddo
     enddo

  enddo ! end of k-loop
  call mpp_update_domains (tracer_mdfl, Dom_mdfl%domain2d, flags=YUPDATE)
  call mpp_update_domains (tracermass_mdfl, Dom_mdfl%domain2d, flags=YUPDATE)
  call mpp_update_domains (mass_mdfl, Dom_mdfl%domain2d, flags=YUPDATE)

  do k=1,nk

     do j=jsc-1,jec
        do i=isc,iec

           Rjp = (tracer_mdfl(i,j+2,k) - tracer_mdfl(i,j+1,k))  &
                 *tmask_mdfl(i,j+2,k)*tmask_mdfl(i,j+1,k)
           Rj  = (tracer_mdfl(i,j+1,k) - tracer_mdfl(i,j,k))    &
                 *tmask_mdfl(i,j+1,k)*tmask_mdfl(i,j,k)
           Rjm = (tracer_mdfl(i,j,k) - tracer_mdfl(i,j-1,k))    &
                *tmask_mdfl(i,j,k)*tmask_mdfl(i,j-1,k)

           massflux = Grd%dxtn(i,j) * Adv_vel%vhrho_nt(i,j,k)

           if (massflux*mass_mdfl(i,j,k)>0.0) then                              
              cfl = abs(massflux) * dtime / mass_mdfl(i,j,k) 
           elseif (massflux*mass_mdfl(i,j+1,k)<0.0) then
              cfl = abs(massflux) * dtime / mass_mdfl(i,j+1,k) 
           else
              cfl = 0.0
           endif

           d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
           d1 = (1.0 - cfl * cfl) * onesixth

           thetaP = Rjm / (sign(1.0e-30,Rj) + Rj)
           thetaM = Rjp / (sign(1.0e-30,Rj) + Rj)

           psiP = d0 + d1 * thetaP
           psiP = psiP*(1.0-sweby_limiter) +                 & 
                  max(0.0, min(min(1.0, d0 + d1 * thetaP),   &
                  (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) ) &
                  *sweby_limiter

           psiM = d0 + d1 * thetaM
           psiM = psiM*(1.0-sweby_limiter) +                 & 
                  max(0.0, min(min(1.0, d0 + d1 * thetaM),   &
                  (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) ) &
                  *sweby_limiter 

           flux_y(i,j,k) =  0.5 * ( ( massflux + abs(massflux) )  &
                * ( tracer_mdfl(i,j,k) + psiP * Rj )              &
                + ( massflux - abs(massflux) )                    &
                * ( tracer_mdfl(i,j+1,k) - psiM * Rj ) )          &
                * tmask_mdfl(i,j,k) * tmask_mdfl(i,j+1,k)

        enddo
     enddo

     ! update tracer
     do j=jsc,jec
        do i=isc,iec

           mass_mdfl(i,j,k) = mass_mdfl(i,j,k) + dtime * (         &
                  Grd%dxtn(i,j-1) * Adv_vel%vhrho_nt(i,j-1,k)      &
                - Grd%dxtn(i,j) * Adv_vel%vhrho_nt(i,j,k) )

           tracermass_mdfl(i,j,k) = tracermass_mdfl(i,j,k) + dtime * ( flux_y(i,j-1,k) - flux_y(i,j,k) )

           if (mass_mdfl(i,j,k)>0.) tracer_mdfl(i,j,k) = tracermass_mdfl(i,j,k) / mass_mdfl(i,j,k)

        enddo
     enddo

    ! calculate the overall tendency
    do j=jsc,jec
      do i=isc,iec

        advect_tracer_mdfl_sweby_test(i,j,k) = tmask_mdfl(i,j,k) * (   &
                    Thickness%rho_dzt(i,j,k,tau) * Tracer_field(i,j,k) &
                     - tracermass_mdfl(i,j,k) * Grd%datr(i,j) ) / dtime

      enddo 
    enddo 

  enddo ! end of k-loop

!aja call tracer_stats(tracer_mdfl(isc:iec,jsc:jec,:),Tmin0,Tmax0,'')

  if(sweby_limiter==1.0) then
     call mpp_clock_end(id_clock_mdfl_sweby_test)
  else 
     call mpp_clock_end(id_clock_dst_linear_test)
  endif 
  
end function advect_tracer_mdfl_sweby_test
! </FUNCTION> NAME="advect_tracer_mdfl_sweby_test"



!#######################################################################
! <FUNCTION NAME="advect_tracer_mdfl_sweby">
!
! <DESCRIPTION>
!
!-------
! NOTE: This is a legacy version of the advect_tracer_mdfl_sweby routine. 
! This version may have problems with limiters that are improperly implemented
! in the case of generalized vertical level coordinates of MOM. Consequently,  
! under certain circumstances, the resulting tracer concentration can 
! exhibit non-monotonic behaviour (i.e., produce extrema).  
! This routine was used in various GFDL model configurations. 
!
! Sept 2011 
! Stephen.Griffies
! Alistair.Adcroft
!-------
!
! Compute tendency due to 3D advection of tracers using a multi-dimensional flux-limited
! method.  This method differs from other methods in the following ways:
!
! 1) Horizontal and vertical advection are combined
! 2) Calculations of the three coordinates (Z, X, Y) are performed sequentially as updates
!        to the tracer, so that the advection components for X and Y depend on Z and Z,Y
!        respectively... This helps limit the flux.
! 3) During the update for each direction, the 3rd order Sweby flux limiter is applied.
! 4) Flux divergence is included within the calculation to also help limit the flux:
!          - During the update for each direction, the divergence in each direction is added.
!          - During the overall tendency calculated, the divergence in all three directions
!                is removed.
! 5) All fluxes are functions of the tracer field at the taum1 time step.  This 
!    means that this method is ideally suited for the twolevel time stepping scheme,
!    in which "taum1=tau", thus enabling twice the tracer time step available for the 
!    threelevel scheme. 
!
! The calculation proceeds as follows:
!
! IMPORTANT NOTE: If this scheme is used at all, it must be used as the option for BOTH
! horizontal and vertical advection.  In the tracer tendency, it is applied as the 
! horizontal term, but applies to vertical as well, for which case the vertical term in
! the tracer tendency equation is set to zero.
!
! This scheme was ported to mom4 from the MIT-GCM by John Dunne and Alistair Adcroft  
! during Summer 2003 
!
! Griffies: 5/27/04
! Optimized by filling 3d arrays prior to sending for mpp_update_domains
! 
! 07/11/2007 by Stephen.Griffies
! When the nonlinear limiters are removed (sweby_limiter=0.0), 
! this scheme reduces to a linear direct space-time method 
! (ADVECT_DST_LINEAR).  
!
! </DESCRIPTION>
!
function advect_tracer_mdfl_sweby(Time, Adv_vel, Thickness, Tracer_field, dtime, sweby_limiter)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field
  real,                         intent(in) :: dtime
  real,                         intent(in) :: sweby_limiter

  real,dimension(isc:iec,jsc:jec)    :: ftp
  real,dimension(isc:iec,jsc:jec)    :: fbt
  real,dimension(isc:iec,jsc:jec)    :: wkm1
  real,dimension(isc:iec,jsc:jec,nk) :: advect_tracer_mdfl_sweby

  integer  :: i, j, k
  integer  :: kp1, kp2, km1
  integer  :: tau, taum1
  real     :: Rjm, Rj, Rjp, cfl, massflux
  real     :: d0, d1, thetaP, psiP 
  real     :: thetaM, psiM

  if(sweby_limiter==1.0) then
     call mpp_clock_begin(id_clock_mdfl_sweby)
  else 
     call mpp_clock_begin(id_clock_dst_linear)
  endif 

  ftp  = 0.0
  fbt  = 0.0
  wkm1 = 0.0
  tracer_mdfl = 0.0
  flux_x      = 0.0
  flux_y      = 0.0

  tau   = Time%tau
  taum1 = Time%taum1

     do k=1,nk

        do j=jsc,jec
           do i=isc,iec
           tracer_mdfl(i,j,k) = Tracer_field(i,j,k)
           enddo
        enddo

        kp1 = min(k+1,nk)
        kp2 = min(k+2,nk)
        km1 = max(k-1,1)

     ! calculate flux at bottom of box and update the tracer
        do j=jsc,jec
           do i=isc,iec

           Rjp = (Tracer_field(i,j,km1) - Tracer_field(i,j,k))    &
                   *tmask_mdfl(i,j,km1)*tmask_mdfl(i,j,k)
           Rj  = (Tracer_field(i,j,k) - Tracer_field(i,j,kp1))    &
                   *tmask_mdfl(i,j,k)*tmask_mdfl(i,j,kp1)
           Rjm = (Tracer_field(i,j,kp1) - Tracer_field(i,j,kp2))  &
                   *tmask_mdfl(i,j,kp1)*tmask_mdfl(i,j,kp2)

              massflux = Grd%dat(i,j) * Adv_vel%wrho_bt(i,j,k)
              cfl = abs(Adv_vel%wrho_bt(i,j,k) * dtime  / Thickness%rho_dzt(i,j,k,tau))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

           psiP = d0 + d1 * thetaP
           psiP = psiP*(1.0-sweby_limiter) +                 &
                  max(0.0, min(min(1.0, d0 + d1 * thetaP),   &
                  (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) ) &
                  *sweby_limiter

           psiM = d0 + d1 * thetaM
           psiM = psiM*(1.0-sweby_limiter) +                 &
                  max(0.0, min(min(1.0, d0 + d1 * thetaM),   &
                  (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) ) &
                  *sweby_limiter  

              fbt(i,j) =  0.5 * ( ( massflux + abs(massflux) )   &
                * ( Tracer_field(i,j,kp1) + psiP * Rj )          &  
                   + ( massflux - abs(massflux) )                &
                * ( Tracer_field(i,j, k ) - psiM * Rj ) )        &
                   * tmask_mdfl(i,j,kp1) * tmask_mdfl(i,j,k)

           tracer_mdfl(i,j,k) = tracer_mdfl(i,j,k)               &
                + dtime / Thickness%rho_dzt(i,j,k,tau) * (       &
                Grd%datr(i,j) * ( fbt(i,j) - ftp(i,j))           &
                + Tracer_field(i,j,k) *                          &
                (wkm1(i,j) - Adv_vel%wrho_bt(i,j,k)) )

              flux_z(i,j,k) = fbt(i,j)
              ftp(i,j)      = fbt(i,j)

           enddo
        enddo

        ! update vertical velocity for next k-level
        do j=jsc,jec
           do i=isc,iec
              wkm1(i,j) = Adv_vel%wrho_bt(i,j,k)
           enddo
        enddo

     enddo ! end of k-loop
  call mpp_update_domains (tracer_mdfl, Dom_mdfl%domain2d, flags=XUPDATE)


  ! calculate flux at the eastern wall of the boxes
     do k=1,nk

        do j=jsc,jec
           do i=isc-1,iec

           Rjp = (tracer_mdfl(i+2,j,k) - tracer_mdfl(i+1,j,k))    &
                   *tmask_mdfl(i+2,j,k)*tmask_mdfl(i+1,j,k)
           Rj  = (tracer_mdfl(i+1,j,k) - tracer_mdfl( i ,j,k) )   &
                   *tmask_mdfl(i+1,j,k)*tmask_mdfl(i,j,k)
           Rjm = (tracer_mdfl(i,j,k) - tracer_mdfl(i-1,j,k))      &
                   *tmask_mdfl(i,j,k)*tmask_mdfl(i-1,j,k)

              massflux = Grd%dyte(i,j) * Adv_vel%uhrho_et(i,j,k)
              cfl = abs(Adv_vel%uhrho_et(i,j,k) * dtime * 2.0    &
                   / ((Thickness%rho_dzt(i,j,k,tau) + Thickness%rho_dzt(i+1,j,k,tau)) * Grd%dxte(i,j)))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

           psiP = d0 + d1 * thetaP
           psiP = psiP*(1.0-sweby_limiter) +                 & 
                  max(0.0, min(min(1.0, d0 + d1 * thetaP),   &
                  (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) ) &
                  *sweby_limiter

           psiM = d0 + d1 * thetaM
           psiM = psiM*(1.0-sweby_limiter) +                 & 
                  max(0.0, min(min(1.0, d0 + d1 * thetaM),   &
                  (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) ) &
                  *sweby_limiter

              flux_x(i,j,k) =  0.5 * ( ( massflux + abs(massflux) ) &
                * ( tracer_mdfl( i ,j,k) + psiP * Rj )              &
                   + ( massflux - abs(massflux) )                   &
                * ( tracer_mdfl(i+1,j,k) - psiM * Rj ) )            &
                   * tmask_mdfl(i,j,k) * tmask_mdfl(i+1,j,k)
           enddo
        enddo

        ! update the tracer
        do j=jsc,jec
           do i=isc,iec
           tracer_mdfl(i,j,k) = tracer_mdfl(i,j,k)                                   &
          + dtime * tmask_mdfl(i,j,k) * Grd%datr(i,j) / Thickness%rho_dzt(i,j,k,tau) &
                   * (flux_x(i-1,j,k) - flux_x(i,j,k)                                &
                + Tracer_field(i,j,k)* (                                             &
                     Grd%dyte(i,j)   * Adv_vel%uhrho_et( i ,j,k)                     &
                   - Grd%dyte(i-1,j) * Adv_vel%uhrho_et(i-1,j,k) ) )
           enddo
        enddo

     enddo ! end of k-loop
  call mpp_update_domains (tracer_mdfl, Dom_mdfl%domain2d, flags=YUPDATE)


  ! calculate flux at the northern wall of the boxes
     do j=jsc,jec
        do i=isc,iec
           wkm1(i,j) = 0.0
        enddo
     enddo

     do k=1,nk

        do j=jsc-1,jec
           do i=isc,iec

           Rjp = (tracer_mdfl(i,j+2,k) - tracer_mdfl(i,j+1,k))  &
                   *tmask_mdfl(i,j+2,k)*tmask_mdfl(i,j+1,k)
           Rj  = (tracer_mdfl(i,j+1,k) - tracer_mdfl(i,j,k))    &
                   *tmask_mdfl(i,j+1,k)*tmask_mdfl(i,j,k)
           Rjm = (tracer_mdfl(i,j,k) - tracer_mdfl(i,j-1,k))    &
                   *tmask_mdfl(i,j,k)*tmask_mdfl(i,j-1,k)

              massflux = Grd%dxtn(i,j) * Adv_vel%vhrho_nt(i,j,k)
              cfl = abs(Adv_vel%vhrho_nt(i,j,k) * dtime * 2.0              &
                   / ((Thickness%rho_dzt(i,j,k,tau) + Thickness%rho_dzt(i,j+1,k,tau)) * Grd%dytn(i,j)))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

           psiP = d0 + d1 * thetaP
           psiP = psiP*(1.0-sweby_limiter) +                 & 
                  max(0.0, min(min(1.0, d0 + d1 * thetaP),   &
                  (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) ) &
                  *sweby_limiter

           psiM = d0 + d1 * thetaM
           psiM = psiM*(1.0-sweby_limiter) +                 & 
                  max(0.0, min(min(1.0, d0 + d1 * thetaM),   &
                  (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) ) &
                  *sweby_limiter 

              flux_y(i,j,k) =  0.5 * ( ( massflux + abs(massflux) ) &
                * ( tracer_mdfl(i,j,k) + psiP * Rj )                &
                   + ( massflux - abs(massflux) )                   &
                * ( tracer_mdfl(i,j+1,k) - psiM * Rj ) )            &
                   * tmask_mdfl(i,j,k) * tmask_mdfl(i,j+1,k)

           enddo
        enddo

     ! update tracer
        do j=jsc,jec
           do i=isc,iec
           tracer_mdfl(i,j,k) = tracer_mdfl(i,j,k)                                          &
                + dtime * tmask_mdfl(i,j,k) * Grd%datr(i,j) / Thickness%rho_dzt(i,j,k,tau)  &
                * (flux_y(i,j-1,k) - flux_y(i,j,k))
           enddo
        enddo

    ! calculate the overall tendency
        do j=jsc,jec
           do i=isc,iec

        tracer_mdfl(i,j,k) = tracer_mdfl(i,j,k)                               &
           + dtime * Tracer_field(i,j,k) / Thickness%rho_dzt(i,j,k,tau) *   ( &
                      Adv_vel%wrho_bt(i,j,k) - wkm1(i,j)                      &
             + Grd%datr(i,j)*(                                                &
                         ( Grd%dyte(i-1,j) * Adv_vel%uhrho_et(i-1,j,k)        &
                         - Grd%dyte( i ,j) * Adv_vel%uhrho_et( i ,j,k))))

        ! advection is applied as a convergence
        advect_tracer_mdfl_sweby(i,j,k) = -Thickness%rho_dzt(i,j,k,tau) *  &
          ( tracer_mdfl(i,j,k) - Tracer_field(i,j,k) )                            &
            / dtime * tmask_mdfl(i,j,k)
          enddo
        enddo

        ! update vertical velocity for next level
        do j=jsc,jec
           do i=isc,iec
              wkm1(i,j) = Adv_vel%wrho_bt(i,j,k)
           enddo
        enddo

     enddo ! end of k-loop


  if(sweby_limiter==1.0) then
     call mpp_clock_end(id_clock_mdfl_sweby)
  else 
     call mpp_clock_end(id_clock_dst_linear)
     endif


end function advect_tracer_mdfl_sweby
! </FUNCTION> NAME="advect_tracer_mdfl_sweby"



!#######################################################################
! <SUBROUTINE NAME="advect_tracer_sweby_all">
!
! <DESCRIPTION>
!
!-------
! NOTE: This is routine suffers from the same problems as the 
! advect_tracer_mdfl_sweby routine. Namely, it has problems  
! with limiters that are improperly implemented in the case of generalized 
! vertical level coordinates of MOM. Consequently,  
! under certain circumstances, the resulting tracer concentration can 
! exhibit non-monotonic behaviour (i.e., produce extrema).  
! This routine was used in various GFDL model configurations. It is retained 
! solely to provide bitwise reproducing those earlier configurations. 
!
! Sept 2011 
! Stephen.Griffies
! Alistair.Adcroft
!-------
!
! Sweby scheme optimized by doing all tracers at once, so 
! can send larger packet to mpp_update_domains.
!
! This scheme is available ONLY when advecting all tracers with the 
! sweby scheme.
!
! This scheme CANNOT be called from compute_adv_diss.  
! 
! Stephen.Griffies
! June 2004
!
! </DESCRIPTION>
!
subroutine advect_tracer_sweby_all(Time, Adv_vel, Dens, T_prog, Thickness, dtime)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_thickness_type),   intent(in)    :: Thickness
  real,                         intent(in)    :: dtime

  real,dimension(isc:iec,jsc:jec)             :: ftp
  real,dimension(isc:iec,jsc:jec)             :: fbt
  real,dimension(isc:iec,jsc:jec)             :: wkm1
  real,dimension(isd:ied,jsd:jed)             :: tmp_flux

  integer                                     :: i, j, k, n
  integer                                     :: kp1, kp2, km1
  integer                                     :: tau, taum1
  real                                        :: Rjm, Rj, Rjp, cfl, massflux
  real                                        :: d0, d1, thetaP, psiP 
  real                                        :: thetaM, psiM
  integer, dimension(num_prog_tracers)        :: id_update

  call mpp_clock_begin(id_clock_mdfl_sweby_all)

  tau     = Time%tau
  taum1   = Time%taum1
  flux_x  = 0.0
  flux_y  = 0.0
  flux_z  = 0.0
  wrk1    = 0.0

  ! calculate flux at bottom face of the T-cells
  do n=1,num_prog_tracers

     ftp  = 0.0
     fbt  = 0.0
     wkm1 = 0.0

         do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              T_prog(n)%wrk1(i,j,k) = 0.0
               enddo
            enddo
         enddo

     do k=1,nk

        do j=jsc,jec
           do i=isc,iec
              tracer_mdfl_all(n)%field(i,j,k) = T_prog(n)%field(i,j,k,taum1)
           enddo
        enddo

        kp1 = min(k+1,nk)
        kp2 = min(k+2,nk)
        km1 = max(k-1,1)
        call mpp_clock_begin(id_clock_mdfl_sweby_cu1)
        do j=jsc,jec
           do i=isc,iec

              Rjp = (T_prog(n)%field(i,j,km1,taum1) - T_prog(n)%field(i,j,k,taum1))    &
                   *tmask_mdfl(i,j,km1)*tmask_mdfl(i,j,k)
              Rj  = (T_prog(n)%field(i,j,k,taum1) - T_prog(n)%field(i,j,kp1,taum1))    &
                   *tmask_mdfl(i,j,k)*tmask_mdfl(i,j,kp1)
              Rjm = (T_prog(n)%field(i,j,kp1,taum1) - T_prog(n)%field(i,j,kp2,taum1))  &
                   *tmask_mdfl(i,j,kp1)*tmask_mdfl(i,j,kp2)

              massflux = Grd%dat(i,j) * Adv_vel%wrho_bt(i,j,k)
              cfl = abs(Adv_vel%wrho_bt(i,j,k) * dtime  / Thickness%rho_dzt(i,j,k,tau))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

              psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),       &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
              psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),       &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

              fbt(i,j) =  0.5 * ( ( massflux + abs(massflux) )         &
                   * ( T_prog(n)%field(i,j,kp1,taum1) + psiP * Rj )    &
                   + ( massflux - abs(massflux) )                      &
                   * ( T_prog(n)%field(i,j, k ,taum1) - psiM * Rj ) )  &
                   * tmask_mdfl(i,j,kp1) * tmask_mdfl(i,j,k)

              wrk1(i,j,k) = Grd%datr(i,j)*( fbt(i,j) - ftp(i,j))     &
                   + T_prog(n)%field(i,j,k,taum1)*(wkm1(i,j) - Adv_vel%wrho_bt(i,j,k)) 

              tracer_mdfl_all(n)%field(i,j,k) = tracer_mdfl_all(n)%field(i,j,k)  &
                   + wrk1(i,j,k) * dtime/Thickness%rho_dzt(i,j,k,tau) 

              flux_z(i,j,k) = fbt(i,j)
              ftp(i,j)      = fbt(i,j)

           enddo
        enddo
        call mpp_clock_end(id_clock_mdfl_sweby_cu1)        

        ! update vertical velocity for next k-level
        do j=jsc,jec
           do i=isc,iec
              wkm1(i,j) = Adv_vel%wrho_bt(i,j,k)
           enddo
        enddo

     enddo ! end of k-loop

     call mpp_clock_begin(id_clock_mdfl_sweby_mpi)
     if (async_domain_update) then        
        id_update(n) = mpp_start_update_domains (tracer_mdfl_all(n)%field(:,:,:), Dom_mdfl%domain2d, &
                                                 flags=XUPDATE)
     else
        id_update(1) = mpp_start_update_domains (tracer_mdfl_all(n)%field(:,:,:), Dom_mdfl%domain2d, &
                                                 flags=XUPDATE, complete=T_prog(n)%complete)
     endif
     call mpp_clock_end(id_clock_mdfl_sweby_mpi)

     call mpp_clock_begin(id_clock_mdfl_sweby_dia)
     if (id_zflux_adv(n) > 0) call diagnose_3d(Time, Grd, id_zflux_adv(n), &
         T_prog(n)%conversion*flux_z(:,:,:))

     if (id_advection_z(n) > 0) call diagnose_3d(Time, Grd, id_advection_z(n), &
         T_prog(n)%conversion*wrk1(:,:,:))
     call mpp_clock_end(id_clock_mdfl_sweby_dia)
    
  enddo ! end of n-loop for tracers  

  if( .not. async_domain_update) then
     call mpp_clock_begin(id_clock_mdfl_sweby_mpi)
     do n=1,num_prog_tracers
        call mpp_complete_update_domains (id_update(1), tracer_mdfl_all(n)%field(:,:,:), Dom_mdfl%domain2d, &
                                          flags=XUPDATE, complete=T_prog(n)%complete)
     enddo
     call mpp_clock_end(id_clock_mdfl_sweby_mpi)
  endif

  ! calculate flux at the eastern face of the T-cells
  do n=1,num_prog_tracers
     if(async_domain_update) then
        call mpp_clock_begin(id_clock_mdfl_sweby_mpi)
        call mpp_complete_update_domains (id_update(n), tracer_mdfl_all(n)%field(:,:,:), Dom_mdfl%domain2d, &
                                          flags=XUPDATE)
        call mpp_clock_end(id_clock_mdfl_sweby_mpi)
     endif
     call mpp_clock_begin(id_clock_mdfl_sweby_cu2)
     do k=1,nk

        do j=jsc,jec
           do i=isc-1,iec

              Rjp = (tracer_mdfl_all(n)%field(i+2,j,k) - tracer_mdfl_all(n)%field(i+1,j,k))    &
                   *tmask_mdfl(i+2,j,k)*tmask_mdfl(i+1,j,k)
              Rj  = (tracer_mdfl_all(n)%field(i+1,j,k) - tracer_mdfl_all(n)%field( i ,j,k) )   &
                   *tmask_mdfl(i+1,j,k)*tmask_mdfl(i,j,k)
              Rjm = (tracer_mdfl_all(n)%field(i,j,k) - tracer_mdfl_all(n)%field(i-1,j,k))      &
                   *tmask_mdfl(i,j,k)*tmask_mdfl(i-1,j,k)

              massflux = Grd%dyte(i,j) * Adv_vel%uhrho_et(i,j,k)
              cfl = abs(Adv_vel%uhrho_et(i,j,k) * dtime * 2.0    &
                   / ((Thickness%rho_dzt(i,j,k,tau) + Thickness%rho_dzt(i+1,j,k,tau)) * Grd%dxte(i,j)))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

              psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),  &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
              psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),  &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

              flux_x(i,j,k) =  0.5 * ( ( massflux + abs(massflux) )      &
                   * ( tracer_mdfl_all(n)%field( i ,j,k) + psiP * Rj )   &
                   + ( massflux - abs(massflux) )                        &
                   * ( tracer_mdfl_all(n)%field(i+1,j,k) - psiM * Rj ) ) &
                   * tmask_mdfl(i,j,k) * tmask_mdfl(i+1,j,k)
           enddo
        enddo

        ! update the tracer
        do j=jsc,jec
           do i=isc,iec
              wrk1(i,j,k) = tmask_mdfl(i,j,k) * Grd%datr(i,j)       &
                   * (flux_x(i-1,j,k) - flux_x(i,j,k)               &
                   + T_prog(n)%field(i,j,k,taum1)* (                &
                     Grd%dyte(i,j)   * Adv_vel%uhrho_et( i ,j,k)    &
                   - Grd%dyte(i-1,j) * Adv_vel%uhrho_et(i-1,j,k) ) )

              tracer_mdfl_all(n)%field(i,j,k) = tracer_mdfl_all(n)%field(i,j,k) &
                   + wrk1(i,j,k)*dtime/Thickness%rho_dzt(i,j,k,tau)
           enddo
        enddo

     enddo ! end of k-loop
     call mpp_clock_end(id_clock_mdfl_sweby_cu2)

     call mpp_clock_begin(id_clock_mdfl_sweby_mpi)
     if (async_domain_update) then
        id_update(n) = mpp_start_update_domains (tracer_mdfl_all(n)%field(:,:,:), Dom_mdfl%domain2d, &
                                        flags=YUPDATE)
     else 
        id_update(1) = mpp_start_update_domains (tracer_mdfl_all(n)%field(:,:,:), Dom_mdfl%domain2d, &
                                 flags=YUPDATE, complete=T_prog(n)%complete)
     endif
     call mpp_clock_end(id_clock_mdfl_sweby_mpi)

     call mpp_clock_begin(id_clock_mdfl_sweby_dia)
     if (id_xflux_adv(n) > 0) call diagnose_3d(Time, Grd, id_xflux_adv(n), T_prog(n)%conversion*flux_x(:,:,:))

     if (id_advection_x(n) > 0) call diagnose_3d(Time, Grd, id_advection_x(n), T_prog(n)%conversion*wrk1(:,:,:))

     if (id_xflux_adv_int_z(n) > 0) then 
         tmp_flux(:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  tmp_flux(i,j) = tmp_flux(i,j) +  flux_x(i,j,k) 
               enddo
            enddo
         enddo
         call diagnose_2d(Time, Grd, id_xflux_adv_int_z(n), T_prog(n)%conversion*tmp_flux(:,:))
     endif
     call mpp_clock_end(id_clock_mdfl_sweby_dia)

     ! this is only needed for diagnostics of the flux through open boundaries.
     if (have_obc) then
        call store_ocean_obc_tracer_flux(Time,T_prog(n),flux_x,n,'z','adv')
     endif
  enddo ! end of n-loop for tracers  

  if(.not. async_domain_update) then
     call mpp_clock_begin(id_clock_mdfl_sweby_mpi)
     do n=1,num_prog_tracers
        call mpp_complete_update_domains (id_update(1), tracer_mdfl_all(n)%field(:,:,:), Dom_mdfl%domain2d, &
                                          flags=YUPDATE, complete=T_prog(n)%complete)
     enddo
     call mpp_clock_end(id_clock_mdfl_sweby_mpi)
  endif


  ! calculate flux at the northern face of the T-cells 
  do n=1,num_prog_tracers 

     if (async_domain_update) then
        call mpp_clock_begin(id_clock_mdfl_sweby_mpi)
        call mpp_complete_update_domains (id_update(n), tracer_mdfl_all(n)%field(:,:,:), Dom_mdfl%domain2d, &
                              flags=YUPDATE)
        call mpp_clock_end(id_clock_mdfl_sweby_mpi)
     endif

     do j=jsc,jec
        do i=isc,iec
           wkm1(i,j) = 0.0
        enddo
     enddo

     do k=1,nk
        call mpp_clock_begin(id_clock_mdfl_sweby_cu3)
        do j=jsc-1,jec
           do i=isc,iec

              Rjp = (tracer_mdfl_all(n)%field(i,j+2,k) - tracer_mdfl_all(n)%field(i,j+1,k))   &
                   *tmask_mdfl(i,j+2,k)*tmask_mdfl(i,j+1,k)
              Rj  = (tracer_mdfl_all(n)%field(i,j+1,k) - tracer_mdfl_all(n)%field(i,j,k))     &
                   *tmask_mdfl(i,j+1,k)*tmask_mdfl(i,j,k)
              Rjm = (tracer_mdfl_all(n)%field(i,j,k)   - tracer_mdfl_all(n)%field(i,j-1,k))   &
                   *tmask_mdfl(i,j,k)*tmask_mdfl(i,j-1,k)

              massflux = Grd%dxtn(i,j) * Adv_vel%vhrho_nt(i,j,k)
              cfl = abs(Adv_vel%vhrho_nt(i,j,k) * dtime * 2.0              &
                   / ((Thickness%rho_dzt(i,j,k,tau) + Thickness%rho_dzt(i,j+1,k,tau)) * Grd%dytn(i,j)))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

              psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),   &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
              psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),   &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

              flux_y(i,j,k) =  0.5 * ( ( massflux + abs(massflux) )       &
                   * ( tracer_mdfl_all(n)%field(i,j,k) + psiP * Rj )      &
                   + ( massflux - abs(massflux) )                         &
                   * ( tracer_mdfl_all(n)%field(i,j+1,k) - psiM * Rj ) )  &
                   * tmask_mdfl(i,j,k) * tmask_mdfl(i,j+1,k)

           enddo
        enddo

        ! calculate the overall tendency and update tracer 
        do j=jsc,jec
           do i=isc,iec

              wrk1(i,j,k) = tmask_mdfl(i,j,k)*Grd%datr(i,j)*(flux_y(i,j-1,k)-flux_y(i,j,k))  &
                            + T_prog(n)%field(i,j,k,taum1) *                                 &
                   (Adv_vel%wrho_bt(i,j,k) - wkm1(i,j)                                       &
                   + Grd%datr(i,j)*                                                          &
                   ( Grd%dyte(i-1,j) * Adv_vel%uhrho_et(i-1,j,k)                             &
                   - Grd%dyte( i ,j) * Adv_vel%uhrho_et( i ,j,k)))

              tracer_mdfl_all(n)%field(i,j,k) = tracer_mdfl_all(n)%field(i,j,k) &
                   + wrk1(i,j,k)*dtime/Thickness%rho_dzt(i,j,k,tau)

              T_prog(n)%wrk1(i,j,k) = &
                Thickness%rho_dzt(i,j,k,tau)*(tracer_mdfl_all(n)%field(i,j,k)-T_prog(n)%field(i,j,k,taum1))/dtime  &
                   *tmask_mdfl(i,j,k) 

           enddo
        enddo
        call mpp_clock_end(id_clock_mdfl_sweby_cu3)

        if (have_obc) call ocean_obc_zero_boundary(T_prog(n)%wrk1(:,:,k), "T")
        do j=jsc,jec
           do i=isc,iec
              T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
                   enddo
                enddo

        ! update vertical velocity for next level
                    do j=jsc,jec
                       do i=isc,iec
              wkm1(i,j) = Adv_vel%wrho_bt(i,j,k)
                       enddo
                    enddo
     enddo ! end of k-loop

     call mpp_clock_begin(id_clock_mdfl_sweby_dia)
     ! this is only needed for diagnostics of the flux through open boundaries.
     if (have_obc) then
        call store_ocean_obc_tracer_flux(Time,T_prog(n),flux_y,n,'m', 'adv')
     endif


     ! diagnostics 

     if (id_yflux_adv(n) > 0) call diagnose_3d(Time, Grd, id_yflux_adv(n), &
                      T_prog(n)%conversion*flux_y(:,:,:))

     if (id_advection_y(n) > 0) call diagnose_3d(Time, Grd, id_advection_y(n), &
                      T_prog(n)%conversion*wrk1(:,:,:))

     if (id_yflux_adv_int_z(n) > 0) then 
         tmp_flux(:,:) = 0.0
                 do k=1,nk
                    do j=jsc,jec
                       do i=isc,iec
                         tmp_flux(i,j) = tmp_flux(i,j) +  flux_y(i,j,k) 
                       enddo
                    enddo
                 enddo
                 call diagnose_2d(Time, Grd, id_yflux_adv_int_z(n), T_prog(n)%conversion*tmp_flux(:,:))
     endif

     if (id_sweby_advect(n) > 0) then 
        call diagnose_3d(Time, Grd, id_sweby_advect(n), T_prog(n)%conversion*T_prog(n)%wrk1(:,:,:))
     endif

     advect_tendency(:,:,:) = 0.0
     do k=1,nk
         do j=jsc,jec
             do i=isc,iec
                 advect_tendency(i,j,k) = T_prog(n)%wrk1(i,j,k)
             enddo
         enddo
     enddo

     ! for watermass_diag 
     if(n==index_temp) then 
         neutral_temp_advect(:,:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                   neutral_temp_advect(i,j,k) = T_prog(n)%wrk1(i,j,k)
                enddo
            enddo
         enddo
     endif
     if(n==index_salt) then 
         neutral_salt_advect(:,:,:) = 0.0
         do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   neutral_salt_advect(i,j,k) = T_prog(n)%wrk1(i,j,k)
                 enddo
              enddo
          enddo
      endif

      ! This is only needed for diagnostics of the flux through open boundaries.
      ! It must be checked, if still correct with the new scheme
      ! for gyre/overturning diagnostics
      if(compute_gyre_overturn_diagnose) then
         call gyre_overturn_diagnose(Time, Adv_vel, T_prog(n), n)
      endif

      call mpp_clock_end(id_clock_mdfl_sweby_dia)
  enddo ! end of tracer n-loop

  call watermass_diag(Time, Dens)
  
  call mpp_clock_end(id_clock_mdfl_sweby_all)


end subroutine advect_tracer_sweby_all
! </SUBROUTINE> NAME="advect_tracer_sweby_all"



!#######################################################################
! <FUNCTION NAME="advect_tracer_psom">
!
! <DESCRIPTION>
!
! Compute advective tendency of dzt*rho*tracer concentration using 
! Prather's second order moment method (SOM): 
!
! Prather, M. J.,"Numerical Advection by Conservation of Second-Order Moments"
! JGR, Vol 91, NO. D6, p 6671-6681, May 20, 1986
!
! Merryfield and Holloway (2003), "Application of an accurate advection 
! algorithm to sea-ice modelling". Ocean Modelling, Vol 5, p 1-15.
!
! MOM 3 code by M. Hofmann and M. A. M. Maqueda
! Ported to mom4p0 by Ronald.Pacanowski
! Ported to mom4p1 by Stephen.Griffies
!
! The preferred limiters are taken from Merryfield and Holloway,
! and generalized by Bill Merryfield for non-constant grid spacing.
!
! IMPORTANT NOTE: This scheme must be used for BOTH horizontal and
! vertical advection.  In the tracer tendency, it is applied as the 
! horizontal term, but applies to vertical as well, for which case the 
! vertical term in the tracer tendency equation is set to zero.
!
! NOTE: When using psom_limit_prather=.true., the tracer has a lower 
! bound of zero.  So this limiter IS NOT appropriate for temperature 
! when measured in degrees C. The preferred, and default, limiter
! is from Merryfield and Holloway.  
!
! </DESCRIPTION>
!
function advect_tracer_psom(Time, Adv_vel, Tracer, Thickness, dtime, ntracer, do_passive_sq)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_prog_tracer_type), intent(inout) :: Tracer
  type(ocean_thickness_type),   intent(in)    :: Thickness
  real,                         intent(in)    :: dtime
  integer,                      intent(in)    :: ntracer 
  logical,                      intent(in)    :: do_passive_sq 

  real, dimension(isc:iec,jsc:jec,nk)         :: advect_tracer_psom
  real, dimension(isc:iec,jsc:jec)            :: ft1

  integer :: i, j, k, lag

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
     '==>Error from ocean_tracer_advect (advect_tracer_psom): ocean_tracer_advect_mod not yet initialized')
  endif 

  if(.not. do_passive_sq) then 
      if(ntracer==index_temp_sq .or. ntracer==index_salt_sq) then 
          advect_tracer_psom=0.0
          return
      endif
  endif

  call mpp_clock_begin(id_clock_psom)

  ! needed for reproducibility across restarts 
  ! and when use different advection schemes. 
  flux_z = 0

  lag = Time%taum1
  do k=1,nk
    do j=jsd,jed
      do i=isd,ied
        sm(i,j,k)        = Thickness%rho_dzt(i,j,k,lag)*Grd%dat(i,j)
        Tracer%s0(i,j,k) = Tracer%field(i,j,k,lag)*sm(i,j,k)*Grd%tmask(i,j,k)
      enddo
    enddo
  enddo

  ! alternate sequential solving in one dimension to improve accuracy in time
  if (mod(Time%itt0,6).eq.0) then
    call psom_x (Adv_vel, Tracer, dtime, lag)
    call psom_y (Adv_vel, Tracer, dtime, lag)
    call psom_z (Adv_vel, Tracer, dtime, lag)
  endif
  if (mod(Time%itt0,6).eq.1) then
    call psom_z (Adv_vel, Tracer, dtime, lag)
    call psom_x (Adv_vel, Tracer, dtime, lag)
    call psom_y (Adv_vel, Tracer, dtime, lag)
  endif
  if (mod(Time%itt0,6).eq.2) then
    call psom_y (Adv_vel, Tracer, dtime, lag)
    call psom_z (Adv_vel, Tracer, dtime, lag)
    call psom_x (Adv_vel, Tracer, dtime, lag)
  endif
  if (mod(Time%itt0,6).eq.3) then
    call psom_x (Adv_vel, Tracer, dtime, lag)
    call psom_z (Adv_vel, Tracer, dtime, lag)
    call psom_y (Adv_vel, Tracer, dtime, lag)
  endif
  if (mod(Time%itt0,6).eq.4) then
    call psom_y (Adv_vel, Tracer, dtime, lag)
    call psom_x (Adv_vel, Tracer, dtime, lag)
    call psom_z (Adv_vel, Tracer, dtime, lag)
  endif
  if (mod(Time%itt0,6).eq.5) then
    call psom_z (Adv_vel, Tracer, dtime, lag)
    call psom_y (Adv_vel, Tracer, dtime, lag)
    call psom_x (Adv_vel, Tracer, dtime, lag)
  endif


  ! to handle redundancies across the bipolar fold 
  if (Grd%tripolar) call mpp_update_domains(flux_x, flux_y, Dom_flux%domain2d, gridtype=CGRID_NE) 

  ft1(:,:) = 0.0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        advect_tracer_psom(i,j,k) = Grd%tmask(i,j,k)*(flux_x(i,j,k) - flux_x(i-1,j,k)&
                                                   +  flux_y(i,j,k) - flux_y(i,j-1,k)&
                                                   +  ft1(i,j)-flux_z(i,j,k))*Grd%datr(i,j)
      enddo
    enddo
    do j=jsc,jec
      do i=isc,iec
        ft1(i,j) = flux_z(i,j,k)
      enddo
    enddo
  enddo

  if(debug_this_module) then 
     call tracer_psom_chksum(Time, Tracer)
  endif 
 
  call mpp_clock_end(id_clock_psom)

end function advect_tracer_psom
! </FUNCTION> NAME="advect_tracer_psom"



!#######################################################################
! <SUBROUTINE NAME="psom_x">
!
! <DESCRIPTION>
! Compute i-advective flux using Prather's SOM. 
! </DESCRIPTION>
!
subroutine psom_x (Adv_vel, Tracer, dtime, lag)

  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_prog_tracer_type), intent(inout) :: Tracer
  real,                         intent(in)    :: dtime
  integer,                      intent(in)    :: lag

  real, dimension(isd:ied,jsd:jed,nk) :: f0, fm, fx, fxx, fy, fyy, fz, fzz, fxy, fxz, fyz
  integer :: i, j, k, ip1
  real    :: slpmax, s1max, s1new, s2new, velocity
  real    :: s0max1, s0min1
  real    :: alf, alfq, alf1, alf1q, temptm
  real    :: dtimer

  call mpp_clock_begin(id_clock_psom_x)
  dtimer = 1.0/dtime

  ! when flux-limiting, limit appropriate moments before transport.
  if (Tracer%psom_limit) then

      if(psom_limit_prather) then 

          ! limiter from Prather
          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   slpmax = 0.0
                   if (Tracer%s0(i,j,k) > 0.0) slpmax = Tracer%s0(i,j,k)
                   s1max  = 1.5*slpmax
                   s1new  = min(s1max,max(-s1max,Tracer%sx(i,j,k)))
                   s2new  = min(slpmax+slpmax-0.3334*abs(s1new), max(abs(s1new)-slpmax,Tracer%sxx(i,j,k)))
                   Tracer%sxy(i,j,k) = min(slpmax,max(-slpmax,Tracer%sxy(i,j,k)))
                   Tracer%sxz(i,j,k) = min(slpmax,max(-slpmax,Tracer%sxz(i,j,k)))
                   Tracer%sx(i,j,k)  = s1new
                   Tracer%sxx(i,j,k) = s2new
                enddo
             enddo
          enddo

      else 

          ! limiter from Merryfield and Holloway, generalized to nonconstant grids
          do k=1,nk
             do j=jsd,jed
                do i=isc,iec

                   s0max1 = Tracer%field(i,j,k,lag)*Grd%tmask(i,j,k)
                   if(Grd%tmask(i-1,j,k) > 0.) s0max1 = max(Tracer%field(i-1,j,k,lag),s0max1)
                   if(Grd%tmask(i+1,j,k) > 0.) s0max1 = max(Tracer%field(i+1,j,k,lag),s0max1) 

                   s0max1 = s0max1*sm(i,j,k)
                   slpmax = -Tracer%s0(i,j,k) + s0max1
                   s1max  = 1.5*slpmax
                   s1new  = min(s1max,max(-s1max,-Tracer%sx(i,j,k)))
                   s2new  = min(slpmax+slpmax-0.3334*abs(s1new), max(abs(s1new)-slpmax,-Tracer%sxx(i,j,k)))
                   Tracer%sxy(i,j,k) = -min(slpmax,max(-slpmax,-Tracer%sxy(i,j,k)))
                   Tracer%sxz(i,j,k) = -min(slpmax,max(-slpmax,-Tracer%sxz(i,j,k)))
                   Tracer%sx(i,j,k)  = -s1new
                   Tracer%sxx(i,j,k) = -s2new

                   s0min1 = Tracer%field(i,j,k,lag)*Grd%tmask(i,j,k)
                   if(Grd%tmask(i-1,j,k) > 0.) s0min1 = min(Tracer%field(i-1,j,k,lag),s0min1)
                   if(Grd%tmask(i+1,j,k) > 0.) s0min1 = min(Tracer%field(i+1,j,k,lag),s0min1)  

                   s0min1 = s0min1*sm(i,j,k)
                   slpmax = Tracer%s0(i,j,k) - s0min1
                   s1max  = 1.5*slpmax
                   s1new  = min(s1max,max(-s1max,Tracer%sx(i,j,k)))
                   s2new  = min(slpmax+slpmax-0.3334*abs(s1new), max(abs(s1new)-slpmax,Tracer%sxx(i,j,k)))
                   Tracer%sxy(i,j,k) = min(slpmax,max(-slpmax,Tracer%sxy(i,j,k)))
                   Tracer%sxz(i,j,k) = min(slpmax,max(-slpmax,Tracer%sxz(i,j,k)))
                   Tracer%sx(i,j,k)  = s1new
                   Tracer%sxx(i,j,k) = s2new

                enddo
             enddo
          enddo

          call mpp_update_domains (Tracer%sx(:,:,:),  Dom_flux%domain2d, complete=.false.)
          call mpp_update_domains (Tracer%sxx(:,:,:), Dom_flux%domain2d, complete=.false.)
          call mpp_update_domains (Tracer%sxy(:,:,:), Dom_flux%domain2d, complete=.false.)
          call mpp_update_domains (Tracer%sxz(:,:,:), Dom_flux%domain2d, complete=.true.)

      endif

  endif

  do k=1,nk
    do j=jsc,jec
      do ip1=isc,iec+1
        i = ip1-1    
        velocity = Adv_vel%uhrho_et(i,j,k)*min(Grd%tmask(i,j,k),Grd%tmask(ip1,j,k))
        if (velocity >= 0.0) then

          ! flux from i to i+1 when u > 0
          ! create temporary moments/masses for partial boxes in transit

          fm(i,j,k)  = velocity*dtime*Grd%dyte(i,j)
          alf        = fm(i,j,k)/(sm(i,j,k)+epsln)
          alfq       = alf*alf
          alf1       = 1.0-alf
          alf1q      = alf1*alf1
          f0(i,j,k)  = alf*(Tracer%s0(i,j,k)+alf1*(Tracer%sx(i,j,k)+(alf1-alf)*Tracer%sxx(i,j,k)))
          fx(i,j,k)  = alfq*(Tracer%sx(i,j,k)+3.0*alf1*Tracer%sxx(i,j,k))
          fxx(i,j,k) = alf*alfq*Tracer%sxx(i,j,k)
          fy(i,j,k)  = alf*(Tracer%sy(i,j,k)+alf1*Tracer%sxy(i,j,k))
          fz(i,j,k)  = alf*(Tracer%sz(i,j,k)+alf1*Tracer%sxz(i,j,k))
          fxy(i,j,k) = alfq*Tracer%sxy(i,j,k)
          fxz(i,j,k) = alfq*Tracer%sxz(i,j,k)
          fyy(i,j,k) = alf*Tracer%syy(i,j,k)
          fzz(i,j,k) = alf*Tracer%szz(i,j,k)
          fyz(i,j,k) = alf*Tracer%syz(i,j,k)
          flux_x(i,j,k) = f0(i,j,k)*dtimer

          ! readjust moments remaining in the box

          sm(i,j,k)         = sm(i,j,k)-fm(i,j,k)
          Tracer%s0(i,j,k)  = Tracer%s0(i,j,k)-f0(i,j,k)
          Tracer%sx(i,j,k)  = alf1q*(Tracer%sx(i,j,k)-3.0*alf*Tracer%sxx(i,j,k))
          Tracer%sxx(i,j,k) = alf1*alf1q*Tracer%sxx(i,j,k)
          Tracer%sy(i,j,k)  = Tracer%sy(i,j,k)-fy(i,j,k)
          Tracer%syy(i,j,k) = Tracer%syy(i,j,k)-fyy(i,j,k)
          Tracer%sz(i,j,k)  = Tracer%sz(i,j,k)-fz(i,j,k)
          Tracer%szz(i,j,k) = Tracer%szz(i,j,k)-fzz(i,j,k)
          Tracer%sxy(i,j,k) = alf1q*Tracer%sxy(i,j,k)
          Tracer%sxz(i,j,k) = alf1q*Tracer%sxz(i,j,k)
          Tracer%syz(i,j,k) = Tracer%syz(i,j,k)-fyz(i,j,k)
        else

          ! flux from i+1 to i when u<0 (i.e., take left side of box ip1)

          fm(i,j,k)  = - velocity*dtime*Grd%dyte(i,j)
          alf        = fm(i,j,k)/(sm(ip1,j,k)+epsln)
          alfq       = alf*alf
          alf1       = 1.0-alf
          alf1q      = alf1*alf1
          f0(i,j,k)  = alf*(Tracer%s0(ip1,j,k)-alf1*(Tracer%sx(ip1,j,k)-(alf1-alf)*Tracer%sxx(ip1,j,k)))
          fx(i,j,k)  = alfq*(Tracer%sx(ip1,j,k)-3.0*alf1*Tracer%sxx(ip1,j,k))
          fxx(i,j,k) = alf*alfq*Tracer%sxx(ip1,j,k)
          fy(i,j,k)  = alf*(Tracer%sy(ip1,j,k)-alf1*Tracer%sxy(ip1,j,k))
          fz(i,j,k)  = alf*(Tracer%sz(ip1,j,k)-alf1*Tracer%sxz(ip1,j,k))
          fxy(i,j,k) = alfq*Tracer%sxy(ip1,j,k)
          fxz(i,j,k) = alfq*Tracer%sxz(ip1,j,k)
          fyy(i,j,k) = alf*Tracer%syy(ip1,j,k)
          fzz(i,j,k) = alf*Tracer%szz(ip1,j,k)
          fyz(i,j,k) = alf*Tracer%syz(ip1,j,k)
          flux_x(i,j,k) = - f0(i,j,k)*dtimer

          ! readjust moments remaining in the box

          sm(ip1,j,k)         = sm(ip1,j,k)-fm(i,j,k)
          Tracer%s0(ip1,j,k)  = Tracer%s0(ip1,j,k)-f0(i,j,k)
          Tracer%sx(ip1,j,k)  = alf1q*(Tracer%sx(ip1,j,k)+3.0*alf*Tracer%sxx(ip1,j,k))
          Tracer%sxx(ip1,j,k) = alf1*alf1q*Tracer%sxx(ip1,j,k)
          Tracer%sy(ip1,j,k)  = Tracer%sy(ip1,j,k)-fy(i,j,k)
          Tracer%syy(ip1,j,k) = Tracer%syy(ip1,j,k)-fyy(i,j,k)
          Tracer%sz(ip1,j,k)  = Tracer%sz(ip1,j,k)-fz(i,j,k)
          Tracer%szz(ip1,j,k) = Tracer%szz(ip1,j,k)-fzz(i,j,k)
          Tracer%sxy(ip1,j,k) = alf1q*Tracer%sxy(ip1,j,k)
          Tracer%sxz(ip1,j,k) = alf1q*Tracer%sxz(ip1,j,k)
          Tracer%syz(ip1,j,k) = Tracer%syz(ip1,j,k)-fyz(i,j,k)
        endif
      enddo
    enddo
  enddo

  call mpp_update_domains (f0(:,:,:),         Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (fx(:,:,:),         Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (fxx(:,:,:),        Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (fy(:,:,:),         Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (fz(:,:,:),         Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (fxy(:,:,:),        Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (fxz(:,:,:),        Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (fyy(:,:,:),        Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (fzz(:,:,:),        Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (fyz(:,:,:),        Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (flux_x(:,:,:),     Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (sm(:,:,:),         Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%s0(:,:,:),  Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%sx(:,:,:),  Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%sxx(:,:,:), Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%sy(:,:,:),  Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%syy(:,:,:), Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%sz(:,:,:),  Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%szz(:,:,:), Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%sxy(:,:,:), Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%sxz(:,:,:), Dom_flux%domain2d, flags=XUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%syz(:,:,:), Dom_flux%domain2d, flags=XUPDATE, complete=.true.)

  ! put the temporary moments fi into appropriate neighboring boxes
  do k = 1, nk
    do j=jsc,jec
      do ip1=isc,iec+1
        i = ip1-1    
        velocity = Adv_vel%uhrho_et(i,j,k)*min(Grd%tmask(i,j,k),Grd%tmask(ip1,j,k))
        if (velocity >= 0.0) then
          ! flux from i to i+1 when u > 0
          sm(ip1,j,k)         = sm(ip1,j,k)+fm(i,j,k)
          alf                 = fm(i,j,k)/(sm(ip1,j,k)+epsln)
          alf1                = 1.0-alf
          temptm              = alf*Tracer%s0(ip1,j,k)-alf1*f0(i,j,k)
          Tracer%s0(ip1,j,k)  = Tracer%s0(ip1,j,k)+f0(i,j,k)
          Tracer%sxx(ip1,j,k) = alf*alf*fxx(i,j,k)+alf1*alf1*Tracer%sxx(ip1,j,k)&
                                +5.0*(alf*alf1*(Tracer%sx(ip1,j,k)-fx(i,j,k))-(alf1-alf)*temptm)
          Tracer%sx(ip1,j,k)  = alf*fx(i,j,k)+alf1*Tracer%sx(ip1,j,k)+3.0*temptm
          Tracer%sxy(ip1,j,k) = alf*fxy(i,j,k)+alf1*Tracer%sxy(ip1,j,k) + 3.0*(-alf1*fy(i,j,k)+alf*Tracer%sy(ip1,j,k))
          Tracer%sxz(ip1,j,k) = alf*fxz(i,j,k)+alf1*Tracer%sxz(ip1,j,k) + 3.0*(-alf1*fz(i,j,k)+alf*Tracer%sz(ip1,j,k))
          Tracer%sy(ip1,j,k)  = Tracer%sy(ip1,j,k)+fy(i,j,k)
          Tracer%syy(ip1,j,k) = Tracer%syy(ip1,j,k)+fyy(i,j,k)
          Tracer%sz(ip1,j,k)  = Tracer%sz(ip1,j,k)+fz(i,j,k)
          Tracer%szz(ip1,j,k) = Tracer%szz(ip1,j,k)+fzz(i,j,k)
          Tracer%syz(ip1,j,k) = Tracer%syz(ip1,j,k)+fyz(i,j,k)
        else
          sm(i,j,k)         = sm(i,j,k)+fm(i,j,k)
          alf               = fm(i,j,k)/(sm(i,j,k)+epsln)
          alf1              = 1.0-alf
          temptm            = -alf*Tracer%s0(i,j,k)+alf1*f0(i,j,k)
          Tracer%s0(i,j,k)  = Tracer%s0(i,j,k)+f0(i,j,k)
          Tracer%sxx(i,j,k) = alf*alf*fxx(i,j,k)+alf1*alf1*Tracer%sxx(i,j,k)&
                              + 5.0*(alf*alf1*(-Tracer%sx(i,j,k)+fx(i,j,k))+(alf1-alf)*temptm)
          Tracer%sx(i,j,k)  = alf*fx(i,j,k)+alf1*Tracer%sx(i,j,k)+3.0*temptm
          Tracer%sxy(i,j,k) = alf*fxy(i,j,k)+alf1*Tracer%sxy(i,j,k) + 3.0*(alf1*fy(i,j,k)-alf*Tracer%sy(i,j,k))
          Tracer%sxz(i,j,k) = alf*fxz(i,j,k)+alf1*Tracer%sxz(i,j,k) + 3.0*(alf1*fz(i,j,k)-alf*Tracer%sz(i,j,k))
          Tracer%sy(i,j,k)  = Tracer%sy(i,j,k)+fy(i,j,k)
          Tracer%syy(i,j,k) = Tracer%syy(i,j,k)+fyy(i,j,k)
          Tracer%sz(i,j,k)  = Tracer%sz(i,j,k)+fz(i,j,k)
          Tracer%szz(i,j,k) = Tracer%szz(i,j,k)+fzz(i,j,k)
          Tracer%syz(i,j,k) = Tracer%syz(i,j,k)+fyz(i,j,k)
        endif
      enddo
    enddo
  enddo
  call mpp_update_domains (sm(:,:,:),         Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%s0(:,:,:),  Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sx(:,:,:),  Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sxx(:,:,:), Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sy(:,:,:),  Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%syy(:,:,:), Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sz(:,:,:),  Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%szz(:,:,:), Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sxy(:,:,:), Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sxz(:,:,:), Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%syz(:,:,:), Dom_flux%domain2d, complete=.true.)

  call mpp_clock_end(id_clock_psom_x)

end subroutine psom_x
! </SUBROUTINE> NAME="psom_x"


!#######################################################################
! <SUBROUTINE NAME="psom_y">
!
! <DESCRIPTION>
! Computes j-advective flux using Prather's SOM. 
! </DESCRIPTION>
!
subroutine psom_y (Adv_vel, Tracer, dtime, lag)

  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_prog_tracer_type), intent(inout) :: Tracer
  real,                         intent(in)    :: dtime
  integer,                      intent(in)    :: lag 

  real, dimension(isd:ied,jsd:jed,nk) :: f0, fm, fx, fxx, fy, fyy, fz, fzz, fxy, fxz, fyz
  integer :: i, j, k, jp1
  real    :: slpmax, s1max, s1new, s2new, velocity
  real    :: s0max1, s0min1
  real    :: alf, alfq, alf1, alf1q, temptm
  real    :: dtimer

  call mpp_clock_begin(id_clock_psom_y)
  dtimer = 1.0/dtime

  if (Tracer%psom_limit) then

      if(psom_limit_prather) then 

          ! limiter from Prather 
          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   slpmax = 0.0
                   if (Tracer%s0(i,j,k) > 0.0) slpmax = Tracer%s0(i,j,k)
                   s1max  = 1.5*slpmax
                   s1new  = min(s1max,max(-s1max,Tracer%sy(i,j,k)))
                   s2new  = min(slpmax+slpmax-0.3334*abs(s1new), max(abs(s1new)-slpmax,Tracer%syy(i,j,k)))
                   Tracer%sxy(i,j,k) = min(slpmax,max(-slpmax,Tracer%sxy(i,j,k)))
                   Tracer%syz(i,j,k) = min(slpmax,max(-slpmax,Tracer%syz(i,j,k)))
                   Tracer%sy(i,j,k)  = s1new
                   Tracer%syy(i,j,k) = s2new
                enddo
             enddo
          enddo

      else 

          ! limiter from Merryfield and Holloway, generalized to nonconstant grids
          do k=1,nk
             do j=jsc,jec
                do i=isd,ied

                   s0max1 = Tracer%field(i,j,k,lag)*Grd%tmask(i,j,k)
                   if(Grd%tmask(i,j-1,k) > 0.) s0max1 = max(Tracer%field(i,j-1,k,lag),s0max1)
                   if(Grd%tmask(i,j+1,k) > 0.) s0max1 = max(Tracer%field(i,j+1,k,lag),s0max1) 

                   s0max1 = s0max1*sm(i,j,k)
                   slpmax = -Tracer%s0(i,j,k) + s0max1
                   s1max  = 1.5*slpmax
                   s1new  = min(s1max,max(-s1max,-Tracer%sy(i,j,k)))
                   s2new  = min(slpmax+slpmax-0.3334*abs(s1new), max(abs(s1new)-slpmax,-Tracer%syy(i,j,k)))
                   Tracer%sxy(i,j,k) = -min(slpmax,max(-slpmax,-Tracer%sxy(i,j,k)))
                   Tracer%syz(i,j,k) = -min(slpmax,max(-slpmax,-Tracer%syz(i,j,k)))
                   Tracer%sy(i,j,k)  = -s1new
                   Tracer%syy(i,j,k) = -s2new

                   s0min1 = Tracer%field(i,j,k,lag)*Grd%tmask(i,j,k)
                   if(Grd%tmask(i,j-1,k) > 0.) s0min1 = min(Tracer%field(i,j-1,k,lag),s0min1)
                   if(Grd%tmask(i,j+1,k) > 0.) s0min1 = min(Tracer%field(i,j+1,k,lag),s0min1)  

                   s0min1 = s0min1*sm(i,j,k)
                   slpmax = Tracer%s0(i,j,k) - s0min1
                   s1max  = 1.5*slpmax
                   s1new  = min(s1max,max(-s1max,Tracer%sy(i,j,k)))
                   s2new  = min(slpmax+slpmax-0.3334*abs(s1new), max(abs(s1new)-slpmax,Tracer%syy(i,j,k)))
                   Tracer%sxy(i,j,k) = min(slpmax,max(-slpmax,Tracer%sxy(i,j,k)))
                   Tracer%syz(i,j,k) = min(slpmax,max(-slpmax,Tracer%syz(i,j,k)))
                   Tracer%sy(i,j,k)  = s1new
                   Tracer%syy(i,j,k) = s2new

                enddo
             enddo
          enddo

          call mpp_update_domains (Tracer%sy(:,:,:),  Dom_flux%domain2d, complete=.false.)
          call mpp_update_domains (Tracer%syy(:,:,:), Dom_flux%domain2d, complete=.false.)
          call mpp_update_domains (Tracer%sxy(:,:,:), Dom_flux%domain2d, complete=.false.)
          call mpp_update_domains (Tracer%syz(:,:,:), Dom_flux%domain2d, complete=.true.)

      endif

  endif

  do k=1,nk
    do jp1=jsc,jec+1
      j = jp1-1    
      do i=isc,iec
        velocity = Adv_vel%vhrho_nt(i,j,k)*min(Grd%tmask(i,j,k),Grd%tmask(i,jp1,k))
        if (velocity >= 0.0) then

          ! flux from j to j+1 when v > 0
          ! create temporary moments/masses for partial boxes in transit

          fm(i,j,k)  = velocity*dtime*Grd%dxtn(i,j)
          alf        = fm(i,j,k)/(sm(i,j,k)+epsln)
          alfq       = alf*alf
          alf1       = 1.0-alf
          alf1q      = alf1*alf1
          f0(i,j,k)  = alf*(Tracer%s0(i,j,k)+alf1*(Tracer%sy(i,j,k)+(alf1-alf)*Tracer%syy(i,j,k)))
          fy(i,j,k)  = alfq*(Tracer%sy(i,j,k)+3.0*alf1*Tracer%syy(i,j,k))
          fyy(i,j,k) = alf*alfq*Tracer%syy(i,j,k)
          fx(i,j,k)  = alf*(Tracer%sx(i,j,k)+alf1*Tracer%sxy(i,j,k))
          fz(i,j,k)  = alf*(Tracer%sz(i,j,k)+alf1*Tracer%syz(i,j,k))
          fxy(i,j,k) = alfq*Tracer%sxy(i,j,k)
          fyz(i,j,k) = alfq*Tracer%syz(i,j,k)
          fxx(i,j,k) = alf*Tracer%sxx(i,j,k)
          fzz(i,j,k) = alf*Tracer%szz(i,j,k)
          fxz(i,j,k) = alf*Tracer%sxz(i,j,k)
          flux_y(i,j,k) = f0(i,j,k)*dtimer

          ! readjust moments remaining in the box

          sm(i,j,k)         = sm(i,j,k)-fm(i,j,k)
          Tracer%s0(i,j,k)  = Tracer%s0(i,j,k)-f0(i,j,k)
          Tracer%sy(i,j,k)  = alf1q*(Tracer%sy(i,j,k)-3.0*alf*Tracer%syy(i,j,k))
          Tracer%syy(i,j,k) = alf1*alf1q*Tracer%syy(i,j,k)
          Tracer%sx(i,j,k)  = Tracer%sx(i,j,k)-fx(i,j,k)
          Tracer%sxx(i,j,k) = Tracer%sxx(i,j,k)-fxx(i,j,k)
          Tracer%sz(i,j,k)  = Tracer%sz(i,j,k)-fz(i,j,k)
          Tracer%szz(i,j,k) = Tracer%szz(i,j,k)-fzz(i,j,k)
          Tracer%sxy(i,j,k) = alf1q*Tracer%sxy(i,j,k)
          Tracer%syz(i,j,k) = alf1q*Tracer%syz(i,j,k)
          Tracer%sxz(i,j,k) = Tracer%sxz(i,j,k)-fxz(i,j,k)
        else

          ! flux from j+1 to j when v<0 (i.e., take south side of box jp1)

          fm(i,j,k)  = - velocity*dtime*Grd%dxtn(i,j)
          alf        = fm(i,j,k)/(sm(i,jp1,k)+epsln)
          alfq       = alf*alf
          alf1       = 1.0-alf
          alf1q      = alf1*alf1
          f0(i,j,k)  = alf*(Tracer%s0(i,jp1,k)-alf1*(Tracer%sy(i,jp1,k)-(alf1-alf)*Tracer%syy(i,jp1,k)))
          fy(i,j,k)  = alfq*(Tracer%sy(i,jp1,k)-3.0*alf1*Tracer%syy(i,jp1,k))
          fyy(i,j,k) = alf*alfq*Tracer%syy(i,jp1,k)
          fx(i,j,k)  = alf*(Tracer%sx(i,jp1,k)-alf1*Tracer%sxy(i,jp1,k))
          fz(i,j,k)  = alf*(Tracer%sz(i,jp1,k)-alf1*Tracer%syz(i,jp1,k))
          fxy(i,j,k) = alfq*Tracer%sxy(i,jp1,k)
          fyz(i,j,k) = alfq*Tracer%syz(i,jp1,k)
          fxx(i,j,k) = alf*Tracer%sxx(i,jp1,k)
          fzz(i,j,k) = alf*Tracer%szz(i,jp1,k)
          fxz(i,j,k) = alf*Tracer%sxz(i,jp1,k)
          flux_y(i,j,k) = - f0(i,j,k)*dtimer

          ! readjust moments remaining in the box

          sm(i,jp1,k)         = sm(i,jp1,k)-fm(i,j,k)
          Tracer%s0(i,jp1,k)  = Tracer%s0(i,jp1,k)-f0(i,j,k)
          Tracer%sy(i,jp1,k)  = alf1q*(Tracer%sy(i,jp1,k)+3.0*alf*Tracer%syy(i,jp1,k))
          Tracer%syy(i,jp1,k) = alf1*alf1q*Tracer%syy(i,jp1,k)
          Tracer%sx(i,jp1,k)  = Tracer%sx(i,jp1,k)-fx(i,j,k)
          Tracer%sxx(i,jp1,k) = Tracer%sxx(i,jp1,k)-fxx(i,j,k)
          Tracer%sz(i,jp1,k)  = Tracer%sz(i,jp1,k)-fz(i,j,k)
          Tracer%szz(i,jp1,k) = Tracer%szz(i,jp1,k)-fzz(i,j,k)
          Tracer%sxy(i,jp1,k) = alf1q*Tracer%sxy(i,jp1,k)
          Tracer%syz(i,jp1,k) = alf1q*Tracer%syz(i,jp1,k)
          Tracer%sxz(i,jp1,k) = Tracer%sxz(i,jp1,k)-fxz(i,j,k)
        endif
      enddo
    enddo
  enddo

  call mpp_update_domains (fm(:,:,:),         Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (f0(:,:,:),         Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (fx(:,:,:),         Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (fxx(:,:,:),        Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (fy(:,:,:),         Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (fz(:,:,:),         Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (fxy(:,:,:),        Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (fxz(:,:,:),        Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (fyy(:,:,:),        Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (fzz(:,:,:),        Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (fyz(:,:,:),        Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (flux_y(:,:,:),     Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (sm(:,:,:),         Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%s0(:,:,:),  Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%sx(:,:,:),  Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%sxx(:,:,:), Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%sy(:,:,:),  Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%syy(:,:,:), Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%sz(:,:,:),  Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%szz(:,:,:), Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%sxy(:,:,:), Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%sxz(:,:,:), Dom_flux%domain2d, flags=YUPDATE, complete=.false.)
  call mpp_update_domains (Tracer%syz(:,:,:), Dom_flux%domain2d, flags=YUPDATE, complete=.true.)

  ! put the temporary moments fj into appropriate neighboring boxes
  do k = 1, nk
    do jp1=jsc,jec+1
      j = jp1-1    
      do i=isc,iec
        velocity = Adv_vel%vhrho_nt(i,j,k)*min(Grd%tmask(i,j,k),Grd%tmask(i,jp1,k))
        if (velocity >= 0.0) then
          ! flux from j to j+1 when v > 0
          sm(i,jp1,k)         = sm(i,jp1,k)+fm(i,j,k)
          alf                 = fm(i,j,k)/(sm(i,jp1,k)+epsln)
          alf1                = 1.0-alf
          temptm              = alf*Tracer%s0(i,jp1,k)-alf1*f0(i,j,k)
          Tracer%s0(i,jp1,k)  = Tracer%s0(i,jp1,k)+f0(i,j,k)
          Tracer%syy(i,jp1,k) = alf*alf*fyy(i,j,k)+alf1*alf1*Tracer%syy(i,jp1,k)&
                                + 5.0*(alf*alf1*(Tracer%sy(i,jp1,k)-fy(i,j,k))-(alf1-alf)*temptm)
          Tracer%sy(i,jp1,k)  = alf*fy(i,j,k)+alf1*Tracer%sy(i,jp1,k)+3.0*temptm
          Tracer%sxy(i,jp1,k) = alf*fxy(i,j,k)+alf1*Tracer%sxy(i,jp1,k) + 3.0*(-alf1*fx(i,j,k)+alf*Tracer%sx(i,jp1,k)) 
          Tracer%syz(i,jp1,k) = alf*fyz(i,j,k)+alf1*Tracer%syz(i,jp1,k) + 3.0*(-alf1*fz(i,j,k)+alf*Tracer%sz(i,jp1,k))
          Tracer%sx(i,jp1,k)  = Tracer%sx(i,jp1,k)+fx(i,j,k)
          Tracer%sxx(i,jp1,k) = Tracer%sxx(i,jp1,k)+fxx(i,j,k)
          Tracer%sz(i,jp1,k)  = Tracer%sz(i,jp1,k)+fz(i,j,k)
          Tracer%szz(i,jp1,k) = Tracer%szz(i,jp1,k)+fzz(i,j,k)
          Tracer%sxz(i,jp1,k) = Tracer%sxz(i,jp1,k)+fxz(i,j,k)
        else
          sm(i,j,k)         = sm(i,j,k)+fm(i,j,k)
          alf               = fm(i,j,k)/(sm(i,j,k)+epsln)
          alf1              = 1.0-alf
          temptm            = -alf*Tracer%s0(i,j,k)+alf1*f0(i,j,k)
          Tracer%s0(i,j,k)  = Tracer%s0(i,j,k)+f0(i,j,k)
          Tracer%syy(i,j,k) = alf*alf*fyy(i,j,k)+alf1*alf1*Tracer%syy(i,j,k)&
                              + 5.0*(alf*alf1*(-Tracer%sy(i,j,k)+fy(i,j,k))+(alf1-alf)*temptm)
          Tracer%sy(i,j,k)  = alf*fy(i,j,k)+alf1*Tracer%sy(i,j,k)+3.0*temptm
          Tracer%sxy(i,j,k) = alf*fxy(i,j,k)+alf1*Tracer%sxy(i,j,k) + 3.0*(alf1*fx(i,j,k)-alf*Tracer%sx(i,j,k))
          Tracer%syz(i,j,k) = alf*fyz(i,j,k)+alf1*Tracer%syz(i,j,k) + 3.0*(alf1*fz(i,j,k)-alf*Tracer%sz(i,j,k))
          Tracer%sx(i,j,k)  = Tracer%sx(i,j,k)+fx(i,j,k)
          Tracer%sxx(i,j,k) = Tracer%sxx(i,j,k)+fxx(i,j,k)
          Tracer%sz(i,j,k)  = Tracer%sz(i,j,k)+fz(i,j,k)
          Tracer%szz(i,j,k) = Tracer%szz(i,j,k)+fzz(i,j,k)
          Tracer%sxz(i,j,k) = Tracer%sxz(i,j,k)+fxz(i,j,k)
        endif
      enddo
    enddo
  enddo
  call mpp_update_domains (sm(:,:,:),         Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%s0(:,:,:),  Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sx(:,:,:),  Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sxx(:,:,:), Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sy(:,:,:),  Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%syy(:,:,:), Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sz(:,:,:),  Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%szz(:,:,:), Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sxy(:,:,:), Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sxz(:,:,:), Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%syz(:,:,:), Dom_flux%domain2d, complete=.true.)

  call mpp_clock_end(id_clock_psom_y)

end subroutine psom_y
! </SUBROUTINE> NAME="psom_y"


!#######################################################################
! <SUBROUTINE NAME="psom_z">
!
! <DESCRIPTION>
! Computes k-advective flux using Prather's SOM. 
! </DESCRIPTION>
!
subroutine psom_z (Adv_vel, Tracer, dtime, lag)

  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_prog_tracer_type), intent(inout) :: Tracer
  real,                         intent(in)    :: dtime
  integer,                      intent(in)    :: lag 

  real, dimension(nk) :: f0, fm, fx, fxx, fy, fyy, fz, fzz, fxy, fxz, fyz
  integer :: i, j, k, km1, kp1
  real    :: slpmax, s1max, s1new, s2new, velocity
  real    :: s0max1, s0min1
  real    :: alf, alfq, alf1, alf1q, temptm
  real    :: dtimer

  call mpp_clock_begin(id_clock_psom_z)
  dtimer = 1.0/dtime

  ! when flux-limiting, limit appropriate moments before transport.
  if (Tracer%psom_limit) then

      if(psom_limit_prather) then 

          ! limiter from Prather
          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   slpmax = 0.0
                   if (Tracer%s0(i,j,k) > 0.0) slpmax = Tracer%s0(i,j,k)
                   s1max  = 1.5*slpmax
                   s1new  = min(s1max,max(-s1max,Tracer%sz(i,j,k)))
                   s2new  = min(slpmax+slpmax-0.3334*abs(s1new), max(abs(s1new)-slpmax,Tracer%szz(i,j,k)))
                   Tracer%sxz(i,j,k) = min(slpmax,max(-slpmax,Tracer%sxz(i,j,k)))
                   Tracer%syz(i,j,k) = min(slpmax,max(-slpmax,Tracer%syz(i,j,k)))
                   Tracer%sz(i,j,k)  = s1new
                   Tracer%szz(i,j,k) = s2new
                enddo
             enddo
          enddo

      else 


          ! limiter from Merryfield and Holloway, generalized to nonconstant grids
          do k=1,nk
             km1=max(1,k-1)
             kp1=min(nk,k+1)
             do j=jsd,jed
                do i=isd,ied

                   s0max1 = Tracer%field(i,j,k,lag)*Grd%tmask(i,j,k)
                   if(Grd%tmask(i,j,km1) > 0.) s0max1 = max(Tracer%field(i,j,km1,lag),s0max1)
                   if(Grd%tmask(i,j,kp1) > 0.) s0max1 = max(Tracer%field(i,j,kp1,lag),s0max1) 

                   s0max1 = s0max1*sm(i,j,k)
                   slpmax = -Tracer%s0(i,j,k) + s0max1
                   s1max  = 1.5*slpmax
                   s1new  = min(s1max,max(-s1max,-Tracer%sz(i,j,k)))
                   s2new  = min(slpmax+slpmax-0.3334*abs(s1new), max(abs(s1new)-slpmax,-Tracer%szz(i,j,k)))
                   Tracer%sxz(i,j,k) = -min(slpmax,max(-slpmax,-Tracer%sxz(i,j,k)))
                   Tracer%syz(i,j,k) = -min(slpmax,max(-slpmax,-Tracer%syz(i,j,k)))
                   Tracer%sz(i,j,k)  = -s1new
                   Tracer%szz(i,j,k) = -s2new

                   s0min1 = Tracer%field(i,j,k,lag)*Grd%tmask(i,j,k)
                   if(Grd%tmask(i,j,km1) > 0.) s0min1 = min(Tracer%field(i,j,km1,lag),s0min1)
                   if(Grd%tmask(i,j,kp1) > 0.) s0min1 = min(Tracer%field(i,j,kp1,lag),s0min1)  

                   s0min1 = s0min1*sm(i,j,k)
                   slpmax = Tracer%s0(i,j,k) - s0min1
                   s1max  = 1.5*slpmax
                   s1new  = min(s1max,max(-s1max,Tracer%sz(i,j,k)))
                   s2new  = min(slpmax+slpmax-0.3334*abs(s1new), max(abs(s1new)-slpmax,Tracer%szz(i,j,k)))
                   Tracer%sxz(i,j,k) = min(slpmax,max(-slpmax,Tracer%sxz(i,j,k)))
                   Tracer%syz(i,j,k) = min(slpmax,max(-slpmax,Tracer%syz(i,j,k)))
                   Tracer%sz(i,j,k)  = s1new
                   Tracer%szz(i,j,k) = s2new

                enddo
             enddo
          enddo

      endif

  endif


  do j=jsc,jec
    do i=isc,iec
      do km1=nk-1,1,-1
        k = km1 + 1
        velocity = Adv_vel%wrho_bt(i,j,km1)*Grd%tmask(i,j,k)
        if (velocity >= 0.0) then

          ! flux from k to k-1 when w > 0
          ! create temporary moments/masses for partial boxes in transit

          fm(k)  = velocity*dtime*Grd%dat(i,j)
          alf    = fm(k)/(sm(i,j,k)+epsln)
          alfq   = alf*alf
          alf1   = 1.0-alf
          alf1q  = alf1*alf1
          f0(k)  = alf*(Tracer%s0(i,j,k)+alf1*(Tracer%sz(i,j,k)+(alf1-alf)*Tracer%szz(i,j,k)))
          fz(k)  = alfq*(Tracer%sz(i,j,k)+3.0*alf1*Tracer%szz(i,j,k))
          fzz(k) = alf*alfq*Tracer%szz(i,j,k)
          fx(k)  = alf*(Tracer%sx(i,j,k)+alf1*Tracer%sxz(i,j,k))
          fy(k)  = alf*(Tracer%sy(i,j,k)+alf1*Tracer%syz(i,j,k))
          fxz(k) = alfq*Tracer%sxz(i,j,k)
          fyz(k) = alfq*Tracer%syz(i,j,k)
          fxx(k) = alf*Tracer%sxx(i,j,k)
          fyy(k) = alf*Tracer%syy(i,j,k)
          fxy(k) = alf*Tracer%sxy(i,j,k)
          flux_z(i,j,km1) = f0(k)*dtimer

          ! readjust moments remaining in the box

          sm(i,j,k)         = sm(i,j,k)-fm(k)
          Tracer%s0(i,j,k)  = Tracer%s0(i,j,k)-f0(k)
          Tracer%sz(i,j,k)  = alf1q*(Tracer%sz(i,j,k)-3.0*alf*Tracer%szz(i,j,k))
          Tracer%szz(i,j,k) = alf1*alf1q*Tracer%szz(i,j,k)
          Tracer%sx(i,j,k)  = Tracer%sx(i,j,k)-fx(k)
          Tracer%sxx(i,j,k) = Tracer%sxx(i,j,k)-fxx(k)
          Tracer%sy(i,j,k)  = Tracer%sy(i,j,k)-fy(k)
          Tracer%syy(i,j,k) = Tracer%syy(i,j,k)-fyy(k)
          Tracer%sxz(i,j,k) = alf1q*Tracer%sxz(i,j,k)
          Tracer%syz(i,j,k) = alf1q*Tracer%syz(i,j,k)
          Tracer%sxy(i,j,k) = Tracer%sxy(i,j,k)-fxy(k)
        else

          ! flux from k-1 to k when w<0 (i.e., take lower side of box k-1)

          fm(k)  = - velocity*dtime*Grd%dat(i,j)
          alf    = fm(k)/(sm(i,j,km1)+epsln)
          alfq   = alf*alf
          alf1   = 1.0-alf
          alf1q  = alf1*alf1
          f0(k)  = alf*(Tracer%s0(i,j,km1)-alf1*(Tracer%sz(i,j,km1)-(alf1-alf)*Tracer%szz(i,j,km1)))
          fz(k)  = alfq*(Tracer%sz(i,j,km1)-3.0*alf1*Tracer%szz(i,j,km1))
          fzz(k) = alf*alfq*Tracer%szz(i,j,km1)
          fx(k)  = alf*(Tracer%sx(i,j,km1)-alf1*Tracer%sxz(i,j,km1))
          fy(k)  = alf*(Tracer%sy(i,j,km1)-alf1*Tracer%syz(i,j,km1))
          fxz(k) = alfq*Tracer%sxz(i,j,km1)
          fyz(k) = alfq*Tracer%syz(i,j,km1)
          fxx(k) = alf*Tracer%sxx(i,j,km1)
          fyy(k) = alf*Tracer%syy(i,j,km1)
          fxy(k) = alf*Tracer%sxy(i,j,km1)
          flux_z(i,j,km1) = - f0(k)*dtimer

          ! readjust moments remaining in the box

          sm(i,j,km1)         = sm(i,j,km1)-fm(k)
          Tracer%s0(i,j,km1)  = Tracer%s0(i,j,km1)-f0(k)
          Tracer%sz(i,j,km1)  = alf1q*(Tracer%sz(i,j,km1)+3.0*alf*Tracer%szz(i,j,km1))
          Tracer%szz(i,j,km1) = alf1*alf1q*Tracer%szz(i,j,km1)
          Tracer%sx(i,j,km1)  = Tracer%sx(i,j,km1)-fx(k)
          Tracer%sxx(i,j,km1) = Tracer%sxx(i,j,km1)-fxx(k)
          Tracer%sy(i,j,km1)  = Tracer%sy(i,j,km1)-fy(k)
          Tracer%syy(i,j,km1) = Tracer%syy(i,j,km1)-fyy(k)
          Tracer%sxz(i,j,km1) = alf1q*Tracer%sxz(i,j,km1)
          Tracer%syz(i,j,km1) = alf1q*Tracer%syz(i,j,km1)
          Tracer%sxy(i,j,km1) = Tracer%sxy(i,j,km1)-fxy(k)
        endif
      enddo

      do km1=nk-1,1,-1
        k = km1 + 1
        velocity = Adv_vel%wrho_bt(i,j,km1)*Grd%tmask(i,j,k)
        if (velocity >= 0.0) then

          ! flux from k to k-1 when w > 0
          sm(i,j,km1)         = sm(i,j,km1)+fm(k)
          alf                 = fm(k)/(sm(i,j,km1)+epsln)
          alf1                = 1.0-alf
          temptm              = alf*Tracer%s0(i,j,km1)-alf1*f0(k)
          Tracer%s0(i,j,km1)  = Tracer%s0(i,j,km1)+f0(k)
          Tracer%szz(i,j,km1) = alf*alf*fzz(k)+alf1*alf1*Tracer%szz(i,j,km1)&
                                + 5.0*(alf*alf1*(Tracer%sz(i,j,km1)-fz(k))-(alf1-alf)*temptm)
          Tracer%sz(i,j,km1)  = alf*fz(k)+alf1*Tracer%sz(i,j,km1)+3.0*temptm
          Tracer%sxz(i,j,km1) = alf*fxz(k)+alf1*Tracer%sxz(i,j,km1) + 3.0*(-alf1*fx(k)+alf*Tracer%sx(i,j,km1))
          Tracer%syz(i,j,km1) = alf*fyz(k)+alf1*Tracer%syz(i,j,km1) + 3.0*(-alf1*fy(k)+alf*Tracer%sy(i,j,km1))
          Tracer%sx(i,j,km1)  = Tracer%sx(i,j,km1)+fx(k)
          Tracer%sxx(i,j,km1) = Tracer%sxx(i,j,km1)+fxx(k)
          Tracer%sy(i,j,km1)  = Tracer%sy(i,j,km1)+fy(k)
          Tracer%syy(i,j,km1) = Tracer%syy(i,j,km1)+fyy(k)
          Tracer%sxy(i,j,km1) = Tracer%sxy(i,j,km1)+fxy(k)
        else
          sm(i,j,k)         = sm(i,j,k)+fm(k)
          alf               = fm(k)/(sm(i,j,k)+epsln)
          alf1              = 1.0-alf
          temptm            = -alf*Tracer%s0(i,j,k)+alf1*f0(k)
          Tracer%s0(i,j,k)  = Tracer%s0(i,j,k)+f0(k)
          Tracer%szz(i,j,k) = alf*alf*fzz(k)+alf1*alf1*Tracer%szz(i,j,k)&
                              + 5.0*(alf*alf1*(-Tracer%sz(i,j,k)+fz(k))+(alf1-alf)*temptm)
          Tracer%sz(i,j,k)  = alf*fz(k)+alf1*Tracer%sz(i,j,k)+3.0*temptm
          Tracer%sxz(i,j,k) = alf*fxz(k)+alf1*Tracer%sxz(i,j,k) + 3.0*(alf1*fx(k)-alf*Tracer%sx(i,j,k))
          Tracer%syz(i,j,k) = alf*fyz(k)+alf1*Tracer%syz(i,j,k) + 3.0*(alf1*fy(k)-alf*Tracer%sy(i,j,k))
          Tracer%sx(i,j,k)  = Tracer%sx(i,j,k)+fx(k)
          Tracer%sxx(i,j,k) = Tracer%sxx(i,j,k)+fxx(k)
          Tracer%sy(i,j,k)  = Tracer%sy(i,j,k)+fy(k)
          Tracer%syy(i,j,k) = Tracer%syy(i,j,k)+fyy(k)
          Tracer%sxy(i,j,k) = Tracer%sxy(i,j,k)+fxy(k)
        endif
      enddo
    enddo
  enddo
  call mpp_update_domains (sm(:,:,:),         Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%s0(:,:,:),  Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sx(:,:,:),  Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sxx(:,:,:), Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sy(:,:,:),  Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%syy(:,:,:), Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sz(:,:,:),  Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%szz(:,:,:), Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sxy(:,:,:), Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%sxz(:,:,:), Dom_flux%domain2d, complete=.false.)
  call mpp_update_domains (Tracer%syz(:,:,:), Dom_flux%domain2d, complete=.true.)

  call mpp_clock_end(id_clock_psom_z)

end subroutine psom_z
! </SUBROUTINE> NAME="psom_z"




!#######################################################################
! <FUNCTION NAME="advect_tracer_mdppm_test">
!
! <DESCRIPTION>
!
! Compute advective tendency of dzt*rho*tracer concentration using 
! multi-dimensional piecewise parabolic method. 
! 
! Controlling parameters:
!   Tracer%ppm_hlimiter=0  -> No limiting in the horizontal
!          ppm_hlimiter=1  -> Full constraint (Colella and Woodward, 1984)
!          ppm_hlimiter=2  -> Improved full constraint (Lin, 2004)
!          ppm_hlimiter=3  -> Huynh limiter (Huynh, 1996)
!
!          ppm_vlimiter=*  -> as for ppm_hlimit but for the vertical
!
! Coded by Alistair.Adcroft
! Jan/Feb 2006
!
! Updated with controlling parameters, April 2006
!
! This is an unsupported test version of the MDPPM scheme. It is not
! meant for general use. 
!
! </DESCRIPTION>
!
function advect_tracer_mdppm_test(Time, Adv_vel, Tracer, Thickness, Tracer_field, dtime)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field
  real,                         intent(in) :: dtime

  real,dimension(isc:iec,jsc:jec,nk)       :: advect_tracer_mdppm_test

  real,dimension(isc:iec,jsc:jec)          :: ftp
  real,dimension(isc:iec,jsc:jec)          :: fbt
  real,dimension(isc:iec,jsc:jec)          :: wkm1

  integer                                  :: i, j, k
  integer                                  :: kp1, kp2, km1, km2
  integer                                  :: tau, taum1

  real                                     :: cfl, massflux, massflux_bt, dMx, dMn
  real                                     :: Sim2,Sim1,Si,Sip1,Sip2

  real                                     :: da2,da4,da3m,da3p
  real                                     :: mskm2,mskm1,msk0,mskp1,mskp2
  real,parameter                           :: oneSixth=1./6., r12=1./12.
  real,parameter                           :: twoThirds=2./3.
  real,dimension(isc-4:iec+4,jsc-4:jec+4)  :: da, aL, aR, a6, d1m, d1p, d1mm, d1pp
  real,dimension(isc-4:iec+4,jsc-4:jec+4,nk) :: dak

  call mpp_clock_begin(id_clock_mdppm_test)

  ftp          = 0.0
  fbt          = 0.0
  wkm1         = 0.0 !aja DELETE?
  tracer_mdppm = 0.0
  tracermass_mdppm = 0.0
  mass_mdppm   = 0.0
  flux_x       = 0.0
  flux_y       = 0.0
  flux_z       = 0.0

  ! initialize arrays with extended halos 
  da           = 0.0
  aL           = 0.0
  aR           = 0.0
  a6           = 0.0
  d1m          = 0.0
  d1p          = 0.0
  d1mm         = 0.0
  d1pp         = 0.0
  dak          = 0.0

  tau   = Time%tau
  taum1 = Time%taum1

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tracer_mdppm(i,j,k) = Tracer_field(i,j,k) * tmask_mdppm(i,j,k)
           mass_mdppm(i,j,k) = Thickness%rho_dzt(i,j,k,tau) * Grd%dat(i,j)
           tracermass_mdppm(i,j,k) = mass_mdppm(i,j,k) * tracer_mdppm(i,j,k)
        enddo
     enddo
  enddo
!
! Calculate vertical flux at the top and bottom of the boxes
!

  ! This block calculates the slope in each cell (aka PLM but
  ! with a fourth-order estimate)
  do k=1,nk
     km2 = max(k-2,1)
     km1 = max(k-1,1)
     kp1 = min(k+1,nk)
     kp2 = min(k+2,nk)
     do j=jsc,jec
        do i=isc,iec

        Sip2 = tracer_mdppm(i,j,kp2)
        Sip1 = tracer_mdppm(i,j,kp1)
        Si   = tracer_mdppm(i,j,k)
        Sim1 = tracer_mdppm(i,j,km1)
        Sim2 = tracer_mdppm(i,j,km2)
        ! Simple unmasked 4th order
        da4 = r12 * ( 8. * ( Sip1 - Sim1 ) + ( Sim2 - Sip2 ) )
        ! Simple 2nd order
        da2 = 0.5 * ( Sip1 - Sim1 )
        ! Simple 3rd order, biased to left
        da3m = r12 * ( 2. * Sim2 - 12. * Sim1 + 6. * Si + 4. *Sip1 )
        ! Simple 3rd order, biased to right
        da3p = r12 * ( - 4. * Sim1 - 6. * Si + 12. *Sip1  - 2. * Sip2 )
        ! Estimate of slope: dak(i,j,k) is the twice the "mismatch" as defined
        ! by Lin which itself is only half of the slope. The factor of two
        ! in the mismatch is confusing so we use simple estimate of slope.
        ! Note: The temperature masks used here are proxies for the advective
        ! flow masks that would be used if they existed. The only functional
        ! difference is that this scheme will not work with thin walls. 
        mskp2 = tmask_mdppm(i,j,kp2) * real(kp2-kp1)
        mskp1 = tmask_mdppm(i,j,kp1) * real(kp1-k)
        msk0  = tmask_mdppm(i,j,k)
        mskm1 = tmask_mdppm(i,j,km1) * real(k-km1)
        mskm2 = tmask_mdppm(i,j,km2) * real(km1-km2)
        dak(i,j,k) = ( (   mskm2 * ( 1. - 0.5 * mskp2 ) * da3m       &
                         + mskp2 * ( 1. - 0.5 * mskm2 ) * da3p )     &
                       + ( 1. - mskm2 ) * ( 1. - mskp2 ) * da2 )
        !!! dak(i,j,k) = da4 ! *** UNMASKED ***
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max( Sip1, Sim1,  Si ) - Si
        dMn = Si - min( Sip1, Sim1, Si )
        dak(i,j,k) = sign(1.,dak(i,j,k)) &
                     * min( abs(dak(i,j,k)) , 2. * min( dMx , dMn ) )  &
                     * mskm1 * mskp1 * msk0
        enddo
     enddo
  enddo

  ! This block estimates the upper (L) and lower (R) values on cell boundaries
  ! in k-direction
  do k=1,nk
     km2 = max(k-2,1)
     km1 = max(k-1,1)
     kp1 = min(k+1,nk)
     kp2 = min(k+2,nk)
     do j=jsc,jec
        do i=isc,iec

        mskp2 = tmask_mdppm(i,j,kp2) * real(kp2-kp1)
        mskp1 = tmask_mdppm(i,j,kp1) * real(kp1-k)
        msk0  = tmask_mdppm(i,j,k)
        mskm1 = tmask_mdppm(i,j,km1) * real(k-km1)
        mskm2 = tmask_mdppm(i,j,km2) * real(km1-km2)
        Sip2 = tracer_mdppm(i,j,kp2)
        Sip1 = tracer_mdppm(i,j,kp1)
        Si   = tracer_mdppm(i,j,k)
        Sim1 = tracer_mdppm(i,j,km1)
        Sim2 = tracer_mdppm(i,j,km2)

        ! Calculate the second difference for use in the Huynh limiter later
        d1m(i,j)  = ( Si - Sim1 ) * mskm1
        d1p(i,j)  = ( Sip1 - Si ) * mskp1
        d1mm(i,j) = ( Sim1 - Sim2 ) * mskm1 * mskm2
        d1pp(i,j) = ( Sip2 - Sip1 ) * mskp1 * mskp2

        Sim1 = Si - d1m(i,j)
        Sip1 = Si + d1p(i,j)
        ! Left edge: Eq. B2 in Lin 1994, MWR (132)
        aL(i,j) = 0.5 * ( Sim1 + Si ) + oneSixth * ( dak(i,j,km1) - dak(i,j,k) )
        ! Right edge: Eq. B2 in Lin 1994, MWR (132)
        aR(i,j) = 0.5 * ( Si + Sip1 ) + oneSixth * ( dak(i,j,k) - dak(i,j,kp1) )        
        enddo
     enddo

     ! Limit the edge values, aL and aR, for monotonicity
     if (Tracer%ppm_hlimiter.eq.1) then
       call ppm_limit_cw84(isc,iec,jsc,jec, tracer_mdppm(isc-4,jsc-4,k), aL, aR)
     elseif (Tracer%ppm_hlimiter.eq.2) then
       call ppm_limit_ifc(isc,iec,jsc,jec, tracer_mdppm(isc-4,jsc-4,k), dak(isc-4,jsc-4,k), aL,aR)
     elseif (Tracer%ppm_hlimiter.eq.3) then
       call ppm_limit_sh(isc,iec,jsc,jec, &
                        tracer_mdppm(isc-4,jsc-4,k), &
                        d1m, d1p, d1mm, d1pp, aL, aR)
     else
       call mpp_error(FATAL,&
       '==>Error from ocean_tracer_advect_mod: must choose ppm_hlimiter=1,2,or 3 in field table.')
     endif

     ! Calculate flux at top or bottom of box depending on sign of flow
     do j=jsc,jec
        do i=isc,iec

           mskp1 = tmask_mdppm(i,j,kp1) * real(kp1-k)
           msk0  = tmask_mdppm(i,j,k)
           mskm1 = tmask_mdppm(i,j,km1) * real(k-km1)
           Si   = tracer_mdppm(i,j,k)

           ! Curvature in k-direction
           a6(i,j) = 6. * Si - 3. * ( aR(i,j) + aL(i,j) )


           massflux = Grd%dat(i,j) * Adv_vel%wrho_bt(i,j,k) * msk0 * mskp1
           if (massflux*Thickness%rho_dzt(i,j,k,tau)<0.0) then
             cfl = abs(Adv_vel%wrho_bt(i,j,k) * dtime / Thickness%rho_dzt(i,j,k,tau))
             flux_z(i,j,k) = massflux                                        &
                           * ( aR(i,j) + 0.5 * cfl * ( ( aL(i,j) - aR(i,j) ) &
                               + ( 1. - twoThirds * cfl ) * a6(i,j) ) )
           else
             flux_z(i,j,k) = 0.0
           endif

           massflux_bt = massflux

           massflux = Grd%dat(i,j) * Adv_vel%wrho_bt(i,j,km1) * msk0 * mskm1
           if (massflux*Thickness%rho_dzt(i,j,k,tau)>0.0) then
               cfl = abs(Adv_vel%wrho_bt(i,j,km1) * dtime / Thickness%rho_dzt(i,j,k,tau))
               flux_z(i,j,km1) = massflux                                      &
                    * ( aL(i,j) + 0.5 * cfl * ( ( aR(i,j) - aL(i,j) ) &
                    + ( 1. - twoThirds * cfl ) * a6(i,j) ) ) 
           endif

           mass_mdppm(i,j,k) = mass_mdppm(i,j,k) + dtime * ( massflux_bt - massflux  )

        enddo
     enddo
  enddo

  do k=1,nk
     km1 = max(k-1,1)
     mskm1 = real(k-km1)
     do j=jsc,jec
        do i=isc,iec
           tracermass_mdppm(i,j,k) = tracermass_mdppm(i,j,k) + dtime * ( flux_z(i,j,k) - mskm1*flux_z(i,j,km1) )

           if (mass_mdppm(i,j,k)>0.) tracer_mdppm(i,j,k) = tracermass_mdppm(i,j,k) / mass_mdppm(i,j,k)

!       tracer_mdppm(i,j,k) = tracer_mdppm(i,j,k)                          &
!            + dtime / Thickness%rho_dzt(i,j,k,tau) * (                    &
!               Grd%datr(i,j) * ( flux_z(i,j,k) - mskm1 * flux_z(i,j,km1)) &
!            + Tracer_field(i,j,k) *                                       &
!              (mskm1 * Adv_vel%wrho_bt(i,j,km1) - Adv_vel%wrho_bt(i,j,k)) )
        enddo
     enddo

     ! update vertical velocity for next k-level
     do j=jsc,jec
        do i=isc,iec
           wkm1(i,j) = Adv_vel%wrho_bt(i,j,k)   !aja DELETE?
        enddo
     enddo

  enddo ! end of k-loop
  call mpp_update_domains (tracer_mdppm, Dom_mdppm%domain2d, flags=XUPDATE)
  call mpp_update_domains (tracermass_mdppm, Dom_mdppm%domain2d, flags=XUPDATE)
  call mpp_update_domains (mass_mdppm, Dom_mdppm%domain2d, flags=XUPDATE)

!
! Calculate flux at the eastern wall of the boxes
!
  do k=1,nk
    do j=jsc,jec
      do i=isc-2,iec+2
        ! This block calculates the slope in each cell (aka PLM but
        ! with a fourth-order estimate)
        Sip2 = tracer_mdppm(i+2,j,k)
        Sip1 = tracer_mdppm(i+1,j,k)
        Si   = tracer_mdppm(i,j,k)
        Sim1 = tracer_mdppm(i-1,j,k)
        Sim2 = tracer_mdppm(i-2,j,k)
        ! Simple unmasked 4th order
        da4 = r12 * ( 8. * ( Sip1 - Sim1 ) + ( Sim2 - Sip2 ) )
        ! Simple 2nd order
        da2 = 0.5 * ( Sip1 - Sim1 )
        ! Simple 3rd order, biased to left
        da3m = r12 * ( 2. * Sim2 - 12. * Sim1 + 6. * Si + 4. *Sip1 )
        ! Simple 3rd order, biased to right
        da3p = r12 * ( - 4. * Sim1 - 6. * Si + 12. *Sip1  - 2. * Sip2 )
        ! Estimate of slope: da(i,j) is the twice the "mismatch" as defined
        ! by Lin which itself is only half of the slope. The factor of two
        ! in the mismatch is confusing so we use simple estimate of slope.
        ! Note: The temperature masks used here are proxies for the advective
        ! flow masks that would be used if they existed. The only functional
        ! difference is that this scheme will not work with thin walls. 
        mskm2 = tmask_mdppm(i-2,j,k)
        mskm1 = tmask_mdppm(i-1,j,k)
        msk0  = tmask_mdppm(i,j,k)
        mskp1 = tmask_mdppm(i+1,j,k)
        mskp2 = tmask_mdppm(i+2,j,k)
        da(i,j) = ( (   mskm2 * ( 1. - 0.5 * mskp2 ) * da3m       &
                      + mskp2 * ( 1. - 0.5 * mskm2 ) * da3p )     &
                    + ( 1. - mskm2 ) * ( 1. - mskp2 ) * da2 )
        !!! da(i,j) = da4 ! *** UNMASKED ***
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max( Sip1, Sim1,  Si ) - Si
        dMn = Si - min( Sip1, Sim1, Si )
        da(i,j) = sign(1.,da(i,j)) * min( abs(da(i,j)) , 2. * min( dMx , dMn ) )  &
                  * mskm1 * mskp1 * msk0
        ! Calculate the second difference for use in the Huynh limiter later
        d1m(i,j)  = ( Si - Sim1 ) * mskm1
        d1p(i,j)  = ( Sip1 - Si ) * mskp1
        d1mm(i,j) = ( Sim1 - Sim2 ) * mskm1 * mskm2
        d1pp(i,j) = ( Sip2 - Sip1 ) * mskp1 * mskp2
      enddo !i

      do i=isc-1,iec+1
        ! This block estimates the left and right values on cell boundaries
        ! in i-direction
        Si   = tracer_mdppm(i,j,k)
        Sim1 = Si - d1m(i,j)
        Sip1 = Si + d1p(i,j)
!       mskm1 = tmask_mdppm(i-1,j,k)
!       mskp1 = tmask_mdppm(i+1,j,k)
        ! Left edge: Eq. B2 in Lin 1994, MWR (132)
        aL(i,j) = 0.5 * ( Sim1 + Si ) + oneSixth * ( da(i-1,j) - da(i,j) )        
        ! ... and masks
!       aL(i,j) = mskm1 * aL(i,j) + ( 1. - mskm1 ) * Si
        ! Right edge: Eq. B2 in Lin 1994, MWR (132)
        aR(i,j) = 0.5 * ( Si + Sip1 ) + oneSixth * ( da(i,j) - da(i+1,j) )        
        ! ... and masks
!       aR(i,j) = mskp1 * aR(i,j) + ( 1. - mskp1 ) * Si
      enddo !i
    enddo !j

    ! Limit the edge values, aL and aR, for monotonicity
    if (Tracer%ppm_hlimiter.eq.1) then
      call ppm_limit_cw84(isc,iec,jsc,jec, tracer_mdppm(isc-4,jsc-4,k), aL, aR)
    elseif (Tracer%ppm_hlimiter.eq.2) then
      call ppm_limit_ifc(isc,iec,jsc,jec, tracer_mdppm(isc-4,jsc-4,k), da, aL, aR)
    elseif (Tracer%ppm_hlimiter.eq.3) then
      call ppm_limit_sh(isc,iec,jsc,jec, &
                       tracer_mdppm(isc-4,jsc-4,k), &
                       d1m, d1p, d1mm, d1pp, aL, aR)
    else
      call mpp_error(FATAL,&
      '==>Error from ocean_tracer_advect_mod: must choose %ppm_hlimiter=1,2, or 3 in field table.')
    endif

    do j=jsc,jec
! NOTE TO SELF(AJA): SHOULD BE ABLE TO ELIMINATE A6
      do i=isc-1,iec+1
        ! Curvature in i-direction                            
        Si   = tracer_mdppm(i,j,k)
        a6(i,j) = 6. * Si - 3. * ( aR(i,j) + aL(i,j) )        
      enddo !i
                                                              
      do i=isc-1,iec                                          
                                                              
        ! Volume flux and advective flux at right edge of cell (i.e. i+1/2) 
        massflux = Grd%dyte(i,j) * Adv_vel%uhrho_et(i,j,k)        

        if (massflux>0.0) then                              
         cfl = abs(massflux) * dtime / mass_mdppm(i,j,k) 
         flux_x(i,j,k) = massflux * tmask_mdppm(i,j,k) * tmask_mdppm(i+1,j,k) &
                   * ( aR(i,j) + 0.5 * cfl * ( ( aL(i,j) - aR(i,j) )          &
                       + ( 1. - twoThirds * cfl ) * a6(i,j) ) )   
        elseif (massflux<0.0) then
         cfl = abs(massflux) * dtime / mass_mdppm(i+1,j,k) 
         flux_x(i,j,k) = massflux * tmask_mdppm(i,j,k) * tmask_mdppm(i+1,j,k) &
                   * ( aL(i+1,j) - 0.5 * cfl * ( ( aR(i+1,j) - aL(i+1,j) )    &
                       + ( 1. + twoThirds * cfl ) * a6(i+1,j) ) ) 
        else
         flux_x(i,j,k) = 0.0
        endif                                                 
                                                              
      enddo !i                                                
    enddo !j                                                  
!                                                             
! Update the tracer                                           
!                                                             
    do j=jsc,jec                                              
      do i=isc,iec                                            

        mass_mdppm(i,j,k) = mass_mdppm(i,j,k) + dtime * (   &
                Grd%dyte(i-1,j) * Adv_vel%uhrho_et(i-1,j,k) &
              - Grd%dyte(i,j) * Adv_vel%uhrho_et(i,j,k) )

        tracermass_mdppm(i,j,k) = tracermass_mdppm(i,j,k) + dtime * ( flux_x(i-1,j,k) - flux_x(i,j,k) )

        if (mass_mdppm(i,j,k)>0.) tracer_mdppm(i,j,k) = tracermass_mdppm(i,j,k) / mass_mdppm(i,j,k)

      enddo !i                                                
    enddo !j                                                  
  enddo !k                                                    
                                                              
! Update the 4-point tracer halo                              
  call mpp_update_domains (tracer_mdppm, Dom_mdppm%domain2d, flags=YUPDATE) 
  call mpp_update_domains (tracermass_mdppm, Dom_mdppm%domain2d, flags=YUPDATE) 
  call mpp_update_domains (mass_mdppm, Dom_mdppm%domain2d, flags=YUPDATE) 
!ajacall stats(mass_mdppm(isc:iec,jsc:jec,:),'M(u)')
!ajacall stats(tracer_mdppm(isc:iec,jsc:jec,:),'T(u)')
!ajacall stats(tracermass_mdppm(isc:iec,jsc:jec,:),'TM(u)')

!aja! No flux at bottom (for final divergence of fluxes below)
!aja do j=jsc,jec
!aja   do i=isc,iec
!aja     wkm1(i,j) = 0.0
!aja   enddo !i
!aja enddo !j

!
! Calculate flux at the northern wall of the boxes
!
  do k=1,nk
    do j=jsc-2,jec+2
      do i=isc,iec
        ! This block calculates the slope in each cell (aka PLM but
        ! with a fourth-order estimate)
        Sip2 = tracer_mdppm(i,j+2,k)
        Sip1 = tracer_mdppm(i,j+1,k)
        Si   = tracer_mdppm(i,j,k)                            
        Sim1 = tracer_mdppm(i,j-1,k)
        Sim2 = tracer_mdppm(i,j-2,k)
        ! Simple unmasked 4th order
        da4 = r12 * ( 8. * ( Sip1 - Sim1 ) + ( Sim2 - Sip2 ) )
        ! Simple 2nd order
        da2 = 0.5 * ( Sip1 - Sim1 )
        ! Simple 3rd order, biased to left
        da3m = r12 * ( 2. * Sim2 - 12. * Sim1 + 6. * Si + 4. *Sip1 )
        ! Simple 3rd order, biased to right
        da3p = r12 * ( - 4. * Sim1 - 6. * Si + 12. *Sip1  - 2. * Sip2 )
        ! Estimate of slope: da(i,j) is the twice the "mismatch" as defined
        ! by Lin which itself is only half of the slope. The factor of two
        ! in the mismatch is confusing so we use simple estimate of slope.
        ! Note: The temperature masks used here are proxies for the advective
        ! flow masks that would be used if they existed. The only functional
        ! difference is that this scheme will not work with thin walls. 
        mskm2 = tmask_mdppm(i,j-2,k)
        mskm1 = tmask_mdppm(i,j-1,k)
        msk0  = tmask_mdppm(i,j,k)
        mskp1 = tmask_mdppm(i,j+1,k)
        mskp2 = tmask_mdppm(i,j+2,k)
        da(i,j) = ( (   mskm2 * ( 1. - 0.5 * mskp2 ) * da3m       &
                      + mskp2 * ( 1. - 0.5 * mskm2 ) * da3p )     &
                    + ( 1. - mskm2 ) * ( 1. - mskp2 ) * da2 )
        !!! da(i,j) = da4 ! *** UNMASKED ***
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max( Sip1, Sim1, Si ) - Si             
        dMn = Si - min( Sip1, Sim1, Si )             
        da(i,j) = sign(1.,da(i,j)) * min( abs(da(i,j)) , 2. * min( dMx , dMn ) )  &
                  * mskm1 * mskp1 * msk0
        ! Calculate the second difference for use in the Huynh limiter later
        d1m(i,j)  = ( Si - Sim1 ) * mskm1
        d1p(i,j)  = ( Sip1 - Si ) * mskp1
        d1mm(i,j) = ( Sim1 - Sim2 ) * mskm1 * mskm2
        d1pp(i,j) = ( Sip2 - Sip1 ) * mskp1 * mskp2
      enddo !i                                                
    enddo !j                                                  
                                                              
    do j=jsc-1,jec+1                                          
      do i=isc,iec                                            
        ! This block estimates the left and right values on cell boundaries
        ! in j-direction
        Si   = tracer_mdppm(i,j,k)                            
        Sim1 = Si - d1m(i,j)
        Sip1 = Si + d1p(i,j)
!       mskm1 = tmask_mdppm(i-1,j,k)
!       mskp1 = tmask_mdppm(i+1,j,k)
        ! Left edge: Eq. B2 in Lin 1994, MWR (132)
        aL(i,j) = 0.5 * ( Sim1 + Si ) + oneSixth * ( da(i,j-1) - da(i,j) )        
        ! ... and masks
!       aL(i,j) = mskm1 * aL(i,j) + ( 1. - mskm1 ) * Si
        ! Right edge: Eq. B2 in Lin 1994, MWR (132)
        aR(i,j) = 0.5 * ( Si + Sip1 ) + oneSixth * ( da(i,j) - da(i,j+1) )        
        ! ... and masks
!       aR(i,j) = mskp1 * aR(i,j) + ( 1. - mskp1 ) * Si
      enddo !i                                                
    enddo !j                                                  

    ! Limit the edge values, aL and aR, for monotonicity
    if (Tracer%ppm_hlimiter.eq.1) then
      call ppm_limit_cw84(isc,iec,jsc,jec, tracer_mdppm(isc-4,jsc-4,k), aL, aR)
    elseif (Tracer%ppm_hlimiter.eq.2) then
      call ppm_limit_ifc(isc,iec,jsc,jec, tracer_mdppm(isc-4,jsc-4,k), da, aL, aR)
    elseif (Tracer%ppm_hlimiter.eq.3) then
      call ppm_limit_sh(isc,iec,jsc,jec, &
                       tracer_mdppm(isc-4,jsc-4,k), &
                       d1m, d1p, d1mm, d1pp, aL, aR)
    else
      call mpp_error(FATAL,&
      '==>Error from ocean_tracer_advect_mod: must choose ppm_hlimiter=1,2, or 3 in field table.')
    endif

! NOTE TO SELF(AJA): SHOULD BE ABLE TO ELIMINATE A6
    do j=jsc-1,jec+1                                          
      do i=isc,iec                                            
        ! Curvature in j-direction                            
        Si   = tracer_mdppm(i,j,k)
        a6(i,j) = 6. * Si - 3. * ( aR(i,j) + aL(i,j) )        
      enddo !i                                                
    enddo !j                                                  
                                                              
    do j=jsc-1,jec                                            
      do i=isc,iec                                            
                                                              
        ! Volume flux and advective flux at right edge of cell (i.e. i+1/2) 
        massflux = Grd%dxtn(i,j) * Adv_vel%vhrho_nt(i,j,k)        

        if (massflux>0.0) then                              
         cfl = abs(massflux) * dtime / mass_mdppm(i,j,k) 
         flux_y(i,j,k) = massflux * tmask_mdppm(i,j,k) * tmask_mdppm(i,j+1,k) &
                   * ( aR(i,j) + 0.5 * cfl * ( ( aL(i,j) - aR(i,j) )          &
                       + ( 1. - 2./3. * cfl ) * a6(i,j) ) )   
        elseif (massflux<0.0) then                              
         cfl = abs(massflux) * dtime / mass_mdppm(i,j+1,k) 
         flux_y(i,j,k) = massflux * tmask_mdppm(i,j,k) * tmask_mdppm(i,j+1,k) &
                   * ( aL(i,j+1) - 0.5 * cfl * ( ( aR(i,j+1) - aL(i,j+1) )    &
                       + ( 1. + 2./3. * cfl ) * a6(i,j+1) ) ) 
        else
         flux_y(i,j,k) = 0.0
        endif                                                 
                                                              
      enddo !i                                                
    enddo !j                                                  

! Update the tracer                                           
    do j=jsc,jec                                              
      do i=isc,iec                                            

        mass_mdppm(i,j,k) = mass_mdppm(i,j,k) + dtime * (   &
                Grd%dxtn(i,j-1) * Adv_vel%vhrho_nt(i,j-1,k) &
              - Grd%dxtn(i,j) * Adv_vel%vhrho_nt(i,j,k) )

        tracermass_mdppm(i,j,k) = tracermass_mdppm(i,j,k) + dtime * ( flux_y(i,j-1,k) - flux_y(i,j,k) )

        if (mass_mdppm(i,j,k)>0.) tracer_mdppm(i,j,k) = tracermass_mdppm(i,j,k) / mass_mdppm(i,j,k)

      enddo !i                                                
    enddo !j                                                  

! Calculate the overall tendency                              
    do j=jsc,jec                                              
      do i=isc,iec                                            
                                                              
! By convention, advection is applied as a negative tendency (i.e., convergence)
        advect_tracer_mdppm_test(i,j,k) = tmask_mdppm(i,j,k) * (       &
                    Thickness%rho_dzt(i,j,k,tau) * Tracer_field(i,j,k) &
                - tracermass_mdppm(i,j,k) * Grd%datr(i,j) ) / dtime

      enddo !i                                                
    enddo !j                                                  
                                                              
                                                              
  enddo !k                                                    

  call mpp_clock_end(id_clock_mdppm_test)

end function advect_tracer_mdppm_test
! </FUNCTION> NAME="advect_tracer_mdppm_test"


!#######################################################################
! <FUNCTION NAME="advect_tracer_mdppm">
!
! <DESCRIPTION>
!
!-------
! NOTE: This is routine may suffer from same problems as the 
! advect_tracer_mdfl_sweby routine. Namely, it has problems  
! with limiters that are improperly implemented in the case of generalized 
! vertical level coordinates of MOM. Consequently,  
! under certain circumstances, the resulting tracer concentration can 
! exhibit non-monotonic behaviour (i.e., produce extrema).  
!
! Sept 2011 
! Stephen.Griffies
! Alistair.Adcroft
!-------
!
! Compute advective tendency of dzt*rho*tracer concentration using 
! multi-dimensional piecewise parabolic method. 
! 
! Controlling parameters:
!   Tracer%ppm_hlimiter=0  -> No limiting in the horizontal
!          ppm_hlimiter=1  -> Full constraint (Colella and Woodward, 1984)
!          ppm_hlimiter=2  -> Improved full constraint (Lin, 2004)
!          ppm_hlimiter=3  -> Huynh limiter (Huynh, 1996)
!
!          ppm_vlimiter=*  -> as for ppm_hlimit but for the vertical
!
! Coded by Alistair.Adcroft
! Jan/Feb 2006
!
! Updated with controlling parameters, April 2006
!
! </DESCRIPTION>
!
function advect_tracer_mdppm(Time, Adv_vel, Tracer, Thickness, Tracer_field, dtime)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field
  real,                         intent(in) :: dtime

  real,dimension(isc:iec,jsc:jec,nk)       :: advect_tracer_mdppm

  real,dimension(isc:iec,jsc:jec)          :: ftp
  real,dimension(isc:iec,jsc:jec)          :: fbt
  real,dimension(isc:iec,jsc:jec)          :: wkm1

  integer                                  :: i, j, k
  integer                                  :: kp1, kp2, km1, km2
  integer                                  :: tau, taum1

  real                                     :: cfl, massflux, dMx, dMn
  real                                     :: Sim2,Sim1,Si,Sip1,Sip2

  real                                     :: da2,da4,da3m,da3p
  real                                     :: mskm2,mskm1,msk0,mskp1,mskp2
  real,parameter                           :: oneSixth=1./6., r12=1./12.
  real,parameter                           :: twoThirds=2./3.
  real,dimension(isc-4:iec+4,jsc-4:jec+4)  :: da, aL, aR, a6, d1m, d1p, d1mm, d1pp
  real,dimension(isc-4:iec+4,jsc-4:jec+4,nk) :: dak

  call mpp_clock_begin(id_clock_mdppm)

  ftp          = 0.0
  fbt          = 0.0
  wkm1         = 0.0
  tracer_mdppm = 0.0
  flux_x       = 0.0
  flux_y       = 0.0

  ! initialize arrays with extended halos 
  da           = 0.0
  aL           = 0.0
  aR           = 0.0
  a6           = 0.0
  d1m          = 0.0
  d1p          = 0.0
  d1mm         = 0.0
  d1pp         = 0.0
  dak          = 0.0

  tau   = Time%tau
  taum1 = Time%taum1

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tracer_mdppm(i,j,k) = Tracer_field(i,j,k)
        enddo
     enddo
  enddo

!
! Calculate vertical flux at the top and bottom of the boxes
!

  ! This block calculates the slope in each cell (aka PLM but
  ! with a fourth-order estimate)
  do k=1,nk
     km2 = max(k-2,1)
     km1 = max(k-1,1)
     kp1 = min(k+1,nk)
     kp2 = min(k+2,nk)
     do j=jsc,jec
        do i=isc,iec

        Sip2 = tracer_mdppm(i,j,kp2)
        Sip1 = tracer_mdppm(i,j,kp1)
        Si   = tracer_mdppm(i,j,k)
        Sim1 = tracer_mdppm(i,j,km1)
        Sim2 = tracer_mdppm(i,j,km2)
        ! Simple unmasked 4th order
        da4 = r12 * ( 8. * ( Sip1 - Sim1 ) + ( Sim2 - Sip2 ) )
        ! Simple 2nd order
        da2 = 0.5 * ( Sip1 - Sim1 )
        ! Simple 3rd order, biased to left
        da3m = r12 * ( 2. * Sim2 - 12. * Sim1 + 6. * Si + 4. *Sip1 )
        ! Simple 3rd order, biased to right
        da3p = r12 * ( - 4. * Sim1 - 6. * Si + 12. *Sip1  - 2. * Sip2 )
        ! Estimate of slope: dak(i,j,k) is the twice the "mismatch" as defined
        ! by Lin which itself is only half of the slope. The factor of two
        ! in the mismatch is confusing so we use simple estimate of slope.
        ! Note: The temperature masks used here are proxies for the advective
        ! flow masks that would be used if they existed. The only functional
        ! difference is that this scheme will not work with thin walls. 
        mskp2 = tmask_mdppm(i,j,kp2) * real(kp2-kp1)
        mskp1 = tmask_mdppm(i,j,kp1) * real(kp1-k)
        msk0  = tmask_mdppm(i,j,k)
        mskm1 = tmask_mdppm(i,j,km1) * real(k-km1)
        mskm2 = tmask_mdppm(i,j,km2) * real(km1-km2)
        dak(i,j,k) = ( (   mskm2 * ( 1. - 0.5 * mskp2 ) * da3m       &
                         + mskp2 * ( 1. - 0.5 * mskm2 ) * da3p )     &
                       + ( 1. - mskm2 ) * ( 1. - mskp2 ) * da2 )
        !!! dak(i,j,k) = da4 ! *** UNMASKED ***
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max( Sip1, Sim1,  Si ) - Si
        dMn = Si - min( Sip1, Sim1, Si )
        dak(i,j,k) = sign(1.,dak(i,j,k)) &
                     * min( abs(dak(i,j,k)) , 2. * min( dMx , dMn ) )  &
                     * mskm1 * mskp1 * msk0
        enddo
     enddo
  enddo

  ! This block estimates the upper (L) and lower (R) values on cell boundaries
  ! in k-direction
  do k=1,nk
     km2 = max(k-2,1)
     km1 = max(k-1,1)
     kp1 = min(k+1,nk)
     kp2 = min(k+2,nk)
     do j=jsc,jec
        do i=isc,iec

        mskp2 = tmask_mdppm(i,j,kp2) * real(kp2-kp1)
        mskp1 = tmask_mdppm(i,j,kp1) * real(kp1-k)
        msk0  = tmask_mdppm(i,j,k)
        mskm1 = tmask_mdppm(i,j,km1) * real(k-km1)
        mskm2 = tmask_mdppm(i,j,km2) * real(km1-km2)
        Sip2 = tracer_mdppm(i,j,kp2)
        Sip1 = tracer_mdppm(i,j,kp1)
        Si   = tracer_mdppm(i,j,k)
        Sim1 = tracer_mdppm(i,j,km1)
        Sim2 = tracer_mdppm(i,j,km2)

        ! Calculate the second difference for use in the Huynh limiter later
        d1m(i,j)  = ( Si - Sim1 ) * mskm1
        d1p(i,j)  = ( Sip1 - Si ) * mskp1
        d1mm(i,j) = ( Sim1 - Sim2 ) * mskm1 * mskm2
        d1pp(i,j) = ( Sip2 - Sip1 ) * mskp1 * mskp2

        Sim1 = Si - d1m(i,j)
        Sip1 = Si + d1p(i,j)
        ! Left edge: Eq. B2 in Lin 1994, MWR (132)
        aL(i,j) = 0.5 * ( Sim1 + Si ) + oneSixth * ( dak(i,j,km1) - dak(i,j,k) )
        ! Right edge: Eq. B2 in Lin 1994, MWR (132)
        aR(i,j) = 0.5 * ( Si + Sip1 ) + oneSixth * ( dak(i,j,k) - dak(i,j,kp1) )        
        enddo
     enddo

     ! Limit the edge values, aL and aR, for monotonicity
     if (Tracer%ppm_hlimiter.eq.1) then
       call ppm_limit_cw84(isc,iec,jsc,jec, tracer_mdppm(isc-4,jsc-4,k), aL, aR)
     elseif (Tracer%ppm_hlimiter.eq.2) then
       call ppm_limit_ifc(isc,iec,jsc,jec, tracer_mdppm(isc-4,jsc-4,k), dak(isc-4,jsc-4,k), aL,aR)
     elseif (Tracer%ppm_hlimiter.eq.3) then
       call ppm_limit_sh(isc,iec,jsc,jec, &
                        tracer_mdppm(isc-4,jsc-4,k), &
                        d1m, d1p, d1mm, d1pp, aL, aR)
     else
       call mpp_error(FATAL,&
       '==>Error from ocean_tracer_advect_mod: must choose ppm_hlimiter=1,2 or 3 in field table.')
     endif

     ! Calculate flux at top or bottom of box depending on sign of flow
     do j=jsc,jec
        do i=isc,iec

           mskp1 = tmask_mdppm(i,j,kp1) * real(kp1-k)
           msk0  = tmask_mdppm(i,j,k)
           mskm1 = tmask_mdppm(i,j,km1) * real(k-km1)
           Si   = tracer_mdppm(i,j,k)

           ! Curvature in k-direction
           a6(i,j) = 6. * Si - 3. * ( aR(i,j) + aL(i,j) )

           massflux = Grd%dat(i,j) * Adv_vel%wrho_bt(i,j,k) * msk0 * mskp1
           if (massflux.le.0.0) then
             cfl = abs(Adv_vel%wrho_bt(i,j,k) * dtime / Thickness%rho_dzt(i,j,k,tau))
             flux_z(i,j,k) = massflux                                        &
                           * ( aR(i,j) + 0.5 * cfl * ( ( aL(i,j) - aR(i,j) ) &
                               + ( 1. - twoThirds * cfl ) * a6(i,j) ) )
           endif

           massflux = Grd%dat(i,j) * Adv_vel%wrho_bt(i,j,km1) * msk0 * mskm1
           if (massflux.gt.0.0) then
               cfl = abs(Adv_vel%wrho_bt(i,j,km1) * dtime / Thickness%rho_dzt(i,j,km1,tau))
               flux_z(i,j,km1) = massflux                                      &
                    * ( aL(i,j) + 0.5 * cfl * ( ( aR(i,j) - aL(i,j) ) &
                    + ( 1. - twoThirds * cfl ) * a6(i,j) ) ) 
           endif

        enddo
     enddo
  enddo

  do k=1,nk
     km1 = max(k-1,1)
     mskm1 = real(k-km1)
     do j=jsc,jec
        do i=isc,iec
        tracer_mdppm(i,j,k) = tracer_mdppm(i,j,k)                          &
             + dtime / Thickness%rho_dzt(i,j,k,tau) * (                    &
                Grd%datr(i,j) * ( flux_z(i,j,k) - mskm1 * flux_z(i,j,km1)) &
             + Tracer_field(i,j,k) *                                       &
               (mskm1 * Adv_vel%wrho_bt(i,j,km1) - Adv_vel%wrho_bt(i,j,k)) )
        enddo
     enddo

     ! update vertical velocity for next k-level
     do j=jsc,jec
        do i=isc,iec
           wkm1(i,j) = Adv_vel%wrho_bt(i,j,k)
        enddo
     enddo

  enddo ! end of k-loop
  call mpp_update_domains (tracer_mdppm, Dom_mdppm%domain2d, flags=XUPDATE)

!
! Calculate flux at the eastern wall of the boxes
!
  do k=1,nk
    do j=jsc,jec
      do i=isc-2,iec+2
        ! This block calculates the slope in each cell (aka PLM but
        ! with a fourth-order estimate)
        Sip2 = tracer_mdppm(i+2,j,k)
        Sip1 = tracer_mdppm(i+1,j,k)
        Si   = tracer_mdppm(i,j,k)
        Sim1 = tracer_mdppm(i-1,j,k)
        Sim2 = tracer_mdppm(i-2,j,k)
        ! Simple unmasked 4th order
        da4 = r12 * ( 8. * ( Sip1 - Sim1 ) + ( Sim2 - Sip2 ) )
        ! Simple 2nd order
        da2 = 0.5 * ( Sip1 - Sim1 )
        ! Simple 3rd order, biased to left
        da3m = r12 * ( 2. * Sim2 - 12. * Sim1 + 6. * Si + 4. *Sip1 )
        ! Simple 3rd order, biased to right
        da3p = r12 * ( - 4. * Sim1 - 6. * Si + 12. *Sip1  - 2. * Sip2 )
        ! Estimate of slope: da(i,j) is the twice the "mismatch" as defined
        ! by Lin which itself is only half of the slope. The factor of two
        ! in the mismatch is confusing so we use simple estimate of slope.
        ! Note: The temperature masks used here are proxies for the advective
        ! flow masks that would be used if they existed. The only functional
        ! difference is that this scheme will not work with thin walls. 
        mskm2 = tmask_mdppm(i-2,j,k)
        mskm1 = tmask_mdppm(i-1,j,k)
        msk0  = tmask_mdppm(i,j,k)
        mskp1 = tmask_mdppm(i+1,j,k)
        mskp2 = tmask_mdppm(i+2,j,k)
        da(i,j) = ( (   mskm2 * ( 1. - 0.5 * mskp2 ) * da3m       &
                      + mskp2 * ( 1. - 0.5 * mskm2 ) * da3p )     &
                    + ( 1. - mskm2 ) * ( 1. - mskp2 ) * da2 )
        !!! da(i,j) = da4 ! *** UNMASKED ***
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max( Sip1, Sim1,  Si ) - Si
        dMn = Si - min( Sip1, Sim1, Si )
        da(i,j) = sign(1.,da(i,j)) * min( abs(da(i,j)) , 2. * min( dMx , dMn ) )  &
                  * mskm1 * mskp1 * msk0
        ! Calculate the second difference for use in the Huynh limiter later
        d1m(i,j)  = ( Si - Sim1 ) * mskm1
        d1p(i,j)  = ( Sip1 - Si ) * mskp1
        d1mm(i,j) = ( Sim1 - Sim2 ) * mskm1 * mskm2
        d1pp(i,j) = ( Sip2 - Sip1 ) * mskp1 * mskp2
      enddo !i

      do i=isc-1,iec+1
        ! This block estimates the left and right values on cell boundaries
        ! in i-direction
        Si   = tracer_mdppm(i,j,k)
        Sim1 = Si - d1m(i,j)
        Sip1 = Si + d1p(i,j)
!       mskm1 = tmask_mdppm(i-1,j,k)
!       mskp1 = tmask_mdppm(i+1,j,k)
        ! Left edge: Eq. B2 in Lin 1994, MWR (132)
        aL(i,j) = 0.5 * ( Sim1 + Si ) + oneSixth * ( da(i-1,j) - da(i,j) )        
        ! ... and masks
!       aL(i,j) = mskm1 * aL(i,j) + ( 1. - mskm1 ) * Si
        ! Right edge: Eq. B2 in Lin 1994, MWR (132)
        aR(i,j) = 0.5 * ( Si + Sip1 ) + oneSixth * ( da(i,j) - da(i+1,j) )        
        ! ... and masks
!       aR(i,j) = mskp1 * aR(i,j) + ( 1. - mskp1 ) * Si
      enddo !i
    enddo !j

    ! Limit the edge values, aL and aR, for monotonicity
    if (Tracer%ppm_hlimiter.eq.1) then
      call ppm_limit_cw84(isc,iec,jsc,jec, tracer_mdppm(isc-4,jsc-4,k), aL, aR)
    elseif (Tracer%ppm_hlimiter.eq.2) then
      call ppm_limit_ifc(isc,iec,jsc,jec, tracer_mdppm(isc-4,jsc-4,k), da, aL, aR)
     elseif (Tracer%ppm_hlimiter.eq.3) then
      call ppm_limit_sh(isc,iec,jsc,jec, &
                        tracer_mdppm(isc-4,jsc-4,k), &
                        d1m, d1p, d1mm, d1pp, aL, aR)
    else
       call mpp_error(FATAL,&
       '==>Error from ocean_tracer_advect_mod: must choose ppm_hlimiter=1,2 or 3 in field table.')
    endif

    do j=jsc,jec
! NOTE TO SELF(AJA): SHOULD BE ABLE TO ELIMINATE A6
      do i=isc-1,iec+1
        ! Curvature in i-direction                            
        Si   = tracer_mdppm(i,j,k)
        a6(i,j) = 6. * Si - 3. * ( aR(i,j) + aL(i,j) )        
      enddo !i
                                                              
      do i=isc-1,iec                                          
                                                              
        ! Volume flux and advective flux at right edge of cell (i.e. i+1/2) 
        massflux = Grd%dyte(i,j) * Adv_vel%uhrho_et(i,j,k)        
        cfl = Adv_vel%uhrho_et(i,j,k) * dtime * 2.0                           &
              / ( (Thickness%rho_dzt(i,j,k,tau) + Thickness%rho_dzt(i+1,j,k,tau) )    &
              * Grd%dxte(i,j))                                
        if (massflux.ge.0.0) then                              
         flux_x(i,j,k) = massflux * tmask_mdppm(i,j,k) * tmask_mdppm(i+1,j,k) &
                   * ( aR(i,j) + 0.5 * cfl * ( ( aL(i,j) - aR(i,j) )          &
                       + ( 1. - twoThirds * cfl ) * a6(i,j) ) )   
        else                                                  
         flux_x(i,j,k) = massflux * tmask_mdppm(i,j,k) * tmask_mdppm(i+1,j,k) &
                   * ( aL(i+1,j) - 0.5 * cfl * ( ( aR(i+1,j) - aL(i+1,j) )    &
                       + ( 1. + twoThirds * cfl ) * a6(i+1,j) ) ) 
        endif                                                 
                                                              
      enddo !i                                                
    enddo !j                                                  
!                                                             
! Update the tracer                                           
!                                                             
    do j=jsc,jec                                              
      do i=isc,iec                                            
        tracer_mdppm(i,j,k) = tracer_mdppm(i,j,k)                             &
            + dtime * tmask_mdppm(i,j,k) * Grd%datr(i,j)                      &
              / Thickness%rho_dzt(i,j,k,tau)                                  &
            * (flux_x(i-1,j,k) - flux_x(i,j,k)                                &
            + Tracer_field(i,j,k)* (                                          &
                      Grd%dyte(i,j) *   Adv_vel%uhrho_et( i ,j,k)             &
                    - Grd%dyte(i-1,j) * Adv_vel%uhrho_et(i-1,j,k) ) )
      enddo !i                                                
    enddo !j                                                  
  enddo !k                                                    
                                                              
! Update the 4-point tracer halo                              
  call mpp_update_domains (tracer_mdppm, Dom_mdppm%domain2d, flags=YUPDATE) 


! No flux at bottom (for final divergence of fluxes below)
  do j=jsc,jec
    do i=isc,iec
      wkm1(i,j) = 0.0
    enddo !i
  enddo !j

!
! Calculate flux at the northern wall of the boxes
!
  do k=1,nk
    do j=jsc-2,jec+2
      do i=isc,iec
        ! This block calculates the slope in each cell (aka PLM but
        ! with a fourth-order estimate)
        Sip2 = tracer_mdppm(i,j+2,k)
        Sip1 = tracer_mdppm(i,j+1,k)
        Si   = tracer_mdppm(i,j,k)                            
        Sim1 = tracer_mdppm(i,j-1,k)
        Sim2 = tracer_mdppm(i,j-2,k)
        ! Simple unmasked 4th order
        da4 = r12 * ( 8. * ( Sip1 - Sim1 ) + ( Sim2 - Sip2 ) )
        ! Simple 2nd order
        da2 = 0.5 * ( Sip1 - Sim1 )
        ! Simple 3rd order, biased to left
        da3m = r12 * ( 2. * Sim2 - 12. * Sim1 + 6. * Si + 4. *Sip1 )
        ! Simple 3rd order, biased to right
        da3p = r12 * ( - 4. * Sim1 - 6. * Si + 12. *Sip1  - 2. * Sip2 )
        ! Estimate of slope: da(i,j) is the twice the "mismatch" as defined
        ! by Lin which itself is only half of the slope. The factor of two
        ! in the mismatch is confusing so we use simple estimate of slope.
        ! Note: The temperature masks used here are proxies for the advective
        ! flow masks that would be used if they existed. The only functional
        ! difference is that this scheme will not work with thin walls. 
        mskm2 = tmask_mdppm(i,j-2,k)
        mskm1 = tmask_mdppm(i,j-1,k)
        msk0  = tmask_mdppm(i,j,k)
        mskp1 = tmask_mdppm(i,j+1,k)
        mskp2 = tmask_mdppm(i,j+2,k)
        da(i,j) = ( (   mskm2 * ( 1. - 0.5 * mskp2 ) * da3m       &
                      + mskp2 * ( 1. - 0.5 * mskm2 ) * da3p )     &
                    + ( 1. - mskm2 ) * ( 1. - mskp2 ) * da2 )
        !!! da(i,j) = da4 ! *** UNMASKED ***
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max( Sip1, Sim1, Si ) - Si             
        dMn = Si - min( Sip1, Sim1, Si )             
        da(i,j) = sign(1.,da(i,j)) * min( abs(da(i,j)) , 2. * min( dMx , dMn ) )  &
                  * mskm1 * mskp1 * msk0
        ! Calculate the second difference for use in the Huynh limiter later
        d1m(i,j)  = ( Si - Sim1 ) * mskm1
        d1p(i,j)  = ( Sip1 - Si ) * mskp1
        d1mm(i,j) = ( Sim1 - Sim2 ) * mskm1 * mskm2
        d1pp(i,j) = ( Sip2 - Sip1 ) * mskp1 * mskp2
      enddo !i                                                
    enddo !j                                                  
                                                              
    do j=jsc-1,jec+1                                          
      do i=isc,iec                                            
        ! This block estimates the left and right values on cell boundaries
        ! in j-direction
        Si   = tracer_mdppm(i,j,k)                            
        Sim1 = Si - d1m(i,j)
        Sip1 = Si + d1p(i,j)
!       mskm1 = tmask_mdppm(i-1,j,k)
!       mskp1 = tmask_mdppm(i+1,j,k)
        ! Left edge: Eq. B2 in Lin 1994, MWR (132)
        aL(i,j) = 0.5 * ( Sim1 + Si ) + oneSixth * ( da(i,j-1) - da(i,j) )        
        ! ... and masks
!       aL(i,j) = mskm1 * aL(i,j) + ( 1. - mskm1 ) * Si
        ! Right edge: Eq. B2 in Lin 1994, MWR (132)
        aR(i,j) = 0.5 * ( Si + Sip1 ) + oneSixth * ( da(i,j) - da(i,j+1) )        
        ! ... and masks
!       aR(i,j) = mskp1 * aR(i,j) + ( 1. - mskp1 ) * Si
      enddo !i                                                
    enddo !j                                                  

    ! Limit the edge values, aL and aR, for monotonicity
    if (Tracer%ppm_hlimiter.eq.1) then
      call ppm_limit_cw84(isc,iec,jsc,jec, tracer_mdppm(isc-4,jsc-4,k), aL, aR)
    elseif (Tracer%ppm_hlimiter.eq.2) then
      call ppm_limit_ifc(isc,iec,jsc,jec, tracer_mdppm(isc-4,jsc-4,k), da, aL, aR)
     elseif (Tracer%ppm_hlimiter.eq.3) then
       call ppm_limit_sh(isc,iec,jsc,jec, &
                        tracer_mdppm(isc-4,jsc-4,k), &
                        d1m, d1p, d1mm, d1pp, aL, aR)
    else
       call mpp_error(FATAL,&
       '==>Error from ocean_tracer_advect_mod: must choose ppm_hlimiter=1,2, or 3 in field table.')
    endif

! NOTE TO SELF(AJA): SHOULD BE ABLE TO ELIMINATE A6
    do j=jsc-1,jec+1                                          
      do i=isc,iec                                            
        ! Curvature in j-direction                            
        Si   = tracer_mdppm(i,j,k)
        a6(i,j) = 6. * Si - 3. * ( aR(i,j) + aL(i,j) )        
      enddo !i                                                
    enddo !j                                                  
                                                              
    do j=jsc-1,jec                                            
      do i=isc,iec                                            
                                                              
        ! Volume flux and advective flux at right edge of cell (i.e. i+1/2) 
        massflux = Grd%dxtn(i,j) * Adv_vel%vhrho_nt(i,j,k)        
        cfl = Adv_vel%vhrho_nt(i,j,k) * dtime * 2.0                              &
              / ( (Thickness%rho_dzt(i,j,k,tau) + Thickness%rho_dzt(i,j+1,k,tau))     &
              * Grd%dytn(i,j))                                
        if (massflux.ge.0.0) then                              
         flux_y(i,j,k) = massflux * tmask_mdppm(i,j,k) * tmask_mdppm(i,j+1,k) &
                   * ( aR(i,j) + 0.5 * cfl * ( ( aL(i,j) - aR(i,j) )          &
                       + ( 1. - 2./3. * cfl ) * a6(i,j) ) )   
        else                                                  
         flux_y(i,j,k) = massflux * tmask_mdppm(i,j,k) * tmask_mdppm(i,j+1,k) &
                   * ( aL(i,j+1) - 0.5 * cfl * ( ( aR(i,j+1) - aL(i,j+1) )    &
                       + ( 1. + 2./3. * cfl ) * a6(i,j+1) ) ) 
        endif                                                 
                                                              
      enddo !i                                                
    enddo !j                                                  

! Update the tracer                                           
    do j=jsc,jec                                              
      do i=isc,iec                                            
        tracer_mdppm(i,j,k) = tracer_mdppm(i,j,k)                             &
            + dtime * Grd%tmask(i,j,k) * Grd%datr(i,j)                        &
              / Thickness%rho_dzt(i,j,k,tau)                                  &
         * (flux_y(i,j-1,k) - flux_y(i,j,k))                              
      enddo !i                                                
    enddo !j                                                  

    ! calculate the overall tendency                              
    do j=jsc,jec                                              
      do i=isc,iec                                            
        tracer_mdppm(i,j,k) = tracer_mdppm(i,j,k)                             &
           + dtime * Tracer_field(i,j,k) / Thickness%rho_dzt(i,j,k,tau)       &
                  * ( Adv_vel%wrho_bt(i,j,k) - wkm1(i,j)                      &
                      + Grd%datr(i,j)*(                                       &
                         ( Grd%dyte(i-1,j) * Adv_vel%uhrho_et(i-1,j,k)        &
                         - Grd%dyte( i ,j) * Adv_vel%uhrho_et( i ,j,k) ) ) )   
                                                              
        ! by convention, advection is applied as a convergence
        advect_tracer_mdppm(i,j,k) = - Thickness%rho_dzt(i,j,k,tau) *  &
          ( tracer_mdppm(i,j,k) - Tracer_field(i,j,k) )                &
            / dtime * tmask_mdppm(i,j,k)                       
      enddo                                                
    enddo                                                  

    ! deep a copy of vertical velocity for next level             
    do j=jsc,jec                                              
      do i=isc,iec                                            
        wkm1(i,j) = Adv_vel%wrho_bt(i,j,k)                       
      enddo                                                 
    enddo                                                   
                                                              
  enddo                                                    

  call mpp_clock_end(id_clock_mdppm)


end function advect_tracer_mdppm
! </FUNCTION> NAME="advect_tracer_mdppm"


!#######################################################################
! <SUBROUTINE NAME="ppm_limit_cw84">
!
! <DESCRIPTION>
!
! Kernel to limit the edge values for PPM following Colella and Woodward, 1984
! 
! Coded by Alistair.Adcroft
! Apr 2006
!
! </DESCRIPTION>
!
subroutine ppm_limit_cw84(isc,iec,jsc,jec, tracer, aL, aR)
implicit none
! Arguments
integer, intent(in) :: isc,iec,jsc,jec
real, dimension(isc-4:iec+4,jsc-4:jec+4), intent(in)    :: tracer
real, dimension(isc-4:iec+4,jsc-4:jec+4), intent(inout) :: aL, aR
! Local
real    :: da2,da3m,da3p,da4,Si
integer :: i,j

    do j=jsc-4,jec+4
      do i=isc-4,iec+4
        ! This block monotonizes the parabola by adjusting the left and right values
        ! Limiter from Colella and Woodward, 1984, Eq. 1.10
        Si   = tracer(i,j)
        if ( ( aR(i,j) - Si ) * ( Si - aL(i,j) ) .le. 0. ) then
           aL(i,j) = Si
           aR(i,j) = Si
        endif
        da2 = aR(i,j) - aL(i,j)            ! Difference of edge values
        da4 = 0.5 * ( aR(i,j) + aL(i,j) )  ! Mean of edge values
        da3m = 6. * da2 * ( Si - da4 )
        da3p = da2 * da2
        if ( da3m .gt.  da3p ) aL(i,j) = 3. * Si - 2. * aR(i,j)
        if ( da3m .lt. -da3p ) aR(i,j) = 3. * Si - 2. * aL(i,j)
      enddo !i
    enddo !j

end subroutine ppm_limit_cw84
! </SUBROUTINE> NAME="ppm_limit_cw84"


!#######################################################################
! <SUBROUTINE NAME="ppm_limit_ifc">
!
! <DESCRIPTION>
!
! Kernel to limit the edge values for PPM using the Improved Full Constraint
! (IFC) of Lin, 2004.
! 
! Coded by Alistair.Adcroft
! Apr 2006
!
! </DESCRIPTION>
!
subroutine ppm_limit_ifc(isc,iec,jsc,jec, tracer, da, aL, aR)
implicit none
! Arguments
integer, intent(in) :: isc,iec,jsc,jec
real, dimension(isc-4:iec+4,jsc-4:jec+4), intent(in)    :: tracer, da
real, dimension(isc-4:iec+4,jsc-4:jec+4), intent(inout) :: aL, aR
! Local
real    :: Si,ada,sda
integer :: i,j

    do j=jsc-4,jec+4
      do i=isc-4,iec+4
        ! This block monotonizes the parabola by adjusting the left and right values
        ! Improved full constraint (IFC) limiter from Lin, 2004
        Si  = tracer(i,j)
        ada = abs(da(i,j))
        sda = sign(1., da(i,j))
        aL(i,j) = Si - sda * min( ada, abs(aL(i,j)-Si) )
        aR(i,j) = Si + sda * min( ada, abs(aR(i,j)-Si) )
      enddo !i
    enddo !j

end subroutine ppm_limit_ifc
! </SUBROUTINE> NAME="ppm_limit_ifc"


!#######################################################################
! <SUBROUTINE NAME="ppm_limit_sh">
!
! <DESCRIPTION>
!
! Kernel to limit the edge values for PPM following the monotonicity-
! preserving approach of Suresh and Huynh, 1997.
! 
! Coded by Alistair.Adcroft
! Apr 2006
!
! </DESCRIPTION>
!
! <NOTE>
! About efficiency: ppm_limit_sh() is not as efficient as it would be if
! we wrote a s/r for each direction. The pre-calculation of d1m, d1p, d1mm and
! d1pp duplicates operations and would be much faster if d1m was replaced
! by d1p(i+1). However, in order to re-use this limiter for the all directions
! (to simplify debugging) I have opted for the less efficient form for now. - AJA
! </NOTE>
subroutine ppm_limit_sh(isc,iec,jsc,jec, tracer, d1m, d1p, d1mm, d1pp, aL, aR)
implicit none
! Arguments
integer, intent(in) :: isc,iec,jsc,jec
real, dimension(isc-4:iec+4,jsc-4:jec+4), intent(in)    :: tracer, d1m, d1p, d1mm, d1pp
real, dimension(isc-4:iec+4,jsc-4:jec+4), intent(inout) :: aL, aR
! Local
real    :: Si,Sim1,Sip1
real    :: x,y,z,w,dM4m,dM4p,qAV,qMP,qUL,qLC,qMD,qMin,qMax
real, parameter :: fourThirds = 4./3.
integer :: i,j

    do j=jsc-1,jec+1
      do i=isc-1,iec+1
        ! This block monotonizes the parabola by adjusting the left and right values
        ! Limiter from Suresh and Huynh, 1997
        Si   = tracer(i,j)
        Sim1 = Si - d1m(i,j)               ! Sim1 = tracer(i-1,j)
        Sip1 = Si + d1p(i,j)               ! Sip1 = tracer(i+1,j)

        z = d1m(i,j) - d1mm(i,j)           ! z = del2(i-1,j)
        w = d1p(i,j) - d1m(i,j)            ! w = del2(i,j)
        x = 4. * w - z
        y = 4. * z - w
        dM4m = max(0., min(x,y,z,w)) + min(0., max(x,y,z,w))
        z = d1pp(i,j) - d1p(i,j)           ! z = del2(i+1,j)
        x = 4. * w - z
        y = 4. * z - w
        dM4p = max(0., min(x,y,z,w)) + min(0., max(x,y,z,w))

        qAV = 0.5 * ( Si + Sip1 )
        x = d1p(i,j)                       ! = Sip1 - Si
        y = 3. * d1m(i,j)                  ! = 3. * ( Si - Sim1 )
        qMP = Si + max(0., min(x, y)) + min(0., max(x,y))
        qUL = Si + y
        qLC = ( Si + 0.5 * d1m(i,j) ) + fourThirds * dM4m
        qMD = qAV - 0.5 * dM4p
        qMin = max( min(qMD,Si,Sip1), min(Si,qUL,qLC) )
        qMax = min( max(qMD,Si,Sip1), max(Si,qUL,qLC) )
        if ( (aR(i,j)-Si)*(aR(i,j)-qMP).gt.1.e-10 )  &
             aR(i,j) = min( max( aR(i,j), qMin ), qMax )

        qAV = 0.5 * ( Si + Sim1 )
        x = -d1m(i,j)                      ! = Sim1 - Si
        y = -3. * d1p(i,j)                 ! = 3. * ( Si - Sip1 )
        qMP = Si + max(0., min(x, y)) + min(0., max(x,y))
        qUL = Si + y
        qLC = ( Si - 0.5 * d1p(i,j) ) + fourThirds * dM4p
        qMD = qAV - 0.5 * dM4m
        qMin = max( min(qMD,Si,Sim1), min(Si,qUL,qLC) )
        qMax = min( max(qMD,Si,Sim1), max(Si,qUL,qLC) )
        if ( (aL(i,j)-Si)*(aL(i,j)-qMP).gt.1.e-10 )  &
             aL(i,j) = min( max( aL(i,j), qMin ), qMax )
      enddo !i
    enddo !j

end subroutine ppm_limit_sh
! </SUBROUTINE> NAME="ppm_limit_sh"


!#######################################################################
! <FUNCTION NAME="advect_tracer_mdmdt_test">
!
! <DESCRIPTION>
!
! Compute advective tendency of dzt*rho*tracer concentration using 
! multi-dimensional piecewise parabolic method. 
! 
! Controlling parameters:
!   Tracer%mdt_scheme=0    -> 7th order FV, unlimited (Daru & Tenaud, 2004)
!          mdt_scheme=1    -> 7th order FV, TVD limiter
!          mdt_scheme=2    -> 7th order FV, MP limiter
!          mdt_scheme=3    -> 7th order MP, 4th order TVD, FCT combination (Adcroft, 2011)
!
! Coded by Alistair.Adcroft
! Aug/Sep 2011
!
! This code is unsupported test code, which is not meant for general use.  
!
! </DESCRIPTION>
!
function advect_tracer_mdmdt_test(Time, Adv_vel, Tracer, Thickness, Tracer_field, dtime)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field
  real,                         intent(in) :: dtime

  real,dimension(isc:iec,jsc:jec,nk)       :: advect_tracer_mdmdt_test

  real,dimension(isc:iec,jsc:jec)          :: ftp
  real,dimension(isc:iec,jsc:jec)          :: fbt
  real,dimension(isc:iec,jsc:jec)          :: wkm1

  integer                                  :: i, j, k
  integer                                  :: kp1, kp2, kp3, kp4, km1, km2, km3
  integer                                  :: tau, taum1

  real                                     :: cfl, massflux

#ifdef DEBUG_OS7MP
  real :: tmin0,tmax0
#endif

  real :: Phi, Qippp, Qipp, Qip, Qi, Qim, Qimm, Qimmm
  real :: MskPPP, MskPP, MskP, Msk, MskM, MskMM, MskMMM, MskC
  real :: Fac, rp1h_cfl
  real :: DelP, DelM, DelPP, DelMM, DelPPP, DelMMM
  real :: Del2, Del2P, Del2M, Del2PP, Del2MM
  real :: Del3P, Del3M, Del3PP, Del3MM
  real :: Del4, Del4P, Del4M
  real :: Del5P, Del5M
  real :: Del6

#undef  DEBUG_OS7MP

  call mpp_clock_begin(id_clock_mdmdt_test)

  ftp         = 0.0
  fbt         = 0.0
  wkm1        = 0.0
  tracer_mdmdt = 0.0
  mass_mdmdt   = 0.0
  tracermass_mdmdt   = 0.0
  flux_x      = 0.0
  flux_y      = 0.0

  tau   = Time%tau
  taum1 = Time%taum1

#ifdef DEBUG_OS7MP
  call get_tracer_stats(Tracer_field(isc:iec,jsc:jec,:),tmask_mdmdt(isc:iec,jsc:jec,:),tmin0,tmax0)
#endif

  do k=1,nk

    do j=jsc,jec
      do i=isc,iec
        tracer_mdmdt(i,j,k) = Tracer_field(i,j,k)
        mass_mdmdt(i,j,k) = Thickness%rho_dzt(i,j,k,tau) * Grd%dat(i,j)
        tracermass_mdmdt(i,j,k) = mass_mdmdt(i,j,k) * Tracer_field(i,j,k)
      enddo
    enddo

    kp1 = min(k+1,nk)
    kp2 = min(k+2,nk)
    kp3 = min(k+3,nk)
    kp4 = min(k+4,nk)
    km1 = max(k-1,1)
    km2 = max(k-2,1)
    km3 = max(k-3,1)

    ! calculate flux at bottom of box and update the tracer
    do j=jsc,jec
      do i=isc,iec

        massflux = Grd%dat(i,j) * Adv_vel%wrho_bt(i,j,k)

        if (massflux==0.0) then                              
          fbt(i,j) = 0.0
        else ! cfl<>0

          ! CFL is based on mass of upwind cell
          if (massflux*Thickness%rho_dzt(i,j,k,tau)<0.0) then                              
            cfl = abs(Adv_vel%wrho_bt(i,j,k)) * dtime / Thickness%rho_dzt(i,j,k,tau)
            Qippp = Tracer_field(i,j,kp3)
            Qipp  = Tracer_field(i,j,kp2)
            Qip   = Tracer_field(i,j,kp1)
            Qi    = Tracer_field(i,j,k)
            Qim   = Tracer_field(i,j,km1)
            Qimm  = Tracer_field(i,j,km2)
            Qimmm = Tracer_field(i,j,km3)
            Msk   = tmask_mdmdt(i,j,k)
            MskP  = tmask_mdmdt(i,j,kp1) * Msk   * real(kp1-k)
            MskPP = tmask_mdmdt(i,j,kp2) * MskP  * real(kp2-kp1)
            MskPPP= tmask_mdmdt(i,j,kp3) * MskPP * real(kp3-kp2)
            MskM  = tmask_mdmdt(i,j,km1) * Msk   * real(k-km1)
            MskMM = tmask_mdmdt(i,j,km2) * MskM  * real(km1-km2)
            MskMMM= tmask_mdmdt(i,j,km3) * MskMM * real(km2-km3)
          elseif (massflux*Thickness%rho_dzt(i,j,kp1,tau)>0.0) then
            cfl = abs(Adv_vel%wrho_bt(i,j,k)) * dtime / Thickness%rho_dzt(i,j,kp1,tau)
            Qippp = Tracer_field(i,j,km2)
            Qipp  = Tracer_field(i,j,km1)
            Qip   = Tracer_field(i,j,k)
            Qi    = Tracer_field(i,j,kp1)
            Qim   = Tracer_field(i,j,kp2)
            Qimm  = Tracer_field(i,j,kp3)
            Qimmm = Tracer_field(i,j,kp4)
            Msk   = tmask_mdmdt(i,j,kp1)         * real(kp1-k)
            MskP  = tmask_mdmdt(i,j,k)   * Msk
            MskPP = tmask_mdmdt(i,j,km1) * MskP  * real(k-km1)
            MskPPP= tmask_mdmdt(i,j,km2) * MskPP * real(km1-km2)
            MskM  = tmask_mdmdt(i,j,kp2) * Msk   * real(kp2-kp1)
            MskMM = tmask_mdmdt(i,j,kp3) * MskM  * real(kp3-kp2)
            MskMMM= tmask_mdmdt(i,j,kp4) * MskMM * real(kp4-kp3)
          else
            stop 'Error A: zero mass in ocean cell in routine advect_tracer_mdmdt_test?'
          endif

#ifdef DEBUG_OS7MP
          if (cfl>0.5) then
            write(0,*) 'w-cfl>1/2'
          endif
          if (cfl>1.0) then
            stop 'w-cfl>1'
          endif
#endif

          ! 1st order method
          Phi = 0.
          MskC = Msk
          ! 2nd order correction [k+1 k]
          Fac = 1.0
          DelP = (Qip - Qi) * MskP
          MskC = MskC * MskP
          Phi = Phi + Fac * DelP * MskC
          ! 3rd order correction [k+1 k k-1]
          Fac = Fac * (cfl + 1.0) / 3.0
          DelM = (Qi - Qim) * MskM
          Del2 = DelP - DelM
          MskC = MskC * MskM
          Phi = Phi - Fac * Del2 * MskC
          ! 4th order correction [k+1 k k-1 k-2]
          Fac = Fac * (cfl - 2.0)/4.0
          DelMM = (Qim - Qimm) * MskMM
          Del2M = DelM - DelMM
          Del3M = Del2 - Del2M
          MskC = MskC * MskMM
          Phi = Phi + Fac * Del3M * MskC
          ! 5th order correction [k+2 k+1 k k-1 k-2]
          Fac = Fac * (cfl - 3.0)/5.0
          DelPP = (Qipp - Qip) * MskPP
          Del2P = DelPP - DelP
          Del3P = Del2P - Del2
          Del4 = Del3P - Del3M
          MskC = MskC * MskPP
          Phi = Phi - Fac * Del4 * MskC
          ! 6th order correction [k+3 k+2 k+1 k k-1 k-2]
          Fac = Fac * (cfl + 2.0)/6.0
          DelPPP = (Qippp - Qipp) * MskPPP
          Del2PP = DelPP - DelP
          Del3PP = Del2PP - Del2P
          Del4P = Del3PP - Del3P
          Del5P = Del4P - Del4
          MskC = MskC * MskPPP
          Phi = Phi + Fac * Del5P * MskC
          ! 7th order correction [k+3 k+2 k+1 k k-1 k-2 k-3]
          Fac = Fac * (cfl + 3.0)/7.0
          DelMMM = (Qimm - Qimmm) * MskMMM
          Del2MM = DelMM - DelMMM
          Del3MM = Del2M - Del2MM
          Del4M = Del3M - Del3MM
          Del5M = Del4 - Del4M
          Del6 = Del5P - Del5M
          MskC = MskC * MskMMM
          Phi = Phi - Fac * Del6 * MskC

          ! Convert Phi into the accuracy function (cf. DT04 eqs. 6 an 7, 10 11)
          if (DelP == 0.0) then
            Phi = 1.0E30
            rp1h_cfl = 1.0E30
          else
            Phi = Phi / DelP
            rp1h_cfl = (DelM / DelP) / cfl
          endif

          ! TVD limiter (DT04 eq. 16)
          if (Tracer%mdt_scheme==1) &
            Phi = max(0.0, min(2.0 / (1.0 - cfl), Phi, 2.0*rp1h_cfl ))

          ! Flux takes form of limited LW flux (DT04, eq. 8)
          Phi = Phi * 0.5 * (1.0 - cfl)
          fbt(i,j) = massflux * (Qi + Phi * DelP)  !!!!! * MskP

        endif ! cfl<>0
 
        mass_mdmdt(i,j,k) = mass_mdmdt(i,j,k) + dtime * Grd%dat(i,j) * ( &
                    Adv_vel%wrho_bt(i,j,k) - wkm1(i,j)  )

#ifdef DEBUG_OS7MP
        if (mass_mdmdt(i,j,k)<0.0) then
          stop 'w: negative mass'
        endif
#endif
 
        tracermass_mdmdt(i,j,k) = tracermass_mdmdt(i,j,k) + dtime * ( fbt(i,j) - ftp(i,j) )

        if (mass_mdmdt(i,j,k)>0.) tracer_mdmdt(i,j,k) = tracermass_mdmdt(i,j,k) / mass_mdmdt(i,j,k)

        flux_z(i,j,k) = fbt(i,j)
        ftp(i,j)      = fbt(i,j)
        wkm1(i,j) = Adv_vel%wrho_bt(i,j,k)

        enddo
     enddo

  enddo ! end of k-loop
  call mpp_update_domains (tracer_mdmdt, Dom_mdmdt%domain2d, flags=XUPDATE)
  call mpp_update_domains (tracermass_mdmdt, Dom_mdmdt%domain2d, flags=XUPDATE)
  call mpp_update_domains (mass_mdmdt, Dom_mdmdt%domain2d, flags=XUPDATE)

#ifdef DEBUG_OS7MP
  call tracer_stats(tracer_mdmdt(isc:iec,jsc:jec,:),Tmin0,Tmax0,'mdmdtafter Z')
#endif

  ! calculate flux at the eastern wall of the boxes
  do k=1,nk

    do j=jsc,jec
      do i=isc-1,iec

        massflux = Grd%dyte(i,j) * Adv_vel%uhrho_et(i,j,k)

        if (massflux==0.0) then                              
          flux_x(i,j,k) = 0.0
        else ! cfl<>0

          if (massflux*mass_mdmdt(i,j,k)>0.0) then                              
            cfl = abs(massflux) * dtime / mass_mdmdt(i,j,k) 
            Qippp = tracer_mdmdt(i+3,j,k)
            Qipp  = tracer_mdmdt(i+2,j,k)
            Qip   = tracer_mdmdt(i+1,j,k)
            Qi    = tracer_mdmdt(i,j,k)
            Qim   = tracer_mdmdt(i-1,j,k)
            Qimm  = tracer_mdmdt(i-2,j,k)
            Qimmm = tracer_mdmdt(i-3,j,k)
            Msk   = tmask_mdmdt(i,j,k)
            MskP  = tmask_mdmdt(i+1,j,k) * Msk
            MskPP = tmask_mdmdt(i+2,j,k) * MskP
            MskPPP= tmask_mdmdt(i+3,j,k) * MskPP
            MskM  = tmask_mdmdt(i-1,j,k) * Msk
            MskMM = tmask_mdmdt(i-2,j,k) * MskM
            MskMMM= tmask_mdmdt(i-3,j,k) * MskMM
          elseif (massflux*mass_mdmdt(i+1,j,k)<0.0) then
            cfl = abs(massflux) * dtime / mass_mdmdt(i+1,j,k) 
            Qippp = tracer_mdmdt(i-2,j,k)
            Qipp  = tracer_mdmdt(i-1,j,k)
            Qip   = tracer_mdmdt(i,j,k)
            Qi    = tracer_mdmdt(i+1,j,k)
            Qim   = tracer_mdmdt(i+2,j,k)
            Qimm  = tracer_mdmdt(i+3,j,k)
            Qimmm = tracer_mdmdt(i+4,j,k)
            Msk   = tmask_mdmdt(i+1,j,k)
            MskP  = tmask_mdmdt(i,j,k) * Msk
            MskPP = tmask_mdmdt(i-1,j,k) * MskP
            MskPPP= tmask_mdmdt(i-2,j,k) * MskPP
            MskM  = tmask_mdmdt(i+2,j,k) * Msk
            MskMM = tmask_mdmdt(i+3,j,k) * MskM
            MskMMM= tmask_mdmdt(i+4,j,k) * MskMM
          else
            stop 'Error B: zero mass in ocean cell in routine advect_tracer_mdmdt_test?'
          endif

#ifdef DEBUG_OS7MP
          if (cfl>0.5) then
            write(0,*) 'u-cfl>1/2'
          endif
          if (cfl>1.0) then
            stop 'u-cfl>1'
          endif
#endif

          ! 1st order method
          Phi = 0.
          MskC = Msk
          ! 2nd order correction [i+1 i]
          Fac = 1.0
          DelP = (Qip - Qi) * MskP
          MskC = MskC * MskP
          Phi = Phi + Fac * DelP * MskC
          ! 3rd order correction [i+1 i i-1]
          Fac = Fac * (cfl + 1.0) / 3.0
          DelM = (Qi - Qim) * MskM
          Del2 = DelP - DelM
          MskC = MskC * MskM
          Phi = Phi - Fac * Del2 * MskC
          ! 4th order correction [i+1 i i-1 i-2]
          Fac = Fac * (cfl - 2.0)/4.0
          DelMM = (Qim - Qimm) * MskMM
          Del2M = DelM - DelMM
          Del3M = Del2 - Del2M
          MskC = MskC * MskMM
          Phi = Phi + Fac * Del3M * MskC
          ! 5th order correction [i+2 i+1 i i-1 i-2]
          Fac = Fac * (cfl - 3.0)/5.0
          DelPP = (Qipp - Qip) * MskPP
          Del2P = DelPP - DelP
          Del3P = Del2P - Del2
          Del4 = Del3P - Del3M
          MskC = MskC * MskPP
          Phi = Phi - Fac * Del4 * MskC
          ! 6th order correction [i+3 i+2 i+1 i i-1 i-2]
          Fac = Fac * (cfl + 2.0)/6.0
          DelPPP = (Qippp - Qipp) * MskPPP
          Del2PP = DelPP - DelP
          Del3PP = Del2PP - Del2P
          Del4P = Del3PP - Del3P
          Del5P = Del4P - Del4
          MskC = MskC * MskPPP
          Phi = Phi + Fac * Del5P * MskC
          ! 7th order correction [i+3 i+2 i+1 i i-1 i-2 i-3]
          Fac = Fac * (cfl + 3.0)/7.0
          DelMMM = (Qimm - Qimmm) * MskMMM
          Del2MM = DelMM - DelMMM
          Del3MM = Del2M - Del2MM
          Del4M = Del3M - Del3MM
          Del5M = Del4 - Del4M
          Del6 = Del5P - Del5M
          MskC = MskC * MskMMM
          Phi = Phi - Fac * Del6 * MskC

          ! Convert Phi into the accuracy function (cf. DT04 eqs. 6 an 7, 10 11)
          if (DelP == 0.0) then
            Phi = 1.0E30
            rp1h_cfl = 1.0E30
          else
            Phi = Phi / DelP
            rp1h_cfl = (DelM / DelP) / cfl
          endif

          ! TVD limiter (DT04 eq. 16)
          if (Tracer%mdt_scheme==1) &
            Phi = max(0.0, min(2.0 / (1.0 - cfl), Phi, 2.0*rp1h_cfl ))

          ! Flux takes form of limited LW flux (DT04, eq. 8)
          Phi = Phi * 0.5 * (1.0 - cfl)
          flux_x(i,j,k) =  massflux * (Qi + Phi * DelP) * MskP

        endif ! cfl<>0
      enddo
    enddo

    ! Update the intermediate tracer-mass and mass and diagnose value of tracer (Easter, 1993)
    do j=jsc,jec
      do i=isc,iec

        mass_mdmdt(i,j,k) = mass_mdmdt(i,j,k) + dtime * (          &
                       Grd%dyte(i-1,j) * Adv_vel%uhrho_et(i-1,j,k) &
                       - Grd%dyte(i,j) * Adv_vel%uhrho_et(i,j,k) )

#ifdef DEBUG_OS7MP
        if (mass_mdmdt(i,j,k)<0.0) then
          stop 'u: negative mass'
        endif
#endif

        tracermass_mdmdt(i,j,k) = tracermass_mdmdt(i,j,k) + dtime * (flux_x(i-1,j,k) - flux_x(i,j,k))

        if (mass_mdmdt(i,j,k)>0.) tracer_mdmdt(i,j,k) = tracermass_mdmdt(i,j,k) / mass_mdmdt(i,j,k)

      enddo
    enddo

  enddo ! end of k-loop
  call mpp_update_domains (tracer_mdmdt, Dom_mdmdt%domain2d, flags=YUPDATE)
  call mpp_update_domains (tracermass_mdmdt, Dom_mdmdt%domain2d, flags=YUPDATE)
  call mpp_update_domains (mass_mdmdt, Dom_mdmdt%domain2d, flags=YUPDATE)

#ifdef DEBUG_OS7MP
  call tracer_stats(tracer_mdmdt(isc:iec,jsc:jec,:),Tmin0,Tmax0,'mdmdtafter U')
#endif

  do k=1,nk

    do j=jsc-1,jec
      do i=isc,iec

        massflux = Grd%dxtn(i,j) * Adv_vel%vhrho_nt(i,j,k)

        if (massflux==0.0) then                              
          flux_y(i,j,k) = 0.0
        else ! cfl<>0

          if (massflux*mass_mdmdt(i,j,k)>0.0) then                              
            cfl = abs(massflux) * dtime / mass_mdmdt(i,j,k) 
            Qippp = tracer_mdmdt(i,j+3,k)
            Qipp  = tracer_mdmdt(i,j+2,k)
            Qip   = tracer_mdmdt(i,j+1,k)
            Qi    = tracer_mdmdt(i,j,k)
            Qim   = tracer_mdmdt(i,j-1,k)
            Qimm  = tracer_mdmdt(i,j-2,k)
            Qimmm = tracer_mdmdt(i,j-3,k)
            Msk   = tmask_mdmdt(i,j,k)
            MskP  = tmask_mdmdt(i,j+1,k) * Msk
            MskPP = tmask_mdmdt(i,j+2,k) * MskP
            MskPPP= tmask_mdmdt(i,j+3,k) * MskPP
            MskM  = tmask_mdmdt(i,j-1,k) * Msk
            MskMM = tmask_mdmdt(i,j-2,k) * MskM
            MskMMM= tmask_mdmdt(i,j-3,k) * MskMM
          elseif (massflux*mass_mdmdt(i,j+1,k)<0.0) then
            cfl = abs(massflux) * dtime / mass_mdmdt(i,j+1,k) 
            Qippp = tracer_mdmdt(i,j-2,k)
            Qipp  = tracer_mdmdt(i,j-1,k)
            Qip   = tracer_mdmdt(i,j,k)
            Qi    = tracer_mdmdt(i,j+1,k)
            Qim   = tracer_mdmdt(i,j+2,k)
            Qimm  = tracer_mdmdt(i,j+3,k)
            Qimmm = tracer_mdmdt(i,j+4,k)
            Msk   = tmask_mdmdt(i,j+1,k)
            MskP  = tmask_mdmdt(i,j,k) * Msk
            MskPP = tmask_mdmdt(i,j-1,k) * MskP
            MskPPP= tmask_mdmdt(i,j-2,k) * MskPP
            MskM  = tmask_mdmdt(i,j+2,k) * Msk
            MskMM = tmask_mdmdt(i,j+3,k) * MskM
            MskMMM= tmask_mdmdt(i,j+4,k) * MskMM
          else
            stop 'Error C: zero mass in ocean cell in routine advect_tracer_mdmdt_test?'
          endif

          ! 1st order method
          Phi = 0.
          MskC = Msk
          ! 2nd order correction [i+1 i]
          Fac = 1.0
          DelP = (Qip - Qi) * MskP
          MskC = MskC * MskP
          Phi = Phi + Fac * DelP * MskC
          ! 3rd order correction [i+1 i i-1]
          Fac = Fac * (cfl + 1.0) / 3.0
          DelM = (Qi - Qim) * MskM
          Del2 = DelP - DelM
          MskC = MskC * MskM
          Phi = Phi - Fac * Del2 * MskC
          ! 4th order correction [i+1 i i-1 i-2]
          Fac = Fac * (cfl - 2.0)/4.0
          DelMM = (Qim - Qimm) * MskMM
          Del2M = DelM - DelMM
          Del3M = Del2 - Del2M
          MskC = MskC * MskMM
          Phi = Phi + Fac * Del3M * MskC
          ! 5th order correction [i+2 i+1 i i-1 i-2]
          Fac = Fac * (cfl - 3.0)/5.0
          DelPP = (Qipp - Qip) * MskPP
          Del2P = DelPP - DelP
          Del3P = Del2P - Del2
          Del4 = Del3P - Del3M
          MskC = MskC * MskPP
          Phi = Phi - Fac * Del4 * MskC
          ! 6th order correction [i+3 i+2 i+1 i i-1 i-2]
          Fac = Fac * (cfl + 2.0)/6.0
          DelPPP = (Qippp - Qipp) * MskPPP
          Del2PP = DelPP - DelP
          Del3PP = Del2PP - Del2P
          Del4P = Del3PP - Del3P
          Del5P = Del4P - Del4
          MskC = MskC * MskPPP
          Phi = Phi + Fac * Del5P * MskC
          ! 7th order correction [i+3 i+2 i+1 i i-1 i-2 i-3]
          Fac = Fac * (cfl + 3.0)/7.0
          DelMMM = (Qimm - Qimmm) * MskMMM
          Del2MM = DelMM - DelMMM
          Del3MM = Del2M - Del2MM
          Del4M = Del3M - Del3MM
          Del5M = Del4 - Del4M
          Del6 = Del5P - Del5M
          MskC = MskC * MskMMM
          Phi = Phi - Fac * Del6 * MskC

          ! Convert Phi into the accuracy function (cf. DT04 eqs. 6 an 7, 10 11)
          if (DelP == 0.0) then
            Phi = 1.0E30
            rp1h_cfl = 1.0E30
          else
            Phi = Phi / DelP
            rp1h_cfl = (DelM / DelP) / cfl
          endif

          ! TVD limiter (DT04 eq. 16)
          if (Tracer%mdt_scheme==1) &
            Phi = max(0.0, min(2.0 / (1.0 - cfl), Phi, 2.0*rp1h_cfl ))

          ! Flux takes form of limited LW flux (DT04, eq. 8)
          Phi = Phi * 0.5 * (1.0 - cfl)
          flux_y(i,j,k) =  massflux * (Qi + Phi * DelP) * MskP

        endif ! cfl<>0
      enddo
    enddo

    ! Update the intermediate tracer-mass and mass and diagnose value of tracer (Easter, 1993)
    do j=jsc,jec
      do i=isc,iec

        mass_mdmdt(i,j,k) = mass_mdmdt(i,j,k) + dtime * (          &
                       Grd%dxtn(i,j-1) * Adv_vel%vhrho_nt(i,j-1,k) &
                       - Grd%dxtn(i,j) * Adv_vel%vhrho_nt(i,j,k) )

#ifdef DEBUG_OS7MP
        if (mass_mdmdt(i,j,k)<0.0) then
          stop 'v: negative mass'
        endif
#endif

        tracermass_mdmdt(i,j,k) = tracermass_mdmdt(i,j,k) + dtime * (flux_y(i,j-1,k) - flux_y(i,j,k))

        if (mass_mdmdt(i,j,k)>0.) tracer_mdmdt(i,j,k) = tracermass_mdmdt(i,j,k) / mass_mdmdt(i,j,k)

      enddo
    enddo

    ! calculate the overall tendency
    do j=jsc,jec
      do i=isc,iec

        advect_tracer_mdmdt_test(i,j,k) = tmask_mdmdt(i,j,k) * ( &
                    Thickness%rho_dzt(i,j,k,tau) * Tracer_field(i,j,k) &
                     - tracermass_mdmdt(i,j,k) * Grd%datr(i,j) ) / dtime

      enddo 
    enddo 

  enddo ! end of k-loop

#ifdef DEBUG_OS7MP
  call tracer_stats(tracer_mdmdt(isc:iec,jsc:jec,:),Tmin0,Tmax0,'mdmdt')
#endif

  call mpp_clock_end(id_clock_mdmdt_test)

end function advect_tracer_mdmdt_test
! </FUNCTION> NAME="advect_tracer_mdmdt_test"



!#######################################################################
! <SUBROUTINE NAME="gyre_overturn_diagnose">
!
! <DESCRIPTION>
! Diagnose tracer transport according to zonal mean and 
! deviations from zonal mean. 
!
! We allow for the use of a basin mask so that the zonal means
! are performed over the global domain and each of five other
! basins.  If no basin mask is read, then assume global domain 
! used to define the zonal means. 
! 
! []=zonal mean; *=deviation from zonal mean
!
! V = rho_dzt dxt T, with rho=rho0 for Boussinesq 
!
! int [V T]   = total advective 
! int [V][T]  = overturning
! int V* T*   = gyre
!
! This routine is much faster than the older routine gyre_overturn_diagnose_old.
!
! In this version, we compute partial zonal sums locally and then perform 
! a "global" sum over a restricted set of pes. 
! This approach results in slightly different answers at the 
! REAL(8) level. However, differences should be negligible for diagnostic 
! purposes at the REAL(4) level. 
! 
! We require no global arrays here, and redundant computation is eliminated.
! 
! original approach 
! Stephen.Griffies
! May 2007
!
! optimization by removing global arrays 
! russell.fiedler@csiro.au
! March 2012
!
! </DESCRIPTION>
  
  subroutine gyre_overturn_diagnose(Time, Adv_vel, Tracer, ntracer)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  integer,                      intent(in) :: ntracer

  integer :: i, j, k, tau, nbasin
  integer :: njnk

  real,dimension(jsc:jec,nk) :: avg_vhrho_nt, avg_tracer, avg_tracer_n
  real,dimension(jsc:jec,nk) :: basin_width       ! width at tracer centre
  real                       :: basin_width_r
  real,dimension(jsc:jec,nk) :: basin_width_jp1   ! width at tracer centre to the "north"
  real                       :: basin_width_jp1_r
  real,dimension(jsc:jec,nk) :: basin_widthn      ! width at tracer northern edger
  real                       :: basin_widthn_r
  real,dimension(jsc:jec)    :: merid_flux_advect_1d
  real,dimension(jsc:jec)    :: merid_flux_gyre_1d
  real,dimension(jsc:jec)    :: merid_flux_over_1d
  real,dimension(jsc:jec,nk) :: merid_flux_advect_2d

  call mpp_clock_begin(id_clock_gyre_over)

  tau = Time%tau
  njnk=nk*(jec-jsc+1)

  do nbasin=0,5

     if(id_merid_flux_advect(nbasin,ntracer) > 0 .or. &
        id_merid_flux_over(nbasin,ntracer)   > 0 .or. & 
        id_merid_flux_gyre(nbasin,ntracer)   > 0) then 

        merid_flux_advect_1d(:)   = 0.0
        merid_flux_advect_2d(:,:) = 0.0
        merid_flux_gyre_1d(:)     = 0.0
        merid_flux_over_1d(:)     = 0.0
        avg_vhrho_nt    = 0.0  
        avg_tracer      = 0.0  
        avg_tracer_n    = 0.0  
        basin_width     = 0.0
        basin_width_jp1 = 0.0
        basin_widthn    = 0.0

        if( do_this_basin(nbasin) ) then
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    basin_width(j,k)  = basin_width(j,k) + Grd%dxt(i,j)*Grd%tmask(i,j,k) &
                      *local_basin_mask(i,j,nbasin)
                    basin_width_jp1(j,k)  = basin_width_jp1(j,k) + Grd%dxt(i,j+1)*Grd%tmask(i,j+1,k) &
                      *local_basin_mask_jp1(i,j,nbasin)
                    basin_widthn(j,k) = basin_widthn(j,k) + Grd%dxtn(i,j)*Grd%tmask(i,j,k) &
                      *local_basin_mask(i,j,nbasin)
                    avg_vhrho_nt(j,k) = avg_vhrho_nt(j,k)  + Grd%dxtn(i,j)*Grd%tmask(i,j,k) &
                      *Grd%tmask(i,j+1,k)*local_basin_mask(i,j,nbasin)*Adv_vel%vhrho_nt(i,j,k)
                    avg_tracer(j,k)   = avg_tracer(j,k) + Grd%dxt(i,j)*Grd%tmask(i,j,k) &
                      *local_basin_mask(i,j,nbasin)* Tracer%field(i,j,k,tau)
                    avg_tracer_n(j,k) = avg_tracer_n(j,k) + Grd%dxt(i,j+1)*Grd%tmask(i,j+1,k) &
                     *local_basin_mask_jp1(i,j,nbasin)* Tracer%field(i,j+1,k,tau)

                     ! this can be done as a sum over 1D, but it gives a result different
                     ! to that using a 2D sum then summing in the vertical.
!                    merid_flux_advect_1d(j) = merid_flux_advect_1d(j) &
!                      + Grd%tmask(i,j,k)*local_basin_mask(i,j,nbasin)*flux_y(i,j,k) 
                    merid_flux_advect_2d(j,k) = merid_flux_advect_2d(j,k) &
                      + Grd%tmask(i,j,k)*local_basin_mask(i,j,nbasin)*flux_y(i,j,k) 
                 enddo
              enddo
           enddo
        endif

        ! these first 3 only need to be done once ever.
        ! do not need jp1 explicitly. Do sum through to jec+1
        call mpp_sum(basin_width ,njnk,pelist_x)

        ! is this really needed? See below.
        call mpp_sum(basin_widthn,njnk,pelist_x)
        call mpp_sum(basin_width_jp1,njnk,pelist_x)

        ! if doing regional sums then we could avoid these by 
        ! checking if any(basin_width(:,1) /= 0.0) ahead of time
        call mpp_sum(avg_vhrho_nt,njnk,pelist_x)
        call mpp_sum(avg_tracer  ,njnk,pelist_x)

        ! no need for this. we have jp1 calculated above.
        call mpp_sum(avg_tracer_n,njnk,pelist_x)

        ! 2D or not 2D, that is the question...
        ! call mpp_sum(merid_flux_advect_1d,jec+1-jsc+1,pelist_x)
        call mpp_sum(merid_flux_advect_2d,njnk,pelist_x)
        merid_flux_advect_1d(:)=sum(merid_flux_advect_2d,dim=2)

        do k=1,nk
           do j=jsc,jec
              if(basin_width_jp1(j,k) > 0.0) then 
                  basin_width_jp1_r = 1.0/basin_width_jp1(j,k)
                  avg_tracer_n(j,k)     = avg_tracer_n(j,k)*basin_width_jp1_r
              endif 
              if(basin_width(j,k) > 0.0) then 
                  basin_width_r     = 1.0/basin_width(j,k)
                  avg_tracer(j,k)     = avg_tracer(j,k)*basin_width_r

! I think the next 3 lines can be replaced by
!                  merid_flux_over_1d(j) = merid_flux_over_1d(j) &
!                  + sum_vhrho_nt(j,k)*0.5*(avg_tracer(j,k)+avg_tracer(j+1,k)) 
! no need to calculate avg_vhrho_nt
!
! i.e. we only need to compute the sum of vhrho_nt
!
                  basin_widthn_r    = 1.0/basin_widthn(j,k)
                  avg_vhrho_nt(j,k)      = avg_vhrho_nt(j,k)*basin_widthn_r
                  merid_flux_over_1d(j) = merid_flux_over_1d(j) &
                  + basin_widthn(j,k)*avg_vhrho_nt(j,k)*0.5*(avg_tracer(j,k)+avg_tracer_n(j,k)) 
              endif
           enddo
        enddo
        do j=jsc,jec
           merid_flux_gyre_1d(j) = merid_flux_advect_1d(j) - merid_flux_over_1d(j)
        enddo

        if (id_merid_flux_advect(nbasin,ntracer) > 0) then
             used = send_data(id_merid_flux_advect(nbasin,ntracer),        &
                    Tracer%conversion*merid_flux_advect_1d,Time%model_time,&
                    is_in=jsc, ie_in=jec)
        endif
        if (id_merid_flux_over(nbasin,ntracer) > 0) then
             used = send_data(id_merid_flux_over(nbasin,ntracer),        &
                    Tracer%conversion*merid_flux_over_1d,Time%model_time,&
                    is_in=jsc, ie_in=jec)
        endif
        if (id_merid_flux_gyre(nbasin,ntracer) > 0) then
             used = send_data(id_merid_flux_gyre(nbasin,ntracer),        &
                    Tracer%conversion*merid_flux_gyre_1d,Time%model_time,&
                    is_in=jsc, ie_in=jec)
        endif

     endif ! endif for id_merid_flux_advect or id_merid_flux_over or id_merid_flux_gyre

  enddo   ! enddo for nbasin 
    

  call mpp_clock_end(id_clock_gyre_over)

end subroutine gyre_overturn_diagnose
! </SUBROUTINE> NAME="gyre_overturn_diagnose"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag">
!
! <DESCRIPTION>
! Compute watermass diagnostics associated with resolved flow
! advection of temperature and salinity.
! </DESCRIPTION>
!
subroutine watermass_diag(Time, Dens)

  type(ocean_time_type),      intent(in)  :: Time
  type(ocean_density_type),   intent(in)  :: Dens 

  integer :: i,j,k,tau

  if(.not. compute_watermass_diag) return 

  tau = Time%tau


  ! both temperature and salinity contributions 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = neutral_temp_advect(i,j,k)*Dens%drhodT(i,j,k) &
                       + neutral_salt_advect(i,j,k)*Dens%drhodS(i,j,k)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k) 
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo
  call diagnose_3d(Time, Grd, id_neut_rho_advect, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_advect, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_advect, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_advect_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_advect_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_advect_on_nrho, wrk4)

  ! temperature contributions 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = neutral_temp_advect(i,j,k)*Dens%drhodT(i,j,k)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k) 
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo
  call diagnose_3d(Time, Grd, id_neut_temp_advect, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_advect, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_temp_advect, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_advect_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_advect_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_advect_on_nrho, wrk4)

  ! salinity contributions 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = neutral_salt_advect(i,j,k)*Dens%drhodS(i,j,k)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k) 
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo
  call diagnose_3d(Time, Grd, id_neut_salt_advect, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_advect, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_salt_advect, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_advect_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_advect_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_advect_on_nrho, wrk4)

end subroutine watermass_diag
! </SUBROUTINE> NAME="watermass_diag"



!#######################################################################
! <SUBROUTINE NAME="compute_adv_diss">
!
! <DESCRIPTION>
!
! Compute the dissipation due to advection truncation errors. 
! This diagnostic requires computation of advection operator acting
! on the squared tracer concentration. 
!
! NOTE: This scheme isolates the dissipation from trucation errors
! in advection ONLY for the following vertical coordinates:
! 1/ geopotential: for all k-levels, except for k=1 
! (due to undulating surface height)
! 2/ pressure: for all k-levels, except for k=kmt
! (due to undulating bottom pressure)
!
! NOTE: For the Quicker advection scheme, we assume the preferred
! two_level time scheme is used here, so that taum1=tau.  
!
! NOTE: If PSOM is used for temp or salt, then we MUST also enable 
! a new passive tracer in the field table, with name 
! passive_temp_sq and passive_salt_sq.  This extra tracer is 
! used for diagnostics alone, and is required due to the extra
! moment fields used for computing the PSOM tendency. 
! If PSOM is used for another tracer besides temp or salt, then
! some extra code needs to be written inside ocean_passive.F90,
! emulating the work done for temp and salt. 
!
! </DESCRIPTION>
!
subroutine compute_adv_diss(Time, Adv_vel, Thickness, T_prog, Tracer, ntracer, dtime)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_prog_tracer_type), intent(inout) :: Tracer
  integer,                      intent(in)    :: ntracer
  real,                         intent(in)    :: dtime

  integer :: tau, taup1
  integer :: i,j,k
  logical :: use_psom=.false.
  real    :: dtimer, term1, term2

  call mpp_clock_begin(id_clock_adv_diss)

  tau    = Time%tau 
  taup1  = Time%taup1
  dtimer = 1.0/dtime 

  wrk1 = 0.0
  wrk2 = 0.0
  wrk3 = 0.0
  wrk4 = 0.0

  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk1(i,j,k) = (Tracer%field(i,j,k,tau))**2
        enddo
     enddo
  enddo

  ! advection operator acting on squared tracer concentration   
  select case (Tracer%horz_advect_scheme)

      case (ADVECT_UPWIND)
          wrk2(isc:iec,jsc:jec,:) =  &
          -horz_advect_tracer_upwind(Adv_vel, wrk1)
      case (ADVECT_2ND_ORDER)
          wrk2(isc:iec,jsc:jec,:) =  &
          -horz_advect_tracer_2nd_order(Adv_vel, wrk1)
      case (ADVECT_4TH_ORDER)
          wrk2(isc:iec,jsc:jec,:) =  &
          -horz_advect_tracer_4th_order(Adv_vel, wrk1)
      case (ADVECT_6TH_ORDER)
          wrk2(isc:iec,jsc:jec,:) =  &
          -horz_advect_tracer_6th_order(Adv_vel, wrk1)
      case (ADVECT_QUICKER)
          wrk2(isc:iec,jsc:jec,:) =  &
          -horz_advect_tracer_quicker(Adv_vel, Tracer, wrk1, wrk1)
      case (ADVECT_QUICKMOM3)
          wrk2(isc:iec,jsc:jec,:) =  &
          -horz_advect_tracer_quickmom3(Adv_vel, Tracer, wrk1, wrk1)
      case (ADVECT_MDFL_SUP_B)
          wrk2(isc:iec,jsc:jec,:) =  &
          -advect_tracer_mdfl_sup_b(Time, Adv_vel, Thickness, wrk1, dtime)
      case (ADVECT_MDFL_SWEBY)
          wrk2(isc:iec,jsc:jec,:) =  &
          -advect_tracer_mdfl_sweby(Time, Adv_vel, Thickness, wrk1, dtime, sweby_limiter=1.0)
      case (ADVECT_DST_LINEAR)
          wrk2(isc:iec,jsc:jec,:) =  &
          -advect_tracer_mdfl_sweby(Time, Adv_vel, Thickness, wrk1, dtime, sweby_limiter=0.0)
      case (ADVECT_DST_LINEAR_TEST)
          wrk2(isc:iec,jsc:jec,:) =  &
          -advect_tracer_mdfl_sweby_test(Time, Adv_vel, Thickness, wrk1, dtime, sweby_limiter=0.0)
      case (ADVECT_MDPPM)
          wrk2(isc:iec,jsc:jec,:) =  &
          -advect_tracer_mdppm(Time, Adv_vel, Tracer, Thickness, wrk1, dtime)
      case (ADVECT_MDMDT_TEST)
          wrk2(isc:iec,jsc:jec,:) =  &
          -advect_tracer_mdmdt_test(Time, Adv_vel, Tracer, Thickness, wrk1, dtime)

      case (ADVECT_PSOM)
          use_psom=.true.  

  end select

      select case (Tracer%vert_advect_scheme)

      case (ADVECT_UPWIND)
          wrk3(isc:iec,jsc:jec,:) = &
          -vert_advect_tracer_upwind(Adv_vel, wrk1)
      case (ADVECT_2ND_ORDER)
          wrk3(isc:iec,jsc:jec,:) = &
          -vert_advect_tracer_2nd_order(Adv_vel, wrk1)
      case (ADVECT_4TH_ORDER)
          wrk3(isc:iec,jsc:jec,:) = &
          -vert_advect_tracer_4th_order(Adv_vel, wrk1)
      case (ADVECT_6TH_ORDER)
          wrk3(isc:iec,jsc:jec,:) = &
          -vert_advect_tracer_6th_order(Adv_vel, wrk1)
      case (ADVECT_QUICKER)
          wrk3(isc:iec,jsc:jec,:) = &
          -vert_advect_tracer_quicker(Adv_vel, Tracer, wrk1, wrk1)
      case (ADVECT_QUICKMOM3)
          wrk3(isc:iec,jsc:jec,:) = &
          -vert_advect_tracer_quickmom3(Adv_vel, Tracer, wrk1, wrk1)

      ! the mdfl, mdppm, dst_linear, and PSOM schemes are three-dimensional, 
      ! with 3-dimensional advection tendency computed in horz_advect_tracer
      case (ADVECT_MDFL_SUP_B)
      case (ADVECT_MDFL_SWEBY)
      case (ADVECT_MDFL_SWEBY_TEST)
      case (ADVECT_DST_LINEAR)
      case (ADVECT_DST_LINEAR_TEST)
      case (ADVECT_MDPPM)
      case (ADVECT_MDPPM_TEST)
      case (ADVECT_PSOM)

  end select

  ! PSOM is a special case, since it requires extra arrays to carry the 
  ! moments of tracer concentration. This then requires an extra tracer
  ! array as well, and these are defined in ocean_passive.F90, and they 
  ! require added settings in the field table for temp_sq and salt_sq. 
  ! The case for temp and salt are the only cases so far supported.
  ! If wish to use this diagnostic for other PSOM advected tracers, then need 
  ! new code in ocean_passive.F90, emulating that for temp_sq and salt_sq. 
  if(use_psom) then 
      wrk3(:,:,:) = 0.0
      if(ntracer==index_temp .and. index_temp_sq > 0) then 
          wrk2(isc:iec,jsc:jec,:) =                                            &
          -advect_tracer_psom(Time, Adv_vel, T_prog(index_temp_sq), Thickness, &
                              dtime, ntracer, do_passive_sq=.true.)
      endif 
      if(ntracer==index_salt .and. index_salt_sq > 0) then 
          wrk2(isc:iec,jsc:jec,:) =                                            &
          -advect_tracer_psom(Time, Adv_vel, T_prog(index_salt_sq), Thickness, &
                              dtime, ntracer, do_passive_sq=.true.)
      endif
  endif

  ! wrk1 is now the total tendency for the squared tracer concentration 
  wrk1(:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = wrk2(i,j,k) + wrk3(i,j,k)
        enddo
     enddo
  enddo

  ! compute the dissipation from advection truncation errors
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           term1 = advect_tendency(i,j,k) &
           *(2.0*Thickness%rho_dzt(i,j,k,tau)*Tracer%field(i,j,k,tau) + dtime*advect_tendency(i,j,k))
           term2 = -Thickness%rho_dzt(i,j,k,taup1)*wrk1(i,j,k)
           wrk4(i,j,k) =  -(Tracer%conversion**2)*dtimer*(term1+term2)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_tracer_adv_diss(ntracer), wrk4(:,:,:))

  if(id_tracer2_advection(ntracer) > 0) then 
     call diagnose_3d(Time, Grd, id_tracer2_advection(ntracer), wrk1(:,:,:)*Tracer%conversion**2)
  endif

  call mpp_clock_end(id_clock_adv_diss)

     
end subroutine compute_adv_diss
! </SUBROUTINE> NAME="compute_adv_diss"


!#######################################################################
! <SUBROUTINE NAME="get_tracer_stats">
!
! <DESCRIPTION>
! Compute the upper/lower values of a 3D field, returning values via arguments
! </DESCRIPTION>
subroutine get_tracer_stats(A,Mask,tmin,tmax)
  real, dimension(isc:iec,jsc:jec,nk), intent(in) :: A
  real, dimension(isc:iec,jsc:jec,nk), intent(in) :: Mask
  real, intent(out) :: tmin,tmax
  integer :: i,j,k
  tmin=1.e30; tmax=-1.e30
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           if (Mask(i,j,k)>0.) then
              tmin=min(tmin,A(i,j,k))
              tmax=max(tmax,A(i,j,k))
           endif
        enddo
     enddo
  enddo
end subroutine get_tracer_stats
! </SUBROUTINE> NAME="get_tracer_stats"


!#######################################################################
! <SUBROUTINE NAME="tracer_stats">
!
! <DESCRIPTION>
! Check the upper/lower values of a 3D field fall between speficied bounds
! reporting points that fall outside. Bring model down uncleanly if bounds are
! exceeded.
!
! NOTE: This is a debugging tool and not for normal use.
!
! </DESCRIPTION>
subroutine tracer_stats(A,tmin0,tmax0,label)
  real, dimension(isc:iec,jsc:jec,nk), intent(in) :: A
  real, intent(in) :: tmin0,tmax0
  character(len=*), intent(in) :: label
  integer :: i,j,k
  real :: tmin,tmax
  tmin=1.e30; tmax=-1.e30

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           if (tmask_mdppm(i,j,k)>0.) then
              tmin=min(tmin,A(i,j,k))
              tmax=max(tmax,A(i,j,k))
              if (j==1.and.A(i,j,k)<tmin0*0.99999999999999) write(0,*) 'under i,k=',i,k,A(i,j,k)
              if (j==1.and.A(i,j,k)>tmax0*1.00000000000001) write(0,*) 'over i,k=',i,k,A(i,j,k)
           endif
        enddo
     enddo
  enddo
  write(0,*) 't_stats: ',label,' min=',tmin,' max=',tmax
  if ((tmin<tmin0*0.99999999999999).or.(tmax>tmax0*1.00000000000001)) then
    write(0,*) 'tmin0=',tmin0,'tmax0=',tmax0
    write(0,*) 'tmin =',tmin, 'tmax =',tmax
!   stop 'Overshoots!'
  endif
  if ((tmin<tmin0*0.9999999999999).or.(tmax>tmax0*1.0000000000001)) then
    stop 'Overshoots!'
  endif
end subroutine
! </SUBROUTINE> NAME="tracer_stats"


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_advect_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_tracer_advect_restart(T_prog)
  type(ocean_prog_tracer_type), intent(in)           :: T_prog(:)

  if(ANY(T_prog(1:num_prog_tracers)%horz_advect_scheme == ADVECT_PSOM) )then
     call save_restart(Adv_restart)
  end if

end subroutine ocean_tracer_advect_restart
! </SUBROUTINE> NAME="ocean_tracer_advect_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_advect_end">
!
! <DESCRIPTION>
! Write the PSOM moments for restarts.  
! </DESCRIPTION>
subroutine ocean_tracer_advect_end(Time, T_prog)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  
  integer :: n
  character(len=128) :: filename

  integer :: stdoutunit 
  stdoutunit=stdout() 

  filename = 'RESTART/ocean_psom_moments.res'

  if(.not. write_a_restart) then
    write(stdoutunit,'(/a)') '==>Warning from ocean_tracer_advect_end: NO restart written for PSOM moments.'
    call mpp_error(WARNING,'==>Warning from ocean_tracer_advect_end: NO restart written for PSOM moments.')
    return
  endif 

  call ocean_tracer_advect_restart(T_prog)

  do n=1,num_prog_tracers
    if (T_prog(n)%horz_advect_scheme == ADVECT_PSOM) then
      call tracer_psom_chksum(Time, T_prog(n))
      write (stdoutunit,*) 'Completed write of psom restart for tracer ', trim(T_prog(n)%name)
    endif

  enddo
  
  return

end subroutine ocean_tracer_advect_end
! </SUBROUTINE> NAME="ocean_tracer_advect_end"



end module ocean_tracer_advect_mod

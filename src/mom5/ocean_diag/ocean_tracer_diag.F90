module ocean_tracer_diag_mod
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Routines for tracer diagnostics 
!</OVERVIEW>
!
!<DESCRIPTION>
! Routines for tracer diagnostics.  Some are printed to ascii output, some are sent 
! to diagnostic manager. 
! </DESCRIPTION>
!
!<NAMELIST NAME="ocean_tracer_diag_nml">
!  <DATA NAME="tracer_conserve_days" UNITS="days" TYPE="real">
!  Number of days between which compute the tracer conservation diagnostics. 
!  </DATA> 
!  <DATA NAME="diag_step" UNITS="dimensionless" TYPE="integer">
!  Number of time steps between which compute the diagnostics.
!  </DATA> 
!
!  <DATA NAME="debug_diagnose_mixingA" TYPE="logical">
!  Set true for help with debugging the diagnostic for mixing.
!  </DATA> 
!  <DATA NAME="debug_diagnose_mixingB" TYPE="logical">
!  Set true for more help with debugging the diagnostic for mixing.
!  Lots of output.
!  </DATA> 
!  <DATA NAME="debug_diagnose_mixingC" TYPE="logical">
!  Set true for more help with debugging the diagnostic for mixing.
!  Lots of output. 
!  </DATA> 
!  <DATA NAME="debug_diagnose_mixingD" TYPE="logical">
!  Set true for more help with debugging the diagnostic for mixing.
!  Lots of output. 
!  </DATA> 
!  <DATA NAME="smooth_kappa_sort" TYPE="integer">
!  Number of 1-2-1 smooths applied to kappa_sort 
!  </DATA> 
!  <DATA NAME="smooth_dzt_rho_sort" TYPE="integer">
!  Number of 1-2-1 smooths applied to rho_sort  
!  </DATA> 
!  <DATA NAME="rho_grad_min" UNITS="kg/m^3/m" TYPE="real">
!  min vertical density gradient (kg/m^3/m) used in computing kappa sorted
!  in the diagnostic mixing sorted. 
!  </DATA> 
!  <DATA NAME="rho_grad_max" UNITS="kg/m^3/m" TYPE="real">
!  max vertical density gradient (kg/m^3/m) used in computing kappa sorted
!  </DATA>
!  <DATA NAME="buoyancy_crit" UNITS="m/s^2" TYPE="real">
!  Critical buoyancy difference relative to surface for computing mixed
!  layer depth. Default buoyancy_crit=0.0003. 
!  </DATA>
!  <DATA NAME="dtheta_crit" UNITS="degC" TYPE="real">
!  Critical temperature difference relative to surface for computing
!  mixed_layer_depth_dtheta . Default dtheta_crit=2.0.
!  </DATA>
!  <DATA NAME="diagnose_mixing_days" UNITS="day" TYPE="real">
!  Days over which time average the thickness weighted density before taking its 
!  time tendency for use in computing the effective diapycnal diffusivity.  
!  </DATA>
!
!  <DATA NAME="psu2ppt" TYPE="real">
!  The preTEOS10 EOS used in MOM requires salinity to 
!  use the Practical Salinity Scale (pss).  This scale is 
!  also known as the Practical Salinity Unit (psu).
!
!  However, salinity as an absolute concentration in 
!  parts per thousand is more convenient to use when 
!  performing budget analyses such as in this module.  
!  Conversion between pss and ppt depends on the precise
!  ratio of ions in the seawater. Hence, the conversion
!  is not constant. However, it is close to a constant,
!  as reported in Jackett etal (2004).  For purposes of 
!  budgets, we take this conversion as a constant.  
!  The conversion is 
!
!  s(ppt) = psu2ppt * s(psu) 
!
!  where again s(psu) is what MOM carries as its 
!  prognostic salinity field when preTEOS10 EOS is used. 
!
!  Jackett etal (2004), correcting a type in equation (53) 
!  of Feistel (2003), report that 
!
!  s(ppt) = 1.004867 * s(psu)
!
!  </DATA>
!
!  <DATA NAME="smooth_mld" TYPE="integer">
!  Smooth the diagnosed mixed layer depth. Default smooth_mld=.false. 
!  </DATA> 
!
!  <DATA NAME="smooth_mld_for_subduction" TYPE="integer">
!  Smooth the diagnosed mixed layer depth to be used for subduction
!  diagnostics. Default smooth_mld_for_subduction=.true. 
!  </DATA> 
!
!  <DATA NAME="do_bitwise_exact_sum" TYPE="logical">
!  Set true to do bitwise exact global sum. When it is false, the global
!  sum will be non-bitwise_exact, but will significantly increase efficiency. 
!  The default value is false.
!  </DATA> 
!
!</NAMELIST>
!
use constants_mod,         only: epsln, c2dbars
use diag_manager_mod,      only: register_diag_field, send_data, need_data
use fms_mod,               only: open_namelist_file, check_nml_error, close_file
use fms_mod,               only: write_version_number, FATAL, WARNING, stdout, stdlog
use mpp_domains_mod,       only: mpp_global_field, mpp_global_sum, mpp_global_max, mpp_global_min
use mpp_domains_mod,       only: mpp_update_domains
use mpp_domains_mod,       only: BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
use mpp_mod,               only: input_nml_file, mpp_error, mpp_max
use mpp_mod,               only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use time_manager_mod,      only: time_type, increment_time

use ocean_density_mod,     only: density_delta_sfc
use ocean_domains_mod,     only: get_local_indices
use ocean_obc_mod,         only: ocean_obc_tracer_flux, ocean_obc_mass_flux
use ocean_operators_mod,   only: FAX, FAY, BAX, BAY, FMX, FMY, FDX_T, FDY_T, S2D
use ocean_parameters_mod,  only: DEPTH_BASED, TWO_LEVEL, THREE_LEVEL, missing_value
use ocean_parameters_mod,  only: rho0, rho0r, sec_in_yr, sec_in_yr_r, grav 
use ocean_parameters_mod,  only: onehalf, onefourth
use ocean_tracer_util_mod, only: tracer_min_max, dzt_min_max, sort_shell_array
use ocean_types_mod,       only: ocean_domain_type, ocean_grid_type, ocean_velocity_type
use ocean_types_mod,       only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,       only: ocean_external_mode_type, ocean_density_type, ocean_adv_vel_type
use ocean_types_mod,       only: ocean_time_type, ocean_time_steps_type, ocean_thickness_type
use ocean_types_mod,       only: ocean_lagrangian_type
use ocean_util_mod,        only: write_timestamp, diagnose_2d
use ocean_tracer_util_mod, only: diagnose_3d_rho
use ocean_workspace_mod,   only: wrk1, wrk2, wrk3, wrk4
use ocean_workspace_mod,   only: wrk1_2d, wrk2_2d, wrk3_2d, wrk4_2d
use ocean_workspace_mod,   only: wrk1_v, wrk2_v , wrk3_v 
use ocean_workspace_mod,   only: wrk1_v2d

implicit none

private

#include <ocean_memory.h>

real    :: dtts
real    :: dteta
real    :: dtime
real    :: dtimer
real    :: aidif 

integer :: num_prog_tracers
integer :: num_diag_tracers
logical :: have_obc 
logical :: use_blobs

integer :: index_temp
integer :: index_salt
integer :: index_frazil

! for area of domain 
real :: cellarea
real :: cellarea_r

! for vertical coordinate 
integer :: vert_coordinate_class=1

! for diagnostics clocks 
integer :: id_compute_subduction 
integer :: id_compute_tracer_mld 
integer :: id_mixed_layer_depth
integer :: id_mixed_layer_depth_dtheta
integer :: id_potrho_mixed_layer
integer :: id_diagnose_depth_of_potrho
integer :: id_diagnose_depth_of_theta
integer :: id_diagnose_tracer_on_rho
integer :: id_diagnose_tracer_zrho_on_rho
integer :: id_tracer_numerical
integer :: id_diagnose_mixing
integer :: id_conservation
integer :: id_total_tracer
integer :: id_total_mass
integer :: id_total_volume
integer :: id_tracer_integrals
integer :: id_tracer_change
integer :: id_tracer_land_cell_check

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

logical :: module_is_initialized = .FALSE.
integer :: tendency=0

character(len=128) :: version=&
     '$Id: ocean_tracer_diag.F90,v 20.0 2013/12/14 00:12:55 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

! for tracer conservation: need to know num_prog_tracers to allocate
real, dimension(:,:,:), allocatable :: tracer_source
real, dimension(:,:,:), allocatable :: tracer_runoff
real, dimension(:,:,:), allocatable :: tracer_calving
real, dimension(:,:,:), allocatable :: tracer_otf
real, dimension(:,:,:), allocatable :: tracer_stf
real, dimension(:,:,:), allocatable :: tracer_btf
real, dimension(:,:,:), allocatable :: tracer_pme
real, dimension(:,:),   allocatable :: tracer_frazil
real, dimension(:,:,:), allocatable :: tracer_tend
real, dimension(:,:,:), allocatable :: tracer_eta_smooth
real, dimension(:,:,:), allocatable :: tracer_pbot_smooth
real, dimension(:,:)  , allocatable :: pme_flux
real, dimension(:,:)  , allocatable :: runoff_flux
real, dimension(:,:)  , allocatable :: source_flux
real, dimension(:,:)  , allocatable :: calving_flux
real, dimension(:,:)  , allocatable :: eta_smooth
real, dimension(:,:)  , allocatable :: pbot_smooth
real, dimension(:,:)  , allocatable :: area_t
real, dimension(:), allocatable     :: tracer_start 
real, dimension(:), allocatable     :: tracer_final 
real                                :: mass_start=0.0
real                                :: mass_final=0.0
real                                :: volume_start=0.0
real                                :: volume_final=0.0
real                                :: tracer_conserve_days=30.0 ! number of days for tracer conservation
integer                             :: itts_tracer, itte_tracer  ! itt starting and ending tracer conservation 
integer                             :: itts_mass, itte_mass      ! itt starting and ending mass conservation

! for saving out the total tracer content and mass
integer, allocatable, dimension(:) :: id_prog_total
integer, allocatable, dimension(:) :: id_prog_surface_ave
integer, allocatable, dimension(:) :: id_prog_surface_area_ave
integer, allocatable, dimension(:) :: id_prog_global_ave
integer :: id_mass_seawater=-1
integer :: id_volume_seawater=-1
integer :: id_press_ave=-1

! for diagnosing mixing between temperature classes. 
logical :: debug_diagnose_mixingA=.false.
logical :: debug_diagnose_mixingB=.false.
logical :: debug_diagnose_mixingC=.false.
logical :: debug_diagnose_mixingD=.false.

real    :: rho_grad_max=1e28   !max vertical density gradient (kg/m^3/m) used in computing kappa
real    :: rho_grad_min=1e-5   !min vertical density gradient (kg/m^3/m) used in computing kappa 
real    :: alpha=1.0           !should=alpha_linear_eos, but choose 1.0 to increase precision 
real    :: thick_column(2)     !thickness of sorted column
real    :: mass_column(2)      !mass of sorted column

integer :: smooth_kappa_sort=0 !for smoothing the diagnosed diffusivity 
integer :: nsortpts            ! number points to be sorted = number wet tracer points

integer, dimension(:),allocatable :: nlayer      ! # wet points per layer when specify this as constant 
real, dimension(:),   allocatable :: wetarea     ! horizontal area as function of depth
real, dimension(:),   allocatable :: normalize   !1.0/(rho0*wetarea(k))
real, dimension(:),   allocatable :: deriv_lay   ! vertical density deriv at bottom of tracer cell
real, dimension(:),   allocatable :: flux_lay    ! dianeutral tracer flux centered at bottom of sorted tracer cell
real, dimension(:),   allocatable :: kappa_lay   ! diagnosed dianeutral diffusivity (m^2/sec) in layer 
real, dimension(:,:), allocatable :: rho_lay     ! sorted density (kg/m^3) averaged onto layers 
real, dimension(:,:), allocatable :: dzt_rho_lay ! dzt*rho (kg/m^2) of a layer 
real, dimension(:,:), allocatable :: dzt_lay     ! thickness (m) of a density layer 
real, dimension(:,:), allocatable :: dzw_lay     ! distance (m) between density centers 
real, dimension(:,:), allocatable :: mass_lay    ! mass (kg) of a layer
real, dimension(:,:), allocatable :: volume_lay  ! volume (m^3) of a layer


! for the spurious mixing diagnostics 
integer :: id_kappa_simple      =-1
integer :: id_kappa_sort        =-1
integer :: id_temp_sort         =-1

! frequency which compute certain ascii diagnostics 
integer :: diag_step = -1

! for diagnostic manager 
logical :: used

! for output
integer :: unit=6

! for mixed layer depth 
integer :: id_mld        =-1
integer :: id_mld_sq     =-1
integer :: id_mld_nrho   =-1

! for subduction diagnostics 
integer :: id_subduction           =-1
integer :: id_subduction_mld       =-1
integer :: id_subduction_dhdt      =-1
integer :: id_subduction_horz      =-1
integer :: id_subduction_vert      =-1
integer :: id_subduction_nrho      =-1
integer :: id_subduction_dhdt_nrho =-1
integer :: id_subduction_horz_nrho =-1
integer :: id_subduction_vert_nrho =-1
integer, dimension(:), allocatable :: id_subduction_tracer
integer, dimension(:), allocatable :: id_subduction_dhdt_tracer
integer, dimension(:), allocatable :: id_subduction_horz_tracer
integer, dimension(:), allocatable :: id_subduction_vert_tracer
integer, dimension(:), allocatable :: id_subduction_mld_zflux_diff
logical :: compute_subduction_diags = .false.


! for tracer integrated over mixed layer
integer, dimension(:), allocatable :: id_tracer_mld
logical :: compute_tracer_mld_diags = .false.


! for mixed layer depth based solely on depth where SST - temp(k) = dtheta
integer :: id_mld_dtheta =-1

! for depth of isopycnal and potential temperature surfaces 
integer :: id_depth_of_potrho=-1
integer :: id_depth_of_theta =-1

! for tracer_on_rho
integer, dimension(:), allocatable :: id_tracer_on_rho
integer, dimension(:), allocatable :: id_tracer_zrho_on_rho

integer, dimension(:), allocatable :: id_global_variance
integer, dimension(:), allocatable :: id_k_variance

! for potential density mixed layer 
integer :: id_potrho_mix_depth=-1
integer :: id_potrho_mix_base=-1
real    :: buoyancy_crit=0.0003 ! (m/s^2) critical buoyancy difference relative to surface 

! for potential temperature based mixed layer 
real :: dtheta_crit=2.0 ! (degC) critical temperature difference relative to surface 

! for salt budgets
real :: psu2ppt=1.004867   ! conversion from psu to ppt concentration 

! for bitwise exact global sums independent of PE number
logical :: do_bitwise_exact_sum = .false.
integer :: global_sum_flag

! for computation of frazil contributions to heat, we need to have 
! frazil_factor equal to that used in the nml from ocean_frazil_mod.
! If tendency=TWO_LEVEL and if use GFDL SIS ocean model, then frazil_factor=1.0
! If tendency=THREE_LEVEL and if use GFDL SIS ocean model, then frazil_factor=.50
real :: frazil_factor=1.0

! for smoothing the diagnosed mld 
logical :: smooth_mld=.false.

! for mld used in subduction diagnostic
logical :: smooth_mld_for_subduction=.true.

public ocean_tracer_diag_init
public ocean_tracer_diagnostics 
public calc_mixed_layer_depth 
public calc_potrho_mixed_layer 
public send_tracer_variance
public diagnose_eta_tend_3dflux

private compute_subduction 
private compute_tracer_mld
private mixed_layer_depth 
private mixed_layer_depth_dtheta 
private potrho_mixed_layer 
private tracer_change
private tracer_integrals 
private tracer_land_cell_check
private tracer_conservation
private mass_conservation
private diagnose_kappa_simple
private diagnose_kappa_sort
private diagnose_depth_of_potrho
private diagnose_depth_of_theta
private diagnose_tracer_on_rho
private diagnose_tracer_zrho_on_rho
private total_tracer
private klevel_total_tracer
private total_mass 
private total_volume
private klevel_total_mass 
private send_total_tracer
private send_global_ave_tracer
private send_global_ave_pressure
private send_surface_ave_tracer
private send_surface_area_ave_tracer
private send_total_mass
private send_total_volume

namelist /ocean_tracer_diag_nml/ tracer_conserve_days , diag_step, psu2ppt,          &
                                 debug_diagnose_mixingA, debug_diagnose_mixingB,     &
                                 debug_diagnose_mixingC, debug_diagnose_mixingD,     &
                                 smooth_kappa_sort, rho_grad_min, rho_grad_max,      &
                                 buoyancy_crit, do_bitwise_exact_sum, frazil_factor, &
                                 smooth_mld, smooth_mld_for_subduction, dtheta_crit

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_diag_init">
!
! <DESCRIPTION>
! Initialize the ocean_tracer_diag module containing subroutines
! diagnosing tracer related properties of the simulation.  These are 
! not terms in the equations, but rather they are diagnosed from 
! terms. 
! </DESCRIPTION>
!
  subroutine ocean_tracer_diag_init(Grid, Domain, Time, Time_steps, Thickness, T_prog, T_diag, Dens, &
                                    ver_coordinate_class, blobs, obc)

    type(ocean_grid_type),        target, intent(in)    :: Grid
    type(ocean_domain_type),      target, intent(in)    :: Domain
    type(ocean_time_type),                intent(in)    :: Time
    type(ocean_time_steps_type),          intent(in)    :: Time_steps 
    type(ocean_thickness_type),           intent(in)    :: Thickness
    type(ocean_prog_tracer_type),         intent(in)    :: T_prog(:)
    type(ocean_diag_tracer_type),         intent(in)    :: T_diag(:)
    type(ocean_density_type),             intent(inout) :: Dens
    integer,                              intent(in)    :: ver_coordinate_class
    logical,                              intent(in)    :: blobs
    logical,                              intent(in)    :: obc

    integer :: n, ioun, io_status, ierr, tau
    integer :: stdoutunit,stdlogunit 
    stdoutunit=stdout();stdlogunit=stdlog() 

    if (module_is_initialized) return

    module_is_initialized = .TRUE.

    num_prog_tracers = size(T_prog,1)
    num_diag_tracers = size(T_diag,1)

    call write_version_number(version, tagname)

#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_tracer_diag_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_tracer_diag_nml')
#else
    ioun = open_namelist_file()
    read(ioun, ocean_tracer_diag_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_tracer_diag_nml')
    call close_file(ioun)
#endif
    write (stdlogunit, ocean_tracer_diag_nml)
    write (stdoutunit,'(/)')
    write (stdoutunit, ocean_tracer_diag_nml)

    if (diag_step == 0) diag_step = 1

    if(do_bitwise_exact_sum) then
       global_sum_flag = BITWISE_EXACT_SUM
    else
       global_sum_flag = NON_BITWISE_EXACT_SUM
    endif

    vert_coordinate_class = ver_coordinate_class 

#ifndef MOM_STATIC_ARRAYS
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    ni = Grid%ni
    nj = Grid%nj
    nk = Grid%nk
#endif

    Dom => Domain
    Grd => Grid

    dtts      = Time_steps%dtts
    dteta     = Time_steps%dteta
    dtime     = Time_steps%dtime_t
    tendency  = Time_steps%tendency
    aidif     = Time_steps%aidif 
    dtimer    = 1.0/(dtime+epsln)
    have_obc  = obc
    use_blobs = blobs
    tau       = Time%tau

    do n=1,num_prog_tracers
       if(T_prog(n)%name == 'temp' ) index_temp = n
       if(T_prog(n)%name == 'salt' ) index_salt = n
    enddo

    index_frazil=-1
    do n=1,num_diag_tracers
       if(T_diag(n)%name == 'frazil' ) index_frazil = n
    enddo

    cellarea   = Grd%tcellsurf
    cellarea_r = 1.0/(epsln + Grd%tcellsurf)

    allocate(area_t(isd:ied,jsd:jed) )    
    area_t(:,:) = Grd%dat(:,:)*Grd%tmask(:,:,1)
    if(have_obc) area_t(:,:) = area_t(:,:)*Grd%obc_tmask(:,:)

    if(tendency==THREE_LEVEL) then
        write (stdoutunit,'(/a)') &
        'Note: tracer and mass/volume conservation tests based on time_tendency==threelevel.'
    elseif(tendency==TWO_LEVEL) then
        write (stdoutunit,'(/a)') &
        'Note: tracer and mass/volume conservation tests based on time_tendency==twolevel.'
    endif

    write (stdoutunit,'(/a,f5.2,a)') &
    ' Note: Set frazil_factor = ',frazil_factor,' for computation of heat diagnostics.'
    write (stdoutunit,'(a)') &
     ' Be sure this agrees with the value set in nml for ocean_frazil_mod'


    ! registers to diagnostic manager 

    id_mass_seawater = register_diag_field ('ocean_model', 'total_mass_seawater', Time%model_time,&
                      'total mass of liquid seawater', 'kg', missing_value=missing_value,         &
                      range=(/-1e1,1e25/), standard_name='sea_water_mass')
    id_volume_seawater = register_diag_field ('ocean_model', 'total_volume_seawater', Time%model_time,&
                      'total volume of liquid seawater', 'm^3', missing_value=missing_value,          &
                      range=(/-1e1,1e25/), standard_name='sea_water_volume')
    
    allocate( id_prog_total(num_prog_tracers) )
    allocate( id_prog_surface_ave(num_prog_tracers) )
    allocate( id_prog_surface_area_ave(num_prog_tracers) )
    allocate( id_prog_global_ave(num_prog_tracers) )
    id_prog_total           = -1
    id_prog_surface_ave     = -1
    id_prog_surface_area_ave= -1
    id_prog_global_ave      = -1
    allocate(id_tracer_on_rho(num_prog_tracers))
    allocate(id_tracer_zrho_on_rho(num_prog_tracers))
    id_tracer_on_rho=-1
    allocate(id_global_variance(num_prog_tracers))
    allocate(id_k_variance(num_prog_tracers))
    id_global_variance(:) = -1
    id_k_variance(:)      = -1

    allocate( id_subduction_tracer(num_prog_tracers) )
    allocate( id_subduction_dhdt_tracer(num_prog_tracers) )
    allocate( id_subduction_horz_tracer(num_prog_tracers) )
    allocate( id_subduction_vert_tracer(num_prog_tracers) )
    allocate( id_subduction_mld_zflux_diff(num_prog_tracers) )
    allocate( id_tracer_mld(num_prog_tracers) )
    id_subduction_tracer(:)         = -1
    id_subduction_dhdt_tracer(:)    = -1
    id_subduction_horz_tracer(:)    = -1
    id_subduction_vert_tracer(:)    = -1
    id_subduction_mld_zflux_diff(:) = -1
    id_tracer_mld(:)                = -1 
    
    do n=1,num_prog_tracers

       if(T_prog(n)%name == 'temp') then 
           id_prog_total(n) = register_diag_field ('ocean_model','total_ocean_heat',                &
                 Time%model_time, 'Total heat in the liquid ocean referenced to 0degC','Joule/1e25',&
                 missing_value=missing_value, range=(/0.0,1e20/))
           id_tracer_on_rho(n) = register_diag_field ('ocean_model', 'temp_on_rho', &
                     Dens%potrho_axes(1:3), Time%model_time,                        &
                     'temperature on potential density surface', 'C',               & 
                     missing_value=missing_value, range=(/0.0,1.e3/))   
           id_tracer_zrho_on_rho(n) = register_diag_field ('ocean_model', 'temp_zrho_on_rho', &
                     Dens%potrho_axes(1:3), Time%model_time,                                  &
                     'temperature*dz/drho on potential density surface', 'degC * m/(kg/m^3)', & 
                     missing_value=missing_value, range=(/-1.e9,1.e9/))   
           id_prog_global_ave(n) = register_diag_field ('ocean_model',trim(T_prog(n)%name)//'_global_ave',&
                 Time%model_time, 'Global mean '//trim(T_prog(n)%name)//' in liquid seawater'             &
                 ,trim(T_prog(n)%units), missing_value=missing_value, range=(/-10.0,1000.0/),             &
                 standard_name='sea_water_potential_temperature')
           id_global_variance(n) = register_diag_field ('ocean_model',trim(T_prog(n)%name)//'_global_variance',&
                 Time%model_time, 'Global '//trim(T_prog(n)%name)//' variance',                                &
                 trim(T_prog(n)%units), missing_value=missing_value, range=(/-1e20,1e20/),                     &
                 standard_name='global_temperature_variance')
           id_k_variance(n) = register_diag_field ('ocean_model',trim(T_prog(n)%name)//'_variance', &
                 Grd%tracer_axes(3:3), Time%model_time, trim(T_prog(n)%name)//' variance',          &
                 trim(T_prog(n)%units), missing_value=missing_value, range=(/-1e20,1e20/),          &
                 standard_name='temperature_variance')

       elseif(T_prog(n)%name =='salt') then 
           id_prog_total(n) = register_diag_field ('ocean_model', 'total_ocean_salt', &
                           Time%model_time, 'total mass of salt in liquid seawater',  &
                           'kg/1e18', missing_value=missing_value, range=(/-1e2,1e10/))
           id_prog_global_ave(n) = register_diag_field ('ocean_model',trim(T_prog(n)%name)//'_global_ave',&
                 Time%model_time, 'Global mean '//trim(T_prog(n)%name)//' in liquid seawater'             &
                 ,trim(T_prog(n)%units), missing_value=missing_value, range=(/-10.0,1000.0/),             &
                 standard_name='sea_water_salinity')
           id_tracer_on_rho(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_on_rho', &
                Dens%potrho_axes(1:3), Time%model_time,                                               &
                trim(T_prog(n)%name)//' on potential density surface', trim(T_prog(n)%units),         &  
                missing_value=missing_value, range=(/0.0,1.e3/))  
           id_tracer_zrho_on_rho(n) = register_diag_field ('ocean_model',       &
                trim(T_prog(n)%name)//'_zrho_on_rho',                           &
                Dens%potrho_axes(1:3), Time%model_time,                         &
                trim(T_prog(n)%name)//' *dz/drho on potential density surface', &
                trim(T_prog(n)%units)//'*m/(kg/m^3)',                           & 
                missing_value=missing_value, range=(/-1.e9,1.e9/))  
           id_global_variance(n) = register_diag_field ('ocean_model',trim(T_prog(n)%name)//'_global_variance',&
                 Time%model_time, 'Global '//trim(T_prog(n)%name)//' variance',                                &
                 trim(T_prog(n)%units), missing_value=missing_value, range=(/-1e20,1e20/),                     &
                 standard_name='global_salinity_variance')
           id_k_variance(n) = register_diag_field ('ocean_model',trim(T_prog(n)%name)//'_variance', &
                 Grd%tracer_axes(3:3), Time%model_time, trim(T_prog(n)%name)//' variance',          &
                 trim(T_prog(n)%units), missing_value=missing_value, range=(/-1e20,1e20/),          &
                 standard_name='salt_variance')

       elseif(T_prog(n)%name(1:3) =='age') then 
           id_prog_total(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_total', &
                           Time%model_time, 'mass integrated '//trim(T_prog(n)%name),             &
                           'yr', missing_value=missing_value, range=(/-1e2,1e10/))
           id_prog_global_ave(n) = register_diag_field ('ocean_model',trim(T_prog(n)%name)//'_global_ave',&
                 Time%model_time, 'Global mean '//trim(T_prog(n)%name)//' in liquid seawater'             &
                 ,trim(T_prog(n)%units), missing_value=missing_value, range=(/-10.0,1000.0/))
           id_tracer_on_rho(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_on_rho', &
                Dens%potrho_axes(1:3), Time%model_time,                                               &
                trim(T_prog(n)%name)//' on potential density surface', trim(T_prog(n)%units),         &  
                missing_value=missing_value, range=(/0.0,1.e3/))  
           id_tracer_zrho_on_rho(n) = register_diag_field ('ocean_model',       &
                trim(T_prog(n)%name)//'_zrho_on_rho',                           &
                Dens%potrho_axes(1:3), Time%model_time,                         &
                trim(T_prog(n)%name)//' *dz/drho on potential density surface', &
                trim(T_prog(n)%units)//'*m/(kg/m^3)',                           & 
                missing_value=missing_value, range=(/-1.e9,1.e9/))  
           id_global_variance(n) = register_diag_field ('ocean_model',trim(T_prog(n)%name)//'_global_variance',&
                 Time%model_time, 'Global '//trim(T_prog(n)%name)//' variance',                                &
                 trim(T_prog(n)%units), missing_value=missing_value, range=(/-1e20,1e20/))
           id_k_variance(n) = register_diag_field ('ocean_model',trim(T_prog(n)%name)//'_variance', &
                 Grd%tracer_axes(3:3), Time%model_time, trim(T_prog(n)%name)//' variance',          &
                 trim(T_prog(n)%units), missing_value=missing_value, range=(/-1e20,1e20/))

       else 
           id_prog_total(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_total', &
                           Time%model_time, 'mass integrated '//trim(T_prog(n)%name),             &
                           'kg/1e18', missing_value=missing_value, range=(/-1e6,1e6/))
           id_prog_global_ave(n) = register_diag_field ('ocean_model',trim(T_prog(n)%name)//'_global_ave', &
            Time%model_time, 'Global mass weighted mean '//trim(T_prog(n)%name)//' in liquid seawater'     &
            ,trim(T_prog(n)%units), missing_value=missing_value, range=(/-10.0,1000.0/))
           id_tracer_on_rho(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_on_rho', &
                Dens%potrho_axes(1:3), Time%model_time,                                               &
                trim(T_prog(n)%name)//' on potential density surface', trim(T_prog(n)%units),         &  
                missing_value=missing_value, range=(/0.0,1.e3/))  
           id_tracer_zrho_on_rho(n) = register_diag_field ('ocean_model',       &
                trim(T_prog(n)%name)//'_zrho_on_rho',                           &
                Dens%potrho_axes(1:3), Time%model_time,                         &
                trim(T_prog(n)%name)//' *dz/drho on potential density surface', &
                trim(T_prog(n)%units)//'*m/(kg/m^3)',                           & 
                missing_value=missing_value, range=(/-1.e9,1.e9/))  
           id_global_variance(n) = register_diag_field ('ocean_model',trim(T_prog(n)%name)//'_global_variance',&
                 Time%model_time, 'Global '//trim(T_prog(n)%name)//' variance',                                &
                 trim(T_prog(n)%units), missing_value=missing_value, range=(/-1e20,1e20/))
           id_k_variance(n) = register_diag_field ('ocean_model',trim(T_prog(n)%name)//'_variance', &
                 Grd%tracer_axes(3:3), Time%model_time, trim(T_prog(n)%name)//' variance',          &
                 trim(T_prog(n)%units), missing_value=missing_value, range=(/-1e20,1e20/))
       endif 

       id_prog_surface_ave(n) = register_diag_field ('ocean_model',trim(T_prog(n)%name)//'_surface_ave',  &
        Time%model_time, 'Global mass weighted mean surface '//trim(T_prog(n)%name)//' in liquid seawater'&
       ,trim(T_prog(n)%units), missing_value=missing_value, range=(/-10.0,1000.0/))

       id_prog_surface_area_ave(n) = register_diag_field ('ocean_model',trim(T_prog(n)%name)//'_surface_area_ave',&
        Time%model_time, 'Global area weighted mean surface '//trim(T_prog(n)%name)//' in liquid seawater'        &
       ,trim(T_prog(n)%units), missing_value=missing_value, range=(/-10.0,1000.0/))


       id_subduction_tracer(n) =                                                 &
        register_diag_field ('ocean_model','subduction_'//trim(T_prog(n)%name),  &
        Grd%tracer_axes(1:2), Time%model_time,                                   &
        'Tracer transport across mld from subduction for '//trim(T_prog(n)%name),&  
        'kg/sec * '//trim(T_prog(n)%units), missing_value=missing_value, range=(/-1.e20,1.e20/))
       if(id_subduction_tracer(n) > 0) compute_subduction_diags = .true.

       id_subduction_dhdt_tracer(n) =                                               &
        register_diag_field ('ocean_model','subduction_dhdt_'//trim(T_prog(n)%name),&
        Grd%tracer_axes(1:2), Time%model_time,                                      &
        'Tracer transport across mld from dhdt for '//trim(T_prog(n)%name),         &
        'kg/sec * '//trim(T_prog(n)%units), missing_value=missing_value, range=(/-1.e20,1.e20/))
       if(id_subduction_dhdt_tracer(n) > 0) compute_subduction_diags = .true.

       id_subduction_horz_tracer(n) =                                                     &
        register_diag_field ('ocean_model','subduction_horz_'//trim(T_prog(n)%name),      &
        Grd%tracer_axes(1:2),Time%model_time,                                             &
        'Tracer transport across mld from from horz advection for '//trim(T_prog(n)%name),&
        'kg/sec * '//trim(T_prog(n)%units), missing_value=missing_value, range=(/-1.e20,1.e20/))
       if(id_subduction_horz_tracer(n) > 0) compute_subduction_diags = .true.

       id_subduction_vert_tracer(n) =                                                    &
        register_diag_field ('ocean_model','subduction_vert_'//trim(T_prog(n)%name),     &
        Grd%tracer_axes(1:2), Time%model_time,                                           &
        'Tracer transport across mld from vertical advection for '//trim(T_prog(n)%name),&
        'kg/sec * '//trim(T_prog(n)%units), missing_value=missing_value, range=(/-1.e20,1.e20/))
       if(id_subduction_vert_tracer(n) > 0) compute_subduction_diags = .true.

       id_subduction_mld_zflux_diff(n) =                                                     &
        register_diag_field ('ocean_model','subduction_mld_zflux_diff'//trim(T_prog(n)%name),&
        Grd%tracer_axes(1:2), Time%model_time,                                               &
        'Tracer transport from vertical diffusion at base of mld for '//trim(T_prog(n)%name),&
        'kg/sec * '//trim(T_prog(n)%units), missing_value=missing_value, range=(/-1.e20,1.e20/))
       if(id_subduction_mld_zflux_diff(n) > 0) compute_subduction_diags = .true.

       id_tracer_mld(n) =                                                                              &
        register_diag_field ('ocean_model',trim(T_prog(n)%name)//'_mld',                               &
        Grd%tracer_axes(1:2), Time%model_time,                                                         &
        'Vertically integrated tracer [sum(rho_dzt*tracer)] in mixed layer for '//trim(T_prog(n)%name),&
        'kg/m^2 * '//trim(T_prog(n)%units), missing_value=missing_value, range=(/-1.e20,1.e20/))
       if(id_tracer_mld(n) > 0) compute_tracer_mld_diags = .true.


    enddo

    id_kappa_sort = register_diag_field ('ocean_model','kappa_sort', &
                    Grd%tracer_axes(3:3), Time%model_time,           &
                    'kappa from sorting', 'm^2/sec',                 &
                    missing_value=missing_value, range=(/-1e6,1e6/))

    id_temp_sort = register_diag_field ('ocean_model','temp_sort', &
                   Grd%tracer_axes(3:3), Time%model_time,          &
                   'sorted temp', 'C',                             &
                   missing_value=missing_value, range=(/-1e6,1e6/))

    id_kappa_simple = register_diag_field ('ocean_model','kappa_simple', &
                      Grd%tracer_axes(3:3), Time%model_time,             &
                      'kappa from horz avg', 'm^2/sec',                  &
                      missing_value=missing_value, range=(/-1e6,1e6/))

    id_mld_dtheta = register_diag_field ('ocean_model', 'mld_dtheta',     &
             Grd%tracer_axes(1:2), Time%model_time,                       &
             'mixed layer depth determined by temperature criteria ', 'm',&
             missing_value=missing_value, range=(/0.0,1.e6/))   

    id_mld = register_diag_field ('ocean_model', 'mld',               &
             Grd%tracer_axes(1:2), Time%model_time,                   &
             'mixed layer depth determined by density criteria ', 'm',&
             missing_value=missing_value, range=(/0.0,1.e6/),         &
             standard_name='ocean_mixed_layer_thickness_defined_by_sigma_t')   

    id_mld_sq = register_diag_field ('ocean_model', 'mld_sq',                  &
             Grd%tracer_axes(1:2), Time%model_time,                            &
             'squared mixed layer depth determined by density criteria', 'm^2',&
             missing_value=missing_value, range=(/0.0,1.e12/),                 &
             standard_name='square_of_ocean_mixed_layer_thickness_defined_by_sigma_t')      

    id_mld_nrho = register_diag_field ('ocean_model', 'mld_nrho',      &
             Grd%tracer_axes(1:2), Time%model_time,                    &
             'neutral density at base of mixed layer depth',  'kg/m^3',&
             missing_value=missing_value, range=(/0.0,1.e6/))   

    id_subduction = register_diag_field ('ocean_model', 'subduction', &
             Grd%tracer_axes(1:2), Time%model_time,                   &
             'rate of mass transferred below the mixed layer base',   &
             'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
    if(id_subduction > 0) compute_subduction_diags = .true.

    id_subduction_mld = register_diag_field ('ocean_model', 'subduction_mld',&
             Grd%tracer_axes(1:2), Time%model_time,                          &
             'mixed layer depth used for subduction diagnostics',            &
             'm', missing_value=missing_value, range=(/-1.e2,1.e20/))
    if(id_subduction_mld > 0) compute_subduction_diags = .true.

    id_subduction_dhdt = register_diag_field ('ocean_model', 'subduction_dhdt',  &
             Grd%tracer_axes(1:2), Time%model_time,                              &
             'rate of mass transferred below the mixed layer base due to dh/dt', &
             'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
    if(id_subduction_dhdt > 0) compute_subduction_diags = .true.

    id_subduction_horz = register_diag_field ('ocean_model', 'subduction_horz',       &
             Grd%tracer_axes(1:2), Time%model_time,                                   &
             'rate of mass transferred below the mixed layer base due to horz advect',&
             'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
    if(id_subduction_horz > 0) compute_subduction_diags = .true.

    id_subduction_vert = register_diag_field ('ocean_model', 'subduction_vert',       &
             Grd%tracer_axes(1:2), Time%model_time,                                   &
             'rate of mass transferred below the mixed layer base due to vert advect',&
             'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
    if(id_subduction_vert > 0) compute_subduction_diags = .true.

    id_subduction_nrho = register_diag_field ('ocean_model', 'subduction_nrho',                     & 
             Dens%neutralrho_axes(1:3), Time%model_time,                                            &
             'rate of mass transferred below the mixed layer base as binned to neutral rho classes',&
             'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
    if(id_subduction_nrho > 0) compute_subduction_diags = .true.

    id_subduction_dhdt_nrho = register_diag_field ('ocean_model', 'subduction_dhdt_nrho',                       & 
             Dens%neutralrho_axes(1:3), Time%model_time,                                                        &
             'rate of mass transferred below the mixed layer base due to dhdt as binned to neutral rho classes',&
             'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
    if(id_subduction_dhdt_nrho > 0) compute_subduction_diags = .true.

    id_subduction_horz_nrho = register_diag_field ('ocean_model', 'subduction_horz_nrho',                                & 
             Dens%neutralrho_axes(1:3), Time%model_time,                                                                 &
             'rate of mass transferred below the mixed layer base due to horz velocity as binned to neutral rho classes',&
             'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
    if(id_subduction_horz_nrho > 0) compute_subduction_diags = .true.

    id_subduction_vert_nrho = register_diag_field ('ocean_model', 'subduction_vert_nrho',                                & 
             Dens%neutralrho_axes(1:3), Time%model_time,                                                                 &
             'rate of mass transferred below the mixed layer base due to vert velocity as binned to neutral rho classes',&
             'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
    if(id_subduction_vert_nrho > 0) compute_subduction_diags = .true.

    id_depth_of_potrho = register_diag_field ('ocean_model', 'depth_of_potrho', &
                         Dens%potrho_axes(1:3), Time%model_time,                &
                         'depth of potential density surface', 'm',             &
                         missing_value=missing_value, range=(/0.0,1.e10/))   

    id_depth_of_theta = register_diag_field ('ocean_model', 'depth_of_theta', &
                        Dens%theta_axes(1:3), Time%model_time,                &
                        'depth of potential temp surface', 'm',               &
                        missing_value=missing_value, range=(/0.0,1.e10/))   

    id_potrho_mix_depth = register_diag_field ('ocean_model','potrho_mix_depth',  &
                          Grd%tracer_axes(1:2), Time%model_time,                  &
                          'Depth of potential density mixed layer','m',           &
                          missing_value=missing_value, range=(/-1e6,1e6/))

    id_potrho_mix_base = register_diag_field ('ocean_model','potrho_mix_base',  &
                         Grd%tracer_axes(1:2), Time%model_time,                 &
                         'Potential density at mixed layer base','kg/m^3',      &
                         missing_value=missing_value, range=(/-1e6,1e6/))

    id_press_ave    = register_diag_field('ocean_model','press_ave',      &
                   Time%model_time, 'global mean absolute ocean pressure, including surface loading',&
                   'dbar', missing_value=missing_value, range=(/-1e1,1e20/))


    ! define clock integers 
    id_compute_subduction          = mpp_clock_id('(Ocean tracer_diag: compute_subduction)',grain=CLOCK_ROUTINE)
    id_compute_tracer_mld          = mpp_clock_id('(Ocean tracer_diag: compute_tracer_mld)',grain=CLOCK_ROUTINE)
    id_mixed_layer_depth_dtheta    = mpp_clock_id('(Ocean tracer_diag: mld_dtheta)'        ,grain=CLOCK_ROUTINE)
    id_mixed_layer_depth           = mpp_clock_id('(Ocean tracer_diag: mld)'               ,grain=CLOCK_ROUTINE)
    id_potrho_mixed_layer          = mpp_clock_id('(Ocean tracer_diag: rho_mld)'           ,grain=CLOCK_ROUTINE)
    id_diagnose_depth_of_potrho    = mpp_clock_id('(Ocean tracer_diag: potrho depth)'      ,grain=CLOCK_ROUTINE)
    id_diagnose_depth_of_theta     = mpp_clock_id('(Ocean tracer_diag: theta depth)'       ,grain=CLOCK_ROUTINE)
    id_diagnose_tracer_on_rho      = mpp_clock_id('(Ocean tracer_diag: tracer on rho)'     ,grain=CLOCK_ROUTINE)
    id_diagnose_tracer_zrho_on_rho = mpp_clock_id('(Ocean tracer_diag: tracer_zrho on rho)',grain=CLOCK_ROUTINE)
    id_tracer_numerical            = mpp_clock_id('(Ocean tracer_diag: numerical)'         ,grain=CLOCK_ROUTINE)
    id_diagnose_mixing             = mpp_clock_id('(Ocean tracer_diag: diag mixing'        ,grain=CLOCK_ROUTINE)
    id_conservation                = mpp_clock_id('(Ocean tracer_diag: conserve)'          ,grain=CLOCK_ROUTINE)
    id_total_mass                  = mpp_clock_id('(Ocean tracer_diag: total water mass)'  ,grain=CLOCK_ROUTINE)
    id_total_volume                = mpp_clock_id('(Ocean tracer_diag: total water volume)',grain=CLOCK_ROUTINE)
    id_total_tracer                = mpp_clock_id('(Ocean tracer_diag: total tracer)'      ,grain=CLOCK_ROUTINE)
    id_tracer_integrals            = mpp_clock_id('(Ocean tracer_diag: integrals)'         ,grain=CLOCK_ROUTINE)
    id_tracer_change               = mpp_clock_id('(Ocean tracer_diag: change)'            ,grain=CLOCK_ROUTINE)
    id_tracer_land_cell_check      = mpp_clock_id('(Ocean tracer_diag: land check)'        ,grain=CLOCK_ROUTINE)


    ! compute tau value for Dens%mld_subduction 
    call calc_mixed_layer_depth(Thickness,                                       &
                                T_prog(index_salt)%field(isd:ied,jsd:jed,:,tau), &
                                T_prog(index_temp)%field(isd:ied,jsd:jed,:,tau), &
                                Dens%rho(isd:ied,jsd:jed,:,tau),                 &
                                Dens%pressure_at_depth(isd:ied,jsd:jed,:),       &
                                Dens%mld_subduction(:,:), smooth_mld_for_subduction) 


    end subroutine ocean_tracer_diag_init
! </SUBROUTINE>  NAME="ocean_tracer_diag_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_diagnostics">
!
! <DESCRIPTION>
! Call diagnostics related to the tracer fields. 
! </DESCRIPTION>

subroutine ocean_tracer_diagnostics(Time, Thickness, T_prog, T_diag, Dens, &
                                    Ext_mode, Velocity, Adv_vel, &
                                    diff_cbt, pme, melt, runoff, calving)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_prog_tracer_type),   intent(in)    :: T_prog(:)
  type(ocean_diag_tracer_type),   intent(in)    :: T_diag(:)
  type(ocean_density_type),       intent(inout) :: Dens
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_velocity_type),      intent(in)    :: Velocity
  type(ocean_adv_vel_type),       intent(in)    :: Adv_vel
  real, dimension(isd:,jsd:,:,:), intent(in)    :: diff_cbt
  real, dimension(isd:,jsd:),     intent(in)    :: pme
  real, dimension(isd:,jsd:),     intent(in)    :: melt
  real, dimension(isd:,jsd:),     intent(in)    :: runoff
  real, dimension(isd:,jsd:),     intent(in)    :: calving 

  type(time_type) :: next_time 

  integer  :: tau, taum1, taup1
  integer  :: n

  if(size(T_prog,1) /= num_prog_tracers) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_tracer_diag_mod (ocean_tracer_diagnostics): T_prog has wrong dimensions')
  endif 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

  next_time = increment_time(Time%model_time, int(dtts), 0)
  
  ! compute mixed layer depth and send to diagnostics manager

  call mpp_clock_begin(id_mixed_layer_depth)
  if(need_data(id_mld, next_time) .or. need_data(id_mld_sq, next_time) &
                                  .or. need_data(id_mld_nrho, next_time)) then 
    call mixed_layer_depth(Thickness,                                                          &
                           T_prog(index_salt)%field(isd:ied,jsd:jed,:,tau),                    &
                           T_prog(index_temp)%field(isd:ied,jsd:jed,:,tau),                    &
                           Dens%rho(isd:ied,jsd:jed,:,tau),Dens%neutralrho(isd:ied,jsd:jed,:), &
                           Dens%pressure_at_depth(isd:ied,jsd:jed,:), Time)
  endif
  call mpp_clock_end(id_mixed_layer_depth)

  call mpp_clock_begin(id_mixed_layer_depth_dtheta)
  if(need_data(id_mld_dtheta, next_time)) then 
    call mixed_layer_depth_dtheta(Thickness,                               &
                           T_prog(index_temp)%field(isd:ied,jsd:jed,:,tau),&
                           Time)
  endif
  call mpp_clock_end(id_mixed_layer_depth_dtheta)

  call mpp_clock_begin(id_potrho_mixed_layer)
  call potrho_mixed_layer(Time, Thickness, Dens)
  call mpp_clock_end(id_potrho_mixed_layer)

  call mpp_clock_begin(id_compute_subduction)
  if(     need_data(id_subduction,      next_time)      .or. need_data(id_subduction_dhdt, next_time)      &
     .or. need_data(id_subduction_horz, next_time)      .or. need_data(id_subduction_vert, next_time)      &
     .or. need_data(id_subduction_nrho, next_time)      .or. need_data(id_subduction_dhdt_nrho, next_time) &
     .or. need_data(id_subduction_horz_nrho, next_time) .or. need_data(id_subduction_vert_nrho, next_time) &
     .or. need_data(id_subduction_mld, next_time)       .or. compute_subduction_diags) then  
    call compute_subduction(Time, Thickness, Velocity, Adv_vel, Dens, T_prog, &
                            T_prog(index_salt)%field(isd:ied,jsd:jed,:,taup1),&
                            T_prog(index_temp)%field(isd:ied,jsd:jed,:,taup1),&
                            diff_cbt)
  endif
  call mpp_clock_end(id_compute_subduction)

  call mpp_clock_begin(id_compute_tracer_mld)
  call compute_tracer_mld(Time, Thickness, Dens, T_prog,                     &
                          T_prog(index_salt)%field(isd:ied,jsd:jed,:,taup1), &
                          T_prog(index_temp)%field(isd:ied,jsd:jed,:,taup1))
  call mpp_clock_end(id_compute_tracer_mld)

  ! compute tracer as function of potential density 
  call mpp_clock_begin(id_diagnose_tracer_on_rho)
  do n=1,num_prog_tracers
     if(need_data(id_tracer_on_rho(n), next_time)) then 
         call diagnose_tracer_on_rho(Time, Dens, T_prog(n), n)
     endif
  enddo
  call mpp_clock_end(id_diagnose_tracer_on_rho)

  ! compute tracer*dz/drho as function of potential density 
  call mpp_clock_begin(id_diagnose_tracer_zrho_on_rho)
  do n=1,num_prog_tracers
     if(need_data(id_tracer_zrho_on_rho(n), next_time)) then 
         call diagnose_tracer_zrho_on_rho(Time, Dens, Thickness, T_prog(n), n)
     endif
  enddo
  call mpp_clock_end(id_diagnose_tracer_zrho_on_rho)


  ! compute depth of isopycnal surfaces and send to diagnostics manager
  call mpp_clock_begin(id_diagnose_depth_of_potrho)
  if(need_data(id_depth_of_potrho, next_time)) then 
    call diagnose_depth_of_potrho(Time, Dens, Thickness)
  endif
  call mpp_clock_end(id_diagnose_depth_of_potrho)

  ! compute depth of potential temp surfaces and send to diagnostics manager
  call mpp_clock_begin(id_diagnose_depth_of_theta)
  if(need_data(id_depth_of_theta, next_time)) then 
    call diagnose_depth_of_theta(Time, Dens, Thickness, T_prog)
  endif
  call mpp_clock_end(id_diagnose_depth_of_theta)

  call mpp_clock_begin(id_diagnose_mixing)

  ! quantify levels of effective mixing using sorting algorithm
  ! method available for tendency==TWO_LEVEL
  if (tendency==TWO_LEVEL) then 
    if(need_data(id_kappa_sort, next_time) .or. need_data(id_temp_sort,next_time)) then 
      call diagnose_kappa_sort(Time, Thickness, T_prog(index_temp))
    endif 
  endif 

  ! quantify levels of mixing using horizontal averaging algorithm 
  if(need_data(id_kappa_simple,next_time)) then
    call diagnose_kappa_simple(Time, T_prog(index_temp))
  endif 

  call mpp_clock_end(id_diagnose_mixing)

  call mpp_clock_begin(id_tracer_numerical)
  if (diag_step > 0) then
      if (mod(Time%itt, diag_step) == 0) then
          call dzt_min_max(Time, Thickness, 'From ocean_tracer_diag, dzt_min_max information')
          do n = 1,num_prog_tracers
             call tracer_min_max(Time, Thickness, T_prog(n))
             call tracer_integrals(Time, Thickness, T_prog(n))
          enddo
          call tracer_change(Time, Thickness, T_prog, T_diag, Ext_mode, &
                             pme, melt, runoff, calving)
          call tracer_land_cell_check (Time, T_prog)
      endif
  endif
  call mpp_clock_end(id_tracer_numerical)
  
  call mpp_clock_begin(id_conservation)
  if (nint(dtts) /= 0) then 
    call mass_conservation   (Time, Thickness, Ext_mode, pme, runoff, calving)
    call tracer_conservation (Time, Thickness, T_prog, T_diag, pme, runoff, calving)
  endif 
  call mpp_clock_end(id_conservation)

  ! send total seawater mass to diag_manager 
  call mpp_clock_begin(id_total_mass)
  if(need_data(id_mass_seawater, next_time)) then 
      call send_total_mass(Time, Thickness)
  endif 
  call mpp_clock_end(id_total_mass)

  ! send total seawater volume to diag_manager 
  call mpp_clock_begin(id_total_volume)
  if(need_data(id_volume_seawater, next_time)) then 
      call send_total_volume(Time, Thickness)
  endif 
  call mpp_clock_end(id_total_volume)

  ! global average pressure 
  if(need_data(id_press_ave, next_time)) then 
     call send_global_ave_pressure(Time, Thickness, Dens)
  endif 

  ! send total tracer to diag_manager 
  call mpp_clock_begin(id_total_tracer)
  do n=1,num_prog_tracers 
    if(need_data(id_prog_total(n), next_time)) then 
      call send_total_tracer(Time, Thickness, T_prog(n), n)
    endif 
    if(need_data(id_prog_global_ave(n), next_time)) then 
      call send_global_ave_tracer(Time, Thickness, T_prog(n), n)
    endif 
    if(need_data(id_prog_surface_ave(n), next_time)) then 
      call send_surface_ave_tracer(Time, Thickness, T_prog(n), n)
    endif 
    if(need_data(id_prog_surface_area_ave(n), next_time)) then 
      call send_surface_area_ave_tracer(Time, T_prog(n), n)
    endif 
  enddo
  call mpp_clock_end(id_total_tracer)

end subroutine ocean_tracer_diagnostics
! </SUBROUTINE>  NAME="ocean_tracer_diagnostics"


!#######################################################################
! <SUBROUTINE NAME="calc_mixed_layer_depth">
!
! <DESCRIPTION>
!
! Calculate the mixed layer depth (m), which is defined as the depth ( > 0 )
! where the buoyancy difference with respect to the surface level is
! equal to buoyancy_crit (m/s^2). 
!
! Note that the mixed layer depth is taken with respect to the ocean surface
! at z=eta_t, so the mixed layer depth is always positive. That is, the mld 
! is here defined as a thickness of water.
!            
! </DESCRIPTION>
!
subroutine calc_mixed_layer_depth(Thickness, salinity, theta, rho, pressure, &
                                  hmxl, smooth_mld_input)

  type(ocean_thickness_type),   intent(in)    :: Thickness
  real, dimension(isd:,jsd:,:), intent(in)    :: salinity
  real, dimension(isd:,jsd:,:), intent(in)    :: theta
  real, dimension(isd:,jsd:,:), intent(in)    :: rho
  real, dimension(isd:,jsd:,:), intent(in)    :: pressure 
  real, dimension(isd:,jsd:),   intent(inout) :: hmxl
  logical, optional,            intent(in)    :: smooth_mld_input

  real, parameter :: epsln=1.0e-20  ! for divisions 
  integer         :: i, j, k, km1, kb
  logical         :: smooth_mld_routine 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag_mod (calc_mixed_layer_depth): module needs initialization ')
  endif 

   if (present(smooth_mld_input)) then
    smooth_mld_routine = smooth_mld_input
  else 
    smooth_mld_routine = smooth_mld
  endif

  hmxl(:,:)   = 0.0
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0

  wrk1(:,:,:) = density_delta_sfc( rho(:,:,:), salinity(:,:,:), theta(:,:,:), pressure(:,:,:))
  do k=2,nk
     do j=jsc,jec
        do i=isc,iec
           wrk2(i,j,k) = -grav*Grd%tmask(i,j,k)*wrk1(i,j,k-1)/(epsln+rho(i,j,k))
        enddo
     enddo
  enddo

  do j=jsc,jec
     do i=isc,iec
        kb=Grd%kmt(i,j)
        if(kb==0) then
            hmxl(i,j) = 0.0
        else
            hmxl(i,j) = Thickness%depth_zwt(i,j,kb)
        endif
     enddo
  enddo

  do k=2,nk
     km1 = k-1  
     do j=jsc,jec
        do i=isc,iec
        kb=Grd%kmt(i,j)
           if (kb == 0) then
               hmxl(i,j) = 0.0
           else
               if ( wrk2(i,j,k) >= buoyancy_crit .and. hmxl(i,j)==Thickness%depth_zwt(i,j,kb)) then
                   hmxl(i,j) = Thickness%depth_zt(i,j,km1)                                &
                             - (Thickness%depth_zt(i,j,km1) - Thickness%depth_zt(i,j,k))  &
                             * (buoyancy_crit-wrk2(i,j,km1))                              &
                             / (wrk2(i,j,k) - wrk2(i,j,km1) + epsln)
               endif
               hmxl(i,j) = hmxl(i,j) * Grd%tmask(i,j,1)
           endif
        enddo
     enddo
  enddo

 ! smooth mld
  if(smooth_mld_routine) then 
      call mpp_update_domains(hmxl(:,:), Dom%domain2d) 
      hmxl(:,:) = S2D(hmxl(:,:))
  endif


end subroutine calc_mixed_layer_depth
! </SUBROUTINE>  NAME="calc_mixed_layer_depth"


!#######################################################################
! <SUBROUTINE NAME="mixed_layer_depth">
!
! <DESCRIPTION>
!
! Diagnose mixed layer depth (m).
! Call calc_mixed_layer_depth to determine the mixed layer depth.
! Also compute neutral density at depth of the mixed layer. 
! </DESCRIPTION>
!
subroutine mixed_layer_depth(Thickness, salinity, theta, rho, neutralrho, pressure, Time)

  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: salinity
  real, dimension(isd:,jsd:,:), intent(in) :: theta
  real, dimension(isd:,jsd:,:), intent(in) :: rho
  real, dimension(isd:,jsd:,:), intent(in) :: neutralrho
  real, dimension(isd:,jsd:,:), intent(in) :: pressure 
  type(ocean_time_type),        intent(in) :: Time

  integer :: i,j,k,kmt
  real, dimension(isd:ied,jsd:jed) :: hmxl

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag_mod (mixed_layer_depth): module needs initialization ')
  endif 

  call calc_mixed_layer_depth(Thickness, salinity, theta, rho, pressure, hmxl)

  call diagnose_2d(Time, Grd, id_mld, hmxl(:,:))
  if(id_mld_sq > 0) then 
     call diagnose_2d(Time, Grd, id_mld_sq, hmxl(:,:)**2)
  endif 

  ! compute neutral density at depth of hmxl 
  if(id_mld_nrho > 0) then
      wrk1_2d(:,:) = 0.0 
      do j=jsc,jec
         do i=isc,iec
            kmt = max(1,Grd%kmt(i,j)) 

            if(hmxl(i,j) < Thickness%depth_zt(i,j,1)) then
                wrk1_2d(i,j) = neutralrho(i,j,1) 
                cycle 
            endif
            if(hmxl(i,j) > Thickness%depth_zt(i,j,kmt)) then
                wrk1_2d(i,j) = neutralrho(i,j,kmt) 
                cycle 
            endif
            do k=2,nk
               if(hmxl(i,j) >  Thickness%depth_zt(i,j,k-1) .and. &
                  hmxl(i,j) <= Thickness%depth_zt(i,j,k)) then 
                   wrk1_2d(i,j) = ( neutralrho(i,j,k-1)*(Thickness%depth_zt(i,j,k)-hmxl(i,j))  &
                                   +neutralrho(i,j,k)  *(hmxl(i,j)-Thickness%depth_zt(i,j,k-1))&
                                   )/Thickness%dzwt(i,j,k-1) 
               endif
            enddo

         enddo
      enddo
      call diagnose_2d(Time, Grd, id_mld_nrho, wrk1_2d(:,:))
  endif


end subroutine mixed_layer_depth
! </SUBROUTINE>  NAME="mixed_layer_depth"



!#######################################################################
! <SUBROUTINE NAME="mixed_layer_depth_dtheta">
!
! <DESCRIPTION>
!
! Calculate the depth required to reach a temperature that is 
! dtheta cooler than the surface temperature. 
!
! Note:

! 1/ mixed_layer_depth_dtheta is taken with respect to the ocean surface
! at z=eta_t.
!            
! 2/ mixed_layer_depth_dtheta is no greater than the ocean depth + eta_t. 
!
! Coded March 2010 by Stephen.Griffies 
!
! </DESCRIPTION>
!
subroutine mixed_layer_depth_dtheta(Thickness, theta, Time)

  type(ocean_thickness_type),   intent(in)  :: Thickness
  real, dimension(isd:,jsd:,:), intent(in)  :: theta
  type(ocean_time_type),        intent(in)  :: Time

  real      :: epsln=1.0e-20  ! for divisions 
  real      :: dtheta 
  integer   :: i, j, k, km1, kb

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag_mod (mixed_layer_depth_dtheta): module needs initialization ')
  endif 

  wrk1(:,:,:)  = 0.0
  wrk1_2d(:,:) = 0.0

  do j=jsc,jec
     do i=isc,iec
        kb=Grd%kmt(i,j)
        if(kb==0) then
            wrk1_2d(i,j) = 0.0
        else
            wrk1_2d(i,j) = Thickness%depth_zwt(i,j,kb)
        endif
     enddo
  enddo

  do k=2,nk
     km1=k-1  
     do j=jsc,jec
        do i=isc,iec
           kb=Grd%kmt(i,j)
           if (kb == 0) then
               wrk1_2d(i,j) = 0.0
           else
               dtheta = theta(i,j,1) - theta(i,j,k)
               if (dtheta >= dtheta_crit .and. wrk1_2d(i,j)==Thickness%depth_zwt(i,j,kb)) then
                   wrk1_2d(i,j) = Thickness%depth_zt(i,j,km1)  &
                                + (dtheta_crit/(dtheta+epsln)) &
                                 *(Thickness%depth_zt(i,j,k) - Thickness%depth_zt(i,j,km1))
               endif
               wrk1_2d(i,j) = wrk1_2d(i,j) * Grd%tmask(i,j,1)
           endif
        enddo
     enddo
  enddo

  call diagnose_2d(Time, Grd, id_mld_dtheta, wrk1_2d(:,:))


end subroutine mixed_layer_depth_dtheta
! </SUBROUTINE>  NAME="mixed_layer_depth_dtheta"



!#######################################################################
! <SUBROUTINE NAME="compute_subduction">
!
! <DESCRIPTION>
! Diagnose subduction rate (kg/sec) based on kinematic method to compute 
! mass transport through base of mixed layer. 
!
! Some approximations made for convenience:
! 1/ use velocity at time tau
!    use tracer at time tau
!    use thickness at time taup1, since all pieces of Thickness
!    have already been updated to taup1.
!
! 2/ horizontally interpolate B-grid to C-grid velocity components,
!    but then vertically interpolate using mld computed at T-points. 
!
!  Stephen.Griffies
!  March 2012 
!  Updated for tracer transport July 2013
! </DESCRIPTION>
!
subroutine compute_subduction(Time, Thickness, Velocity, Adv_vel, Dens, T_prog, salinity, theta, diff_cbt)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness  ! taup1 values 
  type(ocean_velocity_type),    intent(in)    :: Velocity 
  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel    ! tau values 
  type(ocean_density_type),     intent(inout) :: Dens
  type(ocean_prog_tracer_type), intent(in)    :: T_prog(:)  ! for other than temp/saln tracers 
  real, dimension(isd:,jsd:,:), intent(in)    :: salinity   ! taup1 value
  real, dimension(isd:,jsd:,:), intent(in)    :: theta      ! taup1 value
  real, dimension(isd:,jsd:,:,:), intent(in)  :: diff_cbt

  integer :: i,j,k,kmt,n,nmix 
  integer :: tau, taup1 
  real    :: denominator, denominator_r
  real    :: W1, W2 
  real    :: dhdx, dhdy, dhdt
  real, dimension(isd:ied,jsd:jed)   :: mld_tau
  real, dimension(isd:ied,jsd:jed)   :: mld_taup1
  real, dimension(isd:ied,jsd:jed)   :: rho_mld
  real, dimension(isd:ied,jsd:jed)   :: subduction 
  real, dimension(isd:ied,jsd:jed)   :: subduction_dhdt
  real, dimension(isd:ied,jsd:jed)   :: subduction_horz
  real, dimension(isd:ied,jsd:jed)   :: subduction_vert
  real, dimension(isd:ied,jsd:jed,3) :: velocity_mld

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag_mod (compute_subduction): module needs initialization')
  endif 

  tau   = Time%tau
  taup1 = Time%taup1

  ! compute mld_taup1
  ! note Thickness is filled with taup1 thickness fields.   
  call calc_mixed_layer_depth(Thickness, salinity, theta,                &
                              Dens%rho(isd:ied,jsd:jed,:,taup1),         &
                              Dens%pressure_at_depth(isd:ied,jsd:jed,:), &
                              mld_taup1(:,:),           &
                              smooth_mld_for_subduction)

  ! update Dens%mld_subduction 
  do j=jsc,jec
     do i=isc,iec
        mld_tau(i,j)             = Dens%mld_subduction(i,j) 
        Dens%mld_subduction(i,j) = mld_taup1(i,j) 
     enddo
  enddo

  ! interpolate B-grid horizontal velocity to C-grid 
  wrk1_v(:,:,:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           denominator     = Grd%dus(i,j) + Grd%dun(i,j-1) + epsln
           wrk1_v(i,j,k,1) = Grd%tmask(i+1,j,k)*(Velocity%u(i,j,k,1,tau)  *Grd%dun(i,j-1) &
                                               + Velocity%u(i,j-1,k,1,tau)*Grd%dus(i,j))/denominator

           denominator     = Grd%duw(i,j) + Grd%due(i-1,j) + epsln
           wrk1_v(i,j,k,2) = Grd%tmask(i,j+1,k)*(Velocity%u(i-1,j,k,2,tau)*Grd%duw(i,j) &
                                               + Velocity%u(i,j,k,2,tau)  *Grd%due(i-1,j))/denominator
        enddo
     enddo
  enddo

  ! vertical velocity component at T-cell 
  wrk1(:,:,:) = 0.0
  if(vert_coordinate_class==DEPTH_BASED) then
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              wrk1(i,j,k) = rho0r*Adv_vel%wrho_bt(i,j,k)
           enddo
        enddo
     enddo
  else 
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              if(Grd%tmask(i,j,k) > 0.0) then
                  wrk1(i,j,k) = Adv_vel%wrho_bt(i,j,k)/Dens%rho(i,j,k,tau)
              endif
           enddo
        enddo
     enddo
  endif 

  ! vertically interpolate velocity components to mld_tau;
  ! ignore horizontal offset between C-grid velocity and mld.  
  velocity_mld(:,:,:) = 0.0
  rho_mld(:,:)        = 0.0

  ! shallow mld_tau 
  do j=jsc,jec
     do i=isc,iec
        if(mld_tau(i,j) <= Thickness%depth_zt(i,j,1)) then
           velocity_mld(i,j,1) = wrk1_v(i,j,1,1)
           velocity_mld(i,j,2) = wrk1_v(i,j,1,2)
           velocity_mld(i,j,3) = wrk1(i,j,1)
           rho_mld(i,j)        = Dens%rho(i,j,1,tau)
        endif
     enddo
  enddo

  ! intermediate mld_tau 
  do j=jsc,jec
     do i=isc,iec
kloop:  do k=1,nk-1
           if(Thickness%depth_zt(i,j,k) < mld_tau(i,j) .and. mld_tau(i,j) <= Thickness%depth_zt(i,j,k+1)) then 
              if(Grd%tmask(i,j,k+1) > 0) then
                 W1= mld_tau(i,j) - Thickness%depth_zt(i,j,k)
                 W2= Thickness%depth_zt(i,j,k+1) - mld_tau(i,j)
                 denominator_r = 1.0/(W1 + W2 + epsln)
                 velocity_mld(i,j,1) = (wrk1_v(i,j,k,1)*W2 + wrk1_v(i,j,k+1,1)*W1)*denominator_r
                 velocity_mld(i,j,2) = (wrk1_v(i,j,k,2)*W2 + wrk1_v(i,j,k+1,2)*W1)*denominator_r
                 velocity_mld(i,j,3) = (wrk1(i,j,k)    *W2 + wrk1(i,j,k+1)    *W1)*denominator_r
                 rho_mld(i,j)        = (Dens%rho(i,j,k,tau)*W2 + Dens%rho(i,j,k+1,tau)*W1)*denominator_r 
                 exit kloop
              endif
           endif
        enddo kloop
     enddo
  enddo

  ! deep mld_tau 
  do j=jsc,jec
     do i=isc,iec
        kmt = Grd%kmt(i,j)
        if(kmt > 0) then 
           if(mld_tau(i,j) > Thickness%depth_zt(i,j,kmt)) then
              velocity_mld(i,j,1) = wrk1_v(i,j,kmt,1)
              velocity_mld(i,j,2) = wrk1_v(i,j,kmt,2)
              velocity_mld(i,j,3) = wrk1(i,j,kmt)
              rho_mld(i,j)        = Dens%rho(i,j,kmt,tau)
           endif
        endif 
     enddo
  enddo

  ! over-write rho_mld for Boussinesq case
  if(vert_coordinate_class==DEPTH_BASED) then
     do j=jsc,jec
        do i=isc,iec
           rho_mld(i,j) = rho0*Grd%tmask(i,j,1)
        enddo
     enddo
  endif 

  ! save rho*dat 
  wrk1_2d(:,:) = 0.0 
  do j=jsc,jec
     do i=isc,iec
        wrk1_2d(i,j) = rho_mld(i,j)*Grd%dat(i,j)
     enddo
  enddo

  ! compute contributions to subduction 
  call mpp_update_domains(mld_tau(:,:), Dom%domain2d) 
  subduction_dhdt(:,:) = 0.0
  subduction_horz(:,:) = 0.0
  subduction_vert(:,:) = 0.0
  subduction(:,:)      = 0.0
  do j=jsc,jec
     do i=isc,iec
       dhdx                 = (mld_tau(i+1,j)-mld_tau(i,j))*Grd%dxter(i,j) 
       dhdy                 = (mld_tau(i,j+1)-mld_tau(i,j))*Grd%dytnr(i,j)
       dhdt                 = (mld_taup1(i,j) - mld_tau(i,j))*dtimer
       subduction_dhdt(i,j) = -wrk1_2d(i,j)*dhdt
       subduction_horz(i,j) = -wrk1_2d(i,j)*(dhdx*velocity_mld(i,j,1)+dhdy*velocity_mld(i,j,2))
       subduction_vert(i,j) = -wrk1_2d(i,j)*velocity_mld(i,j,3)
       subduction(i,j)      = subduction_dhdt(i,j) + subduction_horz(i,j) + subduction_vert(i,j)  
     enddo
  enddo

  call diagnose_2d(Time, Grd, id_subduction_dhdt, subduction_dhdt(:,:))
  call diagnose_2d(Time, Grd, id_subduction_horz, subduction_horz(:,:))
  call diagnose_2d(Time, Grd, id_subduction_vert, subduction_vert(:,:))
  call diagnose_2d(Time, Grd, id_subduction, subduction(:,:))
  call diagnose_2d(Time, Grd, id_subduction_mld, mld_tau(:,:))

  if(id_subduction_dhdt_nrho > 0) then
     wrk1(:,:,:)      = 0.0
     k=1
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = subduction_dhdt(i,j)
        enddo
     enddo
     call diagnose_3d_rho(Time, Dens, id_subduction_dhdt_nrho, wrk1)
  endif

  if(id_subduction_horz_nrho > 0) then
     wrk1(:,:,:)      = 0.0
     k=1
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = subduction_horz(i,j)
        enddo
     enddo
     call diagnose_3d_rho(Time, Dens, id_subduction_horz_nrho, wrk1)
  endif

  if(id_subduction_vert_nrho > 0) then
     wrk1(:,:,:)      = 0.0
     k=1
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = subduction_vert(i,j)
        enddo
     enddo
     call diagnose_3d_rho(Time, Dens, id_subduction_vert_nrho, wrk1)
  endif

  if(id_subduction_nrho > 0) then
     wrk1(:,:,:)      = 0.0
     k=1
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = subduction(i,j)
        enddo
     enddo
     call diagnose_3d_rho(Time, Dens, id_subduction_nrho, wrk1)
  endif



  !!!!!!! tracer diagnostic calculations !!!!!!!  

  do n=1,num_prog_tracers

     if (n==index_salt) then
        nmix=2
     else
        nmix=1
     endif

     ! compute downgradient vertical diffusive flux (assume time-implicit calculation)
     wrk1(:,:,:) = 0.0
     if(id_subduction_mld_zflux_diff(n) > 0) then   
        do k=1,nk-1
           do j=jsc,jec
              do i=isc,iec
                 wrk1(i,j,k) = -1.0*rho0*Grd%tmask(i,j,k+1)*diff_cbt(i,j,k,nmix)            &
                              *(T_prog(n)%field(i,j,k,taup1)-T_prog(n)%field(i,j,k+1,taup1))&
                              /Thickness%dzwt(i,j,k)
              enddo
           enddo
        enddo
     endif 

     if(    id_subduction_tracer(n)     >0 .or. id_subduction_dhdt_tracer(n)   >0 .or. id_subduction_horz_tracer(n)>0 &
       .or. id_subduction_vert_tracer(n)>0 .or. id_subduction_mld_zflux_diff(n)>0) then   

        ! interpolate tracer concentration and vertical diffusive flux to base of mld 
        wrk1_2d(:,:) = 0.0  ! interpolated tracer 
        wrk2_2d(:,:) = 0.0  ! interpolated zflux 

        ! shallow mld_tau 
        do j=jsc,jec
           do i=isc,iec
              if(mld_tau(i,j) <= Thickness%depth_zt(i,j,1)) then
                 wrk1_2d(i,j) = T_prog(n)%field(i,j,k,tau)
                 wrk2_2d(i,j) = wrk1(i,j,k)
              endif
           enddo
        enddo

        ! intermediate mld_tau 
        do j=jsc,jec
           do i=isc,iec
     kloopB:  do k=1,nk-1
                 if(Thickness%depth_zt(i,j,k) < mld_tau(i,j) .and. mld_tau(i,j) <= Thickness%depth_zt(i,j,k+1)) then 
                    if(Grd%tmask(i,j,k+1) > 0) then
                       W1= mld_tau(i,j) - Thickness%depth_zt(i,j,k)
                       W2= Thickness%depth_zt(i,j,k+1) - mld_tau(i,j)
                       denominator_r = 1.0/(W1 + W2 + epsln)
                       wrk1_2d(i,j)  = (T_prog(n)%field(i,j,k,tau)*W2 + T_prog(n)%field(i,j,k+1,tau)*W1)*denominator_r
                       wrk2_2d(i,j)  = (wrk1(i,j,k)*W2 + wrk1(i,j,k+1)*W1)*denominator_r 
                       exit kloopB
                    endif
                 endif
              enddo kloopB
           enddo
        enddo

        ! deep mld_tau 
        do j=jsc,jec
           do i=isc,iec
              kmt = Grd%kmt(i,j)
              if(kmt > 0) then 
                 if(mld_tau(i,j) > Thickness%depth_zt(i,j,kmt)) then
                    wrk1_2d(i,j) = T_prog(n)%field(i,j,kmt,tau)
                    wrk2_2d(i,j) = wrk1(i,j,kmt)
                 endif
              endif 
           enddo
        enddo

    endif 
 
    ! write diagnostics 
    wrk3_2d(:,:) = 0.0

    if(id_subduction_tracer(n) > 0) then 
       wrk3_2d(:,:) = subduction(:,:)*wrk1_2d(:,:)
       used = send_data (id_subduction_tracer(n), wrk3_2d(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif 
    if(id_subduction_dhdt_tracer(n) > 0) then 
       wrk3_2d(:,:) = subduction_dhdt(:,:)*wrk1_2d(:,:)
       used = send_data (id_subduction_dhdt_tracer(n), wrk3_2d(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif 
    if(id_subduction_horz_tracer(n) > 0) then 
       wrk3_2d(:,:) = subduction_horz(:,:)*wrk1_2d(:,:)
       used = send_data (id_subduction_horz_tracer(n), wrk3_2d(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif 
    if(id_subduction_vert_tracer(n) > 0) then 
       wrk3_2d(:,:) = subduction_vert(:,:)*wrk1_2d(:,:)
       used = send_data (id_subduction_vert_tracer(n), wrk3_2d(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif 
    if(id_subduction_mld_zflux_diff(n) > 0) then 
       wrk3_2d(:,:) = wrk2_2d(:,:)*Grd%dat(:,:)
       used = send_data (id_subduction_mld_zflux_diff(n), wrk3_2d(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif 


  enddo   ! end of n-loop 

end subroutine compute_subduction
! </SUBROUTINE>  NAME="compute_subduction"


!#######################################################################
! <SUBROUTINE NAME="compute_tracer_mld">
!
! <DESCRIPTION>
!
! Diagnose the vertically integrated tracer within the mixed layer.
!
! Work with taup1 values, since Thickness type is filled with taup1
! thickness fields.   
!
! Stephen.Griffies
! July 2013
! </DESCRIPTION>
!
subroutine compute_tracer_mld(Time, Thickness, Dens, T_prog, salinity, theta)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness  ! taup1 values 
  type(ocean_density_type),     intent(inout) :: Dens
  type(ocean_prog_tracer_type), intent(in)    :: T_prog(:) 
  real, dimension(isd:,jsd:,:), intent(in)    :: salinity   ! taup1 value
  real, dimension(isd:,jsd:,:), intent(in)    :: theta      ! taup1 value

  integer :: i,j,k,kp1,n
  integer :: taup1
  real, dimension(isd:ied,jsd:jed) :: mld

  
  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag_mod (compute_tracer_mld): module needs initialization')
  endif

  if(.not. compute_tracer_mld_diags) return 

  taup1 = Time%taup1

  ! compute mld using taup1 values. 
  ! Note the Thickness type is filled with taup1 thickness fields.   
  call calc_mixed_layer_depth(Thickness, salinity, theta,                &
                              Dens%rho(isd:ied,jsd:jed,:,taup1),         &
                              Dens%pressure_at_depth(isd:ied,jsd:jed,:), &
                              mld(:,:), smooth_mld_for_subduction)

 
  ! determine column masking function
  wrk1(:,:,:) = 0.0
  k=1
  do j=jsc,jec
     do i=isc,iec
        if(Grd%tmask(i,j,k)==1.0) then 
            if(Thickness%depth_zwt(i,j,k) >= mld(i,j)) then
                wrk1(i,j,1)    = mld(i,j)/Thickness%depth_zwt(i,j,k)
                wrk1(i,j,2:nk) = 0.0
            endif
        endif
     enddo
  enddo
  
  ! k>1 
  do j=jsc,jec
     do i=isc,iec
        kloopA:    do k=2,nk
           if(Grd%tmask(i,j,k)==1.0) then 
               if(Thickness%depth_zwt(i,j,k)   >= mld(i,j) .and. &
                  Thickness%depth_zwt(i,j,k-1) <  mld(i,j)) then
                   kp1 = min(k+1,nk)
                   wrk1(i,j,1:k-1)  = 1.0
                   wrk1(i,j,k)      = (mld(i,j)-Thickness%depth_zwt(i,j,k-1))/Thickness%dzt(i,j,k)
                   wrk1(i,j,kp1:nk) = 0.0
                   exit kloopA
               endif
           endif
        enddo kloopA
     enddo
  enddo

  k=nk
  do j=jsc,jec
     do i=isc,iec
        if(Grd%tmask(i,j,k)==1.0) then 
            if(Thickness%depth_zwt(i,j,k) <= mld(i,j)) then
                wrk1(i,j,:) = 1.0
            endif
        endif
     enddo
  enddo
 
  ! compute the vertically integrated tracer within the mixed layer   
  do n=1,num_prog_tracers
     if(id_tracer_mld(n) > 0) then
         wrk1_2d(:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1_2d(i,j) = wrk1_2d(i,j) &
                               + wrk1(i,j,k)*Thickness%rho_dzt(i,j,k,taup1)*T_prog(n)%field(i,j,k,taup1)
               enddo
            enddo
         enddo
         used = send_data (id_tracer_mld(n), wrk1_2d(:,:), &
         Time%model_time, rmask=Grd%tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif 
  enddo
  

end subroutine compute_tracer_mld
! </SUBROUTINE>  NAME="compute_tracer_mld"



!#######################################################################
! <SUBROUTINE NAME="tracer_change">
!
! <DESCRIPTION>
!
! Compute change in tracer over a time step and difference between 
! this change and the boundary forcing. 
!
! This routine is very useful for detecting bugs in tracer routines.  
! 
! </DESCRIPTION>
subroutine tracer_change (Time, Thickness, T_prog, T_diag, Ext_mode, &
                          pme, melt, runoff, calving)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_thickness_type),     intent(in) :: Thickness
  type(ocean_prog_tracer_type),   intent(in) :: T_prog(:) 
  type(ocean_diag_tracer_type),   intent(in) :: T_diag(:) 
  type(ocean_external_mode_type), intent(in) :: Ext_mode
  real, dimension(isd:,jsd:),     intent(in) :: pme
  real, dimension(isd:,jsd:),     intent(in) :: melt
  real, dimension(isd:,jsd:),     intent(in) :: runoff
  real, dimension(isd:,jsd:),     intent(in) :: calving

  real, dimension(isd:ied,jsd:jed) :: tmp
  real, dimension(isd:ied,jsd:jed) :: tmpE
  real, dimension(isd:ied,jsd:jed) :: tmpL
  real, dimension(isd:ied,jsd:jed) :: area_k

  real :: mass_error, tracer_error, flux_error
  real :: tracer_input_stf, tracer_input_btf, tracer_input_pme
  real :: tracer_input_otf
  real :: tracer_input_eta_smooth, tracer_input_pbot_smooth
  real :: tracer_input_runoff, tracer_input_calving, tracer_input_total
  real :: pme_input, runoff_input, calving_input, melt_input, obc_input
  real :: eta_t_max, eta_t_min, eta_t_avg
  real :: anompbot_max, anompbot_min
  real :: frazil_heat 
  real :: tracer_th_tend, tracer_src 
  real :: ttracer 
  real :: tracer_L, tracer_E
  real :: dmasstracer
  real :: dmasstracer_E, dmasstracer_L, dmasstracer_L2
  real :: mass_total, mass_source, mass_change
  real :: mass_E, mass_L, mass_T
  real :: mass_changeE, mass_changeL, mass_changeL2
  real :: mass_taup1, mass_taup1_r, mass_taup1_keq1
  real :: mass_eta_smooth, mass_pbot_smooth
  real :: conversion, conv_blob

  integer :: i,j,k,n
  integer :: taum1, tau, taup1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  call mpp_clock_begin(id_tracer_change)

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_tracer_diag_mod (tracer_change): module needs initialization ')
  endif 

  if(size(T_prog,1) /= num_prog_tracers) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_tracer_diag_mod (tracer_change): passing T_prog with wrong dimensions')
  endif 

  write(stdoutunit,*) ' '
  write (stdoutunit,'(//50x,a)') ' Mass and tracer change summary (over one timestep):'
  call write_timestamp(Time%model_time)
  write(stdoutunit,*) ' '


  ! compute total mass and changes over a time step 

  pme_input        = 0.0 ! total mass of evap-precip input to ocean over time step (kg)
  melt_input       = 0.0 ! total mass of ice melt input to ocean over time step (kg)
  runoff_input     = 0.0 ! total mass of runoff water input to ocean over time step (kg)
  calving_input    = 0.0 ! total mass of calving water input to ocean over time step (kg)
  obc_input        = 0.0 ! total mass of water input through open boundaries to ocean over time step (kg)
  eta_t_avg        = 0.0 ! area average of eta_t (m)
  eta_t_max        = 0.0 ! global max of eta_t (m)
  eta_t_min        = 0.0 ! global min of eta_t (m)
  anompbot_max     = 0.0 ! global max of pbot_t-pbot0 (dbar)
  anompbot_min     = 0.0 ! global min of pbot_t-pbot0 (dbar)
  mass_total       = 0.0 ! mass of ocean on tracer cells (kg)
  mass_taup1_keq1  = 0.0 ! mass of ocean in the k=1 cell (kg)
  mass_eta_smooth  = 0.0 ! mass from eta_t surface filtering (kg) 
  mass_pbot_smooth = 0.0 ! mass from pbot_t filtering (kg)
  mass_change      = 0.0 ! change in mass over time step (kg)
  mass_source      = 0.0 ! mass from sources (kg)
 
  mass_L        = 0.0 ! mass of L system
  mass_E        = 0.0 ! mass of E system
  mass_T        = 0.0
  mass_changeL  = 0.0 ! change in mass of L system
  mass_changeE  = 0.0 ! change in mass of E system
  mass_changeL2 = 0.0 ! mass transfer from L to E system
                      ! note: if all is well, mass_changeL = mass_changeL2
  conv_blob     = 0.0 ! tendency due to convergence of blobs (should be close to zero)

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  ! inverse mass of ocean water at time taup1 
  mass_taup1   = total_mass(Thickness, taup1)
  mass_taup1_r = 1/(mass_taup1+epsln)

  ! mass of ocean water at time taup1 and k-level=1 
  mass_taup1_keq1 = klevel_total_mass(Thickness, taup1, 1)

  ! mass of ocean water at time tau
  tmpE(:,:) = 0.0
  do k=1,nk
     tmpE(:,:) = tmpE(:,:) + Grd%tmask(:,:,k)*area_t(:,:)*Thickness%rho_dzt(:,:,k,tau)
  enddo
  mass_E     = mpp_global_sum(Dom%domain2d, tmpE(:,:), global_sum_flag)
  mass_total = mass_E

  ! mass of ocean water at time tau, no with blobs 
  if(use_blobs) then
     tmpL(:,:) = 0.0
     tmp(:,:)  = 0.0
     do k=1,nk
        tmpL(:,:) = tmpL(:,:) + Grd%tmask(:,:,k)*area_t(:,:)*Thickness%rho_dztL(:,:,k,tau)
        tmp(:,:)  = tmp(:,:)  + Grd%tmask(:,:,k)*area_t(:,:)*Thickness%rho_dztT(:,:,k,tau)
     enddo
     mass_L     = mpp_global_sum(Dom%domain2d, tmpL(:,:), global_sum_flag)
     mass_T     = mpp_global_sum(Dom%domain2d, tmp(:,:),  global_sum_flag)
     mass_total = mass_L + mass_E                                                
  endif
 
  ! mass of ice melt input to the ocean 
  tmp(:,:)  = dtime*melt(:,:)*area_t(:,:)
  melt_input = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

  ! mass of pme input to the ocean. subtract melt_input 
  ! since ice melt is already contained in pme. 
  tmp(:,:)  = dtime*pme(:,:)*area_t(:,:)
  pme_input = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag) - melt_input

  ! mass of river runoff water input to the ocean
  tmp(:,:)     = dtime*runoff(:,:)*area_t(:,:)
  runoff_input = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

  ! mass of calving land ice input to the ocean
  tmp(:,:)      = dtime*calving(:,:)*area_t(:,:)
  calving_input = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

  ! area average surface height at time tau
  tmp(:,:)   = area_t(:,:)*Ext_mode%eta_t(:,:,tau)
  eta_t_avg  = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)/Grd%tcellsurf
  
  ! max and min of eta_t
  eta_t_max = mpp_global_max(Dom%domain2d,Ext_mode%eta_t(:,:,tau))
  eta_t_min = mpp_global_min(Dom%domain2d,Ext_mode%eta_t(:,:,tau))

  ! max and min of pbot_t-pbot0
  anompbot_max = c2dbars &
   *mpp_global_max(Dom%domain2d,Ext_mode%pbot_t(:,:,tau)-Thickness%pbot0(:,:))
  anompbot_min = c2dbars &
   *mpp_global_min(Dom%domain2d,Ext_mode%pbot_t(:,:,tau)-Thickness%pbot0(:,:))

  ! change in mass over a time step 
  tmpE(:,:) = 0.0
  do k=1,nk
    tmpE(:,:) = tmpE(:,:) + Grd%tmask(:,:,k)*area_t(:,:) &
                           *(Thickness%rho_dzt(:,:,k,taup1)-Thickness%rho_dzt(:,:,k,taum1))
  enddo
  mass_changeE = mpp_global_sum(Dom%domain2d,tmpE(:,:), global_sum_flag)
  mass_change  = mass_changeE

  ! mass change over a time step with blobs 
  if (use_blobs) then
     tmpL(:,:) = 0.0
     do k=1,nk
        tmpL(:,:) = tmpL(:,:) + Grd%tmask(:,:,k)*area_t(:,:) &
                               *(Thickness%rho_dztL(:,:,k,taup1)-Thickness%rho_dztL(:,:,k,taum1))
     enddo
     mass_changeL = mpp_global_sum(Dom%domain2d,tmpL(:,:), global_sum_flag)
     mass_change  = mass_changeE +  mass_changeL

     tmp(:,:)      = dtime*Grd%tmask(:,:,1)*Grd%dat(:,:)*Thickness%blob_source(:,:)
     mass_changeL2 = mass_changeL2 + mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

     tmp(:,:)  = dtime*area_t(:,:)*Ext_mode%conv_blob(:,:)
     conv_blob = mpp_global_sum(Dom%domain2d, tmp(:,:), global_sum_flag)
  endif

  ! mass input to ocean from sources 
  tmp(:,:) = 0.0
  do k=1,nk    
    tmp(:,:) = tmp(:,:) + dtime*area_t(:,:)*Thickness%mass_source(:,:,k)
  enddo
  mass_source = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

  ! mass input due to smoothing of the surface height (should be roundoff).
  ! note that this contribution is already in mass_source 
  tmp(:,:)        = dtime*area_t(:,:)*Ext_mode%eta_smooth(:,:)
  mass_eta_smooth = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

  ! mass input due to smoothing of the bottom pressure (should be roundoff). 
  ! note that this contribution is already in mass_source 
  tmp(:,:)         = dtime*area_t(:,:)*Ext_mode%pbot_smooth(:,:)
  mass_pbot_smooth = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

  ! compute mass passing across open boundaries 
  if (have_obc) then 
    call ocean_obc_mass_flux(Time, Ext_mode, tmp)
    tmp(:,:)  = tmp(:,:)*dtime
    obc_input = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag) 
  endif 

  mass_error = mass_change-pme_input-melt_input-runoff_input-calving_input-obc_input-mass_source

  write (stdoutunit,'(/a)') ' ----Single time step diagnostics for volume and mass----' 

  write (stdoutunit,'(a,es24.17,a)') ' Global max (pbot_t-pbot0) on tracer cells (tau)   = ',&
                                   anompbot_max,' dbar'
  write (stdoutunit,'(a,es24.17,a)') ' Global min (pbot_t-pbot0) on tracer cells (tau)   = ',&
                                   anompbot_min,' dbar'
  write (stdoutunit,'(a,es24.17,a)') ' Global max surface height on tracer cells (tau)   = ',&
                                   eta_t_max,' m'
  write (stdoutunit,'(a,es24.17,a)') ' Global min surface height on tracer cells (tau)   = ',&
                                   eta_t_min,' m'
  write (stdoutunit,'(a,es24.17,a)') ' Area average surface height on tracer cells (tau) = ',&
                                   eta_t_avg,' m'
  write (stdoutunit,'(a,es24.17,a)') ' Surface area of k=1 tracer cells                  = ',&
                                   Grd%tcellsurf,' m^2'
  write (stdoutunit,'(a,es24.17,a)') ' Area integral of surface height on tracer cells   = ',&
                                   eta_t_avg*Grd%tcellsurf,' m^3'
  if (use_blobs) then
     write (stdoutunit,'(a,es24.17,a)') ' L system mass change (taup1 - taum1)              = ',&
                                     mass_changeL,' kg'
     write (stdoutunit,'(a,es24.17,a)') ' E system mass change (taup1 - taum1)              = ',&
                                     mass_changeE,' kg'
     write (stdoutunit,'(a,es24.17,a)') ' Transfer of mass from L to E (taup1)              = ',&
                                     mass_changeL2,' kg'
  endif
  write (stdoutunit,'(a,es24.17,a)') ' Total Mass change                                 = ',&
                                     mass_change,' kg'  
  write (stdoutunit,'(a,es24.17,a)') ' Mass of ocean tracer cells (tau)                  = ',&
                                   mass_total,' kg'
  if (use_blobs) then
     write (stdoutunit,'(a,es24.17,a)') ' Mass of L system ocean tracer cells (tau)         = ',&
                                     mass_L,' kg'
     write (stdoutunit,'(a,es24.17,a)') ' Mass of E system ocean tracer cells (tau)         = ',&
                                     mass_E,' kg'
     write (stdoutunit,'(a,es24.17,a)') ' Total Mass of ocean tracer cells (tau)            = ',&
                                     mass_total,' kg'
     write (stdoutunit,'(a,es24.17,a)') ' Mismatch between T Mass and L + E mass (T-[L+E])  = ',&
                                     (mass_T - (mass_total)),' kg'
     write (stdoutunit,'(a,es24.17,a)') ' Mass from blob convergence (should be small)      = ',&
                                   conv_blob,' kg'
  endif
  write (stdoutunit,'(a,es24.17,a)') ' Mass from sources                                 = ',&
                                   mass_source,' kg'
  write (stdoutunit,'(a,es24.17,a)') ' Mass from eta smoothing                           = ',&
                                   mass_eta_smooth,' kg'
  write (stdoutunit,'(a,es24.17,a)') ' Mass from pbot smoothing                          = ',&
                                   mass_pbot_smooth,' kg'
  write (stdoutunit,'(a,es24.17,a)') ' Mass of precip-evap input                         = ',&
                                   pme_input,' kg'
  write (stdoutunit,'(a,es24.17,a)') ' Mass of river runoff liquid water input           = ',&
                                   runoff_input,' kg'
  write (stdoutunit,'(a,es24.17,a)') ' Mass of sea ice melt input                        = ',&
                                   melt_input,' kg'
  write (stdoutunit,'(a,es24.17,a)') ' Mass of calving land ice input                    = ',&
                                   calving_input,' kg'
  write (stdoutunit,'(a,es24.17,a)') ' Mass flux through open boundaries                 = ',&
                                   obc_input,' kg'
  write (stdoutunit,'(a,es24.17,a)') ' Mismatch between mass change and water input      = ',&
                                   mass_error,' kg'


  ! compute total tracer and changes over a time step 

  do n=1,num_prog_tracers

     dmasstracer         = 0.0 ! (rho*vol)(taup1)*tracer(taup1)-(rho*vol)(taum1)*tracer(taum1))
     tracer_input_stf    = 0.0 ! tracer quantity input to ocean via stf (>0 means tracer entering ocean)
     tracer_input_btf    = 0.0 ! tracer quantity input to ocean via btf (>0 means tracer leaving ocean)
     tracer_input_otf    = 0.0 ! tracer quantity input through open boundaries (>0 means tracer entering ocean) 
     tracer_input_pme    = 0.0 ! tracer quantity input to ocean via pme flux 
     tracer_input_runoff = 0.0 ! tracer quantity input to ocean via river runoff
     tracer_input_calving= 0.0 ! tracer quantity input to ocean via calving land ice
     tracer_input_total  = 0.0 ! total tracer quantity input to ocean 
     frazil_heat         = 0.0 ! heat (Joule) from frazil formation  
     tracer_th_tend      = 0.0 ! contribution from tracer tendencies including 
                               ! neutral diffusion, sigma diffusion, overflow, xlandmix, 
                               ! sponges, rivers, and other more traditional sources. 
     tracer_src          = 0.0 ! sources that have dimensions tracer concentration per time 
     tracer_input_eta_smooth = 0.0 ! tracer input at surface due to smoothing eta_t (should be roundoff)
     tracer_input_pbot_smooth= 0.0 ! tracer input bottom due to smoothing pbot_t  (should be roundoff)

     tracer_L       = 0.0 ! tracer quantity in the L system
     tracer_E       = 0.0 ! tracer quantity in the E system
     dmasstracer_L  = 0.0 ! change in tracer quantity of the L system
     dmasstracer_L2 = 0.0 ! tracer quantity transferred from the L to E system
     dmasstracer_E  = 0.0 ! change in tracer quantity of the E system

     conversion = T_prog(n)%conversion     

     tmp(:,:) = area_t(:,:)*conversion*dtime*T_prog(n)%stf(:,:)
     tracer_input_stf = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

     tmp(:,:) = area_t(:,:)*conversion*dtime*T_prog(n)%btf(:,:)
     tracer_input_btf = -mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

     tmp(:,:) = area_t(:,:)*conversion*dtime*pme(:,:)*T_prog(n)%tpme(:,:)
     tracer_input_pme = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag) 

     ! runoff can be negative, as for the "geoengineering" trick of siphoning
     ! Arabian Sea water to the Caspian Sea, to keep the Caspian from evaporating
     ! away during a coupled climate simulation.  As negative runoff is ignored  
     ! in the rivermix module (rightly so, as we do not wish to insert tracer 
     ! from negative runoff into the ocean), we also wish to ignore negative 
     ! runoff for the diagnostics here, in order to reflect what the prognostic
     ! equations are doing.      
     tmp(:,:) = 0.0
     do j=jsc,jec
        do i=isc,iec
           if(runoff(i,j) > 0) then        
               tmp(i,j) = area_t(i,j)*conversion*dtime*runoff(i,j)*T_prog(n)%trunoff(i,j)        
           endif
        enddo
     enddo
     tracer_input_runoff = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

     tmp(:,:) = area_t(:,:)*conversion*dtime*calving(:,:)*T_prog(n)%tcalving(:,:)
     tracer_input_calving = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

     tmp(:,:) = area_t(:,:)*conversion*dtime*T_prog(n)%eta_smooth(:,:)
     tracer_input_eta_smooth = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

     tmp(:,:) = area_t(:,:)*conversion*dtime*T_prog(n)%pbot_smooth(:,:)
     tracer_input_pbot_smooth = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

     ! ocean_obc_tracer_flux returns vertically integrated horizontal tracer flux     
     ! at last internal points and zero otherwise     
     if(have_obc) then 
       call ocean_obc_tracer_flux(T_prog(n), tmp) 
       tmp(:,:) = tmp(:,:)*conversion*dtime
       tracer_input_otf = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)
     endif 

     ! river(=calving+runoff) is added via T_prog%th_tendency inside ocean_rivermix_mod
     tracer_input_total =   tracer_input_stf    + tracer_input_btf     + tracer_input_otf &
                          + tracer_input_runoff + tracer_input_calving + tracer_input_pme &
                          + tracer_input_eta_smooth + tracer_input_pbot_smooth

     frazil_heat = 0.0
     if(T_prog(n)%name =='temp' .and. index_frazil > 0) then
       tmp(:,:) = 0.0
       do k=1,nk
         tmp(:,:) = tmp(:,:) + T_diag(index_frazil)%field(:,:,k)*area_t(:,:)*Grd%tmask(:,:,k)/frazil_factor 
       enddo
       frazil_heat        = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag) 
       tracer_input_total = tracer_input_total + frazil_heat
     endif 

     ! tracer changes and total source
     tmp(:,:)     = 0.0    
     tmpE(:,:)    = 0.0    
     wrk1_2d(:,:) = 0.0
     wrk2_2d(:,:) = 0.0
     do k=1,nk
        area_k(:,:) = Grd%dat(:,:)*Grd%tmask(:,:,k)  
        tmpE(:,:)   = tmpE(:,:) + area_k(:,:)*conversion*                           &
                      (Thickness%rho_dzt(:,:,k,taup1)*T_prog(n)%field(:,:,k,taup1)  &
                      -Thickness%rho_dzt(:,:,k,taum1)*T_prog(n)%field(:,:,k,taum1))
        wrk1_2d(:,:) = wrk1_2d(:,:) + area_k(:,:)*conversion*dtime*T_prog(n)%th_tendency(:,:,k)
        wrk2_2d(:,:) = wrk2_2d(:,:) + area_k(:,:)*conversion*dtime*T_prog(n)%source(:,:,k)*Thickness%rho_dzt(:,:,k,taup1)
     enddo
     if(have_obc) then 
        tmpE(:,:)    = tmpE(:,:)*Grd%obc_tmask(:,:)
        wrk1_2d(:,:) = wrk1_2d(:,:)*Grd%obc_tmask(:,:)
        wrk2_2d(:,:) = wrk2_2d(:,:)*Grd%obc_tmask(:,:)
     endif 
     tracer_th_tend = mpp_global_sum(Dom%domain2d, wrk1_2d(:,:), global_sum_flag)
     tracer_src     = mpp_global_sum(Dom%domain2d, wrk2_2d(:,:), global_sum_flag)
     dmasstracer_E  = mpp_global_sum(Dom%domain2d, tmpE(:,:)   , global_sum_flag)
     dmasstracer    = dmasstracer_E

     ! tracer changes and total source, now with blobs 
     if(use_blobs) then 
        tmpL(:,:) = 0.0    
        tmp(:,:)  = 0.0    
        do k=1,nk
           area_k(:,:) = Grd%dat(:,:)*Grd%tmask(:,:,k)  
           tmpL(:,:)   = tmpL(:,:) + conversion*(T_prog(n)%sum_blob(:,:,k,taup1) - T_prog(n)%sum_blob(:,:,k,taum1))
           tmp(:,:)    = tmp(:,:)  + dtime*conversion*area_k(:,:)*(T_prog(n)%tend_blob(:,:,k))           
        enddo
        if(have_obc) then 
           tmpL(:,:) = tmpL(:,:)*Grd%obc_tmask(:,:)
           tmp(:,:)  = tmp(:,:)*Grd%obc_tmask(:,:)
        endif 
        dmasstracer_L  = dmasstracer_L  + mpp_global_sum(Dom%domain2d, tmpL(:,:), global_sum_flag)
        dmasstracer_L2 = dmasstracer_L2 + mpp_global_sum(Dom%domain2d, tmp(:,:),  global_sum_flag)
        dmasstracer    = dmasstracer_L  + dmasstracer_E 
     endif 
 
     ! for cleaner diagnostics, remove contributions from runoff, calving, pme, and smoother,
     ! each of which are added to T_prog%th_tendency inside of ocean_tracer.F90.
     tracer_th_tend = tracer_th_tend -tracer_input_runoff     -tracer_input_calving -tracer_input_pme &
                                     -tracer_input_eta_smooth -tracer_input_pbot_smooth

     ! for aidif=0.0, stf is included in T_prog%th_tendency through vertical diffusion.
     ! for aidif=1.0, it is not in th_tendency; rather, it is included through invtri.
     if(aidif==0.0) then 
       tracer_th_tend = tracer_th_tend - tracer_input_stf
     endif 

     ! tracer in E system
     tmpE(:,:) = 0.0
     do k=1,nk 
        tmpE(:,:) = tmpE(:,:) + conversion*area_t(:,:)*Thickness%rho_dzt(:,:,k,tau)*T_prog(n)%field(:,:,k,tau)
     enddo
     tracer_E = mpp_global_sum(Dom%domain2d, tmpE(:,:), global_sum_flag)
     ttracer  = tracer_E

     ! tracer in L system
     if (use_blobs) then
        tmpL(:,:) = 0.0
        do k=1,nk
           tmpL(:,:) = tmpL(:,:) + conversion*Grd%tmask(:,:,k)*(T_prog(n)%sum_blob(:,:,k,tau))
        enddo
        tracer_L = mpp_global_sum(Dom%domain2d, tmpL(:,:), global_sum_flag)
        ttracer  = tracer_L + tracer_E
     endif

     tracer_error = dmasstracer-tracer_input_total-tracer_th_tend-tracer_src
!     if (use_blobs) tracer_error = tracer_error + dmasstracer_L2
     flux_error   = tracer_error/(dtime*Grd%tcellsurf + epsln)

     if (T_prog(n)%name =='temp') then

         write (stdoutunit,'(/a)') ' ----Single time step diagnostics for tracer '//trim(T_prog(n)%name)//'----' 
         write (stdoutunit,'(a,es24.17,a)') ' Total heat in ocean at (tau) referenced to 0degC   = ',&
                                          ttracer,' J'
         if (use_blobs) then
            write (stdoutunit,'(a,es24.17,a)') ' Total heat change in L system for (taup1-taum1)    = ',&
                                             dmasstracer_L,' J'
            write (stdoutunit,'(a,es24.17,a)') ' Total heat change in E system for (taup1-taum1)    = ',&
                                             dmasstracer_E,' J'
         endif
         write (stdoutunit,'(a,es24.17,a)') ' Total heat change of full system for (taup1-taum1) = ',&
                                          dmasstracer,' J'
         if (use_blobs) then
            write (stdoutunit,'(a,es24.17,a)') ' Total heat transfer from L to E system (tau)       = ',&
                                             dmasstracer_L2,' J'
            write (stdoutunit,'(a,es24.17,a)') ' Total heat in ocean L system at (tau) ref to 0degC = ',&
                                             tracer_L,' J'
            write (stdoutunit,'(a,es24.17,a)') ' Total heat in ocean E system at (tau) ref to 0degC = ',&
                                             tracer_E,' J'
         endif
         write (stdoutunit,'(a,es24.17,a)') ' Total heat input to ocean referenced to 0degC      = ',&
                                          tracer_input_total,' J'
         write (stdoutunit,'(a,es24.17,a)') ' Heat input via surface heat fluxes                 = ',&
                                          tracer_input_stf,' J'
         write (stdoutunit,'(a,es24.17,a)') ' Heat input via bottom heat fluxes                  = ',&
                                          tracer_input_btf,' J'
         write (stdoutunit,'(a,es24.17,a)') ' Heat input via open boundaries                     = ',&
                                          tracer_input_otf,' J'
         write (stdoutunit,'(a,es24.17,a)') ' Heat input via precip-evap+melt                    = ',&
                                          tracer_input_pme,' J'
         write (stdoutunit,'(a,es24.17,a)') ' Heat input via river runoff                        = ',&
                                          tracer_input_runoff,' J'
         write (stdoutunit,'(a,es24.17,a)') ' Heat input via calving land ice                    = ',&
                                          tracer_input_calving,' J'
         write (stdoutunit,'(a,es24.17,a)') ' Heat input via frazil formation                    = ',&
                                          frazil_heat,' J'
         write (stdoutunit,'(a,es24.17,a)') ' Heat input via sources in the source array         = ',&
                                          tracer_src,' J'
         write (stdoutunit,'(a,es24.17,a)') ' Heat input via sources in th_tendency, or errors   = ',&
                                          tracer_th_tend,' J'
         write (stdoutunit,'(a,es24.17,a)') ' Heat input via eta_t smooth                        = ',&
                                          tracer_input_eta_smooth,' J'
         write (stdoutunit,'(a,es24.17,a)') ' Heat input via pbot_t smooth                       = ',&
                                          tracer_input_pbot_smooth,' J'
         write (stdoutunit,'(a,es24.17,a)') ' d(T*rho*dV)                                        = ',&
                                          dmasstracer,' J'
         write (stdoutunit,'(a,es24.17,a)') ' Tracer mismatch: cp*d(rho*dV*T)-input              = ',&
                                          tracer_error,' J'
         write (stdoutunit,'(a,es24.17,a)') ' Mismatch converted to a surface flux               = ',&
                                          flux_error,' W/m^2'

     elseif (T_prog(n)%name(1:3) =='age') then

         tracer_error = dmasstracer-tracer_input_total-dtts*sec_in_yr_r*(mass_taup1-mass_taup1_keq1)
         tracer_error = tracer_error*mass_taup1_r
         flux_error   = sec_in_yr*tracer_error/(dtime*Grd%tcellsurf + epsln)

         write (stdoutunit,'(/a)') ' ----Single time step diagnostics for tracer '//trim(T_prog(n)%name)//'----' 
         write (stdoutunit,'(a,es24.17,a)') ' Total age at (tau) [sum(rho*dV*age)/mass_ocean]     = ',&     
                                          ttracer/mass_total,' yr'
         if (use_blobs) then
            if (mass_L .ne. 0.0) then
               write (stdoutunit,'(a,es24.17,a)') ' Total age in ocean L system at (tau)                = ',&
                                                  tracer_L/mass_L,' yr'
            else
               write (stdoutunit,'(a,es24.17,a)') ' Total age in ocean L system at (tau)                = ',&
                                                  0.0,' yr'
            endif
            write (stdoutunit,'(a,es24.17,a)') ' Total age in ocean E system at (tau)                = ',&
                                               tracer_E/mass_E,' yr'
         endif
         write (stdoutunit,'(a,es24.17,a)') ' dtts*(mass_ocean-mass_k=1)*sec_in_yr_r              = ',&     
                                          dtts*sec_in_yr_r*(mass_taup1-mass_taup1_keq1),' kg*yr'
         write (stdoutunit,'(a,es24.17,a)') ' dtts*(mass_ocean-mass_k=1)*sec_in_yr_r/mass_taup1   = ',&     
                             mass_taup1_r*dtts*sec_in_yr_r*(mass_taup1-mass_taup1_keq1),' yr'
         write (stdoutunit,'(a,es24.17,a)') ' Total age input to ocean                            = ',&
                                          tracer_input_total*mass_taup1_r,' yr'
         write (stdoutunit,'(a,es24.17,a)') ' Age input via surface fluxes                        = ',&
                                          tracer_input_stf*mass_taup1_r,' yr'
         write (stdoutunit,'(a,es24.17,a)') ' Age input via bottom fluxes                         = ',&
                                          tracer_input_btf*mass_taup1_r,' yr'
         write (stdoutunit,'(a,es24.17,a)') ' Age input via open boundaries                       = ',&
                                          tracer_input_otf*mass_taup1_r,' yr'
         write (stdoutunit,'(a,es24.17,a)') ' Age input via precip-evap+melt                      = ',&
                                          tracer_input_pme*mass_taup1_r,' yr'
         write (stdoutunit,'(a,es24.17,a)') ' Age input via river runoff                          = ',&
                                          tracer_input_runoff*mass_taup1_r,' yr'
         write (stdoutunit,'(a,es24.17,a)') ' Age input via calving land ice                      = ',&
                                          tracer_input_calving*mass_taup1_r,' yr'
         write (stdoutunit,'(a,es24.17,a)') ' Age input via sources in source array               = ',&
                                          tracer_src*mass_taup1_r,' yr'
         write (stdoutunit,'(a,es24.17,a)') ' Age input via sources in th_tendency, or errors     = ',&
                                          tracer_th_tend*mass_taup1_r,' yr'
         write (stdoutunit,'(a,es24.17,a)') ' Age input via eta_t smoother                        = ',&
                                          tracer_input_eta_smooth*mass_taup1_r,' yr'
         write (stdoutunit,'(a,es24.17,a)') ' Age input via pbot_t smoother                       = ',&
                                          tracer_input_pbot_smooth*mass_taup1_r,' yr'
         write (stdoutunit,'(a,es24.17,a)') ' d(age*rho*vol)/mass_ocean                           = ',&
                                          dmasstracer*mass_taup1_r,' yr'
         write (stdoutunit,'(a,es24.17,a)') ' Age mismatch: [d(rho*dV*age) - dtts*(Mt-M1)]/Mt     = ',&
                                          tracer_error,' yr'
         write (stdoutunit,'(a,es24.17,a)') ' Mismatch converted to a surface flux                = ',&
                                          flux_error,' 1/m^2'

     else 

         write (stdoutunit,'(/a)') '----Single time step diagnostics for tracer '//trim(T_prog(n)%name)//'----' 
         write (stdoutunit,'(a,es24.17,a)') ' Total tracer in ocean at time  (tau)                = ',&
                                          ttracer,' kg'
         if (use_blobs) then
            write (stdoutunit,'(a,es24.17,a)') ' Total tracer change in L system for (taup1-taum1)   = ',&
                                               dmasstracer_L,' kg'
            write (stdoutunit,'(a,es24.17,a)') ' Total tracer change in E system for (taup1-taum1)   = ',&
                                               dmasstracer_E,' kg'
         endif
         write (stdoutunit,'(a,es24.17,a)') ' Total tracer change in system for (taup1-taum1)     = ',&
                                          dmasstracer,' kg'
         if (use_blobs) then
            write (stdoutunit,'(a,es24.17,a)') ' Total tracer transfer from L to E system (tau)      = ',&
                                               dmasstracer_L2,' kg'
            write (stdoutunit,'(a,es24.17,a)') ' Total tracer in ocean L system at (tau)             = ',&
                                               tracer_L,' kg'
            write (stdoutunit,'(a,es24.17,a)') ' Total tracer in ocean E system at (tau)             = ',&
                                               tracer_E,' kg'
         endif
         write (stdoutunit,'(a,es24.17,a)') ' Total tracer input to ocean                         = ',&
                                          tracer_input_total,' kg'
         write (stdoutunit,'(a,es24.17,a)') ' Tracer input via surface fluxes                     = ',&
                                          tracer_input_stf,' kg'
         write (stdoutunit,'(a,es24.17,a)') ' Tracer input via bottom fluxes                      = ',&
                                          tracer_input_btf,' kg'
         write (stdoutunit,'(a,es24.17,a)') ' Tracer input via open boundaries                    = ',&
                                          tracer_input_otf,' kg'
         write (stdoutunit,'(a,es24.17,a)') ' Tracer input via precip-evap+melt                   = ',&
                                          tracer_input_pme,' kg'
         write (stdoutunit,'(a,es24.17,a)') ' Tracer input via river runoff                       = ',&
                                          tracer_input_runoff,' kg'
         write (stdoutunit,'(a,es24.17,a)') ' Tracer input via calving land ice                   = ',&
                                          tracer_input_calving,' kg'
         write (stdoutunit,'(a,es24.17,a)') ' Tracer input via sources in source array            = ',&
                                          tracer_src,' kg'
         write (stdoutunit,'(a,es24.17,a)') ' Tracer input via sources in th_tendency, or errors  = ',&
                                          tracer_th_tend,' kg'
         write (stdoutunit,'(a,es24.17,a)') ' Tracer input via eta_t smoother                     = ',&
                                          tracer_input_eta_smooth,' kg'
         write (stdoutunit,'(a,es24.17,a)') ' Tracer input via pbot_t smoother                    = ',&
                                          tracer_input_pbot_smooth,' kg'

         if (T_prog(n)%name =='salt') then

         write (stdoutunit,'(a,es24.17,a)') ' .001*d(S*rho*dV)                                    = ', &
                                          psu2ppt*dmasstracer,' kg'
        write (stdoutunit,'(a,es24.17,a)')  ' Tracer mismatch: 0.001*d(rho*dV*T)-input            = ', &
                                          psu2ppt*tracer_error,' kg '
         write (stdoutunit,'(a,es24.17,a)') ' Mismatch converted to a surface flux                = ', &
                                          psu2ppt*flux_error,' kg/(m^2 sec) '

         else  

         write (stdoutunit,'(a,es24.17,a)') ' d(T*rho*dV)                                         = ',dmasstracer,' kg'
         write (stdoutunit,'(a,es24.17,a)') ' Tracer mismatch: d(rho*dV*T)-input                  = ',tracer_error,' kg'
         write (stdoutunit,'(a,es24.17,a)') ' Mismatch converted to a surface flux                = ',flux_error,' kg/(m^2 sec) '

         endif

     endif

  enddo  ! end n-loop

  call mpp_clock_end(id_tracer_change)

end subroutine tracer_change
! </SUBROUTINE>  NAME="tracer_change"

     
!#######################################################################
! <FUNCTION NAME="total_tracer">
!
! <DESCRIPTION>
! Compute integrated tracer in model. 
! </DESCRIPTION>
!
function total_tracer (Tracer, Thickness, index)

  type(ocean_prog_tracer_type), intent(in) :: Tracer
  type(ocean_thickness_type),   intent(in) :: Thickness
  integer,                      intent(in) :: index 
  real                                     :: total_tracer

  real, dimension(isd:ied,jsd:jed) :: tracer_k
  integer                          :: k

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag_mod (total_tracer): module needs initialization ')
  endif 

  total_tracer  = 0.0
  tracer_k(:,:) = 0.0

  if (use_blobs) then
     do k=1,nk
        tracer_k(:,:) = tracer_k(:,:) + Tracer%conversion*Grd%tmask(:,:,k)             &
              *( Grd%dat(:,:)*Thickness%rho_dzt(:,:,k,index)*Tracer%field(:,:,k,index) &
              + Tracer%sum_blob(:,:,k,index))
     enddo
  else 
     do k=1,nk
        tracer_k(:,:) =  tracer_k(:,:) + Tracer%conversion*Grd%tmask(:,:,k)*Grd%dat(:,:) &
                        *Thickness%rho_dzt(:,:,k,index)*Tracer%field(:,:,k,index)
     enddo
  endif 

  if(have_obc) tracer_k(:,:) = tracer_k(:,:)*Grd%obc_tmask(:,:)
  total_tracer  = mpp_global_sum(Dom%domain2d,tracer_k(:,:), global_sum_flag)


end function total_tracer
! </FUNCTION>  NAME="total_tracer"


!#######################################################################
! <FUNCTION NAME="klevel_total_tracer">
!
! <DESCRIPTION>
! Compute integrated tracer on a k-level. 
! </DESCRIPTION>
!
function klevel_total_tracer(Tracer, Thickness, tindex, kindex)

  type(ocean_prog_tracer_type), intent(in) :: Tracer
  type(ocean_thickness_type),   intent(in) :: Thickness
  integer,                      intent(in) :: tindex 
  integer,                      intent(in) :: kindex
  real                                     :: klevel_total_tracer

  real, dimension(isd:ied,jsd:jed) :: tracer_k

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag_mod (total_tracer): module needs initialization ')
  endif 

  klevel_total_tracer=0.0

  if (use_blobs) then
     tracer_k(:,:) =  Tracer%conversion*Grd%tmask(:,:,kindex) &
          *( Grd%dat(:,:)*Thickness%rho_dzt(:,:,kindex,tindex)*Tracer%field(:,:,kindex,tindex) &
          + Tracer %sum_blob(:,:,kindex,tindex) )
  else
     tracer_k(:,:) =  Tracer%conversion*Grd%tmask(:,:,kindex)*Grd%dat(:,:) &
          *Thickness%rho_dzt(:,:,kindex,tindex)*Tracer%field(:,:,kindex,tindex)
  endif
  if(have_obc) tracer_k(:,:) = tracer_k(:,:)*Grd%obc_tmask(:,:)

  klevel_total_tracer  = mpp_global_sum(Dom%domain2d,tracer_k(:,:), global_sum_flag)

end function klevel_total_tracer
! </FUNCTION>  NAME="klevel_total_tracer"


!#######################################################################
! <FUNCTION NAME="total_mass">
!
! <DESCRIPTION>
! Compute total ocean tracer cell mass.  For Boussinesq fluid, 
! mass is determined using rho0 for density.
! </DESCRIPTION>
!
function total_mass (Thickness, index)

  type(ocean_thickness_type), intent(in) :: Thickness
  integer,                    intent(in) :: index

  real, dimension(isd:ied,jsd:jed) :: mass_t
  real                             :: total_mass
  integer                          :: k

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag_mod (total_mass): module needs initialization ')
  endif 

  total_mass  = 0.0
  mass_t(:,:) = 0.0

  if (use_blobs) then
     do k=1,nk
        mass_t(:,:) = mass_t(:,:) + Thickness%rho_dztT(:,:,k,index)*Grd%dat(:,:)*Grd%tmask(:,:,k)
     enddo
  else
     do k=1,nk
        mass_t(:,:) = mass_t(:,:) + Thickness%rho_dzt(:,:,k,index)*Grd%dat(:,:)*Grd%tmask(:,:,k)
     enddo
  endif
  if(have_obc) mass_t(:,:) = mass_t(:,:)*Grd%obc_tmask(:,:)
  total_mass  = mpp_global_sum(Dom%domain2d,mass_t(:,:), global_sum_flag)


end function total_mass
! </FUNCTION>  NAME="total_mass"


!#######################################################################
! <FUNCTION NAME="total_volume">
!
! <DESCRIPTION>
! Compute total ocean tracer cell volume. 
! </DESCRIPTION>
!
function total_volume (Thickness, index)

  type(ocean_thickness_type), intent(in) :: Thickness
  integer,                    intent(in) :: index

  real, dimension(isd:ied,jsd:jed) :: volume_t
  real                             :: total_volume
  integer                          :: k

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag_mod (total_volume): module needs initialization ')
  endif 

  total_volume  = 0.0
  volume_t(:,:) = 0.0

  if (use_blobs) then
     do k=1,nk
        volume_t(:,:) = volume_t(:,:) + Thickness%dztT(:,:,k,index)*Grd%dat(:,:)*Grd%tmask(:,:,k)
     enddo
  else
     do k=1,nk
        volume_t(:,:) = volume_t(:,:) + Thickness%dzt(:,:,k)*Grd%dat(:,:)*Grd%tmask(:,:,k)
     enddo
  endif
  if(have_obc) volume_t(:,:) = volume_t(:,:)*Grd%obc_tmask(:,:)

  total_volume  = mpp_global_sum(Dom%domain2d, volume_t(:,:), global_sum_flag)

end function total_volume
! </FUNCTION>  NAME="total_volume"


!#######################################################################
! <FUNCTION NAME="klevel_total_mass">
!
! <DESCRIPTION>
! Compute ocean tracer cell mass in a k-level.  For Boussinesq fluid, 
! mass is determined using rho0 for density.
! </DESCRIPTION>
!
function klevel_total_mass (Thickness, tindex, kindex)

  type(ocean_thickness_type), intent(in) :: Thickness
  integer,                    intent(in) :: tindex
  integer,                    intent(in) :: kindex

  real, dimension(isd:ied,jsd:jed) :: mass_t
  real                             :: klevel_total_mass

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag_mod (klevel_total_mass): module needs initialization ')
  endif 

  klevel_total_mass = 0.0
  if (use_blobs) then
     mass_t(:,:) = Thickness%rho_dztT(:,:,kindex,tindex)*Grd%dat(:,:)*Grd%tmask(:,:,kindex)
  else
     mass_t(:,:) = Thickness%rho_dzt(:,:,kindex,tindex)*Grd%dat(:,:)*Grd%tmask(:,:,kindex)
  endif
  if(have_obc) mass_t(:,:) = mass_t(:,:)*Grd%obc_tmask(:,:)
  klevel_total_mass = mpp_global_sum(Dom%domain2d,mass_t(:,:), global_sum_flag)

end function klevel_total_mass
! </FUNCTION>  NAME="klevel_total_mass"


!#######################################################################
! <SUBROUTINE NAME="tracer_integrals">
!
! <DESCRIPTION>
! Compute some integrated tracer diagnostics. 
! </DESCRIPTION>
!
subroutine tracer_integrals (Time, Thickness, Tracer)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  
  real, dimension(isd:ied,jsd:jed)   :: mass_t
  real, dimension(isd:ied,jsd:jed,3) :: T_field

  real    :: tbar, tvar, tchg, ttot
  real    :: mass_total 
  integer :: k, tau, taum1, taup1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag_mod (tracer_integrals): module needs initialization ')
  endif 

  call mpp_clock_begin(id_tracer_integrals)

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1
  
  write(stdoutunit,*) ' '
  call write_timestamp(Time%model_time)
  write(stdoutunit,'(/a,a)') 'Tracer integrals for ' , trim(Tracer%name)
  
  tbar = 0.0; tvar = 0.0; tchg = 0.0; ttot = 0.0

  mass_total = total_mass(Thickness, tau)

  wrk1_2d(:,:) = 0.0
  wrk2_2d(:,:) = 0.0
  wrk3_2d(:,:) = 0.0

  if (use_blobs) then
     do k=1,nk
        mass_t(:,:) = Thickness%rho_dzt(:,:,k,tau)*Grd%dat(:,:)*Grd%tmask(:,:,k)
        if(have_obc) mass_t(:,:) = mass_t(:,:) * Grd%obc_tmask(:,:)
        
        wrk1_2d(:,:)  = wrk1_2d(:,:) + Tracer%field(:,:,k,tau)*mass_t(:,:) &
                                     + Tracer%sum_blob(:,:,k,tau)*Grd%tmask(:,:,k)

        T_field(:,:,tau)  = ( Grd%datr(:,:)*Tracer%sum_blob(:,:,k,tau)              &
                             +Tracer%field(:,:,k,tau)*Thickness%rho_dzt(:,:,k,tau) )&
                             /(Thickness%rho_dztT(:,:,k,tau) + epsln)
        wrk2_2d(:,:)  = wrk2_2d(:,:) + T_field(:,:,tau)*T_field(:,:,tau)&
                                       *Thickness%rho_dztT(:,:,k,tau)*Grd%dat(:,:)

        T_field(:,:,taup1) = ( Grd%datr(:,:)*Tracer%sum_blob(:,:,k,taup1)                &
                              +Tracer%field(:,:,k,taup1)*Thickness%rho_dzt(:,:,k,taup1) )&
                              /(Thickness%rho_dztT(:,:,k,taup1) + epsln)
        T_field(:,:,taum1) = ( Grd%datr(:,:)*Tracer%sum_blob(:,:,k,taum1)                &
                              +Tracer%field(:,:,k,taum1)*Thickness%rho_dzt(:,:,k,taum1) )&
                              /(Thickness%rho_dztT(:,:,k,taum1) + epsln)    
        wrk3_2d(:,:)  = wrk3_2d(:,:) + abs(T_field(:,:,taup1)-T_field(:,:,taum1))*dtimer

     enddo

  else
     do k=1,nk
        mass_t(:,:) = Thickness%rho_dzt(:,:,k,tau)*Grd%dat(:,:)*Grd%tmask(:,:,k)
        if(have_obc) mass_t(:,:) = mass_t(:,:) * Grd%obc_tmask(:,:)
        wrk1_2d(:,:) = wrk1_2d(:,:) + Tracer%field(:,:,k,tau)*mass_t(:,:)
        wrk2_2d(:,:) = wrk2_2d(:,:) + Tracer%field(:,:,k,tau)*Tracer%field(:,:,k,tau)*mass_t(:,:)
        wrk3_2d(:,:) = wrk3_2d(:,:) + abs(Tracer%field(:,:,k,taup1)-Tracer%field(:,:,k,taum1))*dtimer
     enddo
  endif
  tbar = mpp_global_sum(Dom%domain2d, wrk1_2d(:,:), global_sum_flag)
  tvar = mpp_global_sum(Dom%domain2d, wrk2_2d(:,:), global_sum_flag)
  tchg = mpp_global_sum(Dom%domain2d, wrk3_2d(:,:), global_sum_flag)


  ttot = total_tracer(Tracer,Thickness,tau)
  tbar = tbar/mass_total
  tvar = tvar/mass_total - tbar**2

  write(stdoutunit, '(/a,a,a,es24.17)') 'Average  ',trim(Tracer%name),' = ', tbar
  write(stdoutunit, '(a,a,a,es24.17)')  'Variance ',trim(Tracer%name),' = ', tvar
  write(stdoutunit, '(a,a,a,es24.17)')  '|dT/dt|  ',trim(Tracer%name),' = ', tchg
  write(stdoutunit, '(a,a,a,es24.17/)') 'Total    ',trim(Tracer%name),' = ', ttot

  call mpp_clock_end(id_tracer_integrals)
  
end subroutine tracer_integrals
! </SUBROUTINE>  NAME="tracer_integrals"



!#######################################################################
! <SUBROUTINE NAME="tracer_land_cell_check">
!
! <DESCRIPTION>
! Check to be sure ocean tracer is zero over land 
! </DESCRIPTION>
!
subroutine tracer_land_cell_check (Time, T_prog)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  integer :: i, j, k, n, num, taup1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  call mpp_clock_begin(id_tracer_land_cell_check)

  if(size(T_prog,1) /= num_prog_tracers) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag_mod (tracer_land_cell_check): passing T_prog with wrong dimensions')
  endif 

  if (.not. module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag_mod (tracer_land_cell_check): module needs initialization ')
  endif 

  taup1 = Time%taup1

  write (stdoutunit,'(1x,/a/)')'Locations (if any) where land cell tracer is non-zero...'
  num = 0
  do n=1,num_prog_tracers 
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              if (Grd%tmask(i,j,k) == 0 .and. (T_prog(n)%field(i,j,k,taup1) /= 0.0) ) then
                  num = num + 1
                  if (num < 10) then
                    write(unit,9100) &
                    i+Dom%ioff,j+Dom%joff,k,Grd%xt(i,j),Grd%yt(i,j),Grd%zt(k),n,T_prog(n)%field(i,j,k,taup1)
                  endif 
              endif
           enddo
        enddo
     enddo
  enddo

  call mpp_max(num)
  if (num > 0) call mpp_error(FATAL)

9100  format(/' " =>Error: Land cell at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f8.3,',',f8.3,',',f12.3,'m) has tracer(n=',i2,') =',e10.3)

  call mpp_clock_end(id_tracer_land_cell_check)


end subroutine tracer_land_cell_check
! </SUBROUTINE>  NAME="tracer_land_cell_check"


!#######################################################################
! <SUBROUTINE NAME="mass_conservation">
!
! <DESCRIPTION>
! Compute change in mass over many time steps, and compare to the 
! input of mass through surface to check for mass conservation.
!
!============================================================
!
! threelevel scheme 
! 
! Here is the logic for the accumulation of the fluxes and 
! comparisons between mass/volumes at the start and the end. 
!
! Consider accumulation over four leap-frog time steps. 
! Ignore time filtering.  
!
! mass(2) = mass(0) + 2dt*F(1)  taup1=2, taum1=0, tau=1 
!
! mass(3) = mass(1) + 2dt*F(2)  taup1=3, taum1=1, tau=2 
!
! mass(4) = mass(2) + 2dt*F(3)  taup1=4, taum1=2, tau=3 
!
! mass(5) = mass(3) + 2dt*F(4)  taup1=5, taum1=3, tau=4 
! 
! Hence,
!
! [mass(4) + mass(5)] = [mass(0) + mass(1)] + 2dt*[F(1)+F(2)+F(3)+F(4)]
! 
! For this example, we have 
!
! itts_mass=1 through itte_mass=4 for accumulating fluxes
!
! itt=itts_mass=1=tau we use taum1=0 and tau=1 to get starting mass
!
! itt=itte_mass=4=tau we use tau=4 and taup1=5 to get the final mass
!
!============================================================
!
! twolevel scheme
! 
! Here is the logic for the accumulation of the fluxes and 
! comparisons between mass/volumes at the start and the end. 
!
! Consider accumulation over four time steps. 
!
! mass(3/2) = mass(1/2) + dt*F(1)   taup1=3/2, taum1=1/2, tau=1 
!
! mass(5/2) = mass(3/2) + dt*F(2)   taup1=5/2, taum1=3/2, tau=2 
!
! mass(7/2) = mass(5/2) + dt*F(3)   taup1=7/2, taum1=5/2, tau=3 
!
! mass(9/2) = mass(7/2) + dt*F(4)   taup1=9/2, taum1=7/2, tau=4 
! 
! Hence,
!
! mass(9/2) = mass(1/2) + dt*[F(1)+F(2)+F(3)+F(4)]
! 
! For this example, we have 
!
! itts_mass=1 through itte_mass=4 for accumulating fluxes
!
! itt=itts_mass=1=tau we use taum1=1/2 to get starting mass
!
! itt=itte_mass=4=tau we use taup1=9/2 to get the final mass
!
! </DESCRIPTION>
!
subroutine mass_conservation (Time, Thickness, Ext_mode, pme, runoff, calving)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_thickness_type),     intent(in) :: Thickness
  type(ocean_external_mode_type), intent(in) :: Ext_mode
  real, dimension(isd:,jsd:),     intent(in) :: pme
  real, dimension(isd:,jsd:),     intent(in) :: runoff
  real, dimension(isd:,jsd:),     intent(in) :: calving

  real :: total_time
  real :: mass_error, mass_error_rate
  real :: mass_input, mass_eta_smooth, mass_pbot_smooth, mass_chg
  real :: pme_input, runoff_input, calving_input, source_input

  integer :: i, j
  integer :: itt, taum1, tau, taup1

  logical, save :: first=.true.  

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag (mass_conservation): module needs initialization ')
  endif 

  itt   = Time%itt
  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  if (first) then
    first = .false.
    allocate (pme_flux(isd:ied,jsd:jed))
    allocate (runoff_flux(isd:ied,jsd:jed))
    allocate (source_flux(isd:ied,jsd:jed))
    allocate (calving_flux(isd:ied,jsd:jed))
    allocate (eta_smooth(isd:ied,jsd:jed))
    allocate (pbot_smooth(isd:ied,jsd:jed))
    pme_flux     = 0.0
    runoff_flux   = 0.0
    source_flux  = 0.0
    calving_flux = 0.0
    eta_smooth   = 0.0
    pbot_smooth  = 0.0
    mass_start   = 0.0
    mass_final   = 0.0
    itts_mass    = itt+2
    itte_mass    = itts_mass-2 + nint(tracer_conserve_days*86400.0/dtts)
    if(itts_mass  >= itte_mass) then 
      write(stdoutunit,*) &
      '==>Warning: increase tracer_conserve_days to properly compute mass and tracer conservation.'
    endif 
  endif

  if(tendency==THREE_LEVEL) then 
      if (itt == itts_mass) then
          mass_start = 0.5*(total_mass(Thickness, taum1)+total_mass(Thickness, tau))    
      endif
      if (itt == itte_mass) then
          mass_final = 0.5*(total_mass(Thickness, tau)+total_mass(Thickness, taup1))
      endif
  elseif(tendency==TWO_LEVEL) then 
      if (itt == itts_mass) then
          mass_start = total_mass(Thickness, taum1)
      endif
      if (itt == itte_mass) then
          mass_final = total_mass(Thickness, taup1)
      endif
  endif

  if(itt >= itts_mass .and. itt <= itte_mass) then 
      do j=jsc,jec
         do i=isc,iec
            pme_flux(i,j)     = pme_flux(i,j)    + dteta*Grd%dat(i,j)*pme(i,j)
            runoff_flux(i,j)  = runoff_flux(i,j) + dteta*Grd%dat(i,j)*runoff(i,j)
            calving_flux(i,j) = calving_flux(i,j)+ dteta*Grd%dat(i,j)*calving(i,j)
            source_flux(i,j)  = source_flux(i,j) + dteta*Grd%dat(i,j)*Ext_mode%source(i,j)
            eta_smooth(i,j)   = eta_smooth(i,j)  + dteta*Grd%tmask(i,j,1)*Grd%dat(i,j)*Ext_mode%eta_smooth(i,j)
            pbot_smooth(i,j)  = pbot_smooth(i,j) + dteta*Grd%tmask(i,j,1)*Grd%dat(i,j)*Ext_mode%pbot_smooth(i,j)
         enddo
      enddo
  endif

  if (itt==itte_mass) then 

      if(have_obc) then
         pme_flux(:,:)     = pme_flux(:,:)*Grd%obc_tmask(:,:)
         runoff_flux(:,:)  = runoff_flux(:,:)*Grd%obc_tmask(:,:)
         calving_flux(:,:) = calving_flux(:,:)*Grd%obc_tmask(:,:)
         source_flux(:,:)  = source_flux(:,:)*Grd%obc_tmask(:,:)
         eta_smooth(:,:)   = eta_smooth(:,:)*Grd%obc_tmask(:,:)
         pbot_smooth(:,:)  = pbot_smooth(:,:)*Grd%obc_tmask(:,:)
         eta_smooth(:,:)   = eta_smooth(:,:)*Grd%obc_tmask(:,:)
         pbot_smooth(:,:)  = pbot_smooth(:,:)*Grd%obc_tmask(:,:)
      endif
      pme_input        = mpp_global_sum(Dom%domain2d,pme_flux(:,:),global_sum_flag)
      runoff_input     = mpp_global_sum(Dom%domain2d,runoff_flux(:,:),global_sum_flag)
      calving_input    = mpp_global_sum(Dom%domain2d,calving_flux(:,:),global_sum_flag)
      source_input     = mpp_global_sum(Dom%domain2d,source_flux(:,:),global_sum_flag)
      mass_eta_smooth  = mpp_global_sum(Dom%domain2d,eta_smooth(:,:), global_sum_flag)
      mass_pbot_smooth = mpp_global_sum(Dom%domain2d,pbot_smooth(:,:),global_sum_flag)
      mass_input       = pme_input + runoff_input + calving_input + source_input + mass_eta_smooth + mass_pbot_smooth

      total_time      = (itte_mass-itts_mass+1)*dtts+epsln
      mass_chg        = mass_final - mass_start
      mass_error      = mass_chg-mass_input
      mass_error_rate = mass_error/total_time 

      write (stdoutunit,'(/50x,a/)')           &
       ' Measures of global integrated mass conservation over multiple time steps'
      write (stdoutunit,'(a,i10,a,es24.17,a)') &
       ' Ocean mass at timestep ',itts_mass,'                = ',mass_start,' kg.'
      write (stdoutunit,'(a,i10,a,es24.17,a)') &
       ' Ocean mass at timestep ',itte_mass,'                = ',mass_final,' kg.'
      write (stdoutunit,'(a,es24.17,a)') &
      ' Change in ocean mass over time interval              = ',mass_chg,' kg.'
      write (stdoutunit,'(a,es24.17,a)') &
      ' Mass input via all contributions over time interval  = ',mass_input,' kg.'
      write (stdoutunit,'(a,es24.17,a)')  &
      ' Error in mass content change                         = ',mass_error,' kg.'
      write (stdoutunit,'(a,es24.17,a)')  &
      ' Error in rate of mass content change                 = ',mass_error_rate,' kg/s.'

      write (stdoutunit,'(a,es24.17,a)') &
      ' Mass input via P-E fluxes over time interval         = ',pme_input,' kg.'
      write (stdoutunit,'(a,es24.17,a)') &
      ' Mass input via runoff fluxes over time interval      = ',runoff_input,' kg.'
      write (stdoutunit,'(a,es24.17,a)') &
      ' Mass input via ice calving fluxes over time interval = ',calving_input,' kg.'
      write (stdoutunit,'(a,es24.17,a)') &
      ' Mass input via source fluxes over time interval      = ',source_input,' kg.'
      write (stdoutunit,'(a,es24.17,a)') &
      ' Mass input via eta_t smoother over time interval     = ',mass_eta_smooth,' kg.'
      write (stdoutunit,'(a,es24.17,a)') &
      ' Mass input via pbot_t smoother over time interval    = ',mass_pbot_smooth,' kg.'
      write (stdoutunit,'(/)')

  endif

end subroutine mass_conservation
! </SUBROUTINE>  NAME="mass_conservation"



!#######################################################################
! <SUBROUTINE NAME="tracer_conservation">
!
! <DESCRIPTION>
! Compute change in global integrated tracer over many time steps,
! and compare to the input of tracer through the boundaries to 
! check for total tracer conservation.
!
! Accumulate fluxes as in the mass_conservation diagnostic. 
!
! </DESCRIPTION>
!
subroutine tracer_conservation (Time, Thickness, T_prog, T_diag, pme, runoff, calving)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_thickness_type),     intent(in) :: Thickness
  type(ocean_prog_tracer_type),   intent(in) :: T_prog(:)
  type(ocean_diag_tracer_type),   intent(in) :: T_diag(:)
  real, dimension(isd:,jsd:),     intent(in) :: pme
  real, dimension(isd:,jsd:),     intent(in) :: runoff
  real, dimension(isd:,jsd:),     intent(in) :: calving
 
  real, dimension(isd:ied,jsd:jed) :: tmp
  real, dimension(isd:ied,jsd:jed) :: dt_con_area
  real, dimension(isd:ied,jsd:jed) :: runoff_mod
  real :: temporary
  real :: tracer_chg
  real :: tracer_error
  real :: tracer_error_rate
  real :: tracer_tend_input
  real :: tracer_calving_input
  real :: tracer_runoff_input
  real :: tracer_pme_input
  real :: tracer_stf_input
  real :: tracer_btf_input
  real :: tracer_otf_input
  real :: tracer_frazil_input
  real :: tracer_total_input
  real :: tracer_source_input
  real :: tracer_eta_smooth_input
  real :: tracer_pbot_smooth_input
  real :: total_time

  integer :: i, j, k, n
  integer :: itt
  integer :: tau, taum1, taup1

  logical, save :: first =.true.  

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag (tracer_conservation): module needs initialization ')
  endif 

  if(size(T_prog,1) /= num_prog_tracers) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag (tracer_conservation): passing T_prog with wrong dimensions')
  endif 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1
  itt   = Time%itt

  if (first) then
    first  = .false.
    allocate (tracer_tend(isd:ied,jsd:jed,num_prog_tracers))
    allocate (tracer_runoff(isd:ied,jsd:jed,num_prog_tracers))
    allocate (tracer_calving(isd:ied,jsd:jed,num_prog_tracers))
    allocate (tracer_pme(isd:ied,jsd:jed,num_prog_tracers))
    allocate (tracer_otf(isd:ied,jsd:jed,num_prog_tracers))
    allocate (tracer_stf(isd:ied,jsd:jed,num_prog_tracers))
    allocate (tracer_btf(isd:ied,jsd:jed,num_prog_tracers))
    allocate (tracer_frazil(isd:ied,jsd:jed))
    allocate (tracer_source(isd:ied,jsd:jed,num_prog_tracers))
    allocate (tracer_eta_smooth(isd:ied,jsd:jed,num_prog_tracers))
    allocate (tracer_pbot_smooth(isd:ied,jsd:jed,num_prog_tracers))
    allocate (tracer_start(num_prog_tracers))
    allocate (tracer_final(num_prog_tracers))
    tracer_source(:,:,:)      = 0.0
    tracer_tend(:,:,:)        = 0.0
    tracer_runoff(:,:,:)      = 0.0
    tracer_calving(:,:,:)     = 0.0
    tracer_pme(:,:,:)         = 0.0
    tracer_otf(:,:,:)         = 0.0
    tracer_stf(:,:,:)         = 0.0
    tracer_btf(:,:,:)         = 0.0
    tracer_frazil(:,:)        = 0.0
    tracer_source(:,:,:)      = 0.0
    tracer_eta_smooth(:,:,:)  = 0.0
    tracer_pbot_smooth(:,:,:) = 0.0
    tracer_start(:)           = 0.0
    tracer_final(:)           = 0.0
    itts_tracer               = itt+4
    itte_tracer               = itts_tracer-4 + nint(tracer_conserve_days*86400.0/dtts)
  endif

  if(tendency==THREE_LEVEL) then 
      if (itt == itts_tracer) then
          do n=1,num_prog_tracers
             tracer_start(n) = 0.5*(  total_tracer(T_prog(n),Thickness,taum1) &
                                    + total_tracer(T_prog(n),Thickness,tau))
          enddo
      endif
      if (itt == itte_tracer) then
          do n=1,num_prog_tracers
             tracer_final(n) = 0.5*(  total_tracer(T_prog(n),Thickness,tau)  &
                                    + total_tracer(T_prog(n),Thickness,taup1))      
          enddo
      endif
  elseif(tendency==TWO_LEVEL) then 
      if (itt == itts_tracer) then
          do n=1,num_prog_tracers
             tracer_start(n) = total_tracer(T_prog(n),Thickness,taum1)    
          enddo
      endif
      if (itt == itte_tracer) then
          do n=1,num_prog_tracers
             tracer_final(n) = total_tracer(T_prog(n),Thickness,taup1)
          enddo
      endif
  endif

  if(itt >= itts_tracer .and. itt <= itte_tracer) then 

      ! runoff can be negative, as for the "geoengineering" trick of siphoning
      ! Arabian Sea water to the Caspian Sea, to keep the Caspian from evaporating
      ! away during a coupled climate simulation.  As negative runoff is ignored  
      ! in the rivermix module (rightly so, as we do not wish to insert tracer 
      ! from negative runoff into the ocean), we also wish to ignore negative 
      ! runoff for the diagnostics here, in order to reflect what the prognostic
      ! equations are doing.      
      runoff_mod(:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec
            if(runoff(i,j) > 0) then        
                runoff_mod(i,j) = runoff(i,j)
            endif
         enddo
      enddo

      do n=1,num_prog_tracers

         ! global integral of th_tendency leaves just the sources 
         ! (and other stuff to be removed below)

         ! define an intermediate array 
         do j=jsc,jec
            do i=isc,iec
               dt_con_area(i,j) = dtts*T_prog(n)%conversion*Grd%dat(i,j)
            enddo
         enddo

         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  tracer_tend(i,j,n) = tracer_tend(i,j,n) +     &
                    dt_con_area(i,j)*Grd%tmask(i,j,k)*T_prog(n)%th_tendency(i,j,k)
                  tracer_source(i,j,n) = tracer_source(i,j,n) + &
                    dt_con_area(i,j)*Grd%tmask(i,j,k)*Thickness%rho_dzt(i,j,k,taup1)*T_prog(n)%source(i,j,k)
               enddo
            enddo
         enddo

         ! btf < 0 means heat is added to ocean; hence the minus sign 
         do j=jsc,jec
            do i=isc,iec        
               temporary                 = dt_con_area(i,j)*Grd%tmask(i,j,1)
               tracer_eta_smooth(i,j,n)  = tracer_eta_smooth(i,j,n)  + temporary*T_prog(n)%eta_smooth(i,j)
               tracer_pbot_smooth(i,j,n) = tracer_pbot_smooth(i,j,n) + temporary*T_prog(n)%pbot_smooth(i,j)
               tracer_runoff(i,j,n)      = tracer_runoff(i,j,n)      + temporary*T_prog(n)%trunoff(i,j)*runoff_mod(i,j)
               tracer_calving(i,j,n)     = tracer_calving(i,j,n)     + temporary*T_prog(n)%tcalving(i,j)*calving(i,j)
               tracer_pme(i,j,n)         = tracer_pme(i,j,n)         + temporary*T_prog(n)%tpme(i,j)*pme(i,j)
               tracer_stf(i,j,n)         = tracer_stf(i,j,n)         + temporary*T_prog(n)%stf(i,j)
               tracer_btf(i,j,n)         = tracer_btf(i,j,n)         - temporary*T_prog(n)%btf(i,j)
            enddo
         enddo
         if(n==index_temp .and. index_frazil > 0) then 
            do k=1,nk
               do j=jsc,jec
                   do i=isc,iec
                      tracer_frazil(i,j) = tracer_frazil(i,j)  &
                         + T_diag(index_frazil)%field(i,j,k)*Grd%dat(i,j)*Grd%tmask(i,j,k)
                   enddo
               enddo
            enddo
         endif

         if (have_obc) then
             call ocean_obc_tracer_flux(T_prog(n), tmp)
             tracer_otf(:,:,n)         = (tracer_otf(:,:,n) + tmp(:,:))*Grd%obc_tmask(:,:)
             tracer_stf(:,:,n)         = tracer_stf(:,:,n)*Grd%obc_tmask(:,:)
             tracer_btf(:,:,n)         = tracer_btf(:,:,n)*Grd%obc_tmask(:,:)
             tracer_runoff(:,:,n)      = tracer_runoff(:,:,n)*Grd%obc_tmask(:,:)
             tracer_calving(:,:,n)     = tracer_calving(:,:,n)*Grd%obc_tmask(:,:)
             tracer_pme(:,:,n)         = tracer_pme(:,:,n)*Grd%obc_tmask(:,:)
             tracer_tend(:,:,n)        = tracer_tend(:,:,n)*Grd%obc_tmask(:,:)
             tracer_source(:,:,n)      = tracer_source(:,:,n)*Grd%obc_tmask(:,:)
             tracer_eta_smooth(:,:,n)  = tracer_eta_smooth(:,:,n)*Grd%obc_tmask(:,:)
             tracer_pbot_smooth(:,:,n) = tracer_pbot_smooth(:,:,n)*Grd%obc_tmask(:,:)
             if (n==index_temp) then 
                tracer_frazil(:,:)     = tracer_frazil(:,:)*Grd%obc_tmask(:,:)
             endif
         endif

      enddo

  endif

  if (itt==itte_tracer) then 

    write (stdoutunit,'(/50x,a/)') &
    ' Measures of global integrated tracer conservation over multiple time steps'

    total_time = (itte_mass-itts_mass+1)*dtts+epsln

    do n=1,num_prog_tracers

      tracer_frazil_input = 0.0
      if (n==index_temp) then 
         tracer_frazil_input   = mpp_global_sum(Dom%domain2d,tracer_frazil(:,:),        global_sum_flag)
      endif 
      tracer_runoff_input      = mpp_global_sum(Dom%domain2d,tracer_runoff(:,:,n),      global_sum_flag)
      tracer_calving_input     = mpp_global_sum(Dom%domain2d,tracer_calving(:,:,n),     global_sum_flag)
      tracer_pme_input         = mpp_global_sum(Dom%domain2d,tracer_pme(:,:,n),         global_sum_flag)
      tracer_otf_input         = mpp_global_sum(Dom%domain2d,tracer_otf(:,:,n),         global_sum_flag)
      tracer_stf_input         = mpp_global_sum(Dom%domain2d,tracer_stf(:,:,n),         global_sum_flag)
      tracer_btf_input         = mpp_global_sum(Dom%domain2d,tracer_btf(:,:,n),         global_sum_flag)
      tracer_tend_input        = mpp_global_sum(Dom%domain2d,tracer_tend(:,:,n),        global_sum_flag)
      tracer_source_input      = mpp_global_sum(Dom%domain2d,tracer_source(:,:,n),      global_sum_flag)
      tracer_eta_smooth_input  = mpp_global_sum(Dom%domain2d,tracer_eta_smooth(:,:,n),  global_sum_flag)
      tracer_pbot_smooth_input = mpp_global_sum(Dom%domain2d,tracer_pbot_smooth(:,:,n), global_sum_flag)

      tracer_total_input =   tracer_stf_input    + tracer_btf_input        + tracer_otf_input &
                           + tracer_runoff_input + tracer_calving_input    + tracer_pme_input &
                           + tracer_frazil_input + tracer_eta_smooth_input + tracer_pbot_smooth_input

      ! runoff, calving, pme, and smooth are added to T_prog(n)%th_tendency inside 
      ! ocean_rivermix_mod (runoff, calving) and ocean_tracer_mod (pme, and smooth).
      ! To avoid double counting their contribution, we remove these terms from 
      ! tracer_tend_input. As a result, tracer_tend_input will only include 
      ! contributions from the global integration of internal flux convergences. 
      ! When integrated globally, these convergences should sum to zero, and such is 
      ! an important check that the contributions to tracer updates are coded properly.  
      tracer_tend_input =  tracer_tend_input - tracer_runoff_input     - tracer_calving_input & 
                          -tracer_pme_input  - tracer_eta_smooth_input - tracer_pbot_smooth_input

      ! stf and btf added to th_tendency when use explicit vertical diffusion (aidif=0)
      if (aidif==0.0) tracer_tend_input = tracer_tend_input - tracer_stf_input - tracer_btf_input

      tracer_chg   = tracer_final(n)-tracer_start(n)
      tracer_error = tracer_chg - tracer_total_input - tracer_tend_input - tracer_source_input
      tracer_error_rate = tracer_error/(total_time*Grd%tcellsurf) 

      if (n==index_temp) then 

        write (stdoutunit,'(a,i10,a,es24.17,a)') ' Ocean heat content at timestep ',itts_tracer,'    = ',tracer_start(n),' J.'
        write (stdoutunit,'(a,i10,a,es24.17,a)') ' Ocean heat content at timestep ',itte_tracer,'    = ',tracer_final(n),' J.'
        write (stdoutunit,'(a,es24.17,a)') ' Change in ocean heat content                         = ',tracer_chg,' J.'
        write (stdoutunit,'(a,es24.17,a)') ' Total Heat input                                     = ',tracer_total_input,' J.'
        write (stdoutunit,'(a,es24.17,a)') ' Error in heat content change                         = ',tracer_error,' J.'
        write (stdoutunit,'(a,es24.17,a)') ' Error in rate of heat content change                 = ',tracer_error_rate,' W/m^2.'
        write (stdoutunit,'(a,es24.17,a)') ' Heat input by tendency terms                         = ',tracer_tend_input,' J.'
        write (stdoutunit,'(a,es24.17,a)') ' Heat input by surface fluxes                         = ',tracer_stf_input,' J.'
        write (stdoutunit,'(a,es24.17,a)') ' Heat input by bottom fluxes                          = ',tracer_btf_input,' J.'
        write (stdoutunit,'(a,es24.17,a)') ' Heat input by runoff and calving                     = ',tracer_runoff_input,' J.'
        write (stdoutunit,'(a,es24.17,a)') ' Heat input by precip-evap+calving                    = ',tracer_pme_input,' J.'
        write (stdoutunit,'(a,es24.17,a)') ' Heat input by open boundaries                        = ',tracer_otf_input,' J.'
        write (stdoutunit,'(a,es24.17,a)') ' Heat input by frazil formation                       = ',tracer_frazil_input,' J.'
        write (stdoutunit,'(a,es24.17,a)') ' Heat input by eta_t smoother                         = ',tracer_eta_smooth_input,' J.'
        write (stdoutunit,'(a,es24.17,a)') ' Heat input by pbot_t smoother                        = ',tracer_pbot_smooth_input,' J.'
        write (stdoutunit,'(a,es24.17,a)') ' Heat input by sources                                = ',tracer_source_input,' J.'
        write (stdoutunit,'(/)')

      elseif(T_prog(n)%name(1:3)=='age') then 

        write (stdoutunit,'(a,i10,a,es24.17,a)') ' Mass weighted '//trim(T_prog(n)%name)// &
                         ' at timestep ',itts_tracer,'  = ',tracer_start(n),' yr*kg.'
        write (stdoutunit,'(a,i10,a,es24.17,a)') ' Mass weighted '//trim(T_prog(n)%name)// &
                         ' at timestep ',itte_tracer,'  = ',tracer_final(n),' yr*kg.'

        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' total tracer input                    = ',tracer_total_input,' kg*yr.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' change in ocean tracer                = ',tracer_chg,' kg*yr.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' error in tracer change                = ',tracer_error,' kg*yr.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' error in rate of tracer change        = ',sec_in_yr*tracer_error_rate,' kg/(m^2).'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by tendency terms               = ',tracer_tend_input,' kg*yr.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by surface fluxes               = ',tracer_stf_input,' kg*yr.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by bottom fluxes                = ',tracer_btf_input,' kg*yr.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by runoff and calving           = ',tracer_runoff_input,' kg*yr.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by precip-evap+calving          = ',tracer_pme_input,' kg*yr.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by open boundaries              = ',tracer_otf_input,' kg*yr.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by eta_t smoother               = ',tracer_eta_smooth_input,' kg*yr.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by pbot_t smoother              = ',tracer_pbot_smooth_input,' kg*yr.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by sources                      = ',tracer_source_input,' kg*yr.'
        write (stdoutunit,'(/)')

      else 

        write (stdoutunit,'(a,i10,a,es24.17,a)') ' Total '//trim(T_prog(n)%name)// &
                         ' mass at timestep ',itts_tracer,'  = ',tracer_start(n),' kg'
        write (stdoutunit,'(a,i10,a,es24.17,a)') ' Total '//trim(T_prog(n)%name)// &
                         ' mass at timestep ',itte_tracer,'  = ',tracer_final(n),' kg'

        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' change in ocean tracer mass                      = ',tracer_chg,' kg.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' error in mass change                             = ',tracer_error,' kg.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' error in rate of mass change                     = ',tracer_error_rate,' kg/(m^2 sec).'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// & 
                         ' input by tendency terms                          = ',tracer_tend_input,' kg.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by surface fluxes                          = ',tracer_stf_input,' kg.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by bottom fluxes                           = ',tracer_btf_input,' kg*yr.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by runoff and calving                      = ',tracer_runoff_input,' kg.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by precip-evap+calving                     = ',tracer_pme_input,' kg.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by open boundaries                         = ',tracer_otf_input,' kg*yr.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by eta_t smoother                          = ',tracer_eta_smooth_input,' kg.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by pbot_t smoother                         = ',tracer_pbot_smooth_input,' kg.'
        write (stdoutunit,'(1x,a,es24.17,a)') trim(T_prog(n)%name)// & 
                         ' input by sources                                 = ',tracer_source_input,' kg.'
        write (stdoutunit,'(/)')

      endif 

    enddo 
  endif

end subroutine tracer_conservation
! </SUBROUTINE>  NAME="tracer_conservation"


  
!#######################################################################
! <SUBROUTINE NAME="diagnose_kappa_sort">
!
! <DESCRIPTION>
! Routine to diagnose the amount of mixing between classes of a 
! particular tracer.  Temperature is used as default.
! Method follows that used in the paper 
!
! Spurious diapycnal mixing associated with advection in a
! z-coordinate ocean model, 2000: S.M. Griffies, R.C.
! Pacanowski, and R.W. Hallberg. Monthly Weather Review, vol 128, 538--564.
!
! This diagnostic is most useful when computing the levels of 
! effective dia-tracer mixing occuring in a model run with 
! zero buoyancy forcing at the boundaries.
!
! Algorithm notes:
!
! -assumes flat ocean bottom--non-flat bottoms loose the precise relation 
!  between sorted depth and true ocean depth.  This is a minor inconvenience.
!  The code is actually written so that the horizontal area of each 
!  layer can be different.  This will allow for this diagnostic to be 
!  used, say, for simple topography, such as bowl or bump.  
!
! -assumes Boussinesq fluid so that consider volume instead of mass of a cell.
!  Also his means that dzt = rho0r*rho_dzt
!
! -assumes area integrated eta_t is zero, so domain volume is static.
!  This is the case when there are no water boundary fluxes. 
! 
! -Results are meaningful only when dxt*dyt*dst of each grid cell
!  is the same.  This is best realized with a beta-plane or f-plane
!  geometry, and with zstar vertical coordinate. 
!
!  My best understanding of this limitation is related to systematic
!  biases in roundoff errors that result when the grid cells have 
!  varying volumes.  
!
!  If choose to use geopotential vertical coordinate, it is best
!  to set linear_free_surface=.true. in ocean_thickness_nml,
!  so that Thickness%rho_dzt = rho0%Grd%dzt.  The sorting model of 
!  mixing has not been generalized to evolving layer thicknesses
!  with geopotential. 
!
!  With zstar, the dst is constant in time, and the sorting method 
!  will sort to a depth in zstar space rather than geopotential space. 
!  This is a trivial distinction in principle, but should help with
!  some roundoff issues in practice. 
!
! -assumes tendency=TWO_LEVEL, which is exploited here to save memory.
!
! Numerical roundoff is a real issue with this diagnostic.
! It is critical that full double precision be used 
! to garner sensible results. 
!
! -defines some global arrays, so requires large memory.
!  This feature can be removed if parallel sort is 
!  implemented. So far, such has not been done.  
!
! -Effective kappa is set to zero at top of top-most cell.
! It is then diagnosed as zero (or roundoff) at bottom 
! of the column if there are zero boundary buoyancy fluxes.  
!
! -when computing density, we do rho=-alpha*theta.
! We drop the rho0 factor in order to reduce roundoff.
! Likewise, we assume alpha=1.0 rather than alpha=alpha_linear_eos
! We use alpha=1.0 to improve precision. 
! With alpha=1.0 and rho=-alpha*theta, the rho variable
! is then just minus theta. 
!
! -minimum vertical density gradient rho_grad_min is necessary to 
! avoid errors with truncation in the division by drho/dz when compute kappa.
! rho_grad_min corresponds roughly to the precision of the computation.
! Physically, with
!      
! N^2 = -(g/rho0)(drho/dz)
!      
! then rho_grad_min sets a minimum N^2 resolved. 
! This corresponds to a frequency f=N/2pi.  The typical 
! period of inertial oscillations in the deep ocean is 6hrs
! (Pickard and Emery, page 55-56).  In the upper ocean, it is
! 10-30 minutes, in pycnocline it is smaller still.  
! So to cover the majority of the ocean's stratification,
! we will want to set rho_grad_min to something smaller than 9e-6. 
!
! To bin the effective diffusivity, it is also useful to have a max
! vertical density gradient.        
!
!-versions: 
!
!  mom4p0 method assumed rigid lid, or zero surface height undulations.
!         Fit the following equation to the model data
!         \partial_{t}(rho_sort) = (F_{k}-F_{k-1})/dzw
! 
!  mom4p1 method fits the following equation to model data 
!         \partial_{t}(dzt_sort*theta_sort) = F_{k} - F_{k-1}
!         Fits this equation assuming two-time level tendency 
!
!  revision: 05/2005
!  revision: 07/2007
!  Stephen.Griffies
!
! </DESCRIPTION>
!
subroutine diagnose_kappa_sort(Time, Thickness, Theta)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: Theta

  real, dimension(:,:,:), allocatable ::  global_tmask  ! for global mask
  real, dimension(:,:,:), allocatable ::  global_dzt    ! for global dzt
  real, dimension(:,:,:), allocatable ::  global_tracer ! for global tracer 
  real, dimension(:,:)  , allocatable ::  global_dat    ! for global area 
  real, dimension(:)    , allocatable ::  rho_sort      ! density (kg/m^3) of sorted parcels
  real, dimension(:)    , allocatable ::  vol_sort      ! volume (m^3) of sorted parcels

  logical, save :: first_enter_routine=.true.
  integer       :: tau, taup1

  integer :: i, j, k, m, n
  integer :: nsort, numzero
 
  real :: kappa_prev, tmp 
  real :: difference

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not.module_is_initialized) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_tracer_diag (diagnose_kappa_sort): module needs initialization ')
  endif

! taum1=tau since tendency=TWO_LEVEL
  tau   = Time%tau
  taup1 = Time%taup1

  allocate (global_tmask(ni,nj,nk)) ; global_tmask=0.0
  call mpp_global_field(Dom%domain2d, Grd%tmask, global_tmask)
  allocate (global_dat(ni,nj))      ; global_dat=0.0
  call mpp_global_field(Dom%domain2d, Grd%dat, global_dat)

  if (first_enter_routine) then
      first_enter_routine = .false.

      nsortpts    = Grd%wet_t_points
      thick_column= 0.0
      mass_column = 0.0

! for diagnosing kappa_sort
      allocate (nlayer(nk))        ; nlayer      = 0 
      allocate (wetarea(nk))       ; wetarea     = 0.0 
      allocate (normalize(nk))     ; normalize   = 0.0 
      allocate (deriv_lay(nk))     ; deriv_lay   = 0.0
      allocate (flux_lay(nk))      ; flux_lay    = 0.0 
      allocate (kappa_lay(nk))     ; kappa_lay   = 0.0
      allocate (rho_lay(nk,2))     ; rho_lay     = 0.0
      allocate (dzt_rho_lay(nk,2)) ; dzt_rho_lay = 0.0
      allocate (dzt_lay(nk,2))     ; dzt_lay     = 0.0
      allocate (dzw_lay(nk,2))     ; dzw_lay     = 0.0
      allocate (mass_lay(nk,2))    ; mass_lay    = 0.0
      allocate (volume_lay(nk,2))  ; volume_lay  = 0.0 

      do k=1,nk      
         nlayer(k) = 0
         wetarea(k)= 0.0
         do j=1,nj
            do i=1,ni
               if(global_tmask(i,j,k) == 1.0) then 
                   nlayer(k)  = nlayer(k)  + 1
                   wetarea(k) = wetarea(k) + global_dat(i,j)
               endif
            enddo
         enddo
         if(nlayer(k) > 0) then 
             normalize(k) = rho0r/wetarea(k)
         else
             normalize(k) = 0.0
         endif 
      enddo

  else 

      do k=1,nk
        rho_lay(k,1)     = rho_lay(k,2)
        dzt_rho_lay(k,1) = dzt_rho_lay(k,2)
        dzt_lay(k,1)     = dzt_lay(k,2)
        dzw_lay(k,1)     = dzw_lay(k,2)
        mass_lay(k,1)    = mass_lay(k,2)
        volume_lay(k,1)  = volume_lay(k,2)
      enddo

  endif

  allocate (global_dzt(ni,nj,nk))    ; global_dzt=0.0
  allocate (global_tracer(ni,nj,nk)) ; global_tracer=0.0
  allocate (vol_sort(nsortpts))      ; vol_sort=0.0
  allocate (rho_sort(nsortpts))      ; rho_sort=0.0

! taup1 time level
! remove rho0 from dzt at a later point
  call mpp_global_field(Dom%domain2d, Theta%field(:,:,:,taup1), global_tracer)
  call mpp_global_field(Dom%domain2d, Thickness%dst(:,:,:), global_dzt)

  nsort=0 
  do k=1,nk      
     do j=1,nj
        do i=1,ni
           if(global_tmask(i,j,k) == 1.0) then 
               nsort=nsort+1
               vol_sort(nsort) = rho0*global_dzt(i,j,k)*global_dat(i,j)
               rho_sort(nsort) = -global_tracer(i,j,k)             ! assume alpha=1.0
           endif
        enddo
     enddo
  enddo
  deallocate(global_tracer)
  deallocate(global_dzt)
  deallocate(global_dat)

  if(debug_diagnose_mixingB) then
      write(stdoutunit,'(/a)') '==========Before sort============' 
      do nsort=1,nsortpts 
         write(stdoutunit,'(a,i7,a,e17.11,a,i7,a,e17.11)') &
              ' rho_sort(',nsort,') = ',rho_sort(nsort), &
              ' vol_sort(',nsort,') = ',vol_sort(nsort)*rho0r 
      enddo
  endif

! reorder rho_sort and carry vol_sort along with it. 
! after the sorting, we will have 
! rho_sort(1)        = smallest rho_sort value, and vol_sort(1)        is its volume
! rho_sort(nsortpts) = largest  rho_sort value, and vol_sort(nsortpts) is its volume
  call sort_shell_array(rho_sort,vol_sort)

  if(debug_diagnose_mixingC) then
      write(stdoutunit,'(/a)') '==========After sort============' 
      do nsort=1,nsortpts 
         write(stdoutunit,'(a,i7,a,e17.11,a,i7,a,e17.11)') &
              ' rho_sort(',nsort,') = ',rho_sort(nsort), &
              ' vol_sort(',nsort,') = ',vol_sort(nsort)*rho0r 
      enddo
  endif

! layer properties, where layers defined by fixed number of points per layer 
  nsort=0
  do k=1,nk
     do n=1,nlayer(k)
        nsort=nsort+1
        dzt_lay(k,2)     = dzt_lay(k,2)     + vol_sort(nsort)
        dzt_rho_lay(k,2) = dzt_rho_lay(k,2) + vol_sort(nsort)*rho_sort(nsort)
     enddo
     dzt_lay(k,2)     = dzt_lay(k,2)*normalize(k)
     dzt_rho_lay(k,2) = dzt_rho_lay(k,2)*normalize(k)
  enddo
  
  deallocate(rho_sort)
  deallocate(vol_sort)

  do k=1,nk
     volume_lay(k,2)  = dzt_lay(k,2)*wetarea(k)
     mass_lay(k,2)    = dzt_rho_lay(k,2)*wetarea(k)
     if(volume_lay(k,2) > 0) then  
         rho_lay(k,2) = mass_lay(k,2)/volume_lay(k,2)
     endif
  enddo

  do k=1,nk-1
     dzw_lay(k,2) = 0.5*(dzt_lay(k,2)+dzt_lay(k+1,2))
  enddo
  dzw_lay(nk,2)   = 0.5*dzt_lay(nk,2)  

  do k=1,nk
     thick_column(2) = thick_column(2) + dzt_lay(k,2)
     mass_column(2)  = mass_column(2)  + mass_lay(k,2)
  enddo


! compute derivatives across averaged layers.
! derivatives are defined at bottom of tracer cell,
! with deriv_lay(k=1) at bottom of the top-most cell,
! and deriv_lay(k=nk) at bottom of the bottom-most cell (and so set to zero). 
  deriv_lay(:) = 0.0
  do k=1,nk-1
     if(dzw_lay(k,1) /= 0.0) then 
         deriv_lay(k) = -(rho_lay(k+1,1)-rho_lay(k,1))/dzw_lay(k,1)
     endif
  enddo

! Compute diapycnal flux of sorted density.
! Flux is defined on the bottom face of the t-cell.  
! Start at column top (k=1), specifying no flux condition.
! Bottom (k=nk) should be diagnosed to have zero flux for case
! where there is no buoyancy forcing.  Finding a zero flux 
! at bottom is a useful check on algorithm integrity.
  flux_lay(:) = 0.0
  k=1
  flux_lay(k) = (dzt_rho_lay(k,2)-dzt_rho_lay(k,1))*dtimer
  do k=2,nk
     flux_lay(k) = flux_lay(k-1) + (dzt_rho_lay(k,2)-dzt_rho_lay(k,1))*dtimer
  enddo


! compute effective diffusivity
  kappa_lay(:) = 0.0
  numzero      = 0  
  do k=1,nk
     if( abs(deriv_lay(k)) >= rho_grad_min  .and.  abs(deriv_lay(k)) <= rho_grad_max ) then
         kappa_lay(k) = -flux_lay(k)/deriv_lay(k)
     else 
         numzero = numzero + 1
     endif
  enddo

  if(smooth_kappa_sort>0) then 
      do m=1,smooth_kappa_sort
         kappa_prev = 0.25*kappa_lay(1)
         do k=2,nk-2
            tmp          =  kappa_lay(k)
            kappa_lay(k) =  kappa_prev + 0.5*kappa_lay(k) + 0.25*kappa_lay(k+1)
            kappa_prev   =  0.25*tmp
         enddo
      enddo
  endif


  if(debug_diagnose_mixingA) then

      write(stdoutunit,'(/a/)')'==========Debug diagnostics for kappa_sort calculation========='

      write(stdoutunit,'(/a)')' Number of parcels within each density layer'
      do k=1,nk
         write(stdoutunit,'(a,i3,a,i7)') ' nlayer(',  k,') = ',nlayer(k)
      enddo

      write(stdoutunit,'(/a)')' Thickness(m)  of density layers'
      do k=1,nk
         difference = dzt_lay(k,2)-dzt_lay(k,1)
         write(stdoutunit,'(a,i3,a,e17.11,a,i3,a,e17.11,a,i3,a,e17.11)') &
              ' dzt_lay(',k,',1)= ',dzt_lay(k,1),                      &
              ' dzt_lay(',k,',2)= ',dzt_lay(k,2),                      &
              ' diff(',k,')= ',difference 
      enddo
      write(stdoutunit,'(/a)') ' Total thickness(m) of sorted column'
      write(stdoutunit,'(a,e17.11,a,e17.11,a,e17.11)')     &
           ' thick_column(tau)= '  ,thick_column(tau)  , &
           ' thick_column(taup1)= ',thick_column(taup1), &
           ' diff= ', thick_column(taup1)-thick_column(tau)


      write(stdoutunit,'(/a)')' Distance(m) between density cell centers '
      do k=1,nk
         difference = dzw_lay(k,2)-dzw_lay(k,1)
         write(stdoutunit,'(a,i3,a,e17.11,a,i3,a,e17.11,a,i3,a,e17.11)') &
              ' dzw_lay(',k,',1)= ',dzw_lay(k,1),                      &
              ' dzw_lay(',k,',2)= ',dzw_lay(k,2),                      &
              ' diff(',k,')= ',difference 
      enddo

      write(stdoutunit,'(/a)') ' Mass (kg) of layers' 
      do k=1,nk 
         difference = mass_lay(k,2)-mass_lay(k,1) 
         write(stdoutunit,'(a,i3,a,e17.11,a,i3,a,e17.11,a,i3,a,e17.11)') &
              ' mass_lay(',k,',1)= ',mass_lay(k,1),                    &
              ' mass_lay(',k,',2)= ',mass_lay(k,2),                    &
              ' diff(',k,')= ',difference 
      enddo
      write(stdoutunit,'(/a)') ' Total mass(kg) of sorted column'
      write(stdoutunit,'(a,e17.11,a,e17.11,a,e17.11)')   &
           ' mass_column(tau)= '  ,mass_column(tau)  , &
           ' mass_column(taup1)= ',mass_column(taup1), &
           ' diff= ', mass_column(taup1)-mass_column(tau)

      write(stdoutunit,'(/a)') ' Volume (m^3) of layers' 
      do k=1,nk 
         difference = volume_lay(k,2)-volume_lay(k,1) 
         write(stdoutunit,'(a,i3,a,e17.11,a,i3,a,e17.11,a,i3,a,e17.11)') &
              ' volume_lay(',k,',1)= ',volume_lay(k,1),                &
              ' volume_lay(',k,',2)= ',volume_lay(k,2),                &
              ' diff(',k,')= ',difference 
      enddo

      write(stdoutunit,'(/a)') ' Density (kg/m^3) of layers' 
      do k=1,nk 
         difference = rho_lay(k,2)-rho_lay(k,1) 
         write(stdoutunit,'(a,i3,a,e17.11,a,i3,a,e17.11,a,i3,a,e17.11)') &
              ' rho_lay(',k,',1)= ',rho_lay(k,1),                      &
              ' rho_lay(',k,',2)= ',rho_lay(k,2),                      &
              ' diff(',k,')= ',difference 
      enddo

      write(stdoutunit,'(/a)') ' Thickness weighted density dzt*rho (kg/m^2)' 
      do k=1,nk
         difference = dzt_rho_lay(k,2)-dzt_rho_lay(k,1) 
         write(stdoutunit,'(a,i3,a,e17.11,a,i3,a,e17.11,a,i3,a,e17.11)') &
              ' dzt_rho_lay(',k,',1)= ',dzt_rho_lay(k,1),              &
              ' dzt_rho_lay(',k,',2)= ',dzt_rho_lay(k,2),              &
              ' diff(',k,')= ',difference
      enddo

      write(stdoutunit,'(/a)')' Terms contributing to effective diffusivity'
      do k=1,nk
         write(stdoutunit,'(a,i3,a,e17.11,a,i3,a,e17.11)') &
              '  flux_lay(',k,')  = ',flux_lay(k),       &
              '  deriv_lay(',k,') = ',deriv_lay(k)
      enddo

      write (stdoutunit,'(/a)') ' Effective diapycnal diffusivity (m^2/sec)'
      do k=1,nk
         write(stdoutunit,'(a,i3,a,e17.11)')' kappa_sort(',k,') = ',kappa_lay(k) 
      enddo
      write(stdoutunit,'(a)') '============================================================='

  endif


  if (id_kappa_sort > 0) then 
      used = send_data (id_kappa_sort, kappa_lay(:), Time%model_time)
  endif
  if (id_temp_sort > 0) then 
      used = send_data (id_temp_sort, -rho_lay(:,1)/alpha, Time%model_time)
  endif


end subroutine diagnose_kappa_sort
! </SUBROUTINE>  NAME="diagnose_kappa_sort"

  
!#######################################################################
! <SUBROUTINE NAME="diagnose_kappa_simple">
!
! <DESCRIPTION>
! Routine to diagnose the amount of mixing between classes of a 
! particular tracer.  Temperature is used as default.
!
! Compute horizontal average of temp to define a stable profile.
! Evolution of this profile defines an effective diffusity. 
!
! This diffusivity is different than the one diagnosed
! from the adiabatic sorting approach.  The sorting approach is 
! more relevant.  The two approaches agree when there 
! is zero baroclinicity, and the present simple scheme is 
! useful ONLY to help debug the sorting routine. 
!
! </DESCRIPTION>
!
subroutine diagnose_kappa_simple(Time, Theta)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_prog_tracer_type), intent(in) :: Theta
  integer                                  :: tau, taum1, taup1

  real :: kappa_simple(nk)  ! diagnosed dianeutral diffusivity (m^2/sec)
  real :: deriv_simple(nk)  ! vertical density gradient centered at top of sorted tracer cell
  real :: flux_simple(nk)   ! dianuetral tracer flux centered at top of sorted tracer cell
  real :: rho_simple(nk,3)  ! density 

  real, dimension(:,:), allocatable     ::  global_dat     ! for global area 
  real, dimension(:,:,:), allocatable   ::  global_tmask   ! for global mask
  real, dimension(:,:,:,:), allocatable ::  global_tracer  ! for global tracer 
  
  integer :: i, j, k
  real    :: cellarea_r 

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag (diagnose_kappa_simple): module needs initialization ')
  endif 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

! NOTE: minimum vertical density gradient rho_grad_min is necessary to 
! avoid errors with truncation in the division by drho/dz when compute kappa.
! rho_grad_min corresponds roughly to the precision of the computation.
! Physically, with
!      
! N^2 = -(g/rho0)(drho/dz)
!      
! then rho_grad_min sets a minimum N^2 resolved. 
! This corresponds to a frequency f=N/2pi.  The typical 
! period of inertial oscillations in the deep ocean is 6hrs
! (Pickard and Emery, page 56).  In the upper ocean, it is
! 10-30 minutes.  So to cover the majority of the ocean's stratification,
! we will want to set rho_grad_min to something smaller than 9e-6. 
!
! To bin the effective diffusivity, it is also useful to have a max
! vertical density gradient.        

  allocate (global_dat(ni,nj)) ; global_dat=0.0
  call mpp_global_field(Dom%domain2d, Grd%dat, global_dat)
  allocate (global_tmask(ni,nj,nk)) ; global_tmask=0.0
  call mpp_global_field(Dom%domain2d, Grd%tmask, global_tmask)

! horizontal tracer area 
  if(debug_diagnose_mixingA) then
      write(stdoutunit,'(a)') ' ' 
      write(stdoutunit,'(a,e17.11)')' in simple, total cross area of domain = ',Grd%tcellsurf
      write(stdoutunit,'(a,e17.11)')' in simple, total volume of domain     = ',Grd%tcellsurf*Grd%zw(nk)
  endif
  cellarea_r = 1.0/(epsln + Grd%tcellsurf)

! horizontally average to get vertical density profile 
  allocate (global_tracer(ni,nj,nk,3)) ; global_tracer=0.0
  call mpp_global_field(Dom%domain2d, Theta%field(:,:,:,:), global_tracer)
  rho_simple(:,:) = 0.0 
  do k=1,nk      
     do j=1,nj
        do i=1,ni
           if(global_tmask(i,j,k) == 1.0) then 
             rho_simple(k,:) = -alpha*global_tracer(i,j,k,:)*global_dat(i,j) + rho_simple(k,:) 
           endif
        enddo
     enddo
     rho_simple(k,:) = rho_simple(k,:)*cellarea_r
  enddo
  deallocate(global_tracer)
  deallocate(global_tmask)
  deallocate(global_dat)

! compute derivatives across averaged layers.
! derivatives are defined at bottom of tracer cell,
! with deriv_simple(k=1) at bottom of top-most cell,
! and deriv_simple(k=nk)=0 at bottom of bottom-most cell. 
  deriv_simple(:) = 0.0
  do k=1,nk-1
     deriv_simple(k) = -(rho_simple(k+1,taup1)-rho_simple(k,taup1))*Grd%dzwr(k-1)
  enddo
  if(debug_diagnose_mixingA) then
    write(stdoutunit,'(/a)') 'density at three time levels' 
    do k=1,nk
      write(stdoutunit,'(a,i3,a,e17.11,a,i3,a,e17.11,a,i3,a,e17.11)') &
       ' rho_simple(taum1,',k,')= ',rho_simple(k,taum1)+rho0, &
       ' rho_simple(tau(',k,')= ',rho_simple(k,tau)+rho0,&
       ' rho_simple(taup1(',k,')= ',rho_simple(k,taup1)+rho0
    enddo 
    write(stdoutunit,'(/a)') 'vertical derivative of density' 
    do k=1,nk
      write(stdoutunit,'(a,i3,a,e17.11)')' deriv_simple(',k,') = ',deriv_simple(k) 
    enddo 
  endif 

! compute effective diffusivity by diagnosing the  
! diapycnal flux. The flux is defined on the bottom
! face of the t-cell. start integration from the 
! ocean top and work down.
  flux_simple(:) = 0.0
  k=1
  flux_simple(k) = (rho_simple(k,taup1)-rho_simple(k,taum1))*Grd%dzt(k)*dtimer
  do k=2,nk
     flux_simple(k) = flux_simple(k-1) + (rho_simple(k,taup1)-rho_simple(k,taum1))*Grd%dzt(k)*dtimer
  enddo

  if(debug_diagnose_mixingA) then
    write(stdoutunit,'(/a)')'vertical flux across a layer' 
    do k=1,nk
      write(stdoutunit,'(a,i3,a,e17.11)')' flux_simple(',k,') = ',flux_simple(k) 
    enddo 
  endif 

  do k=1,nk
     if(abs(deriv_simple(k)) >= rho_grad_min .and. abs(deriv_simple(k)) <= rho_grad_max) then
         kappa_simple(k) = -flux_simple(k)/deriv_simple(k)
     else 
         kappa_simple(k) = 0.0
     endif
  enddo

  if (id_kappa_simple > 0) then 
     used = send_data (id_kappa_simple, kappa_simple(:), Time%model_time)
  endif 

  if(debug_diagnose_mixingA) then
    write (stdoutunit,'(/a)') ' Effective diffusivity (m^2/sec)'
    do k=1,nk
      write(stdoutunit,'(a,i3,a,e17.11)')' kappa_simple(',k,') = ',kappa_simple(k) 
    enddo 
  endif 

end subroutine diagnose_kappa_simple
! </SUBROUTINE>  NAME="diagnose_kappa_simple"


!#######################################################################
! <SUBROUTINE NAME="diagnose_depth_of_potrho">
! <DESCRIPTION>
!
! Diagnose depth (m) of a potential density surface surface relative to 
! the ocean surface at z=eta (not relative to z=0).   
!
! Method uses linear interpolation to find the depth of a 
! potential rho surface.
!
! Scheme currently does not forward (backwards) interpolate if 
! rho surface lies within lowest (uppermost) grid cell.
!
! Diagnostic only makes sense when rho is monotonically
! increasing as go deeper in water column. 
!
! Author: Harper.Simmons
!         Zhi.Liang
! </DESCRIPTION>
!
subroutine diagnose_depth_of_potrho(Time, Dens, Thickness)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_density_type),   intent(in) :: Dens
  type(ocean_thickness_type), intent(in) :: Thickness

  real      :: w1,w2
  integer   :: potrho_nk
  integer   :: i, j, k, n
  real, dimension(isd:ied,jsd:jed,size(Dens%potrho_ref(:))) :: depth_of_potrho 

  potrho_nk = size(Dens%potrho_ref(:))

  if (.not.module_is_initialized) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_tracer_diag (diagnose_depth_of_potrho): module needs initialization ')
  endif

  depth_of_potrho(:,:,:)=-10.0
  do n=1,potrho_nk
     do j=jsc,jec
        do i=isc,iec
kloop:     do k=nk-1,1,-1
              if(    Dens%potrho(i,j,k) < Dens%potrho_ref(n)) then
                  if(Dens%potrho_ref(n) < Dens%potrho(i,j,k+1)) then
                      if(Grd%tmask(i,j,k+1) > 0) then 
                          W1= Dens%potrho_ref(n)   - Dens%potrho(i,j,k)
                          W2= Dens%potrho(i,j,k+1) - Dens%potrho_ref(n)
                          depth_of_potrho(i,j,n) = ( Thickness%depth_zt(i,j,k+1)*W1  &
                                                    +Thickness%depth_zt(i,j,k)  *W2) &
                                                    /(W1 + W2 + epsln)
                          exit kloop 
                      endif
                  endif
              endif
           enddo kloop
        enddo
     enddo
  enddo
  used = send_data (id_depth_of_potrho, depth_of_potrho(:,:,:), &
         Time%model_time, &
         is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=potrho_nk)

end subroutine diagnose_depth_of_potrho
! </SUBROUTINE>  NAME="diagnose_depth_of_potrho"



!#######################################################################
! <SUBROUTINE NAME="diagnose_depth_of_theta">
! <DESCRIPTION>
!
! Diagnose depth (m) of a potential temperature surface relative to 
! the ocean surface at z=eta (not relative to z=0).   
!
! Method uses linear interpolation to find the depth of a 
! potential temp surface.
!
! Scheme currently does not forward (backwards) interpolate if 
! theta surface lies within lowest (uppermost) grid cell.
!
! Diagnostic only makes sense when theta is monotonically
! decreasing as go deeper in water column. 
!
! Based on "diagnose_depth_of_potrho" by Harper.Simmons
!
! Author: Stephen.Griffies
!         Zhi.Liang 
! </DESCRIPTION>
!
subroutine diagnose_depth_of_theta(Time, Dens, Thickness, T_prog)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_density_type),     intent(in) :: Dens
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)

  real      :: w1,w2
  integer   :: theta_nk
  integer   :: i, j, k, n, tau
  real, dimension(isd:ied,jsd:jed,size(Dens%theta_ref(:))) :: depth_of_theta 

  theta_nk = size(Dens%theta_ref(:))
  tau      = Time%tau
 
  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag (diagnose_depth_of_theta): module needs initialization ')
  endif 

  depth_of_theta(:,:,:)=-10.0
  do n=1,theta_nk
     do j=jsc,jec
        do i=isc,iec
kloop:     do k=nk-1,1,-1
              if(    Dens%theta_ref(n) < T_prog(index_temp)%field(i,j,k,tau)  ) then
                  if(T_prog(index_temp)%field(i,j,k+1,tau) < Dens%theta_ref(n)) then
                      if(Grd%tmask(i,j,k+1) > 0) then
                          W1= T_prog(index_temp)%field(i,j,k,tau) - Dens%theta_ref(n) 
                          W2= Dens%theta_ref(n) - T_prog(index_temp)%field(i,j,k+1,tau) 
                          depth_of_theta(i,j,n) = ( Thickness%depth_zt(i,j,k+1)*W1   &
                                                   +Thickness%depth_zt(i,j,k)  *W2 ) &
                                                  /(W1 + W2 + epsln)
                          exit kloop
                      endif
                  endif
              endif
           enddo kloop
        enddo
     enddo
  enddo
  used = send_data (id_depth_of_theta, depth_of_theta(:,:,:), &
         Time%model_time, &
         is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=theta_nk)

end subroutine diagnose_depth_of_theta
! </SUBROUTINE>  NAME="diagnose_depth_of_theta"


!#######################################################################
! <SUBROUTINE NAME="diagnose_tracer_on_rho">
! <DESCRIPTION>
! Diagnose tracer concentration on potential density surface. 
! Method based on diagnose_depth_of_potrho diagnostic.
!
! Author: Stephen.Griffies
!
! Updated Oct 2009 to be more vectorized 
!
! </DESCRIPTION>
!
subroutine diagnose_tracer_on_rho(Time, Dens, Tracer, ntracer)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_density_type),     intent(in) :: Dens
  type(ocean_prog_tracer_type), intent(in) :: Tracer 
  integer,                      intent(in) :: ntracer                       

  real      :: W1, W2
  integer   :: potrho_nk
  integer   :: i, j, k, k_rho, tau
  real, dimension(isd:ied,jsd:jed,size(Dens%potrho_ref(:))) :: tracer_on_rho 

  potrho_nk = size(Dens%potrho_ref(:))
  tau       = Time%tau

  if (.not. module_is_initialized) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_tracer_diag (diagnose_tracer_on_rho): module needs initialization')
  endif

  tracer_on_rho(:,:,:)=0.0

  ! for (i,j) points with potrho_ref < potrho(k=1),   keep tracer_on_rho=0
  ! for (i,j) points with potrho_ref > potrho(k=kmt), keep tracer_on_rho=0
  ! these assumptions mean there is no need to specially handle the endpoints,
  ! since the initial value for tracer_on_rho is 0.

  ! interpolate tracer field onto rho-surface 
  do k_rho=1,potrho_nk
     do k=1,nk-1
        do j=jsc,jec
           do i=isc,iec
              if(     Dens%potrho_ref(k_rho) >  Dens%potrho(i,j,k)  ) then
                  if( Dens%potrho_ref(k_rho) <= Dens%potrho(i,j,k+1)) then 
                      W1= Dens%potrho_ref(k_rho)- Dens%potrho(i,j,k)
                      W2= Dens%potrho(i,j,k+1)  - Dens%potrho_ref(k_rho)
                      tracer_on_rho(i,j,k_rho) = ( Tracer%field(i,j,k+1,tau)*W1  &
                                                  +Tracer%field(i,j,k,tau)  *W2) &
                                                  /(W1 + W2 + epsln)
                  endif
              endif
           enddo
        enddo
     enddo
  enddo

  ! ensure masking is applied to interpolated field 
  do k_rho=1,potrho_nk
     do j=jsc,jec
        do i=isc,iec
           tracer_on_rho(i,j,k_rho) = tracer_on_rho(i,j,k_rho)*Grd%tmask(i,j,1)
        enddo
     enddo
  enddo

  used = send_data (id_tracer_on_rho(ntracer), tracer_on_rho(:,:,:), &
                    Time%model_time,                                 & 
                    is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=potrho_nk)


end subroutine diagnose_tracer_on_rho
! </SUBROUTINE>  NAME="diagnose_tracer_on_rho"


!#######################################################################
! <SUBROUTINE NAME="diagnose_tracer_zrho_on_rho">
! <DESCRIPTION>
! Diagnose tracer concentration * dz/drho on potential density surface.
! This product, when integrated over dx*dy*drho, will yield the same 
! total tracer (to within roundoff) as the usual tracer concentration 
! integrated over dx*dy*dz.  
! 
! compute abs(dz/drho)==dz/drho in order to have tracer_zrho_on_rho
! with same sign as tracer. 
! 
! Method based on diagnose_tracer_on_rho diagnostic.
!
! Author: Stephen.Griffies
! Updated Oct 2009 to be more vectorized 
!
! </DESCRIPTION>
!
subroutine diagnose_tracer_zrho_on_rho(Time, Dens, Thickness, Tracer, ntracer)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_density_type),     intent(in) :: Dens
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: Tracer 
  integer,                      intent(in) :: ntracer                       

  real      :: W1,W2
  real      :: zrho
  integer   :: potrho_nk
  integer   :: i, j, k, n, tau
  real, dimension(isd:ied,jsd:jed,size(Dens%potrho_ref(:))) :: tracer_zrho_on_rho 

  potrho_nk = size(Dens%potrho_ref(:))
  tau       = Time%tau

  if (.not. module_is_initialized) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_tracer_diag (diagnose_tracer_zrho_on_rho): module needs initialization')
  endif

  tracer_zrho_on_rho(:,:,:)=0.0

  ! for (i,j) points with potrho_ref < potrho(k=1),   keep tracer_zrho_on_rho=0
  ! for (i,j) points with potrho_ref > potrho(k=kmt), keep tracer_zrho_on_rho=0
  ! these assumptions mean there is no need to specially handle the endpoints,
  ! since the initial value for tracer_on_rho is 0.

  ! interpolate field onto rho-surface 
  do n=1,potrho_nk
     do k=nk-1,1,-1
        do j=jsc,jec
           do i=isc,iec
              if(    Dens%potrho(i,j,k) < Dens%potrho_ref(n)) then
                  if(Dens%potrho_ref(n) < Dens%potrho(i,j,k+1)) then
                      if(Grd%tmask(i,j,k+1) > 0) then 

                          zrho = Thickness%dzwt(i,j,k)/(-Dens%potrho(i,j,k)+Dens%potrho(i,j,k+1)+epsln)
                          W1= Dens%potrho_ref(n)   - Dens%potrho(i,j,k)
                          W2= Dens%potrho(i,j,k+1) - Dens%potrho_ref(n)
                          tracer_zrho_on_rho(i,j,n) = zrho* &
                               ( Tracer%field(i,j,k+1,tau)*W1  &
                                +Tracer%field(i,j,k,tau)  *W2) &
                               /(W1 + W2 + epsln)              
                      endif
                  endif
              endif
           enddo
        enddo
     enddo
  enddo
  used = send_data (id_tracer_zrho_on_rho(ntracer), tracer_zrho_on_rho(:,:,:), &
                    Time%model_time,                                 & 
                    is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=potrho_nk)


end subroutine diagnose_tracer_zrho_on_rho
! </SUBROUTINE>  NAME="diagnose_tracer_zrho_on_rho"


!#######################################################################
! <SUBROUTINE NAME="calc_potrho_mixed_layer">
!
! <DESCRIPTION>
! Calculate the mixed layer depth and potential density at mixed layer base  
! according to depth at which buoyancy is greater than buoyancy_crit
! relative to the surface. Compute the buoyancy using potential 
! density, rather than the insitu density, since we aim for this 
! diagnostic to be comparable to diagnostics from isopcynal models. 
!
! Note that the mixed layer depth is taken with respect to the ocean surface,
! and so the mixed layer depth is always positive. That is, the mld is here
! defined as a thickness of water.
!           
! </DESCRIPTION>
!
subroutine calc_potrho_mixed_layer (Thickness, Dens, potrho_mix_depth, potrho_mix_base)

  type(ocean_thickness_type), intent(in)            :: Thickness
  type(ocean_density_type),   intent(in)            :: Dens
  real, dimension(isd:,jsd:), intent(out), optional :: potrho_mix_depth
  real, dimension(isd:,jsd:), intent(out), optional :: potrho_mix_base

  integer :: i, j, k
  real    :: depth, crit
  real, dimension(isd:ied,jsd:jed) :: potrho_mix_base_use

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag (calc_potrho_mixed_layer): module needs initialization ')
  endif 

  crit = buoyancy_crit*rho0/grav 
  potrho_mix_base_use = 0.0

  do j=jsc,jec
     do i=isc,iec
        potrho_mix_base_use(i,j) = Dens%potrho(i,j,1) + crit
     enddo
  enddo

  if (present(potrho_mix_base)) then
    potrho_mix_base  = potrho_mix_base_use
  endif

  if (present(potrho_mix_depth)) then
    potrho_mix_depth = 0.0
    do j=jsc,jec
       do i=isc,iec
          depth=max(0.0,Thickness%dzwt(i,j,0))
          do k=2,Grd%kmt(i,j)
             if( Dens%potrho(i,j,k) < potrho_mix_base_use(i,j)) then 
                 depth=depth+Thickness%dzwt(i,j,k)
             else
                 potrho_mix_depth(i,j) = depth + Thickness%dzwt(i,j,k) &
                                     *(potrho_mix_base_use(i,j)-Dens%potrho(i,j,k-1))&
                                     /(Dens%potrho(i,j,k)  -Dens%potrho(i,j,k-1)) 
                 exit  
             endif
          enddo
       enddo
    enddo
  endif

end subroutine calc_potrho_mixed_layer
! </SUBROUTINE> NAME="calc_potrho_mixed_layer"



!#######################################################################
! <SUBROUTINE NAME="potrho_mixed_layer">
!
! <DESCRIPTION>
! Determine mixed layer depth and potential density at mixed layer base  
! according to depth at which buoyancy is greater than buoyancy_crit
! relative to the surface.
! Call calc_potrho_mixed_layer to calculate the quantities.
!           
! </DESCRIPTION>
!
subroutine potrho_mixed_layer (Time, Thickness, Dens)

  type(ocean_time_type),      intent(in) :: Time 
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_density_type),   intent(in) :: Dens

  type(time_type) :: next_time 
  real            :: potrho_mix_depth(isd:ied,jsd:jed)
  real            :: potrho_mix_base(isd:ied,jsd:jed) 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag (potrho_mixed_layer): module needs initialization ')
  endif 

  next_time = increment_time(Time%model_time, int(dtts), 0)
  if (need_data(id_potrho_mix_depth,next_time) .and. id_potrho_mix_base > 0) then
    call calc_potrho_mixed_layer (Thickness, Dens,        &
         potrho_mix_depth = potrho_mix_depth, potrho_mix_base = potrho_mix_base)
  elseif (need_data(id_potrho_mix_depth,next_time)) then
    call calc_potrho_mixed_layer (Thickness, Dens,        &
         potrho_mix_depth = potrho_mix_depth)
  elseif (id_potrho_mix_base > 0) then
    call calc_potrho_mixed_layer (Thickness, Dens,        &
         potrho_mix_base = potrho_mix_base)
  endif

  call diagnose_2d(Time, Grd, id_potrho_mix_depth, potrho_mix_depth(:,:))
  call diagnose_2d(Time, Grd, id_potrho_mix_base, potrho_mix_base(:,:))

end subroutine potrho_mixed_layer
! </SUBROUTINE> NAME="potrho_mixed_layer"


!#######################################################################
! <SUBROUTINE NAME="send_total_mass">
!
! <DESCRIPTION>
! Send total liquid seawater mass to diagnostic manager.
! </DESCRIPTION>
!
subroutine send_total_mass(Time, Thickness)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness

  integer  :: tau
  real     :: tmass

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag (send_total_mass): module needs initialization ')
  endif 

  tau   = Time%tau
  tmass = total_mass(Thickness,tau)
  tmass = tmass
  used  = send_data (id_mass_seawater, tmass, Time%model_time)

end subroutine send_total_mass
! </SUBROUTINE>  NAME="send_total_mass"


!#######################################################################
! <SUBROUTINE NAME="send_total_volume">
!
! <DESCRIPTION>
! Send total liquid seawater mass to diagnostic manager.
! </DESCRIPTION>
!
subroutine send_total_volume(Time, Thickness)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness

  integer  :: tau
  real     :: tvolume

  tau = Time%tau

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag (send_total_volume): module needs initialization ')
  endif 

  tvolume = total_volume(Thickness, tau)
  tvolume = tvolume
  used    = send_data (id_volume_seawater, tvolume, Time%model_time)

end subroutine send_total_volume
! </SUBROUTINE>  NAME="send_total_volume"


!#######################################################################
! <SUBROUTINE NAME="send_total_tracer">
!
! <DESCRIPTION>
! Send total tracer to diagnostic manager.
! </DESCRIPTION>
!
subroutine send_total_tracer(Time, Thickness, Tracer, ntracer)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  integer,                      intent(in) :: ntracer 

  integer  :: tau
  real     :: ttracer

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag (send_total_tracer): module needs initialization ')
  endif 

  tau   = Time%tau

  ttracer = total_tracer(Tracer,Thickness,tau)
  if(Tracer%name=='temp') then 
      ttracer = ttracer*1e-25
  elseif(Tracer%name(1:3)=='age') then 
      ttracer = ttracer
  else 
      ttracer = ttracer*1e-18
  endif
  used = send_data (id_prog_total(ntracer), ttracer, Time%model_time)

end subroutine send_total_tracer
! </SUBROUTINE>  NAME="send_total_tracer"


!#######################################################################
! <SUBROUTINE NAME="send_global_ave_tracer">
!
! <DESCRIPTION>
! Send global averaged  tracer to diagnostic manager.
! </DESCRIPTION>
!
subroutine send_global_ave_tracer(Time, Thickness, Tracer, ntracer)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  integer,                      intent(in) :: ntracer 

  real    :: ave_tracer, totaltracer, totalmass 
  integer :: tau 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag (send_global_ave_tracer): module needs initialization ')
  endif 

  tau=Time%tau 

  totaltracer = total_tracer(Tracer,Thickness,tau)
  totalmass   = total_mass(Thickness,tau)
  ave_tracer  = totaltracer/totalmass/Tracer%conversion

  used = send_data (id_prog_global_ave(ntracer), ave_tracer, Time%model_time)

end subroutine send_global_ave_tracer
! </SUBROUTINE>  NAME="send_global_ave_tracer"


!#######################################################################
! <SUBROUTINE NAME="send_global_ave_pressure">
!
! <DESCRIPTION>
! Send global averaged pressure to diagnostic manager.
! </DESCRIPTION>
!
subroutine send_global_ave_pressure(Time, Thickness, Dens)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_density_type),     intent(in) :: Dens

  real    :: total_press, press_ave, totalmass 
  real, dimension(isd:ied,jsd:jed) :: press_k
  integer :: tau, k

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag (send_global_ave_press): module needs initialization')
  endif 

  tau=Time%tau 

  total_press  = 0.0
  press_k(:,:) = 0.0
  do k=1,nk
    press_k(:,:) =  press_k(:,:) + Grd%tmask(:,:,k)*Grd%dat(:,:) &
                    *Thickness%rho_dzt(:,:,k,tau)*Dens%pressure_at_depth(:,:,k)
  enddo
  if(have_obc) press_k(:,:) = press_k(:,:)*Grd%obc_tmask(:,:)
  total_press = mpp_global_sum(Dom%domain2d, press_k(:,:), global_sum_flag)
  
  totalmass = total_mass(Thickness,tau)
  press_ave = total_press/totalmass

  used = send_data (id_press_ave, press_ave, Time%model_time)

end subroutine send_global_ave_pressure
! </SUBROUTINE>  NAME="send_global_ave_pressure"


!#######################################################################
! <SUBROUTINE NAME="send_surface_ave_tracer">
!
! <DESCRIPTION>
! Send global averaged surface tracer to diagnostic manager.
! Note the presence of a rho_dzt weighting here...
! </DESCRIPTION>
!
subroutine send_surface_ave_tracer(Time, Thickness, Tracer, ntracer)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  integer,                      intent(in) :: ntracer 

  real    :: ave_tracer, totaltracer, totalmass 
  integer :: tau 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag (send_surface_ave_tracer): module needs initialization ')
  endif 

  tau=Time%tau 

  totaltracer = klevel_total_tracer(Tracer,Thickness,tau,1)
  totalmass   = klevel_total_mass(Thickness,tau,1)
  ave_tracer  = totaltracer/totalmass/Tracer%conversion

  used = send_data (id_prog_surface_ave(ntracer), ave_tracer, Time%model_time)

end subroutine send_surface_ave_tracer
! </SUBROUTINE>  NAME="send_surface_ave_tracer"

!#######################################################################
! <SUBROUTINE NAME="send_surface_area_ave_tracer">
!
! <DESCRIPTION>
! Send global area averaged surface tracer to diagnostic manager.
! Note the weigthing is just area, with no thickness nor density. 
! </DESCRIPTION>
!
subroutine send_surface_area_ave_tracer(Time, Tracer, ntracer)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  integer,                      intent(in) :: ntracer 

  real, dimension(isd:ied,jsd:jed) :: tmp
  real    :: ave_tracer
  integer :: i,j,tau 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag (send_surface_area_ave_tracer): module needs initialization')
  endif 

  tau=Time%tau 

  tmp(:,:) = 0.0
  do j=jsc,jec
    do i=isc,iec
       tmp(i,j) = tmp(i,j) + Tracer%field(i,j,1,tau)*area_t(i,j) 
     enddo
  enddo
  ave_tracer = mpp_global_sum(Dom%domain2d, tmp(:,:), global_sum_flag)
  ave_tracer = ave_tracer*cellarea_r

  used = send_data (id_prog_surface_area_ave(ntracer), ave_tracer, Time%model_time)

end subroutine send_surface_area_ave_tracer
! </SUBROUTINE>  NAME="send_surface_area_ave_tracer"


!#######################################################################
! <SUBROUTINE NAME="send_tracer_variance">
!
! <DESCRIPTION>
!
! Compute the global and k-level tracer variance.
!
! </DESCRIPTION>
!
subroutine send_tracer_variance(Time, Tracer, Thickness, n)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  type(ocean_thickness_type),   intent(in) :: Thickness
  integer,                      intent(in) :: n

  integer             :: i,j,k,tau
  real                :: global_variance, tmass, ttracer, ttracer2
  real, dimension(nk) :: k_variance, k_tmass, k_ttracer, k_ttracer2

  tau = Time%tau

  if (id_global_variance(n) > 0 .or. id_k_variance(n) > 0) then
     
     if (use_blobs) then
        wrk1(:,:,:) = 0.0
        wrk2(:,:,:) = 0.0
        wrk3(:,:,:) = 0.0
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 wrk1(i,j,k) = Grd%dat(i,j)*Thickness%rho_dztT(i,j,k,tau)*Grd%tmask(i,j,k) !Global Mass
                 wrk2(i,j,k) = wrk1(i,j,k)*Tracer%fieldT(i,j,k)*Tracer%fieldT(i,j,k)
                 wrk3(i,j,k) = wrk1(i,j,k)*Tracer%fieldT(i,j,k)
               enddo
           enddo
        enddo
        
     else!not use_blobs
        wrk1(:,:,:) = 0.0
        wrk2(:,:,:) = 0.0
        wrk3(:,:,:) = 0.0
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 wrk1(i,j,k) = Grd%dat(i,j)*Thickness%rho_dzt(i,j,k,tau)*Grd%tmask(i,j,k)
                 wrk2(i,j,k) = wrk1(i,j,k)*Tracer%field(i,j,k,tau)*Tracer%field(i,j,k,tau)
                 wrk3(i,j,k) = wrk1(i,j,k)*Tracer%field(i,j,k,tau)
              enddo
           enddo 
        enddo
        
     endif   ! endif for blobs 

     if (id_global_variance(n) > 0) then
        tmass    = mpp_global_sum(Dom%domain2d, wrk1(:,:,:))
        ttracer2 = mpp_global_sum(Dom%domain2d, wrk2(:,:,:))
        ttracer  = mpp_global_sum(Dom%domain2d, wrk3(:,:,:))
        if(tmass > 0.0) then     
           global_variance = ttracer2/tmass - ( ttracer/tmass )**2
        else 
           global_variance = 0.0
        endif 
        used = send_data (id_global_variance(n), global_variance, Time%model_time)
     endif

     if (id_k_variance(n) > 0) then
        do k=1,nk
           k_tmass(k)    = mpp_global_sum(Dom%domain2d, wrk1(:,:,k))
           k_ttracer2(k) = mpp_global_sum(Dom%domain2d, wrk2(:,:,k))
           k_ttracer(k)  = mpp_global_sum(Dom%domain2d, wrk3(:,:,k))
           if(k_tmass(k) > 0.0) then 
              k_variance(k) = k_ttracer2(k)/k_tmass(k) - ( k_ttracer(k)/k_tmass(k) )**2
           else
              k_variance(k) = 0.0
           endif 
        enddo
        used = send_data( id_k_variance(n), k_variance(:), Time%model_time )
     endif

  endif

end subroutine send_tracer_variance
! </SUBROUTINE>  NAME="send_tracer_variance"


!#######################################################################
! <SUBROUTINE NAME="diagnose_eta_tend_3dflux">
!
! <DESCRIPTION>
! Diagnose contribution to global mean sea level evolution arising  
! from a 3d MOM flux computed from a parameterization. 
!
! fluxes are assumed to have the following dimensions:
! flux_x = (dy*dz)*diffusivity*rho*tracer_derivative_x
! flux_y = (dx*dz)*diffusivity*rho*tracer_derivative_y 
! flux_z =         diffusivity*rho*tracer_derivative_z 
!
! Subroutine history: 
! Jan2012 version 1.0: Stephen.Griffies 
! 
! </DESCRIPTION>
!
subroutine diagnose_eta_tend_3dflux (Time, Thickness, Dens,&
           flux_x_temp, flux_y_temp, flux_z_temp,           & 
           flux_x_salt, flux_y_salt, flux_z_salt,           &
           eta_tend, eta_tend_glob)  

  type(ocean_time_type),         intent(in)    :: Time
  type(ocean_thickness_type),    intent(in)    :: Thickness
  type(ocean_density_type),      intent(in)    :: Dens
  real, dimension(isd:,jsd:,:),  intent(in)    :: flux_x_temp
  real, dimension(isd:,jsd:,:),  intent(in)    :: flux_y_temp
  real, dimension(isd:,jsd:,:),  intent(in)    :: flux_z_temp
  real, dimension(isd:,jsd:,:),  intent(in)    :: flux_x_salt
  real, dimension(isd:,jsd:,:),  intent(in)    :: flux_y_salt
  real, dimension(isd:,jsd:,:),  intent(in)    :: flux_z_salt
  real, dimension(isd:,jsd:),    intent(inout) :: eta_tend
  real,                          intent(inout) :: eta_tend_glob 

  integer :: i, j, k, kp1, tau

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_tracer_util_mod (diagnose_eta_tend_3dflux): module needs initialization ')
  endif 

  tau = Time%tau

  wrk1(:,:,:) = 0.0   !  nu_theta    = alpha/rho  = -drhodT/rho**2
  wrk2(:,:,:) = 0.0   !  nu_salinity = -beta/rho  = -drhodS/rho**2
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk1(i,j,k) = -Grd%tmask(i,j,k)*Dens%drhodT(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2) 
           wrk2(i,j,k) = -Grd%tmask(i,j,k)*Dens%drhodS(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2) 
        enddo
     enddo
  enddo

  ! compute partial_x(nu_theta) and partial_x(nu_salinity) at east  face of T-cell 
  ! compute partial_y(nu_theta) and partial_y(nu_salinity) at north face of T-cell 
  wrk1_v(:,:,:,1) = 0.0   !  partial_x(nu_theta)    in 1/(deg C * rho * meter)
  wrk1_v(:,:,:,2) = 0.0   !  partial_x(nu_salinity) in 1/(ppt   * rho * meter)
  wrk2_v(:,:,:,1) = 0.0   !  partial_y(nu_theta)    in 1/(deg C * rho * meter)
  wrk2_v(:,:,:,2) = 0.0   !  partial_y(nu_salinity) in 1/(ppt   * rho * meter)
  do k=1,nk
     wrk3(:,:,1)     = wrk1(:,:,k)
     wrk4(:,:,1)     = wrk2(:,:,k)
     wrk1_v(:,:,k,1) = FDX_T(wrk3(:,:,1))*FMX(Grd%tmask(:,:,k))
     wrk1_v(:,:,k,2) = FDX_T(wrk4(:,:,1))*FMX(Grd%tmask(:,:,k))
     wrk2_v(:,:,k,1) = FDY_T(wrk3(:,:,1))*FMY(Grd%tmask(:,:,k))
     wrk2_v(:,:,k,2) = FDY_T(wrk4(:,:,1))*FMY(Grd%tmask(:,:,k))
  enddo

  ! compute partial_z(nu_theta) and partial_z(nu_salinity) at bottom of T-cell 
  wrk3_v(:,:,:,1) = 0.0   !  partial_z(nu_theta)    in 1/(deg C * rho * meter)
  wrk3_v(:,:,:,2) = 0.0   !  partial_z(nu_salinity) in 1/(ppt   * rho * meter)
  do k=1,nk-1
     kp1=k+1
     do j=jsc,jec
        do i=isc,iec
           wrk3_v(i,j,k,1) = Grd%tmask(i,j,kp1)*(wrk1(i,j,k)-wrk1(i,j,kp1))/Thickness%dzwt(i,j,k)
           wrk3_v(i,j,k,2) = Grd%tmask(i,j,kp1)*(wrk2(i,j,k)-wrk2(i,j,kp1))/Thickness%dzwt(i,j,k)
        enddo
     enddo
  enddo

  ! horizontal flux components have a dzt factor already
  ! included (thickness weighted fluxes).
  ! dytr and dxtr multipliers account for extra length 
  ! factors appearing in the horizontal flux components. 
  ! the vertical flux component is multiplied by dzt
  ! to anticipate the vertical integral, and to bring 
  ! it into consistent dimensions with the horizontal flux
  ! components.  
  wrk3(:,:,:) = 0.0
  do k=1,nk
     kp1=min(nk,k+1)
     do j=jsc,jec
        do i=isc,iec
           wrk3(i,j,k) =   Grd%tmask(i+1,j,k)*Grd%dytr(i,j)                                         &
                         *(flux_x_temp(i,j,k)*wrk1_v(i,j,k,1) + flux_x_salt(i,j,k)*wrk1_v(i,j,k,2)) &
                         + Grd%tmask(i,j+1,k)*Grd%dxtr(i,j)                                         &
                         *(flux_y_temp(i,j,k)*wrk2_v(i,j,k,1) + flux_y_salt(i,j,k)*wrk2_v(i,j,k,2)) &
                         + Grd%tmask(i,j,kp1)*Thickness%dzt(i,j,k)                                  &
                         *(flux_z_temp(i,j,k)*wrk3_v(i,j,k,1) + flux_z_salt(i,j,k)*wrk3_v(i,j,k,2))
        enddo
     enddo
  enddo

  ! sea level tendency via a vertical integral.  
  ! minus sign switches from MOM convention of an "upgradient"
  ! flux, to the desired downgradient form. 
  eta_tend(:,:) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           eta_tend(i,j) = eta_tend(i,j) - wrk3(i,j,k)
        enddo
     enddo
  enddo
  
  ! global integral 
  wrk1_2d(:,:)  = Grd%tmask(:,:,1)*Grd%dat(:,:)*eta_tend(:,:)
  eta_tend_glob = mpp_global_sum(Dom%domain2d, wrk1_2d(:,:), global_sum_flag)*cellarea_r


end subroutine diagnose_eta_tend_3dflux
! </SUBROUTINE> NAME="diagnose_eta_tend_3dflux"


end module ocean_tracer_diag_mod



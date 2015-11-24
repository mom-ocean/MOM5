module ocean_form_drag_mod
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Compute form drag as per Greatbatch and Lamb (1990) and/or 
! Aiki etal (2004).
!</OVERVIEW>
!
!<DESCRIPTION>
! Compute thickness weighted and density weighted tendency
! [velocity*(kg/m^3)*meter/sec] for velocity associated with 
! parameterized form drag.  Code employs the methods described
! in Greatbatch and Lamb (1990) as well as Aiki etal (2004).
!
! This scheme has not been updated for C-grid specific layout.
! However, there is little reason to expect that shifting the 
! viscosity by 1/2 grid point will make sense physically for this
! scheme.  Indeed, one does not expect this scheme to be relevant
! near boundaries anyhow.  So use for both Bgrid and Cgrid should
! be reasonable, though with the caveat that differences will appear
! near boundaries.  
!
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Greatbatch and Lamb, 1990: On parameterizing vertical mixing of 
! momentum in non eddy-resolving ocean models.  Journal of 
! Physical Oceanography. vol. 20, pages 1634-1637.
! </REFERENCE>
!
! <REFERENCE>
! Greatbatch, 1998: Exploring the relationship betwen eddy-induced 
! transport velocity, vertical momentum transfer, and the isopycnal 
! flux of potential vorticity.  Journal of Physical Oceanography,
! vol. 28, pages 422-432.
! </REFERENCE>
!
! <REFERENCE>
! Aiki, Jacobson, and Yamagata, 2004: Parameterizing ocean eddy transports
! from surface to bottom. Geophysical Research Letters, vol. 31.  
! </REFERENCE>
!
! <REFERENCE>
! Ferreira and Marshall, 2006: Formulation and implementation of 
! a residual-mean ocean circulation model. Ocean Modelling,
! vol. 13, pages 86--107.
! </REFERENCE>
!
! <REFERENCE>
! Danabosaglu and Marshall, 2007: Effects of vertical variations of thickness
! diffusivity in an ocean general circulation model.
! Ocean Modelling, vol. 18, pages 122-141.
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_form_drag_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  Needs to be true in order to use this scheme. Default is false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging.
!  </DATA> 
!  <DATA NAME="verbose_init" TYPE="logical">
!  For verbose initialization information.
!  </DATA> 
!
!  <DATA NAME="use_form_drag_gbatch" TYPE="logical">
!  For using the Greatbatch form drag approach, which places  
!  the contributions from "GM" into the momentum equation.
!
!  use_form_drag_gbatch, is the traditional transformed
!  Eulerian mean approach as per Greatbatch and Lamb (1990)
!  and Greatbatch (1998). In this approach, we modify 
!  visc_cbu so that the the "GM-effects" occur through 
!  vertical viscosity.  
!
!  Since thise scheme is experimental, it is not recommended
!  for general use, so we set default 
!  use_form_drag_gbatch=.false.  
!  </DATA> 
!
!  <DATA NAME="form_drag_gbatch_surf_layer" TYPE="logical">
!  Logical to enable the use of a diabatic layer over which the 
!  eddy-induce velocity is assumed to be constant with depth. 
!  Default form_drag_gbatch_surf_layer=.false. 
!  </DATA> 
!  <DATA NAME="ksurf_blayer_min" TYPE="integer">
!  Minimum number of vertical grid points in the surface turbulent 
!  boundary layer for use with form_drag_gbatch_surf_layer=.true.
!  Default ksurf_blayer_min=3. 
!  </DATA> 
!
!  <DATA NAME="form_drag_gbatch_alpha_f2" TYPE="logical">
!  For use of a vertical viscosity with the form_drag_gbatch that is 
!  equal to visc = form_drag_gbatch_alpha * f^2.  This form of the 
!  vertical viscosity is used by Ferreira and Marshall, 2006.
!  Default form_drag_gbatch_alpha_f2=.false.
!  </DATA> 
!  <DATA NAME="form_drag_gbatch_alpha" TYPE="real"  UNITS="m^2*sec">
!  For use of a vertical viscosity with the form_drag_gbatch that is 
!  equal to visc = form_drag_gbatch_alpha * f^2. 
!  Default form_drag_gbatch_alpha=3e8
!  </DATA> 
!
!  <DATA NAME="form_drag_gbatch_f2overN2" TYPE="logical">
!  To compute vertical viscosity according to visc=kappa*(f/N)**2
!  with kappa=gm-diffusivity, f=Coriolis, and N=buoyancy frequency.
!  This is the form suggested by Greatbatch and Lamb (1990).
!  Default form_drag_gbatch_f2overN2=.false.
!  </DATA> 
!
!  <DATA NAME="form_drag_gbatch_f2overNb2" TYPE="logical">
!  To compute vertical viscosity according to visc=kappa*(f/Nb)**2
!  with kappa=gm-diffusivity, f=Coriolis, and Nb=buoyancy frequency
!  just below the diabatic mixed layer. This is the form suggested 
!  by Danabasoglu and Marshall (2007).
!  Default form_drag_gbatch_f2overNb2=.false.
!  </DATA> 
!
!  <DATA NAME="form_drag_gbatch_f2overNo2" TYPE="logical">
!  To compute vertical viscosity according to visc=kappa*(f/No)**2
!  with kappa=gm-diffusivity, f=Coriolis, and No=constant buoyancy
!  frequency set via namelist. 
!  Default form_drag_gbatch_f2overNo2=.false.
!  </DATA> 
!  <DATA NAME="form_drag_gbatch_No" TYPE="real"  UNITS="sec^-1">
!  To compute vertical viscosity according to visc=kappa*(f/No)**2
!  with kappa=gm-diffusivity, f=Coriolis, and No=constant buoyancy
!  frequency.  Default form_drag_gbatch_No=5e-3
!  </DATA> 
!
!  <DATA NAME="form_drag_gbatch_smooth_N2" TYPE="logical">
!  For vertically smoothing the squared buoyancy frequency for 
!  use in computing the vertical viscosity in the form drag 
!  approach. Default form_drag_gbatch_smooth_N2=.false.
!  </DATA> 
!  <DATA NAME="num_121_passes" TYPE="integer">
!  Number of 1-2-1 passes for vertically smoothing the squared
!  buoyancy frequency. Default num_121_passes=1.
!  </DATA> 
!
!  <DATA NAME="visc_cbu_form_drag_max" TYPE="real" UNITS="m^2/sec">
!  Maximum vertical viscosity used for the form drag contribution 
!  to vertical friction from the Greatbatch TEM approach. 
!  Default visc_cbu_form_drag_max=1m^2/sec. 
!  </DATA> 
!  <DATA NAME="vel_form_drag_max" TYPE="real" UNITS="m/sec">
!  For diagnostic purpuses, maximum form drag eddy induced velocity.
!  Default vel_form_drag_max=1m/sec. 
!  </DATA> 
!  <DATA NAME="N_squared_min" TYPE="real" UNITS="1/sec">
!  Minimum sequared buoyancy frequency (N^2) for use in computing 
!  the vertical viscosity from the Greatbatch form drag scheme. 
!  Default N_squared_min=1e-10.
!  </DATA> 
!
!  <DATA NAME="use_form_drag_aiki" TYPE="logical">
!  For using the Aiki form drag approach.
!  </DATA> 
!
!  <DATA NAME="cprime_aiki" TYPE="real"  UNITS="dimensionless">
!  Dimensionless parameters from the Aiki etal scheme. 
!  Default cprime_aiki=0.3. 
!  </DATA> 
!  <DATA NAME="form_drag_aiki_scale_by_gm" TYPE="logical">
! Compute a dimensionless scale function proportional to 
! the GM-diffusivity, for use with the Aiki etal form 
! drag scheme. Default form_drag_aiki_scale_by_gm=.false.
!  </DATA> 
!
!  <DATA NAME="form_drag_aiki_bottom_layer" TYPE="logical">
!  For implementing the Aiki form drag just in a selected number of 
!  bottom layers.  Will still insist that the scheme conserves 
!  momentum, as a form drag scheme should do.  
!  Default form_drag_aiki_bottom_layer=.false. 
!  </DATA> 
!  <DATA NAME="form_drag_aiki_bottom_klevels" TYPE="integer">
!  Number of klevels above the bottom that we choose to close the 
!  Aiki form drag scheme.  Default form_drag_aiki_bottom_klevels=3.
!  This should ideally be less than the minimum number of klevels 
!  in the model.
!  </DATA> 
!  <DATA NAME="form_drag_aiki_scale_by_gradH" TYPE="logical">
!  For scaling the coefficient used by the form drag scheme with 
!  the bottom slope.  
!  </DATA> 
!  <DATA NAME="form_drag_aiki_gradH_power" TYPE="real">
!  For scaling the coefficient used by the form drag scheme with 
!  the bottom slope raised to the power "form_drag_aiki_gradH_power".
!  </DATA> 
!  <DATA NAME="form_drag_aiki_gradH_max" TYPE="logical">
!  For scaling setting the maximum of the bottom slope
!  for use with scaling the form drag coefficient. 
!  Default form_drag_aiki_gradH_max=.05
!  </DATA> 
!
!</NAMELIST>

use constants_mod,     only: epsln
use diag_manager_mod,  only: register_diag_field, register_static_field
use fms_mod,           only: stdout, stdlog, FATAL, WARNING, NOTE
use fms_mod,           only: write_version_number, open_namelist_file, check_nml_error, close_file
use mpp_mod,           only: input_nml_file, mpp_error
use mpp_domains_mod,   only: mpp_global_max, mpp_update_domains

use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: missing_value, onefourth 
use ocean_parameters_mod, only: rho0r, onehalf, grav
use ocean_types_mod,      only: ocean_domain_type, ocean_grid_type
use ocean_types_mod,      only: ocean_thickness_type, ocean_time_type
use ocean_types_mod,      only: ocean_velocity_type, ocean_options_type
use ocean_types_mod,      only: ocean_density_type, ocean_time_steps_type
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3 
use ocean_workspace_mod,  only: wrk1_v, wrk2_v, wrk3_v
use ocean_workspace_mod,  only: wrk1_v2d, wrk2_v2d
use ocean_workspace_mod,  only: wrk1_2d, wrk2_2d, wrk3_2d, wrk4_2d
use ocean_util_mod,       only: diagnose_2d_u, diagnose_3d_u

implicit none

private 

public compute_visc_form_drag
public form_drag_accel
public ocean_form_drag_init

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

real, dimension(:,:),   allocatable :: gradH_scale     ! (dimensionless) scaling according to grad of topog
real, dimension(:,:,:), allocatable :: agm_diffusivity ! (m^2/sec) GM diffusivity 
real, dimension(:,:,:), allocatable :: N_squared       ! (m^2/sec) GM diffusivity 

! for diagnostics 
logical :: used
integer :: id_form_drag_aiki_gradH_scale=-1
integer :: id_form_drag_aiki_u          =-1
integer :: id_form_drag_aiki_v          =-1
integer :: id_form_drag_aiki_power      =-1
integer :: id_form_drag_aiki_coeff      =-1

integer :: id_visc_cbu_form_drag_u    =-1
integer :: id_visc_cbu_form_drag_v    =-1
integer :: id_form_drag_velocity_u    =-1
integer :: id_form_drag_velocity_v    =-1
integer :: id_form_drag_velocity_fu   =-1
integer :: id_form_drag_velocity_fv   =-1
integer :: id_N_squared               =-1
integer :: id_f2overN2                =-1
integer :: id_f2overNb2               =-1
integer :: id_ksurf_blayer_form_drag  =-1
integer :: id_surface_blayer_form_drag=-1


character(len=128) :: version=&
       '$Id: ocean_form_drag.F90,v 20.0 2013/12/14 00:16:34 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk 

real :: onepcprime_aiki2_r
real :: cprime_aiki2

! nml for Aiki form drag scheme 
logical :: use_form_drag_aiki            =.false.
logical :: form_drag_aiki_bottom_layer   =.false.
logical :: form_drag_aiki_scale_by_gradH =.false.
real    :: form_drag_aiki_gradH_max      = 0.05 
real    :: form_drag_aiki_gradH_power    = 1.0
integer :: form_drag_aiki_bottom_klevels = 3
real    :: cprime_aiki                   = 0.3

! nml for Gbatch form drag scheme 
logical :: use_form_drag_gbatch        =.false.
logical :: form_drag_gbatch_surf_layer =.false.
logical :: form_drag_gbatch_alpha_f2   =.false.
logical :: form_drag_gbatch_f2overN2   =.false.
logical :: form_drag_gbatch_f2overNb2  =.false.
logical :: form_drag_gbatch_f2overNo2  =.false.
logical :: form_drag_gbatch_smooth_N2  =.false.
logical :: form_drag_aiki_scale_by_gm  =.false.

integer :: num_121_passes   = 1 
integer :: ksurf_blayer_min = 3

real    :: form_drag_gbatch_alpha     = 3e8  
real    :: form_drag_gbatch_No        = 5e-3
real    :: visc_cbu_form_drag_max     = 1.0
real    :: vel_form_drag_max          = 1.0
real    :: N_squared_min              = 1e-10
real    :: agm_form_drag              = 600.0
real    :: agm_r                      = 0.0016666667

logical :: verbose_init                  = .true.
logical :: use_this_module               = .false.
logical :: debug_this_module             = .false.
logical :: module_is_initialized         = .FALSE.

namelist /ocean_form_drag_nml/ verbose_init, use_this_module, debug_this_module,   &
                  use_form_drag_aiki, cprime_aiki, form_drag_aiki_bottom_layer,    & 
                  form_drag_aiki_bottom_klevels, form_drag_aiki_scale_by_gradH,    &
                  form_drag_aiki_gradH_max, form_drag_aiki_gradH_power,            &
                  form_drag_aiki_scale_by_gm,                                      &
                  use_form_drag_gbatch, visc_cbu_form_drag_max, vel_form_drag_max, &
                  N_squared_min, agm_form_drag,                                    &
                  form_drag_gbatch_surf_layer, ksurf_blayer_min,                   &
                  form_drag_gbatch_alpha_f2, form_drag_gbatch_alpha,               &
                  form_drag_gbatch_f2overN2, form_drag_gbatch_f2overNb2,           &
                  form_drag_gbatch_smooth_N2, num_121_passes,                      &
                  form_drag_gbatch_f2overNo2, form_drag_gbatch_No

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_form_drag_init">
!
! <DESCRIPTION>
! Initial set up for parameterized pressure form drag.
! </DESCRIPTION>
!
subroutine ocean_form_drag_init(Grid, Domain, Time, Time_steps, Ocean_options, debug)

  type(ocean_grid_type),       intent(in), target   :: Grid
  type(ocean_domain_type),     intent(in), target   :: Domain
  type(ocean_time_type),       intent(in)           :: Time
  type(ocean_time_steps_type), intent(in)           :: Time_steps
  type(ocean_options_type),    intent(inout)        :: Ocean_options
  logical,                     intent(in), optional :: debug

  integer :: io_status, ioun, ierr
  integer :: i,j
  integer :: num_visc_schemes 

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if (module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_form_drag_mod (ocean_form_drag_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number(version, tagname)

  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk

  Dom => Domain
  Grd => Grid

  ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_form_drag_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_form_drag_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_form_drag_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_form_drag_nml')
  call close_file (ioun)
#endif
  write (stdlogunit,ocean_form_drag_nml)
  write (stdoutunit,'(/)') 
  write (stdoutunit,ocean_form_drag_nml)

  ! this array needs allocation even if not using module 
  allocate (agm_diffusivity(isd:ied,jsd:jed,nk))
  agm_diffusivity=0.0


  if(use_this_module) then 
      if(use_form_drag_gbatch .and. use_form_drag_aiki) then 
         Ocean_options%form_drag = 'Used form drag from BOTH Gbatch and Aiki: experimental, so beware.'
         write(stdoutunit,'(a)') &
          '==>Note: ocean_form_drag_mod using BOTH Gbatch and Aiki schemes: experimental, so beware!'  
         write(stdoutunit,'(a)') &
         '==>Note: May choose to set "use_gm_skew=.false." in ocean_neutral_physics to avoid double counting.'
      endif 
      if(use_form_drag_gbatch .and. .not. use_form_drag_aiki) then 
         Ocean_options%form_drag = 'Used form drag from Greatbatch: experimental, so beware.'
         write(stdoutunit,'(a)')  &
         '==>Note: ocean_form_drag_mod using Greatbatch scheme: experimental, so beware'  
         write(stdoutunit,'(a)') &
         '==>Note: May choose to set "use_gm_skew=.false." in ocean_neutral_physics to avoid double counting.'
      endif 
      if(.not. use_form_drag_gbatch .and. use_form_drag_aiki) then 
         Ocean_options%form_drag = 'Used form drag from Aiki etal: experimental, so beware.'
         write(stdoutunit,'(a)')  &
         '==>Note: ocean_form_drag_mod using Aiki scheme: experimental, so beware'  
         write(stdoutunit,'(a)') &
         '==>Note: May choose to set "use_gm_skew=.false." in ocean_neutral_physics to avoid double counting.'
      endif 
      if(.not. use_form_drag_gbatch .and. .not. use_form_drag_aiki) then 
         write(stdoutunit,'(a)')  &
         '==>Note: ocean_form_drag_mod: use_this_module=.true. but use_form_drag_gbatch=.false. &use_form_drag_aiki=.false.'  
         Ocean_options%form_drag = 'Did NOT use parameterized form drag.'
      endif 
  else 
      call mpp_error(NOTE, '==>Note from ocean_form_drag_mod: NOT USING this module')
      Ocean_options%form_drag = 'Did NOT use parameterized form drag.'
      return
  endif

  if (PRESENT(debug) .and. .not. debug_this_module) then
      debug_this_module = debug
  endif 
  if(debug_this_module) then 
     write(stdoutunit,'(a)') '==>Note: running ocean_form_drag_mod with debug_this_module=.true.'  
  endif 

  if(Time_steps%aidif /= 1.0) then 
      call mpp_error(FATAL, &
        '==>Error from ocean_form_drag_mod: aidif=1.0 required to handle form drag implicitly.')
  endif


  if(use_form_drag_gbatch) then 
      num_visc_schemes = 0

      if(form_drag_gbatch_alpha_f2) then 
          write(stdoutunit,'(1x,a)')    &
           '==> Note from ocean_form_drag_mod: Vertical viscosity for TEM method is set to constant*f^2'
          num_visc_schemes = 1 + num_visc_schemes
      endif
      if(form_drag_gbatch_f2overN2) then 
          write(stdoutunit,'(1x,a)')    &
           '==> Note from ocean_form_drag_mod: Vertical viscosity for TEM method is set to kappa*(f/N)^2'
          num_visc_schemes = 1 + num_visc_schemes
      endif
      if(form_drag_gbatch_f2overNb2) then 
          write(stdoutunit,'(1x,a)')    &
           '==> Note from ocean_form_drag_mod: Vertical viscosity for TEM method is set to kappa*(f/Nb)^2'
          num_visc_schemes = 1 + num_visc_schemes
      endif
      if(form_drag_gbatch_f2overNo2) then 
          write(stdoutunit,'(1x,a)')    &
           '==> Note from ocean_form_drag_mod: Vertical viscosity for TEM method is set to kappa*(f/No)^2'
          num_visc_schemes = 1 + num_visc_schemes
      endif

      if (num_visc_schemes == 0) then 
          call mpp_error(WARNING, &
           '==>Warning from ocean_form_drag_init: no method chosen for vertical viscosity. Gbatch form drag will do nothing.')
      endif
      if (num_visc_schemes > 1) then 
          call mpp_error(FATAL, &
           '==>Error from ocean_form_drag_init: can choose just one method to compute vertical viscosity')
      endif

      if(form_drag_gbatch_surf_layer) then 
          write(stdoutunit,'(1x,a)')    &
           '==> Note from ocean_form_drag_mod: Vertical viscosity for TEM method modified in surf turb layer'
      endif

  endif

  ! commonly used constants for Aiki scheme 
  onepcprime_aiki2_r = 1.0/(1.0 + cprime_aiki*cprime_aiki)
  cprime_aiki2       = cprime_aiki*cprime_aiki

  ! squared buoyancy frequency 
  allocate (N_squared(isd:ied,jsd:jed,nk))
  N_squared = 0.0

  ! scaling proportional to bottom topog gradient 
  allocate (gradH_scale(isd:ied,jsd:jed))
  gradH_scale(:,:) = 1.0
  if(form_drag_aiki_scale_by_gradH) then 
      do j=jsd,jed
         do i=isd,ied
            gradH_scale(i,j) = &
            (Grd%gradH(i,j)/(epsln+form_drag_aiki_gradH_max))**form_drag_aiki_gradH_power
            gradH_scale(i,j) = min(gradH_scale(i,j),1.0) 
         enddo
      enddo
  endif


  ! Aiki related diagnostics 
  id_form_drag_aiki_gradH_scale = register_static_field ('ocean_model', 'form_drag_aiki_gradH_scale', &
                   Grd%vel_axes_uv(1:2), 'Topography scaling for Aiki form drag', 'dimensionless',    &
                    missing_value=missing_value, range=(/-1e6,1e6/))
  call diagnose_2d_u(Time, Grd, id_form_drag_aiki_gradH_scale, gradH_scale(:,:))

  id_form_drag_aiki_u = register_diag_field ('ocean_model', 'form_drag_aiki_u', &
                   Grd%vel_axes_u(1:3), Time%model_time,                        &
                   'Aiki i-component form stress', 'N/m2',                      &
                    missing_value=missing_value, range=(/-1e6,1e6/))

  id_form_drag_aiki_v = register_diag_field ('ocean_model', 'form_drag_aiki_v', &
                   Grd%vel_axes_v(1:3), Time%model_time,                        &
                   'Aiki j-component form stress', 'N/m2',                      &
                   missing_value=missing_value, range=(/-1e6,1e6/))

  id_form_drag_aiki_power = register_diag_field ('ocean_model', 'form_drag_aiki_power', &
                       Grd%vel_axes_uv(1:3), Time%model_time,                           &
                       'Power sink from Aiki form drag', 'Watt',                        &
                       missing_value=missing_value, range=(/-1e16,1e16/))

  id_form_drag_aiki_coeff= register_diag_field ('ocean_model', 'form_drag_aiki_coeff', &
                           Grd%vel_axes_wu(1:3), Time%model_time,                      &
                           'dimensionless form scaling coefficient for Aiki scheme',   &
                           'dimensionless', missing_value=missing_value, range=(/-10.0,1e6/))


  ! Greatbatch related diagnostics 
  id_N_squared = register_diag_field ('ocean_model', 'N_sqrd_gbatch', Grd%vel_axes_wu(1:3), &
               Time%model_time, 'N^2 on T-cell bottom for Gbatch form drag', 'sec^-2',      &
               missing_value=missing_value, range=(/-1.e2,1e8/))

  id_f2overN2 = register_diag_field ('ocean_model', 'f2overN2', Grd%vel_axes_wu(1:3), &
               Time%model_time, '(f/N)**2 on U-cell bottom', 'dimensionless',         &
               missing_value=missing_value, range=(/-1.0,1e6/))

  id_f2overNb2= register_diag_field ('ocean_model', 'f2overNb2', Grd%vel_axes_wu(1:2),    &
               Time%model_time, '(f/N)**2 on U-cell at base of surf mixed layer', 'sec^2',&
               missing_value=missing_value, range=(/-1.0,1e6/))

  id_visc_cbu_form_drag_u= register_diag_field ('ocean_model', 'visc_cbu_form_drag_u', &
                 Grd%vel_axes_u(1:3), Time%model_time,                                 &
                 'vertical viscosity from form drag on u-velocity component', 'm^2/s', &
                 missing_value=missing_value, range=(/-10.0,1e6/),                     &
                 standard_name='ocean_vertical_momentum_diffusivity_due_to_form_drag')

  id_visc_cbu_form_drag_v= register_diag_field ('ocean_model', 'visc_cbu_form_drag_v', &
                 Grd%vel_axes_v(1:3), Time%model_time,                                 &
                 'vertical viscosity from form drag on v-velocity component', 'm^2/s', &
                 missing_value=missing_value, range=(/-10.0,1e6/))

  id_form_drag_velocity_u= register_diag_field ('ocean_model', 'form_drag_velocity_u', &
                 Grd%vel_axes_u(1:3), Time%model_time,                                 &
                 'eddy induced i-velocity from Gbatch form drag', 'm/s',               &
                 missing_value=missing_value, range=(/-10.0,10.0/))

  id_form_drag_velocity_v= register_diag_field ('ocean_model', 'form_drag_velocity_v', &
                 Grd%vel_axes_v(1:3), Time%model_time,                                 &
                 'eddy induced j-velocity from Gbatch form drag', 'm/s',               &
                 missing_value=missing_value, range=(/-10.0,10.0/))

  id_form_drag_velocity_fu= register_diag_field ('ocean_model', 'form_drag_velocity_fu', &
                 Grd%vel_axes_uv(1:3), Time%model_time,                                  &
                 'f * eddy induced i-velocity from Gbatch form drag', 'm/s',             &
                 missing_value=missing_value, range=(/-10.0,10.0/))

  id_form_drag_velocity_fv= register_diag_field ('ocean_model', 'form_drag_velocity_fv', &
                 Grd%vel_axes_uv(1:3), Time%model_time,                                  &
                 'f * eddy induced j-velocity from Gbatch form drag', 'm/s',             &
                 missing_value=missing_value, range=(/-10.0,10.0/))

  id_ksurf_blayer_form_drag= register_diag_field ('ocean_model', 'ksurf_blayer_form_drag', &
                 Grd%vel_axes_uv(1:2), Time%model_time,                                    &
                 'k-index at the base of the boundary layer for GBatch form drag',         &
                 'unitless', missing_value=missing_value, range=(/-10.0,1e10/))

  id_surface_blayer_form_drag= register_diag_field ('ocean_model', 'surface_blayer_form_drag', &
                 Grd%vel_axes_uv(1:2), Time%model_time,                                        &
                 'U-cell depth of surface boundary layer used in Gbatch form drag',            &
                 'm', missing_value=missing_value, range=(/-10.0,1e10/))


end subroutine ocean_form_drag_init
! </SUBROUTINE> NAME="ocean_form_drag_init"



!#######################################################################
! <SUBROUTINE NAME="compute_visc_form_drag">
! <DESCRIPTION>
! Compute vertical viscosity arising from parameterized 
! form drag according to the methods from the Greatbatch TEM approach.
! 
! Follow the method of Ferreira and Marshall (2006) 
! for transitioning through the diabatic surface mixed layer. 
! We take the surf_blthick from mixed layer schemes (e.g., KPP), and use 
! this to define the region where the eddy induced velocity has 
! zero vertical shear.  Note that we do not introduce an additional
! transitional layer.  
!
! There are four methods to compute the vertical viscosity within
! the interior, with one required to be chosen:
! 1/ form_drag_gbatch_f2overN2
!    visc = kappa*(f/N)^2, with kappa=gm-diffusivity
! 2/ form_drag_gbatch_f2overNb2
!    visc = kappa*(f/Nb)^2, with Nb=buoyancy freq just below surface mixed layer
! 3/ form_drag_gbatch_f2overNo2
!    visc = kappa*(f/No)^2, with No=constant buoyancy freq 
! 4/ form_drag_gbatch_alpha_f2
!    visc = alpha*f^2, with alpha = constant (units of m^2*sec)
!
! </DESCRIPTION>
!
!
subroutine compute_visc_form_drag(Time, Thickness, Velocity, Dens, &
                                  gm_diffusivity, surface_blt, visc_cbu_form_drag)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_velocity_type),      intent(in)    :: Velocity
  type(ocean_density_type),       intent(in)    :: Dens
  real, dimension(isd:,jsd:,:),   intent(in)    :: gm_diffusivity 
  real, dimension(isd:,jsd:),     intent(in)    :: surface_blt
  real, dimension(isd:,jsd:,:,:), intent(inout) :: visc_cbu_form_drag

  integer :: i,j,k,kb,m,n,tau,kbot
  integer :: ksurf_blayer(isd:ied,jsd:jed)
  real    :: surface_blayerT(isd:ied,jsd:jed)
  real    :: surface_blayerU(isd:ied,jsd:jed)
  real    :: N_squaredU, sign_switch
  real    :: tmp,N2_prev, No_r
  logical :: below_surf_blayer

  if(.not. use_this_module .or. .not. use_form_drag_gbatch) then
    visc_cbu_form_drag(:,:,:,:) = 0.0
    agm_diffusivity(:,:,:)      = gm_diffusivity(:,:,:)
    return
  endif 

  tau = Time%tau

  wrk1(:,:,:)                 = 0.0
  wrk2(:,:,:)                 = 0.0
  wrk3_v(:,:,:,:)             = 0.0
  visc_cbu_form_drag(:,:,:,:) = 0.0
  No_r                        = 1.0/(epsln+form_drag_gbatch_No)
  agm_diffusivity(:,:,:)      = gm_diffusivity(:,:,:)

  ! compute N^2 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           N_squared(i,j,k) = -grav*rho0r*Dens%drhodz_zt(i,j,k)
        enddo
     enddo
  enddo
  if(form_drag_gbatch_smooth_N2) then
      do m=1,num_121_passes
         do j=jsd,jed
            do i=isd,ied
               N2_prev = onefourth*N_squared(i,j,1)
               kbot=Grd%kmt(i,j)
               if(kbot>1) then
                   do k=2,kbot-2
                      tmp = N_squared(i,j,k)
                      N_squared(i,j,k) = N2_prev + onehalf*N_squared(i,j,k) &
                                        + onefourth*N_squared(i,j,k+1)
                      N2_prev = onefourth*tmp
                   enddo
               endif
            enddo
         enddo
      enddo
  endif

  ! squared buoyancy frequency on U-cell
  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec
           if(Grd%umask(i,j,k) > 0.0) then 
               N_squaredU  = onefourth*(N_squared(i,j,k)   + N_squared(i+1,j,k)  &
                                      + N_squared(i,j+1,k) + N_squared(i+1,j+1,k))
               wrk2(i,j,k) = abs(N_squaredU) 
               wrk1(i,j,k) = Grd%f(i,j)**2/(wrk2(i,j,k) + N_squared_min)
           endif
        enddo
     enddo
  enddo


  ! viscosity from Greatbatch and Lamb (1990) 
  if(form_drag_gbatch_f2overN2) then
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               if(Grd%umask(i,j,k) > 0.0) then 
                   wrk3_v(i,j,k,1) = Grd%umask(i,j,k+1) &
                        *min(visc_cbu_form_drag_max,wrk1(i,j,k)*agm_diffusivity(i,j,k))
                   wrk3_v(i,j,k,2) = wrk3_v(i,j,k,1)
                   visc_cbu_form_drag(i,j,k,1) = wrk3_v(i,j,k,1) 
                   visc_cbu_form_drag(i,j,k,2) = wrk3_v(i,j,k,2) 
               endif
            enddo
         enddo
      enddo
  endif

  ! viscosity from Greatbatch and Lamb (1990) with constant N 
  if(form_drag_gbatch_f2overNo2) then
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               if(Grd%umask(i,j,k) > 0.0) then 
                   wrk3_v(i,j,k,1) = Grd%umask(i,j,k+1) &
                    *min(visc_cbu_form_drag_max,agm_diffusivity(i,j,k)*(No_r*Grd%f(i,j))**2)
                   wrk3_v(i,j,k,2) = wrk3_v(i,j,k,1)
                   visc_cbu_form_drag(i,j,k,1) = wrk3_v(i,j,k,1) 
                   visc_cbu_form_drag(i,j,k,2) = wrk3_v(i,j,k,2) 
               endif
            enddo
         enddo
      enddo
  endif

  ! viscosity from Ferreira and Marshall (2006)
  if(form_drag_gbatch_alpha_f2) then 
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               if(Grd%umask(i,j,k) > 0.0) then 
                   visc_cbu_form_drag(i,j,k,1) = form_drag_gbatch_alpha*Grd%f(i,j)**2
                   visc_cbu_form_drag(i,j,k,2) = form_drag_gbatch_alpha*Grd%f(i,j)**2
               endif
            enddo
         enddo
      enddo
  endif

  ! determine some properties of the surface diabatic layer 
  if(form_drag_gbatch_surf_layer .or. form_drag_gbatch_f2overNb2) then

      ksurf_blayer(:,:)    = 1
      surface_blayerT(:,:) = 0.0
      surface_blayerU(:,:) = 0.0
      wrk4_2d(:,:)         = 0.0

      ! compute boundary layer depth on U-cell 
      do j=jsc,jec
         do i=isc,iec
            surface_blayerT(i,j) = surface_blt(i,j)
         enddo
      enddo
      call mpp_update_domains (surface_blayerT(:,:), Dom%domain2d)
      surface_blayerU(:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec
            surface_blayerU(i,j) = onefourth*Grd%umask(i,j,k)      &
                 *(surface_blayerT(i,j)  +surface_blayerT(i+1,j)   &
                  +surface_blayerT(i,j+1)+surface_blayerT(i+1,j+1))
         enddo
      enddo

      ! compute k-level at or just below surface boundary layer. 
      ! initialize ksurf_blayer to 1 rather than 0, to avoid 
      ! accessing the zeroth element of an array.   
      wrk3_2d(:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec

            kb               = Grd%kmt(i,j)
            below_surf_blayer= .false.
            ksurf_blayer(i,j)= ksurf_blayer_min
            ksurf_blayer_loop: do k=ksurf_blayer_min,kb-1
               if(Thickness%depth_zwu(i,j,k) >= surface_blayerU(i,j) &
                  .and. .not. below_surf_blayer) then 
                   ksurf_blayer(i,j) = k
                   below_surf_blayer = .true. 
                   exit ksurf_blayer_loop
               endif
            enddo ksurf_blayer_loop

            ! for cases where blayer reaches to the bottom.
            ! need to give ksurf_blayer one above the bottom so 
            ! to be able to compute a nonzero vertical shear. 
            if(.not. below_surf_blayer) then
                ksurf_blayer(i,j) = max(kb-1,ksurf_blayer_min)
            endif

            ! need a real field to send to diagnostic manager
            wrk3_2d(i,j) = float(ksurf_blayer(i,j))

            ! for the form_drag_gbatch_f2overNb2 scheme, save 
            ! (f/Nb)**2, where Nb=buoyancy frequency one level beneath 
            ! diabatic mixed layer. 
            k = min(nk,1+ksurf_blayer(i,j))
            wrk4_2d(i,j) = Grd%f(i,j)**2/(wrk2(i,j,k) + N_squared_min)

         enddo
      enddo

  endif

  ! viscosity from Danabasoglu and Marshall (2007)
  if(form_drag_gbatch_f2overNb2) then
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               if(Grd%umask(i,j,k) > 0.0) then 
                   wrk3_v(i,j,k,1) = Grd%umask(i,j,k+1) &
                        *min(visc_cbu_form_drag_max,wrk4_2d(i,j)*agm_diffusivity(i,j,k))
                   wrk3_v(i,j,k,2) = wrk3_v(i,j,k,1)
                   visc_cbu_form_drag(i,j,k,1) = wrk3_v(i,j,k,1) 
                   visc_cbu_form_drag(i,j,k,2) = wrk3_v(i,j,k,2) 
               endif
            enddo
         enddo
      enddo
  endif 


  ! modify viscosity in surface boundary layer 
  if(form_drag_gbatch_surf_layer) then

      do n=1,2
         do j=jsc,jec
            do i=isc,iec

               ! compute vertical shear at blayer base
               k               = ksurf_blayer(i,j)
               wrk1_v2d(i,j,n) = Grd%umask(i,j,k+1) &
                    *(Velocity%u(i,j,k,n,tau)-Velocity%u(i,j,k+1,n,tau))/Thickness%dzwu(i,j,k) 

               ! modify the vertical viscosity inside blayer.
               ! use linear tapering function, going to zero as approach ocean top.  
               do k=1,ksurf_blayer(i,j)
                  if(surface_blayerU(i,j) > 0.0) then 
                      wrk3_v(i,j,k,n) = wrk3_v(i,j,ksurf_blayer(i,j),n)                       &
                           *(Thickness%depth_zwu(i,j,k)/surface_blayerU(i,j))*wrk1_v2d(i,j,n) &
                           *Thickness%dzwu(i,j,k)/(Velocity%u(i,j,k,n,tau)-Velocity%u(i,j,k+1,n,tau)+epsln)
                      wrk3_v(i,j,k,n) = Grd%umask(i,j,k+1)*min(visc_cbu_form_drag_max,abs(wrk3_v(i,j,k,n)))
                      visc_cbu_form_drag(i,j,k,n) = wrk3_v(i,j,k,n)
                  endif
               enddo

            enddo
         enddo
      enddo

  endif

  ! send some diagnostics 
  call diagnose_3d_u(Time, Grd, id_N_squared, N_squared(:,:,:))
  call diagnose_2d_u(Time, Grd, id_ksurf_blayer_form_drag, wrk3_2d(:,:))
  call diagnose_2d_u(Time, Grd, id_surface_blayer_form_drag, surface_blayerU(:,:))
  call diagnose_3d_u(Time, Grd, id_f2overN2, wrk1(:,:,:))
  call diagnose_2d_u(Time, Grd, id_f2overNb2, wrk4_2d(:,:))
  call diagnose_3d_u(Time, Grd, id_visc_cbu_form_drag_u, visc_cbu_form_drag(:,:,:,1))
  call diagnose_3d_u(Time, Grd, id_visc_cbu_form_drag_v, visc_cbu_form_drag(:,:,:,2))

  ! diagnose the eddy induced velocity 
  if(id_form_drag_velocity_u  > 0 .or. id_form_drag_velocity_v  > 0 .or. &
     id_form_drag_velocity_fu > 0 .or. id_form_drag_velocity_fv > 0) then 

      wrk1_v(:,:,:,:) = 0.0
      wrk2_v(:,:,:,:) = 0.0
      wrk1_2d(:,:)    = 0.0
      wrk2_2d(:,:)    = 0.0

      do n=1,2

         sign_switch = (-1.0)**(n-1) 

         do j=jsc,jec
            do i=isc,iec
               wrk1_2d(i,j) = 0.0
            enddo
         enddo

         do k=1,nk-1
            do j=jsc,jec
               do i=isc,iec
                  wrk2_2d(i,j) =  Grd%umask(i,j,k+1)*wrk3_v(i,j,k,n)        &
                       *(Velocity%u(i,j,k,n,tau)-Velocity%u(i,j,k+1,n,tau)) &
                       /Thickness%dzwu(i,j,k)
                  wrk2_v(i,j,k,n) = sign_switch*Grd%umask(i,j,k)*(wrk1_2d(i,j)-wrk2_2d(i,j)) &
                       /Thickness%dzu(i,j,k)
                  wrk1_v(i,j,k,n) = wrk2_v(i,j,k,n)/(Grd%f(i,j)+epsln)
                  wrk1_2d(i,j)    = wrk2_2d(i,j)
               enddo
            enddo
         enddo
         do j=jsc,jec
            do i=isc,iec
               wrk2_v(i,j,nk,n) = sign_switch*Grd%umask(i,j,nk)*wrk1_2d(i,j) &
                    /Thickness%dzu(i,j,k)
               wrk1_v(i,j,nk,n) = wrk2_v(i,j,nk,n)/(Grd%f(i,j)+epsln)
            enddo
         enddo
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  if(abs(wrk1_v(i,j,k,n)) > vel_form_drag_max) then 
                      wrk1_v(i,j,k,n) = vel_form_drag_max*sign(1.0,wrk1_v(i,j,k,n))
                  endif
               enddo
            enddo
         enddo

      enddo
      
      call diagnose_3d_u(Time, Grd, id_form_drag_velocity_fu, wrk2_v(:,:,:,1))
      call diagnose_3d_u(Time, Grd, id_form_drag_velocity_fv, wrk2_v(:,:,:,2))
      call diagnose_3d_u(Time, Grd, id_form_drag_velocity_u, wrk1_v(:,:,:,1))
      call diagnose_3d_u(Time, Grd, id_form_drag_velocity_v, wrk1_v(:,:,:,2))

  endif


end subroutine compute_visc_form_drag
! </SUBROUTINE> NAME="compute_visc_form_drag"



!#######################################################################
! <SUBROUTINE NAME="form_drag_accel">
!
! <DESCRIPTION>
!  Compute the rho*dz weighted acceleration from the Aiki etal 
!  form drag scheme.   
! </DESCRIPTION>
!
subroutine form_drag_accel(Time, Thickness, Velocity, energy_analysis_step)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_velocity_type),  intent(inout) :: Velocity
  logical,                    intent(in)    :: energy_analysis_step 

  integer :: i,j,k,n
  integer :: tau, taum1
  integer :: kmu
  real    :: agm_ucell, agm_max, agm_max_r

  if(.not. use_this_module .or. .not. use_form_drag_aiki) then
      if(energy_analysis_step) then 
          do n=1,2
             do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      Velocity%wrkv(i,j,k,n) = 0.0
                   enddo
                enddo
             enddo
          enddo
      endif
    return 
  endif 


  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_form_drag_mod (form_drag_accel): module needs to be initialized')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau
  wrk1  = 1.0

  ! Compute a dimensionless scale function proportional to GM-diffusivity
  if(form_drag_aiki_scale_by_gm) then 
      agm_max   = mpp_global_max(Dom%domain2d, agm_diffusivity(:,:,1))
      if(agm_max > 0.0) then  
          agm_max_r = 1.0/agm_max
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   agm_ucell = onefourth*Grd%umask(i,j,k)                   &
                        *(agm_diffusivity(i,j,k)  +agm_diffusivity(i+1,j,k) &
                         +agm_diffusivity(i,j+1,k)+agm_diffusivity(i+1,j+1,k))
                   wrk1(i,j,k) = agm_ucell*agm_r
                enddo
             enddo
          enddo
      endif
  endif

  wrk1_2d(:,:)    = 0.0
  wrk2_2d(:,:)    = 0.0
  wrk1_v2d(:,:,:) = 0.0
  wrk2_v2d(:,:,:) = 0.0
  wrk1_v(:,:,:,:) = 0.0

  ! pre-calculate some 2d coefficients and the inverse mass per area of a column 
  do j=jsc,jec
     do i=isc,iec
        wrk1_v2d(i,j,1) = abs(Grd%f(i,j))*cprime_aiki*onepcprime_aiki2_r
        wrk1_v2d(i,j,2) = Grd%f(i,j)*cprime_aiki2*onepcprime_aiki2_r
        if(Grd%umask(i,j,1) > 0.0) then 
          wrk1_2d(i,j) = 1.0/Thickness%mass_u(i,j,tau)
        endif 
     enddo
  enddo

  ! vertical integral of equation (11) in Aiki etal paper. 
  ! NOTE: wrk1 is a dimensionless coefficient which 
  ! scales the drag according to GM-diffusivity (from ocean_nphysics_mod)
  ! or slope of the topography (default)
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           kmu = Grd%kmu(i,j)
           if(k <= kmu) then  
               wrk1_v(i,j,k,1) = -wrk1_v2d(i,j,1)*(Velocity%u(i,j,k,1,tau)-Velocity%u(i,j,kmu,1,tau)) &
                                 -wrk1_v2d(i,j,2)*(Velocity%u(i,j,k,2,tau)-Velocity%u(i,j,kmu,2,tau)) 
               wrk1_v(i,j,k,2) = -wrk1_v2d(i,j,1)*(Velocity%u(i,j,k,2,tau)-Velocity%u(i,j,kmu,2,tau)) &
                                 +wrk1_v2d(i,j,2)*(Velocity%u(i,j,k,1,tau)-Velocity%u(i,j,kmu,1,tau)) 
               wrk1_v(i,j,k,1) =  wrk1_v(i,j,k,1)*wrk1(i,j,k)*gradH_scale(i,j)
               wrk1_v(i,j,k,2) =  wrk1_v(i,j,k,2)*wrk1(i,j,k)*gradH_scale(i,j)
           endif
        enddo
     enddo
  enddo

  if(form_drag_aiki_bottom_layer) then 

      ! zero those levels not used 
      do j=jsc,jec
         do i=isc,iec
            do k=1,Grd%kmu(i,j)-(form_drag_aiki_bottom_klevels-1)
               wrk1_v(i,j,k,1) = 0.0
               wrk1_v(i,j,k,2) = 0.0
            enddo
         enddo
      enddo

      ! mass per area in the retained levels 
      wrk1_2d(:,:) = 0.0 
      do j=jsc,jec
         do i=isc,iec
            do k=(Grd%kmu(i,j)-form_drag_aiki_bottom_klevels),Grd%kmu(i,j)
               wrk1_2d(i,j) = wrk1_2d(i,j) + Thickness%rho_dzu(i,j,k,tau)  
            enddo
         enddo
      enddo

      do j=jsc,jec
         do i=isc,iec
            if(wrk1_2d(i,j) > 0.0) then 
                wrk1_2d(i,j) = 1.0/wrk1_2d(i,j)
            endif
         enddo
      enddo

  endif


  do n=1,2
     
     ! vertical integral 
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              wrk2_v2d(i,j,n) = wrk2_v2d(i,j,n) + Thickness%rho_dzu(i,j,k,tau)*wrk1_v(i,j,k,n) 
           enddo
        enddo
     enddo

     ! divide by mass per area of the column 
     do j=jsc,jec
        do i=isc,iec
           wrk2_v2d(i,j,n) = wrk2_v2d(i,j,n)*wrk1_2d(i,j)
        enddo
     enddo

  enddo


  ! ensure SGS form drag has zero vertical integral  
  if(form_drag_aiki_bottom_layer) then 

      ! compute form drag which has zero vertical integral over retained levels 
      do n=1,2
         do j=jsc,jec
            do i=isc,iec
               do k=(Grd%kmu(i,j)-form_drag_aiki_bottom_klevels),Grd%kmu(i,j)
                  wrk1_v(i,j,k,n) = Thickness%rho_dzu(i,j,k,tau)*(wrk1_v(i,j,k,n) - wrk2_v2d(i,j,n))
               enddo
            enddo
         enddo
      enddo

  else 

      ! compute form drag which has zero vertical integral over full column
      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1_v(i,j,k,n) = Thickness%rho_dzu(i,j,k,tau)*(wrk1_v(i,j,k,n) - wrk2_v2d(i,j,n))
               enddo
            enddo
         enddo
      enddo

  endif

  ! energetic analysis 
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

  ! add to acceleration 
  else 
      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%accel(i,j,k,n) = Velocity%accel(i,j,k,n) + wrk1_v(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

      call diagnose_3d_u(Time, Grd, id_form_drag_aiki_u, wrk1_v(isc:iec,jsc:jec,:,1))
      call diagnose_3d_u(Time, Grd, id_form_drag_aiki_v, wrk1_v(isc:iec,jsc:jec,:,2))
      if (id_form_drag_aiki_power > 0) then 
          wrk2(:,:,:) = 0.0
          do n=1,2
             do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      wrk2(i,j,k) = wrk2(i,j,k) + Velocity%u(i,j,k,n,tau)*wrk1_v(i,j,k,n)*Grd%dau(i,j)
                   enddo
                enddo
             enddo
          enddo
          call diagnose_3d_u(Time, Grd, id_form_drag_aiki_power, wrk2(isc:iec,jsc:jec,:))
      endif

      if (id_form_drag_aiki_coeff > 0) then
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   wrk2(i,j,k) = wrk1(i,j,k)*gradH_scale(i,j)
                enddo
             enddo
          enddo
          call diagnose_3d_u(Time, Grd, id_form_drag_aiki_coeff, wrk2(:,:,:))
      endif


  endif  ! endif for energy_analysis_step


end subroutine form_drag_accel
! </SUBROUTINE> NAME="form_drag_accel"



end module ocean_form_drag_mod
      
      





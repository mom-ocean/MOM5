module ocean_passive_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Alistair Adcroft 
!</CONTACT>
!
!<OVERVIEW>
! Set up for idealized passive tracers. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module setups various passive tracer configurations.
! Some idealized initial conditions are provided.  
!
! These tracers are used for various idealized purposes, 
! such as tracing water mass transports, testing 
! advection schemes, etc. 
!
! These passive tracers are always initialized within 
! this module with an idealized profile. 
! No preprocessing step is required. To define a common initial condition
! for all passive tracers use the namelist "ocean_passive_nml"
! 
! This setting can be overwritten for each tracer by the namelist feature of the tracer field_table.
! For example
!
!   "namelists","ocean_mod","ocean_passive/patch"
!   restore = f
!   init_condition = temp_sq_init
!
! init_condition is one of :: 'level', 'wall','patch', 'patch_'klevel, with "klevel" an integer
!                                                       for the k-level that will place the patch. 
!                             'exponential', 'shelfbowl', 'rho_surface', 'temp_sq_init', 'salt_sq_init'
! Default is 'patch'
!
! With rho_surface' a density surface can be selected, the density value 
! is defined in the field table namelist, for example 
!
!   init_surface = 1025
!
! However, if an initial file exists for a passive tracer,
! then ocean_tracer_mod will overwrite the passive tracer
! with the tracer concentration in the initial file. In this 
! way, we can, for example, initialize a passive tracer with 
! some profile that is not readily determined via a simple 
! algorithmic procedure. 
! 
! If restoring of a passive tracer to its initial value is enabled by
! setting in the field table
!
!   restore = t
!
! the initial field is used only to find the grid cells where to restore the 
! passive tracer to the initial tracer field. Restoring is done where the tracer 
! concentration exceeds 0.00001. The inital value of passive tracers with restoring
! is always set to '1'. With the field_table namelist 
!
!   init_value = some_real
!
! this can be changed to another but also constant value. 
!
! All passive tracers in this module are dimensionless and are 
! treated the same internally to this module. However, they can 
! generally have different initial conditions and can use 
! different advection schemes. Indeed, one motivation for 
! developing this module is to test advection schemes, with 
! the same initial condition used for each of the tracers,
! but different advection schemes. In this way we can readily 
! determine the difference between advection schemes on various
! profiles within MOM.  
!
! Sample passive tracers are setup here.  The user can modify  
! code in a straightforward manner to change the number of
! passive tracers and/or the initial profiles.  
! 
! </DESCRIPTION>
!
!<NAMELIST NAME="ocean_passive_nml">
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging the module.
!  </DATA>
!
!  <DATA NAME="common_init_condition" TYPE="character">
!  Default for the tracer initial conditions.
!
!  Options are the following:
!  common_init_condition='level'
!  common_init_condition='wall'
!  common_init_condition='patch'
!  common_init_condition='patch_'klevel, with "klevel" an integer
!                         for the k-level that will place the patch. 
!  common_init_condition='exponential'
!  common_init_condition='shelfbowl'
!  common_init_condition='rho_surface'
!  common_init_condition='temp_sq_init'
!  common_init_condition='salt_sq_init'
!  Default=common_init_condition='patch'
!  </DATA>
!
!  <DATA NAME="layer_value" TYPE="real">
!  Value of tracer concentration within the layer. 
!  Default=1.0.
!  </DATA>
!  <DATA NAME="layer_ztop" TYPE="real">
!  Depth at the top of the tracer layer. 
!  </DATA>
!  <DATA NAME="layer_zbot" TYPE="real">
!  Depth at the bottom of the tracer layer. 
!  </DATA>
!  <DATA NAME="patch_init_klevel_gaussian" TYPE="logical">
!  To initialize on the klevel with a gaussian region.  
!  Default=patch_init_klevel_gaussian=.false.
!  </DATA>
!
!  <DATA NAME="wall_value" TYPE="real">
!  Value of tracer concentration within the wall. 
!  Default=1.0.
!  </DATA>
!  <DATA NAME="wall_ratio_south" TYPE="real">
!  Ratio of the full j-range, northward of which 
!  we place the wall. 
!  </DATA>
!  <DATA NAME="wall_ratio_north" TYPE="real">
!  Ratio of the full j-range, southward of which 
!  we place the wall. 
!  </DATA>
!
!  <DATA NAME="patch_value" TYPE="real">
!  Value of the tracer concentration within a patch. 
!  Default=1.0.
!  </DATA>
!  <DATA NAME="patch_ztop" TYPE="real">
!  Depth at the top of the tracer patch. 
!  </DATA>
!  <DATA NAME="patch_zbot" TYPE="real">
!  Depth at the bottom of the tracer patch. 
!  </DATA>
!  <DATA NAME="patch_ratio1" TYPE="real">
!  For setting position of tracer patch.
!  </DATA>
!  <DATA NAME="patch_ratio2" TYPE="real">
!  For setting position of tracer patch.
!  </DATA>
!
!  <DATA NAME="efold_depth" TYPE="real" UNITS="metre">
!  The efolding depth used for exponential tracer profile. 
!  Default=1000.0.
!  </DATA>
!  <DATA NAME="exponential_value" TYPE="real" UNITS="dimensionless">
!  The tracer value at zero depth when choosing the exponential profile.
!  Default=1.0.
!  </DATA>
!
!</NAMELIST>

use constants_mod,     only: pi
use field_manager_mod, only: fm_string_len, fm_path_name_len, fm_field_name_len
use field_manager_mod, only: fm_get_length, fm_get_value, fm_new_value
use fm_util_mod,       only: fm_util_start_namelist, fm_util_end_namelist
use fm_util_mod,       only: fm_util_set_value, fm_util_get_real, fm_util_get_logical
use fm_util_mod,       only: fm_util_get_string_array, fm_util_get_string
use fm_util_mod,       only: fm_util_check_for_bad_fields
use fms_mod,           only: write_version_number 
use fms_mod,           only: open_namelist_file, check_nml_error, close_file
use mpp_mod,           only: input_nml_file, mpp_min, mpp_max, mpp_error, FATAL, WARNING, NOTE, stdout, stdlog

use ocean_domains_mod,    only: get_local_indices
use ocean_tpm_util_mod,   only: otpm_set_prog_tracer, otpm_set_tracer_package, otpm_set_diag_tracer
use ocean_types_mod,      only: ocean_domain_type, ocean_grid_type
use ocean_types_mod,      only: ocean_options_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_diag_tracer_type

implicit none

private

#include <ocean_memory.h>

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()

character(len=128) :: version=&
       '$Id: ocean_passive.F90,v 20.0 2013/12/14 00:17:02 fms Exp $'
character (len=128) :: tagname = &
       '$Name: tikal $'

public ocean_passive_init
public passive_tracer_init
public ocean_passive_tracer_init
public update_tracer_passive

private layer_init
private wall_init
private patch_init
private patch_init_klevel
private exponential_init
private shelfbowl_init

type passive_type
  real, allocatable , dimension(:,:,:) :: mask
  character(len=fm_field_name_len) :: name
  character(len=128)               :: init_condition
  integer                          :: index, diag_index
  real                             :: init_surface
  real                             :: init_value
  logical                          :: restore
end type passive_type
type(passive_type), dimension(:), allocatable, save :: passive

integer :: index_temp =-1
integer :: index_salt =-1  


! defaults for the passive tracer package 

integer, parameter                          :: num_types = 1
character(len=fm_field_name_len), parameter :: package_name     = 'ocean_passive'
character(len=48), parameter                :: mod_name         = 'ocean_passive_mod'
character(len=fm_string_len), parameter     :: default_restart_file = 'ocean_passive.res.nc'

! logical set according to tracer table 
logical :: use_tracer_sq=.false.
logical :: use_tracer_restore=.false.

! set min/max tracer concentration 
! if solution falls outside this range, model will be 
! gracefully brought down.
real :: t_min=-1e6
real :: t_max= 1e6

! min/max tracer concentration beyond which employ upwind 
! advection if set limit_tracer_range_quick=.true.
! (set in ocean_tracer_advect_nml) and/or 
! horizontal diffusion if limit_tracer_range_neutral=.true. 
! (set in ocean_neutral_physics_nml)
real :: t_min_limit=-.10
real :: t_max_limit=1.0

! number of passive tracers 
integer :: instances

! for setting the passive tracer package 
integer :: package_index

logical :: module_is_initialized=.false.
logical :: use_this_module=.false.

! nml settings 
logical :: debug_this_module = .false.

! options are 'patch', 'patch_'klevel, 'layer', 'wall', 'exponential', and 'shelfbowl'
character(len=32) :: common_init_condition='patch'   

real :: layer_value = 1.0
real :: layer_ztop  = 100.0 ! metres
real :: layer_zbot  = 200.0 ! metres

real :: wall_value       = 1.0
real :: wall_ratio_south = 0.33333 
real :: wall_ratio_north = 0.66666 

real :: patch_value  = 1.0    
real :: patch_ztop   = 0.0    ! metres
real :: patch_zbot   = 200.0  ! metres
real :: patch_ratio1 = 0.33333
real :: patch_ratio2 = 0.66666
logical :: patch_init_klevel_gaussian=.false.

real :: shelfbowl_north = 70.0 ! latitude
real :: shelf_value     = 1.0 

real :: efold_depth       = 1000.0  ! metres 
real :: exponential_value = 1.0

namelist /ocean_passive_nml/ debug_this_module, common_init_condition,                        &
                             layer_value, layer_ztop, layer_zbot,                             &   
                             wall_value, wall_ratio_south, wall_ratio_north,                  &   
                             patch_value, patch_ztop, patch_zbot, patch_ratio1, patch_ratio2, &   
                             patch_init_klevel_gaussian, efold_depth, exponential_value,      &
                             shelfbowl_north, shelf_value
contains


!#######################################################################
! <SUBROUTINE NAME="ocean_passive_init">
!
! <DESCRIPTION>
! Initialize the indices for passive tracer fields. 
! This routine is called by ocean_model.F90.
! </DESCRIPTION>
!
subroutine ocean_passive_init(Domain, Grid, Ocean_options, debug)

  type(ocean_domain_type),  intent(in), target   :: Domain
  type(ocean_grid_type),    intent(in), target   :: Grid
  type(ocean_options_type), intent(inout)        :: Ocean_options
  logical,                  intent(in), optional :: debug
  integer                                        :: ioun, io_status, ierr

! -------------------------------------------------------------------------
! start of tracer package material

  ! local parameters 
  character(len=64), parameter :: sub_name = 'ocean_passive_init'
  character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
  character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

  ! local variables 
  integer                                             :: n
  character(len=fm_field_name_len)                    :: name
  character(len=fm_path_name_len)                     :: path_to_names
  character(len=fm_field_name_len+1)                  :: suffix
  character(len=fm_field_name_len+3)                  :: long_suffix
  character(len=256)                                  :: caller_str
  character(len=fm_string_len), pointer, dimension(:) :: good_list
   
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then
      call mpp_error(FATAL, '==>Error in ocean_passive_mod (ocean_passive_init): module initialized.')
  endif
  module_is_initialized=.true.

  ! initialize the ocean passive tracer package and set some defaults. 
  ! each of these defaults can be altered via entries to the field_table. 
  package_index = otpm_set_tracer_package(package_name,                             &
     units = 'dimensionless', min_tracer=t_min, max_tracer=t_max,                   &
     min_tracer_limit=t_min_limit, max_tracer_limit=t_max_limit,                    &
     restart_file=default_restart_file,                                             &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')',                       &
     conversion=1.0, offset=0.0, min_range=-10.0, max_range=100.0,                  &
     flux_units = 'dimensionless', min_flux_range=-1.0e+16, max_flux_range=1.0e+16, &
     vert_adv_scheme='mdppm', horiz_adv_scheme='mdppm', psom_limit=.true.)

  ! check for number of entries in the field_table for passive tracers.  
  path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
  instances = fm_get_length(path_to_names)
  if (instances < 0) then  
    call mpp_error( &
    FATAL, trim(error_header) // ' Could not get number of passive tracer instances.')
  endif  

  ! determine whether to use this module. 
  write (stdoutunit,*) ' '
  if(instances==0) then  
    write (stdoutunit,*) trim(note_header), ' No instances of passive tracers in field_table.'
    use_this_module=.false.
    write(stdoutunit,'(a)') ' '
    write(stdoutunit,'(a)') &
    '==>Note: NOT running with idealized passive tracers.'  
    call mpp_error(NOTE, '==>Note: ocean_passive_mod: NOT using idealized passive tracer module.')
    Ocean_options%passive_tracers = 'Did NOT use the idealized passive tracer module.'
    return 
  else  
    write (stdoutunit,*) trim(note_header), ' ', instances, ' instances of passive tracers in field_table.'
    use_this_module=.true.
    Ocean_options%passive_tracers = 'Ran with idealized passive tracer module.'
  endif  

  ! allocate the passive array
  allocate (passive(instances))

  ! loop over the names of the passive tracers, saving them into the passive array
  do n=1,instances  
     if (fm_get_value(path_to_names, name, index = n)) then  
         passive(n)%name = name
     else 
         write (name,*) n
         call mpp_error(FATAL, trim(error_header) //        &
              ' Bad field name for index ' // trim(name))
     endif
  enddo

  ! determine the tracer name for this instance, 
  ! and set the prog tracer through otpm_set_prog_tracer
  do n=1,instances  
     name = passive(n)%name
     if (name(1:1) .eq. '_') then  
         suffix = ' '
         long_suffix = ' '
     else  
         suffix = '_' // name
         long_suffix = ' (' // trim(name) // ')'
     endif
     passive(n)%index = otpm_set_prog_tracer('passive' // trim(suffix), package_name, &
          longname = 'passive' // trim(long_suffix),                                  &
          caller = trim(mod_name)//'('//trim(sub_name)//')')
  enddo

  ! add the package name to the list of good namelists; used later for consistency check
  if (fm_new_value('/ocean_mod/GOOD/good_namelists', package_name, append = .true.) .le. 0) then 
    call mpp_error(FATAL, trim(error_header) //                           &
        ' Could not add ' // trim(package_name) // ' to "good_namelists" list')
  endif 
 
  ! set defaults for the available namelist parameters 
  caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'
  do n=1,instances  
     call fm_util_start_namelist(package_name, passive(n)%name, caller=caller_str, &
                                 no_overwrite=.true., check=.true.)
     call fm_util_set_value('init_condition', 'patch')
     call fm_util_set_value('init_surface', 1030.0)
     call fm_util_set_value('init_value', 1.0)
     call fm_util_set_value('restore', .false.)
     call fm_util_end_namelist(package_name, passive(n)%name, check = .true., caller = caller_str)
  enddo

  ! check for errors in number of fields in the namelists for this package
  good_list => fm_util_get_string_array('/ocean_mod/GOOD/namelists/' // trim(package_name) // '/good_values',   &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  if (associated(good_list)) then  
      call fm_util_check_for_bad_fields('/ocean_mod/namelists/' // trim(package_name), good_list,       &
           caller = trim(mod_name) // '(' // trim(sub_name) // ')')
      deallocate(good_list)
  else  
      call mpp_error(FATAL,trim(error_header) // ' Empty "' // trim(package_name) // '" list')
  endif

! end of tracer package material
! -------------------------------------------------------------------------

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif

  Dom => Domain
  Grd => Grid

  call write_version_number( version, tagname )

  ! provide for namelist override of defaults
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_passive_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_passive_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_passive_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status, 'ocean_passive_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_passive_nml)
  write (stdlogunit,ocean_passive_nml)

  if (PRESENT(debug) .and. .not. debug_this_module) then
      debug_this_module = debug
  endif
  if(debug_this_module) then 
      write(stdoutunit,'(a)') '==>Note: running ocean_passive with debug_this_module=.true.'  
  endif

  ! set default initial condition according to setting in ocean_passive_nml
  do n=1,instances
     passive(n)%init_condition = trim(common_init_condition)
  enddo

  ! namelist information for each passive tracer from the field_table. 
  ! in particular, can override the common_init_condition so that each 
  ! tracer can, for example, have a distinct initial condition.  
  caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'
  do n=1,instances
     call fm_util_start_namelist(package_name, passive(n)%name, caller = caller_str)
     passive(n)%init_condition =  fm_util_get_string ('init_condition', scalar = .true.)
     passive(n)%init_surface   =  fm_util_get_real('init_surface', scalar = .true.)
     passive(n)%init_value     =  fm_util_get_real('init_value', scalar = .true.)
     passive(n)%restore        =  fm_util_get_logical('restore', scalar = .true.)
     call fm_util_end_namelist(package_name, passive(n)%name, caller = caller_str)
  enddo
  do n=1,instances
     if (passive(n)%restore ) then
        name = passive(n)%name
        if (name(1:1) .eq. '_') then  
           suffix = ' '
           long_suffix = ' '
        else  
           suffix = '_' // name
           long_suffix = ' (' // trim(name) // ')'
        endif
        passive(n)%diag_index = otpm_set_diag_tracer('passive' // trim(suffix) // 'mask',         &
           longname = 'passive' // trim(long_suffix), restart_file='ocean_passive_masks.res.nc',  &
           min_range=-1.0, max_range=2.0, const_init_tracer=.true.,const_init_value=0.0,          &
           caller = caller_str )
        use_tracer_restore=.true.
     endif
  enddo


end subroutine ocean_passive_init
! </SUBROUTINE> NAME="ocean_passive_init"


!#######################################################################
! <SUBROUTINE NAME="passive_tracer_init">
!
! <DESCRIPTION>
! Initialize profiles for the passive tracers.
! This routine is called by ocean_tracer.F90.
! </DESCRIPTION>
!
subroutine passive_tracer_init(time_init, Tracer, initialize_as_a_passive_tracer)

  logical,                      intent(in)    :: time_init
  type(ocean_prog_tracer_type), intent(inout) :: Tracer
  logical,                      intent(inout) :: initialize_as_a_passive_tracer

  integer :: i,j,k,n

  integer :: stdoutunit 
  stdoutunit=stdout() 

  initialize_as_a_passive_tracer=.false.
  if(.not. use_this_module) return 
  if(.not. time_init) return 

  do n=1,instances 

        if(trim(Tracer%name)=='passive_'//trim(passive(n)%name)) then 
            initialize_as_a_passive_tracer=.true.

            if(trim(passive(n)%init_condition)=='layer') then 
                call layer_init(Tracer%field(:,:,:,1))
                write(stdoutunit,'(a)') '==>From passive_tracer_init: initialized tracer '&
                     //trim(passive(n)%name)//' with "layer" profile.' 
            elseif(trim(passive(n)%init_condition)=='wall') then  
                call wall_init(Tracer%field(:,:,:,1))
                write(stdoutunit,'(a)') '==>From passive_tracer_init: initialized tracer '&
                     //trim(passive(n)%name)//' with "wall" profile.' 
            elseif(trim(passive(n)%init_condition)=='patch') then  
                call patch_init(Tracer%field(:,:,:,1))
                write(stdoutunit,'(a)') '==>From passive_tracer_init: initialized tracer '&
                     //trim(passive(n)%name)//' with "patch" profile.' 
            elseif(trim(passive(n)%init_condition)=='patch_klevel') then  
                call patch_init_klevel(Tracer%field(:,:,:,1),trim(Tracer%name))
                write(stdoutunit,'(a)') '==>From passive_tracer_init: initialized tracer '&
                     //trim(passive(n)%name)//' with "patch_klevel" profile, w/ patch on single k-level.' 
            elseif(trim(passive(n)%init_condition)=='exponential') then  
                call exponential_init(Tracer%field(:,:,:,1))
                write(stdoutunit,'(a)') '==>From passive_tracer_init: initialized tracer '&
                     //trim(passive(n)%name)//' with "exponential" profile.' 
            elseif(trim(passive(n)%init_condition)=='shelfbowl') then  
                call shelfbowl_init(Tracer%field(:,:,:,1))
                write(stdoutunit,'(a)') '==>From passive_tracer_init: initialized tracer '&
                     //trim(passive(n)%name)//' with "shelfbowl" profile.' 
            elseif(     trim(passive(n)%init_condition(1:4))=='temp') then 
                write(stdoutunit,'(a)') &
                  '==>From passive_tracer_init: tracer to be initialized after initialize temp.'
            elseif(     trim(passive(n)%init_condition(1:3))=='rho') then 
                write(stdoutunit,'(a)') &
                  '==>From passive_tracer_init: tracer to be initialized after initialize density.'
            elseif(     trim(passive(n)%init_condition(1:12))=='temp_sq_init') then 
                write(stdoutunit,'(a)') &
                  '==>From passive_tracer_init: tracer to be initialized after initialize temp.'
            elseif(     trim(passive(n)%init_condition(1:12))=='salt_sq_init') then 
                write(stdoutunit,'(a)') &
                  '==>From passive_tracer_init: tracer to be initialized after initialize salinity.'
            else
                write(stdoutunit,'(a)') &
                '==>Error in passive_tracer_init: tracer '//trim(Tracer%name)// ' is not initialized.'
                call mpp_error(FATAL, &
                 '==>Error in passive_tracer_init: passive tracer not initialized.')
            endif 
            do k=1,nk
               do j=jsd,jed
                  do i=isd,ied
                     Tracer%field(i,j,k,2) = Tracer%field(i,j,k,1)
                     Tracer%field(i,j,k,3) = Tracer%field(i,j,k,1)
                  enddo
               enddo
            enddo

        endif

  enddo


end subroutine passive_tracer_init
! </SUBROUTINE> NAME="passive_tracer_init"

!#######################################################################
! <SUBROUTINE NAME="ocean_passive_tracer_init">
!
! <DESCRIPTION>
! Initialize restoring to initial values for the passive tracers.
!
! or
! Initialize profiles for the passive tracers which are tagged with 
! specific values of the temperature field or potential density field.
!
! or
! Initialize profiles for the passive tracers which are defined 
! as the temperature or salinity squared. These tracers are used
! for diagnosing the level of spurious mixing associated with 
! PSOM advection.  
! </DESCRIPTION>
!
subroutine ocean_passive_tracer_init(T_prog, T_diag, rho, &
                                     time_init, tau, num_prog_tracers, i_temp, i_salt)

  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_diag_tracer_type), intent(inout) :: T_diag(:)
  real, dimension(isd:,jsd:,:), intent(in)    :: rho
  logical,                      intent(in)    :: time_init
  integer,                      intent(in)    :: tau, num_prog_tracers, i_temp, i_salt

  integer :: i,j,k,n,nt
  real    :: isotherm, isorho

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return 
  if(.not. time_init) return 
  
  index_temp = i_temp
  index_salt = i_salt
  
  do nt=1,num_prog_tracers 
     do n=1,instances 
        
        if(trim(T_prog(nt)%name)=='passive_'//trim(passive(n)%name)) then 

           if( passive(n)%restore ) then 

              T_diag(passive(n)%diag_index)%field = 0.
              do k=1,nk
                do j=jsd,jed
                  do i=isd,ied
                     if (T_prog(nt)%field(i,j,k,tau) > 0.00001) then 
                        ! set a mask to 1 if restoring is required
                        T_diag(passive(n)%diag_index)%field(i,j,k) = 1.
                        T_prog(nt)%field(i,j,k,1) = passive(n)%init_value
                        T_prog(nt)%field(i,j,k,2) = passive(n)%init_value
                        T_prog(nt)%field(i,j,k,3) = passive(n)%init_value
                     endif
                  enddo
                enddo
              enddo
              write(stdoutunit,'(a)') &
                '==>From passive_tracer_init: restore '//trim(T_prog(nt)%name)// ' to initial conditions.'
           else
              write(stdoutunit,'(a)') &
                '==>From passive_tracer_init: do not restore '//trim(T_prog(nt)%name)// ' to initial conditions.'
           endif 


           ! Initialize profiles for the passive tracers which are tagged with 
           ! specific values of the temperature field or potential density field.
           if(trim(passive(n)%init_condition)=='temp_surface') then 
              isotherm=passive(n)%init_surface
              call surface_init(isotherm, 'temp', T_prog(index_temp)%field(:,:,:,tau), & 
                   T_prog(nt)%field(:,:,:,1))
              write(stdoutunit,'(a,a,a,f12.3)')                                        &
                '==>From passive_tracer_surface_init: passive tracer ',T_prog(nt)%name,&
                ' initialized to isotherm (C) = ',isotherm
              write(stdoutunit,'(a)') ' '
              do k=1,nk
                do j=jsd,jed
                   do i=isd,ied
                      T_prog(nt)%field(i,j,k,2) = T_prog(nt)%field(i,j,k,1)
                      T_prog(nt)%field(i,j,k,3) = T_prog(nt)%field(i,j,k,1)
                   enddo
                enddo
             enddo
          endif
           if(trim(passive(n)%init_condition)=='rho_surface') then 
             isorho=passive(n)%init_surface
             call surface_init(isorho, 'rho', rho, T_prog(nt)%field(:,:,:,1))
             write(stdoutunit,'(a,a,a,f12.3)')                                        &
                '==>From passive_tracer_surface_init: passive tracer ',T_prog(nt)%name,&
                ' initialized to neutral rho (kg/m^3) = ',isorho
             write(stdoutunit,'(a)') ' '
             do k=1,nk
               do j=jsd,jed
                  do i=isd,ied
                     T_prog(nt)%field(i,j,k,2) = T_prog(nt)%field(i,j,k,1)
                     T_prog(nt)%field(i,j,k,3) = T_prog(nt)%field(i,j,k,1)
                  enddo
               enddo
             enddo
           endif

           ! Initialize profiles for the passive tracers which are defined 
           ! as the temperature or salinity squared. These tracers are used
           ! for diagnosing the level of spurious mixing associated with 
           ! PSOM advection.  
           if(trim(passive(n)%init_condition)=='temp_sq_init') then 
               do k=1,nk
                  do j=jsd,jed
                     do i=isd,ied
                        T_prog(nt)%field(i,j,k,1)   = T_prog(index_temp)%field(i,j,k,tau)**2
                        T_prog(nt)%field(i,j,k,2)   = T_prog(nt)%field(i,j,k,1)
                        T_prog(nt)%field(i,j,k,3)   = T_prog(nt)%field(i,j,k,1)
                        T_prog(nt)%units            = '(deg C)^2'
                     enddo
                  enddo
               enddo
               write(stdoutunit,'(a,a,a)')                                    &
                    '==>From ocean_tracer_sq_init: tracer ',T_prog(nt)%name,&
                    ' initialized to squared temp (degC^2)'
               write(stdoutunit,'(a)') ' '
               use_tracer_sq=.true.
           endif

           if(trim(passive(n)%init_condition)=='salt_sq_init') then 
               do k=1,nk
                  do j=jsd,jed
                     do i=isd,ied
                        T_prog(nt)%field(i,j,k,1)   = T_prog(index_salt)%field(i,j,k,tau)**2
                        T_prog(nt)%field(i,j,k,2)   = T_prog(nt)%field(i,j,k,1)
                        T_prog(nt)%field(i,j,k,3)   = T_prog(nt)%field(i,j,k,1)
                        T_prog(nt)%units            = '(psu)^2'
                     enddo
                  enddo
               enddo
               write(stdoutunit,'(a,a,a)')                                    &
                    '==>From ocean_tracer_sq_init: tracer ',T_prog(nt)%name,&
                    ' initialized to squared salinity (psu^2)'
               write(stdoutunit,'(a)') ' '
               use_tracer_sq=.true.
           endif

        endif
     enddo    !!n=1,instances 
  enddo       !!nt=1,num_prog_tracers


end subroutine ocean_passive_tracer_init
! </SUBROUTINE> NAME="ocean_passive_tracer_init"


!#######################################################################
! <SUBROUTINE NAME="surface_init">
!
! <DESCRIPTION>
! Initialize passive tracer according to a particular iso-surface.
! </DESCRIPTION>
!
subroutine surface_init(isovalue, surface_name, surface_field, passive_field)

  real, intent(in)                            :: isovalue
  character(len=*), intent(in)                :: surface_name
  real, intent(in)   , dimension(isd:,jsd:,:) :: surface_field
  real, intent(inout), dimension(isd:,jsd:,:) :: passive_field

  integer :: i, j, k

  passive_field(:,:,:) = 0.0

  if(trim(surface_name)=='rho') then 
      do k=1,nk-1
         do j=jsd,jed
            do i=isd,ied
               if(surface_field(i,j,k) <= isovalue .and. surface_field(i,j,k+1) > isovalue) then 
                   passive_field(i,j,k) = Grd%tmask(i,j,k+1)
               endif
            enddo
         enddo
      enddo
  elseif(trim(surface_name)=='temp') then
      do k=1,nk-1
         do j=jsd,jed
            do i=isd,ied
               if(surface_field(i,j,k) >= isovalue .and. surface_field(i,j,k+1) < isovalue) then 
                   passive_field(i,j,k) = Grd%tmask(i,j,k+1)
               endif
            enddo
         enddo
      enddo
  endif

end subroutine surface_init
! </SUBROUTINE> NAME="surface_init"


!#######################################################################
! <SUBROUTINE NAME="layer_init">
!
! <DESCRIPTION>
!
! Initialize tracer inside a depth layer to have a value, and
! zero outside the layer. 
!
! </DESCRIPTION>
!
subroutine layer_init(field)

  real, intent(inout), dimension(isd:,jsd:,:) :: field 

  integer :: i, j, k

  field(:,:,:) = 0.0

  do k=1,nk
     if(Grd%zt(k) >=layer_ztop .and. Grd%zt(k) <= layer_zbot) then 
         do j=jsd,jed
            do i=isd,ied
               field(i,j,k) = layer_value
            enddo
         enddo
     endif
  enddo

end subroutine layer_init
! </SUBROUTINE> NAME="layer_init"


!#######################################################################
! <SUBROUTINE NAME="wall_init">
!
! <DESCRIPTION>
!
! Initialize tracer inside an (i,k) wall to have a value, and
! zero outside the wall. 
!
! </DESCRIPTION>
!
subroutine wall_init(field)

  real, intent(inout), dimension(isd:,jsd:,:) :: field 

  integer :: i, j, k
  integer :: jsouth, jnorth

  jsouth = int(wall_ratio_south*Grd%nj)
  jnorth = int(wall_ratio_north*Grd%nj)

  field(:,:,:) = 0.0

  do j=jsd,jed
     if(j+Dom%joff >= jsouth .and. j+Dom%joff <= jnorth) then                 
         do k=1,nk
            do i=isd,ied
               field(i,j,k) = wall_value
            enddo
         enddo
     endif
  enddo

end subroutine wall_init
! </SUBROUTINE> NAME="wall_init"


!#######################################################################
! <SUBROUTINE NAME="patch_init">
!
! <DESCRIPTION>
!
! Initialize tracer with simple shapes based on vertical layer:
!  Level k=1, square pill box
!  Level k=2, circular pill box (cylinder)
!  Level k=3, circular cone
!  Level k=4, cosine bell
!  Level k=5, Gaussian bell
!  Levels k>5, square patch based on i,j coordinates (original "patch").
!
! </DESCRIPTION>
!
subroutine patch_init(field)

  real, intent(inout), dimension(isd:,jsd:,:) :: field 

  integer :: i, j, k
  integer :: iwest, ieast
  integer :: jsouth, jnorth
  real    :: x_min,x_max,y_min,y_max
  real    :: x_mid,y_mid,xnd,ynd,scl,rr

  integer :: stdoutunit 
  stdoutunit=stdout() 

  x_min = minval( Grd%xu(:,:) )
  call mpp_min(x_min)
  x_max = maxval( Grd%xu(:,:) )
  call mpp_max(x_max)

  y_min = minval( Grd%yu(:,:) )
  call mpp_min(y_min)
  y_max = maxval( Grd%yu(:,:) )
  call mpp_max(y_max)

  scl = 0.5 * ( patch_ratio2 - patch_ratio1 )
  x_mid = x_min + 0.5 * ( patch_ratio1 + patch_ratio2 ) * ( x_max - x_min )
  y_mid = y_min + 0.5 * ( patch_ratio1 + patch_ratio2 ) * ( y_max - y_min )

  field(:,:,:) = 0.0

  iwest  = int(patch_ratio1*Grd%ni)
  ieast  = int(patch_ratio2*Grd%ni)
  jsouth = int(patch_ratio1*Grd%nj)
  jnorth = int(patch_ratio2*Grd%nj)

  ! Discontinuous square box for all k-levels (levels k=1,5 over-written below)
  do k=1,nk
     if(Grd%zt(k) >=patch_ztop .and. Grd%zt(k) <= patch_zbot) then
         do j=jsd,jed
            do i=isd,ied
               if(i+Dom%ioff  >= iwest  .and. i+Dom%ioff <= ieast   .and. &
                  j+Dom%joff  >= jsouth .and. j+Dom%joff <= jnorth ) then
                   field(i,j,k) = patch_value
               endif
            enddo
         enddo
     endif
  enddo

  write(stdoutunit,*) 'passive_init: ',iwest,ieast,jsouth,jnorth
  write(stdoutunit,*) 'passive_init: ',x_min,x_max
  write(stdoutunit,*) 'passive_init: ',y_min,y_max
  write(stdoutunit,*) 'passive_init: ',x_mid,y_mid
  write(stdoutunit,*) 'passive_init: ',scl

  ! k=1: Discontinuous square box
  k=1
  if (nk.ge.k) then
      do j=jsd,jed
         do i=isd,ied
            xnd = 1.-abs(Grd%xt(i,j)-x_mid)/(x_max-x_min)/scl
            ynd = 1.-abs(Grd%yt(i,j)-y_mid)/(y_max-y_min)/scl
            field(i,j,k) = patch_value * max( 0., sign( 1., xnd ) )  &
                                       * max( 0., sign( 1., ynd ) )
         enddo
      enddo
  endif ! k

  ! k=2: Discontinuous circular box (cylinder)
  k=2
  if (nk.ge.k) then
      do j=jsd,jed
         do i=isd,ied
            xnd = (Grd%xt(i,j)-x_mid)/(x_max-x_min)
            ynd = (Grd%yt(i,j)-y_mid)/(y_max-y_min)
            rr  = 1.-sqrt( xnd**2  + ynd**2 )/scl
            field(i,j,k) = patch_value * max( 0., sign( 1., rr ) )
         enddo
      enddo
  endif ! k

  ! k=3: Circular cone
  k=3
  if (nk.ge.k) then
      do j=jsd,jed
         do i=isd,ied
            xnd = (Grd%xt(i,j)-x_mid)/(x_max-x_min)
            ynd = (Grd%yt(i,j)-y_mid)/(y_max-y_min)
            rr  = sqrt( xnd**2  + ynd**2 )
            field(i,j,k) = patch_value * max( 0., 1.-0.5*rr/scl )
         enddo
      enddo
  endif ! k

  ! k=4: Circular cosine bell
  k=4
  if (nk.ge.k) then
      do j=jsd,jed
         do i=isd,ied
            xnd = (Grd%xt(i,j)-x_mid)/(x_max-x_min)
            ynd = (Grd%yt(i,j)-y_mid)/(y_max-y_min)
            rr  = sqrt( xnd**2  + ynd**2 )
            field(i,j,k) = patch_value * 0.5 * ( 1.0+cos(2.*pi*min(0.5,0.25*rr/scl)) )
         enddo
      enddo
  endif ! k

  ! k=5: Circular Gaussian (with 2.5 sigma inside patch)
  k=5
  if (nk.ge.k) then
      do j=jsd,jed
         do i=isd,ied
            xnd = (Grd%xt(i,j)-x_mid)/(x_max-x_min)
            ynd = (Grd%yt(i,j)-y_mid)/(y_max-y_min)
            rr  = sqrt( xnd**2  + ynd**2 )
            field(i,j,k) = patch_value * ( max(0.001, 1.001*exp(-(1.*rr/scl)**2) ) - 0.001 )
         enddo
      enddo
  endif ! k


end subroutine patch_init
! </SUBROUTINE> NAME="patch_init"


!#######################################################################
! <SUBROUTINE NAME="patch_init_klevel">
!
! <DESCRIPTION>
! Initialize tracer with gaussian patch or constant on a level,
! both on a single k-level 
! </DESCRIPTION>
!
subroutine patch_init_klevel(field, name)

  real, dimension(isd:,jsd:,:), intent(inout) :: field 
  character(len=*),             intent(in)    :: name 

  integer :: i, j, n
  integer :: start, strlen, klevel, digit  
  real    :: x_min,x_max,y_min,y_max
  real    :: x_mid,y_mid,xnd,ynd,scl,rr
  logical :: fatal_flag=.false.

  integer :: stdoutunit 
  stdoutunit=stdout() 

  x_min = minval( Grd%xu(:,:) )
  call mpp_min(x_min)
  x_max = maxval( Grd%xu(:,:) )
  call mpp_max(x_max)

  y_min = minval( Grd%yu(:,:) )
  call mpp_min(y_min)
  y_max = maxval( Grd%yu(:,:) )
  call mpp_max(y_max)

  scl = 0.5 * ( patch_ratio2 - patch_ratio1 )
  x_mid = x_min + 0.5 * ( patch_ratio1 + patch_ratio2 ) * ( x_max - x_min )
  y_mid = y_min + 0.5 * ( patch_ratio1 + patch_ratio2 ) * ( y_max - y_min )

  ! determine the k-level that is to be initialized
  ! Algorithm from Zhi.Liang 
  ! start=15 is specific to "passive_patch_" having 14 characters, 
  ! so the 15th character in the name of the tracer begins the klevel.  
  ! Example: passive_patch_123. 
  ! do n=15,17
  ! when n=1, klevel=1       = 1
  ! when n=2, klevel=1*10+2  = 12
  ! when n=3, klevel=12*10+3 = 123  
  start  = 15
  strlen = len_trim(name)
  klevel = 0
  fatal_flag=.false. 
  
  write(stdoutunit,'(a)') ' Initializing passive_patch_klevel for tracer name = ',name 
  do n=start,strlen
     digit = ichar(name(n:n)) - ichar("0")
     if(digit > 9 .or. digit < 0) then 
         write(*,*)'Error: digit =  ',digit 
         fatal_flag=.true.
     endif
     klevel = klevel*10 + digit
  end do

  if(fatal_flag) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_passive_mod (patch_init_klevel): did not find a digit for klevel.')
  endif 

  ! initialize the patch at k=klevel as a circular Gaussian 
  ! (with 2.5 sigma inside patch)
  if(patch_init_klevel_gaussian) then 
      field(:,:,:) = 0.0
      do j=jsd,jed
         do i=isd,ied
            xnd = (Grd%xt(i,j)-x_mid)/(x_max-x_min)
            ynd = (Grd%yt(i,j)-y_mid)/(y_max-y_min)
            rr  = sqrt( xnd**2  + ynd**2 )
            field(i,j,klevel) = Grd%tmask(i,j,klevel) * patch_value &
                                *( max(0.001, 1.001*exp(-(1.*rr/scl)**2) ) - 0.001 )
         enddo
      enddo
  else
      do j=jsd,jed
         do i=isd,ied
            field(i,j,klevel) = patch_value*Grd%tmask(i,j,klevel)
         enddo
      enddo
  endif


end subroutine patch_init_klevel
! </SUBROUTINE> NAME="patch_init_klevel"


!#######################################################################
! <SUBROUTINE NAME="exponential_init">
!
! <DESCRIPTION>
!
! Initialize tracer with a vertical exponential profile. 
!
! </DESCRIPTION>
!
subroutine exponential_init(field)

  real, intent(inout), dimension(isd:,jsd:,:) :: field 

  integer :: i, j, k

  field(:,:,:) = 0.0

  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           field(i,j,k) = exponential_value*exp(-Grd%zt(k)/efold_depth)
        enddo
     enddo
  enddo

end subroutine exponential_init
! </SUBROUTINE> NAME="exponential_init"


!#######################################################################
! <SUBROUTINE NAME="shelfbowl_init">
!
! <DESCRIPTION>
!
! Initialize tracer as in the shelfbowl topography of use for studying
! idealized overflow problems, as in Winton etal (1998). 
!
! </DESCRIPTION>
!
subroutine shelfbowl_init(field)

  real, intent(inout), dimension(isd:,jsd:,:) :: field 

  integer :: i, j, k
  
  field(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           field(i,j,k) = 0.5*shelf_value*(1.0 + 0.5*(1.0-tanh(0.5*(Grd%yt(i,j)-shelfbowl_north))) )
        enddo
     enddo
  enddo

end subroutine shelfbowl_init
! </SUBROUTINE> NAME="shelfbowl_init"


!#######################################################################
! <SUBROUTINE NAME="update_tracer_passive">
!
! <DESCRIPTION>
! Update the squared tracer.
! Restore passive tracers.
! </DESCRIPTION>
!
subroutine update_tracer_passive(num_prog_tracers, T_prog, T_diag, dtts, taup1)

  integer,                      intent(in)    :: num_prog_tracers
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_diag_tracer_type), intent(inout) :: T_diag(:)
  real                        , intent(in)    :: dtts
  integer                     , intent(in)    :: taup1

  integer :: i,j,k,n,nt
  real    :: mask, temp_t

  if(.not. use_this_module) return 

  if(use_tracer_sq) then

     do nt=1,num_prog_tracers 

        if(trim(T_prog(nt)%name)=='passive_temp_sq') then 
            do k=1,nk
               do j=jsd,jed
                  do i=isd,ied
                     T_prog(nt)%field(i,j,k,1) = T_prog(index_temp)%field(i,j,k,1)**2
                     T_prog(nt)%field(i,j,k,2) = T_prog(index_temp)%field(i,j,k,2)**2
                     T_prog(nt)%field(i,j,k,3) = T_prog(index_temp)%field(i,j,k,3)**2
                  enddo
               enddo
            enddo
        endif
        if(trim(T_prog(nt)%name)=='passive_salt_sq') then 
            do k=1,nk
               do j=jsd,jed
                  do i=isd,ied
                     T_prog(nt)%field(i,j,k,1) = T_prog(index_salt)%field(i,j,k,1)**2
                     T_prog(nt)%field(i,j,k,2) = T_prog(index_salt)%field(i,j,k,2)**2
                     T_prog(nt)%field(i,j,k,3) = T_prog(index_salt)%field(i,j,k,3)**2
                  enddo
               enddo
            enddo
        endif 
     enddo
  endif 

  if(use_tracer_restore) then

     do nt=1,num_prog_tracers 
       
        do n=1,instances 

           if( trim(T_prog(nt)%name)=='passive_'//trim(passive(n)%name) .and. passive(n)%restore ) then 
              do k=1,nk
                do j=jsd,jed
                  do i=isd,ied

                     mask = T_diag(passive(n)%diag_index)%field(i,j,k)
                     temp_t = T_prog(nt)%field(i,j,k,taup1)
                     T_prog(nt)%field(i,j,k,taup1) = mask * passive(n)%init_value &
                          + (1.- mask) * T_prog(nt)%field(i,j,k,taup1)

                     ! store the tendency as source for computing tracer conservation diagnostics
                     T_prog(nt)%source(i,j,k) = T_prog(nt)%source(i,j,k) &
                                              + (T_prog(nt)%field(i,j,k,taup1) - temp_t)/dtts

                  enddo
                enddo
              enddo
           endif

        enddo   ! n
     enddo      ! nt

  endif 

end subroutine update_tracer_passive
! </SUBROUTINE> NAME="update_tracer_passive"


end module ocean_passive_mod

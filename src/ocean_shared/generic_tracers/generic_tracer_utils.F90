!----------------------------------------------------------------
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Niki Zadeh 
! </CONTACT>
! 
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> William Cooke
! </REVIEWER>
!
! <OVERVIEW>
! g_tracer_utils module consists of core utility subroutines 
! to be used by all generic tracer modules.
! These include the lowest level functions for adding, 
! allocating memory, and record keeping of individual 
! generic tracers irrespective of their physical/chemical nature.
! </OVERVIEW>
!----------------------------------------------------------------


#include <fms_platform.h>
module g_tracer_utils

  use coupler_types_mod, only: coupler_2d_bc_type, ind_flux, ind_deltap, ind_kw
  use coupler_types_mod, only: ind_alpha, ind_csurf, ind_sc_no
  use FMS_coupler_util,  only: extract_coupler_values, set_coupler_values
  use atmos_ocean_fluxes_mod, only: aof_set_coupler_flux
  use mpp_mod,           only: mpp_error, NOTE, WARNING, FATAL
  use mpp_mod,           only: mpp_pe, mpp_root_pe
  use time_manager_mod,  only: time_type
  use diag_manager_mod,  only: register_diag_field, send_data 


  use field_manager_mod, only: fm_string_len, fm_path_name_len, fm_new_list, fm_change_list, fm_get_value

  use fms_mod,           only: stdout


  implicit none ; private
!-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: generic_tracer_utils.F90,v 20.0 2013/12/14 00:18:12 fms Exp $'
  character(len=128) :: tag = '$Name: tikal $'
!-----------------------------------------------------------------------

  character(len=48), parameter :: mod_name = 'g_tracer_utils'

  ! <DESCRIPTION>
  ! Public types:
  !
  ! Each generic tracer node is an instant of a FORTRAN type with the following member variables.
  ! These member fields are supposed to uniquely define an individual tracer.
  ! One such type shall be instantiated for EACH individual tracer.
  ! <PRE>
  !type g_tracer_type
  !   !A pointer to the next node in the list for the current "linked-list implementation".
  !   type(g_tracer_type), pointer :: next => NULL()  
  !         
  !   !A unique index (for the possible future "array implementation")
  !   integer :: index              
  !
  !   ! Tracer name, descriptive name, package that instantiates it 
  !   character(len=64) :: name, longname, package_name
  ! 
  !   ! Units of measurement for its field and its flux
  !   character(len=64) :: units, flux_units
  !
  !   ! Tracer concentration field in space (and time)
  !   ! MOM keeps the field at 3 time levels, hence 4D.
  !   real, _ALLOCATABLE, dimension(:,:,:,:):: field  _NULL
  !
  !   ! Surface flux, surface gas flux, deltap and kw
  !   real, _ALLOCATABLE, dimension(:,:)    :: stf    _NULL
  ! 
  !   real, _ALLOCATABLE, dimension(:,:)    :: deltap    _NULL
  ! 
  !   real, _ALLOCATABLE, dimension(:,:)    :: kw    _NULL
  ! 
  !   ! Bottom  flux
  !   real, _ALLOCATABLE, dimension(:,:)    :: btf    _NULL
  ! 
  !   ! Bottom  reservoir flux
  !   real, _ALLOCATABLE, dimension(:,:)    :: btm_reservoir    _NULL 
  !
  !   ! Tracer concentration in river runoff
  !   real, _ALLOCATABLE, dimension(:,:)    :: trunoff _NULL 
  !
  !   ! Runoff flux of tracer
  !   real, _ALLOCATABLE, dimension(:,:)    :: runoff_tracer_flux _NULL 
  !
  !   ! Wet deposition flux of tracer
  !   real, _ALLOCATABLE, dimension(:,:)    :: wetdep _NULL 
  !
  !   ! Dry deposition flux of tracer
  !   real, _ALLOCATABLE, dimension(:,:)    :: drydep _NULL 
  !
  !   ! Tracer saturation, alpha, and schmidt number 
  !   real, _ALLOCATABLE, dimension(:,:)    :: csurf  _NULL 
  !
  !   real, _ALLOCATABLE, dimension(:,:)    :: alpha  _NULL 
  !
  !   real, _ALLOCATABLE, dimension(:,:)    :: sc_no  _NULL 
  !
  !   ! An 3D field for vertical movement, esp. for zooplankton, ... 
  !   real, _ALLOCATABLE, dimension(:,:,:)  :: vmove  _NULL
  !   ! An 3D field for random vertical movement, esp. for zooplankton, ... 
  !   real, _ALLOCATABLE, dimension(:,:,:)  :: vdiff  _NULL

  !   ! An 3D field for implicit vertical diffusion
  !   real, _ALLOCATABLE, dimension(:,:,:)  :: vdiffuse_impl  _NULL

  !   ! An auxiliary 3D field for keeping model dependent change tendencies, ... 
  !   real, _ALLOCATABLE, dimension(:,:,:)  :: tendency  _NULL
  !
  !   ! IDs for using diag_manager tools
  !   integer :: diag_id_field=-1, diag_id_stf=-1, diag_id_stf_gas=-1, diag_id_deltap=-1, diag_id_kw=-1, diag_id_trunoff=-1
  !   integer :: diag_id_alpha=-1, diag_id_csurf=-1, diag_id_sc_no=-1, diag_id_aux=-1
  !
  !  ! Tracer Initial concentration if constant everywhere
  !   real    :: const_init_value = 0.0
  !
  !   ! Tracer Sinking rate
  !   real    :: sink_rate   = 0.0
  !
  !   ! Logical switches
  !   logical :: prog        = .false. !Is this a prognostic (.true.) or diagnostic (.false.) tracer?
  !   logical :: move_vertical = .false. ! Enable allocation of fields for active vertical movement
  !   logical :: diff_vertical = .false. ! Enable allocation of fields for random active vertical movement
  !   logical :: flux_gas    = .false. !Is there a gas flux to atmosphere?
  !   logical :: flux_runoff = .false. !Is there a river flux?
  !   logical :: flux_wetdep = .false. !Is there a wet deposition?
  !   logical :: flux_drydep = .false. !Is there a dry deposition?
  !   logical :: flux_bottom = .false. !Is there a flux through bottom?
  !
  !   ! Flux identifiers to be set by aof_set_coupler_flux()
  !   integer :: flux_gas_ind    = -1  
  !   integer :: flux_runoff_ind = -1
  !   integer :: flux_wetdep_ind = -1
  !   integer :: flux_drydep_ind = -1
  !
  !end type g_tracer_type
  !
  ! </PRE>
  ! 
  !
  ! </DESCRIPTION>

  type g_tracer_type
     !A pointer to the next node in the list for the current "linked-list implementation".
     type(g_tracer_type), pointer :: next => NULL()  

     !A unique index (for the possible future "array implementation")
     integer :: index              

     ! Tracer name, descriptive name, package that instantiates it 
     character(len=fm_string_len) :: name, longname, alias, package_name

     ! Tracer molecular wt
     real :: flux_gas_molwt

     ! Tracer flux names recognized by component models (OCN, LND, ICE, ATM) 
     character(len=fm_string_len) :: flux_gas_name, flux_gas_type, flux_runoff_name, flux_wetdep_name, flux_drydep_name
     real, _ALLOCATABLE, dimension(:) :: flux_param, flux_gas_param

     ! IN and OUT (restart) files
     character(len=fm_string_len) :: ice_restart_file, ocean_restart_file
     character(len=fm_string_len) :: flux_gas_restart_file 
     ! Units of measurement for its field and its flux
     character(len=fm_string_len) :: units, flux_units

     ! Tracer concentration field in space (and time)
     ! MOM keeps the prognostic tracer fields at 3 time levels, hence 4D.
     real, pointer, dimension(:,:,:,:):: field  => NULL()
     !The following pointer is intended to point to prognostic tracer field in MOM. Do not allocate!
     real, pointer,      dimension(:,:,:,:):: field4d_ptr => NULL()
     !The following pointer is intended to point to diagnostic tracer field in MOM. Do not allocate!
     real, pointer,      dimension(:,:,:)  :: field3d_ptr => NULL() 
     ! Define a 3-d field pointer so as to retain the lower
     ! and upper bounds for the 3-d version of g_tracer_get_pointer
     ! for the field option
     real, pointer,      dimension(:,:,:)  :: field_3d => NULL()

     ! Surface flux, surface flux of gas, deltap and kw
     real, _ALLOCATABLE, dimension(:,:)    :: stf    _NULL

     real, _ALLOCATABLE, dimension(:,:)    :: stf_gas    _NULL

     real, _ALLOCATABLE, dimension(:,:)    :: deltap    _NULL

     real, _ALLOCATABLE, dimension(:,:)    :: kw    _NULL

     ! Bottom  flux
     real, _ALLOCATABLE, dimension(:,:)    :: btf    _NULL

     ! Bottom  reservoir flux
     real, _ALLOCATABLE, dimension(:,:)    :: btm_reservoir    _NULL 

     ! Tracer concentration in river runoff
     real, _ALLOCATABLE, dimension(:,:)    :: trunoff _NULL 

     ! Runoff flux of tracer
     real, _ALLOCATABLE, dimension(:,:)    :: runoff_tracer_flux _NULL 

     ! Wet deposition flux of tracer
     real, _ALLOCATABLE, dimension(:,:)    :: wetdep _NULL 

     ! Dry deposition flux of tracer
     real, _ALLOCATABLE, dimension(:,:)    :: drydep _NULL 

     ! Tracer saturation, alpha and schmidt number 
     real, _ALLOCATABLE, dimension(:,:)    :: csurf  _NULL 

     real, _ALLOCATABLE, dimension(:,:)    :: alpha  _NULL 

     real, _ALLOCATABLE, dimension(:,:)    :: sc_no  _NULL 

     ! An 3D field for vertical movement, esp. for zooplankton, ... 
     real, _ALLOCATABLE, dimension(:,:,:)  :: vmove  _NULL

     ! An 3D field for random vertical movement, esp. for zooplankton, ... 
     real, _ALLOCATABLE, dimension(:,:,:)  :: vdiff  _NULL

     ! An 3D field for implicit vertical diffusion
     real, _ALLOCATABLE, dimension(:,:,:)  :: vdiffuse_impl  _NULL

     ! An auxiliary 3D field for keeping model dependent change tendencies, ... 
     real, pointer, dimension(:,:,:)  :: tendency  => NULL()

     ! IDs for using diag_manager tools
     integer :: diag_id_field=-1, diag_id_stf=-1, diag_id_stf_gas=-1, diag_id_deltap=-1, diag_id_kw=-1, diag_id_trunoff=-1
     integer :: diag_id_alpha=-1, diag_id_csurf=-1, diag_id_sc_no=-1, diag_id_aux=-1
     integer :: diag_id_btf=-1,diag_id_btm=-1, diag_id_vmove=-1, diag_id_vdiff=-1
     integer :: diag_id_vdiffuse_impl = -1, diag_id_tendency = -1, diag_id_field_taup1 = -1

     ! Tracer Initial concentration if constant everywhere
     real    :: const_init_value = 0.0
     real    :: initial_value = 0.0
     ! Tracer Sinking rate
     real    :: sink_rate   = 0.0

     ! Logical switches
     logical :: prog        = .false. !Is this a prognostic (.true.) or diagnostic (.false.) tracer?
     logical :: move_vertical = .false. ! Enable allocation of fields for active vertical movement
     logical :: diff_vertical = .false. ! Enable allocation of fields for random active vertical movement
     logical :: flux_gas    = .false. !Is there a gas flux to atmosphere?
     logical :: flux_runoff = .false. !Is there a river flux?
     logical :: flux_wetdep = .false. !Is there a wet deposition?
     logical :: flux_drydep = .false. !Is there a dry deposition?
     logical :: flux_bottom = .false. !Is there a flux through bottom?
     logical :: has_btm_reservoir = .false. !Is there a flux bottom reservoir?

     ! Flux identifiers to be set by aof_set_coupler_flux()
     integer :: flux_gas_ind    = -1  
     integer :: flux_runoff_ind = -1
     integer :: flux_wetdep_ind = -1
     integer :: flux_drydep_ind = -1

  end type g_tracer_type


  type g_diag_type
     !A pointer to the next node in the list for the current "linked-list implementation".
     type(g_diag_type), pointer :: next => NULL()
  
     integer :: diag_id = -1
     character(len=fm_string_len) :: name, longname, package_name, units
     !Diagnostic axes 
     integer :: axes(3)
     type(time_type) :: init_time
     real :: missing_value = -1.0e+10
     integer :: Z_diag = 0
     real, pointer, dimension(:,:,:) :: field_ptr 
  end type g_diag_type

  ! <DESCRIPTION>
  ! Public types:
  !
  ! The following type fields are common to ALL generic tracers and hence has to be instantiated only once:
  ! </DESCRIPTION>
  type g_tracer_common
     !Domain extents
     integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk

     !Number of time levels 
     integer :: ntau

     !Diagnostic axes 
     integer :: axes(3)

     !Initial time used for diagnostics (all tracers are instantiated at the same time).
     type(time_type)        :: init_time

     !Grid mask
     real, _ALLOCATABLE, dimension(:,:,:) :: grid_tmask  _NULL !nnz: Make this a pointer, needs to be "target" in models

     !Grid bottom index
     integer, _ALLOCATABLE, dimension(:,:):: grid_kmt    _NULL

     !coast mask
     integer, _ALLOCATABLE, dimension(:,:):: grid_mask_coast    _NULL

     ! IN and OUT (restart) files
     character(len=fm_string_len) :: ice_restart_file, ocean_restart_file
  end type g_tracer_common

  !Keep the state of this common type for ALL tracers
  type(g_tracer_common), target, save :: g_tracer_com



  ! <DESCRIPTION>
  ! Public interfaces:
  ! </DESCRIPTION>
  public :: g_tracer_type
  public :: g_tracer_find
  public :: g_tracer_add
  public :: g_tracer_init
  public :: g_tracer_flux_init
  public :: g_tracer_column_int
  public :: g_tracer_flux_at_depth
  public :: g_tracer_add_param
  public :: g_tracer_set_values
  public :: g_tracer_get_values
  public :: g_tracer_get_pointer
  public :: g_tracer_get_common
  public :: g_tracer_set_common
  public :: g_tracer_set_files
  public :: g_tracer_coupler_set
  public :: g_tracer_coupler_get
  public :: g_tracer_send_diag
  public :: g_tracer_diag
  public :: g_tracer_get_name
  public :: g_tracer_get_alias
  public :: g_tracer_get_next
  public :: g_tracer_register_diag
  public :: g_tracer_is_prog
  public :: g_tracer_vertdiff_G
  public :: g_tracer_vertdiff_M
  public :: g_tracer_start_param_list
  public :: g_tracer_end_param_list
  public :: g_diag_type
  public :: g_diag_field_add
  public :: g_tracer_set_pointer

  ! <INTERFACE NAME="g_tracer_add_param">
  !  <OVERVIEW>
  !   Add a new parameter for the generic tracer package
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine is used to add a new parameter by the calling tracer package.
  !   It provides a mechanism for parameter overwrite through the field_table.
  !For each tracer package there is a field called namelists and there
  !the parameters can be modified from their value set by this method. 
  !E.g., we may have the following in the field_table
  !
  !   "namelists","ocean_mod","generic_topaz"
  !   init = t
  !    /
  !This will overwrite the parameter topaz%init to be .true. at the run time 
  !even though generic_topaz package had in the code 
  !<TT>call g_tracer_add_param('init', topaz%init, .false. )</TT>
  ! 
  !   For the parameters overwrite mechanism to work all calls
  !   for adding new parameters (refer to description for subroutine g_tracer_add_param)
  !   should happen between a <TT>call g_tracer_start_param_list(package_name)</TT>
  !   and a <TT>call g_tracer_end_param_list(package_name)</TT> 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call g_tracer_add_param(param_name, param_variable, param_value )
  !  </TEMPLATE>
  !  <IN NAME="param_name" TYPE="character(len=fm_string_len)">
  !   Name of the  parameter (e.g., "init")
  !  </IN>
  !  <IN NAME="param_variable" TYPE="integer or logical or real">
  !   Variable to contain the  parameter (e.g., "topaz%init")
  !  </IN>
  !  <IN NAME="param_value" TYPE="integer or logical or real">
  !   Value of the  parameter (e.g., ".true.")
  !  </IN>
  ! </INTERFACE>
  interface g_tracer_add_param
     module procedure g_tracer_add_param_real
     module procedure g_tracer_add_param_logical
     module procedure g_tracer_add_param_integer
     module procedure g_tracer_add_param_string      
  end interface

  interface g_tracer_set_pointer
    module procedure g_tracer_set_pointer_3d
    module procedure g_tracer_set_pointer_4d
  end interface g_tracer_set_pointer

  ! <INTERFACE NAME="g_tracer_set_values">
  !  <OVERVIEW>
  !   Set the values of various (array) memebers of the tracer node g_tracer_type
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This function is overloaded to set the values of the following member variables
  !4D arrays:   'field'
  !3D arrays:   'field' , 'tendency'
  !2D arrays:   'alpha','csurf','sc_no','stf','stf_gas','deltap','kw','btf','btm_reservoir','trunoff','runoff_tracer_flux','drydep','wetdep'
  !1D values:   'field','tendency','stf','stf_gas','deltap','kw','btf','btm_reservoir','trunoff','runoff_tracer_flux','drydep','wetdep','btm_reservoir','sink_rate'
  !
  !In case of 1D values all element arrays are set to the particular value.
  !
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call g_tracer_set_values(tracer_list,tracer_name,field_name, array_in ,isd,jsd)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type),    pointer">
  !   Pointer to the head of the generic tracer list.
  !  </IN>
  !  <IN NAME="tracer_name" TYPE="character(len=*)">
  !   Name of the particular tracer.
  !  </IN>
  !  <IN NAME="field_name" TYPE="character(len=*)">
  !   String associated with the member array, one of the following:
  !   'field','tendency','stf','stf_gas','deltap','kw','btf','btm_reservoir','trunoff','runoff_tracer_flux','drydep','wetdep','btm_reservoir','sink_rate'
  !   So the result of this call is tracer%field_name = array_in for the tracer called tracer_name
  !  </IN>
  !  <IN NAME="value" TYPE="real OR real(isd:,jsd:) OR real(isd:,jsd:,:) OR real(isd:,jsd:,:,:)">
  !   Overloaded based on the dimension of argument array_in. 
  !  </IN>
  !  <IN NAME="isd,jsd" TYPE="integer">
  !   Lower bound of the domain for argument array_in
  !  </IN>
  ! </INTERFACE>
  interface g_tracer_set_values
     module procedure g_tracer_set_real
     module procedure g_tracer_set_2D
     module procedure g_tracer_set_3D
     module procedure g_tracer_set_4D
  end interface

  ! <INTERFACE NAME="g_tracer_get_values">
  !  <OVERVIEW>
  !   Reverse of interface g_tracer_set_values for getting the tracer member arrays  in the argument value.
  !  </OVERVIEW>
  !  <TEMPLATE>
  !   call g_tracer_get_values(tracer_list,tracer_name,field_name, array_out ,isd,jsd)
  !  </TEMPLATE>
  !  <DESCRIPTION>
  !   This means "get the values of array  %field_name for tracer tracer_name and put them in argument array_out".
  !  </DESCRIPTION>
  ! </INTERFACE>
  !
  interface g_tracer_get_values
     module procedure g_tracer_get_4D_val
     module procedure g_tracer_get_3D_val
     module procedure g_tracer_get_2D_val
     module procedure g_tracer_get_real
     module procedure g_tracer_get_string
  end interface

  ! <INTERFACE NAME="g_tracer_get_pointer">
  !  <OVERVIEW>
  !   Return the pointer to the requested field of a particular tracer
  !  </OVERVIEW>
  !  <TEMPLATE>
  !       call g_tracer_get_pointer(tracer_list,tracer_name,field_name, array_ptr)
  !  </TEMPLATE>
  !  <DESCRIPTION>
  !   This means "get the pointer of array  %field_name for tracer tracer_name  in argument array_ptr".
  !  </DESCRIPTION>
  ! </INTERFACE>

  interface g_tracer_get_pointer
     module procedure g_tracer_get_4D
     module procedure g_tracer_get_3D
     module procedure g_tracer_get_2D
  end interface

contains

  ! <SUBROUTINE NAME="g_tracer_start_param_list">
  !  <OVERVIEW>
  !   Mark the start of adding new parameters for a package
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   For the parameters override mechanism to work all calls
  !   for adding new parameters (refer to description for subroutine g_tracer_add_param)
  !   should happen between a <TT>call g_tracer_start_param_list(package_name)</TT>
  !   and a <TT>call g_tracer_end_param_list(package_name)</TT> 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call g_tracer_start_param_list(package_name)
  !  </TEMPLATE>
  !  <IN NAME="package_name" TYPE="character(len=fm_string_len)">
  !   Name of the generic tracer package that is adding the parameters (e.g., "generic_cfc")
  !  </IN>
  ! </SUBROUTINE>

  subroutine g_tracer_start_param_list(package_name)
    character(len=fm_string_len), intent(in) :: package_name
    character(len=fm_string_len), parameter  :: sub_name = 'g_tracer_start_param_list'
    character(len=fm_string_len) :: list_path
    integer                      :: list_index

    list_path = '/ocean_mod/namelists/' // trim(package_name) // '/'

    list_index = fm_new_list(list_path)
    if (list_index .le. 0) then  !{
       call mpp_error(FATAL, trim(sub_name) // ' Could not make  the new list' // list_path)
    endif  !}
    
    if (.not. fm_change_list(list_path)) then  !{
       call mpp_error(FATAL, trim(sub_name) // ' Could not change to the new list' // list_path)
    endif  !}

  end subroutine g_tracer_start_param_list

  ! <SUBROUTINE NAME="g_tracer_end_param_list">
  !  <OVERVIEW>
  !   Mark the start of adding new parameters for a package
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   For the parameters override mechanism to work all calls
  !   for adding new parameters (refer to description for subroutine g_tracer_add_param)
  !   should happen between a <TT>call g_tracer_start_param_list(package_name)</TT>
  !   and a <TT>call g_tracer_end_param_list(package_name)</TT> 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call g_tracer_end_param_list(package_name)
  !  </TEMPLATE>
  !  <IN NAME="package_name" TYPE="character(len=fm_string_len)">
  !   Name of the generic tracer package that is adding the parameters (e.g., "generic_cfc")
  !  </IN>
  ! </SUBROUTINE>

  subroutine g_tracer_end_param_list(package_name)
    character(len=fm_string_len) :: package_name

  end subroutine g_tracer_end_param_list

  !Overload interface g_tracer_add_param for real parameter
  subroutine g_tracer_add_param_real(name, var,  value)
    character(len=*), intent(in)  :: name
    real,             intent(in)  :: value
    real,             intent(out) :: var

    real :: x

    ! Need to save "value" since if "var" and "value" are the same
    ! variable, and "name" does not exist, then "var/value" will be
    ! set to 0 in the fm_get_value routine, and "var" cannot then be
    ! set to the supplied default value

    x = value

    if(.NOT. fm_get_value(name, var))  var = x

  end subroutine g_tracer_add_param_real

  !Overload interface g_tracer_add_param for logical parameter
  subroutine g_tracer_add_param_logical(name, var,  value)
    character(len=*), intent(in)  :: name
    logical,          intent(in)  :: value
    logical,          intent(out) :: var

    logical :: x

    ! Need to save "value" since if "var" and "value" are the same
    ! variable, and "name" does not exist, then "var/value" will be
    ! set to false in the fm_get_value routine, and "var" cannot then be
    ! set to the supplied default value

    x = value

    if(.NOT. fm_get_value(name, var))  var = x

  end subroutine g_tracer_add_param_logical

  !Overload interface g_tracer_add_param for integer parameter
  subroutine g_tracer_add_param_integer(name, var,  value)
    character(len=*), intent(in)  :: name
    integer,          intent(in)  :: value
    integer,          intent(out) :: var

    real :: x

    ! Need to save "value" since if "var" and "value" are the same
    ! variable, and "name" does not exist, then "var/value" will be
    ! set to 0 in the fm_get_value routine, and "var" cannot then be
    ! set to the supplied default value

    x = value

    if(.NOT. fm_get_value(name, var))  var = x

  end subroutine g_tracer_add_param_integer

  !Overload interface g_tracer_add_param for string parameter
  subroutine g_tracer_add_param_string(name, var,  value)
    character(len=*), intent(in)  :: name
    character(len=*), intent(in)  :: value
    character(len=*), intent(out) :: var

    character(len=fm_string_len) :: x

    ! Need to save "value" since if "var" and "value" are the same
    ! variable, and "name" does not exist, then "var/value" will be
    ! set to '' in the fm_get_value routine, and "var" cannot then be
    ! set to the supplied default value

    x = value

    if(.NOT. fm_get_value(name, var))  var = x

  end subroutine g_tracer_add_param_string


  ! <SUBROUTINE NAME="g_tracer_add">
  !  <OVERVIEW>
  !   Add a new tracer (node) at the top of the list of generic tracers
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    This subroutine call adds an individual new tracer to the growing list of generic tracers.
  !    It then allocates all the necessary arrays for using this tracer in the Ocean model that requested it.
  !    The information passed into this subroutine should be enough to fully describe the individual tracer
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !  call g_tracer_add(tracer_list,&
  !       package    = 'generic_topaz',&
  !       name       = 'g_dic',               &
  !       longname   = 'g Dissolved Inorganic Carbon', &
  !       units      = 'mol/kg',            &
  !       prog       = .true.,              &
  !       flux_gas       = .true.,                      &
  !       flux_gas_name  = 'co2_flux',                  &
  !       flux_gas_molwt = WTMCO2,                      &
  !       flux_gas_param = (/ 9.36e-07, 9.7561e-06 /),  &
  !       flux_runoff    = .true.,          &
  !       flux_param     = (/12.011e-03  /),  &
  !       flux_bottom    = .true.,          &
  !       ice_restart_file = topaz%ice_restart_file)
  !
  !  </TEMPLATE>
  !  <IN NAME="node_ptr" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head node of the tracer list. This is also going to be the pointer to the node being added after the call.
  !  </IN>
  !  <IN NAME="package" TYPE="character(len=*)">
  !   Name of tracer package adding this node.
  !  </IN>
  !  <IN NAME="name" TYPE="character(len=*)">
  !   Name of this tracer.
  !  </IN>
  !  <IN NAME="longname" TYPE="character(len=*)">
  !   Descriptive name of this tracer..
  !  </IN>
  !  <IN NAME="units" TYPE="character(len=*)">
  !   Concentration units (units of array %field).
  !  </IN>
  !  <IN NAME="prog" TYPE="logical">
  !   .true. for prognastic , .false. for diagnostic tracer.
  !  </IN>
  !
  !  OPTIONAL arguments begin:
  !
  !  <IN NAME="const_init_value" TYPE="real">
  !   Initial value of concenteration if constant.
  !  </IN>
  !  <IN NAME="flux_gas" TYPE="logical">
  !   .true. if there is gas flux exchange with atmos.
  !  </IN>
  !  <IN NAME="flux_gas_name" TYPE="character(len=*)">
  !   Name of the atmospheric tracer to exchange flux with (if flux_gas=.true.).
  !  </IN>
  !  <IN NAME="flux_runoff" TYPE="logical">
  !   .true. if there is runoff flux.
  !  </IN>
  !  <IN NAME="flux_wetdep" TYPE="logical">
  !   .true. if there is wetdep flux.
  !  </IN>
  !  <IN NAME="flux_drydep" TYPE="logical">
  !   .true. if there is drydep flux.
  !  </IN>
  !  <IN NAME="flux_bottom" TYPE="logical">
  !   .true. if there is bottom flux.
  !  </IN>
  !  <IN NAME="btm_reservoir" TYPE="logical">
  !   .true. if there is bottom reservoir.
  !  </IN>
  !  <IN NAME="move_vertical" TYPE="logical">
  !   .true. if there is active vertical movement
  !  </IN>
  !  <IN NAME="diff_vertical" TYPE="logical">
  !   .true. if there is random active vertical movement
  !  </IN>
  !  <IN NAME="flux_gas_molwt" TYPE="real">
  !   Molecular wt of gas defined in constants.F90 (g/mol)
  !  </IN>
  !  <IN NAME="flux_gas_param" TYPE="real, dimension(:)">
  !   Aray of parameters for gas flux (refer to documentation for subroutine aof_set_coupler_flux() ).
  !  </IN>
  !  <IN NAME="flux_param" TYPE="real, dimension(:)">
  !   Aray of parameters for non-gas flux (refer to documentation for subroutine aof_set_coupler_flux() ).
  !  </IN>
  !  <IN NAME="sink_rate" TYPE="real">
  !   Sinking rate if non-zero.
  !  </IN>
  !  <IN NAME="ice_restart_file" TYPE="character(len=*)">
  !   refer to documentation for subroutine aof_set_coupler_flux().
  !  </IN>
  !  <IN NAME="ocean_restart_file" TYPE="character(len=*)">
  !   refer to documentation for subroutine aof_set_coupler_flux().
  !  </IN>
  !  
  ! </SUBROUTINE>

  subroutine g_tracer_add(node_ptr, package, name, longname, units,  prog, const_init_value,init_value,&
       flux_gas, flux_gas_name, flux_runoff, flux_wetdep, flux_drydep, flux_gas_molwt, flux_gas_param, &
       flux_param, flux_bottom, btm_reservoir, move_vertical, diff_vertical, sink_rate, flux_gas_restart_file, flux_gas_type) 

    type(g_tracer_type), pointer :: node_ptr 
    character(len=*),   intent(in) :: package,name,longname,units
    logical,            intent(in) :: prog
    real,               intent(in), optional :: const_init_value
    real,               intent(in), optional :: init_value
    real,               intent(in), optional :: sink_rate
    logical,            intent(in), optional :: flux_gas
    logical,            intent(in), optional :: flux_runoff
    logical,            intent(in), optional :: flux_wetdep
    logical,            intent(in), optional :: flux_drydep
    logical,            intent(in), optional :: flux_bottom
    logical,            intent(in), optional :: btm_reservoir
    logical,            intent(in), optional :: move_vertical
    logical,            intent(in), optional :: diff_vertical
    real,               intent(in), optional :: flux_gas_molwt
    real, dimension(:), intent(in), optional :: flux_gas_param
    real, dimension(:), intent(in), optional :: flux_param
    character(len=*),   intent(in), optional :: flux_gas_name
    character(len=*),   intent(in), optional :: flux_gas_type
    character(len=*),   intent(in), optional :: flux_gas_restart_file

    !
    !       Local parameters
    !

    character(len=fm_string_len), parameter  :: sub_name = 'g_tracer_add'
    character(len=fm_string_len) :: flux_name
    !
    !       Local variables
    !
    type(g_tracer_type), pointer :: g_tracer => NULL()
    integer, save :: index = 0

    !===================================================================
    !Initialize the node
    !===================================================================
    allocate(g_tracer)

    !Specific properties
    index = index + 1
    g_tracer%index        = index
    g_tracer%name         = trim(name)
    g_tracer%longname     = trim(longname)
    g_tracer%package_name = trim(package)
    g_tracer%units        = trim(units)
    g_tracer%prog         = prog 

    !Restart files for tracers 
    g_tracer%ocean_restart_file = trim(g_tracer_com%ocean_restart_file)
    !Restart files for ice fluxes
    g_tracer%ice_restart_file   = trim(g_tracer_com%ice_restart_file)

    !Restart files for csurf, alpha and sc_no for tracers with gas flux default values
    g_tracer%flux_gas_restart_file = trim("ocean_airsea_flux.res.nc")    

    g_tracer%alias        = trim(name)
    !%alias is the global name for this tracer 
    !i.e., the name this tracer is known outside generic modules e.g., to model componenets.
    !Normally, %alias = %name as set above.
    !Sometimes, for debugging purposes the user may want to change the global name of a tracer
    !e.g., to clone a tracer and rrun it under a different name. 
    !This can be done for all tracers by:
    !g_tracer%alias        = trim("g_") // trim(name)


    !===================================================================
    !Allocate and initialize member field arrays
    !===================================================================
    !Note that const_init_value unlike init_value has special meaning in MOM 
    ! and if present the field is not restarted from a file!!
    if(present(const_init_value)) then
       g_tracer%const_init_value = const_init_value
       g_tracer%initial_value = const_init_value
    endif

    if(present(init_value))  g_tracer%initial_value = init_value

    !
    !Determine the fluxes 
    !

    if(present(flux_gas_molwt)) then
       g_tracer%flux_gas_molwt=flux_gas_molwt
    else
       g_tracer%flux_gas_molwt = 0.0
    endif
    if(present(flux_gas_param)) then
       allocate(g_tracer%flux_gas_param(size(flux_gas_param)))
       g_tracer%flux_gas_param=flux_gas_param
    endif
    if(present(flux_param)) then
       allocate(g_tracer%flux_param(size(flux_param)))
       g_tracer%flux_param=flux_param
    endif

    if(present(flux_gas))  g_tracer%flux_gas = flux_gas
    if(g_tracer%flux_gas) then
       g_tracer%flux_gas_name=trim(g_tracer%alias) // trim("_flux")
       if(present(flux_gas_name)) g_tracer%flux_gas_name=flux_gas_name
       if(present(flux_gas_restart_file))  g_tracer%flux_gas_restart_file  = flux_gas_restart_file
       g_tracer%flux_gas_type=trim("air_sea_gas_flux_generic")
       if(present(flux_gas_type)) g_tracer%flux_gas_type=flux_gas_type
    endif

    if(present(flux_runoff))  g_tracer%flux_runoff = flux_runoff
    if(g_tracer%flux_runoff) then
       g_tracer%flux_runoff_name=trim("runoff_") // trim(g_tracer%alias) 
    endif

    if(present(flux_wetdep))  g_tracer%flux_wetdep = flux_wetdep
    if(g_tracer%flux_wetdep) then
       g_tracer%flux_wetdep_name=trim("wet_dep_") // trim(g_tracer%alias)
    endif

    if(present(flux_drydep))  g_tracer%flux_drydep = flux_drydep
    if(g_tracer%flux_drydep) then
       g_tracer%flux_drydep_name=trim("dry_dep_") // trim(g_tracer%alias)
    endif

    if(present(flux_bottom))  g_tracer%flux_bottom = flux_bottom

    if(present(btm_reservoir)) g_tracer%has_btm_reservoir = btm_reservoir

    if(present(move_vertical)) g_tracer%move_vertical = move_vertical

    if(present(diff_vertical)) g_tracer%diff_vertical = diff_vertical

    if(present(sink_rate)) g_tracer%sink_rate = sink_rate

    !===================================================================
    !Reversed Linked List implementation! Make this new node to be the head of the list.
    !===================================================================    

    g_tracer%next => node_ptr 
    node_ptr => g_tracer 


  end subroutine g_tracer_add

  !
  !     Local functiion to remap the bounds of an array
  !     (Thanks to wikipedia for the suggestion)
  !

  function remap_bounds(ilb, jlb, klb, array) result(ptr)

  real, dimension(:,:,:),          pointer              :: ptr

  integer,                                 intent(in)   :: ilb
  integer,                                 intent(in)   :: jlb
  integer,                                 intent(in)   :: klb
  real, dimension(ilb:,jlb:,klb:), target, intent(in)   :: array

  ptr => array

  return
  end function remap_bounds

  subroutine g_tracer_init(g_tracer)
    type(g_tracer_type), pointer :: g_tracer
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed, nk,ntau,axes(3)

    !Get the common values for all tracers
    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes) 

    allocate(g_tracer%field(isd:ied,jsd:jed,nk,ntau));  g_tracer%field(:,:,:,:) = g_tracer%initial_value
    g_tracer%field_3d => remap_bounds(isd, jsd, 1, g_tracer%field(:,:,:,1))

    if(g_tracer%prog) then
       allocate(g_tracer%tendency(isd:ied,jsd:jed,nk)); g_tracer%tendency(:,:,:) = 0.0
       allocate(g_tracer%vdiffuse_impl(isd:ied,jsd:jed,nk))
       g_tracer%vdiffuse_impl(:,:,:) = 0.0
    endif

    if(g_tracer%flux_gas) then
       allocate(g_tracer%alpha(isd:ied,jsd:jed));g_tracer%alpha=0.0
       allocate(g_tracer%csurf(isd:ied,jsd:jed));g_tracer%csurf=0.0
       allocate(g_tracer%stf_gas(isd:ied,jsd:jed)); g_tracer%stf_gas(:,:) = 0.0 
       if(g_tracer%flux_gas_type .eq. 'air_sea_gas_flux_generic') then
          allocate(g_tracer%sc_no(isd:ied,jsd:jed));g_tracer%sc_no=0.0
          allocate(g_tracer%deltap(isd:ied,jsd:jed)); g_tracer%deltap(:,:) = 0.0 
          allocate(g_tracer%kw(isd:ied,jsd:jed)); g_tracer%kw(:,:) = 0.0 
       endif
    endif
    if(g_tracer%flux_runoff) then
       allocate(g_tracer%trunoff(isd:ied,jsd:jed));g_tracer%trunoff(:,:) = 0.0 
       allocate(g_tracer%runoff_tracer_flux(isd:ied,jsd:jed));g_tracer%runoff_tracer_flux(:,:) = 0.0 
    endif

    if(g_tracer%flux_wetdep) then
       allocate(g_tracer%wetdep(isd:ied,jsd:jed));g_tracer%wetdep(:,:) = 0.0 
    endif

    if(g_tracer%flux_drydep) then
       allocate(g_tracer%drydep(isd:ied,jsd:jed));g_tracer%drydep(:,:) = 0.0 
    endif

    if(g_tracer%flux_bottom) then
       allocate(g_tracer%btf(isd:ied,jsd:jed));g_tracer%btf(:,:) = 0.0 
    endif

    if(g_tracer%has_btm_reservoir) then
       allocate(g_tracer%btm_reservoir(isd:ied,jsd:jed));g_tracer%btm_reservoir(:,:) = 0.0 
    endif

    if(g_tracer%move_vertical) then
       allocate(g_tracer%vmove(isd:ied,jsd:jed, nk));g_tracer%vmove(:,:,:) = 0. 
    endif

    if(g_tracer%diff_vertical) then
       allocate(g_tracer%vdiff(isd:ied,jsd:jed, nk));g_tracer%vdiff(:,:,:) = 0. 
    endif
    !Surface flux %stf exists if one of the following fluxes were requested:

    if(g_tracer%flux_gas .or. g_tracer%flux_runoff .or. g_tracer%flux_wetdep .or. g_tracer%flux_drydep) then
       allocate(g_tracer%stf(isd:ied,jsd:jed)); g_tracer%stf(:,:) = 0.0 
    endif
    
  end subroutine g_tracer_init

  subroutine g_tracer_flux_init(g_tracer)
    type(g_tracer_type), pointer :: g_tracer


    !===================================================================
    !Get coupler flux indices
    !===================================================================
    !
    !For each kind of flux allocate the appropriate arrays only if that flux 
    !indicated to exist for the tracer.

    if(g_tracer%flux_gas) then
       g_tracer%flux_gas_ind  = aof_set_coupler_flux(g_tracer%flux_gas_name,                  &
            flux_type         = g_tracer%flux_gas_type,                                       &
            implementation    = 'ocmip2',                                                     &
            mol_wt            = g_tracer%flux_gas_molwt,                                      &
            param             = g_tracer%flux_gas_param,                                      &
            ice_restart_file  = g_tracer%ice_restart_file,                                    &
            ocean_restart_file= g_tracer%flux_gas_restart_file                                &
            )
    endif

    if(g_tracer%flux_runoff) then
       g_tracer%flux_runoff_ind  = aof_set_coupler_flux(g_tracer%flux_runoff_name,            &
            flux_type            = 'land_sea_runoff',                                         &
            implementation       = 'river',                                                   &
            param                = g_tracer%flux_param,                                       &
            ice_restart_file     = g_tracer%ice_restart_file                                  &
            )
    endif

    if(g_tracer%flux_wetdep) then
       g_tracer%flux_wetdep_ind  = aof_set_coupler_flux(g_tracer%flux_wetdep_name,            &
            flux_type            = 'air_sea_deposition',                                      &
            implementation       = 'wet',                                                     &
            param                = g_tracer%flux_param,                                       &
            ice_restart_file       = g_tracer%ice_restart_file                                &
            )
    endif

    if(g_tracer%flux_drydep) then
       g_tracer%flux_drydep_ind  = aof_set_coupler_flux(g_tracer%flux_drydep_name,            &
            flux_type            = 'air_sea_deposition',                                      &
            implementation       = 'dry',                                                     &
            param                = g_tracer%flux_param,                                       &
            ice_restart_file     = g_tracer%ice_restart_file                                  &
            )
    endif
    
    

  end subroutine g_tracer_flux_init


  ! <SUBROUTINE NAME="g_tracer_register_diag">
  !  <OVERVIEW>
  !   Diag-register all the internal fields that were _ALLOCATED for a tracer.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Use diag_manager register_diag_field for each of the field arrays that were _ALLOCATED for a tracer node.
  !   These include %field,  %tendency, %stf, %stf_gas, %deltap, %kw, %btf, %trunoff, %alpha, %csurf, %sc_no, %btm_reservoir. 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call  g_tracer_register_diag(g_tracer)
  !  </TEMPLATE>
  !  <IN NAME="g_tracer" TYPE="type(g_tracer_type), pointer">
  !   Pointer to this tracer node. 
  !  </IN>
  ! </SUBROUTINE>

  subroutine g_tracer_register_diag(g_tracer)
    type(g_tracer_type), pointer :: g_tracer

    character(len=fm_string_len) :: string

    g_tracer%diag_id_field = register_diag_field(g_tracer%package_name, &
         trim(g_tracer%alias),         &
         g_tracer_com%axes(1:3),       &
         g_tracer_com%init_time,       &
         trim(g_tracer%longname),      &
         trim(g_tracer%units),         &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_taup1")
    g_tracer%diag_id_field_taup1 = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:3),       &
         g_tracer_com%init_time,       &
         trim(g_tracer%longname) // ' at taup1',      &
         trim(g_tracer%units),         &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_aux")
    g_tracer%diag_id_aux = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:3),       &
         g_tracer_com%init_time,       &
         trim(string),                 &
         trim(g_tracer%units),         &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_vmove")
    g_tracer%diag_id_vmove = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:3),       &
         g_tracer_com%init_time,       &
         trim('vertical movement'),    &
         trim('m/s'),                  &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_vdiffuse_impl")
    g_tracer%diag_id_vdiffuse_impl = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:3),       &
         g_tracer_com%init_time,       &
         'Implicit vertical diffusion of ' // trim(g_tracer%alias),      &
         trim('mole/m^2/s'),                  &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_tendency")
    g_tracer%diag_id_tendency = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:3),       &
         g_tracer_com%init_time,       &
         'Generic tracer tendency of ' // trim(g_tracer%alias),      &
         trim('mole/m^2/s'),                  &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_vdiff")
    g_tracer%diag_id_vdiff = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:3),       &
         g_tracer_com%init_time,       &
         trim('random movement'),      &
         trim('m/s'),                  &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_stf")
    g_tracer%diag_id_stf = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:2),       &
         g_tracer_com%init_time,       &
         trim('Total flux of ') // trim(g_tracer%alias) // trim(' into Ocean Surface'), &
         trim('mole/m^2/sec'),         &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_stf_gas")
    g_tracer%diag_id_stf_gas = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:2),       &
         g_tracer_com%init_time,       &
         trim('Gas exchange flux of ') // trim(g_tracer%alias) // trim(' into Ocean Surface'), &
         trim('mole/m^2/sec'),         &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_deltap")
    g_tracer%diag_id_deltap = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:2),       &
         g_tracer_com%init_time,       &
         trim('Ocn minus Atm pressure of ') // trim(g_tracer%alias), &
         trim('uatm'),                 &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_kw")
    g_tracer%diag_id_kw = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:2),       &
         g_tracer_com%init_time,       &
         trim('Gas Exchange piston velocity for ') // trim(g_tracer%alias), &
         trim('m/sec'),                &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_btf")
    g_tracer%diag_id_btf = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:2),       &
         g_tracer_com%init_time,       &
         trim('Total flux of ') // trim(g_tracer%alias) // trim(' into Ocean Bottom'), &
         trim('mole/m^2/sec'),         &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_btm_reservoir")
    g_tracer%diag_id_btm = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:2),       &
         g_tracer_com%init_time,       &
         trim('Bottom reservoir of ') // trim(g_tracer%alias), &
         trim(g_tracer%units),         &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_trunoff")
    g_tracer%diag_id_trunoff = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:2),       &
         g_tracer_com%init_time,       &
         trim('River concentration of ') // trim(g_tracer%alias), &
         trim(g_tracer%units),         &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_alpha")
    g_tracer%diag_id_alpha = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:2),       &
         g_tracer_com%init_time,       &
         trim('Atmospheric saturation for ') // trim(g_tracer%alias), &
         trim(g_tracer%units),         &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_csurf")
    g_tracer%diag_id_csurf = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:2),       &
         g_tracer_com%init_time,       &
         trim('Ocean surface gas concentration of ') // trim(g_tracer%alias), &
         trim(g_tracer%units),         &
         missing_value = -1.0e+20)

    string=trim(g_tracer%alias) // trim("_sc_no")
    g_tracer%diag_id_sc_no = register_diag_field(g_tracer%package_name, &
         trim(string),                 &
         g_tracer_com%axes(1:2),       &
         g_tracer_com%init_time,       &
         trim('Ocean surface Schmidt Number for ') // trim(g_tracer%alias), &
         trim(g_tracer%units),         &
         missing_value = -1.0e+20)

  end subroutine g_tracer_register_diag

  ! <SUBROUTINE NAME="g_tracer_coupler_set">
  !  <OVERVIEW>
  !   Set coupler values only for tracers that have _ALLOCATED %alpha, %csurf and %sc_no
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Use coupler_util subroutine set_coupler_values() to set the coupler values
  !   for fluxes to be exchanged with Ice for the requested fluxes.
  !   NOTE:
  !   This is a collective subroutine and will traverese the list of generic tracers and 
  !   set the coupler values for each tracer node accordingly.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call g_tracer_coupler_set(g_tracer_list,IOB_struc,value)
  !  </TEMPLATE>
  !  <IN NAME="g_tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of the generic tracer list.
  !  </IN>
  !  <IN NAME="IOB_struc" TYPE="type(coupler_2d_bc_type)">
  !   The coupler flux IOB structure. 
  !  </IN>
  !  OPTIONAL ARGS:
  !  <IN NAME="value" TYPE="real">
  !   Set the coupler values to a constant (particularly 0) is desired.
  !  </IN>
  ! </SUBROUTINE>

  subroutine g_tracer_coupler_set(g_tracer_list,IOB_struc,value)
    type(g_tracer_type), pointer :: g_tracer_list,g_tracer 
    type(coupler_2d_bc_type), intent(inout) :: IOB_struc
    real, optional :: value

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_coupler_set'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list
    !Go through the list of tracers 
    do  
       !
       !Set coupler values only for tracers that have _ALLOCATED %alpha, %csurf and %sc_no
       !

       if(_ALLOCATED(g_tracer%alpha)) then

          if(present(value)) g_tracer%alpha=value 

          call set_coupler_values(g_tracer%alpha,   &
               BC_struc   = IOB_struc,              &
               BC_index   = g_tracer%flux_gas_ind,  &
               BC_element = ind_alpha,              &
               ilb=g_tracer_com%isd, jlb=g_tracer_com%jsd ,&
               is=g_tracer_com%isc, ie=g_tracer_com%iec,&
               js=g_tracer_com%jsc, je=g_tracer_com%jec &
               )
       endif

       if(_ALLOCATED(g_tracer%csurf)) then

          if(present(value)) g_tracer%csurf=value 

          call set_coupler_values(g_tracer%csurf,   &
               BC_struc   = IOB_struc,              &
               BC_index   = g_tracer%flux_gas_ind,  &
               BC_element = ind_csurf,              &
               ilb=g_tracer_com%isd, jlb=g_tracer_com%jsd ,&
               is=g_tracer_com%isc, ie=g_tracer_com%iec,&
               js=g_tracer_com%jsc, je=g_tracer_com%jec &
               )
       endif

       if(_ALLOCATED(g_tracer%sc_no)) then

          if(present(value)) g_tracer%sc_no=value 

          call set_coupler_values(g_tracer%sc_no,   &
               BC_struc   = IOB_struc,              &
               BC_index   = g_tracer%flux_gas_ind,  &
               BC_element = ind_sc_no,              &
               ilb=g_tracer_com%isd, jlb=g_tracer_com%jsd ,&
               is=g_tracer_com%isc, ie=g_tracer_com%iec,&
               js=g_tracer_com%jsc, je=g_tracer_com%jec &
               )
       endif

       !traverse the linked list till hit NULL
       if(.NOT. associated(g_tracer%next)) exit
       g_tracer => g_tracer%next
    enddo

  end subroutine g_tracer_coupler_set

  ! <SUBROUTINE NAME="g_tracer_coupler_get">
  !  <OVERVIEW>
  !   Get coupler values only for tracers that have _ALLOCATED arrays for the fluxes
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Use coupler_util subroutine extract_coupler_values() to get the coupler values
  !   for fluxes to be exchanged with Ice for the requested fluxes only.
  !   NOTE:
  !   This is a collective subroutine and will traverese the list of generic tracers and 
  !   get the coupler values for each tracer node accordingly.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call g_tracer_coupler_get(g_tracer_list,IOB_struc)
  !  </TEMPLATE>
  !  <IN NAME="g_tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of the generic tracer list.
  !  </IN>
  !  <IN NAME="IOB_struc" TYPE="type(coupler_2d_bc_type)">
  !   The coupler flux IOB structure. 
  !  </IN>
  ! </SUBROUTINE>

  subroutine g_tracer_coupler_get(g_tracer_list,IOB_struc)
    type(g_tracer_type),          pointer :: g_tracer_list, g_tracer 
    type(coupler_2d_bc_type),    intent(in) :: IOB_struc

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_coupler_get'
    real, dimension(:,:), allocatable :: temp_array,stf_array,stf_gas_array,deltap_array, kw_array

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list
    allocate(temp_array(g_tracer_com%isd:g_tracer_com%ied,g_tracer_com%jsd:g_tracer_com%jed));temp_array=0.0
    allocate(stf_array(g_tracer_com%isd:g_tracer_com%ied,g_tracer_com%jsd:g_tracer_com%jed))
    allocate(stf_gas_array(g_tracer_com%isd:g_tracer_com%ied,g_tracer_com%jsd:g_tracer_com%jed))
    allocate(deltap_array(g_tracer_com%isd:g_tracer_com%ied,g_tracer_com%jsd:g_tracer_com%jed))
    allocate(kw_array(g_tracer_com%isd:g_tracer_com%ied,g_tracer_com%jsd:g_tracer_com%jed))


    !Go through the list of tracers 
    do  
       stf_array=0.0

       if(g_tracer%flux_gas) then
          temp_array=0.0
          call extract_coupler_values(BC_struc  =IOB_struc, &
               BC_index  =g_tracer%flux_gas_ind,    & 
               BC_element=ind_flux,                 &
               array_out =temp_array,               &
               conversion=-1.0,                     &
               ilb=g_tracer_com%isd,jlb=g_tracer_com%jsd,& !lower bounds of array_out 
               is=g_tracer_com%isc, ie=g_tracer_com%iec,&
               js=g_tracer_com%jsc, je=g_tracer_com%jec)
          !This does temp_array=conv *BC_struc%bc(flux_gas_ind)%field(ind_flux)%values

          stf_array = stf_array+temp_array !flux_gas contributes to %stf
          stf_gas_array=stf_array


          if(g_tracer%flux_gas_type .eq. 'air_sea_gas_flux_generic') then
             temp_array=0.0
             deltap_array=0.0
             call extract_coupler_values(BC_struc  =IOB_struc, &
                  BC_index  =g_tracer%flux_gas_ind,    & 
                  BC_element=ind_deltap,               &
                  array_out =temp_array,               &
                  conversion=1.0,                      &
                  ilb=g_tracer_com%isd,jlb=g_tracer_com%jsd,& !lower bounds of array_out 
                  is=g_tracer_com%isc, ie=g_tracer_com%iec,&
                  js=g_tracer_com%jsc, je=g_tracer_com%jec)
             !This does temp_array=conv *BC_struc%bc(flux_gas_ind)%field(ind_flux)%values
             
             deltap_array = deltap_array+temp_array
             temp_array=0.0
             kw_array=0.0
             call extract_coupler_values(BC_struc  =IOB_struc, &
                  BC_index  =g_tracer%flux_gas_ind,    & 
                  BC_element=ind_kw,                   &
                  array_out =temp_array,               &
                  conversion=1.0,                      &
                  ilb=g_tracer_com%isd,jlb=g_tracer_com%jsd,& !lower bounds of array_out
                  is=g_tracer_com%isc, ie=g_tracer_com%iec,&
                  js=g_tracer_com%jsc, je=g_tracer_com%jec)
             !This does temp_array=conv *BC_struc%bc(flux_gas_ind)%field(ind_flux)%values
             
             kw_array = kw_array+temp_array
          endif
          
       endif

       if(g_tracer%flux_drydep) then
          temp_array=0.0
          call extract_coupler_values(BC_struc  =IOB_struc, &
               BC_index  =g_tracer%flux_drydep_ind,  &
               BC_element=ind_flux,                  &
               array_out =temp_array,                &
               conversion=-1.0,                      &
               ilb=g_tracer_com%isd,jlb=g_tracer_com%jsd,&
               is=g_tracer_com%isc, ie=g_tracer_com%iec, &
               js=g_tracer_com%jsc, je=g_tracer_com%jec)

          stf_array = stf_array+temp_array !flux_drydep contributes to %stf

          call g_tracer_set_values(g_tracer,g_tracer%name,'drydep',temp_array,&
               g_tracer_com%isd,g_tracer_com%jsd)
       endif

       if(g_tracer%flux_wetdep) then
          temp_array=0.0
          call extract_coupler_values(BC_struc  =IOB_struc, &
               BC_index  =g_tracer%flux_wetdep_ind,  &
               BC_element=ind_flux,                  &
               array_out =temp_array,                &
               conversion=-1.0,                      &
               ilb=g_tracer_com%isd,jlb=g_tracer_com%jsd,&
               is=g_tracer_com%isc, ie=g_tracer_com%iec, &
               js=g_tracer_com%jsc, je=g_tracer_com%jec)

          stf_array = stf_array+temp_array  !flux_wetdep contributes to %stf

          call g_tracer_set_values(g_tracer,g_tracer%name,'wetdep',temp_array,&
               g_tracer_com%isd,g_tracer_com%jsd)
       endif

       if(g_tracer%flux_runoff) then
          temp_array=0.0
          call extract_coupler_values(BC_struc  =IOB_struc, &
               BC_index  =g_tracer%flux_runoff_ind, &
               BC_element=ind_flux,                 &
               array_out =temp_array,               &
               conversion=1.0,                      &
               ilb=g_tracer_com%isd,jlb=g_tracer_com%jsd,&
               is=g_tracer_com%isc, ie=g_tracer_com%iec, &
               js=g_tracer_com%jsc, je=g_tracer_com%jec)

          call g_tracer_set_values(g_tracer,g_tracer%name,'trunoff',temp_array,&
               g_tracer_com%isd,g_tracer_com%jsd)
       endif

       !Any of the following fluxes contribute to %stf
       !gas, wetdep and drydep contribute explicitly here.
       !runoff contributes to %stf in GOLD but not in MOM, 
       !so it will be added later in the model-dependent driver code (GOLD_generic_tracer.F90)

       if(g_tracer%flux_gas .or. g_tracer%flux_drydep .or. g_tracer%flux_wetdep .or. g_tracer%flux_runoff ) then
          call g_tracer_set_values(g_tracer,g_tracer%name,'stf',stf_array,&
               g_tracer_com%isd,g_tracer_com%jsd)
       endif

       if(g_tracer%flux_gas) then
          call g_tracer_set_values(g_tracer,g_tracer%name,'stf_gas',stf_gas_array,&
               g_tracer_com%isd,g_tracer_com%jsd)
          if(g_tracer%flux_gas_type .eq. 'air_sea_gas_flux_generic') then
             call g_tracer_set_values(g_tracer,g_tracer%name,'deltap',deltap_array,&
                  g_tracer_com%isd,g_tracer_com%jsd)
             call g_tracer_set_values(g_tracer,g_tracer%name,'kw',kw_array,&
                  g_tracer_com%isd,g_tracer_com%jsd)
          endif
       endif

       !traverse the linked list till hit NULL
       if(.NOT. associated(g_tracer%next)) exit
       g_tracer => g_tracer%next
    enddo

    deallocate(temp_array, stf_array, stf_gas_array, deltap_array, kw_array)

  end subroutine g_tracer_coupler_get

  ! <SUBROUTINE NAME="g_tracer_set_common">
  !  <OVERVIEW>
  !   Set common values and arrays for ALL generic tracers to share
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   ALL generic tracers share the same properties such as 2D Domain, # of depth levels, # of time steps retained
  !   grid_mask array and initial time.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call g_tracer_set_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes,grid_tmask,grid_kmt,init_time)
  !  </TEMPLATE>
  !  <IN NAME="" TYPE="">
  !   
  !  </IN>
  ! </SUBROUTINE>

  subroutine g_tracer_set_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes,grid_tmask,grid_kmt,init_time)
    integer,                     intent(in) :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes(3)
    real, dimension(isd:,jsd:,:),intent(in) :: grid_tmask
    integer,dimension(isd:,jsd:),intent(in) :: grid_kmt
    type(time_type),             intent(in) :: init_time 

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_set_common'
    integer :: i,j

    !Here we assume that all the tracers in the list have the same following properties

    g_tracer_com%isd=isd
    g_tracer_com%ied=ied
    g_tracer_com%jsd=jsd
    g_tracer_com%jed=jed
    g_tracer_com%isc=isc
    g_tracer_com%iec=iec
    g_tracer_com%jsc=jsc
    g_tracer_com%jec=jec
    g_tracer_com%nk =nk
    g_tracer_com%ntau=ntau
    g_tracer_com%axes=axes
    g_tracer_com%init_time=init_time

    if(.NOT. _ALLOCATED(g_tracer_com%grid_tmask)) allocate(g_tracer_com%grid_tmask(isd:ied,jsd:jed,nk))
    g_tracer_com%grid_tmask=grid_tmask 


    if(.NOT. _ALLOCATED(g_tracer_com%grid_kmt)) allocate(g_tracer_com%grid_kmt(isd:ied,jsd:jed))    
    g_tracer_com%grid_kmt = grid_kmt

    if(.NOT. _ALLOCATED(g_tracer_com%grid_mask_coast)) allocate(g_tracer_com%grid_mask_coast(isd:ied,jsd:jed))

    !Determine the coast line.
    !In order to that grid_tmask must have the proper value on the data domain boundaries isd,ied,jsd,jed
    !so that we can decide if the coast line coinsides with a point on the compute domain boundary

    g_tracer_com%grid_mask_coast(:,:) = 0
    do j =jsc, jec ; do i = isc, iec  
       if (g_tracer_com%grid_tmask(i,j,1) .gt. 0) then 
          if (g_tracer_com%grid_tmask(i-1,j,1) .lt. 1 .or. g_tracer_com%grid_tmask(i,j-1,1) .lt. 1 .or. &
              g_tracer_com%grid_tmask(i+1,j,1) .lt. 1 .or. g_tracer_com%grid_tmask(i,j+1,1) .lt. 1) then !{
             g_tracer_com%grid_mask_coast(i,j) = 1
          endif
       endif
    enddo; enddo     


  end subroutine g_tracer_set_common

  ! <SUBROUTINE NAME="g_tracer_get_common">
  !  <OVERVIEW>
  !   Get common values and arrays for ALL generic tracers to share
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,&
  !     axes,grid_tmask,grid_mask_coast,grid_kmt,init_time)
  !  </TEMPLATE>
  !  <IN NAME="" TYPE="">
  !   
  !  </IN>
  ! </SUBROUTINE>

  subroutine g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,&
       axes,grid_tmask,grid_mask_coast,grid_kmt,init_time)

    integer,               intent(out) :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau
    integer,optional,      intent(out) :: axes(3)
    type(time_type), optional,      intent(out) :: init_time 
    real, optional, dimension(:,:,:),pointer    :: grid_tmask
    integer, optional, dimension(:,:),  pointer :: grid_mask_coast
    integer, optional, dimension(:,:),  pointer :: grid_kmt

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_get_common'

    !Here we assume that all the tracers in the list have the same following properties

    isd=g_tracer_com%isd
    ied=g_tracer_com%ied
    jsd=g_tracer_com%jsd
    jed=g_tracer_com%jed
    isc=g_tracer_com%isc
    iec=g_tracer_com%iec
    jsc=g_tracer_com%jsc
    jec=g_tracer_com%jec
    nk =g_tracer_com%nk
    ntau=g_tracer_com%ntau
    if(present(axes))             axes = g_tracer_com%axes 
    if(present(init_time))        init_time=g_tracer_com%init_time
    if(present(grid_tmask))       grid_tmask => g_tracer_com%grid_tmask
    if(present(grid_mask_coast))  grid_mask_coast=> g_tracer_com%grid_mask_coast
    if(present(grid_kmt))         grid_kmt => g_tracer_com%grid_kmt
!    if(present(ice_restart_file)) ice_restart_file    = g_tracer_com%ice_restart_file
!    if(present(ocean_restart_file)) ocean_restart_file  = g_tracer_com%ocean_restart_file

  end subroutine g_tracer_get_common

  subroutine g_tracer_set_files(ice_restart_file,ocean_restart_file)
    character(len=*),   intent(in) :: ice_restart_file
    character(len=*),   intent(in) :: ocean_restart_file

    g_tracer_com%ice_restart_file    = ice_restart_file
    g_tracer_com%ocean_restart_file  = ocean_restart_file

  end subroutine g_tracer_set_files
    
  !Overload interface g_tracer_get_pointer for 4D fields

  subroutine g_tracer_get_4D(g_tracer_list,name,member,array_ptr)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    real, dimension(:,:,:,:), pointer    :: array_ptr

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_get_4D'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!
    !Find the node which has name=name
    call g_tracer_find(g_tracer,name)
    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list with name="//trim(name))

    select case(member)
    case ('field')
       if(associated(g_tracer%field)) then 
          array_ptr => g_tracer%field
       else
          call mpp_error(FATAL, trim(sub_name)//": Cannot get member variable: "//trim(name)//" % "//trim(member))
       endif
    case default 
       call mpp_error(FATAL, trim(sub_name)//": Not a known member variable: "//trim(name)//" % "//trim(member))   
    end select

  end subroutine g_tracer_get_4D

  !Overload interface g_tracer_get_pointer for 3D fields

  subroutine g_tracer_get_3D(g_tracer_list,name,member,array_ptr)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    real, dimension(:,:,:), pointer    :: array_ptr

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_get_3D'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!
    !Find the node which has name=name
    call g_tracer_find(g_tracer,name)
    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list with name="//trim(name))

    select case(member)
    case ('field') 
       if(associated(g_tracer%field3d_ptr)) then 
          array_ptr => g_tracer%field3d_ptr
       elseif(associated(g_tracer%field_3d)) then
          array_ptr => g_tracer%field_3d
       else
          call mpp_error(FATAL, trim(sub_name)//": Cannot get member variable: "//trim(name)//" % "//trim(member))
       endif
    case ('vmove') 
       array_ptr => g_tracer%vmove
    case ('vdiff') 
       array_ptr => g_tracer%vdiff
    case ('vdiffuse_impl') 
       array_ptr => g_tracer%vdiffuse_impl
    case default 
       call mpp_error(FATAL, trim(sub_name)//": Not a known member variable: "//trim(member))   
    end select

  end subroutine g_tracer_get_3D  

  !Overload interface g_tracer_get_pointer for 2D fields

  subroutine g_tracer_get_2D(g_tracer_list,name,member,array_ptr)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    real, dimension(:,:), pointer    :: array_ptr

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_get_4D'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!
    !Find the node which has name=name
    call g_tracer_find(g_tracer,name)
    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list with name="//trim(name))

    select case(member)
    case ('alpha') 
       array_ptr => g_tracer%alpha
    case ('csurf') 
       array_ptr => g_tracer%csurf
    case ('sc_no') 
       array_ptr => g_tracer%sc_no
    case ('stf') 
       array_ptr => g_tracer%stf
    case ('stf_gas') 
       array_ptr => g_tracer%stf_gas
    case ('deltap') 
       array_ptr => g_tracer%deltap
    case ('kw') 
       array_ptr => g_tracer%kw
    case ('btf') 
       array_ptr => g_tracer%btf
    case ('btm_reservoir') 
       array_ptr => g_tracer%btm_reservoir
    case ('trunoff') 
       array_ptr => g_tracer%trunoff
    case ('runoff_tracer_flux') 
       array_ptr => g_tracer%runoff_tracer_flux
    case ('drydep') 
       array_ptr => g_tracer%drydep
    case ('wetdep') 
       array_ptr => g_tracer%wetdep
    case default 
       call mpp_error(FATAL, trim(sub_name)//": Not a known member variable: "//trim(member))   
    end select

  end subroutine g_tracer_get_2D

  !Overload interface g_tracer_get_values for 4D fields

  subroutine g_tracer_get_4D_val(g_tracer_list,name,member,array,isd,jsd)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer 
    integer,                  intent(in) :: isd,jsd
    real, dimension(isd:,jsd:,:,:), intent(out):: array

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_get_4D_val'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!
    !Find the node which has name=name
    call g_tracer_find(g_tracer,name)
    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list with name="//trim(name))

    select case(member)
    case ('field') 
       if(associated(g_tracer%field)) then 
          array = g_tracer%field
       else
          call mpp_error(FATAL, trim(sub_name)//": Cannot get member variable: "//trim(name)//" % "//trim(member))
       endif
    case default 
       call mpp_error(FATAL, trim(sub_name)//": Not a known member variable: "//trim(name)//" % "//trim(member))   
    end select

  end subroutine g_tracer_get_4D_val

  !Overload interface g_tracer_get_values for 3D fields

  subroutine g_tracer_get_3D_val(g_tracer_list,name,member,array,isd,jsd,ntau,positive)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer 
    integer,                  intent(in) :: isd,jsd
    integer, optional,        intent(in) :: ntau
    logical, optional,        intent(in) :: positive
    real, dimension(isd:,jsd:,:), intent(out):: array
    integer :: tau
    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_get_3D_val'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!
    !Find the node which has name=name
    call g_tracer_find(g_tracer,name)
    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list with name="//trim(name))

    tau = 1
    if(present(ntau)) tau = ntau

    select case(member)
    case ('field') 
       if(associated(g_tracer%field)) then 
          array(:,:,:) = g_tracer%field(:,:,:,tau)
       elseif(associated(g_tracer%field3d_ptr)) then 
          array(:,:,:) = g_tracer%field3d_ptr(:,:,:)
       else
          call mpp_error(FATAL, trim(sub_name)//": Cannot get member variable: "//trim(name)//" % "//trim(member))
       endif
          
       if(present(positive)) array = max(0.0,array)
    case ('tendency') 
       array(:,:,:) = g_tracer%tendency(:,:,:)
    case ('vmove') 
       array(:,:,:) = g_tracer%vmove(:,:,:)
    case ('vdiff') 
       array(:,:,:) = g_tracer%vdiff(:,:,:)
    case ('vdiffuse_impl') 
       array(:,:,:) = g_tracer%vdiffuse_impl(:,:,:)
    case default 
       call mpp_error(FATAL, trim(sub_name)//": Not a known member variable: "//trim(member))   
    end select

  end subroutine g_tracer_get_3D_val

  !Overload interface g_tracer_get_values for 2D fields

  subroutine g_tracer_get_2D_val(g_tracer_list,name,member,array,isd,jsd)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer 
    integer,                  intent(in) :: isd,jsd
    real, dimension(isd:,jsd:), intent(out):: array

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_get_2D_val'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!
    !Find the node which has name=name
    call g_tracer_find(g_tracer,name)
    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list with name="//trim(name))

    select case(member)
    case ('alpha') 
       array = g_tracer%alpha
    case ('csurf') 
       array = g_tracer%csurf
    case ('sc_no') 
       array = g_tracer%sc_no
    case ('stf') 
       array = g_tracer%stf
    case ('stf_gas') 
       array = g_tracer%stf_gas
    case ('deltap') 
       array = g_tracer%deltap
    case ('kw') 
       array = g_tracer%kw
    case ('btf') 
       array = g_tracer%btf
    case ('btm_reservoir') 
       array = g_tracer%btm_reservoir
    case ('trunoff') 
       array = g_tracer%trunoff
    case ('runoff_tracer_flux') 
       array = g_tracer%runoff_tracer_flux
    case ('drydep') 
       array = g_tracer%drydep
    case ('wetdep') 
       array = g_tracer%wetdep
    case default 
       call mpp_error(FATAL, trim(sub_name)//": Not a known member variable: "//trim(member))   
    end select

  end subroutine g_tracer_get_2D_val

  !Overload interface g_tracer_get_values for 1D fields

  subroutine g_tracer_get_real(g_tracer_list,name,member,value)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer 
    real,                     intent(out):: value

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_get_real'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!
    !Find the node which has name=name
    call g_tracer_find(g_tracer,name)
    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list with name="//trim(name))

    select case(member)
    case ('sink_rate') 
       value = g_tracer%sink_rate
    case ('const_init_value') 
       value = g_tracer%const_init_value
    case default 
       call mpp_error(FATAL, trim(sub_name)//": Not a known member variable: "//trim(member))   
    end select

  end subroutine g_tracer_get_real

  !Overload interface g_tracer_get_values for string members

  subroutine g_tracer_get_string(g_tracer_list,name,member,string)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer 
    character(len=fm_string_len), intent(out) :: string
    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_get_string'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!

    !If queried for 'name' return the %name of the head 
    if(member .eq. 'name') then
       string=g_tracer%name
       return
    endif

    !Find the node which has name=name
    call g_tracer_find(g_tracer,name)
    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list with name="//trim(name))

    select case(member)
    case ('longname') 
       string = g_tracer%longname
    case ('alias') 
       string = g_tracer%alias
    case ('units') 
       string = g_tracer%units
    case ('package') 
       string = g_tracer%package_name
    case ('ocean_restart_file') 
       string = g_tracer%ocean_restart_file
    case default 
       call mpp_error(FATAL, trim(sub_name)//": Not a known member variable: "//trim(member))   
    end select

  end subroutine g_tracer_get_string

  !Overload interface g_tracer_set_values for 2D fields

  subroutine g_tracer_set_2D(g_tracer_list,name,member,array,isd,jsd)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    integer,                  intent(in) :: isd,jsd
    real, dimension(isd:,jsd:),intent(in) :: array

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_set_2D'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!

    !Find the node which has name=name
    call g_tracer_find(g_tracer,name)
    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list with name="//trim(name))

    select case(member)
    case ('alpha') 
       g_tracer%alpha  = array 
    case ('csurf')
       g_tracer%csurf  = array
    case ('sc_no')
       g_tracer%sc_no  = array
    case ('stf') 
       g_tracer%stf    = array
    case ('stf_gas') 
       g_tracer%stf_gas= array
    case ('deltap') 
       g_tracer%deltap = array
    case ('kw') 
       g_tracer%kw     = array
    case ('btf') 
       g_tracer%btf    = array
    case ('btm_reservoir') 
       g_tracer%btm_reservoir = array
    case ('trunoff')
       g_tracer%trunoff = array
    case ('runoff_tracer_flux')
       g_tracer%runoff_tracer_flux = array
    case ('drydep')
       g_tracer%drydep = array
    case ('wetdep')
       g_tracer%wetdep = array
    case default 
       call mpp_error(FATAL, trim(sub_name)//": Not a known member variable: "//trim(member))   
    end select

  end subroutine g_tracer_set_2D

  !Overload interface g_tracer_set_values for 3D fields

  subroutine g_tracer_set_3D(g_tracer_list,name,member,array,isd,jsd,ntau)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    integer,                  intent(in) :: isd,jsd
    integer, optional,        intent(in) :: ntau
    real, dimension(isd:,jsd:,:), intent(in)       :: array
    integer :: tau

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_set_3D'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!

    !Find the node which has name=name
    call g_tracer_find(g_tracer,name)
    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list with name="//trim(name))

    tau = 1
    if(present(ntau)) tau = ntau

    select case(member)
    case ('tendency') 
       g_tracer%tendency  = array 
    case ('field') 
       if(associated(g_tracer%field)) then 
          g_tracer%field(:,:,:,tau) = array(:,:,:) 
       elseif(associated(g_tracer%field3d_ptr)) then 
          g_tracer%field3d_ptr(:,:,:) = array(:,:,:) 
       else
          call mpp_error(FATAL, trim(sub_name)//": Cannot set member variable: "//trim(name)//" % "//trim(member))
       endif      
    case ('vmove') 
       g_tracer%vmove  = array 
    case ('vdiff') 
       g_tracer%vdiff  = array 
    case ('vdiffuse_impl') 
       g_tracer%vdiffuse_impl  = array 
    case default 
       call mpp_error(FATAL, trim(sub_name)//": Not a known member variable: "//trim(member))   
    end select

  end subroutine g_tracer_set_3D

  !Overload interface g_tracer_set_values for 4D fields

  subroutine g_tracer_set_4D(g_tracer_list,name,member,array,isd,jsd)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    integer,                  intent(in) :: isd,jsd
    real, dimension(isd:,jsd:,:,:), intent(in)       :: array

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_set_4D'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!

    !Find the node which has name=name
    call g_tracer_find(g_tracer,name)
    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list with name="//trim(name))

    select case(member)
    case ('field')
       if(associated(g_tracer%field)) then
          g_tracer%field = array
       else
          call mpp_error(FATAL, trim(sub_name)//": Cannot set member variable: "//trim(name)//" % "//trim(member))
       endif
    case default 
       call mpp_error(FATAL, trim(sub_name)//": Not a known member variable: "//trim(member))   
    end select

  end subroutine g_tracer_set_4D

  !Overload interface g_tracer_set_values for 1D fields

  subroutine g_tracer_set_real(g_tracer_list,name,member,value)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    real,                     intent(in) :: value

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_set_real'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!

    !Find the node which has name=name
    call g_tracer_find(g_tracer,name)
    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list with name="//trim(name))

    select case(member)
    case ('field') 
       if(associated(g_tracer%field)) then
          g_tracer%field = value !Set all elements to value
       else
          call mpp_error(FATAL, trim(sub_name)//": Cannot set member variable: "//trim(name)//" % "//trim(member))
       endif
    case ('tendency') 
       g_tracer%tendency  = value 
    case ('alpha') 
       g_tracer%alpha     = value 
    case ('csurf') 
       g_tracer%csurf     = value 
    case ('sc_no') 
       g_tracer%sc_no     = value 
    case ('stf') 
       g_tracer%stf       = value 
    case ('stf_gas') 
       g_tracer%stf_gas   = value 
    case ('deltap') 
       g_tracer%deltap    = value 
    case ('kw') 
       g_tracer%kw        = value 
    case ('btf') 
       g_tracer%btf       = value 
    case ('trunoff') 
       g_tracer%trunoff   = value 
    case ('runoff_tracer_flux') 
       g_tracer%runoff_tracer_flux = value 
    case ('btm_reservoir') 
       g_tracer%btm_reservoir = value 
    case ('sink_rate') 
       g_tracer%sink_rate = value 
    case default 
       call mpp_error(FATAL, trim(sub_name)//": Not a known member variable: "//trim(member))   
    end select

  end subroutine g_tracer_set_real

  subroutine g_tracer_set_pointer_4D(g_tracer_list,name,member,array,ilb,jlb)
    character(len=*),               intent(in) :: name
    character(len=*),               intent(in) :: member
    type(g_tracer_type),            pointer    :: g_tracer_list, g_tracer
    integer,                        intent(in) :: ilb,jlb
    real, dimension(ilb:,jlb:,:,:), target, intent(in) :: array

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_set_pointer_4D'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!

    !Find the node which has name=name
    call g_tracer_find(g_tracer,name)
    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list with name="//trim(name))

    select case(member)
    case ('field') 
       if (associated(g_tracer%field )) then
          call mpp_error(NOTE, trim(sub_name) // ": Deallocating generic tracer "//trim(name)//" % "//trim(member))
          deallocate( g_tracer%field )
       endif
       g_tracer%field  => array 
    case default 
       call mpp_error(FATAL, trim(sub_name)//": Not a supported operation for member variable: "//trim(name)//" % "//trim(member))
    end select

  end subroutine g_tracer_set_pointer_4D

  subroutine g_tracer_set_pointer_3D(g_tracer_list,name,member,array,ilb,jlb)
    character(len=*),               intent(in) :: name
    character(len=*),               intent(in) :: member
    type(g_tracer_type),            pointer    :: g_tracer_list, g_tracer
    integer,                        intent(in) :: ilb,jlb
    real, dimension(ilb:,jlb:,:), target, intent(in) :: array
    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_set_pointer_3D'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!

    !Find the node which has name=name
    call g_tracer_find(g_tracer,name)
    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list with name="//trim(name))

    select case(member)
    case ('tendency') 
       if (associated( g_tracer%tendency )) then
          call mpp_error(NOTE, trim(sub_name) // ": Deallocating generic tracer "//trim(name)//" % "//trim(member))
          deallocate( g_tracer%tendency )
       endif
       g_tracer%tendency  => array 
    case ('field') 
       if (associated( g_tracer%field )) then
          call mpp_error(NOTE, trim(sub_name) // ": Deallocating generic tracer "//trim(name)//" % "//trim(member))
          deallocate( g_tracer%field )
       endif
       g_tracer%field3d_ptr  => array 
!       call set_cray_pointer_field(g_tracer%field,array,ilb,jlb)

    case default 
       call mpp_error(FATAL, trim(sub_name)//": Not a supported operation for member variable: "//trim(name)//" % "//trim(member))   
    end select

  end subroutine g_tracer_set_pointer_3D

  !The following does not compile:
  !error #6406: Conflicting attributes or multiple declaration of name.   [FIELD]
  !  pointer(ptr,field)
  !----------------^

!  subroutine set_cray_pointer_field(field,array,ilb,jlb)
!    real, dimension(:,:,:,:), intent(inout)     ::  field
!    integer,                        intent(in) :: ilb,jlb
!    real, dimension(ilb:,jlb:,:), target, intent(in) :: array
!
!    pointer(ptr,field)
!
!    ptr = LOC(array)
!
!  end subroutine set_cray_pointer_field


  ! <SUBROUTINE NAME="g_tracer_find">
  !  <OVERVIEW>
  !   Get the pointer for the named tracer node
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call g_tracer_find(g_tracer,name)
  !  </TEMPLATE>
  !  <IN/OUT NAME="g_tracer" TYPE="type(g_tracer_type),    pointer">
  !   Head of the generic tracer list.
  !   Upon return this will be a pointer to the tracer node called name or NULL if not found 
  !  </IN/OUT>
  !  <IN NAME="name" TYPE="character(len=*)">
  !   Name of a tracer node
  !  </IN>
  ! </SUBROUTINE>


  subroutine g_tracer_find(g_tracer,name)
    character(len=*),         intent(in) :: name
    type(g_tracer_type),    pointer    :: g_tracer

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_find'

    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    !Go through the list of tracers 
    do  
       if(g_tracer%name == name .or. g_tracer%alias == name) exit

       !traverse the linked list till hit NULL
       if(.NOT. associated(g_tracer%next)) then
          g_tracer => NULL() ! name Not found.
          exit
       endif

       g_tracer => g_tracer%next
    enddo
  end subroutine g_tracer_find
  

!#######################################################################
  ! <SUBROUTINE NAME="g_tracer_column_int">
  !  <OVERVIEW>
  !   Calculate the column interval for a given variable
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Calculate the column interval for a given variable
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call g_tracer_column_int(depth, ilb, jlb, var, dzt, rho_dzt, rd, k_level, integral, caller)
  !  </TEMPLATE>
  !  <IN NAME="depth" TYPE="real">
  !   Depth over which to integrate
  !  </IN>
  !  <IN NAME="ilb" TYPE="integer">
  !   Lower bound of 1st dimension of arrays
  !  </IN>
  !  <IN NAME="jlb" TYPE="integer">
  !   Lower bound of 2nd dimension of arrays
  !  </IN>
  !  <IN NAME="var" TYPE="real(:,:,:)">
  !   Variable to integrate
  !  </IN>
  !  <IN NAME="dzt" TYPE="real(:,:,:)">
  !   Layer thicknesses
  !  </IN>
  !  <IN NAME="rho_dzt" TYPE="real(:,:,:)">
  !   Density times layer thicknesses
  !  </IN>
  !  <INOUT NAME="rd" TYPE="real(:,:,:)">
  !   Work array: rho_dzt to be multiplied by var to do the integral (may be used in subsequent calls)
  !  </INOUT>
  !  <INOUT NAME="k_level" TYPE="integer">
  !   K level for maximum depth to perform the integral, if 0 then calculate rd array (may be used in subsequent calls)
  !   If set greater than 0, then the work array can be used in subsequent calls for the same depth to save some
  !   computation. Care should be taken that if k_level is set > 0 that the same depth range is used.
  !  </INOUT>
  !  <OUT NAME="integral" TYPE="real(:,:,:)">
  !   Integral of var over depth
  !  </OUT>
  !  <IN NAME="caller" TYPE="character(len=*), optional">
  !   string indicating caller of this routine, for traceback
  !  </IN>
  ! </SUBROUTINE>

  subroutine g_tracer_column_int(depth, ilb, jlb, var, dzt, rho_dzt, rd, k_level, integral, caller)

    real,                         intent(in)            :: depth
    integer,                      intent(in)            :: ilb
    integer,                      intent(in)            :: jlb
    real, dimension(ilb:,jlb:,:), intent(in)            :: var
    real, dimension(ilb:,jlb:,:), intent(in)            :: dzt
    real, dimension(ilb:,jlb:,:), intent(in)            :: rho_dzt
    real, dimension(ilb:,jlb:,:), intent(inout)         :: rd
    integer,                      intent(inout)         :: k_level
    real, dimension(ilb:,jlb:),   intent(out)           :: integral
    character(len=*),             intent(in), optional  :: caller

!-----------------------------------------------------------------------
!     local parameters

    character(len=fm_string_len), parameter     :: sub_name = 'g_tracer_column_int'

    character(len=256)                          :: caller_str
    character(len=256)                          :: error_header
    character(len=256)                          :: warn_header
    character(len=256)                          :: note_header
    integer                                     :: isc
    integer                                     :: iec
    integer                                     :: jsc
    integer                                     :: jec
    integer                                     :: isd
    integer                                     :: ied
    integer                                     :: jsd
    integer                                     :: jed
    integer                                     :: nk
    integer                                     :: ntau
    real,    dimension(:,:,:), pointer          :: grid_tmask
    integer                                     :: i
    integer                                     :: j
    integer                                     :: k 
    logical                                     :: continue_calc
    real, dimension(:,:), allocatable           :: depth_x
    
    !  Set up the headers for stdout messages.

    if (present(caller)) then
      caller_str = trim(mod_name) // '(' // trim(sub_name) // ')[' // trim(caller) // ']'
    else
      caller_str = trim(mod_name) // '(' // trim(sub_name) // ')[]'
    endif
    error_header = '==> Error from '   // trim(caller_str) // ':'
    warn_header =  '==> Warning from ' // trim(caller_str) // ':'
    note_header =  '==> Note from '    // trim(caller_str) // ':'

    !
    ! Check the depth
    !

    if (depth .le. 0.0) then
      call mpp_error(FATAL, trim(error_header) // ' Depth <= 0,0')
    endif

    !  Set up the module if not already done

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau,  &
         grid_tmask = grid_tmask)

    !
    ! Check the k_level
    !

    if (k_level .gt. nk) then
      call mpp_error(FATAL, trim(error_header) // ' k_level > nk')
    endif

    !
    !   Calculate the integral
    !

    if (k_level .le. 0) then  !{
      allocate (depth_x(isd:ied,jsd:jed))
      depth_x(:,:) = depth
      rd(:,:,:) = 0.0
      do k = 1, nk  !}
        continue_calc = .false.
        do j = jsc, jec  !{
          do i = isc, iec  !{
            if (grid_tmask(i,j,k) .gt. 0.5 .and. depth_x(i,j) .gt. 0.0) then  !{
              k_level = k
              if (depth_x(i,j) .gt. dzt(i,j,k)) then  !{
                continue_calc = .true.
                rd(i,j,k) = rho_dzt(i,j,k)
                depth_x(i,j) = depth_x(i,j) - dzt(i,j,k)
              else  !}{
                rd(i,j,k) = depth_x(i,j) / dzt(i,j,k) * rho_dzt(i,j,k)
                depth_x(i,j) = 0.0
              endif  !}
            endif  !}
          enddo  !} i
        enddo  !} j
        if (.not. continue_calc) then
          exit
        endif
      enddo  !} k
      deallocate (depth_x)
    endif  !}

    integral(:,:) = 0.0
    do k = 1, k_level  !}
      do j = jsc, jec  !{
        do i = isc, iec  !{
          integral(i,j) = integral(i,j) + var(i,j,k) * rd(i,j,k)
        enddo  !} i
      enddo  !} j
    enddo  !} k
            
    return

  end subroutine g_tracer_column_int
  

!#######################################################################
  ! <SUBROUTINE NAME="g_tracer_flux_at_depth">
  !  <OVERVIEW>
  !   Calculate the column interval for a given variable
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Calculate the column interval for a given variable
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call g_tracer_flux_at_depth(depth, ilb, jlb, var, dzt, k_level, frac, initialized, flux, caller)
  !  </TEMPLATE>
  !  <IN NAME="depth" TYPE="real">
  !   Depth over which to integrate
  !  </IN>
  !  <IN NAME="ilb" TYPE="integer">
  !   Lower bound of 1st dimension of arrays
  !  </IN>
  !  <IN NAME="jlb" TYPE="integer">
  !   Lower bound of 2nd dimension of arrays
  !  </IN>
  !  <IN NAME="var" TYPE="real(:,:,:)">
  !   Variable to integrate
  !  </IN>
  !  <IN NAME="dzt" TYPE="real(:,:,:)">
  !   Layer thicknesses
  !  </IN>
  !  <INOUT NAME="k_level" TYPE="real(:,:)">
  !   Work array: array of k level for each grid point at which depth occurs (may be used in future calls)
  !  </INOUT>
  !  <INOUT NAME="frac" TYPE="real(:,:)">
  !   Work array: fraction of level at which depth occurs (may be used in future calls)
  !  </INOUT>
  !  <INOUT NAME="initialized" TYPE="logical">
  !   True if the arrays have been initialized from a previous call, set to true in subroutine.
  !   If true, then the work arrays can be used in subsequent calls for the same depth to save some
  !   computation. Care should be taken that if iniitialized is set to true that the same depth range is used.
  !  </INOUT>
  !  <OUT NAME="flux" TYPE="real(:,:)">
  !   Flux at specified depth
  !  </OUT>
  !  <IN NAME="caller" TYPE="character(len=*), optional">
  !   string indicating caller of this routine, for traceback
  !  </IN>
  ! </SUBROUTINE>

  subroutine g_tracer_flux_at_depth(depth, ilb, jlb, var, dzt, k_level, frac, initialized, flux, caller)

    real,                            intent(in)                 :: depth
    integer,                         intent(in)                 :: ilb
    integer,                         intent(in)                 :: jlb
    real,    dimension(ilb:,jlb:,:), intent(in)                 :: var
    real,    dimension(ilb:,jlb:,:), intent(in)                 :: dzt
    integer, dimension(ilb:,jlb:),   intent(inout)              :: k_level
    real,    dimension(ilb:,jlb:),   intent(inout)              :: frac
    logical,                         intent(inout)              :: initialized
    real,    dimension(ilb:,jlb:),   intent(out)                :: flux
    character(len=*),                intent(in),    optional    :: caller

!-----------------------------------------------------------------------
!     local parameters

    character(len=fm_string_len), parameter     :: sub_name = 'g_tracer_flux_at_depth'

    character(len=256)                          :: caller_str
    character(len=256)                          :: error_header
    character(len=256)                          :: warn_header
    character(len=256)                          :: note_header
    integer                                     :: isc
    integer                                     :: iec
    integer                                     :: jsc
    integer                                     :: jec
    integer                                     :: isd
    integer                                     :: ied
    integer                                     :: jsd
    integer                                     :: jed
    integer                                     :: nk
    integer                                     :: ntau
    real,    dimension(:,:,:), pointer          :: grid_tmask
    integer                                     :: i
    integer                                     :: j
    integer                                     :: k 
    real, dimension(:,:), allocatable           :: depth_x
    logical                                     :: continue_calc
    
    !  Set up the headers for stdout messages.

    if (present(caller)) then
      caller_str = trim(mod_name) // '(' // trim(sub_name) // ')[' // trim(caller) // ']'
    else
      caller_str = trim(mod_name) // '(' // trim(sub_name) // ')[]'
    endif
    error_header = '==> Error from '   // trim(caller_str) // ':'
    warn_header =  '==> Warning from ' // trim(caller_str) // ':'
    note_header =  '==> Note from '    // trim(caller_str) // ':'

    !
    ! Check the depth
    !

    if (depth .le. 0.0) then
      call mpp_error(FATAL, trim(error_header) // ' Depth <= 0,0')
    endif

    !  Set up the module if not already done

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau,  &
         grid_tmask = grid_tmask)

    !
    !   Calculate the flux
    !

    if (.not. initialized) then  !{
      allocate (depth_x(isd:ied,jsd:jed))
      depth_x(:,:) = depth
      frac(:,:) = 0.0
      k_level(:,:) = 0
      do k = 1, nk  !{
        continue_calc = .false.
        do j = jsc, jec  !{
          do i = isc, iec  !{
            if (grid_tmask(i,j,k) .gt. 0.5 .and. depth_x(i,j) .gt. 0.0) then  !{
              if (depth_x(i,j) .gt. dzt(i,j,k)) then  !{
                continue_calc = .true.
                depth_x(i,j) = depth_x(i,j) - dzt(i,j,k)
              else  !}{
                frac(i,j) = depth_x(i,j) / dzt(i,j,k)
                k_level(i,j) = k
                depth_x(i,j) = 0.0
              endif  !}
            endif  !}
          enddo  !} i
        enddo  !} j
        if (.not. continue_calc) then
          exit
        endif
      enddo  !} k
      deallocate (depth_x)
    endif  !}
    initialized = .true.

    flux(:,:) = 0.0
    do j = jsc, jec  !{
      do i = isc, iec  !{
        if (k_level(i,j) .gt. 0) then  !{
          k = k_level(i,j)
          if (k .eq. 1) then
            flux(i,j) = frac(i,j) * var(i,j,k)
          else
            flux(i,j) = (1.0 - frac(i,j)) * var(i,j,k-1) + frac(i,j) * var(i,j,k)
          endif
        endif  !}
      enddo  !} i
    enddo  !} j
            
    return

  end subroutine g_tracer_flux_at_depth


  ! <SUBROUTINE NAME="g_tracer_send_diag">
  !  <OVERVIEW>
  !   Send diagnostics for all registered fields (if in diag_table)
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Collectively sends out the diagnostics for all registered fields of all generic tracers
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call g_tracer_send_diag(g_tracer_list,model_time , tau)
  !  </TEMPLATE>
  !  <IN NAME="g_tracer_list" TYPE="type(g_tracer_type),    pointer">
  !   pointer to the head of the generic tracer list
  !  </IN>
  !  <IN NAME="model_time" TYPE="type(time_type)">
  !   Time that the diagnostics is sent
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   The time step for the %field 4D field to be reported
  !  </IN>
  ! </SUBROUTINE>

  subroutine g_tracer_send_diag(g_tracer_list,model_time,tau)
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    type(time_type),          intent(in) :: model_time
    integer,                  intent(in) :: tau
    integer :: tau_1
    logical :: used

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_send_diag'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!

    !Go through the list of tracers 
    do  
       tau_1=tau
       if (g_tracer%diag_id_field .gt. 0) then
          if(.NOT. g_tracer_is_prog(g_tracer)) tau_1=1

       if(associated(g_tracer%field)) then 
          used = send_data(g_tracer%diag_id_field, g_tracer%field(:,:,:,tau_1), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,:),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc, ks_in=1,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec, ke_in=g_tracer_com%nk)
       elseif(associated(g_tracer%field3d_ptr)) then 
          used = send_data(g_tracer%diag_id_field, g_tracer%field3d_ptr(:,:,:), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,:),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc, ks_in=1,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec, ke_in=g_tracer_com%nk)

       else
          call mpp_error(FATAL, trim(sub_name)//": Cannot send_diag field variable for "//trim(g_tracer%name) )
       endif      


       endif

       if (g_tracer%diag_id_vmove .gt. 0 .and. _ALLOCATED(g_tracer%vmove)) then
          used = send_data(g_tracer%diag_id_vmove, g_tracer%vmove(:,:,:), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,:),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc, ks_in=1,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec, ke_in=g_tracer_com%nk)
       endif

       if (g_tracer%diag_id_vdiff .gt. 0 .and. _ALLOCATED(g_tracer%vdiff)) then
          used = send_data(g_tracer%diag_id_vdiff, g_tracer%vdiff(:,:,:), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,:),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc, ks_in=1,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec, ke_in=g_tracer_com%nk)
       endif

       if (g_tracer%diag_id_aux .gt. 0) then
          used = send_data(g_tracer%diag_id_aux, g_tracer%tendency(:,:,:), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,:),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc, ks_in=1,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec, ke_in=g_tracer_com%nk)
       endif

       if (g_tracer%diag_id_stf .gt. 0 .and. _ALLOCATED(g_tracer%stf)) then
          used = send_data(g_tracer%diag_id_stf, g_tracer%stf(:,:), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,1),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec )
       endif

       if (g_tracer%diag_id_stf_gas .gt. 0 .and. _ALLOCATED(g_tracer%stf_gas)) then
          used = send_data(g_tracer%diag_id_stf_gas, g_tracer%stf_gas(:,:), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,1),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec )
       endif

       if (g_tracer%diag_id_deltap .gt. 0 .and. _ALLOCATED(g_tracer%deltap)) then
          used = send_data(g_tracer%diag_id_deltap, g_tracer%deltap(:,:), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,1),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec )
       endif

       if (g_tracer%diag_id_kw .gt. 0 .and. _ALLOCATED(g_tracer%kw)) then
          used = send_data(g_tracer%diag_id_kw, g_tracer%kw(:,:), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,1),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec )
       endif

       if (g_tracer%diag_id_btf .gt. 0 .and. _ALLOCATED(g_tracer%btf)) then
          used = send_data(g_tracer%diag_id_btf, g_tracer%btf(:,:), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,1),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec )
       endif

       if (g_tracer%diag_id_btm .gt. 0 .and. _ALLOCATED(g_tracer%btm_reservoir)) then
          used = send_data(g_tracer%diag_id_btm, g_tracer%btm_reservoir(:,:), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,1),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec )
       endif

       if (g_tracer%diag_id_trunoff .gt. 0 .and. _ALLOCATED(g_tracer%trunoff)) then
          used = send_data(g_tracer%diag_id_trunoff, g_tracer%trunoff(:,:), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,1),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec )
       endif

       if (g_tracer%diag_id_alpha .gt. 0 .and. _ALLOCATED(g_tracer%alpha)) then
          used = send_data(g_tracer%diag_id_alpha, g_tracer%alpha(:,:), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,1),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec )
       endif

       if (g_tracer%diag_id_csurf .gt. 0 .and. _ALLOCATED(g_tracer%csurf)) then
          used = send_data(g_tracer%diag_id_csurf, g_tracer%csurf(:,:), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,1),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec )
       endif

       if (g_tracer%diag_id_sc_no .gt. 0 .and. _ALLOCATED(g_tracer%sc_no)) then
          used = send_data(g_tracer%diag_id_sc_no, g_tracer%sc_no(:,:), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,1),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec )
       endif

       !traverse the linked list till hit NULL
       if(.NOT. associated(g_tracer%next)) exit
       g_tracer => g_tracer%next
    enddo

  end subroutine g_tracer_send_diag


  ! <SUBROUTINE NAME="g_tracer_diag">
  !  <OVERVIEW>
  !   Send diagnostics for all registered fields at finish (if in diag_table)
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Collectively sends out the diagnostics for all registered fields of all generic tracers
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call g_tracer_diag(g_tracer_list,model_time , tau)
  !  </TEMPLATE>
  !  <IN NAME="g_tracer_list" TYPE="type(g_tracer_type),    pointer">
  !   pointer to the head of the generic tracer list
  !  </IN>
  !  <IN NAME="model_time" TYPE="type(time_type)">
  !   Time that the diagnostics is sent
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   The time step for the %field 4D field to be reported
  !  </IN>
  ! </SUBROUTINE>

  subroutine g_tracer_diag(g_tracer_list, ilb, jlb, rho_dzt_tau, rho_dzt_taup1, model_time, tau, taup1, dtts)
    type(g_tracer_type),    pointer    :: g_tracer_list
    integer,                  intent(in) :: ilb
    integer,                  intent(in) :: jlb
    real, dimension(ilb:,jlb:,:),   intent(in) :: rho_dzt_tau
    real, dimension(ilb:,jlb:,:),   intent(in) :: rho_dzt_taup1
    type(time_type),          intent(in) :: model_time
    integer,                  intent(in) :: tau
    integer,                  intent(in) :: taup1
    real,                     intent(in) :: dtts

    type(g_tracer_type),    pointer    :: g_tracer
    integer :: tau_1
    logical :: used

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_diag'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!

    !Go through the list of tracers 
    do  
       tau_1=taup1
       if (g_tracer%diag_id_field_taup1 .gt. 0) then
          if(.NOT. g_tracer_is_prog(g_tracer)) tau_1=1
          used = send_data(g_tracer%diag_id_field_taup1, g_tracer%field(:,:,:,tau_1), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,:),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc, ks_in=1,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec, ke_in=g_tracer_com%nk)
       endif

       if (g_tracer%diag_id_tendency .gt. 0  .and. g_tracer%prog) then
          used = send_data(g_tracer%diag_id_tendency,&
               (g_tracer%field(:,:,:,taup1)*rho_dzt_taup1(:,:,:) - g_tracer%field(:,:,:,tau)*rho_dzt_tau(:,:,:))/dtts, model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,:),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc, ks_in=1,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec, ke_in=g_tracer_com%nk)
       endif

       if (g_tracer%diag_id_vdiffuse_impl .gt. 0 .and. _ALLOCATED(g_tracer%vdiffuse_impl)) then
          used = send_data(g_tracer%diag_id_vdiffuse_impl, g_tracer%vdiffuse_impl(:,:,:), model_time,&
               rmask = g_tracer_com%grid_tmask(:,:,:),& 
               is_in=g_tracer_com%isc, js_in=g_tracer_com%jsc, ks_in=1,&
               ie_in=g_tracer_com%iec, je_in=g_tracer_com%jec, ke_in=g_tracer_com%nk)
       endif

       !traverse the linked list till hit NULL
       if(.NOT. associated(g_tracer%next)) exit
       g_tracer => g_tracer%next
    enddo

  end subroutine g_tracer_diag


  subroutine g_tracer_traverse(g_tracer_list)
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer 

    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_traverse'

    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer => g_tracer_list !Local pointer. Do not change the input pointer!

    !Go through the list of tracers 
    do  
       !do nothing
       !traverse the linked list till hit NULL
       if(.NOT. associated(g_tracer%next)) exit
       g_tracer => g_tracer%next
    enddo

  end subroutine g_tracer_traverse

  !
  !The following subroutines work with individual tracer nodes
  !
  ! <SUBROUTINE NAME="g_tracer_get_name">
  !  <OVERVIEW>
  !   Get the name of a particular tracer Node
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call g_tracer_get_name(g_tracer,string)
  !  </TEMPLATE>
  !  <IN NAME="g_tracer" TYPE="type(g_tracer_type),    pointer">
  !   Pointer to tracer node
  !  </IN>
  !  <IN NAME="string" TYPE="character(len=*)">
  !   Name of the tracer upon return
  !  </IN>
  ! </SUBROUTINE>

  subroutine g_tracer_get_name(g_tracer,string)
    type(g_tracer_type),    pointer    :: g_tracer 
    character(len=*),        intent(out) :: string
    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_get_name'

    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    string=g_tracer%name
  end subroutine g_tracer_get_name

  subroutine g_tracer_get_alias(g_tracer,string)
    type(g_tracer_type),    pointer    :: g_tracer 
    character(len=*),        intent(out) :: string
    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_get_alias'

    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    string=g_tracer%alias
  end subroutine g_tracer_get_alias

  ! <SUBROUTINE NAME="g_tracer_is_prog">
  !  <OVERVIEW>
  !   Is the tracer prognostic?
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   flag = g_tracer_is_prog(g_tracer)  
  !  </TEMPLATE>
  !  <IN NAME="g_tracer" TYPE="type(g_tracer_type),    pointer">
  !   Pointer to tracer node
  !  </IN>
  !  RETURNS .true. for prognostic tracer, .false. for diagnostic
  ! </SUBROUTINE>
  function g_tracer_is_prog(g_tracer) 
    logical :: g_tracer_is_prog
    type(g_tracer_type),    pointer    :: g_tracer 
    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_is_prog'

    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer_is_prog=g_tracer%prog
  end function g_tracer_is_prog

  ! <SUBROUTINE NAME="g_tracer_get_next">
  !  <OVERVIEW>
  !   get the next tracer in the list
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call g_tracer_get_next(g_tracer,g_tracer_next)
  !  </TEMPLATE>
  !  <IN NAME="g_tracer" TYPE="type(g_tracer_type),    pointer">
  !   Pointer to tracer node
  !  </IN>
  !  <IN NAME="g_tracer_next" TYPE="type(g_tracer_type),    pointer">
  !   Pointer to the next tracer node in the list
  !  </IN>
  ! </SUBROUTINE>
  subroutine g_tracer_get_next(g_tracer,g_tracer_next)
    type(g_tracer_type),    pointer    :: g_tracer,g_tracer_next 
    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_get_next'

    if(.NOT. associated(g_tracer)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    g_tracer_next => g_tracer%next
  end subroutine g_tracer_get_next

  ! <SUBROUTINE NAME="g_tracer_vertdiff_G">
  !  <OVERVIEW>
  !   Vertical Diffusion of a tracer node
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine solves a tridiagonal equation to find and set values of vertically diffused field for a tracer node.
  !   This is ported from GOLD (vertdiff) and simplified
  !   Since the surface flux from the atmosphere (%stf) has the units of mol/m^2/sec the resulting tracer concentration
  !   has units of mol/Kg
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call 
  !  </TEMPLATE>
  !  <IN NAME="" TYPE="">
  !   
  !  </IN>
  ! </SUBROUTINE>

  subroutine g_tracer_vertdiff_G(g_tracer, h_old, ea, eb, dt, kg_m2_to_H, m_to_H, tau, mom)
    type(g_tracer_type),    pointer  :: g_tracer
    real, dimension(g_tracer_com%isd:,g_tracer_com%jsd:,:), intent(in) :: h_old, ea, eb
    real,                   intent(in) :: dt, kg_m2_to_H, m_to_H
    integer,                intent(in) :: tau
    logical,                                                intent(in), optional :: mom

    ! Arguments: h_old -  Layer thickness before entrainment, in m or kg m-2.
    !                     In all the following comments the units of h_old are
    !                     denoted as H.
    !  (in)      ea - The amount of fluid entrained from the layer above, in H.
    !  (in)      eb - The amount of fluid entrained from the layer below, in H.
    !  (in)      dt - The amount of time covered by this call, in s.
    !  (in)      kg_m2_to_H - A conversion factor that translates kg m-2 into
    !                         the units of h_old (H).
    !  (in)      m_to_H - A conversion factor that translates m into the units
    !                     of h_old (H).
    !  (in,opt)  mom - If true, then called from MOM and don't do diagnostic,
    !                  if false or not present, then not from MOM and do diagnostics.

    !   This subroutine solves a tridiagonal equation for the final tracer
    ! concentrations after the dual-entrainments, and possibly sinking or surface
    ! and bottom sources, are applied.  The sinking is implemented with an
    ! fully implicit upwind advection scheme.

    real :: sink_dist(1:g_tracer_com%nk+1)    ! The distance the tracer sinks in a time step, in H.
    real :: sfc_src      ! The time-integrated surface source of the tracer, in
    ! units of H times a concentration.
    real :: btm_src      ! The time-integrated bottom source of the tracer, in
    ! units of H times a concentration.
    real :: b1           ! b1 is used by the tridiagonal solver, in H-1.
    real :: d1           ! d1=1-c1 is used by the tridiagonal solver, nondimensional.
    real :: c1(1:g_tracer_com%nk)     ! c1 is used by the tridiagonal solver, ND.
    real :: h_minus_dsink(1:g_tracer_com%nk)  ! The layer thickness minus the
    ! difference in sinking rates across the layer, in H.
    ! By construction, 0 <= h_minus_dsink < h_old.
    real :: sink(1:g_tracer_com%nk+1) ! The tracer's sinking distances at the
    ! interfaces, limited to prevent characteristics from
    ! crossing within a single timestep, in H.
    real :: b_denom_1    ! The first term in the denominator of b1, in H.
    real :: H_to_kg_m2   ! 1 / kg_m2_to_H.
    integer :: i, j, k, nz
    logical :: do_diagnostic

    !
    !   Save the current state for calculation of the implicit vertical diffusion term
    !

    if (g_tracer%diag_id_vdiffuse_impl .gt. 0) then
      if (present(mom)) then
        do_diagnostic = .not. mom
      else
        do_diagnostic = .false.
      endif
    else
      do_diagnostic = .false.
    endif
    if (do_diagnostic) then
      do j = g_tracer_com%jsc, g_tracer_com%jec
         do i = g_tracer_com%isc, g_tracer_com%iec
            do k = 1, g_tracer_com%nk
               g_tracer%vdiffuse_impl(i,j,k) = g_tracer%field(i,j,k,tau) 
            enddo
         enddo
      enddo
    endif

    d1 = 0.0
    H_to_kg_m2 = 1.0 / kg_m2_to_H

    sink_dist = (dt*g_tracer%sink_rate) * m_to_H

    do j=g_tracer_com%jsc,g_tracer_com%jec ; do i=g_tracer_com%isc,g_tracer_com%iec 

       if (g_tracer_com%grid_tmask(i,j,1) > 0.5) then

          nz=g_tracer_com%grid_kmt(i,j)

          if (g_tracer%move_vertical) then
	    do k=2,nz; sink_dist(k) = (dt*g_tracer%vmove(i,j,k)) * m_to_H; enddo
	  endif
          sfc_src = 0.0 ; btm_src = 0.0 

          ! Find the sinking rates at all interfaces, limiting them if necesary
          ! so that the characteristics do not cross within a timestep.
          !   If a non-constant sinking rate were used, that would be incorprated
          ! here.
          if (_ALLOCATED(g_tracer%btm_reservoir)) then
             do k=2,nz 
                sink(k) = sink_dist(k) ; h_minus_dsink(k) = h_old(i,j,k)
             enddo
             sink(nz+1) = sink_dist(nz+1) 
          else
             sink(nz+1) = 0.0 
             ! Find the limited sinking distance at the interfaces.
             do k=nz,2,-1
                if (sink(k+1) >= sink_dist(k)) then
                   sink(k) = sink_dist(k)
                   h_minus_dsink(k) = h_old(i,j,k) + (sink(k+1) - sink(k))
                elseif (sink(k+1) + h_old(i,j,k) < sink_dist(k)) then
                   sink(k) = sink(k+1) + h_old(i,j,k)
                   h_minus_dsink(k) = 0.0
                else
                   sink(k) = sink_dist(k)
                   h_minus_dsink(k) = (h_old(i,j,k) + sink(k+1)) - sink(k)
                endif
             enddo
          endif

          sink(1) = 0.0 ; h_minus_dsink(1) = (h_old(i,j,1) + sink(2))

          !Avoid sinking tracers with negative concentrations
          do k=2,nz+1
             if(g_tracer%field(i,j,k-1,tau) <= 0.0) sink(k) = 0.0
          enddo

          ! Now solve the tridiagonal equation for the tracer concentrations.

          b_denom_1 = h_minus_dsink(1) + ea(i,j,1)
          b1 = 1.0 / (b_denom_1 + eb(i,j,1))
          d1 = b_denom_1 * b1

          if (_ALLOCATED(g_tracer%stf)) sfc_src = (g_tracer%stf(i,j)*dt)*kg_m2_to_H

          g_tracer%field(i,j,1,tau) = b1*(h_old(i,j,1)*g_tracer%field(i,j,1,tau) + sfc_src)

          do k=2,nz-1 
             c1(k) = eb(i,j,k-1) * b1
             b_denom_1 = h_minus_dsink(k) + d1 * (ea(i,j,k) + sink(k))
             b1 = 1.0 / (b_denom_1 + eb(i,j,k))
             d1 = b_denom_1 * b1

             g_tracer%field(i,j,k,tau) = b1 * (h_old(i,j,k) * g_tracer%field(i,j,k,tau) + &
                  (ea(i,j,k) + sink(k)) * g_tracer%field(i,j,k-1,tau))
          enddo


          c1(nz) = eb(i,j,nz-1) * b1
          b_denom_1 = h_minus_dsink(nz) + d1 * (ea(i,j,nz) + sink(nz))
          b1 = 1.0 / (b_denom_1 + eb(i,j,nz))

          if (_ALLOCATED(g_tracer%btf)) btm_src = (-g_tracer%btf(i,j)*dt)*kg_m2_to_H

          g_tracer%field(i,j,nz,tau) = b1 * ((h_old(i,j,nz) * g_tracer%field(i,j,nz,tau) + btm_src) + &
               (ea(i,j,nz) + sink(nz)) * g_tracer%field(i,j,nz-1,tau))

          if (_ALLOCATED(g_tracer%btm_reservoir)) then 
             g_tracer%btm_reservoir(i,j) = g_tracer%btm_reservoir(i,j) + &
                 (sink(nz+1)*g_tracer%field(i,j,nz,tau))*H_to_kg_m2
          endif

          do k=nz-1,1,-1
             g_tracer%field(i,j,k,tau) = g_tracer%field(i,j,k,tau) + c1(k+1)*g_tracer%field(i,j,k+1,tau)
          enddo

        endif !(g_tracer_com%grid_tmask(i,j,1) > 0.5)

    enddo; enddo ! i,j

    !
    !   Calculate the implicit vertical diffusion term
    !   (Note: not sure if this needs any unit conversion)
    !

    if (do_diagnostic) then
      do j = g_tracer_com%jsc, g_tracer_com%jec
         do i = g_tracer_com%isc, g_tracer_com%iec
            do k = 1, g_tracer_com%nk
               g_tracer%vdiffuse_impl(i,j,k) = g_tracer_com%grid_tmask(i,j,k) *   &
                    (g_tracer%field(i,j,k,tau) - g_tracer%vdiffuse_impl(i,j,k)) / dt
            enddo
         enddo
      enddo
    endif

  end subroutine g_tracer_vertdiff_G

  ! <SUBROUTINE NAME="g_tracer_vertdiff_M">
  !  <OVERVIEW>
  !   Vertical Diffusion of a tracer node
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine solves a tridiagonal equation to find and set values of vertically diffused field for a tracer node.
  !   This is designed to calculate entrainments for MOM tracers and then call the g_tracer_vertdiff_G() above.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call 
  !  </TEMPLATE>
  !  <IN NAME="" TYPE="">
  !   
  !  </IN>
  ! </SUBROUTINE>

  subroutine g_tracer_vertdiff_M(g_tracer,dh, dhw, diff_cbt, dt, rho0,tau)
    type(g_tracer_type),    pointer  :: g_tracer
    real, dimension(g_tracer_com%isd:,g_tracer_com%jsd:,:), intent(in) :: dh, diff_cbt
    real, dimension(g_tracer_com%isd:,g_tracer_com%jsd:,0:), intent(in) :: dhw
    real,                   intent(in) :: dt,rho0
    integer,                intent(in) :: tau

    real, dimension(g_tracer_com%isd:g_tracer_com%ied,0:g_tracer_com%nk) ::  a,b,c, e, f, a1, c1
    real, dimension(g_tracer_com%isd:g_tracer_com%ied,0:g_tracer_com%nk) ::  wposu, wnegu, wposl, wnegl
    real, dimension(g_tracer_com%isd:g_tracer_com%ied,g_tracer_com%jsd:g_tracer_com%jed,0:g_tracer_com%nk) ::  dcb
    real, dimension(g_tracer_com%isd:g_tracer_com%ied) :: bet

    real, dimension(:,:,:), allocatable    :: ea, eb
    integer :: i, j, k, km1, kp1
    real :: eps, factu, factl, wabsu, wabsl, fact1, fact2
    logical :: GOLDtridiag = .true.
    logical :: IOWtridiag  = .false.
    
    GOLDtridiag = .true.
    IOWtridiag  = .false.
    if (g_tracer%move_vertical) then
       GOLDtridiag = .false.
       IOWtridiag  = .true.
    endif
    
    eps = 1.e-30

    !
    !   Save the current state for calculation of the implicit vertical diffusion term
    !

    if (g_tracer%diag_id_vdiffuse_impl .gt. 0) then
      do j = g_tracer_com%jsc, g_tracer_com%jec
         do i = g_tracer_com%isc, g_tracer_com%iec
            do k = 1, g_tracer_com%nk
               g_tracer%vdiffuse_impl(i,j,k) = g_tracer%field(i,j,k,tau) 
            enddo
         enddo
      enddo
    endif

    !
    !Add the contribution of K33_implicit to the diffusivity
    !
    !g_tracer%tendency should contain T_prog(n)%K33_implicit or be zero
    !!nnz: This is not really "tendency". Change name to something else like aux_array!
    do j=g_tracer_com%jsc,g_tracer_com%jec
       do i=g_tracer_com%isc,g_tracer_com%iec
          do k=1,g_tracer_com%nk
             dcb(i,j,k) = rho0*diff_cbt(i,j,k) + g_tracer%tendency(i,j,k) 
          enddo
       enddo
    enddo
    if (g_tracer%diff_vertical) then
      do j=g_tracer_com%jsc,g_tracer_com%jec
        do i=g_tracer_com%isc,g_tracer_com%iec
           do k=1,g_tracer_com%nk
             dcb(i,j,k) = dcb(i,j,k) + rho0*g_tracer%vdiff(i,j,k) 
           enddo
        enddo
      enddo
    endif
    !
    !The following two alternatives for solving the tridiagonal equation produce exact same results.
    !The choice between them should come from performance testing. I have not done this yet.
    !

    if(GOLDtridiag) then

       !===== 1 ===================
       !Via GOLD's vertdiff routine
       !===========================
       !            
       !h_old(i,j,k) = dh(i,j,k) (in kg m-2)

       allocate(   ea(g_tracer_com%isd:g_tracer_com%ied,g_tracer_com%jsd:g_tracer_com%jed,1:g_tracer_com%nk))
       allocate(   eb(g_tracer_com%isd:g_tracer_com%ied,g_tracer_com%jsd:g_tracer_com%jed,1:g_tracer_com%nk))

       do j=g_tracer_com%jsc,g_tracer_com%jec 
          do i=g_tracer_com%isc,g_tracer_com%iec

             ea(i,j,1) = 0.0 
             do k=2,g_tracer_com%nk
                ea(i,j,k) = dt * g_tracer_com%grid_tmask(i,j,k)  *dcb(i,j,k-1) /dhw(i,j,k-1)
             enddo

             do k=1,g_tracer_com%nk-1
                eb(i,j,k) = dt * g_tracer_com%grid_tmask(i,j,k+1)*dcb(i,j,k)   /dhw(i,j,k)
             enddo
             eb(i,j,g_tracer_com%nk) = 0.0 
          enddo
       enddo
       !Note: dh, ea, and eb have units here of kg m-2.
       call g_tracer_vertdiff_G(g_tracer, dh, ea, eb, dt, 1.0, rho0, tau, mom = .true.)

       !Mask out the field over "land" (land under Ocean)
       g_tracer%field(:,:,:,tau) = g_tracer%field(:,:,:,tau) * g_tracer_com%grid_tmask(:,:,:)

       deallocate(ea, eb)

    elseif(IOWtridiag) then

       !
       !OR 
       !
       !===== 2 ===================
       !Via MOM's invtri routine with vertical movement added
       !===========================
       !
       !This is borrowed from MOM invtri
       !     call invtri (T_prog(n)%field(:,:,:,taup1), T_prog(n)%stf, &
       !                  T_prog(n)%btf, wrk2(:,:,:), dtime_t, Grd%kmt,&
       !                  Grd%tmask, Thickness%rho_dzt(:,:,:,taup1), Thickness%dzwt, aidif, nk) 
       !
       do j=g_tracer_com%jsc,g_tracer_com%jec
          do k=1,g_tracer_com%nk
             km1   = max(1,k-1)
             kp1   = min(k+1,g_tracer_com%nk)
             do i=g_tracer_com%isc,g_tracer_com%iec
                fact1  = dt/dh(i,j,k)
		fact2  = rho0*fact1*0.5
		factu  = fact1/dhw(i,j,km1)
                factl  = fact1/dhw(i,j,k)
		wabsu       = abs(g_tracer%vmove(i,j,km1))
		wposu(i,k)  = fact2*(g_tracer%vmove(i,j,km1) + wabsu)*g_tracer_com%grid_tmask(i,j,k)
		wnegu(i,k)  = fact2*(g_tracer%vmove(i,j,km1) - wabsu)*g_tracer_com%grid_tmask(i,j,k)
		wabsl       = abs(g_tracer%vmove(i,j,k))
		wposl(i,k)  = fact2*(g_tracer%vmove(i,j,k  ) + wabsl)*g_tracer_com%grid_tmask(i,j,kp1)
		wnegl(i,k)  = fact2*(g_tracer%vmove(i,j,k  ) - wabsl)*g_tracer_com%grid_tmask(i,j,kp1)
                a1(i,k) = dcb(i,j,km1)*factu*g_tracer_com%grid_tmask(i,j,k)  
                c1(i,k) = dcb(i,j,k)  *factl*g_tracer_com%grid_tmask(i,j,kp1)
                a(i,k) = -(a1(i,k) - wnegu(i,k))
                c(i,k) = -(c1(i,k) + wposl(i,k))
                f(i,k) = g_tracer%field(i,j,k,tau)*g_tracer_com%grid_tmask(i,j,k) 
                b(i,k) = 1.0 + a1(i,k) + c1(i,k) - wnegl(i,k) + wposu(i,k)
             enddo
          enddo

          do i=g_tracer_com%isc,g_tracer_com%iec
             a1(i,1)  = 0.0
	     wnegu(i,1) = 0.0; wposu(i,1) = 0.0
	     a(i,1)  = 0.0
             c1(i,g_tracer_com%nk) = 0.0
	     wposl(i,g_tracer_com%nk) = 0.0; wnegl(i,g_tracer_com%nk) = 0.0 
             c(i,g_tracer_com%nk) = 0.0
             b(i,1)  = 1.0 + a1(i,1) + c1(i,1) - wnegl(i,1) + wposu(i,1)
             b(i,g_tracer_com%nk) = 1.0 + a1(i,g_tracer_com%nk) + c1(i,g_tracer_com%nk) &
	                                - wnegl(i,g_tracer_com%nk) + wposu(i,g_tracer_com%nk)

             ! top and bottom b.c.
             if (_ALLOCATED(g_tracer%stf)) &
                  f(i,1) = g_tracer%field(i,j,1,tau) + g_tracer%stf(i,j)*dt*g_tracer_com%grid_tmask(i,j,1)/dh(i,j,1)
             if (_ALLOCATED(g_tracer%btf)) then
                k = max(2,g_tracer_com%grid_kmt(i,j))
                f(i,k) = g_tracer%field(i,j,k,tau) - g_tracer%btf(i,j)*dt*g_tracer_com%grid_tmask(i,j,k)/dh(i,j,k)
             endif
          enddo

          ! decomposition and forward substitution
          do i=g_tracer_com%isc,g_tracer_com%iec
             bet(i) = g_tracer_com%grid_tmask(i,j,1)/(b(i,1) + eps)
             g_tracer%field(i,j,1,tau) = f(i,1)*bet(i)
          enddo
          do k=2,g_tracer_com%nk
             do i=g_tracer_com%isc,g_tracer_com%iec
                e(i,k) = c(i,k-1)*bet(i)
                bet(i) = g_tracer_com%grid_tmask(i,j,k)/(b(i,k) - a(i,k)*e(i,k) + eps)
                g_tracer%field(i,j,k,tau) = (f(i,k) - a(i,k)*g_tracer%field(i,j,k-1,tau))*bet(i)
             enddo
          enddo

          ! back substitution
          do k=g_tracer_com%nk-1,1,-1
             do i=g_tracer_com%isc,g_tracer_com%iec
                g_tracer%field(i,j,k,tau) = g_tracer%field(i,j,k,tau) - e(i,k+1)*g_tracer%field(i,j,k+1,tau)
             enddo
          enddo
       enddo
    else

       !
       !OR 
       !
       !===== 2 ===================
       !Via MOM's invtri routine. Works only for sink_rate=0
       !===========================
       !
       !This is borrowed from MOM invtri
       !     call invtri (T_prog(n)%field(:,:,:,taup1), T_prog(n)%stf, &
       !                  T_prog(n)%btf, wrk2(:,:,:), dtime_t, Grd%kmt,&
       !                  Grd%tmask, Thickness%rho_dzt(:,:,:,taup1), Thickness%dzwt, aidif, nk) 
       !
       do j=g_tracer_com%jsc,g_tracer_com%jec
          do k=1,g_tracer_com%nk
             km1   = max(1,k-1)
             kp1   = min(k+1,g_tracer_com%nk)
             do i=g_tracer_com%isc,g_tracer_com%iec
                factu  = dt/(dhw(i,j,k-1)*dh(i,j,k))
                factl  = dt/(dhw(i,j,k)*dh(i,j,k))
                a(i,k) = -dcb(i,j,km1)*factu*g_tracer_com%grid_tmask(i,j,k)
                c(i,k) = -dcb(i,j,k)  *factl*g_tracer_com%grid_tmask(i,j,kp1)
                f(i,k) = g_tracer%field(i,j,k,tau)*g_tracer_com%grid_tmask(i,j,k) 
                b(i,k) = 1.0 - a(i,k) - c(i,k)
             enddo
          enddo

          do i=g_tracer_com%isc,g_tracer_com%iec
             a(i,1)  = 0.0
             c(i,g_tracer_com%nk) = 0.0
             b(i,1)  = 1.0 - a(i,1) - c(i,1)
             b(i,g_tracer_com%nk) = 1.0 - a(i,g_tracer_com%nk) - c(i,g_tracer_com%nk)

             ! top and bottom b.c.
             if (_ALLOCATED(g_tracer%stf)) &
                  f(i,1) = g_tracer%field(i,j,1,tau) + g_tracer%stf(i,j)*dt*g_tracer_com%grid_tmask(i,j,1)/dh(i,j,1)
             if (_ALLOCATED(g_tracer%btf)) then
                k = max(2,g_tracer_com%grid_kmt(i,j))
                f(i,k) = g_tracer%field(i,j,k,tau) - g_tracer%btf(i,j)*dt*g_tracer_com%grid_tmask(i,j,k)/dh(i,j,k)
             endif
          enddo

          ! decomposition and forward substitution
          do i=g_tracer_com%isc,g_tracer_com%iec
             bet(i) = g_tracer_com%grid_tmask(i,j,1)/(b(i,1) + eps)
             g_tracer%field(i,j,1,tau) = f(i,1)*bet(i)
          enddo
          do k=2,g_tracer_com%nk
             do i=g_tracer_com%isc,g_tracer_com%iec
                e(i,k) = c(i,k-1)*bet(i)
                bet(i) = g_tracer_com%grid_tmask(i,j,k)/(b(i,k) - a(i,k)*e(i,k) + eps)
                g_tracer%field(i,j,k,tau) = (f(i,k) - a(i,k)*g_tracer%field(i,j,k-1,tau))*bet(i)
             enddo
          enddo

          ! back substitution
          do k=g_tracer_com%nk-1,1,-1
             do i=g_tracer_com%isc,g_tracer_com%iec
                g_tracer%field(i,j,k,tau) = g_tracer%field(i,j,k,tau) - e(i,k+1)*g_tracer%field(i,j,k+1,tau)
             enddo
          enddo
       enddo

    endif

    !
    !   Calculate the implicit vertical diffusion term
    !   (Note: dh = rho_dzt(taup1)
    !

    if (g_tracer%diag_id_vdiffuse_impl .gt. 0) then
      do j = g_tracer_com%jsc, g_tracer_com%jec
         do i = g_tracer_com%isc, g_tracer_com%iec
            do k = 1, g_tracer_com%nk
               g_tracer%vdiffuse_impl(i,j,k) = dh(i,j,k) * g_tracer_com%grid_tmask(i,j,k) *   &
                    (g_tracer%field(i,j,k,tau) - g_tracer%vdiffuse_impl(i,j,k)) / dt
            enddo
         enddo
      enddo
    endif

    return

  end subroutine g_tracer_vertdiff_M


  subroutine g_diag_field_add(node_ptr, diag_id, package_name, name, axes, init_time, longname, units, &
                            missing_value, Z_diag, field_ptr, Zname, Zlongname, Zunits)
    type(g_diag_type), pointer :: node_ptr
    integer, intent(inout) :: diag_id
    CHARACTER(len=*), INTENT(in) :: package_name, name
    INTEGER, INTENT(in) :: axes(:)
    TYPE(time_type), INTENT(in) :: init_time
    CHARACTER(len=*), INTENT(in) :: longname, units
    REAL, OPTIONAL, INTENT(in) :: missing_value
    integer, optional, intent(in) :: Z_diag
    CHARACTER(len=*), optional, INTENT(in) :: Zname, Zlongname, Zunits
    real, optional, pointer :: field_ptr(:,:,:)
    

    type(g_diag_type), pointer :: g_diag => NULL()
    
    !diag register with the original name
    diag_id = register_diag_field(package_name, name, axes, init_time, longname,units, missing_value = missing_value)

    !===================================================================
    !Add this diagnostics to the list that is going to be used later (by GOLD) 
    !===================================================================
    allocate(g_diag)
    
    g_diag%diag_id = diag_id
    g_diag%name         = trim(name)
    g_diag%longname     = trim(longname)
    g_diag%units        = trim(units)
    g_diag%package_name = trim(package_name)
    g_diag%axes         = axes
    g_diag%init_time    = init_time
    g_diag%missing_value= missing_value
    !Is this a Z diag?
    if(present(Z_diag))    g_diag%Z_diag       = Z_diag
    if(present(field_ptr)) g_diag%field_ptr    => field_ptr
    if(present(Zname))     g_diag%name         = trim(Zname)
    if(present(Zlongname)) g_diag%longname     = trim(Zlongname)
    if(present(Zunits))    g_diag%units        = trim(Zunits)


    !===================================================================
    !Reversed Linked List implementation! Make this new node to be the head of the list.
    !===================================================================    

    g_diag%next => node_ptr 
    node_ptr => g_diag 
  end subroutine g_diag_field_add


end module g_tracer_utils
